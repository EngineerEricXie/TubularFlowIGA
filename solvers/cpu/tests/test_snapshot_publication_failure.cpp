// Reuse the physical fixture and its schema/identity regression.
#define main SnapshotPublicationContractMain
#include "test_moving_immersed_flow_snapshot_publisher.cpp"
#undef main

#include <cerrno>
#include <cstring>
#include <iostream>
#include <locale>
#include <map>
#include <new>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <unistd.h>
#include <csignal>

namespace snapshot_read_fault {
const char* suffix = nullptr;
bool injected = false;
}

// Linux test executable only: fail a real filebuf read, leaving the process
// alive. Restrict the fault to files in this test's exclusive output tree.
extern "C" ssize_t read(int descriptor, void* buffer, size_t size)
{
	if (snapshot_read_fault::suffix && !snapshot_read_fault::injected) {
		char link[64], path[4096];
		std::snprintf(link, sizeof(link), "/proc/self/fd/%d", descriptor);
		const auto count = readlink(link, path, sizeof(path)-1);
		if (count >= 0) {
			path[count] = '\0';
			const auto length = std::strlen(snapshot_read_fault::suffix);
			if (std::strstr(path, "snapshot-publication-failure-")
				&& static_cast<std::size_t>(count) >= length
				&& std::strcmp(path+count-length, snapshot_read_fault::suffix) == 0) {
				snapshot_read_fault::injected = true;
				errno = EIO;
				return -1;
			}
		}
	}
	return syscall(SYS_read, descriptor, buffer, size);
}

namespace {
namespace fs = std::filesystem;
using Publisher = iga::MovingImmersedFlowSnapshotPublisher;

// Fault once at the first numeric insertion of the selected stream instance.
// iword is reset by each stream constructor, even when its address is reused.
class FormatterFault : public std::num_put<char> {
public:
	FormatterFault(std::size_t target, int mode) : target_(target), mode_(mode) {}
	mutable std::size_t streams = 0;
	mutable bool injected = false;
protected:
	void Visit(std::ios_base& output) const
	{
		static const int slot = std::ios_base::xalloc();
		if (output.iword(slot)) return;
		output.iword(slot) = 1;
		if (++streams != target_) return;
		injected = true;
		if (mode_ == 0) throw std::runtime_error("injected publication formatter error");
		if (mode_ == 1) throw std::bad_alloc();
		throw 7;
	}
#define SNAPSHOT_PUT(Type) \
	iter_type do_put(iter_type it, std::ios_base& output, char fill, Type value) const override \
	{ Visit(output); return std::num_put<char>::do_put(it, output, fill, value); }
	SNAPSHOT_PUT(bool)
	SNAPSHOT_PUT(long)
	SNAPSHOT_PUT(unsigned long)
	SNAPSHOT_PUT(long long)
	SNAPSHOT_PUT(unsigned long long)
	SNAPSHOT_PUT(double)
	SNAPSHOT_PUT(long double)
#undef SNAPSHOT_PUT
private:
	std::size_t target_;
	int mode_;
};

std::map<std::string, std::string> Files(const fs::path& root)
{
	std::map<std::string, std::string> result;
	for (const auto& item : fs::recursive_directory_iterator(root)) {
		Check(item.path().string().find(".tmp.") == std::string::npos, "temporary publication leaked");
		if (item.is_regular_file()) result.emplace(fs::relative(item.path(), root).string(), Read(item.path()));
	}
	return result;
}

void CheckArrays(const fs::path& path, std::size_t points)
{
	const auto xml = Read(path);
	for (const auto& field : std::vector<std::pair<std::string, int>>{
		{"velocity",3},{"pressure",1},{"vorticity",3},{"q_criterion",1},
		{"enstrophy_density",1},{"speed",1},{"integration_weight",1},
		{"reference_parametric",3},{"background_cell_id",1},{"quadrature_ordinal",1},
		{"cell_kind",1},{"stagnant",1},{"points",3},{"connectivity",1},{"offsets",1},{"types",1}}) {
		const auto name = xml.find("Name=\""+field.first+"\"");
		Check(name != std::string::npos, "VTU array missing");
		const auto begin = xml.find('>', name)+1, end = xml.find('<', begin);
		std::istringstream values(xml.substr(begin, end-begin));
		double value = 0.0;
		std::size_t count = 0;
		while (values >> value) { Check(std::isfinite(value), "nonfinite VTU array value"); ++count; }
		Check(values.eof() && count == points*field.second, "VTU array tuple count differs from NumberOfPoints");
	}
}
}

int main()
{
	const auto root = fs::temp_directory_path()/("snapshot-publication-failure-"+std::to_string(getpid()));
	try {
		Check(fs::create_directory(root), "test output directory must be new");
		const auto fixture = MakeFixture();
		const auto old = Snapshot(fixture, 3, true, -.25), next = Snapshot(fixture, 4, false, 0.);
		const auto stem = root/"field";
		const auto first = Publisher::Publish(old, stem, Identity());
		const auto old_files = Files(root);
		const auto published = Publisher::Publish(next, stem, Identity());
		CheckArrays(published.vtu_path, next.Points().size());
		const auto expected = Files(root);
		std::size_t faults = 0, retries = 0;
		for (int operation = 0; operation < 3; ++operation) {
			auto prepare = [&] {
				if (operation == 0) {
					fs::remove_all(published.epoch_directory);
					Publisher::RebuildCollection(stem);
				}
			};
			auto work = [&] {
				if (operation == 2) Publisher::RebuildCollection(stem);
				else Publisher::Publish(next, stem, Identity());
			};
			std::size_t total = 0;
			for (std::size_t ordinal = 0; ordinal <= total; ++ordinal) {
				for (int mode = 0; mode < (ordinal ? 3 : 1); ++mode) {
					prepare();
					const auto before = Files(root);
					auto* facet = new FormatterFault(ordinal, mode);
					const std::locale locale(std::locale(), facet);
					const auto previous = std::locale::global(locale);
					bool rejected = false;
					try { work(); } catch (...) { rejected = true; }
					std::locale::global(previous);
					if (!ordinal) {
						Check(!rejected, "healthy publication failed");
						total = facet->streams;
						Check(total > 0 && total < 200, "unexpected formatter coverage");
					} else {
						Check(facet->injected && rejected, "publication swallowed formatter failure");
						++faults;
						const auto after = Files(root);
						for (const auto& file : before)
							Check(after.at(file.first) == file.second, "failed publication changed accepted output");
					}
					Publisher::Publish(next, stem, Identity());
					Check(Files(root) == expected, "retry changed epoch or collection bytes");
					++retries;
				}
			}
			std::cout << "snapshot operation=" << operation << " formatter_streams=" << total << '\n';
		}
		for (int operation = 0; operation < 3; ++operation) {
			for (const auto* suffix : {"metrics.json", "fields.vtu"}) {
				if (operation == 0) { fs::remove_all(published.epoch_directory); Publisher::RebuildCollection(stem); }
				const auto before = Files(root);
				snapshot_read_fault::suffix = suffix;
				snapshot_read_fault::injected = false;
				bool rejected = false;
				try {
					if (operation == 2) Publisher::RebuildCollection(stem);
					else Publisher::Publish(next, stem, Identity());
				} catch (...) { rejected = true; }
				snapshot_read_fault::suffix = nullptr;
				Check(snapshot_read_fault::injected && rejected, "publication accepted failed file read");
				Check(Files(root) == before, "failed read altered publication");
				++faults;
				Publisher::Publish(next, stem, Identity());
				Check(Files(root) == expected, "read-failure retry changed output");
				++retries;
			}
		}
		for (const auto fault : {iga::MovingImmersedFlowSnapshotPublicationFault::AfterTemporaryVtuWrite,
			iga::MovingImmersedFlowSnapshotPublicationFault::BeforeDirectoryRename,
			iga::MovingImmersedFlowSnapshotPublicationFault::AfterDirectoryRename}) {
			fs::remove_all(published.epoch_directory); Publisher::RebuildCollection(stem);
			Publisher::SetTestFault(fault);
			bool rejected = false;
			try { Publisher::Publish(next, stem, Identity()); } catch (...) { rejected = true; }
			Check(rejected, "publication interruption missed");
			const auto after = Files(root);
			for (const auto& file : old_files) Check(after.at(file.first) == file.second, "interruption damaged previous epoch");
			Publisher::Publish(next, stem, Identity());
			Check(Files(root) == expected, "interruption retry changed output");
			++faults; ++retries;
		}
		std::cout << "snapshot faults=" << faults << " retries=" << retries << " arrays=16 passed\n";
		fs::remove_all(root);
		return 0;
	} catch (const std::exception& error) {
		snapshot_read_fault::suffix = nullptr;
		std::cerr << error.what() << " evidence=" << root << '\n';
		return 1;
	}
}
