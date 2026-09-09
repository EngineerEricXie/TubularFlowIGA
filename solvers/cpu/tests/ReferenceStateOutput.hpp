#ifndef IGA_TEST_REFERENCE_STATE_OUTPUT_HPP
#define IGA_TEST_REFERENCE_STATE_OUTPUT_HPP

#include "Sha256.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace iga { namespace test {

// Optional, post-acceptance fixture output. This is a numerical reference,
// not a checkpoint: a completed manifest is published only after all fields.
class ReferenceStateOutput {
public:
	ReferenceStateOutput(std::filesystem::path directory, std::string fixture)
		: directory_(std::move(directory)), fixture_(std::move(fixture))
	{
		if (fixture_ != "fsi" && fixture_ != "immersed")
			throw std::invalid_argument("unknown reference fixture");
		if (directory_.empty()) throw std::invalid_argument("empty reference directory");
		if (!directory_.parent_path().empty()) std::filesystem::create_directories(directory_.parent_path());
		if (!std::filesystem::create_directory(directory_))
			throw std::runtime_error("reference output directory already exists");
	}

	template <class Ids, class ValueAt>
	void Add(const std::string& name, const std::string& units, const Ids& ids,
		std::size_t columns, ValueAt value_at)
	{
		if (finished_ || failed_) throw std::logic_error("reference output is closed or failed");
		failed_ = true;
		if (name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos
			|| units.empty() || units.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_/*^- ") != std::string::npos
			|| !names_.insert(name).second || ids.empty() || !columns)
			throw std::invalid_argument("invalid reference field metadata");
		std::ostringstream text;
		text.imbue(std::locale::classic());
		text << std::scientific << std::setprecision(17);
		std::uint64_t previous = 0;
		for (std::size_t row = 0; row < ids.size(); ++row) {
			if constexpr (std::is_signed<typename Ids::value_type>::value)
				if (ids[row] < 0) throw std::invalid_argument("negative reference ID");
			const auto id = static_cast<std::uint64_t>(ids[row]);
			if (row && id <= previous)
				throw std::invalid_argument("reference IDs must be nonnegative and strictly increasing");
			previous = id;
			text << id;
			for (std::size_t column = 0; column < columns; ++column) {
				const double value = value_at(row, column);
				if (!std::isfinite(value)) throw std::invalid_argument("nonfinite reference field");
				text << ' ' << value;
			}
			text << '\n';
		}
		const auto contents = text.str();
		Write(directory_/(name+".txt"), contents);
		Sha256 hash; hash.Append(contents.data(), contents.size());
		fields_.push_back({name, units, hash.Hex(), ids.size(), columns});
		failed_ = false;
	}

	void Finish(double time_s, std::uint64_t step)
	{
		if (finished_ || failed_ || fields_.empty() || !std::isfinite(time_s) || time_s < 0.0)
			throw std::logic_error("cannot publish incomplete reference output");
		failed_ = true;
		std::ostringstream text;
		text.imbue(std::locale::classic());
		text << std::setprecision(17) << "{\n  \"schema_version\": 1,\n  \"kind\": \"fixture_reference_state\",\n"
			<< "  \"fixture\": \"" << fixture_ << "\",\n  \"time_s\": " << time_s
			<< ",\n  \"step\": " << step << ",\n  \"native_gates_passed\": true,\n  \"fields\": [";
		for (std::size_t i = 0; i < fields_.size(); ++i) {
			const auto& field = fields_[i];
			text << (i ? "," : "") << "\n    {\"name\": \"" << field.name << "\", \"file\": \""
				<< field.name << ".txt\", \"units\": \"" << field.units << "\", \"rows\": " << field.rows
				<< ", \"columns\": " << field.columns << ", \"sha256\": \"" << field.sha256 << "\"}";
		}
		text << "\n  ]\n}\n";
		Write(directory_/"manifest.tmp", text.str());
		std::filesystem::rename(directory_/"manifest.tmp", directory_/"manifest.json");
		finished_ = true;
		failed_ = false;
	}

private:
	struct Field { std::string name, units, sha256; std::size_t rows, columns; };
	static void Write(const std::filesystem::path& path, const std::string& contents)
	{
		std::ofstream output(path, std::ios::binary);
		if (!output) throw std::runtime_error("cannot open reference output");
		output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
		output.close();
		if (!output) throw std::runtime_error("cannot complete reference output");
	}
	std::filesystem::path directory_;
	std::string fixture_;
	std::set<std::string> names_;
	std::vector<Field> fields_;
	bool failed_ = false, finished_ = false;
};

inline std::unique_ptr<ReferenceStateOutput> ReferenceOutputFromArguments(
	int argc, char** argv, const std::string& fixture, int ranks)
{
	std::string directory;
	bool found = false;
	for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--reference-output") {
		if (found || ++i == argc || std::string(argv[i]).empty() || std::string(argv[i]).front() == '-')
			throw std::invalid_argument("--reference-output requires one unique directory argument");
		found = true; directory = argv[i];
	}
	if (!found) return nullptr;
	if (ranks != 1) throw std::invalid_argument("fixture reference export currently requires one MPI rank");
	return std::unique_ptr<ReferenceStateOutput>(new ReferenceStateOutput(directory, fixture));
}

} } // namespace iga::test

#endif
