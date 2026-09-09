#ifndef IGA_COLLECTIVE_PETSC_OPTIONS_HPP
#define IGA_COLLECTIVE_PETSC_OPTIONS_HPP

#include "CollectiveFailure.hpp"
#include <petscoptions.h>

#include <map>
#include <optional>
#include <set>
#include <string>

namespace iga {
namespace petsc_options_detail {

inline void Check(PetscErrorCode code, const char* operation)
{
	if (code) throw std::runtime_error(std::string(operation)+" failed");
}

struct Buffers {
	explicit Buffers(PetscOptions database) : options(database) {}
	Buffers(const Buffers&) = delete;
	Buffers& operator=(const Buffers&) = delete;
	~Buffers()
	{
		if (have_unused) PetscOptionsLeftRestore(options, &count, &names, &values);
		PetscFree(all);
	}
	PetscOptions options;
	char* all = nullptr;
	PetscInt count = 0;
	char** names = nullptr;
	char** values = nullptr;
	bool have_unused = false;
};

} // namespace petsc_options_detail

// Local inspection of visible option entries. Preserve unused-option tracking:
// LeftGet supplies unused entries; FindPair is called only for already-used
// entries. Do not clone via GetAll/InsertString: GetAll's display string is
// ambiguous when a value contains spaces or text resembling another option.
using PetscOptionEntries = std::map<std::string, std::optional<std::string>>;

inline PetscOptionEntries CapturePetscOptionEntries(PetscOptions options,
	const std::set<std::string>& excluded = {})
{
	using petsc_options_detail::Check;
	petsc_options_detail::Buffers buffers(options);
	Check(PetscOptionsGetAll(options, &buffers.all), "PetscOptionsGetAll");
	Check(PetscOptionsLeftGet(options, &buffers.count, &buffers.names, &buffers.values),
		"PetscOptionsLeftGet");
	buffers.have_unused = true;
	std::map<std::string, const char*> unused;
	for (PetscInt i = 0; i < buffers.count; ++i)
		unused.emplace("-"+std::string(buffers.names[i]), buffers.values[i]);
	const std::string_view display(buffers.all ? buffers.all : "");
	PetscOptionEntries entries;
	std::size_t offset = 0;
	while (offset < display.size()) {
		const auto end = display.find(' ', offset);
		if (end == std::string_view::npos || end == offset || display[offset] != '-')
			throw std::runtime_error("unsupported PETSc option display encoding");
		const std::string name(display.substr(offset, end-offset));
		const char* value = nullptr;
		const auto found = unused.find(name);
		if (found != unused.end()) value = found->second;
		else {
			PetscBool present = PETSC_FALSE;
			Check(PetscOptionsFindPair(options, nullptr, name.c_str(), &value, &present),
				"PetscOptionsFindPair");
			if (!present) throw std::runtime_error("PETSc option enumeration changed during capture");
		}
		offset = end+1;
		if (value) {
			const std::string_view exact(value);
			if (display.substr(offset, exact.size()) != exact
				|| offset+exact.size() >= display.size() || display[offset+exact.size()] != ' ')
				throw std::runtime_error("unsupported PETSc option value encoding");
			offset += exact.size()+1;
		}
		if (!excluded.count(name)
			&& !entries.emplace(name, value ? std::optional<std::string>(value) : std::nullopt).second)
			throw std::runtime_error("duplicate PETSc option during capture");
	}
	return entries;
}

inline std::string SerializePetscOptionEntries(const PetscOptionEntries& entries)
{
	// Length framing distinguishes an option-like substring inside one value
	// from a separate option, including embedded spaces, newlines and quotes.
	std::string result;
	for (const auto& entry : entries) {
		result += std::to_string(entry.first.size())+":"+entry.first;
		result += entry.second ? "V"+std::to_string(entry.second->size())+":"+*entry.second : "F";
	}
	return result;
}

inline std::string CapturePetscOptions(PetscOptions options,
	const std::set<std::string>& excluded = {})
{
	return SerializePetscOptionEntries(CapturePetscOptionEntries(options, excluded));
}

inline void RequireCollectivePetscOptions(MPI_Comm communicator, PetscOptions options = nullptr,
	const std::set<std::string>& excluded = {})
{
	std::string snapshot;
	CollectiveLocalStage(communicator, "PETSc option capture", [&] {
		snapshot = CapturePetscOptions(options, excluded);
	});
	RequireCollectiveSameText(communicator, "PETSc option agreement", snapshot);
}

} // namespace iga

#endif
