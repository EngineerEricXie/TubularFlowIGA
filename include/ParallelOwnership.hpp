#ifndef IGA_PARALLEL_OWNERSHIP_HPP
#define IGA_PARALLEL_OWNERSHIP_HPP

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {

// Global node IDs and local owned-node offsets are deliberately separate.
struct OwnedNodeRange {
	std::uint64_t begin = 0, end = 0, global_nodes = 0;

	void Validate() const
	{
		if (begin > end || end > global_nodes)
			throw std::invalid_argument("invalid owned node range");
	}
	bool Owns(std::uint64_t global_node) const
	{
		return global_node >= begin && global_node < end;
	}
	std::uint64_t LocalNode(std::uint64_t global_node) const
	{
		Validate();
		if (!Owns(global_node)) throw std::out_of_range("node is not locally owned");
		return global_node-begin;
	}
	std::uint64_t GlobalNode(std::uint64_t local_node) const
	{
		Validate();
		if (local_node >= end-begin) throw std::out_of_range("local owned node offset");
		return begin+local_node;
	}
};

inline std::uint64_t CheckedFieldRowCount(std::uint64_t nodes,
	std::int64_t fields, std::uint64_t index_limit)
{
	if (fields <= 0) throw std::invalid_argument("fields must be positive");
	const auto count = static_cast<std::uint64_t>(fields);
	if (nodes > index_limit/count)
		throw std::overflow_error("field row count exceeds index range");
	return nodes*count;
}

inline std::uint64_t CheckedFieldRow(std::int64_t node, std::int64_t field,
	std::uint64_t global_nodes, std::int64_t fields, std::uint64_t index_limit)
{
	CheckedFieldRowCount(global_nodes, fields, index_limit);
	if (node < 0 || static_cast<std::uint64_t>(node) >= global_nodes
		|| field < 0 || field >= fields)
		throw std::out_of_range("global node or field index");
	return static_cast<std::uint64_t>(node)*static_cast<std::uint64_t>(fields)
		+static_cast<std::uint64_t>(field);
}

// A required-state buffer may contain owned nodes as well as remote nodes.
// Its array positions are not global IDs. Reject duplicates instead of making
// a second, ambiguous mapping for the same physical node.
inline void ValidateRequiredNodeIds(const std::vector<std::int32_t>& nodes,
	std::uint64_t global_nodes)
{
	std::int64_t previous = -1;
	for (auto node : nodes) {
		if (node < 0 || static_cast<std::uint64_t>(node) >= global_nodes || node <= previous)
			throw std::invalid_argument("required nodes must be valid, sorted and unique global IDs");
		previous = node;
	}
}

} // namespace iga

#endif
