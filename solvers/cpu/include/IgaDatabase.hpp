#ifndef IGA_DATABASE_HPP
#define IGA_DATABASE_HPP

#include <array>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

constexpr std::array<char, 8> kMagic{{'N', 'T', 'I', 'G', 'A', 'D', 'B', '2'}};
constexpr std::uint32_t kVersion = 3;
constexpr std::uint32_t kBezierPointCount = 64;

struct Header {
	std::array<char, 8> magic{};
	std::uint32_t version = 0;
	std::uint32_t ranks = 0;
	std::uint64_t elements = 0;
	std::uint64_t nodes = 0;
	std::uint32_t bezier_points = 0;
	std::uint32_t reserved = 0;
	std::uint64_t rank_index_offset = 0;
};

struct Element {
	std::uint64_t id = 0;
	std::int32_t type = 0;
	std::int32_t owner = 0;
	std::vector<std::int32_t> connectivity;
	std::vector<std::array<double, 64>> extraction;
	std::array<std::array<double, 3>, 64> bezier_points{};
};

template <class T>
inline void Write(std::ostream& out, const T& value)
{
	out.write(reinterpret_cast<const char*>(&value), sizeof(T));
	if (!out) throw std::runtime_error("binary database write failed");
}

template <class T>
inline void Read(std::istream& in, T& value)
{
	in.read(reinterpret_cast<char*>(&value), sizeof(T));
	if (!in) throw std::runtime_error("binary database is truncated");
}

inline Header ReadHeader(std::istream& in)
{
	Header h;
	in.read(h.magic.data(), h.magic.size());
	Read(in, h.version);
	Read(in, h.ranks);
	Read(in, h.elements);
	Read(in, h.nodes);
	Read(in, h.bezier_points);
	Read(in, h.reserved);
	Read(in, h.rank_index_offset);
	if (h.magic != kMagic) throw std::runtime_error("not a supported .ntiga database");
	if (h.version != kVersion) throw std::runtime_error("unsupported database version");
	if (h.bezier_points != kBezierPointCount) throw std::runtime_error("unsupported Bezier point count");
	if (h.nodes == 0 || h.rank_index_offset == 0) throw std::runtime_error("invalid database ownership index");
	return h;
}

class Database {
public:
	explicit Database(const std::string& path) : input_(path, std::ios::binary)
	{
		if (!input_) throw std::runtime_error("cannot open database: " + path);
		header_ = ReadHeader(input_);
		offsets_.resize(CheckedSize(header_.elements + 1));
		owners_.resize(CheckedSize(header_.elements));
		input_.read(reinterpret_cast<char*>(offsets_.data()),
			static_cast<std::streamsize>(offsets_.size() * sizeof(std::uint64_t)));
		input_.read(reinterpret_cast<char*>(owners_.data()),
			static_cast<std::streamsize>(owners_.size() * sizeof(std::int32_t)));
		if (!input_) throw std::runtime_error("database index is truncated");
		if (offsets_.empty() || offsets_.front() >= offsets_.back())
			throw std::runtime_error("invalid database offsets");
		input_.seekg(static_cast<std::streamoff>(header_.rank_index_offset));
		rank_offsets_.resize(static_cast<std::size_t>(header_.ranks) + 1);
		input_.read(reinterpret_cast<char*>(rank_offsets_.data()),
			static_cast<std::streamsize>(rank_offsets_.size() * sizeof(std::uint64_t)));
		if (!input_ || rank_offsets_.empty()) throw std::runtime_error("rank index is truncated");
		rank_elements_.resize(CheckedSize(rank_offsets_.back()));
		input_.read(reinterpret_cast<char*>(rank_elements_.data()),
			static_cast<std::streamsize>(rank_elements_.size() * sizeof(std::uint64_t)));
		if (!input_) throw std::runtime_error("rank element index is truncated");
	}

	const Header& header() const { return header_; }
	const std::vector<std::int32_t>& owners() const { return owners_; }

	Element Load(std::uint64_t index)
	{
		if (index >= header_.elements) throw std::out_of_range("element index");
		input_.clear();
		input_.seekg(static_cast<std::streamoff>(offsets_[CheckedSize(index)]));
		if (!input_) throw std::runtime_error("cannot seek to element record");
		Element e;
		std::uint32_t nen = 0;
		Read(input_, e.id);
		Read(input_, e.type);
		Read(input_, e.owner);
		Read(input_, nen);
		if (nen == 0 || nen > 4096) throw std::runtime_error("invalid element basis count");
		e.connectivity.resize(nen);
		e.extraction.resize(nen);
		input_.read(reinterpret_cast<char*>(e.connectivity.data()),
			static_cast<std::streamsize>(nen * sizeof(std::int32_t)));
		for (auto& row : e.extraction) {
			row.fill(0.0);
			std::uint8_t nonzeros = 0;
			Read(input_, nonzeros);
			for (std::uint8_t j = 0; j < nonzeros; ++j) {
				std::uint8_t column = 0;
				double value = 0.0;
				Read(input_, column);
				Read(input_, value);
				if (column >= kBezierPointCount) throw std::runtime_error("invalid extraction column");
				row[column] = value;
			}
		}
		input_.read(reinterpret_cast<char*>(e.bezier_points.data()),
			static_cast<std::streamsize>(sizeof(e.bezier_points)));
		if (!input_) throw std::runtime_error("element record is truncated");
		return e;
	}

	std::vector<Element> LoadOwned(std::int32_t rank)
	{
		if (rank < 0 || static_cast<std::uint32_t>(rank) >= header_.ranks)
			throw std::out_of_range("rank");
		std::vector<Element> result;
		for (std::uint64_t i = 0; i < header_.elements; ++i)
			if (owners_[CheckedSize(i)] == rank) result.push_back(Load(i));
		return result;
	}

	std::vector<std::uint64_t> RequiredElementIndices(std::int32_t rank) const
	{
		if (rank < 0 || static_cast<std::uint32_t>(rank) >= header_.ranks)
			throw std::out_of_range("rank");
		const auto begin = rank_offsets_[static_cast<std::size_t>(rank)];
		const auto end = rank_offsets_[static_cast<std::size_t>(rank) + 1];
		return {rank_elements_.begin() + static_cast<std::ptrdiff_t>(begin),
			rank_elements_.begin() + static_cast<std::ptrdiff_t>(end)};
	}

	std::vector<Element> LoadRequired(std::int32_t rank)
	{
		auto indices = RequiredElementIndices(rank);
		std::vector<Element> result;
		result.reserve(indices.size());
		for (auto index : indices) result.push_back(Load(index));
		return result;
	}

	std::pair<std::uint64_t, std::uint64_t> NodeRange(std::int32_t rank) const
	{
		if (rank < 0 || static_cast<std::uint32_t>(rank) >= header_.ranks)
			throw std::out_of_range("rank");
		const auto q = header_.nodes / header_.ranks;
		const auto rem = header_.nodes % header_.ranks;
		const auto r = static_cast<std::uint64_t>(rank);
		const auto begin = r * q + std::min(r, rem);
		const auto end = begin + q + (r < rem ? 1 : 0);
		return {begin, end};
	}

private:
	static std::size_t CheckedSize(std::uint64_t n)
	{
		if (n > std::numeric_limits<std::size_t>::max()) throw std::overflow_error("database too large");
		return static_cast<std::size_t>(n);
	}

	std::ifstream input_;
	Header header_;
	std::vector<std::uint64_t> offsets_;
	std::vector<std::int32_t> owners_;
	std::vector<std::uint64_t> rank_offsets_;
	std::vector<std::uint64_t> rank_elements_;
};

} // namespace iga

#endif
