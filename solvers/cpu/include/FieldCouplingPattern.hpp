#ifndef IGA_FIELD_COUPLING_PATTERN_HPP
#define IGA_FIELD_COUPLING_PATTERN_HPP

#include <petscsys.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

class FieldCouplingPattern {
public:
	explicit FieldCouplingPattern(std::size_t fields = 0)
		: fields_(fields), pair_index_(fields*fields, -1) {}

	static FieldCouplingPattern Dense(std::size_t fields)
	{
		FieldCouplingPattern result(fields);
		for (std::size_t equation = 0; equation < fields; ++equation)
			for (std::size_t trial = 0; trial < fields; ++trial)
				result.Add(equation, trial);
		return result;
	}

	void Add(std::size_t equation, std::size_t trial)
	{
		CheckField(equation);
		CheckField(trial);
		auto& index = pair_index_[equation*fields_+trial];
		if (index >= 0) return;
		index = static_cast<int>(pairs_.size());
		pairs_.push_back({equation, trial});
	}

	std::size_t fields() const { return fields_; }
	std::size_t pairs() const { return pairs_.size(); }
	const std::vector<std::pair<std::size_t, std::size_t>>& active_pairs() const { return pairs_; }

	bool Active(std::size_t equation, std::size_t trial) const
	{
		CheckField(equation);
		CheckField(trial);
		return pair_index_[equation*fields_+trial] >= 0;
	}

	std::size_t PairIndex(std::size_t equation, std::size_t trial) const
	{
		CheckField(equation);
		CheckField(trial);
		const auto index = pair_index_[equation*fields_+trial];
		if (index < 0) throw std::invalid_argument("inactive field-coupling pair");
		return static_cast<std::size_t>(index);
	}

	std::size_t ActiveTrials(std::size_t equation) const
	{
		CheckField(equation);
		return static_cast<std::size_t>(std::count_if(pairs_.begin(), pairs_.end(),
			[equation](const auto& pair) { return pair.first == equation; }));
	}

private:
	void CheckField(std::size_t field) const
	{
		if (field >= fields_) throw std::out_of_range("field-coupling index");
	}

	std::size_t fields_ = 0;
	std::vector<int> pair_index_;
	std::vector<std::pair<std::size_t, std::size_t>> pairs_;
};

class FieldBlockElementMatrix {
public:
	FieldBlockElementMatrix() = default;
	explicit FieldBlockElementMatrix(FieldCouplingPattern pattern)
		: pattern_(std::move(pattern)) {}

	void Reset(std::size_t nodes)
	{
		nodes_ = nodes;
		values_.resize(pattern_.pairs()*nodes_*nodes_);
		std::fill(values_.begin(), values_.end(), PetscScalar{0.0});
	}

	PetscScalar& At(std::size_t equation, std::size_t trial,
		std::size_t row_node, std::size_t column_node)
	{
		return values_.at(Offset(equation, trial, row_node, column_node));
	}

	const PetscScalar& At(std::size_t equation, std::size_t trial,
		std::size_t row_node, std::size_t column_node) const
	{
		return values_.at(Offset(equation, trial, row_node, column_node));
	}

	const PetscScalar* Block(std::size_t pair) const
	{
		if (pair >= pattern_.pairs()) throw std::out_of_range("element field block");
		return values_.data()+static_cast<std::ptrdiff_t>(pair*nodes_*nodes_);
	}

	PetscScalar* Block(std::size_t equation, std::size_t trial)
	{
		return values_.data()+static_cast<std::ptrdiff_t>(
			pattern_.PairIndex(equation, trial)*nodes_*nodes_);
	}

	std::size_t nodes() const { return nodes_; }
	const FieldCouplingPattern& pattern() const { return pattern_; }
	const std::vector<PetscScalar>& values() const { return values_; }

private:
	std::size_t Offset(std::size_t equation, std::size_t trial,
		std::size_t row_node, std::size_t column_node) const
	{
		if (row_node >= nodes_ || column_node >= nodes_)
			throw std::out_of_range("element block node");
		return pattern_.PairIndex(equation, trial)*nodes_*nodes_+row_node*nodes_+column_node;
	}

	FieldCouplingPattern pattern_;
	std::size_t nodes_ = 0;
	std::vector<PetscScalar> values_;
};

} // namespace iga

#endif
