#ifndef IGA_EXACT_DYADIC_PREDICATES_HPP
#define IGA_EXACT_DYADIC_PREDICATES_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {
namespace exact_dyadic {

// A finite binary64 value is an integer mantissa times a power of two.  This
// small signed-integer implementation keeps the predicate fallback exact
// without adding a multiprecision library dependency.
struct Unsigned {
	std::vector<std::uint64_t> limbs;

	bool IsZero() const { return limbs.empty(); }
	void Trim() { while (!limbs.empty() && limbs.back() == 0) limbs.pop_back(); }
};

inline void CheckSize(std::size_t size)
{
	if (size > 128) throw std::overflow_error("exact dyadic predicate magnitude overflows");
}
inline int Compare(const Unsigned& left, const Unsigned& right)
{
	if (left.limbs.size() != right.limbs.size()) return left.limbs.size() < right.limbs.size() ? -1 : 1;
	for (std::size_t i = left.limbs.size(); i != 0; --i)
		if (left.limbs[i-1] != right.limbs[i-1]) return left.limbs[i-1] < right.limbs[i-1] ? -1 : 1;
	return 0;
}
inline Unsigned AddMagnitude(const Unsigned& left, const Unsigned& right)
{
	const std::size_t largest = left.limbs.size() > right.limbs.size() ? left.limbs.size() : right.limbs.size();
	if (largest >= 128) throw std::overflow_error("exact dyadic predicate magnitude overflows");
	Unsigned result; const std::size_t size = largest+1;
	CheckSize(size); result.limbs.resize(size, 0);
	std::uint64_t carry = 0;
	for (std::size_t i = 0; i+1 < size; ++i) {
		const std::uint64_t first = (i < left.limbs.size() ? left.limbs[i] : 0), second = (i < right.limbs.size() ? right.limbs[i] : 0);
		const std::uint64_t partial = first+second; const bool first_carry = partial < first;
		result.limbs[i] = partial+carry; carry = (first_carry || result.limbs[i] < partial) ? 1 : 0;
	}
	result.limbs.back() = carry; result.Trim(); return result;
}
inline Unsigned SubtractMagnitude(const Unsigned& left, const Unsigned& right)
{
	if (Compare(left, right) < 0) throw std::logic_error("exact dyadic magnitude subtraction underflows");
	Unsigned result; result.limbs.resize(left.limbs.size(), 0); std::uint64_t borrow = 0;
	for (std::size_t i = 0; i < left.limbs.size(); ++i) {
		const std::uint64_t value = i < right.limbs.size() ? right.limbs[i] : 0;
		const std::uint64_t middle = left.limbs[i]-value;
		const bool first_borrow = left.limbs[i] < value;
		result.limbs[i] = middle-borrow;
		borrow = (first_borrow || middle < borrow) ? 1 : 0;
	}
	if (borrow != 0) throw std::logic_error("exact dyadic magnitude borrow remains");
	result.Trim(); return result;
}
inline Unsigned ShiftLeft(const Unsigned& value, std::size_t bits)
{
	if (value.IsZero()) return {};
	const std::size_t words = bits/64, remainder = bits%64;
	if (words >= 128 || value.limbs.size() > 127-words) throw std::overflow_error("exact dyadic predicate magnitude overflows");
	const std::size_t size = value.limbs.size()+words+1;
	Unsigned result; result.limbs.assign(words, 0); result.limbs.resize(size, 0);
	std::uint64_t carry = 0;
	for (std::size_t i = 0; i < value.limbs.size(); ++i) {
		result.limbs[i+words] = remainder == 0 ? value.limbs[i] : (value.limbs[i] << remainder)|carry;
		carry = remainder == 0 ? 0 : value.limbs[i] >> (64-remainder);
	}
	result.limbs[value.limbs.size()+words] = carry; result.Trim(); return result;
}
inline Unsigned MultiplyMagnitude(const Unsigned& left, const Unsigned& right)
{
	if (left.IsZero() || right.IsZero()) return {};
	if (left.limbs.size() > 128 || right.limbs.size() > 128 || left.limbs.size() > 128-right.limbs.size()) throw std::overflow_error("exact dyadic predicate magnitude overflows");
	Unsigned result;
	for (std::size_t word = 0; word < right.limbs.size(); ++word) for (std::size_t bit = 0; bit < 64; ++bit)
		if ((right.limbs[word] & (UINT64_C(1) << bit)) != 0) {
			if (word > (std::numeric_limits<std::size_t>::max()-bit)/64) throw std::overflow_error("exact dyadic predicate shift overflows");
			result = AddMagnitude(result, ShiftLeft(left, 64*word+bit));
		}
	return result;
}

struct Number {
	int sign = 0;
	int exponent = 0;
	Unsigned magnitude;
};
inline Number FromDouble(double value)
{
	static_assert(std::numeric_limits<double>::is_iec559 && sizeof(double) == sizeof(std::uint64_t), "exact dyadic predicates require binary64 IEEE-754 doubles");
	std::uint64_t bits = 0; std::memcpy(&bits, &value, sizeof(bits));
	const std::uint64_t fraction = bits&UINT64_C(0x000fffffffffffff); const unsigned exponent_bits = static_cast<unsigned>((bits >> 52)&UINT64_C(0x7ff));
	if (exponent_bits == 0x7ff) throw std::invalid_argument("exact dyadic predicate received nonfinite coordinate");
	if (exponent_bits == 0 && fraction == 0) return {};
	Number result; result.sign = (bits >> 63) == 0 ? 1 : -1; result.exponent = exponent_bits == 0 ? -1074 : static_cast<int>(exponent_bits)-1075;
	result.magnitude.limbs.push_back(exponent_bits == 0 ? fraction : (fraction|UINT64_C(0x0010000000000000))); return result;
}
inline Number Negate(Number value) { value.sign = -value.sign; return value; }
inline Number Add(const Number& left, const Number& right)
{
	if (left.sign == 0) return right;
	if (right.sign == 0) return left;
	const int exponent = left.exponent < right.exponent ? left.exponent : right.exponent;
	const long long left_delta = static_cast<long long>(left.exponent)-static_cast<long long>(exponent), right_delta = static_cast<long long>(right.exponent)-static_cast<long long>(exponent);
	if (left_delta < 0 || right_delta < 0 || static_cast<unsigned long long>(left_delta) > std::numeric_limits<std::size_t>::max() || static_cast<unsigned long long>(right_delta) > std::numeric_limits<std::size_t>::max()) throw std::overflow_error("exact dyadic predicate exponent difference overflows");
	const std::size_t left_shift = static_cast<std::size_t>(left_delta), right_shift = static_cast<std::size_t>(right_delta);
	const auto left_magnitude = ShiftLeft(left.magnitude, left_shift), right_magnitude = ShiftLeft(right.magnitude, right_shift);
	Number result; result.exponent = exponent;
	if (left.sign == right.sign) { result.sign = left.sign; result.magnitude = AddMagnitude(left_magnitude, right_magnitude); }
	else { const int comparison = Compare(left_magnitude, right_magnitude); if (comparison == 0) return {}; result.sign = comparison > 0 ? left.sign : right.sign; result.magnitude = comparison > 0 ? SubtractMagnitude(left_magnitude, right_magnitude) : SubtractMagnitude(right_magnitude, left_magnitude); }
	if (result.magnitude.IsZero()) return {};
	return result;
}
inline Number Subtract(const Number& left, const Number& right) { return Add(left, Negate(right)); }
inline Number Multiply(const Number& left, const Number& right)
{
	if (left.sign == 0 || right.sign == 0) return {};
	Number result; result.sign = left.sign*right.sign;
	const long long exponent = static_cast<long long>(left.exponent)+static_cast<long long>(right.exponent);
	if (exponent > std::numeric_limits<int>::max() || exponent < std::numeric_limits<int>::min()) throw std::overflow_error("exact dyadic predicate exponent overflows");
	result.exponent = static_cast<int>(exponent); result.magnitude = MultiplyMagnitude(left.magnitude, right.magnitude); return result;
}
inline int Sign(const Number& value) { return value.sign; }
inline int CompareAbsolute(const Number& left, const Number& right)
{
	if (left.sign == 0 || right.sign == 0) return left.sign == 0 ? (right.sign == 0 ? 0 : -1) : 1;
	auto highest_bit = [](const Number& value) {
		std::uint64_t limb = value.magnitude.limbs.back(); std::size_t bits = 0;
		while (limb != 0) { ++bits; limb >>= 1; }
		return static_cast<long long>(value.exponent)+static_cast<long long>((value.magnitude.limbs.size()-1)*64+bits);
	};
	const long long left_high = highest_bit(left), right_high = highest_bit(right);
	if (left_high != right_high) return left_high < right_high ? -1 : 1;
	const int exponent = left.exponent < right.exponent ? left.exponent : right.exponent;
	const auto left_magnitude = ShiftLeft(left.magnitude, static_cast<std::size_t>(static_cast<long long>(left.exponent)-exponent));
	const auto right_magnitude = ShiftLeft(right.magnitude, static_cast<std::size_t>(static_cast<long long>(right.exponent)-exponent));
	return Compare(left_magnitude, right_magnitude);
}
using Point = std::array<Number, 3>;
inline Point PointFromDouble(const std::array<double, 3>& point) { return {{FromDouble(point[0]), FromDouble(point[1]), FromDouble(point[2])}}; }
inline Point SubtractPoint(const Point& left, const Point& right) { return {{Subtract(left[0], right[0]), Subtract(left[1], right[1]), Subtract(left[2], right[2])}}; }
inline Point Cross(const Point& left, const Point& right) { return {{Subtract(Multiply(left[1], right[2]), Multiply(left[2], right[1])), Subtract(Multiply(left[2], right[0]), Multiply(left[0], right[2])), Subtract(Multiply(left[0], right[1]), Multiply(left[1], right[0]))}}; }
inline Number Dot(const Point& left, const Point& right) { return Add(Add(Multiply(left[0], right[0]), Multiply(left[1], right[1])), Multiply(left[2], right[2])); }
inline int CoordinateDifferenceSign(double left, double right) { return Sign(Subtract(FromDouble(left), FromDouble(right))); }
enum class CollinearIntervalContact { none, point, overlap };
inline CollinearIntervalContact ClassifyCollinearIntervals(const std::array<double, 3>& a, const std::array<double, 3>& b,
	const std::array<double, 3>& c, const std::array<double, 3>& d)
{
	std::array<Number, 3> direction{{Subtract(FromDouble(b[0]), FromDouble(a[0])), Subtract(FromDouble(b[1]), FromDouble(a[1])), Subtract(FromDouble(b[2]), FromDouble(a[2]))}};
	std::size_t axis = 0; while (axis < 3 && Sign(direction[axis]) == 0) ++axis;
	if (axis == 3) throw std::invalid_argument("collinear interval has a zero first segment");
	for (std::size_t candidate = axis+1; candidate < 3; ++candidate)
		if (Sign(direction[candidate]) != 0 && CompareAbsolute(direction[candidate], direction[axis]) > 0) axis = candidate;
	auto less = [axis](double left, double right) { return CoordinateDifferenceSign(left, right) < 0; };
	const double low_a = less(a[axis], b[axis]) ? a[axis] : b[axis], high_a = low_a == a[axis] ? b[axis] : a[axis];
	const double low_c = less(c[axis], d[axis]) ? c[axis] : d[axis], high_c = low_c == c[axis] ? d[axis] : c[axis];
	const double low = less(low_a, low_c) ? low_c : low_a, high = less(high_a, high_c) ? high_a : high_c;
	const int relation = CoordinateDifferenceSign(high, low);
	return relation < 0 ? CollinearIntervalContact::none : (relation == 0 ? CollinearIntervalContact::point : CollinearIntervalContact::overlap);
}
inline int Orient2D(const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c, std::size_t x, std::size_t y)
{
	const auto ae = PointFromDouble(a), be = PointFromDouble(b), ce = PointFromDouble(c); const auto u = SubtractPoint(be, ae), v = SubtractPoint(ce, ae); return Sign(Subtract(Multiply(u[x], v[y]), Multiply(u[y], v[x])));
}
inline int Orient3D(const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c, const std::array<double, 3>& p)
{
	const auto ae = PointFromDouble(a), be = PointFromDouble(b), ce = PointFromDouble(c), pe = PointFromDouble(p); return Sign(Dot(Cross(SubtractPoint(be, ae), SubtractPoint(ce, ae)), SubtractPoint(pe, ae)));
}
inline int DirectionCrossComponent(const std::array<double, 3>& a, const std::array<double, 3>& b,
	const std::array<double, 3>& c, const std::array<double, 3>& d, std::size_t x, std::size_t y)
{
	const auto ae = PointFromDouble(a), be = PointFromDouble(b), ce = PointFromDouble(c), de = PointFromDouble(d);
	const auto u = SubtractPoint(be, ae), v = SubtractPoint(de, ce); return Sign(Subtract(Multiply(u[x], v[y]), Multiply(u[y], v[x])));
}
inline std::array<Number, 4> SegmentTriangleValues(const std::array<double, 3>& start, const std::array<double, 3>& end,
	const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c)
{
	const auto se = PointFromDouble(start), ee = PointFromDouble(end), ae = PointFromDouble(a), be = PointFromDouble(b), ce = PointFromDouble(c);
	const auto direction = SubtractPoint(ee, se), e1 = SubtractPoint(be, ae), e2 = SubtractPoint(ce, ae), tvec = SubtractPoint(se, ae);
	const auto p = Cross(direction, e2), q = Cross(tvec, e1);
	return {{Dot(e1, p), Dot(tvec, p), Dot(direction, q), Dot(e2, q)}};
}

} // namespace exact_dyadic
} // namespace iga

#endif
