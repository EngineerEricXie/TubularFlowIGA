#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace tubular {

struct Vec3
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	Vec3& operator+=(const Vec3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	Vec3& operator-=(const Vec3& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
	Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }
};

inline Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
inline Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
inline Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(Vec3 a, double s) { return a *= s; }
inline Vec3 operator*(double s, Vec3 a) { return a *= s; }
inline Vec3 operator/(Vec3 a, double s) { return a /= s; }

inline double Dot(const Vec3& a, const Vec3& b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Vec3 Cross(const Vec3& a, const Vec3& b)
{
	return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

inline double NormSquared(const Vec3& a) { return Dot(a, a); }
inline double Norm(const Vec3& a) { return std::sqrt(NormSquared(a)); }
inline bool IsFinite(const Vec3& a) { return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z); }

inline Vec3 Normalized(const Vec3& a, const std::string& context)
{
	const double n = Norm(a);
	if (!std::isfinite(n) || n <= 1.0e-14)
		throw std::runtime_error(context + ": expected a finite nonzero vector");
	return a/n;
}

inline double ClampUnit(double value)
{
	return std::max(-1.0, std::min(1.0, value));
}

inline Vec3 RotateAroundAxis(const Vec3& point, const Vec3& axis, double angle)
{
	const Vec3 unit = Normalized(axis, "RotateAroundAxis");
	const double c = std::cos(angle);
	const double s = std::sin(angle);
	return point*c + Cross(unit, point)*s + unit*Dot(unit, point)*(1.0-c);
}

inline Vec3 RotateSurface(const Vec3& point, const Vec3& old_normal, const Vec3& new_normal)
{
	const Vec3 from = Normalized(old_normal, "RotateSurface old normal");
	const Vec3 to = Normalized(new_normal, "RotateSurface new normal");
	const double cosine = ClampUnit(Dot(from, to));
	if (cosine >= 1.0-1.0e-14) return point;
	Vec3 axis;
	double angle;
	if (cosine <= -1.0+1.0e-14) {
		const std::array<double,3> magnitudes{std::abs(from.x), std::abs(from.y), std::abs(from.z)};
		const auto least = static_cast<int>(std::min_element(magnitudes.begin(), magnitudes.end())-magnitudes.begin());
		Vec3 helper;
		if (least == 0) helper.x = 1.0;
		else if (least == 1) helper.y = 1.0;
		else helper.z = 1.0;
		axis = Cross(from, helper);
		angle = std::acos(-1.0);
	} else {
		axis = Cross(from, to);
		angle = std::acos(cosine);
	}
	return RotateAroundAxis(point, axis, angle);
}

inline Vec3 ProjectToLine(const Vec3& a, const Vec3& b, const Vec3& point)
{
	const Vec3 direction = a-b;
	const double length_squared = NormSquared(direction);
	if (!std::isfinite(length_squared) || length_squared <= 1.0e-28)
		throw std::runtime_error("ProjectToLine: coincident line points");
	return b + direction*(Dot(direction, point-b)/length_squared);
}

inline Vec3 ProjectToPlane(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& point)
{
	const Vec3 normal = Normalized(Cross(b-a, c-a), "ProjectToPlane");
	return point-normal*Dot(point-a, normal);
}

inline double Determinant(const std::array<Vec3,3>& columns)
{
	return Dot(columns[0], Cross(columns[1], columns[2]));
}

} // namespace tubular
