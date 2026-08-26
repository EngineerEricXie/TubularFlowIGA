#include "BoundarySupport.hpp"

#include <array>
#include <iostream>
#include <set>
#include <stdexcept>

namespace {

void Require(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

std::set<int> Set(const std::array<int, 16>& values)
{
	return {values.begin(), values.end()};
}

} // namespace

int main()
{
	try {
		for (int face = 0; face < 6; ++face)
			Require(Set(iga::FaceBezierColumns(face)).size() == 16,
				"Bezier face columns are not unique");
		Require(Set(iga::FaceBezierColumns(0))
			== std::set<int>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}),
			"w=0 mapping failed");
		Require(Set(iga::FaceBezierColumns(5))
			== std::set<int>({48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}),
			"w=1 mapping failed");
		Require(Set(iga::FaceBezierColumns(1))
			== std::set<int>({0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51}),
			"v=0 mapping failed");
		Require(Set(iga::FaceBezierColumns(2))
			== std::set<int>({3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63}),
			"u=1 mapping failed");
		Require(Set(iga::FaceBezierColumns(3))
			== std::set<int>({12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63}),
			"v=1 mapping failed");
		Require(Set(iga::FaceBezierColumns(4))
			== std::set<int>({0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60}),
			"u=0 mapping failed");
		std::cout << "boundary support tests passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
