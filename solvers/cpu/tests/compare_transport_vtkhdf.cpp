#include <hdf5.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

std::vector<double> Read(hid_t file, const char* name, std::vector<hsize_t>& dimensions)
{
	const auto dataset = H5Dopen2(file, name, H5P_DEFAULT);
	if (dataset < 0) throw std::runtime_error("missing dataset");
	const auto space = H5Dget_space(dataset);
	const int rank = H5Sget_simple_extent_ndims(space);
	if (rank < 1) throw std::runtime_error("bad dataset rank");
	dimensions.resize(rank);
	if (H5Sget_simple_extent_dims(space, dimensions.data(), nullptr) < 0) throw std::runtime_error("bad dimensions");
	std::size_t count = 1;
	for (auto d : dimensions) count *= d;
	std::vector<double> values(count);
	if (H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0)
		throw std::runtime_error("cannot read dataset");
	H5Sclose(space); H5Dclose(dataset);
	return values;
}
int main(int argc, char** argv)
{
	if (argc != 3) return 2;
	try {
		const auto first = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
		const auto second = H5Fopen(argv[2], H5F_ACC_RDONLY, H5P_DEFAULT);
		if (first < 0 || second < 0) throw std::runtime_error("cannot open HDF output");
		for (const char* name : {"/VTKHDF/Steps/Values", "/VTKHDF/Points", "/VTKHDF/Connectivity",
			"/VTKHDF/Offsets", "/VTKHDF/PointData/three_red", "/VTKHDF/PointData/three_blue"}) {
			std::vector<hsize_t> da, db;
			const auto a = Read(first, name, da), b = Read(second, name, db);
			if (da != db || a.empty()) throw std::runtime_error("shape differs or empty dataset");
			double norm = 0, delta = 0;
			for (std::size_t i=0; i<a.size(); ++i) {
				if (!std::isfinite(a[i]) || !std::isfinite(b[i])) throw std::runtime_error("nonfinite HDF value");
				norm = std::hypot(norm,a[i]); delta = std::hypot(delta,b[i]-a[i]);
			}
			if (norm ? delta/norm > 1e-6 : delta > 1e-12) throw std::runtime_error("HDF numeric gate failed");
			std::cout.precision(17);
			std::cout << name << " values=" << a.size() << " relative_l2=" << (norm ? delta/norm : delta) << '\n';
		}
		H5Fclose(second); H5Fclose(first);
		std::cout << "HDF comparison passed\n";
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
