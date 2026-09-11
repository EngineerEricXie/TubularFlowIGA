#define main Phase9FixtureReferenceMain
#include "test_phase9_multiscale_closure.cpp"
#undef main

int main(int argc, char** argv)
{
	try {
		if (argc==3 && std::string(argv[1])=="--validate") {
			const auto metrics = Validate(argv[2], 2, 0.01);
			std::cout << "duct_validation=passed pressure_residual=" << metrics.max_pressure_residual
				<< " flow_residual=" << metrics.max_flow_residual
				<< " mass_imbalance=" << metrics.max_three_d_balance << '\n';
			return 0;
		}
		Require(argc==5, "usage: hpc_duct_solver_fixture ROOT TRANSVERSE AXIAL RANKS | --validate OUTPUT");
		auto positive = [](const char* text, unsigned long maximum, const char* name) {
			std::size_t used = 0; const auto value = std::stoul(text, &used);
			Require(used==std::string(text).size() && value>=1 && value<=maximum,
				std::string("fixture ")+name+" must be in [1,"+std::to_string(maximum)+"]");
			return static_cast<int>(value);
		};
		const int transverse = positive(argv[2],128,"transverse dimension");
		const int axial = positive(argv[3],128,"axial dimension");
		const int ranks = positive(argv[4],4096,"rank count");
		Require(ranks<=transverse*transverse*axial, "fixture requires at least one element per rank");
		const fs::path root = argv[1];
		Require(fs::create_directories(root), "fixture output already exists");
		PrepareCase(root, 0.01, 2, ranks, "duct.ntiga");
		iga::test::WriteC2SquareDuctDatabase(root/"duct.ntiga", transverse, ranks, axial);
		WriteThreeDCase(root/"root3d", transverse, axial, 0.01, 2);
		std::cout << "duct_fixture transverse=" << transverse << " axial=" << axial
			<< " elements=" << transverse*transverse*axial
			<< " nodes=" << (transverse+3)*(transverse+3)*(axial+3)
			<< " ranks=" << ranks << '\n';
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
