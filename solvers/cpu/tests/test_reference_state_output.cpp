#include "ReferenceStateOutput.hpp"

#include <iostream>
#include <limits>

namespace {
int checks = 0;
void Require(bool condition, const char* message)
{
	++checks;
	if (!condition) throw std::runtime_error(message);
}
template <class Action> void Reject(Action action)
{
	bool rejected = false;
	try { action(); } catch (const std::exception&) { rejected = true; }
	Require(rejected, "expected rejection");
}
}

int main(int argc, char** argv)
{
	try {
		if (argc != 2) throw std::invalid_argument("usage: reference_state_output_test NEW_DIRECTORY");
		const std::filesystem::path root(argv[1]);
		if (!std::filesystem::create_directory(root)) throw std::runtime_error("test directory exists");
		const std::vector<std::uint64_t> ids{9007199254740993ULL, 9007199254740995ULL};
		iga::test::ReferenceStateOutput output(root/"valid", "fsi");
		const auto value = [](std::size_t row, std::size_t column) { return (row+1)*0.1+column; };
		output.Add("velocity", "m/s", ids, 3, value);
		Require(!std::filesystem::exists(root/"valid/manifest.json"), "premature completion manifest");
		output.Finish(1.05, 1);
		Require(std::filesystem::exists(root/"valid/manifest.json"), "missing completion manifest");
		std::ifstream input(root/"valid/velocity.txt");
		for (std::size_t row = 0; row < ids.size(); ++row) {
			std::uint64_t id = 0; input >> id;
			Require(id == ids[row], "64-bit material ID rounded");
			for (std::size_t column = 0; column < 3; ++column) {
				double actual = 0; input >> actual;
				Require(actual == value(row, column), "field value did not round-trip exactly");
			}
		}
		Reject([&] { output.Finish(1.05, 1); });
		Reject([&] { output.Add("another", "m", ids, 1, value); });
		Reject([&] { iga::test::ReferenceStateOutput overwrite(root/"valid", "fsi"); });
		iga::test::ReferenceStateOutput bad(root/"bad", "immersed");
		Reject([&] { bad.Add("pressure", "pa", ids, 1, [](std::size_t, std::size_t) {
			return std::numeric_limits<double>::quiet_NaN(); }); });
		Reject([&] { bad.Finish(0.0, 0); });
		Require(!std::filesystem::exists(root/"bad/manifest.json"), "failed output published manifest");
		iga::test::ReferenceStateOutput duplicate(root/"duplicate", "immersed");
		Reject([&] { duplicate.Add("pressure", "pa", std::vector<int>{1,1}, 1, value); });
		iga::test::ReferenceStateOutput negative(root/"negative", "immersed");
		Reject([&] { negative.Add("pressure", "pa", std::vector<int>{-1}, 1, value); });
		iga::test::ReferenceStateOutput unsafe(root/"unsafe", "immersed");
		Reject([&] { unsafe.Add("../field", "pa", ids, 1, value); });
		iga::test::ReferenceStateOutput io(root/"io", "immersed");
		std::filesystem::create_directory(root/"io/pressure.txt");
		Reject([&] { io.Add("pressure", "pa", ids, 1, value); });
		Reject([&] { io.Finish(0.0, 0); });
		Require(!std::filesystem::exists(root/"io/manifest.json"), "I/O failure published manifest");
		std::cout << "reference_state_output checks=" << checks << " passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "reference_state_output: " << error.what() << '\n';
		return 1;
	}
}
