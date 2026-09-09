#include "ParallelElementBatch.hpp"

#include <iostream>
#include <string>
#include <thread>

namespace {
int checks = 0;
void Require(bool condition, const char* message)
{
	++checks;
	if (!condition) throw std::runtime_error(message);
}

void Check(int threads)
{
	const auto caller = std::this_thread::get_id();
	std::vector<std::size_t> consumed;
	std::size_t prepared_count = 0;
	const auto prepare = [&](std::size_t index) {
		Require(std::this_thread::get_id() == caller, "prepare entered a worker");
		Require(++prepared_count-consumed.size() <= 8, "unbounded prepared batch");
		return std::vector<double>(index%7+1, static_cast<double>(index));
	};
	const auto compute = [](const std::vector<double>& values, std::size_t) {
		std::vector<double> local = values;
		for (auto& value : local) value = value*value+1;
		return local;
	};
	const auto consume = [&](const auto&, const std::vector<double>& values, std::size_t index) {
		Require(std::this_thread::get_id() == caller, "consumer entered a worker");
		Require(index == consumed.size(), "insertion order changed");
		Require(values.size() == index%7+1, "worker scratch size contaminated");
		for (const auto value : values) Require(value == index*index+1, "worker scratch contaminated");
		consumed.push_back(index);
	};
	const auto statistics = iga::ForEachElementBatch(35, {threads,8}, prepare, compute, consume);
	Require(consumed.size() == 35 && statistics.items == 35, "missing batch work");
	Require(statistics.batches == 5 && statistics.maximum_resident_items == 8, "invalid batch bounds");
	Require(statistics.maximum_team_size == threads, "requested team was not exercised");
	consumed.clear(); prepared_count = 0;
	bool caught = false;
	try {
		iga::ForEachElementBatch(35, {threads,8}, prepare,
			[&](const auto& values, std::size_t index) {
				if (index == 10 || index == 12) throw std::runtime_error("worker-"+std::to_string(index));
				return compute(values,index);
			}, consume);
	} catch (const std::runtime_error& error) { caught = std::string(error.what()) == "worker-10"; }
	Require(caught, "worker exception escaped or selected nondeterministically");
	Require(consumed.size() == 8 && prepared_count == 16, "failed batch was partially consumed");
}
}

int main()
{
	try {
#ifdef _OPENMP
		omp_set_dynamic(0);
#endif
		Check(1);
#ifdef _OPENMP
		Check(2); Check(4);
#else
		bool rejected = false;
		try { iga::ForEachElementBatch(1,{2,1},[](std::size_t i){return i;},
			[](std::size_t i,std::size_t){return i;},[](auto,auto,auto){}); }
		catch (const std::invalid_argument&) { rejected = true; }
		Require(rejected, "serial build silently accepted multiple threads");
#endif
		std::cout << "parallel_element_batch checks=" << checks << " passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "parallel_element_batch: " << error.what() << '\n';
		return 1;
	}
}
