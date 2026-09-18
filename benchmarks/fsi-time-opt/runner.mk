include $(PETSC_DIR)/lib/petsc/conf/variables
CXX=mpicxx
FLAGS=-O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp -UNDEBUG -I$(SOURCE)/solvers/cpu/include -I$(SOURCE)/include -I$(SOURCE)/solvers/cpu/tests $(PETSC_CC_INCLUDES)
.PHONY: benchmark volume
benchmark:
	$(CXX) $(FLAGS) $(WORKLOAD_SOURCE)/benchmarks/fsi-time-opt/benchmark.cpp -o $(BUILD)/benchmark $(PETSC_LIB)
volume:
	$(CXX) $(FLAGS) -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING '-DIGA_WRAPPED_TEST="test_parallel_immersed_volume.cpp"' $(WORKLOAD_SOURCE)/benchmarks/fsi-time-opt/mpi_test_wrapper.cpp -o $(BUILD)/volume $(PETSC_LIB)
