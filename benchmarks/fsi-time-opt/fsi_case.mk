include $(PETSC_DIR)/lib/petsc/conf/variables
CXX=mpicxx
FLAGS=-O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp -UNDEBUG -I$(SOURCE)/solvers/cpu/include -I$(SOURCE)/include -I$(SOURCE)/solvers/cpu/tests $(PETSC_CC_INCLUDES)
.PHONY: fsi_case
fsi_case:
	$(CXX) $(FLAGS) -DIGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING '-DIGA_WRAPPED_TEST="$(CASE_SOURCE)"' $(SOURCE)/benchmarks/fsi-time-opt/mpi_test_wrapper.cpp -o $(BUILD)/$(CASE_NAME) $(PETSC_LIB)
