include $(PETSC_DIR)/lib/petsc/conf/variables
CXX=mpicxx
FLAGS=-O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp -UNDEBUG -I$(SOURCE)/solvers/cpu/include -I$(SOURCE)/include -I$(WORKLOAD_SOURCE)/solvers/cpu/include $(PETSC_CC_INCLUDES)
.PHONY: heartbeat native_reference_test
heartbeat:
	$(CXX) $(FLAGS) $(WORKLOAD_SOURCE)/solvers/cpu/src/bifurcation_fsi.cpp -o $(BUILD)/heartbeat $(PETSC_LIB)
native_reference_test:
	$(CXX) $(FLAGS) $(WORKLOAD_SOURCE)/solvers/cpu/tests/test_native_fsi_reference_output.cpp -o $(BUILD)/native_reference_test
