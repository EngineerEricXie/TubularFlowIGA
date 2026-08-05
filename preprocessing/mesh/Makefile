CXX ?= g++
CPPFLAGS ?= -Iinclude
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS ?= -lstdc++fs

TARGET := tubular_mesh
CORE := src/SwcGraph.cpp src/BSpline.cpp src/HexMesh.cpp src/MeshGenerator.cpp
OBJECTS := $(CORE:.cpp=.o)
TEST_TARGET := mesh_core_test

.PHONY: all test clean

all: $(TARGET)

$(TARGET): $(OBJECTS) src/main.o
	$(CXX) $^ -o $@ $(LDLIBS)

$(TEST_TARGET): $(OBJECTS) tests/test_core.o
	$(CXX) $^ -o $@ $(LDLIBS)

test: $(TEST_TARGET)
	./$(TEST_TARGET) ../../meshgeneration/template

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(TEST_TARGET) $(OBJECTS) src/main.o tests/test_core.o
