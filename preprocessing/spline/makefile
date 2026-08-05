CXX ?= g++
CPPFLAGS ?= -I/usr/include/eigen3
CXXFLAGS ?= -O3 -std=c++11 -Wall -Wextra -fopenmp
LDLIBS ?= -fopenmp -lpthread -lm -ldl

TARGET := spline
OBJECTS := BasicDataStructure.o kernel.o KnotInsertion.o main.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(OBJECTS)
