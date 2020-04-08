# Build simulation

FDTD_PATH := ../../include/OtFdtd/

CXX				:= g++
CXXFLAGS	:= -Wall -O3 --std=c++14 -I$(FDTD_PATH)

TARGET		:= forceCalc

all: $(TARGET)

$(TARGET): sim.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run:
	./$(TARGET)

.PHONY: all run

