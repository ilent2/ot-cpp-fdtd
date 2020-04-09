# Build simulation

FDTD_PATH := ../../include/OtFdtd/
OTT_PATH := ../../include/OttToolbox

CXX				:= g++
CXXFLAGS	:= -Wall -O3 --std=c++14 -I$(FDTD_PATH) -I$(OTT_PATH)

TARGET		:= deflection

all: $(TARGET)

$(TARGET): sim.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run:
	./$(TARGET)

.PHONY: all run

