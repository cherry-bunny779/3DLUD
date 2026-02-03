# Block LU Decomposition Simulator Makefile

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g
INCLUDES = -I./include

SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
TARGET = block_lu_sim

.PHONY: all clean run test verbose

all: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Run with default settings
run: $(TARGET)
	./$(TARGET)

# Run with small matrix for testing
test: $(TARGET)
	./$(TARGET) -n 8 -p 4 -b 4 -mac 1 -div 10

# Run with verbose output
verbose: $(TARGET)
	./$(TARGET) -n 8 -p 4 -b 4 -v

# Run with larger matrix
large: $(TARGET)
	./$(TARGET) -n 64 -p 8 -b 8 -mac 1 -div 10

# Header dependencies
$(BUILD_DIR)/main.o: $(INC_DIR)/block_lu_simulator.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/block_lu_simulator.o: $(INC_DIR)/block_lu_simulator.hpp $(INC_DIR)/pe_array.hpp $(INC_DIR)/memory.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/pe_array.o: $(INC_DIR)/pe_array.hpp $(INC_DIR)/pe.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/pe.o: $(INC_DIR)/pe.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/memory.o: $(INC_DIR)/memory.hpp $(INC_DIR)/config.hpp
