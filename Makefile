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

.PHONY: all clean run test verbose test3d compare large

all: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Run with default settings (2D)
run: $(TARGET)
	./$(TARGET)

# Run with small matrix for testing (2D)
test: $(TARGET)
	./$(TARGET) -n 8 -p 4 -b 4 -mac 1 -div 10

# Run with verbose output (2D)
verbose: $(TARGET)
	./$(TARGET) -n 8 -p 4 -b 4 -v

# Run with larger matrix (2D)
large: $(TARGET)
	./$(TARGET) -n 64 -p 8 -b 8 -mac 1 -div 10

# Run 3D simulation
test3d: $(TARGET)
	./$(TARGET) -n 8 -p 4 -b 4 -3d -z 4 -tsv 1

# Run 3D simulation with verbose output
verbose3d: $(TARGET)
	./$(TARGET) -n 8 -p 4 -b 4 -3d -z 4 -v

# Run comparison mode (2D vs 3D)
compare: $(TARGET)
	./$(TARGET) -n 16 -p 4 -b 4 -compare -z 4

# Run larger comparison
compare_large: $(TARGET)
	./$(TARGET) -n 32 -p 8 -b 8 -compare -z 4

# Test with different layer counts
compare_layers: $(TARGET)
	@echo "=== 2 Layers ===" && ./$(TARGET) -n 16 -p 4 -b 4 -compare -z 2
	@echo "\n=== 4 Layers ===" && ./$(TARGET) -n 16 -p 4 -b 4 -compare -z 4
	@echo "\n=== 8 Layers ===" && ./$(TARGET) -n 16 -p 4 -b 4 -compare -z 8

# Header dependencies - 2D
$(BUILD_DIR)/main.o: $(INC_DIR)/block_lu_simulator.hpp $(INC_DIR)/block_lu_simulator_3d.hpp $(INC_DIR)/config.hpp $(INC_DIR)/config_3d.hpp
$(BUILD_DIR)/block_lu_simulator.o: $(INC_DIR)/block_lu_simulator.hpp $(INC_DIR)/pe_array.hpp $(INC_DIR)/memory.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/pe_array.o: $(INC_DIR)/pe_array.hpp $(INC_DIR)/pe.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/pe.o: $(INC_DIR)/pe.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/memory.o: $(INC_DIR)/memory.hpp $(INC_DIR)/config.hpp

# Header dependencies - 3D
$(BUILD_DIR)/config_3d.o: $(INC_DIR)/config_3d.hpp $(INC_DIR)/config.hpp
$(BUILD_DIR)/pe_3d.o: $(INC_DIR)/pe_3d.hpp $(INC_DIR)/pe.hpp
$(BUILD_DIR)/pe_array_3d.o: $(INC_DIR)/pe_array_3d.hpp $(INC_DIR)/pe_3d.hpp $(INC_DIR)/config_3d.hpp
$(BUILD_DIR)/block_lu_simulator_3d.o: $(INC_DIR)/block_lu_simulator_3d.hpp $(INC_DIR)/pe_array_3d.hpp $(INC_DIR)/memory.hpp $(INC_DIR)/config_3d.hpp
