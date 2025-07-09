# Compiler
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11

# Target executable name
TARGET = 3DLUD_model

# Source files
SRC = 3DLUD_model.cpp

# Build target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean up build artifacts
clean:
	rm -f $(TARGET)