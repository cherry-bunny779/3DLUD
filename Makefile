# Compiler
CXX = g++
CC  = gcc

CXXFLAGS = -O2 -Wall -std=c++11
CFLAGS   = -O2 -Wall

# Target executable name
TARGET = 3DLUD_model

# Source files
CPP_SRC = 3DLUD_model.cpp LU_decomp_blocked.cpp

# Object files
OBJ = $(CPP_SRC:.cpp=.o)

# Build target
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Compile C++ source
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(TARGET) *.o

.PHONY: all clean