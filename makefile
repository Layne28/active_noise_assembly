#Some organization borrowed from: https://github.com/tscott8706/cpp-csv-col-replacer/blob/master

# Build executable with:
# % make
# Delete object files and executable with:
# % make clean
# Rebuild all objects and executable with:
# % make -B

SRC_DIR := src
OBJ_DIR := build
BIN_DIR := bin
TEST_SRC_DIR := test/src
TEST_OBJ_DIR := test/build

EXECUTABLE := assembly
TEST_EXECUTABLE := test_assembly
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
TEST_SOURCES := $(wildcard $(TEST_SRC_DIR)/*.cpp)
HEADERS := $(wildcard $(SRC_DIR)/*.hpp) $(wildcard ../Common-Tools/*.h)
TEST_HEADERS := $(wildcard $(TEST_SRC_DIR)/*.hpp)

#CXX := g++
CXX := h5c++

SHELL = /bin/sh

# Flags to pass to the compiler; per the recommendations of the GNU Scientific Library
CXXFLAGS:= -std=c++17 -Wextra -pedantic -Wall -W -Wmissing-declarations -Wuninitialized -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -fshort-enums -fno-common -m64 -fopenmp -I$(HOME)/.local/include -I../Common-Tools

# Compiler flags controling optimization levels. Use -O3 for full optimization,
# but make sure your results are consistent
# -g includes debugging information. You can also add -pg here for profiling 
PROFILE=-g -pg
OPTFLAGS:=$(PROFILE) -O2

# Flags to pass to the linker; -lm links in the standard c math library
#LDFLAGS:= -fopenmp -lm -lgsl -lgslcblas -llapack -lblas -larmadillo -langen -lfftw3 $(PROFILE) -L$(HOME)/.local/lib 
LDFLAGS:= -fopenmp -lm -lgsl -lgslcblas -lopenblas -larmadillo -lstdc++fs -langen -lfftw3 -lhdf5 -lhdf5_cpp $(PROFILE) -L$(HOME)/.local/lib 

# Variable to compose names of object files from the names of sources
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
OBJECTS_NO_MAIN = $(filter-out $(OBJ_DIR)/main.o,$(OBJECTS))

#When compiling tests, include all objects in actual program except for main
#(there's a main function in the test folder)
TEST_OBJECTS := $(patsubst $(TEST_SRC_DIR)/%.cpp,$(TEST_OBJ_DIR)/%.o,$(TEST_SOURCES))
TEST_OBJECTS += $(OBJECTS_NO_MAIN)

# Default target depends on sources and headers to detect changes
all: $(SOURCES) $(HEADERS)  $(BIN_DIR)/$(EXECUTABLE)
install: 
	install bin/* $(HOME)/.local/bin/
test: $(TEST_SOURCES) $(TEST_HEADERS) $(BIN_DIR)/$(TEST_EXECUTABLE)

# Rule to compile a source file to object code
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(OPTFLAGS) $< -o $@
#$(TEST_OBJECTS): $(TEST_OBJ_DIR)/%.o : $(TEST_SRC_DIR)/%.cpp
$(TEST_OBJ_DIR)/%.o : $(TEST_SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

# Build the executable by linking all objects
$(BIN_DIR)/$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@
$(BIN_DIR)/$(TEST_EXECUTABLE): $(TEST_OBJECTS)
	$(CXX) $(TEST_OBJECTS) $(LDFLAGS) -o $@

# clean up so we can start over (removes executable!)
clean:
	rm -f $(OBJ_DIR)/*.o $(TEST_OBJ_DIR)/*.o $(BIN_DIR)/$(EXECUTABLE) $(BIN_DIR)/$(TEST_EXECUTABLE)
