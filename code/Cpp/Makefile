EXE=CFSModel

SRC_DIR = src
OBJ_DIR = obj
OUTPUT_DIR = outputs

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OUTPUT = $(wildcard $(OUTPUT_DIR)/*.txt)

CXX = g++
CXXFLAGS = -Wall -g -Iinclude

.Phony: all clean

all: $(EXE)

$(EXE): $(OBJ)
				$(CXX) $(CXXFLAGS) $^ -o $@
				@echo "\n ****** CFSModel program built! ******\n"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ) $(EXE) $(OUTPUT)
