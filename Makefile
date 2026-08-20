CXX = g++
INCLUDE = -I ./src/include
SRC = ./src
OBJ = ./obj
BIN = ./bin
TEST = ./test

CXXFLAGS = -std=c++11 -O3 $(shell pkg-config --cflags pinocchio) -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR -DPINOCCHIO_WITH_URDFDOM -lurdfdom_model -lboost_system -lpthread $(INCLUDE)

LDFLAGS = $(shell pkg-config --libs pinocchio)

all: $(BIN)/acrobot $(BIN)/iiwa_14

$(OBJ)/iLQRSolver.o: $(SRC)/iLQRSolver.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/iLQRSolver.cpp -o $(OBJ)/iLQRSolver.o

$(BIN)/acrobot: $(TEST)/acrobot.cpp $(OBJ)/iLQRSolver.o
	$(CXX) $(CXXFLAGS) $(TEST)/acrobot.cpp $(OBJ)/iLQRSolver.o -o $(BIN)/acrobot $(LDFLAGS)

$(BIN)/iiwa_14: $(TEST)/iiwa_14.cpp $(OBJ)/iLQRSolver.o
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14.cpp $(OBJ)/iLQRSolver.o -o $(BIN)/iiwa_14 $(LDFLAGS)
