CXX = g++ -g3 -fno-omit-frame-pointer
INCLUDE = -I ./src/include
SRC = ./src
OBJ = ./obj
BIN = ./bin
TEST = ./test

CXXFLAGS = -std=c++11 -O3 $(shell pkg-config --cflags pinocchio) -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR -DPINOCCHIO_WITH_URDFDOM -lurdfdom_model -lboost_system -lpthread $(INCLUDE)

LDFLAGS = $(shell pkg-config --libs pinocchio)

all: $(BIN)/acrobot $(BIN)/iiwa_14 $(BIN)/iiwa_14_fixed $(BIN)/iiwa_14_printless $(BIN)/iiwa_14_fixed_printless

$(OBJ)/iLQRSolver.o: $(SRC)/iLQRSolver.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/iLQRSolver.cpp -o $(OBJ)/iLQRSolver.o

$(OBJ)/iLQRSolverFixed.o: $(SRC)/iLQRSolverFixed.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/iLQRSolverFixed.cpp -o $(OBJ)/iLQRSolverFixed.o

$(BIN)/acrobot: $(TEST)/acrobot.cpp $(OBJ)/iLQRSolver.o
	$(CXX) $(CXXFLAGS) $(TEST)/acrobot.cpp $(OBJ)/iLQRSolver.o -o $(BIN)/acrobot $(LDFLAGS)

$(BIN)/iiwa_14: $(TEST)/iiwa_14.cpp $(OBJ)/iLQRSolver.o
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14.cpp $(OBJ)/iLQRSolver.o -o $(BIN)/iiwa_14 $(LDFLAGS)

$(BIN)/iiwa_14_fixed: $(TEST)/iiwa_14_fixed.cpp $(OBJ)/iLQRSolverFixed.o
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14_fixed.cpp $(OBJ)/iLQRSolverFixed.o -o $(BIN)/iiwa_14_fixed $(LDFLAGS)

$(BIN)/iiwa_14_printless: $(TEST)/iiwa_14_printless.cpp $(OBJ)/iLQRSolver.o
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14_printless.cpp $(OBJ)/iLQRSolver.o -o $(BIN)/iiwa_14_printless $(LDFLAGS)

$(BIN)/iiwa_14_fixed_printless: $(TEST)/iiwa_14_fixed_printless.cpp $(OBJ)/iLQRSolverFixed.o
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14_fixed_printless.cpp $(OBJ)/iLQRSolverFixed.o -o $(BIN)/iiwa_14_fixed_printless $(LDFLAGS)
