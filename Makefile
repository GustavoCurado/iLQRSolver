CXX = g++ -g3 -fno-omit-frame-pointer -fopenmp
INCLUDE = -I ./include
SRC = ./src
OBJ = ./obj
BIN = ./bin
TEST = ./test

CXXFLAGS = -std=c++11 -O3 $(shell pkg-config --cflags pinocchio) -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR -DPINOCCHIO_WITH_URDFDOM -lurdfdom_model -lboost_system -lpthread $(INCLUDE)

LDFLAGS = $(shell pkg-config --libs pinocchio)

all: $(BIN)/acrobot $(BIN)/iiwa_14 $(BIN)/iiwa_14_mpc $(BIN)/iiwa_14_printless $(BIN)/opt_iiwa_14 $(BIN)/opt_iiwa_14_printless $(BIN)/opt_iiwa_14_mpc

$(OBJ)/iLQRSolverBase.o: $(SRC)/iLQRSolverBase.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/iLQRSolverBase.cpp -o $(OBJ)/iLQRSolverBase.o

$(OBJ)/iLQRSolverPinocchio.o: $(SRC)/iLQRSolverPinocchio.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/iLQRSolverPinocchio.cpp -o $(OBJ)/iLQRSolverPinocchio.o

$(OBJ)/iLQRSolverOptimized.o: $(SRC)/iLQRSolverOptimized.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/iLQRSolverOptimized.cpp -o $(OBJ)/iLQRSolverOptimized.o

SOLVER_OBJS = $(OBJ)/iLQRSolverBase.o $(OBJ)/iLQRSolverPinocchio.o $(OBJ)/iLQRSolverOptimized.o

$(BIN)/acrobot: $(TEST)/acrobot.cpp $(SOLVER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST)/acrobot.cpp $(SOLVER_OBJS) -o $(BIN)/acrobot $(LDFLAGS)

$(BIN)/iiwa_14: $(TEST)/iiwa_14.cpp $(SOLVER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14.cpp $(SOLVER_OBJS) -o $(BIN)/iiwa_14 $(LDFLAGS)

$(BIN)/iiwa_14_mpc: $(TEST)/iiwa_14_mpc.cpp $(SOLVER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14_mpc.cpp $(SOLVER_OBJS) -o $(BIN)/iiwa_14_mpc $(LDFLAGS)

$(BIN)/iiwa_14_printless: $(TEST)/iiwa_14_printless.cpp $(SOLVER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST)/iiwa_14_printless.cpp $(SOLVER_OBJS) -o $(BIN)/iiwa_14_printless $(LDFLAGS)

$(BIN)/opt_iiwa_14: $(TEST)/opt_iiwa_14.cpp $(SOLVER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST)/opt_iiwa_14.cpp $(SOLVER_OBJS) -o $(BIN)/opt_iiwa_14 $(LDFLAGS)

$(BIN)/opt_iiwa_14_printless: $(TEST)/opt_iiwa_14_printless.cpp $(SOLVER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST)/opt_iiwa_14_printless.cpp $(SOLVER_OBJS) -o $(BIN)/opt_iiwa_14_printless $(LDFLAGS)

$(BIN)/opt_iiwa_14_mpc: $(TEST)/opt_iiwa_14_mpc.cpp $(SOLVER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST)/opt_iiwa_14_mpc.cpp $(SOLVER_OBJS) -o $(BIN)/opt_iiwa_14_mpc $(LDFLAGS)
