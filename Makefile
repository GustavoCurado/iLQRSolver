CXX = g++
INCLUDE = -I ./include
SRC = ./src
OBJ = ./obj
BIN = ./bin
TEST = ./test

# Flags common to all executables
COMMON_FLAGS = -std=c++11 -O3 -fno-omit-frame-pointer -fopenmp $(INCLUDE)  #Maybe I could add -DNDBUG

# Flags for Pinocchio-based executables
PIN_CXXFLAGS = $(COMMON_FLAGS) $(shell pkg-config --cflags pinocchio) -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR -DPINOCCHIO_WITH_URDFDOM
PIN_LDFLAGS  = $(shell pkg-config --libs pinocchio) -lurdfdom_model -lboost_system -lpthread

# Flags for Optimized executables
OPT_CXXFLAGS = $(COMMON_FLAGS) $(shell pkg-config --cflags eigen3) -march=native -mavx
OPT_LDFLAGS  = -lpthread

# Targets
PIN_BINS = $(BIN)/acrobot $(BIN)/iiwa_14 $(BIN)/iiwa_14_mpc $(BIN)/iiwa_14_printless $(BIN)/benchmark_1 $(BIN)/benchmark_2
OPT_BINS = $(BIN)/opt_iiwa_14 $(BIN)/opt_iiwa_14_mpc $(BIN)/opt_iiwa_14_printless $(BIN)/opt_benchmark_1 $(BIN)/opt_benchmark_2

all: $(PIN_BINS) $(OPT_BINS)

$(OBJ)/iLQRSolverBase_pin.o: $(SRC)/iLQRSolverBase.cpp
	$(CXX) $(PIN_CXXFLAGS) -c $< -o $@

$(OBJ)/iLQRSolverPinocchio.o: $(SRC)/iLQRSolverPinocchio.cpp
	$(CXX) $(PIN_CXXFLAGS) -c $< -o $@

PIN_OBJS = $(OBJ)/iLQRSolverBase_pin.o $(OBJ)/iLQRSolverPinocchio.o

$(OBJ)/iLQRSolverBase_opt.o: $(SRC)/iLQRSolverBase.cpp
	$(CXX) $(OPT_CXXFLAGS) -c $< -o $@

$(OBJ)/iLQRSolverOptimized.o: $(SRC)/iLQRSolverOptimized.cpp
	$(CXX) $(OPT_CXXFLAGS) -c $< -o $@

OPT_OBJS = $(OBJ)/iLQRSolverBase_opt.o $(OBJ)/iLQRSolverOptimized.o

$(BIN)/acrobot: $(TEST)/acrobot.cpp $(PIN_OBJS)
	$(CXX) $(PIN_CXXFLAGS) $< $(PIN_OBJS) -o $@ $(PIN_LDFLAGS)

$(BIN)/iiwa_14: $(TEST)/iiwa_14.cpp $(PIN_OBJS)
	$(CXX) $(PIN_CXXFLAGS) $< $(PIN_OBJS) -o $@ $(PIN_LDFLAGS)

$(BIN)/iiwa_14_mpc: $(TEST)/iiwa_14_mpc.cpp $(PIN_OBJS)
	$(CXX) $(PIN_CXXFLAGS) $< $(PIN_OBJS) -o $@ $(PIN_LDFLAGS)

$(BIN)/iiwa_14_printless: $(TEST)/iiwa_14_printless.cpp $(PIN_OBJS)
	$(CXX) $(PIN_CXXFLAGS) $< $(PIN_OBJS) -o $@ $(PIN_LDFLAGS)

$(BIN)/opt_iiwa_14: $(TEST)/opt_iiwa_14.cpp $(OPT_OBJS)
	$(CXX) $(OPT_CXXFLAGS) $< $(OPT_OBJS) -o $@ $(OPT_LDFLAGS)

$(BIN)/opt_iiwa_14_mpc: $(TEST)/opt_iiwa_14_mpc.cpp $(OPT_OBJS)
	$(CXX) $(OPT_CXXFLAGS) $< $(OPT_OBJS) -o $@ $(OPT_LDFLAGS)

$(BIN)/opt_iiwa_14_printless: $(TEST)/opt_iiwa_14_printless.cpp $(OPT_OBJS)
	$(CXX) $(OPT_CXXFLAGS) $< $(OPT_OBJS) -o $@ $(OPT_LDFLAGS)

$(BIN)/benchmark_1: $(TEST)/benchmark_1.cpp $(PIN_OBJS)
	$(CXX) $(PIN_CXXFLAGS) $< $(PIN_OBJS) -o $@ $(PIN_LDFLAGS)

$(BIN)/opt_benchmark_1: $(TEST)/benchmark_1.cpp $(OPT_OBJS)
	$(CXX) $(OPT_CXXFLAGS) -DOPTIMIZED_BACKEND $< $(OPT_OBJS) -o $@ $(OPT_LDFLAGS)

$(BIN)/benchmark_2: $(TEST)/benchmark_2.cpp $(PIN_OBJS)
	$(CXX) $(PIN_CXXFLAGS) $< $(PIN_OBJS) -o $@ $(PIN_LDFLAGS)

$(BIN)/opt_benchmark_2: $(TEST)/benchmark_2.cpp $(OPT_OBJS)
	$(CXX) $(OPT_CXXFLAGS) -DOPTIMIZED_BACKEND $< $(OPT_OBJS) -o $@ $(OPT_LDFLAGS)

clean:
	rm -f $(OBJ)/*.o $(BIN)/*

.PHONY: all clean
