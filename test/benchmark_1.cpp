#ifdef OPTIMIZED_BACKEND
	#include "iLQRSolverOptimized.hpp"
	using SolverType = iLQRSolverOptimized;
	const std::string BACKEND_NAME = "Optimized";
#else
	#include "iLQRSolverPinocchio.hpp"
	using SolverType = iLQRSolverPinocchio;
	const std::string BACKEND_NAME = "Pinocchio";
#endif

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <time.h>
#include <omp.h>


// Macro from Plancher's code that converts timespec into microsseconds
#define time_delta_us(start, end) \
    (1e6 * static_cast<double>(end.tv_sec - start.tv_sec) + \
     1e-3 * static_cast<double>(end.tv_nsec - start.tv_nsec))


void printStats(const std::string& label, std::vector<double>& times)
{
	std::sort(times.begin(), times.end());

	size_t n = times.size();
	double sum = std::accumulate(times.begin(), times.end(), 0.0);
	double mean = sum / n;

	double sq_sum = 0.0;
	for (double t : times) 
	{
		sq_sum += (t - mean) * (t - mean);
	}

	double stdev = std::sqrt(sq_sum / n);
	double min_val = times.front();
	double median  = times[n / 2];
	double p99     = times[(n * 99) / 100];

	printf("%-38s | Mediana: %8.2f us | Media: %8.2f us | Std: %6.2f us | Min: %8.2f us | P99: %8.2f us\n", label.c_str(), median, mean, stdev, min_val, p99);
}


void runSingleComputationTest(SolverType& solver, int single_iters = 100000)
{
	std::cout << "\n========================================================================================\n";
	std::cout << "  [" << BACKEND_NAME << "] SINGLE COMPUTATION LATENCY (N = 1, Average of " << single_iters << " iterations)\n";
	std::cout << "========================================================================================\n";

	int Nx = 14, Nu = 7;

	// Generates a random (x, u) state
	std::default_random_engine rng(42);
	std::normal_distribution<double> dist(0.0, 1.0);

	Eigen::VectorXd x(Nx), u(Nu);
	for (int i = 0; i < Nx; i++) x[i] = dist(rng);
	for (int i = 0; i < Nu; i++) u[i] = dist(rng);

	Eigen::VectorXd xdot(Nx);
	Eigen::MatrixXd Jx(Nx, Nx), Ju(Nx, Nu);
	Eigen::MatrixXd A(Nx, Nx),  B(Nx, Nu);

	struct timespec start, end;

	// Testing continuousDynamics()
	clock_gettime(CLOCK_MONOTONIC, &start);
	for (int i = 0; i < single_iters; i++)
	{
		xdot = solver.continuousDynamics(x, u);
	}

	clock_gettime(CLOCK_MONOTONIC, &end);
	double time_dyn = time_delta_us(start, end) / single_iters;

	// Testing getContinuousJacs()
	clock_gettime(CLOCK_MONOTONIC, &start);
	for (int i = 0; i < single_iters; i++)
	{
		solver.getContinuousJacs(x, u, xdot, Jx, Ju);
	}

	clock_gettime(CLOCK_MONOTONIC, &end);
	double time_jac = time_delta_us(start, end) / single_iters;

	// Testing computeDerivatives()
	clock_gettime(CLOCK_MONOTONIC, &start);
	for (int i = 0; i < single_iters; i++)
	{
		solver.computeDerivatives(x, u, A, B);
	}

	clock_gettime(CLOCK_MONOTONIC, &end);
	double time_compute = time_delta_us(start, end) / single_iters;


	printf("continuousDynamics (1 call) | %s: %6.2f us\n", BACKEND_NAME.c_str(), time_dyn);
	printf("getContinuousJacs  (1 call) | %s: %6.2f us\n", BACKEND_NAME.c_str(), time_jac);
	printf("computeDerivatives (1 call) | %s: %6.2f us\n", BACKEND_NAME.c_str(), time_compute);
}


void runHorizonBatchTest(SolverType& solver, int N_points, int test_iters = 5000)
{
	std::cout << "\n========================================================================================\n";
	std::cout << "         [" << BACKEND_NAME << "] HORIZON N = " << N_points << " POINTS (" << test_iters << " repetitions)\n";
	std::cout << "========================================================================================\n";

	int Nx = 14, Nu = 7;

	// Generates random trajectory of size N
	std::default_random_engine rng(123 + N_points);
	std::normal_distribution<double> dist(0.0, 1.0);
	std::vector<Eigen::VectorXd> xs(N_points, Eigen::VectorXd(Nx));
	std::vector<Eigen::VectorXd> us(N_points, Eigen::VectorXd(Nu));

	for (int k = 0; k < N_points; k++)
	{
		for (int i = 0; i < Nx; i++) xs[k][i] = dist(rng);
		for (int i = 0; i < Nu; i++) us[k][i] = dist(rng);
	}

	std::vector<Eigen::MatrixXd> Atraj(N_points, Eigen::MatrixXd::Zero(Nx, Nx));
	std::vector<Eigen::MatrixXd> Btraj(N_points, Eigen::MatrixXd::Zero(Nx, Nu));

	struct timespec start, end;

	std::vector<double> times_1t(test_iters);
	std::vector<double> times_mt(test_iters);

	// 1 thread
	for (int i = 0; i < test_iters; i++)
	{
		clock_gettime(CLOCK_MONOTONIC, &start);
		for (int k = 0; k < N_points; k++)
		{
			solver.computeDerivatives(xs[k], us[k], Atraj[k], Btraj[k]);
		}

		clock_gettime(CLOCK_MONOTONIC, &end);
		times_1t[i] = time_delta_us(start, end);
	}

	// Multithread
	for (int i = 0; i < test_iters; i++)
	{
		clock_gettime(CLOCK_MONOTONIC, &start);

		#pragma omp parallel for num_threads(N_THREADS)
		for (int k = 0; k < N_points; k++)
		{
			solver.computeDerivatives(xs[k], us[k], Atraj[k], Btraj[k]);
		}

		clock_gettime(CLOCK_MONOTONIC, &end);
		times_mt[i] = time_delta_us(start, end);
	}

	printStats(BACKEND_NAME + " (1 thread)", times_1t);
	printStats(BACKEND_NAME + " (multithread)", times_mt);

	double med_1t = times_1t[test_iters / 2];
	double med_mt = times_mt[test_iters / 2];

	printf("  -> Speedup (MT vs 1T): %.2fx\n", med_1t / med_mt);
}

int main()
{
	int Nt = 64;
	double h = 0.01;

	#ifdef OPTIMIZED_BACKEND
		iLQRSolverOptimized solver(h, Nt);
	#else
		iLQRSolverPinocchio solver("utils/iiwa_14.urdf", h, Nt);
	#endif

	solver.setGravity(0.0);

	// Single call latency test
	runSingleComputationTest(solver, 100000);

	// Size N horizon latency test
	for (int N : {16, 32, 64, 128, 500})
	{
		int iters = (N >= 500) ? 1000 : 5000;
		runHorizonBatchTest(solver, N, iters);
	}

	return 0;
}
