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
#include <algorithm>
#include <numeric>
#include <cmath>


void printBreakdown(const std::string& title, int iters_count, double avg_deriv_us, double avg_riccati_us, double avg_rollout_us, double avg_total_us)
{
	double pct_deriv   = 100.0 * avg_deriv_us / avg_total_us;
	double pct_riccati = 100.0 * avg_riccati_us / avg_total_us;
	double pct_rollout = 100.0 * avg_rollout_us / avg_total_us;

	printf("\n--- %s (Repetitions: %d) ---\n", title.c_str(), iters_count);
	printf("  1. Derivatives                : %8.2f ms  (%5.1f%%)\n", avg_deriv_us * 1e-3, pct_deriv);
	printf("  2. Backward Pass              : %8.2f ms  (%5.1f%%)\n", avg_riccati_us * 1e-3, pct_riccati);
	printf("  3. Forward Rollout            : %8.2f ms  (%5.1f%%)\n", avg_rollout_us * 1e-3, pct_rollout);
	printf("  ----------------------------------------------------\n");
	printf("  TOTAL solve()                 : %8.2f ms (100.0%%)\n", avg_total_us * 1e-3);
}


void runFullTrajectoryBenchmark(bool parallel_derivatives, int num_trials = 30)
{
	double h = 0.01;
	double Tfinal = 5.0;
	int Nt = (int)(Tfinal / h) + 1;

	#ifdef OPTIMIZED_BACKEND
		iLQRSolverOptimized ilqr(h, Nt);
	#else
		iLQRSolverPinocchio ilqr("utils/iiwa_14.urdf", h, Nt);
	#endif

	ilqr.setParallelDerivatives(parallel_derivatives);

	int Nx = ilqr.getNx();
	int Nu = ilqr.getNu();

	Eigen::Matrix<double, 14, 1> x0 = Eigen::Matrix<double, 14, 1>::Zero();
	Eigen::Matrix<double, 14, 1> xgoal;
	xgoal << 2.00, 1.20, -0.30, -1.50, 0.40, 1.00, 0.50,
	         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	Eigen::Matrix<double, 14, 1> Q_diag;
	Q_diag << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
	          0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001;

	Eigen::MatrixXd Q = Q_diag.asDiagonal();
	Eigen::MatrixXd R = 0.001 * Eigen::MatrixXd::Identity(Nu, Nu);
	Eigen::MatrixXd Qn = 1000.0 * Eigen::MatrixXd::Identity(Nx, Nx);

	ilqr.setCosts(Q, R, Qn);
	ilqr.setBoundaryConditions(x0, xgoal);
	ilqr.setGravity(0.0);

	double sum_deriv = 0.0, sum_riccati = 0.0, sum_rollout = 0.0, sum_total = 0.0;
	int iters = 0;

	for (int t = 0; t < num_trials; t++)
	{
		ilqr.resetTrajectories();
		ilqr.solve(10000);

		sum_deriv   += ilqr.getTimeDerivativesUs();
		sum_riccati += ilqr.getTimeRiccatiUs();
		sum_rollout += ilqr.getTimeRolloutUs();
		sum_total   += ilqr.getTimeTotalUs();
		iters        = ilqr.getLastIters();
	}

	std::string mode = parallel_derivatives ? "Multithread" : "Single-thread";
	printBreakdown(BACKEND_NAME + " - Full Trajectory [Nt=501] - " + mode,
	               num_trials,
	               sum_deriv / num_trials,
	               sum_riccati / num_trials,
	               sum_rollout / num_trials,
	               sum_total / num_trials);
}


int main()
{
	std::cout << "====================================================================\n";
	std::cout << "         END-TO-END iLQR BREAKDOWN (" << BACKEND_NAME << ")\n";
	std::cout << "====================================================================\n";

	runFullTrajectoryBenchmark(false, 30);
	runFullTrajectoryBenchmark(true, 30);

	return 0;
}
