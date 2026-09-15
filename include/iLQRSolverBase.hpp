#pragma once

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <cmath>
#include <vector>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::VectorXd)

class iLQRSolverBase
{
	protected:
		int Nx;
		int Nu;
		double h; //Time step
		int Nt;   //Horizon size
		double Tfinal;

		bool use_parallel_derivatives = false;

		//Cost matrixes
		Eigen::MatrixXd Q;
		Eigen::MatrixXd R;
		Eigen::MatrixXd Qn;

		Eigen::VectorXd x0;
		Eigen::VectorXd xgoal;
		std::vector<Eigen::VectorXd> xtraj; 
		std::vector<Eigen::VectorXd> utraj;
		std::vector<Eigen::VectorXd> xtraj_new; 
		std::vector<Eigen::VectorXd> utraj_new;

		std::vector<Eigen::VectorXd> d; 
		std::vector<Eigen::MatrixXd> K; 
		std::vector<Eigen::VectorXd> p; 
		std::vector<Eigen::MatrixXd> P; 

		Eigen::LDLT<Eigen::MatrixXd> ldlt;

		void initializeMemory();

		virtual Eigen::VectorXd continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u) = 0; // x' = f(x, u)
		virtual void getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju) = 0;

		Eigen::VectorXd discreteDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u);   // xk+1 = f(xk, uk)

		//Calculates A = df/dx and B = df/du; f being discreteDynamics
		void computeDerivatives(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& A, Eigen::MatrixXd& B);

		double stageCost(const Eigen::VectorXd& x, const Eigen::VectorXd& u);
		double terminalCost(const Eigen::VectorXd& x);
		double calculateTotalCost(const std::vector<Eigen::VectorXd>& x_traj, const std::vector<Eigen::VectorXd>& u_traj);

		bool backwardPass(double& delta_J, double& max_d);
		double forwardRollout(double alpha);

	public:
		iLQRSolverBase(double h, int Nt);
		virtual ~iLQRSolverBase() = default;

		void setCosts(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R, const Eigen::MatrixXd& Qn);
		void setBoundaryConditions(const Eigen::VectorXd& x0, const Eigen::VectorXd& xgoal);

		bool solve(int max_iters = 1000, double tol = 1e-3);

		int getNx() const { return Nx; }
		int getNu() const { return Nu; }
		std::vector<Eigen::VectorXd> getXTraj() const { return xtraj; }
		std::vector<Eigen::VectorXd> getUTraj() const { return utraj; }
		std::vector<Eigen::MatrixXd> getKTraj() const { return K; }
};
