#pragma once

#include "iLQRSolverBase.hpp"

class iLQRSolverGPU : public iLQRSolverBase
{
	private:
		void initializeMemory();

		Eigen::VectorXd continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u) override;
		void getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju) override;

	public:
		iLQRSolverGPU(double h, int Nt);
};
