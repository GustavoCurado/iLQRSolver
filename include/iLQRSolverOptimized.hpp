#pragma once

#include "iLQRSolverBase.hpp"

class iLQRSolverOptimized : public iLQRSolverBase
{
	private:
		// Arrays fixos para o robô Kuka LBR iiwa 14 (7 juntas) exigidos pelo código do Plancher
		double d_T[36 * 7]; // Transformações
		double d_I[36 * 7]; // Inércias

		void initializeMemory();

		Eigen::VectorXd continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u) override;
		void getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju) override;

	public:
		iLQRSolverOptimized(double h, int Nt);
		void setGravity(double g_z);
};
