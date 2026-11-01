#include "iLQRSolverOptimized.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>

Eigen::Matrix<double, 14, 1> getMovingGoal(double t)
{
	Eigen::Matrix<double, 14, 1> goal = Eigen::Matrix<double, 14, 1>::Zero();

	goal[0] =  ( 0.6 * std::sin(0.2 * M_PI * t));
	goal[1] =  ( 0.5 * std::sin(0.4 * M_PI * t) + 0.3);
	goal[2] =  (-0.4 * std::sin(0.2 * M_PI * t));
	goal[3] =  (-0.8 * std::sin(0.4 * M_PI * t) - 0.5);
	goal[4] =  ( 0.4 * std::sin(0.2 * M_PI * t));
	goal[5] =  ( 0.6 * std::sin(0.4 * M_PI * t) + 0.2);
	goal[6] =  ( 0.5 * std::sin(0.2 * M_PI * t));

	return goal;
}

int main()
{
	double h = 0.01;
	int Nt = 64;

	iLQRSolverOptimized ilqr(h, Nt);
	ilqr.setParallelDerivatives(true);

	int Nx = ilqr.getNx();
	int Nu = ilqr.getNu();

	Eigen::Matrix<double, 14, 1> x0 = Eigen::Matrix<double, 14, 1>::Zero();
	Eigen::Matrix<double, 14, 1> xgoal = getMovingGoal(0.0);

	Eigen::Matrix<double, 14, 1> Q_diag, Qn_diag;
	Q_diag << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
	          0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
	Eigen::MatrixXd Q = Q_diag.asDiagonal();
	Eigen::MatrixXd R = 0.0001 * Eigen::MatrixXd::Identity(Nu, Nu);
	Qn_diag << 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
		   10.0,   10.0,   10.0,   10.0,   10.0,   10.0,   10.0;
	Eigen::MatrixXd Qn = Qn_diag.asDiagonal();

	ilqr.setCosts(Q, R, Qn);
	ilqr.setGravity(0.0); //No gravity
	ilqr.setBoundaryConditions(x0, xgoal);

	std::vector<Eigen::VectorXd> traj_real_x;
	std::vector<Eigen::VectorXd> traj_real_u;
	std::vector<Eigen::MatrixXd> traj_real_K;

	Eigen::VectorXd estado_atual = x0;
	traj_real_x.push_back(estado_atual);

	double tempo_simulacao = 10.0;
	double time_limit_ms = 2.0;
	int mpc_steps = (int)(tempo_simulacao / h);

	std::cout << "Rodando simulacao MPC com " << mpc_steps << " passos..." << std::endl;

	double tempo_total_mpc = 0.0;
	int total_iters = 0;
	double erro_medio_pos = 0.0;

	for (int i = 0; i < mpc_steps; ++i)
	{
		double t_atual = i * h;
		xgoal = getMovingGoal(t_atual);
		ilqr.setGoal(xgoal);

		if (i > 0)
			ilqr.advanceMPCStep(estado_atual);

		auto inicio_mpc = std::chrono::high_resolution_clock::now();
		ilqr.solve(100, 1e-4, time_limit_ms);
		auto fim_mpc = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double, std::milli> tempo_passo = fim_mpc - inicio_mpc;
		tempo_total_mpc += tempo_passo.count();
		total_iters += ilqr.getLastIters();

		auto utraj = ilqr.getUTraj();
		auto Ktraj = ilqr.getKTraj();

		Eigen::VectorXd controle = utraj[0];
		Eigen::MatrixXd ganho = Ktraj[0];

		estado_atual = ilqr.simulatePlant(estado_atual, h, 100);

		erro_medio_pos += (xgoal.head(7) - estado_atual.head(7)).norm();

		traj_real_u.push_back(controle);
		traj_real_K.push_back(ganho);
		traj_real_x.push_back(estado_atual);
	}

	std::cout << "\nEstado final alcancado:\n" << traj_real_x.back().transpose() << std::endl;
	std::cout << "\nAlvo final:\n" << xgoal.transpose() << std::endl;
	std::cout << "\nErro final |xgoal - xtraj|:\n" << (xgoal.transpose() - traj_real_x.back().transpose()).cwiseAbs() << std::endl << std::endl;
	std::cout << "Media de iteracoes por passo do MPC: " << (double)total_iters / mpc_steps << std::endl;
	std::cout << "Tempo medio por passo do MPC: " << tempo_total_mpc / mpc_steps << " ms" << std::endl;

	std::ofstream arquivo("trajetoria_iiwa.csv");

	for (int i = 1; i <= Nx; ++i) arquivo << "x" << i << ",";
	for (int i = 1; i <= Nu; ++i) arquivo << "u" << i << ",";
	for (int i = 1; i <= Nu * Nx; ++i) arquivo << (i < Nu * Nx ? "K" + std::to_string(i) + "," : "K" + std::to_string(i) + "\n");

	for (size_t i = 0; i < traj_real_u.size(); ++i)
	{
		for (int j = 0; j < Nx; ++j) arquivo << traj_real_x[i][j] << ",";

		for (int j = 0; j < Nu; ++j) arquivo << traj_real_u[i][j] << ",";

		for (int row = 0; row < Nu; ++row)
		{
			for (int col = 0; col < Nx; ++col)
			{
				arquivo << traj_real_K[i](row, col);
				if (row == Nu - 1 && col == Nx - 1) arquivo << "\n";
				else arquivo << ",";
			}
		}
	}

	for (int j = 0; j < Nx; ++j) arquivo << traj_real_x.back()[j] << ",";
	for (int j = 0; j < Nu; ++j) arquivo << "0.0,";
	for (int i = 0; i < Nu * Nx - 1; ++i) arquivo << "0.0,";
	arquivo << "0.0\n";

	arquivo.close();
	std::cout << "Trajetoria MPC salva com sucesso em 'trajetoria_iiwa.csv'!" << std::endl;

	return 0;
}
