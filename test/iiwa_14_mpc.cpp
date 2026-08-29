#include "iLQRSolver.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

int main()
{
	std::string urdf_filename = "utils/iiwa_14.urdf";
	
	double h = 0.01;         
	double Tfinal = 1.0;     
	int Nt = (int)(Tfinal/h) + 1;

	iLQRSolver ilqr(urdf_filename, h, Nt);

	int Nx = ilqr.getNx();
	int Nu = ilqr.getNu();

	Eigen::Matrix<double, 14, 1> x0;
	x0 << 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
	      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	
	Eigen::Matrix<double, 14, 1> xgoal;
	xgoal << 2.00, 1.20, -0.30, -1.50, 0.40, 1.00, 0.50,
	         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	Eigen::Matrix<double, 14, 1> Q_diag;                                   
	Q_diag << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
	Eigen::MatrixXd Q = Q_diag.asDiagonal();
	Eigen::MatrixXd R = 0.1 * Eigen::MatrixXd::Identity(Nu, Nu);
	Eigen::MatrixXd Qn = 10000.0 * Eigen::MatrixXd::Identity(Nx, Nx);

	ilqr.setCosts(Q, R, Qn);
	ilqr.setGravity(0.0);  //No gravity

	ilqr.setBoundaryConditions(x0, xgoal);
	ilqr.solve(100, 1e-4);

	std::vector<Eigen::VectorXd> traj_real_x;
	std::vector<Eigen::VectorXd> traj_real_u;
	
	Eigen::VectorXd estado_atual = x0;
	traj_real_x.push_back(estado_atual);

	double tempo_simulacao = 5.0;
	int mpc_steps = (int)(tempo_simulacao / h);

	std::cout << "Rodando simulacao MPC..." << std::endl;

	for (int i = 0; i < mpc_steps; ++i)
	{
		ilqr.setBoundaryConditions(estado_atual, xgoal);
		
		ilqr.solve(1);

		auto utraj = ilqr.getUTraj();
		auto xtraj = ilqr.getXTraj();

		Eigen::VectorXd controle = utraj[0];
		estado_atual = xtraj[1]; 

		traj_real_u.push_back(controle);
		traj_real_x.push_back(estado_atual);
	}

	std::ofstream arquivo("trajetoria_iiwa.csv");

	for (int i = 1; i <= Nx; ++i) arquivo << "x" << i << ",";
	for (int i = 1; i < Nu; ++i) arquivo << "u" << i << ",";
	arquivo << "u" << Nu << "\n";

	for (size_t i = 0; i < traj_real_u.size(); ++i)
	{    
		for (int j = 0; j < Nx; ++j) arquivo << traj_real_x[i][j] << ",";
		for (int j = 0; j < Nu - 1; ++j) arquivo << traj_real_u[i][j] << ",";
		arquivo << traj_real_u[i][Nu - 1] << "\n";
	}

	arquivo.close();
	std::cout << "Trajetoria MPC salva com sucesso em 'trajetoria_iiwa.csv'!" << std::endl;

	return 0;
}
