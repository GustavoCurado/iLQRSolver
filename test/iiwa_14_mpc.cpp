#include "iLQRSolverPinocchio.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>

int main()
{
	std::string urdf_filename = "utils/iiwa_14.urdf";
	
	double h = 0.01;         
	double Tfinal = 0.5;     
	//int Nt = (int)(Tfinal/h) + 1;
	int Nt = 64;

	iLQRSolverPinocchio ilqr(urdf_filename, h, Nt);

	int Nx = ilqr.getNx();
	int Nu = ilqr.getNu();

	Eigen::Matrix<double, 14, 1> x0;
	x0 << 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
	      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	
	Eigen::Matrix<double, 14, 1> xgoal;
	xgoal << 2.00, 1.20, -0.30, -1.50, 0.40, 1.00, 0.50,
	         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	Eigen::Matrix<double, 14, 1> Q_diag, Qn_diag;                                   
	Q_diag << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
		  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
	Eigen::MatrixXd Q = Q_diag.asDiagonal();
	Eigen::MatrixXd R = 0.01 * Eigen::MatrixXd::Identity(Nu, Nu);
	Qn_diag << 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
		   10.0,   10.0,   10.0,   10.0,   10.0,   10.0,   10.0;
	Eigen::MatrixXd Qn = Qn_diag.asDiagonal();

	ilqr.setCosts(Q, R, Qn);
	ilqr.setGravity(0.0);  //No gravity

	ilqr.setBoundaryConditions(x0, xgoal);

//	auto inicio_cold = std::chrono::high_resolution_clock::now();
//	ilqr.solve(100, 1e-4);
//	auto fim_cold = std::chrono::high_resolution_clock::now();
	
//	std::chrono::duration<double, std::milli> tempo_cold = fim_cold - inicio_cold;

	std::vector<Eigen::VectorXd> traj_real_x;
	std::vector<Eigen::VectorXd> traj_real_u;
	std::vector<Eigen::MatrixXd> traj_real_K;
	
	Eigen::VectorXd estado_atual = x0;
	traj_real_x.push_back(estado_atual);

	double tempo_simulacao = 5.0;
	int mpc_steps = (int)(tempo_simulacao / h);

	std::cout << "Rodando simulacao MPC com " << mpc_steps << " passos..." << std::endl;

	double tempo_total_mpc = 0.0;

	for (int i = 0; i < mpc_steps; ++i)
	{
		if (i > 0) {
			ilqr.advanceMPCStep(estado_atual);
		}
		
		auto inicio_rti = std::chrono::high_resolution_clock::now();
		ilqr.solve(1);
		auto fim_rti = std::chrono::high_resolution_clock::now();
		
		std::chrono::duration<double, std::milli> tempo_rti = fim_rti - inicio_rti;
		tempo_total_mpc += tempo_rti.count();

		auto utraj = ilqr.getUTraj();
		auto xtraj = ilqr.getXTraj();
		auto Ktraj = ilqr.getKTraj();

		Eigen::VectorXd controle = utraj[0];
		Eigen::MatrixXd ganho = Ktraj[0];
		estado_atual = xtraj[1]; 

		traj_real_u.push_back(controle);
		traj_real_K.push_back(ganho);
		traj_real_x.push_back(estado_atual);
	}
	std::cout << "\nEstado final alcancado:\n" << traj_real_x.back().transpose() << std::endl;
	std::cout << std::endl << "Erro |xgoal - xtraj|:\n" << (xgoal.transpose() - traj_real_x.back().transpose()).cwiseAbs() << std::endl << std::endl;

	std::cout << "Tempo medio por passo do MPC (Warm-Start): " << tempo_total_mpc / mpc_steps << " ms" << std::endl;

	std::ofstream arquivo("trajetoria_iiwa.csv");

	for (int i = 1; i <= Nx; ++i) arquivo << "x" << i << ",";
	for (int i = 1; i <= Nu; ++i) arquivo << "u" << i << ",";
	for (int i = 1; i <= Nu * Nx; ++i) arquivo << (i < Nu * Nx ? "K" + std::to_string(i) + "," : "K" + std::to_string(i) + "\n");

	for (size_t i = 0; i < traj_real_u.size(); ++i)
	{    
		for (int j = 0; j < Nx; ++j) arquivo << traj_real_x[i][j] << ",";
		
		for (int j = 0; j < Nu; ++j) arquivo << traj_real_u[i][j] << ",";
		
		for (int row = 0; row < Nu; ++row) {
			for (int col = 0; col < Nx; ++col) {
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
	std::cout << "Trajetoria MPC salva com sucesso em 'trajetoria_iiwa.csv'!" << std::endl;;

	return 0;
}
