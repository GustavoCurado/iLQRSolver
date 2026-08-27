#include "iLQRSolver.hpp"

int main()
{
	std::string urdf_filename = "utils/iiwa_14.urdf";

	double h = 0.01;
	double Tfinal = 5.0;
	int Nt = (int)(Tfinal/h) + 1;

	iLQRSolver ilqr(urdf_filename, h, Nt);

	int Nx = ilqr.getNx();
	int Nu = ilqr.getNu();

/*	Eigen::VectorXd x0(Nx);
	x0 << 0.01, 0.10, 0.01, -0.10, 0.01, 0.10, 0.01,            //Esse funciona
	      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        Eigen::VectorXd xgoal(Nx);
	xgoal << 02.00, 1.20, -0.30, -1.50, 0.40, 1.00, 0.50,
	         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; */

	Eigen::Matrix<double, 14, 1> x0;
	x0 << 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
	      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	
	Eigen::Matrix<double, 14, 1> xgoal;
	xgoal << 02.00, 1.20, -0.30, -1.50, 0.40, 1.00, 0.50,
	         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	Eigen::Matrix<double, 14, 1> Q_diag;                                   
	Q_diag << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
	
	Eigen::Matrix<double, 14, 14> Q = Q_diag.asDiagonal();
	Eigen::Matrix<double, 7, 7> R = 0.1 * Eigen::Matrix<double, 7, 7>::Identity();
	Eigen::Matrix<double, 14, 14> Qn = 10000.0 * Eigen::Matrix<double, 14, 14>::Identity();

/*	Eigen::VectorXd Q_diag(Nx);                                        //Converge, mas estoura limite
	Q_diag << 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
		      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
	Eigen::MatrixXd Q = Q_diag.asDiagonal();
	Eigen::MatrixXd R = 0.001 * Eigen::MatrixXd::Identity(Nu, Nu);
	Eigen::VectorXd Qn_diag(Nx);
	Qn_diag << 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
	           100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0;
	Eigen::MatrixXd Qn = Qn_diag.asDiagonal(); */

	ilqr.setCosts(Q, R, Qn);
	ilqr.setBoundaryConditions(x0, xgoal);
	ilqr.setGravity(0.0);                    //Gravidade tá zerada

	auto start_time = std::chrono::high_resolution_clock::now();
	ilqr.solve();
	auto end_time = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> duration = end_time - start_time;

	auto xtraj = ilqr.getXTraj();
	auto utraj = ilqr.getUTraj();

	std::cout << "Estado final alcancado:\n" << xtraj.back().transpose() << std::endl;
	std::cout << std::endl << "Erro |xgoal - xtraj|:\n" << (xgoal.transpose() - xtraj.back().transpose()).cwiseAbs() << std::endl << std::endl;

	std::cout << "Tempo de execucao do iLQR: " << duration.count() << " ms" << std::endl;

/*	std::ofstream arquivo("trajetoria_iiwa.csv");

	for (int i = 1; i <= Nu; ++i) arquivo << "theta" << i << ",";     //.csv sem as matrizes K
	for (int i = 1; i <= Nu; ++i) arquivo << "omega" << i << ",";
	for (int i = 1; i < Nu; ++i) arquivo << "torque" << i << ",";
	arquivo << "torque" << Nu << "\n";

	for (int i = 0; i < Nt - 1; ++i)
	{    
		for (int j = 0; j < Nx; ++j) 
		    arquivo << xtraj[i][j] << ",";

		for (int j = 0; j < Nu - 1; ++j)
		    arquivo << utraj[i][j] << ",";

		arquivo << utraj[i][Nu - 1] << "\n";
	}

	for (int j = 0; j < Nx; ++j)
	{
		arquivo << xtraj.back()[j] << ",";
	}

	for (int j = 0; j < Nu - 1; ++j)
	{
		arquivo << "0.0,";
	}

	arquivo << "0.0\n";

	arquivo.close();
	std::cout << "Trajetoria salva com sucesso em 'trajetoria_iiwa.csv'!" << std::endl; */ 

	auto Ktraj = ilqr.getKTraj();

	std::ofstream arquivo("trajetoria_iiwa.csv");

	for (int i = 1; i <= Nx; ++i) arquivo << "x" << i << ",";
	for (int i = 1; i <= Nu; ++i) arquivo << "u" << i << ",";
	for (int i = 1; i <= Nu * Nx; ++i) arquivo << (i < Nu * Nx ? "K" + std::to_string(i) + "," : "K" + std::to_string(i) + "\n");

	for (int i = 0; i < Nt - 1; ++i)
	{    
		for (int j = 0; j < Nx; ++j) arquivo << xtraj[i][j] << ",";

		for (int j = 0; j < Nu; ++j) arquivo << utraj[i][j] << ",";

		for (int row = 0; row < Nu; ++row) {
			for (int col = 0; col < Nx; ++col) {
				arquivo << Ktraj[i](row, col);
				if (row == Nu - 1 && col == Nx - 1) arquivo << "\n";
				else arquivo << ",";
			}
		}
	}

	for (int j = 0; j < Nx; ++j) arquivo << xtraj.back()[j] << ",";
	for (int j = 0; j < Nu; ++j) arquivo << "0.0,";
	for (int i = 0; i < Nu * Nx - 1; ++i) arquivo << "0.0,";
	arquivo << "0.0\n";

	arquivo.close();
	std::cout << "Trajetoria salva com sucesso em 'trajetoria_iiwa.csv'!" << std::endl;

	return 0;
}


