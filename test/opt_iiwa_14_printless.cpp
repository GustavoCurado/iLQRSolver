#include "iLQRSolverOptimized.hpp"

int main()
{
	double h = 0.01;
	double Tfinal = 5.0;
	int Nt = (int)(Tfinal/h) + 1;
//	int Nt = 64;

	iLQRSolverOptimized ilqr(h, Nt);

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
	Q_diag << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
		  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001;
	Eigen::MatrixXd Q = Q_diag.asDiagonal();
	Eigen::MatrixXd R = 0.001 * Eigen::MatrixXd::Identity(Nu, Nu);
	Eigen::MatrixXd Qn = 1000.0 * Eigen::MatrixXd::Identity(Nx, Nx);

	ilqr.setCosts(Q, R, Qn);
	ilqr.setBoundaryConditions(x0, xgoal);
	ilqr.setGravity(0.0);                    //Gravidade tá zerada

	ilqr.solve();

	return 0;
}


