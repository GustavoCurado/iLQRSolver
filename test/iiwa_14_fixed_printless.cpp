#include "iLQRSolverFixed.hpp"

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

	ilqr.setCosts(Q, R, Qn);
	ilqr.setBoundaryConditions(x0, xgoal);
	ilqr.setGravity(0.0);                    //Gravidade tá zerada

	ilqr.solve();
	ilqr.solve();
	ilqr.solve();
	ilqr.solve();
	ilqr.solve();
	ilqr.solve();
	ilqr.solve();
	ilqr.solve();
	ilqr.solve();
	ilqr.solve();

	return 0;
}


