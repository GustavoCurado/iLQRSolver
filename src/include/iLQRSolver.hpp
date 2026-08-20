#pragma once

//#include "fast-rbd-gradients/utils/experiment_helpers.h" // include constants and other experiment consistency helpers

#define BOOST_BIND_GLOBAL_PLACEHOLDERS // fixes a warning I get

#include "pinocchio/spatial/fwd.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/center-of-mass.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/container/aligned-vector.hpp"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <cmath>
#include <vector>
#include <chrono>

#define EIGEN_DEFAULT_TO_ROW_MAJOR
#include <Eigen/Dense>
#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::VectorXd)

//LEMBRA DE MODIFICAR ESSA MONTANHA DE INCLUDE; PODE TER COISA DESNECESSÁRIA AQUI

class iLQRSolver
{
	private:
		pinocchio::Model model;
		pinocchio::Data data;

		int Nx;
		int Nu;
		double h; //Time step
		int Nt;   //Horizon size
		double Tfinal;
		
		//Actuator selection matrix
		Eigen::MatrixXd S;

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

		Eigen::VectorXd continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u); // x' = f(x, u)
		Eigen::VectorXd discreteDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u);   // xk+1 = f(xk, uk)

		//Calculates A = df/dx and B = df/du; f being discreteDynamics
		void computeDerivatives(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& A, Eigen::MatrixXd& B);

		double stageCost(const Eigen::VectorXd& x, const Eigen::VectorXd& u);
		double terminalCost(const Eigen::VectorXd& x);
		double calculateTotalCost(const std::vector<Eigen::VectorXd>& x_traj, const std::vector<Eigen::VectorXd>& u_traj);

		bool backwardPass(double& delta_J, double& max_d);
		double forwardRollout(double alpha);

	public:
		iLQRSolver(const pinocchio::Model& model, double h, int Nt);
		iLQRSolver(const std::string& urdf_path, double h, int Nt);
		iLQRSolver(const pinocchio::Model& model, const Eigen::MatrixXd& S, double h, int Nt);
		iLQRSolver(const std::string& urdf_path, const Eigen::MatrixXd& S, double h, int Nt);

		void setSelectionMatrix(const Eigen::MatrixXd& S);
		void setCosts(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R, const Eigen::MatrixXd& Qn);
		void setBoundaryConditions(const Eigen::VectorXd& x0, const Eigen::VectorXd& xgoal);
		void setViscousDamping(const Eigen::VectorXd& c);
		void setGravity(double g_z);

		bool solve(int max_iters = 1000, double tol = 1e-3);

		int getNx() const { return Nx; }
		int getNu() const { return Nu; }
		std::vector<Eigen::VectorXd> getXTraj() const { return xtraj; }
		std::vector<Eigen::VectorXd> getUTraj() const { return utraj; }
		std::vector<Eigen::MatrixXd> getKTraj() const { return K; }
};

