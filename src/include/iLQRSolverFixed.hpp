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
		Eigen::Matrix<double, 7, 7> S;

		//Cost matrixes
		Eigen::Matrix<double, 14, 14> Q;
		Eigen::Matrix<double, 7, 7> R;
		Eigen::Matrix<double, 14, 14> Qn;

		Eigen::Matrix<double, 14, 1> x0;
		Eigen::Matrix<double, 14, 1> xgoal;
		std::vector<Eigen::Matrix<double, 14, 1>> xtraj; 
		std::vector<Eigen::Matrix<double, 7, 1>> utraj;
		std::vector<Eigen::Matrix<double, 14, 1>> xtraj_new; 
		std::vector<Eigen::Matrix<double, 7, 1>> utraj_new;

		std::vector<Eigen::Matrix<double, 7, 1>> d; 
		std::vector<Eigen::Matrix<double, 7, 14>> K; 
		std::vector<Eigen::Matrix<double, 14, 1>> p; 
		std::vector<Eigen::Matrix<double, 14, 14>> P; 

		Eigen::LDLT<Eigen::Matrix<double, 7, 7>> ldlt;

		void initializeMemory();

		Eigen::Matrix<double, 14, 1> continuousDynamics(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u); // x' = f(x, u)
		Eigen::Matrix<double, 14, 1> discreteDynamics(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u);   // xk+1 = f(xk, uk)

		//Calculates A = df/dx and B = df/du; f being discreteDynamics
		void computeDerivatives(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u, Eigen::Matrix<double, 14, 14>& A, Eigen::Matrix<double, 14, 7>& B);

		double stageCost(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u);
		double terminalCost(const Eigen::Matrix<double, 14, 1>& x);
		double calculateTotalCost(const std::vector<Eigen::Matrix<double, 14, 1>>& x_traj, const std::vector<Eigen::Matrix<double, 7, 1>>& u_traj);

		bool backwardPass(double& delta_J, double& max_d);
		double forwardRollout(double alpha);

	public:
		iLQRSolver(const pinocchio::Model& model, double h, int Nt);
		iLQRSolver(const std::string& urdf_path, double h, int Nt);
		iLQRSolver(const pinocchio::Model& model, const Eigen::Matrix<double, 7, 7>& S, double h, int Nt);
		iLQRSolver(const std::string& urdf_path, const Eigen::Matrix<double, 7, 7>& S, double h, int Nt);

		void setSelectionMatrix(const Eigen::Matrix<double, 7, 7>& S);
		void setCosts(const Eigen::Matrix<double, 14, 14>& Q, const Eigen::Matrix<double, 7, 7>& R, const Eigen::Matrix<double, 14, 14>& Qn);
		void setBoundaryConditions(const Eigen::Matrix<double, 14, 1>& x0, const Eigen::Matrix<double, 14, 1>& xgoal);
		void setViscousDamping(const Eigen::Matrix<double, 7, 1>& c);
		void setGravity(double g_z);

		bool solve(int max_iters = 1000, double tol = 1e-3);

		int getNx() const { return Nx; }
		int getNu() const { return Nu; }
		std::vector<Eigen::Matrix<double, 14, 1>> getXTraj() const { return xtraj; }
		std::vector<Eigen::Matrix<double, 7, 1>> getUTraj() const { return utraj; }
		std::vector<Eigen::Matrix<double, 7, 14>> getKTraj() const { return K; }
};

