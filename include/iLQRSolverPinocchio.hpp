#pragma once
#include "iLQRSolverBase.hpp"
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

//LEMBRA DE MODIFICAR ESSA MONTANHA DE INCLUDE; PODE TER COISA DESNECESSÁRIA AQUI

class iLQRSolverPinocchio : public iLQRSolverBase
{
	private:
		pinocchio::Model model;
		pinocchio::Data data;
		
		//Actuator selection matrix
		Eigen::MatrixXd S;

		void initializeMemory();

		Eigen::VectorXd continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u) override; // x' = f(x, u)
		void getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju) override;

	public:
		iLQRSolverPinocchio(const pinocchio::Model& model, double h, int Nt);
		iLQRSolverPinocchio(const std::string& urdf_path, double h, int Nt);
		iLQRSolverPinocchio(const pinocchio::Model& model, const Eigen::MatrixXd& S, double h, int Nt);
		iLQRSolverPinocchio(const std::string& urdf_path, const Eigen::MatrixXd& S, double h, int Nt);

		void setSelectionMatrix(const Eigen::MatrixXd& S);
		void setViscousDamping(const Eigen::VectorXd& c);
		void setGravity(double g_z);
};

