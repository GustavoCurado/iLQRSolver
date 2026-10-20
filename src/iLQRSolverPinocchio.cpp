#include "iLQRSolverPinocchio.hpp"

#include <omp.h>

iLQRSolverPinocchio::iLQRSolverPinocchio(const pinocchio::Model& model, double h, int Nt) : iLQRSolverBase(h, Nt)
{
	this->model = model;
	this->S = Eigen::MatrixXd::Identity(this->model.nv, this->model.nv);
	this->initializeMemory();
}

iLQRSolverPinocchio::iLQRSolverPinocchio(const std::string& urdf_path, double h, int Nt) : iLQRSolverBase(h, Nt)
{
	pinocchio::urdf::buildModel(urdf_path, this->model);
	this->S = Eigen::MatrixXd::Identity(this->model.nv, this->model.nv);
	this->initializeMemory();
}

iLQRSolverPinocchio::iLQRSolverPinocchio(const pinocchio::Model& model, const Eigen::MatrixXd& S, double h, int Nt) : iLQRSolverBase(h, Nt)
{
	this->model = model;
	this->S = S;
	this->initializeMemory();
}

iLQRSolverPinocchio::iLQRSolverPinocchio(const std::string& urdf_path, const Eigen::MatrixXd& S, double h, int Nt) : iLQRSolverBase(h, Nt)
{
	pinocchio::urdf::buildModel(urdf_path, this->model);
	this->S = S;
	this->initializeMemory();
}

void iLQRSolverPinocchio::initializeMemory()
{
	Nx = model.nq + model.nv;

	this->datas.assign(N_THREADS, pinocchio::Data(this->model));

	if (S.size() > 0)
	{
		Nu = S.cols();
	}
	else
	{
		Nu = model.nv;
		S = Eigen::MatrixXd::Identity(model.nv, model.nv);
	}

	iLQRSolverBase::initializeMemory();
}

Eigen::VectorXd iLQRSolverPinocchio::continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	pinocchio::Data& local_data = this->datas[omp_get_thread_num()];

	Eigen::VectorXd q = x.head(model.nq);
	Eigen::VectorXd v = x.tail(model.nv);

	Eigen::VectorXd tau = S * u;

	//Subtracts viscous damping from the torques
	tau -= model.damping.cwiseProduct(v);

	//Calculates joint acceleration and saves it into data
	pinocchio::aba(model, local_data, q, v, tau);

	Eigen::VectorXd xdot(Nx);
	xdot.head(model.nv) = v;
	xdot.tail(model.nv) = local_data.ddq;

	return xdot;
}

void iLQRSolverPinocchio::getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& xdot, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju)
{
	pinocchio::Data& local_data = this->datas[omp_get_thread_num()];

	Eigen::VectorXd q = x.head(model.nq);
	Eigen::VectorXd v = x.tail(model.nv);

	Eigen::VectorXd tau = S * u - model.damping.cwiseProduct(v);

	//Calculates dqdd/dq and dqdd/dv and saves it in data
	pinocchio::computeABADerivatives(model, local_data, q, v, tau);

	xdot.resize(Nx);
	xdot.head(model.nv) = v;
	xdot.tail(model.nv) = local_data.ddq;

	Jx.setZero(Nx, Nx);
	Jx.topRightCorner(model.nq, model.nv).setIdentity(); 
	Jx.bottomLeftCorner(model.nv, model.nq) = local_data.ddq_dq;
	Jx.bottomRightCorner(model.nv, model.nv) = local_data.ddq_dv - local_data.Minv * model.damping.asDiagonal();

	Ju.setZero(Nx, Nu);
	Ju.bottomRows(model.nv) = local_data.Minv * S;  // dqdd/du = M^-1
}

void iLQRSolverPinocchio::setSelectionMatrix(const Eigen::MatrixXd& S)
{
	this->S = S;

	if (this->Nu != S.cols())
	{
		this->Nu = S.cols();

		this->R = Eigen::MatrixXd::Zero(this->Nu, this->Nu);

		this->utraj.assign(Nt - 1, Eigen::VectorXd::Zero(this->Nu));
		this->utraj_new.assign(Nt - 1, Eigen::VectorXd::Zero(this->Nu));

		this->d.assign(Nt - 1, Eigen::VectorXd::Zero(this->Nu));
		this->K.assign(Nt - 1, Eigen::MatrixXd::Zero(this->Nu, this->Nx));
	}
}

void iLQRSolverPinocchio::setViscousDamping(const Eigen::VectorXd& c)
{
	this->model.damping = c;
}

void iLQRSolverPinocchio::setGravity(double g_z)
{
	this->model.gravity.linear() = Eigen::Vector3d(0.0, 0.0, g_z);
}
