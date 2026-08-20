#include "iLQRSolver.hpp"


iLQRSolver::iLQRSolver(const pinocchio::Model& model, double h, int Nt)
{
	this->model = model;
	this->data = pinocchio::Data(this->model);
	this->h = h;
	this->Nt = Nt;
	this->initializeMemory();
}

iLQRSolver::iLQRSolver(const std::string& urdf_path, double h, int Nt)
{
	pinocchio::urdf::buildModel(urdf_path, this->model);
	this->data = pinocchio::Data(this->model);
	this->h = h;
	this->Nt = Nt;
	this->initializeMemory();
}

iLQRSolver::iLQRSolver(const pinocchio::Model& model, const Eigen::MatrixXd& S, double h, int Nt)
{
	this->model = model;
	this->data = pinocchio::Data(this->model);
	this->S = S;
	this->h = h;
	this->Nt = Nt;
	this->initializeMemory();
}

iLQRSolver::iLQRSolver(const std::string& urdf_path, const Eigen::MatrixXd& S, double h, int Nt)
{
	pinocchio::urdf::buildModel(urdf_path, this->model);
	this->data = pinocchio::Data(this->model);
	this->S = S;
	this->h = h;
	this->Nt = Nt;
	this->initializeMemory();
}

void iLQRSolver::initializeMemory()
{
	Nx = model.nq + model.nv;

	if (S.size() > 0)
	{
		Nu = S.cols();
	}
	else
	{
		Nu = model.nv;
		S = Eigen::MatrixXd::Identity(model.nv, model.nv);
	}

	Q = Eigen::MatrixXd::Zero(Nx, Nx);
	R = Eigen::MatrixXd::Zero(Nu, Nu);
	Qn = Eigen::MatrixXd::Zero(Nx, Nx);

	x0 = Eigen::VectorXd::Zero(Nx);
	xgoal = Eigen::VectorXd::Zero(Nx);

	xtraj.assign(Nt, Eigen::VectorXd::Zero(Nx));
	utraj.assign(Nt - 1, Eigen::VectorXd::Zero(Nu));

	xtraj_new.assign(Nt, Eigen::VectorXd::Zero(Nx));
	utraj_new.assign(Nt - 1, Eigen::VectorXd::Zero(Nu));

	d.assign(Nt - 1, Eigen::VectorXd::Zero(Nu));
	K.assign(Nt - 1, Eigen::MatrixXd::Zero(Nu, Nx));

	p.assign(Nt, Eigen::VectorXd::Zero(Nx));
	P.assign(Nt, Eigen::MatrixXd::Zero(Nx, Nx));
}

Eigen::VectorXd iLQRSolver::continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	Eigen::VectorXd q = x.head(model.nq);
	Eigen::VectorXd v = x.tail(model.nv);

	Eigen::VectorXd tau = S * u;

	//Subtracts viscous damping from the torques
	tau -= model.damping.cwiseProduct(v);

	//Calculates joint acceleration and saves it into data
	pinocchio::aba(model, data, q, v, tau);

	Eigen::VectorXd xdot(Nx);
	xdot.head(model.nv) = v;
	xdot.tail(model.nv) = data.ddq;

	return xdot;
}

Eigen::VectorXd iLQRSolver::discreteDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	//The 4th-order Runge-Kutta method is used for dynamics discretization here

	Eigen::VectorXd k1 = continuousDynamics(x, u);
	Eigen::VectorXd k2 = continuousDynamics(x + 0.5 * h * k1, u);
	Eigen::VectorXd k3 = continuousDynamics(x + 0.5 * h * k2, u);
	Eigen::VectorXd k4 = continuousDynamics(x + h * k3, u);

	return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void iLQRSolver::computeDerivatives(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& A, Eigen::MatrixXd& B)
{
	// Calculates Jx = [0  I ; dqdd/dq dqdd/dv] and Ju = [0 ; dqdd/du]
	auto getContinuousJacs = [&](const Eigen::VectorXd& state, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju) 
	{
		Eigen::VectorXd q = state.head(model.nq);
		Eigen::VectorXd v = state.tail(model.nv);

		Eigen::VectorXd tau = S * u - model.damping.cwiseProduct(v);

		//Calculates dqdd/dq and dqdd/dv and saves it in data
		pinocchio::computeABADerivatives(model, data, q, v, tau);

		Jx.setZero(Nx, Nx);
		Jx.topRightCorner(model.nq, model.nv).setIdentity(); 
		Jx.bottomLeftCorner(model.nv, model.nq) = data.ddq_dq;
		Jx.bottomRightCorner(model.nv, model.nv) = data.ddq_dv - data.Minv * model.damping.asDiagonal();

		Ju.setZero(Nx, Nu);
		Ju.bottomRows(model.nv) = data.Minv * S;  // dqdd/du = M^-1
	};

	Eigen::VectorXd k1 = continuousDynamics(x, u);
	Eigen::VectorXd k2 = continuousDynamics(x + 0.5 * h * k1, u);
	Eigen::VectorXd k3 = continuousDynamics(x + 0.5 * h * k2, u);

	Eigen::MatrixXd Jx1, Jx2, Jx3, Jx4;
	Eigen::MatrixXd Ju1, Ju2, Ju3, Ju4;

	getContinuousJacs(x, Jx1, Ju1);
	getContinuousJacs(x + 0.5 * h * k1, Jx2, Ju2);
	getContinuousJacs(x + 0.5 * h * k2, Jx3, Ju3);
	getContinuousJacs(x + h * k3, Jx4, Ju4);

	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(Nx, Nx);

	Eigen::MatrixXd dx_k1 = Jx1;
	Eigen::MatrixXd dx_k2 = Jx2 * (I + 0.5 * h * dx_k1);
	Eigen::MatrixXd dx_k3 = Jx3 * (I + 0.5 * h * dx_k2);
	Eigen::MatrixXd dx_k4 = Jx4 * (I + h * dx_k3);

	A = I + (h / 6.0) * (dx_k1 + 2.0 * dx_k2 + 2.0 * dx_k3 + dx_k4);

	Eigen::MatrixXd du_k1 = Ju1;
	Eigen::MatrixXd du_k2 = Jx2 * (0.5 * h * du_k1) + Ju2;
	Eigen::MatrixXd du_k3 = Jx3 * (0.5 * h * du_k2) + Ju3;
	Eigen::MatrixXd du_k4 = Jx4 * (h * du_k3) + Ju4;

	B = (h / 6.0) * (du_k1 + 2.0 * du_k2 + 2.0 * du_k3 + du_k4);
}

double iLQRSolver::stageCost(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	Eigen::VectorXd error = x - xgoal;
	return 0.5 * (error.dot(Q * error) + u.dot(R * u));
}

double iLQRSolver::terminalCost(const Eigen::VectorXd& x)
{
	Eigen::VectorXd error = x - xgoal;
	return 0.5 * error.dot(Qn * error);
}

double iLQRSolver::calculateTotalCost(const std::vector<Eigen::VectorXd>& x_traj, const std::vector<Eigen::VectorXd>& u_traj)
{
	double total_cost = 0.0;

	for (int k = 0; k < Nt - 1; ++k)
	{
		total_cost += stageCost(x_traj[k], u_traj[k]);
	}

	total_cost += terminalCost(x_traj.back());

	return total_cost;
}

void iLQRSolver::setSelectionMatrix(const Eigen::MatrixXd& S)
{
	this->S = S;

	//The problem with this function is that, if S changes after the class constructor has been called, we are forced to change the dimensions of many variables, including R
	//Because of that, calling setCosts() after this one would be mandatory
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

void iLQRSolver::setCosts(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R, const Eigen::MatrixXd& Qn)
{
	this->Q = Q;
	this->R = R;
	this->Qn = Qn;
}

void iLQRSolver::setBoundaryConditions(const Eigen::VectorXd& x0, const Eigen::VectorXd& xgoal)
{
	this->x0 = x0;
	this->xgoal = xgoal;

	if (!this->xtraj.empty())
	{
		this->xtraj[0] = this->x0;
		this->xtraj_new[0] = this->x0;
	}
}

void iLQRSolver::setViscousDamping(const Eigen::VectorXd& c)
{
	this->model.damping = c;
}

void iLQRSolver::setGravity(double g_z)
{
	this->model.gravity.linear() = Eigen::Vector3d(0.0, 0.0, g_z);
}

bool iLQRSolver::backwardPass(double& delta_J, double& max_d)
{
	p.back() = Qn * (xtraj.back() - xgoal);
	P.back() = Qn;

	delta_J = 0.0;
	max_d = 0.0;

	Eigen::MatrixXd A(Nx, Nx);
	Eigen::MatrixXd B(Nx, Nu);

	for (int k = Nt - 2; k >= 0; k--)
	{
		Eigen::VectorXd q = Q * (xtraj[k] - xgoal);
		Eigen::VectorXd r = R * utraj[k];    

		computeDerivatives(xtraj[k], utraj[k], A, B);

		Eigen::VectorXd gx = q + A.transpose() * p[k+1];
		Eigen::VectorXd gu = r + B.transpose() * p[k+1];

		Eigen::MatrixXd Gxx = Q + A.transpose() * P[k+1] * A;
		Eigen::MatrixXd Guu = R + B.transpose() * P[k+1] * B;
		Eigen::MatrixXd Gxu = A.transpose() * P[k+1] * B;      
		Eigen::MatrixXd Gux = B.transpose() * P[k+1] * A;

		double beta = 1e-5; 
		ldlt.compute(Guu);

		//Guu regularization
		while (!(ldlt.info() == Eigen::Success && ldlt.isPositive()))
		{
		    Guu.diagonal().array() += beta;
		    beta *= 2.0;
		    
		    if (beta > 1e4)
			return false; 

		    ldlt.compute(Guu);
		}

		d[k] = ldlt.solve(gu);
		K[k] = ldlt.solve(Gux);

		p[k] = gx - K[k].transpose() * gu + K[k].transpose() * Guu * d[k] - Gxu * d[k];
		P[k] = Gxx + K[k].transpose() * Guu * K[k] - Gxu * K[k] - K[k].transpose() * Gux;

		//Estimated cost reduction (used in line search)
		delta_J += gu.dot(d[k]); 

		double current_max_d = d[k].cwiseAbs().maxCoeff();
		if (current_max_d > max_d)
		{
		    max_d = current_max_d;
		}
	}

	return true;
}

double iLQRSolver::forwardRollout(double alpha)
{
	xtraj_new[0] = x0;

	for (int k = 0; k < Nt - 1; ++k)
	{
		Eigen::VectorXd delta_x = xtraj_new[k] - xtraj[k];

		utraj_new[k] = utraj[k] - alpha * d[k] - K[k] * delta_x;
		xtraj_new[k+1] = discreteDynamics(xtraj_new[k], utraj_new[k]);
	}

	return calculateTotalCost(xtraj_new, utraj_new);
}

bool iLQRSolver::solve(int max_iters, double tol)
{
	//Initial rollout
	xtraj[0] = x0;
	for (int k = 0; k < Nt - 1; ++k)
	{
		xtraj[k + 1] = discreteDynamics(xtraj[k], utraj[k]);
	}

	double J = calculateTotalCost(xtraj, utraj);
	double max_d = tol + 1.0; //This initialization value is here just to force entry inside loop

	int iter = 0;

	while (max_d > tol && iter < max_iters)
	{
		iter++;

		max_d = 0.0;
		double delta_J = 0.0;

		bool success = backwardPass(delta_J, max_d);
		if (!success)
		    return false;

		double alpha = 1.0;
		double J_new;

		//Line search on the forward rollout
		do
		{
			J_new = forwardRollout(alpha);

			if (J_new <= J - 1e-4 * alpha * delta_J)
				break;

			alpha *= 0.5;
		    
		} while (alpha > 1e-8); //This tolerance is here to avoid underflow, which would cause an infinite loop

		xtraj = xtraj_new;
		utraj = utraj_new;
		J = J_new;
	}

	return true;
}
