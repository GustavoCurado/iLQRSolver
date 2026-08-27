#include "iLQRSolverFixed.hpp"


iLQRSolver::iLQRSolver(const pinocchio::Model& model, double h, int Nt)
{
	this->model = model;
	this->data = pinocchio::Data(this->model);
	this->h = h;
	this->Nt = Nt;
	this->S = Eigen::Matrix<double, 7, 7>::Identity();
	this->initializeMemory();
}

iLQRSolver::iLQRSolver(const std::string& urdf_path, double h, int Nt)
{
	pinocchio::urdf::buildModel(urdf_path, this->model);
	this->data = pinocchio::Data(this->model);
	this->h = h;
	this->Nt = Nt;
	this->S = Eigen::Matrix<double, 7, 7>::Identity();
	this->initializeMemory();
}

iLQRSolver::iLQRSolver(const pinocchio::Model& model, const Eigen::Matrix<double, 7, 7>& S, double h, int Nt)
{
	this->model = model;
	this->data = pinocchio::Data(this->model);
	this->S = S;
	this->h = h;
	this->Nt = Nt;
	this->initializeMemory();
}

iLQRSolver::iLQRSolver(const std::string& urdf_path, const Eigen::Matrix<double, 7, 7>& S, double h, int Nt)
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
	Nu = 7;

	Q = Eigen::Matrix<double, 14, 14>::Zero();
	R = Eigen::Matrix<double, 7, 7>::Zero();
	Qn = Eigen::Matrix<double, 14, 14>::Zero();

	x0 = Eigen::Matrix<double, 14, 1>::Zero();
	xgoal = Eigen::Matrix<double, 14, 1>::Zero();

	xtraj.assign(Nt, Eigen::Matrix<double, 14, 1>::Zero());
	utraj.assign(Nt - 1, Eigen::Matrix<double, 7, 1>::Zero());

	xtraj_new.assign(Nt, Eigen::Matrix<double, 14, 1>::Zero());
	utraj_new.assign(Nt - 1, Eigen::Matrix<double, 7, 1>::Zero());

	d.assign(Nt - 1, Eigen::Matrix<double, 7, 1>::Zero());
	K.assign(Nt - 1, Eigen::Matrix<double, 7, 14>::Zero());

	p.assign(Nt, Eigen::Matrix<double, 14, 1>::Zero());
	P.assign(Nt, Eigen::Matrix<double, 14, 14>::Zero());
}

Eigen::Matrix<double, 14, 1> iLQRSolver::continuousDynamics(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u)
{
	Eigen::Matrix<double, 7, 1> q = x.head(model.nq);
	Eigen::Matrix<double, 7, 1> v = x.tail(model.nv);

	Eigen::Matrix<double, 7, 1> tau = S * u;

	//Subtracts viscous damping from the torques
	tau -= model.damping.cwiseProduct(v);

	//Calculates joint acceleration and saves it into data
	pinocchio::aba(model, data, q, v, tau);

	Eigen::Matrix<double, 14, 1> xdot;
	xdot.head(model.nv) = v;
	xdot.tail(model.nv) = data.ddq;

	return xdot;
}

Eigen::Matrix<double, 14, 1> iLQRSolver::discreteDynamics(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u)
{
	//The 4th-order Runge-Kutta method is used for dynamics discretization here

	Eigen::Matrix<double, 14, 1> k1 = continuousDynamics(x, u);
	Eigen::Matrix<double, 14, 1> k2 = continuousDynamics(x + 0.5 * h * k1, u);
	Eigen::Matrix<double, 14, 1> k3 = continuousDynamics(x + 0.5 * h * k2, u);
	Eigen::Matrix<double, 14, 1> k4 = continuousDynamics(x + h * k3, u);

	return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void iLQRSolver::computeDerivatives(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u, Eigen::Matrix<double, 14, 14>& A, Eigen::Matrix<double, 14, 7>& B)
{
	// Calculates Jx = [0  I ; dqdd/dq dqdd/dv] and Ju = [0 ; dqdd/du]
	auto getContinuousJacs = [&](const Eigen::Matrix<double, 14, 1>& state, Eigen::Matrix<double, 14, 14>& Jx, Eigen::Matrix<double, 14, 7>& Ju) 
	{
		Eigen::Matrix<double, 7, 1> q = state.head(model.nq);
		Eigen::Matrix<double, 7, 1> v = state.tail(model.nv);

		Eigen::Matrix<double, 7, 1> tau = S * u - model.damping.cwiseProduct(v);

		//Calculates dqdd/dq and dqdd/dv and saves it in data
		pinocchio::computeABADerivatives(model, data, q, v, tau);

		Jx.setZero();
		Jx.topRightCorner(model.nq, model.nv).setIdentity(); 
		Jx.bottomLeftCorner(model.nv, model.nq) = data.ddq_dq;
		Jx.bottomRightCorner(model.nv, model.nv) = data.ddq_dv - data.Minv * model.damping.asDiagonal();

		Ju.setZero();
		Ju.bottomRows(model.nv) = data.Minv * S;  // dqdd/du = M^-1
	};

	Eigen::Matrix<double, 14, 1> k1 = continuousDynamics(x, u);
	Eigen::Matrix<double, 14, 1> k2 = continuousDynamics(x + 0.5 * h * k1, u);
	Eigen::Matrix<double, 14, 1> k3 = continuousDynamics(x + 0.5 * h * k2, u);

	Eigen::Matrix<double, 14, 14> Jx1, Jx2, Jx3, Jx4;
	Eigen::Matrix<double, 14, 7> Ju1, Ju2, Ju3, Ju4;

	getContinuousJacs(x, Jx1, Ju1);
	getContinuousJacs(x + 0.5 * h * k1, Jx2, Ju2);
	getContinuousJacs(x + 0.5 * h * k2, Jx3, Ju3);
	getContinuousJacs(x + h * k3, Jx4, Ju4);

	Eigen::Matrix<double, 14, 14> I = Eigen::Matrix<double, 14, 14>::Identity();

	Eigen::Matrix<double, 14, 14> dx_k1 = Jx1;
	Eigen::Matrix<double, 14, 14> dx_k2 = Jx2 * (I + 0.5 * h * dx_k1);
	Eigen::Matrix<double, 14, 14> dx_k3 = Jx3 * (I + 0.5 * h * dx_k2);
	Eigen::Matrix<double, 14, 14> dx_k4 = Jx4 * (I + h * dx_k3);

	A = I + (h / 6.0) * (dx_k1 + 2.0 * dx_k2 + 2.0 * dx_k3 + dx_k4);

	Eigen::Matrix<double, 14, 7> du_k1 = Ju1;
	Eigen::Matrix<double, 14, 7> du_k2 = Jx2 * (0.5 * h * du_k1) + Ju2;
	Eigen::Matrix<double, 14, 7> du_k3 = Jx3 * (0.5 * h * du_k2) + Ju3;
	Eigen::Matrix<double, 14, 7> du_k4 = Jx4 * (h * du_k3) + Ju4;

	B = (h / 6.0) * (du_k1 + 2.0 * du_k2 + 2.0 * du_k3 + du_k4);
}

double iLQRSolver::stageCost(const Eigen::Matrix<double, 14, 1>& x, const Eigen::Matrix<double, 7, 1>& u)
{
	Eigen::Matrix<double, 14, 1> error = x - xgoal;
	return 0.5 * (error.dot(Q * error) + u.dot(R * u));
}

double iLQRSolver::terminalCost(const Eigen::Matrix<double, 14, 1>& x)
{
	Eigen::Matrix<double, 14, 1> error = x - xgoal;
	return 0.5 * error.dot(Qn * error);
}

double iLQRSolver::calculateTotalCost(const std::vector<Eigen::Matrix<double, 14, 1>>& x_traj, const std::vector<Eigen::Matrix<double, 7, 1>>& u_traj)
{
	double total_cost = 0.0;

	for (int k = 0; k < Nt - 1; ++k)
	{
		total_cost += stageCost(x_traj[k], u_traj[k]);
	}

	total_cost += terminalCost(x_traj.back());

	return total_cost;
}

void iLQRSolver::setSelectionMatrix(const Eigen::Matrix<double, 7, 7>& S)
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

void iLQRSolver::setCosts(const Eigen::Matrix<double, 14, 14>& Q, const Eigen::Matrix<double, 7, 7>& R, const Eigen::Matrix<double, 14, 14>& Qn)
{
	this->Q = Q;
	this->R = R;
	this->Qn = Qn;
}

void iLQRSolver::setBoundaryConditions(const Eigen::Matrix<double, 14, 1>& x0, const Eigen::Matrix<double, 14, 1>& xgoal)
{
	this->x0 = x0;
	this->xgoal = xgoal;

	if (!this->xtraj.empty())
	{
		this->xtraj[0] = this->x0;
		this->xtraj_new[0] = this->x0;
	}
}

void iLQRSolver::setViscousDamping(const Eigen::Matrix<double, 7, 1>& c)
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

	Eigen::Matrix<double, 14, 14> A;
	Eigen::Matrix<double, 14, 7> B;

	for (int k = Nt - 2; k >= 0; k--)
	{
		Eigen::Matrix<double, 14, 1> q = Q * (xtraj[k] - xgoal);
		Eigen::Matrix<double, 7, 1> r = R * utraj[k];    

		computeDerivatives(xtraj[k], utraj[k], A, B);

		Eigen::Matrix<double, 14, 1> gx = q + A.transpose() * p[k+1];
		Eigen::Matrix<double, 7, 1> gu = r + B.transpose() * p[k+1];

		Eigen::Matrix<double, 14, 14> Gxx = Q + A.transpose() * P[k+1] * A;
		Eigen::Matrix<double, 7, 7> Guu = R + B.transpose() * P[k+1] * B;
		Eigen::Matrix<double, 14, 7> Gxu = A.transpose() * P[k+1] * B;      
		Eigen::Matrix<double, 7, 14> Gux = B.transpose() * P[k+1] * A;

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
		Eigen::Matrix<double, 14, 1> delta_x = xtraj_new[k] - xtraj[k];

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
