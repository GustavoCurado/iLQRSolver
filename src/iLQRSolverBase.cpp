#include "iLQRSolverBase.hpp"

iLQRSolverBase::iLQRSolverBase(double h, int Nt)
{
	this->h = h;
	this->Nt = Nt;
    // Note: Nx and Nu must be initialized by the derived class before initializeMemory is called.
}

void iLQRSolverBase::initializeMemory()
{
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

Eigen::VectorXd iLQRSolverBase::discreteDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	//The 4th-order Runge-Kutta method is used for dynamics discretization here

	Eigen::VectorXd k1 = continuousDynamics(x, u);
	Eigen::VectorXd k2 = continuousDynamics(x + 0.5 * h * k1, u);
	Eigen::VectorXd k3 = continuousDynamics(x + 0.5 * h * k2, u);
	Eigen::VectorXd k4 = continuousDynamics(x + h * k3, u);

	return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void iLQRSolverBase::computeDerivatives(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& A, Eigen::MatrixXd& B)
{
	Eigen::VectorXd k1(Nx), k2(Nx), k3(Nx), k4(Nx);
	Eigen::MatrixXd Jx1(Nx, Nx), Jx2(Nx, Nx), Jx3(Nx, Nx), Jx4(Nx, Nx);
	Eigen::MatrixXd Ju1(Nx, Nu), Ju2(Nx, Nu), Ju3(Nx, Nu), Ju4(Nx, Nu);

	getContinuousJacs(x, u, k1, Jx1, Ju1);
	getContinuousJacs(x + 0.5 * h * k1, u, k2, Jx2, Ju2);
	getContinuousJacs(x + 0.5 * h * k2, u, k3, Jx3, Ju3);
	getContinuousJacs(x + h * k3, u, k4, Jx4, Ju4);

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

double iLQRSolverBase::stageCost(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
	Eigen::VectorXd error = x - xgoal;
	return 0.5 * (error.dot(Q * error) + u.dot(R * u));
}

double iLQRSolverBase::terminalCost(const Eigen::VectorXd& x)
{
	Eigen::VectorXd error = x - xgoal;
	return 0.5 * error.dot(Qn * error);
}

double iLQRSolverBase::calculateTotalCost(const std::vector<Eigen::VectorXd>& x_traj, const std::vector<Eigen::VectorXd>& u_traj)
{
	double total_cost = 0.0;

	for (int k = 0; k < Nt - 1; ++k)
	{
		total_cost += stageCost(x_traj[k], u_traj[k]);
	}

	total_cost += terminalCost(x_traj.back());

	return total_cost;
}

void iLQRSolverBase::setCosts(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R, const Eigen::MatrixXd& Qn)
{
	this->Q = Q;
	this->R = R;
	this->Qn = Qn;
}

void iLQRSolverBase::setBoundaryConditions(const Eigen::VectorXd& x0, const Eigen::VectorXd& xgoal)
{
	this->x0 = x0;
	this->xgoal = xgoal;

	if (!this->xtraj.empty())
	{
		this->xtraj[0] = this->x0;
		this->xtraj_new[0] = this->x0;
	}
}

void iLQRSolverBase::setGoal(const Eigen::VectorXd& xgoal)
{
	this->xgoal = xgoal;
}

bool iLQRSolverBase::backwardPass(double& delta_J, double& max_d)
{
	p.back() = Qn * (xtraj.back() - xgoal);
	P.back() = Qn;

	delta_J = 0.0;
	max_d = 0.0;

	Eigen::MatrixXd A(Nx, Nx);
	Eigen::MatrixXd B(Nx, Nu);
	std::vector<Eigen::MatrixXd> Atraj(Nt-1), Btraj(Nt-1);

	struct timespec t_start, t_end;
	clock_gettime(CLOCK_MONOTONIC, &t_start);

	#pragma omp parallel for num_threads(N_THREADS) if(use_parallel_derivatives)
	for (int k = 0; k < Nt - 1; k++)
	{
		Eigen::MatrixXd local_A(Nx, Nx);
		Eigen::MatrixXd local_B(Nx, Nu);
		computeDerivatives(xtraj[k], utraj[k], local_A, local_B);
		Atraj[k] = local_A;
		Btraj[k] = local_B;
	}

	clock_gettime(CLOCK_MONOTONIC, &t_end);
	time_derivatives_us += time_delta_us(t_start, t_end);

	clock_gettime(CLOCK_MONOTONIC, &t_start);

	for (int k = Nt - 2; k >= 0; k--)
	{
		Eigen::VectorXd q = Q * (xtraj[k] - xgoal);
		Eigen::VectorXd r = R * utraj[k];    

		A = Atraj[k];
		B = Btraj[k];

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

	clock_gettime(CLOCK_MONOTONIC, &t_end);
	time_riccati_us += time_delta_us(t_start, t_end);

	return true;
}

double iLQRSolverBase::forwardRollout(double alpha)
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

bool iLQRSolverBase::solve(int max_iters, double tol, double time_limit_ms)
{
	time_derivatives_us = 0.0;
	time_riccati_us = 0.0;
	time_rollout_us = 0.0;
	time_total_us = 0.0;

	struct timespec solve_start, solve_end, roll_start, roll_end, iter_now;
	clock_gettime(CLOCK_MONOTONIC, &solve_start);

	//Initial rollout
	clock_gettime(CLOCK_MONOTONIC, &roll_start);
	xtraj[0] = x0;
	for (int k = 0; k < Nt - 1; ++k)
	{
		xtraj[k + 1] = discreteDynamics(xtraj[k], utraj[k]);
	}

	double J = calculateTotalCost(xtraj, utraj);

	clock_gettime(CLOCK_MONOTONIC, &roll_end);
	time_rollout_us += time_delta_us(roll_start, roll_end);

	double max_d = tol + 1.0; //This initialization value is here just to force entry inside loop

	int iter = 0;

	while (max_d > tol && iter < max_iters)
	{
		if (time_limit_ms > 0.0 && iter > 0)
		{
			clock_gettime(CLOCK_MONOTONIC, &iter_now);
			if (time_delta_us(solve_start, iter_now) >= time_limit_ms * 1000.0)
				break;
		}

		iter++;

		max_d = 0.0;
		double delta_J = 0.0;

		bool success = backwardPass(delta_J, max_d);
		if (!success)
			return false;

		double alpha = 1.0;
		double J_new;

		//Line search on the forward rollout
		clock_gettime(CLOCK_MONOTONIC, &roll_start);
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

		clock_gettime(CLOCK_MONOTONIC, &roll_end);
		time_rollout_us += time_delta_us(roll_start, roll_end);
	}

	clock_gettime(CLOCK_MONOTONIC, &solve_end);
	time_total_us = time_delta_us(solve_start, solve_end);
	last_iters = iter;

	return true;
}

void iLQRSolverBase::advanceMPCStep(const Eigen::VectorXd& estado_atual)
{
	this->x0 = estado_atual;
	if (!this->xtraj.empty()) {
		this->xtraj[0] = this->x0;
		this->xtraj_new[0] = this->x0;
	}

	for (int k = 0; k < Nt - 2; ++k)
	{
		utraj[k] = utraj[k+1];
	}
}

void iLQRSolverBase::resetTrajectories()
{
	xtraj.assign(Nt, Eigen::VectorXd::Zero(Nx));
	utraj.assign(Nt - 1, Eigen::VectorXd::Zero(Nu));
	xtraj[0] = x0;
}

Eigen::VectorXd iLQRSolverBase::simulatePlant(const Eigen::VectorXd& x_inicial, double dt_total, int substeps)
{
	double dt = dt_total / substeps;
	Eigen::VectorXd x = x_inicial;

	for (int s = 0; s < substeps; s++)
	{
		double t_decorrido = s * dt;

		int k = static_cast<int>(t_decorrido / h);
		if (k >= Nt - 1)
			k = Nt - 2;

		Eigen::VectorXd u = utraj[k] - K[k] * (x - xtraj[k]);

		Eigen::VectorXd k1 = continuousDynamics(x, u);
		Eigen::VectorXd k2 = continuousDynamics(x + 0.5 * dt * k1, u);
		Eigen::VectorXd k3 = continuousDynamics(x + 0.5 * dt * k2, u);
		Eigen::VectorXd k4 = continuousDynamics(x + dt * k3, u);

		x = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}

	return x;
}
