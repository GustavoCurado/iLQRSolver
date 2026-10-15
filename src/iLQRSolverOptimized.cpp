#include "iLQRSolverOptimized.hpp"

// Variável para controlar a gravidade
static double g_plancher_gravity = 9.81;

// Basicamente isso está aqui para enganar o código do Plancher para usar a variável desse código
#define GRAVITY g_plancher_gravity

#include "../utils/dynamicsGradient.h"

iLQRSolverOptimized::iLQRSolverOptimized(double h, int Nt) : iLQRSolverBase(h, Nt)
{
    this->use_parallel_derivatives = true;
    this->initializeMemory();
}

void iLQRSolverOptimized::initializeMemory()
{
    // Kuka LBR iiwa 14 tem 7 juntas
    this->Nx = 14;
    this->Nu = 7;
    
    iLQRSolverBase::initializeMemory(); // Chama a inicialização base (matrizes Q, R, xtraj, etc.)

    // Inicializa as transformações e inércias hardcoded do Plancher
    initTransforms<double>(this->d_T);
    initInertiaTensors<double>(this->d_I);
}

Eigen::VectorXd iLQRSolverOptimized::continuousDynamics(const Eigen::VectorXd& x, const Eigen::VectorXd& u)
{
    double q[7];
    double v[7];
    double tau[7];
    
    for(int i = 0; i < 7; i++) {
        q[i] = x[i];
        v[i] = x[i + 7];
        tau[i] = u[i]; // Assumindo S = Identidade
    }
    
    double qdd[7] = {0.0};
    double Minv[49] = {0.0};
    
    // Chamada da função do Plancher (Templates: T, MPC_MODE=false, VEL_DAMPING=true, DENSE_MINV=false)
    forwardDynamics<double, false, true, false>(qdd, tau, v, q, Minv);
    
    Eigen::VectorXd xdot = Eigen::VectorXd::Zero(14);
    xdot.head<7>() = x.tail<7>();
    
    xdot.tail<7>() = Eigen::Map<Eigen::Matrix<double, 7, 1>>(qdd);
    
    return xdot;
}

void iLQRSolverOptimized::getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju)
{
    // Inicialização de currKnot
    knot<double> currKnot = {0};
    for(int i = 0; i < 7; i++) {
        currKnot.q[i] = x[i];
        currKnot.qd[i] = x[i + 7];
        currKnot.u[i] = u[i];
    }
    
    double dqdd_array[98];
    
    // Chamada da função de gradiente do Plancher (Templates: T, QDD_MINV_PASSED_IN=false, MPC_MODE=false, VEL_DAMPING=true, EXTERNAL_TI=false)
    forwardDynamicsGradient<double, false, false, true, false>(dqdd_array, &currKnot, this->d_T, this->d_I);
    
    Jx.setZero(14, 14);
    Jx.topRightCorner<7, 7>().setIdentity();
    
    Jx.bottomLeftCorner<7, 7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(&dqdd_array[0]);
    Jx.bottomRightCorner<7, 7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(&dqdd_array[49]);
    
    // O código do Plancher gera o M inverso como uma matriz triangular superior
    // Dá bug se o código só usasse direto o Minv dele
    triangularToDense<double, 7, true>(currKnot.Minv);
    
    Ju.setZero(14, 7);
    Ju.bottomRows<7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(currKnot.Minv);
}

void iLQRSolverOptimized::setGravity(double g_z)
{
    g_plancher_gravity = g_z;
}
