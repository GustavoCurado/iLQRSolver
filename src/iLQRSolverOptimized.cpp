#include "iLQRSolverOptimized.hpp"

// Variável estática (isolada neste arquivo) para controlar a gravidade
static double g_plancher_gravity = 9.81;

// Enganamos o código do Plancher para usar a nossa variável dinâmica
#define GRAVITY g_plancher_gravity

#include "../../fast-rbd-gradients/helpers_CPU/dynamicsGradient.h"

iLQRSolverOptimized::iLQRSolverOptimized(double h, int Nt) : iLQRSolverBase(h, Nt)
{
    this->initializeMemory();
}

void iLQRSolverOptimized::initializeMemory()
{
    // Kuka LBR iiwa 14 tem 7 juntas
    this->Nx = 14;
    this->Nu = 7;
    
    iLQRSolverBase::initializeMemory(); // Chama a inicialização base (matrizes Q, R, xtraj, etc.)

    // Inicializa as transformações e inércias "hardcoded" do Plancher
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
    
    double qdd[7];
    double Minv[49];
    
    // Chamada do Plancher (Templates: T, MPC_MODE=false, VEL_DAMPING=false, DENSE_MINV=false)
    forwardDynamics<double, false, false, false>(qdd, tau, v, q, Minv);
    
    Eigen::VectorXd xdot = Eigen::VectorXd::Zero(14);
    xdot.head<7>() = x.tail<7>();
    
    // O Map não copia memória, ele lê o array direto pro Eigen
    xdot.tail<7>() = Eigen::Map<Eigen::Matrix<double, 7, 1>>(qdd);
    
    return xdot;
}

void iLQRSolverOptimized::getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju)
{
    knot<double> currKnot;
    for(int i = 0; i < 7; i++) {
        currKnot.q[i] = x[i];
        currKnot.qd[i] = x[i + 7];
        currKnot.u[i] = u[i];
    }
    
    double dqdd_array[98];
    
    // Chamada do gradiente do Plancher (Templates: T, QDD_MINV_PASSED_IN=false, MPC_MODE=false, VEL_DAMPING=false, EXTERNAL_TI=false)
    forwardDynamicsGradient<double, false, false, false, false>(dqdd_array, &currKnot, this->d_T, this->d_I);
    
    Jx.setZero(14, 14);
    Jx.topRightCorner<7, 7>().setIdentity();
    
    // O código do Plancher e o math.h lidam com as matrizes no formato Column-Major (que é o padrão C++ para simulações assim).
    // Como no seu iLQRSolverBase.hpp o Eigen está setado para Row-Major (#define EIGEN_DEFAULT_TO_ROW_MAJOR),
    // nós precisamos usar 'Eigen::ColMajor' no Map para que a tradução entre as memórias seja perfeita!
    Jx.bottomLeftCorner<7, 7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(&dqdd_array[0]);
    Jx.bottomRightCorner<7, 7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(&dqdd_array[49]);
    
    Ju.setZero(14, 7);
    Ju.bottomRows<7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(currKnot.Minv);
}

void iLQRSolverOptimized::setGravity(double g_z)
{
    // Apenas alteramos a variável estática, e o próximo frame matemático
    // executado pelo Plancher irá ler este novo valor!
    g_plancher_gravity = g_z;
}
