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
    
    double qdd[7] = {0.0};
    double Minv[49] = {0.0};
    
    // Chamada do Plancher (Templates: T, MPC_MODE=false, VEL_DAMPING=true, DENSE_MINV=false)
    forwardDynamics<double, false, true, false>(qdd, tau, v, q, Minv);
    
    Eigen::VectorXd xdot = Eigen::VectorXd::Zero(14);
    xdot.head<7>() = x.tail<7>();
    
    // O Map não copia memória, ele lê o array direto pro Eigen
    xdot.tail<7>() = Eigen::Map<Eigen::Matrix<double, 7, 1>>(qdd);
    
    return xdot;
}

void iLQRSolverOptimized::getContinuousJacs(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& Jx, Eigen::MatrixXd& Ju)
{
    // É extremamente importante inicializar com zeros, pois o Plancher lê qdd
    // internamente durante o Forward Pass para calcular o bias. Se não inicializarmos,
    // ele vai ler lixo de memória da stack e explodir a simulação inteira!
    knot<double> currKnot = {0};
    for(int i = 0; i < 7; i++) {
        currKnot.q[i] = x[i];
        currKnot.qd[i] = x[i + 7];
        currKnot.u[i] = u[i];
    }
    
    double dqdd_array[98];
    
    // Chamada do gradiente do Plancher (Templates: T, QDD_MINV_PASSED_IN=false, MPC_MODE=false, VEL_DAMPING=true, EXTERNAL_TI=false)
    forwardDynamicsGradient<double, false, false, true, false>(dqdd_array, &currKnot, this->d_T, this->d_I);
    
    Jx.setZero(14, 14);
    Jx.topRightCorner<7, 7>().setIdentity();
    
    // O código do Plancher gera o array em formato Column-Major nativamente.
    Jx.bottomLeftCorner<7, 7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(&dqdd_array[0]);
    Jx.bottomRightCorner<7, 7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(&dqdd_array[49]);
    
    // O código do Plancher gera o M inverso apenas como uma matriz triangular superior
    // para economizar cálculos. Precisamos espelhar a matriz para preencher a parte de baixo!
    triangularToDense<double, 7, true>(currKnot.Minv);
    
    Ju.setZero(14, 7);
    Ju.bottomRows<7>() = Eigen::Map<Eigen::Matrix<double, 7, 7, Eigen::ColMajor>>(currKnot.Minv);
}

void iLQRSolverOptimized::setGravity(double g_z)
{
    // Apenas alteramos a variável estática, e o próximo frame matemático
    // executado pelo Plancher irá ler este novo valor!
    g_plancher_gravity = g_z;
}
