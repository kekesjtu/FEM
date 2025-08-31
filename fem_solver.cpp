#include "fem_solver.h"
#include "gauss_quadrature.h"
#include "shape_functions.h"
#include "problem_definition.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

// --- 全局变量定义 ---
int N, M;
const int n = 2; // 线性单元
std::vector<double> K_global;
std::vector<double> b;
std::vector<double> u;
std::vector<std::vector<int>> T;
std::vector<double> P;
std::vector<BoundaryCondition> boundary_conditions;

// --- 内部辅助函数声明 (仅在此文件内使用) ---
double calculateStiffnessEntry(int e, int alpha, int beta, int numGaussPoints);
double calculateLoadEntry(int e, int beta, int numGaussPoints);
void computeGeometricQuantities_1D(
    const std::vector<double>& element_coords, double xi,
    const std::vector<double>& dN_dxi_trial, const std::vector<double>& dN_dxi_test,
    double& jacobianDet, std::vector<double>& dN_dx_trial, std::vector<double>& dN_dx_test);

// --- FEM 函数实现 ---

void preprocess() {
    defineProblem(N, M, boundary_conditions);
    K_global.assign(N * N, 0.0);
    b.assign(N, 0.0);
    u.assign(N, 0.0);
    T.assign(M, std::vector<int>(n));
    P.assign(N, 0.0);
    const double domain_length = 1.0;
    for (int i = 0; i < N; ++i) P[i] = domain_length * i / (N - 1);
    for (int e = 0; e < M; ++e) { T[e][0] = e; T[e][1] = e + 1; }
}

void assemble() {
    for (int e = 0; e < M; ++e) {
        const int numGaussPoints_e = 3;
        for (int alpha = 0; alpha < n; ++alpha) {
            for (int beta = 0; beta < n; ++beta) {
                double K_e_val = calculateStiffnessEntry(e, alpha, beta, numGaussPoints_e);
                int global_row = T[e][beta];
                int global_col = T[e][alpha];
                K_global[global_row * N + global_col] += K_e_val;
            }
        }
        for (int beta = 0; beta < n; ++beta) {
            double b_e_val = calculateLoadEntry(e, beta, numGaussPoints_e);
            int global_row = T[e][beta];
            b[global_row] += b_e_val;
        }
    }
}

void applyBoundaryConditions() {
    for (const auto& bc : boundary_conditions) {
        int idx = bc.node_index;
        double val = bc.value;
        for (int j = 0; j < N; ++j) {
            if (j != idx) {
                b[j] -= K_global[j * N + idx] * val;
            }
        }
        for (int j = 0; j < N; ++j) {
            K_global[idx * N + j] = 0.0;
            K_global[j * N + idx] = 0.0;
        }
        K_global[idx * N + idx] = 1.0;
        b[idx] = val;
    }
}

void solveLinearSystem() {
    Eigen::Map<Eigen::MatrixXd> K_eigen(K_global.data(), N, N);
    Eigen::Map<Eigen::VectorXd> b_eigen(b.data(), N);
    Eigen::VectorXd u_eigen = K_eigen.lu().solve(b_eigen);
    for(int i = 0; i < N; ++i) u[i] = u_eigen(i);
}

void postprocess() {
    using namespace std;
    cout << "\n--- 计算结果 ---" << endl;
    cout << left << setw(10) << "节点ID"
         << left << setw(15) << "坐标 (x)"
         << left << setw(20) << "FEM解 (u)"
         << left << setw(20) << "精确解"
         << left << setw(20) << "误差" << endl;
    cout << string(85, '-') << endl;

    for (int i = 0; i < N; ++i) {
        double x = P[i];
        double exact_sol = exact_solution_u(x);
        double error = abs(u[i] - exact_sol);
        
        cout << left << setw(10) << i
             << left << setw(15) << fixed << setprecision(4) << x
             << left << setw(20) << scientific << setprecision(5) << u[i]
             << left << setw(20) << scientific << setprecision(5) << exact_sol
             << left << setw(20) << scientific << setprecision(5) << error
             << endl;
    }
}

// --- 内部辅助函数定义 ---

double calculateStiffnessEntry(int e, int alpha, int beta, int numGaussPoints) {
    std::vector<double> element_coords_e(n);
    for (int i = 0; i < n; ++i) element_coords_e[i] = P[T[e][i]];
    
    std::vector<double> gaussPoints, gaussWeights;
    getGaussPoints(numGaussPoints, gaussPoints, gaussWeights);
    
    double entryValue = 0.0;
    for (int gp = 0; gp < numGaussPoints; ++gp) {
        double xi = gaussPoints[gp], weight = gaussWeights[gp];
        
        std::vector<double> dN_dxi_trial_gp(n), dN_dxi_test_gp(n);
        for (int i = 0; i < n; ++i) {
            dN_dxi_trial_gp[i] = shapeFunctionDerivative_trial(i, xi);
            dN_dxi_test_gp[i] = shapeFunctionDerivative_test(i, xi);
        }
        
        double jacobianDet_e;
        std::vector<double> dN_dx_trial_gp(n), dN_dx_test_gp(n);
        computeGeometricQuantities_1D(element_coords_e, xi, dN_dxi_trial_gp, dN_dxi_test_gp, jacobianDet_e, dN_dx_trial_gp, dN_dx_test_gp);
        
        double dN_dx_alpha = dN_dx_trial_gp[alpha];
        double dN_dx_beta = dN_dx_test_gp[beta];

        double x_gp = 0.0;
        for (int i = 0; i < n; ++i) x_gp += shapeFunction_trial(i, xi) * element_coords_e[i];
        
        const double c_x = coefficient_c(x_gp);
        
        double integrand = c_x * dN_dx_beta * dN_dx_alpha;
        entryValue += integrand * jacobianDet_e * weight;
    }
    return entryValue;
}

double calculateLoadEntry(int e, int beta, int numGaussPoints) {
    std::vector<double> element_coords_e(n);
    for(int i=0; i<n; ++i) element_coords_e[i] = P[T[e][i]];

    std::vector<double> gaussPoints, gaussWeights;
    getGaussPoints(numGaussPoints, gaussPoints, gaussWeights);
    
    double entryValue = 0.0;
    for (int gp = 0; gp < numGaussPoints; ++gp) {
        double xi = gaussPoints[gp], weight = gaussWeights[gp];
        double N_beta = shapeFunction_test(beta, xi);
        double jacobianDet_e = (element_coords_e[1] - element_coords_e[0]) / 2.0;
        
        double x_gp = 0.0;
        for (int i = 0; i < n; ++i) x_gp += shapeFunction_trial(i, xi) * element_coords_e[i];
        
        const double f_x = source_term_f(x_gp);
        
        double integrand = f_x * N_beta;
        entryValue += integrand * jacobianDet_e * weight;
    }
    return entryValue;
}

void computeGeometricQuantities_1D(
    const std::vector<double>& element_coords, double xi,
    const std::vector<double>& dN_dxi_trial, const std::vector<double>& dN_dxi_test,
    double& jacobianDet, std::vector<double>& dN_dx_trial, std::vector<double>& dN_dx_test)
{
    jacobianDet = 0.0;
    for (int i = 0; i < n; ++i) jacobianDet += dN_dxi_trial[i] * element_coords[i];
    for (int i = 0; i < n; ++i) {
        dN_dx_trial[i] = dN_dxi_trial[i] / jacobianDet;
        dN_dx_test[i] = dN_dxi_test[i] / jacobianDet;
    }
}
