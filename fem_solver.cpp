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
std::vector<double> K_global;//全局刚度矩阵，为了兼容库Eigen，没有使用二维数组
std::vector<double> b;//载荷向量
std::vector<double> u;//位移向量
std::vector<std::vector<int>> T;//单元编号、局部编号
std::vector<double> P;//全局编号与对应节点坐标
std::vector<BoundaryCondition> boundary_conditions;

// --- 内部辅助函数声明 (仅在此文件内使用) ---
double calculateStiffnessEntry(int e, int alpha, int beta, int numGaussPoints);//alpha是试探函数的索引，beta是测试（加权）函数的索引
double calculateLoadEntry(int e, int beta, int numGaussPoints);
void computeGeometricQuantities_1D(
    const std::vector<double>& element_coords, double x_ref,
    const std::vector<double>& dN_dx_ref_trial, const std::vector<double>& dN_dx_ref_test,
    double& jacobianDet, std::vector<double>& dN_dx_trial, std::vector<double>& dN_dx_test,
    double& x_physical);

// --- FEM 函数实现 ---

//生成或者导入网格
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

//组装全局刚度矩阵和载荷向量
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

//施加边界条件（在全局刚度矩阵和载荷向量上操作）
void applyBoundaryConditions() {
    // --- 第一步：处理所有罗宾和诺伊曼边界条件 ---
    // 这些条件只是对原始矩阵进行“加法”修正
    for (const auto& bc : boundary_conditions) {
        if (bc.K_bc == 1) { // 只处理 K_bc != 0 的情况
            int idx = bc.node_index;
            double s = 1.0; // 符号位

            // 判断是左端点还是右端点
            if (idx == 0) {
                s = -1.0; // 左端点，符号为负
            } else if (idx == N - 1) {
                s = 1.0;  // 右端点，符号为正
            } else {
                // 对于内部节点，理论上不应该有自然边界条件，但此处作为容错
                continue; 
            }

            K_global[idx * N + idx] += s * bc.L_bc;
            b[idx] += s * bc.q_bc;
        }
    }

    // --- 第二步：处理所有狄利克雷边界条件 ---
    // 狄利克雷条件会“覆盖”矩阵的某些部分，因此必须在最后处理
    for (const auto& bc : boundary_conditions) {
        if (bc.K_bc==0) { // 只处理 K_bc == 0 的情况
            int idx = bc.node_index;

            if (std::fabs(bc.L_bc) < 1e-9) {
                continue; // 无效条件
            }
            
            double val = bc.q_bc / bc.L_bc;
            
            // 修正载荷向量 b
            for (int j = 0; j < N; ++j) {
                if (j != idx) {
                    b[j] -= K_global[j * N + idx] * val;
                }
            }
            
            // 修正刚度矩阵 K (划零法)
            for (int j = 0; j < N; ++j) {
                K_global[idx * N + j] = 0.0;
                K_global[j * N + idx] = 0.0;
            }
            K_global[idx * N + idx] = 1.0;
            b[idx] = val;
        }
    }
}

//求解线性方程组
void solveLinearSystem() {
    Eigen::Map<Eigen::MatrixXd> K_eigen(K_global.data(), N, N);
    Eigen::Map<Eigen::VectorXd> b_eigen(b.data(), N);
    Eigen::VectorXd u_eigen = K_eigen.lu().solve(b_eigen);
    for(int i = 0; i < N; ++i) u[i] = u_eigen(i);
}

//输出结果
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
        double x_ref = gaussPoints[gp], weight = gaussWeights[gp];//这个地方得到的点是在【-1, 1】参考区间内的

        std::vector<double> dN_dx_ref_trial_gp(n), dN_dx_ref_test_gp(n);
        for (int i = 0; i < n; ++i) {
            dN_dx_ref_trial_gp[i] = shapeFunctionDerivative_trial(i, x_ref);//得到参考试探函数alpha=i的导数在x_ref处的值
            dN_dx_ref_test_gp[i] = shapeFunctionDerivative_test(i, x_ref);//得到参考测试函数beta=i的导数在x_ref处的值
        }
        
        double jacobianDet_e;
        std::vector<double> dN_dx_trial_gp(n), dN_dx_test_gp(n);
        double x_gp;
        computeGeometricQuantities_1D(element_coords_e, x_ref, dN_dx_ref_trial_gp, dN_dx_ref_test_gp, jacobianDet_e, dN_dx_trial_gp, dN_dx_test_gp, x_gp);
        
        const double c_x = coefficient_c(x_gp);
        
        double integrand = c_x * dN_dx_test_gp[beta] * dN_dx_trial_gp[alpha];
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
        double x_ref = gaussPoints[gp], weight = gaussWeights[gp];
        
        // 获取参考坐标系的导数（虽然载荷向量不需要导数，但为了使用统一的几何变换函数）
        std::vector<double> dN_dx_ref_trial_gp(n), dN_dx_ref_test_gp(n);
        for (int i = 0; i < n; ++i) {
            dN_dx_ref_trial_gp[i] = shapeFunctionDerivative_trial(i, x_ref);
            dN_dx_ref_test_gp[i] = shapeFunctionDerivative_test(i, x_ref);
        }
        
        // 使用统一的几何变换函数
        double jacobianDet_e;
        std::vector<double> dN_dx_trial_gp(n), dN_dx_test_gp(n);
        double x_gp;
        computeGeometricQuantities_1D(element_coords_e, x_ref, dN_dx_ref_trial_gp, dN_dx_ref_test_gp, jacobianDet_e, dN_dx_trial_gp, dN_dx_test_gp, x_gp);

        double N_ref_test_beta = shapeFunction_test(beta, x_ref);  // 参考测试函数在参考坐标处的值,和真实测试函数在物理坐标处的值一样
        const double f_x = source_term_f(x_gp);              // 源项在物理坐标处的值
        double integrand = f_x * N_ref_test_beta; 
        entryValue += integrand * jacobianDet_e * weight;
    }
    return entryValue;
}

//进行仿射变换
void computeGeometricQuantities_1D(
    const std::vector<double>& element_coords, double x_ref,
    const std::vector<double>& dN_dx_ref_trial, const std::vector<double>& dN_dx_ref_test,
    double& jacobianDet, std::vector<double>& dN_dx_trial, std::vector<double>& dN_dx_test,
    double& x_physical)
{
    // 计算雅可比行列式
    jacobianDet = 0.0;
    for (int i = 0; i < n; ++i) jacobianDet += dN_dx_ref_trial[i] * element_coords[i];
    
    // 计算物理坐标系下的导数
    for (int i = 0; i < n; ++i) {
        dN_dx_trial[i] = dN_dx_ref_trial[i] / jacobianDet;
        dN_dx_test[i] = dN_dx_ref_test[i] / jacobianDet;
    }
    
    // 计算物理坐标：参考坐标映射到物理坐标
    x_physical = 0.0;
    for (int i = 0; i < n; ++i) {
        x_physical += shapeFunction_trial(i, x_ref) * element_coords[i];
    }
}
