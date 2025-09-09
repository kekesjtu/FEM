#include "fem_solver.h"
#include "gauss_quadrature.h"
#include "shape_functions.h"
#include "problem_definition.h"
#include "error_analysis.h"
#include "geometry_mapping.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

// --- 全局变量定义 ---
int N, M;
const int n = 2; // 线性单元
Eigen::SparseMatrix<double> K_global;//全局刚度矩阵，使用稀疏矩阵存储
Eigen::VectorXd b;//载荷向量
Eigen::VectorXd u;//解向量
std::vector<std::vector<int>> T;//单元编号、局部编号
std::vector<double> P;//全局编号与对应节点坐标
std::vector<BoundaryCondition> boundary_conditions;

// --- 内部辅助函数声明 (仅在此文件内使用) ---
double calculateStiffnessEntry(int e, int alpha, int beta, int numGaussPoints);//alpha是试探函数的索引，beta是测试（加权）函数的索引
double calculateLoadEntry(int e, int beta, int numGaussPoints);

// --- FEM 函数实现 ---

//生成或者导入网格
void preprocess() {
    defineProblem(N, M, boundary_conditions);
    K_global.resize(N, N);  // 创建NxN的稀疏矩阵
    K_global.reserve(Eigen::VectorXi::Constant(N, 3));  // 每行预计有3个非零元素（1D情况）
    b = Eigen::VectorXd::Zero(N);  // 零向量
    u = Eigen::VectorXd::Zero(N);  // 零向量
    T.assign(M, std::vector<int>(n));
    P.assign(N, 0.0);
    const double domain_length = 1.0;
    for (int i = 0; i < N; ++i) P[i] = domain_length * i / (N - 1);
    for (int e = 0; e < M; ++e) { T[e][0] = e; T[e][1] = e + 1; }
}

//组装全局刚度矩阵和载荷向量
void assemble() {
    // 使用三元组列表收集所有元素
    typedef Eigen::Triplet<double> T_entry;
    std::vector<T_entry> tripletList;
    tripletList.reserve(3 * N); // 预留空间
    
    for (int e = 0; e < M; ++e) {
        const int numGaussPoints = 3;
        for (int alpha = 0; alpha < n; ++alpha) {
            for (int beta = 0; beta < n; ++beta) {
                double K_e_val = calculateStiffnessEntry(e, alpha, beta, numGaussPoints);
                int global_row = T[e][beta];
                int global_col = T[e][alpha];
                tripletList.push_back(T_entry(global_row, global_col, K_e_val));
            }
        }
        for (int beta = 0; beta < n; ++beta) {
            double b_e_val = calculateLoadEntry(e, beta, numGaussPoints);
            int global_row = T[e][beta];
            b(global_row) += b_e_val;
        }
    }
    
    // 从三元组列表构建稀疏矩阵
    K_global.setFromTriplets(tripletList.begin(), tripletList.end());
    // 对具有相同索引的元素进行求和
    K_global.makeCompressed();
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

            // 对于稀疏矩阵，使用coeffRef来修改元素
            K_global.coeffRef(idx, idx) += s * bc.L_bc;
            b(idx) += s * bc.q_bc;
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
                    b(j) -= K_global.coeff(j, idx) * val;
                }
            }
            
            // 转换为压缩格式以便修改
            K_global.makeCompressed();
            
            // 将整行和整列设置为零 (稀疏矩阵需要逐个元素处理)
            for (int j = 0; j < N; ++j) {
                if (j != idx) {
                    K_global.coeffRef(idx, j) = 0.0;
                    K_global.coeffRef(j, idx) = 0.0;
                }
            }
            K_global.coeffRef(idx, idx) = 1.0;
            b(idx) = val;
        }
    }
}

//求解线性方程组
void solveLinearSystem() {
    // 使用SparseLU求解器求解稀疏线性系统
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(K_global);
    solver.factorize(K_global);
    if(solver.info() != Eigen::Success) {
        std::cerr << "分解失败！" << std::endl;
        return;
    }
    u = solver.solve(b);
    if(solver.info() != Eigen::Success) {
        std::cerr << "求解失败！" << std::endl;
        return;
    }
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
        double error = abs(u(i) - exact_sol);
        
        cout << left << setw(10) << i
             << left << setw(15) << fixed << setprecision(4) << x
             << left << setw(20) << scientific << setprecision(5) << u(i)
             << left << setw(20) << scientific << setprecision(5) << exact_sol
             << left << setw(20) << scientific << setprecision(5) << error
             << endl;
    }
    
    // --- 误差分析 ---
    cout << "\n--- 误差分析 ---" << endl;
    
    // 计算各种误差范数
    double max_error = computeMaxError(u, exact_solution_u);
    double l2_error = computeL2Error(u, exact_solution_u);
    double h1_error = computeH1Error(u, exact_solution_du_dx);
    
    cout << "误差范数:" << endl;
    cout << "  最大误差 (L∞范数): " << scientific << setprecision(6) << max_error << endl;
    cout << "  L2范数误差:        " << scientific << setprecision(6) << l2_error << endl;
    cout << "  H1范数误差:        " << scientific << setprecision(6) << h1_error << endl;
}

// --- 内部辅助函数定义 ---

double calculateStiffnessEntry(int e, int alpha, int beta, int numGaussPoints) {
    // 创建几何映射对象
    GeometryMapping1D mapping(e);  // 直接使用单元索引构造
    
    std::vector<double> gaussPoints, gaussWeights;
    getGaussPoints(numGaussPoints, gaussPoints, gaussWeights);
    
    double entryValue = 0.0;
    for (int gp = 0; gp < numGaussPoints; ++gp) {
        double x_ref_gp = gaussPoints[gp];
        double weight = gaussWeights[gp];
        
        // 获取物理坐标和雅可比
        double x_gp = mapping.mapToPhysical(x_ref_gp);
        double jacobian = mapping.getJacobian();
        
        // 计算形函数导数
        double dN_dx_ref_trial = shapeFunctionDerivative_trial(alpha, x_ref_gp);
        double dN_dx_ref_test = shapeFunctionDerivative_test(beta, x_ref_gp);
        
        // 转换为物理坐标系下的导数
        double dN_dx_trial = mapping.transformDerivative(dN_dx_ref_trial);
        double dN_dx_test = mapping.transformDerivative(dN_dx_ref_test);
        
        // 计算积分被积函数
        const double c_x = coefficient_c(x_gp);
        double integrand = c_x * dN_dx_test * dN_dx_trial;

        entryValue += integrand * jacobian * weight;
    }
    return entryValue;
}

double calculateLoadEntry(int e, int beta, int numGaussPoints) {
    // 创建几何映射对象
    GeometryMapping1D mapping(e);  // 直接使用单元索引构造

    std::vector<double> gaussPoints, gaussWeights;
    getGaussPoints(numGaussPoints, gaussPoints, gaussWeights);
    
    double entryValue = 0.0;
    for (int gp = 0; gp < numGaussPoints; ++gp) {
        double x_ref_gp = gaussPoints[gp];
        double weight = gaussWeights[gp];
        
        // 获取物理坐标和雅可比
        double x_gp = mapping.mapToPhysical(x_ref_gp);
        double jacobian = mapping.getJacobian();

        // 计算形函数值（载荷向量只需要形函数值，不需要导数）,而参考坐标系下的形函数值与物理坐标系下的形函数值相同
        double N_test_beta = shapeFunction_test(beta, x_ref_gp);
        
        // 计算源项
        const double f_x = source_term_f(x_gp);
        
        double integrand = f_x * N_test_beta; 
        entryValue += integrand * jacobian * weight;
    }
    return entryValue;
}
