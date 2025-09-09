#include "fem_solver_2d.h"
#include "gauss_quadrature_2d.h"
#include "shape_functions_2d.h"
#include "geometry_mapping_2d.h"
#include "problem_definition.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

// --- 全局变量定义 ---
int N, M;
const int n = 3; // 三角形单元有3个节点
Eigen::SparseMatrix<double> K_global;
Eigen::VectorXd b;
Eigen::VectorXd u;
std::vector<std::vector<int>> T;
std::vector<double> P;
std::vector<std::vector<int>> boundary_edges;
std::vector<BoundaryCondition> boundary_conditions;

// --- 内部辅助函数声明 ---
double calculateStiffnessEntry2D(int e, int alpha, int beta, int numGaussPoints);
double calculateLoadEntry2D(int e, int beta, int numGaussPoints);

// --- 网格生成函数实现 ---
void generateMesh2D(const Domain& domain, int N1, int N2, const BoundaryTypes& bc_types) {
    // 清空旧数据
    P.clear();
    T.clear();
    boundary_edges.clear();

    if (N1 < 1 || N2 < 1) {
        std::cerr << "Interval counts N1 and N2 must be at least 1." << std::endl;
        return;
    }

    // 生成节点坐标 (P vector)
    const int Nx = N1 + 1;
    const int Ny = N2 + 1;
    N = Nx * Ny;  // 设置全局节点数
    P.reserve(N * 2);
    const double hx = (domain.right - domain.left) / N1;
    const double hy = (domain.top - domain.bottom) / N2;

    // Column-Major order: 外循环为列(x)，内循环为行(y)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            P.push_back(domain.left + i * hx);
            P.push_back(domain.bottom + j * hy);
        }
    }

    // 生成三角形单元 (T vector)
    M = 2 * N1 * N2;  // 设置全局单元数
    T.reserve(M);
    auto node_2d_to_1d = [&](int i, int j) { return i * Ny + j; };

    // Column-Major order for elements
    for (int ie = 0; ie < N1; ++ie) {
        for (int je = 0; je < N2; ++je) {
            int idx_BL = node_2d_to_1d(ie, je);
            int idx_BR = node_2d_to_1d(ie + 1, je);
            int idx_TL = node_2d_to_1d(ie, je + 1);
            int idx_TR = node_2d_to_1d(ie + 1, je + 1);

            // Lower triangle (直角在左下角): 第一个局部编号在左下角 BL
            T.push_back({idx_BL, idx_BR, idx_TL}); 
            // Upper triangle (直角在右上角): 第一个局部编号在左上角 TL
            T.push_back({idx_TL, idx_BR, idx_TR}); 
        }
    }

    // 生成边界边（简化版本，不完整实现）
    boundary_edges.reserve(2 * N1 + 2 * N2);
    auto element_2d_to_1d = [&](int ie, int je, bool is_upper) {
        return ie * (2 * N2) + 2 * je + (is_upper ? 1 : 0);
    };

    // 这里可以根据需要添加完整的边界边生成逻辑
    // 目前简化处理
}

// --- FEM 函数实现 ---

void preprocess2D() {
    // 定义计算域和网格参数
    Domain domain = {0.0, 1.0, 0.0, 1.0}; // 单位正方形
    int N1 = 4, N2 = 4; // 网格划分
    BoundaryTypes bc_types = {-1, -2, -1, -2}; // 边界条件类型
    
    // 生成网格
    generateMesh2D(domain, N1, N2, bc_types);
    
    // 初始化矩阵和向量
    K_global.resize(N, N);
    K_global.reserve(Eigen::VectorXi::Constant(N, 9)); // 每行预计有9个非零元素（2D情况）
    b = Eigen::VectorXd::Zero(N);
    u = Eigen::VectorXd::Zero(N);
    
    std::cout << "网格生成完成：" << N << " 个节点，" << M << " 个单元" << std::endl;
}

void assemble2D() {
    // 使用三元组列表收集所有元素
    typedef Eigen::Triplet<double> T_entry;
    std::vector<T_entry> tripletList;
    tripletList.reserve(9 * M); // 每个单元贡献9个矩阵元素
    
    for (int e = 0; e < M; ++e) {
        const int numGaussPoints = 3; // 使用3点高斯积分
        for (int alpha = 0; alpha < n; ++alpha) {
            for (int beta = 0; beta < n; ++beta) {
                double K_e_val = calculateStiffnessEntry2D(e, alpha, beta, numGaussPoints);
                int global_row = T[e][beta];
                int global_col = T[e][alpha];
                tripletList.push_back(T_entry(global_row, global_col, K_e_val));
            }
        }
        for (int beta = 0; beta < n; ++beta) {
            double b_e_val = calculateLoadEntry2D(e, beta, numGaussPoints);
            int global_row = T[e][beta];
            b(global_row) += b_e_val;
        }
    }
    
    // 从三元组列表构建稀疏矩阵
    K_global.setFromTriplets(tripletList.begin(), tripletList.end());
    K_global.makeCompressed();
}

void applyBoundaryConditions2D() {
    // 简化的边界条件处理
    // 对于单位正方形，假设边界上的节点具有零值狄利克雷边界条件
    
    for (int i = 0; i < N; ++i) {
        double x = P[i * 2];
        double y = P[i * 2 + 1];
        
        // 检查是否在边界上
        bool is_boundary = (std::abs(x) < 1e-10) || (std::abs(x - 1.0) < 1e-10) || 
                          (std::abs(y) < 1e-10) || (std::abs(y - 1.0) < 1e-10);
        
        if (is_boundary) {
            // 应用零值狄利克雷边界条件
            // 修正载荷向量
            for (int j = 0; j < N; ++j) {
                if (j != i) {
                    b(j) -= K_global.coeff(j, i) * 0.0; // 边界值为0
                }
            }
            
            // 修改矩阵行和列
            K_global.makeCompressed();
            for (int j = 0; j < N; ++j) {
                if (j != i) {
                    K_global.coeffRef(i, j) = 0.0;
                    K_global.coeffRef(j, i) = 0.0;
                }
            }
            K_global.coeffRef(i, i) = 1.0;
            b(i) = 0.0;
        }
    }
}

void solveLinearSystem2D() {
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

void postprocess2D() {
    using namespace std;
    cout << "\n--- 二维有限元计算结果 ---" << endl;
    cout << left << setw(10) << "节点ID"
         << left << setw(15) << "坐标 (x)"
         << left << setw(15) << "坐标 (y)"
         << left << setw(20) << "FEM解 (u)" << endl;
    cout << string(60, '-') << endl;

    for (int i = 0; i < N; ++i) {
        double x = P[i * 2];
        double y = P[i * 2 + 1];
        
        cout << left << setw(10) << i
             << left << setw(15) << fixed << setprecision(4) << x
             << left << setw(15) << fixed << setprecision(4) << y
             << left << setw(20) << scientific << setprecision(5) << u(i)
             << endl;
    }
}

// --- 内部辅助函数定义 ---

double calculateStiffnessEntry2D(int e, int alpha, int beta, int numGaussPoints) {
    // 创建几何映射对象
    GeometryMapping2D mapping(e);
    
    std::vector<double> gaussPoints, gaussWeights;
    getGaussPointsTriangle(numGaussPoints, gaussPoints, gaussWeights);
    
    double entryValue = 0.0;
    for (int gp = 0; gp < numGaussPoints; ++gp) {
        double x_ref = gaussPoints[gp * 2];
        double y_ref = gaussPoints[gp * 2 + 1];
        double weight = gaussWeights[gp];
        
        // 获取物理坐标
        double x_gp, y_gp;
        mapping.mapToPhysical(x_ref, y_ref, x_gp, y_gp);
        
        // 计算形函数在参考坐标系下的导数
        double dN_dx_ref_trial = shapeFunctionDerivativeXRef_trial(alpha, x_ref, y_ref);
        double dN_dy_ref_trial = shapeFunctionDerivativeYRef_trial(alpha, x_ref, y_ref);
        double dN_dx_ref_test = shapeFunctionDerivativeXRef_test(beta, x_ref, y_ref);
        double dN_dy_ref_test = shapeFunctionDerivativeYRef_test(beta, x_ref, y_ref);
        
        // 转换为物理坐标系下的导数
        double dN_dx_trial, dN_dy_trial, dN_dx_test, dN_dy_test;
        mapping.transformGradient(dN_dx_ref_trial, dN_dy_ref_trial, dN_dx_trial, dN_dy_trial);
        mapping.transformGradient(dN_dx_ref_test, dN_dy_ref_test, dN_dx_test, dN_dy_test);
        
        // 计算积分被积函数（假设扩散系数为1）
        double integrand = (dN_dx_test * dN_dx_trial + dN_dy_test * dN_dy_trial);
        
        entryValue += integrand * mapping.getJacobianDeterminant() * weight;
    }
    return entryValue;
}

double calculateLoadEntry2D(int e, int beta, int numGaussPoints) {
    // 创建几何映射对象
    GeometryMapping2D mapping(e);
    
    std::vector<double> gaussPoints, gaussWeights;
    getGaussPointsTriangle(numGaussPoints, gaussPoints, gaussWeights);
    
    double entryValue = 0.0;
    for (int gp = 0; gp < numGaussPoints; ++gp) {
        double x_ref = gaussPoints[gp * 2];
        double y_ref = gaussPoints[gp * 2 + 1];
        double weight = gaussWeights[gp];
        
        // 获取物理坐标
        double x_gp, y_gp;
        mapping.mapToPhysical(x_ref, y_ref, x_gp, y_gp);
        
        // 计算形函数值
        double N_test_beta = shapeFunction2D_test(beta, x_ref, y_ref);
        
        // 计算源项（假设 f = 1.0）
        const double f_xy = 1.0;
        
        double integrand = f_xy * N_test_beta; 
        entryValue += integrand * mapping.getJacobianDeterminant() * weight;
    }
    return entryValue;
}
