#ifndef FEM_SOLVER_2D_H
#define FEM_SOLVER_2D_H

#include <Eigen/Sparse>
#include <vector>
#include "config.h"

// --- 全局变量与数据结构 ---
// 使用 extern 关键字声明全局变量，定义将在 fem_solver_2d.cpp 中
extern int N, M;                              // N: 节点数, M: 单元数
extern const int n;                           // 每个单元的节点数 (三角形单元为3)
extern Eigen::SparseMatrix<double> K_global;  // 全局刚度矩阵
extern Eigen::VectorXd b;                     // 载荷向量
extern Eigen::VectorXd u;                     // 解向量
extern std::vector<std::vector<int>> T;       // 单元连接矩阵 [element][local_node] = global_node
extern std::vector<double> P;                 // 节点坐标 [x0, y0, x1, y1, ...]
extern std::vector<Config::Boundary> boundary_edges;  // 边界边信息

// --- 核心 FEM 函数声明 ---

/**
 * @brief 预处理步骤：生成网格、初始化矩阵和向量
 */
void preprocess2D();

/**
 * @brief 组装全局刚度矩阵和载荷向量
 */
void assemble2D();

/**
 * @brief 施加边界条件
 */
void applyBoundaryConditions2D();

/**
 * @brief 求解线性方程组 K_global * u = b
 */
void solveLinearSystem2D(const std::string& solver_type = "CG",
                         const std::string& preconditioner_type = "DiagonalPreconditioner",
                         double tol = 1e-8, int max_iter = 1000, bool verbose = true);

/**
 * @brief 后处理步骤：输出结果并计算误差
 */
void postprocess2D();

#endif  // FEM_SOLVER_2D_H
