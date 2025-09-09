#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <vector>
#include "problem_definition.h"
#include <Eigen/Sparse>

// --- 全局变量与数据结构 ---
// 使用 extern 关键字声明全局变量，定义将在 fem_solver.cpp 中
extern int N, M;
extern const int n; // 线性单元
extern Eigen::SparseMatrix<double> K_global; // 使用稀疏矩阵代替稠密矩阵
extern Eigen::VectorXd b; // 使用Eigen向量
extern Eigen::VectorXd u; // 使用Eigen向量
extern std::vector<std::vector<int>> T;
extern std::vector<double> P;
extern std::vector<BoundaryCondition> boundary_conditions;

// --- 核心 FEM 函数声明 ---

/**
 * @brief 预处理步骤：初始化网格、矩阵和向量
 */
void preprocess();

/**
 * @brief 组装全局刚度矩阵和载荷向量
 */
void assemble();

/**
 * @brief 施加狄利克雷边界条件
 */
void applyBoundaryConditions();

/**
 * @brief 求解线性方程组 K_global * u = b
 */
void solveLinearSystem();

/**
 * @brief 后处理步骤：输出结果并计算误差
 */
void postprocess();

#endif // FEM_SOLVER_H
