#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <Eigen/Sparse>
#include <memory>
#include <vector>
#include "config.h"

// --- 全局变量与数据结构 ---
// 使用 extern 关键字声明全局变量，网格和问题定义通过Config对象管理
extern int N, M;                              // N: 节点数, M: 单元数
extern int n;                                 // 每个单元的节点数 (根据单元类型动态确定)
extern Eigen::SparseMatrix<double> K_global;  // 全局刚度矩阵
extern Eigen::VectorXd b;                     // 载荷向量
extern Eigen::VectorXd u;                     // 解向量
extern std::shared_ptr<Config> config;        // 配置对象，包含网格和问题定义

// --- 核心 FEM 函数声明 ---

/**
 * @brief 预处理步骤：生成网格、初始化矩阵和向量
 */
void preprocess();

/**
 * @brief 组装全局刚度矩阵和载荷向量
 */
void assemble();

/**
 * @brief 施加边界条件
 */
void applyBoundaryConditions();

/**
 * @brief 求解线性方程组 K_global * u = b
 */
void solveLinearSystem(const std::string& solver_type = "CG",
                       const std::string& preconditioner_type = "DiagonalPreconditioner",
                       double tol = 1e-8, int max_iter = 1000, bool verbose = true);

/**
 * @brief 后处理步骤：输出结果并计算误差
 */
void postprocess();

#endif  // FEM_SOLVER_H
