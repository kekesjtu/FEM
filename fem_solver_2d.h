#ifndef FEM_SOLVER_2D_H
#define FEM_SOLVER_2D_H

#include <vector>
#include "problem_definition.h"
#include <Eigen/Sparse>

// --- 全局变量与数据结构 ---
// 使用 extern 关键字声明全局变量，定义将在 fem_solver_2d.cpp 中
extern int N, M;                              // N: 节点数, M: 单元数
extern const int n;                           // 每个单元的节点数 (三角形单元为3)
extern Eigen::SparseMatrix<double> K_global; // 全局刚度矩阵
extern Eigen::VectorXd b;                     // 载荷向量
extern Eigen::VectorXd u;                     // 解向量
extern std::vector<std::vector<int>> T;       // 单元连接矩阵 [element][local_node] = global_node
extern std::vector<double> P;                 // 节点坐标 [x0, y0, x1, y1, ...]
extern std::vector<std::vector<int>> boundary_edges; // 边界边信息
extern std::vector<BoundaryCondition> boundary_conditions;

// 网格生成参数
struct Domain {
    double left, right, bottom, top;
};

struct BoundaryTypes {
    int bottom; // e.g., -1 for Dirichlet
    int right;  // e.g., -2 for Neumann
    int top;
    int left;
};

// --- 核心 FEM 函数声明 ---

/**
 * @brief 预处理步骤：生成网格、初始化矩阵和向量
 */
void preprocess2D();

/**
 * @brief 生成二维三角形网格
 * @param domain 计算域
 * @param N1 x方向单元数
 * @param N2 y方向单元数
 * @param bc_types 边界条件类型
 */
void generateMesh2D(const Domain& domain, int N1, int N2, const BoundaryTypes& bc_types);

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
void solveLinearSystem2D();

/**
 * @brief 后处理步骤：输出结果并计算误差
 */
void postprocess2D();

#endif // FEM_SOLVER_2D_H
