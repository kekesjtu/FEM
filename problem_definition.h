#ifndef PROBLEM_DEFINITION_H
#define PROBLEM_DEFINITION_H

#include <vector>

// 结构体：定义边界条件
struct BoundaryCondition {
    int node_index; // 节点索引
    double value;   // 边界条件的值
};

/**
 * @brief 定义方程中的物理系数 c(x)
 *        对应于方程 -d/dx( c(x) * du/dx ) = f(x)
 * @param x 坐标
 * @return double 系数c(x)的值
 */
double coefficient_c(double x);

/**
 * @brief 定义方程中的源项 f(x)
 *        对应于方程 -d/dx( c(x) * du/dx ) = f(x)
 * @param x 坐标
 * @return double 源项f(x)的值
 */
double source_term_f(double x);

/**
 * @brief 定义问题的精确解 u(x) (用于后处理中计算误差)
 * @param x 坐标
 * @return double 精确解u(x)的值
 */
double exact_solution_u(double x);

/**
 * @brief 定义问题的网格参数和边界条件
 * @param N_out 输出参数，全局总节点数
 * @param M_out 输出参数，单元总数
 * @param bcs_out 输出参数，边界条件向量
 */
void defineProblem(int& N_out, int& M_out, std::vector<BoundaryCondition>& bcs_out);

#endif // PROBLEM_DEFINITION_H
