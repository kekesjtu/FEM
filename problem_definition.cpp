#include "problem_definition.h"
#include <cmath>

/**
 * @brief 定义方程中的物理系数 c(x)
 *        对应于方程 -d/dx( c(x) * du/dx ) = f(x)
 */
double coefficient_c(double x) {
    // 对于新问题，c(x) = e^x
    return exp(x);
}

/**
 * @brief 定义方程中的源项 f(x)
 *        对应于方程 -d/dx( c(x) * du/dx ) = f(x)
 */
double source_term_f(double x) {
    // 对于新问题，f(x) = -e^x * [cos(x) - 2sin(x) - x*cos(x) - x*sin(x)]
    return -exp(x) * (cos(x) - 2 * sin(x) - x * cos(x) - x * sin(x));
}

/**
 * @brief 定义问题的精确解 u(x) (用于后处理中计算误差)
 */
double exact_solution_u(double x) {
    // 对于新问题，u(x) = x * cos(x)
    return x * cos(x);
}

/**
 * @brief 定义网格和边界条件
 */
void defineProblem(int& N_out, int& M_out, std::vector<BoundaryCondition>& bcs_out) {
    // 1. 网格参数
    N_out = 4; // 全局总节点数
    M_out = N_out - 1; // 单元总数

    // 2. 边界条件
    bcs_out.clear();
    bcs_out.push_back({0, 0.0});         // 节点 0 (x=0) 处, u=0
    bcs_out.push_back({N_out - 1, cos(1.0)}); // 节点 N-1 (x=1) 处, u=cos(1)
}
