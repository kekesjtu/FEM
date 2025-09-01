#ifndef PROBLEM_DEFINITION_H
#define PROBLEM_DEFINITION_H

#include <vector>

// 结构体：定义通用边界条件
// 表示形式为: K * (c * du/dx) + L * u = q
struct BoundaryCondition {
    int node_index;     // 作用的节点索引
    int K_bc;        // 导数项的系数
    double L_bc;        // 值项的系数
    double q_bc;        // 右侧项

    // 辅助构造函数
    // 构造狄利克雷边界条件: u = val
    // 对应 K=0, L=1, q=val
    static BoundaryCondition Dirichlet(int index, double val) {
        return {index, 0, 1.0, val};
    }
    // 构造罗宾边界条件: c*du/dx + h*u = g
    // 对应 K=1, L=h, q=g
    static BoundaryCondition Robin(int index, double h, double g) {
        return {index, 1, h, g};
    }
    // 构造诺曼边界条件: c*du/dx = q_flux
    // 对应 K=1, L=0, q=q_flux
    static BoundaryCondition Neumann(int index, double q_flux) {
        return {index, 1, 0.0, q_flux};
    }
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
