#include "problem_definition.h"
#include <cmath>

/**
 * @brief 定义方程中的物理系数 c(x)
 *        对应于方程 -d/dx( c(x) * du/dx ) = f(x)
 */
double coefficient_c(double x) {

    return exp(x);
}

/**
 * @brief 定义方程中的源项 f(x)
 *        对应于方程 -d/dx( c(x) * du/dx ) = f(x)
 */
double source_term_f(double x) {
    
    return -exp(x) * (cos(x) - 2 * sin(x) - x * cos(x) - x * sin(x));
}

/**
 * @brief 定义问题的精确解 u(x) (用于后处理中计算误差)
 */
double exact_solution_u(double x) {
  
    return x * cos(x);
}

/**
 * @brief 定义网格和边界条件
 */
void defineProblem(int& N_out, int& M_out, std::vector<BoundaryCondition>& bcs_out) {
    // 1. 网格参数
    N_out = 4; // 全局总节点数
    M_out = N_out - 1; // 单元总数

    // 2. 边界条件 (使用新的辅助构造函数)
    bcs_out.clear();
    // 节点 0 (x=0) 处, u=0 (狄利克雷)
    //bcs_out.push_back(BoundaryCondition::Dirichlet(0, 0.0)); 
    // 节点 N-1 (x=1) 处, u=cos(1) (狄利克雷)
    //bcs_out.push_back(BoundaryCondition::Dirichlet(N_out - 1, cos(1.0))); 

    // ---- 其他边界条件示例 (如果需要，可以取消注释) ----
    // 示例1: 罗宾边界条件
    bcs_out.push_back(BoundaryCondition::Robin(0, 1.0, 1.0));
    bcs_out.push_back(BoundaryCondition::Robin(N_out - 1, 1.0, 2.718*(cos(1.0)-sin(1.0))+cos(1.0)));
    // 示例2: 诺曼边界条件
    //bcs_out.push_back(BoundaryCondition::Neumann(N_out - 1, 2.718*(cos(1.0)-sin(1.0))));

}
