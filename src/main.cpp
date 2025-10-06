#include <iostream>
#include "error_analysis.h"
#include "fem_solver.h"

/**
 * @brief 有限元求解器的主函数
 *
 * @return int 退出代码 (0 表示成功)
 */
int main()
{
    preprocess();
    std::cout << "步骤 1: 网格生成和预处理完成。" << std::endl;

    assemble();
    std::cout << "步骤 2: 全局矩阵组装完成。" << std::endl;

    applyBoundaryConditions();
    std::cout << "步骤 3: 边界条件施加完成。" << std::endl;

    solveLinearSystem();  // 默认参数为"CG", "DiagonalPreconditioner", 1e-8, 1000, true
    std::cout << "步骤 4: 线性方程组求解完成。" << std::endl;

    postprocess();
    std::cout << "步骤 5: 后处理完成。" << std::endl;

    return 0;
}
