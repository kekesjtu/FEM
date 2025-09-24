#include <iostream>
#include "error_analysis_2d.h"
#include "fem_solver_2d.h"

/**
 * @brief 二维有限元求解器的主函数
 *
 * @return int 退出代码 (0 表示成功)
 */
int main()
{
    preprocess2D();
    std::cout << "步骤 1: 二维网格生成和预处理完成。" << std::endl;

    assemble2D();
    std::cout << "步骤 2: 二维全局矩阵组装完成。" << std::endl;

    applyBoundaryConditions2D();
    std::cout << "步骤 3: 二维边界条件施加完成。" << std::endl;

    solveLinearSystem2D();  // 默认参数为"CG", "DiagonalPreconditioner", 1e-8, 1000, true
    std::cout << "步骤 4: 二维线性方程组求解完成。" << std::endl;

    postprocess2D();
    std::cout << "步骤 5: 二维后处理完成。" << std::endl;

    // 清理全局配置
    FEMConfig::cleanup();

    return 0;
}
