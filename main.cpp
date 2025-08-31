#include "fem_solver.h"
#include <iostream>

/**
 * @brief 主函数，程序的入口点
 * 
 * @return int 退出代码 (0 表示成功)
 */
int main() {
    preprocess();
    std::cout << "步骤 1: 预处理完成。" << std::endl;
    
    assemble();
    std::cout << "步骤 2: 全局矩阵组装完成。" << std::endl;
    
    applyBoundaryConditions();
    std::cout << "步骤 3: 边界条件施加完成。" << std::endl;
    
    solveLinearSystem();
    std::cout << "步骤 4: 线性方程组求解完成。" << std::endl;
    
    postprocess();
    std::cout << "步骤 5: 后处理完成。" << std::endl;
    
    return 0;
}
