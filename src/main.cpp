#include <iostream>
#include "electrothermal_solver.h"
#include "fem_solver.h"
#include "problem_definitions.h"

/**
 * @brief 单场求解示例：Robin边界条件验证
 */
void solveSingleField()
{
    std::cout << "\n================================================" << std::endl;
    std::cout << "     单场求解" << std::endl;
    std::cout << "================================================\n" << std::endl;

    // 创建配置对象（只负责网格）
    auto config = std::make_shared<Config>();

    std::string problem_type =
        "CubeUniformCharge";  // 可选: "RobinTest", "ElectricField", "ThermalField",
                              // "CubeUniformCharge", "Custom"
    auto problem = ProblemLibrary::ProblemFactory::create(problem_type);

    std::cout << "求解问题: " << problem->getName() << std::endl;

    // 创建FEM求解器
    FEMSolver solver(config, problem);

    // 一键求解（边界条件、系数、源项已在 problem 中定义）
    solver.solveComplete();
}

/**
 * @brief 电热耦合求解示例
 *
 * 使用预定义的问题库设置电场和热场的边界条件
 */
void solveElectrothermal()
{
    std::cout << "\n================================================" << std::endl;
    std::cout << "        电热耦合仿真" << std::endl;
    std::cout << "================================================\n" << std::endl;

    // 创建网格配置（共享同一个网格）
    auto config = std::make_shared<Config>();

    // 🎯 定义物理参数并创建电热耦合问题配置
    ProblemLibrary::ElectrothermalParams params;
    params.sigma0 = 1.0;   // 参考电导率 [S/m]
    params.alpha = 0.003;  // 温度系数 [1/K]
    params.T0 = 300.0;     // 参考温度 [K]
    // 热导率使用 createThermalField() 中的默认值（1.0）

    // 🎯 使用工厂模式创建电热耦合问题
    auto et_problem = ProblemLibrary::ProblemFactory::createElectrothermal(params);

    std::cout << "电场边界条件: " << et_problem.electric_problem->getName() << std::endl;
    std::cout << "  - 左半圆边界 (x<0): V = 1.0 V (施加电压)" << std::endl;
    std::cout << "  - 右半圆边界 (x>0): V = 0.0 V (接地)" << std::endl;

    std::cout << "热场边界条件: " << et_problem.thermal_problem->getName() << std::endl;
    std::cout << "  - 所有边界: ΔT = 0 (保持参考温度 " << params.T0 << " K)" << std::endl;

    // 创建电热耦合求解器
    ElectrothermalSolver et_solver(config, et_problem,
                                   1e-6,  // 收敛容差
                                   100,   // 最大迭代次数
                                   true,  // 详细输出
                                   0.7);  // 松弛因子：0.5表示新旧解各占50%（可调整0.3-0.7）

    // 执行耦合求解
    et_solver.solve();

    // 输出结果
    et_solver.outputResults("results/electrothermal");

    std::cout << "\n电热耦合求解完成！" << std::endl;
    std::cout << "请使用ParaView查看以下文件：" << std::endl;
    std::cout << "  - results/electrothermal_V.vtu (电势场)" << std::endl;
    std::cout << "  - results/electrothermal_T.vtu (温度场)" << std::endl;
    std::cout << "  - results/electrothermal_Q.vtu (焦耳热密度)" << std::endl;
}

/**
 * @brief 主函数
 */
int main()
{
    // 选择求解模式
    int mode = 1;  // 1: 单场求解(3D正方体均匀电荷), 2: 电热耦合

    std::cout << "================================================" << std::endl;
    std::cout << "          有限元求解器" << std::endl;
    std::cout << "================================================" << std::endl;
    std::cout << "模式: " << (mode == 1 ? "单场求解(3D)" : "电热耦合") << std::endl;

    if (mode == 1)
    {
        solveSingleField();
    }
    else if (mode == 2)
    {
        solveElectrothermal();
    }

    std::cout << "\n================================================" << std::endl;
    std::cout << "          求解完成" << std::endl;
    std::cout << "================================================" << std::endl;

    return 0;
}
