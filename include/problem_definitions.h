#ifndef PROBLEM_DEFINITIONS_H
#define PROBLEM_DEFINITIONS_H

#include <memory>
#include "config.h"
#include "problem_setup.h"

/**
 * @brief 预定义问题库
 *
 * 包含常用的物理问题设置，用户可以直接调用这些函数创建问题定义
 * 也可以作为参考，自定义自己的问题设置
 */
namespace ProblemLibrary
{

/**
 * @brief 电热耦合问题的物理参数配置
 */
struct ElectrothermalParams
{
    double sigma0;  // 参考电导率 [S/m]
    double alpha;   // 电导率温度系数 [1/K]
    double T0;      // 参考温度 [K]

    ElectrothermalParams(double sigma0_ = 1.0, double alpha_ = 0.003, double T0_ = 300.0)
        : sigma0(sigma0_), alpha(alpha_), T0(T0_)
    {
    }
};

/**
 * @brief 电热耦合问题配置（包含两个场和物理参数）
 */
struct ElectrothermalProblem
{
    std::shared_ptr<ProblemSetup> electric_problem;
    std::shared_ptr<ProblemSetup> thermal_problem;
    ElectrothermalParams params;
};

/**
 * @brief Robin 边界条件验证算例
 *
 * 控制方程：-∇²u = -6 在单位圆域
 * 边界条件：∂u/∂n + u = 6 - 3x² 在单位圆边界
 * 解析解：u(x,y) = x² + 2y²
 */
inline std::shared_ptr<ProblemSetup> createRobinTest()
{
    auto problem = std::make_shared<ProblemSetup>("Robin Test");

    // 设置边界条件 - 直接操作 boundarys 数组
    problem->setBoundarySetupFunction(
        [](std::shared_ptr<Config> config)
        {
            auto& boundarys = config->getBoundaryMutable();
            const auto& coords = config->getNodeCoordinates();
            int dim = config->getDimension();

            for (auto& boundary : boundarys)
            {
                // 计算边界中点坐标
                double x = 0.0;
                for (int node_idx : boundary.global_node_indices_in_element)
                {
                    x += coords[node_idx * dim + 0];
                }
                x /= boundary.global_node_indices_in_element.size();

                // Robin 边界条件：∂u/∂n + u = 6 - 3x²
                double g = 6.0 - 3.0 * x * x;
                boundary.bc = Config::BoundaryCondition::Robin(1.0, g);
            }
        });

    // 设置物理问题：控制方程 -∇²u = -6
    problem->setCoefficient([](const std::vector<double>&) -> double { return 1.0; });
    problem->setSource([](const std::vector<double>&) -> double { return -6.0; });

    // 设置精确解（用于误差分析）
    problem->setExactSolutionU(
        [](const std::vector<double>& coords) -> double
        {
            double x = coords[0];
            double y = coords[1];
            return x * x + 2.0 * y * y;  // u = x² + 2y²
        });

    problem->setExactSolutionGradients(
        [](const std::vector<double>& coords, std::vector<double>& gradients) -> double
        {
            double x = coords[0];
            double y = coords[1];
            gradients[0] = 2.0 * x;      // ∂u/∂x = 2x
            gradients[1] = 4.0 * y;      // ∂u/∂y = 4y
            return x * x + 2.0 * y * y;  // u = x² + 2y²
        });

    return problem;
}

/**
 * @brief 电热耦合 - 电场设置
 *
 * 物理背景：
 * - 边界条件：根据位置设置不同电压，产生电流
 * - 控制方程：-∇·(σ(T)∇V) = 0
 *
 * 典型设置（实心圆盘）：
 * - 左半圆边界：高电压（V = 1.0 V）
 * - 右半圆边界：接地（V = 0.0 V）
 * - 这样会产生横向电流，从而产生焦耳热
 */
inline std::shared_ptr<ProblemSetup> createElectricField()
{
    auto problem = std::make_shared<ProblemSetup>("Electric Field");

    // 设置边界条件 - 直接操作 boundarys 数组
    problem->setBoundarySetupFunction(
        [](std::shared_ptr<Config> config)
        {
            auto& boundarys = config->getBoundaryMutable();
            const auto& coords = config->getNodeCoordinates();
            int dim = config->getDimension();

            const double epsilon = 0.5e-1;  // 容差：约0.1%的半径

            // 方案：为每条边界边的所有节点单独判断，避免冲突
            // 这样即使边界边跨越交界线，每个节点也能得到正确的值
            for (auto& boundary : boundarys)
            {
                // 检查这条边的所有节点的x坐标
                double x_min = 1e10, x_max = -1e10;
                for (int node_idx : boundary.global_node_indices_in_element)
                {
                    double x = coords[node_idx * dim + 0];
                    x_min = std::min(x_min, x);
                    x_max = std::max(x_max, x);
                }

                // 判断边界边的位置
                if (x_max < -epsilon)  // 完全在左侧
                {
                    boundary.bc = Config::BoundaryCondition::Dirichlet(10.0);  // 10V
                }
                else if (x_min > epsilon)  // 完全在右侧
                {
                    boundary.bc = Config::BoundaryCondition::Dirichlet(0.0);  // 接地
                }
                else  // 跨越交界或在交界处 (x_min <= ε 且 x_max >= -ε)
                {
                    // 对于跨越交界的边，使用Robin边界条件避免奇异性
                    // Robin: k·∂T/∂n + h·T = g, 这里简化为自然边界条件
                    // 或者使用线性插值的Dirichlet
                    double x_center = (x_min + x_max) / 2.0;
                    // 线性插值：V(x) = 5 - 5*x/ε 在[-ε, ε]区间
                    double V_interp = 5.0 * (1.0 - x_center / epsilon);
                    V_interp = std::max(0.0, std::min(10.0, V_interp));  // 限制在[0,10]
                    boundary.bc = Config::BoundaryCondition::Dirichlet(V_interp);
                }
            }
        });

    // 设置物理问题：-∇·(σ∇V) = 0
    // 系数将在电热耦合求解器中动态处理（温度依赖）
    // 这里设置为占位符，实际不会被直接使用
    problem->setCoefficient([](const std::vector<double>&) -> double { return 1.0; });
    problem->setSource([](const std::vector<double>&) -> double { return 0.0; });

    // 电场问题通常没有解析解，不设置精确解

    return problem;
}

/**
 * @brief 电热耦合 - 热场设置
 *
 * 物理背景：
 * - 边界温度固定在参考温度 T0 - Dirichlet 边界条件
 * - 控制方程：-∇·(k∇T) = Q(x,y)，其中 Q 为焦耳热
 *
 * 注意：求解温度增量场 ΔT，边界条件 ΔT = 0 意味着边界保持在 T = T0
 */
inline std::shared_ptr<ProblemSetup> createThermalField()
{
    auto problem = std::make_shared<ProblemSetup>("Thermal Field");

    // 设置边界条件 - 直接操作 boundarys 数组
    problem->setBoundarySetupFunction(
        [](std::shared_ptr<Config> config)
        {
            auto& boundarys = config->getBoundaryMutable();

            // 所有边界温度增量设置为 ΔT = 0（对应绝对温度 T = T0）
            for (auto& boundary : boundarys)
            {
                boundary.bc = Config::BoundaryCondition::Dirichlet(0.0);  // ΔT = 0
            }
        });

    // 设置物理问题：-∇·(k∇T) = Q
    // 热导率（将在电热耦合求解器中使用）
    // 源项将在电热耦合求解器中动态处理（焦耳热）
    problem->setCoefficient([](const std::vector<double>&) -> double { return 1.0; });
    problem->setSource([](const std::vector<double>&) -> double { return 0.0; });

    // 热场问题通常没有解析解，不设置精确解

    return problem;
}

/**
 * @brief 自定义问题模板
 *
 * 用户可以复制此函数并修改，创建自己的问题定义
 */
inline std::shared_ptr<ProblemSetup> createCustomProblem()
{
    auto problem = std::make_shared<ProblemSetup>("Custom Problem");

    // 🎯 1. 定义边界条件
    problem->setBoundarySetupFunction(
        [](std::shared_ptr<Config> config)
        {
            auto& boundarys = config->getBoundaryMutable();
            // const auto& coords = config->getNodeCoordinates();
            // int dim = config->getDimension();

            // 示例：对所有边界设置 Dirichlet 边界条件
            for (auto& boundary : boundarys)
            {
                // TODO: 用户根据边界位置自定义边界条件
                boundary.bc = Config::BoundaryCondition::Dirichlet(0.0);
            }
        });

    // 🎯 2. 定义扩散系数和源项
    problem->setCoefficient(
        [](const std::vector<double>&) -> double
        {
            // TODO: 用户自定义扩散系数 c(x,y)
            return 1.0;
        });

    problem->setSource(
        [](const std::vector<double>&) -> double
        {
            // TODO: 用户自定义源项 f(x,y)
            return 0.0;
        });

    // 🎯 3. 如果有解析解，用户可以在这里定义（可选，用于误差分析）
    problem->setExactSolutionU(
        [](const std::vector<double>&) -> double
        {
            // TODO: 用户自定义精确解 u(x,y)
            return 0.0;
        });

    return problem;
}

/**
 * @brief 创建电热耦合问题配置
 *
 * @param params 物理参数（sigma0, alpha, T0）
 * @return 电热耦合问题配置
 *
 * 物理背景：
 * - 电场方程：-∇·(σ(T)∇V) = 0，其中 σ(T) = σ0[1 + α(T-T0)]
 * - 热场方程：-∇·(k∇T) = Q(V,T)，其中 Q = σ(T)|∇V|²
 * - 边界条件：复用 createElectricField() 和 createThermalField()
 * - 热导率 k：使用 createThermalField() 中的默认值（1.0）
 */
inline ElectrothermalProblem createElectrothermalProblem(
    const ElectrothermalParams& params = ElectrothermalParams())
{
    ElectrothermalProblem et_problem;
    et_problem.params = params;

    // 🎯 直接复用已定义的电场和热场问题
    et_problem.electric_problem = createElectricField();
    et_problem.thermal_problem = createThermalField();

    // 热导率使用 createThermalField() 中的默认值
    // 电场的电导率和热场的源项将在 ElectrothermalSolver 中动态处理

    return et_problem;
}

}  // namespace ProblemLibrary

#endif  // PROBLEM_DEFINITIONS_H
