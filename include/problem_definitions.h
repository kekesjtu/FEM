#ifndef PROBLEM_DEFINITIONS_H
#define PROBLEM_DEFINITIONS_H

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "boundary_condition.h"
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
class RobinTestProblem : public ProblemSetup
{
  public:
    RobinTestProblem() : ProblemSetup("Robin Test")
    {
        // 设置物理问题：控制方程 -∇²u = -6
        setCoefficient([](const std::vector<double>&) -> double { return 1.0; });
        setSource([](const std::vector<double>&) -> double { return -6.0; });

        // 设置精确解（用于误差分析）
        setExactSolutionU(
            [](const std::vector<double>& coords) -> double
            {
                double x = coords[0];
                double y = coords[1];
                return x * x + 2.0 * y * y;  // u = x² + 2y²
            });

        setExactSolutionGradients(
            [](const std::vector<double>& coords, std::vector<double>& gradients) -> double
            {
                double x = coords[0];
                double y = coords[1];
                gradients[0] = 2.0 * x;      // ∂u/∂x = 2x
                gradients[1] = 4.0 * y;      // ∂u/∂y = 4y
                return x * x + 2.0 * y * y;  // u = x² + 2y²
            });
    }

    void setupBoundaryConditions(std::shared_ptr<Config> config) override
    {
        const auto& boundaries = config->getBoundary();
        const auto& coords = config->getNodeCoordinates();
        int dim = config->getDimension();

        boundary_conditions_.resize(boundaries.size());

        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            const auto& boundary = boundaries[i];

            // 计算边界中点坐标
            double x = 0.0;
            for (int node_idx : boundary.global_node_indices_in_element)
            {
                x += coords[node_idx * dim + 0];
            }
            x /= boundary.global_node_indices_in_element.size();

            // Robin 边界条件：∂u/∂n + u = 6 - 3x²
            double g = 6.0 - 3.0 * x * x;
            boundary_conditions_[i] = BoundaryCondition::Robin(1.0, g);
        }
    }
};

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
class ElectricFieldProblem : public ProblemSetup
{
  public:
    ElectricFieldProblem() : ProblemSetup("Electric Field")
    {
        // 设置物理问题：-∇·(σ∇V) = 0
        // 系数将在电热耦合求解器中动态处理（温度依赖）
        // 这里设置为占位符，实际不会被直接使用
        setCoefficient([](const std::vector<double>&) -> double { return 1.0; });
        setSource([](const std::vector<double>&) -> double { return 0.0; });

        // 电场问题通常没有解析解，不设置精确解
    }

    void setupBoundaryConditions(std::shared_ptr<Config> config) override
    {
        const auto& boundaries = config->getBoundary();
        const auto& coords = config->getNodeCoordinates();
        const auto& boundary_entities = config->getBoundaryGeometricEntities();
        int dim = config->getDimension();

        boundary_conditions_.resize(boundaries.size());

        // 使用几何实体来施加边界条件
        // 根据网格分析:
        // Entity 0, 1: X范围[-1.000, 0.000], X中心范围[-0.997, -0.052] → 左半圆 → V = 10.0V
        // Entity 2, 3: X范围[0.000, 1.000], X中心范围[0.052, 0.997] → 右半圆 → V = 0.0V
        //
        // 注意: 所有实体都包含x=0轴上的节点(因为边界单元连接两个节点)
        //       但边界单元的中心坐标明确区分了左右半圆
        //       Entity 0,1 的边界单元中心都在 x<0 区域
        //       Entity 2,3 的边界单元中心都在 x>0 区域

        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            const auto& boundary = boundaries[i];
            int entity_id = boundary_entities[i];

            // 使用几何实体直接判断 - 网格已经正确划分了左右半圆
            if (entity_id == 0 || entity_id == 1)
            {
                // Entity 0 和 Entity 1 都在左半圆(x<0) → V = 10.0V
                boundary_conditions_[i] = BoundaryCondition::Dirichlet(10.0);
            }
            else  // entity_id == 2 或 3
            {
                // Entity 2 和 Entity 3 都在右半圆(x>0) → V = 0.0V (接地)
                boundary_conditions_[i] = BoundaryCondition::Dirichlet(0.0);
            }
        }
    }
};

/**
 * @brief 电热耦合 - 热场设置
 *
 * 物理背景：
 * - 边界温度固定在参考温度 T0 - Dirichlet 边界条件
 * - 控制方程：-∇·(k∇T) = Q(x,y)，其中 Q 为焦耳热
 *
 * 注意：求解温度增量场 ΔT，边界条件 ΔT = 0 意味着边界保持在 T = T0
 */
class ThermalFieldProblem : public ProblemSetup
{
  public:
    ThermalFieldProblem() : ProblemSetup("Thermal Field")
    {
        // 设置物理问题：-∇·(k∇T) = Q
        // 热导率（将在电热耦合求解器中使用）
        // 源项将在电热耦合求解器中动态处理（焦耳热）
        setCoefficient([](const std::vector<double>&) -> double { return 1.0; });
        setSource([](const std::vector<double>&) -> double { return 0.0; });

        // 热场问题通常没有解析解，不设置精确解
    }

    void setupBoundaryConditions(std::shared_ptr<Config> config) override
    {
        const auto& boundaries = config->getBoundary();
        boundary_conditions_.resize(boundaries.size());

        // 所有边界温度增量设置为 ΔT = 0（对应绝对温度 T = T0）
        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            boundary_conditions_[i] = BoundaryCondition::Dirichlet(0.0);  // ΔT = 0
        }
    }
};

/**
 * @brief 3D泊松方程 - 正方体均匀电荷问题
 *
 * 物理背景：
 * - 控制方程：-∇·(ε∇V) = ρ，其中 ε = 1.0（介电常数），ρ = 10.0（电荷密度）
 * - 边界条件：六个面全部接地（V = 0）- Dirichlet边界条件
 * - 几何：单位正方体 [0,1]³
 *
 * 注意：这是一个经典的3D泊松方程算例，用于验证3D有限元求解器
 */
class CubeUniformChargeProblem : public ProblemSetup
{
  public:
    CubeUniformChargeProblem() : ProblemSetup("Cube Uniform Charge")
    {
        // 设置物理问题：-∇·(ε∇V) = ρ
        // 介电常数 ε = 1.0
        setCoefficient([](const std::vector<double>&) -> double { return 1.0; });

        // 均匀电荷密度 ρ = 10.0
        // 注意：泊松方程 -∇²V = ρ/ε，这里源项为 f = ρ/ε = 10.0/1.0 = 10.0
        setSource([](const std::vector<double>&) -> double { return 10.0; });

        // 该问题没有简单的解析解，不设置精确解
        // 可以通过网格细化收敛性测试来验证求解器正确性
    }

    void setupBoundaryConditions(std::shared_ptr<Config> config) override
    {
        const auto& boundaries = config->getBoundary();
        boundary_conditions_.resize(boundaries.size());

        // 所有边界面（正方体的六个面）都设置为接地：V = 0
        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            boundary_conditions_[i] = BoundaryCondition::Dirichlet(0.0);
        }
    }
};

/**
 * @brief 3D泊松方程 - 单位球域均匀源项问题（带解析解）
 *
 * 物理背景：
 * - 控制方程：-∇²u = 10 在单位球域 Ω = {(x,y,z) : x²+y²+z² < 1}
 * - 边界条件：u = 0 在球面边界 ∂Ω
 * - 几何：单位球域（半径 R = 1）
 * 该问题用于验证3D有限元求解器在球域上的精度
 */
class UnitSphereUniformSourceProblem : public ProblemSetup
{
  public:
    UnitSphereUniformSourceProblem() : ProblemSetup("Unit Sphere Uniform Source")
    {
        // 设置物理问题：-∇²u = 10
        setCoefficient([](const std::vector<double>&) -> double { return 1.0; });

        // 均匀源项 f = 10.0
        setSource([](const std::vector<double>&) -> double { return 10.0; });

        // 设置解析解：u(r) = (5/3)(1 - r²)
        setExactSolutionU(
            [](const std::vector<double>& coords) -> double
            {
                double x = coords[0];
                double y = coords[1];
                double z = coords[2];
                double r2 = x * x + y * y + z * z;  // r² = x² + y² + z²
                return (5.0 / 3.0) * (1.0 - r2);    // u = (5/3)(1 - r²)
            });

        // 设置解析解的梯度：∇u = -(10/3)(x, y, z)
        setExactSolutionGradients(
            [](const std::vector<double>& coords, std::vector<double>& gradients) -> double
            {
                double x = coords[0];
                double y = coords[1];
                double z = coords[2];
                double r2 = x * x + y * y + z * z;

                // ∇u = du/dr · ∇r = -2(5/3)r · (x,y,z)/r = -(10/3)(x,y,z)
                gradients[0] = -(10.0 / 3.0) * x;  // ∂u/∂x
                gradients[1] = -(10.0 / 3.0) * y;  // ∂u/∂y
                gradients[2] = -(10.0 / 3.0) * z;  // ∂u/∂z

                return 0.0;
            });
    }

    void setupBoundaryConditions(std::shared_ptr<Config> config) override
    {
        const auto& boundaries = config->getBoundary();
        boundary_conditions_.resize(boundaries.size());

        // 所有边界面（球面）都设置为零 Dirichlet 边界条件：u = 0
        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            boundary_conditions_[i] = BoundaryCondition::Dirichlet(0.0);
        }
    }
};

/**
 * @brief 自定义问题模板
 *
 * 用户可以复制此类并修改，创建自己的问题定义
 */
class CustomProblem : public ProblemSetup
{
  public:
    CustomProblem() : ProblemSetup("Custom Problem")
    {
        //  定义扩散系数和源项
        setCoefficient(
            [](const std::vector<double>&) -> double
            {
                // TODO: 用户自定义扩散系数 c(x,y)
                return 1.0;
            });

        setSource(
            [](const std::vector<double>&) -> double
            {
                // TODO: 用户自定义源项 f(x,y)
                return 0.0;
            });

        //  如果有解析解，用户可以在这里定义（可选，用于误差分析）
        setExactSolutionU(
            [](const std::vector<double>&) -> double
            {
                // TODO: 用户自定义精确解 u(x,y)
                return 0.0;
            });
    }

    void setupBoundaryConditions(std::shared_ptr<Config> config) override
    {
        // 1. 定义边界条件
        const auto& boundaries = config->getBoundary();
        boundary_conditions_.resize(boundaries.size());

        // 示例：对所有边界设置 Dirichlet 边界条件
        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            // TODO: 用户根据边界位置自定义边界条件
            boundary_conditions_[i] = BoundaryCondition::Dirichlet(0.0);
        }
    }
};

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
 * - 热导率 k：使用 ThermalFieldProblem 中的默认值（1.0）
 */
inline ElectrothermalProblem createElectrothermalProblem(
    const ElectrothermalParams& params = ElectrothermalParams())
{
    ElectrothermalProblem et_problem;
    et_problem.params = params;

    // 直接创建电场和热场问题
    et_problem.electric_problem = std::make_shared<ElectricFieldProblem>();
    et_problem.thermal_problem = std::make_shared<ThermalFieldProblem>();

    // 热导率使用 ThermalFieldProblem 中的默认值
    // 电场的电导率和热场的源项将在 ElectrothermalSolver 中动态处理

    return et_problem;
}

// ========================================
// 问题工厂 - 统一接口
// ========================================

/**
 * @brief 问题工厂类 - 通过字符串名称创建问题实例
 *
 * 使用工厂模式,提供统一接口来创建不同类型的问题
 *
 * 使用示例:
 * ```cpp
 * auto problem = ProblemFactory::createProblem("RobinTest");
 * auto problem = ProblemFactory::createProblem("ElectricField");
 * auto problem = ProblemFactory::createProblem("ThermalField");
 * auto problem = ProblemFactory::createProblem("Custom");
 * ```
 */
class ProblemFactory
{
  public:
    /**
     * @brief 创建问题实例
     * @param problem_name 问题名称（不区分大小写）
     * @return 问题实例指针
     * @throws std::invalid_argument 如果问题名称未知
     *
     * 支持的问题名称:
     * - "RobinTest" 或 "Robin" - Robin边界条件验证算例
     * - "ElectricField" 或 "Electric" - 电场问题
     * - "ThermalField" 或 "Thermal" - 热场问题
     * - "CubeUniformCharge" 或 "Cube" - 正方体均匀电荷问题
     * - "UnitSphereUniformSource" 或 "Sphere" - 单位球域均匀源项问题（带解析解）
     * - "Custom" - 自定义问题模板
     */
    static std::shared_ptr<ProblemSetup> createProblem(const std::string& problem_name)
    {
        // 转换为小写以实现不区分大小写
        std::string name_lower = toLower(problem_name);

        if (name_lower == "robintest" || name_lower == "robin")
        {
            return std::make_shared<RobinTestProblem>();
        }
        else if (name_lower == "electricfield" || name_lower == "electric")
        {
            return std::make_shared<ElectricFieldProblem>();
        }
        else if (name_lower == "thermalfield" || name_lower == "thermal")
        {
            return std::make_shared<ThermalFieldProblem>();
        }
        else if (name_lower == "custom")
        {
            return std::make_shared<CustomProblem>();
        }
        else if (name_lower == "cubeuniformcharge" || name_lower == "cube")
        {
            return std::make_shared<CubeUniformChargeProblem>();
        }
        else if (name_lower == "unitsphereuniformsource" || name_lower == "sphere")
        {
            return std::make_shared<UnitSphereUniformSourceProblem>();
        }
        else
        {
            throw std::invalid_argument("未知的问题类型: " + problem_name +
                                        "\n支持的类型: RobinTest, ElectricField, ThermalField, "
                                        "CubeUniformCharge, UnitSphereUniformSource, Custom");
        }
    }

    /**
     * @brief 创建电热耦合问题配置
     * @param params 物理参数
     * @return 电热耦合问题配置
     */
    static ElectrothermalProblem createElectrothermal(
        const ElectrothermalParams& params = ElectrothermalParams())
    {
        return createElectrothermalProblem(params);
    }

    /**
     * @brief 获取所有可用的问题类型列表
     * @return 问题类型名称列表
     */
    static std::vector<std::string> getAvailableProblems()
    {
        return {"RobinTest",         "ElectricField",           "ThermalField",
                "CubeUniformCharge", "UnitSphereUniformSource", "Custom"};
    }

  private:
    /**
     * @brief 将字符串转换为小写
     */
    static std::string toLower(const std::string& str)
    {
        std::string result = str;
        for (char& c : result)
        {
            c = std::tolower(static_cast<unsigned char>(c));
        }
        return result;
    }
};

}  // namespace ProblemLibrary

#endif  // PROBLEM_DEFINITIONS_H
