/**
 * @file test_simple_3d_poisson.cpp
 * @brief 测试简单的3D泊松方程求解
 * 
 * 问题设置:
 * -∇²u = 1, 在单位立方体 [0,1]³
 * u = 0, 在所有边界
 * 
 * 这个问题有近似解析解(对于正方体中心):
 * u(0.5, 0.5, 0.5) ≈ 0.0234 (使用级数解的前几项估计)
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../include/fem_solver.h"
#include "../include/config.h"
#include "../include/problem_setup.h"

class SimpleCubeProblem : public ProblemSetup
{
public:
    SimpleCubeProblem() : ProblemSetup("Simple Cube f=1")
    {
        // -∇²u = 1
        setCoefficient([](const std::vector<double>&) -> double { return 1.0; });
        setSource([](const std::vector<double>&) -> double { return 1.0; });
    }
    
    void setupBoundaryConditions(std::shared_ptr<Config> config) override
    {
        const auto& boundaries = config->getBoundary();
        boundary_conditions_.resize(boundaries.size());
        
        // 所有边界 u = 0
        for (size_t i = 0; i < boundaries.size(); ++i)
        {
            boundary_conditions_[i] = BoundaryCondition::Dirichlet(0.0);
        }
    }
};

int main()
{
    std::cout << "======================================" << std::endl;
    std::cout << "  简单3D泊松方程测试" << std::endl;
    std::cout << "  -∇²u = 1" << std::endl;
    std::cout << "  u|∂Ω = 0" << std::endl;
    std::cout << "======================================\n" << std::endl;
    
    // 创建配置
    auto config = std::make_shared<Config>();
    std::cout << "加载网格: " << config->getMeshFilename() << std::endl;
    std::cout << "维度: " << config->getDimension() << "D" << std::endl;
    std::cout << "节点数: " << config->getNodesNum() << std::endl;
    std::cout << "单元数: " << config->getElementsNum() << std::endl;
    std::cout << "边界单元数: " << config->getBoundary().size() << std::endl;
    
    // 创建问题
    auto problem = std::make_shared<SimpleCubeProblem>();
    
    // 创建求解器
    FEMSolver solver(config, problem);
    
    // 求解
    solver.solveComplete();
    
    // 分析结果
    const auto& u = solver.getSolution();
    const auto& coords = config->getNodeCoordinates();
    int dim = config->getDimension();
    
    // 找到最大解和中心点的解
    double u_max = 0.0;
    double u_center = 0.0;
    int center_node = -1;
    
    for (int i = 0; i < config->getNodesNum(); ++i)
    {
        double x = coords[i * dim + 0];
        double y = coords[i * dim + 1];
        double z = coords[i * dim + 2];
        
        if (std::abs(u[i]) > std::abs(u_max))
        {
            u_max = u[i];
        }
        
        // 寻找最接近中心 (0.5, 0.5, 0.5) 的节点
        double dist = std::sqrt(std::pow(x - 0.5, 2) + 
                               std::pow(y - 0.5, 2) + 
                               std::pow(z - 0.5, 2));
        if (dist < 0.05 && center_node == -1)  // 在中心附近5%范围内
        {
            center_node = i;
            u_center = u[i];
            std::cout << "\n找到接近中心的节点 #" << i << std::endl;
            std::cout << "  坐标: (" << x << ", " << y << ", " << z << ")" << std::endl;
            std::cout << "  解值: " << u_center << std::endl;
        }
    }
    
    std::cout << "\n======================================" << std::endl;
    std::cout << "  结果分析" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "最大解值: " << u_max << std::endl;
    std::cout << "中心点解值: " << u_center << std::endl;
    std::cout << "\n粗略理论估计:" << std::endl;
    std::cout << "对于 -∇²u = 1, u|∂Ω = 0" << std::endl;
    std::cout << "中心点解约为: 0.02 - 0.03" << std::endl;
    std::cout << "(根据有限差分方法的粗略估计)" << std::endl;
    
    // 检查解的合理性
    if (u_max > 0.01 && u_max < 0.05)
    {
        std::cout << "\n✓ 解的数量级合理" << std::endl;
    }
    else
    {
        std::cout << "\n✗ 解的数量级可能有误!" << std::endl;
    }
    
    return 0;
}
