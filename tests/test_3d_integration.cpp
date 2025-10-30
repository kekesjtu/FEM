/**
 * @file test_3d_integration.cpp
 * @brief 测试3D有限元的数值积分正确性
 * 
 * 验证:
 * 1. 雅可比行列式计算
 * 2. 高斯积分权重
 * 3. 形函数梯度变换
 * 4. 刚度矩阵单元积分
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../include/shape_functions.h"
#include "../include/geometry_mapping.h"
#include "../include/gauss_quadrature.h"
#include "../include/config.h"

void testUnitTetrahedronVolume()
{
    std::cout << "\n=== 测试1: 单位四面体体积计算 ===" << std::endl;
    
    // 标准参考四面体的顶点
    // 节点0: (0,0,0), 节点1: (1,0,0), 节点2: (0,1,0), 节点3: (0,0,1)
    std::vector<double> coords = {
        0.0, 0.0, 0.0,  // 节点0
        1.0, 0.0, 0.0,  // 节点1
        0.0, 1.0, 0.0,  // 节点2
        0.0, 0.0, 1.0   // 节点3
    };
    
    TetrahedronLinearMapping mapping(coords);
    double volume = mapping.getElementVolume();
    double expected_volume = 1.0 / 6.0;  // 标准四面体体积
    
    std::cout << "计算体积: " << volume << std::endl;
    std::cout << "理论体积: " << expected_volume << std::endl;
    std::cout << "相对误差: " << std::abs(volume - expected_volume) / expected_volume << std::endl;
    
    if (std::abs(volume - expected_volume) < 1e-10)
    {
        std::cout << "✓ 体积计算正确" << std::endl;
    }
    else
    {
        std::cout << "✗ 体积计算错误!" << std::endl;
    }
}

void testGaussIntegrationWeights()
{
    std::cout << "\n=== 测试2: 高斯积分权重和 ===" << std::endl;
    
    // 对于参考单元,所有权重之和应该等于参考单元的体积
    TetrahedronGauss1Point gauss1;
    TetrahedronGauss4Point gauss4;
    
    double sum1 = 0.0;
    for (int i = 0; i < gauss1.getNumPoints(); ++i)
    {
        sum1 += gauss1.getWeights()[i];
    }
    
    double sum4 = 0.0;
    for (int i = 0; i < gauss4.getNumPoints(); ++i)
    {
        sum4 += gauss4.getWeights()[i];
    }
    
    double expected_sum = 1.0 / 6.0;  // 参考四面体体积
    
    std::cout << "1点高斯积分权重和: " << sum1 << std::endl;
    std::cout << "4点高斯积分权重和: " << sum4 << std::endl;
    std::cout << "期望值: " << expected_sum << std::endl;
    
    if (std::abs(sum1 - expected_sum) < 1e-10 && std::abs(sum4 - expected_sum) < 1e-10)
    {
        std::cout << "✓ 高斯积分权重正确" << std::endl;
    }
    else
    {
        std::cout << "✗ 高斯积分权重错误!" << std::endl;
    }
}

void testShapeFunctionPartitionOfUnity()
{
    std::cout << "\n=== 测试3: 形函数单位分解性质 ===" << std::endl;
    
    TetrahedronLinearShapeFunction shape;
    
    // 在参考单元内任意一点测试
    std::vector<double> test_point = {0.2, 0.3, 0.1};
    
    double sum = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        sum += shape.computeTrialFunction(i, test_point);
    }
    
    std::cout << "形函数值之和: " << sum << std::endl;
    std::cout << "期望值: 1.0" << std::endl;
    
    if (std::abs(sum - 1.0) < 1e-10)
    {
        std::cout << "✓ 形函数满足单位分解性质" << std::endl;
    }
    else
    {
        std::cout << "✗ 形函数不满足单位分解性质!" << std::endl;
    }
}

void testGradientTransformation()
{
    std::cout << "\n=== 测试4: 梯度变换正确性 ===" << std::endl;
    
    // 使用单位四面体
    std::vector<double> coords = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    };
    
    TetrahedronLinearMapping mapping(coords);
    TetrahedronLinearShapeFunction shape;
    
    std::vector<double> ref_point = {0.25, 0.25, 0.25};
    
    // 检查梯度的散度(应该为零,因为线性形函数)
    double div_sum = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        std::vector<double> grad_ref = shape.computeTrialGradients(i, ref_point);
        std::vector<double> grad_phys;
        mapping.transformGradient(grad_ref, grad_phys, ref_point);
        
        std::cout << "节点" << i << " 物理梯度: (" 
                  << grad_phys[0] << ", " 
                  << grad_phys[1] << ", " 
                  << grad_phys[2] << ")" << std::endl;
        
        // 线性形函数的梯度和应该是零向量
        div_sum += grad_phys[0] + grad_phys[1] + grad_phys[2];
    }
    
    std::cout << "所有节点梯度之和的散度: " << div_sum << std::endl;
    
    if (std::abs(div_sum) < 1e-10)
    {
        std::cout << "✓ 梯度变换正确" << std::endl;
    }
    else
    {
        std::cout << "✗ 梯度变换可能有误!" << std::endl;
    }
}

void testIntegrationOfConstantFunction()
{
    std::cout << "\n=== 测试5: 积分常函数 f(x,y,z) = 1 ===" << std::endl;
    
    // 使用一个任意四面体
    std::vector<double> coords = {
        0.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        0.0, 3.0, 0.0,
        0.0, 0.0, 4.0
    };
    
    TetrahedronLinearMapping mapping(coords);
    TetrahedronGauss4Point gauss;
    
    double integral = 0.0;
    const auto& points = gauss.getPoints();
    const auto& weights = gauss.getWeights();
    
    for (int gp = 0; gp < gauss.getNumPoints(); ++gp)
    {
        std::vector<double> ref_coord(3);
        ref_coord[0] = points[gp * 3 + 0];
        ref_coord[1] = points[gp * 3 + 1];
        ref_coord[2] = points[gp * 3 + 2];
        
        double jac_det = mapping.getJacobianDet(ref_coord);
        integral += 1.0 * jac_det * weights[gp];
    }
    
    double expected_volume = mapping.getElementVolume();
    
    std::cout << "数值积分结果: " << integral << std::endl;
    std::cout << "四面体体积: " << expected_volume << std::endl;
    std::cout << "理论体积: " << (2.0 * 3.0 * 4.0) / 6.0 << std::endl;
    
    if (std::abs(integral - expected_volume) < 1e-10)
    {
        std::cout << "✓ 常函数积分正确" << std::endl;
    }
    else
    {
        std::cout << "✗ 常函数积分错误!" << std::endl;
    }
}

void testStiffnessMatrixSymmetry()
{
    std::cout << "\n=== 测试6: 刚度矩阵对称性 ===" << std::endl;
    
    // 使用单位四面体测试
    std::vector<double> coords = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    };
    
    TetrahedronLinearMapping mapping(coords);
    TetrahedronLinearShapeFunction shape;
    TetrahedronGauss4Point gauss;
    
    // 计算单元刚度矩阵 K_ij = ∫ ∇φ_i · ∇φ_j dV
    double K[4][4] = {{0}};
    
    const auto& points = gauss.getPoints();
    const auto& weights = gauss.getWeights();
    
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            double entry = 0.0;
            
            for (int gp = 0; gp < gauss.getNumPoints(); ++gp)
            {
                std::vector<double> ref_coord(3);
                ref_coord[0] = points[gp * 3 + 0];
                ref_coord[1] = points[gp * 3 + 1];
                ref_coord[2] = points[gp * 3 + 2];
                
                std::vector<double> grad_i_ref = shape.computeTrialGradients(i, ref_coord);
                std::vector<double> grad_j_ref = shape.computeTrialGradients(j, ref_coord);
                
                std::vector<double> grad_i_phys, grad_j_phys;
                mapping.transformGradient(grad_i_ref, grad_i_phys, ref_coord);
                mapping.transformGradient(grad_j_ref, grad_j_phys, ref_coord);
                
                double dot_product = 0.0;
                for (int d = 0; d < 3; ++d)
                {
                    dot_product += grad_i_phys[d] * grad_j_phys[d];
                }
                
                double jac_det = mapping.getJacobianDet(ref_coord);
                entry += dot_product * jac_det * weights[gp];
            }
            
            K[i][j] = entry;
        }
    }
    
    // 输出刚度矩阵
    std::cout << "刚度矩阵 K:" << std::endl;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            std::cout << std::setw(10) << std::setprecision(4) << K[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
    // 检查对称性
    double max_asymmetry = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = i + 1; j < 4; ++j)
        {
            double diff = std::abs(K[i][j] - K[j][i]);
            max_asymmetry = std::max(max_asymmetry, diff);
        }
    }
    
    std::cout << "最大非对称性: " << max_asymmetry << std::endl;
    
    if (max_asymmetry < 1e-10)
    {
        std::cout << "✓ 刚度矩阵对称" << std::endl;
    }
    else
    {
        std::cout << "✗ 刚度矩阵不对称!" << std::endl;
    }
    
    // 检查行和(应该接近零,因为常数场应该在零能模式中)
    std::cout << "\n行和(应接近零):" << std::endl;
    for (int i = 0; i < 4; ++i)
    {
        double row_sum = 0.0;
        for (int j = 0; j < 4; ++j)
        {
            row_sum += K[i][j];
        }
        std::cout << "行" << i << ": " << row_sum << std::endl;
    }
}

int main()
{
    std::cout << "======================================" << std::endl;
    std::cout << "  3D有限元数值积分正确性测试" << std::endl;
    std::cout << "======================================" << std::endl;
    
    testUnitTetrahedronVolume();
    testGaussIntegrationWeights();
    testShapeFunctionPartitionOfUnity();
    testGradientTransformation();
    testIntegrationOfConstantFunction();
    testStiffnessMatrixSymmetry();
    
    std::cout << "\n======================================" << std::endl;
    std::cout << "  测试完成" << std::endl;
    std::cout << "======================================" << std::endl;
    
    return 0;
}
