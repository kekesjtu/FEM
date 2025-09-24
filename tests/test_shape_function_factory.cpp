#include <iostream>
#include <memory>
#include "shape_functions_2d.h"

/**
 * @brief 测试形函数工厂的简化接口
 */
int main()
{
    try
    {
        std::cout << "=== 形函数工厂模式测试 ===" << std::endl;

        // ================ 方式1：使用统一的createShapeFunction接口 ================
        std::cout << "\n1. 使用createShapeFunction创建形函数对象：" << std::endl;

        // 创建三角形线性形函数（默认线性阶数）
        auto triangleShapeFunc =
            ShapeFunctionFactory::createShapeFunction(ShapeFunctionFactory::ElementType::Triangle);

        // 也可以显式指定阶数
        auto triangleLinearFunc = ShapeFunctionFactory::createShapeFunction(
            ShapeFunctionFactory::ElementType::Triangle, ShapeFunctionFactory::Order::Linear);

        std::cout << "   成功创建三角形形函数，节点数: " << triangleShapeFunc->getNumNodes()
                  << std::endl;

        // ================ 方式2：使用单例模式的getShapeFunction接口 ================
        std::cout << "\n2. 使用getShapeFunction获取共享实例：" << std::endl;

        // 获取共享的三角形形函数实例（避免重复创建）
        auto* sharedTriangleFunc =
            ShapeFunctionFactory::getShapeFunction(ShapeFunctionFactory::ElementType::Triangle);

        std::cout << "   获取共享三角形形函数实例，节点数: " << sharedTriangleFunc->getNumNodes()
                  << std::endl;

        // ================ 测试形函数计算 ================
        std::cout << "\n3. 测试形函数计算：" << std::endl;

        double xi = 0.2, eta = 0.3;
        std::vector<double> coords = {xi, eta};

        for (int i = 0; i < 3; ++i)
        {
            double value = triangleShapeFunc->computeTrialFunction(i, coords);
            auto derivatives = triangleShapeFunc->computeTrialDerivatives(i, coords);

            std::cout << "   节点" << i << ": N = " << value << ", dN/dxi = " << derivatives[0]
                      << ", dN/deta = " << derivatives[1] << std::endl;
        }

        // ================ 验证新工厂接口与直接方法调用的一致性 ================
        std::cout << "\n4. 验证新工厂接口的一致性：" << std::endl;

        // 使用工厂接口创建的对象
        auto directTriangleFunc =
            ShapeFunctionFactory::createShapeFunction(ShapeFunctionFactory::ElementType::Triangle);
        auto* directTrianglePtr = dynamic_cast<TriangleShapeFunction*>(directTriangleFunc.get());

        for (int i = 0; i < 3; ++i)
        {
            double direct_value = directTrianglePtr->computeTrialFunction2D(i, xi, eta);
            double direct_dxi = directTrianglePtr->computeTrialDerivativeXi(i, xi, eta);
            double direct_deta = directTrianglePtr->computeTrialDerivativeEta(i, xi, eta);

            std::cout << "   节点" << i << " (直接调用): N = " << direct_value
                      << ", dN/dxi = " << direct_dxi << ", dN/deta = " << direct_deta << std::endl;
        }

        // ================ 测试不支持的类型 ================
        std::cout << "\n5. 测试异常处理：" << std::endl;
        try
        {
            auto quadFunc = ShapeFunctionFactory::createShapeFunction(
                ShapeFunctionFactory::ElementType::Quadrilateral);
        }
        catch (const std::invalid_argument& e)
        {
            std::cout << "   预期异常: " << e.what() << std::endl;
        }

        std::cout << "\n=== 所有测试完成! ===" << std::endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "测试失败: " << e.what() << std::endl;
        return 1;
    }
}