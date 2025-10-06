#include "shape_functions.h"
#include "config.h"
#include <stdexcept>

// ================ 面向对象实现 ================

/**
 * @brief 线性三角形单元形函数值计算
 *
 * 对于标准参考三角形，三个顶点为：
 * 节点0: (0, 0)
 * 节点1: (1, 0)
 * 节点2: (0, 1)
 *
 * 形函数为：
 * N0(xi, eta) = 1 - xi - eta
 * N1(xi, eta) = xi
 * N2(xi, eta) = eta
 */
double TriangleLinearShapeFunction::computeTrialFunction(int element_node_index,
                                                   const std::vector<double>& coords) const
{
    if (coords.size() != 2)
    {
        throw std::invalid_argument("Reference coordinates must have 2 components (xi, eta).");
    }
    double xi = coords[0];
    double eta = coords[1];
    switch (element_node_index)
    {
        case 0:
            return 1.0 - xi - eta;
        case 1:
            return xi;
        case 2:
            return eta;
        default:
            return 0.0;  // 错误情况
    }
}

double TriangleLinearShapeFunction::computeTestFunction(int element_node_index,
                                                  const std::vector<double>& coords) const
{
    // 对于伽辽金法，试探函数和检验函数相同
    return computeTrialFunction(element_node_index, coords);
}

/**
 * @brief 线性三角形单元形函数对xi的偏导数
 */
std::vector<double> TriangleLinearShapeFunction::computeTrialGradients(int element_node_index,
                                                   const std::vector<double>& coords) const
{
    switch (element_node_index)
    {
        case 0:
            return {-1.0, -1.0};
        case 1:
            return {1.0, 0.0};
        case 2:
            return {0.0, 1.0};
        default:
            return {0.0, 0.0};  // 错误情况
    }
}

/**
 * @brief 线性三角形单元形函数对eta的偏导数
 */
std::vector<double> TriangleLinearShapeFunction::computeTestGradients(int element_node_index,
                                                         const std::vector<double>& coords) const
{
    switch (element_node_index)
    {
        case 0:
            return {-1.0, -1.0};
        case 1:
            return {1.0, 0.0};
        case 2:
            return {0.0, 1.0};
        default:
            return {0.0, 0.0};  // 错误情况
    }
}

// ================ 工厂模式实现 ================

std::unique_ptr<ShapeFunction> ShapeFunctionFactory::createShapeFunction(
    std::shared_ptr<Config> config_)
{
    switch (config_->getElementType())
    {
        case Config::ElementType::TRIANGLE:
            switch (config_->getOrder())
            {
                case 1:
                    return std::make_unique<TriangleLinearShapeFunction>();
                default:
                    throw std::invalid_argument("Unsupported order for Triangle");
            }

        case Config::ElementType::QUADRILATERAL:
            throw std::invalid_argument(
                "Requested ElementType is not implemented in ShapeFunctionFactory");
        default:
            throw std::invalid_argument("Unknown ElementType");
    }
}