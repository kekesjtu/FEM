#include "shape_functions_2d.h"
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
double TriangleShapeFunction::computeTrialFunction2D(int node_index, double xi, double eta) const
{
    switch (node_index)
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

/**
 * @brief 线性三角形单元形函数对xi的偏导数
 */
double TriangleShapeFunction::computeTrialDerivativeXi(int node_index, double xi, double eta) const
{
    switch (node_index)
    {
        case 0:
            return -1.0;
        case 1:
            return 1.0;
        case 2:
            return 0.0;
        default:
            return 0.0;  // 错误情况
    }
}

/**
 * @brief 线性三角形单元形函数对eta的偏导数
 */
double TriangleShapeFunction::computeTrialDerivativeEta(int node_index, double xi, double eta) const
{
    switch (node_index)
    {
        case 0:
            return -1.0;
        case 1:
            return 0.0;
        case 2:
            return 1.0;
        default:
            return 0.0;  // 错误情况
    }
}

// ================ 工厂模式实现 ================

// ShapeFunctionFactory 静态成员初始化
std::unique_ptr<TriangleShapeFunction> ShapeFunctionFactory::triangleShapeFunction_ = nullptr;

std::unique_ptr<ShapeFunction> ShapeFunctionFactory::createShapeFunction(ElementType elementType,
                                                                         Order order)
{
    switch (elementType)
    {
        case ElementType::Triangle:
            return createTriangleShapeFunction(order);
        case ElementType::Quadrilateral:
            // 预留：四边形单元支持
            throw std::invalid_argument("Quadrilateral shape functions not implemented yet");
        case ElementType::Tetrahedron:
            // 预留：四面体单元支持
            throw std::invalid_argument("Tetrahedron shape functions not implemented yet");
        case ElementType::Hexahedron:
            // 预留：六面体单元支持
            throw std::invalid_argument("Hexahedron shape functions not implemented yet");
        default:
            throw std::invalid_argument("Unsupported element type");
    }
}

ShapeFunction* ShapeFunctionFactory::getShapeFunction(ElementType elementType, Order order)
{
    switch (elementType)
    {
        case ElementType::Triangle:
            return getTriangleShapeFunction(order);
        case ElementType::Quadrilateral:
            // 预留：四边形单元支持
            throw std::invalid_argument("Quadrilateral shape functions not implemented yet");
        case ElementType::Tetrahedron:
            // 预留：四面体单元支持
            throw std::invalid_argument("Tetrahedron shape functions not implemented yet");
        case ElementType::Hexahedron:
            // 预留：六面体单元支持
            throw std::invalid_argument("Hexahedron shape functions not implemented yet");
        default:
            throw std::invalid_argument("Unsupported element type");
    }
}

std::unique_ptr<ShapeFunction2D> ShapeFunctionFactory::createShapeFunction2D(
    ElementType elementType, Order order)
{
    switch (elementType)
    {
        case ElementType::Triangle:
            return createTriangleShapeFunction(order);
        case ElementType::Quadrilateral:
            // 预留：四边形单元支持
            throw std::invalid_argument("Quadrilateral 2D shape functions not implemented yet");
        default:
            throw std::invalid_argument("Unsupported 2D element type");
    }
}

std::unique_ptr<TriangleShapeFunction> ShapeFunctionFactory::createTriangleShapeFunction(
    Order order)
{
    switch (order)
    {
        case Order::Linear:
            return std::unique_ptr<TriangleShapeFunction>(new TriangleShapeFunction());
        case Order::Quadratic:
            // 预留：二次形函数支持
            throw std::invalid_argument("Quadratic triangle shape functions not implemented yet");
        case Order::Cubic:
            // 预留：三次形函数支持
            throw std::invalid_argument("Cubic triangle shape functions not implemented yet");
        default:
            throw std::invalid_argument("Unsupported shape function order");
    }
}

TriangleShapeFunction* ShapeFunctionFactory::getTriangleShapeFunction(Order order)
{
    if (order != Order::Linear)
    {
        throw std::invalid_argument("Only linear triangle shape functions are currently supported");
    }

    if (!triangleShapeFunction_)
    {
        triangleShapeFunction_ = createTriangleShapeFunction(order);
    }
    return triangleShapeFunction_.get();
}
