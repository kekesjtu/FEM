#include "gauss_quadrature_2d.h"
#include <stdexcept>

// ============================================================================
// GaussPoint 基类实现
// ============================================================================

double GaussPoint::getWeight(int index) const
{
    if (index < 0 || index >= numPoints_)
    {
        throw std::out_of_range("Index out of range for Gauss quadrature");
    }
    return weights_[index];
}

// ============================================================================
// TriangleGaussPoint 实现
// ============================================================================

void TriangleGaussPoint::getPoint(int index, double& x_ref, double& y_ref) const
{
    if (index < 0 || index >= numPoints_)
    {
        throw std::out_of_range("Index out of range for triangle Gauss quadrature");
    }
    x_ref = points_[index * 2];
    y_ref = points_[index * 2 + 1];
}

// ============================================================================
// 具体三角形高斯点类实现
// ============================================================================

// TriangleGauss1Point 实现
TriangleGauss1Point::TriangleGauss1Point() : TriangleGaussPoint(1)
{
    initializeGaussPoints();
}

void TriangleGauss1Point::initializeGaussPoints()
{
    // 一点积分：三角形中心
    points_.resize(2);
    weights_.resize(1);

    points_[0] = 1.0 / 3.0;  // x_ref
    points_[1] = 1.0 / 3.0;  // y_ref
    weights_[0] = 0.5;       // 权重（对于单位面积参考三角形）
}

// TriangleGauss3Point 实现
TriangleGauss3Point::TriangleGauss3Point() : TriangleGaussPoint(3)
{
    initializeGaussPoints();
}

void TriangleGauss3Point::initializeGaussPoints()
{
    // 三点积分：三个边的中点
    points_.resize(6);
    weights_.resize(3);

    double weight = 1.0 / 6.0;

    // 点1: 边(0,1)的中点
    points_[0] = 0.5;  // x_ref
    points_[1] = 0.0;  // y_ref
    weights_[0] = weight;

    // 点2: 边(1,2)的中点
    points_[2] = 0.5;  // x_ref
    points_[3] = 0.5;  // y_ref
    weights_[1] = weight;

    // 点3: 边(0,2)的中点
    points_[4] = 0.0;  // x_ref
    points_[5] = 0.5;  // y_ref
    weights_[2] = weight;
}

// TriangleGauss4Point 实现
TriangleGauss4Point::TriangleGauss4Point() : TriangleGaussPoint(4)
{
    initializeGaussPoints();
}

void TriangleGauss4Point::initializeGaussPoints()
{
    // 四点积分：三角形中心 + 三个顶点附近
    points_.resize(8);
    weights_.resize(4);

    // 点1: 三角形中心
    points_[0] = 1.0 / 3.0;  // x_ref
    points_[1] = 1.0 / 3.0;  // y_ref
    weights_[0] = -27.0 / 96.0;

    // 其他三个点（对称位置）
    double a = 0.6;
    double b = 0.2;
    double weight = 25.0 / 96.0;

    // 点2
    points_[2] = a;  // x_ref
    points_[3] = b;  // y_ref
    weights_[1] = weight;

    // 点3
    points_[4] = b;  // x_ref
    points_[5] = a;  // y_ref
    weights_[2] = weight;

    // 点4
    points_[6] = b;  // x_ref
    points_[7] = b;  // y_ref
    weights_[3] = weight;
}

// ============================================================================
// GaussPointFactory 实现
// ============================================================================

namespace
{
/**
 * @brief 创建三角形高斯积分对象的辅助函数
 * @tparam T 返回的智能指针类型
 * @param numPoints 积分点数量
 * @return 三角形高斯积分对象的智能指针
 */
template <typename T>
std::unique_ptr<T> createTriangleGaussPointHelper(int numPoints)
{
    switch (numPoints)
    {
        case 1:
            return std::unique_ptr<T>(new TriangleGauss1Point());
        case 3:
            return std::unique_ptr<T>(new TriangleGauss3Point());
        case 4:
            return std::unique_ptr<T>(new TriangleGauss4Point());
        default:
            throw std::invalid_argument("Unsupported number of Gauss points for triangle");
    }
}
}  // namespace

std::unique_ptr<GaussPoint> GaussPointFactory::createGaussPoint(ElementType elementType,
                                                                int numPoints)
{
    switch (elementType)
    {
        case ElementType::Triangle:
            return createTriangleGaussPointHelper<GaussPoint>(numPoints);
        case ElementType::Quadrilateral:
            throw std::invalid_argument("Quadrilateral elements not implemented yet");
        case ElementType::Tetrahedron:
            throw std::invalid_argument("Tetrahedron elements not implemented yet");
        case ElementType::Hexahedron:
            throw std::invalid_argument("Hexahedron elements not implemented yet");
        default:
            throw std::invalid_argument("Unknown element type");
    }
}

std::unique_ptr<GaussPoint2D> GaussPointFactory::createGaussPoint2D(ElementType elementType,
                                                                    int numPoints)
{
    switch (elementType)
    {
        case ElementType::Triangle:
            return createTriangleGaussPointHelper<GaussPoint2D>(numPoints);
        case ElementType::Quadrilateral:
            throw std::invalid_argument("Quadrilateral elements not implemented yet");
        default:
            throw std::invalid_argument("Unsupported 2D element type");
    }
}

std::unique_ptr<TriangleGaussPoint> GaussPointFactory::createTriangleGaussPoint(int numPoints)
{
    return createTriangleGaussPointHelper<TriangleGaussPoint>(numPoints);
}
