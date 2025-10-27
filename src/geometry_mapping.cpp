#include "geometry_mapping.h"
#include <cmath>  // For std::abs and std::sqrt

// =========================================================================
// ==              LineLinearMapping Implementation                       ==
// =========================================================================

/**
 * @brief 构造函数
 * 线段嵌入在 embed_dim 维空间中（1D/2D/3D）
 */
LineLinearMapping::LineLinearMapping(const std::vector<double>& element_coords, int embed_dim)
    : GeometryMapping1D(element_coords, 1, 2, embed_dim)  // order=1 (线性), nodes=2, embedding_dim
{
    // 坐标验证已在基类中完成
    computeJacobian();
}

/**
 * @brief 计算线段的雅可比（边长的一半）
 * 对于线性线段，雅可比是常数
 */
void LineLinearMapping::computeJacobian()
{
    // 线段的参考坐标是 t ∈ [-1, 1]
    // 物理坐标映射：x(t) = (1-t)/2 * x0 + (1+t)/2 * x1
    //
    // dx/dt = (-x0 + x1) / 2 = (x1 - x0) / 2
    //
    // 雅可比行列式就是 |dx/dt| = sqrt(sum((x1_i - x0_i)^2)) / 2
    //
    // 坐标存储格式：[x0, y0, x1, y1] 或 [x0, y0, z0, x1, y1, z1]

    double length_squared = 0.0;
    for (int d = 0; d < embedding_dimension_; ++d)
    {
        double coord0 = element_coords_[d];
        double coord1 = element_coords_[embedding_dimension_ + d];
        double dx = coord1 - coord0;
        length_squared += dx * dx;
    }

    double length = std::sqrt(length_squared);
    jacobian_det_ = length / 2.0;  // 雅可比 = 边长 / 2

    if (jacobian_det_ < 1e-15)
    {
        throw std::runtime_error("Degenerate line element: length is zero.");
    }

    // 对于1D单元，"逆雅可比"就是 1/jacobian_det_
    jacobian_inv_.resize(1, std::vector<double>(1));
    jacobian_inv_[0][0] = 1.0 / jacobian_det_;
}

void LineLinearMapping::mapToPhysical(const std::vector<double>& coord_ref,
                                      std::vector<double>& coord_phys) const
{
    if (coord_ref.size() != 1)
    {
        throw std::invalid_argument("Reference coordinate must have 1 component (t).");
    }

    double t = coord_ref[0];
    coord_phys.resize(embedding_dimension_);

    // 线性形函数：N0 = (1-t)/2, N1 = (1+t)/2
    double N0 = (1.0 - t) / 2.0;
    double N1 = (1.0 + t) / 2.0;

    for (int d = 0; d < embedding_dimension_; ++d)
    {
        coord_phys[d] = N0 * element_coords_[d] + N1 * element_coords_[embedding_dimension_ + d];
    }
}

double LineLinearMapping::getJacobianDet(const std::vector<double>& coord_ref) const
{
    // 对于线性线段，雅可比是常数
    (void)coord_ref;  // 避免未使用参数警告
    return jacobian_det_;
}

void LineLinearMapping::transformGradient(const std::vector<double>& gradient_ref,
                                          std::vector<double>& gradient_phys,
                                          const std::vector<double>& coord_ref) const
{
    // 对于1D单元嵌入在高维空间中的情况，
    // 参考坐标系下的梯度 dN/dt 需要转换到物理空间
    // gradient_phys = (dN/dt) * (dt/ds)，其中 s 是弧长参数
    // dt/ds = 1/jacobian_det_

    (void)coord_ref;  // 线性单元，雅可比是常数

    if (gradient_ref.size() != 1)
    {
        throw std::invalid_argument("Gradient in reference coordinates must have 1 component.");
    }

    // 对于1D边界单元，梯度仍然是1维的（沿着边界的导数）
    gradient_phys.resize(1);
    gradient_phys[0] = gradient_ref[0] * jacobian_inv_[0][0];
}

double LineLinearMapping::getElementLength() const
{
    // 边长 = 2 * jacobian_det_
    return 2.0 * jacobian_det_;
}

// =========================================================================
// ==              TriangleLinearMapping Implementation                   ==
// =========================================================================

/**
 * @brief 计算该单元的常数雅可比矩阵、其行列式和逆矩阵
 * 对于线性三角形，这些量在整个单元内都是不变的。
 */
void TriangleLinearMapping::computeJacobian()
{
    // 从基类的 element_coords 中提取节点坐标以便于计算
    // [x0, y0, x1, y1, x2, y2]
    const double x0 = element_coords_[0], y0 = element_coords_[1];
    const double x1 = element_coords_[2], y1 = element_coords_[3];
    const double x2 = element_coords_[4], y2 = element_coords_[5];

    // 雅可比矩阵 J = [[dx/dxi,  dx/deta],
    //                [dy/dxi,  dy/deta]]
    //
    // 对于线性三角形，形函数为:
    // N0 = 1 - xi - eta
    // N1 = xi
    // N2 = eta
    //
    // 坐标映射为:
    // x(xi, eta) = N0*x0 + N1*x1 + N2*x2
    // y(xi, eta) = N0*y0 + N1*y1 + N2*y2
    //
    // 求导可得:
    // dx/dxi  = x1 - x0
    // dx/deta = x2 - x0
    // dy/dxi  = y1 - y0
    // dy/deta = y2 - y0
    const double dx_dxi = x1 - x0;
    const double dx_deta = x2 - x0;
    const double dy_dxi = y1 - y0;
    const double dy_deta = y2 - y0;

    // 计算雅可比行列式
    jacobian_det_ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    if (std::abs(jacobian_det_) < 1e-15)
    {
        throw std::runtime_error("Degenerate triangle element: jacobian determinant is zero.");
    }

    // 计算逆雅可比矩阵 J^-1 = [[dxi/dx,  dxi/dy],
    //                          [deta/dx, deta/dy]]
    const double inv_det = 1.0 / jacobian_det_;
    jacobian_inv_.resize(2, std::vector<double>(2));
    jacobian_inv_[0][0] = dy_deta * inv_det;   // dxi/dx
    jacobian_inv_[0][1] = -dx_deta * inv_det;  // dxi/dy
    jacobian_inv_[1][0] = -dy_dxi * inv_det;   // deta/dx
    jacobian_inv_[1][1] = dx_dxi * inv_det;    // deta/dy

    // 这边改成每个jacobian指向一个函数，以适应高次单元
}

void TriangleLinearMapping::mapToPhysical(const std::vector<double>& coord_ref,
                                          std::vector<double>& coord_phys) const
{
    if (coord_ref.size() != 2)
    {
        throw std::invalid_argument("Reference coordinate must have 2 components (xi, eta).");
    }
    coord_phys.resize(2);

    const double xi = coord_ref[0];
    const double eta = coord_ref[1];

    // 使用形函数进行线性插值
    const double N0 = 1.0 - xi - eta;
    const double N1 = xi;
    const double N2 = eta;

    // x = N0*x0 + N1*x1 + N2*x2
    coord_phys[0] = N0 * element_coords_[0] + N1 * element_coords_[2] + N2 * element_coords_[4];
    // y = N0*y0 + N1*y1 + N2*y2
    coord_phys[1] = N0 * element_coords_[1] + N1 * element_coords_[3] + N2 * element_coords_[5];
}

double TriangleLinearMapping::getJacobianDet(const std::vector<double>& coord_ref) const
{
    // 对于线性三角形，雅可比行列式是常数，因此忽略输入参数 coord_ref
    (void)coord_ref;  // 避免 "unused parameter" 警告
    return jacobian_det_;
}

void TriangleLinearMapping::transformGradient(const std::vector<double>& gradient_ref,
                                              std::vector<double>& gradient_phys,
                                              const std::vector<double>& coord_ref) const
{
    // 对于线性三角形，逆雅可比矩阵是常数，因此忽略输入参数 coord_ref
    (void)coord_ref;  // 避免 "unused parameter" 警告

    if (gradient_ref.size() != 2)
    {
        throw std::invalid_argument("Reference gradient must have 2 components (dN/dxi, dN/deta).");
    }
    gradient_phys.resize(2);

    const double dN_dxi = gradient_ref[0];
    const double dN_deta = gradient_ref[1];

    // [dN/dx] = [dxi/dx  deta/dx] [dN/dxi ]
    // [dN/dy]   [dxi/dy  deta/dy] [dN/deta]
    // 注意：这里用的是 J_inv 的转置 (J_inv^T)
    // dN/dx = (dxi/dx) * dN/dxi + (deta/dx) * dN/deta
    gradient_phys[0] = jacobian_inv_[0][0] * dN_dxi + jacobian_inv_[1][0] * dN_deta;
    // dN/dy = (dxi/dy) * dN/dxi + (deta/dy) * dN/deta
    gradient_phys[1] = jacobian_inv_[0][1] * dN_dxi + jacobian_inv_[1][1] * dN_deta;
}

double TriangleLinearMapping::getElementArea() const
{
    // 三角形面积 = 0.5 * |det(J)|
    return 0.5 * std::abs(jacobian_det_);
}

std::unique_ptr<GeometryMapping> GeometryMappingFactory::createMapping(
    const std::vector<double>& element_coords, std::shared_ptr<Config> config_, bool is_boundary)
{
    // 根据 is_boundary 选择体单元类型或边界单元类型
    int element_type = is_boundary ? config_->getBoundaryElementType() : config_->getElementType();
    int order = config_->getOrder();
    int embed_dim = config_->getDimension();

    switch (element_type)
    {
        case Config::ElementType::LINE:
            switch (order)
            {
                case 1:
                    return std::make_unique<LineLinearMapping>(element_coords, embed_dim);
                default:
                    throw std::invalid_argument("Unsupported order for Line element");
            }

        case Config::ElementType::TRIANGLE:
            switch (order)
            {
                case 1:
                    return std::make_unique<TriangleLinearMapping>(element_coords);
                default:
                    throw std::invalid_argument("Unsupported order for Triangle");
            }

        case Config::ElementType::QUADRILATERAL:
            throw std::invalid_argument(
                "Requested ElementType is not implemented in GeometryMappingFactory");
        default:
            throw std::invalid_argument("Unknown ElementType");
    }
}