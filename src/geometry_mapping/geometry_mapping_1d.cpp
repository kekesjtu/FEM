#include <cmath>  // For std::sqrt
#include "geometry_mapping.h"


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
