#include <cmath>  // For std::abs
#include "geometry_mapping.h"


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
