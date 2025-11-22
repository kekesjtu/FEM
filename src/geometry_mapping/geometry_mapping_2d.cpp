#include <cmath>  // For std::abs, std::sqrt
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
    // 从基类的 element_coords 中提取节点坐标
    // 2D: [x0, y0, x1, y1, x2, y2] (6个元素)
    // 3D边界: [x0, y0, z0, x1, y1, z1, x2, y2, z2] (9个元素)

    const int dim = embedding_dimension_;  // 嵌入空间维度

    // 提取三个顶点的坐标
    const double x0 = element_coords_[0 * dim + 0], y0 = element_coords_[0 * dim + 1];
    const double x1 = element_coords_[1 * dim + 0], y1 = element_coords_[1 * dim + 1];
    const double x2 = element_coords_[2 * dim + 0], y2 = element_coords_[2 * dim + 1];

    // 雅可比矩阵的计算
    // 对于2D三角形（dim=2）：J 是 2×2 矩阵
    // 对于3D曲面三角形（dim=3）：J 是 3×2 矩阵，需要用不同方法计算行列式

    const double dx_dxi = x1 - x0;
    const double dx_deta = x2 - x0;
    const double dy_dxi = y1 - y0;
    const double dy_deta = y2 - y0;

    if (dim == 2)
    {
        // 2D情况：标准的2×2雅可比行列式
        jacobian_det_ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    }
    else if (dim == 3)
    {
        // 3D曲面情况：需要考虑z坐标
        const double z0 = element_coords_[0 * dim + 2];
        const double z1 = element_coords_[1 * dim + 2];
        const double z2 = element_coords_[2 * dim + 2];

        const double dz_dxi = z1 - z0;
        const double dz_deta = z2 - z0;

        // 雅可比矩阵 J = [dx/dxi  dx/deta ]
        //                [dy/dxi  dy/deta ]
        //                [dz/dxi  dz/deta ]
        //
        // 对于3×2矩阵，"行列式"实际是曲面法向量的模长
        // |J| = ||∂r/∂ξ × ∂r/∂η||
        // 其中 ∂r/∂ξ = (dx/dxi, dy/dxi, dz/dxi)
        //      ∂r/∂η = (dx/deta, dy/deta, dz/deta)

        // 叉积: n = ∂r/∂ξ × ∂r/∂η
        double nx = dy_dxi * dz_deta - dz_dxi * dy_deta;
        double ny = dz_dxi * dx_deta - dx_dxi * dz_deta;
        double nz = dx_dxi * dy_deta - dy_dxi * dx_deta;

        // 雅可比行列式 = ||n||
        jacobian_det_ = std::sqrt(nx * nx + ny * ny + nz * nz);
    }
    else
    {
        throw std::runtime_error("TriangleLinearMapping only supports 2D or 3D embedding");
    }

    if (std::abs(jacobian_det_) < 1e-15)
    {
        throw std::runtime_error("Degenerate triangle element: jacobian determinant is zero.");
    }

    // 计算逆雅可比矩阵（仅对2D情况）
    if (dim == 2)
    {
        const double inv_det = 1.0 / jacobian_det_;
        jacobian_inv_.resize(2, std::vector<double>(2));
        jacobian_inv_[0][0] = dy_deta * inv_det;   // dxi/dx
        jacobian_inv_[0][1] = -dx_deta * inv_det;  // dxi/dy
        jacobian_inv_[1][0] = -dy_dxi * inv_det;   // deta/dx
        jacobian_inv_[1][1] = dx_dxi * inv_det;    // deta/dy
    }
    // 对于3D曲面，暂不计算逆矩阵（当前边界积分不需要梯度）
}

void TriangleLinearMapping::mapToPhysical(const std::vector<double>& coord_ref,
                                          std::vector<double>& coord_phys) const
{
    if (coord_ref.size() != 2)
    {
        throw std::invalid_argument("Reference coordinate must have 2 components (xi, eta).");
    }
    coord_phys.resize(embedding_dimension_);

    const double xi = coord_ref[0];
    const double eta = coord_ref[1];

    // 使用形函数进行线性插值
    const double N0 = 1.0 - xi - eta;
    const double N1 = xi;
    const double N2 = eta;

    const int dim = embedding_dimension_;

    // 对每个坐标维度进行插值
    for (int d = 0; d < dim; ++d)
    {
        coord_phys[d] = N0 * element_coords_[0 * dim + d] + N1 * element_coords_[1 * dim + d] +
                        N2 * element_coords_[2 * dim + d];
    }
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
