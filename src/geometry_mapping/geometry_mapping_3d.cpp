#include <cmath>  // For std::abs
#include "geometry_mapping.h"


// =========================================================================
// ==              TetrahedronLinearMapping Implementation                ==
// =========================================================================

/**
 * @brief 计算四面体单元的常数雅可比矩阵、其行列式和逆矩阵
 * 对于线性四面体，这些量在整个单元内都是不变的。
 */
void TetrahedronLinearMapping::computeJacobian()
{
    // 从基类的 element_coords 中提取节点坐标
    // [x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3]
    const double x0 = element_coords_[0], y0 = element_coords_[1], z0 = element_coords_[2];
    const double x1 = element_coords_[3], y1 = element_coords_[4], z1 = element_coords_[5];
    const double x2 = element_coords_[6], y2 = element_coords_[7], z2 = element_coords_[8];
    const double x3 = element_coords_[9], y3 = element_coords_[10], z3 = element_coords_[11];

    // 雅可比矩阵 J = [[dx/dxi,  dx/deta,  dx/dzeta],
    //                [dy/dxi,  dy/deta,  dy/dzeta],
    //                [dz/dxi,  dz/deta,  dz/dzeta]]
    //
    // 对于线性四面体，形函数为:
    // N0 = 1 - xi - eta - zeta
    // N1 = xi
    // N2 = eta
    // N3 = zeta
    //
    // 坐标映射为:
    // x(xi, eta, zeta) = N0*x0 + N1*x1 + N2*x2 + N3*x3
    // y(xi, eta, zeta) = N0*y0 + N1*y1 + N2*y2 + N3*y3
    // z(xi, eta, zeta) = N0*z0 + N1*z1 + N2*z2 + N3*z3
    //
    // 求导可得:
    // dx/dxi   = x1 - x0,  dx/deta  = x2 - x0,  dx/dzeta = x3 - x0
    // dy/dxi   = y1 - y0,  dy/deta  = y2 - y0,  dy/dzeta = y3 - y0
    // dz/dxi   = z1 - z0,  dz/deta  = z2 - z0,  dz/dzeta = z3 - z0

    const double dx_dxi = x1 - x0, dx_deta = x2 - x0, dx_dzeta = x3 - x0;
    const double dy_dxi = y1 - y0, dy_deta = y2 - y0, dy_dzeta = y3 - y0;
    const double dz_dxi = z1 - z0, dz_deta = z2 - z0, dz_dzeta = z3 - z0;

    // 计算雅可比行列式 (3×3矩阵的行列式)
    // det(J) = dx/dxi * (dy/deta * dz/dzeta - dy/dzeta * dz/deta)
    //        - dx/deta * (dy/dxi * dz/dzeta - dy/dzeta * dz/dxi)
    //        + dx/dzeta * (dy/dxi * dz/deta - dy/deta * dz/dxi)
    jacobian_det_ = dx_dxi * (dy_deta * dz_dzeta - dy_dzeta * dz_deta) -
                    dx_deta * (dy_dxi * dz_dzeta - dy_dzeta * dz_dxi) +
                    dx_dzeta * (dy_dxi * dz_deta - dy_deta * dz_dxi);

    if (std::abs(jacobian_det_) < 1e-15)
    {
        throw std::runtime_error("Degenerate tetrahedron element: jacobian determinant is zero.");
    }

    // 计算逆雅可比矩阵 J^-1 (使用代数余子式方法)
    const double inv_det = 1.0 / jacobian_det_;
    jacobian_inv_.resize(3, std::vector<double>(3));

    // J^-1 的第一行 [dxi/dx, dxi/dy, dxi/dz]
    jacobian_inv_[0][0] = (dy_deta * dz_dzeta - dy_dzeta * dz_deta) * inv_det;
    jacobian_inv_[0][1] = -(dx_deta * dz_dzeta - dx_dzeta * dz_deta) * inv_det;
    jacobian_inv_[0][2] = (dx_deta * dy_dzeta - dx_dzeta * dy_deta) * inv_det;

    // J^-1 的第二行 [deta/dx, deta/dy, deta/dz]
    jacobian_inv_[1][0] = -(dy_dxi * dz_dzeta - dy_dzeta * dz_dxi) * inv_det;
    jacobian_inv_[1][1] = (dx_dxi * dz_dzeta - dx_dzeta * dz_dxi) * inv_det;
    jacobian_inv_[1][2] = -(dx_dxi * dy_dzeta - dx_dzeta * dy_dxi) * inv_det;

    // J^-1 的第三行 [dzeta/dx, dzeta/dy, dzeta/dz]
    jacobian_inv_[2][0] = (dy_dxi * dz_deta - dy_deta * dz_dxi) * inv_det;
    jacobian_inv_[2][1] = -(dx_dxi * dz_deta - dx_deta * dz_dxi) * inv_det;
    jacobian_inv_[2][2] = (dx_dxi * dy_deta - dx_deta * dy_dxi) * inv_det;
}

void TetrahedronLinearMapping::mapToPhysical(const std::vector<double>& coord_ref,
                                             std::vector<double>& coord_phys) const
{
    if (coord_ref.size() != 3)
    {
        throw std::invalid_argument("Reference coordinate must have 3 components (xi, eta, zeta).");
    }
    coord_phys.resize(3);

    const double xi = coord_ref[0];
    const double eta = coord_ref[1];
    const double zeta = coord_ref[2];

    // 使用形函数进行线性插值
    const double N0 = 1.0 - xi - eta - zeta;
    const double N1 = xi;
    const double N2 = eta;
    const double N3 = zeta;

    // x = N0*x0 + N1*x1 + N2*x2 + N3*x3
    coord_phys[0] = N0 * element_coords_[0] + N1 * element_coords_[3] + N2 * element_coords_[6] +
                    N3 * element_coords_[9];
    // y = N0*y0 + N1*y1 + N2*y2 + N3*y3
    coord_phys[1] = N0 * element_coords_[1] + N1 * element_coords_[4] + N2 * element_coords_[7] +
                    N3 * element_coords_[10];
    // z = N0*z0 + N1*z1 + N2*z2 + N3*z3
    coord_phys[2] = N0 * element_coords_[2] + N1 * element_coords_[5] + N2 * element_coords_[8] +
                    N3 * element_coords_[11];
}

double TetrahedronLinearMapping::getJacobianDet(const std::vector<double>& coord_ref) const
{
    // 对于线性四面体，雅可比行列式是常数，因此忽略输入参数 coord_ref
    (void)coord_ref;  // 避免 "unused parameter" 警告
    return jacobian_det_;
}

void TetrahedronLinearMapping::transformGradient(const std::vector<double>& gradient_ref,
                                                 std::vector<double>& gradient_phys,
                                                 const std::vector<double>& coord_ref) const
{
    // 对于线性四面体，逆雅可比矩阵是常数，因此忽略输入参数 coord_ref
    (void)coord_ref;  // 避免 "unused parameter" 警告

    if (gradient_ref.size() != 3)
    {
        throw std::invalid_argument(
            "Reference gradient must have 3 components (dN/dxi, dN/deta, dN/dzeta).");
    }
    gradient_phys.resize(3);

    const double dN_dxi = gradient_ref[0];
    const double dN_deta = gradient_ref[1];
    const double dN_dzeta = gradient_ref[2];

    // [dN/dx] = [dxi/dx  deta/dx  dzeta/dx] [dN/dxi  ]
    // [dN/dy]   [dxi/dy  deta/dy  dzeta/dy] [dN/deta ]
    // [dN/dz]   [dxi/dz  deta/dz  dzeta/dz] [dN/dzeta]
    // 注意：这里用的是 J_inv 的转置 (J_inv^T)
    gradient_phys[0] = jacobian_inv_[0][0] * dN_dxi + jacobian_inv_[1][0] * dN_deta +
                       jacobian_inv_[2][0] * dN_dzeta;
    gradient_phys[1] = jacobian_inv_[0][1] * dN_dxi + jacobian_inv_[1][1] * dN_deta +
                       jacobian_inv_[2][1] * dN_dzeta;
    gradient_phys[2] = jacobian_inv_[0][2] * dN_dxi + jacobian_inv_[1][2] * dN_deta +
                       jacobian_inv_[2][2] * dN_dzeta;
}

double TetrahedronLinearMapping::getElementVolume() const
{
    // 四面体体积 = (1/6) * |det(J)|
    return (1.0 / 6.0) * std::abs(jacobian_det_);
}
