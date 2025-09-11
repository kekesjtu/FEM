#include "geometry_mapping_2d.h"
#include <cmath>
#include <stdexcept>
#include "fem_solver_2d.h"  // 为了访问全局变量P和T


GeometryMapping2D::GeometryMapping2D(int element_index)
{
    // 从全局T和P矩阵获取单元的三个顶点坐标
    if (element_index >= static_cast<int>(T.size()) || element_index < 0)
    {
        throw std::invalid_argument("Invalid element index");
    }

    x_coords.resize(3);
    y_coords.resize(3);

    for (int i = 0; i < 3; ++i)
    {
        int node_idx = T[element_index][i];
        x_coords[i] = P[node_idx * 2];      // x坐标
        y_coords[i] = P[node_idx * 2 + 1];  // y坐标
    }

    computeJacobian();
}

GeometryMapping2D::GeometryMapping2D(const std::vector<double>& coords)
{
    if (coords.size() != 6)
    {
        throw std::invalid_argument("coords must have 6 elements: [x0, y0, x1, y1, x2, y2]");
    }

    x_coords.resize(3);
    y_coords.resize(3);

    for (int i = 0; i < 3; ++i)
    {
        x_coords[i] = coords[i * 2];
        y_coords[i] = coords[i * 2 + 1];
    }

    computeJacobian();
}

void GeometryMapping2D::computeJacobian()
{
    // 雅可比矩阵 J = [dx/dx_ref   dx/dy_ref]
    //               [dy/dx_ref   dy/dy_ref]
    //
    // 对于线性三角形单元：
    // x(x_ref, y_ref) = x0*(1-x_ref-y_ref) + x1*x_ref + x2*y_ref
    // y(x_ref, y_ref) = y0*(1-x_ref-y_ref) + y1*x_ref + y2*y_ref
    //
    // 所以：
    // dx/dx_ref  = x1 - x0
    // dx/dy_ref = x2 - x0
    // dy/dx_ref  = y1 - y0
    // dy/dy_ref = y2 - y0

    double dx_dx_ref = x_coords[1] - x_coords[0];
    double dx_dy_ref = x_coords[2] - x_coords[0];
    double dy_dx_ref = y_coords[1] - y_coords[0];
    double dy_dy_ref = y_coords[2] - y_coords[0];

    // 计算雅可比行列式
    jacobian_det = dx_dx_ref * dy_dy_ref - dx_dy_ref * dy_dx_ref;

    if (std::abs(jacobian_det) < 1e-15)
    {
        throw std::runtime_error("Degenerate element: jacobian determinant is zero");
    }

    // 计算逆雅可比矩阵
    jacobian_inv.resize(2, std::vector<double>(2));
    jacobian_inv[0][0] = dy_dy_ref / jacobian_det;   // dx_ref/dx
    jacobian_inv[0][1] = -dx_dy_ref / jacobian_det;  // dx_ref/dy
    jacobian_inv[1][0] = -dy_dx_ref / jacobian_det;  // dy_ref/dx
    jacobian_inv[1][1] = dx_dx_ref / jacobian_det;   // dy_ref/dy
}

void GeometryMapping2D::mapToPhysical(double x_ref, double y_ref, double& x, double& y) const
{
    // 线性插值
    x = x_coords[0] * (1.0 - x_ref - y_ref) + x_coords[1] * x_ref + x_coords[2] * y_ref;
    y = y_coords[0] * (1.0 - x_ref - y_ref) + y_coords[1] * x_ref + y_coords[2] * y_ref;
}

void GeometryMapping2D::mapToReference(double x, double y, double& x_ref, double& y_ref) const
{
    // 求解线性方程组：
    // x = x0 + (x1-x0)*x_ref + (x2-x0)*y_ref
    // y = y0 + (y1-y0)*x_ref + (y2-y0)*y_ref
    //
    // 重写为：
    // x - x0 = (x1-x0)*x_ref + (x2-x0)*y_ref
    // y - y0 = (y1-y0)*x_ref + (y2-y0)*y_ref

    double dx = x - x_coords[0];
    double dy = y - y_coords[0];

    // 使用逆雅可比矩阵
    double dx_dx_ref = x_coords[1] - x_coords[0];
    double dx_dy_ref = x_coords[2] - x_coords[0];
    double dy_dx_ref = y_coords[1] - y_coords[0];
    double dy_dy_ref = y_coords[2] - y_coords[0];

    x_ref = (dy_dy_ref * dx - dx_dy_ref * dy) / jacobian_det;
    y_ref = (-dy_dx_ref * dx + dx_dx_ref * dy) / jacobian_det;
}

double GeometryMapping2D::getJacobianDet() const
{
    return jacobian_det;
}

void GeometryMapping2D::transformGradient(double dN_dx_ref, double dN_dy_ref, double& dN_dx,
                                          double& dN_dy) const
{
    // [dN/dx] = [dx_ref/dx   dy_ref/dx ] [dN/dx_ref ]
    // [dN/dy]   [dx_ref/dy   dy_ref/dy ] [dN/dy_ref]

    dN_dx = jacobian_inv[0][0] * dN_dx_ref + jacobian_inv[1][0] * dN_dy_ref;
    dN_dy = jacobian_inv[0][1] * dN_dx_ref + jacobian_inv[1][1] * dN_dy_ref;
}

double GeometryMapping2D::getElementArea() const
{
    return 0.5 * std::abs(jacobian_det);
}
