#include <stdexcept>
#include "shape_functions.h"


// ============================================================================
// 三维线性四面体单元形函数实现
// ============================================================================

/**
 * @brief 线性四面体单元形函数值计算
 *
 * 对于标准参考四面体，四个顶点为：
 * 节点0: (0, 0, 0)
 * 节点1: (1, 0, 0)
 * 节点2: (0, 1, 0)
 * 节点3: (0, 0, 1)
 *
 * 形函数为：
 * N0(xi, eta, zeta) = 1 - xi - eta - zeta
 * N1(xi, eta, zeta) = xi
 * N2(xi, eta, zeta) = eta
 * N3(xi, eta, zeta) = zeta
 */
double TetrahedronLinearShapeFunction::computeTrialFunction(int element_node_index,
                                                            const std::vector<double>& coords) const
{
    if (coords.size() != 3)
    {
        throw std::invalid_argument(
            "Reference coordinates must have 3 components (xi, eta, zeta).");
    }
    double xi = coords[0];
    double eta = coords[1];
    double zeta = coords[2];
    switch (element_node_index)
    {
        case 0:
            return 1.0 - xi - eta - zeta;
        case 1:
            return xi;
        case 2:
            return eta;
        case 3:
            return zeta;
        default:
            return 0.0;  // 错误情况
    }
}

double TetrahedronLinearShapeFunction::computeTestFunction(int element_node_index,
                                                           const std::vector<double>& coords) const
{
    // 对于伽辽金法，试探函数和检验函数相同
    return computeTrialFunction(element_node_index, coords);
}

/**
 * @brief 线性四面体单元形函数梯度计算
 *
 * 梯度向量 = [dN/dxi, dN/deta, dN/dzeta]^T
 *
 * dN0/dxi = -1,   dN0/deta = -1,   dN0/dzeta = -1
 * dN1/dxi =  1,   dN1/deta =  0,   dN1/dzeta =  0
 * dN2/dxi =  0,   dN2/deta =  1,   dN2/dzeta =  0
 * dN3/dxi =  0,   dN3/deta =  0,   dN3/dzeta =  1
 */
std::vector<double> TetrahedronLinearShapeFunction::computeTrialGradients(
    int element_node_index, const std::vector<double>& coords) const
{
    (void)coords;  // 线性形函数的导数是常数，不依赖于坐标
    switch (element_node_index)
    {
        case 0:
            return {-1.0, -1.0, -1.0};
        case 1:
            return {1.0, 0.0, 0.0};
        case 2:
            return {0.0, 1.0, 0.0};
        case 3:
            return {0.0, 0.0, 1.0};
        default:
            return {0.0, 0.0, 0.0};  // 错误情况
    }
}

std::vector<double> TetrahedronLinearShapeFunction::computeTestGradients(
    int element_node_index, const std::vector<double>& coords) const
{
    // 对于伽辽金法，试探函数和检验函数相同
    return computeTrialGradients(element_node_index, coords);
}
