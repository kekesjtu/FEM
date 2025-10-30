#include <stdexcept>
#include "shape_functions.h"


// ============================================================================
// 一维线性线段单元形函数实现
// ============================================================================

/**
 * @brief 线性线段单元形函数值计算
 *
 * 对于标准参考线段，两个端点为：
 * 节点0: t = -1
 * 节点1: t = +1
 *
 * 形函数为：
 * N0(t) = (1 - t) / 2
 * N1(t) = (1 + t) / 2
 */
double LineLinearShapeFunction::computeTrialFunction(int element_node_index,
                                                     const std::vector<double>& coords) const
{
    if (coords.size() != 1)
    {
        throw std::invalid_argument("Reference coordinates must have 1 component (t).");
    }
    double t = coords[0];
    switch (element_node_index)
    {
        case 0:
            return (1.0 - t) / 2.0;
        case 1:
            return (1.0 + t) / 2.0;
        default:
            return 0.0;  // 错误情况
    }
}

double LineLinearShapeFunction::computeTestFunction(int element_node_index,
                                                    const std::vector<double>& coords) const
{
    // 对于伽辽金法，试探函数和检验函数相同
    return computeTrialFunction(element_node_index, coords);
}

/**
 * @brief 线性线段单元形函数对t的偏导数
 * dN0/dt = -1/2
 * dN1/dt = +1/2
 */
std::vector<double> LineLinearShapeFunction::computeTrialGradients(
    int element_node_index, const std::vector<double>& coords) const
{
    (void)coords;  // 线性形函数的导数是常数，不依赖于坐标
    switch (element_node_index)
    {
        case 0:
            return {-0.5};
        case 1:
            return {0.5};
        default:
            return {0.0};  // 错误情况
    }
}

std::vector<double> LineLinearShapeFunction::computeTestGradients(
    int element_node_index, const std::vector<double>& coords) const
{
    // 对于伽辽金法，试探函数和检验函数相同
    return computeTrialGradients(element_node_index, coords);
}
