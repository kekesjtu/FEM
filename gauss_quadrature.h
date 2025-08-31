#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H

#include <vector>

/**
 * @brief 获取一维高斯-勒让德积分点和权重
 * 
 * @param numPoints 需要的积分点数量 (支持 1 到 5 点)
 * @param points 输出参数，存储积分点的局部坐标 (xi)
 * @param weights 输出参数，存储对应的权重
 */
void getGaussPoints(int numPoints, std::vector<double>& points, std::vector<double>& weights);

#endif // GAUSS_QUADRATURE_H
