#ifndef GAUSS_QUADRATURE_2D_H
#define GAUSS_QUADRATURE_2D_H

#include <vector>

/**
 * @brief 获取三角形单元上的高斯积分点和权重
 *
 * @param numPoints 积分点数量 (1, 3, 4, 或 7)
 * @param gaussPoints 输出的积分点坐标 [x_ref0, y_ref0, x_ref1, y_ref1, ...]
 * @param gaussWeights 输出的积分权重
 */
void getGaussPointsTriangle(int numPoints, std::vector<double>& gaussPoints,
                            std::vector<double>& gaussWeights);

/**
 * @brief 获取三角形单元上一点积分的积分点和权重（中心点）
 *
 * @param gaussPoints 输出的积分点坐标 [x_ref, y_ref]
 * @param gaussWeights 输出的积分权重
 */
void getGaussPointsTriangle1(std::vector<double>& gaussPoints, std::vector<double>& gaussWeights);

/**
 * @brief 获取三角形单元上三点积分的积分点和权重
 *
 * @param gaussPoints 输出的积分点坐标 [x_ref0, y_ref0, x_ref1, y_ref1, x_ref2, y_ref2]
 * @param gaussWeights 输出的积分权重
 */
void getGaussPointsTriangle3(std::vector<double>& gaussPoints, std::vector<double>& gaussWeights);

/**
 * @brief 获取三角形单元上四点积分的积分点和权重
 *
 * @param gaussPoints 输出的积分点坐标 [x_ref0, y_ref0, x_ref1, y_ref1, x_ref2, y_ref2, x_ref3,
 * y_ref3]
 * @param gaussWeights 输出的积分权重
 */
void getGaussPointsTriangle4(std::vector<double>& gaussPoints, std::vector<double>& gaussWeights);

#endif  // GAUSS_QUADRATURE_2D_H
