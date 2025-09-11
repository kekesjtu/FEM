#ifndef SHAPE_FUNCTIONS_2D_H
#define SHAPE_FUNCTIONS_2D_H

/**
 * @brief 计算二维线性三角形单元的试探函数值 N_alpha(x_ref, y_ref)
 *
 * @param alpha 形函数的局部索引 (0, 1, 或 2)
 * @param x_ref 单元内的第一个参考坐标
 * @param y_ref 单元内的第二个参考坐标
 * @return double 形函数的值
 */
double shapeFunction2D_trial(int alpha, double x_ref, double y_ref);

/**
 * @brief 计算二维线性三角形单元的试探函数对x_ref的偏导数 dN_alpha/dx_ref
 *
 * @param alpha 形函数的局部索引 (0, 1, 或 2)
 * @param x_ref 单元内的第一个参考坐标
 * @param y_ref 单元内的第二个参考坐标
 * @return double 形函数对x_ref的偏导数值
 */
double shapeFunctionDerivativeXRef_trial(int alpha, double x_ref, double y_ref);

/**
 * @brief 计算二维线性三角形单元的试探函数对y_ref的偏导数 dN_alpha/dy_ref
 *
 * @param alpha 形函数的局部索引 (0, 1, 或 2)
 * @param x_ref 单元内的第一个参考坐标
 * @param y_ref 单元内的第二个参考坐标
 * @return double 形函数对y_ref的偏导数值
 */
double shapeFunctionDerivativeYRef_trial(int alpha, double x_ref, double y_ref);

/**
 * @brief 计算二维线性三角形单元的检验函数值 N_beta(x_ref, y_ref)
 *        (在伽辽金法中，检验函数与试探函数相同)
 *
 * @param beta 形函数的局部索引 (0, 1, 或 2)
 * @param x_ref 单元内的第一个参考坐标
 * @param y_ref 单元内的第二个参考坐标
 * @return double 形函数的值
 */
double shapeFunction2D_test(int beta, double x_ref, double y_ref);

/**
 * @brief 计算二维线性三角形单元的检验函数对x_ref的偏导数 dN_beta/dx_ref
 *        (在伽辽金法中，检验函数与试探函数相同)
 *
 * @param beta 形函数的局部索引 (0, 1, 或 2)
 * @param x_ref 单元内的第一个参考坐标
 * @param y_ref 单元内的第二个参考坐标
 * @return double 形函数对x_ref的偏导数值
 */
double shapeFunctionDerivativeXRef_test(int beta, double x_ref, double y_ref);

/**
 * @brief 计算二维线性三角形单元的检验函数对y_ref的偏导数 dN_beta/dy_ref
 *        (在伽辽金法中，检验函数与试探函数相同)
 *
 * @param beta 形函数的局部索引 (0, 1, 或 2)
 * @param x_ref 单元内的第一个参考坐标
 * @param y_ref 单元内的第二个参考坐标
 * @return double 形函数对y_ref的偏导数值
 */
double shapeFunctionDerivativeYRef_test(int beta, double x_ref, double y_ref);

#endif  // SHAPE_FUNCTIONS_2D_H
