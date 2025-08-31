#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTIONS_H

/**
 * @brief 计算一维线性单元的试探函数值 N_alpha(xi)
 * 
 * @param alpha 形函数的局部索引 (0 或 1)
 * @param xi 单元内的局部坐标 (-1 到 1)
 * @return double 形函数的值
 */
double shapeFunction_trial(int alpha, double xi);

/**
 * @brief 计算一维线性单元的试探函数导数 dN_alpha/dxi
 * 
 * @param alpha 形函数的局部索引 (0 或 1)
 * @param xi 单元内的局部坐标 (-1 到 1)
 * @return double 形函数的导数值
 */
double shapeFunctionDerivative_trial(int alpha, double xi);

/**
 * @brief 计算一维线性单元的检验函数值 N_beta(xi)
 *        (在伽辽金法中，检验函数与试探函数相同)
 * 
 * @param beta 形函数的局部索引 (0 或 1)
 * @param xi 单元内的局部坐标 (-1 到 1)
 * @return double 形函数的值
 */
double shapeFunction_test(int beta, double xi);

/**
 * @brief 计算一维线性单元的检验函数导数 dN_beta/dxi
 *        (在伽辽金法中，检验函数与试探函数相同)
 * 
 * @param beta 形函数的局部索引 (0 或 1)
 * @param xi 单元内的局部坐标 (-1 到 1)
 * @return double 形函数的导数值
 */
double shapeFunctionDerivative_test(int beta, double xi);


#endif // SHAPE_FUNCTIONS_H
