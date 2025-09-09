#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H

#include <Eigen/Dense>

/**
 * @brief 计算单元内某点的数值解
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param xi 单元内的局部坐标 (-1 到 1)
 * @return double 数值解的值
 */
double numerical_solution_in_element(int element_index, const Eigen::VectorXd& u, double xi);

/**
 * @brief 计算单元内某点的数值解导数
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param xi 单元内的局部坐标 (-1 到 1)
 * @return double 数值解导数的值
 */
double numerical_solution_Derivative_in_element(int element_index, const Eigen::VectorXd& u, double xi);

/**
 * @brief 计算最大误差（无穷范数）
 * @param u 全局解向量
 * @param exact_solution 精确解函数指针
 * @return double 最大误差
 */
double computeMaxError(const Eigen::VectorXd& u, double (*exact_solution)(double));

/**
 * @brief 计算L2范数误差
 * @param u 全局解向量
 * @param exact_solution 精确解函数指针
 * @return double L2范数误差
 */
double computeL2Error(const Eigen::VectorXd& u, double (*exact_solution)(double));

/**
 * @brief 计算H1范数误差
 * @param u 全局解向量
 * @param exact_solution_derivative 精确解导数函数指针
 * @return double H1范数误差
 */
double computeH1Error(const Eigen::VectorXd& u, double (*exact_solution_derivative)(double));

#endif // ERROR_ANALYSIS_H
