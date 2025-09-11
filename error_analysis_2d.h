#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H

#include <Eigen/Dense>

/**
 * @brief 计算单元内某点的数值解
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param x_ref 单元内的局部坐标x (-1 到 1)
 * @param y_ref 单元内的局部坐标y (-1 到 1)
 * @return double 数值解的值
 */
double numerical_solution_in_element(int element_index, const Eigen::VectorXd& u, double x_ref,
                                     double y_ref);

/**
 * @brief 计算单元内某点的数值解梯度
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param x_ref 单元内的局部坐标x (-1 到 1)
 * @param y_ref 单元内的局部坐标y (-1 到 1)
 * @param grad_x 输出：x方向导数
 * @param grad_y 输出：y方向导数
 */
void numerical_solution_gradient_in_element(int element_index, const Eigen::VectorXd& u,
                                            double x_ref, double y_ref, double& grad_x,
                                            double& grad_y);

/**
 * @brief 计算最大误差（无穷范数）
 * @param u 全局解向量
 * @param exact_solution 精确解函数指针
 * @return double 最大误差
 */
double computeMaxError(const Eigen::VectorXd& u, double (*exact_solution)(double, double));

/**
 * @brief 计算L2范数误差
 * @param u 全局解向量
 * @param exact_solution 精确解函数指针
 * @return double L2范数误差
 */
double computeL2Error(const Eigen::VectorXd& u, double (*exact_solution)(double, double));

/**
 * @brief 计算H1半范数误差（只包含梯度项）
 * @param u 全局解向量
 * @param exact_solution_du_dx 精确解x方向导数函数指针
 * @param exact_solution_du_dy 精确解y方向导数函数指针
 * @return double H1半范数误差
 */
double computeH1Error(const Eigen::VectorXd& u, double (*exact_solution_du_dx)(double, double),
                      double (*exact_solution_du_dy)(double, double));

#endif  // ERROR_ANALYSIS_H
