#include "error_analysis_2d.h"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "fem_solver_2d.h"
#include "gauss_quadrature_2d.h"
#include "geometry_mapping_2d.h"
#include "problem_definition.h"
#include "shape_functions_2d.h"

/**
 * @brief 计算单元内某点的数值解
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param x_ref 单元内的局部坐标x (-1 到 1)
 * @param y_ref 单元内的局部坐标y (-1 到 1)
 * @return double 数值解的值
 */
double numerical_solution_in_element(int element_index, const Eigen::VectorXd& u, double x_ref,
                                     double y_ref)
{
    double numerical_solution = 0.0;
    for (int alpha = 0; alpha < n; ++alpha)  // 遍历单元内所有节点
    {
        int global_node_index = T[element_index][alpha];
        numerical_solution += u(global_node_index) * shapeFunction2D_trial(alpha, x_ref, y_ref);
    }
    return numerical_solution;
}

/**
 * @brief 计算单元内某点的数值解梯度
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param x_ref 单元内的局部坐标x (-1 到 1)
 * @param y_ref 单元内的局部坐标y (-1 到 1)
 * @param numerical_du_dx 输出：x方向导数
 * @param numerical_du_dy 输出：y方向导数
 */
void numerical_solution_gradient_in_element(int element_index, const Eigen::VectorXd& u,
                                            double x_ref, double y_ref, double& numerical_du_dx,
                                            double& numerical_du_dy)
{
    numerical_du_dx = 0.0;
    numerical_du_dy = 0.0;

    // 使用几何映射类计算导数变换
    GeometryMapping2D mapping(element_index);  // 直接使用单元索引构造

    for (int alpha = 0; alpha < n; ++alpha)
    {
        int global_node_index = T[element_index][alpha];
        double dN_dx_ref = shapeFunctionDerivativeXRef_trial(alpha, x_ref, y_ref);
        double dN_dy_ref = shapeFunctionDerivativeYRef_trial(alpha, x_ref, y_ref);
        double dN_dx, dN_dy;
        mapping.transformGradient(dN_dx_ref, dN_dy_ref, dN_dx, dN_dy);
        numerical_du_dx += u(global_node_index) * dN_dx;
        numerical_du_dy += u(global_node_index) * dN_dy;
    }
}

/**
 * @brief 计算最大误差（无穷范数）
 * @param u 全局解向量
 * @param exact_solution 精确解函数指针
 * @return double 最大误差
 */
double computeMaxError(const Eigen::VectorXd& u, double (*exact_solution)(double, double))
{
    double max_error = 0.0;

    // 遍历每个单元
    for (int element_index = 0; element_index < M; ++element_index)
    {
        // 在每个单元内取多个点进行误差计算
        // 使用三角形单元的合理采样点
        int num_points_per_side = 5;  // 每边取5个点
        GeometryMapping2D mapping(element_index);

        for (int i = 0; i < num_points_per_side; ++i)
        {
            for (int j = 0; j < num_points_per_side - i; ++j)
            {
                // 在标准三角形参考单元内采样
                double x_ref = double(i) / (num_points_per_side - 1);
                double y_ref = double(j) / (num_points_per_side - 1);

                // 确保点在三角形内 (x_ref + y_ref <= 1)
                if (x_ref + y_ref <= 1.0)
                {
                    double x, y;
                    mapping.mapToPhysical(x_ref, y_ref, x, y);  // 转换为全局坐标

                    double numerical_val =
                        numerical_solution_in_element(element_index, u, x_ref, y_ref);
                    double exact_val = exact_solution(x, y);
                    double error = std::abs(exact_val - numerical_val);
                    if (error > max_error)
                    {
                        max_error = error;
                    }
                }
            }
        }
    }
    return max_error;
}

/**
 * @brief 计算L2范数误差
 * @param u 全局解向量
 * @param exact_solution 精确解函数指针
 * @return double L2范数误差
 */
double computeL2Error(const Eigen::VectorXd& u, double (*exact_solution)(double, double))
{
    double l2_error = 0.0;

    // 使用高斯积分计算L2误差
    int numGaussPoints = 3;  // 使用3点高斯积分
    std::vector<double> gauss_points, gauss_weights;
    getGaussPointsTriangle(numGaussPoints, gauss_points, gauss_weights);

    // 遍历每个单元
    for (int element_index = 0; element_index < M; ++element_index)
    {
        GeometryMapping2D mapping(element_index);    // 直接使用单元索引构造
        double jacobian = mapping.getJacobianDet();  // 雅可比行列式

        // 在每个单元上进行高斯积分
        for (int gp = 0; gp < numGaussPoints; ++gp)
        {
            double x_ref = gauss_points[2 * gp];      // 局部坐标
            double y_ref = gauss_points[2 * gp + 1];  // 局部坐标
            double weight = gauss_weights[gp];        // 权重
            double x, y;
            mapping.mapToPhysical(x_ref, y_ref, x, y);  // 全局坐标

            double numerical_val = numerical_solution_in_element(element_index, u, x_ref, y_ref);
            double exact_val = exact_solution(x, y);
            double error_at_point = exact_val - numerical_val;
            l2_error += weight * jacobian * error_at_point * error_at_point;
        }
    }
    return std::sqrt(l2_error);
}

/**
 * @brief 计算H1半范数误差
 * @param u 全局解向量
 * @param exact_solution_du_dx 精确解x方向导数函数指针
 * @param exact_solution_du_dy 精确解y方向导数函数指针
 * @return double H1半范数误差
 */
double computeH1Error(const Eigen::VectorXd& u, double (*exact_solution_du_dx)(double, double),
                      double (*exact_solution_du_dy)(double, double))
{
    double h1_error = 0.0;

    // 使用高斯积分计算H1半范数误差
    int numGaussPoints = 3;  // 使用3点高斯积分
    std::vector<double> gauss_points, gauss_weights;
    getGaussPointsTriangle(numGaussPoints, gauss_points, gauss_weights);

    // 遍历每个单元
    for (int element_index = 0; element_index < M; ++element_index)
    {
        GeometryMapping2D mapping(element_index);    // 直接使用单元索引构造
        double jacobian = mapping.getJacobianDet();  // 雅可比行列式

        // 在每个单元上进行高斯积分
        for (int gp = 0; gp < numGaussPoints; ++gp)
        {
            double x_ref = gauss_points[2 * gp];      // 局部坐标
            double y_ref = gauss_points[2 * gp + 1];  // 局部坐标
            double weight = gauss_weights[gp];        // 权重
            double x, y;
            mapping.mapToPhysical(x_ref, y_ref, x, y);  // 全局坐标

            // 计算数值解的梯度
            double numerical_du_dx, numerical_du_dy;
            numerical_solution_gradient_in_element(element_index, u, x_ref, y_ref, numerical_du_dx,
                                                   numerical_du_dy);

            double du_dx_error = exact_solution_du_dx(x, y) - numerical_du_dx;
            double du_dy_error = exact_solution_du_dy(x, y) - numerical_du_dy;

            // H1半范数只包含梯度项
            double integrand = du_dx_error * du_dx_error + du_dy_error * du_dy_error;
            h1_error += weight * jacobian * integrand;
        }
    }

    return std::sqrt(h1_error);
}
