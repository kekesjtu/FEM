#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fem_solver.h"
#include "shape_functions.h"
#include "problem_definition.h"
#include "gauss_quadrature.h"
#include "geometry_mapping.h"
#include "error_analysis.h"

//还没有更改成二维！！！！
/**
 * @brief 计算单元内某点的数值解
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param x_ref 单元内的局部坐标 (-1 到 1)
 * @return double 数值解的值
 */
double numerical_solution_in_element(int element_index, const Eigen::VectorXd& u, double x_ref) {
    double numerical_solution = 0.0;
    for (int alpha = 0; alpha < n; ++alpha) {
        int global_node = T[element_index][alpha];
        numerical_solution += u(global_node) * shapeFunction_trial(alpha, x_ref);
    }
    return numerical_solution;
}

/**
 * @brief 计算单元内某点的数值解导数
 * @param element_index 单元索引
 * @param u 全局解向量
 * @param x_ref 单元内的局部坐标 (-1 到 1)
 * @return double 数值解导数的值
 */
double numerical_solution_Derivative_in_element(int element_index, const Eigen::VectorXd& u, double x_ref) {
    double numerical_solution_Derivative = 0.0;
    
    // 使用几何映射类计算导数变换
    GeometryMapping1D mapping(element_index);  // 直接使用单元索引构造
    
    for (int alpha = 0; alpha < n; ++alpha) {
        int global_node = T[element_index][alpha];
        double dN_dx_ref = shapeFunctionDerivative_trial(alpha, x_ref);
        double dN_dx = mapping.transformDerivative(dN_dx_ref);
        numerical_solution_Derivative += u(global_node) * dN_dx;
    }
    return numerical_solution_Derivative;
}

/**
 * @brief 计算最大误差（无穷范数）
 * @param u 全局解向量
 * @param exact_solution 精确解函数指针
 * @return double 最大误差
 */
double computeMaxError(const Eigen::VectorXd& u, double (*exact_solution)(double)) {
    double max_error = 0.0;
    
    // 遍历每个单元
    for (int element_index = 0; element_index < M; ++element_index) {
        // 在每个单元内取多个点进行误差计算
        int num_points = 11; // 在每个单元内取11个点
        for (int i = 0; i < num_points; ++i) {
            double x_ref = -1.0 + 2.0 * i / (num_points - 1); // 局部坐标从-1到1
            GeometryMapping1D mapping(element_index);
            double x = mapping.mapToPhysical(x_ref); // 转换为全局坐标
            
            double numerical_val = numerical_solution_in_element(element_index, u, x_ref);
            double exact_val = exact_solution(x);
            double error = std::abs(exact_val - numerical_val);
            
            if (error > max_error) {
                max_error = error;
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
double computeL2Error(const Eigen::VectorXd& u, double (*exact_solution)(double)) {
    double l2_error = 0.0;
    
    // 使用高斯积分计算L2误差
    int numGaussPoints = 3; // 使用3点高斯积分
    std::vector<double> gauss_points, gauss_weights;
    getGaussPoints(numGaussPoints, gauss_points, gauss_weights);
    
    // 遍历每个单元
    for (int element_index = 0; element_index < M; ++element_index) {
        GeometryMapping1D mapping(element_index);  // 直接使用单元索引构造
        double jacobian = mapping.getJacobian(); // 雅可比行列式
        
        // 在每个单元上进行高斯积分
        for (int gp = 0; gp < numGaussPoints; ++gp) {
            double x_ref = gauss_points[gp]; // 局部坐标
            double weight = gauss_weights[gp]; // 权重
            double x = mapping.mapToPhysical(x_ref); // 全局坐标
            
            double numerical_val = numerical_solution_in_element(element_index, u, x_ref);
            double exact_val = exact_solution(x);
            double error_at_point = exact_val - numerical_val;
            
            l2_error += weight * jacobian * error_at_point * error_at_point;
        }
    }
    
    return std::sqrt(l2_error);
}

/**
 * @brief 计算H1范数误差
 * @param u 全局解向量
 * @param exact_solution_derivative 精确解导数函数指针
 * @return double H1范数误差
 */
double computeH1Error(const Eigen::VectorXd& u, double (*exact_solution_derivative)(double)) {
    double h1_error = 0.0;
    
    // 使用高斯积分计算H1误差
    int numGaussPoints = 3; // 使用3点高斯积分
    std::vector<double> gauss_points, gauss_weights;
    getGaussPoints(numGaussPoints, gauss_points, gauss_weights);
    
    // 遍历每个单元
    for (int element_index = 0; element_index < M; ++element_index) {
        GeometryMapping1D mapping(element_index);  // 直接使用单元索引构造
        double jacobian = mapping.getJacobian(); // 雅可比行列式
        
        // 在每个单元上进行高斯积分
        for (int gp = 0; gp < numGaussPoints; ++gp) {
            double x_ref = gauss_points[gp]; // 局部坐标
            double weight = gauss_weights[gp]; // 权重
            double x = mapping.mapToPhysical(x_ref); // 全局坐标
            
            double numerical_derivative = numerical_solution_Derivative_in_element(element_index, u, x_ref);
            double exact_derivative = exact_solution_derivative(x);
            double error_at_point = exact_derivative - numerical_derivative;
            
            h1_error += weight * jacobian * error_at_point * error_at_point;
        }
    }
    
    return std::sqrt(h1_error);
}
