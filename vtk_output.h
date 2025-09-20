#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include <Eigen/Dense>
#include <string>

/**
 * @brief 输出VTK格式文件用于ParaView可视化
 *
 * 输出包含三角形网格和数值解的VTU文件，可以在ParaView中打开进行可视化
 *
 * @param filename 输出文件名 (不包含.vtu扩展名)
 * @param solution 求解向量，包含每个节点的数值解
 */
void outputVTKNumericalSolution(const std::string& filename, const Eigen::VectorXd& solution);

/**
 * @brief 输出包含解析解的VTK文件
 *
 * @param filename 输出文件名 (不包含.vtu扩展名)
 * @param exact_func 精确解函数指针
 */
void outputVTKExactSolution(const std::string& filename, double (*exact_func)(double, double));

/**
 * @brief 输出包含数值解、解析解和误差的对比VTK文件
 *
 * @param filename 输出文件名 (不包含.vtu扩展名)
 * @param solution 数值解向量
 * @param exact_func 精确解函数指针
 */
void outputVTKComparison(const std::string& filename, const Eigen::VectorXd& solution,
                         double (*exact_func)(double, double));

/**
 * @brief 输出加密采样网格的VTK文件（使用与maxerror相同的采样策略）
 *
 * 在每个原始单元内部按照computeMaxError的采样方式生成密集的采样点，
 * 然后为每个采样点计算数值解和精确解，以节点方式输出到VTK文件中
 *
 * @param filename 输出文件名 (不包含.vtu扩展名)
 * @param solution 数值解向量
 * @param exact_func 精确解函数指针
 * @param num_points_per_side 每边采样点数（默认5，与computeMaxError一致）
 */
void outputVTKDenseSamplingError(const std::string& filename, const Eigen::VectorXd& solution,
                                 double (*exact_func)(double, double), int num_points_per_side = 5);

/**
 * @brief 确保输出目录存在
 *
 * @param directory 目录路径
 * @return bool 目录创建是否成功
 */
bool ensureDirectoryExists(const std::string& directory);

#endif  // VTK_OUTPUT_H