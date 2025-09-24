#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include <Eigen/Core>  // 只需要VectorXd，不需要整个Dense模块
#include <memory>
#include <string>
#include "mesh_hierarchy.h"

/**
 * @brief VTK输出基类
 */
class VTKOutput
{
  public:
    VTKOutput() = default;
    virtual ~VTKOutput() = default;

    /**
     * @brief 输出VTK格式文件用于ParaView可视化（纯虚函数）
     *
     * @param filename 输出文件名 (不包含.vtu扩展名)
     * @param solution 求解向量，包含每个节点的数值解
     */
    virtual void outputNumericalSolution(const std::string& filename,
                                         const Eigen::VectorXd& solution) = 0;

    /**
     * @brief 输出包含解析解的VTK文件（纯虚函数）
     *
     * @param filename 输出文件名 (不包含.vtu扩展名)
     * @param exact_func 精确解函数指针
     */
    virtual void outputExactSolution(const std::string& filename,
                                     double (*exact_func)(double, double)) = 0;

    /**
     * @brief 输出包含数值解、解析解和误差的对比VTK文件（纯虚函数）
     *
     * @param filename 输出文件名 (不包含.vtu扩展名)
     * @param solution 数值解向量
     * @param exact_func 精确解函数指针
     */
    virtual void outputComparison(const std::string& filename, const Eigen::VectorXd& solution,
                                  double (*exact_func)(double, double)) = 0;
};

/**
 * @brief 二维VTK输出类，注入Mesh2D
 */
class VTKOutput2D : public VTKOutput
{
  private:
    std::shared_ptr<Mesh2D> mesh2D_;

  public:
    explicit VTKOutput2D(std::shared_ptr<Mesh2D> mesh2D) : mesh2D_(mesh2D)
    {
    }

    /**
     * @brief 输出二维VTK格式文件用于ParaView可视化
     *
     * 输出包含三角形网格和数值解的VTU文件，可以在ParaView中打开进行可视化
     *
     * @param filename 输出文件名 (不包含.vtu扩展名)
     * @param solution 求解向量，包含每个节点的数值解
     */
    void outputNumericalSolution(const std::string& filename,
                                 const Eigen::VectorXd& solution) override;

    /**
     * @brief 输出包含解析解的二维VTK文件
     *
     * @param filename 输出文件名 (不包含.vtu扩展名)
     * @param exact_func 精确解函数指针
     */
    void outputExactSolution(const std::string& filename,
                             double (*exact_func)(double, double)) override;

    /**
     * @brief 输出包含数值解、解析解和误差的对比二维VTK文件
     *
     * @param filename 输出文件名 (不包含.vtu扩展名)
     * @param solution 数值解向量
     * @param exact_func 精确解函数指针
     */
    void outputComparison(const std::string& filename, const Eigen::VectorXd& solution,
                          double (*exact_func)(double, double)) override;

    /**
     * @brief 输出加密采样网格的二维VTK文件
     *
     * @param filename 输出文件名 (不包含.vtu扩展名)
     * @param solution 数值解向量
     * @param exact_func 精确解函数指针
     * @param num_points_per_side 每边采样点数（默认5）
     */
    void outputDenseSamplingError(const std::string& filename, const Eigen::VectorXd& solution,
                                  double (*exact_func)(double, double),
                                  int num_points_per_side = 5);

    /**
     * @brief 获取注入的Mesh2D对象
     *
     * @return std::shared_ptr<Mesh2D> 注入的网格对象
     */
    std::shared_ptr<Mesh2D> getMesh2D() const
    {
        return mesh2D_;
    }
};

/**
 * @brief 确保输出目录存在
 *
 * @param directory 目录路径
 * @return bool 目录创建是否成功
 */
bool ensureDirectoryExists(const std::string& directory);

#endif  // VTK_OUTPUT_H