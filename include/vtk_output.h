#ifndef VTK_OUTPUT_H
#define VTK_OUTPUT_H

#include <Eigen/Core>
#include <functional>
#include <memory>
#include <string>
#include "config.h"
#include "problem_setup.h"

/**
 * @brief VTK输出基类
 *
 * 负责将有限元计算的结果（数值解、精确解、误差）输出为VTK格式文件，
 * 以便在ParaView等后处理软件中进行可视化。
 *
 * 设计原则：此类存储从Config对象中提取的、用于可视化的核心网格数据副本。
 * 这种设计使其在数据层面自包含，便于单独测试和调试。
 * 对于需要复杂计算（如密集采样）的方法，将Config对象作为参数传入，
 * 以便在方法内部访问工厂来创建所需组件（如GeometryMapping）。
 */
class VTKOutput
{
  protected:
    // 保存全局配置，所有网格与问题定义从此处读取
    std::shared_ptr<Config> config_;
    // 数值解向量（与网格节点一一对应）
    const Eigen::VectorXd& solution_;
    // 每个单元边上采样点数（用于加密采样）
    int num_points_per_side_ = 5;

  public:
    VTKOutput(std::shared_ptr<Config> config, const Eigen::VectorXd& solution,
              int num_points_per_side)
        : config_(std::move(config)), solution_(solution), num_points_per_side_(num_points_per_side)
    {
    }

    virtual ~VTKOutput() = default;

    /**
     * @brief 输出数值解的VTK文件
     * @param filename 输出文件名 (不含.vtu)
     */
    virtual void outputNumericalSolution(const std::string& filename) = 0;

    /**
     * @brief 输出精确解的VTK文件
     * @param filename 输出文件名 (不含.vtu)
     * @param exact_func 用于计算精确解的函数
     */
    virtual void outputExactSolution(
        const std::string& filename,
        const std::function<double(const std::vector<double>&)>& exact_func) = 0;

    /**
     * @brief 输出加密采样的误差VTK文件
     * @param filename 输出文件名 (不含.vtu)
     * @param config 配置对象，用于访问工厂
     * @param problem 问题设置对象，用于获取精确解
     */
    virtual void outputDenseSamplingError(const std::string& filename,
                                          std::shared_ptr<Config> config,
                                          std::shared_ptr<ProblemSetup> problem) = 0;

  protected:
    bool ensureDirectoryExists(const std::string& directory);
};

/**
 * @brief 二维VTK输出类
 */
class VTKOutput2D : public VTKOutput
{
  public:
    VTKOutput2D(std::shared_ptr<Config> config, const Eigen::VectorXd& solution,
                int num_points_per_side)
        : VTKOutput(std::move(config), solution, num_points_per_side)
    {
    }
};

/**
 * @brief 三维VTK输出类
 */
class VTKOutput3D : public VTKOutput
{
  public:
    VTKOutput3D(std::shared_ptr<Config> config, const Eigen::VectorXd& solution,
                int num_points_per_side)
        : VTKOutput(std::move(config), solution, num_points_per_side)
    {
    }
};

/**
 * @brief 二维三角形单元VTK输出类
 */
class TriangleVTKOutput2D : public VTKOutput2D
{
  public:
    TriangleVTKOutput2D(std::shared_ptr<Config> config, const Eigen::VectorXd& solution,
                        int num_points_per_side)
        : VTKOutput2D(std::move(config), solution, num_points_per_side)
    {
    }

    void outputNumericalSolution(const std::string& filename) override;

    void outputExactSolution(
        const std::string& filename,
        const std::function<double(const std::vector<double>&)>& exact_func) override;

    void outputDenseSamplingError(const std::string& filename, std::shared_ptr<Config> config,
                                  std::shared_ptr<ProblemSetup> problem) override;

  private:
    // 辅助函数，计算单元内任意参考坐标点的数值解
    double evaluateNumericalSolutionAt(int element_index, const std::vector<double>& coords_ref,
                                       std::shared_ptr<Config> config) const;
};

/**
 * @brief 三维四面体单元VTK输出类
 */
class TetrahedronVTKOutput3D : public VTKOutput3D
{
  public:
    TetrahedronVTKOutput3D(std::shared_ptr<Config> config, const Eigen::VectorXd& solution,
                           int num_points_per_side)
        : VTKOutput3D(std::move(config), solution, num_points_per_side)
    {
    }

    void outputNumericalSolution(const std::string& filename) override;

    void outputExactSolution(
        const std::string& filename,
        const std::function<double(const std::vector<double>&)>& exact_func) override;

    void outputDenseSamplingError(const std::string& filename, std::shared_ptr<Config> config,
                                  std::shared_ptr<ProblemSetup> problem) override;

  private:
    // 辅助函数，计算单元内任意参考坐标点的数值解
    double evaluateNumericalSolutionAt(int element_index, const std::vector<double>& coords_ref,
                                       std::shared_ptr<Config> config) const;
};

/**
 * @brief VTK输出工厂类
 */
class VTKOutputFactory
{
  public:
    /**
     * @brief 根据Config创建VTK输出对象
     * @param config 配置对象
     * @param solution 数值解向量
     * @return 指向创建的VTK输出对象的 unique_ptr
     */
    static std::unique_ptr<VTKOutput> createVTKOutput(std::shared_ptr<Config> config,
                                                      const Eigen::VectorXd& solution);
};

#endif  // VTK_OUTPUT_H