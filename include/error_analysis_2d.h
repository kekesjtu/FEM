#ifndef ERROR_ANALYSIS_REDESIGNED_H
#define ERROR_ANALYSIS_REDESIGNED_H

#include <Eigen/Dense>
#include <map>
#include <memory>
#include <string>
#include "mesh_hierarchy.h"

/**
 * @brief 误差分析基类 - 重新设计
 *
 * 核心理念：
 * 1. 误差范数的数学定义与维度无关，只是积分域不同
 * 2. 维度信息应该在构造时确定，不需要每个类都管理
 * 3. 提供统一的误差计算接口
 */
class ErrorAnalysis
{
  public:
    /**
     * @brief 误差范数类型
     */
    enum class NormType
    {
        L_INFINITY,  // L∞范数误差（最大误差）
        L2,          // L2范数误差
        H1_SEMINORM  // H1半范数误差（只包含导数项）
    };

    /**
     * @brief 误差结果容器
     * 使用 map 存储不同范数的误差值
     */
    using ErrorResults = std::map<NormType, double>;

  protected:
    ErrorResults results_;  // 误差结果
    const int dimension_;   // 空间维数（构造时确定，不可变）

  public:
    explicit ErrorAnalysis(int dim) : dimension_(dim)
    {
    }
    virtual ~ErrorAnalysis() = default;

    // 获取基本信息
    int getDimension() const
    {
        return dimension_;
    }
    const ErrorResults& getResults() const
    {
        return results_;
    }

    // 获取特定范数误差
    double getError(NormType norm) const
    {
        auto it = results_.find(norm);
        return (it != results_.end()) ? it->second : 0.0;
    }

    // 清空结果
    void clearResults()
    {
        results_.clear();
    }

    // 核心接口 - 计算指定范数误差
    virtual double computeNormError(NormType norm, const Eigen::VectorXd& solution) = 0;

    // 批量计算所有支持的范数误差
    virtual void computeAllNormErrors(const Eigen::VectorXd& solution) = 0;

    // 输出误差摘要
    virtual void printErrorSummary() const;

    // 检查是否支持某种范数
    virtual bool supportsNorm(NormType norm) const = 0;
};

/**
 * @brief 二维误差分析类
 *
 * 专门处理二维问题的误差计算
 * 重点：不再管理维度，专注于二维特定的计算逻辑
 */
class ErrorAnalysis2D : public ErrorAnalysis
{
  public:
    // 精确解函数类型
    using ExactSolutionFunc = double (*)(double x, double y);
    using ExactGradientFunc = double (*)(double x, double y);

  private:
    // 网格依赖（解耦全局变量）
    std::shared_ptr<Mesh2D> mesh_;

    // 精确解相关
    ExactSolutionFunc exact_solution_;
    ExactGradientFunc exact_du_dx_;
    ExactGradientFunc exact_du_dy_;

    // 计算参数
    int max_error_sampling_points_;  // L∞误差采样点数
    int gauss_integration_points_;   // 数值积分点数

  public:
    /**
     * @brief 构造函数
     * 维度固定为2，注入mesh依赖，专注于精确解设置
     */
    ErrorAnalysis2D(std::shared_ptr<Mesh2D> mesh, ExactSolutionFunc exact_sol = nullptr,
                    ExactGradientFunc exact_dx = nullptr, ExactGradientFunc exact_dy = nullptr);

    virtual ~ErrorAnalysis2D() = default;

    // 设置精确解
    void setExactSolution(ExactSolutionFunc func)
    {
        exact_solution_ = func;
    }
    void setExactGradients(ExactGradientFunc du_dx, ExactGradientFunc du_dy)
    {
        exact_du_dx_ = du_dx;
        exact_du_dy_ = du_dy;
    }

    // 设置计算参数
    void setSamplingPoints(int points)
    {
        max_error_sampling_points_ = points;
    }
    void setGaussPoints(int points)
    {
        gauss_integration_points_ = points;
    }

    // 实现基类接口
    double computeNormError(NormType norm, const Eigen::VectorXd& solution) override;
    void computeAllNormErrors(const Eigen::VectorXd& solution) override;
    bool supportsNorm(NormType norm) const override;

    // 二维特定的误差计算方法
    double computeLInfinityError(const Eigen::VectorXd& solution);
    double computeL2Error(const Eigen::VectorXd& solution);
    double computeH1SeminormError(const Eigen::VectorXd& solution);

    // 辅助计算函数
    double evaluateNumericalSolution(int element_idx, const Eigen::VectorXd& solution, double xi,
                                     double eta) const;
    void evaluateNumericalGradient(int element_idx, const Eigen::VectorXd& solution, double xi,
                                   double eta, double& grad_x, double& grad_y) const;

    // 网格信息访问（避免直接暴露mesh指针）
    int getNumNodes() const
    {
        return mesh_->getNumNodes();
    }
    int getNumElements() const
    {
        return mesh_->getNumElements();
    }
    int getNodesPerElement() const
    {
        return mesh_->getNodesPerElement();
    }
    std::vector<double> getElementNodes(int element_id) const
    {
        return mesh_->getElementNodes(element_id);
    }
    const std::vector<std::vector<int>>& getConnectivity() const
    {
        return mesh_->getElementConnectivity();
    }
    const std::vector<double>& getNodeCoordinates() const
    {
        return mesh_->getNodeCoordinates();
    }

    // 输出函数
    void printErrorSummary() const override;
    void printDetailedNodeErrors(const Eigen::VectorXd& solution, int max_nodes = 10) const;

  private:
    bool hasExactSolution() const
    {
        return exact_solution_ != nullptr;
    }
    bool hasExactGradients() const
    {
        return exact_du_dx_ != nullptr && exact_du_dy_ != nullptr;
    }
};

/**
 * @brief 误差分析工厂
 *
 * 简化对象创建，避免直接调用构造函数
 */
class ErrorAnalysisFactory
{
  public:
    /**
     * @brief 创建二维误差分析器
     */
    static std::unique_ptr<ErrorAnalysis2D> create2D(
        std::shared_ptr<Mesh2D> mesh, ErrorAnalysis2D::ExactSolutionFunc exact_sol = nullptr,
        ErrorAnalysis2D::ExactGradientFunc exact_dx = nullptr,
        ErrorAnalysis2D::ExactGradientFunc exact_dy = nullptr);

    /**
     * @brief 创建带完整精确解的二维误差分析器
     */
    static std::unique_ptr<ErrorAnalysis2D> create2DWithFullSolution(
        std::shared_ptr<Mesh2D> mesh, ErrorAnalysis2D::ExactSolutionFunc exact_sol,
        ErrorAnalysis2D::ExactGradientFunc exact_dx, ErrorAnalysis2D::ExactGradientFunc exact_dy);
};

/**
 * @brief 全局配置类（解决维度管理问题）
 *
 * 如果多个类都需要维度信息，可以使用这个单例类统一管理
 */
class FEMConfig
{
  private:
    static FEMConfig* instance_;
    int dimension_;

    FEMConfig(int dim) : dimension_(dim)
    {
    }

  public:
    static void initialize(int dimension)
    {
        if (!instance_)
        {
            instance_ = new FEMConfig(dimension);
        }
    }

    static FEMConfig& getInstance()
    {
        if (!instance_)
        {
            throw std::runtime_error("FEMConfig not initialized!");
        }
        return *instance_;
    }

    int getDimension() const
    {
        return dimension_;
    }

    static void cleanup()
    {
        delete instance_;
        instance_ = nullptr;
    }
};

#endif  // ERROR_ANALYSIS_REDESIGNED_H