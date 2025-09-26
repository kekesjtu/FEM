#ifndef ERROR_ANALYSIS_REDESIGNED_H
#define ERROR_ANALYSIS_REDESIGNED_H

#include <Eigen/Dense>
#include <map>
#include <memory>
#include <string>
#include "config.h"

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

  protected:
    std::shared_ptr<Config> config_;

    const Eigen::VectorXd& solution;

  public:
    // 参数化构造函数
    ErrorAnalysis(std::shared_ptr<Config> config, const Eigen::VectorXd& solution)
        : config_(config), solution(solution)
    {
    }

    virtual ~ErrorAnalysis() = default;

    // 输出误差摘要
    virtual void printErrorSummary() const = 0;
    // 输出函数
    void printErrorSummary() const;
    void printDetailedNodeErrors(const Eigen::VectorXd& solution, int max_nodes = 10) const;

  private:
    bool hasExactSolution() const
    {
        return exact_solution_ != nullptr;  // 完善报错
    }
    bool hasExactGradients() const
    {
        return exact_gradients_ != nullptr;  // 完善报错
    }
    virtual double calculateNumericalSolution(int element_idx, const Eigen::VectorXd& solution,
                                              std::vector<double> coords) const;
    virtual void calculateNumericalGradient(int element_idx, const Eigen::VectorXd& solution,
                                            std::vector<double> coords,
                                            std::vector<double>& gradients) const;
    virtual bool supportsNorm(NormType norm) const;
    virtual double computeNormError(NormType norm, const Eigen::VectorXd& solution) const = 0;
    virtual double computeLInfinityError(const Eigen::VectorXd& solution) const = 0;
    virtual double computeL2Error(const Eigen::VectorXd& solution) const = 0;
    virtual double computeH1SeminormError(const Eigen::VectorXd& solution) const = 0;
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
    /**
     * @brief 构造函数
     * 维度固定为2，注入mesh依赖，专注于精确解设置
     */
    ErrorAnalysis2D(std::shared_ptr<Mesh2D> mesh, std::shared_ptr<ComputeConfig> config,
                    std::shared_ptr<ProblemDef> problem_def, const Eigen::VectorXd& solution)
        : ErrorAnalysis(mesh, config, problem_def, solution)
    {
    }

    virtual ~ErrorAnalysis2D() = default;

    void printErrorSummary() const override;

  private:
    double calculateNumericalSolution(int element_idx, const Eigen::VectorXd& solution,
                                      std::vector<double> coords) const;
    void calculateNumericalGradient(int element_idx, const Eigen::VectorXd& solution,
                                    std::vector<double> coords,
                                    std::vector<double>& gradients) const;
    bool supportsNorm(NormType norm) const;
    double computeNormError(NormType norm, const Eigen::VectorXd& solution) const = 0;
    double computeLInfinityError(const Eigen::VectorXd& solution) const;
    double computeL2Error(const Eigen::VectorXd& solution) const;
    double computeH1SeminormError(const Eigen::VectorXd& solution) const;
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
    static std::unique_ptr<ErrorAnalysis2D> createAnalyzer(
        std::shared_ptr<Mesh2D> mesh, ErrorAnalysis2D::ExactSolutionFunc exact_sol = nullptr,
        ErrorAnalysis2D::ExactGradientFunc exact_dx = nullptr,
        ErrorAnalysis2D::ExactGradientFunc exact_dy = nullptr);
};