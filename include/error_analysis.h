#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>
#include "config.h"

class ErrorAnalysis
{
  public:
    enum class NormType
    {
        L_INFINITY,  // L∞范数误差
        L2,          // L2范数误差
        H1_SEMINORM  // H1半范数误差
    };

    // 构造函数，只依赖Config和数值解
    ErrorAnalysis(std::shared_ptr<Config> config, const Eigen::VectorXd& solution);

    ~ErrorAnalysis() = default;

    // 输出误差摘要
    void printErrorSummary() const;
    // 打印详细的节点误差
    void printDetailedNodeErrors(int max_nodes = 10) const;

    // 计算范数误差
    double computeNormError(NormType norm) const;
    double computeLInfinityError() const;
    double computeL2Error() const;
    double computeH1SeminormError() const;

  private:
    // 内部辅助函数，检查Config中是否定义了精确解
    bool hasExactSolution() const;
    bool hasExactGradients() const;

    double calculateNumericalSolution(int element_idx, const std::vector<double>& coords) const;
    std::vector<double> calculateNumericalGradient(int element_idx,
                                                   const std::vector<double>& coords) const;

  protected:
    std::shared_ptr<Config> config_;
    const Eigen::VectorXd& solution_;
};

#endif  // ERROR_ANALYSIS_H