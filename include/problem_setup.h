#ifndef PROBLEM_SETUP_H
#define PROBLEM_SETUP_H

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "boundary_condition.h"

// 前向声明以避免循环包含
class Config;

/**
 * @brief 物理问题设置类
 *
 * 负责定义具体物理场的边界条件、精确解等与物理问题相关的设置
 * 与 Config 类分离，Config 只负责网格、积分点等通用配置
 */
class ProblemSetup
{
  public:
    /**
     * @brief 构造函数
     * @param name 问题名称（用于标识，如 "Electric Field", "Thermal Field"）
     */
    explicit ProblemSetup(const std::string& name = "Unnamed Problem");

    /**
     * @brief 设置边界条件（应该在求解前调用）
     * @param config 配置对象，用于访问边界几何信息
     *
     * 子类必须实现此方法，根据边界几何信息设置物理边界条件
     * 边界条件将存储在 boundary_conditions_ 中，索引对应 config->getBoundary()
     */
    virtual void setupBoundaryConditions(std::shared_ptr<Config> config) = 0;

    /**
     * @brief 设置边界条件设置函数（Lambda方式，供便捷使用）
     * @param func 边界条件设置函数，接收 Config 对象
     *
     * 示例：
     * problem->setBoundarySetupFunction([](std::shared_ptr<Config> config) {
     *     auto& boundary_conditions = problem->getBoundaryConditionsMutable();
     *     const auto& boundaries = config->getBoundary();
     *     const auto& coords = config->getNodeCoordinates();
     *     int dim = config->getDimension();
     *
     *     boundary_conditions.resize(boundaries.size());
     *     for (size_t i = 0; i < boundaries.size(); ++i) {
     *         const auto& boundary = boundaries[i];
     *         // 计算边界中点
     *         double x = 0.0, y = 0.0;
     *         for (int node_idx : boundary.global_node_indices_in_element) {
     *             x += coords[node_idx * dim + 0];
     *             y += coords[node_idx * dim + 1];
     *         }
     *         x /= boundary.global_node_indices_in_element.size();
     *         y /= boundary.global_node_indices_in_element.size();
     *
     *         // 根据位置设置边界条件
     *         if (x < 0) boundary_conditions[i] = BoundaryCondition::Dirichlet(0.0);
     *         else boundary_conditions[i] = BoundaryCondition::Dirichlet(1.0);
     *     }
     * });
     */
    void setBoundarySetupFunction(std::function<void(std::shared_ptr<Config>)> func);

    /**
     * @brief 设置扩散系数函数 c(x,y)
     * @param func 系数函数，输入物理坐标，返回系数值
     */
    void setCoefficient(std::function<double(const std::vector<double>&)> func);

    /**
     * @brief 设置源项函数 f(x,y)
     * @param func 源项函数，输入物理坐标，返回源项值
     */
    void setSource(std::function<double(const std::vector<double>&)> func);

    /**
     * @brief 设置精确解函数（用于误差分析，可选）
     * @param u_func 精确解函数 u(x,y,...)
     */
    void setExactSolutionU(std::function<double(const std::vector<double>&)> func);

    /**
     * @brief 设置精确解梯度函数（用于误差分析，可选）
     * @param grad_func 精确解梯度函数 ∇u(x,y,...)
     */
    void setExactSolutionGradients(
        std::function<double(const std::vector<double>&, std::vector<double>&)> func);

    /**
     * @brief 计算扩散系数 c(x,y)
     * @param coords 点坐标
     * @return 扩散系数值
     */
    double coefficient(const std::vector<double>& coords) const;

    /**
     * @brief 计算源项 f(x,y)
     * @param coords 点坐标
     * @return 源项值
     */
    double source(const std::vector<double>& coords) const;

    /**
     * @brief 计算精确解（如果已设置）
     * @param coords 点坐标
     * @return 精确解值
     */
    double exactSolutionU(const std::vector<double>& coords) const;

    /**
     * @brief 计算精确解梯度（如果已设置）
     * @param coords 点坐标
     * @param gradients 输出梯度向量
     * @return 精确解值
     */
    double exactSolutionGradients(const std::vector<double>& coords,
                                  std::vector<double>& gradients) const;

    /**
     * @brief 获取边界条件列表（供 FEMSolver 使用）
     * @return 边界条件列表，索引对应 Config::getBoundary()
     */
    const std::vector<BoundaryCondition>& getBoundaryConditions() const;

    /**
     * @brief 获取指定边界的边界条件
     * @param boundary_index 边界索引（对应 Config::getBoundary() 中的索引）
     * @return 边界条件
     */
    const BoundaryCondition& getBoundaryCondition(size_t boundary_index) const;

    /**
     * @brief 获取可修改的边界条件列表（供子类使用）
     * @return 边界条件列表的可修改引用
     */
    std::vector<BoundaryCondition>& getBoundaryConditionsMutable();

    /**
     * @brief 检查是否设置了精确解
     */
    bool hasExactSolution() const;

    /**
     * @brief 获取问题名称
     */
    std::string getName() const
    {
        return name_;
    }

  protected:
    std::string name_;  ///< 问题名称

    /**
     * @brief 边界条件列表（索引对应 Config 中的 Boundary）
     *
     * boundary_conditions_[i] 对应 config->getBoundary()[i]
     * 由子类在 setupBoundaryConditions() 中设置
     */
    std::vector<BoundaryCondition> boundary_conditions_;

    // 边界条件设置函数（Lambda方式，可选）
    std::function<void(std::shared_ptr<Config>)> boundary_setup_func_;

  private:
    // 物理问题函数
    std::function<double(const std::vector<double>&)> coefficient_func_;  // 扩散系数 c(x,y)
    std::function<double(const std::vector<double>&)> source_func_;       // 源项 f(x,y)

    // 精确解函数（可选，用于误差分析）
    std::function<double(const std::vector<double>&)> exact_solution_u_func_;
    std::function<double(const std::vector<double>&, std::vector<double>&)>
        exact_solution_gradients_func_;
};

#endif  // PROBLEM_SETUP_H
