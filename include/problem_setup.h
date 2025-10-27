#ifndef PROBLEM_SETUP_H
#define PROBLEM_SETUP_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

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
     * @param config 配置对象，用于访问和修改边界条件
     *
     * 此方法调用存储的边界条件设置函数，直接操作 config 的 boundarys 数组
     */
    void setupBoundaryConditions(std::shared_ptr<Config> config);

    /**
     * @brief 设置边界条件设置函数
     * @param func 边界条件设置函数，接收 Config 对象，直接操作其 boundarys 数组
     *
     * 示例：
     * problem->setBoundarySetupFunction([](std::shared_ptr<Config> config) {
     *     auto& boundarys = config->getBoundaryMutable();
     *     const auto& coords = config->getNodeCoordinates();
     *     int dim = config->getDimension();
     *
     *     for (auto& boundary : boundarys) {
     *         // 计算边界中点
     *         double x = 0.0, y = 0.0;
     *         for (int node_idx : boundary.global_node_indices_in_element) {
     *             x += coords[node_idx * dim + 0];
     *             y += coords[node_idx * dim + 1];
     *         }
     *         x /= boundary.global_node_indices_in_element.size();
     *         y /= boundary.global_node_indices_in_element.size();
     *         double r = std::sqrt(x*x + y*y);
     *
     *         // 根据位置设置边界条件
     *         if (r < 0.5) boundary.bc = Config::BoundaryCondition::Dirichlet(0.0);
     *         else boundary.bc = Config::BoundaryCondition::Dirichlet(1.0);
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

  private:
    std::string name_;  // 问题名称

    // 边界条件设置函数（操作 boundarys 数组）
    std::function<void(std::shared_ptr<Config>)> boundary_setup_func_;

    // 物理问题函数
    std::function<double(const std::vector<double>&)> coefficient_func_;  // 扩散系数 c(x,y)
    std::function<double(const std::vector<double>&)> source_func_;       // 源项 f(x,y)

    // 精确解函数（可选，用于误差分析）
    std::function<double(const std::vector<double>&)> exact_solution_u_func_;
    std::function<double(const std::vector<double>&, std::vector<double>&)>
        exact_solution_gradients_func_;
};

#endif  // PROBLEM_SETUP_H
