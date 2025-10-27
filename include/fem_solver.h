#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <Eigen/Sparse>
#include <memory>
#include <vector>
#include "config.h"
#include "problem_setup.h"

/**
 * @brief 有限元求解器类
 *
 * 封装了完整的有限元求解流程，支持自定义系数和源项函数。
 * 适用于求解形式为 -∇·(c∇u) = f 的偏微分方程。
 */
class FEMSolver
{
  private:
    // 配置对象
    std::shared_ptr<Config> config_;
    std::shared_ptr<ProblemSetup> problem_;  // 物理问题设置

    // 网格信息
    int N_;  // 节点数
    int M_;  // 单元数
    int n_;  // 每个单元的节点数

    // 有限元矩阵和向量
    Eigen::SparseMatrix<double> K_global_;  // 全局刚度矩阵
    Eigen::VectorXd b_;                     // 载荷向量
    Eigen::VectorXd u_;                     // 解向量

    // 内部辅助函数
    double calculateStiffnessEntry(int e, int alpha, int beta);
    double calculateLoadEntry(int e, int beta);
    double calculateBoundaryStiffness(const Config::Boundary& boundary, int local_i, int local_j);
    double calculateBoundaryLoad(const Config::Boundary& boundary, int local_i);

  public:
    /**
     * @brief 构造函数
     * @param config 配置对象，包含网格信息
     * @param problem 问题设置，包含边界条件、系数、源项和精确解
     */
    FEMSolver(std::shared_ptr<Config> config, std::shared_ptr<ProblemSetup> problem);

    // === 高级接口：一键求解 ===

    /**
     * @brief 完整求解流程（推荐用于单场问题）
     *
     * 依次执行：preprocess -> assemble -> applyBoundaryConditions -> solve -> postprocess
     */
    void solveComplete();

    // === 低级接口：分步求解（供高级用户/耦合问题使用） ===

    /**
     * @brief 预处理步骤：初始化矩阵和向量
     */
    void preprocess();

    /**
     * @brief 组装全局刚度矩阵和载荷向量
     */
    void assemble();

    /**
     * @brief 设置外部组装的刚度矩阵（用于耦合问题）
     * @param K_external 外部组装的刚度矩阵
     */
    void setStiffnessMatrix(const Eigen::SparseMatrix<double>& K_external);

    /**
     * @brief 设置外部组装的载荷向量（用于耦合问题）
     * @param b_external 外部组装的载荷向量
     */
    void setLoadVector(const Eigen::VectorXd& b_external);

    /**
     * @brief 施加边界条件
     */
    void applyBoundaryConditions();

    /**
     * @brief 求解线性方程组 K * u = b
     */
    void solve();

    /**
     * @brief 后处理：输出结果并计算误差
     */
    void postprocess();

    // === 访问器 ===

    /**
     * @brief 获取求解结果向量
     */
    const Eigen::VectorXd& getSolution() const;

    /**
     * @brief 获取节点数
     */
    int getNodesNum() const;

    /**
     * @brief 获取单元数
     */
    int getElementsNum() const;

    /**
     * @brief 获取配置对象
     */
    std::shared_ptr<Config> getConfig() const;

    /**
     * @brief 获取问题设置对象
     */
    std::shared_ptr<ProblemSetup> getProblem() const;

    // === 梯度计算（用于电热耦合） ===

    /**
     * @brief 计算指定点处的解梯度范数平方 |∇u|²
     *
     * 用于电热耦合中计算焦耳热密度 Q = σ|∇V|²
     *
     * @param element_idx 单元索引
     * @param ref_coords 参考坐标系中的点坐标
     * @return 梯度范数的平方 (∂u/∂x)² + (∂u/∂y)²
     */
    double computeGradientNormSquared(int element_idx, const std::vector<double>& ref_coords);
};

#endif  // FEM_SOLVER_H
