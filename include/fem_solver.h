#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include <Eigen/Sparse>
#include <memory>
#include <string>
#include <vector>
#include "config.h"
#include "material.h"
#include "problem_setup.h"

/**
 * @brief 有限元求解器类 (V2版本 - JSON驱动)
 *
 * 封装了完整的有限元求解流程,支持通过JSON配置文件定义材料属性、边界条件和源项。
 * 适用于求解形式为 -∇·(c∇u) = f 的偏微分方程。
 */
class FEMSolver
{
  private:
    // 配置对象
    std::shared_ptr<Config> config_;
    std::shared_ptr<ProblemSetup> problem_;       // 问题配置
    std::shared_ptr<MaterialLibrary> materials_;  // 材料库
    std::string field_name_;                      // 当前求解的场名称 ("electric", "thermal"等)

    // 网格信息
    int N_;  // 节点数
    int M_;  // 单元数
    int n_;  // 每个单元的节点数

    // 有限元矩阵和向量
    Eigen::SparseMatrix<double> K_global_;  // 全局刚度矩阵
    Eigen::VectorXd b_;                     // 载荷向量
    Eigen::VectorXd u_;                     // 解向量

    // 三元组列表（用于延迟矩阵构建，优化性能）
    typedef Eigen::Triplet<double> T_entry;
    std::vector<T_entry> triplet_list_;

    // 内部辅助函数
    double calculateStiffnessEntry(int e, int alpha, int beta);
    double calculateLoadEntry(int e, int beta);
    double calculateBoundaryStiffness(const Config::Boundary& boundary, int local_i, int local_j,
                                      size_t boundary_idx);
    double calculateBoundaryLoad(const Config::Boundary& boundary, int local_i,
                                 size_t boundary_idx);

  public:
    /**
     * @brief 构造函数 (V2 API)
     * @param config 配置对象，包含网格信息
     * @param problem 问题配置(来自JSON)
     * @param materials 材料库(来自JSON)
     * @param field_name 要求解的物理场名称 (如 "electric", "thermal")
     */
    FEMSolver(std::shared_ptr<Config> config, std::shared_ptr<ProblemSetup> problem,
              std::shared_ptr<MaterialLibrary> materials, const std::string& field_name);

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
     *
     * 此函数依次调用 assembleStiffnessMatrix() 和 assembleLoadVector()
     */
    void assemble();

    /**
     * @brief 仅组装全局刚度矩阵（不处理载荷向量）
     *
     * 用于耦合问题中需要单独组装刚度矩阵的情况
     */
    void assembleStiffnessMatrix();

    /**
     * @brief 仅组装载荷向量（不处理刚度矩阵）
     *
     * 用于耦合问题中需要单独组装载荷向量的情况
     */
    void assembleLoadVector();

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
     * @brief 设置三元组列表（用于耦合问题的自定义组装）
     * @param triplets 外部组装的三元组列表
     *
     * 此方法允许耦合求解器直接设置三元组列表，
     * 使得 applyBoundaryConditions() 可以正确处理边界条件
     */
    void setTripletList(const std::vector<T_entry>& triplets);

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
