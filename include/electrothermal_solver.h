#ifndef ELECTROTHERMAL_SOLVER_H
#define ELECTROTHERMAL_SOLVER_H

#include <Eigen/Dense>
#include <memory>
#include <string>
#include "config.h"
#include "fem_solver.h"
#include "problem_definitions.h"

/**
 * @brief 电热耦合求解器类
 *
 * 求解双向耦合的电-热问题：
 * - 电场方程：-∇·(σ(T)∇V) = 0
 * - 热场方程：-∇·(k∇T) = Q(V,T)
 *
 * 其中：
 * - σ(T) = σ0[1 + α(T-T0)] : 温度依赖的电导率
 * - Q(V,T) = σ(T)|∇V|² : 焦耳热源
 *
 * 使用固定点迭代求解耦合问题。
 */
class ElectrothermalSolver
{
  private:
    // 两个场的FEM求解器
    std::shared_ptr<FEMSolver> electric_solver_;  // 电场求解器
    std::shared_ptr<FEMSolver> thermal_solver_;   // 热场求解器

    // 物理参数
    double sigma0_;  // 参考电导率 [S/m]
    double alpha_;   // 电导率温度系数 [1/K]
    double T0_;      // 参考温度 [K]

    // 迭代参数
    double tol_;                // 收敛容差（相对误差）
    int max_iter_;              // 最大迭代次数
    bool verbose_;              // 是否输出详细信息
    double relaxation_factor_;  // 松弛因子 (0,1]，用于稳定迭代

    // 当前场解
    Eigen::VectorXd T_current_;  // 当前温度增量场 ΔT [K]（绝对温度 = T0 + ΔT）
    Eigen::VectorXd V_current_;  // 当前电势场 [V]

    int iterations_;  // 实际迭代次数

    // 辅助函数

    /**
     * @brief 计算指定单元指定点的焦耳热密度
     * @param element_idx 单元索引
     * @param ref_coords 参考坐标
     * @return 焦耳热密度 Q = σ(T)|∇V|² [W/m³]
     */
    double computeJouleHeat(int element_idx, const std::vector<double>& ref_coords);

    /**
     * @brief 计算指定单元指定点的温度增量 ΔT（用于计算指定点电导率）
     * @param element_idx 单元索引
     * @param ref_coords 参考坐标
     * @return 温度增量 ΔT [K]（绝对温度 = T0 + ΔT）
     */
    double interpolateTemperature(int element_idx, const std::vector<double>& ref_coords);

    /**
     * @brief 计算温度依赖的电导率
     * @param delta_T 温度增量 ΔT [K]
     * @return 电导率 σ(T) [S/m]，其中 T = T0 + ΔT
     */
    double computeConductivity(double delta_T) const;

    /**
     * @brief 组装电场刚度矩阵（考虑温度依赖的电导率）
     *
     * 在每个积分点处插值温度，计算该点的电导率，精确组装刚度矩阵。
     * 这是方案1的简化实现：直接在ElectrothermalSolver中重写电场组装。
     *
     * @return 组装好的刚度矩阵
     */
    Eigen::SparseMatrix<double> assembleElectricFieldCoupled();

    /**
     * @brief 计算单元刚度矩阵的一个元素（电场，考虑温度耦合）
     * @param e 单元索引
     * @param alpha 试探函数索引
     * @param beta 检验函数索引
     * @return 刚度矩阵元素值
     */
    double calculateElectricStiffnessEntry(int e, int alpha, int beta);

    /**
     * @brief 组装热场载荷向量（考虑焦耳热源）
     *
     * 在每个积分点处计算焦耳热 Q = σ(T)|∇V|²，精确组装载荷向量。
     *
     * @return 组装好的载荷向量
     */
    Eigen::VectorXd assembleThermalLoadCoupled();

    /**
     * @brief 计算单元载荷向量的一个元素（热场，考虑焦耳热源）
     * @param e 单元索引
     * @param beta 检验函数索引
     * @return 载荷向量元素值
     */
    double calculateThermalLoadEntry(int e, int beta);

  public:
    /**
     * @brief 构造函数
     *
     * @param config 网格配置对象（电场和热场共享同一个网格）
     * @param et_problem 电热耦合问题配置（包含两个场的ProblemSetup和物理参数）
     * @param tol 收敛容差（默认1e-6）
     * @param max_iter 最大迭代次数（默认100）
     * @param verbose 是否输出详细信息（默认true）
     * @param relaxation 松弛因子（默认0.5），范围(0,1]，值越小迭代越稳定但可能更慢
     */
    ElectrothermalSolver(std::shared_ptr<Config> config,
                         const ProblemLibrary::ElectrothermalProblem& et_problem, double tol = 1e-6,
                         int max_iter = 100, bool verbose = true, double relaxation = 0.5);

    /**
     * @brief 执行电热耦合迭代求解
     *
     * 使用固定点迭代：
     * 1. 初始化 T^0 = T0
     * 2. 循环：
     *    a. 计算 σ^n = σ0[1 + α(T^(n-1) - T0)]
     *    b. 求解电场：-∇·(σ^n∇V) = 0 → V^n
     *    c. 计算焦耳热：Q^n = σ^n|∇V^n|²
     *    d. 求解热场：-∇·(k∇T) = Q^n → T^n
     *    e. 检查收敛：||T^n - T^(n-1)||/||T^n|| < tol
     */
    void solve();

    /**
     * @brief 输出结果到VTK文件
     *
     * @param prefix 文件前缀（默认"results/electrothermal"）
     *
     * 输出三个文件：
     * - {prefix}_V.vtu : 电势分布
     * - {prefix}_T.vtu : 温度分布
     * - {prefix}_Q.vtu : 焦耳热密度分布
     */
    void outputResults(const std::string& prefix = "results/electrothermal");

    // 访问器

    /**
     * @brief 获取电势场解
     */
    const Eigen::VectorXd& getElectricPotential() const;

    /**
     * @brief 获取温度增量场解 ΔT（绝对温度 = T0 + ΔT）
     */
    const Eigen::VectorXd& getTemperature() const;

    /**
     * @brief 获取迭代次数（solve()后有效）
     */
    int getIterations() const;
};

#endif  // ELECTROTHERMAL_SOLVER_H
