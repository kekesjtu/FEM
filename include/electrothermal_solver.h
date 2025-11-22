#ifndef ELECTROTHERMAL_SOLVER_H
#define ELECTROTHERMAL_SOLVER_H

#include <Eigen/Dense>
#include <memory>
#include <string>
#include "config.h"
#include "fem_solver.h"
#include "material.h"
#include "problem_setup.h"

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
    // 配置和材料
    std::shared_ptr<Config> config_;
    std::shared_ptr<ProblemSetup> problem_;
    std::shared_ptr<MaterialLibrary> materials_;

    // 两个场的FEM求解器
    std::shared_ptr<FEMSolver> electric_solver_;  // 电场求解器
    std::shared_ptr<FEMSolver> thermal_solver_;   // 热场求解器

    // 迭代参数(从problem_的coupling配置读取)
    double tol_;                // 收敛容差（相对误差）
    int max_iter_;              // 最大迭代次数
    double relaxation_factor_;  // 松弛因子 (0,1]，用于稳定迭代
    double T0_;                 // 参考温度 [K]

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

    /**
     * @brief 计算节点焦耳热密度场（使用单元平均）
     *
     * 在每个单元中心计算焦耳热密度，然后分配到单元节点并平均。
     *
     * @return 节点焦耳热密度向量 [W/m³]
     */
    Eigen::VectorXd computeJouleHeatField();

  public:
    /**
     * @brief 构造函数 (V2 API - JSON驱动)
     *
     * @param config 网格配置对象（电场和热场共享同一个网格）
     * @param problem 问题配置对象(包含electric和thermal两个场)
     * @param materials 材料库对象
     * @param verbose 是否输出详细信息（默认true）
     */
    ElectrothermalSolver(std::shared_ptr<Config> config, std::shared_ptr<ProblemSetup> problem,
                         std::shared_ptr<MaterialLibrary> materials);

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
};

#endif  // ELECTROTHERMAL_SOLVER_H
