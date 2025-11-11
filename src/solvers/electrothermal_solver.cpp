#include "electrothermal_solver.h"
#include <iostream>
#include "gauss_quadrature.h"
#include "geometry_mapping.h"
#include "shape_functions.h"
#include "vtk_output.h"

// 构造函数
ElectrothermalSolver::ElectrothermalSolver(std::shared_ptr<Config> config,
                                           const ProblemLibrary::ElectrothermalProblem& et_problem,
                                           double tol, int max_iter, bool verbose,
                                           double relaxation)
    : sigma0_(et_problem.params.sigma0),
      alpha_(et_problem.params.alpha),
      T0_(et_problem.params.T0),
      tol_(tol),
      max_iter_(max_iter),
      verbose_(verbose),
      relaxation_factor_(relaxation),
      iterations_(0)
{
    // 创建电场和热场的FEM求解器（共享同一个网格配置）
    electric_solver_ = std::make_shared<FEMSolver>(config, et_problem.electric_problem);
    thermal_solver_ = std::make_shared<FEMSolver>(config, et_problem.thermal_problem);

    // 初始化温度增量场为0（绝对温度 = T0 + ΔT = T0）
    int N = config->getNodesNum();
    T_current_ = Eigen::VectorXd::Zero(N);  // ΔT = 0，对应绝对温度 T0
    V_current_ = Eigen::VectorXd::Zero(N);

    if (verbose_)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "    电热耦合求解器初始化完成" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "电场问题: " << et_problem.electric_problem->getName() << std::endl;
        std::cout << "热场问题: " << et_problem.thermal_problem->getName() << std::endl;
        std::cout << "\n物理参数：" << std::endl;
        std::cout << "  参考电导率 σ0 = " << sigma0_ << " S/m" << std::endl;
        std::cout << "  温度系数   α  = " << alpha_ << " /K" << std::endl;
        std::cout << "  参考温度   T0 = " << T0_ << " K" << std::endl;
        std::cout << "  热导率     k  = 1.0 W/(m·K) (来自热场问题定义)" << std::endl;
        std::cout << "\n迭代参数：" << std::endl;
        std::cout << "  收敛容差       = " << tol_ << std::endl;
        std::cout << "  最大迭代次数   = " << max_iter_ << std::endl;
        std::cout << "  松弛因子       = " << relaxation_factor_ << " (用于稳定收敛)" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }
}

// 计算温度依赖的电导率
// 注意：输入参数是温度增量 ΔT，绝对温度 T_abs = T0 + ΔT
// σ(T) = σ0[1 + α(T - T0)] = σ0[1 + α·ΔT]
double ElectrothermalSolver::computeConductivity(double delta_T) const
{
    return sigma0_ * (1.0 + alpha_ * delta_T);
}

// 插值计算指定点的温度增量 ΔT
double ElectrothermalSolver::interpolateTemperature(int element_idx,
                                                    const std::vector<double>& ref_coords)
{
    auto config = thermal_solver_->getConfig();
    const auto& connectivity = config->getElementConnectivity();
    int n = connectivity[0].size();  // 每个单元的节点数

    // 创建形函数对象
    auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config);

    // 插值：ΔT(ξ,η) = Σ ΔT_i * N_i(ξ,η)
    double delta_T = 0.0;
    for (int i = 0; i < n; ++i)
    {
        int node_idx = connectivity[element_idx][i];
        double N_i = shapeFunction->computeTrialFunction(i, ref_coords);
        delta_T += T_current_(node_idx) * N_i;
    }

    return delta_T;
}

// 计算焦耳热密度
double ElectrothermalSolver::computeJouleHeat(int element_idx,
                                              const std::vector<double>& ref_coords)
{
    // Q = σ(T) * |∇V|²，其中 T = T0 + ΔT

    // 1. 插值得到该点的温度增量 ΔT
    double delta_T = interpolateTemperature(element_idx, ref_coords);

    // 2. 计算温度依赖的电导率 σ(T) = σ0[1 + α·ΔT]
    double sigma = computeConductivity(delta_T);

    // 3. 计算电势梯度的范数平方 |∇V|²
    double grad_V_squared = electric_solver_->computeGradientNormSquared(element_idx, ref_coords);

    // 4. 计算焦耳热密度
    return sigma * grad_V_squared;
}

// 计算单元刚度矩阵元素（电场，考虑温度耦合）
double ElectrothermalSolver::calculateElectricStiffnessEntry(int e, int alpha, int beta)
{
    auto config = electric_solver_->getConfig();
    const auto& connectivity = config->getElementConnectivity();
    const auto& coordinates = config->getNodeCoordinates();
    int dimension = config->getDimension();
    int n = connectivity[e].size();  // 每个单元的节点数

    // 提取单元节点坐标
    std::vector<double> element_coords(n * dimension);
    for (int i = 0; i < n; ++i)
    {
        int node_idx = connectivity[e][i];
        for (int d = 0; d < dimension; ++d)
        {
            element_coords[i * dimension + d] = coordinates[node_idx * dimension + d];
        }
    }

    // 创建几何映射对象
    auto mapping = GeometryMappingFactory::createMapping(element_coords, config);

    // 创建高斯积分对象
    auto gaussPoint = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config->getElementType()),
        config->getAssembleGaussPointsNum());

    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    double entryValue = 0.0;

    // 在每个积分点上计算
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        // 参考坐标
        std::vector<double> pointsRef(dimension);
        for (int d = 0; d < dimension; ++d)
        {
            pointsRef[d] = points[gpIndex * dimension + d];
        }

        // ✅ 关键：在积分点处插值温度增量 ΔT
        double delta_T_gp = interpolateTemperature(e, pointsRef);

        // ✅ 计算该点的电导率 σ(T) = σ0[1 + α·ΔT]
        double sigma_gp = computeConductivity(delta_T_gp);

        // 计算形函数梯度
        auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config);
        std::vector<double> trialGradientsRef =
            shapeFunction->computeTrialGradients(alpha, pointsRef);
        std::vector<double> testGradientsRef = shapeFunction->computeTestGradients(beta, pointsRef);

        // 转换为物理坐标系下的梯度
        std::vector<double> trialGradientsPhys;
        std::vector<double> testGradientsPhys;
        mapping->transformGradient(trialGradientsRef, trialGradientsPhys, pointsRef);
        mapping->transformGradient(testGradientsRef, testGradientsPhys, pointsRef);

        // 计算积分被积函数：σ(T) * ∇φ_alpha · ∇φ_beta
        double integrand = 0.0;
        for (int d = 0; d < dimension; ++d)
        {
            integrand += trialGradientsPhys[d] * testGradientsPhys[d];
        }
        integrand *= sigma_gp;

        // 累加积分
        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }

    return entryValue;
}

// 组装电场刚度矩阵（考虑温度耦合）
Eigen::SparseMatrix<double> ElectrothermalSolver::assembleElectricFieldCoupled()
{
    auto config = electric_solver_->getConfig();
    int N = config->getNodesNum();
    int M = config->getElementsNum();
    const auto& connectivity = config->getElementConnectivity();
    int n = connectivity[0].size();

    // 使用三元组列表收集所有元素
    typedef Eigen::Triplet<double> T_entry;
    std::vector<T_entry> tripletList;
    tripletList.reserve(n * n * M);

    // 遍历所有单元
    for (int e = 0; e < M; ++e)
    {
        for (int alpha = 0; alpha < n; ++alpha)
        {
            for (int beta = 0; beta < n; ++beta)
            {
                // ✅ 使用新的计算函数，在积分点处考虑温度
                double K_e_val = calculateElectricStiffnessEntry(e, alpha, beta);
                int global_row = connectivity[e][beta];
                int global_col = connectivity[e][alpha];
                tripletList.push_back(T_entry(global_row, global_col, K_e_val));
            }
        }
    }

    // ✅ 关键修复：直接设置 electric_solver_ 的三元组列表
    // 这样 applyBoundaryConditions() 就可以正确处理边界条件
    electric_solver_->setTripletList(tripletList);

    // 仍然构建并返回稀疏矩阵（用于调试/验证）
    Eigen::SparseMatrix<double> K_global(N, N);
    K_global.setFromTriplets(tripletList.begin(), tripletList.end());
    K_global.makeCompressed();

    return K_global;
}

// 计算单元载荷向量元素（热场，考虑焦耳热源）
double ElectrothermalSolver::calculateThermalLoadEntry(int e, int beta)
{
    auto config = thermal_solver_->getConfig();
    const auto& connectivity = config->getElementConnectivity();
    const auto& coordinates = config->getNodeCoordinates();
    int dimension = config->getDimension();
    int n = connectivity[e].size();

    // 提取单元节点坐标
    std::vector<double> element_coords(n * dimension);
    for (int i = 0; i < n; ++i)
    {
        int node_idx = connectivity[e][i];
        for (int d = 0; d < dimension; ++d)
        {
            element_coords[i * dimension + d] = coordinates[node_idx * dimension + d];
        }
    }

    // 创建几何映射对象
    auto mapping = GeometryMappingFactory::createMapping(element_coords, config);

    // 创建高斯积分对象
    auto gaussPoint = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config->getElementType()),
        config->getAssembleGaussPointsNum());

    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    double entryValue = 0.0;

    // 在每个积分点上计算
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        // 参考坐标
        std::vector<double> pointsRef(dimension);
        for (int d = 0; d < dimension; ++d)
        {
            pointsRef[d] = points[gpIndex * dimension + d];
        }

        // ✅ 关键：在积分点计算焦耳热密度
        double Q_gp = computeJouleHeat(e, pointsRef);

        // 计算形函数值
        auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config);
        double N_beta = shapeFunction->computeTestFunction(beta, pointsRef);

        // 计算积分被积函数：Q * N_beta
        double integrand = Q_gp * N_beta;

        // 累加积分
        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }

    return entryValue;
}

// 组装热场载荷向量（考虑焦耳热源）
Eigen::VectorXd ElectrothermalSolver::assembleThermalLoadCoupled()
{
    auto config = thermal_solver_->getConfig();
    int N = config->getNodesNum();
    int M = config->getElementsNum();
    const auto& connectivity = config->getElementConnectivity();
    int n = connectivity[0].size();

    // 初始化载荷向量
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);

    // 遍历所有单元
    for (int e = 0; e < M; ++e)
    {
        for (int beta = 0; beta < n; ++beta)
        {
            // ✅ 使用新的计算函数，在积分点处精确计算焦耳热
            double b_e_val = calculateThermalLoadEntry(e, beta);
            int global_row = connectivity[e][beta];
            b(global_row) += b_e_val;
        }
    }

    return b;
}

// 主求解函数
void ElectrothermalSolver::solve()
{
    if (verbose_)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "      开始电热耦合迭代求解" << std::endl;
        std::cout << "      (使用精确的积分点耦合)" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    // 预处理两个求解器(边界条件在此阶段设置)
    electric_solver_->preprocess();
    thermal_solver_->preprocess();

    // 固定点迭代
    double final_error = 0.0;  // 记录最终误差
    for (iterations_ = 1; iterations_ <= max_iter_; ++iterations_)
    {
        if (verbose_)
        {
            std::cout << "\n--- 第 " << iterations_ << " 次迭代 ---" << std::endl;
        }

        // ============================================
        // 步骤 1: 求解电场（使用上一步的温度场）
        // ============================================

        if (verbose_)
        {
            std::cout << "\n[电场求解] (积分点精确耦合)" << std::endl;
        }

        // ✅ 使用精确的积分点耦合组装刚度矩阵
        if (verbose_)
        {
            std::cout << "  在每个积分点插值温度，计算σ(T)..." << std::endl;
        }

        // ✅ 组装刚度矩阵并设置三元组列表
        // 注意：assembleElectricFieldCoupled() 内部会调用 electric_solver_->setTripletList()
        Eigen::SparseMatrix<double> K_electric_coupled = assembleElectricFieldCoupled();

        if (verbose_)
        {
            std::cout << "  组装完成，非零元素数: " << K_electric_coupled.nonZeros() << std::endl;
        }

        // ✅ 重置载荷向量为0（电场方程无源项）
        Eigen::VectorXd b_electric =
            Eigen::VectorXd::Zero(electric_solver_->getConfig()->getNodesNum());
        electric_solver_->setLoadVector(b_electric);

        // ✅ 施加边界条件（现在可以从三元组正确处理）
        electric_solver_->applyBoundaryConditions();
        electric_solver_->solve();

        // 获取电场解
        V_current_ = electric_solver_->getSolution();

        if (verbose_)
        {
            std::cout << "电势范围: [" << V_current_.minCoeff() << ", " << V_current_.maxCoeff()
                      << "] V" << std::endl;
        }

        // ============================================
        // 步骤 2: 求解热场（使用刚求得的电场）
        // ============================================

        if (verbose_)
        {
            std::cout << "\n[热场求解] (积分点精确焦耳热)" << std::endl;
        }

        // 保存上一步的温度场用于收敛判断
        Eigen::VectorXd T_old = T_current_;

        // ✅ 使用精确的积分点焦耳热组装载荷向量
        // 注意：热导率由 thermal_problem 定义（默认为 1.0）
        if (verbose_)
        {
            std::cout << "  在每个积分点计算焦耳热 Q = σ(T)|∇V|²..." << std::endl;
        }

        Eigen::VectorXd b_thermal_coupled = assembleThermalLoadCoupled();

        if (verbose_)
        {
            std::cout << "  载荷向量组装完成" << std::endl;
            double Q_total = b_thermal_coupled.sum();
            std::cout << "  总焦耳热功率: " << Q_total << " W" << std::endl;
        }

        // 组装热场刚度矩阵（需要热导率 k）
        // 我们需要调用 assemble() 来组装刚度矩阵
        // 但需要确保 thermal_solver_ 的 problem_ 有正确的系数和源项
        thermal_solver_->assemble();

        // ✅ 使用我们精确组装的载荷向量
        // 注意：thermal_solver_->assemble() 已经组装了b（源项为0），
        //      所以直接设置为我们的焦耳热载荷向量
        thermal_solver_->setLoadVector(b_thermal_coupled);

        // 施加边界条件和求解
        thermal_solver_->applyBoundaryConditions();
        thermal_solver_->solve();

        // 获取温度场解
        Eigen::VectorXd T_new = thermal_solver_->getSolution();

        // ⚠️ 应用欠松弛以提高稳定性
        // T_current = ω * T_new + (1-ω) * T_old
        // 松弛因子ω越小，迭代越稳定但可能更慢
        if (iterations_ > 1)
        {
            T_current_ = relaxation_factor_ * T_new + (1.0 - relaxation_factor_) * T_old;
        }
        else
        {
            T_current_ = T_new;  // 第一次迭代不应用松弛
        }

        if (verbose_)
        {
            std::cout << "温度增量 ΔT 范围: [" << T_current_.minCoeff() << ", "
                      << T_current_.maxCoeff() << "] K" << std::endl;
            std::cout << "绝对温度范围: [" << (T0_ + T_current_.minCoeff()) << ", "
                      << (T0_ + T_current_.maxCoeff()) << "] K" << std::endl;
            if (iterations_ > 1)
            {
                std::cout << "（已应用松弛因子 " << relaxation_factor_ << "）" << std::endl;
            }
        }

        // ============================================
        // 步骤 3: 检查收敛
        // ============================================

        double T_norm = T_current_.norm();
        double T_diff_norm = (T_current_ - T_old).norm();
        double relative_error = (T_norm > 1e-10) ? (T_diff_norm / T_norm) : T_diff_norm;
        final_error = relative_error;  // 记录每次迭代的误差

        if (verbose_)
        {
            std::cout << "\n[收敛检查]" << std::endl;
            std::cout << "温度变化相对误差: " << relative_error << std::endl;
        }

        if (relative_error < tol_)
        {
            if (verbose_)
            {
                std::cout << "\n========================================" << std::endl;
                std::cout << "    收敛！迭代次数: " << iterations_ << std::endl;
                std::cout << "========================================\n" << std::endl;
            }
            return;
        }
    }

    // 未收敛
    std::cout << "\n警告：达到最大迭代次数 (" << max_iter_ << ")，未完全收敛。" << std::endl;
    std::cout << "最终相对误差: " << final_error << " (容差要求: " << tol_ << ")" << std::endl;
    std::cout << "提示：可能需要调整容差、增加迭代次数，或检查边界条件设置" << std::endl;
}

// 输出结果
void ElectrothermalSolver::outputResults(const std::string& prefix)
{
    if (verbose_)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "         输出结果到VTK文件" << std::endl;
        std::cout << "========================================" << std::endl;
    }

    // 输出电势场
    auto vtkElectric = VTKOutputFactory::createVTKOutput(electric_solver_->getConfig(), V_current_);
    vtkElectric->outputNumericalSolution(prefix + "_V");
    std::cout << "已输出: " << prefix << "_V.vtu (电势场)" << std::endl;

    // 输出温度场（转换为绝对温度）
    Eigen::VectorXd T_absolute = T_current_.array() + T0_;  // T = T0 + ΔT
    auto vtkThermal = VTKOutputFactory::createVTKOutput(thermal_solver_->getConfig(), T_absolute);
    vtkThermal->outputNumericalSolution(prefix + "_T");
    std::cout << "已输出: " << prefix << "_T.vtu (温度场 - 绝对温度)" << std::endl;

    // 计算并输出焦耳热密度场
    // 在每个节点计算焦耳热密度（使用单元平均）
    int N = T_current_.size();
    Eigen::VectorXd Q_nodes = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd node_count = Eigen::VectorXd::Zero(N);

    auto config = electric_solver_->getConfig();
    const auto& connectivity = config->getElementConnectivity();
    int M = config->getElementsNum();
    int dimension = config->getDimension();

    for (int e = 0; e < M; ++e)
    {
        // 在单元中心计算焦耳热
        std::vector<double> center_ref(dimension, 1.0 / 3.0);  // 三角形中心
        double Q_center = computeJouleHeat(e, center_ref);

        // 分配到单元节点
        int n = connectivity[e].size();
        for (int i = 0; i < n; ++i)
        {
            int node_idx = connectivity[e][i];
            Q_nodes(node_idx) += Q_center;
            node_count(node_idx) += 1.0;
        }
    }

    // 节点平均
    for (int i = 0; i < N; ++i)
    {
        if (node_count(i) > 0)
        {
            Q_nodes(i) /= node_count(i);
        }
    }

    auto vtkHeat = VTKOutputFactory::createVTKOutput(config, Q_nodes);
    vtkHeat->outputNumericalSolution(prefix + "_Q");
    std::cout << "已输出: " << prefix << "_Q.vtu (焦耳热密度场)" << std::endl;

    std::cout << "\n焦耳热密度范围: [" << Q_nodes.minCoeff() << ", " << Q_nodes.maxCoeff()
              << "] W/m³" << std::endl;

    if (verbose_)
    {
        std::cout << "========================================\n" << std::endl;
    }
}

// 访问器实现
const Eigen::VectorXd& ElectrothermalSolver::getElectricPotential() const
{
    return V_current_;
}

const Eigen::VectorXd& ElectrothermalSolver::getTemperature() const
{
    return T_current_;
}

int ElectrothermalSolver::getIterations() const
{
    return iterations_;
}
