#include "electrothermal_solver.h"
#include <iostream>
#include "gauss_quadrature.h"
#include "geometry_mapping.h"
#include "logger.h"
#include "material.h"
#include "shape_functions.h"
#include "vtk_output.h"

// 构造函数
ElectrothermalSolver::ElectrothermalSolver(std::shared_ptr<Config> config,
                                           std::shared_ptr<ProblemSetup> problem,
                                           std::shared_ptr<MaterialLibrary> materials)
    : config_(config), problem_(problem), materials_(materials), iterations_(0)
{
    // 验证参数
    if (!config || !problem || !materials)
    {
        throw std::runtime_error("ElectrothermalSolver: 配置参数不能为空");
    }

    // 验证必须有electric和thermal两个场
    if (!problem->hasField("electric") || !problem->hasField("thermal"))
    {
        throw std::runtime_error("ElectrothermalSolver: 问题必须定义electric和thermal两个场");
    }

    // 验证必须有耦合配置
    if (!problem->has_coupling)
    {
        throw std::runtime_error("ElectrothermalSolver: 问题必须定义coupling配置");
    }

    // 从耦合配置读取参数
    const auto& coupling = problem->coupling;
    tol_ = coupling.tolerance;
    max_iter_ = coupling.max_iterations;
    relaxation_factor_ = coupling.relaxation_factor;
    T0_ = coupling.reference_temperature;

    // 创建电场和热场的FEM求解器（共享同一个网格配置）
    electric_solver_ = std::make_shared<FEMSolver>(config, problem, materials, "electric");
    thermal_solver_ = std::make_shared<FEMSolver>(config, problem, materials, "thermal");

    // 初始化温度增量场为0（绝对温度 = T0 + ΔT = T0）
    int N = config->getNodesNum();
    T_current_ = Eigen::VectorXd::Zero(N);  // ΔT = 0，对应绝对温度 T0
    V_current_ = Eigen::VectorXd::Zero(N);

    LOG_HEADER("电热耦合求解器初始化完成");
    LOG_INFO("问题名称: " + problem->name);
    LOG_INFO("问题描述: " + problem->description);
    LOG_INFO("");
    LOG_INFO("物理场：");
    LOG_INFO("  - Electric Field (电场)");
    LOG_INFO("  - Thermal Field (热场)");
    LOG_INFO("");
    LOG_INFO("耦合参数：");
    LOG_INFO("  耦合类型       = " + coupling.type);
    LOG_INFO("  参考温度   T0  = " + std::to_string(T0_) + " K");
    LOG_INFO("  收敛容差       = " + std::to_string(tol_));
    LOG_INFO("  最大迭代次数   = " + std::to_string(max_iter_));
    LOG_INFO("  松弛因子       = " + std::to_string(relaxation_factor_) + " (用于稳定收敛)");
}

// 计算温度依赖的电导率
// 注意：输入参数是温度增量 ΔT，绝对温度 T_abs = T0 + ΔT
// σ(T) 从材料属性中评估
double ElectrothermalSolver::computeConductivity(double delta_T) const
{
    // 获取第一个域的材料(假设整个域使用相同材料)
    const auto& element_entities = config_->getElementGeometricEntities();
    int entity_id = element_entities[0];
    std::string material_name = problem_->getMaterialForDomain(entity_id);
    auto material = materials_->getMaterial(material_name);

    // 构造评估上下文,设置温度
    EvaluationContext ctx;
    ctx.T = T0_ + delta_T;  // 绝对温度

    // 从材料属性评估电导率
    return material->electrical_conductivity.evaluate(ctx);
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

// 计算节点焦耳热密度场
Eigen::VectorXd ElectrothermalSolver::computeJouleHeatField()
{
    int N = T_current_.size();
    Eigen::VectorXd Q_nodes = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd node_count = Eigen::VectorXd::Zero(N);

    auto config = electric_solver_->getConfig();
    const auto& connectivity = config->getElementConnectivity();
    int M = config->getElementsNum();
    int dimension = config->getDimension();

    // 在每个单元中心计算焦耳热
    for (int e = 0; e < M; ++e)
    {
        // 单元中心参考坐标
        std::vector<double> center_ref;
        if (dimension == 2)
        {
            center_ref = {1.0 / 3.0, 1.0 / 3.0};  // 三角形中心
        }
        else if (dimension == 3)
        {
            center_ref = {0.25, 0.25, 0.25};  // 四面体中心
        }
        else
        {
            center_ref = {0.5};  // 1D线段中心
        }

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

    return Q_nodes;
}

// 主求解函数
void ElectrothermalSolver::solve()
{
    TIMER_SCOPE("电热耦合总求解");

    LOG_HEADER("开始电热耦合迭代求解");

    // 预处理两个求解器(边界条件在此阶段设置)
    {
        TIMER_SCOPE_DEBUG("求解器预处理");
        electric_solver_->preprocess();
        thermal_solver_->preprocess();
    }

    // 固定点迭代
    double final_error = 0.0;  // 记录最终误差
    for (iterations_ = 1; iterations_ <= max_iter_; ++iterations_)
    {
        TIMER_SCOPE_DEBUG("第 " + std::to_string(iterations_) + " 次迭代");

        LOG_INFO("");
        LOG_INFO("--- 第 " + std::to_string(iterations_) + " 次迭代 ---");

        // ============================================
        // 步骤 1: 求解电场（使用刚求得的温度场）
        // ============================================

        Eigen::SparseMatrix<double> K_electric_coupled;  // 在外部定义，供后续日志使用

        {
            TIMER_SCOPE("  [电场] 完整求解");

            LOG_INFO("[电场求解]");

            // 组装刚度矩阵并设置三元组列表
            // 注意：assembleElectricFieldCoupled() 内部会调用 electric_solver_->setTripletList()
            {
                TIMER_SCOPE_DEBUG("    - 组装温度相关刚度矩阵");
                K_electric_coupled = assembleElectricFieldCoupled();
                LOG_DEBUG("  组装完成，非零元素数: " +
                          std::to_string(K_electric_coupled.nonZeros()));
            }

            // 组装载荷向量
            {
                TIMER_SCOPE_DEBUG("    - 组装载荷向量");
                electric_solver_->assembleLoadVector();
            }

            // 应用边界条件
            {
                TIMER_SCOPE_DEBUG("    - 应用边界条件");
                electric_solver_->applyBoundaryConditions();
            }

            // 求解线性系统
            {
                TIMER_SCOPE_DEBUG("    - CG求解");
                electric_solver_->solve();
            }

            // 更新电势场
            V_current_ = electric_solver_->getSolution();

            LOG_INFO("电势范围: [" + std::to_string(V_current_.minCoeff()) + ", " +
                     std::to_string(V_current_.maxCoeff()) + "] V");
        }

        // ============================================
        // 步骤 2: 求解热场（使用刚求得的电场）
        // ============================================

        {
            TIMER_SCOPE("  [热场] 完整求解");

            LOG_INFO("[热场求解]");

            // 保存上一步的温度场用于收敛判断
            Eigen::VectorXd T_old = T_current_;

            // 组装热场刚度矩阵（需要热导率 k）
            {
                TIMER_SCOPE_DEBUG("    - 组装热导率矩阵");
                thermal_solver_->assembleStiffnessMatrix();
            }

            // 组装热场载荷向量（包含焦耳热源）
            {
                TIMER_SCOPE_DEBUG("    - 组装焦耳热载荷");
                Eigen::VectorXd b_thermal_coupled = assembleThermalLoadCoupled();
                thermal_solver_->setLoadVector(b_thermal_coupled);
            }

            // 施加边界条件和求解
            {
                TIMER_SCOPE_DEBUG("    - 应用边界条件");
                thermal_solver_->applyBoundaryConditions();
            }

            {
                TIMER_SCOPE_DEBUG("    - CG求解");
                thermal_solver_->solve();
            }

            // 获取温度场解
            Eigen::VectorXd T_new = thermal_solver_->getSolution();

            // 应用欠松弛以提高稳定性
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

            LOG_DEBUG("（已应用松弛因子 " + std::to_string(relaxation_factor_) + "）");
            LOG_INFO("温度增量 ΔT 范围: [" + std::to_string(T_current_.minCoeff()) + ", " +
                     std::to_string(T_current_.maxCoeff()) + "] K");
            LOG_INFO("绝对温度范围: [" + std::to_string(T0_ + T_current_.minCoeff()) + ", " +
                     std::to_string(T0_ + T_current_.maxCoeff()) + "] K");

            // ============================================
            // 步骤 3: 检查收敛
            // ============================================

            double T_norm = T_current_.norm();
            double T_diff_norm = (T_current_ - T_old).norm();
            double relative_error = (T_norm > 1e-10) ? (T_diff_norm / T_norm) : T_diff_norm;
            final_error = relative_error;  // 记录每次迭代的误差

            LOG_INFO("[收敛检查] 温度变化相对误差: " + std::to_string(relative_error));

            if (relative_error < tol_)
            {
                LOG_INFO("");
                LOG_INFO("收敛！迭代次数: " + std::to_string(iterations_) +
                         ", 最终误差: " + std::to_string(relative_error));
                return;
            }
        }
    }

    // 循环结束后仍未收敛
    LOG_INFO("");
    LOG_WARNING("达到最大迭代次数 (" + std::to_string(max_iter_) + ")，未完全收敛。");
    LOG_WARNING("最终相对误差: " + std::to_string(final_error) +
                " (容差要求: " + std::to_string(tol_) + ")");
    LOG_WARNING("警告：可能需要调整容差、增加迭代次数，或检查边界条件设置");
}

// 输出结果
void ElectrothermalSolver::outputResults(const std::string& prefix)
{
    LOG_HEADER("输出结果到VTK文件");

    // 输出电势场
    auto vtkElectric = VTKOutputFactory::createVTKOutput(electric_solver_->getConfig(), V_current_);
    vtkElectric->outputNumericalSolution(prefix + "_V");
    LOG_INFO("已输出: " + prefix + "_V.vtu (电势场)");

    // 输出温度场（转换为绝对温度）
    Eigen::VectorXd T_absolute = T_current_.array() + T0_;  // T = T0 + ΔT
    auto vtkThermal = VTKOutputFactory::createVTKOutput(thermal_solver_->getConfig(), T_absolute);
    vtkThermal->outputNumericalSolution(prefix + "_T");
    LOG_INFO("已输出: " + prefix + "_T.vtu (温度场 - 绝对温度)");

    // 计算并输出焦耳热密度场
    Eigen::VectorXd Q_nodes = computeJouleHeatField();

    auto config = electric_solver_->getConfig();
    auto vtkHeat = VTKOutputFactory::createVTKOutput(config, Q_nodes);
    vtkHeat->outputNumericalSolution(prefix + "_Q");
    LOG_INFO("已输出: " + prefix + "_Q.vtu (焦耳热密度场)");

    LOG_INFO("");
    LOG_INFO("焦耳热密度范围: [" + std::to_string(Q_nodes.minCoeff()) + ", " +
             std::to_string(Q_nodes.maxCoeff()) + "] W/m³");
}