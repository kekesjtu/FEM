#include "fem_solver.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <ctime>
#include <iostream>
#include "error_analysis.h"
#include "gauss_quadrature.h"
#include "geometry_mapping.h"
#include "logger.h"
#include "material.h"
#include "shape_functions.h"
#include "vtk_output.h"

// --- FEMSolver 类实现 (V2版本 - JSON驱动) ---

// 构造函数
FEMSolver::FEMSolver(std::shared_ptr<Config> config, std::shared_ptr<ProblemSetup> problem,
                     std::shared_ptr<MaterialLibrary> materials, const std::string& field_name)
    : config_(config),
      problem_(problem),
      materials_(materials),
      field_name_(field_name),
      N_(0),
      M_(0),
      n_(0)
{
    // 验证参数
    if (!config)
    {
        throw std::runtime_error("Config对象不能为空");
    }
    if (!problem)
    {
        throw std::runtime_error("ProblemSetup对象不能为空");
    }
    if (!materials)
    {
        throw std::runtime_error("MaterialLibrary对象不能为空");
    }
    if (!problem->hasField(field_name))
    {
        throw std::runtime_error("问题配置中不存在场: " + field_name);
    }
}

// 高级接口：完整求解
void FEMSolver::solveComplete()
{
    TIMER_SCOPE("FEM完整求解 (" + field_name_ + ")");

    LOG_SUBHEADER("求解物理场: " + field_name_);

    {
        TIMER_SCOPE("1.网格生成和预处理完成");
        preprocess();
    }

    {
        TIMER_SCOPE("2.全局矩阵组装完成");
        assemble();
    }

    {
        TIMER_SCOPE("3.边界条件施加完成");
        applyBoundaryConditions();
    }

    {
        TIMER_SCOPE("4: 线性方程组求解完成");
        solve();
    }

    {
        TIMER_SCOPE_DEBUG("5.后处理完成");
        postprocess();
    }
}

// 预处理
void FEMSolver::preprocess()
{
    // 从Config获取网格信息
    N_ = config_->getNodesNum();
    M_ = config_->getElementsNum();
    n_ = config_->getNodesNumPerElement();

    // 初始化矩阵和向量
    K_global_.resize(N_, N_);
    // 动态计算稀疏矩阵预留空间每个节点大约连接n²个其他节点
    K_global_.reserve(Eigen::VectorXi::Constant(N_, n_ * n_));
    b_ = Eigen::VectorXd::Zero(N_);
    u_ = Eigen::VectorXd::Zero(N_);

    LOG_DEBUG("初始化场 [" + field_name_ + "]: " + std::to_string(N_) + " 个节点, " +
              std::to_string(M_) + " 个单元");
}

void FEMSolver::assemble()
{
    assembleStiffnessMatrix();
    assembleLoadVector();
}

void FEMSolver::assembleStiffnessMatrix()
{
    // 清空并预留空间
    triplet_list_.clear();
    triplet_list_.reserve(n_ * n_ * M_);  // 每个单元贡献n²个矩阵元素

    const auto& connectivity = config_->getElementConnectivity();

    // 只收集三元组，不立即构建矩阵
    for (int e = 0; e < M_; ++e)
    {
        for (int alpha = 0; alpha < n_; ++alpha)
        {
            for (int beta = 0; beta < n_; ++beta)
            {
                double K_e_val = calculateStiffnessEntry(e, alpha, beta);
                int global_row = connectivity[e][beta];
                int global_col = connectivity[e][alpha];
                triplet_list_.push_back(T_entry(global_row, global_col, K_e_val));
            }
        }
    }

    LOG_DEBUG("组装刚度矩阵：收集了 " + std::to_string(triplet_list_.size()) + " 个矩阵元素");
    // 注意：矩阵构建延迟到 applyBoundaryConditions() 中进行
}

void FEMSolver::assembleLoadVector()
{
    // 确保 b_ 已经初始化（在 preprocess 中完成）
    if (b_.size() != N_)
    {
        b_ = Eigen::VectorXd::Zero(N_);
    }
    else
    {
        b_.setZero();  // 重置为0
    }

    const auto& connectivity = config_->getElementConnectivity();

    // 组装载荷向量
    for (int e = 0; e < M_; ++e)
    {
        for (int beta = 0; beta < n_; ++beta)
        {
            double b_e_val = calculateLoadEntry(e, beta);
            int global_row = connectivity[e][beta];
            b_(global_row) += b_e_val;
        }
    }
}

void FEMSolver::setStiffnessMatrix(const Eigen::SparseMatrix<double>& K_external)
{
    K_global_ = K_external;
    K_global_.makeCompressed();
}

void FEMSolver::setLoadVector(const Eigen::VectorXd& b_external)
{
    b_ = b_external;
}

void FEMSolver::setTripletList(const std::vector<T_entry>& triplets)
{
    triplet_list_ = triplets;
}

void FEMSolver::applyBoundaryConditions()
{
    // 支持通用边界条件形式：K * (c * du/dn) + L * u = q
    // 三种类型：
    // 1. Dirichlet (K=0, L≠0): u = q/L
    // 2. Neumann (K≠0, L=0): c*du/dn = q/K
    // 3. Robin (K≠0, L≠0): c*du/dn + (L/K)*u = q/K

    const auto& boundarys = config_->getBoundary();
    const auto& field = problem_->getField(field_name_);

    LOG_DEBUG("处理边界条件，共 " + std::to_string(boundarys.size()) + " 个边界单元");

    // 统计不同类型的边界条件
    int dirichlet_count = 0;
    int neumann_count = 0;
    int robin_count = 0;

    // 标记 Dirichlet 边界节点
    std::vector<bool> is_dirichlet_node(N_, false);
    std::vector<double> dirichlet_values(N_, 0.0);

    // 第一遍：处理 Neumann 和 Robin 边界条件，直接添加到三元组列表
    for (size_t boundary_idx = 0; boundary_idx < boundarys.size(); ++boundary_idx)
    {
        const auto& boundary = boundarys[boundary_idx];

        // 通过几何实体ID查找边界条件
        int entity_id = config_->getBoundaryGeometricEntity(boundary_idx);

        // 检查该实体是否有边界条件定义
        // 首先查找具体entity_id，如果没有则查找"all"(-1)
        auto it_specific = field.entity_to_bc.find(entity_id);
        auto it_all = field.entity_to_bc.find(-1);

        if (it_specific == field.entity_to_bc.end() && it_all == field.entity_to_bc.end())
        {
            // 没有定义边界条件,跳过(默认为自然边界条件)
            continue;
        }

        // 优先使用具体的entity_id，其次使用"all"
        const auto& bc =
            (it_specific != field.entity_to_bc.end()) ? it_specific->second : it_all->second;

        // 判断边界条件类型
        if (bc.K_bc == 0)
        {
            // Dirichlet 边界条件：只标记，稍后处理
            dirichlet_count++;
            for (int node_idx : boundary.global_node_indices_in_element)
            {
                if (node_idx >= 0 && node_idx < N_)
                {
                    is_dirichlet_node[node_idx] = true;
                    dirichlet_values[node_idx] = bc.q_bc / bc.L_bc;
                }
            }
        }
        else if (bc.L_bc == 0.0)
        {
            // Neumann 边界条件：只修改载荷向量
            neumann_count++;
            const auto& edge_nodes = boundary.global_node_indices_in_element;

            for (size_t local_i = 0; local_i < edge_nodes.size(); ++local_i)
            {
                int global_i = edge_nodes[local_i];
                double load_contribution = calculateBoundaryLoad(boundary, local_i, boundary_idx);
                b_(global_i) += load_contribution;
            }
        }
        else
        {
            // Robin 边界条件：直接添加到三元组列表
            robin_count++;
            const auto& edge_nodes = boundary.global_node_indices_in_element;

            for (size_t local_i = 0; local_i < edge_nodes.size(); ++local_i)
            {
                int global_i = edge_nodes[local_i];

                for (size_t local_j = 0; local_j < edge_nodes.size(); ++local_j)
                {
                    int global_j = edge_nodes[local_j];
                    double stiffness_contribution =
                        calculateBoundaryStiffness(boundary, local_i, local_j, boundary_idx);
                    triplet_list_.push_back(T_entry(global_i, global_j, stiffness_contribution));
                }

                double load_contribution = calculateBoundaryLoad(boundary, local_i, boundary_idx);
                b_(global_i) += load_contribution;
            }
        }
    }

    LOG_DEBUG("边界条件统计:");
    LOG_DEBUG("  Dirichlet 边界单元数: " + std::to_string(dirichlet_count));
    LOG_DEBUG("  Neumann 边界单元数: " + std::to_string(neumann_count));
    LOG_DEBUG("  Robin 边界单元数: " + std::to_string(robin_count));

    // 统计 Dirichlet 节点
    std::vector<int> dirichlet_nodes;
    double max_dirichlet_val = -1e10;
    double min_dirichlet_val = 1e10;
    for (int i = 0; i < N_; ++i)
    {
        if (is_dirichlet_node[i])
        {
            dirichlet_nodes.push_back(i);
            max_dirichlet_val = std::max(max_dirichlet_val, dirichlet_values[i]);
            min_dirichlet_val = std::min(min_dirichlet_val, dirichlet_values[i]);
        }
    }
    LOG_DEBUG("  Dirichlet 边界节点数: " + std::to_string(dirichlet_nodes.size()));
    LOG_DEBUG("  Dirichlet 值范围: [" + std::to_string(min_dirichlet_val) + ", " +
              std::to_string(max_dirichlet_val) + "]");

    // 第二遍：处理 Dirichlet 边界条件，从三元组中过滤
    if (!dirichlet_nodes.empty())
    {
        // 第一步:先修改右端项(需要用到原始矩阵元素)
        for (const auto& triplet : triplet_list_)
        {
            int row = triplet.row();
            int col = triplet.col();
            double value = triplet.value();

            // 如果列是Dirichlet节点但行不是,修改右端项
            if (is_dirichlet_node[col] && !is_dirichlet_node[row])
            {
                b_(row) -= value * dirichlet_values[col];
            }
        }

        // 第二步:过滤三元组,只保留内部节点间的耦合和Dirichlet节点的对角元素
        std::vector<T_entry> filtered_triplets;
        filtered_triplets.reserve(triplet_list_.size());

        // 标记已经添加过对角元素的Dirichlet节点
        std::vector<bool> dirichlet_diag_added(N_, false);

        for (const auto& triplet : triplet_list_)
        {
            int row = triplet.row();
            int col = triplet.col();

            // 如果是Dirichlet节点的对角元素,只添加一次并设为1
            if (row == col && is_dirichlet_node[row])
            {
                if (!dirichlet_diag_added[row])
                {
                    filtered_triplets.push_back(T_entry(row, col, 1.0));
                    dirichlet_diag_added[row] = true;
                }
                // 否则跳过,避免重复添加
            }
            // 如果行和列都不是Dirichlet节点,保留原值
            else if (!is_dirichlet_node[row] && !is_dirichlet_node[col])
            {
                filtered_triplets.push_back(triplet);
            }
            // 其他情况(行或列是Dirichlet节点但不是对角元素):删除,保持矩阵对称
        }

        // 替换三元组列表
        triplet_list_ = std::move(filtered_triplets);

        // 设置Dirichlet节点的右端项
        for (int i : dirichlet_nodes)
        {
            b_(i) = dirichlet_values[i];
        }
    }

    // 最后：一次性从三元组构建稀疏矩阵
    LOG_DEBUG("  从 " + std::to_string(triplet_list_.size()) + " 个三元组构建稀疏矩阵...");
    K_global_.setFromTriplets(triplet_list_.begin(), triplet_list_.end());
    K_global_.makeCompressed();
    LOG_DEBUG("  矩阵构建完成，非零元素数: " + std::to_string(K_global_.nonZeros()));

    // 调试：检查载荷向量的范围
    LOG_DEBUG("  载荷向量统计:");
    LOG_DEBUG("    最小值: " + std::to_string(b_.minCoeff()));
    LOG_DEBUG("    最大值: " + std::to_string(b_.maxCoeff()));
    LOG_DEBUG("    平均值: " + std::to_string(b_.mean()));

    // 清空三元组列表释放内存
    triplet_list_.clear();
    triplet_list_.shrink_to_fit();
}

void FEMSolver::solve()
{
    // 从ProblemSetup获取求解器参数
    const std::string& solver_type = problem_->solver.solver_type;
    const std::string& preconditioner_type = problem_->solver.preconditioner_type;
    double tol = problem_->solver.tolerance;
    int max_iter = problem_->solver.max_iterations;

    clock_t start = clock();

    int iterations = 0;
    double error = 0.0;

    if (solver_type == "SparseLU")
    {
        // 保留原有的直接求解方法
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(K_global_);
        solver.factorize(K_global_);
        if (solver.info() != Eigen::Success)
        {
            LOG_ERROR("分解失败！");
            return;
        }
        u_ = solver.solve(b_);
        if (solver.info() != Eigen::Success)
        {
            LOG_ERROR("求解失败！");
            return;
        }
    }
    else if (solver_type == "CG")
    {
        // 共轭梯度法 - 使用成员变量 K_global_, b_, u_
        if (preconditioner_type == "DiagonalPreconditioner")
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                                     Eigen::DiagonalPreconditioner<double>>
                solver;
            solver.setMaxIterations(max_iter);
            solver.setTolerance(tol);
            solver.compute(K_global_);  // 直接使用成员变量 K_global_
            u_ = solver.solve(b_);      // 直接使用成员变量 b_，结果存储到成员变量 u_

            iterations = solver.iterations();
            error = solver.error();
        }
        else
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
            solver.setMaxIterations(max_iter);
            solver.setTolerance(tol);
            solver.compute(K_global_);
            u_ = solver.solve(b_);

            iterations = solver.iterations();
            error = solver.error();
        }

        // 输出CG求解器的迭代信息
        LOG_DEBUG("CG求解完成: 迭代 " + std::to_string(iterations) + " 次, 误差 " +
                  std::to_string(error));
    }
    // 可以添加其他迭代求解器...

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    // 输出求解信息
    LOG_DEBUG("求解耗时: " + std::to_string(time_spent) + " 秒");
    double residual_norm = (K_global_ * u_ - b_).norm() / b_.norm();
    LOG_DEBUG("相对残差: " + std::to_string(residual_norm));
}

void FEMSolver::postprocess()
{
    // === 输出求解结果统计信息 ===
    LOG_INFO("求解结果统计 [" + field_name_ + "]:");
    LOG_INFO("  最小值: " + std::to_string(u_.minCoeff()));
    LOG_INFO("  最大值: " + std::to_string(u_.maxCoeff()));
    LOG_INFO("  平均值: " + std::to_string(u_.mean()));

    // 找出异常值
    int count_negative = 0;
    int count_very_large = 0;
    double threshold_large = 1e6;  // 超过1e6认为异常大
    for (int i = 0; i < u_.size(); ++i)
    {
        if (u_(i) < 0)
            count_negative++;
        if (u_(i) > threshold_large)
            count_very_large++;
    }
    if (count_negative > 0)
    {
        LOG_WARNING("  发现 " + std::to_string(count_negative) + " 个负值节点");
    }
    if (count_very_large > 0)
    {
        LOG_WARNING("  发现 " + std::to_string(count_very_large) + " 个异常大值节点 (>" +
                    std::to_string(threshold_large) + ")");
    }

    // 创建误差分析器
    auto errorAnalyzer = std::make_shared<ErrorAnalysis>(config_, problem_, u_, field_name_);

    // 输出误差分析摘要
    errorAnalyzer->printErrorSummary();

    // 输出节点误差详情（前10个节点）
    errorAnalyzer->printDetailedNodeErrors(10);

    // 创建VTK输出对象
    auto vtkOutput = VTKOutputFactory::createVTKOutput(config_, u_);

    // 输出数值解
    vtkOutput->outputNumericalSolution("results/numerical_solution");

    // 检查是否有精确解(V2版本使用FieldConfiguration)
    const auto& field = problem_->getField(field_name_);
    bool has_exact = field.has_exact_solution;

    // 输出精确解
    if (has_exact)
    {
        // 创建Lambda函数包装MaterialProperty的evaluate
        auto exact_func = [&field](const std::vector<double>& coords) -> double
        {
            EvaluationContext ctx;
            ctx.x = (coords.size() >= 1) ? coords[0] : 0.0;
            ctx.y = (coords.size() >= 2) ? coords[1] : 0.0;
            ctx.z = (coords.size() >= 3) ? coords[2] : 0.0;
            return field.exact_solution_u.evaluate(ctx);
        };
        vtkOutput->outputExactSolution("results/exact_solution", exact_func);
        LOG_INFO("2. exact_solution.vtu - 查看精确解分布");
    }

    LOG_INFO(" ParaView可视化指南:");
    LOG_INFO("1. numerical_solution.vtu - 查看数值解分布");
    if (has_exact)
    {
        LOG_INFO("2. exact_solution.vtu - 查看精确解分布");
    }
    else
    {
        LOG_INFO("注意：该问题未定义解析解，仅输出数值解文件");
    }
}

// === 访问器实现 ===

const Eigen::VectorXd& FEMSolver::getSolution() const
{
    return u_;
}

int FEMSolver::getNodesNum() const
{
    return N_;
}

int FEMSolver::getElementsNum() const
{
    return M_;
}

std::shared_ptr<Config> FEMSolver::getConfig() const
{
    return config_;
}

std::shared_ptr<ProblemSetup> FEMSolver::getProblem() const
{
    return problem_;
}

// --- 内部辅助函数定义 ---

double FEMSolver::calculateStiffnessEntry(int e, int alpha, int beta)
{
    // 从Config获取单元节点坐标
    const auto& connectivity = config_->getElementConnectivity();
    const auto& coordinates = config_->getNodeCoordinates();
    const auto& element_entities = config_->getElementGeometricEntities();
    int dimension = config_->getDimension();

    std::vector<double> element_coords(n_ * dimension);
    for (int i = 0; i < n_; ++i)
    {
        int node_idx = connectivity[e][i];
        for (int d = 0; d < dimension; ++d)
        {
            element_coords[i * dimension + d] = coordinates[node_idx * dimension + d];
        }
    }

    // 获取该单元的材料
    int entity_id = element_entities[e];
    std::string material_name = problem_->getMaterialForDomain(entity_id);
    auto material = materials_->getMaterial(material_name);

    // 创建几何映射对象
    auto mapping = GeometryMappingFactory::createMapping(element_coords, config_);

    // 使用工厂创建体单元的高斯积分对象
    auto gaussPoint = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config_->getElementType()),
        config_->getAssembleGaussPointsNum());

    // 获取积分点坐标和权重向量
    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    double entryValue = 0.0;
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        std::vector<double> pointsRef(dimension);
        for (int d = 0; d < dimension; ++d)
        {
            pointsRef[d] = points[gpIndex * dimension + d];
        }

        // 计算形函数在参考坐标系下的梯度
        auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config_);
        std::vector<double> trialGradientsRef =
            shapeFunction->computeTrialGradients(alpha, pointsRef);
        std::vector<double> testGradientsRef = shapeFunction->computeTestGradients(beta, pointsRef);

        // 转换为物理坐标系下的导数
        std::vector<double> trialGradientsPhys;
        std::vector<double> testGradientsPhys;
        mapping->transformGradient(trialGradientsRef, trialGradientsPhys, pointsRef);
        mapping->transformGradient(testGradientsRef, testGradientsPhys, pointsRef);

        // 获取物理坐标（用于计算变系数）
        std::vector<double> pointsPhys;
        mapping->mapToPhysical(pointsRef, pointsPhys);

        // 构造评估上下文
        EvaluationContext ctx;
        ctx.x = (dimension >= 1) ? pointsPhys[0] : 0.0;
        ctx.y = (dimension >= 2) ? pointsPhys[1] : 0.0;
        ctx.z = (dimension >= 3) ? pointsPhys[2] : 0.0;

        // 从材料获取扩散系数 (对于电场是sigma,对于热场是k)
        double coeff = 1.0;
        if (field_name_ == "electric")
        {
            coeff = material->electrical_conductivity.evaluate(ctx);
        }
        else if (field_name_ == "thermal")
        {
            coeff = material->thermal_conductivity.evaluate(ctx);
        }
        else
        {
            throw std::runtime_error("未知的场类型: " + field_name_);
        }

        // 计算积分被积函数：c(x,y) * ∇φ_alpha · ∇φ_beta
        double integrand = 0.0;
        for (int d = 0; d < dimension; ++d)
        {
            integrand += trialGradientsPhys[d] * testGradientsPhys[d];
        }
        integrand *= coeff;

        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }
    return entryValue;
}

double FEMSolver::calculateLoadEntry(int e, int beta)
{
    // 从Config获取单元节点坐标
    const auto& connectivity = config_->getElementConnectivity();
    const auto& coordinates = config_->getNodeCoordinates();
    const auto& element_entities = config_->getElementGeometricEntities();  // 获取单元的几何实体ID
    int dimension = config_->getDimension();

    std::vector<double> element_coords(n_ * dimension);
    for (int i = 0; i < n_; ++i)
    {
        int node_idx = connectivity[e][i];
        for (int d = 0; d < dimension; ++d)
        {
            element_coords[i * dimension + d] = coordinates[node_idx * dimension + d];
        }
    }

    // 创建几何映射对象
    auto mapping = GeometryMappingFactory::createMapping(element_coords, config_);

    // 使用工厂创建体单元的高斯积分对象
    auto gaussPoint = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config_->getElementType()),
        config_->getAssembleGaussPointsNum());

    // 获取积分点坐标和权重向量
    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    // 获取当前场的配置
    const auto& field = problem_->getField(field_name_);

    // 获取该单元的几何实体ID，用于查找对应的源项
    int entity_id = element_entities[e];

    double entryValue = 0.0;
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        std::vector<double> pointsRef(dimension);
        for (int d = 0; d < dimension; ++d)
        {
            pointsRef[d] = points[gpIndex * dimension + d];
        }

        // 获取物理坐标
        std::vector<double> pointsPhys;
        mapping->mapToPhysical(pointsRef, pointsPhys);

        // 计算形函数值
        auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config_);
        double N_test_beta = shapeFunction->computeTestFunction(beta, pointsRef);

        // 构造评估上下文
        EvaluationContext ctx;
        ctx.x = (dimension >= 1) ? pointsPhys[0] : 0.0;
        ctx.y = (dimension >= 2) ? pointsPhys[1] : 0.0;
        ctx.z = (dimension >= 3) ? pointsPhys[2] : 0.0;

        // 计算源项（从 FieldConfiguration 根据单元的几何实体ID获取）
        double f_xy = 0.0;
        if (!field.source_computed)  // 如果不是由求解器计算的源项
        {
            // 首先查找具体entity_id的源项，如果没有则使用默认源项(-1)
            auto it_specific = field.domain_to_source.find(entity_id);
            auto it_default = field.domain_to_source.find(-1);

            if (it_specific != field.domain_to_source.end())
            {
                f_xy = it_specific->second.evaluate(ctx);
            }
            else if (it_default != field.domain_to_source.end())
            {
                f_xy = it_default->second.evaluate(ctx);
            }
            // 否则源项为0（未指定）
        }
        // 否则源项为0,将由耦合求解器在外部设置

        double integrand = f_xy * N_test_beta;
        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }
    return entryValue;
}

/**
 * @brief 计算边界单元上的刚度矩阵贡献（用于 Robin 边界条件）
 * @note 边界单元在不同维度: 1D点, 2D边, 3D面
 */
double FEMSolver::calculateBoundaryStiffness(const Config::Boundary& boundary, int local_i,
                                             int local_j, size_t boundary_idx)
{
    // 获取边界单元的节点索引
    const auto& boundary_nodes = boundary.global_node_indices_in_element;
    int boundary_nodes_num = config_->getBoundaryNodesPerElement();
    if (static_cast<int>(boundary_nodes.size()) != boundary_nodes_num)
    {
        throw std::runtime_error("边界节点数不匹配");
    }

    // 获取节点坐标
    const auto& coordinates = config_->getNodeCoordinates();
    int embed_dim = config_->getDimension();

    // 构造边界单元的坐标数组
    std::vector<double> boundary_coords(boundary_nodes_num * embed_dim);
    for (int i = 0; i < boundary_nodes_num; ++i)
    {
        int node_idx = boundary_nodes[i];
        for (int d = 0; d < embed_dim; ++d)
        {
            boundary_coords[i * embed_dim + d] = coordinates[node_idx * embed_dim + d];
        }
    }

    // 使用工厂创建边界单元的几何映射对象
    auto mapping = GeometryMappingFactory::createMapping(boundary_coords, config_, true);

    // 使用工厂创建边界单元的形函数对象
    auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config_, true);

    // 使用工厂创建边界单元的高斯积分对象
    auto gaussPoint = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config_->getBoundaryElementType()),
        config_->getBoundaryGaussPointsNum());
    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    // 通过几何实体ID获取边界条件
    int entity_id = config_->getBoundaryGeometricEntity(boundary_idx);
    const auto& field = problem_->getField(field_name_);

    // 检查该实体是否有边界条件定义
    auto it_specific = field.entity_to_bc.find(entity_id);
    auto it_all = field.entity_to_bc.find(-1);

    if (it_specific == field.entity_to_bc.end() && it_all == field.entity_to_bc.end())
    {
        return 0.0;  // 没有边界条件,返回0
    }
    const auto& bc =
        (it_specific != field.entity_to_bc.end()) ? it_specific->second : it_all->second;

    // 边界条件系数
    double L_over_K = bc.L_bc / bc.K_bc;

    // 调试：输出第一个Robin边界的参数
    static bool first_robin_logged = false;
    if (!first_robin_logged && bc.K_bc != 0 && bc.L_bc != 0.0)
    {
        LOG_DEBUG("Robin边界参数: K=" + std::to_string(bc.K_bc) + ", L=" + std::to_string(bc.L_bc) +
                  ", q=" + std::to_string(bc.q_bc));
        LOG_DEBUG("  L/K = " + std::to_string(L_over_K));
        first_robin_logged = true;
    }

    // 边界参考坐标维度 = 体维度 - 1
    int boundary_dim = config_->getDimension() - 1;

    double entryValue = 0.0;
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        // 动态构造边界参考坐标
        std::vector<double> coord_ref(boundary_dim);
        for (int d = 0; d < boundary_dim; ++d)
        {
            coord_ref[d] = points[gpIndex * boundary_dim + d];
        }

        // 使用形函数类计算形函数值
        double N_i = shapeFunction->computeTestFunction(local_i, coord_ref);
        double N_j = shapeFunction->computeTrialFunction(local_j, coord_ref);

        // 使用几何映射获取雅可比
        double jacobian = mapping->getJacobianDet(coord_ref);

        // 累加积分
        entryValue += L_over_K * N_i * N_j * jacobian * weights[gpIndex];
    }

    return entryValue;
}

/**
 * @brief 计算边界单元上的载荷向量贡献（用于 Neumann 和 Robin 边界条件）
 * @note 边界单元在不同维度: 1D点, 2D边, 3D面
 */
double FEMSolver::calculateBoundaryLoad(const Config::Boundary& boundary, int local_i,
                                        size_t boundary_idx)
{
    // 获取边界单元的节点索引
    const auto& boundary_nodes = boundary.global_node_indices_in_element;
    int boundary_nodes_num = config_->getBoundaryNodesPerElement();
    if (static_cast<int>(boundary_nodes.size()) != boundary_nodes_num)
    {
        throw std::runtime_error("边界节点数不匹配");
    }

    // 获取节点坐标
    const auto& coordinates = config_->getNodeCoordinates();
    int embed_dim = config_->getDimension();

    // 构造边界单元的坐标数组
    std::vector<double> boundary_coords(boundary_nodes_num * embed_dim);
    for (int i = 0; i < boundary_nodes_num; ++i)
    {
        int node_idx = boundary_nodes[i];
        for (int d = 0; d < embed_dim; ++d)
        {
            boundary_coords[i * embed_dim + d] = coordinates[node_idx * embed_dim + d];
        }
    }

    // 使用工厂创建边界单元的几何映射对象
    auto mapping = GeometryMappingFactory::createMapping(boundary_coords, config_, true);

    // 使用工厂创建边界单元的形函数对象
    auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config_, true);

    // 使用工厂创建边界单元的高斯积分对象
    auto gaussPoint = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config_->getBoundaryElementType()),
        config_->getBoundaryGaussPointsNum());
    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    // 通过几何实体ID获取边界条件
    int entity_id = config_->getBoundaryGeometricEntity(boundary_idx);
    const auto& field = problem_->getField(field_name_);

    // 检查该实体是否有边界条件定义
    auto it_specific = field.entity_to_bc.find(entity_id);
    auto it_all = field.entity_to_bc.find(-1);

    if (it_specific == field.entity_to_bc.end() && it_all == field.entity_to_bc.end())
    {
        return 0.0;  // 没有边界条件,返回0
    }
    const auto& bc =
        (it_specific != field.entity_to_bc.end()) ? it_specific->second : it_all->second;

    // 边界参考坐标维度 = 体维度 - 1
    int boundary_dim = config_->getDimension() - 1;

    // 调试：输出第一个Robin边界载荷参数
    static bool first_robin_load_logged = false;
    if (!first_robin_load_logged && bc.K_bc != 0)
    {
        double q_over_K = bc.q_bc / bc.K_bc;
        LOG_DEBUG("Robin边界载荷: q=" + std::to_string(bc.q_bc) + ", K=" + std::to_string(bc.K_bc));
        LOG_DEBUG("  q/K = " + std::to_string(q_over_K));
        first_robin_load_logged = true;
    }

    double entryValue = 0.0;
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        // 动态构造边界参考坐标
        std::vector<double> coord_ref(boundary_dim);
        for (int d = 0; d < boundary_dim; ++d)
        {
            coord_ref[d] = points[gpIndex * boundary_dim + d];
        }

        // 映射到物理坐标
        std::vector<double> coord_phys;
        mapping->mapToPhysical(coord_ref, coord_phys);

        // 计算边界条件右端项的值
        double q_value = bc.q_bc;
        double q_over_K = q_value / bc.K_bc;

        // 使用形函数类计算形函数值
        double N_i = shapeFunction->computeTrialFunction(local_i, coord_ref);

        // 使用几何映射获取雅可比
        double jacobian = mapping->getJacobianDet(coord_ref);

        // 累加积分
        entryValue += q_over_K * N_i * jacobian * weights[gpIndex];
    }

    return entryValue;
}

// === 梯度计算实现（用于电热耦合） ===

double FEMSolver::computeGradientNormSquared(int element_idx, const std::vector<double>& ref_coords)
{
    // 从Config获取单元节点坐标
    const auto& connectivity = config_->getElementConnectivity();
    const auto& coordinates = config_->getNodeCoordinates();
    int dimension = config_->getDimension();

    // 检查单元索引有效性
    if (element_idx < 0 || element_idx >= M_)
    {
        throw std::runtime_error("Invalid element index in computeGradientNormSquared");
    }

    // 提取单元坐标
    std::vector<double> element_coords(n_ * dimension);
    for (int i = 0; i < n_; ++i)
    {
        int node_idx = connectivity[element_idx][i];
        for (int d = 0; d < dimension; ++d)
        {
            element_coords[i * dimension + d] = coordinates[node_idx * dimension + d];
        }
    }

    // 创建几何映射对象
    auto mapping = GeometryMappingFactory::createMapping(element_coords, config_);

    // 创建形函数对象
    auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config_);

    // 计算该点处的解梯度
    std::vector<double> grad_u_phys(dimension, 0.0);

    for (int i = 0; i < n_; ++i)
    {
        // 获取节点i的解值
        int node_idx = connectivity[element_idx][i];
        double u_i = u_(node_idx);

        // 计算形函数在参考坐标下的梯度
        std::vector<double> grad_N_ref = shapeFunction->computeTrialGradients(i, ref_coords);

        // 转换为物理坐标下的梯度
        std::vector<double> grad_N_phys;
        mapping->transformGradient(grad_N_ref, grad_N_phys, ref_coords);

        // 累加到解的梯度：∇u = Σ u_i * ∇N_i
        for (int d = 0; d < dimension; ++d)
        {
            grad_u_phys[d] += u_i * grad_N_phys[d];
        }
    }

    // 计算梯度范数的平方：|∇u|² = (∂u/∂x)² + (∂u/∂y)²
    double norm_squared = 0.0;
    for (int d = 0; d < dimension; ++d)
    {
        norm_squared += grad_u_phys[d] * grad_u_phys[d];
    }

    return norm_squared;
}
