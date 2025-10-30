#include "fem_solver.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <ctime>
#include <iostream>
#include "error_analysis.h"
#include "gauss_quadrature.h"
#include "geometry_mapping.h"
#include "shape_functions.h"
#include "vtk_output.h"

// --- FEMSolver 类实现 ---

// 构造函数
FEMSolver::FEMSolver(std::shared_ptr<Config> config, std::shared_ptr<ProblemSetup> problem)
    : config_(config), problem_(problem), N_(0), M_(0), n_(0)
{
    // 初始化时不做具体操作，等待 preprocess() 调用
    // 所有物理问题设置（边界条件、系数、源项）应在 ProblemSetup 中完成
}

// 高级接口：完整求解
void FEMSolver::solveComplete()
{
    preprocess();
    std::cout << "步骤 1: 网格生成和预处理完成。" << std::endl;

    assemble();
    std::cout << "步骤 2: 全局矩阵组装完成。" << std::endl;

    applyBoundaryConditions();
    std::cout << "步骤 3: 边界条件施加完成。" << std::endl;

    solve();
    std::cout << "步骤 4: 线性方程组求解完成。" << std::endl;

    postprocess();
    std::cout << "步骤 5: 后处理完成。" << std::endl;
}

// 预处理
void FEMSolver::preprocess()
{
    // 从Config获取网格信息
    N_ = config_->getNodesNum();
    M_ = config_->getElementsNum();
    n_ = config_->getNodesNumPerElement();

    int dimension = config_->getDimension();

    // 在预处理阶段设置边界条件（直接操作 config 的 boundarys 数组）
    if (problem_)
    {
        problem_->setupBoundaryConditions(config_);
    }

    // 初始化矩阵和向量
    K_global_.resize(N_, N_);
    // 动态计算稀疏矩阵预留空间每个节点大约连接n²个其他节点
    K_global_.reserve(Eigen::VectorXi::Constant(N_, n_ * n_));
    b_ = Eigen::VectorXd::Zero(N_);
    u_ = Eigen::VectorXd::Zero(N_);

    std::cout << "网格生成完成：" << N_ << " 个节点，" << M_ << " 个单元" << std::endl;
    std::cout << "维度: " << dimension << "D" << std::endl;
    std::cout << "单元类型: 每个单元" << n_ << "个节点" << std::endl;
}

void FEMSolver::assemble()
{
    // 使用三元组列表收集所有元素
    typedef Eigen::Triplet<double> T_entry;
    std::vector<T_entry> tripletList;
    tripletList.reserve(n_ * n_ * M_);  // 每个单元贡献n²个矩阵元素

    const auto& connectivity = config_->getElementConnectivity();

    for (int e = 0; e < M_; ++e)
    {
        for (int alpha = 0; alpha < n_; ++alpha)
        {
            for (int beta = 0; beta < n_; ++beta)
            {
                double K_e_val = calculateStiffnessEntry(e, alpha, beta);
                int global_row = connectivity[e][beta];
                int global_col = connectivity[e][alpha];
                tripletList.push_back(T_entry(global_row, global_col, K_e_val));
            }
        }
        for (int beta = 0; beta < n_; ++beta)
        {
            double b_e_val = calculateLoadEntry(e, beta);
            int global_row = connectivity[e][beta];
            b_(global_row) += b_e_val;
        }
    }

    // 从三元组列表构建稀疏矩阵
    K_global_.setFromTriplets(tripletList.begin(), tripletList.end());
    K_global_.makeCompressed();
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

void FEMSolver::applyBoundaryConditions()
{
    // 支持通用边界条件形式：K * (c * du/dn) + L * u = q
    // 三种类型：
    // 1. Dirichlet (K=0, L≠0): u = q/L
    // 2. Neumann (K≠0, L=0): c*du/dn = q/K
    // 3. Robin (K≠0, L≠0): c*du/dn + (L/K)*u = q/K

    const auto& boundarys = config_->getBoundary();
    const auto& boundary_conditions = problem_->getBoundaryConditions();

    if (boundarys.size() != boundary_conditions.size())
    {
        throw std::runtime_error("边界几何数量与边界条件数量不匹配");
    }

    std::cout << "处理边界条件，共 " << boundarys.size() << " 个边界单元" << std::endl;

    // 统计不同类型的边界条件
    int dirichlet_count = 0;
    int neumann_count = 0;
    int robin_count = 0;

    // 标记 Dirichlet 边界节点
    std::vector<bool> is_dirichlet_node(N_, false);
    std::vector<double> dirichlet_values(N_, 0.0);

    // 首先处理 Neumann 和 Robin 边界条件（修改刚度矩阵和载荷向量）
    for (size_t boundary_idx = 0; boundary_idx < boundarys.size(); ++boundary_idx)
    {
        const auto& boundary = boundarys[boundary_idx];
        const auto& bc = boundary_conditions[boundary_idx];

        // 判断边界条件类型
        if (bc.K_bc == 0)
        {
            // Dirichlet 边界条件：稍后处理
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
            // Robin 边界条件：修改刚度矩阵和载荷向量
            robin_count++;
            const auto& edge_nodes = boundary.global_node_indices_in_element;

            // 修改刚度矩阵
            for (size_t local_i = 0; local_i < edge_nodes.size(); ++local_i)
            {
                int global_i = edge_nodes[local_i];

                for (size_t local_j = 0; local_j < edge_nodes.size(); ++local_j)
                {
                    int global_j = edge_nodes[local_j];
                    double stiffness_contribution =
                        calculateBoundaryStiffness(boundary, local_i, local_j, boundary_idx);
                    K_global_.coeffRef(global_i, global_j) += stiffness_contribution;
                }

                // 修改载荷向量
                double load_contribution = calculateBoundaryLoad(boundary, local_i, boundary_idx);
                b_(global_i) += load_contribution;
            }
        }
    }

    std::cout << "边界条件统计：" << std::endl;
    std::cout << "  Dirichlet 边界单元数: " << dirichlet_count << std::endl;
    std::cout << "  Neumann 边界单元数: " << neumann_count << std::endl;
    std::cout << "  Robin 边界单元数: " << robin_count << std::endl;

    // 确保矩阵已压缩
    K_global_.makeCompressed();

    // 最后处理 Dirichlet 边界条件（强加条件，修改矩阵结构）
    std::vector<int> dirichlet_nodes;
    for (int i = 0; i < N_; ++i)
    {
        if (is_dirichlet_node[i])
        {
            dirichlet_nodes.push_back(i);
        }
    }

    std::cout << "  Dirichlet 边界节点数: " << dirichlet_nodes.size() << std::endl;

    // 批量处理 Dirichlet 边界条件
    for (int i : dirichlet_nodes)
    {
        double boundary_value = dirichlet_values[i];

        // 遍历第i列的非0元素
        for (Eigen::SparseMatrix<double>::InnerIterator it(K_global_, i); it; ++it)
        {
            int j = it.row();  // 取当前遍历到的行号
            if (j != i)
            {
                // ⚠️ 重要：必须先保存it.value()的值，再修改矩阵元素
                // 原因：K_global_.coeffRef(j,i)=0会修改迭代器当前指向的元素
                //      如果先修改再读取it.value()，可能读到已被修改的值（0）
                //      这会导致边界条件的贡献无法正确传递到内部节点
                double contribution = it.value() * boundary_value;

                K_global_.coeffRef(j, i) = 0.0;
                b_(j) -= contribution;
            }
        }

        // 遍历第i行的非0元素，outerSize()返回列数
        for (int k = 0; k < K_global_.outerSize(); ++k)
        {
            if (k != i)
            {
                K_global_.coeffRef(i, k) = 0.0;
            }
        }
        K_global_.prune(0.0);  // 移除所有值为0的元素，保持矩阵稀疏性

        // 设置对角元素和右端项
        K_global_.coeffRef(i, i) = 1.0;
        b_(i) = boundary_value;
    }
}

void FEMSolver::solve()
{
    // 从Config获取求解器参数
    const std::string& solver_type = config_->getSolverType();
    const std::string& preconditioner_type = config_->getPreconditionerType();
    double tol = config_->getSolverTolerance();
    int max_iter = config_->getSolverMaxIterations();

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
            std::cerr << "分解失败！" << std::endl;
            return;
        }
        u_ = solver.solve(b_);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "求解失败！" << std::endl;
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

        // 默认输出CG求解器的迭代信息
        std::cout << "CG求解完成：" << std::endl;
        std::cout << "  迭代次数: " << iterations << std::endl;
        std::cout << "  估计误差: " << error << std::endl;
    }
    // 可以添加其他迭代求解器...

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    // 默认输出求解信息
    std::cout << "求解耗时: " << time_spent << " 秒" << std::endl;
    double residual_norm = (K_global_ * u_ - b_).norm() / b_.norm();
    std::cout << "相对残差: " << residual_norm << std::endl;
}

void FEMSolver::postprocess()
{
    using namespace std;

    int dimension = config_->getDimension();
    cout << "网格信息: " << N_ << " 个节点，" << M_ << " 个单元" << endl;
    cout << "维度: " << dimension << "D" << endl;
    cout << "单元类型: 每个单元" << n_ << "个节点" << endl;

    // 创建误差分析器
    auto errorAnalyzer = std::make_shared<ErrorAnalysis>(config_, problem_, u_);

    // 输出误差分析摘要
    errorAnalyzer->printErrorSummary();

    // 输出节点误差详情（前10个节点）
    errorAnalyzer->printDetailedNodeErrors(10);

    // 输出VTK文件用于ParaView可视化
    cout << "\n--- " << dimension << "维VTK文件输出 ---" << endl;

    // 创建VTK输出对象
    auto vtkOutput = VTKOutputFactory::createVTKOutput(config_, u_);

    // 输出数值解
    vtkOutput->outputNumericalSolution("results/numerical_solution");

    // 输出精确解（如果 ProblemSetup 中有定义）
    bool has_exact = false;
    if (problem_ && problem_->hasExactSolution())
    {
        auto exact_func = [&](const std::vector<double>& coords) -> double
        { return problem_->exactSolutionU(coords); };
        vtkOutput->outputExactSolution("results/exact_solution", exact_func);
        has_exact = true;
    }

    // 输出加密采样误差文件（如果有精确解）
    if (problem_ && problem_->hasExactSolution())
    {
        vtkOutput->outputDenseSamplingError("results/comparison", config_, problem_);
    }

    cout << "\n ParaView可视化指南:" << endl;
    cout << "1. numerical_solution.vtu - 查看数值解分布" << endl;
    if (has_exact)
    {
        cout << "2. exact_solution.vtu     - 查看解析解分布" << endl;
        cout << "3. comparison.vtu         - 误差分析" << endl;
    }
    else
    {
        cout << "注意：该问题未定义解析解，仅输出数值解文件" << endl;
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

        // 计算扩散系数（从 ProblemSetup 获取）
        double coeff = 1.0;
        if (problem_)
        {
            coeff = problem_->coefficient(pointsPhys);
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

        // 计算源项（从 ProblemSetup 获取）
        double f_xy = 0.0;
        if (problem_)
        {
            f_xy = problem_->source(pointsPhys);
        }

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

    // 从 ProblemSetup 获取边界条件
    const auto& bc = problem_->getBoundaryCondition(boundary_idx);

    // 边界条件系数
    double L_over_K = bc.L_bc / bc.K_bc;

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

    // 从 ProblemSetup 获取边界条件
    const auto& bc = problem_->getBoundaryCondition(boundary_idx);

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
