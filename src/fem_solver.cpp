#include "fem_solver.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include "error_analysis.h"
#include "gauss_quadrature.h"
#include "geometry_mapping.h"
#include "shape_functions.h"
#include "vtk_output.h"

// --- 全局变量定义 ---
int N, M;
int n;  // 每个单元的节点数，根据单元类型动态确定
Eigen::SparseMatrix<double> K_global;
Eigen::VectorXd b;
Eigen::VectorXd u;

std::shared_ptr<Config> config = std::make_shared<Config>();  // 配置对象，包含网格和问题定义

// --- 内部辅助函数声明 ---
double calculateStiffnessEntry(int e, int alpha, int beta, int numGaussPoints);
double calculateLoadEntry(int e, int beta, int numGaussPoints);

// --- FEM 函数实现 ---

void preprocess()
{
    // 从Config获取网格信息
    N = config->getNodesNum();
    M = config->getElementsNum();

    // 动态确定单元节点数（根据网格类型）
    const auto& connectivity = config->getElementConnectivity();
    if (!connectivity.empty())
    {
        n = connectivity[0].size();  // 第一个单元的节点数
    }
    else
    {
        n = 3;  // 默认三角形单元
    }

    // 动态计算稀疏矩阵预留空间
    int dimension = config->getDimension();
    int estimated_nnz_per_row = n * n;  // 每个节点大约连接n²个其他节点

    // 初始化矩阵和向量
    K_global.resize(N, N);
    K_global.reserve(Eigen::VectorXi::Constant(N, estimated_nnz_per_row));
    b = Eigen::VectorXd::Zero(N);
    u = Eigen::VectorXd::Zero(N);

    std::cout << "网格生成完成：" << N << " 个节点，" << M << " 个单元" << std::endl;
    std::cout << "维度: " << dimension << "D" << std::endl;
    std::cout << "单元类型: 每个单元" << n << "个节点" << std::endl;
}

void assemble()
{
    // 使用三元组列表收集所有元素
    typedef Eigen::Triplet<double> T_entry;
    std::vector<T_entry> tripletList;
    tripletList.reserve(n * n * M);  // 每个单元贡献n²个矩阵元素

    const auto& connectivity = config->getElementConnectivity();

    // 使用config中设置的高斯积分点数
    int numGaussPoints = config->getAssembleGaussPoints();

    for (int e = 0; e < M; ++e)
    {
        for (int alpha = 0; alpha < n; ++alpha)
        {
            for (int beta = 0; beta < n; ++beta)
            {
                double K_e_val = calculateStiffnessEntry(e, alpha, beta, numGaussPoints);
                int global_row = connectivity[e][beta];
                int global_col = connectivity[e][alpha];
                tripletList.push_back(T_entry(global_row, global_col, K_e_val));
            }
        }
        for (int beta = 0; beta < n; ++beta)
        {
            double b_e_val = calculateLoadEntry(e, beta, numGaussPoints);
            int global_row = connectivity[e][beta];
            b(global_row) += b_e_val;
        }
    }

    // 从三元组列表构建稀疏矩阵
    K_global.setFromTriplets(tripletList.begin(), tripletList.end());
    K_global.makeCompressed();
}

void applyBoundaryConditions()
{
    std::cout << "开始应用边界条件..." << std::endl;

    // 支持通用边界条件形式：K * (c * du/dn) + L * u = q

    // 应用默认狄利克雷边界条件（u=0）到所有边界节点
    // 注意：实际的边界条件信息需要从COMSOL网格中提取或手动指定
    std::vector<bool> is_boundary_node(N, false);

    // 简化处理：将所有边界节点设为齐次狄利克雷条件
    // TODO: 根据实际问题设置合适的边界条件
    const auto& coords = config->getNodeCoordinates();
    int dimension = config->getDimension();

    for (int i = 0; i < N; ++i)
    {
        // 通用的边界检测（假设单位球/圆边界）
        double distance_squared = 0.0;
        for (int d = 0; d < dimension; ++d)
        {
            double coord = coords[i * dimension + d];
            distance_squared += coord * coord;
        }

        // 检查是否在边界上
        if (std::abs(distance_squared - 1.0) < 1e-6)
        {
            is_boundary_node[i] = true;
        }
    }

    // 确保矩阵已压缩
    K_global.makeCompressed();

    // 应用齐次狄利克雷边界条件 u = 0
    std::vector<int> dirichlet_nodes;
    for (int i = 0; i < N; ++i)
    {
        if (is_boundary_node[i])
        {
            dirichlet_nodes.push_back(i);
        }
    }

    std::cout << "狄利克雷边界条件节点数量: " << dirichlet_nodes.size() << std::endl;

    // 批量处理狄利克雷边界条件
    for (int i : dirichlet_nodes)
    {
        double boundary_value = 0.0;  // 齐次边界条件

        // 遍历第i列的非0元素
        for (Eigen::SparseMatrix<double>::InnerIterator it(K_global, i); it; ++it)
        {
            int j = it.row();  // 取当前遍历到的行号
            if (j != i)
            {
                K_global.coeffRef(j, i) = 0.0;
                b(j) -= it.value() * boundary_value;
            }
        }

        // 遍历第i行的非0元素，outerSize()返回列数
        for (int k = 0; k < K_global.outerSize(); ++k)
        {
            if (k != i)
            {
                K_global.coeffRef(i, k) = 0.0;
            }
        }
        K_global.prune(0.0);  // 移除所有值为0的元素，保持矩阵稀疏性

        // 设置对角元素和右端项
        K_global.coeffRef(i, i) = 1.0;
        b(i) = boundary_value;
    }
}

void solveLinearSystem(const std::string& solver_type, const std::string& preconditioner_type,
                       double tol, int max_iter, bool verbose)
{
    clock_t start = clock();

    int iterations = 0;
    double error = 0.0;

    if (solver_type == "SparseLU")
    {
        // 保留原有的直接求解方法
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(K_global);
        solver.factorize(K_global);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "分解失败！" << std::endl;
            return;
        }
        u = solver.solve(b);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "求解失败！" << std::endl;
            return;
        }
    }
    else if (solver_type == "CG")
    {
        // 共轭梯度法 - 使用现有的K_global, b, u
        if (preconditioner_type == "DiagonalPreconditioner")
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                                     Eigen::DiagonalPreconditioner<double>>
                solver;
            solver.setMaxIterations(max_iter);
            solver.setTolerance(tol);
            solver.compute(K_global);  // 直接使用现有的K_global
            u = solver.solve(b);       // 直接使用现有的b，结果存储到现有的u

            iterations = solver.iterations();
            error = solver.error();
        }
        else
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
            solver.setMaxIterations(max_iter);
            solver.setTolerance(tol);
            solver.compute(K_global);
            u = solver.solve(b);

            iterations = solver.iterations();
            error = solver.error();
        }

        if (verbose)  // 通过verbose判断是否输出迭代信息
        {
            std::cout << "CG求解完成：" << std::endl;
            std::cout << "  迭代次数: " << iterations << std::endl;
            std::cout << "  估计误差: " << error << std::endl;
        }
    }
    // 可以添加其他迭代求解器...

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    if (verbose)
    {
        std::cout << "求解耗时: " << time_spent << " 秒" << std::endl;
        double residual_norm = (K_global * u - b).norm() / b.norm();
        std::cout << "相对残差: " << residual_norm << std::endl;
    }
}

void postprocess()
{
    using namespace std;

    int dimension = config->getDimension();
    cout << "网格信息: " << N << " 个节点，" << M << " 个单元" << endl;
    cout << "维度: " << dimension << "D" << endl;
    cout << "单元类型: 每个单元" << n << "个节点" << endl;

    // 创建误差分析器
    auto errorAnalyzer = std::make_shared<ErrorAnalysis>(config, u);

    // 输出误差分析摘要
    errorAnalyzer->printErrorSummary();

    // 输出节点误差详情（前10个节点）
    errorAnalyzer->printDetailedNodeErrors(10);

    // 输出VTK文件用于ParaView可视化
    cout << "\n--- " << dimension << "维VTK文件输出 ---" << endl;

    // 创建VTK输出对象
    auto vtkOutput = VTKOutputFactory::createVTKOutput(config, u);

    // 输出数值解
    vtkOutput->outputNumericalSolution("results/numerical_solution");

    // 输出精确解
    auto exact_func = [&](const std::vector<double>& coords) -> double
    { return config->exact_solution_u(coords); };
    vtkOutput->outputExactSolution("results/exact_solution", exact_func);

    // 输出加密采样误差文件
    vtkOutput->outputDenseSamplingError("results/comparison", config);

    cout << "\n ParaView可视化指南:" << endl;
    cout << "1. numerical_solution.vtu - 查看数值解分布" << endl;
    cout << "2. exact_solution.vtu     - 查看解析解分布" << endl;
    cout << "3. comparison.vtu         - 误差分析" << endl;
}

// --- 内部辅助函数定义 ---

double calculateStiffnessEntry(int e, int alpha, int beta, int numGaussPoints)
{
    // 从Config获取单元节点坐标
    const auto& connectivity = config->getElementConnectivity();
    const auto& coordinates = config->getNodeCoordinates();
    int dimension = config->getDimension();

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

    // 使用新的高斯点类
    auto gaussPoint = GaussPointFactory::createGaussPoint(GaussPointFactory::ElementType::Triangle,
                                                          numGaussPoints);

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
        auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config);
        std::vector<double> trialGradientsRef =
            shapeFunction->computeTrialGradients(alpha, pointsRef);
        std::vector<double> testGradientsRef = shapeFunction->computeTestGradients(beta, pointsRef);

        // 转换为物理坐标系下的导数
        std::vector<double> trialGradientsPhys;
        std::vector<double> testGradientsPhys;
        mapping->transformGradient(trialGradientsRef, trialGradientsPhys, pointsRef);
        mapping->transformGradient(testGradientsRef, testGradientsPhys, pointsRef);

        // 计算积分被积函数（假设扩散系数为1）
        double integrand = 0.0;
        for (int d = 0; d < dimension; ++d)
        {
            integrand += trialGradientsPhys[d] * testGradientsPhys[d];
        }

        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }
    return entryValue;
}

double calculateLoadEntry(int e, int beta, int numGaussPoints)
{
    // 从Config获取单元节点坐标
    const auto& connectivity = config->getElementConnectivity();
    const auto& coordinates = config->getNodeCoordinates();
    int dimension = config->getDimension();

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

    // 使用新的高斯点类
    auto gaussPoint = GaussPointFactory::createGaussPoint(GaussPointFactory::ElementType::Triangle,
                                                          numGaussPoints);

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
        auto shapeFunction = ShapeFunctionFactory::createShapeFunction(config);
        double N_test_beta = shapeFunction->computeTestFunction(beta, pointsRef);

        // 计算源项
        const double f_xy = config->source_term_f(pointsPhys);
        double integrand = f_xy * N_test_beta;
        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }
    return entryValue;
}
