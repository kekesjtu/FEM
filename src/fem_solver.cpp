#include "fem_solver_2d.h"
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
const int n = 3;  // 一阶三角形单元有3个节点
Eigen::SparseMatrix<double> K_global;
Eigen::VectorXd b;
Eigen::VectorXd u;
std::vector<std::vector<int>> T;
std::vector<double> P;
std::vector<BoundaryEdge> boundary_edges;

Config config;  // 配置对象，包含网格和问题定义
auto shapeFunction = ShapeFunctionFactory::createShapeFunction(std::make_shared<Config>(config));

// --- 内部辅助函数声明 ---
double calculateStiffnessEntry2D(int e, int alpha, int beta, int numGaussPoints);
double calculateLoadEntry2D(int e, int beta, int numGaussPoints);

// --- FEM 函数实现 ---

void preprocess2D()
{
    // 调用defineProblem函数获取问题定义，获得N,M,生成网格，P,T,boundary_edges
    defineProblem_by_mesh_importer();

    // 设置全局变量N和M
    N = static_cast<int>(P.size()) / 2;  // 节点数 = 坐标数 / 2
    M = static_cast<int>(T.size());      // 单元数

    // 初始化矩阵和向量
    K_global.resize(N, N);
    K_global.reserve(Eigen::VectorXi::Constant(N, 9));  // 每行预计有9个非零元素（2D情况）
    b = Eigen::VectorXd::Zero(N);
    u = Eigen::VectorXd::Zero(N);

    std::cout << "网格生成完成：" << N << " 个节点，" << M << " 个单元" << std::endl;
    std::cout << "边界边数量：" << boundary_edges.size() << std::endl;
}

void assemble2D()
{
    // 使用三元组列表收集所有元素
    typedef Eigen::Triplet<double> T_entry;
    std::vector<T_entry> tripletList;
    tripletList.reserve(9 * M);  // 每个单元贡献9个矩阵元素

    for (int e = 0; e < M; ++e)
    {
        const int numGaussPoints = 3;  // 使用3点高斯积分
        for (int alpha = 0; alpha < n; ++alpha)
        {
            for (int beta = 0; beta < n; ++beta)
            {
                double K_e_val = calculateStiffnessEntry2D(e, alpha, beta, numGaussPoints);
                int global_row = T[e][beta];
                int global_col = T[e][alpha];
                tripletList.push_back(T_entry(global_row, global_col, K_e_val));
            }
        }
        for (int beta = 0; beta < n; ++beta)
        {
            double b_e_val = calculateLoadEntry2D(e, beta, numGaussPoints);
            int global_row = T[e][beta];
            b(global_row) += b_e_val;
        }
    }

    // 从三元组列表构建稀疏矩阵
    K_global.setFromTriplets(tripletList.begin(), tripletList.end());
    K_global.makeCompressed();
}

void applyBoundaryConditions2D()
{
    std::cout << "开始应用边界条件..." << std::endl;

    // 支持通用边界条件形式：K * (c * du/dn) + L * u = q

    // 首先收集所有边界节点及其边界条件
    std::vector<bool> is_boundary_node(N, false);
    std::vector<BoundaryCondition> node_bc(N);

    for (const auto& edge : boundary_edges)
    {
        int node1 = edge.global_node_index1;
        int node2 = edge.global_node_index2;

        // 标记边界节点并设置边界条件
        is_boundary_node[node1] = true;
        is_boundary_node[node2] = true;
        node_bc[node1] = edge.bc;
        node_bc[node2] = edge.bc;
    }

    // 确保矩阵已压缩
    K_global.makeCompressed();

    // 收集所有狄利克雷边界条件节点，批量处理
    std::vector<int> dirichlet_nodes;
    std::vector<double> dirichlet_values;

    for (int i = 0; i < N; ++i)
    {
        if (is_boundary_node[i])
        {
            const BoundaryCondition& bc = node_bc[i];

            if (bc.K_bc == 0)
            {
                // 狄利克雷边界条件
                double boundary_value = bc.q_bc / bc.L_bc;
                dirichlet_nodes.push_back(i);
                dirichlet_values.push_back(boundary_value);
            }
            else
            {
                // 诺曼或罗宾边界条件
                if (bc.L_bc == 0.0)
                {
                    b(i) += bc.q_bc;
                }
                else
                {
                    K_global.coeffRef(i, i) += bc.L_bc;
                    b(i) += bc.q_bc;
                }
            }
        }
    }

    std::cout << "狄利克雷边界条件节点数量: " << dirichlet_nodes.size() << std::endl;

    // 批量处理狄利克雷边界条件
    for (size_t idx = 0; idx < dirichlet_nodes.size(); ++idx)
    {
        int i = dirichlet_nodes[idx];  // i存储迪利克雷节点的全局编号
        double boundary_value = dirichlet_values[idx];

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

void solveLinearSystem2D(const std::string& solver_type, const std::string& preconditioner_type,
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

void postprocess2D()
{
    using namespace std;

    cout << "网格信息: " << N << " 个节点，" << M << " 个单元" << endl;
    cout << "边界边数量: " << boundary_edges.size() << endl;

    // 初始化全局配置
    FEMConfig::initialize(2);

    // 创建临时mesh对象（从全局变量转换）
    auto mesh = std::make_shared<TriangleMesh2D>();
    mesh->setNodes(N, P);
    mesh->setElements(M, T);

    // 创建重构后的误差分析器
    auto errorAnalyzer =
        ErrorAnalysisFactory::create2DWithFullSolution(mesh,                  // mesh对象
                                                       exact_solution_u,      // 精确解
                                                       exact_solution_du_dx,  // x方向导数
                                                       exact_solution_du_dy   // y方向导数
        );

    // 设置计算参数
    errorAnalyzer->setSamplingPoints(5);  // L∞误差采样点数
    errorAnalyzer->setGaussPoints(3);     // 数值积分点数

    // 批量计算所有范数误差
    errorAnalyzer->computeAllNormErrors(u);

    // 输出详细的误差分析摘要
    errorAnalyzer->printErrorSummary();

    // 输出节点误差详情（前10个节点）
    errorAnalyzer->printDetailedNodeErrors(u, 10);

    // 输出VTK文件用于ParaView可视化
    cout << "\n--- 二维VTK文件输出 ---" << endl;

    // 创建VTKOutput2D对象，注入mesh2D对象
    auto vtkOutput2D = std::make_shared<VTKOutput2D>(mesh);

    // 输出二维数值解
    vtkOutput2D->outputNumericalSolution("results/numerical_solution", u);

    // 输出二维解析解
    vtkOutput2D->outputExactSolution("results/exact_solution", exact_solution_u);

    // 输出二维对比文件
    vtkOutput2D->outputDenseSamplingError("results/comparison", u, exact_solution_u);

    cout << "\n ParaView可视化指南:" << endl;
    cout << "1. numerical_solution.vtu - 查看数值解分布" << endl;
    cout << "2. exact_solution.vtu     - 查看解析解分布" << endl;
    cout << "3. comparison_dense.vtu   - 误差分析" << endl;
}

// --- 内部辅助函数定义 ---

double calculateStiffnessEntry2D(int e, int alpha, int beta, int numGaussPoints)
{
    // 获取单元节点坐标（使用全局变量作为临时方案）
    std::vector<double> element_coords(6);
    for (int i = 0; i < 3; ++i)
    {
        int node_idx = T[e][i];
        element_coords[i * 2] = P[node_idx * 2];          // x坐标
        element_coords[i * 2 + 1] = P[node_idx * 2 + 1];  // y坐标
    }

    // 创建几何映射对象
    auto mapping = std::make_unique<GeometryMapping2D>(element_coords);

    // 使用新的高斯点类
    auto gaussPoint = GaussPointFactory::createGaussPoint(GaussPointFactory::ElementType::Triangle,
                                                          numGaussPoints);

    // 获取积分点坐标和权重向量
    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    double entryValue = 0.0;
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        std::vector<double> pointsRef = {points[gpIndex * 2], points[gpIndex * 2 + 1]};

        // 计算形函数在参考坐标系下的梯度
        std::vector<double> trialGradientsRef;
        std::vector<double> testGradientsRef;
        trialGradientsRef = shapeFunction->computeTrialGradients(alpha, pointsRef);
        testGradientsRef = shapeFunction->computeTestGradients(beta, pointsRef);

        // 转换为物理坐标系下的导数
        std::vector<double> trialGradientsPhys;
        std::vector<double> testGradientsPhys;
        mapping->transformGradient(trialGradientsRef, trialGradientsPhys, pointsRef);
        mapping->transformGradient(testGradientsRef, testGradientsPhys, pointsRef);

        // 计算积分被积函数（假设扩散系数为1）
        double integrand = (trialGradientsPhys[0] * testGradientsPhys[0] +
                            trialGradientsPhys[1] * testGradientsPhys[1]);

        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }
    return entryValue;
}

double calculateLoadEntry2D(int e, int beta, int numGaussPoints)
{
    // 获取单元节点坐标（使用全局变量作为临时方案）
    std::vector<double> element_coords(6);
    for (int i = 0; i < 3; ++i)
    {
        int node_idx = T[e][i];
        element_coords[i * 2] = P[node_idx * 2];          // x坐标
        element_coords[i * 2 + 1] = P[node_idx * 2 + 1];  // y坐标
    }

    // 创建几何映射对象
    auto mapping = std::make_unique<GeometryMapping2D>(element_coords);

    // 使用新的高斯点类
    auto gaussPoint = GaussPointFactory::createGaussPoint(GaussPointFactory::ElementType::Triangle,
                                                          numGaussPoints);

    // 获取积分点坐标和权重向量
    const auto& points = gaussPoint->getPoints();
    const auto& weights = gaussPoint->getWeights();

    double entryValue = 0.0;
    for (int gpIndex = 0; gpIndex < gaussPoint->getNumPoints(); ++gpIndex)
    {
        std::vector<double> pointsRef = {points[gpIndex * 2], points[gpIndex * 2 + 1]};

        // 获取物理坐标
        std::vector<double> pointsPhys;
        mapping->mapToPhysical(pointsRef, pointsPhys);

        // 计算形函数值
        double N_test_beta = shapeFunction->computeTestFunction(beta, pointsRef);

        // 计算电荷密度（源项）
        const double f_xy = source_term_f(pointsPhys[0], pointsPhys[1]);
        double integrand = f_xy * N_test_beta;
        entryValue += integrand * mapping->getJacobianDet(pointsRef) * weights[gpIndex];
    }
    return entryValue;
}
