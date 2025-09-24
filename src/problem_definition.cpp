#include "problem_definition.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "comsol_mesh_importer.h"
#include "fem_solver_2d.h"

/**
 * @brief 定义二维泊松方程中的扩散系数 c(x,y)
 *        对应于方程 -∇·(c(x,y)∇u) = f(x,y)
 *        在静电场问题中，c表示介电常数，这里取为常数1
 */
double coefficient_c(double x, double y)
{
    return 1.0;  // 真空中的相对介电常数
}

/**
 * @brief 定义二维泊松方程中的源项 f(x,y)
 *        对应于方程 -∇²u = f(x,y)
 *        在静电场问题中，f表示电荷密度分布
 */
double source_term_f(double x, double y)
{
    const double pi = 3.14159265358979323846;
    return 10.0;
}

/**
 * @brief 定义二维泊松方程问题的精确解 u(x,y) (用于后处理中计算误差)
 *        在静电场问题中，u表示静电势φ(x,y)
 */
double exact_solution_u(double x, double y)
{
    const double pi = 3.14159265358979323846;
    return 10.0 / 4 * (1 * 1 - x * x - y * y);
}

/**
 * @brief 定义二维泊松方程问题的精确解导数 du/dx (用于后处理中计算H1误差)
 */
double exact_solution_du_dx(double x, double y)
{
    const double pi = 3.14159265358979323846;
    return -10.0 / 2 * x;
}

/**
 * @brief 定义二维泊松方程问题的精确解导数 du/dy (用于后处理中计算H1误差)
 */
double exact_solution_du_dy(double x, double y)
{
    const double pi = 3.14159265358979323846;
    return -10.0 / 2 * y;
}

// --- 网格生成函数实现 ---
void generateMesh2D_rectangle(const Domain &domain, int N1, int N2,
                              std::vector<BoundaryEdge> &boundary_edges)
{
    // --- 1. Clear old data ---
    P.clear();
    T.clear();
    boundary_edges.clear();

    if (N1 < 1 || N2 < 1)
    {
        std::cerr << "Interval counts N1 and N2 must be at least 1." << std::endl;
        return;
    }

    // --- 2. Generate Node Coordinates (P vector) ---
    const int Nx = N1 + 1;
    const int Ny = N2 + 1;
    P.reserve(Nx * Ny * 2);
    const double hx = (domain.right - domain.left) / N1;
    const double hy = (domain.top - domain.bottom) / N2;

    // Column-Major order: Outer loop for columns (x), inner for rows (y)
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            P.push_back(domain.left + i * hx);
            P.push_back(domain.bottom + j * hy);
        }
    }

    // --- 3. Generate Triangular Elements (T vector) ---
    T.reserve(2 * N1 * N2);
    auto node_2d_to_1d = [&](int i, int j) { return i * Ny + j; };

    // Column-Major order for elements
    for (int ie = 0; ie < N1; ++ie)
    {
        for (int je = 0; je < N2; ++je)
        {
            int idx_BL = node_2d_to_1d(ie, je);
            int idx_BR = node_2d_to_1d(ie + 1, je);
            int idx_TL = node_2d_to_1d(ie, je + 1);
            int idx_TR = node_2d_to_1d(ie + 1, je + 1);

            // Lower triangle (直角在左下角): 第一个局部编号在左下角 BL
            T.push_back({idx_BL, idx_BR, idx_TL});
            // Upper triangle (直角在右上角): 第一个局部编号在左上角 TL
            T.push_back({idx_TL, idx_BR, idx_TR});
        }
    }

    // --- 4. Generate Boundary Edges ---
    boundary_edges.reserve(2 * N1 + 2 * N2);
    auto element_2d_to_1d = [&](int ie, int je, bool is_upper)
    { return ie * (2 * N2) + 2 * je + (is_upper ? 1 : 0); };

    // BOTTOM boundary (left to right)
    for (int ie = 0; ie < N1; ++ie)
    {
        int node1 = node_2d_to_1d(ie, 0);
        int node2 = node_2d_to_1d(ie + 1, 0);
        int elem_idx = element_2d_to_1d(ie, 0, false);  // Belongs to the lower triangle
        // 默认先设置成Dirichlet边界条件
        boundary_edges.push_back({BoundaryCondition::Dirichlet(0.0), elem_idx, node1, node2});
    }

    // RIGHT boundary (bottom to top)
    for (int je = 0; je < N2; ++je)
    {
        int node1 = node_2d_to_1d(N1, je);
        int node2 = node_2d_to_1d(N1, je + 1);
        int elem_idx = element_2d_to_1d(
            N1 - 1, je, true);  // Belongs to the upper triangle of the last column of cells
        boundary_edges.push_back({BoundaryCondition::Dirichlet(0.0), elem_idx, node1, node2});
    }

    // TOP boundary (right to left)
    for (int ie = N1 - 1; ie >= 0; --ie)
    {
        int node1 = node_2d_to_1d(ie + 1, N2);
        int node2 = node_2d_to_1d(ie, N2);
        int elem_idx = element_2d_to_1d(
            ie, N2 - 1, true);  // Belongs to the upper triangle of the top row of cells
        boundary_edges.push_back({BoundaryCondition::Dirichlet(0.0), elem_idx, node1, node2});
    }

    // LEFT boundary (top to bottom)
    for (int je = N2 - 1; je >= 0; --je)
    {
        int node1 = node_2d_to_1d(0, je + 1);
        int node2 = node_2d_to_1d(0, je);
        int elem_idx = element_2d_to_1d(
            0, je, false);  // Belongs to the lower triangle of the first column of cells
        boundary_edges.push_back({BoundaryCondition::Dirichlet(0.0), elem_idx, node1, node2});
    }
}

/**
 * @brief 定义二维泊松方程问题的网格和边界条件
 *        静电场问题：平行板电容器模型
 *        计算域：单位正方形[0,1]×[0,1]
 *        控制方程：-∇²φ = ρ(x,y)
 */
void defineProblem()
{
    Domain domain;  // 矩形区域
    int N1, N2;

    // 1. 定义计算域
    domain = {0.0, 1.0, 0.0, 1.0};  // 单位正方形

    // 2. 定义网格参数
    N1 = 80;  // x方向格子数
    N2 = 80;  // y方向格子数

    // 生成网格，得到boundary_edges信息
    generateMesh2D_rectangle(domain, N1, N2, boundary_edges);

    // 3. 定义边界条件
    // 在boundary_edges中手动输入边界条件，boundary_edges的索引是按照逆时针方向从左下角单元的底边开始
    for (int i = 0; i < static_cast<int>(boundary_edges.size()); ++i)
    {
        boundary_edges[i].bc = BoundaryCondition::Dirichlet(0.0);
    }
}

// 有了直接导入网格的程序之后，将不需要原来的defineProblem函数。
void defineProblem_by_mesh_importer()
{
    // 1. 创建COMSOL网格导入器
    ComsolMeshImporter importer;

    // 2. 导入网格文件
    std::string mesh_file = "circle_mesh2.mphtxt";
    if (!importer.importMesh(mesh_file))
    {
        std::cerr << "网格导入失败！" << std::endl;
        return;
    }

    // 3. 填充全局数组P、T、boundary_edges（边界条件为默认值）
    importer.populateGlobalArrays();

    // 4. 在defineProblem中显式设置边界条件
    // 对于圆域问题，设置齐次Dirichlet边界条件
    std::cout << "设置边界条件..." << std::endl;
    for (auto &boundary_edge : boundary_edges)
    {
        boundary_edge.bc = BoundaryCondition::Dirichlet(0.0);
    }

    std::cout << "边界条件设置完成，共 " << boundary_edges.size() << " 条边界边" << std::endl;
}