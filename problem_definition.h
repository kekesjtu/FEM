#ifndef PROBLEM_DEFINITION_H
#define PROBLEM_DEFINITION_H

#include <vector>

// 结构体：定义通用边界条件
// 表示形式为: K * (c * du/dn) + L * u = q
struct BoundaryCondition
{
    int K_bc;     // 导数项的系数
    double L_bc;  // 值项的系数
    double q_bc;  // 右侧项

    // 辅助构造函数
    // 构造狄利克雷边界条件: u = val
    // 对应 K=0, L=1, q=val
    static BoundaryCondition Dirichlet(double val)
    {
        return {0, 1.0, val};
    }
    // 构造罗宾边界条件: c*du/dn + h*u = g
    // 对应 K=1, L=h, q=g
    static BoundaryCondition Robin(double h, double g)
    {
        return {1, h, g};
    }
    // 构造诺曼边界条件: c*du/dn = q_flux
    // 对应 K=1, L=0, q=q_flux
    static BoundaryCondition Neumann(double q_flux)
    {
        return {1, 0.0, q_flux};
    }
};

// 网格生成参数
struct Domain
{
    double left, right, bottom, top;
};

// 结构体：定义边界边信息
struct BoundaryEdge
{
    BoundaryCondition bc;    // 边界条件类型
    int element_index;       // 边界边所属单元索引
    int global_node_index1;  // 边界边的第一个全局节点索引
    int global_node_index2;  // 边界边的第二个全局节点索引
};

/**
 * @brief 定义二维泊松方程中的扩散系数 c(x,y)
 *        对应于方程 -∇·(c(x,y)∇u) = f(x,y)
 *        在静电场问题中，c表示介电常数
 * @param x x坐标（在一维兼容接口中使用）
 * @return double 系数c的值
 */
double coefficient_c(double x, double y);

/**
 * @brief 定义二维泊松方程中的源项 f(x,y)
 *        对应于方程 -∇²u = f(x,y)
 *        在静电场问题中，f表示电荷密度分布
 * @param x x坐标（在一维兼容接口中使用）
 * @return double 源项f的值
 */
double source_term_f(double x, double y);

/**
 * @brief 定义二维泊松方程问题的精确解 u(x,y) (用于后处理中计算误差)
 *        在静电场问题中，u表示静电势φ(x,y)
 * @param x x坐标（在一维兼容接口中使用）
 * @return double 精确解u的值
 */
double exact_solution_u(double x, double y);

/**
 * @brief 定义二维泊松方程问题的精确解导数 du/dx (用于后处理中计算H1误差)
 * @param x x坐标（在一维兼容接口中使用）
 * @return double 精确解导数du/dx的值
 */
double exact_solution_du_dx(double x, double y);

/**
 * @brief 定义二维泊松方程问题的精确解导数 du/dy (用于后处理中计算H1误差)
 * @param x x坐标（在一维兼容接口中使用）
 * @return double 精确解导数du/dy的值
 */
double exact_solution_du_dy(double x, double y);

/**
 * @brief 生成二维三角形网格
 * @param domain 计算域
 * @param N1 x方向单元数
 * @param N2 y方向单元数
 * @param boundary_edges 边界条件信息（输出参数）
 */
void generateMesh2D_rectangle(const Domain& domain, int N1, int N2,
                              std::vector<BoundaryEdge>& boundary_edges);

/**
 * @brief 定义二维泊松方程问题的网格参数，调用generateMesh2D_rectangle，并指定边界条件
 */
void defineProblem();

/**
 * @brief 定义二维泊松方程问题的网格参数，调用ComsolMeshImporter导入网格，并指定边界条件
 */
void defineProblem_by_mesh_importer();

#endif  // PROBLEM_DEFINITION_H
