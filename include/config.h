#ifndef CONFIG_H
#define CONFIG_H

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

// 简单的配置类，用于统一管理与误差分析/积分相关的计算参数。
// 该头同时提供别名 ComputeConfig = Config，以便与旧代码兼容。
class Config
{
  public:
    enum ElementType
    {
        TRIANGLE,
        QUADRILATERAL,
        // 三维还没有实现
    };

    // 问题定义相关结构体
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
        // 对应 K=1, L=0, q_flux
        static BoundaryCondition Neumann(double q_flux)
        {
            return {1, 0.0, q_flux};
        }
    };

    // 结构体：定义边界边信息
    struct Boundary
    {
        BoundaryCondition bc;                             // 边界条件类型
        int element_index;                                // 边界边所属单元索引
        std::vector<int> global_node_indices_in_element;  // 边界边的全局节点索引列表
    };

  private:
    // 网格相关
    std::string mesh_filename_ = "circle_mesh2.mphtxt";  // 网格文件名，用户可直接在此修改
    int dimension_;
    int nodes_num_;
    int elements_num_;
    int nodes_num_per_element_;
    ElementType element_type_;
    std::vector<double> node_coordinates_;  // 节点坐标矩阵 [node_id * dimension + coord_id]
    std::vector<std::vector<int>>
        element_connectivity_;        // 单元连接矩阵 [element_id][local_node_id] = global_node_id
    std::vector<Boundary> boundarys;  // 边界边信息

    // 矩阵组装相关
    int gauss_assemble_points_num_;

    // 形函数阶数
    int shape_function_order_ = 1;  // 形函数阶数默认为线性，用户可直接在此修改

    // 误差分析相关
    int max_error_sampling_points_num_ = 5;  // 最大误差采样点数，用户可直接在此修改
    int gauss_error_points_num_ = 3;         // 高斯误差积分点数，用户可直接在此修改

    // 求解器相关参数
    std::string solver_type_ = "CG";  // 求解器类型："SparseLU"或"CG"，用户可直接在此修改
    std::string preconditioner_type_ =
        "DiagonalPreconditioner";       // 预条件子类型，用户可直接在此修改
    double solver_tolerance_ = 1e-6;    // 求解器收敛容差，用户可直接在此修改
    int solver_max_iterations_ = 1000;  // 求解器最大迭代次数，用户可直接在此修改

    // 文件中将会规定dimension、node_coordinates、element_connectivity、
    // int nodes_num_,int elements_num_,int nodes_num_per_element_;
    void loadMeshFromFile(const std::string& filename);

    // 这里要根据维度和单元类型设置gauss_assemble_points_.n个高斯点能精确积分2n+1阶多项式
    void setGaussAssemblePoints();

    // 使用std::function，更灵活
    std::function<double(const std::vector<double>&)> coefficient_c_func;
    std::function<double(const std::vector<double>&)> source_term_f_func;
    std::function<double(const std::vector<double>&)> exact_solution_u_func;
    std::function<double(const std::vector<double>&, std::vector<double>&)>
        exact_solution_gradients_func;

    void setCoefficientC(std::function<double(const std::vector<double>&)> func);
    void setSourceTermF(std::function<double(const std::vector<double>&)> func);
    void setExactSolutionU(std::function<double(const std::vector<double>&)> func);
    void setExactSolutionGradients(
        std::function<double(const std::vector<double>&, std::vector<double>&)> func);

    // 用户手动配置问题的函数 - 在这里定义具体问题
    void setProblem()
    {
        int dimension = dimension_;
        // 示例1：泊松方程
        setCoefficientC([](const std::vector<double>& coords) -> double { return 1.0; });

        setSourceTermF(
            [dimension](const std::vector<double>& coords) -> double
            {
                if (coords.size() != dimension)
                {
                    throw std::invalid_argument("坐标维度不匹配");
                }
                else
                {
                    return 10.0;
                }
            });

        setExactSolutionU(
            [dimension](const std::vector<double>& coords) -> double
            {
                if (coords.size() != dimension)
                {
                    throw std::invalid_argument("坐标维度不匹配");
                }
                else
                {
                    return 10.0 / 4 * (1 * 1 - coords[0] * coords[0] - coords[1] * coords[1]);
                }
            });

        setExactSolutionGradients(
            [dimension](const std::vector<double>& coords, std::vector<double>& gradients) -> double
            {
                if (coords.size() != dimension || gradients.size() != dimension)
                {
                    throw std::invalid_argument("坐标或梯度维度不匹配");
                }
                else
                {
                    gradients[0] = -10.0 / 2 * coords[0];  // du/dx
                    gradients[1] = -10.0 / 2 * coords[1];  // du/dy
                    return 0.0;                            // 返回值未使用
                }
            });
    }

  public:
    Config();  // 默认构造函数，使用默认参数

    // 访问器
    const std::string& getMeshFilename() const;
    void setMeshFilename(const std::string& filename);
    int getDimension() const;
    int getSamplingPointsNum() const;
    int getErrorGaussPointsNum() const;
    int getAssembleGaussPointsNum() const;
    int getNodesNum() const;
    int getElementsNum() const;
    int getNodesPerElement() const;
    int getOrder() const;
    int getElementType() const;
    const std::vector<double>& getNodeCoordinates() const;
    const std::vector<std::vector<int>>& getElementConnectivity() const;
    bool hasExactSolution() const;
    bool hasExactGradients() const;

    // 求解器参数访问器
    const std::string& getSolverType() const;
    const std::string& getPreconditionerType() const;
    double getSolverTolerance() const;
    int getSolverMaxIterations() const;

    // 调用函数
    double coefficient_c(const std::vector<double>& coords) const;
    double source_term_f(const std::vector<double>& coords) const;
    double exact_solution_u(const std::vector<double>& coords) const;
    double exact_solution_gradients(const std::vector<double>& coords,
                                    std::vector<double>& gradients) const;
};

#endif  // CONFIG_H
