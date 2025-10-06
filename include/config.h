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
        // 三维的四面体
    };

  private:
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
        // 对应 K=1, L=0, q=q_flux
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
    int gauss_assemble_points_;

    // 形函数阶数
    const int order_ = 1;  // 形函数阶数，默认为线性

    // 误差分析相关
    const int max_error_sampling_points_ = 5;
    const int gauss_error_points_ = 3;

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

    void setCoefficientC(std::function<double(const std::vector<double>&)> func)
    {
        coefficient_c_func = func;
    }

    void setSourceTermF(std::function<double(const std::vector<double>&)> func)
    {
        source_term_f_func = func;
    }

    void setExactSolutionU(std::function<double(const std::vector<double>&)> func)
    {
        exact_solution_u_func = func;
    }

    void setExactSolutionGradients(
        std::function<double(const std::vector<double>&, std::vector<double>&)> func)
    {
        exact_solution_gradients_func = func;
    }

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
    Config(int max_error_sampling_points, int gauss_error_points, int order_)
        : max_error_sampling_points_(max_error_sampling_points),
          gauss_error_points_(gauss_error_points),
          order_(order_)
    {
        loadMeshFromFile("circle_mesh2.mphtxt");  // 网格相关变量初始化完毕
        // setGaussAssemblePoints();                 // 矩阵组装计算参数初始化完毕
        setProblem();  // 问题定义初始化完毕
    }

    // 访问器
    int getDimension() const
    {
        return dimension_;
    }
    int getSamplingPoints() const
    {
        return max_error_sampling_points_;
    }
    int getErrorGaussPoints() const
    {
        return gauss_error_points_;
    }

    int getAssembleGaussPoints() const
    {
        return gauss_assemble_points_;
    }

    int getNodesNum() const
    {
        return nodes_num_;
    }

    int getElementsNum() const
    {
        return elements_num_;
    }

    int getNodesPerElement() const
    {
        return nodes_num_per_element_;
    }

    int getOrder() const
    {
        return order_;
    }

    int getElementType() const
    {
        return element_type_;
    }

    const std::vector<double>& getNodeCoordinates() const
    {
        return node_coordinates_;
    }

    const std::vector<std::vector<int>>& getElementConnectivity() const
    {
        return element_connectivity_;
    }

    bool hasExactSolution() const
    {
        return exact_solution_u_func != nullptr;
    }

    bool hasExactGradients() const
    {
        return exact_solution_gradients_func != nullptr;
    }

    // 调用函数
    double coefficient_c(const std::vector<double>& coords) const
    {
        if (coefficient_c_func)
        {
            return coefficient_c_func(coords);
        }
        throw std::runtime_error("coefficient_c function not set");
    }

    double source_term_f(const std::vector<double>& coords) const
    {
        if (source_term_f_func)
        {
            return source_term_f_func(coords);
        }
        throw std::runtime_error("source_term_f function not set");
    }

    double exact_solution_u(const std::vector<double>& coords) const
    {
        if (exact_solution_u_func)
        {
            return exact_solution_u_func(coords);
        }
        throw std::runtime_error("exact_solution_u function not set");
    }

    double exact_solution_gradients(const std::vector<double>& coords,
                                    std::vector<double>& gradients) const
    {
        if (exact_solution_gradients_func)
        {
            return exact_solution_gradients_func(coords, gradients);
        }
        throw std::runtime_error("exact_solution_gradients function not set");
    }
};

#endif  // CONFIG_H
