#include "config.h"
#include "comsol_mesh_importer.h"

// 构造函数实现
Config::Config()
{
    loadMeshFromFile(mesh_filename_);  // 网格相关变量初始化完毕
    setGaussAssemblePoints();          // 矩阵组装计算参数初始化完毕
    setProblem();                      // 问题定义初始化完毕
}

// 高斯积分点数设置
void Config::setGaussAssemblePoints()
{
    // 对于形函数阶数order_的单元，刚度矩阵积分涉及2*order_阶多项式
    // n个高斯点能精确积分2n-1阶多项式，所以需要至少(2*shape_function_order_+1+1)/2 = shape_function_order_+1个点
    // 但实际中通常取稍多一些以保证精度

    if (dimension_ == 2)
    {
        if (element_type_ == TRIANGLE)
        {
            // 三角形单元：
            // 1点：精确积分1阶多项式（常数）
            // 3点：精确积分2阶多项式
            // 4点：精确积分3阶多项式
            // 7点：精确积分5阶多项式
            if (shape_function_order_ == 1)
            {
                gauss_assemble_points_num_ = 3;  // 线性单元用3点高斯积分
            }
            else if (shape_function_order_ == 2)
            {
                gauss_assemble_points_num_ = 4;  // 二次单元用4点高斯积分
            }
            else
            {
                gauss_assemble_points_num_ = 7;  // 高阶单元用7点高斯积分
            }
        }
    }
    else
    {
        // 一维情况或其他未实现的维度
        gauss_assemble_points_num_ = shape_function_order_ + 1;
    }

    std::cout << "设置高斯积分点数: " << gauss_assemble_points_num_ << " (维度=" << dimension_
              << ", 单元类型=" << element_type_ << ", 阶数=" << shape_function_order_ << ")" << std::endl;
}

// Setter方法实现
void Config::setCoefficientC(std::function<double(const std::vector<double>&)> func)
{
    coefficient_c_func = func;
}

void Config::setSourceTermF(std::function<double(const std::vector<double>&)> func)
{
    source_term_f_func = func;
}

void Config::setExactSolutionU(std::function<double(const std::vector<double>&)> func)
{
    exact_solution_u_func = func;
}

void Config::setExactSolutionGradients(
    std::function<double(const std::vector<double>&, std::vector<double>&)> func)
{
    exact_solution_gradients_func = func;
}

// 访问器实现
const std::string& Config::getMeshFilename() const
{
    return mesh_filename_;
}

void Config::setMeshFilename(const std::string& filename)
{
    mesh_filename_ = filename;
    // 重新加载网格文件
    loadMeshFromFile(mesh_filename_);
    // 重新设置高斯积分点
    setGaussAssemblePoints();
}

int Config::getDimension() const
{
    return dimension_;
}

int Config::getSamplingPointsNum() const
{
    return max_error_sampling_points_num_;
}

int Config::getErrorGaussPointsNum() const
{
    return gauss_error_points_num_;
}

int Config::getAssembleGaussPointsNum() const
{
    return gauss_assemble_points_num_;
}

int Config::getNodesNum() const
{
    return nodes_num_;
}

int Config::getElementsNum() const
{
    return elements_num_;
}

int Config::getNodesPerElement() const
{
    return nodes_num_per_element_;
}

int Config::getOrder() const
{
    return shape_function_order_;
}

int Config::getElementType() const
{
    return element_type_;
}

const std::vector<double>& Config::getNodeCoordinates() const
{
    return node_coordinates_;
}

const std::vector<std::vector<int>>& Config::getElementConnectivity() const
{
    return element_connectivity_;
}

bool Config::hasExactSolution() const
{
    return exact_solution_u_func != nullptr;
}

bool Config::hasExactGradients() const
{
    return exact_solution_gradients_func != nullptr;
}

// 求解器参数访问器实现
const std::string& Config::getSolverType() const
{
    return solver_type_;
}

const std::string& Config::getPreconditionerType() const
{
    return preconditioner_type_;
}

double Config::getSolverTolerance() const
{
    return solver_tolerance_;
}

int Config::getSolverMaxIterations() const
{
    return solver_max_iterations_;
}

// 调用函数实现
double Config::coefficient_c(const std::vector<double>& coords) const
{
    if (coefficient_c_func)
    {
        return coefficient_c_func(coords);
    }
    throw std::runtime_error("coefficient_c function not set");
}

double Config::source_term_f(const std::vector<double>& coords) const
{
    if (source_term_f_func)
    {
        return source_term_f_func(coords);
    }
    throw std::runtime_error("source_term_f function not set");
}

double Config::exact_solution_u(const std::vector<double>& coords) const
{
    if (exact_solution_u_func)
    {
        return exact_solution_u_func(coords);
    }
    throw std::runtime_error("exact_solution_u function not set");
}

double Config::exact_solution_gradients(const std::vector<double>& coords,
                                        std::vector<double>& gradients) const
{
    if (exact_solution_gradients_func)
    {
        return exact_solution_gradients_func(coords, gradients);
    }
    throw std::runtime_error("exact_solution_gradients function not set");
}

// 网格加载实现

void Config::loadMeshFromFile(const std::string& filename)
{
    // 使用COMSOL网格导入器导入网格
    ComsolMeshImporter importer;
    if (importer.importMesh(filename))
    {
        dimension_ = importer.getDimension();
        nodes_num_ = importer.getNodesNum();
        elements_num_ = importer.getElementsNum();
        nodes_num_per_element_ = importer.getNodesPerElement();
        element_connectivity_ = importer.getTriangularElements();

        // 转换节点坐标格式
        const auto& nodes = importer.getNodes();
        node_coordinates_.reserve(nodes.size() * dimension_);
        for (const auto& node : nodes)
        {
            node_coordinates_.push_back(node.first);   // x
            node_coordinates_.push_back(node.second);  // y
        }

        // TODO: 更完善的单元类型判断
        if (nodes_num_per_element_ == 3)
        {
            element_type_ = TRIANGLE;
        }
    }
    else
    {
        std::cerr << "在 Config 中加载网格失败: " << filename << std::endl;
        throw std::runtime_error("无法加载网格文件");
    }
}