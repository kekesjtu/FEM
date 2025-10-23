#include "config.h"
#include <iostream>
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
    // n个高斯点能精确积分2n-1阶多项式，所以需要至少(2*shape_function_order_+1+1)/2 =
    // shape_function_order_+1个点 但实际中通常取稍多一些以保证精度

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
              << ", 单元类型=" << element_type_ << ", 阶数=" << shape_function_order_ << ")"
              << std::endl;
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

const std::vector<Config::Boundary>& Config::getBoundary() const
{
    return boundarys;
}

bool Config::hasExactSolution() const
{
    return exact_solution_u_func != nullptr;
}

bool Config::hasExactGradients() const
{
    return exact_solution_gradients_func != nullptr;
}

// 边界单元相关访问器实现
int Config::getBoundaryElementType() const
{
    // 对于2D问题，边界是一维线段
    // 对于3D问题（未来），边界是二维面
    if (dimension_ == 2)
    {
        return LINE;
    }
    else
    {
        throw std::runtime_error("Boundary element type for dimension > 2 not implemented");
    }
}

int Config::getBoundaryNodesPerElement() const
{
    // 对于线性边界单元
    if (shape_function_order_ == 1)
    {
        if (dimension_ == 2)
        {
            return 2;  // 线段有2个节点
        }
        else
        {
            throw std::runtime_error(
                "Boundary nodes per element for dimension > 2 not implemented");
        }
    }
    else
    {
        throw std::runtime_error("Higher order boundary elements not implemented");
    }
}

int Config::getBoundaryGaussPointsNum() const
{
    // 对于线性边界单元，使用2点或3点高斯积分
    // 根据形函数阶数选择合适的积分点数
    if (shape_function_order_ == 1)
    {
        return 2;  // 线性边界单元用2点高斯积分
    }
    else if (shape_function_order_ == 2)
    {
        return 3;  // 二次边界单元用3点高斯积分
    }
    else
    {
        throw std::runtime_error("Higher order boundary elements not implemented");
    }
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

        // 初始化边界边信息
        initializeBoundarys(importer.getEdgeElements());
    }
    else
    {
        std::cerr << "在 Config 中加载网格失败: " << filename << std::endl;
        throw std::runtime_error("无法加载网格文件");
    }
}

void Config::initializeBoundarys(const std::vector<std::vector<int>>& edge_elements)
{
    boundarys.clear();
    boundarys.reserve(edge_elements.size());

    for (const auto& edge : edge_elements)
    {
        Boundary boundary;

        // ============================================================
        // 用户可在此处自定义边界条件
        // ============================================================

        // 选项 1: 齐次 Dirichlet 边界条件（默认）
        boundary.bc = BoundaryCondition::Dirichlet(0.0);

        // 选项 2: 非齐次 Dirichlet 边界条件
        // boundary.bc = BoundaryCondition::Dirichlet(1.0);  // u = 1.0

        // 选项 3: Neumann 边界条件
        // boundary.bc = BoundaryCondition::Neumann(2.0);  // c*du/dn = 2.0

        // 选项 4: Robin 边界条件
        // boundary.bc = BoundaryCondition::Robin(2.0, 3.0);  // c*du/dn + 2.0*u = 3.0

        // ============================================================

        // 找到包含这条边的单元
        boundary.element_index = findElementContainingEdge(edge);

        // 设置边界边的全局节点索引
        boundary.global_node_indices_in_element = edge;

        boundarys.push_back(boundary);
    }

    std::cout << "初始化边界边信息: " << boundarys.size() << " 条边界边" << std::endl;
    std::cout << "提示：可在 config.cpp 的 initializeBoundarys() 函数中自定义边界条件类型"
              << std::endl;
}

int Config::findElementContainingEdge(const std::vector<int>& edge_nodes) const
{
    // 简化实现：返回第一个包含所有边界节点的单元
    // 实际应用中可能需要更复杂的逻辑
    for (int i = 0; i < elements_num_; ++i)
    {
        const auto& element = element_connectivity_[i];
        bool contains_all_nodes = true;

        for (int edge_node : edge_nodes)
        {
            bool found = false;
            for (int elem_node : element)
            {
                if (elem_node == edge_node)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                contains_all_nodes = false;
                break;
            }
        }

        if (contains_all_nodes)
        {
            return i;
        }
    }

    // 如果找不到，返回-1表示错误
    return -1;
}