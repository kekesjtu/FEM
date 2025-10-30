#include "config.h"
#include <iostream>
#include "comsol_mesh_importer.h"

// 构造函数实现
Config::Config()
{
    loadMeshFromFile(mesh_filename_);  // 网格相关变量初始化完毕
    setGaussAssemblePoints();          // 矩阵组装计算参数初始化完毕
}

// 高斯积分点数设置
void Config::setGaussAssemblePoints()
{
    // 对于形函数阶数order_的单元，刚度矩阵积分涉及2*order_阶多项式
    // n个高斯点能精确积分2n-1阶多项式，所以需要至少(2*order_+1+1)/2 =
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
            if (order_ == 1)
            {
                gauss_assemble_points_num_ = 3;  // 线性单元用3点高斯积分
            }
            else if (order_ == 2)
            {
                gauss_assemble_points_num_ = 4;  // 二次单元用4点高斯积分
            }
            else
            {
                gauss_assemble_points_num_ = 7;  // 高阶单元用7点高斯积分
            }
        }
    }
    else if (dimension_ == 3)
    {
        if (element_type_ == TETRAHEDRON)
        {
            // 四面体单元：
            // 1点：精确积分1阶多项式（常数）
            // 4点：精确积分2阶多项式
            // 5点：精确积分3阶多项式
            // 11点：精确积分4阶多项式
            if (order_ == 1)
            {
                gauss_assemble_points_num_ = 4;  // 线性单元用4点高斯积分
            }
            else if (order_ == 2)
            {
                gauss_assemble_points_num_ = 5;  // 二次单元用5点高斯积分
            }
            else
            {
                gauss_assemble_points_num_ = 11;  // 高阶单元用11点高斯积分
            }
        }
    }
    else
    {
        // 一维情况或其他未实现的维度
        gauss_assemble_points_num_ = order_ + 1;
    }

    std::cout << "设置高斯积分点数: " << gauss_assemble_points_num_ << " (维度=" << dimension_
              << ", 单元类型=" << element_type_ << ", 阶数=" << order_ << ")" << std::endl;
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

int Config::getNodesNumPerElement() const
{
    return nodes_num_per_element_;
}

int Config::getOrder() const
{
    return order_;
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

// 边界单元相关访问器实现
int Config::getBoundaryElementType() const
{
    // 对于2D问题，边界是一维线段
    // 对于3D问题，边界是二维三角形（对于四面体网格）
    if (dimension_ == 2)
    {
        return LINE;
    }
    else if (dimension_ == 3)
    {
        // 对于四面体网格，边界是三角形
        if (element_type_ == TETRAHEDRON)
        {
            return TRIANGLE;
        }
        else
        {
            throw std::runtime_error("3D boundary element type only supports tetrahedron mesh");
        }
    }
    else
    {
        throw std::runtime_error("Unsupported dimension for boundary element type");
    }
}

int Config::getBoundaryNodesPerElement() const
{
    // 对于线性边界单元
    if (order_ == 1)
    {
        if (dimension_ == 2)
        {
            return 2;  // 2D边界：线段有2个节点
        }
        else if (dimension_ == 3)
        {
            return 3;  // 3D边界：三角形有3个节点
        }
        else
        {
            throw std::runtime_error("Boundary nodes per element for unsupported dimension");
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
    // 根据维度和形函数阶数选择合适的积分点数
    if (order_ == 1)
    {
        if (dimension_ == 2)
        {
            return 2;  // 2D边界（线段）：线性边界单元用2点高斯积分
        }
        else if (dimension_ == 3)
        {
            return 3;  // 3D边界（三角形）：线性边界单元用3点高斯积分
        }
        else
        {
            throw std::runtime_error("Unsupported dimension for boundary Gauss points");
        }
    }
    else if (order_ == 2)
    {
        if (dimension_ == 2)
        {
            return 3;  // 2D边界：二次边界单元用3点高斯积分
        }
        else if (dimension_ == 3)
        {
            return 4;  // 3D边界：二次边界单元用4点高斯积分
        }
        else
        {
            throw std::runtime_error("Unsupported dimension for boundary Gauss points");
        }
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

// 网格加载实现

void Config::loadMeshFromFile(const std::string& filename)
{
    // 使用COMSOL网格导入器导入网格
    ComsolMeshImporter importer;
    if (importer.importMesh(filename))
    {
        // 获取网格基本信息
        dimension_ = importer.getDimension();
        nodes_num_ = importer.getNodesNum();
        elements_num_ = importer.getElementsNum();
        nodes_num_per_element_ = importer.getNodesNumPerElement();

        // 获取网格数据(已经是Config需要的格式)
        node_coordinates_ = importer.getNodeCoordinates();
        element_connectivity_ = importer.getElementConnectivity();

        // 根据节点数判断单元类型
        if (dimension_ == 2 && nodes_num_per_element_ == 3)
        {
            element_type_ = TRIANGLE;
        }
        else if (dimension_ == 3 && nodes_num_per_element_ == 4)
        {
            element_type_ = TETRAHEDRON;
        }

        // 初始化边界信息
        initializeBoundarys(importer.getBoundaryElements());
    }
    else
    {
        std::cerr << "在 Config 中加载网格失败: " << filename << std::endl;
        throw std::runtime_error("无法加载网格文件");
    }
}

void Config::initializeBoundarys(const std::vector<std::vector<int>>& boundary_elements)
{
    boundarys.clear();
    boundarys.reserve(boundary_elements.size());

    for (const auto& boundary_element : boundary_elements)
    {
        Boundary boundary;

        // 找到包含这个边界单元的体单元
        // (1D: 包含边界点的线段, 2D: 包含边界边的三角形, 3D: 包含边界面的四面体)
        boundary.element_index = findElementContainingEdge(boundary_element);

        // 设置边界单元的全局节点索引
        boundary.global_node_indices_in_element = boundary_element;

        boundarys.push_back(boundary);
    }

    // 根据维度输出相应的边界单元类型
    std::string boundary_type;
    if (dimension_ == 1)
    {
        boundary_type = "个边界点";
    }
    else if (dimension_ == 2)
    {
        boundary_type = "条边界边";
    }
    else if (dimension_ == 3)
    {
        boundary_type = "个边界面";
    }
    else
    {
        boundary_type = "个边界单元";
    }

    std::cout << "初始化边界信息: " << boundarys.size() << " " << boundary_type << std::endl;
    std::cout << "提示：边界条件将由 ProblemSetup 在求解时动态设置" << std::endl;
}

int Config::findElementContainingEdge(const std::vector<int>& edge_nodes) const
{
    // 查找包含所有边界节点的体单元
    // 1D: 查找包含边界点的线段单元
    // 2D: 查找包含边界边两个节点的三角形单元
    // 3D: 查找包含边界面所有节点的四面体单元
    for (int i = 0; i < elements_num_; ++i)
    {
        const auto& element = element_connectivity_[i];
        bool contains_all_nodes = true;

        for (int boundary_node : edge_nodes)
        {
            bool found = false;
            for (int elem_node : element)
            {
                if (elem_node == boundary_node)
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