#include "mesh_hierarchy.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

std::vector<double> Mesh::getElementNodes(int element_id) const
{
    if (element_id < 0 || element_id >= num_elements_)
    {
        throw std::out_of_range("单元ID超出范围");
    }

    // 根据维度确定坐标数量
    int coords_per_node = static_cast<int>(dimension_);
    const auto& connectivity = element_connectivity_[element_id];

    std::vector<double> element_coords(connectivity.size() * coords_per_node);

    for (size_t i = 0; i < connectivity.size(); ++i)
    {
        int node_idx = connectivity[i];
        if (node_idx < 0 || node_idx >= num_nodes_)
        {
            throw std::out_of_range("节点ID超出范围");
        }

        for (int coord = 0; coord < coords_per_node; ++coord)
        {
            element_coords[i * coords_per_node + coord] =
                node_coordinates_[node_idx * coords_per_node + coord];
        }
    }

    return element_coords;
}

// ================================
// Mesh2D 类实现
// ================================

void Mesh2D::printMeshInfo() const
{
    std::cout << "\n=== 二维网格信息 ===" << std::endl;
    std::cout << "节点数量: " << num_nodes_ << std::endl;
    std::cout << "单元数量: " << num_elements_ << std::endl;
    std::cout << "边界边数量: " << boundary_edges_.size() << std::endl;
    std::cout << "网格维度: " << static_cast<int>(dimension_) << "D" << std::endl;
    std::cout << "每个单元节点数: " << getNodesPerElement() << std::endl;
}

void Mesh2D::clear()
{
    num_nodes_ = 0;
    num_elements_ = 0;
    node_coordinates_.clear();
    element_connectivity_.clear();
    boundary_edges_.clear();
}

// ================================
// TriangleMesh2D 类实现
// ================================

bool TriangleMesh2D::isValid() const
{
    if (!Mesh2D::isValid())
    {
        return false;
    }

    // 检查连接表大小
    if (element_connectivity_.size() != static_cast<size_t>(num_elements_))
    {
        return false;
    }

    // 检查每个三角形的节点数
    for (const auto& element : element_connectivity_)
    {
        if (element.size() != NODES_PER_TRIANGLE)
        {
            return false;
        }

        // 检查节点ID是否有效
        for (int node_id : element)
        {
            if (node_id < 0 || node_id >= num_nodes_)
            {
                return false;
            }
        }
    }

    return true;
}

void TriangleMesh2D::computeBoundaryEdges()
{
    // 简化版边界边计算（实际实现可能需要更复杂的算法）
    // 这里只是为了保持接口完整性
    boundary_edges_.clear();

    // TODO: 实现边界边检测算法
    // 目前先创建一个空的边界边列表
    std::cout << "边界边计算完成，找到 " << boundary_edges_.size() << " 条边界边" << std::endl;
}

// ================================
// MeshFactory 类实现
// ================================

std::unique_ptr<Mesh> MeshFactory::createMesh(MeshType type)
{
   //先不需要实现
   return nullptr;
}

std::unique_ptr<Mesh> MeshFactory::createMeshFromFile(const std::string& filename)
{
    //模仿comsol_mesh_importer.cpp的逻辑进行导入，读取信息创建网格对象
    //首先进行实现
    return nullptr;
}