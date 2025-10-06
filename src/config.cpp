#include "config.h"
#include "comsol_mesh_importer.h"

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