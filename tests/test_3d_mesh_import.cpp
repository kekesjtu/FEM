/**
 * @file test_3d_mesh_import.cpp
 * @brief 测试3D四面体网格导入功能
 *
 * 测试ComsolMeshImporter是否能正确解析3D网格:
 * - 识别空间维度为3
 * - 解析节点坐标(3D)
 * - 识别四面体单元(tet, 4节点)
 * - 识别三角形边界单元(tri作为边界)
 * - 忽略低维单元(edg在3D中, vtx)
 */

#include <iostream>
#include <string>
#include <vector>
#include "../include/comsol_mesh_importer.h"

int main()
{
    std::cout << "================================================" << std::endl;
    std::cout << "     3D四面体网格导入测试" << std::endl;
    std::cout << "================================================\n" << std::endl;

    // 创建导入器
    ComsolMeshImporter importer;

    // 导入3D四面体网格
    std::string mesh_file = "cube_tet_mesh.mphtxt";
    std::cout << "开始导入COMSOL 3D网格文件: " << mesh_file << std::endl;

    if (!importer.importMesh(mesh_file))
    {
        std::cerr << "错误: 网格导入失败!" << std::endl;
        return 1;
    }

    std::cout << "\n✓ 网格导入成功!\n" << std::endl;

    // 获取网格数据
    int dimension = importer.getDimension();
    int nodes_num = importer.getNodesNum();
    const auto& node_coords = importer.getNodeCoordinates();
    const auto& elements = importer.getElementConnectivity();
    const auto& boundaries = importer.getBoundaryElements();

    // 显示统计信息
    std::cout << "=== 网格统计信息 ===" << std::endl;
    std::cout << "空间维度: " << dimension << "D" << std::endl;
    std::cout << "节点总数: " << nodes_num << std::endl;
    std::cout << "体单元数(四面体): " << elements.size() << std::endl;
    std::cout << "边界单元数(三角形): " << boundaries.size() << std::endl;

    if (!elements.empty())
    {
        std::cout << "每个体单元节点数: " << elements[0].size() << std::endl;
    }
    if (!boundaries.empty())
    {
        std::cout << "每个边界单元节点数: " << boundaries[0].size() << std::endl;
    }
    std::cout << "节点坐标数组大小: " << node_coords.size() << " (应为 " << nodes_num << " × "
              << dimension << " = " << nodes_num * dimension << ")" << std::endl;

    // 验证数据一致性
    std::cout << "\n=== 数据一致性检查 ===" << std::endl;
    bool all_passed = true;

    // 检查1: 维度应该是3
    if (dimension == 3)
    {
        std::cout << "✓ 维度检查通过: 3D" << std::endl;
    }
    else
    {
        std::cerr << "✗ 维度检查失败: 期望3D, 实际" << dimension << "D" << std::endl;
        all_passed = false;
    }

    // 检查2: 节点坐标数组大小
    if (node_coords.size() == static_cast<size_t>(nodes_num * dimension))
    {
        std::cout << "✓ 节点坐标数组大小正确" << std::endl;
    }
    else
    {
        std::cerr << "✗ 节点坐标数组大小错误" << std::endl;
        all_passed = false;
    }

    // 检查3: 四面体单元应该有4个节点
    if (!elements.empty() && elements[0].size() == 4)
    {
        std::cout << "✓ 四面体单元节点数正确(4)" << std::endl;
    }
    else if (!elements.empty())
    {
        std::cerr << "✗ 四面体单元节点数错误: " << elements[0].size() << std::endl;
        all_passed = false;
    }

    // 检查4: 三角形边界应该有3个节点
    if (!boundaries.empty() && boundaries[0].size() == 3)
    {
        std::cout << "✓ 三角形边界节点数正确(3)" << std::endl;
    }
    else if (!boundaries.empty())
    {
        std::cerr << "✗ 三角形边界节点数错误: " << boundaries[0].size() << std::endl;
        all_passed = false;
    }

    // 检查5: 应该有体单元
    if (elements.size() > 0)
    {
        std::cout << "✓ 成功解析体单元(四面体)" << std::endl;
    }
    else
    {
        std::cerr << "✗ 未解析到体单元" << std::endl;
        all_passed = false;
    }

    // 检查6: 应该有边界单元
    if (boundaries.size() > 0)
    {
        std::cout << "✓ 成功解析边界单元(三角形)" << std::endl;
    }
    else
    {
        std::cerr << "✗ 未解析到边界单元" << std::endl;
        all_passed = false;
    }

    // 显示部分数据样本
    std::cout << "\n=== 数据样本 ===" << std::endl;

    // 显示前3个节点的坐标
    std::cout << "前3个节点坐标(x, y, z):" << std::endl;
    for (int i = 0; i < std::min(3, nodes_num); ++i)
    {
        std::cout << "  节点" << i << ": (" << node_coords[i * dimension + 0] << ", "
                  << node_coords[i * dimension + 1] << ", " << node_coords[i * dimension + 2] << ")"
                  << std::endl;
    }

    // 显示前3个四面体单元
    std::cout << "\n前3个四面体单元(4个节点索引):" << std::endl;
    for (size_t i = 0; i < std::min(size_t(3), elements.size()); ++i)
    {
        std::cout << "  单元" << i << ": [";
        for (size_t j = 0; j < elements[i].size(); ++j)
        {
            std::cout << elements[i][j];
            if (j < elements[i].size() - 1)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }

    // 显示前3个三角形边界
    std::cout << "\n前3个三角形边界(3个节点索引):" << std::endl;
    for (size_t i = 0; i < std::min(size_t(3), boundaries.size()); ++i)
    {
        std::cout << "  边界" << i << ": [";
        for (size_t j = 0; j < boundaries[i].size(); ++j)
        {
            std::cout << boundaries[i][j];
            if (j < boundaries[i].size() - 1)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }

    // 最终结果
    std::cout << "\n================================================" << std::endl;
    if (all_passed)
    {
        std::cout << "✓✓✓ 所有测试通过! ✓✓✓" << std::endl;
        std::cout << "3D四面体网格导入功能正常" << std::endl;
    }
    else
    {
        std::cout << "✗✗✗ 部分测试失败 ✗✗✗" << std::endl;
    }
    std::cout << "================================================" << std::endl;

    return all_passed ? 0 : 1;
}
