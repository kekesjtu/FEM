#ifndef COMSOL_MESH_IMPORTER_H
#define COMSOL_MESH_IMPORTER_H

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "problem_definition.h"

/**
 * @brief COMSOL网格导入器类
 *
 * 专门用于解析COMSOL Multiphysics生成的.mphtxt格式网格文件
 * 并将其转换为FEM求解器所需的P、T、boundary_edges格式
 */
class ComsolMeshImporter
{
  public:
    /**
     * @brief 构造函数
     */
    ComsolMeshImporter();

    /**
     * @brief 析构函数
     */
    ~ComsolMeshImporter();

    /**
     * @brief 导入COMSOL网格文件
     * @param filename 网格文件路径
     * @return bool 导入是否成功
     */
    bool importMesh(const std::string& filename);

    /**
     * @brief 填充全局数组P、T、boundary_edges
     * 必须在importMesh成功后调用
     * 注意：边界条件将设置为默认值，需要在defineProblem中重新设置
     */
    void populateGlobalArrays();

    /**
     * @brief 获取导入统计信息
     */
    void printImportStatistics() const;

  private:
    // 解析状态
    bool is_imported_;
    // 导入的文件名
    std::string filename_;

    // 网格数据（简化存储）
    std::vector<std::pair<double, double>> nodes_;       // 节点坐标，索引就是节点编号
    std::vector<std::vector<int>> triangular_elements_;  // 三角形单元连接
    std::vector<std::vector<int>> edge_elements_;        // 边界边连接

    // 解析方法
    bool parseFile(const std::string& filename);
    bool parseNodes(std::ifstream& file, std::string& line);
    bool parseElements(std::ifstream& file, std::string& line);
    bool parseElementType(std::ifstream& file, std::string& line);
    bool skipToNextSection(std::ifstream& file, std::string& line);

    // 数据处理方法
    int findElementContainingEdge(int node1, int node2) const;

    // 辅助方法
    void clearData();
    void trimString(std::string& str) const;
};

#endif  // COMSOL_MESH_IMPORTER_H
