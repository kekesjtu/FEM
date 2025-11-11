#ifndef COMSOL_MESH_IMPORTER_H
#define COMSOL_MESH_IMPORTER_H

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "config.h"

/**
 * @brief COMSOL网格导入器类
 *
 * 专门用于解析COMSOL Multiphysics生成的.mphtxt格式网格文件
 * 支持1D/2D/3D网格的通用解析:
 * - 1D: 线段单元(体), 点(边界)
 * - 2D: 三角形/四边形(体), 边(边界)
 * - 3D: 四面体/六面体等(体), 三角形/四边形面(边界)
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
     * @brief 打印导入统计信息
     */
    void printImportStatistics() const;

    // --- 通用数据访问接口 (对应Config类成员) ---

    /**
     * @brief 获取空间维度
     * @return 空间维度 (1D/2D/3D)
     */
    int getDimension() const;

    /**
     * @brief 获取节点总数
     * @return 节点数量
     */
    int getNodesNum() const;

    /**
     * @brief 获取体单元总数
     * @return 单元数量
     */
    int getElementsNum() const;

    /**
     * @brief 获取每个体单元的节点数
     * @return 节点数量
     */
    int getNodesNumPerElement() const;

    /**
     * @brief 获取节点坐标 (展平数组)
     * @return 节点坐标向量 [node_id * dimension + coord_id]
     * @note 对应 Config::node_coordinates_
     */
    const std::vector<double>& getNodeCoordinates() const;

    /**
     * @brief 获取体单元连接关系
     * @return 单元连接矩阵 [element_id][local_node_id] = global_node_id
     * @note 对应 Config::element_connectivity_
     */
    const std::vector<std::vector<int>>& getElementConnectivity() const;

    /**
     * @brief 获取边界单元连接关系
     * @return 边界单元连接矩阵 [boundary_element_id][local_node_id] = global_node_id
     * @note 用于初始化 Config::boundarys
     *       - 1D: 边界点 (1个节点)
     *       - 2D: 边界边 (2个节点)
     *       - 3D: 边界面 (3或4个节点,三角形或四边形)
     */
    const std::vector<std::vector<int>>& getBoundaryElements() const;

    /**
     * @brief 获取边界单元几何实体编码
     * @return 边界几何实体向量 [boundary_id] = entity_id
     * @note 对应 Config::boundary_geometric_entities_
     */
    const std::vector<int>& getBoundaryGeometricEntities() const;

    /**
     * @brief 获取体单元几何实体编码
     * @return 体单元几何实体向量 [element_id] = entity_id
     * @note 对应 Config::element_geometric_entities_
     */
    const std::vector<int>& getElementGeometricEntities() const;

  private:
    // 单元分类枚举
    enum class ElementClassification
    {
        VOLUME,    ///< 体单元 (dimension维)
        BOUNDARY,  ///< 边界单元 (dimension-1维)
        IGNORED    ///< 忽略的单元 (dimension-2维或更低)
    };

    // 解析状态
    bool is_imported_;
    std::string filename_;

    // 网格基本信息 (对应 Config 成员)
    int dimension_;              ///< 空间维度 (1D/2D/3D)
    int nodes_num_;              ///< 节点总数
    int elements_num_;           ///< 体单元总数
    int nodes_num_per_element_;  ///< 每个体单元的节点数

    // 网格数据 (通用格式,对应 Config 成员)
    std::vector<double> node_coordinates_;  ///< 节点坐标 [node_id * dimension + coord_id]
    std::vector<std::vector<int>>
        element_connectivity_;  ///< 体单元连接 [element_id][local_node_id]
    std::vector<std::vector<int>>
        boundary_elements_;  ///< 边界单元连接 [boundary_id][local_node_id]

    // 几何实体编码 (新增)
    std::vector<int>
        boundary_geometric_entities_;  ///< 边界单元几何实体编码 [boundary_id] = entity_id
    std::vector<int> element_geometric_entities_;  ///< 体单元几何实体编码 [element_id] = entity_id

    // 主解析方法
    bool parseFile(const std::string& filename);

    // 顶层解析方法
    bool parseDimension(std::ifstream& file);
    bool parseNodes(std::ifstream& file);
    bool parseElements(std::ifstream& file);

    // 单元解析方法(并列关系)
    bool parseElementType(std::ifstream& file, std::string& element_type);
    bool parseNodesPerElement(std::ifstream& file, int& nodes_per_element);
    bool parseNumElements(std::ifstream& file, int& num_elements);
    bool parseElementConnectivity(std::ifstream& file, int num_elements, int nodes_per_element,
                                  const std::string& element_type);
    bool parseGeometricEntityIndices(std::ifstream& file, int num_elements,
                                     std::vector<int>& geometric_entities);

    // 单元分类判断模块
    ElementClassification classifyElement(const std::string& element_type,
                                          int nodes_per_element) const;
    // 辅助解析方法
    bool skipToNextSection(std::ifstream& file, std::string& line);
    void skipEmptyLines(std::ifstream& file, std::string& line);  // 工具方法
    void clearData();
    void trimString(std::string& str) const;
};

#endif  // COMSOL_MESH_IMPORTER_H
