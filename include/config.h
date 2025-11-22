#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>

// 简单的配置类，用于统一管理与误差分析/积分相关的计算参数。
// 该头同时提供别名 ComputeConfig = Config，以便与旧代码兼容。
class Config
{
  public:
    enum ElementType
    {
        LINE,           // 线段单元（一维，用于边界）
        TRIANGLE,       // 三角形单元（二维）
        QUADRILATERAL,  // 四边形单元（二维）
        TETRAHEDRON,    // 四面体单元（三维）
        // 其他3D单元类型可以在这里添加
    };

    /**
     * @brief 边界几何信息结构体
     *
     * 只包含边界的几何信息,不包含物理边界条件
     * 物理边界条件由 ProblemSetup 管理,以支持多物理场耦合
     *
     * 边界单元在不同维度下的含义:
     * - 1D空间: 边界点 (0维, 1个节点)
     * - 2D空间: 边界边 (1维, 2个节点)
     * - 3D空间: 边界面 (2维, 3或4个节点,三角形或四边形)
     */
    struct Boundary
    {
        int element_index;  ///< 边界单元所属的体单元索引(用于查找相邻体单元信息)
        std::vector<int> global_node_indices_in_element;  ///< 边界单元的全局节点索引列表
    };

  private:
    // 网格相关
    std::string mesh_filename_;
    int dimension_;
    int nodes_num_;
    int elements_num_;
    int nodes_num_per_element_;
    ElementType element_type_;
    std::vector<double> node_coordinates_;  // 节点坐标矩阵 [node_id * dimension + coord_id]
    std::vector<std::vector<int>>
        element_connectivity_;        // 单元连接矩阵 [element_id][local_node_id] = global_node_id
    std::vector<Boundary> boundarys;  // 边界单元信息 (1D:点, 2D:边, 3D:面)

    // 几何实体编码 (新增)
    std::vector<int>
        boundary_geometric_entities_;  ///< 边界单元几何实体编码 [boundary_id] = entity_id
    std::vector<int> element_geometric_entities_;  ///< 体单元几何实体编码 [element_id] = entity_id

    // 矩阵组装相关
    int gauss_assemble_points_num_;

    // 形函数阶数
    int order_ = 1;  // 形函数阶数默认为线性，用户可直接在此修改

    // 误差分析相关
    int max_error_sampling_points_num_ = 5;  // 最大误差采样点数，用户可直接在此修改
    int gauss_error_points_num_ = 3;         // 高斯误差积分点数，用户可直接在此修改

    // 网格单位缩放因子
    double mesh_unit_scale_ = 1.0;  // 默认1.0（不缩放），设为0.001可将mm转换为m

    // 文件中将会规定dimension、node_coordinates、element_connectivity、
    // int nodes_num_,int elements_num_,int nodes_num_per_element_;
    void loadMeshFromFile(const std::string& filename);

    // 这里要根据维度和单元类型设置gauss_assemble_points_.n个高斯点能精确积分2n+1阶多项式
    void setGaussAssemblePoints();

    /**
     * @brief 初始化边界单元信息
     * @param boundary_elements 边界单元连接数组 (来自ComsolMeshImporter)
     *                          1D: 点(1节点), 2D: 边(2节点), 3D: 面(3或4节点)
     */
    void initializeBoundarys(const std::vector<std::vector<int>>& boundary_elements);

    // 查找包含指定边的单元
    int findElementContainingEdge(const std::vector<int>& edge_nodes) const;

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
    int getNodesNumPerElement() const;
    int getOrder() const;
    int getElementType() const;
    const std::vector<double>& getNodeCoordinates() const;
    const std::vector<std::vector<int>>& getElementConnectivity() const;
    const std::vector<Boundary>& getBoundary() const;

    // 几何实体编码访问器 (新增)
    /**
     * @brief 获取边界单元几何实体编码
     * @return 边界几何实体向量 [boundary_id] = entity_id
     */
    const std::vector<int>& getBoundaryGeometricEntities() const;

    /**
     * @brief 获取单个边界单元的几何实体ID
     * @param boundary_idx 边界单元索引
     * @return 几何实体ID
     */
    int getBoundaryGeometricEntity(int boundary_idx) const;

    /**
     * @brief 获取体单元几何实体编码
     * @return 体单元几何实体向量 [element_id] = entity_id
     */
    const std::vector<int>& getElementGeometricEntities() const;

    // 边界单元相关访问器
    int getBoundaryElementType() const;      // 返回边界单元类型（LINE）
    int getBoundaryNodesPerElement() const;  // 返回边界单元节点数（2）
    int getBoundaryGaussPointsNum() const;   // 返回边界积分使用的高斯点数

    // 网格单位缩放
    void setMeshUnitScale(double scale);  // 设置网格单位缩放因子（例如0.001将mm转为m）
    double getMeshUnitScale() const;
};

#endif  // CONFIG_H
