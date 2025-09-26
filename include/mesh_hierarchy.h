// 网格层次结构设计方案

#ifndef MESH_HIERARCHY_H
#define MESH_HIERARCHY_H

#include <memory>
#include <string>
#include <vector>

/**
 * @brief 网格基类
 *
 * 提供所有网格类型的通用接口和数据
 */
class Mesh
{
  public:
    /**
     * @brief 网格维度枚举（可用 int 替代，枚举更安全）
     */

  protected:
    int num_nodes_;                         // 节点总数
    int num_elements_;                      // 单元总数
    int dimension_;                         // 网格维度
    std::vector<double> node_coordinates_;  // 节点坐标矩阵 [node_id * dimension + coord_id]
    std::vector<std::vector<int>>
        element_connectivity_;  // 单元连接矩阵 [element_id][local_node_id] = global_node_id

  public:
    explicit Mesh(int dim) : dimension_(dim), num_nodes_(0), num_elements_(0)
    {
    }
    virtual ~Mesh() = default;

    // 基本信息访问
    int getNumNodes() const
    {
        return num_nodes_;
    }
    int getNumElements() const
    {
        return num_elements_;
    }

    int getDimension() const
    {
        return dimension_;
    }

    // 坐标访问方法
    const std::vector<double>& getNodeCoordinates() const
    {
        return node_coordinates_;
    }

    // 连接矩阵访问方法
    const std::vector<std::vector<int>>& getElementConnectivity() const
    {
        return element_connectivity_;
    }

    // 获取单元节点坐标的便捷方法
    std::vector<double> getElementNodes(int element_id) const;

    // 设置方法（用于从P、T初始化mesh）
    void setNodes(int num_nodes, const std::vector<double>& coordinates)
    {
        num_nodes_ = num_nodes;
        node_coordinates_ = coordinates;
    }

    void setElements(int num_elements, const std::vector<std::vector<int>>& connectivity)
    {
        num_elements_ = num_elements;
        element_connectivity_ = connectivity;
    }

    // 纯虚函数 - 子类必须实现
    virtual void printMeshInfo() const = 0;
    virtual void clear() = 0;

    // 虚函数 - 子类可选择重写
    virtual bool isValid() const
    {
        return num_nodes_ > 0 && num_elements_ > 0;
    }
};

/**
 * @brief 二维网格类
 *
 * 专门处理二维网格的通用功能
 */
class Mesh2D : public Mesh
{
  public:
    /**
     * @brief 边界条件结构体
     */
    struct BoundaryEdge
    {
        int global_node_index1;
        int global_node_index2;
        // 可以包含边界条件信息...
    };

  protected:
    std::vector<BoundaryEdge> boundary_edges_;  // 边界边信息

  public:
    Mesh2D() : Mesh(2)  // 2D网格
    {
    }
    virtual ~Mesh2D() = default;

    // 边界相关
    int getNumBoundaryEdges() const
    {
        return static_cast<int>(boundary_edges_.size());
    }

    const std::vector<BoundaryEdge>& getBoundaryEdges() const
    {
        return boundary_edges_;
    }
    // 实现基类虚函数
    void printMeshInfo() const override;
    void clear() override;

    // 二维特定的虚函数
    virtual int getNodesPerElement() const = 0;
};

/**
 * @brief 三角形网格类
 *
 * 专门处理三角形网格
 */
class TriangleMesh2D : public Mesh2D
{
  public:
    static constexpr int NODES_PER_TRIANGLE = 3;

  public:
    TriangleMesh2D() = default;
    virtual ~TriangleMesh2D() = default;

    int getNodesPerElement() const override
    {
        return NODES_PER_TRIANGLE;
    }

    // 验证网格有效性
    bool isValid() const override;

  protected:
    // 内部辅助函数
    void computeBoundaryEdges();
};

/**
 * @brief 网格工厂类
 */
class MeshFactory
{
    /**
     * @brief 网格类型枚举
     */
    enum class MeshType
    {
        TRIANGLE_2D,
        QUADRILATERAL_2D,
        // 可以扩展其他类型...
    };

  public:
    /**
     * @brief 创建指定类型的网格
     */
    static std::unique_ptr<Mesh> createMesh(MeshType type);

    /**
     * @brief 从文件创建网格（自动检测类型）
     */
    static std::unique_ptr<Mesh> createMeshFromFile(const std::string& filename);
};

#endif  // MESH_HIERARCHY_H