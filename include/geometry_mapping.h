#ifndef GEOMETRY_MAPPING_H
#define GEOMETRY_MAPPING_H

#include <memory>
#include <stdexcept>  // 用于抛出异常
#include <vector>
#include "config.h"

/**
 * @brief 几何映射的抽象基类
 *
 * 定义了所有单元几何映射的通用接口。
 * 这个基类的设计考虑了雅可比矩阵可能在单元内变化的情况（例如高阶单元），
 * 因此核心函数都接受一个参考坐标作为参数。
 */
class GeometryMapping
{
  protected:
    std::vector<double> element_coords_;             // 单元节点坐标
    int nodes_num_per_element_;                      // 单元节点数量
    int order_;                                      // 单元的阶数
    int dimension_;                                  // 参考坐标维度 (1D/2D/3D)
    int embedding_dimension_;                        // 嵌入空间维度 (1D/2D/3D)
    double jacobian_det_;                            // 雅可比行列式 (常数)
    std::vector<std::vector<double>> jacobian_inv_;  // 逆雅可比矩阵 (常数)

  private:
    virtual void computeJacobian() = 0;

  public:
    /**
     * @brief 构造函数 - 通过节点坐标创建几何映射
     * @param element_coords 单元节点坐标 (e.g., 2D空间中的线段: [x0,y0,x1,y1])
     * @param nodes 节点数量
     * @param order 单元阶数
     * @param ref_dim 参考坐标维度 (线段=1, 三角形=2, 四面体=3)
     * @param embed_dim 嵌入空间维度 (坐标数组使用此维度验证)
     */
    GeometryMapping(const std::vector<double>& element_coords, int nodes, int order, int ref_dim,
                    int embed_dim)
        : element_coords_(element_coords),
          nodes_num_per_element_(nodes),
          order_(order),
          dimension_(ref_dim),
          embedding_dimension_(embed_dim)
    {
        // 坐标数组大小验证：节点数 × 嵌入维度
        if (element_coords_.size() != static_cast<size_t>(nodes * embed_dim))
        {
            throw std::invalid_argument(
                "Coordinate vector size does not match nodes and embedding dimension.");
        }
    }

    virtual ~GeometryMapping() = default;

    /**
     * @brief 从参考坐标映射到物理坐标
     * @param coord_ref 参考坐标 (e.g., 2D: [xi,eta])
     * @param coord_phys 输出的物理坐标 (e.g., 2D: [x,y])
     */
    virtual void mapToPhysical(const std::vector<double>& coord_ref,
                               std::vector<double>& coord_phys) const = 0;

    /**
     * @brief 在指定的参考点计算雅可比行列式
     * @param coord_ref 计算雅可比行列式的参考坐标点
     * @return double 该点的雅可比行列式
     */
    virtual double getJacobianDet(const std::vector<double>& coord_ref) const = 0;

    /**
     * @brief 在指定的参考点，将参考坐标系下的梯度转换为物理坐标系下的梯度
     * @param gradient_ref 参考坐标系下的梯度 (e.g., 2D: [dN/dxi, dN/deta])
     * @param gradient_phys 输出的物理坐标系下的梯度 (e.g., 2D: [dN/dx, dN/dy])
     * @param coord_ref 计算变换矩阵所在的参考坐标点
     */
    virtual void transformGradient(const std::vector<double>& gradient_ref,
                                   std::vector<double>& gradient_phys,
                                   const std::vector<double>& coord_ref) const = 0;
};

/**
 * @brief 一维几何映射的中间抽象基类
 *
 * 继承自通用的 GeometryMapping，并为所有 1D 单元（线段）添加通用接口，如计算长度。
 * 线段单元可以嵌入在1D、2D或3D空间中（例如2D问题的边界或3D问题的边）。
 */
class GeometryMapping1D : public GeometryMapping
{
  public:
    /**
     * @brief 构造函数
     * @param element_coords 单元节点坐标（嵌入在embed_dim维空间中）
     * @param order 插值阶数
     * @param nodes 节点数量
     * @param embed_dim 嵌入空间维度（1/2/3）
     */
    GeometryMapping1D(const std::vector<double>& element_coords, int order, int nodes,
                      int embed_dim)
        : GeometryMapping(element_coords, nodes, order,
                          1,          // ref_dim: 1D参考单元
                          embed_dim)  // embed_dim: 嵌入空间维度
    {
    }

    /**
     * @brief 获取单元的长度
     */
    virtual double getElementLength() const = 0;
};

/**
 * @brief 二维几何映射的中间抽象基类
 *
 * 继承自通用的 GeometryMapping，并为所有 2D 单元添加通用接口，如计算面积。
 * 二维单元可以嵌入在2D或3D空间中（例如3D问题中的曲面边界）。
 */
class GeometryMapping2D : public GeometryMapping
{
  public:
    /**
     * @brief 构造函数
     * @param element_coords 单元节点坐标
     * @param order 插值阶数
     * @param nodes 节点数量
     * @param embed_dim 嵌入空间维度（默认2，可以是3用于3D中的曲面）
     */
    GeometryMapping2D(const std::vector<double>& element_coords, int order, int nodes,
                      int embed_dim = 2)
        : GeometryMapping(element_coords, nodes, order,
                          2,          // ref_dim: 2D参考单元
                          embed_dim)  // embed_dim: 嵌入空间维度
    {
    }

    /**
     * @brief 获取单元的面积
     * @note 这通常通过对雅可比行列式在参考单元上积分得到。
     */
    virtual double getElementArea() const = 0;
};

/**
 * @brief 三维几何映射的中间抽象基类
 *
 * 继承自通用的 GeometryMapping，并为所有 3D 单元添加通用接口，如计算体积。
 * 三维单元通常嵌入在3D空间中。
 */
class GeometryMapping3D : public GeometryMapping
{
  public:
    /**
     * @brief 构造函数
     * @param element_coords 单元节点坐标
     * @param order 插值阶数
     * @param nodes 节点数量
     * @param embed_dim 嵌入空间维度（通常为3）
     */
    GeometryMapping3D(const std::vector<double>& element_coords, int order, int nodes,
                      int embed_dim = 3)
        : GeometryMapping(element_coords, nodes, order,
                          3,          // ref_dim: 3D参考单元
                          embed_dim)  // embed_dim: 嵌入空间维度
    {
    }

    /**
     * @brief 获取单元的体积
     * @note 这通常通过对雅可比行列式在参考单元上积分得到。
     */
    virtual double getElementVolume() const = 0;
};

// =========================================================================
// ==              1D 单元的具体实现 (Concrete Implementations)         ==
// =========================================================================

/**
 * @brief 一维线性线段单元的几何映射类 (Line Linear)
 *
 * 对于线性线段，雅可比（即边长的一半）在整个单元内是常数。
 * 线段可以嵌入在1D、2D或3D空间中。
 *
 * @par 参考坐标系
 * - 维度: 1D
 * - 参数范围: t ∈ [-1, 1]
 *
 * @par 物理坐标映射
 * x(t) = N0(t)*x0 + N1(t)*x1, 其中 N0=(1-t)/2, N1=(1+t)/2
 *
 * @par 雅可比
 * dx/dt 是长度为 embed_dim 的向量
 * jacobian_det = ||dx/dt|| = edge_length / 2
 */
class LineLinearMapping : public GeometryMapping1D
{
  private:
    void computeJacobian() override;

  public:
    /**
     * @brief 构造函数 - 创建1D线性线段映射
     * @param element_coords 两个端点的坐标
     *        1D: [x0, x1]
     *        2D: [x0, y0, x1, y1]
     *        3D: [x0, y0, z0, x1, y1, z1]
     * @param embed_dim 嵌入空间的维度（1/2/3）
     */
    LineLinearMapping(const std::vector<double>& element_coords, int embed_dim);

    void mapToPhysical(const std::vector<double>& coord_ref,
                       std::vector<double>& coord_phys) const override;

    double getJacobianDet(const std::vector<double>& coord_ref) const override;

    void transformGradient(const std::vector<double>& gradient_ref,
                           std::vector<double>& gradient_phys,
                           const std::vector<double>& coord_ref) const override;

    double getElementLength() const override;
};

// =========================================================================
// ==              2D 单元的具体实现 (Concrete Implementations)         ==
// =========================================================================

/**
 * @brief 二维线性三角形单元的几何映射类 (Triangle Linear)
 *
 * 对于线性三角形，雅可比矩阵在整个单元内是常数。
 * 三角形可以嵌入在2D或3D空间中（3D中的三角形可作为四面体的边界）。
 *
 * @par 参考坐标系
 * - 维度: 2D
 * - 参数范围: ξ,η ≥ 0, ξ+η ≤ 1
 *
 * @par 物理坐标映射
 * x(ξ,η) = N0*x0 + N1*x1 + N2*x2
 * 其中 N0=1-ξ-η, N1=ξ, N2=η
 */
class TriangleLinearMapping : public GeometryMapping2D
{
  private:
    void computeJacobian() override;

  public:
    /**
     * @brief 构造函数 - 创建2D线性三角形映射
     * @param element_coords 三个顶点的坐标
     *        2D: [x0, y0, x1, y1, x2, y2]
     *        3D: [x0, y0, z0, x1, y1, z1, x2, y2, z2]
     * @param embed_dim 嵌入空间维度（默认2）
     */
    TriangleLinearMapping(const std::vector<double>& element_coords, int embed_dim = 2)
        : GeometryMapping2D(element_coords, 1, 3, embed_dim)  // order=1, nodes=3
    {
        computeJacobian();  // 计算常数雅可比矩阵
    }

    void mapToPhysical(const std::vector<double>& coord_ref,
                       std::vector<double>& coord_phys) const override;

    double getJacobianDet(const std::vector<double>& coord_ref) const override;

    void transformGradient(const std::vector<double>& gradient_ref,
                           std::vector<double>& gradient_phys,
                           const std::vector<double>& coord_ref) const override;

    double getElementArea() const override;
};

// =========================================================================
// ==              3D 单元的具体实现 (Concrete Implementations)         ==
// =========================================================================

/**
 * @brief 三维线性四面体单元的几何映射类 (Tetrahedron Linear)
 *
 * 对于线性四面体，雅可比矩阵在整个单元内是常数。
 * 四面体嵌入在3D空间中。
 *
 * @par 参考坐标系
 * - 维度: 3D
 * - 参数范围: ξ,η,ζ ≥ 0, ξ+η+ζ ≤ 1
 *
 * @par 物理坐标映射
 * x(ξ,η,ζ) = N0*x0 + N1*x1 + N2*x2 + N3*x3
 * 其中 N0=1-ξ-η-ζ, N1=ξ, N2=η, N3=ζ
 *
 * @par 雅可比矩阵
 * J = [dx/dξ, dx/dη, dx/dζ] 是3×3矩阵
 * J = [x1-x0, x2-x0, x3-x0]
 */
class TetrahedronLinearMapping : public GeometryMapping3D
{
  private:
    void computeJacobian() override;

  public:
    /**
     * @brief 构造函数 - 创建3D线性四面体映射
     * @param element_coords 四个顶点的坐标 [x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3]
     * @param embed_dim 嵌入空间维度（默认3）
     */
    TetrahedronLinearMapping(const std::vector<double>& element_coords, int embed_dim = 3)
        : GeometryMapping3D(element_coords, 1, 4, embed_dim)  // order=1, nodes=4
    {
        computeJacobian();  // 计算常数雅可比矩阵
    }

    void mapToPhysical(const std::vector<double>& coord_ref,
                       std::vector<double>& coord_phys) const override;

    double getJacobianDet(const std::vector<double>& coord_ref) const override;

    void transformGradient(const std::vector<double>& gradient_ref,
                           std::vector<double>& gradient_phys,
                           const std::vector<double>& coord_ref) const override;

    double getElementVolume() const override;
};

/**
 * @brief 几何映射工厂类
 *
 * 根据单元类型和阶次，创建相应的几何映射对象。
 * 这种工厂模式将对象的创建逻辑与使用逻辑解耦，使代码更易于维护和扩展。
 */
class GeometryMappingFactory
{
  public:
    /**
     * @brief 统一的几何映射创建接口，所有信息从 Config 获取
     *
     * @param element_coords 单元节点的坐标
     * @param config_ 配置对象，包含单元类型、阶数和维度信息
     * @param is_boundary 是否创建边界单元的几何映射（默认 false，创建体单元）
     * @return 指向创建的几何映射对象的 unique_ptr
     */
    static std::unique_ptr<GeometryMapping> createMapping(const std::vector<double>& element_coords,
                                                          std::shared_ptr<Config> config_,
                                                          bool is_boundary = false);
};

#endif  // GEOMETRY_MAPPING_H