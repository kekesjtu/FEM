#ifndef GEOMETRY_MAPPING_H
#define GEOMETRY_MAPPING_H

#include <stdexcept>  // 用于抛出异常
#include <vector>
#include <memory>

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
    std::vector<double> element_coords;  // 单元节点坐标
    int num_nodes;                       // 单元节点数量
    int dimension;                       // 空间维度 (2D/3D)

  public:
    /**
     * @brief 构造函数 - 通过节点坐标创建几何映射
     * @param coords 单元节点坐标 (e.g., 2D: [x0,y0,x1,y1,...])
     * @param nodes 节点数量
     * @param dim 空间维度
     */
    GeometryMapping(const std::vector<double>& coords, int nodes, int dim)
        : element_coords(coords), num_nodes(nodes), dimension(dim)
    {
        if (coords.size() != static_cast<size_t>(nodes * dim))
        {
            throw std::invalid_argument(
                "Coordinate vector size does not match nodes and dimension.");
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

    // 获取基本信息
    int getNumNodes() const
    {
        return num_nodes;
    }
    int getDimension() const
    {
        return dimension;
    }
};

/**
 * @brief 二维几何映射的中间抽象基类
 *
 * 继承自通用的 GeometryMapping，并为所有 2D 单元添加通用接口，如计算面积。
 */
class GeometryMapping2D : public GeometryMapping
{
  public:
    /**
     * @brief 构造函数
     * @param coords 单元节点坐标
     * @param nodes 节点数量
     */
    GeometryMapping2D(const std::vector<double>& coords, int nodes)
        : GeometryMapping(coords, nodes, 2)  // 维度固定为 2
    {
    }

    /**
     * @brief 获取单元的面积
     * @note 这通常通过对雅可比行列式在参考单元上积分得到。
     */
    virtual double getElementArea() const = 0;
};

// =========================================================================
// ==              2D 单元的具体实现 (Concrete Implementations)         ==
// =========================================================================

/**
 * @brief 二维线性三角形单元的几何映射类 (Triangle Linear)
 *
 * 对于线性三角形，雅可比矩阵在整个单元内是常数。
 */
class TriangleLinearMapping : public GeometryMapping2D
{
  private:
    double jacobian_det;                            // 雅可比行列式 (常数)
    std::vector<std::vector<double>> jacobian_inv;  // 逆雅可比矩阵 (常数)

    /**
     * @brief 计算常数雅可比矩阵及其相关量
     */
    void computeConstantJacobian();

  public:
    /**
     * @brief 构造函数 - 创建2D线性三角形映射
     * @param coords 三个顶点的坐标 [x0, y0, x1, y1, x2, y2]
     */
    TriangleLinearMapping(const std::vector<double>& coords);

    void mapToPhysical(const std::vector<double>& coord_ref,
                       std::vector<double>& coord_phys) const override;

    double getJacobianDet(const std::vector<double>& coord_ref) const override;

    void transformGradient(const std::vector<double>& gradient_ref,
                           std::vector<double>& gradient_phys,
                           const std::vector<double>& coord_ref) const override;

    double getElementArea() const override;
};

//-----------------------------------------------------------------------------
// 工厂类
// =========================================================================

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
     * @brief 单元类型枚举
     */
    enum class ElementType
    {
        Triangle,      // 三角形单元
        Quadrilateral  // 四边形单元
    };

    /**
     * @brief 创建几何映射对象
     *
     * @param type 单元的几何类型 (e.g., Triangle)
     * @param order 单元的插值阶次 (e.g., Linear)
     * @param coords 该单元所有节点的坐标
     * @return 指向创建的几何映射对象的 unique_ptr。返回基类指针以支持多态。
     */
    static std::unique_ptr<GeometryMapping> createMapping(ElementType type, int order,
                                                          const std::vector<double>& element_coords);
};

#endif  // GEOMETRY_MAPPING_H