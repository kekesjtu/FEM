#ifndef SHAPE_FUNCTIONS_2D_H
#define SHAPE_FUNCTIONS_2D_H

#include <memory>
#include <vector>

/**
 * @brief 形函数基类 - 定义形函数的通用接口
 */
class ShapeFunction
{
  public:
    virtual ~ShapeFunction() = default;

    /**
     * @brief 获取形函数的节点数量
     */
    virtual int getNumNodes() const = 0;

    /**
     * @brief 计算试探函数值
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数值
     */
    virtual double computeTrialFunction(int node_index,
                                        const std::vector<double>& coords) const = 0;

    /**
     * @brief 计算检验函数值（伽辽金法中与试探函数相同）
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数值
     */
    virtual double computeTestFunction(int node_index, const std::vector<double>& coords) const
    {
        return computeTrialFunction(node_index, coords);
    }

    /**
     * @brief 计算试探函数导数
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数导数向量
     */
    virtual std::vector<double> computeTrialDerivatives(
        int node_index, const std::vector<double>& coords) const = 0;

    /**
     * @brief 计算检验函数导数（伽辽金法中与试探函数相同）
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数导数向量
     */
    virtual std::vector<double> computeTestDerivatives(int node_index,
                                                       const std::vector<double>& coords) const
    {
        return computeTrialDerivatives(node_index, coords);
    }
};

/**
 * @brief 二维形函数基类
 */
class ShapeFunction2D : public ShapeFunction
{
  public:
    virtual ~ShapeFunction2D() = default;

    /**
     * @brief 获取参考坐标维度（二维）
     */
    virtual int getDimension() const
    {
        return 2;
    }

    /**
     * @brief 计算试探函数对第一个参考坐标的导数
     */
    virtual double computeTrialDerivativeXi(int node_index, double xi, double eta) const = 0;

    /**
     * @brief 计算试探函数对第二个参考坐标的导数
     */
    virtual double computeTrialDerivativeEta(int node_index, double xi, double eta) const = 0;

    /**
     * @brief 计算检验函数对第一个参考坐标的导数
     */
    virtual double computeTestDerivativeXi(int node_index, double xi, double eta) const
    {
        return computeTrialDerivativeXi(node_index, xi, eta);
    }

    /**
     * @brief 计算检验函数对第二个参考坐标的导数
     */
    virtual double computeTestDerivativeEta(int node_index, double xi, double eta) const
    {
        return computeTrialDerivativeEta(node_index, xi, eta);
    }

    // 重写基类的方法，提供二维坐标接口
    double computeTrialFunction(int node_index, const std::vector<double>& coords) const override
    {
        return computeTrialFunction2D(node_index, coords[0], coords[1]);
    }

    std::vector<double> computeTrialDerivatives(int node_index,
                                                const std::vector<double>& coords) const override
    {
        return {computeTrialDerivativeXi(node_index, coords[0], coords[1]),
                computeTrialDerivativeEta(node_index, coords[0], coords[1])};
    }

  protected:
    /**
     * @brief 二维形函数值计算
     */
    virtual double computeTrialFunction2D(int node_index, double xi, double eta) const = 0;
};

/**
 * @brief 线性三角形单元形函数类
 */
class TriangleShapeFunction : public ShapeFunction2D
{
  public:
    TriangleShapeFunction() = default;
    virtual ~TriangleShapeFunction() = default;

    int getNumNodes() const override
    {
        return 3;
    }

    double computeTrialDerivativeXi(int node_index, double xi, double eta) const override;
    double computeTrialDerivativeEta(int node_index, double xi, double eta) const override;

    // 为了向后兼容，将此方法设为public
    double computeTrialFunction2D(int node_index, double xi, double eta) const override;
};

/**
 * @brief 形函数工厂类，根据单元类型和阶数创建相应的形函数对象
 */
class ShapeFunctionFactory
{
  public:
    /**
     * @brief 单元类型枚举
     */
    enum class ElementType
    {
        Triangle,       // 三角形单元
        Quadrilateral,  // 四边形单元 (预留)
        Tetrahedron,    // 四面体单元 (预留)
        Hexahedron      // 六面体单元 (预留)
    };

    /**
     * @brief 形函数阶数枚举
     */
    enum class Order
    {
        Linear,     // 线性形函数
        Quadratic,  // 二次形函数 (预留)
        Cubic       // 三次形函数 (预留)
    };

    /**
     * @brief 统一的形函数创建接口 - 用户只需调用这一个函数
     * @param elementType 单元类型
     * @param order 形函数阶数 (默认线性)
     * @return 形函数对象的智能指针，自动返回对应的具体类型
     *
     * 使用示例：
     * auto shapeFunc = ShapeFunctionFactory::createShapeFunction(
     *     ShapeFunctionFactory::ElementType::Triangle,
     *     ShapeFunctionFactory::Order::Linear
     * );
     *
     * // 简化版本（使用默认线性阶数）：
     * auto shapeFunc = ShapeFunctionFactory::createShapeFunction(
     *     ShapeFunctionFactory::ElementType::Triangle
     * );
     */
    static std::unique_ptr<ShapeFunction> createShapeFunction(ElementType elementType,
                                                              Order order = Order::Linear);

    /**
     * @brief 获取共享的形函数实例 (单例模式) - 避免重复创建
     * @param elementType 单元类型
     * @param order 形函数阶数 (默认线性)
     * @return 形函数对象指针
     *
     * 使用示例：
     * auto* shapeFunc = ShapeFunctionFactory::getShapeFunction(
     *     ShapeFunctionFactory::ElementType::Triangle
     * );
     */
    static ShapeFunction* getShapeFunction(ElementType elementType, Order order = Order::Linear);

  private:
    // ================ 内部实现方法，用户无需直接调用 ================

    /**
     * @brief 创建2D形函数对象
     */
    static std::unique_ptr<ShapeFunction2D> createShapeFunction2D(ElementType elementType,
                                                                  Order order);

    /**
     * @brief 创建三角形形函数对象
     */
    static std::unique_ptr<TriangleShapeFunction> createTriangleShapeFunction(
        Order order = Order::Linear);

    /**
     * @brief 获取三角形形函数实例
     */
    static TriangleShapeFunction* getTriangleShapeFunction(Order order = Order::Linear);

    // 存储单例实例
    static std::unique_ptr<TriangleShapeFunction> triangleShapeFunction_;
};

#endif  // SHAPE_FUNCTIONS_2D_H
