/**
 * @file improved_factory_solution.h
 * @brief 改进当前工厂模式设计，提供更优雅的派生类功能访问方式
 */

#ifndef IMPROVED_FACTORY_SOLUTION_H
#define IMPROVED_FACTORY_SOLUTION_H

#include <memory>
#include <stdexcept>
#include <vector>


/**
 * @brief 形函数基类 - 增强版本
 */
class ShapeFunction
{
  public:
    virtual ~ShapeFunction() = default;

    // 基本接口
    virtual int getNumNodes() const = 0;
    virtual double computeTrialFunction(int node_index,
                                        const std::vector<double>& coords) const = 0;
    virtual std::vector<double> computeTrialDerivatives(
        int node_index, const std::vector<double>& coords) const = 0;

    // 类型查询接口 - 避免直接使用dynamic_cast
    virtual bool isType(const std::string& type) const = 0;

    // 模板方法：安全的类型转换
    template <typename T>
    T* as()
    {
        return dynamic_cast<T*>(this);
    }

    template <typename T>
    const T* as() const
    {
        return dynamic_cast<const T*>(this);
    }

    // 模板方法：检查是否为指定类型
    template <typename T>
    bool is() const
    {
        return dynamic_cast<const T*>(this) != nullptr;
    }
};

/**
 * @brief 二维形函数基类
 */
class ShapeFunction2D : public ShapeFunction
{
  public:
    virtual ~ShapeFunction2D() = default;

    virtual int getDimension() const
    {
        return 2;
    }

    // 二维特有的接口
    virtual double computeTrialFunction2D(int node_index, double xi, double eta) const = 0;
    virtual double computeTrialDerivativeXi(int node_index, double xi, double eta) const = 0;
    virtual double computeTrialDerivativeEta(int node_index, double xi, double eta) const = 0;

    // 实现基类方法
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

    bool isType(const std::string& type) const override
    {
        return type == "ShapeFunction2D" || type == "ShapeFunction";
    }
};

/**
 * @brief 三角形形函数类
 */
class TriangleShapeFunction : public ShapeFunction2D
{
  public:
    int getNumNodes() const override
    {
        return 3;
    }

    double computeTrialFunction2D(int node_index, double xi, double eta) const override
    {
        switch (node_index)
        {
            case 0:
                return 1.0 - xi - eta;
            case 1:
                return xi;
            case 2:
                return eta;
            default:
                return 0.0;
        }
    }

    double computeTrialDerivativeXi(int node_index, double xi, double eta) const override
    {
        switch (node_index)
        {
            case 0:
                return -1.0;
            case 1:
                return 1.0;
            case 2:
                return 0.0;
            default:
                return 0.0;
        }
    }

    double computeTrialDerivativeEta(int node_index, double xi, double eta) const override
    {
        switch (node_index)
        {
            case 0:
                return -1.0;
            case 1:
                return 0.0;
            case 2:
                return 1.0;
            default:
                return 0.0;
        }
    }

    bool isType(const std::string& type) const override
    {
        return type == "TriangleShapeFunction" || ShapeFunction2D::isType(type);
    }

    // 三角形特有的方法
    double computeArea(const std::vector<double>& coords) const
    {
        // coords 包含6个元素：x1,y1,x2,y2,x3,y3
        double x1 = coords[0], y1 = coords[1];
        double x2 = coords[2], y2 = coords[3];
        double x3 = coords[4], y3 = coords[5];
        return 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    }
};

/**
 * @brief 工厂类 - 支持多种返回类型
 */
class ShapeFunctionFactory
{
  public:
    enum class ElementType
    {
        Triangle,
        Quadrilateral
    };
    enum class Order
    {
        Linear,
        Quadratic
    };

    // 基础工厂方法 - 返回基类指针
    static std::unique_ptr<ShapeFunction> createShapeFunction(ElementType type,
                                                              Order order = Order::Linear)
    {
        switch (type)
        {
            case ElementType::Triangle:
                return createTriangleShapeFunction(order);
            default:
                throw std::invalid_argument("Unsupported element type");
        }
    }

    // 类型特化的工厂方法
    static std::unique_ptr<ShapeFunction2D> createShapeFunction2D(ElementType type,
                                                                  Order order = Order::Linear)
    {
        auto basePtr = createShapeFunction(type, order);
        auto derivedPtr =
            std::unique_ptr<ShapeFunction2D>(dynamic_cast<ShapeFunction2D*>(basePtr.release()));
        if (!derivedPtr)
        {
            throw std::runtime_error("Created shape function is not 2D");
        }
        return derivedPtr;
    }

    static std::unique_ptr<TriangleShapeFunction> createTriangleShapeFunction(
        Order order = Order::Linear)
    {
        switch (order)
        {
            case Order::Linear:
                return std::make_unique<TriangleShapeFunction>();
            default:
                throw std::invalid_argument("Unsupported order for triangle shape function");
        }
    }

    // 单例模式的获取方法（保持向后兼容）
    static ShapeFunction* getShapeFunction(ElementType type, Order order = Order::Linear)
    {
        static std::unordered_map<std::pair<ElementType, Order>, std::unique_ptr<ShapeFunction>,
                                  std::hash<std::underlying_type_t<ElementType>>>
            instances;

        auto key = std::make_pair(type, order);
        if (instances.find(key) == instances.end())
        {
            instances[key] = createShapeFunction(type, order);
        }
        return instances[key].get();
    }
};

/**
 * @brief 辅助函数 - 安全的类型特化调用
 */
template <typename T, typename Func>
auto withShapeFunction(ShapeFunction* basePtr, Func&& func) -> decltype(func(std::declval<T*>()))
{
    if (auto* derivedPtr = basePtr->as<T>())
    {
        return func(derivedPtr);
    }
    throw std::runtime_error("Shape function is not of the expected type");
}

/**
 * @brief 使用示例 - 展示各种使用方式
 */
class ImprovedFactoryExample
{
  public:
    void demonstrateUsage()
    {
        // 方式1：直接创建具体类型
        auto triangleFunc = ShapeFunctionFactory::createTriangleShapeFunction();
        double area = triangleFunc->computeArea({0, 0, 1, 0, 0, 1});

        // 方式2：基类指针 + 安全转换
        auto baseFunc =
            ShapeFunctionFactory::createShapeFunction(ShapeFunctionFactory::ElementType::Triangle);

        // 使用模板方法进行安全转换
        if (baseFunc->is<TriangleShapeFunction>())
        {
            auto* trianglePtr = baseFunc->as<TriangleShapeFunction>();
            double dxi = trianglePtr->computeTrialDerivativeXi(0, 0.2, 0.3);
        }

        // 方式3：使用辅助函数
        double result = withShapeFunction<TriangleShapeFunction>(
            baseFunc.get(), [](TriangleShapeFunction* triangle)
            { return triangle->computeTrialDerivativeXi(0, 0.2, 0.3); });

        // 方式4：类型查询
        if (baseFunc->isType("TriangleShapeFunction"))
        {
            // 安全地使用特定功能
        }
    }

    // 重构你的计算函数
    double calculateStiffnessEntry2D_improved(int e, int alpha, int beta, ShapeFunction* shapeFunc)
    {
        // 直接使用基类方法计算导数
        double xi = 0.2, eta = 0.3;  // 高斯点坐标
        std::vector<double> coords = {xi, eta};

        auto derivatives_alpha = shapeFunc->computeTrialDerivatives(alpha, coords);
        auto derivatives_beta = shapeFunc->computeTrialDerivatives(beta, coords);

        // 对于2D情况，derivatives[0]是dN/dxi, derivatives[1]是dN/deta
        return derivatives_alpha[0] * derivatives_beta[0] +
               derivatives_alpha[1] * derivatives_beta[1];
    }

    // 或者，如果需要特定功能
    double calculateStiffnessEntry2D_specialized(int e, int alpha, int beta,
                                                 ShapeFunction* shapeFunc)
    {
        return withShapeFunction<TriangleShapeFunction>(
            shapeFunc,
            [=](TriangleShapeFunction* triangle)
            {
                double xi = 0.2, eta = 0.3;
                double dN_alpha_xi = triangle->computeTrialDerivativeXi(alpha, xi, eta);
                double dN_alpha_eta = triangle->computeTrialDerivativeEta(alpha, xi, eta);
                double dN_beta_xi = triangle->computeTrialDerivativeXi(beta, xi, eta);
                double dN_beta_eta = triangle->computeTrialDerivativeEta(beta, xi, eta);
                return dN_alpha_xi * dN_beta_xi + dN_alpha_eta * dN_beta_eta;
            });
    }
};

#endif  // IMPROVED_FACTORY_SOLUTION_H