/**
 * @file visitor_pattern_solution.h
 * @brief 使用访问者模式解决工厂模式中基类指针访问派生类功能的问题
 */

#ifndef VISITOR_PATTERN_SOLUTION_H
#define VISITOR_PATTERN_SOLUTION_H

#include <memory>
#include <vector>

// 前向声明访问者
class ShapeFunctionVisitor;

/**
 * @brief 形函数基类 - 支持访问者模式
 */
class ShapeFunction
{
  public:
    virtual ~ShapeFunction() = default;

    // 基本接口保持不变
    virtual int getNumNodes() const = 0;
    virtual double computeTrialFunction(int node_index,
                                        const std::vector<double>& coords) const = 0;

    // 新增：接受访问者
    virtual void accept(ShapeFunctionVisitor& visitor) = 0;
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

    // 抽象的二维特有接口
    virtual double computeTrialFunction2D(int node_index, double xi, double eta) const = 0;
    virtual double computeTrialDerivativeXi(int node_index, double xi, double eta) const = 0;
    virtual double computeTrialDerivativeEta(int node_index, double xi, double eta) const = 0;
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

    double computeTrialFunction(int node_index, const std::vector<double>& coords) const override
    {
        return computeTrialFunction2D(node_index, coords[0], coords[1]);
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

    void accept(ShapeFunctionVisitor& visitor) override;
};

/**
 * @brief 访问者基类 - 定义对不同形函数类型的操作
 */
class ShapeFunctionVisitor
{
  public:
    virtual ~ShapeFunctionVisitor() = default;
    virtual void visitTriangleShapeFunction(TriangleShapeFunction& triangle) = 0;
    // 可扩展：添加对其他形函数类型的访问
    // virtual void visitQuadShapeFunction(QuadShapeFunction& quad) = 0;
};

/**
 * @brief 具体访问者：计算刚度矩阵元素
 */
class StiffnessCalculationVisitor : public ShapeFunctionVisitor
{
  private:
    int alpha_, beta_;
    double xi_, eta_;
    double result_;

  public:
    StiffnessCalculationVisitor(int alpha, int beta, double xi, double eta)
        : alpha_(alpha), beta_(beta), xi_(xi), eta_(eta), result_(0.0)
    {
    }

    void visitTriangleShapeFunction(TriangleShapeFunction& triangle) override
    {
        // 使用三角形特有的方法
        double dN_alpha_xi = triangle.computeTrialDerivativeXi(alpha_, xi_, eta_);
        double dN_alpha_eta = triangle.computeTrialDerivativeEta(alpha_, xi_, eta_);
        double dN_beta_xi = triangle.computeTrialDerivativeXi(beta_, xi_, eta_);
        double dN_beta_eta = triangle.computeTrialDerivativeEta(beta_, xi_, eta_);

        result_ = dN_alpha_xi * dN_beta_xi + dN_alpha_eta * dN_beta_eta;
    }

    double getResult() const
    {
        return result_;
    }
};

// 实现accept方法
inline void TriangleShapeFunction::accept(ShapeFunctionVisitor& visitor)
{
    visitor.visitTriangleShapeFunction(*this);
}

/**
 * @brief 工厂类保持不变，返回基类指针
 */
class ShapeFunctionFactory
{
  public:
    enum class ElementType
    {
        Triangle,
        Quadrilateral
    };

    static std::unique_ptr<ShapeFunction> createShapeFunction(ElementType type)
    {
        switch (type)
        {
            case ElementType::Triangle:
                return std::make_unique<TriangleShapeFunction>();
            default:
                throw std::invalid_argument("Unsupported element type");
        }
    }
};

/**
 * @brief 使用示例
 */
class VisitorPatternExample
{
  public:
    void demonstrateUsage()
    {
        // 通过工厂创建，得到基类指针
        auto shapeFunction =
            ShapeFunctionFactory::createShapeFunction(ShapeFunctionFactory::ElementType::Triangle);

        // 创建访问者来执行特定操作
        StiffnessCalculationVisitor visitor(0, 1, 0.2, 0.3);

        // 通过访问者模式调用派生类特有功能
        shapeFunction->accept(visitor);

        double stiffnessValue = visitor.getResult();
        // 使用计算结果...
    }
};

#endif  // VISITOR_PATTERN_SOLUTION_H