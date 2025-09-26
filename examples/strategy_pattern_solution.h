/**
 * @file strategy_pattern_solution.h
 * @brief 使用策略模式解决工厂模式中基类指针访问派生类功能的问题
 */

#ifndef STRATEGY_PATTERN_SOLUTION_H
#define STRATEGY_PATTERN_SOLUTION_H

#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>


/**
 * @brief 形函数操作策略基类
 */
class ShapeFunctionStrategy
{
  public:
    virtual ~ShapeFunctionStrategy() = default;

    // 定义各种操作的接口
    virtual double computeDerivativeXi(int node_index, double xi, double eta) const = 0;
    virtual double computeDerivativeEta(int node_index, double xi, double eta) const = 0;
    virtual double computeStiffnessContribution(int alpha, int beta, double xi,
                                                double eta) const = 0;
};

/**
 * @brief 三角形形函数策略实现
 */
class TriangleShapeFunctionStrategy : public ShapeFunctionStrategy
{
  public:
    double computeDerivativeXi(int node_index, double xi, double eta) const override
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

    double computeDerivativeEta(int node_index, double xi, double eta) const override
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

    double computeStiffnessContribution(int alpha, int beta, double xi, double eta) const override
    {
        double dN_alpha_xi = computeDerivativeXi(alpha, xi, eta);
        double dN_alpha_eta = computeDerivativeEta(alpha, xi, eta);
        double dN_beta_xi = computeDerivativeXi(beta, xi, eta);
        double dN_beta_eta = computeDerivativeEta(beta, xi, eta);

        return dN_alpha_xi * dN_beta_xi + dN_alpha_eta * dN_beta_eta;
    }
};

/**
 * @brief 改进的形函数基类 - 包含策略
 */
class ShapeFunction
{
  protected:
    std::unique_ptr<ShapeFunctionStrategy> strategy_;

  public:
    ShapeFunction(std::unique_ptr<ShapeFunctionStrategy> strategy) : strategy_(std::move(strategy))
    {
    }

    virtual ~ShapeFunction() = default;

    // 基本接口
    virtual int getNumNodes() const = 0;
    virtual double computeTrialFunction(int node_index,
                                        const std::vector<double>& coords) const = 0;

    // 通过策略提供派生类特有功能
    double computeDerivativeXi(int node_index, double xi, double eta) const
    {
        return strategy_->computeDerivativeXi(node_index, xi, eta);
    }

    double computeDerivativeEta(int node_index, double xi, double eta) const
    {
        return strategy_->computeDerivativeEta(node_index, xi, eta);
    }

    double computeStiffnessContribution(int alpha, int beta, double xi, double eta) const
    {
        return strategy_->computeStiffnessContribution(alpha, beta, xi, eta);
    }
};

/**
 * @brief 三角形形函数实现
 */
class TriangleShapeFunction : public ShapeFunction
{
  public:
    TriangleShapeFunction() : ShapeFunction(std::make_unique<TriangleShapeFunctionStrategy>())
    {
    }

    int getNumNodes() const override
    {
        return 3;
    }

    double computeTrialFunction(int node_index, const std::vector<double>& coords) const override
    {
        double xi = coords[0], eta = coords[1];
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
};

/**
 * @brief 工厂类 - 返回基类指针
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
class StrategyPatternExample
{
  public:
    void demonstrateUsage()
    {
        // 通过工厂创建，得到基类指针
        auto shapeFunction =
            ShapeFunctionFactory::createShapeFunction(ShapeFunctionFactory::ElementType::Triangle);

        // 直接通过基类接口使用派生类特有功能
        double dxi = shapeFunction->computeDerivativeXi(0, 0.2, 0.3);
        double deta = shapeFunction->computeDerivativeEta(0, 0.2, 0.3);
        double stiffness = shapeFunction->computeStiffnessContribution(0, 1, 0.2, 0.3);

        // 无需类型转换，直接使用
    }
};

#endif  // STRATEGY_PATTERN_SOLUTION_H