/**
 * @file type_erasure_solution.h
 * @brief 使用类型擦除技术解决工厂模式中基类指针访问派生类功能的问题
 */

#ifndef TYPE_ERASURE_SOLUTION_H
#define TYPE_ERASURE_SOLUTION_H

#include <functional>
#include <memory>
#include <vector>


/**
 * @brief 形函数概念接口（类型擦除实现）
 */
class ShapeFunction
{
  private:
    // 内部接口类
    struct ShapeFunctionConcept
    {
        virtual ~ShapeFunctionConcept() = default;
        virtual int getNumNodes() const = 0;
        virtual double computeTrialFunction(int node_index,
                                            const std::vector<double>& coords) const = 0;
        virtual double computeDerivativeXi(int node_index, double xi, double eta) const = 0;
        virtual double computeDerivativeEta(int node_index, double xi, double eta) const = 0;
        virtual double computeStiffnessContribution(int alpha, int beta, double xi,
                                                    double eta) const = 0;
        virtual std::unique_ptr<ShapeFunctionConcept> clone() const = 0;
    };

    // 模板实现类
    template <typename T>
    struct ShapeFunctionModel : public ShapeFunctionConcept
    {
        T data_;

        ShapeFunctionModel(T data) : data_(std::move(data))
        {
        }

        int getNumNodes() const override
        {
            return data_.getNumNodes();
        }

        double computeTrialFunction(int node_index,
                                    const std::vector<double>& coords) const override
        {
            return data_.computeTrialFunction(node_index, coords);
        }

        double computeDerivativeXi(int node_index, double xi, double eta) const override
        {
            // 使用 if constexpr (C++17) 或 SFINAE 检查方法是否存在
            if constexpr (requires { data_.computeDerivativeXi(node_index, xi, eta); })
            {
                return data_.computeDerivativeXi(node_index, xi, eta);
            }
            else
            {
                // 默认实现或抛出异常
                throw std::runtime_error(
                    "computeDerivativeXi not supported for this shape function type");
            }
        }

        double computeDerivativeEta(int node_index, double xi, double eta) const override
        {
            if constexpr (requires { data_.computeDerivativeEta(node_index, xi, eta); })
            {
                return data_.computeDerivativeEta(node_index, xi, eta);
            }
            else
            {
                throw std::runtime_error(
                    "computeDerivativeEta not supported for this shape function type");
            }
        }

        double computeStiffnessContribution(int alpha, int beta, double xi,
                                            double eta) const override
        {
            if constexpr (requires { data_.computeStiffnessContribution(alpha, beta, xi, eta); })
            {
                return data_.computeStiffnessContribution(alpha, beta, xi, eta);
            }
            else
            {
                // 默认实现
                double dN_alpha_xi = computeDerivativeXi(alpha, xi, eta);
                double dN_alpha_eta = computeDerivativeEta(alpha, xi, eta);
                double dN_beta_xi = computeDerivativeXi(beta, xi, eta);
                double dN_beta_eta = computeDerivativeEta(beta, xi, eta);
                return dN_alpha_xi * dN_beta_xi + dN_alpha_eta * dN_beta_eta;
            }
        }

        std::unique_ptr<ShapeFunctionConcept> clone() const override
        {
            return std::make_unique<ShapeFunctionModel>(data_);
        }
    };

    std::unique_ptr<ShapeFunctionConcept> pimpl_;

  public:
    // 构造函数模板 - 可以接受任何类型
    template <typename T>
    ShapeFunction(T data) : pimpl_(std::make_unique<ShapeFunctionModel<T>>(std::move(data)))
    {
    }

    // 拷贝构造
    ShapeFunction(const ShapeFunction& other) : pimpl_(other.pimpl_->clone())
    {
    }

    // 移动构造
    ShapeFunction(ShapeFunction&&) = default;

    // 赋值操作
    ShapeFunction& operator=(const ShapeFunction& other)
    {
        if (this != &other)
        {
            pimpl_ = other.pimpl_->clone();
        }
        return *this;
    }

    ShapeFunction& operator=(ShapeFunction&&) = default;

    // 公共接口 - 统一的接口，无需知道具体类型
    int getNumNodes() const
    {
        return pimpl_->getNumNodes();
    }

    double computeTrialFunction(int node_index, const std::vector<double>& coords) const
    {
        return pimpl_->computeTrialFunction(node_index, coords);
    }

    double computeDerivativeXi(int node_index, double xi, double eta) const
    {
        return pimpl_->computeDerivativeXi(node_index, xi, eta);
    }

    double computeDerivativeEta(int node_index, double xi, double eta) const
    {
        return pimpl_->computeDerivativeEta(node_index, xi, eta);
    }

    double computeStiffnessContribution(int alpha, int beta, double xi, double eta) const
    {
        return pimpl_->computeStiffnessContribution(alpha, beta, xi, eta);
    }
};

/**
 * @brief 具体的形函数实现类 - 不需要继承任何基类
 */
class TriangleShapeFunction
{
  public:
    int getNumNodes() const
    {
        return 3;
    }

    double computeTrialFunction(int node_index, const std::vector<double>& coords) const
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

    // 三角形特有的方法
    double computeDerivativeXi(int node_index, double xi, double eta) const
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

    double computeDerivativeEta(int node_index, double xi, double eta) const
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

    double computeStiffnessContribution(int alpha, int beta, double xi, double eta) const
    {
        double dN_alpha_xi = computeDerivativeXi(alpha, xi, eta);
        double dN_alpha_eta = computeDerivativeEta(alpha, xi, eta);
        double dN_beta_xi = computeDerivativeXi(beta, xi, eta);
        double dN_beta_eta = computeDerivativeEta(beta, xi, eta);

        return dN_alpha_xi * dN_beta_xi + dN_alpha_eta * dN_beta_eta;
    }
};

/**
 * @brief 工厂类 - 返回类型擦除后的统一对象
 */
class ShapeFunctionFactory
{
  public:
    enum class ElementType
    {
        Triangle,
        Quadrilateral
    };

    static ShapeFunction createShapeFunction(ElementType type)
    {
        switch (type)
        {
            case ElementType::Triangle:
                return ShapeFunction(TriangleShapeFunction{});
            default:
                throw std::invalid_argument("Unsupported element type");
        }
    }
};

/**
 * @brief 使用示例
 */
class TypeErasureExample
{
  public:
    void demonstrateUsage()
    {
        // 通过工厂创建，得到统一的ShapeFunction对象
        auto shapeFunction =
            ShapeFunctionFactory::createShapeFunction(ShapeFunctionFactory::ElementType::Triangle);

        // 直接调用所有方法，无需类型转换
        double value = shapeFunction.computeTrialFunction(0, {0.2, 0.3});
        double dxi = shapeFunction.computeDerivativeXi(0, 0.2, 0.3);
        double deta = shapeFunction.computeDerivativeEta(0, 0.2, 0.3);
        double stiffness = shapeFunction.computeStiffnessContribution(0, 1, 0.2, 0.3);

        // 完全不需要继承层次结构，但提供统一接口
    }
};

#endif  // TYPE_ERASURE_SOLUTION_H