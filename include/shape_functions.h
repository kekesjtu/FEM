#ifndef SHAPE_FUNCTIONS_2D_H
#define SHAPE_FUNCTIONS_2D_H

#include <memory>
#include <vector>
#include "config.h"

/**
 * @brief 形函数基类 - 定义形函数的通用接口
 */
class ShapeFunction
{
  public:
    virtual ~ShapeFunction() = default;

    /**
     * @brief 计算试探函数值
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数值
     */
    virtual double computeTrialFunction(int element_node_index,
                                        const std::vector<double>& coords) const = 0;

    /**
     * @brief 计算检验函数值（伽辽金法中与试探函数相同）
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数值
     */
    virtual double computeTestFunction(int element_node_index,
                                       const std::vector<double>& coords) const = 0;

    /**
     * @brief 计算试探函数导数
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数导数向量
     */
    virtual std::vector<double> computeTrialGradients(int element_node_index,
                                                      const std::vector<double>& coords) const = 0;

    /**
     * @brief 计算检验函数导数（伽辽金法中与试探函数相同）
     * @param node_index 节点索引
     * @param coords 参考坐标
     * @return 形函数导数向量
     */
    virtual std::vector<double> computeTestGradients(int element_node_index,
                                                     const std::vector<double>& coords) const = 0;
};

/**
 * @brief 一维形函数基类
 */
class ShapeFunction1D : public ShapeFunction
{
  public:
    virtual ~ShapeFunction1D() = default;
};

/**
 * @brief 二维形函数基类
 */
class ShapeFunction2D : public ShapeFunction
{
  public:
    virtual ~ShapeFunction2D() = default;
};

/**
 * @brief 三维形函数基类
 */
class ShapeFunction3D : public ShapeFunction
{
  public:
    virtual ~ShapeFunction3D() = default;
};

/**
 * @brief 线性线段单元形函数类（一维）
 * 用于边界积分计算
 */
class LineLinearShapeFunction : public ShapeFunction1D
{
  public:
    LineLinearShapeFunction() = default;
    virtual ~LineLinearShapeFunction() = default;

    double computeTrialFunction(int element_node_index,
                                const std::vector<double>& coords) const override;
    double computeTestFunction(int element_node_index,
                               const std::vector<double>& coords) const override;
    std::vector<double> computeTrialGradients(int element_node_index,
                                              const std::vector<double>& coords) const override;
    std::vector<double> computeTestGradients(int element_node_index,
                                             const std::vector<double>& coords) const override;
};

/**
 * @brief 线性三角形单元形函数类
 */
class TriangleLinearShapeFunction : public ShapeFunction2D
{
  public:
    TriangleLinearShapeFunction() = default;
    virtual ~TriangleLinearShapeFunction() = default;

    double computeTrialFunction(int element_node_index,
                                const std::vector<double>& coords) const override;
    double computeTestFunction(int element_node_index,
                               const std::vector<double>& coords) const override;
    std::vector<double> computeTrialGradients(int element_node_index,
                                              const std::vector<double>& coords) const override;
    std::vector<double> computeTestGradients(int element_node_index,
                                             const std::vector<double>& coords) const override;
};

/**
 * @brief 线性四面体单元形函数类
 *
 * 对于标准参考四面体，四个顶点为：
 * 节点0: (0, 0, 0)
 * 节点1: (1, 0, 0)
 * 节点2: (0, 1, 0)
 * 节点3: (0, 0, 1)
 *
 * 形函数为：
 * N0(xi, eta, zeta) = 1 - xi - eta - zeta
 * N1(xi, eta, zeta) = xi
 * N2(xi, eta, zeta) = eta
 * N3(xi, eta, zeta) = zeta
 */
class TetrahedronLinearShapeFunction : public ShapeFunction3D
{
  public:
    TetrahedronLinearShapeFunction() = default;
    virtual ~TetrahedronLinearShapeFunction() = default;

    double computeTrialFunction(int element_node_index,
                                const std::vector<double>& coords) const override;
    double computeTestFunction(int element_node_index,
                               const std::vector<double>& coords) const override;
    std::vector<double> computeTrialGradients(int element_node_index,
                                              const std::vector<double>& coords) const override;
    std::vector<double> computeTestGradients(int element_node_index,
                                             const std::vector<double>& coords) const override;
};

/**
 * @brief 形函数工厂类，根据单元类型和阶数创建相应的形函数对象
 */
class ShapeFunctionFactory
{
  public:
    /**
     * @brief 统一的形函数创建接口，所有信息从 Config 获取
     * @param config_ 配置对象，包含单元类型和阶数信息
     * @param is_boundary 是否创建边界单元的形函数（默认 false，创建体单元）
     * @return 形函数对象的智能指针，自动返回对应的具体类型
     */
    static std::unique_ptr<ShapeFunction> createShapeFunction(std::shared_ptr<Config> config_,
                                                              bool is_boundary = false);
};

#endif  // SHAPE_FUNCTIONS_2D_H
