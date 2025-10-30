#ifndef GAUSS_QUADRATURE_2D_H
#define GAUSS_QUADRATURE_2D_H

#include <memory>
#include <vector>
#include "config.h"  // 使用 Config 的 ElementType

/**
 * @brief 高斯点基类，负责生成和存储高斯点信息
 * 定义所有高斯积分点的通用接口，并提供通用实现
 */
class GaussPoint
{
  protected:
    std::vector<double> points_;
    std::vector<double> weights_;
    int numPoints_;

  public:
    GaussPoint(int numPoints) : numPoints_(numPoints)
    {
    }
    virtual ~GaussPoint() = default;

    /**
     * @brief 获取积分点数量
     * @return 积分点数量
     */
    int getNumPoints() const
    {
        return numPoints_;
    }

    /**
     * @brief 获取积分点坐标向量 (格式由子类定义)
     * @return 积分点坐标向量
     */
    const std::vector<double>& getPoints() const
    {
        return points_;
    }

    /**
     * @brief 获取积分权重向量
     * @return 积分权重向量
     */
    const std::vector<double>& getWeights() const
    {
        return weights_;
    }

    /**
     * @brief 获取维度
     * @return 空间维度 (1D, 2D, 3D)
     */
    virtual int getDimension() const = 0;
};

/**
 * @brief 一维高斯点基类
 * 为一维高斯积分（线段）添加特定的接口
 */
class GaussPoint1D : public GaussPoint
{
  public:
    GaussPoint1D(int numPoints) : GaussPoint(numPoints)
    {
    }
    virtual ~GaussPoint1D() = default;

    /**
     * @brief 获取维度 (固定为1D)
     * @return 1
     */
    int getDimension() const override final
    {
        return 1;
    }
};

/**
 * @brief 二维高斯点基类
 * 为二维高斯积分添加特定的接口，并提供积分点存储
 */
class GaussPoint2D : public GaussPoint
{
  public:
    GaussPoint2D(int numPoints) : GaussPoint(numPoints)
    {
    }
    virtual ~GaussPoint2D() = default;

    /**
     * @brief 获取维度 (固定为2D)
     * @return 2
     */
    int getDimension() const override final
    {
        return 2;
    }
};

/**
 * @brief 三维高斯点基类
 * 为三维高斯积分添加特定的接口
 */
class GaussPoint3D : public GaussPoint
{
  public:
    GaussPoint3D(int numPoints) : GaussPoint(numPoints)
    {
    }
    virtual ~GaussPoint3D() = default;

    /**
     * @brief 获取维度 (固定为3D)
     * @return 3
     */
    int getDimension() const override final
    {
        return 3;
    }
};

/**
 * @brief 线段单元高斯点基类
 * 提供线段单元高斯积分的通用实现
 */
class LineGaussPoint : public GaussPoint1D
{
  public:
    LineGaussPoint(int numPoints) : GaussPoint1D(numPoints)
    {
    }
    virtual ~LineGaussPoint() = default;
};

/**
 * @brief 三角形单元高斯点基类
 * 提供三角形单元高斯积分的通用实现
 */
class TriangleGaussPoint : public GaussPoint2D
{
  public:
    TriangleGaussPoint(int numPoints) : GaussPoint2D(numPoints)
    {
    }
    virtual ~TriangleGaussPoint() = default;
};

/**
 * @brief 四面体单元高斯点基类
 * 提供四面体单元高斯积分的通用实现
 */
class TetrahedronGaussPoint : public GaussPoint3D
{
  public:
    TetrahedronGaussPoint(int numPoints) : GaussPoint3D(numPoints)
    {
    }
    virtual ~TetrahedronGaussPoint() = default;
};

/**
 * @brief 线段单元二点积分
 */
class LineGauss2Point : public LineGaussPoint
{
  public:
    LineGauss2Point();
};

/**
 * @brief 线段单元三点积分
 */
class LineGauss3Point : public LineGaussPoint
{
  public:
    LineGauss3Point();
};

/**
 * @brief 三角形单元一点积分
 */
class TriangleGauss1Point : public TriangleGaussPoint
{
  public:
    TriangleGauss1Point();
};

/**
 * @brief 三角形单元三点积分
 */
class TriangleGauss3Point : public TriangleGaussPoint
{
  public:
    TriangleGauss3Point();
};

/**
 * @brief 三角形单元四点积分
 */
class TriangleGauss4Point : public TriangleGaussPoint
{
  public:
    TriangleGauss4Point();
};

/**
 * @brief 四面体单元一点积分
 * 对应于四面体的中心点，精度为一阶
 */
class TetrahedronGauss1Point : public TetrahedronGaussPoint
{
  public:
    TetrahedronGauss1Point();
};

/**
 * @brief 四面体单元四点积分
 * 对应于四面体的四个顶点，精度为二阶
 */
class TetrahedronGauss4Point : public TetrahedronGaussPoint
{
  public:
    TetrahedronGauss4Point();
};

/**
 * @brief 高斯点工厂类，根据单元类型和积分点数创建相应的高斯积分对象
 */
class GaussPointFactory
{
  public:
    /**
     * @brief 根据单元类型和积分点数创建高斯积分对象
     * @param elementType 单元类型（使用 Config::ElementType）
     * @param numPoints 积分点数量
     * @return 高斯积分对象的智能指针
     */
    static std::unique_ptr<GaussPoint> createGaussPoint(Config::ElementType elementType,
                                                        int numPoints);
};

#endif  // GAUSS_QUADRATURE_2D_H
