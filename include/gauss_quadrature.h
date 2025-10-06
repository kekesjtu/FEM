#ifndef GAUSS_QUADRATURE_2D_H
#define GAUSS_QUADRATURE_2D_H

#include <memory>
#include <vector>

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

  protected:
    /**
     * @brief 初始化高斯点和权重 - 由派生类实现
     */
    virtual void initializeGaussPoints() = 0;
};

/**
 * @brief 三角形单元一点积分
 */
class TriangleGauss1Point : public TriangleGaussPoint
{
  public:
    TriangleGauss1Point();

  protected:
    void initializeGaussPoints() override;
};

/**
 * @brief 三角形单元三点积分
 */
class TriangleGauss3Point : public TriangleGaussPoint
{
  public:
    TriangleGauss3Point();

  protected:
    void initializeGaussPoints() override;
};

/**
 * @brief 三角形单元四点积分
 */
class TriangleGauss4Point : public TriangleGaussPoint
{
  public:
    TriangleGauss4Point();

  protected:
    void initializeGaussPoints() override;
};

/**
 * @brief 高斯点工厂类，根据单元类型和积分点数创建相应的高斯积分对象
 */
class GaussPointFactory
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
     * @brief 创建高斯积分对象 (通用接口)
     * @param elementType 单元类型
     * @param numPoints 积分点数量
     * @return 高斯积分对象的智能指针
     */
    static std::unique_ptr<GaussPoint> createGaussPoint(ElementType elementType, int numPoints);
};

#endif  // GAUSS_QUADRATURE_2D_H
