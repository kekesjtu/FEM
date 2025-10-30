#include "gauss_quadrature.h"

// ============================================================================
// 具体四面体高斯点类实现
// ============================================================================

// TetrahedronGauss1Point 实现
TetrahedronGauss1Point::TetrahedronGauss1Point() : TetrahedronGaussPoint(1)
{
    // 直接在构造函数中初始化,避免调用虚函数
    // 一点积分：四面体中心
    // 标准参考四面体的中心点: (1/4, 1/4, 1/4)
    points_.resize(3);
    weights_.resize(1);

    points_[0] = 0.25;        // xi
    points_[1] = 0.25;        // eta
    points_[2] = 0.25;        // zeta
    weights_[0] = 1.0 / 6.0;  // 权重（对于单位体积参考四面体）
}

// TetrahedronGauss4Point 实现
TetrahedronGauss4Point::TetrahedronGauss4Point() : TetrahedronGaussPoint(4)
{
    // 直接在构造函数中初始化,避免调用虚函数
    // 四点积分：对称分布的四个点
    // 这是标准的四点积分公式，可精确积分二阶多项式
    points_.resize(12);
    weights_.resize(4);

    // 对称坐标值
    double a = 0.5854101966249685;  // (5 + sqrt(5))/20
    double b = 0.1381966011250105;  // (5 - sqrt(5))/20
    double weight = 1.0 / 24.0;     // 每个点的权重

    // 点1: (b, b, b)
    points_[0] = b;
    points_[1] = b;
    points_[2] = b;
    weights_[0] = weight;

    // 点2: (a, b, b)
    points_[3] = a;
    points_[4] = b;
    points_[5] = b;
    weights_[1] = weight;

    // 点3: (b, a, b)
    points_[6] = b;
    points_[7] = a;
    points_[8] = b;
    weights_[2] = weight;

    // 点4: (b, b, a)
    points_[9] = b;
    points_[10] = b;
    points_[11] = a;
    weights_[3] = weight;
}
