#include <cmath>
#include "gauss_quadrature.h"

// ============================================================================
// 具体线段高斯点类实现
// ============================================================================

// LineGauss2Point 实现
LineGauss2Point::LineGauss2Point() : LineGaussPoint(2)
{
    // 直接在构造函数中初始化,避免调用虚函数
    // 二点高斯积分：标准区间 [-1, 1]
    points_.resize(2);
    weights_.resize(2);

    // 高斯点位置
    double a = 1.0 / std::sqrt(3.0);  // 约 0.577350269

    // 点1: t = -a
    points_[0] = -a;
    weights_[0] = 1.0;

    // 点2: t = a
    points_[1] = a;
    weights_[1] = 1.0;
}

// LineGauss3Point 实现
LineGauss3Point::LineGauss3Point() : LineGaussPoint(3)
{
    // 直接在构造函数中初始化,避免调用虚函数
    // 三点高斯积分：标准区间 [-1, 1]
    points_.resize(3);
    weights_.resize(3);

    double a = std::sqrt(3.0 / 5.0);  // 约 0.774596669

    // 点1: t = -a
    points_[0] = -a;
    weights_[0] = 5.0 / 9.0;

    // 点2: t = 0
    points_[1] = 0.0;
    weights_[1] = 8.0 / 9.0;

    // 点3: t = a
    points_[2] = a;
    weights_[2] = 5.0 / 9.0;
}
