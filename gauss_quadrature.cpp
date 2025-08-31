#include "gauss_quadrature.h"

void getGaussPoints(int numPoints, std::vector<double>& points, std::vector<double>& weights) {
    points.resize(numPoints);
    weights.resize(numPoints);

    switch (numPoints) {
        case 1:
            points[0] = 0.0;
            weights[0] = 2.0;
            break;
        case 2:
            points[0] = -0.5773502691896257;
            points[1] = 0.5773502691896257;
            weights[0] = 1.0;
            weights[1] = 1.0;
            break;
        case 3:
            points[0] = -0.7745966692414834;
            points[1] = 0.0;
            points[2] = 0.7745966692414834;
            weights[0] = 0.5555555555555556;
            weights[1] = 0.8888888888888888;
            weights[2] = 0.5555555555555556;
            break;
        case 4:
            points[0] = -0.8611363115940526;
            points[1] = -0.3399810435848563;
            points[2] = 0.3399810435848563;
            points[3] = 0.8611363115940526;
            weights[0] = 0.3478548451374538;
            weights[1] = 0.6521451548625461;
            weights[2] = 0.6521451548625461;
            weights[3] = 0.3478548451374538;
            break;
        case 5:
            points[0] = -0.9061798459386640;
            points[1] = -0.5384693101056831;
            points[2] = 0.0;
            points[3] = 0.5384693101056831;
            points[4] = 0.9061798459386640;
            weights[0] = 0.2369268850561891;
            weights[1] = 0.4786286704993665;
            weights[2] = 0.5688888888888889;
            weights[3] = 0.4786286704993665;
            weights[4] = 0.2369268850561891;
            break;
        default:
            // 默认或错误处理：使用1点积分
            points.resize(1);
            weights.resize(1);
            points[0] = 0.0;
            weights[0] = 2.0;
            break;
    }
}
