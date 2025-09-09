#include "gauss_quadrature_2d.h"
#include <stdexcept>

void getGaussPointsTriangle(int numPoints, 
                           std::vector<double>& gaussPoints, 
                           std::vector<double>& gaussWeights) {
    switch (numPoints) {
        case 1:
            getGaussPointsTriangle1(gaussPoints, gaussWeights);
            break;
        case 3:
            getGaussPointsTriangle3(gaussPoints, gaussWeights);
            break;
        case 4:
            getGaussPointsTriangle4(gaussPoints, gaussWeights);
            break;
        default:
            throw std::invalid_argument("Unsupported number of Gauss points for triangle");
    }
}

void getGaussPointsTriangle1(std::vector<double>& gaussPoints, 
                            std::vector<double>& gaussWeights) {
    // 一点积分：三角形中心
    gaussPoints.clear();
    gaussWeights.clear();
    
    gaussPoints.push_back(1.0/3.0);  // x_ref
    gaussPoints.push_back(1.0/3.0);  // y_ref
    
    gaussWeights.push_back(0.5);     // 权重（对于单位面积参考三角形）
}

void getGaussPointsTriangle3(std::vector<double>& gaussPoints, 
                            std::vector<double>& gaussWeights) {
    // 三点积分：三个边的中点
    gaussPoints.clear();
    gaussWeights.clear();
    
    // 权重（每个点相等）
    double weight = 1.0/6.0;
    
    // 点1: 边(0,1)的中点
    gaussPoints.push_back(0.5);      // x_ref
    gaussPoints.push_back(0.0);      // y_ref
    gaussWeights.push_back(weight);
    
    // 点2: 边(1,2)的中点
    gaussPoints.push_back(0.5);      // x_ref
    gaussPoints.push_back(0.5);      // y_ref
    gaussWeights.push_back(weight);
    
    // 点3: 边(0,2)的中点
    gaussPoints.push_back(0.0);      // x_ref
    gaussPoints.push_back(0.5);      // y_ref
    gaussWeights.push_back(weight);
}

void getGaussPointsTriangle4(std::vector<double>& gaussPoints, 
                            std::vector<double>& gaussWeights) {
    // 四点积分：三角形中心 + 三个顶点附近
    gaussPoints.clear();
    gaussWeights.clear();
    
    // 点1: 三角形中心
    gaussPoints.push_back(1.0/3.0);  // x_ref
    gaussPoints.push_back(1.0/3.0);  // y_ref
    gaussWeights.push_back(-27.0/96.0);
    
    // 其他三个点（对称位置）
    double a = 0.6;
    double b = 0.2;
    double weight = 25.0/96.0;
    
    // 点2
    gaussPoints.push_back(a);        // x_ref
    gaussPoints.push_back(b);        // y_ref
    gaussWeights.push_back(weight);
    
    // 点3
    gaussPoints.push_back(b);        // x_ref
    gaussPoints.push_back(a);        // y_ref
    gaussWeights.push_back(weight);
    
    // 点4
    gaussPoints.push_back(b);        // x_ref
    gaussPoints.push_back(b);        // y_ref
    gaussWeights.push_back(weight);
}
