#include "geometry_mapping.h"
#include "fem_solver.h"
#include "shape_functions.h"
#include <stdexcept>
#include <cmath>

// --- GeometryMapping1D 类实现 ---

GeometryMapping1D::GeometryMapping1D(int element_index) {
    // 通过单元索引获取节点坐标
    x_left = P[T[element_index][0]];   // 左端点
    x_right = P[T[element_index][1]];  // 右端点
    
    if (x_right <= x_left) {
        throw std::invalid_argument("右端点坐标必须大于左端点坐标");
    }
    
    // 计算雅可比行列式：dx/dx_ref = (x_right - x_left) / 2
    jacobian = (x_right - x_left) / 2.0;
}

GeometryMapping1D::GeometryMapping1D(const std::vector<double>& global_element_node) {
    if (global_element_node.size() != 2) {
        throw std::invalid_argument("线性单元必须有2个节点坐标");
    }
    
    x_left = global_element_node[0];
    x_right = global_element_node[1];
    
    if (x_right <= x_left) {
        throw std::invalid_argument("右端点坐标必须大于左端点坐标");
    }
    
    // 计算雅可比行列式：dx/dx_ref = (x_right - x_left) / 2
    jacobian = (x_right - x_left) / 2.0;
}

double GeometryMapping1D::mapToPhysical(double x_ref) const {
    // x = (x_left + x_right)/2 + (x_right - x_left)/2 * x_ref
    return 0.5 * (x_left + x_right) + 0.5 * (x_right - x_left) * x_ref;
}

double GeometryMapping1D::mapToReference(double x) const {
    // x_ref = 2*(x - (x_left + x_right)/2) / (x_right - x_left)
    return 2.0 * (x - 0.5 * (x_left + x_right)) / (x_right - x_left);
}

double GeometryMapping1D::getJacobian() const {
    return jacobian;
}

double GeometryMapping1D::getInverseJacobian() const {
    return 1.0 / jacobian;
}

double GeometryMapping1D::transformDerivative(double dN_dx_ref) const {
    // dN/dx = dN/dxi * dxi/dx = dN/dxi / jacobian
    return dN_dx_ref * getInverseJacobian();
}

double GeometryMapping1D::getElementLength() const {
    return x_right - x_left;
}
