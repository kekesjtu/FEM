#ifndef GEOMETRY_MAPPING_H
#define GEOMETRY_MAPPING_H

#include <vector>

/**
 * @brief 一维线性单元的几何映射类
 * 处理参考单元[-1,1]与物理单元[x_left, x_right]之间的映射
 */
class GeometryMapping1D {
private:
    double x_left;   // 物理单元左端点
    double x_right;  // 物理单元右端点
    double jacobian; // 雅可比行列式 dx/dx_ref
    
public:
    /**
     * @brief 构造函数 - 通过单元索引创建几何映射
     * @param element_index 单元索引
     */
    GeometryMapping1D(int element_index);
    
    /**
     * @brief 构造函数 - 通过节点坐标创建几何映射
     * @param global_element_node 单元节点坐标向量 [x_left, x_right]
     */
    GeometryMapping1D(const std::vector<double>& global_element_node);
    
    /**
     * @brief 从参考坐标映射到物理坐标
     * @param x_ref 参考坐标 [-1, 1]
     * @return double 物理坐标
     */
    double mapToPhysical(double x_ref) const;
    
    /**
     * @brief 从物理坐标映射到参考坐标
     * @param x 物理坐标
     * @return double 参考坐标 [-1, 1]
     */
    double mapToReference(double x) const;
    
    /**
     * @brief 获取雅可比行列式 dx/dxi
     * @return double 雅可比行列式
     */
    double getJacobian() const;
    
    /**
     * @brief 获取逆雅可比 dxi/dx
     * @return double 逆雅可比
     */
    double getInverseJacobian() const;
    
    /**
     * @brief 将参考坐标系下的导数转换为物理坐标系下的导数
     * @param dN_dx_ref 参考坐标系下的导数
     * @return double 物理坐标系下的导数
     */
    double transformDerivative(double dN_dx_ref) const;
    
    /**
     * @brief 获取单元长度
     * @return double 单元长度
     */
    double getElementLength() const;
};

#endif // GEOMETRY_MAPPING_H
