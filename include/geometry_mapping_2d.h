#ifndef GEOMETRY_MAPPING_2D_H
#define GEOMETRY_MAPPING_2D_H

#include <vector>

/**
 * @brief 二维线性三角形单元的几何映射类
 * 处理参考三角形与物理三角形之间的映射
 */
class GeometryMapping2D
{
  private:
    std::vector<double> x_coords;                   // 三个顶点的x坐标 [x0, x1, x2]
    std::vector<double> y_coords;                   // 三个顶点的y坐标 [y0, y1, y2]
    double jacobian_det;                            // 雅可比行列式
    std::vector<std::vector<double>> jacobian_inv;  // 逆雅可比矩阵 2x2

    void computeJacobian();

  public:
    /**
     * @brief 构造函数 - 通过节点坐标创建几何映射
     * @param coords 单元三个顶点的坐标 [x0, y0, x1, y1, x2, y2]
     */
    GeometryMapping2D(const std::vector<double>& coords);

    /**
     * @brief 从参考坐标映射到物理坐标
     * @param x_ref 参考坐标第一分量
     * @param y_ref 参考坐标第二分量
     * @param x 输出物理坐标x
     * @param y 输出物理坐标y
     */
    void mapToPhysical(double x_ref, double y_ref, double& x, double& y) const;

    /**
     * @brief 从物理坐标映射到参考坐标
     * @param x 物理坐标x
     * @param y 物理坐标y
     * @param x_ref 输出参考坐标第一分量
     * @param y_ref 输出参考坐标第二分量
     */
    void mapToReference(double x, double y, double& x_ref, double& y_ref) const;

    /**
     * @brief 获取雅可比行列式
     * @return double 雅可比行列式
     */
    double getJacobianDet() const;

    /**
     * @brief 将参考坐标系下的梯度转换为物理坐标系下的梯度
     * @param dN_dx_ref 参考坐标系下对x_ref的偏导数
     * @param dN_dy_ref 参考坐标系下对y_ref的偏导数
     * @param dN_dx 输出物理坐标系下对x的偏导数
     * @param dN_dy 输出物理坐标系下对y的偏导数
     */
    void transformGradient(double dN_dx_ref, double dN_dy_ref, double& dN_dx, double& dN_dy) const;

    /**
     * @brief 获取单元面积
     * @return double 单元面积
     */
    double getElementArea() const;
};

#endif  // GEOMETRY_MAPPING_2D_H
