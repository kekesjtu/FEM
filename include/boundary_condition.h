#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

/**
 * @file boundary_condition.h
 * @brief Boundary condition definition (decoupled from mesh geometry)
 *
 * This file defines the BoundaryCondition struct which is used by ProblemSetup
 * to specify physics-dependent boundary conditions, independently of mesh geometry.
 */

/**
 * @brief 边界条件定义
 *
 * 通用边界条件形式: K * c * du/dn + L * u = q
 * - K=0, L=1: Dirichlet 边界条件 (u = q)
 * - K=1, L=0: Neumann 边界条件 (c * du/dn = q)
 * - K=1, L≠0: Robin 边界条件 (c * du/dn + L * u = q)
 */
struct BoundaryCondition
{
    int K_bc;     ///< 导数项的系数 (0 or 1)
    double L_bc;  ///< 值项的系数
    double q_bc;  ///< 右侧项

    // ========================================
    // Helper constructor functions
    // ========================================

    /**
     * @brief 构造狄利克雷边界条件: u = val
     * @param val 边界值
     * @return Dirichlet 边界条件 (K=0, L=1, q=val)
     */
    static BoundaryCondition Dirichlet(double val)
    {
        return {0, 1.0, val};
    }

    /**
     * @brief 构造罗宾边界条件: c*du/dn + h*u = g
     * @param h 值项系数
     * @param g 右侧项值
     * @return Robin 边界条件 (K=1, L=h, q=g)
     */
    static BoundaryCondition Robin(double h, double g)
    {
        return {1, h, g};
    }

    /**
     * @brief 构造诺曼边界条件: c*du/dn = q_flux
     * @param q_flux 法向通量值
     * @return Neumann 边界条件 (K=1, L=0, q=q_flux)
     */
    static BoundaryCondition Neumann(double q_flux)
    {
        return {1, 0.0, q_flux};
    }
};

#endif  // BOUNDARY_CONDITION_H
