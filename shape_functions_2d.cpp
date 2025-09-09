#include "shape_functions_2d.h"

/**
 * @brief 计算二维线性三角形单元的试探函数值
 * 
 * 对于标准参考三角形，三个顶点为：
 * 节点0: (0, 0)
 * 节点1: (1, 0) 
 * 节点2: (0, 1)
 * 
 * 形函数为：
 * N0(x_ref, y_ref) = 1 - x_ref - y_ref
 * N1(x_ref, y_ref) = x_ref
 * N2(x_ref, y_ref) = y_ref
 */
double shapeFunction2D_trial(int alpha, double x_ref, double y_ref) {
    switch (alpha) {
        case 0:
            return 1.0 - x_ref - y_ref;
        case 1:
            return x_ref;
        case 2:
            return y_ref;
        default:
            return 0.0; // 错误情况
    }
}

/**
 * @brief 计算二维线性三角形单元的试探函数对x_ref的偏导数
 */
double shapeFunctionDerivativeXRef_trial(int alpha, double x_ref, double y_ref) {
    switch (alpha) {
        case 0:
            return -1.0;
        case 1:
            return 1.0;
        case 2:
            return 0.0;
        default:
            return 0.0; // 错误情况
    }
}

/**
 * @brief 计算二维线性三角形单元的试探函数对y_ref的偏导数
 */
double shapeFunctionDerivativeYRef_trial(int alpha, double x_ref, double y_ref) {
    switch (alpha) {
        case 0:
            return -1.0;
        case 1:
            return 0.0;
        case 2:
            return 1.0;
        default:
            return 0.0; // 错误情况
    }
}

/**
 * @brief 检验函数与试探函数相同（伽辽金法）
 */
double shapeFunction2D_test(int beta, double x_ref, double y_ref) {
    return shapeFunction2D_trial(beta, x_ref, y_ref);
}

double shapeFunctionDerivativeXRef_test(int beta, double x_ref, double y_ref) {
    return shapeFunctionDerivativeXRef_trial(beta, x_ref, y_ref);
}

double shapeFunctionDerivativeYRef_test(int beta, double x_ref, double y_ref) {
    return shapeFunctionDerivativeYRef_trial(beta, x_ref, y_ref);
}
