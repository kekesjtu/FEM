#include "shape_functions.h"

//以下均为局部参考函数，定义在【-1，1】单元上
double shapeFunction_trial(int alpha, double x_ref) {
    if (alpha == 0) return 0.5 * (1.0 - x_ref);
    if (alpha == 1) return 0.5 * (1.0 + x_ref);
    return 0.0; // 如果索引无效，返回0
}

double shapeFunctionDerivative_trial(int alpha, double x_ref) {
    if (alpha == 0) return -0.5;
    if (alpha == 1) return 0.5;
    return 0.0; // 如果索引无效，返回0
}

double shapeFunction_test(int beta, double x_ref) {
    // 伽辽金法：检验函数与试探函数相同
    return shapeFunction_trial(beta, x_ref);
}

double shapeFunctionDerivative_test(int beta, double x_ref) {
    // 伽辽金法：检验函数导数与试探函数导数相同
    return shapeFunctionDerivative_trial(beta, x_ref);
}
