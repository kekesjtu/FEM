#include "shape_functions.h"

double shapeFunction_trial(int alpha, double xi) {
    if (alpha == 0) return 0.5 * (1.0 - xi);
    if (alpha == 1) return 0.5 * (1.0 + xi);
    return 0.0; // 如果索引无效，返回0
}

double shapeFunctionDerivative_trial(int alpha, double xi) {
    if (alpha == 0) return -0.5;
    if (alpha == 1) return 0.5;
    return 0.0; // 如果索引无效，返回0
}

double shapeFunction_test(int beta, double xi) {
    // 伽辽金法：检验函数与试探函数相同
    return shapeFunction_trial(beta, xi);
}

double shapeFunctionDerivative_test(int beta, double xi) {
    // 伽辽金法：检验函数导数与试探函数导数相同
    return shapeFunctionDerivative_trial(beta, xi);
}
