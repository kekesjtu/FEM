// 测试 muParser 库
#include <iostream>
#include <string>
#include "third_party/muparser/muParser.h"

int main()
{
    std::cout << "测试 muParser..." << std::endl;

    try
    {
        // 创建解析器
        mu::Parser parser;

        // 定义变量
        double T = 350.0;  // 温度 [K]
        double sigma0 = 5.96e7;
        double alpha = 0.00393;
        double T0 = 293.15;

        // 设置表达式
        std::string expression = "sigma0 / (1 + alpha * (T - T0))";

        // 注册变量
        parser.DefineVar("T", &T);
        parser.DefineVar("sigma0", &sigma0);
        parser.DefineVar("alpha", &alpha);
        parser.DefineVar("T0", &T0);

        // 设置表达式
        parser.SetExpr(expression);

        // 计算结果
        double result = parser.Eval();

        std::cout << "\n表达式: " << expression << std::endl;
        std::cout << "变量值:" << std::endl;
        std::cout << "  T = " << T << " K" << std::endl;
        std::cout << "  sigma0 = " << sigma0 << " S/m" << std::endl;
        std::cout << "  alpha = " << alpha << " 1/K" << std::endl;
        std::cout << "  T0 = " << T0 << " K" << std::endl;
        std::cout << "\n结果: σ(T) = " << result << " S/m" << std::endl;

        // 测试改变温度
        std::cout << "\n测试温度变化:" << std::endl;
        for (double temp : {300.0, 350.0, 400.0, 450.0})
        {
            T = temp;
            double sigma = parser.Eval();
            std::cout << "  T = " << temp << " K  =>  σ = " << sigma << " S/m" << std::endl;
        }

        std::cout << "\n✅ muParser测试成功!" << std::endl;
    }
    catch (mu::Parser::exception_type &e)
    {
        std::cout << "❌ 错误: " << e.GetMsg() << std::endl;
        return 1;
    }

    return 0;
}
