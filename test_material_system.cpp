#include <iomanip>
#include <iostream>
#include "material_v2.h"


void test_constant_property()
{
    std::cout << "\n=== 测试常数属性 ===" << std::endl;

    MaterialProperty k_const(401.0);

    EvaluationContext ctx;
    ctx.T = 300.0;
    ctx.x = 0.5;
    ctx.y = 0.5;
    ctx.z = 0.5;

    double value = k_const.evaluate(ctx);
    std::cout << "常数属性值: " << value << std::endl;

    if (std::abs(value - 401.0) < 1e-10)
    {
        std::cout << "✓ 常数属性测试通过" << std::endl;
    }
    else
    {
        std::cout << "✗ 常数属性测试失败" << std::endl;
    }
}

void test_expression_property()
{
    std::cout << "\n=== 测试表达式属性 ===" << std::endl;

    // 创建温度相关的电导率: σ(T) = σ0 / (1 + α(T - T0))
    std::string formula = "sigma0 / (1 + alpha * (T - T0))";
    std::vector<std::string> variables = {"T"};
    std::map<std::string, double> parameters = {
        {"sigma0", 5.96e7}, {"alpha", 0.00393}, {"T0", 293.15}};

    MaterialProperty sigma_expr(formula, variables, parameters);

    // 在不同温度下测试
    std::vector<double> temperatures = {293.15, 323.15, 373.15, 423.15};

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n温度(K)   电导率(S/m)" << std::endl;
    std::cout << "------------------------" << std::endl;

    for (double T : temperatures)
    {
        EvaluationContext ctx;
        ctx.T = T;

        double sigma = sigma_expr.evaluate(ctx);
        std::cout << std::setw(7) << T << "   " << std::scientific << std::setprecision(3) << sigma
                  << std::endl;
    }

    // 验证在T0时σ = σ0
    EvaluationContext ctx;
    ctx.T = 293.15;
    double sigma_at_T0 = sigma_expr.evaluate(ctx);

    if (std::abs(sigma_at_T0 - 5.96e7) < 1e5)
    {
        std::cout << "\n✓ 表达式属性测试通过" << std::endl;
    }
    else
    {
        std::cout << "\n✗ 表达式属性测试失败" << std::endl;
    }
}

void test_json_loading()
{
    std::cout << "\n=== 测试JSON加载 ===" << std::endl;

    try
    {
        // 创建简单的JSON字符串
        std::string json_str = R"({
            "materials": {
                "test_material": {
                    "description": "测试材料",
                    "thermal_conductivity": {
                        "type": "constant",
                        "value": 100.0
                    },
                    "electrical_conductivity": {
                        "type": "expression",
                        "formula": "sigma0 * T / T0",
                        "variables": ["T"],
                        "parameters": {
                            "sigma0": 1e6,
                            "T0": 300.0
                        }
                    },
                    "density": {
                        "type": "constant",
                        "value": 1000.0
                    },
                    "specific_heat": {
                        "type": "constant",
                        "value": 500.0
                    }
                }
            },
            "domain_materials": {
                "1": "test_material"
            }
        })";

        MaterialLibrary lib;
        lib.loadFromJSONString(json_str);

        // 测试获取材料
        auto material = lib.getMaterial("test_material");
        std::cout << "材料名称: " << material->name << std::endl;
        std::cout << "描述: " << material->description << std::endl;

        // 测试属性计算
        EvaluationContext ctx;
        ctx.T = 300.0;

        double k = material->thermal_conductivity.evaluate(ctx);
        double sigma = material->electrical_conductivity.evaluate(ctx);
        double rho = material->density.evaluate(ctx);
        double cp = material->specific_heat.evaluate(ctx);

        std::cout << std::fixed << std::setprecision(1);
        std::cout << "热导率: " << k << " W/(m·K)" << std::endl;
        std::cout << "电导率: " << std::scientific << sigma << " S/m" << std::endl;
        std::cout << "密度: " << std::fixed << rho << " kg/m³" << std::endl;
        std::cout << "比热: " << cp << " J/(kg·K)" << std::endl;

        // 测试域映射
        auto domain_material = lib.getMaterialForDomain(1);
        if (domain_material->name == "test_material")
        {
            std::cout << "\n✓ JSON加载和域映射测试通过" << std::endl;
        }
        else
        {
            std::cout << "\n✗ 域映射测试失败" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cout << "✗ JSON加载测试失败: " << e.what() << std::endl;
    }
}

void test_json_file()
{
    std::cout << "\n=== 测试JSON文件加载 ===" << std::endl;

    try
    {
        MaterialLibrary lib;
        lib.loadFromJSON("materials_example.json");

        std::cout << "已加载材料列表:" << std::endl;
        for (const auto& name : lib.getMaterialNames())
        {
            std::cout << "  - " << name << std::endl;
        }

        // 测试铜材料的温度相关电导率
        auto copper = lib.getMaterial("copper");
        std::cout << "\n铜材料电导率随温度变化:" << std::endl;
        std::cout << "温度(K)   电导率(S/m)" << std::endl;
        std::cout << "------------------------" << std::endl;

        std::vector<double> temps = {273.15, 293.15, 323.15, 373.15, 423.15};
        for (double T : temps)
        {
            EvaluationContext ctx;
            ctx.T = T;
            double sigma = copper->electrical_conductivity.evaluate(ctx);
            std::cout << std::fixed << std::setprecision(2) << std::setw(7) << T << "   "
                      << std::scientific << std::setprecision(3) << sigma << std::endl;
        }

        // 测试硅材料
        auto silicon = lib.getMaterial("silicon");
        std::cout << "\n硅材料在不同温度下的属性:" << std::endl;

        EvaluationContext ctx;
        ctx.T = 300.0;
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "T = " << ctx.T << " K:" << std::endl;
        std::cout << "  热导率: " << silicon->thermal_conductivity.evaluate(ctx) << " W/(m·K)"
                  << std::endl;
        std::cout << "  电导率: " << std::scientific
                  << silicon->electrical_conductivity.evaluate(ctx) << " S/m" << std::endl;

        ctx.T = 400.0;
        std::cout << std::fixed << "T = " << ctx.T << " K:" << std::endl;
        std::cout << "  热导率: " << silicon->thermal_conductivity.evaluate(ctx) << " W/(m·K)"
                  << std::endl;
        std::cout << "  电导率: " << std::scientific
                  << silicon->electrical_conductivity.evaluate(ctx) << " S/m" << std::endl;

        std::cout << "\n✓ JSON文件加载测试通过" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << "✗ JSON文件加载测试失败: " << e.what() << std::endl;
    }
}

void test_evaluation_context()
{
    std::cout << "\n=== 测试求值上下文 ===" << std::endl;

    // 创建依赖多个变量的属性
    std::string formula = "k0 * (1 + alpha * (T - T0)) * (1 + 0.1 * x)";
    std::vector<std::string> variables = {"T", "x"};
    std::map<std::string, double> parameters = {{"k0", 100.0}, {"alpha", 0.001}, {"T0", 300.0}};

    MaterialProperty k_complex(formula, variables, parameters);

    EvaluationContext ctx;
    ctx.T = 350.0;
    ctx.x = 0.5;

    double k = k_complex.evaluate(ctx);

    std::cout << "复杂表达式测试:" << std::endl;
    std::cout << "公式: " << formula << std::endl;
    std::cout << "T = " << ctx.T << ", x = " << ctx.x << std::endl;
    std::cout << "计算结果: " << k << std::endl;

    // 手动验证: k = 100 * (1 + 0.001*(350-300)) * (1 + 0.1*0.5)
    //            = 100 * 1.05 * 1.05 = 110.25
    double expected = 100.0 * 1.05 * 1.05;

    if (std::abs(k - expected) < 1e-6)
    {
        std::cout << "✓ 复杂求值上下文测试通过" << std::endl;
    }
    else
    {
        std::cout << "✗ 复杂求值上下文测试失败 (期望: " << expected << ")" << std::endl;
    }
}

int main()
{
    std::cout << "======================================" << std::endl;
    std::cout << "材料系统测试程序" << std::endl;
    std::cout << "======================================" << std::endl;

    try
    {
        test_constant_property();
        test_expression_property();
        test_evaluation_context();
        test_json_loading();
        test_json_file();

        std::cout << "\n======================================" << std::endl;
        std::cout << "所有测试完成!" << std::endl;
        std::cout << "======================================" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\n异常: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
