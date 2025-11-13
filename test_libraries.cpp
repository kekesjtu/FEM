// 测试 nlohmann/json 和 ExprTk 库是否可用
#include <iostream>
#include "third_party/exprtk/exprtk.hpp"
#include "third_party/nlohmann/json.hpp"


using json = nlohmann::json;

int main()
{
    // 测试 nlohmann/json
    std::cout << "测试 nlohmann/json..." << std::endl;
    json j = {{"name", "test"}, {"value", 42}, {"active", true}};
    std::cout << "  JSON: " << j.dump(2) << std::endl;

    // 测试 ExprTk
    std::cout << "\n测试 ExprTk..." << std::endl;
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    double x = 5.0;
    std::string expression_string = "3 * x + 10";

    symbol_table_t symbol_table;
    symbol_table.add_variable("x", x);

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    if (parser.compile(expression_string, expression))
    {
        double result = expression.value();
        std::cout << "  表达式: " << expression_string << std::endl;
        std::cout << "  x = " << x << std::endl;
        std::cout << "  结果: " << result << std::endl;
    }
    else
    {
        std::cout << "  表达式解析失败!" << std::endl;
    }

    std::cout << "\n✅ 所有库测试成功!" << std::endl;
    return 0;
}
