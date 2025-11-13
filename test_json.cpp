// 测试 nlohmann/json 库
#include <iostream>
#include "third_party/nlohmann/json.hpp"

using json = nlohmann::json;

int main()
{
    std::cout << "测试 nlohmann/json..." << std::endl;

    // 创建JSON对象
    json material_config = {
        {"name", "copper"},
        {"thermal_conductivity", 401.0},
        {"electrical_conductivity",
         {{"type", "expression"},
          {"formula", "sigma0 / (1 + alpha * (T - T0))"},
          {"parameters", {{"sigma0", 5.96e7}, {"alpha", 0.00393}, {"T0", 293.15}}}}}};

    // 打印JSON
    std::cout << "\n材料配置:" << std::endl;
    std::cout << material_config.dump(2) << std::endl;

    // 读取数据
    std::string name = material_config["name"];
    double k = material_config["thermal_conductivity"];

    std::cout << "\n读取数据:" << std::endl;
    std::cout << "  材料名: " << name << std::endl;
    std::cout << "  热导率: " << k << " W/(m·K)" << std::endl;

    std::cout << "\n✅ JSON库测试成功!" << std::endl;
    return 0;
}
