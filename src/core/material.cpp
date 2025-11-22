#include "material.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "logger.h"

// ============================================================================
// 辅助函数：安全地解析域ID
// ============================================================================

/**
 * @brief 安全地将字符串转换为域ID
 * @param domain_str 域ID字符串（必须是数字）
 * @return 域ID
 * @throw std::invalid_argument 如果字符串不是有效的数字
 */
static int parseDomainId(const std::string& domain_str)
{
    try
    {
        return std::stoi(domain_str);
    }
    catch (const std::exception& e)
    {
        throw std::invalid_argument("Invalid domain ID: '" + domain_str + "'. Expected integer.");
    }
}

// ============================================================================
// MaterialProperty 实现
// ============================================================================

MaterialProperty::MaterialProperty()
    : type_(Type::CONSTANT), constant_value_(1.0), parser_initialized_(false)
{
}

MaterialProperty::MaterialProperty(double value)
    : type_(Type::CONSTANT), constant_value_(value), parser_initialized_(false)
{
}

MaterialProperty::MaterialProperty(const std::string& formula,
                                   const std::vector<std::string>& variables,
                                   const std::map<std::string, double>& parameters)
    : type_(Type::EXPRESSION),
      constant_value_(0.0),
      formula_(formula),
      variables_(variables),
      parameters_(parameters),
      parser_initialized_(false)
{
}

MaterialProperty MaterialProperty::fromJSON(const json& j)
{
    if (!j.contains("type"))
    {
        throw std::runtime_error("Material property JSON must contain 'type' field");
    }

    std::string type_str = j["type"];

    if (type_str == "constant")
    {
        if (!j.contains("value"))
        {
            throw std::runtime_error("Constant material property must have 'value' field");
        }
        return MaterialProperty(j["value"].get<double>());
    }
    else if (type_str == "expression")
    {
        if (!j.contains("formula"))
        {
            throw std::runtime_error("Expression material property must have 'formula' field");
        }

        std::string formula = j["formula"];
        std::vector<std::string> variables;
        std::map<std::string, double> parameters;

        // 读取依赖的变量
        if (j.contains("variables"))
        {
            variables = j["variables"].get<std::vector<std::string>>();
        }

        // 读取参数
        if (j.contains("parameters"))
        {
            parameters = j["parameters"].get<std::map<std::string, double>>();
        }

        return MaterialProperty(formula, variables, parameters);
    }
    else
    {
        throw std::runtime_error("Unknown material property type: " + type_str);
    }
}

void MaterialProperty::initializeParser(const EvaluationContext& ctx) const
{
    if (parser_initialized_)
    {
        return;
    }

    parser_ = std::make_shared<mu::Parser>();

    // 定义公式中的参数(常数)
    for (const auto& [name, value] : parameters_)
    {
        parser_->DefineConst(name, value);
    }

    // 为每个变量创建存储空间,并定义到parser中
    for (const auto& var_name : variables_)
    {
        // 创建变量存储并初始化为0
        variable_values_[var_name] = 0.0;
        // 让parser使用这个存储的指针
        parser_->DefineVar(var_name, &variable_values_[var_name]);
    }

    // 设置表达式
    try
    {
        parser_->SetExpr(formula_);
        parser_initialized_ = true;
    }
    catch (mu::Parser::exception_type& e)
    {
        std::string msg = "Failed to parse material property formula: ";
        msg += e.GetMsg();
        throw std::runtime_error(msg);
    }
}

double MaterialProperty::evaluate(const EvaluationContext& ctx) const
{
    if (type_ == Type::CONSTANT)
    {
        return constant_value_;
    }

    // 表达式类型
    if (!parser_initialized_)
    {
        initializeParser(ctx);
    }

    // 更新变量值
    auto var_map = const_cast<EvaluationContext&>(ctx).getVariableMap();
    for (const auto& var_name : variables_)
    {
        if (var_map.find(var_name) != var_map.end())
        {
            variable_values_[var_name] = *var_map[var_name];
        }
        else
        {
            throw std::runtime_error("Unknown variable in material property formula: " + var_name);
        }
    }

    try
    {
        return parser_->Eval();
    }
    catch (mu::Parser::exception_type& e)
    {
        std::string msg = "Error evaluating material property: ";
        msg += e.GetMsg();
        throw std::runtime_error(msg);
    }
}

// ============================================================================
// Material 实现
// ============================================================================

Material::Material()
    : name("DefaultMaterial"),
      description("Default material with unit properties"),
      thermal_conductivity(1.0),
      electrical_conductivity(1.0),
      density(1.0),
      specific_heat(1.0)
{
}

Material::Material(const std::string& name)
    : name(name),
      description(""),
      thermal_conductivity(1.0),
      electrical_conductivity(1.0),
      density(1.0),
      specific_heat(1.0)
{
}

std::shared_ptr<Material> Material::fromJSON(const std::string& material_name, const json& j)
{
    auto material = std::make_shared<Material>(material_name);

    // 读取描述
    if (j.contains("description"))
    {
        material->description = j["description"];
    }

    // 读取各个属性
    if (j.contains("thermal_conductivity"))
    {
        material->thermal_conductivity = MaterialProperty::fromJSON(j["thermal_conductivity"]);
        LOG_DEBUG("    - 读取thermal_conductivity属性");
    }

    if (j.contains("electrical_conductivity"))
    {
        material->electrical_conductivity =
            MaterialProperty::fromJSON(j["electrical_conductivity"]);
        LOG_DEBUG("    - 读取electrical_conductivity属性");
    }

    if (j.contains("density"))
    {
        material->density = MaterialProperty::fromJSON(j["density"]);
        LOG_DEBUG("    - 读取density属性");
    }

    if (j.contains("specific_heat"))
    {
        material->specific_heat = MaterialProperty::fromJSON(j["specific_heat"]);
        LOG_DEBUG("    - 读取specific_heat属性");
    }

    return material;
}

json Material::toJSON() const
{
    json j;
    j["name"] = name;
    j["description"] = description;

    // 属性转换辅助函数
    auto property_to_json = [](const MaterialProperty& prop) -> json
    {
        json p;
        if (prop.getType() == MaterialProperty::Type::CONSTANT)
        {
            p["type"] = "constant";
            p["value"] = prop.getConstantValue();
        }
        else
        {
            p["type"] = "expression";
            p["formula"] = prop.getFormula();
            if (!prop.getVariables().empty())
            {
                p["variables"] = prop.getVariables();
            }
            if (!prop.getParameters().empty())
            {
                p["parameters"] = prop.getParameters();
            }
        }
        return p;
    };

    j["thermal_conductivity"] = property_to_json(thermal_conductivity);
    j["electrical_conductivity"] = property_to_json(electrical_conductivity);
    j["density"] = property_to_json(density);
    j["specific_heat"] = property_to_json(specific_heat);

    return j;
}

// ============================================================================
// MaterialLibrary 实现
// ============================================================================

MaterialLibrary::MaterialLibrary()
{
}

void MaterialLibrary::loadFromJSON(const std::string& json_file)
{
    std::ifstream file(json_file);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open material library file: " + json_file);
    }

    json j;
    try
    {
        file >> j;
    }
    catch (json::parse_error& e)
    {
        throw std::runtime_error("Failed to parse JSON file: " + std::string(e.what()));
    }

    loadFromJSONString(j.dump());
}

void MaterialLibrary::loadFromJSONString(const std::string& json_str)
{
    json j;
    try
    {
        j = json::parse(json_str);
    }
    catch (json::parse_error& e)
    {
        throw std::runtime_error("Failed to parse JSON string: " + std::string(e.what()));
    }

    // 加载材料
    if (j.contains("materials"))
    {
        for (auto& [material_name, material_data] : j["materials"].items())
        {
            auto material = Material::fromJSON(material_name, material_data);
            addMaterial(material_name, material);
            LOG_DEBUG("  加载材料: '" + material_name + "'");
        }
    }

    // 加载域到材料的映射
    if (j.contains("domain_materials"))
    {
        for (auto& [domain_id_str, material_name] : j["domain_materials"].items())
        {
            int domain_id = parseDomainId(domain_id_str);
            assignMaterialToDomain(domain_id, material_name);
        }
    }

    // 加载边界域到材料的映射
    if (j.contains("boundary_domain_materials"))
    {
        for (auto& [domain_id_str, material_name] : j["boundary_domain_materials"].items())
        {
            int domain_id = parseDomainId(domain_id_str);
            assignMaterialToBoundaryDomain(domain_id, material_name);
        }
    }
}

void MaterialLibrary::addMaterial(const std::string& name, std::shared_ptr<Material> material)
{
    materials_[name] = material;
}

std::shared_ptr<Material> MaterialLibrary::getMaterial(const std::string& name) const
{
    auto it = materials_.find(name);
    if (it == materials_.end())
    {
        throw std::runtime_error("Material not found: " + name);
    }
    return it->second;
}

void MaterialLibrary::assignMaterialToDomain(int domain_id, const std::string& material_name)
{
    // 检查材料是否存在
    if (materials_.find(material_name) == materials_.end())
    {
        throw std::runtime_error("Cannot assign unknown material: " + material_name);
    }
    domain_to_material_[domain_id] = material_name;
}

void MaterialLibrary::assignMaterialToBoundaryDomain(int boundary_domain_id,
                                                     const std::string& material_name)
{
    if (materials_.find(material_name) == materials_.end())
    {
        throw std::runtime_error("Cannot assign unknown material: " + material_name);
    }
    boundary_domain_to_material_[boundary_domain_id] = material_name;
}

std::shared_ptr<Material> MaterialLibrary::getMaterialForDomain(int domain_id) const
{
    auto it = domain_to_material_.find(domain_id);
    if (it == domain_to_material_.end())
    {
        throw std::runtime_error("No material assigned to domain: " + std::to_string(domain_id));
    }
    return getMaterial(it->second);
}

std::shared_ptr<Material> MaterialLibrary::getMaterialForBoundaryDomain(
    int boundary_domain_id) const
{
    auto it = boundary_domain_to_material_.find(boundary_domain_id);
    if (it == boundary_domain_to_material_.end())
    {
        throw std::runtime_error("No material assigned to boundary domain: " +
                                 std::to_string(boundary_domain_id));
    }
    return getMaterial(it->second);
}

std::vector<std::string> MaterialLibrary::getMaterialNames() const
{
    std::vector<std::string> names;
    for (const auto& [name, material] : materials_)
    {
        names.push_back(name);
    }
    return names;
}

void MaterialLibrary::saveToJSON(const std::string& json_file) const
{
    json j;

    // 保存材料
    json materials_json;
    for (const auto& [name, material] : materials_)
    {
        materials_json[name] = material->toJSON();
    }
    j["materials"] = materials_json;

    // 保存域映射
    json domain_materials_json;
    for (const auto& [domain_id, material_name] : domain_to_material_)
    {
        domain_materials_json[std::to_string(domain_id)] = material_name;
    }
    j["domain_materials"] = domain_materials_json;

    // 保存边界域映射
    json boundary_domain_materials_json;
    for (const auto& [domain_id, material_name] : boundary_domain_to_material_)
    {
        boundary_domain_materials_json[std::to_string(domain_id)] = material_name;
    }
    j["boundary_domain_materials"] = boundary_domain_materials_json;

    // 写入文件
    std::ofstream file(json_file);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file for writing: " + json_file);
    }
    file << j.dump(2);  // 2个空格缩进
}
