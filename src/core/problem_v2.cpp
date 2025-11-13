#include "problem_v2.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

// ============================================================================
// SolverConfiguration 实现
// ============================================================================

SolverConfiguration SolverConfiguration::fromJSON(const json& j)
{
    SolverConfiguration solver;

    if (j.contains("solver_type"))
    {
        solver.solver_type = j["solver_type"];
    }
    if (j.contains("preconditioner_type"))
    {
        solver.preconditioner_type = j["preconditioner_type"];
    }
    if (j.contains("tolerance"))
    {
        solver.tolerance = j["tolerance"];
    }
    if (j.contains("max_iterations"))
    {
        solver.max_iterations = j["max_iterations"];
    }

    return solver;
}

// ============================================================================
// 辅助函数：从JSON创建BoundaryCondition
// ============================================================================

static BoundaryCondition parseBoundaryConditionFromJSON(const json& j)
{
    if (!j.contains("type"))
    {
        throw std::runtime_error("Boundary condition must have 'type' field");
    }

    std::string type_str = j["type"];

    if (type_str == "Dirichlet")
    {
        double val = 0.0;
        if (j.contains("value"))
        {
            val = j["value"].get<double>();
        }
        return BoundaryCondition::Dirichlet(val);
    }
    else if (type_str == "Neumann")
    {
        double q_flux = 0.0;
        if (j.contains("value"))
        {
            q_flux = j["value"].get<double>();
        }
        return BoundaryCondition::Neumann(q_flux);
    }
    else if (type_str == "Robin")
    {
        double h = 1.0;
        double g = 0.0;
        if (j.contains("coefficient"))
        {
            h = j["coefficient"].get<double>();
        }
        if (j.contains("value"))
        {
            g = j["value"].get<double>();
        }
        return BoundaryCondition::Robin(h, g);
    }
    else
    {
        throw std::runtime_error("Unknown boundary condition type: " + type_str);
    }
}

// ============================================================================
// FieldConfiguration 实现
// ============================================================================

FieldConfiguration FieldConfiguration::fromJSON(const std::string& field_name, const json& j)
{
    FieldConfiguration field;
    field.name = field_name;

    // 解析源项
    if (j.contains("source"))
    {
        const auto& source_json = j["source"];

        if (source_json.contains("type") && source_json["type"] == "computed")
        {
            field.source_computed = true;
            field.source = MaterialProperty(0.0);  // 占位符
        }
        else
        {
            field.source_computed = false;
            if (source_json.is_number())
            {
                field.source = MaterialProperty(source_json.get<double>());
            }
            else if (source_json.is_object())
            {
                field.source = MaterialProperty::fromJSON(source_json);
            }
        }
    }

    // 解析边界条件
    if (j.contains("boundary_conditions"))
    {
        for (const auto& bc_item : j["boundary_conditions"])
        {
            if (!bc_item.contains("entities"))
            {
                throw std::runtime_error("Boundary condition must specify 'entities'");
            }

            BoundaryCondition bc = parseBoundaryConditionFromJSON(bc_item);

            // 处理实体列表
            const auto& entities = bc_item["entities"];
            if (entities.is_array())
            {
                for (const auto& entity_id : entities)
                {
                    field.entity_to_bc[entity_id.get<int>()] = bc;
                }
            }
            else if (entities.is_string() && entities.get<std::string>() == "all")
            {
                // "all" 标记会在后续处理
                // 这里先存储到一个特殊的实体ID（-1）
                field.entity_to_bc[-1] = bc;
            }
        }
    }

    // 解析精确解
    if (j.contains("exact_solution"))
    {
        field.has_exact_solution = true;
        const auto& exact = j["exact_solution"];

        if (exact.contains("u"))
        {
            field.exact_solution_u = MaterialProperty::fromJSON(exact["u"]);
        }
        if (exact.contains("grad_x"))
        {
            field.exact_solution_grad_x = MaterialProperty::fromJSON(exact["grad_x"]);
        }
        if (exact.contains("grad_y"))
        {
            field.exact_solution_grad_y = MaterialProperty::fromJSON(exact["grad_y"]);
        }
        if (exact.contains("grad_z"))
        {
            field.exact_solution_grad_z = MaterialProperty::fromJSON(exact["grad_z"]);
        }
    }

    return field;
}

// ============================================================================
// CouplingConfiguration 实现
// ============================================================================

CouplingConfiguration CouplingConfiguration::fromJSON(const json& j)
{
    CouplingConfiguration coupling;

    if (j.contains("type"))
    {
        coupling.type = j["type"];
    }
    if (j.contains("max_iterations"))
    {
        coupling.max_iterations = j["max_iterations"];
    }
    if (j.contains("tolerance"))
    {
        coupling.tolerance = j["tolerance"];
    }
    if (j.contains("relaxation_factor"))
    {
        coupling.relaxation_factor = j["relaxation_factor"];
    }
    if (j.contains("reference_temperature"))
    {
        coupling.reference_temperature = j["reference_temperature"];
    }

    return coupling;
}

// ============================================================================
// ProblemSetupV2 实现
// ============================================================================

ProblemSetupV2::ProblemSetupV2()
    : name("unnamed_problem"),
      description(""),
      dimension(2),
      mesh_file(""),
      material_library_file(""),
      has_coupling(false)
{
}

void ProblemSetupV2::loadFromJSON(const std::string& json_file)
{
    std::ifstream file(json_file);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open problem file: " + json_file);
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

    parseBasicInfo(j);
    parseDomainMaterials(j);
    parseFields(j);
    parseCoupling(j);
}

void ProblemSetupV2::loadFromJSONString(const std::string& json_str)
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

    parseBasicInfo(j);
    parseDomainMaterials(j);
    parseFields(j);
    parseSolver(j);
    parseCoupling(j);
}

void ProblemSetupV2::parseBasicInfo(const json& j)
{
    if (j.contains("problem_name"))
    {
        name = j["problem_name"];
    }
    if (j.contains("description"))
    {
        description = j["description"];
    }
    if (j.contains("dimension"))
    {
        dimension = j["dimension"];
    }
    if (j.contains("mesh_file"))
    {
        mesh_file = j["mesh_file"];
    }
    if (j.contains("material_library_file"))
    {
        material_library_file = j["material_library_file"];
    }
}

void ProblemSetupV2::parseDomainMaterials(const json& j)
{
    if (j.contains("domain_materials"))
    {
        for (auto& [domain_id_str, material_name] : j["domain_materials"].items())
        {
            int domain_id = std::stoi(domain_id_str);
            domain_to_material[domain_id] = material_name;
        }
    }
}

void ProblemSetupV2::parseFields(const json& j)
{
    if (j.contains("fields"))
    {
        for (auto& [field_name, field_json] : j["fields"].items())
        {
            fields[field_name] = FieldConfiguration::fromJSON(field_name, field_json);
        }
    }
}

void ProblemSetupV2::parseSolver(const json& j)
{
    if (j.contains("solver"))
    {
        solver = SolverConfiguration::fromJSON(j["solver"]);
    }
    // 如果没有solver字段，使用默认值
}

void ProblemSetupV2::parseCoupling(const json& j)
{
    if (j.contains("coupling"))
    {
        has_coupling = true;
        coupling = CouplingConfiguration::fromJSON(j["coupling"]);
    }
}

const FieldConfiguration& ProblemSetupV2::getField(const std::string& field_name) const
{
    auto it = fields.find(field_name);
    if (it == fields.end())
    {
        throw std::runtime_error("Field not found: " + field_name);
    }
    return it->second;
}

bool ProblemSetupV2::hasField(const std::string& field_name) const
{
    return fields.find(field_name) != fields.end();
}

std::string ProblemSetupV2::getMaterialForDomain(int domain_id) const
{
    auto it = domain_to_material.find(domain_id);
    if (it == domain_to_material.end())
    {
        throw std::runtime_error("No material assigned to domain: " + std::to_string(domain_id));
    }
    return it->second;
}

BoundaryCondition ProblemSetupV2::getBoundaryConditionForEntity(const std::string& field_name,
                                                                int entity_id) const
{
    const auto& field = getField(field_name);

    // 首先检查是否有"all"边界条件
    auto it_all = field.entity_to_bc.find(-1);

    // 然后查找具体的实体
    auto it = field.entity_to_bc.find(entity_id);

    if (it != field.entity_to_bc.end())
    {
        return it->second;
    }
    else if (it_all != field.entity_to_bc.end())
    {
        return it_all->second;
    }
    else
    {
        throw std::runtime_error("No boundary condition found for entity " +
                                 std::to_string(entity_id) + " in field '" + field_name + "'");
    }
}

void ProblemSetupV2::validate(bool warn_only) const
{
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    // 检查基本信息
    if (mesh_file.empty())
    {
        errors.push_back("Mesh file not specified");
    }

    if (material_library_file.empty())
    {
        warnings.push_back("Material library file not specified");
    }

    if (domain_to_material.empty())
    {
        errors.push_back("No domain materials assigned");
    }

    if (fields.empty())
    {
        errors.push_back("No fields defined");
    }

    // 检查每个场的配置
    for (const auto& [field_name, field] : fields)
    {
        if (!field.has_exact_solution)
        {
            warnings.push_back("Field '" + field_name +
                               "' has no exact solution (error analysis disabled)");
        }

        if (field.entity_to_bc.empty())
        {
            errors.push_back("Field '" + field_name + "' has no boundary conditions");
        }
    }

    // 输出警告
    for (const auto& warning : warnings)
    {
        std::cout << "⚠️  警告: " << warning << std::endl;
    }

    // 处理错误
    if (!errors.empty())
    {
        std::string error_msg = "Problem configuration validation failed:\n";
        for (const auto& error : errors)
        {
            error_msg += "  ❌ " + error + "\n";
        }

        if (warn_only)
        {
            std::cout << error_msg;
        }
        else
        {
            throw std::runtime_error(error_msg);
        }
    }
}

void ProblemSetupV2::saveToJSON(const std::string& json_file) const
{
    json j;

    // 基本信息
    j["problem_name"] = name;
    j["description"] = description;
    j["dimension"] = dimension;
    j["mesh_file"] = mesh_file;
    j["material_library_file"] = material_library_file;

    // 域材料映射
    json domain_materials_json;
    for (const auto& [domain_id, material_name] : domain_to_material)
    {
        domain_materials_json[std::to_string(domain_id)] = material_name;
    }
    j["domain_materials"] = domain_materials_json;

    // 字段配置（简化版本，完整版本需要更多实现）
    json fields_json;
    for (const auto& [field_name, field] : fields)
    {
        json field_json;
        // TODO: 完整实现字段序列化
        fields_json[field_name] = field_json;
    }
    j["fields"] = fields_json;

    // 耦合配置
    if (has_coupling)
    {
        json coupling_json;
        coupling_json["type"] = coupling.type;
        coupling_json["max_iterations"] = coupling.max_iterations;
        coupling_json["tolerance"] = coupling.tolerance;
        coupling_json["relaxation_factor"] = coupling.relaxation_factor;
        coupling_json["reference_temperature"] = coupling.reference_temperature;
        j["coupling"] = coupling_json;
    }

    // 写入文件
    std::ofstream file(json_file);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file for writing: " + json_file);
    }
    file << j.dump(2);
}

// ============================================================================
// ProblemLibrary 实现
// ============================================================================

std::map<std::string, std::string> ProblemLibrary::predefined_problems_;

std::shared_ptr<ProblemSetupV2> ProblemLibrary::createFromJSON(const std::string& json_file)
{
    auto problem = std::make_shared<ProblemSetupV2>();
    problem->loadFromJSON(json_file);
    return problem;
}

std::shared_ptr<ProblemSetupV2> ProblemLibrary::create(const std::string& problem_name)
{
    auto it = predefined_problems_.find(problem_name);
    if (it == predefined_problems_.end())
    {
        throw std::runtime_error("Unknown predefined problem: " + problem_name);
    }

    return createFromJSON(it->second);
}

std::vector<std::string> ProblemLibrary::getAvailableProblems()
{
    std::vector<std::string> names;
    for (const auto& [name, path] : predefined_problems_)
    {
        names.push_back(name);
    }
    return names;
}

void ProblemLibrary::registerProblem(const std::string& name, const std::string& json_file_path)
{
    predefined_problems_[name] = json_file_path;
}
