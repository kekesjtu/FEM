#include "problem_setup.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "logger.h"

// ============================================================================
// 辅助函数：安全地解析域ID
// ============================================================================

/**
 * @brief 安全地将字符串转换为域ID
 * @param domain_str 域ID字符串（可以是数字或"all"）
 * @return 域ID（-1表示"all"）
 * @throw std::invalid_argument 如果字符串既不是数字也不是"all"
 */
static int parseDomainId(const std::string& domain_str)
{
    if (domain_str == "all")
    {
        return -1;
    }
    try
    {
        return std::stoi(domain_str);
    }
    catch (const std::exception& e)
    {
        throw std::invalid_argument("Invalid domain ID: '" + domain_str +
                                    "'. Expected integer or 'all'.");
    }
}

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

    // 统一使用 K, L, q 字段表示边界条件
    // 通用形式: K * c * du/dn + L * u = q
    // - Dirichlet (K=0, L≠0): u = q/L
    // - Neumann (K≠0, L=0): c*du/dn = q/K
    // - Robin (K≠0, L≠0): c*du/dn + (L/K)*u = q/K

    if (type_str == "Dirichlet")
    {
        double val = 0.0;
        if (j.contains("value"))
        {
            val = j["value"].get<double>();
        }
        LOG_DEBUG("解析Dirichlet边界条件: value=" + std::to_string(val));
        return BoundaryCondition::Dirichlet(val);
    }
    else if (type_str == "Neumann")
    {
        // Neumann: K=1, L=0, q=flux
        int K = 1;
        double q = 0.0;

        if (j.contains("K"))
        {
            K = j["K"].get<double>();
        }
        if (j.contains("q"))
        {
            q = j["q"].get<double>();
        }
        LOG_DEBUG("解析Neumann边界条件: K=" + std::to_string(K) + ", q=" + std::to_string(q));
        return BoundaryCondition::Neumann(q);
    }
    else if (type_str == "Robin")
    {
        // Robin: K≠0, L≠0
        int K = 1;
        double L = 1.0;
        double q = 0.0;

        if (j.contains("K"))
        {
            K = j["K"].get<double>();
        }
        if (j.contains("L"))
        {
            L = j["L"].get<double>();
        }
        if (j.contains("q"))
        {
            q = j["q"].get<double>();
        }
        LOG_DEBUG("解析Robin边界条件: K=" + std::to_string(K) + ", L=" + std::to_string(L) +
                  ", q=" + std::to_string(q));
        return BoundaryCondition::Robin(L, q);
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
            // computed 类型：源项由其他场计算得到
            field.source_computed = true;
            // domain_to_source 保持为空
        }
        else if (source_json.contains("domains"))
        {
            // 新格式：支持每个域有不同的源项
            field.source_computed = false;
            for (const auto& [domain_str, domain_source] : source_json["domains"].items())
            {
                int domain_id = parseDomainId(domain_str);
                if (domain_source.is_number())
                {
                    field.domain_to_source[domain_id] =
                        MaterialProperty(domain_source.get<double>());
                }
                else if (domain_source.is_object())
                {
                    field.domain_to_source[domain_id] = MaterialProperty::fromJSON(domain_source);
                }
            }
        }
        else
        {
            // 旧格式：单一源项，映射到默认域 (-1)
            field.source_computed = false;
            if (source_json.is_number())
            {
                field.domain_to_source[-1] = MaterialProperty(source_json.get<double>());
            }
            else if (source_json.is_object())
            {
                field.domain_to_source[-1] = MaterialProperty::fromJSON(source_json);
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
                    int eid = entity_id.get<int>();
                    field.entity_to_bc[eid] = bc;
                    LOG_DEBUG("  绑定边界条件到实体 " + std::to_string(eid) +
                              ": K=" + std::to_string(bc.K_bc) + ", L=" + std::to_string(bc.L_bc) +
                              ", q=" + std::to_string(bc.q_bc));
                }
            }
            else if (entities.is_string() && entities.get<std::string>() == "all")
            {
                // "all" 标记会在后续处理
                // 这里先存储到一个特殊的实体ID（-1）
                field.entity_to_bc[-1] = bc;
                LOG_DEBUG("  绑定边界条件到所有实体(all): K=" + std::to_string(bc.K_bc) +
                          ", L=" + std::to_string(bc.L_bc) + ", q=" + std::to_string(bc.q_bc));
            }
        }
    }

    // 解析精确�?
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
// LogConfiguration 实现
// ============================================================================

LogConfiguration LogConfiguration::fromJSON(const json& j)
{
    LogConfiguration log;

    if (j.contains("log_level"))
    {
        log.log_level = j["log_level"];
    }
    if (j.contains("log_file"))
    {
        log.log_file = j["log_file"];
    }

    return log;
}

// ============================================================================
// ProblemSetup 实现
// ============================================================================

ProblemSetup::ProblemSetup()
    : name("unnamed_problem"),
      description(""),
      dimension(2),
      mesh_file(""),
      material_library_file(""),
      mesh_unit_scale(1.0),
      has_coupling(false)
{
}

void ProblemSetup::loadFromJSON(const std::string& json_file)
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

    // 先解析日志配置并立即应用，避免后续解析过程输出不必要的DEBUG日志
    parseLog(j);
    applyLogConfiguration();

    parseBasicInfo(j);
    parseDomainMaterials(j);
    parseFields(j);
    parseSolver(j);
    parseCoupling(j);
}

void ProblemSetup::parseBasicInfo(const json& j)
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
    if (j.contains("mesh_unit_scale"))
    {
        mesh_unit_scale = j["mesh_unit_scale"];
        LOG_DEBUG("从JSON读取网格单位缩放因子: " + std::to_string(mesh_unit_scale));
    }
}

void ProblemSetup::parseDomainMaterials(const json& j)
{
    if (j.contains("domain_materials"))
    {
        const auto& dm = j["domain_materials"];

        // 支持两种格式：
        // 格式1（传统）：{"1": "材料名", "2": "材料名", ...}
        // 格式2（批量）：{"材料名": [1, 2, 3], "材料名2": [4, 5], ...}

        // 检测格式：如果第一个元素的值是数组，则为格式2
        bool is_batch_format = false;
        if (!dm.empty())
        {
            auto first_item = dm.items().begin();
            if (first_item.value().is_array())
            {
                is_batch_format = true;
            }
        }

        if (is_batch_format)
        {
            // 格式2：{"材料名": [域ID列表]}
            for (auto& [material_name, domain_ids] : dm.items())
            {
                if (domain_ids.is_array())
                {
                    for (const auto& domain_id : domain_ids)
                    {
                        int did = domain_id.get<int>();
                        domain_to_material[did] = material_name;
                        LOG_DEBUG("  域 " + std::to_string(did) + " -> 材料 '" + material_name +
                                  "'");
                    }
                }
            }
        }
        else
        {
            // 格式1：{"域ID": "材料名"}
            for (auto& [domain_id_str, material_name] : dm.items())
            {
                int domain_id = parseDomainId(domain_id_str);
                std::string mat_name = material_name;
                domain_to_material[domain_id] = mat_name;
                LOG_DEBUG("  域 " + std::to_string(domain_id) + " -> 材料 '" + mat_name + "'");
            }
        }
    }
}

void ProblemSetup::parseFields(const json& j)
{
    if (j.contains("fields"))
    {
        for (auto& [field_name, field_json] : j["fields"].items())
        {
            fields[field_name] = FieldConfiguration::fromJSON(field_name, field_json);
        }
    }
}

void ProblemSetup::parseSolver(const json& j)
{
    if (j.contains("solver"))
    {
        solver = SolverConfiguration::fromJSON(j["solver"]);
    }
    // 如果没有solver字段，使用默认值
}

void ProblemSetup::parseLog(const json& j)
{
    if (j.contains("log"))
    {
        log = LogConfiguration::fromJSON(j["log"]);
    }
    // 如果没有log字段，使用默认值
}

void ProblemSetup::applyLogConfiguration()
{
    // 应用日志配置到全局Logger
    if (!log.log_level.empty())
    {
        Logger::getInstance().setLogLevel(stringToLogLevel(log.log_level));
        LOG_DEBUG("已应用日志级别: " + log.log_level);
    }
    // 暂不支持log_file，可以后续添加
}

void ProblemSetup::parseCoupling(const json& j)
{
    if (j.contains("coupling"))
    {
        has_coupling = true;
        coupling = CouplingConfiguration::fromJSON(j["coupling"]);
    }
}

const FieldConfiguration& ProblemSetup::getField(const std::string& field_name) const
{
    auto it = fields.find(field_name);
    if (it == fields.end())
    {
        throw std::runtime_error("Field not found: " + field_name);
    }
    return it->second;
}

bool ProblemSetup::hasField(const std::string& field_name) const
{
    return fields.find(field_name) != fields.end();
}

std::string ProblemSetup::getMaterialForDomain(int domain_id) const
{
    auto it = domain_to_material.find(domain_id);
    if (it == domain_to_material.end())
    {
        throw std::runtime_error("No material assigned to domain: " + std::to_string(domain_id));
    }
    return it->second;
}

BoundaryCondition ProblemSetup::getBoundaryConditionForEntity(const std::string& field_name,
                                                              int entity_id) const
{
    const auto& field = getField(field_name);

    // 首先检查是否有"all"边界条件
    auto it_all = field.entity_to_bc.find(-1);

    // 然后查找具体的实�?
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

void ProblemSetup::validate(bool warn_only) const
{
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    // 检查基本信�?
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

    // 检查每个场的配�?
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
            error_msg += "  �?" + error + "\n";
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

void ProblemSetup::saveToJSON(const std::string& json_file) const
{
    json j;

    // 基本信息
    j["problem_name"] = name;
    j["description"] = description;
    j["dimension"] = dimension;
    j["mesh_file"] = mesh_file;
    j["material_library_file"] = material_library_file;

    // 域材料映�?
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
        // TODO: 完整实现字段序列�?
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

std::shared_ptr<ProblemSetup> ProblemLibrary::createFromJSON(const std::string& json_file)
{
    auto problem = std::make_shared<ProblemSetup>();
    problem->loadFromJSON(json_file);
    return problem;
}

std::shared_ptr<ProblemSetup> ProblemLibrary::create(const std::string& problem_name)
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
