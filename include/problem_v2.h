#ifndef PROBLEM_V2_H
#define PROBLEM_V2_H

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "boundary_condition.h"
#include "material_v2.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

/**
 * @brief 物理场定义（电场、热场等）
 */
struct FieldConfiguration
{
    std::string name;

    // 源项
    MaterialProperty source;
    bool source_computed = false;  // 是否由求解器计算（如焦耳热）

    // 边界条件映射：几何实体ID -> 边界条件
    std::map<int, BoundaryCondition> entity_to_bc;

    // 精确解（可选）
    bool has_exact_solution = false;
    MaterialProperty exact_solution_u;
    MaterialProperty exact_solution_grad_x;
    MaterialProperty exact_solution_grad_y;
    MaterialProperty exact_solution_grad_z;

    static FieldConfiguration fromJSON(const std::string& field_name, const json& j);
};

/**
 * @brief 求解器配置
 */
struct SolverConfiguration
{
    std::string solver_type = "CG";                              // "SparseLU" 或 "CG"
    std::string preconditioner_type = "DiagonalPreconditioner";  // 预条件子类型
    double tolerance = 1e-8;                                     // 收敛容差
    int max_iterations = 1000;                                   // 最大迭代次数

    static SolverConfiguration fromJSON(const json& j);
};

/**
 * @brief 耦合问题配置（如电热耦合）
 */
struct CouplingConfiguration
{
    std::string type;
    int max_iterations = 100;
    double tolerance = 1e-6;
    double relaxation_factor = 0.7;
    double reference_temperature = 300.0;

    static CouplingConfiguration fromJSON(const json& j);
};

/**
 * @brief 问题配置类 V2
 *
 * 负责：
 * - 指定网格文件
 * - 指定几何实体使用的材料
 * - 设置源项
 * - 设置边界条件（通过几何实体）
 * - 提供精确解（可选）
 */
class ProblemSetupV2
{
  public:
    // 基本信息
    std::string name;
    std::string description;
    int dimension = 2;
    std::string mesh_file;
    std::string material_library_file;

    // 几何实体 -> 材料名映射
    std::map<int, std::string> domain_to_material;

    // 物理场配置（支持多场问题）
    std::map<std::string, FieldConfiguration> fields;  // "electric", "thermal"等

    // 求解器配置
    SolverConfiguration solver;

    // 耦合配置（可选）
    bool has_coupling = false;
    CouplingConfiguration coupling;

    /**
     * @brief 默认构造
     */
    ProblemSetupV2();

    /**
     * @brief 从JSON文件加载问题配置
     */
    void loadFromJSON(const std::string& json_file);

    /**
     * @brief 从JSON字符串加载
     */
    void loadFromJSONString(const std::string& json_str);

    /**
     * @brief 获取指定场的配置
     */
    const FieldConfiguration& getField(const std::string& field_name) const;

    /**
     * @brief 检查是否有指定的场
     */
    bool hasField(const std::string& field_name) const;

    /**
     * @brief 获取域对应的材料名
     */
    std::string getMaterialForDomain(int domain_id) const;

    /**
     * @brief 获取边界实体对应的边界条件
     */
    BoundaryCondition getBoundaryConditionForEntity(const std::string& field_name,
                                                    int entity_id) const;

    /**
     * @brief 验证配置完整性
     *
     * 检查：
     * - 网格文件是否存在
     * - 材料库文件是否存在
     * - 所有域是否都分配了材料
     * - 边界条件是否完整
     *
     * @param warn_only 如果为true，只发出警告；如果为false，抛出异常
     */
    void validate(bool warn_only = false) const;

    /**
     * @brief 保存为JSON文件
     */
    void saveToJSON(const std::string& json_file) const;

  private:
    void parseBasicInfo(const json& j);
    void parseDomainMaterials(const json& j);
    void parseFields(const json& j);
    void parseSolver(const json& j);
    void parseCoupling(const json& j);
};

/**
 * @brief 问题库工厂类
 *
 * 用于创建和管理预定义问题
 */
class ProblemLibrary
{
  public:
    /**
     * @brief 从JSON文件创建问题
     */
    static std::shared_ptr<ProblemSetupV2> createFromJSON(const std::string& json_file);

    /**
     * @brief 从预定义问题名创建
     */
    static std::shared_ptr<ProblemSetupV2> create(const std::string& problem_name);

    /**
     * @brief 列出所有可用的预定义问题
     */
    static std::vector<std::string> getAvailableProblems();

    /**
     * @brief 注册预定义问题
     */
    static void registerProblem(const std::string& name, const std::string& json_file_path);

  private:
    static std::map<std::string, std::string> predefined_problems_;
};

#endif  // PROBLEM_V2_H
