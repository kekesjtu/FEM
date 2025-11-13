#ifndef MATERIAL_V2_H
#define MATERIAL_V2_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "muparser/muParser.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

/**
 * @brief 材料属性求值的上下文
 *
 * 包含计算材料属性所需的所有物理场信息
 */
struct EvaluationContext
{
    // 空间坐标
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    // 物理场
    double T = 300.0;  // 温度 [K]
    double V = 0.0;    // 电势 [V]

    // 场梯度(可选)
    double grad_T_x = 0.0;
    double grad_T_y = 0.0;
    double grad_T_z = 0.0;
    double grad_V_x = 0.0;
    double grad_V_y = 0.0;
    double grad_V_z = 0.0;

    // 时间(瞬态问题)
    double t = 0.0;

    /**
     * @brief 转换为变量映射,供muParser使用
     */
    std::map<std::string, double*> getVariableMap()
    {
        // 注意:返回指针映射,因为muParser需要指针来动态更新值
        return {{"x", &x},
                {"y", &y},
                {"z", &z},
                {"T", &T},
                {"V", &V},
                {"t", &t},
                {"grad_T_x", &grad_T_x},
                {"grad_T_y", &grad_T_y},
                {"grad_T_z", &grad_T_z},
                {"grad_V_x", &grad_V_x},
                {"grad_V_y", &grad_V_y},
                {"grad_V_z", &grad_V_z}};
    }
};

/**
 * @brief 材料属性类
 *
 * 支持两种类型:
 * 1. 常数 - 直接存储数值
 * 2. 表达式 - 使用muParser解析和计算
 */
class MaterialProperty
{
  public:
    enum class Type
    {
        CONSTANT,   // 常数
        EXPRESSION  // 表达式
    };

    MaterialProperty();

    /**
     * @brief 从常数值构造
     */
    explicit MaterialProperty(double value);

    /**
     * @brief 从表达式构造
     * @param formula 数学表达式字符串
     * @param variables 依赖的变量列表
     * @param parameters 公式中的常数参数
     */
    MaterialProperty(const std::string& formula, const std::vector<std::string>& variables = {},
                     const std::map<std::string, double>& parameters = {});

    /**
     * @brief 从JSON加载
     */
    static MaterialProperty fromJSON(const json& j);

    /**
     * @brief 计算属性值
     */
    double evaluate(const EvaluationContext& ctx) const;

    /**
     * @brief 获取类型
     */
    Type getType() const
    {
        return type_;
    }

    /**
     * @brief 获取常数值(仅当类型为CONSTANT时有效)
     */
    double getConstantValue() const
    {
        return constant_value_;
    }

    /**
     * @brief 获取表达式字符串(如果是表达式类型)
     */
    const std::string& getFormula() const
    {
        return formula_;
    }

    /**
     * @brief 获取变量列表
     */
    const std::vector<std::string>& getVariables() const
    {
        return variables_;
    }

    /**
     * @brief 获取参数映射
     */
    const std::map<std::string, double>& getParameters() const
    {
        return parameters_;
    }

  private:
    Type type_;
    double constant_value_;
    std::string formula_;
    std::vector<std::string> variables_;
    std::map<std::string, double> parameters_;

    // muParser解析器(使用mutable允许在const函数中修改)
    mutable std::shared_ptr<mu::Parser> parser_;
    mutable bool parser_initialized_;

    /**
     * @brief 初始化表达式解析器
     */
    void initializeParser(const EvaluationContext& ctx) const;
};

/**
 * @brief 材料类
 *
 * 包含材料的所有物理属性
 */
class Material
{
  public:
    std::string name;
    std::string description;

    // 材料属性
    MaterialProperty thermal_conductivity;     // 热导率 k [W/(m·K)]
    MaterialProperty electrical_conductivity;  // 电导率 σ [S/m]
    MaterialProperty density;                  // 密度 ρ [kg/m³]
    MaterialProperty specific_heat;            // 比热容 cp [J/(kg·K)]

    /**
     * @brief 默认构造函数
     */
    Material();

    /**
     * @brief 构造函数
     */
    explicit Material(const std::string& name);

    /**
     * @brief 从JSON加载材料
     */
    static std::shared_ptr<Material> fromJSON(const std::string& material_name, const json& j);

    /**
     * @brief 转换为JSON
     */
    json toJSON() const;
};

/**
 * @brief 材料库类
 *
 * 管理所有材料和材料与几何域的映射关系
 */
class MaterialLibrary
{
  public:
    MaterialLibrary();

    /**
     * @brief 从JSON文件加载材料库
     */
    void loadFromJSON(const std::string& json_file);

    /**
     * @brief 从JSON字符串加载
     */
    void loadFromJSONString(const std::string& json_str);

    /**
     * @brief 添加材料
     */
    void addMaterial(const std::string& name, std::shared_ptr<Material> material);

    /**
     * @brief 获取材料
     */
    std::shared_ptr<Material> getMaterial(const std::string& name) const;

    /**
     * @brief 设置域到材料的映射
     */
    void assignMaterialToDomain(int domain_id, const std::string& material_name);

    /**
     * @brief 设置边界域到材料的映射
     */
    void assignMaterialToBoundaryDomain(int boundary_domain_id, const std::string& material_name);

    /**
     * @brief 根据域ID获取材料
     */
    std::shared_ptr<Material> getMaterialForDomain(int domain_id) const;

    /**
     * @brief 根据边界域ID获取材料
     */
    std::shared_ptr<Material> getMaterialForBoundaryDomain(int boundary_domain_id) const;

    /**
     * @brief 获取所有材料名称
     */
    std::vector<std::string> getMaterialNames() const;

    /**
     * @brief 保存到JSON文件
     */
    void saveToJSON(const std::string& json_file) const;

  private:
    std::map<std::string, std::shared_ptr<Material>> materials_;
    std::map<int, std::string> domain_to_material_;
    std::map<int, std::string> boundary_domain_to_material_;
};

#endif  // MATERIAL_V2_H
