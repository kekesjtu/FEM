#ifndef MATERIAL_H
#define MATERIAL_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "muparser/muParser.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

/**
 * @brief Evaluation context for material properties
 *
 * Contains all physical field information needed to compute material properties
 */
struct EvaluationContext
{
    // Spatial coordinates
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    // Physical fields
    double T = 300.0;  // Temperature [K]
    double V = 0.0;    // Electric potential [V]

    // Field gradients (optional)
    double grad_T_x = 0.0;
    double grad_T_y = 0.0;
    double grad_T_z = 0.0;
    double grad_V_x = 0.0;
    double grad_V_y = 0.0;
    double grad_V_z = 0.0;

    // Time (for transient problems)
    double t = 0.0;

    /**
     * @brief Convert to variable map for muParser
     */
    std::map<std::string, double*> getVariableMap()
    {
        // Note: Returns pointer map because muParser needs pointers for dynamic updates
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
 * @brief Material property class
 *
 * Supports two types:
 * 1. Constant - stores numerical value directly
 * 2. Expression - uses muParser to parse and evaluate
 */
class MaterialProperty
{
  public:
    enum class Type
    {
        CONSTANT,   // Constant value
        EXPRESSION  // Mathematical expression
    };

    MaterialProperty();

    /**
     * @brief Construct from constant value
     */
    explicit MaterialProperty(double value);

    /**
     * @brief Construct from expression
     * @param formula Mathematical expression string
     * @param variables List of dependent variables
     * @param parameters Constant parameters in the formula
     */
    MaterialProperty(const std::string& formula, const std::vector<std::string>& variables = {},
                     const std::map<std::string, double>& parameters = {});

    /**
     * @brief Load from JSON
     */
    static MaterialProperty fromJSON(const json& j);

    /**
     * @brief Evaluate property value
     */
    double evaluate(const EvaluationContext& ctx) const;

    /**
     * @brief Get type
     */
    Type getType() const
    {
        return type_;
    }

    /**
     * @brief Get constant value (only valid when type is CONSTANT)
     */
    double getConstantValue() const
    {
        return constant_value_;
    }

    /**
     * @brief Get formula string (if type is EXPRESSION)
     */
    const std::string& getFormula() const
    {
        return formula_;
    }

    /**
     * @brief Get variable list
     */
    const std::vector<std::string>& getVariables() const
    {
        return variables_;
    }

    /**
     * @brief Get parameter map
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

    // muParser parser (mutable to allow modification in const functions)
    mutable std::shared_ptr<mu::Parser> parser_;
    mutable bool parser_initialized_;

    // Storage for variable values (muParser needs pointers to these)
    mutable std::map<std::string, double> variable_values_;

    /**
     * @brief Initialize expression parser
     */
    void initializeParser(const EvaluationContext& ctx) const;
};

/**
 * @brief Material class
 *
 * Contains all physical properties of a material
 */
class Material
{
  public:
    std::string name;
    std::string description;

    // Material properties
    MaterialProperty thermal_conductivity;     // k [W/(m·K)]
    MaterialProperty electrical_conductivity;  // sigma [S/m]
    MaterialProperty density;                  // rho [kg/m³]
    MaterialProperty specific_heat;            // cp [J/(kg·K)]

    /**
     * @brief Default constructor
     */
    Material();

    /**
     * @brief Constructor
     */
    explicit Material(const std::string& name);

    /**
     * @brief Load material from JSON
     */
    static std::shared_ptr<Material> fromJSON(const std::string& material_name, const json& j);

    /**
     * @brief Convert to JSON
     */
    json toJSON() const;
};

/**
 * @brief Material library class
 *
 * Manages materials and their mappings to domains/boundaries
 */
class MaterialLibrary
{
  public:
    MaterialLibrary();

    /**
     * @brief Load material library from JSON file
     */
    void loadFromJSON(const std::string& json_file);

    /**
     * @brief Load from JSON string
     */
    void loadFromJSONString(const std::string& json_str);

    /**
     * @brief Add material
     */
    void addMaterial(const std::string& name, std::shared_ptr<Material> material);

    /**
     * @brief Get material
     */
    std::shared_ptr<Material> getMaterial(const std::string& name) const;

    /**
     * @brief Set domain to material mapping
     */
    void assignMaterialToDomain(int domain_id, const std::string& material_name);

    /**
     * @brief Set boundary domain to material mapping
     */
    void assignMaterialToBoundaryDomain(int boundary_domain_id, const std::string& material_name);

    /**
     * @brief Get material for domain
     */
    std::shared_ptr<Material> getMaterialForDomain(int domain_id) const;

    /**
     * @brief Get material for boundary domain
     */
    std::shared_ptr<Material> getMaterialForBoundaryDomain(int boundary_domain_id) const;

    /**
     * @brief Get all material names
     */
    std::vector<std::string> getMaterialNames() const;

    /**
     * @brief Save to JSON file
     */
    void saveToJSON(const std::string& json_file) const;

  private:
    std::map<std::string, std::shared_ptr<Material>> materials_;
    std::map<int, std::string> domain_to_material_;
    std::map<int, std::string> boundary_domain_to_material_;
};

#endif  // MATERIAL_H
