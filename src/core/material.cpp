#include "material.h"
#include <iostream>
#include "config.h"

// ============================================================================
// Material 类实现
// ============================================================================

Material::Material(const std::string& name) : name_(name)
{
}

void Material::setCoefficient(std::function<double(const std::vector<double>&)> func)
{
    coefficient_func_ = func;
    temperature_dependent_coefficient_func_ = nullptr;  // 清除温度依赖
}

void Material::setTemperatureDependentCoefficient(
    std::function<double(const std::vector<double>&, double)> func)
{
    temperature_dependent_coefficient_func_ = func;
    coefficient_func_ = nullptr;  // 清除普通系数
}

void Material::setSource(std::function<double(const std::vector<double>&)> func)
{
    source_func_ = func;
}

double Material::coefficient(const std::vector<double>& coords) const
{
    if (temperature_dependent_coefficient_func_)
    {
        throw std::runtime_error("材料 '" + name_ +
                                 "' 的系数依赖温度，请使用 coefficient(coords, temperature)");
    }

    if (!coefficient_func_)
    {
        throw std::runtime_error("材料 '" + name_ + "' 未设置系数函数");
    }

    return coefficient_func_(coords);
}

double Material::coefficient(const std::vector<double>& coords, double temperature) const
{
    // 优先使用温度依赖函数
    if (temperature_dependent_coefficient_func_)
    {
        return temperature_dependent_coefficient_func_(coords, temperature);
    }

    // 回退到不依赖温度的函数
    if (coefficient_func_)
    {
        return coefficient_func_(coords);
    }

    throw std::runtime_error("材料 '" + name_ + "' 未设置系数函数");
}

double Material::source(const std::vector<double>& coords) const
{
    if (!source_func_)
    {
        return 0.0;  // 默认无源项
    }
    return source_func_(coords);
}

// ============================================================================
// MaterialLibrary 类实现
// ============================================================================

int MaterialLibrary::addMaterial(std::shared_ptr<Material> material)
{
    int material_id = static_cast<int>(materials_.size());
    materials_.push_back(material);
    material_name_to_id_[material->getName()] = material_id;
    return material_id;
}

std::shared_ptr<Material> MaterialLibrary::getMaterialByName(const std::string& name) const
{
    auto it = material_name_to_id_.find(name);
    if (it == material_name_to_id_.end())
    {
        throw std::runtime_error("未找到材料: " + name);
    }
    return materials_[it->second];
}

std::shared_ptr<Material> MaterialLibrary::getMaterialByID(int material_id) const
{
    if (material_id < 0 || material_id >= static_cast<int>(materials_.size()))
    {
        throw std::out_of_range("材料ID超出范围: " + std::to_string(material_id));
    }
    return materials_[material_id];
}

void MaterialLibrary::assignMaterialToDomain(int domain_id, int material_id)
{
    if (material_id < 0 || material_id >= static_cast<int>(materials_.size()))
    {
        throw std::out_of_range("材料ID超出范围: " + std::to_string(material_id));
    }
    domain_to_material_[domain_id] = material_id;
}

void MaterialLibrary::assignMaterialToDomain(int domain_id, const std::string& material_name)
{
    auto it = material_name_to_id_.find(material_name);
    if (it == material_name_to_id_.end())
    {
        throw std::runtime_error("未找到材料: " + material_name);
    }
    assignMaterialToDomain(domain_id, it->second);
}

void MaterialLibrary::assignMaterialToBoundaryDomain(int boundary_domain_id, int material_id)
{
    if (material_id < 0 || material_id >= static_cast<int>(materials_.size()))
    {
        throw std::out_of_range("材料ID超出范围: " + std::to_string(material_id));
    }
    boundary_domain_to_material_[boundary_domain_id] = material_id;
}

void MaterialLibrary::assignMaterialToBoundaryDomain(int boundary_domain_id,
                                                     const std::string& material_name)
{
    auto it = material_name_to_id_.find(material_name);
    if (it == material_name_to_id_.end())
    {
        throw std::runtime_error("未找到材料: " + material_name);
    }
    assignMaterialToBoundaryDomain(boundary_domain_id, it->second);
}

std::shared_ptr<Material> MaterialLibrary::getMaterialForDomain(int domain_id) const
{
    auto it = domain_to_material_.find(domain_id);
    if (it == domain_to_material_.end())
    {
        throw std::runtime_error("域 " + std::to_string(domain_id) + " 未分配材料");
    }
    return materials_[it->second];
}

std::shared_ptr<Material> MaterialLibrary::getMaterialForElement(int element_id) const
{
    if (element_to_domain_.empty())
    {
        throw std::runtime_error(
            "单元-域映射未设置。请先调用 setElementToDomainMap() 或 loadFromConfig()");
    }

    if (element_id < 0 || element_id >= static_cast<int>(element_to_domain_.size()))
    {
        throw std::out_of_range("单元ID超出范围: " + std::to_string(element_id));
    }

    int domain_id = element_to_domain_[element_id];
    return getMaterialForDomain(domain_id);
}

std::shared_ptr<Material> MaterialLibrary::getMaterialForBoundary(
    int boundary_id, std::shared_ptr<Config> config) const
{
    if (!config)
    {
        throw std::runtime_error("Config 对象为空");
    }

    // 从 Config 获取边界几何实体编码
    const auto& boundary_geometric_entities = config->getBoundaryGeometricEntities();

    if (boundary_id < 0 || boundary_id >= static_cast<int>(boundary_geometric_entities.size()))
    {
        throw std::out_of_range("边界ID超出范围: " + std::to_string(boundary_id));
    }

    // 获取边界的几何实体 ID
    int entity_id = boundary_geometric_entities[boundary_id];

    // 在边界域映射中查找对应的材料
    auto it = boundary_domain_to_material_.find(entity_id);
    if (it == boundary_domain_to_material_.end())
    {
        throw std::runtime_error("边界几何域 " + std::to_string(entity_id) + " 未分配材料");
    }

    return materials_[it->second];
}

void MaterialLibrary::setElementToDomainMap(const std::vector<int>& element_to_domain)
{
    element_to_domain_ = element_to_domain;
}

void MaterialLibrary::loadFromConfig(std::shared_ptr<Config> config)
{
    if (!config)
    {
        throw std::runtime_error("Config 对象为空");
    }

    // 从 Config 读取体单元几何实体映射
    const auto& element_geometric_entities = config->getElementGeometricEntities();

    if (!element_geometric_entities.empty())
    {
        element_to_domain_ = element_geometric_entities;
        std::cout << "从 Config 加载了 " << element_to_domain_.size() << " 个单元的几何实体映射"
                  << std::endl;
    }
    else
    {
        std::cerr << "警告: Config 中没有体单元几何实体信息" << std::endl;
    }
}

std::vector<std::string> MaterialLibrary::getAllMaterialNames() const
{
    std::vector<std::string> names;
    names.reserve(materials_.size());
    for (const auto& material : materials_)
    {
        names.push_back(material->getName());
    }
    return names;
}
