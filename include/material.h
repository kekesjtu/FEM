#ifndef MATERIAL_H
#define MATERIAL_H

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// 前置声明
class Config;

/**
 * @brief 材料属性设置类
 *
 * 管理材料的物理属性（系数、源项等），与 ProblemSetup 类似的设计模式
 *
 * 设计理念：
 * - 每个材料对应一个几何域（domain）
 * - 通过 Config 读取单元-域映射关系
 * - 支持常数、空间依赖、温度依赖的系数
 *
 * 使用流程：
 * 1. 创建 Material 对象并设置系数函数
 * 2. 创建 MaterialLibrary 管理多个材料
 * 3. 通过 Config 读取单元-域映射（将来实现）
 * 4. 在求解器中通过单元编号获取对应材料的系数
 */
class Material
{
  public:
    /**
     * @brief 构造函数
     * @param name 材料名称
     */
    Material(const std::string& name = "DefaultMaterial");

    /**
     * @brief 获取材料名称
     */
    const std::string& getName() const
    {
        return name_;
    }

    // ============================================================================
    // 设置函数（类似 ProblemSetup 的接口风格）
    // ============================================================================

    /**
     * @brief 设置材料系数函数 c(x, y, z)
     * @param func 系数函数，可以是常数或空间依赖
     *
     * 示例：
     * - 常数: setCoefficient([](const std::vector<double>&) { return 400.0; });
     * - 空间依赖: setCoefficient([](const std::vector<double>& coords) {
     *              return 400.0 + 50.0 * coords[0];
     *            });
     */
    void setCoefficient(std::function<double(const std::vector<double>&)> func);

    /**
     * @brief 设置温度依赖的系数函数 c(x, y, z, T)（用于非线性问题）
     * @param func 系数函数，依赖于位置和温度
     *
     * 示例（电导率的温度依赖）：
     * setTemperatureDependentCoefficient([](const std::vector<double>&, double T) {
     *     double sigma0 = 5.96e7;
     *     double alpha = 0.00393;
     *     double T0 = 300.0;
     *     return sigma0 / (1.0 + alpha * (T - T0));
     * });
     */
    void setTemperatureDependentCoefficient(
        std::function<double(const std::vector<double>&, double)> func);

    /**
     * @brief 设置源项函数 f(x, y, z)（可选，某些材料可能有内部热源等）
     * @param func 源项函数
     */
    void setSource(std::function<double(const std::vector<double>&)> func);

    // ============================================================================
    // 访问函数（类似 ProblemSetup 的接口风格）
    // ============================================================================

    /**
     * @brief 计算材料系数 c(x, y, z)
     * @param coords 空间坐标
     * @return 系数值
     */
    double coefficient(const std::vector<double>& coords) const;

    /**
     * @brief 计算材料系数 c(x, y, z, T)（温度依赖）
     * @param coords 空间坐标
     * @param temperature 温度
     * @return 系数值
     */
    double coefficient(const std::vector<double>& coords, double temperature) const;

    /**
     * @brief 计算源项 f(x, y, z)
     * @param coords 空间坐标
     * @return 源项值
     */
    double source(const std::vector<double>& coords) const;

    /**
     * @brief 检查是否设置了系数函数
     */
    bool hasCoefficient() const
    {
        return coefficient_func_ != nullptr;
    }

    /**
     * @brief 检查是否设置了温度依赖的系数函数
     */
    bool isTemperatureDependent() const
    {
        return temperature_dependent_coefficient_func_ != nullptr;
    }

    /**
     * @brief 检查是否设置了源项函数
     */
    bool hasSource() const
    {
        return source_func_ != nullptr;
    }

  private:
    std::string name_;

    // 系数函数（空间依赖）
    std::function<double(const std::vector<double>&)> coefficient_func_;

    // 系数函数（温度依赖）
    std::function<double(const std::vector<double>&, double)>
        temperature_dependent_coefficient_func_;

    // 源项函数
    std::function<double(const std::vector<double>&)> source_func_;
};

/**
 * @brief 材料库类
 *
 * 管理多个材料，并通过几何域（domain）编号关联到单元
 *
 * 工作流程：
 * 1. 创建并添加材料到材料库
 * 2. 将材料分配给不同的几何域
 * 3. 从 Config 读取单元-域映射（将来实现）
 * 4. 通过单元编号查询对应的材料
 */
class MaterialLibrary
{
  public:
    /**
     * @brief 添加材料到库中
     * @param material 材料对象
     * @return 材料ID（从0开始）
     */
    int addMaterial(std::shared_ptr<Material> material);

    /**
     * @brief 通过名称获取材料
     * @param name 材料名称
     * @return 材料对象指针
     */
    std::shared_ptr<Material> getMaterialByName(const std::string& name) const;

    /**
     * @brief 通过ID获取材料
     * @param material_id 材料ID
     * @return 材料对象指针
     */
    std::shared_ptr<Material> getMaterialByID(int material_id) const;

    /**
     * @brief 将材料分配给几何域
     * @param domain_id 几何域编号（从0开始）
     * @param material_id 材料ID
     *
     * 注意：domain_id 对应网格文件中的几何实体标签
     */
    void assignMaterialToDomain(int domain_id, int material_id);

    /**
     * @brief 将材料分配给几何域（通过材料名称）
     * @param domain_id 几何域编号
     * @param material_name 材料名称
     */
    void assignMaterialToDomain(int domain_id, const std::string& material_name);

    /**
     * @brief 将材料分配给边界几何域
     * @param boundary_domain_id 边界几何域编号（entity_id）
     * @param material_id 材料ID
     */
    void assignMaterialToBoundaryDomain(int boundary_domain_id, int material_id);

    /**
     * @brief 将材料分配给边界几何域（通过材料名称）
     * @param boundary_domain_id 边界几何域编号
     * @param material_name 材料名称
     */
    void assignMaterialToBoundaryDomain(int boundary_domain_id, const std::string& material_name);

    /**
     * @brief 获取指定域的材料
     * @param domain_id 几何域编号
     * @return 材料对象指针
     */
    std::shared_ptr<Material> getMaterialForDomain(int domain_id) const;

    /**
     * @brief 获取指定单元的材料
     * @param element_id 单元编号
     * @return 材料对象指针
     *
     * 注意：需要先通过 setElementToDomainMap() 设置映射关系
     *      或者通过 Config 读取（将来实现）
     */
    std::shared_ptr<Material> getMaterialForElement(int element_id) const;

    /**
     * @brief 获取指定边界的材料
     * @param boundary_id 边界单元编号（在 Config::boundarys 中的索引）
     * @param config Config 对象，用于查询几何实体映射
     * @return 材料对象指针
     */
    std::shared_ptr<Material> getMaterialForBoundary(int boundary_id,
                                                     std::shared_ptr<Config> config) const;

    /**
     * @brief 设置单元到域的映射
     * @param element_to_domain element_to_domain[e] = domain_id
     *
     * 注意：这是临时接口，将来会从 Config 中读取
     */
    void setElementToDomainMap(const std::vector<int>& element_to_domain);

    /**
     * @brief 从 Config 读取单元-域映射（将来实现）
     * @param config Config 对象
     *
     * 这个函数将来会读取 Config 中存储的几何实体信息
     */
    void loadFromConfig(std::shared_ptr<Config> config);

    /**
     * @brief 获取材料数量
     */
    int getMaterialCount() const
    {
        return static_cast<int>(materials_.size());
    }

    /**
     * @brief 获取所有材料名称
     */
    std::vector<std::string> getAllMaterialNames() const;

  private:
    // 材料存储
    std::vector<std::shared_ptr<Material>> materials_;

    // 材料名称到ID的映射
    std::map<std::string, int> material_name_to_id_;

    // 域到材料的映射（domain_to_material_[domain_id] = material_id）
    std::map<int, int> domain_to_material_;

    // 边界域到材料的映射（boundary_domain_to_material_[boundary_domain_id] = material_id）
    std::map<int, int> boundary_domain_to_material_;

    // 单元到域的映射（element_to_domain_[element_id] = domain_id）
    // 将来会从 Config 中读取，现在作为临时存储
    std::vector<int> element_to_domain_;
};

#endif  // MATERIAL_H
