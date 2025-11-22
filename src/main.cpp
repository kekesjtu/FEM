#include <iostream>
#include <stdexcept>
#include "config.h"
#include "electrothermal_solver.h"
#include "fem_solver.h"
#include "logger.h"
#include "material.h"
#include "problem_setup.h"

/**
 * @brief 单场求解示例（JSON驱动）
 * @param json_file 问题配置文件路径
 */
void solveSingleFieldFromJSON(const std::string& json_file)
{
    TIMER_SCOPE("单场求解总耗时");

    LOG_HEADER("单场求解 (JSON驱动)");

    try
    {
        std::shared_ptr<ProblemSetup> problem;
        std::shared_ptr<MaterialLibrary> materials;
        std::shared_ptr<Config> config;
        std::string field_name;

        // 1. 从JSON加载问题配置（内部会自动应用日志配置）
        {
            TIMER_SCOPE_DEBUG("加载问题配置");
            problem = ProblemLibrary::createFromJSON(json_file);
            LOG_INFO("问题名称: " + problem->name);
            LOG_INFO("问题描述: " + problem->description);
        }

        // 2. 从JSON加载材料库
        {
            TIMER_SCOPE_DEBUG("加载材料库");
            materials = std::make_shared<MaterialLibrary>();
            materials->loadFromJSON(problem->material_library_file);
            LOG_INFO("材料库: " + problem->material_library_file);
        }

        // 3. 创建配置对象并加载网格
        {
            TIMER_SCOPE("加载网格");
            config = std::make_shared<Config>();
            // 从problem读取网格单位缩放因子
            config->setMeshUnitScale(problem->mesh_unit_scale);
            config->setMeshFilename(problem->mesh_file);
            LOG_INFO("网格文件: " + problem->mesh_file);
        }

        // 4. 确定要求解的场（假设是第一个场）
        field_name = problem->fields.begin()->first;

        // 5. 创建FEM求解器
        FEMSolver solver(config, problem, materials, field_name);

        // 6. 一键求解（内部已有详细计时）
        solver.solveComplete();

        LOG_INFO("\n单场求解完成！");
    }
    catch (const std::exception& e)
    {
        LOG_ERROR(e.what());
        throw;
    }
}

/**
 * @brief 电热耦合求解示例（JSON驱动）
 * @param json_file 问题配置文件路径
 */
void solveElectrothermalFromJSON(const std::string& json_file)
{
    LOG_HEADER("电热耦合仿真 (JSON驱动)");

    try
    {
        // 1. 从JSON加载问题配置（内部会自动应用日志配置）
        auto problem = ProblemLibrary::createFromJSON(json_file);
        LOG_INFO("问题名称: " + problem->name);
        LOG_INFO("问题描述: " + problem->description);

        // 验证是否为电热耦合问题
        if (!problem->has_coupling)
        {
            throw std::runtime_error("配置文件不包含耦合配置");
        }
        if (!problem->hasField("electric") || !problem->hasField("thermal"))
        {
            throw std::runtime_error("电热耦合问题必须定义electric和thermal两个场");
        }

        // 2. 从JSON加载材料库
        auto materials = std::make_shared<MaterialLibrary>();
        materials->loadFromJSON(problem->material_library_file);
        LOG_INFO("材料库: " + problem->material_library_file);

        // 3. 创建配置对象并加载网格
        auto config = std::make_shared<Config>();
        // 从problem读取网格单位缩放因子
        config->setMeshUnitScale(problem->mesh_unit_scale);
        config->setMeshFilename(problem->mesh_file);

        // 4. 创建电热耦合求解器
        ElectrothermalSolver et_solver(config, problem, materials);

        // 5. 求解
        et_solver.solve();

        // 6. 输出结果
        et_solver.outputResults("results/electrothermal");

        LOG_INFO("\n电热耦合求解完成！");
        LOG_INFO("请使用ParaView查看以下文件：");
        LOG_INFO("  - results/electrothermal_V.vtu (电势场)");
        LOG_INFO("  - results/electrothermal_T.vtu (温度场)");
        LOG_INFO("  - results/electrothermal_Q.vtu (焦耳热密度)");
    }
    catch (const std::exception& e)
    {
        LOG_ERROR(e.what());
        throw;
    }
}

/**
 * @brief 主函数（V2版本 - 完全JSON驱动）
 */
int main(int argc, char* argv[])
{
    // 初始化日志系统（默认INFO级别，由JSON配置文件控制）
    Logger::getInstance().setLogLevel(LogLevel::INFO);

    LOG_HEADER("有限元求解器");

    try
    {
        std::string json_file = "problem_complex_thermal.json";

        // 支持命令行参数指定JSON文件
        if (argc > 1)
        {
            json_file = argv[1];
        }

        LOG_INFO("配置文件: " + json_file);

        // 根据文件名判断是单场还是耦合问题
        if (json_file.find("electrothermal") != std::string::npos)
        {
            solveElectrothermalFromJSON(json_file);
        }
        else
        {
            solveSingleFieldFromJSON(json_file);
        }

        LOG_HEADER("求解完成");

        return 0;
    }
    catch (const std::exception& e)
    {
        LOG_ERROR("\n程序异常终止: " + std::string(e.what()));
        return 1;
    }
}
