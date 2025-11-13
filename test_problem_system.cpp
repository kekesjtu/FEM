#include <iomanip>
#include <iostream>
#include "material_v2.h"
#include "problem_v2.h"

void test_problem_loading()
{
    std::cout << "\n=== 测试问题配置加载 ===" << std::endl;

    try
    {
        // 加载电热耦合问题
        auto problem = std::make_shared<ProblemSetupV2>();
        problem->loadFromJSON("problem_electrothermal_2d.json");

        std::cout << "问题名称: " << problem->name << std::endl;
        std::cout << "描述: " << problem->description << std::endl;
        std::cout << "维度: " << problem->dimension << std::endl;
        std::cout << "网格文件: " << problem->mesh_file << std::endl;
        std::cout << "材料库文件: " << problem->material_library_file << std::endl;

        // 检查域材料映射
        std::cout << "\n域材料映射:" << std::endl;
        for (const auto& [domain_id, material_name] : problem->domain_to_material)
        {
            std::cout << "  域 " << domain_id << " -> " << material_name << std::endl;
        }

        // 检查场配置
        std::cout << "\n配置的场:" << std::endl;
        for (const auto& [field_name, field] : problem->fields)
        {
            std::cout << "  - " << field_name << std::endl;
            std::cout << "    源项计算: " << (field.source_computed ? "是" : "否") << std::endl;
            std::cout << "    边界条件数: " << field.entity_to_bc.size() << std::endl;
            std::cout << "    精确解: " << (field.has_exact_solution ? "是" : "否") << std::endl;
        }

        // 检查耦合配置
        if (problem->has_coupling)
        {
            std::cout << "\n耦合配置:" << std::endl;
            std::cout << "  类型: " << problem->coupling.type << std::endl;
            std::cout << "  最大迭代次数: " << problem->coupling.max_iterations << std::endl;
            std::cout << "  收敛容差: " << problem->coupling.tolerance << std::endl;
            std::cout << "  松弛因子: " << problem->coupling.relaxation_factor << std::endl;
            std::cout << "  参考温度: " << problem->coupling.reference_temperature << " K"
                      << std::endl;
        }

        // 检查求解器配置
        std::cout << "\n求解器配置:" << std::endl;
        std::cout << "  求解器类型: " << problem->solver.solver_type << std::endl;
        std::cout << "  预条件子: " << problem->solver.preconditioner_type << std::endl;
        std::cout << "  收敛容差: " << problem->solver.tolerance << std::endl;
        std::cout << "  最大迭代次数: " << problem->solver.max_iterations << std::endl;

        std::cout << "\n✓ 电热耦合问题加载测试通过" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << "✗ 问题加载失败: " << e.what() << std::endl;
    }
}

void test_boundary_conditions()
{
    std::cout << "\n=== 测试边界条件查询 ===" << std::endl;

    try
    {
        auto problem = std::make_shared<ProblemSetupV2>();
        problem->loadFromJSON("problem_electrothermal_2d.json");

        // 测试电场边界条件
        std::cout << "\n电场边界条件:" << std::endl;
        for (int entity_id = 0; entity_id <= 3; ++entity_id)
        {
            auto bc = problem->getBoundaryConditionForEntity("electric", entity_id);
            std::string bc_type;
            if (bc.K_bc == 0 && bc.L_bc == 1.0)
            {
                bc_type = "Dirichlet";
            }
            else if (bc.K_bc == 1 && bc.L_bc == 0.0)
            {
                bc_type = "Neumann";
            }
            else
            {
                bc_type = "Robin";
            }

            std::cout << "  实体 " << entity_id << ": " << bc_type << ", 值 = " << bc.q_bc
                      << std::endl;
        }

        std::cout << "\n✓ 边界条件查询测试通过" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << "✗ 边界条件查询失败: " << e.what() << std::endl;
    }
}

void test_sphere_problem()
{
    std::cout << "\n=== 测试球域问题（带精确解） ===" << std::endl;

    try
    {
        auto problem = std::make_shared<ProblemSetupV2>();
        problem->loadFromJSON("problem_sphere_3d.json");

        std::cout << "问题名称: " << problem->name << std::endl;
        std::cout << "维度: " << problem->dimension << std::endl;
        std::cout << "网格文件: " << problem->mesh_file << std::endl;

        // 检查精确解
        if (problem->hasField("thermal"))
        {
            const auto& field = problem->getField("thermal");
            if (field.has_exact_solution)
            {
                std::cout << "\n精确解已定义，测试几个点:" << std::endl;

                std::vector<std::tuple<double, double, double>> test_points = {
                    {0.0, 0.0, 0.0},  // 中心点
                    {0.5, 0.0, 0.0},
                    {0.0, 0.5, 0.0},
                    {0.0, 0.0, 0.5}};

                std::cout << std::fixed << std::setprecision(4);
                std::cout << "\n点坐标         u(x,y,z)   ∂u/∂x    ∂u/∂y    ∂u/∂z" << std::endl;
                std::cout << "--------------------------------------------------------"
                          << std::endl;

                for (const auto& [x, y, z] : test_points)
                {
                    EvaluationContext ctx;
                    ctx.x = x;
                    ctx.y = y;
                    ctx.z = z;

                    double u = field.exact_solution_u.evaluate(ctx);
                    double grad_x = field.exact_solution_grad_x.evaluate(ctx);
                    double grad_y = field.exact_solution_grad_y.evaluate(ctx);
                    double grad_z = field.exact_solution_grad_z.evaluate(ctx);

                    std::cout << "(" << std::setw(4) << x << "," << std::setw(4) << y << ","
                              << std::setw(4) << z << ")  " << std::setw(7) << u << "  "
                              << std::setw(7) << grad_x << "  " << std::setw(7) << grad_y << "  "
                              << std::setw(7) << grad_z << std::endl;
                }
            }
        }

        std::cout << "\n✓ 球域问题加载测试通过" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << "✗ 球域问题加载失败: " << e.what() << std::endl;
    }
}

void test_validation()
{
    std::cout << "\n=== 测试配置验证 ===" << std::endl;

    try
    {
        auto problem = std::make_shared<ProblemSetupV2>();
        problem->loadFromJSON("problem_electrothermal_2d.json");

        std::cout << "\n验证配置完整性..." << std::endl;
        problem->validate(false);  // 严格模式

        std::cout << "✓ 配置验证通过" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << "✗ 配置验证失败: " << e.what() << std::endl;
    }
}

void test_integration_with_materials()
{
    std::cout << "\n=== 测试问题与材料库集成 ===" << std::endl;

    try
    {
        // 加载材料库
        auto material_lib = std::make_shared<MaterialLibrary>();
        material_lib->loadFromJSON("materials_example.json");

        // 加载问题
        auto problem = std::make_shared<ProblemSetupV2>();
        problem->loadFromJSON("problem_electrothermal_2d.json");

        std::cout << "\n模拟求解流程:" << std::endl;

        // 模拟遍历域
        for (const auto& [domain_id, material_name] : problem->domain_to_material)
        {
            std::cout << "\n域 " << domain_id << ":" << std::endl;
            std::cout << "  材料: " << material_name << std::endl;

            // 从材料库获取材料
            auto material = material_lib->getMaterial(material_name);

            // 模拟计算材料属性
            EvaluationContext ctx;
            ctx.T = 300.0;
            ctx.x = 0.5;
            ctx.y = 0.5;

            double k = material->thermal_conductivity.evaluate(ctx);
            double sigma = material->electrical_conductivity.evaluate(ctx);

            std::cout << "  热导率: " << k << " W/(m·K)" << std::endl;
            std::cout << "  电导率: " << std::scientific << sigma << " S/m" << std::endl;
        }

        // 模拟边界条件应用
        std::cout << "\n电场边界条件应用:" << std::endl;
        const auto& electric_field = problem->getField("electric");
        for (const auto& [entity_id, bc] : electric_field.entity_to_bc)
        {
            if (entity_id >= 0)
            {  // 跳过"all"标记
                std::cout << "  实体 " << entity_id << ": V = " << bc.q_bc << " V" << std::endl;
            }
        }

        std::cout << "\n✓ 问题与材料库集成测试通过" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cout << "✗ 集成测试失败: " << e.what() << std::endl;
    }
}

int main()
{
    std::cout << "======================================" << std::endl;
    std::cout << "问题配置系统测试程序" << std::endl;
    std::cout << "======================================" << std::endl;

    try
    {
        test_problem_loading();
        test_boundary_conditions();
        test_sphere_problem();
        test_validation();
        test_integration_with_materials();

        std::cout << "\n======================================" << std::endl;
        std::cout << "所有测试完成!" << std::endl;
        std::cout << "======================================" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\n异常: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
