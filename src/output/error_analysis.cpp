#include "error_analysis.h"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include "config.h"
#include "gauss_quadrature.h"
#include "geometry_mapping.h"
#include "material.h"
#include "shape_functions.h"

ErrorAnalysis::ErrorAnalysis(std::shared_ptr<Config> config, std::shared_ptr<ProblemSetup> problem,
                             const Eigen::VectorXd& solution, const std::string& field_name)
    : config_(config), problem_(problem), solution_(solution), field_name_(field_name)
{
}

bool ErrorAnalysis::hasExactSolution() const
{
    if (!problem_ || !problem_->hasField(field_name_))
    {
        return false;
    }
    const auto& field = problem_->getField(field_name_);
    return field.has_exact_solution;
}

bool ErrorAnalysis::hasExactGradients() const
{
    if (!problem_ || !problem_->hasField(field_name_))
    {
        return false;
    }
    const auto& field = problem_->getField(field_name_);
    return field.has_exact_solution;  // V2系统中精确解和梯度一起定义
}

double ErrorAnalysis::computeNormError(NormType norm) const
{
    // 检查是否有解析解
    if (!hasExactSolution())
    {
        return std::nan("");  // 返回 NaN 表示无法计算
    }

    switch (norm)
    {
        case NormType::L_INFINITY:
            return computeLInfinityError();
        case NormType::L2:
            return computeL2Error();
        case NormType::H1_SEMINORM:
            return computeH1SeminormError();
    }
    return 0.0;
}

void ErrorAnalysis::printErrorSummary() const
{
    std::cout << "\n=== Error Analysis Summary ===" << std::endl;
    std::cout << "Space Dimension: " << config_->getDimension() << "D" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    // 检查是否有解析解
    if (!hasExactSolution())
    {
        std::cout << "警告: 该问题没有定义解析解，无法进行误差分析" << std::endl;
        std::cout << std::flush;
        return;
    }

    const std::vector<std::pair<NormType, std::string>> norm_names = {
        {NormType::L_INFINITY, "L-infinity norm error"},
        {NormType::L2, "L2 norm error"},
        {NormType::H1_SEMINORM, "H1 seminorm error"}};

    for (const auto& [norm, name] : norm_names)
    {
        std::cout << name << ": " << std::scientific << std::setprecision(6)
                  << computeNormError(norm) << std::endl;
    }
    std::cout << std::flush;
}

void ErrorAnalysis::printDetailedNodeErrors(int max_nodes) const
{
    // 检查是否有解析解
    if (!hasExactSolution())
    {
        std::cout << "\n=== Detailed Node Errors ===" << std::endl;
        std::cout << "警告: 该问题没有定义解析解，无法进行详细节点误差分析" << std::endl;
        return;
    }

    const int sdim = config_->getDimension();
    std::cout << "\n=== Detailed Node Errors ===" << std::endl;
    std::cout << std::left << std::setw(8) << "Node";
    for (int i = 0; i < sdim; ++i)
        std::cout << std::setw(15) << "coord_" + std::to_string(i);
    std::cout << std::setw(18) << "Numerical" << std::setw(18) << "Exact" << std::setw(18)
              << "Absolute Error" << std::endl;
    std::cout << std::string(8 + sdim * 15 + 3 * 18, '-') << std::endl;

    const int num_nodes = config_->getNodesNum();
    const int display_count = std::min(max_nodes, num_nodes);
    const auto& nodes = config_->getNodeCoordinates();

    for (int i = 0; i < display_count; ++i)
    {
        std::cout << std::left << std::setw(8) << i;
        std::vector<double> coords(sdim);
        for (int d = 0; d < sdim; ++d)
        {
            coords[d] = nodes[i * sdim + d];
            std::cout << std::fixed << std::setprecision(6) << std::setw(15) << coords[d];
        }

        double numerical = solution_(i);

        // 使用MaterialProperty评估精确解
        const auto& field = problem_->getField(field_name_);
        EvaluationContext ctx;
        ctx.x = (sdim >= 1) ? coords[0] : 0.0;
        ctx.y = (sdim >= 2) ? coords[1] : 0.0;
        ctx.z = (sdim >= 3) ? coords[2] : 0.0;
        double exact = field.exact_solution_u.evaluate(ctx);

        double error = std::abs(exact - numerical);

        std::cout << std::scientific << std::setprecision(6) << std::setw(18) << numerical
                  << std::setw(18) << exact << std::setw(18) << error << std::endl;
    }

    if (display_count < num_nodes)
    {
        std::cout << "... (omitting " << (num_nodes - display_count) << " other nodes)"
                  << std::endl;
    }
}

double ErrorAnalysis::calculateNumericalSolution(int element_idx,
                                                 const std::vector<double>& coords_ref) const
{
    double numerical_value = 0.0;
    auto shape_func = ShapeFunctionFactory::createShapeFunction(config_);
    const auto& elements = config_->getElementConnectivity();
    const int nodes_per_element = config_->getNodesNumPerElement();

    for (int i = 0; i < nodes_per_element; ++i)
    {
        int global_node_idx = elements[element_idx][i];
        double shape_value = shape_func->computeTrialFunction(i, coords_ref);
        numerical_value += solution_(global_node_idx) * shape_value;
    }
    return numerical_value;
}

std::vector<double> ErrorAnalysis::calculateNumericalGradient(
    int element_idx, const std::vector<double>& coords_ref) const
{
    const int sdim = config_->getDimension();
    std::vector<double> numerical_gradient(sdim, 0.0);

    auto shape_func = ShapeFunctionFactory::createShapeFunction(config_);
    const auto& elements = config_->getElementConnectivity();
    const int nodes_per_element = config_->getNodesNumPerElement();
    const auto& all_nodes = config_->getNodeCoordinates();

    std::vector<double> element_nodes_coords;
    element_nodes_coords.reserve(nodes_per_element * sdim);
    for (int i = 0; i < nodes_per_element; ++i)
    {
        int global_node_idx = elements[element_idx][i];
        for (int d = 0; d < sdim; ++d)
        {
            element_nodes_coords.push_back(all_nodes[global_node_idx * sdim + d]);
        }
    }

    auto mapping = GeometryMappingFactory::createMapping(element_nodes_coords, config_);

    for (int i = 0; i < nodes_per_element; ++i)
    {
        int global_node_idx = elements[element_idx][i];
        auto grad_ref = shape_func->computeTrialGradients(i, coords_ref);
        std::vector<double> grad_phys(sdim);
        mapping->transformGradient(grad_ref, grad_phys, coords_ref);

        for (int d = 0; d < sdim; ++d)
        {
            numerical_gradient[d] += solution_(global_node_idx) * grad_phys[d];
        }
    }
    return numerical_gradient;
}

double ErrorAnalysis::computeLInfinityError() const
{
    // 检查是否有解析解
    if (!hasExactSolution())
    {
        return std::nan("");  // 返回 NaN 表示无法计算
    }

    double max_error = 0.0;
    const int num_elements = config_->getElementsNum();
    const int sdim = config_->getDimension();
    const int nodes_per_element = config_->getNodesNumPerElement();
    const auto& elements = config_->getElementConnectivity();
    const auto& all_nodes = config_->getNodeCoordinates();

    // 获取采样点数配置
    int num_points_per_side = config_->getSamplingPointsNum();

    // 遍历每个单元进行采样
    for (int element_index = 0; element_index < num_elements; ++element_index)
    {
        // 获取单元节点坐标
        std::vector<double> element_coords;
        element_coords.reserve(nodes_per_element * sdim);
        for (int j = 0; j < nodes_per_element; ++j)
        {
            int global_node_idx = elements[element_index][j];
            for (int d = 0; d < sdim; ++d)
            {
                element_coords.push_back(all_nodes[global_node_idx * sdim + d]);
            }
        }

        // 创建几何映射对象
        auto mapping = GeometryMappingFactory::createMapping(element_coords, config_);

        // 在单元内进行采样
        std::vector<int> indices(sdim, 0);
        bool done = false;

        while (!done)
        {
            std::vector<double> coord_ref(sdim);

            // 生成参考坐标
            for (int d = 0; d < sdim; ++d)
            {
                coord_ref[d] = double(indices[d]) / (num_points_per_side - 1);
            }

            // 检查点是否在单元内（对于三角形：所有坐标之和 <= 1）
            double sum = 0.0;
            for (int d = 0; d < sdim; ++d)
            {
                sum += coord_ref[d];
            }

            if (sum <= 1.0)
            {
                std::vector<double> coord_phys;

                // 将参考坐标转换为物理坐标
                mapping->mapToPhysical(coord_ref, coord_phys);

                // 计算该点的数值解
                double numerical_val = calculateNumericalSolution(element_index, coord_ref);

                // 计算该点的精确解 - 使用MaterialProperty
                const auto& field = problem_->getField(field_name_);
                EvaluationContext ctx;
                ctx.x = (sdim >= 1) ? coord_phys[0] : 0.0;
                ctx.y = (sdim >= 2) ? coord_phys[1] : 0.0;
                ctx.z = (sdim >= 3) ? coord_phys[2] : 0.0;
                double exact_val = field.exact_solution_u.evaluate(ctx);

                // 计算误差
                double error = std::abs(exact_val - numerical_val);
                max_error = std::max(max_error, error);
            }

            // 更新索引，用于进位法生成采样点
            int carry = 1;
            for (int d = sdim - 1; d >= 0 && carry > 0; --d)
            {
                indices[d] += carry;
                if (indices[d] >= num_points_per_side)
                {
                    indices[d] = 0;
                    carry = 1;
                }
                else
                {
                    carry = 0;
                }
            }

            if (carry > 0)
            {
                done = true;
            }
        }
    }
    return max_error;
}

double ErrorAnalysis::computeL2Error() const
{
    // 检查是否有解析解
    if (!hasExactSolution())
    {
        return std::nan("");  // 返回 NaN 表示无法计算
    }

    double l2_error_squared = 0.0;
    const int num_elements = config_->getElementsNum();
    const int sdim = config_->getDimension();
    const int nodes_per_element = config_->getNodesNumPerElement();
    const auto& elements = config_->getElementConnectivity();
    const auto& all_nodes = config_->getNodeCoordinates();

    auto gauss_points = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config_->getElementType()),
        config_->getErrorGaussPointsNum());

    for (int i = 0; i < num_elements; ++i)
    {
        std::vector<double> element_nodes_coords;
        element_nodes_coords.reserve(nodes_per_element * sdim);
        for (int j = 0; j < nodes_per_element; ++j)
        {
            int global_node_idx = elements[i][j];
            for (int d = 0; d < sdim; ++d)
            {
                element_nodes_coords.push_back(all_nodes[global_node_idx * sdim + d]);
            }
        }

        auto mapping = GeometryMappingFactory::createMapping(element_nodes_coords, config_);

        for (int gp = 0; gp < gauss_points->getNumPoints(); ++gp)
        {
            std::vector<double> gp_coords_ref;
            const auto& all_gp_coords = gauss_points->getPoints();
            gp_coords_ref.assign(all_gp_coords.begin() + gp * sdim,
                                 all_gp_coords.begin() + (gp + 1) * sdim);

            double weight = gauss_points->getWeights()[gp];
            std::vector<double> phys_coords(sdim);
            mapping->mapToPhysical(gp_coords_ref, phys_coords);
            double jacobian = mapping->getJacobianDet(gp_coords_ref);

            double numerical = calculateNumericalSolution(i, gp_coords_ref);

            // 使用MaterialProperty评估精确解
            const auto& field = problem_->getField(field_name_);
            EvaluationContext ctx;
            ctx.x = (sdim >= 1) ? phys_coords[0] : 0.0;
            ctx.y = (sdim >= 2) ? phys_coords[1] : 0.0;
            ctx.z = (sdim >= 3) ? phys_coords[2] : 0.0;
            double exact = field.exact_solution_u.evaluate(ctx);

            double error = exact - numerical;

            l2_error_squared += error * error * weight * jacobian;
        }
    }

    return std::sqrt(l2_error_squared);
}

double ErrorAnalysis::computeH1SeminormError() const
{
    // 检查是否有解析解
    if (!hasExactSolution())
    {
        return std::nan("");  // 返回 NaN 表示无法计算
    }

    double h1_error_squared = 0.0;
    const int num_elements = config_->getElementsNum();
    const int sdim = config_->getDimension();
    const int nodes_per_element = config_->getNodesNumPerElement();
    const auto& elements = config_->getElementConnectivity();
    const auto& all_nodes = config_->getNodeCoordinates();

    auto gauss_points = GaussPointFactory::createGaussPoint(
        static_cast<Config::ElementType>(config_->getElementType()),
        config_->getErrorGaussPointsNum());

    for (int i = 0; i < num_elements; ++i)
    {
        std::vector<double> element_nodes_coords;
        element_nodes_coords.reserve(nodes_per_element * sdim);
        for (int j = 0; j < nodes_per_element; ++j)
        {
            int global_node_idx = elements[i][j];
            for (int d = 0; d < sdim; ++d)
            {
                element_nodes_coords.push_back(all_nodes[global_node_idx * sdim + d]);
            }
        }

        auto mapping = GeometryMappingFactory::createMapping(element_nodes_coords, config_);

        for (int gp = 0; gp < gauss_points->getNumPoints(); ++gp)
        {
            std::vector<double> gp_coords_ref;
            const auto& all_gp_coords = gauss_points->getPoints();
            gp_coords_ref.assign(all_gp_coords.begin() + gp * sdim,
                                 all_gp_coords.begin() + (gp + 1) * sdim);

            double weight = gauss_points->getWeights()[gp];
            std::vector<double> phys_coords(sdim);
            mapping->mapToPhysical(gp_coords_ref, phys_coords);
            double jacobian = mapping->getJacobianDet(gp_coords_ref);

            auto num_grad = calculateNumericalGradient(i, gp_coords_ref);

            // 使用MaterialProperty评估精确梯度
            const auto& field = problem_->getField(field_name_);
            EvaluationContext ctx;
            ctx.x = (sdim >= 1) ? phys_coords[0] : 0.0;
            ctx.y = (sdim >= 2) ? phys_coords[1] : 0.0;
            ctx.z = (sdim >= 3) ? phys_coords[2] : 0.0;

            std::vector<double> exact_grad(sdim);
            if (sdim >= 1)
                exact_grad[0] = field.exact_solution_grad_x.evaluate(ctx);
            if (sdim >= 2)
                exact_grad[1] = field.exact_solution_grad_y.evaluate(ctx);
            if (sdim >= 3)
                exact_grad[2] = field.exact_solution_grad_z.evaluate(ctx);

            double grad_error_squared = 0;
            for (int d = 0; d < sdim; ++d)
            {
                double error_d = exact_grad[d] - num_grad[d];
                grad_error_squared += error_d * error_d;
            }

            h1_error_squared += grad_error_squared * weight * jacobian;
        }
    }

    return std::sqrt(h1_error_squared);
}
