#include "error_analysis_2d.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include "gauss_quadrature_2d.h"
#include "geometry_mapping_2d.h"
#include "shape_functions_2d.h"

// ================================
// ErrorAnalysis 基类实现
// ================================

void ErrorAnalysis::printErrorSummary() const
{
    std::cout << "\n=== 误差分析摘要 ===" << std::endl;
    std::cout << "空间维数: " << dimension_ << "D" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    // 按顺序输出各种范数误差
    const std::vector<std::pair<NormType, std::string>> norm_names = {
        {NormType::L_INFINITY, "L∞范数误差"},
        {NormType::L2, "L2范数误差"},
        {NormType::H1_SEMINORM, "H1半范数误差"}};

    for (const auto [norm, name] : norm_names)
    {
        std::cout << name << ": " << std::scientific << std::setprecision(6)
                  << computeNormError(norm, solution) << std::endl;
    }
}

// ================================
// ErrorAnalysis2D 类实现
// ================================

bool ErrorAnalysis::supportsNorm(NormType norm) const
{
    switch (norm)
    {
        case NormType::L_INFINITY:
            return hasExactSolution();
        case NormType::L2:
            return hasExactSolution();
        case NormType::H1_SEMINORM:
            return hasExactGradients();
        default:
            return false;
    }
}

double ErrorAnalysis::computeNormError(NormType norm, const Eigen::VectorXd& solution) const
{
    if (!supportsNorm(norm))
    {
        std::cerr << "警告: 缺少精确解或梯度" << std::endl;
        return 0.0;
    }

    double error = 0.0;
    switch (norm)
    {
        case NormType::L_INFINITY:
            error = computeLInfinityError(solution);
            break;
        case NormType::L2:
            error = computeL2Error(solution);
            break;
        case NormType::H1_SEMINORM:
            error = computeH1SeminormError(solution);
            break;
    }
    return error;
}

double ErrorAnalysis::computeLInfinityError(const Eigen::VectorXd& solution) const
{
    double max_error = 0.0;

    // 在每个单元内密集采样求最大误差
    for (int elem = 0; elem < mesh_->getNumElements(); ++elem)
    {
        // 获取单元节点坐标（使用 mesh 方法）
        auto element_coords = mesh_->getElementNodes(elem);

        auto mapping = GeometryMappingFactory::createMapping(
            GeometryMappingFactory::ElementType::Triangle, 1, element_coords);

        // 在三角形参考单元内采样
        for (int i = 0; i < max_error_sampling_points_; ++i)
        {
            for (int j = 0; j < max_error_sampling_points_ - i; ++j)
            {
                double xi = double(i) / (max_error_sampling_points_ - 1);
                double eta = double(j) / (max_error_sampling_points_ - 1);

                // 确保在参考三角形内
                if (xi + eta <= 1.0)
                {
                    // 转换到物理坐标
                    double x, y;
                    mapping.mapToPhysical(xi, eta, x, y);

                    // 计算误差
                    double numerical = calculateNumericalSolution(elem, solution, xi, eta);
                    double exact = exact_solution_(x, y);
                    double error = std::abs(exact - numerical);

                    max_error = std::max(max_error, error);
                }
            }
        }
    }

    return max_error;
}

double ErrorAnalysis::computeL2Error(const Eigen::VectorXd& solution) const
{
    if (!hasExactSolution())
    {
        return 0.0;
    }

    double l2_error_squared = 0.0;

    // 使用高斯积分计算 L2 误差
    auto gauss_points = GaussPointFactory::createGaussPoint(
        GaussPointFactory::ElementType::Triangle, gauss_integration_points_);

    for (int elem = 0; elem < mesh_->getNumElements(); ++elem)
    {
        // 获取单元节点坐标（使用 mesh 方法）
        auto element_coords = mesh_->getElementNodes(elem);

        GeometryMapping2D mapping(element_coords);
        double jacobian = mapping.getJacobianDet();

        // 获取积分点坐标和权重向量
        const auto& points = gauss_points->getPoints();
        const auto& weights = gauss_points->getWeights();

        // 高斯积分
        for (int gp = 0; gp < gauss_points->getNumPoints(); ++gp)
        {
            double xi = points[gp * 2];
            double eta = points[gp * 2 + 1];
            double weight = weights[gp];

            // 物理坐标
            double x, y;
            mapping.mapToPhysical(xi, eta, x, y);

            // 计算误差
            double numerical = calculateNumericalSolution(elem, solution, xi, eta);
            double exact = exact_solution_(x, y);
            double error = exact - numerical;

            l2_error_squared += weight * jacobian * error * error;
        }
    }

    return std::sqrt(l2_error_squared);
}

double ErrorAnalysis::computeH1SeminormError(const Eigen::VectorXd& solution) const
{
    if (!hasExactGradients())
    {
        return 0.0;
    }

    double h1_error_squared = 0.0;

    // 使用高斯积分计算 H1 半范数误差
    auto gauss_points = GaussPointFactory::createGaussPoint(
        GaussPointFactory::ElementType::Triangle, gauss_integration_points_);

    for (int elem = 0; elem < mesh_->getNumElements(); ++elem)
    {
        // 获取单元节点坐标（使用 mesh 方法）
        auto element_coords = mesh_->getElementNodes(elem);

        GeometryMapping2D mapping(element_coords);
        double jacobian = mapping.getJacobianDet();

        // 获取积分点坐标和权重向量
        const auto& points = gauss_points->getPoints();
        const auto& weights = gauss_points->getWeights();

        // 高斯积分
        for (int gp = 0; gp < gauss_points->getNumPoints(); ++gp)
        {
            double xi = points[gp * 2];
            double eta = points[gp * 2 + 1];
            double weight = weights[gp];

            // 物理坐标
            double x, y;
            mapping.mapToPhysical(xi, eta, x, y);

            // 计算梯度误差
            double num_grad_x, num_grad_y;
            calculateNumericalGradient(elem, solution, xi, eta, num_grad_x, num_grad_y);

            double exact_grad_x = exact_du_dx_(x, y);
            double exact_grad_y = exact_du_dy_(x, y);

            double grad_error_x = exact_grad_x - num_grad_x;
            double grad_error_y = exact_grad_y - num_grad_y;

            // H1 半范数只包含梯度项
            double integrand = grad_error_x * grad_error_x + grad_error_y * grad_error_y;
            h1_error_squared += weight * jacobian * integrand;
        }
    }

    return std::sqrt(h1_error_squared);
}

double ErrorAnalysis::calculateNumericalSolution(int element_idx, const Eigen::VectorXd& solution,
                                                 std::vector<double> coords) const
{
    double numerical_value = 0.0;
    const auto& connectivity = mesh_->getElementConnectivity();

    for (int alpha = 0; alpha < mesh_->getNodesPerElement(); ++alpha)
    {
        int global_node = connectivity[element_idx][alpha];
        double shape_value = g_triangleShapeFunction->computeTrialFunction2D(alpha, coords[0], coords[1]);
        numerical_value += solution(global_node) * shape_value;
    }

    return numerical_value;
}

void ErrorAnalysis::calculateNumericalGradient(int element_idx, const Eigen::VectorXd& solution,
                                               std::vector<double> coords,
                                               std::vector<double>& gradients) const
{
    gradients.resize(2);
    double& grad_x = gradients[0];
    double& grad_y = gradients[1];

    // 获取单元节点坐标（使用 mesh 方法）
    auto element_coords = mesh_->getElementNodes(element_idx);
    const auto& connectivity = mesh_->getElementConnectivity();

    GeometryMapping2D mapping(element_coords);

    for (int alpha = 0; alpha < mesh_->getNodesPerElement(); ++alpha)
    {
        int global_node = connectivity[element_idx][alpha];

        // 参考单元内的形函数导数
        double dphi_dxi = g_triangleShapeFunction->computeTrialDerivativeXi(alpha, xi, eta);
        double dphi_deta = g_triangleShapeFunction->computeTrialDerivativeEta(alpha, xi, eta);

        // 变换到物理空间
        double dphi_dx, dphi_dy;
        mapping.transformGradient(dphi_dxi, dphi_deta, dphi_dx, dphi_dy);

        // 累加
        grad_x += solution(global_node) * dphi_dx;
        grad_y += solution(global_node) * dphi_dy;
    }
}

void ErrorAnalysis2D::printDetailedNodeErrors(const Eigen::VectorXd& solution, int max_nodes) const
{
    if (!hasExactSolution())
    {
        std::cout << "无精确解，无法显示节点误差详情" << std::endl;
        return;
    }

    std::cout << "\n=== 节点误差详情 ===" << std::endl;
    std::cout << std::left << std::setw(6) << "节点" << std::setw(12) << "x坐标" << std::setw(12)
              << "y坐标" << std::setw(15) << "数值解" << std::setw(15) << "精确解" << std::setw(15)
              << "绝对误差" << std::endl;
    std::cout << std::string(75, '-') << std::endl;

    int display_count = std::min(max_nodes, mesh_->getNumNodes());
    const auto& coordinates = mesh_->getNodeCoordinates();

    for (int i = 0; i < display_count; ++i)
    {
        double x = coordinates[2 * i];
        double y = coordinates[2 * i + 1];
        double numerical = solution(i);
        double exact = exact_solution_(x, y);
        double error = std::abs(exact - numerical);

        std::cout << std::left << std::setw(6) << i << std::setw(12) << std::fixed
                  << std::setprecision(6) << x << std::setw(12) << std::fixed
                  << std::setprecision(6) << y << std::setw(15) << std::scientific
                  << std::setprecision(6) << numerical << std::setw(15) << std::scientific
                  << std::setprecision(6) << exact << std::setw(15) << std::scientific
                  << std::setprecision(6) << error << std::endl;
    }

    if (display_count < mesh_->getNumNodes())
    {
        std::cout << "... (省略剩余 " << (mesh_->getNumNodes() - display_count) << " 个节点)"
                  << std::endl;
    }
}

// ================================
// ErrorAnalysisFactory 工厂实现
// ================================

std::unique_ptr<ErrorAnalysis2D> ErrorAnalysisFactory::create2D(
    std::shared_ptr<Mesh2D> mesh, ErrorAnalysis2D::ExactSolutionFunc exact_sol,
    ErrorAnalysis2D::ExactGradientFunc exact_dx, ErrorAnalysis2D::ExactGradientFunc exact_dy)
{
    return std::make_unique<ErrorAnalysis2D>(mesh, exact_sol, exact_dx, exact_dy);
}

std::unique_ptr<ErrorAnalysis2D> ErrorAnalysisFactory::create2DWithFullSolution(
    std::shared_ptr<Mesh2D> mesh, ErrorAnalysis2D::ExactSolutionFunc exact_sol,
    ErrorAnalysis2D::ExactGradientFunc exact_dx, ErrorAnalysis2D::ExactGradientFunc exact_dy)
{
    if (!exact_sol || !exact_dx || !exact_dy)
    {
        throw std::invalid_argument("完整解析解需要提供解函数和两个方向的导数函数");
    }

    return std::make_unique<ErrorAnalysis2D>(mesh, exact_sol, exact_dx, exact_dy);
}