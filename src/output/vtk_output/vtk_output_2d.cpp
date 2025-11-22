#include <Eigen/Core>
#include <cmath>     // 用于std::abs
#include <fstream>   // 用于文件输出
#include <iomanip>   // 用于输出格式控制
#include <iostream>  // 用于控制台输出
#include "config.h"
#include "geometry_mapping.h"  // 用于几何映射
#include "shape_functions.h"   // 用于形函数计算
#include "vtk_output.h"


// ============================================================================
// TriangleVTKOutput2D 类的实现
// ============================================================================

void TriangleVTKOutput2D::outputNumericalSolution(const std::string& filename)
{
    ensureDirectoryExists("results");

    std::string fullname = filename + ".vtu";
    std::ofstream file(fullname);

    if (!file.is_open())
    {
        std::cerr << "无法创建VTK文件: " << fullname << std::endl;
        return;
    }

    std::cout << "正在写入二维VTK文件: " << fullname << std::endl;

    const auto& coordinates = config_->getNodeCoordinates();
    const auto& connectivity = config_->getElementConnectivity();
    int nodes_num = config_->getNodesNum();
    int elements_num = config_->getElementsNum();
    int dimension = config_->getDimension();

    file << std::fixed << std::setprecision(6);

    // VTK文件头
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << nodes_num << "\" NumberOfCells=\"" << elements_num
         << "\">\n";

    // 输出节点坐标
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        // 2D网格扩展为3D，z坐标设为0
        double x = coordinates[i * dimension];
        double y = coordinates[i * dimension + 1];
        file << "          " << x << " " << y << " 0.0\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // 输出单元连接信息
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    // 单元类型 5=VTK_TRIANGLE
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          5\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // 输出节点数据
    file << "      <PointData Scalars=\"numerical_solution\">\n";
    file << "        <DataArray type=\"Float64\" Name=\"numerical_solution\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << solution_(i) << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"x_coordinate\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << coordinates[i * dimension] << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"y_coordinate\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << coordinates[i * dimension + 1] << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </PointData>\n";

    // 输出单元数据
    file << "      <CellData>\n";
    file << "        <DataArray type=\"Int32\" Name=\"element_id\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
    std::cout << "数值解VTK文件写入完成!" << std::endl;
}

void TriangleVTKOutput2D::outputExactSolution(
    const std::string& filename,
    const std::function<double(const std::vector<double>&)>& exact_func)
{
    ensureDirectoryExists("results");

    std::string fullname = filename + ".vtu";
    std::ofstream file(fullname);

    if (!file.is_open())
    {
        std::cerr << "无法创建VTK文件: " << fullname << std::endl;
        return;
    }

    std::cout << "正在写入二维解析解VTK文件: " << fullname << std::endl;

    const auto& coordinates = config_->getNodeCoordinates();
    const auto& connectivity = config_->getElementConnectivity();
    int nodes_num = config_->getNodesNum();
    int elements_num = config_->getElementsNum();

    file << std::fixed << std::setprecision(6);

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << nodes_num << "\" NumberOfCells=\"" << elements_num
         << "\">\n";

    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << coordinates[i * 2] << " " << coordinates[i * 2 + 1] << " 0.0\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          5\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    file << "      <PointData Scalars=\"exact_solution\">\n";
    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        std::vector<double> coords = {coordinates[i * 2], coordinates[i * 2 + 1]};
        file << "          " << exact_func(coords) << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </PointData>\n";

    file << "      <CellData>\n";
    file << "        <DataArray type=\"Int32\" Name=\"element_id\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
    std::cout << "解析解VTK文件写入完成!" << std::endl;
}

void TriangleVTKOutput2D::outputDenseSamplingError(const std::string& filename,
                                                   std::shared_ptr<Config> config,
                                                   std::shared_ptr<ProblemSetup> problem,
                                                   const std::string& field_name)
{
    if (!problem || !problem->hasField(field_name))
    {
        std::cout << "警告：未设置精确解或不存在该场，跳过误差输出" << std::endl;
        return;
    }

    const auto& field = problem->getField(field_name);
    if (!field.has_exact_solution)
    {
        std::cout << "警告：该场没有定义解析解，跳过误差输出" << std::endl;
        return;
    }

    ensureDirectoryExists("results");

    std::string fullname = filename + ".vtu";
    std::ofstream file(fullname);

    if (!file.is_open())
    {
        std::cerr << "无法创建加密采样VTK文件: " << fullname << std::endl;
        return;
    }

    std::cout << "正在写入二维加密采样VTK文件: " << fullname << std::endl;
    std::cout << "采样密度: 每个三角形单元 " << num_points_per_side_ << "×" << num_points_per_side_
              << " 采样点" << std::endl;

    const auto& coordinates = config_->getNodeCoordinates();
    const auto& connectivity = config_->getElementConnectivity();
    int nodes_num = config_->getNodesNum();
    int elements_num = config_->getElementsNum();
    int dimension = config_->getDimension();

    std::vector<double> dense_points;
    std::vector<double> dense_numerical;
    std::vector<double> dense_exact;

    auto exact_func = [problem, field_name](const std::vector<double>& coords) -> double
    {
        const auto& field_local = problem->getField(field_name);
        EvaluationContext ctx;
        int dim = static_cast<int>(coords.size());
        ctx.x = (dim >= 1) ? coords[0] : 0.0;
        ctx.y = (dim >= 2) ? coords[1] : 0.0;
        ctx.z = (dim >= 3) ? coords[2] : 0.0;
        return field_local.exact_solution_u.evaluate(ctx);
    };

    // 遍历每个原始单元进行采样
    for (int e = 0; e < elements_num; ++e)
    {
        for (int i = 0; i < num_points_per_side_; ++i)
        {
            for (int j = 0; j < num_points_per_side_ - i; ++j)
            {
                double xi = double(i) / (num_points_per_side_ - 1);
                double eta = double(j) / (num_points_per_side_ - 1);

                if (xi + eta <= 1.0)
                {
                    std::vector<double> coord_ref = {xi, eta};
                    std::vector<double> coord_phys;

                    std::vector<double> element_coords;
                    for (int local_node = 0; local_node < config_->getNodesNumPerElement();
                         ++local_node)
                    {
                        int global_node = connectivity[e][local_node];
                        element_coords.push_back(coordinates[global_node * dimension]);
                        element_coords.push_back(coordinates[global_node * dimension + 1]);
                    }

                    auto mapping = GeometryMappingFactory::createMapping(element_coords, config);
                    mapping->mapToPhysical(coord_ref, coord_phys);

                    dense_points.push_back(coord_phys[0]);
                    dense_points.push_back(coord_phys[1]);
                    dense_points.push_back(0.0);

                    double numerical_val = evaluateNumericalSolutionAt(e, coord_ref, config);
                    dense_numerical.push_back(numerical_val);

                    double exact_val = exact_func(coord_phys);
                    dense_exact.push_back(exact_val);
                }
            }
        }
    }

    int total_points = dense_points.size() / 3;
    std::cout << "总采样点数: " << total_points << " (原始节点数: " << nodes_num << ")"
              << std::endl;

    file << std::fixed << std::setprecision(6);

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << total_points << "\" NumberOfCells=\"" << total_points
         << "\">\n";

    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_points.size(); i += 3)
    {
        file << "          " << dense_points[i] << " " << dense_points[i + 1] << " "
             << dense_points[i + 2] << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < total_points; ++i)
    {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < total_points; ++i)
    {
        file << "          " << i + 1 << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < total_points; ++i)
    {
        file << "          1\n";  // VTK_VERTEX
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    file << "      <PointData Scalars=\"numerical_solution\">\n";

    file << "        <DataArray type=\"Float64\" Name=\"numerical_solution\" format=\"ascii\">\n";
    for (double val : dense_numerical)
    {
        file << "          " << val << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (double val : dense_exact)
    {
        file << "          " << val << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"error\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_numerical.size(); ++i)
    {
        file << "          " << (dense_numerical[i] - dense_exact[i]) << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"absolute_error\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_numerical.size(); ++i)
    {
        file << "          " << std::abs(dense_numerical[i] - dense_exact[i]) << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"relative_error\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_numerical.size(); ++i)
    {
        double rel_error = 0.0;
        if (std::abs(dense_exact[i]) > 1e-14)
        {
            rel_error = std::abs(dense_numerical[i] - dense_exact[i]) / std::abs(dense_exact[i]);
        }
        file << "          " << rel_error << "\n";
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";
    file << "      <CellData>\n";
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
    std::cout << "加密采样误差分析VTK文件写入完成!" << std::endl;
}

double TriangleVTKOutput2D::evaluateNumericalSolutionAt(int element_index,
                                                        const std::vector<double>& coords_ref,
                                                        std::shared_ptr<Config> config) const
{
    double numerical_solution = 0.0;
    const auto& connectivity = config_->getElementConnectivity();
    int n = connectivity[element_index].size();

    auto shape_func = ShapeFunctionFactory::createShapeFunction(config);

    for (int alpha = 0; alpha < n; ++alpha)
    {
        int global_node_index = connectivity[element_index][alpha];
        double shape_val = shape_func->computeTrialFunction(alpha, coords_ref);
        numerical_solution += solution_(global_node_index) * shape_val;
    }
    return numerical_solution;
}
