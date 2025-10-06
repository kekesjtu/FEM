#include "vtk_output.h"
#include <direct.h>  // 用于Windows的_mkdir
#include <sys/stat.h>
#include <Eigen/Core>
#include <cmath>     // 用于std::abs
#include <fstream>   // 用于文件输出
#include <iomanip>   // 用于输出格式控制
#include <iostream>  // 用于控制台输出
#include "config.h"
#include "geometry_mapping.h"  // 用于几何映射
#include "shape_functions.h"   // 用于形函数计算

bool VTKOutput::ensureDirectoryExists(const std::string& directory)
{
    struct stat info;
    if (stat(directory.c_str(), &info) != 0)
    {
        // 目录不存在，尝试创建
#ifdef _WIN32
        int result = _mkdir(directory.c_str());
#else
        int result = mkdir(directory.c_str(), 0755);
#endif
        if (result == 0)
        {
            std::cout << "创建输出目录: " << directory << std::endl;
            return true;
        }
        else
        {
            std::cerr << "无法创建目录: " << directory << std::endl;
            return false;
        }
    }
    return true;  // 目录已存在
}

// VTKOutput2D类的实现
void TriangleVTKOutput2D::outputNumericalSolution(const std::string& filename)
{
    // 确保results目录存在
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

    // 设置输出精度
    file << std::fixed << std::setprecision(6);

    // VTK文件头
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";

    // Piece信息
    file << "    <Piece NumberOfPoints=\"" << nodes_num << "\" NumberOfCells=\"" << elements_num
         << "\">\n";

    // === 输出节点坐标 ===
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

    // === 输出单元连接信息 ===
    file << "      <Cells>\n";

    // 连接信息
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << "\n";
    }
    file << "        </DataArray>\n";

    // 每个单元的偏移量（对于三角形，每个单元有3个节点）
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    // 单元类型（5 = 三角形）
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          5\n";  // VTK_TRIANGLE
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // === 输出节点数据 ===
    file << "      <PointData Scalars=\"numerical_solution\">\n";

    // 数值解
    file << "        <DataArray type=\"Float64\" Name=\"numerical_solution\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << solution_(i) << "\n";
    }
    file << "        </DataArray>\n";

    // 节点坐标（方便在ParaView中使用）
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

    // === 输出单元数据（可选） ===
    file << "      <CellData>\n";

    // 单元ID
    file << "        <DataArray type=\"Int32\" Name=\"element_id\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";

    file << "      </CellData>\n";

    // 文件结尾
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
    // 确保results目录存在
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
    int dimension = config_->getDimension();

    // 设置输出精度
    file << std::fixed << std::setprecision(6);

    // VTK文件头
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";

    // Piece信息
    file << "    <Piece NumberOfPoints=\"" << nodes_num << "\" NumberOfCells=\"" << elements_num
         << "\">\n";

    // === 输出节点坐标 ===
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << coordinates[i * 2] << " " << coordinates[i * 2 + 1] << " 0.0\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // === 输出单元连接信息 ===
    file << "      <Cells>\n";

    // 连接信息
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << "\n";
    }
    file << "        </DataArray>\n";

    // 每个单元的偏移量
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    // 单元类型（5 = 三角形）
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          5\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // === 输出节点数据 ===
    file << "      <PointData Scalars=\"exact_solution\">\n";

    // 解析解
    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        std::vector<double> coords = {coordinates[i * 2], coordinates[i * 2 + 1]};
        double exact_val = exact_func(coords);
        file << "          " << exact_val << "\n";
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";

    // === 输出单元数据 ===
    file << "      <CellData>\n";
    file << "        <DataArray type=\"Int32\" Name=\"element_id\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << i << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </CellData>\n";

    // 文件结尾
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();

    std::cout << "解析解VTK文件写入完成!" << std::endl;
}

void TriangleVTKOutput2D::outputDenseSamplingError(const std::string& filename,
                                                   std::shared_ptr<Config> config)
{
    // 确保results目录存在
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

    // 计算加密后的总点数和单元数
    std::vector<double> dense_points;     // 存储所有采样点的坐标
    std::vector<double> dense_numerical;  // 存储数值解
    std::vector<double> dense_exact;      // 存储精确解

    auto exact_func = [config](const std::vector<double>& coords) -> double
    { return config->exact_solution_u(coords); };

    // 遍历每个原始单元进行采样
    for (int e = 0; e < elements_num; ++e)
    {
        // 在每个三角形内部采样
        for (int i = 0; i < num_points_per_side_; ++i)
        {
            for (int j = 0; j < num_points_per_side_ - i; ++j)
            {
                // 在标准三角形参考单元内采样
                double xi = double(i) / (num_points_per_side_ - 1);
                double eta = double(j) / (num_points_per_side_ - 1);

                // 确保点在三角形内 (xi + eta <= 1)
                if (xi + eta <= 1.0)
                {
                    std::vector<double> coord_ref = {xi, eta};
                    std::vector<double> coord_phys;

                    // 获取单元节点坐标
                    std::vector<double> element_coords;
                    for (int local_node = 0; local_node < config_->getNodesPerElement();
                         ++local_node)
                    {
                        int global_node = connectivity[e][local_node];
                        element_coords.push_back(coordinates[global_node * dimension]);
                        element_coords.push_back(coordinates[global_node * dimension + 1]);
                    }

                    // 将参考坐标转换为物理坐标
                    auto mapping = GeometryMappingFactory::createMapping(element_coords, config);
                    mapping->mapToPhysical(coord_ref, coord_phys);

                    // 存储点坐标
                    dense_points.push_back(coord_phys[0]);
                    dense_points.push_back(coord_phys[1]);
                    dense_points.push_back(0.0);

                    // 计算该点的数值解
                    double numerical_val = evaluateNumericalSolutionAt(e, coord_ref, config);
                    dense_numerical.push_back(numerical_val);

                    // 计算该点的精确解
                    double exact_val = exact_func(coord_phys);
                    dense_exact.push_back(exact_val);
                }
            }
        }
    }

    int total_points = dense_points.size() / 3;  // 每个点有3个坐标分量
    std::cout << "总采样点数: " << total_points << " (原始节点数: " << nodes_num << ")"
              << std::endl;

    // 设置输出精度
    file << std::fixed << std::setprecision(6);

    // VTK文件头 - 使用点云格式，不生成单元连接
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" "
            "version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << total_points << "\" NumberOfCells=\"" << total_points
         << "\">\n";

    // === 输出节点坐标 ===
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" "
            "NumberOfComponents=\"3\" "
            "format=\"ascii\">\n";
    for (size_t i = 0; i < dense_points.size(); i += 3)
    {
        file << "          " << dense_points[i] << " " << dense_points[i + 1] << " "
             << dense_points[i + 2] << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // === 输出单元信息（每个点作为一个顶点单元）===
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" "
            "Name=\"connectivity\" "
            "format=\"ascii\">\n";
    for (int i = 0; i < total_points; ++i)
    {
        file << "          " << i << "\n";  // 每个点自成一个顶点单元
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" "
            "Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < total_points; ++i)
    {
        file << "          " << i + 1 << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"UInt8\" "
            "Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < total_points; ++i)
    {
        file << "          1\n";  // VTK_VERTEX = 1
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // === 输出节点数据 ===
    file << "      <PointData "
            "Scalars=\"numerical_solution\">\n";

    // 数值解
    file << "        <DataArray type=\"Float64\" "
            "Name=\"numerical_solution\" "
            "format=\"ascii\">\n";
    for (double val : dense_numerical)
    {
        file << "          " << val << "\n";
    }
    file << "        </DataArray>\n";

    // 精确解
    file << "        <DataArray type=\"Float64\" "
            "Name=\"exact_solution\" "
            "format=\"ascii\">\n";
    for (double val : dense_exact)
    {
        file << "          " << val << "\n";
    }
    file << "        </DataArray>\n";

    // 误差
    file << "        <DataArray type=\"Float64\" "
            "Name=\"error\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_numerical.size(); ++i)
    {
        double error = dense_numerical[i] - dense_exact[i];
        file << "          " << error << "\n";
    }
    file << "        </DataArray>\n";

    // 绝对误差
    file << "        <DataArray type=\"Float64\" "
            "Name=\"absolute_error\" "
            "format=\"ascii\">\n";
    for (size_t i = 0; i < dense_numerical.size(); ++i)
    {
        double abs_error = std::abs(dense_numerical[i] - dense_exact[i]);
        file << "          " << abs_error << "\n";
    }
    file << "        </DataArray>\n";

    // 相对误差
    file << "        <DataArray type=\"Float64\" "
            "Name=\"relative_error\" "
            "format=\"ascii\">\n";
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

    // 文件结尾
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
    int n = connectivity[element_index].size();  // 每个单元的节点数

    auto shape_func = ShapeFunctionFactory::createShapeFunction(config);

    for (int alpha = 0; alpha < n; ++alpha)  // 遍历单元内所有节点
    {
        int global_node_index = connectivity[element_index][alpha];
        double shape_val = shape_func->computeTrialFunction(alpha, coords_ref);
        numerical_solution += solution_(global_node_index) * shape_val;
    }
    return numerical_solution;
}

std::unique_ptr<VTKOutput> VTKOutputFactory::createVTKOutput(std::shared_ptr<Config> config,
                                                             const Eigen::VectorXd& solution)
{
    switch (config->getElementType())
    {
        case Config::ElementType::TRIANGLE:
            return std::make_unique<TriangleVTKOutput2D>(config, solution,
                                                         config->getSamplingPoints());
        case Config::ElementType::QUADRILATERAL:
            throw std::invalid_argument("Requested ElementType is not implemented in VTKOutput");
        default:
            throw std::invalid_argument("Unknown ElementType");
    }
}