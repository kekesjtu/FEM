#include "vtk_output.h"
#include <direct.h>  // 用于Windows的_mkdir
#include <sys/stat.h>
#include <cmath>                  // 用于std::abs
#include <fstream>                // 用于文件输出
#include <iomanip>                // 用于输出格式控制
#include <iostream>               // 用于控制台输出
#include "geometry_mapping_2d.h"  // 用于几何映射
#include "shape_functions_2d.h"   // 用于形函数计算

// 获取全局形函数实例
static auto* g_shapeFunction =
    ShapeFunctionFactory::getShapeFunction(ShapeFunctionFactory::ElementType::Triangle);
static auto* g_triangleShapeFunction = dynamic_cast<TriangleShapeFunction*>(g_shapeFunction);

/**
 * @brief 本地实现：计算单元内某点的数值解（使用mesh对象）
 */
static double evaluateNumericalSolutionWithMesh(int element_index, const Eigen::VectorXd& solution,
                                                double x_ref, double y_ref,
                                                std::shared_ptr<Mesh2D> mesh)
{
    double numerical_solution = 0.0;
    const auto& connectivity = mesh->getElementConnectivity();
    int n = connectivity[element_index].size();  // 每个单元的节点数

    for (int alpha = 0; alpha < n; ++alpha)  // 遍历单元内所有节点
    {
        int global_node_index = connectivity[element_index][alpha];
        numerical_solution += solution(global_node_index) *
                              g_triangleShapeFunction->computeTrialFunction2D(alpha, x_ref, y_ref);
    }
    return numerical_solution;
}

bool ensureDirectoryExists(const std::string& directory)
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
void VTKOutput2D::outputNumericalSolution(const std::string& filename,
                                          const Eigen::VectorXd& solution)
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

    const auto& coordinates = mesh2D_->getNodeCoordinates();
    const auto& connectivity = mesh2D_->getElementConnectivity();
    int num_nodes = mesh2D_->getNumNodes();
    int num_elements = mesh2D_->getNumElements();
    int spatial_dim = mesh2D_->getDimension();  // 直接获取维度整数值

    // 设置输出精度
    file << std::fixed << std::setprecision(6);

    // VTK文件头
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";

    // Piece信息
    file << "    <Piece NumberOfPoints=\"" << num_nodes << "\" NumberOfCells=\"" << num_elements
         << "\">\n";

    // === 输出节点坐标 ===
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        // 2D网格扩展为3D，z坐标设为0
        double x = coordinates[i * spatial_dim];
        double y = coordinates[i * spatial_dim + 1];
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
    for (int i = 0; i < num_elements; ++i)
    {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    // 单元类型（5 = 三角形）
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < num_elements; ++i)
    {
        file << "          5\n";  // VTK_TRIANGLE
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // === 输出节点数据 ===
    file << "      <PointData Scalars=\"numerical_solution\">\n";

    // 数值解
    file << "        <DataArray type=\"Float64\" Name=\"numerical_solution\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        file << "          " << solution(i) << "\n";
    }
    file << "        </DataArray>\n";

    // 节点坐标（方便在ParaView中使用）
    file << "        <DataArray type=\"Float64\" Name=\"x_coordinate\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        file << "          " << coordinates[i * spatial_dim] << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"y_coordinate\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        file << "          " << coordinates[i * spatial_dim + 1] << "\n";
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";

    // === 输出单元数据（可选） ===
    file << "      <CellData>\n";

    // 单元ID
    file << "        <DataArray type=\"Int32\" Name=\"element_id\" format=\"ascii\">\n";
    for (int i = 0; i < num_elements; ++i)
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

void VTKOutput2D::outputExactSolution(const std::string& filename,
                                      double (*exact_func)(double, double))
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

    const auto& coordinates = mesh2D_->getNodeCoordinates();
    const auto& connectivity = mesh2D_->getElementConnectivity();
    int num_nodes = mesh2D_->getNumNodes();
    int num_elements = mesh2D_->getNumElements();
    int spatial_dim = static_cast<int>(mesh2D_->getDimension());  // 获取空间维度

    // 设置输出精度
    file << std::fixed << std::setprecision(6);

    // VTK文件头
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";

    // Piece信息
    file << "    <Piece NumberOfPoints=\"" << num_nodes << "\" NumberOfCells=\"" << num_elements
         << "\">\n";

    // === 输出节点坐标 ===
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
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
    for (int i = 0; i < num_elements; ++i)
    {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    // 单元类型（5 = 三角形）
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < num_elements; ++i)
    {
        file << "          5\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // === 输出节点数据 ===
    file << "      <PointData Scalars=\"exact_solution\">\n";

    // 解析解
    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        double x = coordinates[i * 2];
        double y = coordinates[i * 2 + 1];
        double exact_val = exact_func(x, y);
        file << "          " << exact_val << "\n";
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";

    // === 输出单元数据 ===
    file << "      <CellData>\n";
    file << "        <DataArray type=\"Int32\" Name=\"element_id\" format=\"ascii\">\n";
    for (int i = 0; i < num_elements; ++i)
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

void VTKOutput2D::outputComparison(const std::string& filename, const Eigen::VectorXd& solution,
                                   double (*exact_func)(double, double))
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

    std::cout << "正在写入二维对比VTK文件: " << fullname << std::endl;

    const auto& coordinates = mesh2D_->getNodeCoordinates();
    const auto& connectivity = mesh2D_->getElementConnectivity();
    int num_nodes = mesh2D_->getNumNodes();
    int num_elements = mesh2D_->getNumElements();
    int spatial_dim = mesh2D_->getDimension();  // 直接获取维度整数值

    // 设置输出精度
    file << std::fixed << std::setprecision(6);

    // VTK文件头
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";

    // Piece信息
    file << "    <Piece NumberOfPoints=\"" << num_nodes << "\" NumberOfCells=\"" << num_elements
         << "\">\n";

    // === 输出节点坐标 ===
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        file << "          " << coordinates[i * 2] << " " << coordinates[i * 2 + 1] << " 0.0\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // === 输出单元连接信息 ===
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < num_elements; ++i)
    {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < num_elements; ++i)
    {
        file << "          5\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // === 输出节点数据 ===
    file << "      <PointData Scalars=\"numerical_solution\">\n";

    // 数值解
    file << "        <DataArray type=\"Float64\" Name=\"numerical_solution\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        file << "          " << solution(i) << "\n";
    }
    file << "        </DataArray>\n";

    // 精确解
    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        double x = coordinates[i * spatial_dim];
        double y = coordinates[i * spatial_dim + 1];
        double exact_val = exact_func(x, y);
        file << "          " << exact_val << "\n";
    }
    file << "        </DataArray>\n";

    // 误差
    file << "        <DataArray type=\"Float64\" Name=\"error\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        double x = coordinates[i * spatial_dim];
        double y = coordinates[i * spatial_dim + 1];
        double exact_val = exact_func(x, y);
        double error = solution(i) - exact_val;
        file << "          " << error << "\n";
    }
    file << "        </DataArray>\n";

    // 绝对误差
    file << "        <DataArray type=\"Float64\" Name=\"absolute_error\" format=\"ascii\">\n";
    for (int i = 0; i < num_nodes; ++i)
    {
        double x = coordinates[i * spatial_dim];
        double y = coordinates[i * spatial_dim + 1];
        double exact_val = exact_func(x, y);
        double abs_error = std::abs(solution(i) - exact_val);
        file << "          " << abs_error << "\n";
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";

    // === 输出单元数据 ===
    file << "      <CellData>\n";
    file << "        <DataArray type=\"Int32\" Name=\"element_id\" format=\"ascii\">\n";
    for (int i = 0; i < num_elements; ++i)
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

    std::cout << "对比VTK文件写入完成!" << std::endl;
}

void VTKOutput2D::outputDenseSamplingError(const std::string& filename,
                                           const Eigen::VectorXd& solution,
                                           double (*exact_func)(double, double),
                                           int num_points_per_side)
{
    // 确保results目录存在
    ensureDirectoryExists("results");

    std::string fullname = filename + "_dense.vtu";
    std::ofstream file(fullname);

    if (!file.is_open())
    {
        std::cerr << "无法创建加密采样VTK文件: " << fullname << std::endl;
        return;
    }

    std::cout << "正在写入二维加密采样VTK文件: " << fullname << std::endl;
    std::cout << "采样密度: 每个三角形单元 " << num_points_per_side << "×" << num_points_per_side
              << " 采样点" << std::endl;

    const auto& coordinates = mesh2D_->getNodeCoordinates();
    const auto& connectivity = mesh2D_->getElementConnectivity();
    int num_elements = mesh2D_->getNumElements();

    // 计算加密后的总点数和单元数
    std::vector<double> dense_points;     // 存储所有采样点的坐标 [x1,y1,0, x2,y2,0, ...]
    std::vector<double> dense_numerical;  // 存储数值解
    std::vector<double> dense_exact;      // 存储精确解

    // 遍历每个原始单元进行采样
    for (int e = 0; e < num_elements; ++e)
    {
        // 使用mesh2D的getElementNodes方法获取单元节点坐标
        std::vector<double> element_coords = mesh2D_->getElementNodes(e);

        GeometryMapping2D mapping(element_coords);

        // 在每个三角形内部采样
        for (int i = 0; i < num_points_per_side; ++i)
        {
            for (int j = 0; j < num_points_per_side - i; ++j)
            {
                // 在标准三角形参考单元内采样
                double x_ref = double(i) / (num_points_per_side - 1);
                double y_ref = double(j) / (num_points_per_side - 1);

                // 确保点在三角形内 (x_ref + y_ref <= 1)
                if (x_ref + y_ref <= 1.0)
                {
                    double x, y;
                    mapping.mapToPhysical(x_ref, y_ref, x, y);  // 转换为全局坐标

                    // 存储点坐标
                    dense_points.push_back(x);
                    dense_points.push_back(y);
                    dense_points.push_back(0.0);

                    // 计算该点的数值解
                    double numerical_val =
                        evaluateNumericalSolutionWithMesh(e, solution, x_ref, y_ref, mesh2D_);
                    dense_numerical.push_back(numerical_val);

                    // 计算该点的精确解
                    double exact_val = exact_func(x, y);
                    dense_exact.push_back(exact_val);
                }
            }
        }
    }

    int total_points = dense_points.size() / 3;  // 每个点有3个坐标分量
    std::cout << "总采样点数: " << total_points << " (原始节点数: " << mesh2D_->getNumNodes() << ")"
              << std::endl;

    // 设置输出精度
    file << std::fixed << std::setprecision(6);

    // VTK文件头 - 使用点云格式，不生成单元连接
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << total_points << "\" NumberOfCells=\"" << total_points
         << "\">\n";

    // === 输出节点坐标 ===
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_points.size(); i += 3)
    {
        file << "          " << dense_points[i] << " " << dense_points[i + 1] << " "
             << dense_points[i + 2] << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // === 输出单元信息（每个点作为一个顶点单元）===
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < total_points; ++i)
    {
        file << "          " << i << "\n";  // 每个点自成一个顶点单元
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
        file << "          1\n";  // VTK_VERTEX = 1
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // === 输出节点数据 ===
    file << "      <PointData Scalars=\"numerical_solution\">\n";

    // 数值解
    file << "        <DataArray type=\"Float64\" Name=\"numerical_solution\" format=\"ascii\">\n";
    for (double val : dense_numerical)
    {
        file << "          " << val << "\n";
    }
    file << "        </DataArray>\n";

    // 精确解
    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (double val : dense_exact)
    {
        file << "          " << val << "\n";
    }
    file << "        </DataArray>\n";

    // 误差
    file << "        <DataArray type=\"Float64\" Name=\"error\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_numerical.size(); ++i)
    {
        double error = dense_numerical[i] - dense_exact[i];
        file << "          " << error << "\n";
    }
    file << "        </DataArray>\n";

    // 绝对误差
    file << "        <DataArray type=\"Float64\" Name=\"absolute_error\" format=\"ascii\">\n";
    for (size_t i = 0; i < dense_numerical.size(); ++i)
    {
        double abs_error = std::abs(dense_numerical[i] - dense_exact[i]);
        file << "          " << abs_error << "\n";
    }
    file << "        </DataArray>\n";

    // 相对误差
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

    // 文件结尾
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();

    std::cout << "加密采样误差分析VTK文件写入完成!" << std::endl;
}