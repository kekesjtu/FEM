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

// ============================================================================
// TetrahedronVTKOutput3D 类的实现
// ============================================================================

void TetrahedronVTKOutput3D::outputNumericalSolution(const std::string& filename)
{
    ensureDirectoryExists("results");
    
    std::string fullname = filename + ".vtu";
    std::ofstream file(fullname);
    
    if (!file.is_open())
    {
        std::cerr << "无法创建VTK文件: " << fullname << std::endl;
        return;
    }
    
    std::cout << "正在写入三维VTK文件: " << fullname << std::endl;
    
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
    file << "    <Piece NumberOfPoints=\"" << nodes_num << "\" NumberOfCells=\"" << elements_num << "\">\n";
    
    // 输出节点坐标
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        double x = coordinates[i * dimension];
        double y = coordinates[i * dimension + 1];
        double z = coordinates[i * dimension + 2];
        file << "          " << x << " " << y << " " << z << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    
    // 输出单元连接信息
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << " " << element[3] << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << (i + 1) * 4 << "\n";
    }
    file << "        </DataArray>\n";
    
    // 单元类型 10=VTK_TETRA
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          10\n";
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
    
    file << "        <DataArray type=\"Float64\" Name=\"z_coordinate\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << coordinates[i * dimension + 2] << "\n";
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

void TetrahedronVTKOutput3D::outputExactSolution(
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
    
    std::cout << "正在写入三维解析解VTK文件: " << fullname << std::endl;
    
    const auto& coordinates = config_->getNodeCoordinates();
    const auto& connectivity = config_->getElementConnectivity();
    int nodes_num = config_->getNodesNum();
    int elements_num = config_->getElementsNum();
    int dimension = config_->getDimension();
    
    file << std::fixed << std::setprecision(6);
    
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << nodes_num << "\" NumberOfCells=\"" << elements_num << "\">\n";
    
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        double x = coordinates[i * dimension];
        double y = coordinates[i * dimension + 1];
        double z = coordinates[i * dimension + 2];
        file << "          " << x << " " << y << " " << z << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << " " << element[3] << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << (i + 1) * 4 << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          10\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    
    file << "      <PointData Scalars=\"exact_solution\">\n";
    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        std::vector<double> coords = {coordinates[i * dimension], 
                                      coordinates[i * dimension + 1],
                                      coordinates[i * dimension + 2]};
        file << "          " << exact_func(coords) << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </PointData>\n";
    
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
    
    file.close();
    std::cout << "解析解VTK文件写入完成!" << std::endl;
}

void TetrahedronVTKOutput3D::outputDenseSamplingError(const std::string& filename,
                                                      std::shared_ptr<Config> config,
                                                      std::shared_ptr<ProblemSetup> problem)
{
    ensureDirectoryExists("results");
    
    std::string fullname = filename + ".vtu";
    std::ofstream file(fullname);
    
    if (!file.is_open())
    {
        std::cerr << "无法创建VTK文件: " << fullname << std::endl;
        return;
    }
    
    auto exact_func = [problem](const std::vector<double>& coords) -> double
    { return problem->exactSolutionU(coords); };
    
    if (!problem->hasExactSolution())
    {
        std::cerr << "警告: 没有提供精确解，无法输出误差文件" << std::endl;
        file.close();
        return;
    }
    
    std::cout << "正在写入三维误差VTK文件: " << fullname << std::endl;
    std::cout << "注意: 3D加密采样功能尚未完全实现，仅在节点处评估误差" << std::endl;
    
    const auto& coordinates = config_->getNodeCoordinates();
    const auto& connectivity = config_->getElementConnectivity();
    int nodes_num = config_->getNodesNum();
    int elements_num = config_->getElementsNum();
    int dimension = config_->getDimension();
    
    file << std::fixed << std::setprecision(6);
    
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << nodes_num << "\" NumberOfCells=\"" << elements_num << "\">\n";
    
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        double x = coordinates[i * dimension];
        double y = coordinates[i * dimension + 1];
        double z = coordinates[i * dimension + 2];
        file << "          " << x << " " << y << " " << z << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : connectivity)
    {
        file << "          " << element[0] << " " << element[1] << " " << element[2] << " " << element[3] << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          " << (i + 1) * 4 << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < elements_num; ++i)
    {
        file << "          10\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    
    file << "      <PointData Scalars=\"error\">\n";
    
    file << "        <DataArray type=\"Float64\" Name=\"numerical_solution\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        file << "          " << solution_(i) << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        std::vector<double> coords = {coordinates[i * dimension], 
                                      coordinates[i * dimension + 1],
                                      coordinates[i * dimension + 2]};
        file << "          " << exact_func(coords) << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"Float64\" Name=\"error\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        std::vector<double> coords = {coordinates[i * dimension], 
                                      coordinates[i * dimension + 1],
                                      coordinates[i * dimension + 2]};
        double error = solution_(i) - exact_func(coords);
        file << "          " << error << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "        <DataArray type=\"Float64\" Name=\"absolute_error\" format=\"ascii\">\n";
    for (int i = 0; i < nodes_num; ++i)
    {
        std::vector<double> coords = {coordinates[i * dimension], 
                                      coordinates[i * dimension + 1],
                                      coordinates[i * dimension + 2]};
        double abs_error = std::abs(solution_(i) - exact_func(coords));
        file << "          " << abs_error << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "      </PointData>\n";
    
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
    
    file.close();
    std::cout << "误差VTK文件写入完成!" << std::endl;
}

double TetrahedronVTKOutput3D::evaluateNumericalSolutionAt(
    int element_index, const std::vector<double>& coords_ref, std::shared_ptr<Config> config) const
{
    // TODO: 实现四面体单元内插值
    // 目前outputDenseSamplingError仅在节点处评估，所以暂不实现
    return 0.0;
}
