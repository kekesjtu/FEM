#include "vtk_output.h"
#include <direct.h>  // 用于Windows的_mkdir
#include <sys/stat.h>
#include <iostream>

// ============================================================================
// VTKOutput 基类实现
// ============================================================================

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

// ============================================================================
// VTKOutputFactory 工厂类实现
// ============================================================================

std::unique_ptr<VTKOutput> VTKOutputFactory::createVTKOutput(std::shared_ptr<Config> config,
                                                             const Eigen::VectorXd& solution)
{
    switch (config->getElementType())
    {
        case Config::ElementType::TRIANGLE:
            return std::make_unique<TriangleVTKOutput2D>(config, solution,
                                                         config->getSamplingPointsNum());
        case Config::ElementType::TETRAHEDRON:
            return std::make_unique<TetrahedronVTKOutput3D>(config, solution,
                                                            config->getSamplingPointsNum());
        case Config::ElementType::QUADRILATERAL:
            throw std::invalid_argument("Requested ElementType is not implemented in VTKOutput");
        default:
            throw std::invalid_argument("Unknown ElementType");
    }
}
