#include "comsol_mesh_importer.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

ComsolMeshImporter::ComsolMeshImporter() : is_imported_(false), dimension_(0)
{
}

ComsolMeshImporter::~ComsolMeshImporter()
{
    clearData();
}

bool ComsolMeshImporter::importMesh(const std::string& filename)
{
    filename_ = filename;
    clearData();

    std::cout << "开始导入COMSOL网格文件: " << filename << std::endl;

    if (!parseFile(filename))
    {
        std::cerr << "解析网格文件失败" << std::endl;
        return false;
    }

    // 由于节点索引是连续的，无需重新映射
    is_imported_ = true;
    printImportStatistics();

    return true;
}

void ComsolMeshImporter::populateGlobalArrays()
{
    if (!is_imported_)
    {
        std::cerr << "错误：必须先成功导入网格才能填充全局数组" << std::endl;
        return;
    }
    // 此函数为旧版兼容代码，在新设计中应被废弃
    // 仅为示例保留
}

bool ComsolMeshImporter::parseFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }

    std::string line;

    // 优先解析维度
    while (std::getline(file, line))
    {
        if (line.find("sdim") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> dimension_; // 直接读取行首的数字
            break;
        }
    }
    
    // 重置文件流以重新搜索
    file.clear();
    file.seekg(0, std::ios::beg);

    // 跳过文件头，找到节点坐标部分
    while (std::getline(file, line))
    {
        if (line.find("# Mesh vertex coordinates") != std::string::npos)
        {
            break;
        }
    }

    if (file.eof())
    {
        std::cerr << "未找到节点坐标部分" << std::endl;
        return false;
    }

    // 解析节点坐标
    if (!parseNodes(file, line))
    {
        return false;
    }

    // 解析单元信息
    if (!parseElements(file, line))
    {
        return false;
    }

    file.close();
    return true;
}

bool ComsolMeshImporter::parseNodes(std::ifstream& file, std::string& line)
{
    while (std::getline(file, line))
    {
        trimString(line);

        // 跳过空行和注释
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        // 检查是否到达单元类型部分
        if (line.find("# number of element types") != std::string::npos ||
            line.find("element types") != std::string::npos)
        {
            break;
        }

        // 尝试解析坐标
        std::istringstream iss(line);
        double x, y, z; // 为3D做准备
        if (dimension_ == 2 && (iss >> x >> y))
        {
            nodes_.push_back({x, y});
        }
        else if (dimension_ == 3 && (iss >> x >> y >> z))
        {
            // 当前数据结构只支持2D，但为未来扩展留下接口
            nodes_.push_back({x, y}); // 暂时只存xy
        }
        else
        {
            // 如果不是坐标行，可能已经到达下一部分
            break;
        }
    }

    std::cout << "解析得到 " << nodes_.size() << " 个节点" << std::endl;
    return !nodes_.empty();
}
bool ComsolMeshImporter::parseElements(std::ifstream& file, std::string& line)
{
    std::cout << "开始解析单元信息..." << std::endl;

    // 寻找单元类型信息
    while (std::getline(file, line))
    {
        if (line.find("# Type #") != std::string::npos)
        {
            std::cout << "找到单元类型: " << line << std::endl;
            if (!parseElementType(file, line))
            {
                return false;
            }
        }
        // 检查是否结束
        if (file.eof())
        {
            break;
        }
    }

    std::cout << "单元解析完成，共解析到 " << triangular_elements_.size() << " 个三角形单元"
              << std::endl;
    return true;
}

bool ComsolMeshImporter::parseElementType(std::ifstream& file, std::string& line)
{
    std::cout << "parseElementType: 开始解析单元类型..." << std::endl;

    // 读取单元类型名称，跳过空行
    std::string type_line;
    do
    {
        if (!std::getline(file, type_line))
        {
            return false;
        }
        trimString(type_line);
    } while (type_line.empty());  // 跳过空行

    std::string element_type;
    if (type_line.find("tri") != std::string::npos)
    {
        element_type = "tri";
        std::cout << "parseElementType: 识别为三角形单元" << std::endl;
    }
    else if (type_line.find("edg") != std::string::npos)
    {
        element_type = "edg";
        std::cout << "parseElementType: 识别为边单元" << std::endl;
    }
    else if (type_line.find("vtx") != std::string::npos)
    {
        element_type = "vtx";
        std::cout << "parseElementType: 识别为顶点单元" << std::endl;
    }
    else
    {
        std::cout << "parseElementType: 未知类型: " << type_line << "，跳过此类型" << std::endl;
        // 跳过未知类型
        return skipToNextSection(file, line);
    }

    std::cout << "parseElementType: 开始解析 " << element_type << " 类型单元" << std::endl;

    // 读取每个单元的节点数，跳过空行
    int nodes_per_element;
    std::string nodes_line;
    do
    {
        if (!std::getline(file, nodes_line))
        {
            std::cout << "parseElementType: 无法读取节点数行" << std::endl;
            return false;
        }
        trimString(nodes_line);
    } while (nodes_line.empty());

    std::cout << "parseElementType: 节点数行: '" << nodes_line << "'" << std::endl;

    if (!(std::istringstream(nodes_line) >> nodes_per_element))
    {
        std::cout << "parseElementType: 无法解析节点数" << std::endl;
        return false;
    }
    std::cout << "parseElementType: 每个单元节点数: " << nodes_per_element << std::endl;

    // 读取单元数量，跳过空行
    int num_elements;
    std::string elements_line;
    do
    {
        if (!std::getline(file, elements_line))
        {
            std::cout << "parseElementType: 无法读取单元数量行" << std::endl;
            return false;
        }
        trimString(elements_line);
    } while (elements_line.empty());

    std::cout << "parseElementType: 单元数量行: '" << elements_line << "'" << std::endl;

    if (!(std::istringstream(elements_line) >> num_elements))
    {
        std::cout << "parseElementType: 无法解析单元数量" << std::endl;
        return false;
    }
    std::cout << "parseElementType: 单元数量: " << num_elements << std::endl;

    // 跳过 "# Elements" 行，可能包含空行
    std::string elements_header;
    do
    {
        if (!std::getline(file, elements_header))
        {
            std::cout << "parseElementType: 无法读取Elements标题行" << std::endl;
            return false;
        }
        trimString(elements_header);
    } while (elements_header.empty());

    std::cout << "parseElementType: Elements标题行: " << elements_header << std::endl;

    // 读取单元连接信息
    for (int i = 0; i < num_elements; ++i)
    {
        if (!std::getline(file, line))
        {
            std::cout << "parseElementType: 无法读取第 " << i << " 个单元" << std::endl;
            return false;
        }

        std::istringstream iss(line);
        std::vector<int> node_indices;
        int node_id;

        while (iss >> node_id)
        {
            node_indices.push_back(node_id);
        }

        if (static_cast<int>(node_indices.size()) != nodes_per_element)
        {
            std::cout << "parseElementType: 单元 " << i << " 节点数不匹配，期望 "
                      << nodes_per_element << "，实际 " << node_indices.size() << std::endl;
            continue;
        }

        // 直接存储节点连接，不需要复杂的结构体
        if (element_type == "tri")
        {
            triangular_elements_.push_back(node_indices);
        }
        else if (element_type == "edg")
        {
            edge_elements_.push_back(node_indices);
        }
        // 忽略点单元(vtx)
    }

    // 跳过几何实体索引部分
    std::string geo_line;
    // 先跳过空行
    while (std::getline(file, geo_line))
    {
        trimString(geo_line);
        if (!geo_line.empty())
            break;
    }

    // 如果找到"# number of geometric entity indices"行，跳过整个几何实体索引部分
    if (geo_line.find("# number of geometric entity indices") != std::string::npos)
    {
        // 跳过"# Geometric entity indices"标题行
        while (std::getline(file, geo_line))
        {
            trimString(geo_line);
            if (geo_line.find("# Geometric entity indices") != std::string::npos)
                break;
        }

        // 跳过所有几何实体索引数据（num_elements行）
        for (int i = 0; i < num_elements; ++i)
        {
            std::getline(file, geo_line);
        }
    }

    std::cout << "parseElementType: 成功解析 " << element_type << " 类型单元 " << num_elements
              << " 个" << std::endl;
    return true;
}

int ComsolMeshImporter::findElementContainingEdge(int node1, int node2) const
{
    // 在三角形单元中查找包含指定边的单元
    for (size_t i = 0; i < triangular_elements_.size(); ++i)
    {
        const auto& nodes = triangular_elements_[i];

        // 检查三角形的三条边
        if ((nodes[0] == node1 && nodes[1] == node2) || (nodes[0] == node2 && nodes[1] == node1) ||
            (nodes[1] == node1 && nodes[2] == node2) || (nodes[1] == node2 && nodes[2] == node1) ||
            (nodes[2] == node1 && nodes[0] == node2) || (nodes[2] == node2 && nodes[0] == node1))
        {
            return static_cast<int>(i);
        }
    }

    return 0;  // 如果未找到，返回第一个单元（作为默认值）
}

void ComsolMeshImporter::printImportStatistics() const
{
    std::cout << "\n=== 网格导入统计信息 ===" << std::endl;
    std::cout << "文件名: " << filename_ << std::endl;
    std::cout << "网格维度: " << dimension_ << std::endl;
    std::cout << "节点总数: " << nodes_.size() << std::endl;
    std::cout << "三角形单元数: " << triangular_elements_.size() << std::endl;
    std::cout << "边界边数: " << edge_elements_.size() << std::endl;
    std::cout << "========================\n" << std::endl;
}

bool ComsolMeshImporter::skipToNextSection(std::ifstream& file, std::string& line)
{
    // 跳过当前部分，寻找下一个 "# Type #" 或文件结束
    while (std::getline(file, line))
    {
        if (line.find("# Type #") != std::string::npos)
        {
            return true;
        }
    }
    return true;  // 到达文件末尾
}

void ComsolMeshImporter::clearData()
{
    nodes_.clear();
    triangular_elements_.clear();
    edge_elements_.clear();
    is_imported_ = false;
    dimension_ = 0;
}

void ComsolMeshImporter::trimString(std::string& str) const
{
    // 移除前导和尾随空格
    str.erase(str.begin(), std::find_if(str.begin(), str.end(),
                                        [](unsigned char ch) { return !std::isspace(ch); }));

    str.erase(
        std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) { return !std::isspace(ch); })
            .base(),
        str.end());
}

// --- 数据访问接口实现 ---

int ComsolMeshImporter::getDimension() const
{
    return dimension_;
}

const std::vector<std::pair<double, double>>& ComsolMeshImporter::getNodes() const
{
    return nodes_;
}

const std::vector<std::vector<int>>& ComsolMeshImporter::getTriangularElements() const
{
    return triangular_elements_;
}

int ComsolMeshImporter::getNodesNum() const
{
    return nodes_.size();
}

int ComsolMeshImporter::getElementsNum() const
{
    return triangular_elements_.size();
}

int ComsolMeshImporter::getNodesPerElement() const
{
    if (!triangular_elements_.empty())
    {
        return triangular_elements_[0].size();
    }
    return 0;
}
