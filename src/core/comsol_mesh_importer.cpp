#include "comsol_mesh_importer.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

ComsolMeshImporter::ComsolMeshImporter()
    : is_imported_(false), dimension_(0), nodes_num_(0), elements_num_(0), nodes_num_per_element_(0)
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

bool ComsolMeshImporter::parseFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }

    // 1. 解析空间维度
    if (!parseDimension(file))
    {
        std::cerr << "解析维度信息失败" << std::endl;
        return false;
    }
    std::cout << "空间维度: " << dimension_ << "D" << std::endl;

    // 2. 解析节点坐标
    if (!parseNodes(file))
    {
        std::cerr << "解析节点坐标失败" << std::endl;
        return false;
    }

    // 3. 解析单元信息
    if (!parseElements(file))
    {
        std::cerr << "解析单元信息失败" << std::endl;
        return false;
    }

    file.close();
    return true;
}

bool ComsolMeshImporter::parseDimension(std::ifstream& file)
{
    std::string line;
    while (std::getline(file, line))
    {
        if (line.find("sdim") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> dimension_;
            return dimension_ > 0 && dimension_ <= 3;
        }
    }
    return false;
}

bool ComsolMeshImporter::parseNodes(std::ifstream& file)
{
    std::string line;

    // 找到节点坐标部分
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

    std::vector<std::vector<double>> temp_nodes;  // 临时存储 [node_id][coords]

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

        // 尝试解析坐标 - 支持1D/2D/3D
        std::istringstream iss(line);
        std::vector<double> coord;
        double value;

        while (iss >> value)
        {
            coord.push_back(value);
        }

        // 验证坐标维度是否匹配
        if (!coord.empty() && static_cast<int>(coord.size()) == dimension_)
        {
            temp_nodes.push_back(coord);
        }
        else if (!coord.empty())
        {
            // 如果不是坐标行，可能已经到达下一部分
            break;
        }
    }

    // 转换为展平数组 [node_id * dimension + coord_id]
    nodes_num_ = static_cast<int>(temp_nodes.size());
    node_coordinates_.clear();
    node_coordinates_.reserve(nodes_num_ * dimension_);

    for (const auto& node : temp_nodes)
    {
        for (double coord : node)
        {
            node_coordinates_.push_back(coord);
        }
    }

    std::cout << "解析得到 " << nodes_num_ << " 个 " << dimension_ << "D 节点" << std::endl;
    return nodes_num_ > 0;
}

bool ComsolMeshImporter::parseElements(std::ifstream& file)
{
    std::cout << "开始解析单元信息..." << std::endl;

    std::string line;
    // 寻找单元类型信息
    while (std::getline(file, line))
    {
        if (line.find("# Type #") != std::string::npos)
        {
            std::cout << "找到单元类型: " << line << std::endl;

            // 1. 解析单元类型(tri/tet/edg/vtx)
            std::string element_type;
            if (!parseElementType(file, element_type))
            {
                // 需要跳过的类型(vtx, 3D中的edg, 未知类型)
                // 读取并跳过: 节点数行、单元数行、连接矩阵
                int nodes_per_element, num_elements;
                if (parseNodesPerElement(file, nodes_per_element) &&
                    parseNumElements(file, num_elements))
                {
                    // 跳过连接矩阵
                    for (int i = 0; i < num_elements; ++i)
                    {
                        std::getline(file, line);
                    }
                }
                continue;
            }

            // 2. 解析每个单元的节点数
            int nodes_per_element;
            if (!parseNodesPerElement(file, nodes_per_element))
            {
                return false;
            }

            // 3. 解析单元数量
            int num_elements;
            if (!parseNumElements(file, num_elements))
            {
                return false;
            }

            // 4. 解析单元连接矩阵
            if (!parseElementConnectivity(file, num_elements, nodes_per_element, element_type))
            {
                return false;
            }

            std::cout << "成功解析 " << element_type << " 类型单元 " << num_elements << " 个"
                      << std::endl;
        }

        // 检查是否结束
        if (file.eof())
        {
            break;
        }
    }

    // 设置体单元统计信息
    elements_num_ = static_cast<int>(element_connectivity_.size());
    if (elements_num_ > 0)
    {
        nodes_num_per_element_ = static_cast<int>(element_connectivity_[0].size());
    }

    std::cout << "单元解析完成，共解析到 " << elements_num_ << " 个体单元，每个 "
              << nodes_num_per_element_ << " 个节点" << std::endl;
    std::cout << "边界单元数: " << boundary_elements_.size() << std::endl;

    return elements_num_ > 0;
}

bool ComsolMeshImporter::parseElementType(std::ifstream& file, std::string& element_type)
{
    // 读取单元类型名称，跳过空行
    std::string type_line;
    do
    {
        if (!std::getline(file, type_line))
        {
            return false;
        }
        trimString(type_line);
    } while (type_line.empty());

    // 识别单元类型
    // 支持的类型: tri(三角形), tet(四面体), quad(四边形), hex(六面体),
    //           prism(棱柱), pyramid(金字塔), edg(边), vtx(顶点)

    if (type_line.find("tri") != std::string::npos)
    {
        element_type = "tri";
        std::cout << "识别为三角形单元" << std::endl;
        return true;
    }
    else if (type_line.find("tet") != std::string::npos)
    {
        element_type = "tet";
        std::cout << "识别为四面体单元" << std::endl;
        return true;
    }
    else if (type_line.find("quad") != std::string::npos)
    {
        element_type = "quad";
        std::cout << "识别为四边形单元" << std::endl;
        return true;
    }
    else if (type_line.find("hex") != std::string::npos)
    {
        element_type = "hex";
        std::cout << "识别为六面体单元" << std::endl;
        return true;
    }
    else if (type_line.find("prism") != std::string::npos)
    {
        element_type = "prism";
        std::cout << "识别为棱柱单元" << std::endl;
        return true;
    }
    else if (type_line.find("pyramid") != std::string::npos ||
             type_line.find("pyr") != std::string::npos)
    {
        element_type = "pyramid";
        std::cout << "识别为金字塔单元" << std::endl;
        return true;
    }
    else if (type_line.find("edg") != std::string::npos)
    {
        element_type = "edg";
        // edg在不同维度有不同意义:
        // 1D: 体单元(线段)
        // 2D: 边界单元(边)
        // 3D: 忽略(低两维)
        if (dimension_ == 3)
        {
            std::cout << "3D网格中忽略边单元(低两维)" << std::endl;
            return false;  // 跳过
        }
        else if (dimension_ == 2)
        {
            std::cout << "识别为边单元(2D边界)" << std::endl;
        }
        else  // dimension_ == 1
        {
            std::cout << "识别为边单元(1D体单元)" << std::endl;
        }
        return true;
    }
    else if (type_line.find("vtx") != std::string::npos)
    {
        element_type = "vtx";
        // vtx在不同维度有不同意义:
        // 1D: 边界单元(端点)
        // 2D: 忽略(低两维)
        // 3D: 忽略(低三维)
        if (dimension_ == 1)
        {
            std::cout << "识别为顶点单元(1D边界)" << std::endl;
            return true;
        }
        else
        {
            std::cout << "忽略顶点单元(低两维或更多)" << std::endl;
            return false;  // 跳过
        }
    }
    else
    {
        std::cout << "未知类型: " << type_line << "，跳过" << std::endl;
        return false;
    }
}

void ComsolMeshImporter::printImportStatistics() const
{
    std::cout << "\n=== 网格导入统计信息 ===" << std::endl;
    std::cout << "文件名: " << filename_ << std::endl;
    std::cout << "空间维度: " << dimension_ << "D" << std::endl;
    std::cout << "节点总数: " << nodes_num_ << std::endl;
    std::cout << "体单元数: " << elements_num_ << std::endl;
    std::cout << "每个体单元节点数: " << nodes_num_per_element_ << std::endl;
    std::cout << "边界单元数: " << boundary_elements_.size() << std::endl;

    if (!boundary_elements_.empty())
    {
        std::cout << "每个边界单元节点数: " << boundary_elements_[0].size() << std::endl;
    }

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
    node_coordinates_.clear();
    element_connectivity_.clear();
    boundary_elements_.clear();
    is_imported_ = false;
    dimension_ = 0;
    nodes_num_ = 0;
    elements_num_ = 0;
    nodes_num_per_element_ = 0;
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

// --- 通用数据访问接口实现 (对应Config类成员) ---

int ComsolMeshImporter::getDimension() const
{
    return dimension_;
}

int ComsolMeshImporter::getNodesNum() const
{
    return nodes_num_;
}

int ComsolMeshImporter::getElementsNum() const
{
    return elements_num_;
}

int ComsolMeshImporter::getNodesNumPerElement() const
{
    return nodes_num_per_element_;
}

const std::vector<double>& ComsolMeshImporter::getNodeCoordinates() const
{
    return node_coordinates_;
}

const std::vector<std::vector<int>>& ComsolMeshImporter::getElementConnectivity() const
{
    return element_connectivity_;
}

const std::vector<std::vector<int>>& ComsolMeshImporter::getBoundaryElements() const
{
    return boundary_elements_;
}

// --- 辅助解析方法实现 ---

bool ComsolMeshImporter::parseNodesPerElement(std::ifstream& file, int& nodes_per_element)
{
    std::string line;

    // 跳过空行
    do
    {
        if (!std::getline(file, line))
        {
            std::cerr << "无法读取每个单元的节点数" << std::endl;
            return false;
        }
        trimString(line);
    } while (line.empty());

    // 解析节点数
    if (!(std::istringstream(line) >> nodes_per_element))
    {
        std::cerr << "无法解析节点数: '" << line << "'" << std::endl;
        return false;
    }

    std::cout << "每个单元节点数: " << nodes_per_element << std::endl;
    return true;
}

bool ComsolMeshImporter::parseNumElements(std::ifstream& file, int& num_elements)
{
    std::string line;

    // 跳过空行
    do
    {
        if (!std::getline(file, line))
        {
            std::cerr << "无法读取单元数量" << std::endl;
            return false;
        }
        trimString(line);
    } while (line.empty());

    // 解析单元数量
    if (!(std::istringstream(line) >> num_elements))
    {
        std::cerr << "无法解析单元数量: '" << line << "'" << std::endl;
        return false;
    }

    std::cout << "单元数量: " << num_elements << std::endl;
    return true;
}

bool ComsolMeshImporter::parseElementConnectivity(std::ifstream& file, int num_elements,
                                                  int nodes_per_element,
                                                  const std::string& element_type)
{
    std::string line;

    // 跳过 "# Elements" 标题行
    skipEmptyLines(file, line);

    // 读取单元连接矩阵
    for (int i = 0; i < num_elements; ++i)
    {
        if (!std::getline(file, line))
        {
            std::cerr << "无法读取第 " << i << " 个单元" << std::endl;
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
            std::cerr << "单元 " << i << " 节点数不匹配，期望 " << nodes_per_element << "，实际 "
                      << node_indices.size() << std::endl;
            continue;
        }

        // 使用单元分类判断模块
        if (isVolumeElement(element_type, nodes_per_element))
        {
            element_connectivity_.push_back(node_indices);
        }
        else if (isBoundaryElement(element_type, nodes_per_element))
        {
            boundary_elements_.push_back(node_indices);
        }
        // 否则忽略 (已在parseElementType中处理)
    }

    // 跳过几何实体索引部分
    std::string geo_line;
    while (std::getline(file, geo_line))
    {
        trimString(geo_line);
        if (!geo_line.empty())
            break;
    }

    // 如果找到几何实体索引部分，跳过它
    if (geo_line.find("# number of geometric entity indices") != std::string::npos)
    {
        // 跳过"# Geometric entity indices"标题行
        while (std::getline(file, geo_line))
        {
            trimString(geo_line);
            if (geo_line.find("# Geometric entity indices") != std::string::npos)
                break;
        }

        // 跳过所有几何实体索引数据
        for (int i = 0; i < num_elements; ++i)
        {
            std::getline(file, geo_line);
        }
    }

    return true;
}

void ComsolMeshImporter::skipEmptyLines(std::ifstream& file, std::string& line)
{
    do
    {
        if (!std::getline(file, line))
        {
            return;
        }
        trimString(line);
    } while (line.empty());
}

// ==================== 单元分类判断模块 ====================

/**
 * @brief 对单元进行分类:体单元、边界单元或忽略
 *
 * 分类规则:
 * - 体单元(dimension维): 与空间维度相同
 *   * 1D: 线单元 (2节点)
 *   * 2D: 三角形/四边形 (3/4节点)
 *   * 3D: 四面体/六面体 (4/8节点)
 *
 * - 边界单元(dimension-1维): 比空间维度低一维
 *   * 1D: 点 (1节点)
 *   * 2D: 边 (2节点)
 *   * 3D: 三角形/四边形面 (3/4节点)
 *
 * - 忽略单元(dimension-2维或更低): 比空间维度低两维或更多
 *   * 2D: 点单元 (0维)
 *   * 3D: 边单元 (1维), 点单元 (0维)
 *
 * @param element_type 单元类型字符串 ("tri", "tet", "edg", "vtx"等)
 * @param nodes_per_element 每个单元的节点数
 * @return ElementClassification 单元分类
 */
ComsolMeshImporter::ElementClassification ComsolMeshImporter::classifyElement(
    const std::string& element_type, int nodes_per_element) const
{
    // 1D空间
    if (dimension_ == 1)
    {
        if (element_type == "edg" && nodes_per_element == 2)
        {
            return ElementClassification::VOLUME;  // 1D体单元是线段
        }
        else if (element_type == "vtx" && nodes_per_element == 1)
        {
            return ElementClassification::BOUNDARY;  // 1D边界是点
        }
        else
        {
            return ElementClassification::IGNORED;
        }
    }
    // 2D空间
    else if (dimension_ == 2)
    {
        if ((element_type == "tri" && nodes_per_element == 3) ||
            (element_type == "quad" && nodes_per_element == 4))
        {
            return ElementClassification::VOLUME;  // 2D体单元是三角形或四边形
        }
        else if (element_type == "edg" && nodes_per_element == 2)
        {
            return ElementClassification::BOUNDARY;  // 2D边界是边
        }
        else if (element_type == "vtx" && nodes_per_element == 1)
        {
            return ElementClassification::IGNORED;  // 2D中点是低两维,忽略
        }
        else
        {
            return ElementClassification::IGNORED;
        }
    }
    // 3D空间
    else if (dimension_ == 3)
    {
        if ((element_type == "tet" && nodes_per_element == 4) ||
            (element_type == "hex" && nodes_per_element == 8) ||
            (element_type == "prism" && nodes_per_element == 6) ||
            (element_type == "pyramid" && nodes_per_element == 5))
        {
            return ElementClassification::VOLUME;  // 3D体单元是四面体、六面体等
        }
        else if ((element_type == "tri" && nodes_per_element == 3) ||
                 (element_type == "quad" && nodes_per_element == 4))
        {
            return ElementClassification::BOUNDARY;  // 3D边界是三角形或四边形面
        }
        else if (element_type == "edg" && nodes_per_element == 2)
        {
            return ElementClassification::IGNORED;  // 3D中边是低两维,忽略
        }
        else if (element_type == "vtx" && nodes_per_element == 1)
        {
            return ElementClassification::IGNORED;  // 3D中点是低三维,忽略
        }
        else
        {
            return ElementClassification::IGNORED;
        }
    }
    else
    {
        return ElementClassification::IGNORED;  // 未知维度
    }
}

/**
 * @brief 判断是否为体单元
 */
bool ComsolMeshImporter::isVolumeElement(const std::string& element_type,
                                         int nodes_per_element) const
{
    return classifyElement(element_type, nodes_per_element) == ElementClassification::VOLUME;
}

/**
 * @brief 判断是否为边界单元
 */
bool ComsolMeshImporter::isBoundaryElement(const std::string& element_type,
                                           int nodes_per_element) const
{
    return classifyElement(element_type, nodes_per_element) == ElementClassification::BOUNDARY;
}
