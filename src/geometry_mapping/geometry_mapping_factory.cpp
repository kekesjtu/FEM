#include <stdexcept>
#include "config.h"
#include "geometry_mapping.h"


// ================ 工厂模式实现 ================

std::unique_ptr<GeometryMapping> GeometryMappingFactory::createMapping(
    const std::vector<double>& element_coords, std::shared_ptr<Config> config_, bool is_boundary)
{
    // 根据 is_boundary 选择体单元类型或边界单元类型
    int element_type = is_boundary ? config_->getBoundaryElementType() : config_->getElementType();
    int order = config_->getOrder();
    int embed_dim = config_->getDimension();

    switch (element_type)
    {
        case Config::ElementType::LINE:
            switch (order)
            {
                case 1:
                    return std::make_unique<LineLinearMapping>(element_coords, embed_dim);
                default:
                    throw std::invalid_argument("Unsupported order for Line element");
            }

        case Config::ElementType::TRIANGLE:
            switch (order)
            {
                case 1:
                    return std::make_unique<TriangleLinearMapping>(element_coords);
                default:
                    throw std::invalid_argument("Unsupported order for Triangle");
            }

        case Config::ElementType::QUADRILATERAL:
            throw std::invalid_argument(
                "Requested ElementType is not implemented in GeometryMappingFactory");

        case Config::ElementType::TETRAHEDRON:
            switch (order)
            {
                case 1:
                    return std::make_unique<TetrahedronLinearMapping>(element_coords, embed_dim);
                default:
                    throw std::invalid_argument("Unsupported order for Tetrahedron");
            }

        default:
            throw std::invalid_argument("Unknown ElementType");
    }
}
