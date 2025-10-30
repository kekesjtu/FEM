#include <stdexcept>
#include "config.h"
#include "shape_functions.h"


// ================ 工厂模式实现 ================

std::unique_ptr<ShapeFunction> ShapeFunctionFactory::createShapeFunction(
    std::shared_ptr<Config> config_, bool is_boundary)
{
    // 根据 is_boundary 选择体单元类型或边界单元类型
    int element_type = is_boundary ? config_->getBoundaryElementType() : config_->getElementType();
    int order = config_->getOrder();

    switch (element_type)
    {
        case Config::ElementType::LINE:
            switch (order)
            {
                case 1:
                    return std::make_unique<LineLinearShapeFunction>();
                default:
                    throw std::invalid_argument("Unsupported order for Line element");
            }

        case Config::ElementType::TRIANGLE:
            switch (order)
            {
                case 1:
                    return std::make_unique<TriangleLinearShapeFunction>();
                default:
                    throw std::invalid_argument("Unsupported order for Triangle");
            }

        case Config::ElementType::QUADRILATERAL:
            throw std::invalid_argument(
                "Requested ElementType is not implemented in ShapeFunctionFactory");

        case Config::ElementType::TETRAHEDRON:
            switch (order)
            {
                case 1:
                    return std::make_unique<TetrahedronLinearShapeFunction>();
                default:
                    throw std::invalid_argument("Unsupported order for Tetrahedron");
            }

        default:
            throw std::invalid_argument("Unknown ElementType");
    }
}
