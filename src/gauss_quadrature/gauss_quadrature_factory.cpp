#include <stdexcept>
#include "config.h"
#include "gauss_quadrature.h"


// ============================================================================
// GaussPointFactory 实现
// ============================================================================

std::unique_ptr<GaussPoint> GaussPointFactory::createGaussPoint(Config::ElementType elementType,
                                                                int numPoints)
{
    switch (elementType)
    {
        case Config::ElementType::LINE:
            switch (numPoints)
            {
                case 2:
                    return std::make_unique<LineGauss2Point>();
                case 3:
                    return std::make_unique<LineGauss3Point>();
                default:
                    throw std::invalid_argument("Unsupported number of Gauss points for Line");
            }

        case Config::ElementType::TRIANGLE:
            switch (numPoints)
            {
                case 1:
                    return std::make_unique<TriangleGauss1Point>();
                case 3:
                    return std::make_unique<TriangleGauss3Point>();
                case 4:
                    return std::make_unique<TriangleGauss4Point>();
                default:
                    throw std::invalid_argument("Unsupported number of Gauss points for Triangle");
            }

        case Config::ElementType::QUADRILATERAL:
            // 目前未实现其他单元类型的高斯点
            throw std::invalid_argument(
                "Requested ElementType is not implemented in GaussPointFactory");

        case Config::ElementType::TETRAHEDRON:
            switch (numPoints)
            {
                case 1:
                    return std::make_unique<TetrahedronGauss1Point>();
                case 4:
                    return std::make_unique<TetrahedronGauss4Point>();
                default:
                    throw std::invalid_argument(
                        "Unsupported number of Gauss points for Tetrahedron");
            }

        default:
            throw std::invalid_argument("Unknown ElementType");
    }
}
