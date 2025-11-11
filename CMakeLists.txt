cmake_minimum_required(VERSION 3.16)

# 项目设置
project(FEM 
    VERSION 1.0.0
    DESCRIPTION "Finite Element Method Solver"
    LANGUAGES CXX
)

# 生成 compile_commands.json 用于 IntelliSense
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置构建类型默认为Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# 编译选项
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -fdiagnostics-color=auto)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
endif()

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 包含目录
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/third_party
)

# 查找Eigen库（如果使用find_package方式）
# find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 按模块收集源文件 - 更清晰的组织结构
set(CORE_SOURCES
    "${CMAKE_SOURCE_DIR}/src/core/config.cpp"
    "${CMAKE_SOURCE_DIR}/src/core/comsol_mesh_importer.cpp"
    "${CMAKE_SOURCE_DIR}/src/core/problem_setup.cpp"
    "${CMAKE_SOURCE_DIR}/src/core/material.cpp"
)

set(SHAPE_FUNCTION_SOURCES
    "${CMAKE_SOURCE_DIR}/src/shape_functions/shape_functions_1d.cpp"
    "${CMAKE_SOURCE_DIR}/src/shape_functions/shape_functions_2d.cpp"
    "${CMAKE_SOURCE_DIR}/src/shape_functions/shape_functions_3d.cpp"
    "${CMAKE_SOURCE_DIR}/src/shape_functions/shape_functions_factory.cpp"
)

set(GEOMETRY_MAPPING_SOURCES
    "${CMAKE_SOURCE_DIR}/src/geometry_mapping/geometry_mapping_1d.cpp"
    "${CMAKE_SOURCE_DIR}/src/geometry_mapping/geometry_mapping_2d.cpp"
    "${CMAKE_SOURCE_DIR}/src/geometry_mapping/geometry_mapping_3d.cpp"
    "${CMAKE_SOURCE_DIR}/src/geometry_mapping/geometry_mapping_factory.cpp"
)

set(GAUSS_QUADRATURE_SOURCES
    "${CMAKE_SOURCE_DIR}/src/gauss_quadrature/gauss_quadrature_1d.cpp"
    "${CMAKE_SOURCE_DIR}/src/gauss_quadrature/gauss_quadrature_2d.cpp"
    "${CMAKE_SOURCE_DIR}/src/gauss_quadrature/gauss_quadrature_3d.cpp"
    "${CMAKE_SOURCE_DIR}/src/gauss_quadrature/gauss_quadrature_factory.cpp"
)

set(SOLVER_SOURCES
    "${CMAKE_SOURCE_DIR}/src/solvers/fem_solver.cpp"
    "${CMAKE_SOURCE_DIR}/src/solvers/electrothermal_solver.cpp"
)

set(VTK_OUTPUT_SOURCES
    "${CMAKE_SOURCE_DIR}/src/output/vtk_output/vtk_output_base.cpp"
    "${CMAKE_SOURCE_DIR}/src/output/vtk_output/vtk_output_2d.cpp"
    "${CMAKE_SOURCE_DIR}/src/output/vtk_output/vtk_output_3d.cpp"
)

set(OUTPUT_SOURCES
    ${VTK_OUTPUT_SOURCES}
    "${CMAKE_SOURCE_DIR}/src/output/error_analysis.cpp"
)

# 主程序源文件
set(MAIN_SOURCE
    "${CMAKE_SOURCE_DIR}/src/main.cpp"
)

# 组合所有库源文件
set(FEM_LIB_SOURCES
    ${CORE_SOURCES}
    ${SHAPE_FUNCTION_SOURCES}
    ${GEOMETRY_MAPPING_SOURCES}
    ${GAUSS_QUADRATURE_SOURCES}
    ${SOLVER_SOURCES}
    ${OUTPUT_SOURCES}
)

# 可执行文件源文件 = 库源文件 + main
set(FEM_SOURCES
    ${FEM_LIB_SOURCES}
    ${MAIN_SOURCE}
)

# 收集头文件
file(GLOB_RECURSE FEM_HEADERS
    "${CMAKE_SOURCE_DIR}/include/*.h"
)

# 创建可执行文件
add_executable(FEM_solution ${FEM_SOURCES} ${FEM_HEADERS})

# 链接Eigen库（如果使用find_package方式）
# target_link_libraries(FEM_solution Eigen3::Eigen)

# 如果使用本地Eigen（当前方式）
target_include_directories(FEM_solution PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/third_party
)

# 设置目标属性
set_target_properties(FEM_solution PROPERTIES
    OUTPUT_NAME "FEM_solution"
    # DEBUG_POSTFIX "_debug"  # 注释掉以保持统一的文件名
)

# 安装规则
install(TARGETS FEM_solution
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# 可选：创建一个库版本（如果您想将FEM求解器作为库使用）
# 注意：FEM_LIB_SOURCES 已经在上面定义，不包含 main.cpp
add_library(FEM_lib STATIC ${FEM_LIB_SOURCES})
target_include_directories(FEM_lib PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/third_party
)

# 为库设置输出名称
set_target_properties(FEM_lib PROPERTIES
    OUTPUT_NAME "fem"
    # ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# 如果您有测试
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# 打印一些配置信息
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")
message(STATUS "Binary directory: ${CMAKE_BINARY_DIR}")