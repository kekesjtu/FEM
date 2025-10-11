# 测试CMakeLists.txt
cmake_minimum_required(VERSION 3.16)

# 查找测试框架（如果使用Google Test）
# find_package(GTest REQUIRED)

# 如果您想创建简单的测试
file(GLOB TEST_SOURCES "*.cpp")

foreach(test_source ${TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_executable(${test_name} ${test_source})
    
    # 链接到FEM库
    target_link_libraries(${test_name} FEM_lib)
    
    # 如果使用Google Test
    # target_link_libraries(${test_name} GTest::GTest GTest::Main)
    
    # 添加测试
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()