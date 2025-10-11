# CMake工具链文件 for MSYS2/MinGW64
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

# 设置编译器路径
set(CMAKE_C_COMPILER "E:/msys64/ucrt64/bin/gcc.exe")
set(CMAKE_CXX_COMPILER "E:/msys64/ucrt64/bin/g++.exe")

# 设置查找模式
set(CMAKE_FIND_ROOT_PATH "E:/msys64/ucrt64")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# 设置编译器标志
set(CMAKE_C_FLAGS_INIT "-Wall")
set(CMAKE_CXX_FLAGS_INIT "-Wall")

# 确保使用正确的运行时库
set(CMAKE_EXE_LINKER_FLAGS_INIT "-static-libgcc -static-libstdc++")
set(CMAKE_SHARED_LINKER_FLAGS_INIT "-static-libgcc -static-libstdc++")