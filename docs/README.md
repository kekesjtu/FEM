# 有限元求解器 (FEM)

基于C++和Eigen库的高性能二维有限元求解器，专门用于求解椭圆型偏微分方程。采用线性三角形单元，支持复杂几何域、多种边界条件和精确的误差分析。

## 核心功能

- 🧮 **椭圆型PDE求解**：泊松方程 `-∇·(c∇u) = f`
- 🔺 **线性三角形单元**：标准有限元离散化
- 📐 **网格支持**：COMSOL网格导入器
- 🎯 **边界条件**：Dirichlet、Neumann、Robin边界条件
- 📊 **误差分析**：L∞、L2、H1范数计算
- 🎨 **可视化**：VTK输出，支持ParaView
- ⚙️ **配置管理**：统一的配置类管理计算参数
- 🔧 **CMake构建**：现代化的构建系统

## 项目结构

```
FEM/
├── 核心模块/
│   ├── main.cpp                      # 主程序入口
│   ├── fem_solver.h/.cpp             # 核心FEM求解器
│   ├── config.h/.cpp                 # 配置管理类
│   └── error_analysis.h/.cpp         # 误差分析模块
├── 数值计算/
│   ├── shape_functions.h/.cpp        # 三角形形函数
│   ├── gauss_quadrature.h/.cpp       # 高斯积分
│   └── geometry_mapping.h/.cpp       # 几何坐标映射
├── 网格与IO/
│   ├── comsol_mesh_importer.h/.cpp   # COMSOL网格导入
│   ├── vtk_output.h/.cpp             # VTK可视化输出
│   ├── circle_mesh.mphtxt            # 网格数据文件
│   └── circle_mesh2.mphtxt
├── 构建系统/
│   ├── CMakeLists.txt                # CMake构建配置
│   ├── CMakePresets.json             # CMake预设
│   └── cmake/mingw-toolchain.cmake   # MinGW工具链配置
├── 依赖库/
│   └── third_party/Eigen/            # 线性代数库
├── 输出结果/
│   └── results/                      # VTK可视化文件和结果图像
└── 文档/
    └── docs/                         # 技术文档和使用教程
```

## 算法流程

1. **预处理** - 网格导入/生成，数据结构初始化
2. **矩阵组装** - 遍历单元计算局部矩阵，组装全局系统
3. **边界条件** - 识别边界节点，修改系统矩阵  
4. **求解** - 稀疏LU分解求解线性系统
5. **后处理** - 误差计算，VTK文件输出

## 快速开始

### 环境要求

- C++17 编译器 (g++/MSVC)
- CMake 3.16+
- Eigen 3.x 线性代数库 (已包含在third_party目录)

### 编译运行

**VS Code环境**（推荐）:
1. 按 `Ctrl+Shift+P` → "Tasks: Run Task" → "CMake Build"
2. 运行: "Run CMake FEM Solution"
3. 或使用一键构建运行: "Build and Run"

**命令行编译**:
```bash
# 配置项目
cmake -B build -S . -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug

# 构建项目
cmake --build build --config Debug

# 运行程序
./build/bin/FEM_solution.exe
```

### 项目配置

项目包含完整的VS Code配置文件：
- `tasks.json`: 构建和运行任务配置
- `CMakePresets.json`: CMake预设配置
- `FEM.code-workspace`: 工作区配置

详细的VS Code使用教程请参考 `docs/vscode使用教程.md`

### 输出示例

```
步骤 1: 网格生成和预处理完成。
步骤 2: 全局矩阵组装完成。
步骤 3: 边界条件施加完成。
步骤 4: 线性方程组求解完成。
步骤 5: 后处理完成。

=== 详细误差分析 ===
L∞范数:  2.456e-03
L2范数:   8.234e-04  
H1半范数: 3.167e-02

✓ VTK文件已生成 (results/)
```

## ParaView可视化

输出的VTK文件位于 `results/` 目录：
- `comparison.vtu` - 完整对比数据（推荐）
- `numerical_solution.vtu` - 数值解
- `exact_solution.vtu` - 精确解  
- `solution.vtu` - 基本解数据

可视化步骤：
1. 打开ParaView，加载 `results/comparison.vtu`
2. 选择显示字段: `numerical_solution`、`exact_solution`、`error`
3. 调整色标和等值线显示
4. 对比分析数值解与精确解的差异

项目还生成PNG格式的结果图像，方便快速查看。

## 技术特性

- **高性能**: 稀疏矩阵存储，优化的矩阵组装算法
- **模块化**: 清晰的代码结构，便于扩展和维护  
- **工业级**: 支持COMSOL网格导入，复杂几何处理
- **精确验证**: 完整误差分析，收敛性验证
- **可视化**: VTK输出，支持ParaView专业可视化
- **现代化构建**: CMake构建系统，VS Code集成
- **配置灵活**: 统一的Config类管理问题参数和边界条件

## 主要模块

### 核心求解器 (`fem_solver`)
- 全局矩阵组装
- 边界条件处理
- 线性系统求解
- 预处理和后处理

### 配置管理 (`config`)
- 问题参数配置
- 边界条件定义（Dirichlet、Neumann、Robin）
- 几何和数值参数管理

### 误差分析 (`error_analysis`)
- L∞、L2、H1范数计算
- 收敛性分析
- 精确解比较

### 网格处理 (`comsol_mesh_importer`)
- COMSOL `.mphtxt` 格式导入
- 节点和单元数据解析
- 边界识别

### 可视化输出 (`vtk_output`)
- VTK格式文件生成
- 数值解和精确解输出
- 误差分布可视化

## 依赖要求

- C++17 编译器 (g++/MSVC)
- CMake 3.16+
- Eigen 3.x 线性代数库 (已包含)
- 标准库支持

## 应用领域

- 静电场分析
- 热传导问题  
- 结构分析
- 流体压力场
- 偏微分方程数值解

---

*现代化C++有限元求解器，适用于科学计算、工程分析和数值方法研究*

## 开发状态

- **当前分支**: `config类` - 重构配置管理系统
- **主要版本**: v1.0.0
- **构建系统**: CMake 3.16+
- **编程语言**: C++17
