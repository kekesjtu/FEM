# 有限元求解器 (FEM)

基于C++和Eigen库的高性能**1D/2D/3D**有限元求解器，专门用于求解椭圆型偏微分方程。支持线性单元（线段、三角形、四面体），复杂几何域、多种边界条件和精确的误差分析。

## 核心功能

- 🧮 **椭圆型PDE求解**：泊松方程 `-∇·(c∇u) = f`
- 🔺 **多维度支持**：1D线段、2D三角形、3D四面体单元
- 📐 **网格支持**：通用COMSOL网格导入器（自动识别维度）
- 🎯 **边界条件**：Dirichlet、Neumann、Robin边界条件
- 📊 **误差分析**：L∞、L2、H1范数计算
- 🎨 **可视化**：VTK输出，支持ParaView（2D/3D）
- ⚙️ **配置管理**：统一的配置类管理计算参数
- 🔧 **CMake构建**：现代化的构建系统
- 🔬 **问题库**：预定义问题（Robin边界验证、电热耦合等）

## 项目结构

```
FEM/
├── 核心模块/
│   ├── src/solvers/
│   │   ├── fem_solver.h/.cpp            # 核心FEM求解器
│   │   └── electrothermal_solver.h/.cpp # 电热耦合求解器
│   ├── src/core/
│   │   ├── config.h/.cpp                # 配置管理类
│   │   └── problem_setup.h/.cpp         # 问题设置基类
│   └── main.cpp                         # 主程序入口
├── 数值计算/
│   ├── src/shape_functions/             # 形函数模块
│   │   ├── shape_functions_1d.cpp       # 1D线性形函数
│   │   ├── shape_functions_2d.cpp       # 2D三角形形函数
│   │   ├── shape_functions_3d.cpp       # 3D四面体形函数
│   │   └── shape_functions_factory.cpp  # 形函数工厂
│   ├── src/gauss_quadrature/            # 高斯积分
│   │   ├── gauss_quadrature_1d.cpp      # 1D高斯积分
│   │   ├── gauss_quadrature_2d.cpp      # 2D高斯积分
│   │   ├── gauss_quadrature_3d.cpp      # 3D高斯积分
│   │   └── gauss_quadrature_factory.cpp # 积分工厂
│   └── src/geometry_mapping/            # 几何映射
│       ├── geometry_mapping_1d.cpp      # 1D线段映射
│       ├── geometry_mapping_2d.cpp      # 2D三角形映射
│       ├── geometry_mapping_3d.cpp      # 3D四面体映射
│       └── geometry_mapping_factory.cpp # 映射工厂
├── 网格与IO/
│   ├── src/core/
│   │   ├── comsol_mesh_importer.h/.cpp  # 通用COMSOL网格导入器
│   │   └── config.h/.cpp                # 配置类（自动从mesh_data加载网格）
│   ├── src/output/vtk_output/           # VTK输出模块
│   │   ├── vtk_output_base.cpp          # VTK输出基类
│   │   ├── vtk_output_2d.cpp            # 2D VTK输出
│   │   └── vtk_output_3d.cpp            # 3D VTK输出
│   └── mesh_data/                       # 网格文件目录（所有网格统一放在此文件夹）
│       ├── circle_mesh.mphtxt           # 2D圆形网格
│       ├── circle_mesh2.mphtxt          # 2D细化网格
│       ├── circle_mesh3.mphtxt          # 2D细化网格
│       ├── cube_tet_mesh.mphtxt         # 3D立方体网格
│       └── sphere_tet_mesh.mphtxt       # 3D球体网格
├── 问题库/
│   └── include/problem_definitions.h    # 预定义问题
│       ├── RobinTestProblem             # Robin边界验证
│       ├── ElectricFieldProblem         # 电场问题
│       ├── ThermalFieldProblem          # 热场问题
│       ├── CubeUniformChargeProblem     # 3D均匀电荷问题
│       └── CustomProblem                # 自定义问题模板
├── 测试/
│   └── tests/
│       ├── test_3d_mesh_import.cpp      # 3D网格导入测试
│       ├── test_3d_integration.cpp      # 3D数值积分验证
│       └── test_simple_3d_poisson.cpp   # 3D泊松方程测试
├── 构建系统/
│   ├── CMakeLists.txt                   # 主CMake配置
│   ├── CMakePresets.json                # CMake预设
│   └── cmake/mingw-toolchain.cmake      # MinGW工具链
├── 依赖库/
│   └── third_party/Eigen/               # 线性代数库
├── 输出结果/
│   └── results/                         # VTK可视化文件
└── 文档/
    └── docs/
        ├── README.md                    # 本文档
        ├── 3D求解器验证报告.md         # 3D算法验证
        ├── 电热耦合实现总结.md         # 电热耦合文档
        ├── 三维扩展架构设计方案.md     # 3D扩展设计
        └── vscode使用教程.md           # VS Code使用指南
```

## 算法流程

1. **预处理** - 网格导入（自动识别1D/2D/3D），数据结构初始化
2. **矩阵组装** - 遍历单元计算局部矩阵，组装全局稀疏系统
3. **边界条件** - 识别边界节点/边/面，施加边界条件  
4. **求解** - 稀疏LU分解或共轭梯度法求解线性系统
5. **后处理** - 误差计算，VTK文件输出（支持2D/3D）

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

### 网格文件配置

所有网格文件统一存放在 `mesh_data/` 文件夹下：

```
mesh_data/
├── circle_mesh.mphtxt        # 2D圆形网格
├── circle_mesh2.mphtxt       # 2D细化网格  
├── circle_mesh3.mphtxt       # 2D超细网格
├── cube_tet_mesh.mphtxt      # 3D立方体四面体网格
└── sphere_tet_mesh.mphtxt    # 3D球体四面体网格
```

**切换网格文件**：在 `include/config.h` 中修改：

```cpp
std::string mesh_filename_ = "cube_tet_mesh.mphtxt";  // 修改为所需的网格文件名
```

**注意**：
- 只需提供文件名，不需要包含 `mesh_data/` 路径前缀
- 程序会自动从 `mesh_data/` 文件夹加载网格
- 支持的格式：COMSOL `.mphtxt` 格式
- 网格维度和单元类型会自动识别

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

**2D可视化**：
- `comparison.vtu` - 完整对比数据（推荐）
- `numerical_solution.vtu` - 数值解
- `exact_solution.vtu` - 精确解  
- `solution.vtu` - 基本解数据

**3D可视化**：
- `numerical_solution.vtu` - 3D数值解（四面体网格）
- 支持体渲染、切片、等值面显示

**电热耦合可视化**：
- `electrothermal_V.vtu` - 电势场分布
- `electrothermal_T.vtu` - 温度场分布
- `electrothermal_Q.vtu` - 焦耳热密度分布

可视化步骤：
1. 打开ParaView，加载相应的 `.vtu` 文件
2. 选择显示字段: `numerical_solution`、`exact_solution`、`error`
3. **3D模式**: 使用 Slice/Clip 查看内部场分布
4. 调整色标和等值线显示
5. 对比分析数值解与精确解的差异

项目还生成PNG格式的结果图像（2D），方便快速查看。

## 技术特性

- **多维度支持**: 完整的1D/2D/3D有限元框架
- **高性能**: 稀疏矩阵存储，优化的矩阵组装算法
- **模块化架构**: 
  - 工厂模式创建形函数、几何映射、高斯积分
  - 清晰的代码结构，便于扩展和维护
- **通用网格导入器**: 
  - 自动识别网格维度和单元类型
  - 智能分类体单元和边界单元
- **工业级**: 支持COMSOL网格导入，复杂几何处理
- **精确验证**: 
  - 完整误差分析（L∞/L2/H1范数）
  - 线性响应验证（通过测试）
  - 数值积分正确性验证
- **专业可视化**: 
  - VTK输出，支持ParaView
  - 2D/3D数据可视化
  - 多场耦合结果展示
- **现代化构建**: CMake构建系统，VS Code完整集成
- **配置灵活**: 统一的Config类，问题库系统

## 主要模块

### 核心求解器 (`fem_solver`)
- 全局刚度矩阵组装（支持1D/2D/3D）
- 多种边界条件处理（Dirichlet/Neumann/Robin）
- 线性系统求解（稀疏LU、共轭梯度法）
- 预处理和后处理流程

### 电热耦合求解器 (`electrothermal_solver`)
- 电场-热场双向耦合
- 温度依赖电导率 σ(T) = σ₀[1 + α(T-T₀)]
- 焦耳热生成 Q = σ(T)|∇V|²
- 迭代求解收敛算法

### 配置管理 (`config`)
- 自动网格维度识别
- 问题参数配置
- 边界条件定义（Dirichlet、Neumann、Robin）
- 高斯积分点自动设置
- 几何和数值参数管理

### 问题库 (`problem_definitions`)
- **RobinTestProblem**: Robin边界条件验证
- **ElectricFieldProblem**: 2D电场问题
- **ThermalFieldProblem**: 2D热传导问题
- **CubeUniformChargeProblem**: 3D均匀电荷分布
- **CustomProblem**: 自定义问题模板

### 形函数模块 (`shape_functions`)
- **1D**: 线性线段单元
- **2D**: 线性三角形单元
- **3D**: 线性四面体单元
- 工厂模式统一创建接口

### 几何映射 (`geometry_mapping`)
- 参考坐标↔物理坐标映射
- 雅可比矩阵计算
- 梯度坐标变换
- 支持1D线段、2D三角形、3D四面体

### 高斯积分 (`gauss_quadrature`)
- 1D: 1-5点高斯积分
- 2D三角形: 1点、3点、7点积分
- 3D四面体: 1点、4点积分
- 自动选择积分点数量

### 误差分析 (`error_analysis`)
- L∞、L2、H1范数计算
- 支持2D和3D误差分析
- 收敛性分析
- 精确解比较

### 网格处理 (`comsol_mesh_importer`)
- 通用COMSOL `.mphtxt` 格式导入
- 自动识别1D/2D/3D网格
- 智能单元分类（体单元/边界单元）
- 节点和单元数据解析
- 边界自动识别

### 可视化输出 (`vtk_output`)
- VTK UnstructuredGrid格式
- 2D三角形网格输出
- 3D四面体网格输出
- 数值解、精确解、误差分布
- 工厂模式支持多维度输出

## 依赖要求

- C++17 编译器 (g++/MSVC)
- CMake 3.16+
- Eigen 3.x 线性代数库 (已包含)
- 标准库支持

## 应用领域

- **静电场分析**: 电势分布、电荷密度计算（支持3D）
- **热传导问题**: 稳态热传导、热源分布
- **电热耦合**: 温度依赖电导率的耦合问题
- **结构分析**: 位移场计算
- **流体压力场**: 不可压缩流动压力求解
- **偏微分方程数值解**: 通用椭圆型PDE求解器

## 3D求解器验证

3D四面体有限元求解器已通过完整验证：

### ✅ 数值积分测试
- 体积计算精度: < 1e-10
- 高斯积分权重: 精确
- 形函数单位分解: 满足
- 梯度变换: 正确
- 刚度矩阵: 对称且行和为零

### ✅ 线性响应测试
- 源项 f=1: u_max ≈ 0.056
- 源项 f=10: u_max ≈ 0.56
- 比例关系: **10.0:1.0** (精确)
- **结论**: 求解器满足线性叠加原理

### ✅ 解的数量级验证
对于单位立方体 [0,1]³，边界u=0，源项f=1:
- 理论估计: u_center ~ 0.01-0.1
- 数值结果: u_max = 0.056
- **结论**: 解的数量级合理

详细验证报告见 `docs/3D求解器验证报告.md`

---

*现代化C++ 1D/2D/3D有限元求解器，适用于科学计算、工程分析和数值方法研究*

## 开发状态

- **当前分支**: `config类` - 重构配置管理系统，实现3D扩展
- **主要版本**: v2.0.0 (3D支持)
- **构建系统**: CMake 3.16+
- **编程语言**: C++17
- **核心特性**: 
  - ✅ 1D/2D/3D完整支持
  - ✅ 通用网格导入器
  - ✅ 模块化架构（工厂模式）
  - ✅ 电热耦合求解
  - ✅ 问题库系统
  - ✅ 完整验证测试

## 贡献与反馈

欢迎提交Issue和Pull Request！

## 许可证

本项目遵循MIT许可证。
