# 二维有限元求解器 (FEM_2D)

这是一个功能完整的二维有限元方法（FEM）求解器，专门用于求解二维泊松方程等椭圆型偏微分方程。采用线性三角形单元，支持复杂几何域和各种边界条件类型。现已支持高密度采样点云VTK输出，便于在ParaView中实现与有限元插值一致的高精度可视化。
## 项目概述

求解器可处理形式为 `-∇·(c(x,y)∇u) = f(x,y)` 的椭圆型偏微分方程，具有以下特点：

- **模块化架构**：清晰的代码组织，便于维护和扩展
- **通用性强**：支持用户自定义问题参数和几何域
- **专业级精度**：完整的误差分析和收敛性验证
- **高性能计算**：针对大规模稀疏矩阵优化
- **工业级网格支持**：集成COMSOL网格导入器
- **高密度采样点云输出**：支持单元内部高密度采样，VTK点云文件便于插值一致性可视化
### 📊 误差分析
- **L∞范数误差**（最大误差，支持单元内部高密度采样）
- **L2范数误差**（积分误差）
- **H1半范数误差**（梯度误差）
- **节点误差字段简化**：VTK输出中的节点误差已简化为节点处误差，采样点云用于高精度误差分布分析
### 🎨 可视化功能
- **VTK文件输出** - 支持ParaView可视化
- **数值解可视化** - 彩色等值线图
- **解析解对比** - 精确解参考图
- **误差分布图** - 绝对和相对误差可视化
- **多数据对比** - 一个文件包含所有分析数据
- **高密度采样点云输出** - 单元内部采样点云（outputVTKDenseSampling），用于插值一致性和高精度误差分布分析
## 项目结构

```
FEM_2D/
├── main_2d.cpp                     # 程序主入口
├── fem_solver_2d.h/.cpp            # 核心求解器
├── problem_definition.h/.cpp       # 问题定义和网格生成
├── comsol_mesh_importer.h/.cpp     # COMSOL网格导入器
├── vtk_output.h/.cpp               # VTK可视化文件输出（含高密度采样点云输出）
├── shape_functions_2d.h/.cpp       # 三角形形函数实现
├── gauss_quadrature_2d.h/.cpp      # 高斯积分点和权重
├── geometry_mapping_2d.h/.cpp      # 几何坐标映射
├── error_analysis_2d.h/.cpp        # 误差分析和计算
├── circle_mesh.mphtxt               # 圆形域网格文件
├── circle_mesh2.mphtxt              # 精细圆形域网格文件
├── Eigen/                           # Eigen线性代数库
├── FEM_2D_solution.exe              # 编译后的可执行文件
├── results/                         # VTK输出文件目录
│   ├── numerical_solution.vtu      # 数值解可视化
│   ├── exact_solution.vtu          # 解析解可视化
│   ├── comparison.vtu              # 对比分析文件
│   └── dense_sampling_points.vtu   # 高密度采样点云（outputVTKDenseSampling）
└── README.md                        # 项目文档
```
## 算法流程

求解器采用标准的有限元方法流程：

### 1. 预处理 (`preprocess2D()`)
- 导入网格或生成结构化网格
- 初始化全局刚度矩阵和载荷向量
- 设置数据结构

### 2. 矩阵组装 (`assemble2D()`)
- 遍历所有三角形单元
- 计算局部刚度矩阵和载荷向量
- 组装全局系统矩阵

### 3. 边界条件施加 (`applyBoundaryConditions2D()`)
- 识别边界节点
- 根据边界条件类型修改系统矩阵
- 优化的批量处理算法

### 4. 线性系统求解 (`solveLinearSystem2D()`)
- 使用Eigen SparseLU求解器
- 稀疏矩阵分解和回代

### 5. 后处理 (`postprocess2D()`)
- 计算误差范数（支持单元内部高密度采样）
- 输出数值解和精确解对比
- 生成分析报告
## ParaView可视化功能 🎨

### VTK文件输出

程序会自动生成以下VTK文件用于可视化分析：

1. **`numerical_solution.vtu`** - 数值解可视化
   - 显示有限元计算结果
   - 包含网格结构和节点解值

2. **`exact_solution.vtu`** - 解析解可视化  
   - 显示理论精确解
   - 用于对比验证

3. **`comparison.vtu`** - 综合对比分析 (🌟推荐)
   - 数值解 (`numerical_solution`)
   - 解析解 (`exact_solution`) 
   - 误差 (`error`)
   - 绝对误差 (`absolute_error`)
   - 相对误差 (`relative_error`)

4. **`dense_sampling_points.vtu`** - 高密度采样点云
   - 单元内部高密度采样点（outputVTKDenseSampling函数生成）
   - 采样点包含数值解、解析解、误差等字段
   - 用于插值一致性分析和高精度误差分布可视化

### ParaView操作步骤

1. **打开文件**: File → Open → 选择`.vtu`文件（可同时加载网格和点云文件）
2. **应用数据**: 点击"Apply"按钮  
3. **选择数据**: 在"Coloring"下拉菜单中选择要显示的字段（如solution、exact_solution、error等）
4. **调整显示**: 
   - 网格文件选择"Surface With Edges"显示结构
   - 点云文件选择"Points"或"Surface"，可用"Glyph"放大采样点
   - 调整色标范围和颜色映射
   - 添加等值线: Filters → Contour
5. **对比分析**: 
   - 打开多个窗口同时显示数值解、解析解和采样点云误差分布
   - 使用误差字段识别高误差区域

### 可视化效果预览

- **🌈 等值线图**: 清晰显示解的空间分布
- **🔷 网格结构**: 三角形单元的几何形状
- **📊 误差分布**: 高误差区域的可视化识别（支持采样点云）
- **📈 3D表面**: 将2D解拉伸为立体显示
- **🟡 高密度采样点云**: 采样点分布与maxerror一致，插值效果与有限元一致
## 版本历史

### 当前版本特性

- ✅ 完整的二维有限元求解器
- ✅ 线性三角形单元实现  
- ✅ 自动三角形网格生成
- ✅ 完整的误差分析模块（L∞、L2、H1 范数，支持单元内部高密度采样）
- ✅ 高性能边界条件处理
- ✅ 稀疏矩阵支持
- ✅ 多种边界条件支持（Dirichlet、Neumann、Robin）
- ✅ 模块化代码结构
- ✅ VS Code 开发环境配置
- ✅ 详细的性能优化文档
- ✅ 高密度采样点云VTK输出，插值一致性可视化
- ✅ 节点误差字段简化，采样点云用于高精度误差分布分析

该求解器为学习和研究二维有限元方法提供了一个完整、高效的参考实现。适用于求解椭圆型偏微分方程，在静电场、热传导、结构分析等领域有广泛应用。

## 依赖库

- **Eigen 3.x**: 线性代数库，用于稀疏矩阵运算
- **C++11**: 支持现代C++特性
- **标准库**: iostream, vector, cmath等

## 编译和运行

### 系统要求

- C++11 或更高版本的编译器 (推荐使用 g++ 或 MSVC)
- Eigen 库（已包含在项目中）
- Windows 系统（当前配置）

### 使用 VS Code 编译和运行

项目配置了 VS Code 任务：

1. **编译项目**: 按 `Ctrl+Shift+P`，选择 "Tasks: Run Task"，然后选择 "Build FEM 2D Project"
2. **运行程序**: 选择 "Run FEM 2D Solution" 任务

### 命令行编译

```bash
g++ -std=c++11 -I . main_2d.cpp fem_solver_2d.cpp problem_definition.cpp \
    shape_functions_2d.cpp gauss_quadrature_2d.cpp geometry_mapping_2d.cpp \
    error_analysis_2d.cpp comsol_mesh_importer.cpp vtk_output.cpp \
    -o FEM_2D_solution.exe
```

### 运行

```bash
./FEM_2D_solution.exe
```

## 输出示例

程序会输出详细的网格信息、计算结果和误差分析：

```
步骤 1: 二维网格生成和预处理完成。
网格生成完成：1681 个节点，3200 个单元
边界边数量：160

步骤 2: 二维全局矩阵组装完成。

步骤 3: 二维边界条件施加完成。
开始应用边界条件...
狄利克雷边界条件节点数量: 160
边界条件应用完成！

步骤 4: 二维线性方程组求解完成。

步骤 5: 二维后处理完成。
网格信息: 1681 个节点，3200 个单元
边界边数量: 160

=== 详细误差分析 ===
--- 误差范数结果 ---
L∞范数（h^2）:  2.456e-03
L2范数（h^2）:        8.234e-04  
H1半范数（h^1）:      3.167e-02

--- VTK文件输出 ---
✓ 数值解VTK文件已生成
✓ 解析解VTK文件已生成
✓ 对比VTK文件已生成

🎯 ParaView可视化指南:
1. numerical_solution.vtu - 查看数值解分布
2. exact_solution.vtu     - 查看解析解分布
3. comparison.vtu         - 对比分析(推荐)
```

## ParaView可视化功能 🎨

### VTK文件输出

程序会自动生成三个VTK文件用于可视化分析：

1. **`numerical_solution.vtu`** - 数值解可视化
   - 显示有限元计算结果
   - 包含网格结构和节点解值

2. **`exact_solution.vtu`** - 解析解可视化  
   - 显示理论精确解
   - 用于对比验证

3. **`comparison.vtu`** - 综合对比分析 (🌟推荐)
   - 数值解 (`numerical_solution`)
   - 解析解 (`exact_solution`) 
   - 误差 (`error`)
   - 绝对误差 (`absolute_error`)
   - 相对误差 (`relative_error`)

### ParaView操作步骤

1. **打开文件**: File → Open → 选择`.vtu`文件
2. **应用数据**: 点击"Apply"按钮  
3. **选择数据**: 在"Coloring"下拉菜单中选择要显示的字段
4. **调整显示**: 
   - 选择"Surface With Edges"显示网格
   - 调整色标范围和颜色映射
   - 添加等值线: Filters → Contour
5. **对比分析**: 
   - 打开多个窗口同时显示数值解和解析解
   - 使用误差字段识别计算精度分布

### 可视化效果预览

- **🌈 等值线图**: 清晰显示解的空间分布
- **🔷 网格结构**: 三角形单元的几何形状
- **📊 误差分布**: 高误差区域的可视化识别
- **📈 3D表面**: 将2D解拉伸为立体显示

## 如何解决新问题

### 步骤 1: 修改问题定义

打开 `problem_definition.cpp` 文件，修改以下函数：

```cpp
// 1. 定义方程系数 c(x,y)
double coefficient_c(double x, double y) {
    return 1.0;  // 修改为你的 c(x,y)
}

// 2. 定义源项 f(x,y)
double source_term_f(double x, double y) {
    const double pi = 3.14159265358979323846;
    return 10.0;  // 修改为你的 f(x,y)
}

// 3. 定义精确解（如果已知）
double exact_solution_u(double x, double y) {
    const double pi = 3.14159265358979323846;
    return 10.0 / 4 * (1 - x * x - y * y);  // 修改为你的精确解
}

// 4. 定义精确解的梯度（如果已知）
double exact_solution_du_dx(double x, double y) {
    {
    return -10.0 / 2 * y;  // 修改为你的精确解y方向导数
}
```

### 步骤 2: 设置网格和边界条件

选择合适的网格生成方式：

**方法1: 使用结构化网格**

在 `defineProblem()` 函数中：

```cpp
void defineProblem() {
    Domain domain = {0.0, 1.0, 0.0, 1.0};  // 计算域 [x_min, x_max, y_min, y_max]
    
    // 设置网格密度
    int N1 = 40;  // x方向单元数
    int N2 = 40;  // y方向单元数
    
    // 生成网格
    generateMesh2D_rectangle(domain, N1, N2, boundary_edges);
    
    // 设置边界条件
    for (auto& edge : boundary_edges) {
        edge.bc = BoundaryCondition::Dirichlet(0.0);  // 设置为你需要的边界条件
        // 或者使用:
        // edge.bc = BoundaryCondition::Neumann(flux_value);
        // edge.bc = BoundaryCondition::Robin(h, g);
    }
}
```

**方法2: 使用COMSOL网格文件**

在 `defineProblem_by_mesh_importer()` 函数中修改网格文件名。

### 步骤 3: 重新编译和运行

重新编译程序并运行即可看到新问题的求解结果和误差分析。

## 边界条件类型

项目支持三种边界条件，后两种还没有实现：

1. **Dirichlet边界条件**: `u = value`
   ```cpp
   BoundaryCondition::Dirichlet(value)
   ```

2. **Neumann边界条件**: `c * ∇u·n = flux`  (n为外法向量)
   ```cpp
   BoundaryCondition::Neumann(flux)
   ```

3. **Robin边界条件**: `c * ∇u·n + h * u = g`
   ```cpp
   BoundaryCondition::Robin(h, g)
   ```

## 应用示例

### 静电场问题

当前默认实现的是一个静电场问题：
- **控制方程**: `-∇²φ = ρ(x,y)/ε₀`
- **物理意义**: φ为静电势，ρ为电荷密度
- **边界条件**: φ=0 (接地条件)
- **计算域**: 圆形域（通过COMSOL网格）

### 热传导问题

可通过修改系数函数求解热传导问题：
- **控制方程**: `-∇·(k∇T) = Q`
- **物理意义**: T为温度，k为热导率，Q为热源
- **边界条件**: 可设置温度或热流边界条件

## 性能优化

项目针对大规模网格进行了重要的性能优化：

### 边界条件处理优化

- **优化前**: 40×40网格的边界条件处理需要几分钟
- **优化后**: 相同规模的处理时间缩短到几秒钟

主要优化措施：
1. **批量处理**: 按边界条件类型分组处理
2. **稀疏矩阵优化访问**: 只访问非零元素，避免全矩阵遍历  
3. **减少矩阵压缩**: 将多次压缩操作合并为一次


## 扩展说明

该框架设计时考虑了扩展性：

- **高阶单元**: 可通过扩展 `shape_functions_2d` 模块支持高阶三角形或四边形单元
- **自适应网格**: 可集成网格细化算法
- **非线性问题**: 可在 `fem_solver_2d` 模块中添加牛顿迭代求解器
- **三维扩展**: 可扩展到四面体单元求解三维问题
- **并行计算**: 可利用 Eigen 的并行支持或 OpenMP 加速

### 编译环境

项目配置为Windows环境下的C++开发：

- **编译器**: g++ (MSYS2/MinGW-w64)
- **路径**: `E:\msys64\ucrt64\bin\g++.exe`
- **标准**: C++11
- **优化**: 调试版本 (`-g` 标志)

### VS Code任务配置

项目包含完整的VS Code任务配置：

1. **Build FEM 2D Project**: 编译整个项目
2. **Run FEM 2D Solution**: 运行求解器（依赖编译任务）
3. **C/C++: g++.exe 生成活动文件**: 编译当前文件

### 网格文件支持

项目支持COMSOL Multiphysics网格导入：
- `circle_mesh.mphtxt`: 圆形域粗网格
- `circle_mesh2.mphtxt`: 圆形域精细网格
- 自动解析.mphtxt格式，提取节点坐标、单元连接和边界信息

## 数学背景

### 控制方程

项目求解二维椭圆型偏微分方程：

```
-∇·(c(x,y)∇u) = f(x,y)  in Ω
```

其中：
- `Ω` 是二维计算域
- `u(x,y)` 是未知函数
- `c(x,y)` 是扩散系数
- `f(x,y)` 是源项

### 弱形式

采用Galerkin加权余量法，弱形式为：

```
∫_Ω c(x,y) ∇u·∇v dΩ = ∫_Ω f(x,y)v dΩ + ∫_Γ_N g_N v dΓ
```

其中`v`为试探函数。

### 有限元离散

- **单元类型**: 线性三角形单元 (P1)
- **形函数**: 线性基函数
- **积分**: 1点、3点、4点或7点高斯积分
- **矩阵**: 稀疏存储格式

## 技术亮点

### 1. 通用边界条件框架

统一的边界条件表示形式：`K*(c*∂u/∂n) + L*u = q`

- K=0, L=1: Dirichlet条件 `u = q`
- K=1, L=0: Neumann条件 `c*∂u/∂n = q`
- K=1, L≠0: Robin条件 `c*∂u/∂n + L*u = q`

### 2. 几何映射系统

完整的参考单元到物理单元的映射：
- 正向映射：参考坐标 → 物理坐标
- 反向映射：物理坐标 → 参考坐标  
- 梯度变换：参考梯度 → 物理梯度
- 雅可比行列式：积分变换

### 3. 多重误差度量

- **L∞误差**: `max|u_h - u_exact|`（点误差）
- **L2误差**: `√∫(u_h - u_exact)² dΩ`（积分误差）
- **H1误差**: `√∫|∇(u_h - u_exact)|² dΩ`（梯度误差）

### 4. 高性能稀疏求解

- Eigen::SparseMatrix存储
- SparseLU直接求解器
- 三元组方式高效组装
- 针对边界条件的批量优化

## 代码质量

- **模块化设计**: 清晰的功能分离
- **完整注释**: 详细的函数和变量说明
- **错误处理**: 边界情况检查和错误报告
- **性能监控**: 详细的计算统计信息
- **可扩展性**: 支持用户自定义扩展

*注: 本文档详细描述了FEM_2D项目的使用方法和技术细节，建议结合代码注释一起阅读以获得最佳理解效果。*

### VS Code 配置

项目包含完整的 VS Code 开发环境配置：

- **编译任务**: `Build FEM 2D Project` 
- **运行任务**: `Run FEM 2D Solution`
- **调试配置**: 支持 gdb 调试
- **智能感知**: 配置了 C++ 路径和包含目录

### 依赖库

- **Eigen**: 高性能线性代数库，用于稀疏矩阵运算
- **C++11 标准库**: 使用现代 C++ 特性

## 版本历史

### 当前版本特性

- ✅ 完整的二维有限元求解器
- ✅ 线性三角形单元实现  
- ✅ 自动三角形网格生成
- ✅ 完整的误差分析模块（L∞、L2、H1 范数）
- ✅ 高性能边界条件处理
- ✅ 稀疏矩阵支持
- ✅ 多种边界条件支持（Dirichlet、Neumann、Robin）
- ✅ 模块化代码结构
- ✅ VS Code 开发环境配置
- ✅ 详细的性能优化文档

该求解器为学习和研究二维有限元方法提供了一个完整、高效的参考实现。适用于求解椭圆型偏微分方程，在静电场、热传导、结构分析等领域有广泛应用。
