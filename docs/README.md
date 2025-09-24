# 二维有限元求解器 (FEM_2D)

基于C++和Eigen库的高性能二维有限元求解器，专门用于求解椭圆型偏微分方程。采用线性三角形单元，支持复杂几何域、多种边界条件和精确的误差分析。

## 核心功能

- 🧮 **椭圆型PDE求解**：泊松方程 `-∇·(c∇u) = f`
- 🔺 **线性三角形单元**：标准有限元离散化
- 📐 **网格支持**：内置生成器 + COMSOL导入器
- 🎯 **边界条件**：Dirichlet、Neumann、Robin
- 📊 **误差分析**：L∞、L2、H1范数计算
- 🎨 **可视化**：VTK输出，支持ParaView

## 项目结构

```
FEM_2D/
├── 核心模块/
│   ├── main_2d.cpp              # 主程序入口
│   ├── fem_solver_2d.h/.cpp     # 核心FEM求解器
│   ├── problem_definition.h/.cpp # 问题定义(PDE参数、边界条件)
│   └── error_analysis_2d.h/.cpp  # 误差分析模块
├── 数值计算/
│   ├── shape_functions_2d.h/.cpp     # 三角形形函数
│   ├── gauss_quadrature_2d.h/.cpp    # 高斯积分
│   └── geometry_mapping_2d.h/.cpp    # 几何坐标映射
├── 网格与IO/
│   ├── comsol_mesh_importer.h/.cpp   # COMSOL网格导入
│   ├── vtk_output.h/.cpp             # VTK可视化输出
│   ├── circle_mesh.mphtxt            # 网格数据文件
│   └── circle_mesh2.mphtxt
├── 依赖库/
│   └── Eigen/                        # 线性代数库
├── 输出结果/
│   └── results/                      # VTK可视化文件
└── 文档/
    └── md/                           # 技术文档
```

## 算法流程

1. **预处理** - 网格导入/生成，数据结构初始化
2. **矩阵组装** - 遍历单元计算局部矩阵，组装全局系统
3. **边界条件** - 识别边界节点，修改系统矩阵  
4. **求解** - 稀疏LU分解求解线性系统
5. **后处理** - 误差计算，VTK文件输出

## 快速开始

### 编译运行

**VS Code环境**（推荐）:
- 按 `Ctrl+Shift+P` → "Tasks: Run Task" → "Build FEM 2D Project"
- 运行: "Run FEM 2D Solution"

**命令行编译**:
```bash
g++ -std=c++17 -I . main_2d.cpp fem_solver_2d.cpp problem_definition.cpp \
    shape_functions_2d.cpp gauss_quadrature_2d.cpp geometry_mapping_2d.cpp \
    error_analysis_2d.cpp comsol_mesh_importer.cpp vtk_output.cpp \
    mesh_hierarchy.cpp -o FEM_2D_solution.exe
```

### 输出示例

```
步骤 1: 二维网格生成和预处理完成。
网格生成完成：1681 个节点，3200 个单元

步骤 2-4: 矩阵组装、边界条件、求解完成。

步骤 5: 二维后处理完成。
=== 详细误差分析 ===
L∞范数（h^2）:  2.456e-03
L2范数（h^2）:   8.234e-04  
H1半范数（h^1）: 3.167e-02

✓ VTK文件已生成 (results/)
```

## ParaView可视化

1. 打开 `results/comparison.vtu` (推荐 - 包含所有数据)
2. 选择显示字段: `numerical_solution`、`exact_solution`、`error`
3. 调整色标和等值线显示
4. 对比分析数值解与精确解

## 技术特性

- **高性能**: 稀疏矩阵存储，优化的矩阵组装算法
- **模块化**: 清晰的代码结构，便于扩展和维护  
- **工业级**: 支持COMSOL网格导入，复杂几何处理
- **精确验证**: 完整误差分析，收敛性验证
- **可视化**: VTK输出，支持ParaView专业可视化

## 依赖要求

- C++17 编译器 (g++/MSVC)
- Eigen 3.x 线性代数库 (已包含)
- 标准库支持

## 应用领域

- 静电场分析
- 热传导问题  
- 结构分析
- 流体压力场
- 偏微分方程数值解

---

*二维有限元方法的完整C++实现，适用于学习研究和工程应用*
