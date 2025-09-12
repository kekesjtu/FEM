# 二维有限元求解器 (FEM_2D)

这是一个用于求解二维泊松方程的通用有限元方法（FEM）求解器，采用线性三角形单元。求解器可处理形式为 `-∇·(c(x,y)∇u) = f(x,y)` 的椭圆型偏微分方程。

该求解器具有模块化设计，支持完整的误差分析，并针对大规模网格进行了性能优化。

## 项目特色

- **二维三角形网格**：自动生成结构化三角形网格
- **模块化设计**：清晰的代码组织，便于维护和扩展  
- **通用性强**：只需修改问题定义模块即可解决不同的二维物理问题
- **完整的误差分析**：自动计算 L∞、L2 和 H1 范数误差
- **高性能优化**：针对大规模网格优化了边界条件处理
- **稀疏矩阵支持**：使用 Eigen 稀疏矩阵提高计算效率
- **灵活的边界条件**：支持 Dirichlet、Neumann 和 Robin 边界条件

## 项目结构

```
FEM_2D/
├── Eigen/                          # Eigen库 (线性代数)
├── main_2d.cpp                     # 程序主入口
├── fem_solver_2d.h/.cpp            # 二维有限元求解器核心
├── problem_definition.h/.cpp       # 【用户定义】二维问题描述
├── shape_functions_2d.h/.cpp       # 二维三角形形函数库
├── gauss_quadrature_2d.h/.cpp      # 二维高斯积分
├── geometry_mapping_2d.h/.cpp      # 二维几何映射模块
├── error_analysis_2d.h/.cpp        # 二维误差分析模块
├── boundary_conditions_optimization.md # 边界条件性能优化文档
├── FEM_2D_solution.exe             # 编译后的可执行文件
└── README.md                       # 本文档
```

### 模块详细说明

#### 核心求解器模块

- **`main_2d.cpp`**: 程序入口，按顺序调用二维 FEM 分析的各个步骤
- **`fem_solver_2d` (.h/.cpp)**: 求解器核心，包含全局变量声明和二维 FEM 流程实现
- **`gauss_quadrature_2d` (.h/.cpp)**: 提供二维三角形单元的高斯积分点和权重
- **`shape_functions_2d` (.h/.cpp)**: 二维线性三角形形函数及其梯度的实现

#### 功能模块

- **`geometry_mapping_2d` (.h/.cpp)**: 
  - 二维几何映射接口
  - 支持参考单元到物理三角形单元的坐标变换
  - 提供雅可比矩阵和梯度变换功能

- **`error_analysis_2d` (.h/.cpp)**:
  - 完整的二维误差分析功能
  - 支持 L∞ 范数（最大误差）计算
  - 支持 L2 范数误差计算（使用高斯积分）
  - 支持 H1 半范数误差计算（梯度误差）

#### 用户定义模块

- **`problem_definition` (.h/.cpp)**:
  - **这是解决新问题时主要需要修改的模块**
  - 定义方程系数 `coefficient_c(x,y)`
  - 定义源项 `source_term_f(x,y)`
  - 定义精确解 `exact_solution_u(x,y)` 和梯度 `exact_solution_du_dx(x,y)`, `exact_solution_du_dy(x,y)`
  - 设置网格参数和边界条件
  - 自动生成结构化三角形网格

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
    error_analysis_2d.cpp -o FEM_2D_solution.exe
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
处理 160 个边界节点的狄利克雷边界条件
边界条件应用完成

步骤 4: 二维线性方程组求解完成。
线性系统求解完成

步骤 5: 二维后处理完成。

--- 二维误差分析结果 ---
误差范数：
  最大误差 (L∞范数): 2.456e-03
  L2范数误差:        8.234e-04  
  H1半范数误差:      3.167e-02

网格信息：
  单元数量: 3200
  节点数量: 1681
  网格规模: 40×40 (x方向×y方向)
```
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
    return 2.0 * pi * pi * sin(pi * x) * sin(pi * y);  // 修改为你的 f(x,y)
}

// 3. 定义精确解（如果已知）
double exact_solution_u(double x, double y) {
    const double pi = 3.14159265358979323846;
    return sin(pi * x) * sin(pi * y);  // 修改为你的精确解
}

// 4. 定义精确解的梯度（如果已知）
double exact_solution_du_dx(double x, double y) {
    const double pi = 3.14159265358979323846;
    return pi * cos(pi * x) * sin(pi * y);  // 修改为你的精确解x方向导数
}

double exact_solution_du_dy(double x, double y) {
    const double pi = 3.14159265358979323846;
    return pi * sin(pi * x) * cos(pi * y);  // 修改为你的精确解y方向导数
}
```

### 步骤 2: 设置网格和边界条件

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

### 步骤 3: 重新编译和运行

重新编译程序并运行即可看到新问题的求解结果和误差分析。

## 边界条件类型

项目支持三种边界条件：

1. **Dirichlet 边界条件**: `u = value`
   ```cpp
   BoundaryCondition::Dirichlet(value)
   ```

2. **Neumann 边界条件**: `c * ∇u·n = flux`  (n为外法向量)
   ```cpp
   BoundaryCondition::Neumann(flux)
   ```

3. **Robin 边界条件**: `c * ∇u·n + h * u = g`
   ```cpp
   BoundaryCondition::Robin(h, g)
   ```

## 应用示例

### 静电场问题

当前默认实现的是一个静电场问题：
- **控制方程**: `-∇²φ = ρ(x,y)/ε₀`
- **物理意义**: φ为静电势，ρ为电荷密度
- **边界条件**: φ=0 (接地条件)
- **计算域**: 单位正方形 [0,1]×[0,1]

### 热传导问题

可通过修改系数函数求解热传导问题：
- **控制方程**: `-∇·(k∇T) = Q`
- **物理意义**: T为温度，k为热导率，Q为热源
- **边界条件**: 可设置温度或热流边界条件

## 性能优化

### 边界条件处理优化

项目针对大规模网格进行了重要的性能优化，详见 `boundary_conditions_optimization.md`：

- **优化前**: 40×40网格的边界条件处理需要几分钟
- **优化后**: 相同规模的处理时间缩短到几秒钟

主要优化措施：
1. **批量处理**: 按边界条件类型分组处理
2. **稀疏矩阵优化访问**: 只访问非零元素，避免全矩阵遍历  
3. **减少矩阵压缩**: 将多次压缩操作合并为一次

### 网格规模建议

- **小规模测试**: 10×10 网格 (121个节点)
- **中等规模**: 20×20 网格 (441个节点)  
- **大规模计算**: 40×40 网格 (1681个节点)
- **超大规模**: 80×80 网格 (6561个节点，需要足够内存)

## 扩展说明

该框架设计时考虑了扩展性：

- **高阶单元**: 可通过扩展 `shape_functions_2d` 模块支持高阶三角形或四边形单元
- **自适应网格**: 可集成网格细化算法
- **非线性问题**: 可在 `fem_solver_2d` 模块中添加牛顿迭代求解器
- **三维扩展**: 可扩展到四面体单元求解三维问题
- **并行计算**: 可利用 Eigen 的并行支持或 OpenMP 加速

## 开发环境

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
