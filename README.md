# 通用一维有限元求解器

这是一个用于求解一般形式为 `-d/dx(c(x) * du/dx) = f(x)` 的二阶常微分方程的通用一维有限元方法（FEM）求解器。

该求解器具有模块化设计，支持自动误差分析，并提供了统一的几何映射接口。

## 项目特色

- **模块化设计**：清晰的代码组织，便于维护和扩展
- **通用性强**：只需修改问题定义模块即可解决新的物理问题
- **完整的误差分析**：自动计算 L∞、L2 和 H1 范数误差
- **统一几何映射**：简洁的参考单元到物理单元映射接口
- **稀疏矩阵支持**：使用 Eigen 稀疏矩阵提高计算效率
- **灵活的边界条件**：支持 Dirichlet、Neumann 和 Robin 边界条件

## 项目结构

```
E:/code/
├── Eigen/                   # Eigen库 (线性代数)
├── main.cpp                 # 程序主入口
├── fem_solver.h/.cpp        # 通用有限元求解器核心
├── problem_definition.h/.cpp # 【用户定义】问题描述
├── shape_functions.h/.cpp   # 形函数库
├── gauss_quadrature.h/.cpp  # 高斯积分
├── geometry_mapping.h/.cpp  # 几何映射模块
├── error_analysis.h/.cpp    # 误差分析模块
└── README.md                # 本文档
```

### 模块详细说明

#### 核心求解器模块

- **`main.cpp`**: 程序入口，按顺序调用 FEM 分析的各个步骤
- **`fem_solver` (.h/.cpp)**: 求解器核心，包含全局变量声明和 FEM 流程实现
- **`gauss_quadrature` (.h/.cpp)**: 提供高斯-勒让德积分点和权重
- **`shape_functions` (.h/.cpp)**: 一维线性形函数及其导数的实现

#### 新增功能模块

- **`geometry_mapping` (.h/.cpp)**: 
  - 统一的几何映射接口
  - 支持参考单元 [-1,1] 到物理单元的坐标变换
  - 提供雅可比变换和导数变换功能

- **`error_analysis` (.h/.cpp)**:
  - 完整的误差分析功能
  - 支持 L∞ 范数（最大误差）计算
  - 支持 L2 范数误差计算（使用高斯积分）
  - 支持 H1 范数误差计算（导数误差）

#### 用户定义模块

- **`problem_definition` (.h/.cpp)**:
  - **这是解决新问题时主要需要修改的模块**
  - 定义方程系数 `coefficient_c(x)`
  - 定义源项 `source_term_f(x)`
  - 定义精确解 `exact_solution_u(x)` 和导数 `exact_solution_du_dx(x)`
  - 设置网格参数和边界条件

## 编译和运行

### 系统要求

- C++11 或更高版本的编译器
- Eigen 库（已包含在项目中）

### 使用 VS Code 编译和运行

项目配置了 VS Code 任务：

1. **编译项目**: 按 `Ctrl+Shift+P`，选择 "Tasks: Run Task"，然后选择 "Build FEM Project"
2. **运行程序**: 选择 "Run FEM Solution" 任务

### 命令行编译

```bash
g++ -std=c++11 -I . main.cpp fem_solver.cpp problem_definition.cpp \
    shape_functions.cpp gauss_quadrature.cpp geometry_mapping.cpp \
    error_analysis.cpp -o FEM_solution.exe
```

### 运行

```bash
./FEM_solution.exe
```

## 输出示例

程序会输出详细的计算结果和误差分析：

```
--- 计算结果 ---
节点ID  坐标 (x)     FEM解 (u)          精确解           误差
-------------------------------------------------------------------------------------
0         0.0000         -1.94583e-02        0.00000e+00         1.94583e-02
1         0.3333         2.99454e-01         3.14986e-01         1.55313e-02
2         0.6667         5.08010e-01         5.23925e-01         1.59146e-02
3         1.0000         5.20929e-01         5.40302e-01         1.93733e-02

--- 误差分析 ---
误差范数:
  最大误差 (L∞范数): 4.587392e-02
  L2范数误差:        3.024610e-02
  H1范数误差:        1.401112e-01

网格信息:
  单元数量: 3
  节点数量: 4
  平均单元长度: 0.333333
```

## 如何解决新问题

### 步骤 1: 修改问题定义

打开 `problem_definition.cpp` 文件，修改以下函数：

```cpp
// 1. 定义方程系数 c(x)
double coefficient_c(double x) {
    return exp(x);  // 修改为你的 c(x)
}

// 2. 定义源项 f(x)
double source_term_f(double x) {
    return -exp(x) * (cos(x) - 2 * sin(x) - x * cos(x) - x * sin(x));  // 修改为你的 f(x)
}

// 3. 定义精确解（如果已知）
double exact_solution_u(double x) {
    return x * cos(x);  // 修改为你的精确解
}

// 4. 定义精确解的导数（如果已知）
double exact_solution_du_dx(double x) {
    return cos(x) - x * sin(x);  // 修改为你的精确解导数
}
```

### 步骤 2: 设置网格和边界条件

在 `defineProblem()` 函数中：

```cpp
void defineProblem(int& N_out, int& M_out, std::vector<BoundaryCondition>& bcs_out) {
    // 设置网格
    N_out = 11;  // 节点数
    M_out = N_out - 1;  // 单元数
    
    // 设置边界条件
    bcs_out.clear();
    bcs_out.push_back(BoundaryCondition::Dirichlet(0, 0.0));      // 左端点: u=0
    bcs_out.push_back(BoundaryCondition::Dirichlet(N_out-1, 1.0)); // 右端点: u=1
    // 或者设置 Neumann/Robin 边界条件
    // bcs_out.push_back(BoundaryCondition::Neumann(0, flux_value));
    // bcs_out.push_back(BoundaryCondition::Robin(N_out-1, h, g));
}
```

### 步骤 3: 重新编译和运行

重新编译程序并运行即可看到新问题的求解结果和误差分析。

## 边界条件类型

项目支持三种边界条件：

1. **Dirichlet 边界条件**: `u = value`
   ```cpp
   BoundaryCondition::Dirichlet(node_index, value)
   ```

2. **Neumann 边界条件**: `c * du/dx = flux`
   ```cpp
   BoundaryCondition::Neumann(node_index, flux)
   ```

3. **Robin 边界条件**: `c * du/dx + h * u = g`
   ```cpp
   BoundaryCondition::Robin(node_index, h, g)
   ```

## 扩展说明

该框架设计时考虑了扩展性：

- **高阶单元**: 可通过扩展 `shape_functions` 模块支持高阶单元
- **多维问题**: 可扩展 `geometry_mapping` 模块支持二维/三维问题  
- **非线性问题**: 可在 `fem_solver` 模块中添加迭代求解器
- **其他物理场**: 通过修改 `problem_definition` 模块可求解不同的物理问题

## 版本历史

### 最新版本特性

- ✅ 完整的误差分析模块（L∞、L2、H1 范数）
- ✅ 统一的几何映射接口
- ✅ 稀疏矩阵支持
- ✅ 多种边界条件支持
- ✅ 模块化代码结构
- ✅ VS Code 开发环境配置

该求解器为学习和研究有限元方法提供了一个清晰、完整的参考实现。
