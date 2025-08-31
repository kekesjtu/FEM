# 通用一维有限元求解器

这是一个用于求解一般形式为 `-d/dx(c(x) * du/dx) = f(x)` 的二阶常微分方程的通用一维有限元方法（FEM）求解器。

此版本采用纯函数式和全局变量的结构，以匹配原始代码的风格。

## 项目结构

代码库被清晰地划分为几个独立的模块，以实现高度的通用性和可维护性。核心思想是：**用户只需修改 `problem_definition` 模块即可解决新的物理问题，而无需触及求解器的核心代码。**

```
E:/code/
├── Eigen/              # Eigen库 (线性代数)
├── main.cpp            # 程序主入口
├── fem_solver.h        # 通用有限元求解器 (头文件, 全局变量声明)
├── fem_solver.cpp      # 通用有限元求解器 (实现)
├── problem_definition.h # 【用户定义】问题描述 (头文件)
├── problem_definition.cpp # 【用户定义】问题描述 (实现)
├── shape_functions.h   # 形函数库 (头文件)
├── shape_functions.cpp # 形函数库 (实现)
├── gauss_quadrature.h  # 高斯积分 (头文件)
├── gauss_quadrature.cpp# 高斯积分 (实现)
└── README.md           # 本文档
```

### 模块说明

- **`main.cpp`**: 
  - 程序的入口点。
  - 按顺序调用有限元分析的各个步骤函数 (`preprocess`, `assemble`, `solve`, `postprocess`)。

- **`fem_solver` (.h/.cpp)**:
  - 这是求解器的核心，`.h` 文件中声明了所有全局变量和核心FEM函数。
  - `.cpp` 文件中定义了全局变量并实现了FEM流程函数。
  - **通常情况下，你不需要修改此模块。**

- **`problem_definition` (.h/.cpp)**:
  - **这是解决新问题时唯一需要修改的模块。**
  - 你可以在这里定义：
    - `coefficient_c(x)`: 方程中的物理系数。
    - `source_term_f(x)`: 方程中的源项。
    - `exact_solution_u(x)`: 问题的精确解（用于误差分析）。
    - `defineProblem()`: 网格参数（节点数）和边界条件。

- **`shape_functions` (.h/.cpp)**:
  - 定义了有限元方法中使用的基函数（形函数）及其导数。
  - 当前实现为一维线性形函数。

- **`gauss_quadrature` (.h/.cpp)**:
  - 提供了一维高斯-勒让德积分所需的积分点和权重，用于在计算刚度矩阵和载荷向量时进行数值积分。

## 如何编译和运行

你需要一个支持 C++11 或更高版本的编译器（如 G++）以及 Eigen 库。

### 编译命令示例 (使用 G++)

```bash
# 确保你的编译器可以找到 Eigen 库的头文件
# -I E:/code/Eigen 是一个例子，请根据你的实际路径修改
g++ -std=c++11 -I E:/code/Eigen main.cpp fem_solver.cpp problem_definition.cpp shape_functions.cpp gauss_quadrature.cpp -o FEM_solution_new.exe
```

### 运行

```bash
./FEM_solution_new.exe
```

## 如何解决一个新问题

1.  **打开 `problem_definition.cpp` 文件。**
2.  **修改 `coefficient_c(x)` 函数** 以匹配你新问题中的 `c(x)`。
3.  **修改 `source_term_f(x)` 函数** 以匹配你新问题中的 `f(x)`。
4.  **(可选) 修改 `exact_solution_u(x)` 函数**，如果你的问题有已知的精确解，填入此处可以方便地进行误差评估。
5.  **修改 `defineProblem()` 函数**：
    - 设置 `N_out`（总节点数）。
    - 在 `bcs_out` 向量中定义你的狄利克雷边界条件（指定节点索引和值）。
6.  **重新编译并运行程序。**
