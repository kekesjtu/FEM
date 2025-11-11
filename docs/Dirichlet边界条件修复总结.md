# Dirichlet边界条件修复总结

## 问题描述

在电热耦合求解器中发现一个严重的bug:**Dirichlet边界条件未被精确满足**。

### 症状
- 边界条件设置:`V = 10.0V`(左半圆),`V = 0.0V`(右半圆)
- 求解结果:电势范围 `[0, 9.6781]V` 而不是预期的 `[0, 10]V`
- 即使设置了60个边界节点为10V,求解后最大值仍然小于10V

## 根本原因分析

### 问题1:三元组重复累加导致对角元素错误

在应用Dirichlet边界条件时,代码遍历所有三元组,对每个Dirichlet节点的对角元素添加 `T_entry(row, col, 1.0)`:

```cpp
// 原始错误代码
for (const auto& triplet : triplet_list_)
{
    int row = triplet.row();
    int col = triplet.col();
    
    if (row == col && is_dirichlet_node[row])
    {
        filtered_triplets.push_back(T_entry(row, col, 1.0));  // ❌ 重复添加!
    }
}
```

**问题**:一个节点的对角元素 `K(i,i)` 可能出现在多个三元组中(来自不同单元的组装)。如果一个边界节点被3个单元共享,则对角元素会被设置3次,`setFromTriplets()` 会累加这些值,导致 `K(i,i) = 3` 而不是 `1.0`。

**验证输出**:
```
[调试] Dirichlet节点0: b=10, K对角=3, 该行非零元素数=1
```

**后果**:对于Dirichlet节点 i,方程变成 `3.0 * u_i = 10.0`,求解得 `u_i = 3.33V`,而不是预期的 `10V`。

### 问题2:矩阵对称性处理

在应用Dirichlet边界条件时,需要正确处理矩阵的对称性:

- **Dirichlet节点 i 的行**: 只保留 `K(i,i) = 1.0`,删除 `K(i,j)` (j≠i)
- **Dirichlet节点 i 的列**: 删除 `K(j,i)` (j≠i),并修改右端项 `b(j) -= K(j,i) * u_i`

这样可以保证:
1. 矩阵保持对称性(CG求解器要求)
2. Dirichlet条件被精确满足:`1.0 * u_i = u_i`

## 解决方案

### 修复代码

```cpp
// 第二遍：处理 Dirichlet 边界条件，从三元组中过滤
if (!dirichlet_nodes.empty())
{
    std::cout << "  应用 Dirichlet 边界条件，过滤三元组..." << std::endl;

    // 第一步:先修改右端项(需要用到原始矩阵元素)
    for (const auto& triplet : triplet_list_)
    {
        int row = triplet.row();
        int col = triplet.col();
        double value = triplet.value();

        // 如果列是Dirichlet节点但行不是,修改右端项
        if (is_dirichlet_node[col] && !is_dirichlet_node[row])
        {
            b_(row) -= value * dirichlet_values[col];
        }
    }

    // 第二步:过滤三元组,只保留内部节点间的耦合和Dirichlet节点的对角元素
    std::vector<T_entry> filtered_triplets;
    filtered_triplets.reserve(triplet_list_.size());
    
    // ✅ 关键修复:标记已经添加过对角元素的Dirichlet节点
    std::vector<bool> dirichlet_diag_added(N_, false);

    for (const auto& triplet : triplet_list_)
    {
        int row = triplet.row();
        int col = triplet.col();

        // 如果是Dirichlet节点的对角元素,只添加一次并设为1
        if (row == col && is_dirichlet_node[row])
        {
            if (!dirichlet_diag_added[row])
            {
                filtered_triplets.push_back(T_entry(row, col, 1.0));
                dirichlet_diag_added[row] = true;  // ✅ 标记已添加
            }
            // 否则跳过,避免重复添加
        }
        // 如果行和列都不是Dirichlet节点,保留原值
        else if (!is_dirichlet_node[row] && !is_dirichlet_node[col])
        {
            filtered_triplets.push_back(triplet);
        }
        // 其他情况(行或列是Dirichlet节点但不是对角元素):删除,保持矩阵对称
    }

    // 替换三元组列表
    triplet_list_ = std::move(filtered_triplets);

    // 设置Dirichlet节点的右端项
    for (int i : dirichlet_nodes)
    {
        b_(i) = dirichlet_values[i];
    }
}
```

### 关键改进点

1. **使用标记数组 `dirichlet_diag_added`**: 确保每个Dirichlet节点的对角元素只添加一次
2. **两步处理**: 先修改右端项(需要原始矩阵元素),再过滤三元组
3. **保持对称性**: 对称地删除Dirichlet节点的行和列(除对角元素外)

## 验证结果

修复后的输出:
```
Dirichlet 值范围: [0, 10]          ✓ 边界条件设置正确
[调试] Dirichlet节点0: b=10, K对角=1, 该行非零元素数=1  ✓ 对角元素正确
解向量u_范围: [0, 10]              ✓ 求解结果正确
电势范围: [0, 10] V                ✓ Dirichlet条件被精确满足
```

## 影响范围

这个bug影响**所有使用Dirichlet边界条件的问题**:
- ✅ 电热耦合问题(电场和热场)
- ✅ Robin边界条件验证算例
- ✅ 3D泊松方程问题
- ✅ 所有标准FEM问题

## 经验教训

1. **稀疏矩阵三元组累加**: `setFromTriplets()` 会自动累加相同位置的元素,必须避免重复添加
2. **边界条件应用顺序**: 必须先修改右端项(基于原始矩阵),再修改矩阵结构
3. **CG求解器要求**: 共轭梯度法要求矩阵对称正定,边界条件应用必须保持对称性
4. **调试技巧**: 检查关键节点的矩阵行结构和右端项,可以快速定位问题

## 相关文件

- `src/solvers/fem_solver.cpp`: 主要修复位置(`applyBoundaryConditions()`)
- `include/problem_definitions.h`: 边界条件设置(使用几何实体)
- `src/solvers/electrothermal_solver.cpp`: 电热耦合求解器

---

**修复日期**: 2025年11月12日  
**影响版本**: v1.0 - v1.2  
**修复版本**: v1.3+
