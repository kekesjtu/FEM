# muParser 指针生命周期问题修复

## 问题描述

在使用 muParser 库评估材料属性表达式时，发现**电导率计算结果为负数**（`σ = -3.919e+08`），导致焦耳热也为负值。这是一个严重的物理错误，因为电导率和焦耳热在物理上必须为正。

### 症状表现

```
[调试详细] ΔT = 0, σ = -3.919e+08, |∇V|² = 21.2487
[调试] 单元 0 积分点 0: Q = -8.32738e+09 W/m³
```

- 温度增量 `ΔT = 0` (正确)
- 电势梯度平方 `|∇V|² = 21.2487` (正数，正确)
- **电导率 `σ = -3.919e+08`** ❌ (负数，错误！)
- 焦耳热 `Q = σ × |∇V|² = -8.32738e+09` ❌ (负数，错误！)

### 材料配置

```json
"electrical_conductivity": {
  "type": "expression",
  "formula": "sigma0 / (1 + alpha * (T - T0))",
  "variables": ["T"],
  "parameters": {
    "sigma0": 5.96e7,
    "alpha": 0.00393,
    "T0": 293.15
  }
}
```

期望结果：当 `T = 300K` 时，
```
σ = 5.96e7 / (1 + 0.00393 × (300 - 293.15))
  = 5.96e7 / (1 + 0.02689)
  = 5.96e7 / 1.02689
  ≈ 5.80e7 S/m  (正数！)
```

## 根本原因分析

### 问题根源：悬空指针（Dangling Pointer）

在 `MaterialProperty` 类中，muParser 需要变量的**指针**来动态更新值。原始实现存在严重的指针生命周期问题：

#### 原始错误实现

**头文件 (`material.h`)**:
```cpp
class MaterialProperty {
private:
    mutable std::shared_ptr<mu::Parser> parser_;
    mutable bool parser_initialized_;
    // ❌ 缺少存储变量值的成员
};
```

**实现文件 (`material.cpp`)** - 错误的初始化逻辑:
```cpp
void MaterialProperty::initializeParser(const EvaluationContext& ctx) const
{
    parser_ = std::make_shared<mu::Parser>();
    
    // 获取临时对象 ctx 的变量指针映射
    auto var_map = const_cast<EvaluationContext&>(ctx).getVariableMap();
    
    // ❌ 关键错误：使用临时对象的指针
    for (const auto& var_name : variables_) {
        if (var_map.find(var_name) != var_map.end()) {
            parser_->DefineVar(var_name.c_str(), var_map[var_name]);
            // var_map[var_name] 指向 ctx.T, ctx.V 等成员
            // 但 ctx 是临时对象，函数返回后就被销毁！
        }
    }
}

double MaterialProperty::evaluate(const EvaluationContext& ctx) const
{
    if (!parser_initialized_) {
        initializeParser(ctx);  // ctx 在这里传入
    }
    
    // ❌ 此时 parser 中的指针指向已销毁的对象
    return parser_->Eval();  // 读取野指针，返回随机值！
}
```

#### 问题时序图

```
时刻 1: evaluate(ctx_temp) 被调用
   ├─> initializeParser(ctx_temp)
   │   ├─> parser_->DefineVar("T", &ctx_temp.T)
   │   │   // parser 记录了指针: T_ptr = 0x1234 (指向 ctx_temp.T)
   │   └─> return
   └─> ctx_temp 对象被销毁，地址 0x1234 失效

时刻 2: 第二次 evaluate(ctx_new) 被调用
   ├─> parser 已初始化，跳过 initializeParser
   ├─> parser_->Eval()
   │   └─> 读取 T_ptr (0x1234) ❌ 野指针！
   │       返回随机垃圾值，可能是负数
   └─> 返回错误结果
```

### 为什么 `ctx_temp` 会被销毁？

这是 C++ 函数参数传递的核心规则：

#### 1. 按值传递创建临时副本

```cpp
// 函数签名
double MaterialProperty::evaluate(const EvaluationContext& ctx) const;

// 调用方式
EvaluationContext ctx;
ctx.T = 300.0;
double result = property.evaluate(ctx);  // ctx 被传递
```

尽管参数是**引用**（`const EvaluationContext&`），但在某些调用场景下，编译器会创建**临时对象**：

```cpp
// 场景 1: 直接传入临时对象
property.evaluate(EvaluationContext{300.0});  
// ⚠️ 花括号创建临时对象，函数返回后立即销毁

// 场景 2: 隐式类型转换
property.evaluate({300.0});  
// ⚠️ 列表初始化创建临时对象

// 场景 3: 在循环中重复创建
for (int i = 0; i < n; ++i) {
    EvaluationContext ctx_local;  // 栈上局部对象
    ctx_local.T = temperature[i];
    property.evaluate(ctx_local);  
    // ⚠️ ctx_local 在每次循环结束时销毁
}
```

#### 2. 对象生命周期的三个关键时刻

```cpp
void someFunction() {
    {  // <-- 作用域开始
        EvaluationContext ctx;  // ① 构造：在栈上分配内存
        ctx.T = 300.0;
        
        // ctx.T 的地址，比如 0x7fff1234
        std::cout << "Address: " << &ctx.T << std::endl;
        
        property.evaluate(ctx);  // ② 使用：传递引用（地址）
        
    }  // <-- ③ 销毁：作用域结束，栈内存被回收
    
    // ❌ ctx 已不存在，地址 0x7fff1234 的内容不确定
}
```

#### 3. `initializeParser(ctx_temp)` 之后发生了什么

```cpp
// 实际的调用链
double computeConductivity(double delta_T) {
    EvaluationContext ctx;      // ① 在栈上创建 ctx
    ctx.T = T0_ + delta_T;      //    地址: 0x7fff5678
    
    // ② 第一次调用 evaluate
    return material->electrical_conductivity.evaluate(ctx);
    //     └─> initializeParser(ctx)
    //         └─> parser_->DefineVar("T", &ctx.T)  
    //             // 存储指针: 0x7fff5678
    //         └─> return
    //     └─> parser_->Eval()  // ✅ 此时 ctx 还活着
    //     └─> return 5.96e7
    
}  // ③ 函数结束，ctx 被销毁，栈空间回收

// 下一次调用（可能在不同的积分点）
double computeConductivity(double delta_T) {
    EvaluationContext ctx;      // ④ 在栈上创建新的 ctx
    ctx.T = T0_ + delta_T;      //    地址可能是: 0x7fff9abc (不同地址！)
    
    // ⑤ 第二次调用 evaluate
    return material->electrical_conductivity.evaluate(ctx);
    //     └─> parser 已初始化，跳过 initializeParser
    //     └─> parser_->Eval()
    //         └─> 读取 parser 中存储的指针: 0x7fff5678
    //             ❌ 这个地址已失效！返回垃圾值
}
```

#### 4. 栈内存重用导致的问题

```
第一次调用后的栈状态:
┌──────────────┬─────────┐
│ 地址         │ 内容    │
├──────────────┼─────────┤
│ 0x7fff5678   │ 300.0   │ ← ctx.T (第一次调用)
└──────────────┴─────────┘
parser 记录: T_ptr = 0x7fff5678 ✅

函数返回，栈空间被回收:
┌──────────────┬─────────┐
│ 地址         │ 内容    │
├──────────────┼─────────┤
│ 0x7fff5678   │ ???     │ ← 内存已释放，内容不确定
└──────────────┴─────────┘
parser 仍然记录: T_ptr = 0x7fff5678 ❌ 悬空！

第二次调用，栈可能被新数据覆盖:
┌──────────────┬─────────┐
│ 地址         │ 内容    │
├──────────────┼─────────┤
│ 0x7fff5678   │-3.9e8   │ ← 被其他变量占用的垃圾数据
└──────────────┴─────────┘
parser 读取: *T_ptr = -3.9e8 ❌ 错误的值！
```

### 关键认知误区：muParser 的指针绑定机制

#### ❌ 常见错误理解

很多人（包括最初的实现者）误以为：

> "每次调用 `evaluate(ctx)` 时，muParser 都会接收新的 `ctx`，自动更新变量值"

这是**完全错误**的理解！

#### ✅ muParser 的实际工作方式

muParser 使用**指针绑定**机制，分为两个阶段：

**阶段 1: 初始化（只执行一次）**
```cpp
parser->DefineVar("T", &some_variable);
// muParser 内部存储：
// variables_["T"] = 0x7fff5678  （指针地址）
// 
// ⚠️ 重点：muParser 只记录了指针，不会复制值！
```

**阶段 2: 每次评估**
```cpp
parser->Eval();
// muParser 内部执行：
// 1. 找到公式中的变量 "T"
// 2. 查找之前存储的指针: 0x7fff5678
// 3. 解引用读取值: value = *(0x7fff5678)
// 4. 计算公式: sigma0 / (1 + alpha * (value - T0))
//
// ⚠️ 重点：muParser 不知道 ctx 的存在，
//          它只知道那个地址 0x7fff5678！
```

#### 这意味着什么？

```cpp
// 第一次调用
EvaluationContext ctx1;
ctx1.T = 300.0;  // 假设地址 0x1000
parser->DefineVar("T", &ctx1.T);  // parser 记录: T_ptr = 0x1000
parser->Eval();  // 读取 *(0x1000) = 300.0 ✅

// 第二次调用 - 传入新的 ctx
EvaluationContext ctx2;
ctx2.T = 400.0;  // 假设地址 0x2000

parser->Eval();  
// ❌ 问题：parser 仍然读取 *(0x1000)
// 而不是 *(0x2000)！
// 
// 因为 DefineVar 只调用了一次，
// parser 不知道有新的 ctx2 存在！
```

#### 正确的使用方式有两种

**方式 1: 每次评估前更新指针指向的值**（我们采用的方案）
```cpp
class MaterialProperty {
    mutable std::map<std::string, double> variable_values_;
    // ✅ 变量值存储在对象内部，生命周期由对象管理
};

void initializeParser() {
    variable_values_["T"] = 0.0;
    parser->DefineVar("T", &variable_values_["T"]);
    // parser 记录: T_ptr = &(this->variable_values_["T"])
    // ✅ 这个指针在整个对象生命周期内有效
}

double evaluate(const EvaluationContext& ctx) {
    // ✅ 每次评估前，更新成员变量的值
    variable_values_["T"] = ctx.T;  // 300.0 -> variable_values_["T"]
    return parser->Eval();  // 读取 &variable_values_["T"] ✅
}
```

**方式 2: 使用持久化的上下文对象**（不推荐，耦合度高）
```cpp
class MaterialProperty {
    mutable EvaluationContext persistent_ctx_;  // 持久化的上下文
};

void initializeParser() {
    parser->DefineVar("T", &persistent_ctx_.T);
    // ✅ 指针指向成员，生命周期正确
}

double evaluate(const EvaluationContext& ctx) {
    // 复制整个上下文
    persistent_ctx_ = ctx;
    return parser->Eval();
}
```

### 为什么会返回负数？

野指针读取的内存可能被：
1. 新的栈帧覆盖
2. 堆分配重用
3. 包含垃圾数据（未初始化内存）

在本例中，读取到的值恰好是负数（`-3.919e+08`），这完全是**未定义行为**的结果。

## 调试方法

### 1. 添加详细日志

在 `computeJouleHeat` 中添加调试输出：

```cpp
double ElectrothermalSolver::computeJouleHeat(int element_idx,
                                              const std::vector<double>& ref_coords)
{
    double delta_T = interpolateTemperature(element_idx, ref_coords);
    double sigma = computeConductivity(delta_T);
    double grad_V_squared = electric_solver_->computeGradientNormSquared(element_idx, ref_coords);

    // ✅ 调试输出：检查各个物理量
    if (element_idx == 0 && iterations_ == 1 && verbose_) {
        std::cout << "  [调试详细] ΔT = " << delta_T 
                  << ", σ = " << sigma 
                  << ", |∇V|² = " << grad_V_squared << std::endl;
    }

    return sigma * grad_V_squared;
}
```

### 2. 检查 muParser 变量绑定

在 `MaterialProperty::evaluate` 中添加：

```cpp
double MaterialProperty::evaluate(const EvaluationContext& ctx) const
{
    if (!parser_initialized_) {
        initializeParser(ctx);
        
        // ✅ 调试：打印 parser 绑定的变量地址
        std::cout << "[Parser Debug] Variables bound:" << std::endl;
        for (const auto& var : variables_) {
            std::cout << "  " << var << " @ " 
                      << parser_->GetVar()[var] << std::endl;
        }
    }
    
    // ✅ 调试：打印评估前的变量值
    std::cout << "[Eval] T = " << ctx.T << std::endl;
    
    return parser_->Eval();
}
```

### 3. 使用 Valgrind/AddressSanitizer 检测内存错误

编译时启用地址消毒器：

```bash
cmake -DCMAKE_CXX_FLAGS="-fsanitize=address -g" ..
```

这会检测到：
```
ERROR: AddressSanitizer: stack-use-after-scope
READ of size 8 at 0x7fff12345678
```

### 4. 逐步验证公式

手动计算期望值并与实际输出对比：

```cpp
// 期望: σ(T=300) = 5.96e7 / (1 + 0.00393 × (300-293.15)) ≈ 5.80e7
double expected = 5.96e7 / (1 + 0.00393 * (300.0 - 293.15));
double actual = material->electrical_conductivity.evaluate(ctx);
std::cout << "Expected: " << expected << ", Actual: " << actual << std::endl;
```

## 解决方案

### 核心思路

**MaterialProperty 必须拥有变量值的存储空间**，而不是依赖外部临时对象。

### 修复后的实现

**头文件修改**:
```cpp
class MaterialProperty {
private:
    mutable std::shared_ptr<mu::Parser> parser_;
    mutable bool parser_initialized_;
    
    // ✅ 新增：存储变量值的副本
    mutable std::map<std::string, double> variable_values_;
};
```

**初始化修复**:
```cpp
void MaterialProperty::initializeParser(const EvaluationContext& ctx) const
{
    if (parser_initialized_) return;
    
    parser_ = std::make_shared<mu::Parser>();
    
    // 定义常数参数
    for (const auto& [name, value] : parameters_) {
        parser_->DefineConst(name, value);
    }
    
    // ✅ 关键修复：为每个变量创建自己的存储空间
    for (const auto& var_name : variables_) {
        variable_values_[var_name] = 0.0;  // 初始化
        parser_->DefineVar(var_name, &variable_values_[var_name]);
        // 现在指针指向 MaterialProperty 的成员，生命周期正确！
    }
    
    parser_->SetExpr(formula_);
    parser_initialized_ = true;
}
```

**评估修复**:
```cpp
double MaterialProperty::evaluate(const EvaluationContext& ctx) const
{
    if (type_ == Type::CONSTANT) {
        return constant_value_;
    }
    
    if (!parser_initialized_) {
        initializeParser(ctx);
    }
    
    // ✅ 关键修复：每次评估前更新变量值
    auto var_map = const_cast<EvaluationContext&>(ctx).getVariableMap();
    for (const auto& var_name : variables_) {
        if (var_map.find(var_name) != var_map.end()) {
            variable_values_[var_name] = *var_map[var_name];
            // 更新存储在 MaterialProperty 中的值
        }
    }
    
    return parser_->Eval();  // ✅ 现在读取的是有效指针
}
```

### 修复后的指针生命周期

```
时刻 1: evaluate(ctx_temp) 被调用
   ├─> initializeParser(ctx_temp)
   │   ├─> variable_values_["T"] = 0.0
   │   ├─> parser_->DefineVar("T", &variable_values_["T"])
   │   │   // parser 记录: T_ptr = &(this->variable_values_["T"])
   │   │   // 这个指针始终有效，因为指向 MaterialProperty 成员
   │   └─> return
   └─> ctx_temp 被销毁（没关系，我们没用它的指针）

时刻 2: evaluate(ctx_new) 被调用
   ├─> variable_values_["T"] = *ctx_new.getVariableMap()["T"]
   │   // 从新上下文复制值到成员变量
   ├─> parser_->Eval()
   │   └─> 读取 T_ptr (&variable_values_["T"]) ✅ 有效指针！
   │       返回正确计算结果
   └─> 返回 5.80e7 (正确！)
```

## 验证修复

### 运行结果对比

**修复前**:
```
[调试详细] ΔT = 0, σ = -3.919e+08, |∇V|² = 21.2487  ❌
焦耳热 Q = -8.32738e+09 W/m³  ❌
```

**修复后**:
```
[调试详细] ΔT = 0, σ = 5.96e+07, |∇V|² = 21.2487  ✅
焦耳热 Q = 1.266e+09 W/m³  ✅
```

### 物理合理性检查

- ✅ 电导率为正：`σ = 5.96e7 S/m`
- ✅ 焦耳热为正：`Q > 0`
- ✅ 温度依赖正确：σ 随温度升高而降低（α > 0）

## 经验教训

### 1. muParser 的使用陷阱

muParser 通过**指针**访问变量，要求：
- 指针在整个 parser 生命周期内有效
- 变量值通过指针更新，而非重新绑定

### 2. C++ 对象生命周期

- ⚠️ 临时对象在表达式结束后立即销毁
- ⚠️ 函数参数是临时副本（除非是引用）
- ✅ 类成员的生命周期与对象相同

### 3. 调试指针问题的工具

| 工具             | 用途                 |
| ---------------- | -------------------- |
| AddressSanitizer | 运行时检测内存错误   |
| Valgrind         | 内存泄漏和野指针检测 |
| 详细日志         | 追踪变量值和指针地址 |
| 断言检查         | 验证物理量的合理性   |

### 4. 防御性编程

```cpp
// ✅ 在关键位置添加合理性检查
double sigma = computeConductivity(delta_T);
assert(sigma > 0 && "Electrical conductivity must be positive!");

double Q = sigma * grad_V_squared;
assert(Q >= 0 && "Joule heat must be non-negative!");
```

## 相关问题排查清单

遇到类似"物理量符号错误"的问题时，按此顺序检查：

1. ✅ **物理公式正确性**：公式本身是否正确？
2. ✅ **输入数据合理性**：参数值是否在合理范围？
3. ✅ **表达式解析**：muParser 是否正确解析公式？
4. ✅ **变量绑定**：变量指针是否有效？（本次问题所在）
5. ✅ **数值稳定性**：计算过程是否溢出/下溢？
6. ✅ **单位一致性**：所有量是否使用一致单位？

## 参考资料

- [muParser 官方文档](http://beltoforion.de/en/muparser/)
- [C++ 对象生命周期](https://en.cppreference.com/w/cpp/language/lifetime)
- [AddressSanitizer 使用指南](https://github.com/google/sanitizers/wiki/AddressSanitizer)

## 修改文件清单

| 文件                                    | 修改内容                                  |
| --------------------------------------- | ----------------------------------------- |
| `include/material.h`                    | 添加 `variable_values_` 成员              |
| `src/core/material.cpp`                 | 修改 `initializeParser()` 和 `evaluate()` |
| `src/solvers/electrothermal_solver.cpp` | 添加调试输出                              |

---

**修复日期**: 2025-11-14  
**问题严重性**: 高（导致物理结果完全错误）  
**解决状态**: ✅ 已完全解决
