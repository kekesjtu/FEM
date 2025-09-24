setting:
系统默认设置（不可修改）-用户设置-工作区设置-文件夹设置

后者的设置会覆盖前者的设置，若没有设置某一项，将继续使用前者的设置。



我们可以这样理解此层次

用户设置即全局设置，用户自行设定好后，每次打开VSCode即使用的此设定，若某项无设定即使用默认设置。

工作区设置即工作环境设置，可对不同的工作环境是用不同的工作环境，若某项无设定，即使用上一层设置。

文件夹设置即为项目设置，将一个文件夹当成一个项目，对同一个工作环境下的不同项目，使用不同的设置，若某项无设定，即使用上一层设置。

## VS Code 配置文件详解

### .vscode 文件夹
VS Code会在项目根目录创建`.vscode`文件夹，用于存储项目特定的配置文件。主要包含以下配置文件：

### tasks.json - 任务配置文件
`tasks.json`文件定义了可以在VS Code中运行的自动化任务，如编译、构建、测试等。

**主要作用：**
- 定义编译命令和参数
- 配置构建过程
- 设置自动化脚本
- 配置问题匹配器（识别编译错误）

**常见配置项：**
```json
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "编译项目",           // 任务名称
            "type": "shell",             // 任务类型：shell、process
            "command": "g++",            // 要执行的命令
            "args": ["-g", "main.cpp"],  // 命令参数
            "group": "build",            // 任务组：build、test
            "problemMatcher": ["$gcc"],  // 问题匹配器
            "detail": "编译C++项目"       // 任务描述
        }
    ]
}
```

### launch.json - 调试配置文件
`launch.json`文件配置调试器的启动参数，定义如何运行和调试程序。

**主要作用：**
- 配置调试器类型和参数
- 设置程序启动方式
- 定义环境变量
- 指定预启动任务（如编译）

**常见配置项：**
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "调试程序",                    // 配置名称
            "type": "cppdbg",                     // 调试器类型
            "request": "launch",                  // 启动类型：launch、attach
            "program": "${workspaceFolder}/app.exe", // 程序路径
            "args": [],                           // 程序参数
            "stopAtEntry": false,                 // 是否在入口处停止
            "cwd": "${workspaceFolder}",          // 工作目录
            "environment": [],                    // 环境变量
            "MIMode": "gdb",                      // 调试器模式
            "preLaunchTask": "编译项目"            // 启动前执行的任务
        }
    ]
}
```

### 两者的关系

**简单理解：**
- **tasks.json** 决定exe文件怎么生成（编译），或者简单运行
- **launch.json** 决定exe文件怎么被运行（调试运行。可配置断点）

**详细说明：**
- **tasks.json** 定义"怎么做"：如何编译、如何构建
- **launch.json** 定义"怎么运行"：如何启动程序、如何调试
- **协作关系**：launch.json可以通过`preLaunchTask`属性调用tasks.json中定义的任务

### 实际应用场景

**编译和运行流程：**
1. 按F5或点击运行按钮
2. VS Code检查launch.json中的`preLaunchTask`
3. 自动执行tasks.json中对应的编译任务
4. 编译成功后，启动程序进行调试

**常用快捷键：**
- `F5` - 启动调试
- `Ctrl+F5` - 运行而不调试
- `Ctrl+Shift+P` → "Tasks: Run Task" - 手动运行任务

**右上角运行按钮：**
- **▶️ 运行**：等价于 `Ctrl+F5`，快速运行程序（无调试）
- **🐛 运行和调试**：等价于 `F5`，启动调试模式
- 两个按钮都会调用 launch.json 配置和预启动任务

## 深入理解配置文件参数

### tasks.json 参数详解

#### 基本结构
```json
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "任务名称",
            "type": "shell",  // 或 "cppbuild"
            "command": "编译器路径",
            "args": ["参数数组"]
        }
    ]
}
```

#### 核心参数说明

**`command` + `args` = 完整终端命令**
- `command`: 要执行的程序路径
- `args`: 参数数组，每个元素是一个参数

例如终端命令：
```bash
g++ -g -std=c++11 -I /path main.cpp -o program.exe
```

在tasks.json中表示为：
```json
{
    "command": "g++",
    "args": [
        "-g",
        "-std=c++11", 
        "-I", "/path",      // -I 只影响紧跟的参数
        "main.cpp",
        "-o", "program.exe"
    ]
}
```

**任务类型 (`type`)**
- `"shell"`: 通用Shell任务，直接在系统终端执行
- `"cppbuild"`: C++专用构建任务，提供更好的错误解析和IDE集成
- `"process"`: 直接启动进程，不通过shell

**任务组 (`group`)**
- `"build"`: 构建组，可通过Ctrl+Shift+B快速运行
- `"test"`: 测试组，用于运行测试
- `{"kind": "build", "isDefault": true}`: 设置为默认构建任务

**问题匹配器 (`problemMatcher`)**
- `["$gcc"]`: 识别GCC编译器的错误格式
- 自动将编译错误显示在VS Code的问题面板
- 提供错误高亮和快速跳转功能
- `cppbuild`类型会自动配置，`shell`类型需要手动指定

**依赖关系 (`dependsOn`)**
```json
{
    "label": "运行程序",
    "command": "./program.exe",
    "dependsOn": "编译程序"  // 运行前先执行编译任务
}
```

**运行选项 (`options`)**
```json
{
    "options": {
        "cwd": "${workspaceFolder}",  // 工作目录
        "env": {                      // 环境变量
            "PATH": "/custom/path"
        },
        "shell": {                    // Shell配置
            "executable": "cmd.exe",
            "args": ["/d", "/c", "chcp 65001 &&"]  // 解决中文乱码
        }
    }
}
```

### launch.json 参数详解

#### 调试器配置

**调试器类型 (`type`)**
- `"cppdbg"`: C++调试器，用于调试C/C++程序

**启动类型 (`request`)**
- `"launch"`: 启动新进程进行调试
- `"attach"`: 附加到已运行的进程

**目标程序 (`program`)**
- 指定要调试的可执行文件的完整路径

**程序参数 (`args`)**
- 传递给程序的命令行参数数组
- 例：`["input.txt", "--verbose"]`

#### 调试行为控制

**入口停止 (`stopAtEntry`)**
- `false`: 程序直接运行到断点或结束
- `true`: 程序在main函数第一行暂停

**工作目录 (`cwd`)**
- 设置程序运行时的当前工作目录
- 影响相对路径解析和文件操作

**调试器模式 (`MIMode`)**
- `"gdb"`: 使用GDB调试器（Linux/MSYS2）
- `"lldb"`: 使用LLDB调试器（macOS）

**调试器路径 (`miDebuggerPath`)**
- 指定调试器可执行文件的完整路径
- 确保VS Code能找到正确的调试器

**环境变量 (`environment`)**
```json
"environment": [
    {"name": "DEBUG_MODE", "value": "1"},
    {"name": "PATH", "value": "/custom/path"}
]
```

**预启动任务 (`preLaunchTask`)**
- 调试前自动执行的任务名称
- 通常用于编译程序

## 运行方式对比

### 两种运行方式的区别

| 运行方式 | 通过tasks.json | 通过launch.json |
|----------|----------------|------------------|
| **触发方式** | 手动选择任务 | F5, Ctrl+F5 |
| **运行类型** | 直接命令行执行 | 调试器启动 |
| **适用场景** | 简单运行、测试脚本 | 开发调试 |
| **功能特点** | 基础运行 | 断点、变量监视 |
| **错误处理** | 终端输出 | 智能错误分析 |

### 字符编码问题

**问题现象：**
- tasks.json运行可能出现中文乱码
- launch.json调试运行通常没有乱码

**原因分析：**
- tasks.json直接使用系统终端（可能是GBK编码）
- launch.json通过调试器运行（有更好的编码处理）

**解决方案：**
```json
{
    "options": {
        "shell": {
            "executable": "cmd.exe",
            "args": ["/d", "/c", "chcp 65001 &&"]  // 设置UTF-8编码
        }
    }
}
```

## 配置最佳实践

### 项目配置建议

1. **使用cppbuild类型进行编译**
   - 更好的错误解析
   - 自动配置问题匹配器
   - 与VS Code深度集成

2. **设置合理的任务依赖**
   - 运行任务依赖编译任务
   - 确保运行的是最新版本

3. **配置正确的路径和环境**
   - 使用绝对路径指定编译器
   - 设置正确的工作目录和包含路径

4. **处理字符编码问题**
   - 在shell选项中设置UTF-8编码
   - 确保中文输出正常显示

### 常见变量

- `${workspaceFolder}`: 工作区根目录
- `${file}`: 当前打开的文件
- `${fileBasenameNoExtension}`: 不带扩展名的文件名
- `${fileDirname}`: 当前文件所在目录

这些配置文件让VS Code成为一个强大的C++开发环境，通过合理配置可以实现一键编译、调试和运行！