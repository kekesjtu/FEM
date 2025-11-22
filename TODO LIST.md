# TODO LIST

## 1. 目前硬编码为载荷向量为0，后续修改✅已完成
electrothermal_solver.cpp  
```cpp
 // 重置载荷向量为0（电场方程无源项）
        Eigen::VectorXd b_electric =
            Eigen::VectorXd::Zero(electric_solver_->getConfig()->getNodesNum());
        electric_solver_->setLoadVector(b_electric);
```


## 2.删除electrothermal_solver.cpp  中的所有`if (verbose_)`判断，直接打印日志。✅已完成

## 3. 后续将单场求解的assemble函数拆分成assembleStiffnessMatrix和assembleLoadVector两个函数，方便调用。✅已完成
electrothermal_solver.cpp 
```cpp
    // 组装热场刚度矩阵（需要热导率 k）
    // 我们需要调用 assemble() 来组装刚度矩阵
    // 但需要确保 thermal_solver_ 的 problem_ 有正确的系数和源项
    thermal_solver_->assemble();

    Eigen::VectorXd b_thermal_coupled = assembleThermalLoadCoupled();
    // 使用我们精确组装的载荷向量
    // 注意：thermal_solver_->assemble() 已经组装了b（源项为0），
    //      所以直接设置为我们的焦耳热载荷向量
    thermal_solver_->setLoadVector(b_thermal_coupled);
```

## 4.由温度场计算焦耳热密度场之后单独用一个函数封装✅已完成
electrothermal_solver.cpp
```cpp
  // 计算并输出焦耳热密度场
    // 在每个节点计算焦耳热密度（使用单元平均）
    int N = T_current_.size();
    Eigen::VectorXd Q_nodes = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd node_count = Eigen::VectorXd::Zero(N);

    auto config = electric_solver_->getConfig();
    const auto& connectivity = config->getElementConnectivity();
    int M = config->getElementsNum();
    int dimension = config->getDimension();

    for (int e = 0; e < M; ++e)
    {
        // 在单元中心计算焦耳热
        std::vector<double> center_ref(dimension, 1.0 / 3.0);  // 三角形中心
        double Q_center = computeJouleHeat(e, center_ref);

        // 分配到单元节点
        int n = connectivity[e].size();
        for (int i = 0; i < n; ++i)
        {
            int node_idx = connectivity[e][i];
            Q_nodes(node_idx) += Q_center;
            node_count(node_idx) += 1.0;
        }
    }

    // 节点平均
    for (int i = 0; i < N; ++i)
    {
        if (node_count(i) > 0)
        {
            Q_nodes(i) /= node_count(i);
        }
    }

    auto vtkHeat = VTKOutputFactory::createVTKOutput(config, Q_nodes);
    vtkHeat->outputNumericalSolution(prefix + "_Q");
    LOG_INFO("已输出: " + prefix + "_Q.vtu (焦耳热密度场)");

    LOG_INFO("");
    LOG_INFO("焦耳热密度范围: [" + std::to_string(Q_nodes.minCoeff()) + ", " +
             std::to_string(Q_nodes.maxCoeff()) + "] W/m³");
```
已经完成

## 5.将部分不由config类管理的东西移出config类，移到problem_setup ✅已完成
求解器参数（solver_type, preconditioner_type, tolerance, max_iterations）和日志参数（log_level, log_file）已移至ProblemSetup类管理。
fem_solver.cpp已修改为从`problem_->solver`和`problem_->log`读取参数。

## 6. 输出精确解功能需要更新VTK输出类以支持MaterialProperty ✅已完成
vtk_output_2d.cpp和vtk_output_3d.cpp的`outputExactSolution`方法已支持接收函数对象。
fem_solver.cpp使用lambda包装MaterialProperty：`[&field](coords) { return field.exact_solution_u.evaluate(ctx); }`
精确解VTK输出已经可以正常工作。

## 7.缩放因子支持不同单位的坐标数据 ✅已完成
实现了mesh_unit_scale机制：
- ProblemSetup添加`mesh_unit_scale`字段，从JSON读取
- Config类添加`setMeshUnitScale()`/`getMeshUnitScale()`方法
- ComsolMeshImporter在解析坐标时应用缩放：`node_coordinates_.push_back(coord * mesh_unit_scale_)`
- main.cpp从problem读取并传递给config：`config->setMeshUnitScale(problem->mesh_unit_scale)`
- JSON配置示例：`"mesh_unit_scale": 0.001` (mm→m转换)


## 8.json文件中，all、-1、`[1,2,3,4]`等施加方式需要统一。另外，有的地方实体使用entity字段，有的地方不使用，需要统一。