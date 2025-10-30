#include "problem_setup.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "config.h"

ProblemSetup::ProblemSetup(const std::string& name) : name_(name)
{
    // 默认构造函数
}

void ProblemSetup::setBoundarySetupFunction(std::function<void(std::shared_ptr<Config>)> func)
{
    boundary_setup_func_ = func;
}

const std::vector<BoundaryCondition>& ProblemSetup::getBoundaryConditions() const
{
    return boundary_conditions_;
}

const BoundaryCondition& ProblemSetup::getBoundaryCondition(size_t boundary_index) const
{
    if (boundary_index >= boundary_conditions_.size())
    {
        throw std::out_of_range("边界索引 " + std::to_string(boundary_index) +
                                " 超出范围 (问题: " + name_ + ")");
    }
    return boundary_conditions_[boundary_index];
}

std::vector<BoundaryCondition>& ProblemSetup::getBoundaryConditionsMutable()
{
    return boundary_conditions_;
}

void ProblemSetup::setCoefficient(std::function<double(const std::vector<double>&)> func)
{
    coefficient_func_ = func;
}

void ProblemSetup::setSource(std::function<double(const std::vector<double>&)> func)
{
    source_func_ = func;
}

double ProblemSetup::coefficient(const std::vector<double>& coords) const
{
    if (!coefficient_func_)
    {
        throw std::runtime_error("扩散系数函数未设置 (问题: " + name_ + ")");
    }
    return coefficient_func_(coords);
}

double ProblemSetup::source(const std::vector<double>& coords) const
{
    if (!source_func_)
    {
        throw std::runtime_error("源项函数未设置 (问题: " + name_ + ")");
    }
    return source_func_(coords);
}

void ProblemSetup::setExactSolutionU(std::function<double(const std::vector<double>&)> func)
{
    exact_solution_u_func_ = func;
}

void ProblemSetup::setExactSolutionGradients(
    std::function<double(const std::vector<double>&, std::vector<double>&)> func)
{
    exact_solution_gradients_func_ = func;
}

bool ProblemSetup::hasExactSolution() const
{
    return exact_solution_u_func_ != nullptr;
}
double ProblemSetup::exactSolutionU(const std::vector<double>& coords) const
{
    if (!exact_solution_u_func_)
    {
        static bool warning_shown = false;
        if (!warning_shown)
        {
            std::cerr << "警告: 精确解函数未设置 (问题: " << name_ << ")，返回 NaN" << std::endl;
            warning_shown = true;
        }
        return std::nan("");
    }
    return exact_solution_u_func_(coords);
}

double ProblemSetup::exactSolutionGradients(const std::vector<double>& coords,
                                            std::vector<double>& gradients) const
{
    if (!exact_solution_gradients_func_)
    {
        static bool warning_shown = false;
        if (!warning_shown)
        {
            std::cerr << "警告: 精确解梯度函数未设置 (问题: " << name_ << ")，返回 NaN"
                      << std::endl;
            warning_shown = true;
        }
        // 填充梯度向量为 NaN
        for (auto& grad : gradients)
        {
            grad = std::nan("");
        }
        return std::nan("");
    }
    return exact_solution_gradients_func_(coords, gradients);
}