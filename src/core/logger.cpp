#include "logger.h"
#include <algorithm>
#include <iostream>
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"


Logger::Logger() : current_level_(LogLevel::INFO)
{
    // 创建彩色控制台 sink
    console_sink_ = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

    // 创建 logger（初始只有控制台输出）
    spdlog_logger_ = std::make_shared<spdlog::logger>("fem_logger", console_sink_);

    // 设置默认格式（不显示时间戳，保持简洁）
    spdlog_logger_->set_pattern("%v");

    // 设置为默认logger
    spdlog::set_default_logger(spdlog_logger_);

    // 设置默认日志级别
    setLogLevel(LogLevel::INFO);
}

Logger::~Logger()
{
    spdlog::shutdown();
}

Logger& Logger::getInstance()
{
    static Logger instance;
    return instance;
}

void Logger::setLogLevel(LogLevel level)
{
    current_level_ = level;

    // 转换为 spdlog 的日志级别
    switch (level)
    {
        case LogLevel::DEBUG:
            spdlog_logger_->set_level(spdlog::level::debug);
            break;
        case LogLevel::INFO:
            spdlog_logger_->set_level(spdlog::level::info);
            break;
        case LogLevel::WARNING:
            spdlog_logger_->set_level(spdlog::level::warn);
            break;
        case LogLevel::ERR:
            spdlog_logger_->set_level(spdlog::level::err);
            break;
        case LogLevel::NONE:
            spdlog_logger_->set_level(spdlog::level::off);
            break;
    }
}

LogLevel Logger::getLogLevel() const
{
    return current_level_;
}

void Logger::setLogLevelFromString(const std::string& level_str)
{
    std::string upper_str = level_str;
    std::transform(upper_str.begin(), upper_str.end(), upper_str.begin(), ::toupper);

    if (upper_str == "DEBUG")
    {
        setLogLevel(LogLevel::DEBUG);
    }
    else if (upper_str == "INFO")
    {
        setLogLevel(LogLevel::INFO);
    }
    else if (upper_str == "WARNING")
    {
        setLogLevel(LogLevel::WARNING);
    }
    else if (upper_str == "ERROR" || upper_str == "ERR")
    {
        setLogLevel(LogLevel::ERR);
    }
    else if (upper_str == "NONE")
    {
        setLogLevel(LogLevel::NONE);
    }
    else
    {
        spdlog::warn("未知的日志级别 '{}', 使用默认级别 INFO", level_str);
        setLogLevel(LogLevel::INFO);
    }
}

void Logger::enableFileOutput(const std::string& filename)
{
    if (filename.empty())
    {
        disableFileOutput();
        return;
    }

    try
    {
        // 创建文件 sink
        file_sink_ = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, true);

        // 重新创建 logger，同时输出到控制台和文件
        spdlog::sinks_init_list sinks{console_sink_, file_sink_};
        spdlog_logger_ = std::make_shared<spdlog::logger>("fem_logger", sinks);

        // 恢复格式和日志级别
        spdlog_logger_->set_pattern("%v");
        setLogLevel(current_level_);

        // 设置为默认logger
        spdlog::set_default_logger(spdlog_logger_);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        spdlog::error("无法打开日志文件 {}: {}", filename, ex.what());
    }
}

void Logger::disableFileOutput()
{
    if (file_sink_)
    {
        // 重新创建 logger，只输出到控制台
        spdlog_logger_ = std::make_shared<spdlog::logger>("fem_logger", console_sink_);
        spdlog_logger_->set_pattern("%v");
        setLogLevel(current_level_);
        spdlog::set_default_logger(spdlog_logger_);

        file_sink_.reset();
    }
}

void Logger::debug(const std::string& message)
{
    spdlog_logger_->debug("[DEBUG] {}", message);
}

void Logger::info(const std::string& message)
{
    spdlog_logger_->info("{}", message);
}

void Logger::warning(const std::string& message)
{
    spdlog_logger_->warn("[WARNING] {}", message);
}

void Logger::error(const std::string& message)
{
    spdlog_logger_->error("[ERROR] {}", message);
}

void Logger::header(const std::string& title, LogLevel level)
{
    if (level < current_level_)
    {
        return;
    }

    std::string line(50, '=');
    std::string output;

    if (title.empty())
    {
        output = line;
    }
    else
    {
        // 居中显示标题
        size_t padding = (50 - title.length()) / 2;
        if (padding < 2)
            padding = 2;

        output = "\n" + line + "\n" + std::string(padding, ' ') + title + "\n" + line;
    }

    spdlog_logger_->info("{}", output);
}

void Logger::subheader(const std::string& title, LogLevel level)
{
    if (level < current_level_)
    {
        return;
    }

    std::string line(40, '-');
    std::string output;

    if (title.empty())
    {
        output = line;
    }
    else
    {
        output = "\n" + title + "\n" + line;
    }

    spdlog_logger_->info("{}", output);
}

void Logger::separator(LogLevel level)
{
    if (level < current_level_)
    {
        return;
    }

    std::string output(40, '-');
    spdlog_logger_->info("{}", output);
}

void Logger::flush()
{
    spdlog_logger_->flush();
}
