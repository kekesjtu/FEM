#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

// 前向声明 spdlog 类型
namespace spdlog
{
class logger;
namespace sinks
{
class sink;
}
}  // namespace spdlog

/**
 * @brief 日志级别枚举
 */
enum class LogLevel
{
    DEBUG = 0,    // 详细的调试信息
    INFO = 1,     // 一般信息
    WARNING = 2,  // 警告信息
    ERR = 3,      // 错误信息
    NONE = 4      // 不输出任何日志
};

/**
 * @brief 将字符串转换为LogLevel
 */
inline LogLevel stringToLogLevel(const std::string& level_str)
{
    if (level_str == "DEBUG")
        return LogLevel::DEBUG;
    if (level_str == "INFO")
        return LogLevel::INFO;
    if (level_str == "WARNING")
        return LogLevel::WARNING;
    if (level_str == "ERROR")
        return LogLevel::ERR;
    if (level_str == "NONE")
        return LogLevel::NONE;
    return LogLevel::INFO;  // 默认
}

/**
 * @brief 统一的日志输出工具类（基于 spdlog）
 *
 * 提供分级日志输出功能，支持：
 * - 四个日志级别：DEBUG, INFO, WARNING, ERROR
 * - 彩色控制台输出
 * - 文件输出
 * - 高性能异步日志
 * - 线程安全（单例模式）
 */
class Logger
{
  public:
    /**
     * @brief 获取Logger单例
     */
    static Logger& getInstance();

    /**
     * @brief 设置日志级别
     * @param level 日志级别
     */
    void setLogLevel(LogLevel level);

    /**
     * @brief 获取当前日志级别
     */
    LogLevel getLogLevel() const;

    /**
     * @brief 设置日志级别（从字符串）
     * @param level_str 日志级别字符串："DEBUG", "INFO", "WARNING", "ERROR", "NONE"
     */
    void setLogLevelFromString(const std::string& level_str);

    /**
     * @brief 启用/禁用文件输出
     * @param filename 日志文件名（空字符串表示禁用）
     */
    void enableFileOutput(const std::string& filename = "");

    /**
     * @brief 禁用文件输出
     */
    void disableFileOutput();

    /**
     * @brief 输出调试信息
     */
    void debug(const std::string& message);

    /**
     * @brief 输出一般信息
     */
    void info(const std::string& message);

    /**
     * @brief 输出警告信息
     */
    void warning(const std::string& message);

    /**
     * @brief 输出错误信息
     */
    void error(const std::string& message);

    /**
     * @brief 输出标题分隔线（主标题，50字符宽）
     * @param title 标题文本（可选）
     * @param level 日志级别（默认INFO）
     */
    void header(const std::string& title = "", LogLevel level = LogLevel::INFO);

    /**
     * @brief 输出次级标题分隔线（40字符宽）
     * @param title 标题文本（可选）
     * @param level 日志级别（默认INFO）
     */
    void subheader(const std::string& title = "", LogLevel level = LogLevel::INFO);

    /**
     * @brief 输出简单分隔线
     * @param level 日志级别（默认DEBUG）
     */
    void separator(LogLevel level = LogLevel::DEBUG);

    /**
     * @brief 刷新输出缓冲
     */
    void flush();

    // 禁止拷贝和赋值
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

  private:
    Logger();
    ~Logger();

    LogLevel current_level_;  // 当前日志级别

    // spdlog 相关
    std::shared_ptr<spdlog::logger> spdlog_logger_;
    std::shared_ptr<spdlog::sinks::sink> console_sink_;
    std::shared_ptr<spdlog::sinks::sink> file_sink_;
};

// ============================================================================
// 便捷的全局日志函数
// ============================================================================

/**
 * @brief 全局日志输出函数（简化调用）
 */
inline void LOG_DEBUG(const std::string& message)
{
    Logger::getInstance().debug(message);
}

inline void LOG_INFO(const std::string& message)
{
    Logger::getInstance().info(message);
}

inline void LOG_WARNING(const std::string& message)
{
    Logger::getInstance().warning(message);
}

inline void LOG_ERROR(const std::string& message)
{
    Logger::getInstance().error(message);
}

inline void LOG_HEADER(const std::string& title = "")
{
    Logger::getInstance().header(title);
}

inline void LOG_SUBHEADER(const std::string& title = "")
{
    Logger::getInstance().subheader(title);
}

inline void LOG_SEPARATOR()
{
    Logger::getInstance().separator();
}

// ============================================================================
// 性能计时工具
// ============================================================================

/**
 * @brief RAII 风格的作用域计时器
 *
 * 使用示例:
 *   {
 *       TIMER_SCOPE("矩阵组装");
 *       // ... 需要计时的代码
 *   }  // 自动输出耗时
 */
class ScopedTimer
{
  public:
    /**
     * @brief 构造函数，记录起始时间
     * @param name 计时任务名称
     * @param level 输出日志级别（默认INFO）
     */
    ScopedTimer(const std::string& name, LogLevel level = LogLevel::INFO)
        : name_(name), level_(level), start_(std::chrono::high_resolution_clock::now())
    {
    }

    /**
     * @brief 析构函数，自动计算并输出耗时
     */
    ~ScopedTimer()
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<double, std::milli>(end - start_).count();

        std::ostringstream oss;
        oss << name_ << " 耗时: " << std::fixed << std::setprecision(2) << duration << " ms";

        // 根据级别输出
        Logger& logger = Logger::getInstance();
        if (level_ == LogLevel::DEBUG)
        {
            logger.debug(oss.str());
        }
        else if (level_ == LogLevel::INFO)
        {
            logger.info(oss.str());
        }
        else if (level_ == LogLevel::WARNING)
        {
            logger.warning(oss.str());
        }
        else if (level_ == LogLevel::ERR)
        {
            logger.error(oss.str());
        }
    }

    // 禁止拷贝
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

  private:
    std::string name_;
    LogLevel level_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
};

// 便捷宏
#define TIMER_SCOPE(name) ScopedTimer __timer_##__LINE__(name, LogLevel::INFO)
#define TIMER_SCOPE_DEBUG(name) ScopedTimer __timer_##__LINE__(name, LogLevel::DEBUG)

#endif  // LOGGER_H
