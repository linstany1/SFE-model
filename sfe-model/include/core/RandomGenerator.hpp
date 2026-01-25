/**
 * @file RandomGenerator.hpp
 * @brief 随机数生成器封装类，确保模拟的确定性可复现
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 关键要求：
 * - 使用 std::mt19937 作为核心引擎
 * - 所有随机过程使用统一 RNG
 * - 相同 seed 必须产生完全一致的输出
 * 
 * 提供的接口：
 * - getUniform(min, max): 均匀分布随机数
 * - getBernoulli(p): 伯努利试验
 * - getInt(min, max): 整数随机数
 * - getGaussian(mean, stddev): 高斯分布
 */

#ifndef SFE_RANDOM_GENERATOR_HPP
#define SFE_RANDOM_GENERATOR_HPP

#include <random>
#include <cstdint>
#include <mutex>

namespace sfe {

/**
 * @class RandomGenerator
 * @brief 线程安全的随机数生成器单例类
 * 
 * 使用 Mersenne Twister (mt19937) 算法，提供高质量的伪随机数。
 * 通过固定种子确保模拟结果的可复现性。
 */
class RandomGenerator {
public:
    /**
     * @brief 获取单例实例
     * @return RandomGenerator 引用
     */
    static RandomGenerator& getInstance() {
        static RandomGenerator instance;
        return instance;
    }
    
    // 禁止拷贝和移动
    RandomGenerator(const RandomGenerator&) = delete;
    RandomGenerator& operator=(const RandomGenerator&) = delete;
    RandomGenerator(RandomGenerator&&) = delete;
    RandomGenerator& operator=(RandomGenerator&&) = delete;
    
    // ========================================================================
    // 初始化
    // ========================================================================
    
    /**
     * @brief 设置随机数种子
     * @param seed 种子值
     * 
     * 必须在模拟开始前调用，确保可复现性
     */
    void setSeed(uint64_t seed) {
        std::lock_guard<std::mutex> lock(mutex_);
        engine_.seed(static_cast<std::mt19937::result_type>(seed));
        seed_ = seed;
    }
    
    /**
     * @brief 获取当前种子
     * @return 种子值
     */
    uint64_t getSeed() const {
        return seed_;
    }
    
    // ========================================================================
    // 均匀分布
    // ========================================================================
    
    /**
     * @brief 生成 [0, 1) 范围内的均匀分布随机数
     * @return double 随机数
     */
    double getUniform01() {
        std::lock_guard<std::mutex> lock(mutex_);
        return uniform_01_(engine_);
    }
    
    /**
     * @brief 生成 [min, max) 范围内的均匀分布随机数
     * @param min 最小值（包含）
     * @param max 最大值（不包含）
     * @return double 随机数
     */
    double getUniform(double min, double max) {
        std::lock_guard<std::mutex> lock(mutex_);
        std::uniform_real_distribution<double> dist(min, max);
        return dist(engine_);
    }
    
    // ========================================================================
    // 整数分布
    // ========================================================================
    
    /**
     * @brief 生成 [min, max] 范围内的均匀分布整数
     * @param min 最小值（包含）
     * @param max 最大值（包含）
     * @return int 随机整数
     */
    int getInt(int min, int max) {
        std::lock_guard<std::mutex> lock(mutex_);
        std::uniform_int_distribution<int> dist(min, max);
        return dist(engine_);
    }
    
    // ========================================================================
    // 伯努利试验
    // ========================================================================
    
    /**
     * @brief 执行伯努利试验
     * @param probability 成功概率 [0, 1]
     * @return true 如果成功，false 如果失败
     * 
     * 用于死亡判定、定居判定等概率事件
     */
    bool getBernoulli(double probability) {
        if (probability <= 0.0) return false;
        if (probability >= 1.0) return true;
        
        std::lock_guard<std::mutex> lock(mutex_);
        return uniform_01_(engine_) < probability;
    }
    
    /**
     * @brief 执行伯努利试验（别名，语义更清晰）
     * @param probability 成功概率 [0, 1]
     * @return true 如果事件发生
     */
    bool trial(double probability) {
        return getBernoulli(probability);
    }
    
    // ========================================================================
    // 高斯（正态）分布
    // ========================================================================
    
    /**
     * @brief 生成高斯分布随机数
     * @param mean 均值
     * @param stddev 标准差
     * @return double 随机数
     */
    double getGaussian(double mean, double stddev) {
        std::lock_guard<std::mutex> lock(mutex_);
        std::normal_distribution<double> dist(mean, stddev);
        return dist(engine_);
    }
    
    // ========================================================================
    // 泊松分布
    // ========================================================================
    
    /**
     * @brief 生成泊松分布随机数
     * @param lambda 期望值
     * @return int 随机数
     */
    int getPoisson(double lambda) {
        if (lambda <= 0.0) return 0;
        
        std::lock_guard<std::mutex> lock(mutex_);
        std::poisson_distribution<int> dist(lambda);
        return dist(engine_);
    }
    
    // ========================================================================
    // 指数分布
    // ========================================================================
    
    /**
     * @brief 生成指数分布随机数
     * @param lambda 速率参数 (1/mean)
     * @return double 随机数
     */
    double getExponential(double lambda) {
        if (lambda <= 0.0) return 0.0;
        
        std::lock_guard<std::mutex> lock(mutex_);
        std::exponential_distribution<double> dist(lambda);
        return dist(engine_);
    }
    
    // ========================================================================
    // 工具方法
    // ========================================================================
    
    /**
     * @brief 在 1m x 1m 网格内生成随机坐标
     * @param cell_x 网格 x 索引
     * @param cell_y 网格 y 索引
     * @param[out] x 生成的 x 坐标
     * @param[out] y 生成的 y 坐标
     * 
     * 用于幼苗定居时分配随机位置
     */
    void getRandomPositionInCell(int cell_x, int cell_y, double& x, double& y) {
        std::lock_guard<std::mutex> lock(mutex_);
        x = cell_x + uniform_01_(engine_);
        y = cell_y + uniform_01_(engine_);
    }
    
    /**
     * @brief 在指定范围内生成不重叠的随机坐标
     * @param cell_x 网格 x 索引
     * @param cell_y 网格 y 索引
     * @param existing_positions 已存在的位置列表 [(x1,y1), (x2,y2), ...]
     * @param min_distance 最小间距
     * @param[out] x 生成的 x 坐标
     * @param[out] y 生成的 y 坐标
     * @param max_attempts 最大尝试次数
     * @return true 如果成功找到位置
     */
    bool getRandomPositionWithMinDistance(
        int cell_x, int cell_y,
        const std::vector<std::pair<double, double>>& existing_positions,
        double min_distance,
        double& x, double& y,
        int max_attempts = 100
    ) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        double min_dist_sq = min_distance * min_distance;
        
        for (int attempt = 0; attempt < max_attempts; ++attempt) {
            double test_x = cell_x + uniform_01_(engine_);
            double test_y = cell_y + uniform_01_(engine_);
            
            bool valid = true;
            for (const auto& pos : existing_positions) {
                double dx = test_x - pos.first;
                double dy = test_y - pos.second;
                if (dx * dx + dy * dy < min_dist_sq) {
                    valid = false;
                    break;
                }
            }
            
            if (valid) {
                x = test_x;
                y = test_y;
                return true;
            }
        }
        
        // 如果找不到合适位置，返回网格中心
        x = cell_x + 0.5;
        y = cell_y + 0.5;
        return false;
    }
    
    /**
     * @brief 获取底层引擎（仅用于高级用途）
     * @return mt19937 引擎引用
     * 
     * 注意：直接使用引擎可能破坏线程安全性
     */
    std::mt19937& getEngine() {
        return engine_;
    }

private:
    RandomGenerator() : seed_(42), engine_(42), uniform_01_(0.0, 1.0) {}
    
    uint64_t seed_;                              ///< 当前种子
    std::mt19937 engine_;                        ///< Mersenne Twister 引擎
    std::uniform_real_distribution<double> uniform_01_; ///< [0,1) 分布缓存
    std::mutex mutex_;                           ///< 线程安全互斥锁
};

// ============================================================================
// 便捷全局函数（直接调用单例）
// ============================================================================

/**
 * @brief 便捷函数：获取 [0, 1) 均匀分布随机数
 */
inline double rngUniform01() {
    return RandomGenerator::getInstance().getUniform01();
}

/**
 * @brief 便捷函数：获取 [min, max) 均匀分布随机数
 */
inline double rngUniform(double min, double max) {
    return RandomGenerator::getInstance().getUniform(min, max);
}

/**
 * @brief 便捷函数：执行伯努利试验
 */
inline bool rngBernoulli(double p) {
    return RandomGenerator::getInstance().getBernoulli(p);
}

/**
 * @brief 便捷函数：获取 [min, max] 均匀分布整数
 */
inline int rngInt(int min, int max) {
    return RandomGenerator::getInstance().getInt(min, max);
}

} // namespace sfe

#endif // SFE_RANDOM_GENERATOR_HPP
