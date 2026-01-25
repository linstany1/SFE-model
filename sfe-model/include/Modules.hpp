/**
 * @file Modules.hpp
 * @brief 过程模块统一包含头文件
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * 包含所有过程模块：
 * - ClimateData: 气候数据和降尺度
 * - WaterBalance: 水分平衡计算
 * - SpatialHash: 空间哈希索引
 * - LightEngine: 光照计算引擎
 * - SeedDispersal: 种子散布
 * - Establishment: 幼苗定居
 * - Mortality: 死亡率计算
 * - Growth: 生长计算
 */

#ifndef SFE_MODULES_HPP
#define SFE_MODULES_HPP

#include "modules/ClimateData.hpp"
#include "modules/WaterBalance.hpp"
#include "modules/SpatialHash.hpp"
#include "modules/LightEngine.hpp"
#include "modules/SeedDispersal.hpp"
#include "modules/Establishment.hpp"
#include "modules/Mortality.hpp"
#include "modules/Growth.hpp"

#endif // SFE_MODULES_HPP
