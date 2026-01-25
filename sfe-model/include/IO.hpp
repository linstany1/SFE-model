/**
 * @file IO.hpp
 * @brief I/O 模块统一包含头文件
 * 
 * 亚高山森林演替模型 (SFE-Model)
 * 
 * Phase 2: 数据结构与 I/O 解析
 * 
 * 包含所有文件解析器：
 * - FileParser: 通用文件解析工具
 * - ClimateParser: 气候数据解析 (climate.txt)
 * - SpeciesParser: 物种参数解析 (species_params.txt)
 * - SiteParser: 站点信息解析 (site.txt)
 * - SeedSourceParser: 种源组成解析 (seed_source_composition.txt)
 * - DistanceMatrixParser: 距离矩阵解析 (distance_matrix.txt)
 */

#ifndef SFE_IO_HPP
#define SFE_IO_HPP

#include "io/FileParser.hpp"
#include "io/SpeciesProfile.hpp"
#include "io/ClimateParser.hpp"
#include "io/SpeciesParser.hpp"
#include "io/SiteParser.hpp"
#include "io/DistanceMatrixParser.hpp"

#endif // SFE_IO_HPP
