// Optional KD-Tree based wall distance calculation.
// If k-d tree library is not available (HAVE_KDTREE undefined), a brute force
// fallback will be used. This removes the previous hard-coded absolute path.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#ifdef HAVE_KDTREE
#include <kdtree.h>
#include <memory>
#endif

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

// 壁距離の点 (壁の点・評価点)。ビルドの型に依らず double (plan architecture-float-state-double-geometry §4.6)。
struct Point {
  double x, y, z;
  Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

// 壁距離の評価に使う位置 (double)。nullptr の配列は msh の geom_float の値を広げて使う。
// 変換器は幾何の正本 (gmshReader::geo64) の配列を渡し、/VALUE/wall_dist を戻り値から直接書く (§4.6)。
struct WallDistPositions64 {
  const double *nodeCoord = nullptr; // [3*msh.nodes.size()]  節点座標
  const double *cellCent  = nullptr; // [3*msh.cells.size()]  セル (CV) 重心
  const double *planeCent = nullptr; // [3*msh.planes.size()] 面重心
};

// 壁距離の本体: 壁の点と評価点を選んで (選び方は node/cell・wallDistExtraPhysIDs とも従来どおり)、
// double の点群から double の配列 [msh.cells.size()] を返す。壁が無ければ [msh.nCells_all] の 0。
// HAVE_KDTREE の有無・FORGE_WALLDIST_BRUTE のどの経路も同じく double で計算する。
std::vector<double> calcWallDistance64(const solverConfig &cfg, const mesh &msh, const WallDistPositions64 &pos);

// ソルバ用の写し: calcWallDistance64 (msh の geom_float の位置) の結果を var.c["wall_dist"] (flow_float) に丸めて入れ、
// デバイスへ送る。
void calcWallDistance_kdtree(solverConfig &cfg, mesh &msh, variables &var);
