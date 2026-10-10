#pragma once

// node の壁関数 (SST wallTreatmentSST 1、ransWallFunction_d.cu) と WMLES の取得層 (wmlesWallModel_d.cu) が使う
// 代表内点と壁法線距離 (plans/active/architecture-float-state-double-geometry.md §4.2 3.、段 ②)。
//
// 境界面 ib ごとに、壁ノード W の入射内部双対面の相手のうち内部ノード (wall_flag = 0) で、壁内向き法線 −n̂ との
// cos が最大のもの I を選ぶ (SU2 Normal_Neighbor 流。|d| ≤ kSmall と d·(−n̂) ≤ 0 は候補にしない、同率は CSR で先のもの)。
// d = x_I − x_W、dn = −d·n̂。n̂ はデバイスの面ベクトル (sx..ss、軸対称 method 0 では半径の重みを掛けた後) から作る。
// 選び方・候補なしの扱い・同率の扱いは段 ① までのカーネル内の選択と同じで、変わるのは d の出どころだけ:
// 座標の差を double の値の位置 (mesh::cc64) で取り、dn・|d|・cos を double で求めて 1 回だけ flow_float に丸める。
// 幾何は実行中に変わらないので、読み込み後に 1 回だけ作る (従来は毎 step カーネルが選び直していた)。
// 角では同じ壁ノードでも面ごとに n̂ が違うので、壁ノードごとでなく境界面ごとに持つ (conjugateWall の合算した法線は使わない)。

#include <vector>

#include "cuda_forge/cudaConfig.cuh"
#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

struct WallRepPoints {
    geom_int nb = 0;
    geom_int*   irep_d = nullptr;   // [nb] 代表内点 I (候補なしは -1)
    flow_float* dn_d   = nullptr;   // [nb] 壁内向き距離 −(x_I − x_W)·n̂ (候補なしは kSmall)。カーネルは max(dn, kSmall) を y に使う
    flow_float* dist_d = nullptr;   // [nb] |x_I − x_W| (診断 rep_dist・rep_toff 用。候補なしは kSmall)
    flow_float* cos_d  = nullptr;   // [nb] dn/|d| (診断 rep_cos 用。候補なしは -2)
};

// 本番の代表点。node (wall_flag_d あり) の壁 bcond (wall / wall_isothermal) について、初回の呼び出しで
// すべての壁 bcond の分を作ってデバイスに置く (double の値の位置を一時的にデバイスへ上げ、作ったら解放する)。
// 壁でない bcond・cell では nb = 0 の空を返す。
const WallRepPoints& wallRepPoints(cudaConfig& cuda_cfg, mesh& msh, variables& var, const bcond& bc);

// 診断用 (FORGE_DIAG_GEOM_STAGE2_DUMP、plan §6.3): 座標の出どころを選んで作り直し、ホストに返す (本番の配列は触らない)。
//   legacy = false: 本番と同じ (double の値の位置)。true: 段 ① までのカーネルと同じ (デバイスの ccx..ccz を flow_float のまま引く)。
void wallRepPointsHost(cudaConfig& cuda_cfg, mesh& msh, variables& var, const bcond& bc, bool legacy,
                       std::vector<geom_int>& irep, std::vector<flow_float>& dn,
                       std::vector<flow_float>& dist, std::vector<flow_float>& cosv);
