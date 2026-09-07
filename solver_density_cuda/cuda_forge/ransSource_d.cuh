#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

void ransSource_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// SST ブレンド関数 F1 を現在状態 (k, ω, ∇k·∇ω, wall_dist) から `sstF1` に書く前処理。ransGradient の後・
// ransTransport (σ ブレンド拡散) の前に呼び、拡散と生成が同一 step の F1 を使うようにする (ラグ除去)。
void ransBlendF1_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);