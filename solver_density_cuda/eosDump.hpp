#pragma once
// 診断 (env `FORGE_DUMP_EOS_STEP=<k>` と `FORGE_DUMP_EOS_FILE=<path.h5>`、既定 off・出力専用で数値は変えない)。
// 凍結入力に対する EOS 1 回の全出力のビット比較 (plans/active/tooling-sern-te-wake-grid.md §5.1 #2 の追加の受入れ試験) 用。
//
// step k (1 起点 = res_<k> と同じ番号) の最初の assembleResidualPre の中で、境界ピン (no-slip・軸・等温壁) と
// 化学種・受動種の原始量の後、dependentVariables (EOS) の直前と直後に、var.c_d の全 cell 配列
// (壁・ghost を除外せず nCells_all 全体、型そのまま float32) を `/pre/<名前>`・`/post/<名前>` に書く。
// 物性 DB (device の SpeciesThermo / SpeciesThermoF を詰め物なしでバイト列に詰めたもの) を `/db/*`、
// EOS が読む設定のスカラー・経路・solverConfig.yaml の本文を `/meta` の属性に書き、ファイルを閉じて終了する (exit 0)。
// 比較は solver_density_cuda/tools/compare_eos_dump.py。env が無いときは何もしない (同期も出力も入れない)。

#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaConfig.cuh"

// EOS の直前 (step が FORGE_DUMP_EOS_STEP と一致する最初の呼び出しのみ)。
void eosDumpBefore(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step);
// EOS の直後。eosDumpBefore が書いた step なら全配列を書いて終了する (戻らない)。それ以外は何もしない。
void eosDumpAfterAndExit(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step);
