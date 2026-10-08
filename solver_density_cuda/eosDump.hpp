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
//
// 診断専用の再生 (env `FORGE_EOS_REPLAY_FILE=<h5>`、上の 2 つと併用、既定 off)。**状態を書き換えるので出力専用ではない**:
// この env を立てた run は EOS の比較 (§5.1 #2 の追加の試験 (3)、床の下の入力の直接比較) 以外に使わない。
// step k の EOS の直前 (境界ピンの後) に、その h5 の `/pre/<名前>` (上のダンプと同じ形式) を var.c_d の同名配列へ読み込み、
// 再ピンや別の EOS を挟まずに EOS を 1 回通して `/pre`・`/post` をダンプし終了する (ダンプの属性に replay = 1・replay_file)。
// 読み込む配列の名前の集合 (var.c_d の null でない配列と過不足なし)・長さ (nCells_all)・dtype (float32) が一致しなければ、
// 最初の時間ステップの EOS より前に終了する (exit 2)。env が無いときは何もしない。

#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaConfig.cuh"

// EOS の直前 (step が FORGE_DUMP_EOS_STEP と一致する最初の呼び出しのみ)。再生が有効なら、ここで凍結状態を読み込んでから書く。
void eosDumpBefore(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step);
// EOS の直後。eosDumpBefore が書いた step なら全配列を書いて終了する (戻らない)。それ以外は何もしない。
void eosDumpAfterAndExit(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step);
