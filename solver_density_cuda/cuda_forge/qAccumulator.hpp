// qAccumulator.hpp — FP64 影アキュムレータの host 側インタフェース
//   (plans/active/time_integration-fp64-accumulator.md §4.3)
// device 側の commit / reconcile は qAccumulator_d.cuh (__device__ inline) にある。
// こちらは **host C++ (variables.cpp) から呼べる**よう、__device__ コードを含めない。
#pragma once

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"

// 現在の FP32 の Q から FP64 正本を初期化する (restart・初期化直後に 1 回)。
void qaccInitFromQ(double* const qacc[5], flow_float* const q[5], geom_int nCells);

// 採用カウンタを 0 に戻す (monitor 区間ごと)。
void qaccResetAdoptCounter(int* adopt_d);

// 採用カウンタを読み出す (D2H・monitorInterval ごとにログへ)。
int qaccReadAdoptCounter(const int* adopt_d);
