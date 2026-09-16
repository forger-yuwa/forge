#pragma once
// 周期 node group の max/min gather 用の浮動小数 atomicMax / atomicMin と gather カーネル
// (plans/active/species-passive-scalar-unification.md §4.8, codex result-2 M3)。
//
// CUDA 組込みの atomicMax/atomicMin は整数のみなので、compare-and-swap ループで浮動小数の max/min を組む。
// 値の比較は浮動小数のまま行う (順序保存整数写像 i = bits>=0 ? bits : bits^0x7fffffff と同じ全順序を、encode/decode
// パス無しで得る) ので、負値 (-2 < -1)・±0・±Inf・大きさによらず正しい。group は高々 8 member (面 2・辺 4・角 8) なので
// CAS の競合は無視できる。
// 規約: ±0 は同値扱い (max(-0,+0) はどちらが残ってもよい = 値として等しい)。NaN の member は no-op (group を汚さない;
// 明示的に skip)。root が NaN のときは root が残る (その CV の場が既に NaN なので detectNaN 側で捕まえる)。
// float / double (flowFormat.hpp の flow_float 切替) の両方を持つ。
// __global__ カーネルを持つので、ライブラリ側は periodicNode_d.cu からだけ include する (単体試験 TU は別実行体)。
#include "flowFormat.hpp"

__device__ inline void atomicMaxFloat(float* addr, float val)
{
    if (isnan(val)) return;   // NaN member は skip
    int* ai = reinterpret_cast<int*>(addr);
    int old = *ai;
    while (__int_as_float(old) < val) {
        const int assumed = old;
        old = atomicCAS(ai, assumed, __float_as_int(val));
        if (old == assumed) break;
    }
}
__device__ inline void atomicMinFloat(float* addr, float val)
{
    if (isnan(val)) return;
    int* ai = reinterpret_cast<int*>(addr);
    int old = *ai;
    while (__int_as_float(old) > val) {
        const int assumed = old;
        old = atomicCAS(ai, assumed, __float_as_int(val));
        if (old == assumed) break;
    }
}
__device__ inline void atomicMaxFloat(double* addr, double val)
{
    if (isnan(val)) return;
    unsigned long long* ai = reinterpret_cast<unsigned long long*>(addr);
    unsigned long long old = *ai;
    while (__longlong_as_double(old) < val) {
        const unsigned long long assumed = old;
        old = atomicCAS(ai, assumed, __double_as_longlong(val));
        if (old == assumed) break;
    }
}
__device__ inline void atomicMinFloat(double* addr, double val)
{
    if (isnan(val)) return;
    unsigned long long* ai = reinterpret_cast<unsigned long long*>(addr);
    unsigned long long old = *ai;
    while (__longlong_as_double(old) > val) {
        const unsigned long long assumed = old;
        old = atomicCAS(ai, assumed, __double_as_longlong(val));
        if (old == assumed) break;
    }
}

// member (root != self) の値を root へ max / min で集約する (broadcast は periodicBroadcast1FromRoot_d)。
__global__ void periodicGatherMax1ToRoot_d(geom_int nCells, geom_int* root, flow_float* a)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) { const geom_int r = root[ic]; if (r != ic) atomicMaxFloat(&a[r], a[ic]); }
}
__global__ void periodicGatherMin1ToRoot_d(geom_int nCells, geom_int* root, flow_float* a)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) { const geom_int r = root[ic]; if (r != ic) atomicMinFloat(&a[r], a[ic]); }
}
