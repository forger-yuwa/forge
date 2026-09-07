#include <vector>
#include <string>
#include "periodicNode_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

// node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5)。
//
// periodicRoot[ic] は ic の所属する周期 group の root (master) CV index。root/非周期では root==ic。
// 継ぎ目で割れた 1 つの CV の両側部分が同じ group に入る。継ぎ目には双対面を作らず (周期半割面は移流
// ループから除外)、両側の内部双対面が各 partial CV の res を組む。これらを group で足し合わせ、全員へ
// 同じ和を書き戻すことで「1 CV」の残差にする。合併体積 (buildPeriodicNodeGroups) と合わせ、両者が
// 同 res・同 vol で更新され同期する。free-stream 保存は closure と継ぎ目半割面の相殺で成立 (plans §4.5.3)。

// 第1パス: slave (root != self) の res を root へ atomicAdd で集約。
__global__ void periodicGatherToRoot_d
(
    geom_int nCells, geom_int* root,
    flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const geom_int r = root[ic];
        if (r != ic) {
            atomicAdd(&res_ro[r],   res_ro[ic]);
            atomicAdd(&res_roUx[r], res_roUx[ic]);
            atomicAdd(&res_roUy[r], res_roUy[ic]);
            atomicAdd(&res_roUz[r], res_roUz[ic]);
            atomicAdd(&res_roe[r],  res_roe[ic]);
        }
    }
}

// 第2パス: root の集約済み和を slave へ broadcast (group 全員が同じ和を持つ)。
__global__ void periodicBroadcastFromRoot_d
(
    geom_int nCells, geom_int* root,
    flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const geom_int r = root[ic];
        if (r != ic) {
            res_ro[ic]   = res_ro[r];
            res_roUx[ic] = res_roUx[r];
            res_roUy[ic] = res_roUy[r];
            res_roUz[ic] = res_roUz[r];
            res_roe[ic]  = res_roe[r];
        }
    }
}

// 単一配列版 gather/broadcast (k/ω 残差・勾配など可変本数の配列用)。
__global__ void periodicGather1ToRoot_d(geom_int nCells, geom_int* root, flow_float* a)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) { const geom_int r = root[ic]; if (r != ic) atomicAdd(&a[r], a[ic]); }
}
__global__ void periodicBroadcast1FromRoot_d(geom_int nCells, geom_int* root, flow_float* a)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) { const geom_int r = root[ic]; if (r != ic) a[ic] = a[r]; }
}

void periodicNodeGather_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || msh.periodicRoot_d == nullptr || msh.nPeriodicMembers == 0) return;

    periodicGatherToRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.periodicRoot_d,
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    periodicBroadcastFromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.periodicRoot_d,
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    // RANS SST の k/ω 残差、多成分の化学種残差 res_roY{s}、凝縮モーメント残差も merged CV で合算する
    // (周期境界の DOF 同一視)。化学種を合算しないと seam で移流・反応ソースが部分体積分しか効かず
    // (面 1/2・辺 1/4・角 1/8)、エネルギー (合算済) と不整合になる (case/35 run_0049 で発覚, 2026-09-04)。
    std::vector<std::string> extra;
    if (cfg.LESorRANS == 2 && cfg.RANSmodel == 1) { extra.push_back("res_roK"); extra.push_back("res_roOmega"); }
    for (const auto& nm : var.speciesVarNames)     extra.push_back("res_" + nm);
    for (const auto& nm : var.condMomentConsNames) extra.push_back("res_" + nm);
    for (const auto& k : extra) {
        auto it = var.c_d.find(k);
        if (it == var.c_d.end() || it->second == nullptr) continue;
        periodicGather1ToRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(msh.nCells, msh.periodicRoot_d, it->second);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        periodicBroadcast1FromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(msh.nCells, msh.periodicRoot_d, it->second);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    }
}

// NS 保存量 (ro,roUx,roUy,roUz,roe) を周期 group root から member へミラー (§4.5.9)。
// 残差 gather だけでは初期 desync・丸めで master/slave の保存量が drift する (実測: Uz が seam で ~15 m/s 食い違い)。
// 継ぎ目隣接面が master/slave で別 state を読みフラックス不整合→seam 圧力欠陥を生むため、保存量更新直後に slave=master を強制する。
void periodicMirrorNSState_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || msh.periodicRoot_d == nullptr || msh.nPeriodicMembers == 0) return;
    periodicBroadcastFromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.periodicRoot_d,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    // 化学種・凝縮モーメントの保存量も root→member でミラー (残差 gather と対にして drift を防ぐ)。
    std::vector<std::string> extra;
    for (const auto& nm : var.speciesVarNames)     extra.push_back(nm);
    for (const auto& nm : var.condMomentConsNames) extra.push_back(nm);
    for (const auto& k : extra) {
        auto it = var.c_d.find(k);
        if (it == var.c_d.end() || it->second == nullptr) continue;
        periodicBroadcast1FromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(msh.nCells, msh.periodicRoot_d, it->second);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    }
}

// RANS SST: k/ω 状態 (roK, roOmega) を周期 group の root へ member からミラー (§4.5)。
// point-implicit SST (applySSTPointImplicit) は per-cell 対角で更新するため master/slave が別値になり
// drift する。更新直後に master 値を共有させて同一視を保つ (dq ミラーの k/ω 版)。SST 更新の直後に呼ぶ。
void periodicMirrorScalarState_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || msh.periodicRoot_d == nullptr || msh.nPeriodicMembers == 0) return;
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) return;
    for (const char* k : {"roK", "roOmega"}) {
        auto it = var.c_d.find(k);
        if (it == var.c_d.end() || it->second == nullptr) continue;
        periodicBroadcast1FromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(msh.nCells, msh.periodicRoot_d, it->second);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    }
}

// 勾配の periodic gather (§4.5 拡張)。Green-Gauss 勾配 ∇φ[ic]=(1/V[ic])Σφ_f·S_f は、合併体積
// (buildPeriodicNodeGroups で V[master]=V[slave]=V合併) のおかげで ∇φ_master+∇φ_slave=(1/V合併)Σ_両側 =
// 真の合併勾配になる。よって残差と同じ「和→broadcast」で boundary periodic node の片側勾配を厳密合併に直す。
// 前提: calcGradient_b_d で periodic 半割面を除外し勾配=内部のみにしておくこと (calcGradient_d.cu)。
// 非軸対称限定 (grad_volume=volume)。軸対称 (A_planar) は §4.5.8 回転で別途。cell/非周期では no-op。
void periodicGradientGather_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || msh.periodicRoot_d == nullptr || msh.nPeriodicMembers == 0) return;
    if (cfg.isAxisymmetric == 1) return; // 軸対称は grad_volume=A_planar で sum-trick が成り立たない (回転周期は §4.5.8)

    const char* keys[] = {
        "dUxdx","dUxdy","dUxdz", "dUydx","dUydy","dUydz", "dUzdx","dUzdy","dUzdz",
        "drodx","drody","drodz", "dPdx","dPdy","dPdz", "dTdx","dTdy","dTdz", "divU",
        // RANS SST 勾配 (未割当ならスキップ): k/ω の seam 拡散・生産の片側勾配を合併
        "dKdx","dKdy","dKdz", "dOmegadx","dOmegady","dOmegadz"
    };
    for (const char* k : keys) {
        auto it = var.c_d.find(k);
        if (it == var.c_d.end() || it->second == nullptr) continue;
        flow_float* a = it->second;
        periodicGather1ToRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(msh.nCells, msh.periodicRoot_d, a);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        periodicBroadcast1FromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(msh.nCells, msh.periodicRoot_d, a);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    }
}

// block-DPLUR sweep 後の dq ミラー (§4.5.7)。root の補正 dq を member へ broadcast し、
// 周期同一視ノードの dq drift (多値化→発散) を防ぐ。periodicBroadcastFromRoot_d を dq buffer に流用。
void periodicMirrorDq_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.discretization != "node" || msh.periodicRoot_d == nullptr || msh.nPeriodicMembers == 0) return;

    if (cfg.blockDPLUR == 1) {
        periodicBroadcastFromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            msh.nCells, msh.periodicRoot_d,
            var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"],
            var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"]
        );
    } else {
        periodicBroadcastFromRoot_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            msh.nCells, msh.periodicRoot_d,
            var.c_d["dq_ro_old"], var.c_d["dq_roUx_old"], var.c_d["dq_roUy_old"],
            var.c_d["dq_roUz_old"], var.c_d["dq_roe_old"]
        );
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
