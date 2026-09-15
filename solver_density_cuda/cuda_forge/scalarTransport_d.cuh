#pragma once

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

// 汎用スカラ輸送コア。特定の物理 (RANS / species など) には依存せず、
// 1 変数ぶんの保存量 ρφ について移流・拡散の残差組み立てと時間積分を提供する。
// RANS k/ω への適用は ransTransport_d.* が ScalarTransportDesc を構築して呼び出す。
struct ScalarTransportDesc {
    flow_float* phi;
    // φ の勾配 (node 壁の ghostless 弱形式拡散 ∇φ·S で使用; 内部面・cell では未使用で nullptr 可)。
    flow_float* dphidx;
    flow_float* dphidy;
    flow_float* dphidz;
    flow_float* rho_phi;
    flow_float* rho_phi_N;
    flow_float* rho_phi_M;
    flow_float* res_rho_phi;
    flow_float* res_rho_phi_m;
    flow_float* src_jac;        // 源項消散ヤコビアン対角 (∂D/∂(ρφ)≥0)。陽解法 RK の point-implicit 減衰に使用
    flow_float* transport_diag; // 輸送(移流+拡散)ヤコビアン対角 [m³/s]。advection/diffusion kernel で集計
    flow_float floor;           // realizability 下限。陰解法 point-implicit と整合
    flow_float sigma;           // 拡散係数スケール (有効粘性 = vis_lam + sigma·vis_turb)。F1 ブレンド時は k-ω 側 (F1=1) の値
    int diffusion;              // 1: 汎用拡散 (μ ベース, RANS k/ω) を使う / 0: 移流のみ
                                //   (化学種は Fick 拡散を speciesTransport 側で別途扱うため 0)
    flow_float sigma2 = static_cast<flow_float>(0.0);   // F1 ブレンド時の k-ε 側 (F1=0) σ (sstSigmaBlend=1 のときのみ使用)
    flow_float* F1 = nullptr;                             // ブレンド関数 (nullptr なら sigma 定数)
};

// 複数スカラーの移流 (+汎用拡散) を 1 回の面ループで積む (plan performance-3d-node-sst-speedup: k/ω・化学種の
// 面ループ融合)。面幾何・massflux・ρ・μ を 1 度だけ読み、各スカラーは φ の gather と atomicAdd だけ増える。
// 面ごとの演算は単一版 (scalar_advection_first_order_d / scalar_diffusion_first_order_d) と同一順序。
#define SCALAR_MULTI_MAX 4
void scalarTransportResidualMulti_d(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                                    const ScalarTransportDesc* descs, int n);

// 1 変数ぶんの移流 + (任意) 拡散残差を組み立てる。
// res_rho_phi / transport_diag は呼び出し側でゼロ初期化済みであること。同期は行わない。
void scalarTransportResidual_d(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                               const ScalarTransportDesc& desc);

// 1 変数ぶんの時間積分 (RK4 / point-implicit RK)。同期は行わない。
void scalarTimeIntegration_d(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                             const ScalarTransportDesc& desc);
