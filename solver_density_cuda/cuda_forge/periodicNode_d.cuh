#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

#include <array>
#include <string>
#include <vector>

// node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5)。
// setPeriodicPartner + buildPeriodicNodeGroups で構築した周期ノード group (periodicRoot) に対し、
// 保存量残差 res_* を group 全員で足し合わせ、全員へ同じ和を書き戻す (gather + broadcast)。
// 合併体積 (buildPeriodicNodeGroups で vol を group 合算) と組み合わせると、両側部分 CV が
// 同 res・同 vol で更新され「1 つの CV」として bit 一致同期する。assembleResidual の末尾
// (全 flux/source 積算 + 壁/軸射影の後) に毎反復呼ぶ。cell モード / 非周期では no-op。
void periodicNodeGather_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// block-DPLUR の各 sweep 後に呼ぶ (§4.5.7)。周期 group の補正 dq を root(master) から member(slave) へ
// ミラー (slave=master) し、同一視ノードが多値化 (drift) して発散するのを防ぐ。master/slave は別 plane を
// 持つため別 dq が出るが、両者は同じ DOF なので master 解を共有させる。blockDPLUR=1 は dq_block_old_*、
// =0 は dq_*_old を対象にする。cell/非周期では no-op。
void periodicMirrorDq_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 勾配の periodic gather (§4.5 拡張)。calcGradient 後に呼び、合併体積を利用して boundary periodic node の
// Green-Gauss 勾配を「和→broadcast」で厳密合併に直す (片側勾配 → 両側)。2次再構成・粘性の精度に効く。
// 前提: calcGradient_b_d で periodic 半割面を除外しておく。非軸対称限定。cell/非周期では no-op。
void periodicGradientGather_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// RANS SST: k/ω 状態 (roK, roOmega) を周期 group root から member へミラー (§4.5)。point-implicit SST 更新の
// 直後に呼び、周期同一視ノードの k/ω drift を防ぐ。非 SST / cell / 非周期では no-op。
void periodicMirrorScalarState_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);
// 遷移モデル (γ–Re_θt) の保存量 roGamma/roReth を root→member ミラー (transition: none では no-op)。
void periodicMirrorTransitionState_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// NS 保存量 (ro,roUx,roUy,roUz,roe) を周期 group root から member へミラー (§4.5.9)。残差 gather は「同 res を
// 異なる state に足す」ため初期 desync (例: 非周期な seed 摂動) や丸めで master/slave の保存量が drift し、
// 継ぎ目隣接面が master/slave で別 state を読んでフラックス不整合 (seam 圧力欠陥) を生む。各 RK stage の保存量
// 更新直後・初期化時に呼んで slave=master を強制し DOF を真に 1 個にする。cell/非周期では no-op。
void periodicMirrorNSState_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 汎用 1 配列版 (化学種・受動種の gather/mirror 用; plans/active/species-passive-scalar-unification.md §4.1-5)。
// node 周期 DOF 同一視が有効か (cell / 非周期では false)。
bool periodicNodeActive(const solverConfig& cfg, const mesh& msh);
bool periodicSeamMergeActive(const solverConfig& cfg, const mesh& msh);
// a を周期 group で「和→broadcast」(残差・輸送対角・勾配の合併)。非有効なら no-op。
void periodicGatherArray_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , flow_float* a);
// a を root→member でミラー (状態・dq の同一視)。非有効なら no-op。
void periodicBroadcastArray_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , flow_float* a);
// a を周期 group で「max→broadcast」/「min→broadcast」(リミッタの極値 Q_max/Q_min と ψ の合併;
// plans/active/species-passive-scalar-unification.md §4.8)。全符号で正しい (periodicAtomic_d.cuh の CAS 版 atomic)。非有効なら no-op。
void periodicGatherMaxArray_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , flow_float* a);
void periodicGatherMinArray_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , flow_float* a);
// 化学種の保存量 roY{s} を root→member でミラー (化学種更新の直後に呼ぶ; §4.1-5)。
void periodicMirrorSpeciesState_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 診断 (env `FORGE_DUMP_PREGATHER=<path>`、既定 off・出力専用で数値は変えない。plan gradient-scalar-lsq-unification §5.1 #4a)。
// 周期 gather の**直前**の勾配の局所配列を、tag ごとに最初の 1 回だけ raw float [nVar][nCells][3] (成分 x,y,z が最内) で
// `<path>.<tag>` に、変数名を `<path>.<tag>.names` (1 行目 "nVar nCells"、以降 1 行 1 名) に書く。D2H は非 atomic の読み出しのみ。
// env が無いときは何もしない (同期も入れない)。既存の FORGE_DUMP_SCALARGRAD (GG の面寄与) とは別物。
// env が有効か (初回に getenv して固定)。呼び出し側はこれが false なら引数の vector も組まない。
bool preGatherDumpEnabled();
void preGatherDump(const std::string& tag, geom_int nCells, const std::vector<std::string>& names,
                   const std::vector<std::array<const flow_float*, 3>>& grads);
// main の periodicGradientGather 直前 (初期化 "init"・ループ初回 "loop1") 用: NS 6 量 (ro,Ux,Uy,Uz,P,T) × 3 成分と、
// gg 経路で同じ gather に登録される dY{s}・受動種 (ξ・凝縮モーメント) の勾配。scalarGradient: lsq では dY・受動種は
// 各 wrapper 内で合併済みなのでここには含めない (wrapper 側で gather 直前に別 tag で書く)。
void preGatherDumpMain(solverConfig& cfg , mesh& msh , variables& var , const char* tag);
