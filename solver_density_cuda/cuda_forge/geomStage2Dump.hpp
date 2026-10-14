#pragma once

// 段 ② の診断 (plans/active/architecture-float-state-double-geometry.md §6.3 の 1・2 と closure の記録)。既定 off。
//
//   FORGE_DIAG_GEOM_STAGE2_DUMP=<h5>
//     起動の段取り (読み込み・周期の相手・構造量・LSQ の係数・ラインの構築) が済んだ時点で、読み込み時に作った量と接続を
//     <h5> に書いて終了する (時間更新・res_0 出力なし)。float のビルドと FP64 のビルドの出力を、面・incidence・境界面の
//     番号で突き合わせて比べるためのもの。各量について、本番の値 (new = 段 ②: double の写しから) に加えて、
//     段 ① までの作り方 (legacy = geom_float / flow_float の座標から) で作り直した値も書く (旧版との不一致を数えるため)。
//     legacy の作り直しは本番の配列を触らない。
//
// 書く量 (h5 のパス。添字の並びは各データセットの属性 "layout" にも書く):
//   /mesh/plane_cells [2*nPlanes]、/mesh/cell_planes_index [nCells+1]、/mesh/cell_planes [nInc]
//   /lines/{new,legacy}/{offsets,cells,prev,next}   ラインの接続 (lineImplicit の有無によらず作る)。属性 lines_built_in_run、
//                                                    lines_device_equal_new (本番のデバイスの接続と new の一致、未構築は -1)
//   /periodic/<physID>/{partnerPlnID,partnerCellID}  本番 (bint)。/periodic/<physID>/partnerPlnID_legacy は作り直し
//   /periodic/periodicRoot [nCells]                  node の周期 group の root (本番)
//   /lsq/{new,legacy}/cInt [3*nInc]                 gradLSQ=2 の係数 (new は本番の配列 lsq_coef_view)
//   /lsq/{new,legacy}/seam_group, seam_class [nInc] 継ぎ目の同値類 (継ぎ目でない incidence は -1)
//   /lsq/{new,legacy}/M6 [6*nCells]                 gradLSQ=1 の M (Mxx,Mxy,Mxz,Myy,Myz,Mzz を nCells ずつ)
//   /wall/<physID>/{iPlanes,iCells,irep,dn,dist,cos}(_legacy)  node の壁の代表内点 (境界面ごと)。y = max(dn, 1e-12)
//   /axisym/{A_closure_x,A_closure_y,A_planar}      デバイスの値 (軸対称 method 0)。診断として closure_sum64_{x,y}
//                                                    (Σ±S_device を double で足した丸める前の値)、closure_abs64_{x,y} (Σ|S_device|)、
//                                                    closure_raw64_{x,y} (生の double の面ベクトルと半径から作った closure)
//   /d1d2/<physID>/{d1,d2,jdof,jdof2,ok}(_legacy)   第一・第二内部点 (弱形式の等温壁・CHT・interfaceDiag が使う。全壁 bcond)
//   /delta_les [nCells_all]、/delta_les_legacy      SST-DES の Δmax

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaConfig.cuh"

namespace geomStage2Dump {

// <path> に書いて main の終了コードを返す。
int run(const char* path, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

}  // namespace geomStage2Dump
