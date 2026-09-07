#pragma once

#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "cuda_forge/cudaConfig.cuh"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>



class variables {
public:
    std::map<std::string, std::vector<flow_float>> c; // host cell variables
    std::map<std::string, std::vector<flow_float>> p; // host plane variables
    std::map<std::string, flow_float*> c_d; // device cell variables
    std::map<std::string, flow_float*> p_d; // device plane variables


    // 化学種輸送 (M2): registerSpecies() で 1 化学種ごとの保存量/派生量を末尾に追加するため
    // 非 const とする。registerSpecies は変数構築後 allocVariables 前に 1 度だけ呼ぶ。
    std::list<std::string> cellValNames =
    {
        "Ux"   , "Uy"    , "Uz"    , "T"      , "P" , "deltaP" ,
        "k"    , "omega" ,
        "UxN"  , "UyN"   , "UzN"   , "TN"     ,
        "roN"  , "roUxN" , "roUyN" , "roUzN"  , "roeN" , "roKN" , "roOmegaN" ,
        "roNN" , "roUxNN", "roUyNN", "roUzNN" , "roeNN" , "roKNN" , "roOmegaNN" ,
        "ro"   , "roUx"  , "roUy"  , "roUz"   , "roe"  , "roK" , "roOmega" ,
        "dTdx" , "dTdy"  , "dTdz",
        "dUxdx", "dUxdy" , "dUxdz",
        "dUydx", "dUydy" , "dUydz",
        "dUzdx", "dUzdy" , "dUzdz",
        "dPdx" , "dPdy"  , "dPdz",
        "drodx", "drody" , "drodz",
        "dKdx" , "dKdy"  , "dKdz",
        "dOmegadx", "dOmegady", "dOmegadz",

        "divU*vol" , "divU" , "divU_star",

        "ducros" , 

        "limiter_Ux" , 
        "limiter_Uy" , 
        "limiter_Uz" , 
        "limiter_ro" , 
        "limiter_P" , 
        "limiter_T" , 
        "vis_lam" ,
        "vis_turb" , "wall_dist" ,
        // node-centered: 壁ノードは壁面上 (wall_dist=0) で SST 壁 omega BC 6ν/βy² が overflow する。
        // 壁ノードの y を内部隣接ノード (壁面ノード除く) の wall_dist 平均に置換した量 (node のみ算出; cell=wall_dist)。
        "wall_y_eff" ,
        "thermCond" ,  // 熱伝導率 λ [W/(m·K)] セル毎 (M3 kinetic theory / 定数)

        "gamma" , "cp" , "Rmix" ,  // Rmix: 混合比気体定数 R [J/(kg·K)] セル毎 (M6 TP 流束キャッシュ用; CPG は (γ-1)cp/γ)

        "cfl"   , "cfl_pseudo", "dt_local",
        "res_ro"    , "res_roUx"   , "res_roUy"   , "res_roUz"   , "res_roe", "res_roK", "res_roOmega",
        "res_ro_m"  , "res_roUx_m" , "res_roUy_m" , "res_roUz_m" , "res_roe_m", "res_roK_m", "res_roOmega_m",
        "dq_ro_old" , "dq_roUx_old", "dq_roUy_old", "dq_roUz_old", "dq_roe_old",
        "dq_ro_new" , "dq_roUx_new", "dq_roUy_new", "dq_roUz_new", "dq_roe_new",
        "dq_block_old_0", "dq_block_old_1", "dq_block_old_2", "dq_block_old_3", "dq_block_old_4",
        "dq_block_new_0", "dq_block_new_1", "dq_block_new_2", "dq_block_new_3", "dq_block_new_4",
        "rhs_block_0", "rhs_block_1", "rhs_block_2", "rhs_block_3", "rhs_block_4",
        "diag_block_00", "diag_block_01", "diag_block_02", "diag_block_03", "diag_block_04",
        "diag_block_10", "diag_block_11", "diag_block_12", "diag_block_13", "diag_block_14",
        "diag_block_20", "diag_block_21", "diag_block_22", "diag_block_23", "diag_block_24",
        "diag_block_30", "diag_block_31", "diag_block_32", "diag_block_33", "diag_block_34",
        "diag_block_40", "diag_block_41", "diag_block_42", "diag_block_43", "diag_block_44",
        "sonic"   , "Ht" , "ros", // s: entropy

        "roM"  , "roUxM", "roUyM", "roUzM", "roeM" , "roKM", "roOmegaM",

        // Mesh Structure
        "volume" , "ccx" , "ccy" , "ccz",
        "A_planar",  // 軸対称: planar 面積 (= V を r で割った値)。非軸対称時は 0。勾配計算の分母。
        // axisRFloor>0 のみ: 床適用後の離散閉性 Σ_f S_f (outward)。hoop ソース/Jacobian の面積を
        // これに差し替えると任意の床でも一様圧力が厳密に力ゼロ (全面床の CV は 0 = ソース自然消滅)。
        // A_planar (勾配分母) とは役割が別なので分離する。
        "A_closure_x", "A_closure_y",
        "axisym_uy_over_r",
        "axisym_divU",
        // axisymMethod==1 (SU2 流 planar+ソース) 用: AuxVar A0=μ_tot·v/y, A1=A0·v, A2=A0·u と
        // その GG 勾配のうちソースが使う 4 成分 (∂A0/∂x, ∂A0/∂y, ∂A1/∂y, ∂A2/∂x)。
        // method 0 では未使用 (0 のまま)。
        "axisym_aux0", "axisym_aux1", "axisym_aux2",
        "dAux0dx", "dAux0dy", "dAux1dy", "dAux2dx",
        // SST 陰解法 (segregated point-implicit) 用: 消散項のヤコビアン対角
        //   src_jac_k     = ∂Dk/∂(ρk)   = β* ω
        //   src_jac_omega = ∂Dω/∂(ρω)   = 2 β ω
        "src_jac_k",
        "src_jac_omega",
        // SST automatic wall treatment (wallTreatmentSST==1) 用: wall-adjacent セルの
        // wall-function 生産 P_k = ρu_τ⁴/ν·g(1-g) (g=du⁺/dy⁺) を格納。非 wall セルは -1
        // (inactive)。ransSource が >=0 のセルで標準 P_k をこの値に置換する (methods/turbulence §6.5(d))。
        "wf_pk",
        // SST automatic wall treatment (node, wallTreatmentSST==1) 用: 第一内層ノードの **k Dirichlet 値**
        // (保存量 roK_wf = ρ·k_wf, k_wf=ω_w·μ_t,wall/ρ = SU2 SetTurbVars_WF 流)。非対象セルは -1 (inactive)。
        // >=0 のノードで k をこの値にハード Dirichlet し近壁 k 蓄積 (再付着 μ_t ピーク) を断つ。cell では常に -1。
        "roK_wf",
        // node 入口 (Dirichlet スカラー境界) ピン識別: rans_dirichlet_scalar_boundary_d が入口ノードで 1 に
        // セットし、ransSource が ==1 のセルで res_roK/res_roOmega/src_jac_* を 0 化 (保存量整合+残差除外)。
        // 既定 0 (非入口)。cell では常に 0 (node 分岐を通らない)。
        "scalarDirichletPin",
        // 診断 (k 収支調査用, node/cell 近壁 SST 比較): ransSource が確定した最終 k 生産項 P_k
        // (wall-function 置換後)。res_roK へ加える値そのもの。non-RANS では未使用。
        "Pk_diag",
        // 診断 (2026-08-11): Crocco 型 断熱回復温度モデル T_aw = T_rep + r·U_t²/(2cp), r=Pr_lam^{1/3}
        // (SU2 壁関数の熱的閉包と同式)。壁関数対象セル (cell=第一セル / node=壁ノード) に格納、
        // 非対象は -1。現状は**出力のみ** (壁状態へは未適用 — 適用は保存整合の plan 化後)。
        "Taw_diag",
        // SST 断熱壁 SU2 式熱結合 (experimental, sstThermalWallFunction==2, plan
        // turbulence-sst-su2-taw-coupling): 勾配計算後の working primitive overlay。壁関数対象の
        // 断熱壁ノードに Taw_diag と同値を格納 (非対象は -1)。viscousFlux_d が W 入射内部辺の
        // 端点温度としてのみ参照 (SU2 SetTemperature 相当)。状態 T/roe・GG/LSQ 勾配・境界 bvar
        // には一切影響しない。mode 0/1 では常に -1 (ビット不変)。
        "Taw_Prim_Overlay",
        // SST 断熱壁 defect-flux 閉包 (experimental, sstThermalWallFunction==3, theory §6.5(f) mode 3):
        // 壁ノードに H_T·m̂ (H_T=ρ_rep·cp·u_τ/T⁺ [W/m²K], m̂=内向き壁単位法線) を格納。
        // viscousFlux_d が片端のみ非零の W-I 辺で、エネルギー行の (伝導+仕事) を
        // F=(S·H⃗)(Taw−T_W) の ± ペアに置換する。非壁/他モードは 0 (inactive)。
        "Taw_HTnx", "Taw_HTny", "Taw_HTnz",
        // 診断 (2026-08-12, plan turbulence-node-wf-omega-source §4.2): W-I 双対面で **実際に
        // 運動量残差へ加えられた** 粘性 traction を壁ノードに集計する。狙った τ_w が本当に
        // 作用しているかの確認用 (twall 出力はモデル目標値へ再スケール済みで独立検証にならない)。
        //   wi_ftan     : AddTauWall 再スケール **後** の接線力の大きさ [N] (面積込み)
        //   wi_fnrm     : 同じく法線成分の **符号付き** 和 [N] (相殺しうるので単独で判断しない)
        //   wi_fnrm_abs : 同じく法線成分の **絶対値** 和 [N] (相殺なし。|Fn|/|Ft| はこちらで見る)
        //   wi_ftan_res : 再スケール **前** の解像接線力 [N]
        // 毎ステップ 0 クリアして atomicAdd で積む。壁ノード以外は 0。
        "wi_ftan", "wi_fnrm", "wi_fnrm_abs", "wi_ftan_res",
        // 診断 (2026-08-13, plan turbulence-node-wf-omega-source §4.1): omega 方程式の項別収支。
        // res_roOmega に加える前後を分解して、平衡がどの項で決まっているかを直接見る。
        // 単位はすべて [kg/(m·s²)] 相当 = res_roOmega と同じ (体積込み)。
        //   omg_prod  : 生産 P_ω·V   (現行 = α ρ S²·V)
        //   omg_dest  : 消滅 D_ω·V   (= β ρ ω²·V)
        //   omg_cross : 交差拡散 CD_ω·V
        //   omg_trans : 輸送 (対流+拡散) = ソース加算 **前** の res_roOmega
        // FORGE_OMEGA_BUDGET=1 のときだけ確保・出力される (OFF ならメモリも出力も増えない)。
        "omg_prod", "omg_dest", "omg_cross", "omg_trans",
        //   omg_axisym: 軸対称 (axisymMethod==1) の 1/y ソース·V。omg_trans はこれが加わる **前** の
        //   残差なので、軸対称ケース (case/40 等) で収支を取るには本項が要る。
        "omg_axisym",
        // E3 (plan turbulence-node-wf-omega-source §5): 壁関数の第一内層ノードに、壁法則整合の
        // ひずみ二乗 S_wf^2 = (u_tau^2·g/nu)^2 を格納 (非対象は -1 = inactive)。
        // ransSource が omega 生産にのみ使う (k 生産は既に wf_pk で置換済み)。
        // FORGE_WF_OMEGA_SOURCE=1 のときだけ ransSource が参照する (既定はビット不変)。
        "wf_sprod",
        "wf_g",          // 診断: 壁関数が採用した g=du+/dy+ (FORGE_WF_CLOSURE_DIAG)
        // 診断 (2026-08-13, plan turbulence-node-wf-representative-point §3.1): node 壁関数の
        // 代表点 (Normal_Neighbor) 選択の幾何。**壁ノード**に格納 (非壁は -1)。
        //   rep_id      : 選ばれた irep のインデックス (float 化。取り違え検出用)
        //   rep_y       : 壁関数が使う法線距離 y=(x_I-x_W)·(-n)
        //   rep_dist    : 実距離 |x_I-x_W|
        //   rep_cos     : best_cos (内向き法線との cos。現行の唯一の選択基準)
        //   rep_toff    : 接線方向オフセット |(x_I-x_W) - y·(-n)|
        //   rep_wdratio : wall_dist[irep]/y (壁距離場との整合)
        // FORGE_WF_REP_DIAG=1 のときだけ確保・出力される。
        //   rep_nx/ny/nz: 壁の **外向き** 単位法線 (後処理で法線上補間をするのに要る。
        //                  代表点方向で代用すると曲面で 15-20 度ずれる)
        "rep_id", "rep_y", "rep_dist", "rep_cos", "rep_toff", "rep_wdratio",
        "rep_nx", "rep_ny", "rep_nz",
        // 診断/分離実験 (2026-08-12, plan §3 の 2×2 factorial): node 壁関数の **第一内層ノード
        // (irep)** を 1 でマークする (壁ノード・非対象は 0)。`wf_pk>=0` は壁ノードにも立つため
        // 「第一内層だけ」を選ぶマスクとしては使えない。cell では常に 0 (= cell ビット不変が構造的に保証)。
        "wf_irep_flag",
        // SST automatic wall treatment (node, wallTreatmentSST==1) / WMLES (node) 用: 壁ノードの
        // モデル壁せん断応力 τ_w=ρu_τ² [Pa] を格納 (非壁ノードは -1 = inactive)。viscousFlux_d が片端のみ
        // 壁ノードの内部双対面 (W-I エッジ) の接線粘性応力を τ_w に再スケールする (SU2 AddTauWall 相当)。
        // cell モードでは常に -1 (cell は viscousFlux_wall_d の壁面 modeled τ_w を内部セルに課すため不要)。
        "Tau_Wall",
        // node k/ω Dirichlet (SU2 SetTurbVars_WF 流): 第一内層ノードの ω 固定値 ρω_w (Menter ブレンド)。
        // 非対象ノードは -1。roK_wf と対で explicit/implicit 両経路の状態ピンに使う (2026-07-20:
        // ω を壁ノードのみピンすると第一内点がせん断駆動 P_ω で暴騰し νt が立たない — case/38)。
        "roOmega_wf",
        // WMLES 等温壁 (node) 用: 壁ノードのモデル壁熱流束 q_w [W/m²] (壁→流体正, 非対象は -1)。
        // viscousFlux_d の AddQWall (W-I 熱流束置換) と zeroWallMomentumResidual の res_roe ピン
        // マーカを兼ねる (methods/turbulence §10.4)。断熱 WMLES / cell では常に -1。
        "Qw_Wall",
        // SST 陰解法用: k/ω 輸送項（移流+拡散）のヤコビアン対角（point-implicit）[m³/s]。
        // 1次風上移流の Σ_f max(±ṁ,0)/ρ と拡散 Σ_f (μ_face/ρ)(|δ|/dcc) を面ループで集計し、
        // D_φ = V/Δτ + V·src_jac + transport_diag に加える。壁近傍の細セルで陽的輸送 stiff 性を
        // 陰化し安定 cfl_pseudo を一桁以上引き上げる。defect-correction なので定常解は不変。
        "transport_diag_k",
        "transport_diag_omega",
        // SST-DES (DESmode>0) 用 (methods/turbulence §8)。delta_les は幾何セットアップ時に 1 回計算する
        // 静的量、他 3 つは毎 step の診断量。float32 で f_d が 0/1 へ張り付くのは NaN より危険なため
        // r_d / f_d は必ず出力する (plan §4.2)。
        "delta_les",   // per-cell グリッドスケール Δmax (隣接重心距離の最大)
        "l_des",       // per-step l_DDES / l_IDDES (DESmode 1/2)
        "fd_shield",   // per-step シールド関数 (DESmode 1: f_d [0=RANS,1=LES], 2: f̃_d [1=RANS,0=LES] ★向き逆)
        "rd_des",      // per-step r_d (DESmode 1: (νt+ν) 版, 2: r_dt=νt 版)。clamp[0,10] 上限張り付き=飽和診断
        "fe_iddes"     // per-step IDDES log-layer mismatch 補正 f_e (DESmode==2 のみ非零。theory §8.6)
    };

    const std::list<std::string> planeValNames = 
    {
        "cfl_pln" , "cfl_pseudo_pln",
        //"US" , "T" ,  
        //"USN", "Ux" , "Uy" , "Uz", "ro" , "P",

        // Mesh Structure
        "sx"  , "sy" , "sz" , "ss" ,
        "sx_planar" , "sy_planar" , "sz_planar" , "ss_planar" ,
        "pcx" , "pcy" , "pcz" ,
        "fx"  , 
        "massflux" ,
        "dcc"   // dcc: distance between two cell centers
    };

    // 化学種 (M2): registerSpecies() が Y{s} を出力対象に追加するため非 const。
    std::list<std::string> output_cellValNames =
    {
        "ro"    , "Ux"    , "Uy"    , "Uz"  , "T" , "P" ,
        "roUx"  , "roUy"  , "roUz"  , "roe" , "roK" , "roOmega" , "k" , "omega" ,
        "cfl"   , "volume", "sonic" , 
        "dUxdx" , "dUxdy" , "dUxdz" , 
        "dUydx" , "dUydy" , "dUydz" , 
        "dUzdx" , "dUzdy" , "dUzdz" ,
        "drodx" , "drody" , "drodz" ,
        "dPdx"  , "dPdy"  , "dPdz" ,
        "dKdx"  , "dKdy"  , "dKdz" ,
        "dOmegadx", "dOmegady", "dOmegadz",

        "ducros",

        "limiter_ro" , 
        "limiter_Ux" , 
        "limiter_Uy" , 
        "limiter_Uz" , 
        "limiter_P" , 
        //"limiter_Ht" , 

        "dt_local",

        "axisym_divU",

        "vis_lam" , "vis_turb" , "wall_dist" , "thermCond" ,

        // k 収支調査 (node/cell 近壁 SST 比較, 一時診断): 壁関数生産 wf_pk(>=0=壁関数処理/-1=未処理),
        // k Dirichlet 値 roK_wf, 最終 k 生産 Pk_diag, 消散ヤコビ src_jac_k, k 輸送対角 transport_diag_k
        "wf_pk" , "roK_wf" , "Pk_diag" , "Taw_diag" , "src_jac_k" , "transport_diag_k" ,

        // omega 残差収支調査 (入口×壁コーナー残差プラトーの局在, 一時診断): 残差 res_roOmega,
        // 源項ヤコビ src_jac_omega, 輸送対角 transport_diag_omega
        "res_roOmega" , "src_jac_omega" , "transport_diag_omega" , "res_roK" ,

        // W-I 実力診断 (plan turbulence-node-wf-omega-source §4.2)
        "wi_ftan" , "wi_fnrm" , "wi_fnrm_abs" , "wi_ftan_res" , "wf_irep_flag" ,
        "omg_prod" , "omg_dest" , "omg_cross" , "omg_trans" , "omg_axisym" , "wf_sprod" , "wf_g" ,
        "rep_id" , "rep_y" , "rep_dist" , "rep_cos" , "rep_toff" , "rep_wdratio" ,
        "rep_nx" , "rep_ny" , "rep_nz" ,

        // SST-DES 診断 (DESmode>0 で意味を持つ。DESmode==0 でも 0 埋めで出力されるだけ・無害)
        "delta_les" , "l_des" , "fd_shield" , "rd_des" , "fe_iddes"
    };

    //not yet implemented
    const std::list<std::string> read_cellValNames =
    {
        "ro"  , "roUx"  , "roUy"  , "roUz"  , "roe" , "wall_dist" , "roK" , "roOmega"
    };

    // 化学種輸送 (M2)。registerSpecies() で実際に登録された化学種数 (>=2 で機構が起動)。
    // 単成分 (0 または 1) のときは化学種変数を一切登録せず M1 と同一経路を保つ。
    int nSpeciesRegistered = 0;

    // 化学種名 (config 由来, index 順)。s 番目の保存質量分率は "roY{s}"。
    std::vector<std::string> speciesVarNames; // ["roY0","roY1",...]

    // 非平衡凝縮 (4 モーメント方程式)。registerCondensation() で登録された凝縮種数。
    // 0 のときは凝縮変数を一切登録せず従来経路を保つ。methods/condensation/ 参照。
    int nCondSpeciesRegistered = 0;
    // 凝縮モーメントの保存量名 (登録順)。1 凝縮種あたり 4 本: "rog_{s}","roQ2_{s}","roQ1_{s}","roQ0_{s}"。
    // 原始量は ro を外した名前 ("g_{s}" 等)、N/M/res/src_jac/transport_diag も命名規約で導出する
    // (condensationTransport_d.cu / update_d.cu / main.cpp が参照)。
    std::vector<std::string> condMomentConsNames;

    variables();

    //variables(const int& ,  mesh&);
    ~variables();

    // 化学種 (M2): cellValNames / output_cellValNames へ 1 化学種ごとの変数を追加し、
    // c / c_d マップにも空エントリを作る。allocVariables より前に 1 度だけ呼ぶこと。
    // nSpecies <= 1 のときは何もしない (単成分は M1 経路)。
    void registerSpecies(int nSpecies, int chemistry = 0);

    // 非平衡凝縮 (Phase 1): 凝縮種ごとに 4 モーメント (ρg,ρQ2,ρQ1,ρQ0) の保存量・原始量・RK ステージ・
    // 残差・point-implicit 対角を cellValNames / c / c_d へ追加する。allocVariables より前に
    // 1 度だけ呼ぶ。nCondSpecies <= 0 のときは何もしない (従来経路)。
    void registerCondensation(int nCondSpecies);

    void allocVariables(const int &useGPU , mesh& msh);


    void copyVariables_cell_plane_H2D_all();
    void copyVariables_cell_H2D (std::list<std::string>);
    void copyVariables_plane_H2D(std::list<std::string>);

    void copyVariables_cell_plane_D2H_all();
    void copyVariables_cell_D2H (std::list<std::string>);
    void copyVariables_plane_D2H(std::list<std::string>);

    void setStructuralVariables(solverConfig& cfg, cudaConfig& cuda_cfg , mesh& msh);
    void setStructuralVariables_d(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh);

    void readValueHDF5(std::string fname , mesh& msh,
                       flow_float kInit = 0.0, flow_float omegaInit = 0.0);
};