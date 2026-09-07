#pragma once

#include "flowFormat.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "yaml-cpp/yaml.h"
#include <stdexcept>


class solverConfig
{
private:
    std::string solConfigFileName;

public:
    std::string meshFormat; 
    std::string meshFileName; 
    std::string valueFileName; 

    int gpu;

    std::string solver;

    //Time
    int endTimeControl; // 0: use dt , 1: cfl
    int nStepOuter;
    int outStepInterval;
    int outStepStart;

    int dtControl; // 0: use dt , 1: cfl
    flow_float totalTime=0.0;
    flow_float dt;
    flow_float dt_pseudo;
    flow_float cfl;
    flow_float cfl_pseudo;
    flow_float implicitRelax = 1.0;
    flow_float implicitRelaxSST = -1.0; // -1: implicitRelax に倒置 (既定動作不変)
    // 陰的更新の正値性ガード (commit 時の局所 under-relax)。0=OFF (既定・ビット同一迂回)。
    // >0 で「1 step で ro・内部エネルギーが alpha 倍未満に落ちる」セルの Δq を半減列で縮小。
    // plans/active/time_integration-update-positivity-guard.md
    flow_float updateGuardAlpha = 0.0;
    // line-implicit (壁法線ライン block-Thomas を DPLUR に埋め込む)。0=OFF (既定)。
    // 1 で高 AR 積層方向の結合を直接解に昇格し cfl_pseudo 上限を引き上げる。
    // blockDPLUR==1 専用・lowMachPrecond>=2 とは併用不可 (config 検証で拒否)。
    // plans/active/time_integration-line-implicit.md
    int lineImplicit = 0;
    // line-implicit v2 試作 (plans/active/time_integration-line-implicit-viscous-v2.md)。lineImplicit==1 専用。
    int lineKFreeze = 0;              // 1: dual-time サブ反復間で K/diag/LU を凍結 (subiter 0 のみ構築)
    int lineViscCoupling = 0;         // 1: line 面にスカラー粘性結合 K+=α·I (対角 2α→α)
    flow_float lineViscousDtRelief = 0.0;  // θ: on-line セルの擬似 dt 粘性項を (1−θ) 倍
    int lineDtDirectional = 0;        // 1: 方向別 dt — line 面の λ (音響込み) を CFL の max から除外
    int lineDtWallRelief = 0;         // 1: (診断) wall 種境界半割面の λ も on-line セルの CFL max から除外
    // 軸対称 near-axis 安定化: 擬似時間スペクトル半径に軸項 λ_axis=β·(|u_r|+c)·A_planar を加える。
    // 近軸 (r→0) で Δτ∝CFL·r/(|u_r|+c) を自然に与え半径運動量不安定を抑える。0=不変 (既定)。
    flow_float axisTimestepBeta = 0.0;
    int blockDPLUR = 0;
    // 多成分 face 整合再構成: 0 (既定・ビット不変, 組成は owner セル 1 次=mixed-order)、
    // 1: Y_s を ρ と同じ勾配+limiter で face へ 2 次再構成し thermo/species 流束で同一 face 組成を使う。
    int speciesFaceReconstruction = 0;
    // multispeciesRhoYCommonLimiter: opt-in 診断。0 (既定・ビット不変)、
    // 1: ρ と全 species に共通リミタ ψ_ρY=min(ψ_ρ, min_s ψ_Y_s) を適用し ρ_f=ρ(Y_f) 整合だけを切り分ける
    //    (p・速度は各自のリミタのまま)。nSpecies>1 かつ speciesFaceReconstruction>=1 で有効。
    int multispeciesRhoYCommonLimiter = 0;
    int speciesImplicitCoupling = 0; // 多成分 TP 陰解法 (timeIntegration==11, nSpecies>=2) の化学種更新方式。
                                     // 0: 従来 segregated 点陰的 forward-Euler (既定・ビット不変)。
                                     // 1: 緩和整合 scalar-DPLUR (流れ block と同一 dt_local/implicitRelax/nStepInner
                                     //    sweep で ρY_s を緩和し、要因2 の擬似時間緩和ミスマッチを解消)。
                                     // 詳細: plans/accepted/thermophysics-species-implicit-coupling.md。
    int implicitSolvePrecision = 0; // block-DPLUR 線形 solve の内部精度。0: float (既定・高速), 1: double。
                                    // 残差/状態は float のまま、Jacobian 構築+5×5 solve のみ double 化する混合精度
                                    // (iterative refinement)。軸対称 近軸の float 陰解固着 (Uy が −15 でなく
                                    // −0.6 固着) を根治するが double は遅い (RTX で ~×2.6)。詳細:
                                    // plans/archived/precision-mixed-axisym.md。blockDPLUR==1・lowMachPrecond 0/1 でのみ有効。
    // Ducros センサによる衝撃近傍の MUSCL リミタ強制 1 次化。0: off (既定・1 次化を使わない;
    // 衝撃近傍も MUSCL 2 次のまま)、1: on (衝撃近傍でリミタを強制 1 次化する従来挙動)。
    // off のときは ducros センサ配列を 0 に潰すため duc≤0.8 で apply_ducros_limiter がリミタを素通しする。
    // 検証 (case36/SLAU): on/off で場の差<0.1% (SLAU 自身の散逸が支配し 1 次化はほぼ冗長)。
    int ducrosLimiter = 0;
    int detectNaN = 0;             // 0: off (既定), 1: 保存量+P の非有限値を検査し、
                                   // 見つけたら res_nan_<step>.h5 をダンプして即停止する診断モード。
                                   // off のときは検査を一切行わないため通常実行はビット不変。
                                   // 検査は fused 1 カーネルで device フラグへ集約し、detectNaNInterval
                                   // ステップごとにのみ host へ読み出す (per-step 同期を避ける)。
    int detectNaNInterval = 1;     // detectNaN==1 時にフラグを host 読み出しする間隔 [step] (既定 1=毎ステップ)。
                                   // 大きくすると per-step 同期が減るが NaN 検知が最大 (interval-1) step 遅れる。
    // モニタリング出力 (= 残差 CSV の device→host flush と、`max cfl`/`dt` の console 出力) を行う共通間隔 [step]。
    // 既定 1=毎ステップ (従来挙動)。>1 で per-step の host 同期 (残差 reduction の D2H・max cfl の thrust::max_element)
    // をまとめて間引く。残差は毎ステップ device バッファに記録され、この間隔ごとに一括 D2H+CSV 書き出しされる
    // (毎ステップの値・行構成は不変)。`max cfl`/`dt` は console にこの間隔ごとに表示される。
    // 注意: dt 適応 (dtControl==1) はこの間隔とは独立 (下記)。詳細:
    // plans/accepted/architecture-residual-monitor-async.md, architecture-perphase-profiling-hotspot.md。
    int monitorInterval = 1;
    int lowMachPrecond = 0;        // Weiss-Smith 低マッハ前処理。0: off (従来・ビット不変),
                                   // 1: フラックス散逸前処理のみ (RHS, c_hat→c'。収束解を低マッハ域で変更),
                                   // 2: RHS(c' 散逸)+LHS 完全 Γ⁻¹A 前処理 (block-DPLUR 擬似時間項+setDT 拡大),
                                   // 3: LHS 完全前処理のみ (RHS 散逸は c_hat のまま=lowMachPrecond 0 と
                                   //    収束解ビット一致・保存性不変。純粋に収束加速だけを狙う LHS 操作)。
                                   // 2/3 は timeIntegration==11 && blockDPLUR==1 必須。
    flow_float precondEps = 0.15;  // 低マッハ前処理の停留点フロア ε (Ur=min(c,max(|u|,ε·c)))。
                                   // ε 小ほど低マッハ振動を強く減衰するが ε≲0.1 は発散 (ε=0.05 で NaN)。
                                   // ε=0.15: M4 ノズルで limit-cycle 振幅 −32% (検証済), ε=0.3: −17%。
    int lowMachThornber = 0;       // 0: off (従来), 1: Thornber 型再構成補正 (SLAU の L/R 速度ジャンプを
                                   // z=min(M,1) で縮約)。lowMachPrecond と直交・併用可。SLAU 経路のみ。
    // (旧 keepDissipation は廃止。以下 keepDissType は別設計: KEEP 中心流束は不変のまま独立な散逸レイヤを加算する)
    int keepDissType = 0;          // KEEP 用 opt-in 散逸レイヤ (plans/accepted/convection-keep-es-dissipation.md)。
                                   // 0: off (純粋 KEEP・散逸なし・ビット不変, 既定)
                                   // 1: scalar ES 散逸 F -= 0.5*σ*λ'*ΔU (全成分同一 λ' → 渦も食う)
                                   // 2: matrix ES 散逸 F -= 0.5*σ*R|Λ'|S Rᵀ Δw (音響のみ c' 込み → 渦を守る)
                                   //    (λ'=|Un|+c', lowMachPrecond>=1 なら c'=lowMachCprime, else c。
                                   //     LLF 型は Δw·ΔU>=0 で ES。KE 非増加は非保証 → σ は TGV L2 で較正)
    flow_float keepDissCoeff = 0.05; // 上記 σ。0 で実質 off。既定 0.05 (L1: 市松~4桁減衰/400step, L2: TGV KE cost 2.7%)。過大は解像乱流を殺す (σ=0.15 で 8.4%)。
    int keepDissCprime = 1;        // KEEP 散逸レイヤの音響波速に c'=lowMachCprime を使うか (lowMachPrecond から独立)。
    int keepDissPrecond = 0;       // 1: matrix CPG 枝の音響対散逸を Turkel 前処理 2×2 に置換 (既定 0=従来・ビット不変)。
                                   //    Δp 散逸 ∝ c²/Ur 増強 (低マッハ市松キラー)・ΔUn 散逸 ∝ Ur 縮小 (渦保護)、
                                   //    M>=1 で標準 Roe |A| に復帰。plans/active/convection-keep-diss-lowmach-precond.md
    // f_d 駆動 σ ブレンド (turbulence-iddes-sst §4.8, 既定 0=ビット不変)。1 で ES 散逸の σ を
    // face ごと σ_f=max(keepDissCoeff, (1-f̃_d)·keepDissCoeffMax) に置換 (DESmode>0 必須)。
    int keepDissFdBlend = 0;
    flow_float keepDissCoeffMax = 1.0;
    // 高周波圧力欠陥駆動 mass-flux 補正 (Rhie–Chow 型市松キラー, matrix CPG 枝のみ)。
    // δṁ=−½C_cb·S·δp^HF/Ur を [1,ū,H̄t] で全保存量へ整合配布 (sign gate で面ごと ES 保証)。
    // keepDissJump>=1 必須。plans/active/convection-keep-cb-pressure-correction.md
    flow_float keepDissCbCoeff = 0.0;   // C_cb (0=off・ビット不変)。σ_min×precond の Δp 項への加算
    flow_float keepDissCbEps   = 0.10;  // ε_cb: Ur=min(c,max(|u|,ε·c)) のカットオフ Mach
    // opBlend (fdblend×precond) の標準 Roe 側に生ジャンプを使う (RANS/grey 帯=完全風上相当)。
    // jump2 の再構成フィルタは 2-4Δ を保護するため、σ→1 でも 2-3Δ 定在モードへの実効ジャンプが
    // 中央値 6% に落ちて散逸ゾーニングが無効化される (case/39 丘頂鋸歯)。生ジャンプは一様場で
    // 厳密 0 → free-stream 安全。keepDissJump==2 かつ opBlend 有効時のみ意味を持つ。
    int keepDissOpBlendRaw = 0;
    int keepDissJump = 0;          // 散逸ジャンプ: 0=生 (既定・ビット不変) / 1=再構成後ジャンプ /
                                   // 2=再構成+sign-property クリップ (TeCNO 型・証明付き ES, 推奨)。
                                   // matrix CPG 枝のみ。市松 (2Δ) は中心勾配ゼロで再構成が効かずフル減衰のまま、
                                   // 解像スケールのドレインを桁減 (plans/accepted/convection-keep-diss-recon-jump.md)。
                                   // c' の実体は「σ の面ごとマッハ自動調整」。単一マッハ領域なら 0 + σ~0.015 で
                                   // 同等 Pareto 点 (市松 3.9e-8/KE 1.10% ≈ c'+σ0.05 の 7.1e-8/1.36%) に届く。
                                   // 1 (既定): マッハ混在流 (チャンバー低M+プルーム超音速) でグローバル σ が
                                   //   両立できないケース向け。0: フル c — 単一領域で σ を手動較正する運用。
    flow_float dt_max;
    flow_float dt_min;

    int unsteady; // steady , unsteady
    int dualTime; // 0: off , 1: on
    int timeIntegration; // 1: Euler explicit, 3: 3rd Runge explicit, 4: 4th RK, 11: implicit defect-correction

    int isImplicit; // 0: exp , 1:imp;

    // for inner loop
    int nStage;
    int nStepInner;
    int nSubIterDualTime = 20; // dual-time: 物理ステップあたりの擬似時間サブ反復数
    int bdfOrder = 2;          // dual-time: 物理時間 BDF 次数 (1 or 2、初回ステップは BDF1)
    flow_float unsteadyDiagCoef = 0.0; // dual-time: 陰解法対角へ加える物理時間項係数 a/Δt（定常は 0）。driver が毎ステップ設定
    std::vector<flow_float> coef_N;
    std::vector<flow_float> coef_M; 
    std::vector<flow_float> coef_Res;

    std::vector<flow_float> coef_DT_4thRunge;
    std::vector<flow_float> coef_Res_4thRunge;

    flow_float  convMethod;
    int limiter;    // 0: off, 1: Barth-Jespersen, 2: Venkata, -1: legacy

    // free-stream 保存: 対流流束の圧力項を (p_tilde - pRef)*s で組み、非直交メッシュで
    // 大きな p*s を float32 加算する際の桁落ち(metric closure 由来の偽運動量源)を抑える。
    // 既定 0.0 で従来挙動(ビット不変)。一様基準(動作/フリーストリーム)静圧を入れる。
    flow_float  pRef = 0.0;

    // U∞≠0 の動く一様流保存: KEEP (CPG) の移流項を参照一様流 (roRef, uRef) の流束 F∞(s) との
    // 差分形因数分解で組む (plans/active/convection-freestream-preserving-flux.md §8)。
    // F∞ は s 線形の定数流束なので保存厳密・セル和は解析ゼロのゲージ。既定 roRef=0.0 = off (ビット不変)。
    flow_float  roRef = 0.0;
    flow_float  uRefX = 0.0;
    flow_float  uRefY = 0.0;
    flow_float  uRefZ = 0.0;

    int LESorRANS; // 0:no 1:LES 2:RANS
    int LESmodel; // 1:WALE
    int RANSmodel = 0; // 0:none 1:SST
    int scalarDiffusion = 1; // 0:advection-only 1:advection+diffusion
    int dilatationCorrection = 2; // SST生産項の圧縮性補正 0:off 1:deviatoric(A) 2:deviatoric+isotropic(A+B) 既定:2
    // SST ω 交差拡散の point-implicit Jacobian (plans/active/turbulence-sst-omega-crossdiff-jacobian.md)。
    // 1 で CD_ω>0 のとき対角へ +CD_ω/(ρω) を加える (収束解不変・サブ反復収束の改善)。既定 0 = ビット不変。
    int sstCrossDiffJac = 0;
    int katoLaunder = 0; // SST生産項 Kato-Launder 補正 0:標準 mu_t S^2 1:mu_t S Omega 既定:0 (methods/turbulence §7.5)
    int wallTreatmentSST = 1; // SST壁処理 0:low-Re壁解像(60ν/β₁y²) 1:automatic(y⁺非依存) 既定:1 (methods/turbulence §6.5)
    int sstThermalWallFunction = 0; // SST壁関数の熱的閉包 (§6.5(f))。0=OFF / 1=output-only (Tsb=Taw_diag, 場非介入, 生産baseline) / 2=experimental SU2 coupled (Taw primitive overlay, 未採用) / 3=experimental defect-flux (保存的壁層エネルギー流束 H_T(Taw−T_W), T[W]→Taw を残差の解として実現)。wallTreatmentSST==1時のみ有効
    int sstEnergyWallFunction = 0;  // SST壁関数のエネルギー流束置換: 等温壁のKader q_w (§6.5(g))。wallTreatmentSST==1×wall_isothermalのみ有効。既定0=OFF
                              // 注: 既定を 0→1 に変更 (2026-06-28, user 指示)。cell 含む全ケースが automatic 壁関数に
                              //     なるため、低Re 前提の検証 (case/26 Cf, flat-plate 回帰 等) は要再検証。
    // node-centered の第一内層ノードで k をハード Dirichlet (k=ω_w·μ_t,wall/ρ, SU2 SetTurbVars_WF 流) にする。
    // 既定 1=ON (near-wall k 蓄積=再付着 μ_t ピークを除去)。0 で k を解く prod-fix (x_R がやや短め=cell/SU2 寄り)。
    // 副作用: 非平衡再付着で u_τ→0 → k_wf→0 と乱流を抑え x_R が伸びる (case/18.backstep: 7.63→8.67)。
    // wallTreatmentSST==1 かつ node のときだけ有効。cell では無視。
    int nodeKwfDirichlet = 1;
    // node ω Dirichlet (第一内層ノード, SU2 SetTurbVars_WF 対応)。壁ノード ω アンカは維持したまま
    // 第一内層にも ρω_w をピンする (νt(第一内層)=ν_t,wall を構造保証)。凹コーナーで複数壁が同一
    // 第一内層ノードを共有すると値 race の懸念 (methods/turbulence §3.7) があるため既定 OFF。
    // 周期チャネル等の角なし WMLES/EWT ケースで opt-in する (case/38)。
    int nodeOmegaWfDirichlet = 0;
    // SST 初期乱流 (IC に roK/roOmega が無いときの初期値)。既定 0 = 従来動作 (k=ω=0, ビット不変)。
    // ただし ω=0 は mu_t=k/ω が ill-posed で cold start から発散しやすい (特に node wt=1)。
    // 非 SST IC (層流場/メッシュ変換直後) から SST を始める場合は freestream 値 (例 inlet と同じ
    // k, ω) を設定して 0 初期化を避ける。読み込みは variables::readValueHDF5。
    flow_float kInit = 0.0;     // 初期 k [m²/s²] (roK 欠落時 roK=ro*kInit)
    flow_float omegaInit = 0.0; // 初期 ω [1/s]   (roOmega 欠落時 roOmega=ro*omegaInit)
    flow_float turbulentPrandtl = 0.85; // 乱流プラントル数 Pr_t。乱流熱伝導 k_t=cp*mu_t/Pr_t に使用 (既定 0.85)
    // 層流熱伝導のモデル: 0=一定 (cfg.thermCond, 従来・既定) / 1=constant-Pr (k=μ(T)·cp/prandtlLam)。
    // Sutherland 粘性 + 一定 k は高温で分子 Pr が 1.5-1.9 に漂い非物理 (2026-08-11 レビュー指摘)。
    int thermCondMethod = 0;
    flow_float prandtlLam = 0.72;   // 層流プラントル数 (thermCondMethod=1 と壁関数 T_aw 回復係数 r=Pr^{1/3} に使用)

    // WMLES 代数壁応力モデル (Reichardt+Kader, methods/turbulence §10)。有効化は bcondConfig の
    // 壁単位 ints: wallModelLES=1 (LESorRANS!=2 のみ)。以下はモデルパラメータ (turbulence セクション)。
    flow_float wmlesNewtonTol = 1.0e-6; // u_τ Newton の相対許容誤差
    int        wmlesNewtonMaxIt = 20;   // u_τ Newton の最大反復 (warm start 時は 1-3 回で収束想定)
    flow_float wmlesPrt = 0.9;          // Kader 温度壁法則の乱流プラントル数 (SGS 側 turbulentPrandtl とは独立)

    // 一様体積力 [N/m³] (周期チャネル駆動等, wmles plan §5-7 / methods/time_integration)。
    // 運動量に f_i·V、エネルギーに (f·u)·V を residual へ加算する。既定 0 = off (ビット不変)。
    // 軸対称は未対応 (config 読込時に拒否)。陰解法では明示ソース扱い (Jacobian 無し。f·u の
    // 対角寄与は微小のため省略、必要になれば追加)。
    flow_float bodyForceX = 0.0;
    flow_float bodyForceY = 0.0;
    flow_float bodyForceZ = 0.0;

    // 質量流量一定制御 (bodyForceCtrl=1): 体積平均 ρu_x を bodyForceCtrlTarget に保つよう
    // 毎物理ステップ bodyForceX を deadbeat 更新する (Benocci & Pinelli 1990 の圧縮性版,
    // methods/time_integration/implementation.md §一様体積力)。unsteady=1 必須。x 成分のみ。
    int        bodyForceCtrl = 0;
    flow_float bodyForceCtrlTarget = 0.0; // 目標 <ρu_x>_V [kg/m²/s]
    flow_float bodyForceCtrlRelax = 1.0;  // 更新則の緩和係数 γ

    // SST-DES (Detached-Eddy Simulation) length scale 修正。SST k 消滅項を
    //   D_k = β* ρ k ω = ρ k^{3/2}/l_RANS  →  ρ k^{3/2}/l_DDES
    // へ切り替え、剥離域のみ LES、付着 BL は RANS のままにする (methods/turbulence §8)。
    //   0: RANS (既定・既存 SST とビット不変)
    //   1: DDES (Spalart 2006 シールド f_d で BL を保護)
    //   2: IDDES (Phase 2・未実装、現状は範囲のみ受理)
    int DESmode = 0;
    // DES グリッドスケール係数 C_DES (l_LES = C_DES·Δmax)。SST は k-ω/k-ε ブレンドで
    //   C_DES = F1·C_DES_kw + (1-F1)·C_DES_ke (Strelets 2001)。
    flow_float C_DES_kw = 0.78;
    flow_float C_DES_ke = 0.61;

    int isCompressible;
    int isAxisymmetric = 0;

    // 軸対称の定式化: 0 = r 重み幾何 (B 流儀, 従来・既定・ビット不変) / 1 = SU2 流
    // (planar 幾何 + 1/y ソース項。軸ノードは通常 DOF として解く。plan
    // axisymmetric-su2-source-formulation)。
    int axisymMethod = 0;

    // 軸対称 r 床 (axisymMethod==0, SU2 の y<EPS ガードの r 重み版): r < axisRFloor の面・セルは
    // r 重みを床値で打ち切り (幾何が消えない)、その帯のセルは hoop ソースも Jacobian も
    // 課さない (uy_over_r も 0)。node 軸行の体積消失→真空化への対症。0 (既定) = 従来
    // (床 1e-20 = 実質ゼロ面積)。値はメッシュ依存: 軸行重心 < axisRFloor < 第一内点行重心。
    flow_float axisRFloor = 0.0;

    // 軸対称 hoop ソース面積を「解析 A_planar」でなく「離散閉性 Σ_f r_f S_f (床適用後の面ベクトル和)」に
    // する (0:既定 従来・ビット不変 / 1:離散閉性)。一様圧の半径運動量は
    //   res_roUy = −Σ_f p r_f S_f,y + p·A
    // なので A=Σ_f r_f S_f,y なら**丸め誤差に依らず厳密に打ち消す** (自由流保持 / GCL 相当)。
    // float32 メトリック (geom_float=float) では高 AR 壁 CV で Σ r_f S_f,y と A_planar が最大数十 % ずれ
    // (case/43: 生産 NS y+1.5 で 59%、実害は壁 P で ~3 Pa = 4e-5)、その偽半径力を根本から消す。
    // axisRFloor>0 のときは従来どおり (床の閉性補正が同じ配列を使うため) 常に有効。
    int hoopAreaFromClosure = 0;

    // 離散化レイアウト。"cell": cell-centered FVM (既定・従来)、"node": node-centered
    // (中点双対 median-dual) FVM。methods/discretization/ 参照。
    std::string discretization = "cell";

    // 境界隣接 CV の 2 次再構成を 1 次に落とす (0:off 既定, 1:on)。node-centered の壁近傍
    // 高マッハ発散の原因切り分け診断 / ロバス化。methods/discretization.md §7。
    int bndFirstOrder = 0;

    // node-centered 軸対称: 軸上 CV (R=0) の cell 中心に CV 面積加重重心を使うか (1:既定, ゼロ回転体積回避)。
    // 0 で node 座標 (R=0) を使う。0 は軸ソース OFF と併用前提 (converter で消費)。
    int axisCentroidShift = 1;

    // node-centered 壁 no-slip Dirichlet (既定 1=ON)。1 で 壁ノード u=0 (init+運動量残差射影) +
    // block-DPLUR の壁運動量3行 in-Jacobian decouple (SU2 DeleteValsRowi 相当, plan
    // discretization-node-wall-implicit-dirichlet)。これが無いと壁速度が再循環域でドリフトする。
    // 非 node (cell) / explicit では no-op。0 で旧挙動 (弱形式半割面のみ)。
    int nodeWallDirichlet = 1;

    // (撤去 2026-08-16) nodeAxisDirichlet: 軸ノードを第一内点コピーで置換する対症。保存を破り軸を 1 次化するため
    // 削除。軸ノードは通常 DOF として解き、u_r=0 は壁 no-slip と同じ三点セット (状態ピン+残差 0+Jacobian 行) で課す
    // (node × 軸対称で常時 ON)。plan architecture-node-option-consolidation。

    // node-centered の値の位置は**常にノード座標** (2026-08-16 に旧「双対 CV 重心」規約を撤去)。solver 読込時に
    // centCoords ← ノード座標、双対重心の半径 r̄ は mesh::rEff に退避して軸対称の回転体積 V=r̄·A_planar だけが使う。
    // 根拠 = case/43 (軸半 CV と伸縮メッシュ壁行の O(1) 不整合が消える)。converter 側の書式変更は plan Step 3。

    // node × 軸対称: 軸ノードの u_r=0 は常に三点セット (毎ステージ状態ピン roUy=0 / res_roUy=0 / block-DPLUR の
    // roUy 行単位化) で課す (2026-08-16 に常時 ON 化。残差 0 化だけでは連成 solve で u_r が漏れる, case/43)。

    // node-centered の対流 2 次再構成の目標点は**常にエッジ中点** ½(x_A+x_B) (SU2 と同じ; 2026-08-16 常時 ON 化)。
    // 双対面重心を目標にすると伸縮メッシュで壁法線勾配が接線面に混入し高 AR で発散する (case/43)。

    // 勾配を最小二乗 (LSQ) で計算する (0:既定 Green-Gauss, 1:LSQ 毎ステップ solve,
    // 2:LSQ 係数事前計算 = setup 1回 double 固有分解→擬似逆→float32 係数 gather。推奨は 2)。
    // node-centered の median-dual で GG 面勾配が非一様メッシュの線形場を再現できない (近壁で O(1) 誤差)
    // ための対策。境界は bvar 境界値を LSQ 点として含め勾配を閉じる。methods/discretization.md §7.3/§7.3.1。
    int gradLSQ = 0;

    // gradLSQ=2 の退化判定閾値: M=Σ d̂d̂ᵀ の固有値が λ_max×この値 未満のモードを落とし
    // 擬似逆行列にする (退化方向の勾配 1 次化)。近傍方向が共線/共面なノードの LSQ 発散対策。
    double gradLSQDegenThresh = 1.0e-2;

    // (撤去 2026-08-16) nodeMidpointFx: 値=ノード座標では幾何 fx が自動的に中点相当になるため不要。

    // node-centered 壁摩擦応力 (twall) を「壁ノードに接続する内部双対面 (壁ノード↔内部ノード) の粘性力
    // 集約」で算出する (1:既定, 0:旧=壁ノード勾配 ∇U[W]·S・退化ミラーゴースト dcc ベース)。node の壁ノードは
    // 壁面上に乗り wall_dist=0/dcc 退化のため旧法の twall は近壁で特異スパイク・偶奇振動する。新法は
    // 壁端を no-slip 値 Uxb=0、面勾配を 1/2(∇[W]+∇[I]) として viscousFlux_d と同形の over-relaxed で
    // 各内部双対面の粘性運動量 flux を壁ノード CV に集約し、壁半割面積で割って twall とする (= 流体が
    // 壁 CV に及ぼす粘性 traction; magnitude=τ_w)。運動量残差・y+ には触れない (内部ノード運動量は
    // viscousFlux_d が担うため二重計上回避; twall は出力専用で場は不変)。cell モードは無関係 (従来どおり
    // viscousFlux_wall_d が res+twall を算出)。methods/diffusion.md, plans/diffusion-node-wall-viscous-distance.md §11。
    int nodeWallStressEdgeKernel = 1;
    int thermalMethod;   // 0: calorically perfect (定数 cp/γ), 2: 多成分 thermally-perfect (NASA-9)
    int viscMethod;      // 0: 定数, 1: Sutherland, 2: kinetic theory (Chapman-Enskog)

    flow_float ro;
    flow_float visc;
    flow_float thermCond;
    flow_float cp;
    flow_float gamma;

    // EOS 正値化フロア (dependentVariables で適用)。膨張領域の ro→0,P→0 による速度爆発を防ぐ
    // 安全弁だが、無次元・低圧ケース (例: Taylor-Green は P0=1/γ≈0.71 Pa) では既定 1.0 Pa が
    // 場全体をクランプして解を破壊する。ケースの圧力/密度スケールに合わせ下げられるよう config 化。
    // 既定値は従来ハードコード値と同一なので未指定なら従来挙動 (ビット不変)。
    flow_float pMin  = 1.0;     // 圧力フロア [Pa]
    flow_float roMin = 1.0e-4;  // 密度フロア [kg/m³]
    flow_float tMin  = 1.0e-4;  // 温度フロア [K]

    // 多成分 thermally-perfect gas (thermalMethod==2)。calorically-perfect 経路では未使用。
    int nSpecies = 1;                          // 化学種数 (既定 1 = 単成分)
    std::vector<std::string> speciesNames;     // 混合を構成する化学種名。順序が index s を定義
    std::string speciesDBFile = "";            // 任意: NASA-9/LJ 係数の外部 DB (yaml)。空なら内蔵 DB
    int speciesDiffusionMethod = 1;            // 0: 定数 Schmidt, 1: kinetic theory 混合平均拡散
    double thermoHrefTemp = 0.0;               // >0: 各化学種のエンタルピー基準を h_s(thermoHrefTemp)=0 へ
                                               // オフセット (sensible-enthalpy datum)。非反応流では物理不変
                                               // だが、種ごとに桁違いの生成エンタルピー (H2O≈-13.4MJ/kg) を
                                               // 除いて多成分 implicit の roe/roY 緩和ミスマッチ起因の T ジャンプ
                                               // を抑え安定化する。0 (既定) で従来 NASA 絶対基準・ビット不変。
    // 有限速度化学 (physProp.chemistry, methods/chemistry.md)。既定 off で全経路ビット不変。
    int chemEnabled = 0;                       // 1: 反応ソース項 ω_s・反応熱 Q̇ を有効化 (thermalMethod==2, nSpecies>=2 必須)
    std::string chemMechanismFile = "";        // 反応機構 (Cantera YAML サブセット)
    double chemTmaxReaction = 6000.0;          // 速度式評価の温度上限 [K]
    double chemFreezeBelowT = 0.0;             // この温度未満で反応を凍結 (ω=0) [K]
    int chemJacobianMode = 1;                  // 0: 陽ソースのみ, 1: 対角 point-implicit (src_jac_Y), 2: 全ブロック (Phase 2)
    flow_float Sc = 0.7;                       // 定数 Schmidt 数 (speciesDiffusionMethod==0)
    flow_float Sc_t = 0.7;                     // 乱流 Schmidt 数 (D_t=mu_t/(ro*Sc_t))。
                                               // turbulence.turbulentSchmidt でも設定可 (physProp.Sc_t は後方互換、turbulence 優先)

    // 非平衡凝縮 (4 モーメント方程式 ρg,ρQ2,ρQ1,ρQ0)。methods/condensation/ 参照。
    // Phase 1 はモーメントを受動スカラー (ソース=0) として輸送するのみ。既定 off で従来経路ビット不変。
    int condensation = 0;    // 0: off (既定), 1: on
    int nCondSpecies = 0;    // 凝縮種数。condensation==1 のとき >=1。当面 1 (N2)
    int condModel = 0;       // 凝縮種の物性/核生成/成長モデル。0: N2 (CNT_Iland+Goodheart, CPG),
                             //   1: H2O (Murphy-Koop, CNT+Kantrowitz+Hertz-Knudsen, carrier+TP)
    int condGasSpecies = -1; // carrier+condensible: 凝縮する気相化学種の index (roY{s})。
                             //   -1: pure-condensible (気相=凝縮種, N2 Arthur)。>=0: H2O 等の希薄凝縮 (Wyslouzil)
    int condKantrowitz = 0;  // 核生成の Kantrowitz 非等温補正。0: off (等温 CNT, 既定), 1: on
    int condGrowthModel = 0; // 成長則。0: 既定 (H2O=Hertz-Knudsen, N2=Goodheart), 1: Gyarmathy(熱伝導律速)
    double condGyarmathyC = 3.18; // Gyarmathy の Knudsen 補正係数 1/(1+C·Kn)。標準 3.18 (感度評価用に可変)
    int condTwoTemp = 0;     // 液滴温度 T_d 考慮 (Hertz-Knudsen 経路, 準定常 Hill バランス)。0: 一温度(既定), 1: 二温度
    // 蒸発 (plans/accepted/condensation-evaporation.md)。S=p_v/p_sat<=1 のセルで負成長 (蒸発) と液滴消滅を許可する。
    int    condEvaporation = 1;   // 1: S<=1 で蒸発 (既定, 2026-08-18 検証後に既定化), 0: 蒸発なし (旧挙動, dr/dt<0 を 0 クランプ)
    double condEvapRmin    = 1.0e-9; // 完全蒸発 (4 モーメント 0 化) とみなす体積平均半径 r30 [m]。r30<2*rmin で消滅
    int    condEvapKelvin  = 0;   // 蒸発駆動力に Kelvin 項 p_d=p_sat exp(K_e/r) を含める (0: 平面 p_sat, 既定; 正帰還回避)
    // 平衡凝縮 (plans/accepted/condensation-equilibrium.md): 核生成・成長を経ず各セルで p_v=p_sat(T) の g_eq へ緩和。
    int    condEquilibrium = 0;   // 0: 非平衡 (既定) / 1: 平衡凝縮・緩和形 (ソース S_g=αρΔ/dt, モーメント Q0-Q2 ソース 0)
                                  // / 2: 平衡凝縮・EOS 拘束形 (dependentVariables で (T,g) 同時反転し rog へ射影、rog 輸送は凍結;
                                  //      plans/accepted/condensation-equilibrium-eos.md)
    double condEqRelax     = 1.0; // 平衡緩和係数 α: S_g = α ρ (g_eq - g)/dt_loc
    double condEqDTmax     = 10.0; // 平衡分岐の 1 step 潜熱 ΔT 上限 [K] (非平衡の 1 K では急膨張で平衡に追いつかない: case/44 va2 M5 で S 67 残留)
    double condEqDgMax     = 0.05; // 平衡分岐の 1 step Δg 上限
      //          int isCompressible = physProp["isCompressible"].as<int>();
      //          if (isCompressible == 0) flow_float ro = physProp["isCompressible"]["ro"].as<flow_float>();
      //          flow_float visc = physProp["visc"].as<flow_float>();
      //          flow_float thermCond = physProp["thermCond"].as<flow_float>();
      //          flow_float cp = physProp["cp"].as<flow_float>();

    std::string initial;

    flow_float Pref = 1.0/1.4;
    flow_float Tref = 1.0;

    solverConfig();

    void read(std::string);
    void initTimeIntegrationScheme(int timeIntegration);
    int mainLoopCount() const;
    int perStepIterationCount() const;
    const char* perStepIterationLabel() const;
};

