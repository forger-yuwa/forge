# time_integration: line-implicit v2 試作 (粘性結合 + K 凍結) — 上限確認の限定実験

- status: done
- 起票: 2026-09-02
- related_docs: [methods/time_integration/implementation.md](../../methods/time_integration/implementation.md)
- 先行: [plans/accepted/time_integration-line-implicit.md](../accepted/time_integration-line-implicit.md) (v1: 対流のみ、DDES A/B で不採用)

## 目的 (製品化ではない)

v1 の DDES A/B (case/39 ny160, run_diag_lineimp_*) で確定したのは「**現行 v1 は割に合わない**」
(サブ反復収束 +8% / step 単価 2.44 倍) であって line-implicit の到達点ではない。本 plan は
**上限確認のための一段だけ**を区切って実施し、「現実装が遅かった」のか「この離散化・ケースでは
line-implicit の上限自体が低い」のかを決着させる。

損益分岐: 現収束性能 (nSub 20→17 相当) のままなら line のサブ反復単価を ctrl の **~1.18 倍以下**に
落とす必要がある。K 凍結だけで届かなければ、粘性結合による大幅な nSub 削減が必須。

## 実験項目

1. **K 凍結 (`lineKFreeze`)**: K ブロック抽出 (単位ベクトル×5 の Jacobian 列展開 = v1 の主コスト仮説)
   を物理 step の最初のサブ反復のみ実行し、以後のサブ反復で再利用。純粋な計算コスト下限を測る。
   近似の質: defect-correction の LHS はもともと近似であり、サブ反復間の状態変化は O(ΔQ)。
2. **壁法線 line へのスカラー粘性結合**: line 面に限り K_prev/K_next へ −α_f·I を加算
   (α_f = ν_eff·δ/dcc, 既存 `viscous_diag`=2α と同じメトリック)。対角は既存の 2α を維持
   (行和優位を保つ安全側。真の対称形 [−α, 2α, −α] が line 内で完成する)。
3. **粘性 CFL 割引 (`lineViscousDtRelief` θ∈[0,1])**: on-line セルに限り `setDT` の粘性
   スペクトル半径項 2ν_eff/(ρ·dx_min) を (1−θ) 倍。θ=0 (現状) → 0.5 → 0.8 → 1.0 で
   安定限界と必要 nSub を測る。dx_min は壁近傍で壁法線=line 方向なので方向整合は近似的に成立。
4. **同一プロトコル比較**: case/39 ny160・発達場 (run_0023 res_16000)・300 step 窓・
   総 wall time で ctrl (point) と比較。測るもの: (a) K 凍結後のサブ反復単価比、
   (b) θ を上げたときの安定限界、(c) 同等残差に必要な nSub、(d) 総時間が ctrl を下回るか。

## 判定

- **勝ち筋**: (K凍結後単価) × (必要 nSub / 20) < 1 × ctrl 単価 → line-implicit は本実装
  (完全粘性 5×5 Jacobian・adaptive nSub) を検討する価値がある。
- **負け筋**: θ を上げても nSub が削れない/不安定 → 「対流 defect-correction が律速で、
  粘性 line 陰化の上限は低い」と結論し、v2 は記録して閉じる。

## 変更ログ

- 2026-09-02: 起票。
- 2026-09-02: 実装・検証完了 (**勝ち筋で決着 — 「現実装が遅かった」が正**)。case/39 ny160 DDES
  発達場・300 step 窓 (`run_diag_lineimp2_*`, IC=run_0023 res_16000, ctrl=point 345 s):
  - **コスト分解**: v1 モノリシック 843 s (2.44×) の主犯は **Thomas の毎 sweep 再分解**
    (rhs 非依存の LU + W + Kprev·W 625 積を 5 回/subiter 重複)。factor/solve 分離 (厳密) で
    547 s (1.59×)、`lineKFreeze` (サブ反復間凍結) で **456 s (1.32×) = コスト下限**。
    split/kfreeze は v1 と収束軌道が一致 (分離の厳密性・凍結の無害性を同時検証)。
  - **粘性結合・dt 割引はこのケースでは僅差** (roe +0.02 桁)。理由: 壁第一セル (Δn≈3.5e-5 m)
    の λ_visc/λ_acoustic = 2ν/(Δn·c) ≈ **0.02** — 圧縮性 pseudo-dt の壁法線律速は
    **音響 (c/Δn) であり v1 の対流・音響 line が既に陰化している**。粘性 Δn² 律速が主役に
    なるのは Δn < 2ν/c (本ケース y+≈0.02 相当) の超極薄セルのみ。θ=1.0 でも安定。
  - **本命は pseudo-CFL 引き上げ**: point は cfl_pseudo 2 で発散 (step 99)、**line は 2/4/8 全て
    300 step 安定**。cp4 で ctrl@20 品質に subiter 12-13 到達 (ctrl は 19)、cp8+nSub20 は
    roUx 17 倍深い・ωバースト最大 10.8→8.5 に低減。
  - **実測の勝ち**: `cp4 + nSub13` = **302 s (ctrl の 0.88 倍) で ctrl 同水準の step 終端残差**
    (rms_roUx 5.6e-7 vs 8.0e-7)。同時間 (nSub15) なら ctrl@20 より深い。
  - **判定**: M6 定常ノズル (streamwise 対流律速 → line 無効) と DDES dual-time (壁法線音響律速
    → line 有効) で**律速モードが違う**。DDES 生産推奨 = `lineImplicit:1 + lineKFreeze:1 +
    cfl_pseudo 4 + nSub 13-15` (粘性結合/割引は任意)。
  - **未検証の注意**: 300 step 窓のみ — nSub 削減のバースト余裕 (run_0013 の教訓: バースト最悪値は
    乱流発達と共に成長) は長時間 run で要検証。cp8×nSub 削減の組合せ、adaptive early-exit は今後。
- 2026-09-03: **方向別 dt (`lineDtDirectional: 1`) 追試** — ユーザー指摘「陰化した方向の λ は
  Δτ を縛らなくてよいはず」の実証。line 面の λ (音響込み) を CFL の max から完全除外し、
  壁セルの Δτ を off-line (streamwise) 基準 (×AR≈74, BDF 物理項が対角支配する Newton 的
  レジーム) にする。結果 (`run_diag_lineimp2_dir_cp4*`, 300 step):
  - **cp4+nSub20 で安定** (478 s)。収束は line_cp4 と同等〜微改善 (roUx 3.13→3.37 桁)、
    **ωバースト最大 9.5→6.15 に低減** (壁 Δτ 拡大が k-ω 隔離更新の突発を均す)。
  - **nSub8 は品質不足** (0.77/1.35 桁 < ctrl@20 の 1.53/2.17) — 方向別 dt でもサブ反復の
    収縮率は上がらない。**収縮律速はもはや壁擬似時間でなく off-line lag / SST segregated 結合**。
    ctrl@20 品質に必要な nSub は 12-13 のまま。
  - 運用示唆: 同品質の最速点は cp4+nSub13 (~0.87×) で不変だが、**directional はバースト余裕を
    +50% 積み増す**ので nSub 削減時の安全マージンとして併用推奨。
- 2026-09-03: **記録の補正 (Codex レビュー 3 点)**:
  1. **「壁セル Δτ が AR≈74 倍」は現実装では不成立** — `setDT` の directional 除外は
     `line_prev/next` に一致する**内部 line 面だけ**で、壁ノードの**境界半割面は CFL の max に
     残る**。壁 CV は内部 line 面と境界面が同じ V/S (実測 1.408e-5 m) なので、境界面が同じ音響
     制約を残す。**AR 倍化は壁の 1 つ内側以降の line CV のみ**。壁境界面 λ の除外可否
     (Dirichlet 行 decouple 済みなら安全かもしれない) は未検証の将来項目。
  2. ωバースト半減の帰属は「directional dt により半減 (実測)」に留める — 「壁 Δτ 拡大が均した」
     は未分離の推論 (1. の通り壁セル自身の Δτ は伸びていない)。
  3. 収縮律速の候補は **3 者**: off-line lag / segregated SST / **2次 KEEP RHS×1次 FVS LHS の
     defect-correction 不整合** (M6 の cfl×relax 飽和と同根の可能性)。
- 2026-09-03: **推奨構成の直接 A/B (directional cp4 × nSub13/15)** — 結果は下記追記。
- 2026-09-03: **推奨構成の直接 A/B 完了** (`run_diag_lineimp2_dir_cp4_nsub13/15`, 300 step):
  - `directional cp4 + nSub13`: **317 s = ctrl の 0.92 倍**、step 終端品質は同等〜良
    (roUx 4.4e-7 vs ctrl 8.0e-7, roe 2.45e-4 vs 2.13e-4)、ωバースト最大 6.15 (ctrl 10.8)。
  - `directional cp4 + nSub15`: 363 s (1.05 倍) で全量 ctrl より良い (roUx 2.0e-7, roe 1.33e-4)。
  - 非 directional nSub13 (302 s, 0.88 倍) より directional は +5% 遅いが、バースト余裕 −35%
    (9.5→6.15) を買う取引。**生産候補 = directional cp4+nSub13 (速度優先) / nSub15 (品質・余裕優先)**。
  - 残作業 (本採用の条件): 生産再開時に nSub15 で長時間 (数万 step) のバースト余裕検証を 1 本。
- 2026-09-03: **壁境界半割面 Λ の除外実験 (`lineDtWallRelief`, opt-in 診断)** — Codex 提案の
  3 分岐実験。wall 種 bcond の境界面のみ明示フラグ (19642 面) で setDT の max から除外
  (inlet/outlet/periodic は残す。LHS の境界 A⁺ 対角は不変)。結果は**分岐③: 発散**
  (`run_diag_lineimp2_dirwall_cp4`, step ~80-100 で ro が非有限 → detectNaN 停止。
  対照 `dir_cp4` 再走は完走)。最初に落ちたのは **ro** (種の位置は dump 削除により未特定)。
  **解釈の補正 (Codex レビュー)**: node 強壁の RHS 壁対流流束は F_b=(0, pn, 0) で壁質量流束は
  恒等的にゼロ → **連続行に「欠けている物理的境界 Jacobian」はほぼ無い** (境界面の A⁺S 対角は
  もともと人工的な対角強化 [gpu-implicit-plan.md 参照])。今回証明されたのは
  **「壁面 Λ 由来の pseudo-time 質量項 (V/Δτ) が現行反復の壁 CV 減衰として必要」**まで。
  除外時の発散原因は境界 Jacobian 欠落とは確定できず、LHS/RHS 不整合・K 凍結・非線形更新過大を
  含む複数候補のまま。将来壁 Δτ を広げたい場合の手順: 発散直前のセル位置・δρ/ρ・V/Δτ を確認
  → 局所緩和 / 正値性ガード / K 凍結解除で原因を分ける (境界 Jacobian 追加はその後)。
- 2026-09-03: **dt_local 壁距離プロファイル (dir_cp4, `dt_profile.csv`)** が Codex 指摘を定量確認:
  壁ノード帯 (wd<5e-5) の dt_med **1.9e-6** vs 内側 **1.4-1.7e-5 (~8 倍差)** — directional の
  AR 倍恩恵は第一内点以降のみで、壁 CV 自身は境界面 Λ で絞られたまま。内側は ~7.5·dt で
  BDF 物理項支配レジーム到達済み。⇒ **「壁擬似時間律速の完全除去」は未達のまま**が正確で、
  収縮律速 3 候補 (off-line lag / segregated SST / defect-correction 不整合) の切り分けは、
  壁端点を除去できない以上、別経路 (例: SST 連成陰化 or LHS 2 次化) からになる。
- 2026-09-03: **収縮律速の消去法完了 (診断 4 腕, 各 100 step 窓 20-95, 対照 dir_cp4)**:
  ① `FORGE_FREEZE_TURB=1`: roe/roUx 収縮 **2.07/3.37 桁で対照と完全一致** → segregated SST は
  平均流収縮を妨げていない (注: μt は凍結 k,ω から毎回再計算されるので「ほぼ固定」)。
  ② `implicitRelaxSST` 0.7/1.0: 平均流不変。**ω 収縮はむしろ悪化** (0.39→0.26→0.05 桁) =
  ω は relax 0.5 で縁辺安定。SST 2×2 line 化は ω 品質・バースト向けで平均流 nSub は減らない。
  ③ `nStepInner` 5→10: **全チェックポイントで収縮一致** → 線形系は 5 sweep で解き切れている
  (off-line lag は近似 LHS の線形解の中で処理済み)。
  ⇒ **平均流のサブ反復収縮律速 = defect-correction 不整合 (1次 FVS LHS × 2次 KEEP+ES RHS) に確定**。
  根治候補: line K の FD 化 (RHS 整合をライン方向だけ入れる, v2 機構流用可) / JFNK
  (DPLUR+line を前処理に) / adaptive early-exit は律速に依らず ~20-30% 得。
- 2026-09-03: **訂正 (Codex レビュー・重大)**: 上記「消去法完了」の freeze A/B は**無効だった** —
  `FORGE_FREEZE_TURB` のゲートは定常経路 (`implicitNonlinearUpdate`) 専用で、dual-time の
  `advanceImplicitDualTime` は `applySSTPointImplicit` を無条件に呼んでいた (同一経路の 2 回走で
  一致は自明)。**「defect-correction 不整合が単独原因として確定」は撤回** (commit d3b55d6d の
  結論を本追補で差し替え)。dual-time にもゲート+起動ログを追加して再走 (`run_diag_lineimp2_frz2`):
  - 凍結の実証: res_100 の roK/roOmega は**内部 (wd>1e-4) で IC とビット一致**。変化は壁ピン帯
    (wd≤3.2e-5, roOmega の 2.5% ノード = `nodeOmegaWfDirichlet` 等の壁 Dirichlet — 両腕共通) のみ。
  - **結果: 平均流収縮は active と実質一致** (roe@20 2.07/2.07, roUx 3.37/3.36 — 0.01 差が
    「今回は本当に別経路」の傍証)。
  ⇒ 支持される結論は「**この窓では SST フィードバックは平均流収縮の主律速でない**」まで。
  nStepInner 一致の読みも訂正: 「線形系を解き切った」でなく「**同じ近似 LHS をさらに解いても
  非線形収縮は改善しない**」(sweep 増は無価値、の実用結論は不変)。
  残候補: defect-correction 不整合 / implicitRelax=0.5 の緩和上限 / lineKFreeze の古い LHS /
  off-line Jacobian の忠実度 / 壁 CV の pseudo-time 制約 → ir0.7・kfreeze-off の追試で切り分け続行。
- 2026-09-03: **残候補 probe 2 本 (各 100 step, 対照 dir_cp4)**:
  - `implicitRelax` 0.5→0.7 (`run_diag_lineimp2_ir07`): **発散** (NaN)。relax 0.5 は任意の減衰でなく
    **安定必須** — 近似 LHS の defect-correction が要求する緩和で、収縮率の上限 (毎 subiter ≤半歩)
    を課すが外せない。⇒ この候補は「defect-correction 不整合」に吸収される (LHS が忠実なら
    relax→1 で Newton 級収縮が許されるはず)。
  - `lineKFreeze` off = K/LU 毎 subiter 再構築 (`run_diag_lineimp2_nokfrz`): **収縮完全一致**
    (roe@20 2.07 / roUx@20 3.37)。**「K 凍結の古い LHS」候補は消去** — 凍結はコスト削減のみで
    収縮に無害と直接実証。
  ⇒ 総括: 安く動かせる要素 (SST/sweep 数/K 鮮度/relax) を全て振っても収縮は 2.07/3.37 に固着。
  **残るのは近似 LHS の忠実度ファミリー (1次 FVS×2次 KEEP の defect-correction 不整合 +
  off-line Jacobian 忠実度 + 壁 CV pseudo-time 制約)** で、単独犯の特定はこの窓では不能。
  打ち手は変わらず: adaptive early-exit (無条件) / line K の FD 化 / JFNK。
- 2026-09-03: **CFL-nSub マップ + 残差局在の測定 (Codex 提案の最安切り分け, run_diag_cflmap_*)**:
  - **Phase 1 (cp マップ, directional+nSub20 固定, 各 100 step)**: cp8 OK / cp12 OK / **cp14 発散
    @step3 / cp16 @step2 / cp24 @step1** → 上限は [12,14)。**収縮は cp8 で飽和** (roe@10:
    cp4 1.22 / cp8 1.27 / cp12 1.28) — BDF 項が対角支配に入った後は Δτ 増が効かない構図。
    振動兆候ゼロ (roe 増加 subiter 率 0.000)。min ro=1.16/min P=1.02e5 で EOS 床とは無縁
    (M6 の床洗浄とは終端症状が別物)。**運用最良点 = cp8** (上限へ 1.5 倍の余裕)。
  - **発散モードの指紋 (cp14 の inner 履歴)**: (ro, roUy, roe) が **subiter 間で単調増幅**
    (step0 subiter11 から成長・×1.6/3subiter)、roUz/roK は無傷 → 壁法線 (y) の音響/質量-
    エネルギー系の反復モード。NaN 拡散は 1 step で 86-100% に及ぶため種セル特定には
    per-subiter 検査が必要 (将来項目)。
  - **Phase 3 (残差局在, cp8 + `FORGE_OUT_RESIDUALS=1` [main.cpp に env 追加, res_*/dq_block_new_*
    を h5 出力へ])**: 相対補正要求 η=|dq|/max(|Q|,Q_ref) は**丘頂直前のせん断層 (x/h≈8.6,
    wd≈0.09-0.11h) に全 5 式が同座で局在**。壁ノード帯は p99 で BL 上部より小さく **壁 CV は
    支配的でない**。max η ~4e-5 と絶対値も健全。⇒ 判定木では **off-line (streamwise/spanwise)
    枝** — line K の忠実化だけでは動かず、次に効かせるなら streamwise 第 2 ライン族 (x 周期 →
    cyclic Thomas) か、物理活性域の要求として受容 (現状 cp8+nSub13 で実用充分) の二択。
- 2026-09-03: **診断修正 + 局所収縮率 g + cp8 nSub 直接 A/B (Codex 補正 3 点への対応)**:
  - 診断修正: `FORGE_OUT_RESIDUALS` の dq 出力を **dq_block_old_*** (swap 後の最終補正) に修正
    (旧 dq_block_new_* は 1 sweep 前 — 定性結論は不変だが正した)。`FORGE_RESID_SNAP=1` を追加
    (subiter0 の res/dq を res_*_m / dq_*_new スロットへ退避し出力 → g を場で測れる)。
  - **局所収縮率 g=|dq_final|/|dq_sub0| (cp8, step100, run_diag_cflmap_cp8_g)**: 中央値 0.000-0.004、
    p99 ≤0.039、**増幅 (g>1) ノードはゼロ**。帯別も壁ノード帯 0.010 vs コア 0.004 の緩い勾配のみ。
    ⇒ **局所的な数値律速は存在しない — 収縮率は全域一様**。判定木は「off-line 枝」でなく
    **「全域一様 → 1次 FVS LHS × 2次 KEEP RHS の演算子不整合が本命」** に確定 (先の off-line 判定は
    η 最大位置の誤読で撤回済み)。streamwise 第 2 ライン族は根拠を失い保留。
  - **cp8 × nSub 直接 A/B (300 step, 同一 IC)**: cp8+nSub13 = **315 s (ctrl の 0.91 倍) で
    全指標 ctrl 同等以上** (ro 6.6e-10 / roUx 3.3e-7 / roe 2.13e-4, ωバースト 6.02)。
    cp8+nSub12 = 291 s (0.84 倍) だが roe 3.0e-4 (+42%) の妥協。cp4+nSub13 とはほぼ等価。
  - **生産推奨 (最終)**: `lineImplicit1 + lineKFreeze1 + lineDtDirectional1 + cfl_pseudo 8 + nSub 13`
    (同品質 0.91×・バースト −44%)。速度優先なら nSub12 (0.84×, roe 微妥協)。
    さらなる短縮は演算子不整合ファミリー (JFNK / LHS 散逸整合) か adaptive early-exit。
- 2026-09-03: **訂正 (Codex レビュー): g 解析のマスク欠陥 — 「全域一様・FVS-KEEP 確定」を撤回**。
  先の「増幅ゼロ」は分母フィルタ (ini>p50) が「初期 dq が小さく後から成長する点」を定義から
  除外していた。マスクなし再計算 (`analyze_g.py`, マスク明記の再現スクリプトとして case に残置):
  **roe g>1 = 15,242 ノード / roUy 3,972** (最終側マスク |dq_fin|>1%max でも roe 7,613 / roUy 1,135)。
  成長点の空間分布が示唆的: **丘頂域 (x/h 0-1, 7-9) に集中し、roUy 成長点は壁至近
  (wd/h med 0.004)** — **cp14 発散モードと同じ (ro,roUy,roe) ファミリー・同じ場所**。
  ⇒ 支持される言明: 「subiter0 で活発な領域の大部分は良く収縮する。一方、遅れて補正が成長する
  点が丘頂・近壁に実在し、その方向・原因は未判定」。off-line 確定でも FVS-KEEP 確定でもない。
  `FORGE_RESID_SNAP=<m>` を subiter 指定可能に拡張し、m=0/5/10/15/19 の履歴で成長点の
  軌跡を追う (run_diag_cflmap_cp8_g{,5,10,15})。
- 2026-09-03: **長回し検証の条件を訂正**: 最終採用を nSub13 とするなら**長回しも nSub13 で行う**
  (nSub15 の長回しでは nSub13 のバースト耐性は証明できない — run_0013 の nSub12 が step3839 で
  破綻した前歴があるため重要)。
- 2026-09-03: **subiter 軌跡解析 (m=0/5/10/15/19, run_diag_cflmap_cp8_g{,5,10,15}) — 「成長点」の正体は床**:
  同一 run 内の頑健な比較 (m0 と final は同ファイル) で、roe の g>1 集合は **subiter0 で対照の
  1/60 (8.6e-3 vs 5.1e-1) = 最初からほぼ収束済みの点**であり、subiter19 には両集合が**同じ dq 床
  (~1.6e-2) に合流**する (roUy も同型: 床 ~2.5e-5)。⇒ **cp8 に局所増幅モードは存在しない**。
  サブ反復は ~15 回で全域共通の dq 床に到達して止まる — これが収縮の実像で、nSub13 前後が
  最適という運用結果とも整合。**床の起源 (LHS/RHS 不整合のリミットサイクル / DDES 乱流の
  物理揺らぎ / relax 上限) は未分離** — ここから先は凍結場 (物理揺らぎ排除) での床測定などが要る。
  注意: run 間の最終 dq 再現性は中央値 60% 差 (深収束 dq はジッタ支配) — 集合統計のみ有効、
  点単位のクロス run 比較は不可 (`analyze_g.py` に明記)。
- 2026-09-03: **定常 M6 ノズルへの v2 適用可否 (run_0015_cflsweep/tp_cfl6r07_line2 等, warm 2000 step)**:
  **不可 — v2 でも定常側は損のまま**。① line v2 (directional なし): 末尾残差 −10% だが
  **step 単価 1.92 倍** (15.9→30.6 s) — 2D 小メッシュはライン 1250 本=1250 スレッドで占有率が
  壊滅し、DDES (9821 本, 1.32×) と逆に v1 比でも重い。② **directional は cfl6/8 とも発散 @step25**:
  壁法線音響で絞られていた BL セルの小 Δτ が実は安定マージンだった (dual-time の BDF 対角の
  ような保護が定常には無い)。⇒ **v2 は dual-time DDES 専用の武器**と結論。定常 RANS の生産は
  従来どおり point-DPLUR + cfl6/8+relax0.7。
- 2026-09-03: **訂正 (Codex レビュー) + back-to-back 正式 A/B (同一バイナリ・同一 IC,
  run_0015_cflsweep/tp_cfl6r07_{point2,line2b,linemono,line_ni3,line_ni2,linedir2})**:
  前項の「v2 でも定常側は損」「DDES 専用」は言い過ぎだった。正式値:
  | 腕 (cfl6+r0.7, warm 2000 step) | Time | 末尾 roe |
  | point (ni5) | 15.99 s | 0.626 |
  | **split line (ni5)** | **28.97 s (1.81×)** | 0.556 (−11%) |
  | mono line (ni5, FORGE_LINE_MONO=1) | 43.90 s (2.75×) | 0.556 |
  | split line ni3 | 25.08 s (1.57×) | 0.577 |
  | split line ni2 | 発散 @785 | — |
  - **factor/solve 分離は定常でも有効** (2.75×→1.81×, −34%) — 「効かない」は誤り。ただし最良
    チューニング (ni3) でも 1.57× コスト vs roe −8% で**この case では採算不成立**。ni2 発散 =
    off-line lag に sweep ≥3 が必要。「1250 line の占有率壊滅」は未プロファイルの仮説に格下げ。
  - **directional 発散の帰属を訂正**: 種は **x/r_t 71-89 の下流域・全断面** (壁至近 3%・近軸 1%,
    res_nan_26.h5 保持) — 「BL セルの露出」でなく**既知の streamwise 内部モード**が、directional で
    Δτ の上がった下流域セルで先に点火したもの (cp6/cp8 同 step の理由も状態依存で説明がつく)。
    定常律速 = streamwise 対流、の既存結論を補強。
  - **結論の正確な形**: case/45 M6 定常の現設定では non-directional v2 は採算不成立・directional は
    使用不可。**他の定常高 AR 問題への一般化はしない** (3D 定常や別ケースは未検証)。
