# 35. uniform periodic box (Taylor-Green)

全 6 面 Cartesian periodic の 3D 構造化 hex box (32³, `mesh/tg.msh`)。一様流/Taylor-Green 初期値で
周期境界・3D 経路の最小検証に使う。3D median-dual (M4, [`discretization-median-dual-3d`](../../plans/active/discretization-median-dual-3d.md))
の双対メッシュ生成 (体積総和・closure) 検証ケースでもある。

## 計算 run 一覧

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_uniform` | cell モード一様流ベースライン (convMethod=0) | `res_*.h5`, `residual_history.csv` | ref |
| `run_0002_node_dualcheck` | node モード (`discretization: node`) で 3D median-dual 生成の幾何検証 | 体積 dual=248.05=primal (relErr 3.5e-8)・closure normalized 2.5e-6・負体積 0 (`convert.log`)。solver run は未実施 | active |
| `run_0003_node_periodic_smoke` | node 3D で全面 periodic の forge smoke (20 step, §4.5 実装前) | **step 4 で DIVERGED**。半割面を主面接続する旧挙動は幾何が誤りで発散 → §4.5 実装の動機 | 削除済み(2026-08-31) |
| `run_0004_node_periodic_freestream` | **§4.5 実装後**: 構造化 periodic 立方体(8³ hex, 全方向 periodic)+ 一様 +x 流で free-stream 保存検証 | **PASS**: 全 50 step で残差**機械ゼロ**。union-find `217 slave merged`=729−512(triperiodic 一意 CV)で完全一致。`conv main planes periodic 0`(移流から periodic 除外) | ref (検証エビデンス) |
| `run_0005_node_periodic_tg` | node 全面 periodic + Taylor-Green IC (M0.1) explicit RK3 | step4 DIVERGED。union-find `3169 merged`=35937−32768 正。発散は periodic 起因でなく node+低マッハ+explicit 既知不安定(run_0006 で実証) | 削除済み(2026-08-31) |
| `run_0006_node_slip_tg_isolation` | **切り分け**: run_0005 の periodic を slip に置換(periodic 無し)、他同一 | step4 で**同様に DIVERGED** → TG 発散は periodic 実装でなく node+低マッハ+explicit RK3 の不安定と確定(periodic 無罪) | ref (切り分けエビデンス) |
| `run_0007_node_periodic_tg_lmp2` | run_0005 + lowMachPrecond=2 (explicit) | step1 DIVERGED。explicit では低マッハ安定化せず(node 低マッハは implicit 必須)。非自明 periodic は implicit(§4.5.7)待ち | 削除済み(2026-08-31) |
| `run_0008_node_periodic_tg_implicit` | node 全面 periodic + TG + implicit block-DPLUR(lowMachPrecond=2, cfl_pseudo=1, nStepInner=5) | step64 で **NaN 発散**(explicit step4 より大幅に粘るが発散) | 削除済み(2026-08-31) |
| `run_0009_node_slip_tg_implicit_iso` | **切り分け**: run_0008 の periodic を slip に置換、他同一 | step99 まで **NaN なし**(plateau/roe rising だが発散せず)。**implicit では periodic だけ NaN 発散** → §4.5.7 Jacobian 行 fold が implicit periodic 安定化に必要と確定(RHS gather だけでは LHS 不整合) | ref (§4.5.7 必要性の根拠) |
| `run_0010_node_periodic_tg_implicit_mirror` | **§4.5.7 実装後**: run_0008 + sweep ごと dq ミラー(slave=master) | **step100 完走・NaN なし・場健全**(ρ 0.83–1.02, P/T/roe 有界, 渦発達)。run_0009(slip)と同等プラトーになり periodic 固有発散消失。**implicit periodic 安定化達成** | ref (§4.5.7 検証エビデンス) |
| `run_0011_node_periodic_tg_gradgather` | **勾配 periodic gather 実装後**: run_0010 同設定で勾配 gather 有効 | step100 安定・場健全(ρ 0.825–1.017)で run_0010 と実質不変、seam viscous が片側→両側へ微修正(回帰なし)。free-stream も機械ゼロ維持 | ref (勾配 gather 検証) |
| `run_0012_gradcheck` | **勾配 gather 定量検証**: TG 解析勾配 ∂ux/∂x=M0cosx cosy cosz と照合 | 境界(合併)勾配 mean\|err\|=5.3e-4(\|g\|max 0.398 vs 解析0.4)で内部5.8e-4と同等、master==slave一致(2.9e-7)。periodic境界勾配が内部同等に正しいと確認 | ref (勾配検証エビデンス) |
| `run_0014_cell_keep_cbd_pure` | **L1a/b 市松診断** (ラダー L1, [survey §10](../../notes/investigations/convection-central-scheme-oscillation-control.md)): cell 32³, 一様 M0.1 流+圧力のみ市松摂動 ε=1e-3 (`make_checkerboard_ic.py`), pure KEEP (keepDissType=0) | **null-mode 実証**: 全 rms 残差が厳密ゼロ (中心平均が市松を相殺し原理的に見えない)・A_cb=1e-3 が 400step 後も**不変** | ref (L1a null-mode エビデンス) |
| `run_0015_cell_keep_cbd_diss` | 同 IC + **keepDissType=1 (σ=0.15, lowMachPrecond=1)** | rms_roe のみ非ゼロ (散逸が Δp 経由でのみ市松を見る=設計通り)。**A_cb 1e-3→8.0e-10 (6桁減衰/400step)** | ref (L1b 減衰エビデンス) |
| `run_0016_cell_slau_cbd_ref` | 同 IC + SLAU (対照) | SLAU は市松を ro/roUx/roe 全部で見る (rms 即非ゼロ)。A_cb 1e-3→1.2e-8 | ref (対照) |
| `run_0017_cell_keep_cbd_diss005` | 同 IC + keepDissType=1 **σ=0.05** (較正) | A_cb 1e-3→1.3e-7 (~4桁減衰/400step)。TGV KE cost 2.7% (case/09 run_0024) とのバランスで **σ=0.05 を既定に採用** | ref (σ 較正) |
| `run_0018_cell_keep_cbd_dissmat005` | 同 IC + **keepDissType=2 (matrix ES, σ=0.05)** | A_cb 1e-3→**7.1e-8** (scalar σ=0.05 の 1.3e-7 より良) — 音響のみ c' 込み固有値の選択的減衰が市松に有効 | ref (Step 2 matrix L1b) |
| `run_0019_cell_keep_cbd_dissmat005_fullc` | matrix σ=0.05 を**フル c** (lowMachPrecond=0, 旧結線) で | A_cb 1e-3→**9.2e-9** (最速)。ただし TGV KE cost 4.35% (case/09 run_0026) と3.2倍悪化 → c' の価値の実測 | ref (c' vs フル c 切り分け) |
| `run_0020_cell_keep_cbd_dissmat005_knobchk` | **keepDissCprime ノブ分離後**の動作確認 (lowMachPrecond=0 + 既定 keepDissCprime=1) | A_cb→7.5e-8 ≈ run_0018 (c' 挙動を lowMachPrecond と独立に再現、差は atomicAdd ノイズ級) | ref (ノブ分離検証) |
| `run_0021_cell_keep_cbd_dissmat_fullc_s0015` | フル c で σ 弱化の代替可否: σ=0.015 | A_cb→**3.9e-8** = c' σ=0.05 (7.1e-8) と同等以上。**単一マッハ領域なら σ 弱化で c' 代替可** (KE も 1.10% = case/09 run_0027) | ref (σ掃引・c'代替検証) |
| `run_0022_cell_keep_cbd_dissmat_fullc_s004` | 同 σ=0.04 | A_cb→3.1e-10。**σ=0.05 (9.2e-9) より速い非単調** = 2Δ 定在音響波の過減衰領域を観測 | ref (σ非単調の記録) |
| `run_0023_node_keep_cbd_pure` | **node L1**: node 変換 (median-dual 35937CV) + 市松 IC (COORD パリティ・17969/17968)、pure KEEP | **node でも null-mode 実証** (step0 全 rms ~1e-10 = 機械ゼロ、A_cb 400step 不変)。explicit RK3 でも発散せず | ref (node L1a) |
| `run_0024_node_keep_cbd_dissmat005` | node + matrix σ=0.05 (既定) | **[訂正済]** 初報の「継ぎ目 3.3e-4 停滞」は**測定バイアス** (境界サブセットの sum(par)=2≠0 × 平均圧の DC 項 = 予測 3.25e-4 と一致)。共分散射影では **all/interior/boundary が 6.8e-6 に完全一致 = 継ぎ目弱化なし・一様減衰**。cell (7.1e-8) との桁差は **node dt が半分 (4.46e-2 vs 8.92e-2, 周期 half-CV 体積の CFL 律速) なだけ**で、**単位時間あたり減衰率は node 0.116 vs cell 0.122 桁/時間で一致**。plot_checkerboard_decay.py は共分散射影に修正済 | ref (node L1b・バイアス訂正) |
| `run_0025_cell_keep_tp_cbd_pure` | **Step 3 (TP 単成分) L1**: thermalMethod=2 [N2], 次元付き一様場 (P_eff=238kPa, T=705K, u=10m/s) + 市松 ΔP (modify-p-only IC) + pure KEEP | A_cb 9.58e-4 が 400step **不変** = **null-mode は TP でも成立**・発散なし | ref (TP L1a) |
| `run_0026_cell_keep_tp_cbd_dissmat005` | 同 + **matrix ES (TP 版: w=(g-½q)/T 系, S=ρ/(2γR)/ρ/cp/ρT)** σ=0.05 | A_cb → **1.18e-7 (~3.9桁/400step)**・NaN なし。TP matrix 散逸の動作実証 | ref (TP L1b matrix) |
| `run_0027_cell_keep_tp_cbd_scalar005` | 同 + scalar ES (TP: Δroe=保存量ジャンプ・c=sonic) σ=0.05 | A_cb → 1.10e-7・NaN なし | ref (TP L1b scalar) |
| `run_0028_cell_keep_mc_cbd_pure` | **Step 4 (多成分) L1**: [N2,O2] Y=0.7/0.3 一様混合 (P=230kPa/T=709K) + 市松, pure KEEP | 一様混合 400step 安定 (プローブ)・市松に対し **null-mode 不変** (多成分でも成立) | ref (MC L1a) |
| `run_0029_cell_keep_mc_cbd_dissmat005` | 同 + matrix ES σ=0.05 | **★未解決バグ**: 最初の 50step は正常減衰後、**2.3e-4 でプラトー** (単成分は 1.2e-7 まで減衰) | ref (MC matrix バグ再現) |
| `run_0030_cell_keep_mc_comp_smoke` | 組成半割 (N2\|O2) + matrix σ=0.05 (安定性 smoke) | 400step NaN なし・場有界 (P 209-244kPa, T 681-751K)。ゼロ濃度種の正則性 (Y ln X ガード) は動作 | ref (MC 安定 smoke) |
| `run_0031_cell_keep_tp_rerun_newbin` | **bisect①**: 単成分 TP を新バイナリ (_mix 経由 ns=1) で再実行 | 9.6e-4→**1.2e-7** = 旧バイナリと一致 → **ns=1 経路に退行なし** | ref (bisect) |
| `run_0049_node_h2_ignition` | **有限速度化学 Phase 1 の GPU 検証** ([plan](../../plans/active/chemistry-finite-rate-h2.md)): node 8³ 全面 periodic = 定積反応器, TP 9 種 (Jachimowski 9sp20r), 量論 H₂-air 1200 K 1 atm, explicit RK3 dt 5e-9, `chemistry.enabled: 1`。IC は `make_h2_ignition_ic.py`, 図は `plot_h2_ignition.py` | **Cantera と一致**: 着火遅れ 32.0 µs (参照 32.2), 平衡 T 2948 K (参照 2944, +0.15 %), 全ノード一様 (T 差 0.004 K)。`h2_ignition_vs_cantera.png`。**副産物**: node 周期の化学種残差 gather 欠落を発見・修正 (`periodicNode_d.cu`) | ref (化学検証エビデンス) |
| `run_0050_regr_node_periodic_freestream` | 化学 Phase 1 実装後の回帰: run_0004 (node 周期一様流) を化学 OFF で再実行 | 全残差 機械ゼロ (旧と同一) — `periodicNodeGather` 拡張の副作用なし | ref (回帰) |
| `run_0051_regr_cell_keep_mc_cbd_pure` | 同: run_0028 (cell TP 2 成分 KEEP null-mode) を再実行 | 全残差 機械ゼロ (旧と同一) — TP 多成分経路ビット不変 | ref (回帰) |
| `run_0032_cell_keep_mc_y10_bisect` | **bisect②**: 2種だが Y=[1,0] (物理的に純 N2) + matrix | **プラトー再現 (2.2e-4)** → バグは実在混合でなく **ns≥2 コード経路**。組成は厳密に [1,0] 維持・res_0 は ns=1 と 7桁一致なのに **step0 の rms_roUx のみ 0.07% 差** (rms_ro/rms_roe は一致) — 原因未特定 | ref (bisect・バグ切り分け) |
| `run_0033_cell_keep_mc_y10_scalar_bisect` | **bisect③**: 同 Y=[1,0] で **scalar** (type 1) | **1.1e-7 までクリーン減衰** = 多成分ソルバ機構 (種輸送/renormalize/依存変数) は無罪。**バグは matrix×ns≥2 に限定** | ref (bisect・scalar 多成分 OK) |
| `run_0034_cell_keep_mc_cbd_dissmat005_fix` | **バグ根治後**: 実混合 N2/O2 0.7/0.3 + matrix σ=0.05 (run_0029 と同 IC) | A_cb 9.77e-4→**4.79e-8 (完全減衰)**。プラトー根治を実混合で確認 | ref (Step 4 matrix 根治検証) |
| `run_0035_cell_keep_tp_regr_fix` | 根治後の単成分回帰 (run_0031 と同 IC) | 1.18e-7 で不変 (回帰なし) | ref (回帰確認) |
| `run_0036_cell_keep_mc_comp_smoke_fix` | 根治後の組成半割 smoke | NaN なし・場有界 (P 212-242kPa, T 684-741K) | ref (smoke 再確認) |
| `run_0037_node_keep_cbd_regr_advgaugebin` | **回帰**: 移流基準差分 (advGauge, [plan §8](../../plans/active/convection-freestream-preserving-flux.md)) 実装後バイナリで run_0023 (node KEEP 市松 pure) を再実行 | step0 全 rms ~1e-10 機械ゼロ・400step 安定 = run_0023 と同挙動 → **node KEEP 経路に advGauge 実装の影響なし** (node は roRef 未対応で config が fail-fast) | ref (advGauge node 回帰) |
| `run_0038_cell_keep_velpert_dbg` | **切り分け**: M0.1+速度擾乱 a=1e-2 (case/33 発散の対照; Cartesian・全面 periodic・cell KEEP pure)。サブ `steady_variant` は unsteady:0 版 | unsteady:1 は 50step 安定 (物理減衰) だが **unsteady:0 (定常局所dt) は Cartesian でも step11 発散** → case/33 run_0008-0010 の発散は「音響過渡×定常局所 dt」の一般限界との交絡と確定 ([plan §8.4](../../plans/active/convection-freestream-preserving-flux.md)) | ref (定常モード限界の切り分け) |
| `run_0039_node_keep_cbd_dissmat005_jump` | **L1 node × keepDissJump=1** ([recon-jump plan](../../plans/accepted/convection-keep-diss-recon-jump.md)): run_0024 と同一 IC + 再構成ジャンプ散逸 | A_cb 1e-3 → **6.73e-6**/400step = 生ジャンプ (6.78e-6) と同一 → **市松減衰は再構成で無傷** (2Δ の中心勾配ゼロ性の実機確認) | ref (recon-jump L1 node) |
| `run_0040_cell_keep_cbd_dissmat005_jump` | L1 **cell** × keepDissJump=1 (run_0018 と同一 IC) | A_cb → 6.29e-8 = 生 (7.09e-8) と同等以上 → cell 経路も健全 | ref (recon-jump L1 cell) |
| `run_0041_node_keep_cbd_regr_jumpbin` | 回帰: keepDissJump 実装後バイナリで run_0024 (既定 off) 再実行 | A_cb → 6.781e-6 = run_0024 (6.782e-6) と 4 桁一致 → **既定 off に回帰なし** | ref (回帰エビデンス) |
| `run_0042_node_keep_cbd_dissmat005_jump2` | L1 node × **keepDissJump=2 (sign-property クリップ)** | A_cb → 6.99e-6 = 生 (6.78e-6)/jump=1 (6.73e-6) と同等 → **証明付き ES 化でも市松減衰維持** | ref (jump=2 L1) |
| `run_0043_node_keep_cbd_dissmat005_precond` | **L1 node × `keepDissPrecond=1`** (Turkel 前処理音響散逸, [plan](../../plans/active/convection-keep-diss-lowmach-precond.md)): run_0024 と同一 IC + 前処理 2×2 (Δp 散逸 ∝c²/Ur 増強) | A_cb 1e-3 → **4.78e-8 (4.3桁/400step)** = c' 版 run_0024 (6.78e-6) の **142 分の 1** → 真の市松への設計効果を実証 | ref (precond L1 node) |
| `run_0045_cell_keep_cbd_precond_jump2` / `run_0047_node_...` | **cb 補正の L1 基準** (precond+jump2, C_cb=0): cell A_cb 1e-3→1.7e-7@50step, node →2.0e-5 プラトー (周期 seam 残留モード, cb 対象外) | ref (cb L1 基準) |
| `run_0046_cell_..._cb010` / `run_0048_node_..._cb010` / `dbg_*` | **L1 × `keepDissCbCoeff`** ([plan](../../plans/active/convection-keep-cb-pressure-correction.md)): **cell C_cb=0.02 で 7.1e-4→1.0e-6 を 5 step** (基準は 50 step で 1.7e-7)・node 0.1 で早期減衰 10 倍速 = 符号/設計どおり市松狙い撃ち。**cell C_cb=0.1 (実効 C/ε=1.0) は explicit RK3 CFL0.5 で市松モード overshoot 発散** (×3.7/step) → explicit では C/ε≲0.2 に制限、implicit は本番掃引で確認 | ref (cb L1 検証) |
| `run_0044_cell_keep_cbd_dissmat005_precond` | 同 **cell** (run_0018 対照) | A_cb → **1.51e-8** vs run_0018 の 7.09e-8 (**4.7×**) = cell でも改善 | ref (precond L1 cell) |

**matrix×多成分プラトーの根治記録 (2026-07-19 深掘り)**: 真因 = **ΣρY_k ≠ ρ の共通モードノイズ** (ρY と ρ は別カーネル/別 atomicAdd 順で更新されるため ~1e-7 の不一致が発生)。未正規化の Y=ρY_k/ρ がこのノイズを持ち込み、エントロピー変数の s⁰ 項 (×~10) と ln X (対数微分特異) が w 空間ノイズに増幅、matrix 散逸の弱減衰 (エントロピー/せん断) 方向に注入され続けて市松減衰がプラトー化していた (scalar は全モード強減衰のため隠蔽)。**修正 = カーネル内で Y を正規化 (Σ=1 強制)**。切り分けの一時 run (dbg_*: 面単位ビット比較で KEEP_d 散逸は同一と証明→pure/scalar 対照→ns=1 強制 bisect で経路特定→Y 正規化で根治確認) は削除済み。
