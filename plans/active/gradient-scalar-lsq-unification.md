# node のスカラー勾配を LSQ に統一する (k/ω・化学種・受動種・凝縮モーメント)

## メタ

- **area**: `gradient`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/gradient.md`](../../methods/gradient.md) (「境界寄与」node × 周期の継ぎ目、スカラー勾配)
  - [`methods/discretization.md`](../../methods/discretization.md) §7.3 / §7.3.1 (「LSQ は NS だけ、スカラーは GG」の記述を更新する)
- **related_plans**:
  - [`boundary-node-periodic-gradient-fix.md`](../accepted/boundary-node-periodic-gradient-fix.md) (前提。§4.1 の合併 LSQ 係数は「変数に依らない形」で作ってある。§5.1 #7 の継続課題 (1)(3)(4) を本 plan に移す)
  - [`boundary-node-rotational-periodic.md`](boundary-node-rotational-periodic.md) (回転周期。本 plan では扱わない)
- **created**: `2026-09-26`
- **owner**: `sano`

## 1. 目的

node の勾配は、NS の原始量 ($\rho, u, P, T$) だけが LSQ (`gradLSQ: 2`、事前計算係数) で、スカラー ($k,\omega$、化学種 $Y_s$、受動種 $\xi$、凝縮モーメント) は Green–Gauss (GG) である。
2026-09-26 のユーザ決定 (「そろえるべき」) により、node のスカラー勾配も同じ LSQ に揃える。

根拠 (前提 plan の測定):

- GG は float32 で部分和 $\Sigma\phi_f S_f$ が $\sim\phi/h$ の打ち消しになり、**壁半割面込みの CV では定数場が閉じない** (channel 実測、壁節点で約 $5.3\,\varepsilon|\phi|/h$、継ぎ目に依らない。前提 plan §5.1 #7(3))。LSQ は定数場で厳密に 0。
- NS の LSQ は継ぎ目の合併係数 (`lsqPre_mergePeriodic`) を持ち、係数は変数に依らないので、スカラーも**同じ係数を流用**できる (並進周期の継ぎ目の扱いが NS と自動的に揃う)。
- 同じ場に NS は LSQ、スカラーは GG と作用素が混在していると、境界の扱い (LSQ は内部隣接のみ、GG は owner 値の半割面) が変数ごとに違う。

## 2. スコープ

- **やる**: node のスカラー勾配 4 系統 (`ransGradient` の $k,\omega$、`speciesGradient` の $Y_s$、`passiveGradient` の $\xi$・凝縮モーメント) を、NS と同じ事前計算 LSQ 係数 (`cInt`) による gather に置き換える。周期 gather の登録、出力配列名は現行どおり。検証ハーネス (`case/09.Taylor-Green/_g0_lsq_seam/`) の GG 参照を LSQ 参照 (`lsq_merged_ref`) に切り替える。前提 plan §5.1 #7 (1)(3)(4) の引き取り。
- **やらない**: cell モード (GG のまま。ユーザ方針で cell は使わない)。回転周期 (別 plan)。軸対称×周期の継ぎ目合併 (前提 plan #7(5)。LSQ 化で半割面閉包の問題は消えるが、合併は NS と同じく片側のまま)。スカラーのリミタの変更。遷移モデル (スカラー勾配を持たない)。

## 3. 関連 docs と前提

- 現状の事実 (2026-09-26 調査、`solver_density_cuda/`):
  - LSQ 係数 `cInt[3*ilp]` は `cell_planes` CSR の incidence ごとに float32 で 3 個、`calcGradient_d_wrapper` 内の static ローカル (`calcGradient_d.cu:959`)。境界 incidence (`ip >= nNormalPlanes`) は 0 (`:822`)、stencil は実在の内部隣接ノードだけ (ghost・bvar・境界点なし)。M⁺ は倍精度、打ち切り 1e-2。
  - 適用は `lsqPreGrad_internal_d` (`:850-894`) のノード並列 gather (atomic なし)。
  - スカラー GG: $k,\omega$ は `calc_scalar_gradient_face_d` (`ransTransport_d.cu:16-75`、境界面は owner 値 k[ic0]・ω[ic0])、$Y_s$・$\xi$・モーメントは `species_gradient_d` (`speciesTransport_d.cu:641-677`、境界面は実質 owner 値)。軸対称は `A_planar`・`s*_planar`。
  - 使用先: $k,\omega$ 勾配は SST の $CD_{k\omega}$・$F_1$ (`ransSource_d.cu:344-350`、`rans_sst_blend_f1_d`)、`axisymMethod: 1` の半径方向ソース (`ransSource_d.cu:257-271`、既定 0・使用 run 0 件)、それらを通じた拡散係数 σ(F1)μt の面補間 (`scalarTransport_d.cu:116-124`)。移流 (1 次) と粘性流束は勾配を使わない。~~$CD_{k\omega}$・$F_1$ だけ~~ (codex plan M1 で訂正)。$Y_s$ は SLAU の面再構成 (`Yd_recon`、Venkat ψ_Y は nSpecies ≥ 2)。受動種は SLAU S3 の再構成 (`dPdx_recon`、無次元化 Venkat) と、それを通じて FCT の Pface。
  - 周期 gather: `periodicGradientGather_d_wrapper` に NS 18 本・divU・dY・受動種を登録 (k/ω は `ransGradient` 直後の専用 gather)。合併係数の部分和は、group 内で状態が同値 (root→member ミラー、スカラーもミラー済み) なら和で合併 LSQ になる。
  - 軸対称の NS LSQ は係数が純幾何 (planar LSQ) で、planar GG と同じ量を推定する。

## 4. 設計方針 (2026-09-26 `diagnostician` 確定、同日 codex plan (GO-with-changes C0/M7/m2) を全件採用して改訂)

**前提の事実**: 既定の設定で生きているスカラー勾配は $k,\omega$ だけ (`speciesFaceReconstruction` の既定 0、`solverConfig.hpp:110`、`main.cpp:1453`)。ただし $Y$・$\xi$・モーメントの opt-in 経路には case/16 に使用実績 (`run_0476` 5 種 + トレーサ SFR 2、`run_0482` 凝縮 S3) があるので、検証は縮小しない (codex M6)。

### 4.1 係数の共有

- `cInt` を static ローカルから外に出し、**読み取り専用 (const ポインタ) のアクセサ**で公開する。アクセサは構築時の (nCells, nInc, mesh アドレス) も返し、スカラー wrapper は一致を assert する。単一メッシュ・単一プロセス前提を明記。NS と同じ配列を使う (合併済み係数を含む)。係数の再計算・変数別の係数は作らない。

### 4.2 適用カーネル

- 汎用の多変数 LSQ gather カーネルを 1 本足す: 入力 `flow_float**` (N 変数)、出力 3N 配列。ノード並列 (atomic なし)、アキュムレータは最大 4 変数ずつのチャンク (CSR と `cInt` を読み直す)。REG・spill は `cuobjdump -res-usage` で記録。
- **差分形** $\nabla\phi_i=\sum_j c_{ij}(\phi_j-\phi_i)$ (NS と同じ)。定数場の勾配は厳密に 0 (前提 plan #7(3) の GG 壁閉包 5.3ε は消える)。
- 境界の扱いは NS と完全に同一: 内部隣接のみ、境界 incidence は係数 0 を掛けるのでなく `ip >= nNormalPlanes` で **skip** (NS `calcGradient_d.cu:869` と同じ)、疑似点なし。ghost 出力のゼロ初期化の契約は維持 (`speciesTransport_d.cu:1265-1268`)。根拠: accepted の node 境界勾配の原則 (DOF-only) と同じで、変数ごとに境界の扱いが違う現状の方が不整合。GG の内部面値が使う幾何 fx (高 AR 曲面壁で 0.07–0.96 に振れる既知欠陥) も使わない。
- 軸対称: 係数は planar LSQ のまま (NS と同じ)。GG の `A_planar` 除算は不要。
- 壁節点の $k,\omega$: 低 Re + `sstNodeWallKPin` (既定 ON) では壁ノードは k=0・ω=ω_w にピン (`ransBoundary_d.cu:63-79`) され、残差・対角も 0 (`ransSource_d.cu:303-307`)。**壁ノードの F1 は k=0 により arg1_a = arg1_c = 0 で厳密に 0** (`ransSource_d.cu:343-350`、勾配作用素に依らない。~~1 に張り付く~~ は codex plan M1 で訂正)。作用素の差が入るのは (i) 第一内層ノードの $CD_{k\omega}$・$F_1$、(ii) それを通じた W–I 面の σ(F1) 補間、(iii) `axisymMethod: 1` の半径方向ソース (既定 0、使用 run 0 件)。S2 case/48 で第一内層 F1 の L∞ 差を記録する。
- **起動時に roK/roOmega も root→member ミラーする** (`periodicMirrorScalarState` を `main.cpp:1257` の直後に追加。現状の初期ミラーは NS・化学種・roXi のみで、roK/roOmega は更新後 `main.cpp:1722,1947` だけ)。同期入力ではビット不変 (S1 で確認)。

### 4.3 周期

- lsq 経路のスカラー gather (和 → broadcast) は係数合併と同じ述語 `periodicSeamMergeActive` (`periodicNode_d.cu:215-224`: node ∧ 並進 type 0 ∧ 非軸対称) を条件にする。k/ω の専用 gather (`ransTransport_d.cu:246`) は既にこの述語。回転周期・軸対称×周期は片側 LSQ (NS と同じ)。GG 経路の登録 (`periodicGradientGather`) は現行どおり。NS gather の述語統一は回転 plan #0a。
- 壁∩継ぎ目の ω_w (`wall_y_eff` は部分 stencil の最短距離、`ransBoundary_d.cu:153-203`) の group 同値性は **本 plan の必要条件ではない既存問題** (BC 値の問題で作用素と無関係、GG も同じ露出。整合格子では一致する見込み)。S0-e の診断で確認し、非零なら §5.1 に F 項目 (`wall_y_eff` を group-min に)。
- node + 回転周期 (type 1) は起動時に警告 (`[config] node 回転周期は未対応 (plan boundary-node-rotational-periodic)`)。エラー化は回転周期 plan の #0a。

### 4.4 切り替えと provenance

- **2 段ゲート**: Phase 1 = opt-in (`mesh.scalarGradient: gg|lsq`、既定 `gg`) の実装と S0–S3 → codex result 1 回目 → Phase 2 = 既定を `lsq` に切り替え (S4、設計 DB、docs) → codex result 2 回目。**既定化 (#6) の前提**: Phase 1 の result が GO、`tooling-stage-manifest-launch-binding` (別 plan、#2g) の完了、S3 ≤ 5 %。cell は GG 固定。
  opt-in から入る理由: 共有ワークツリーで並行セッションの node SST run と設計 DB が検証途中で黙って変わるのを避ける。
- 起動エコー `'scalarGradient' effective: <値> (default|explicit)`、`RUN_PROVENANCE` と起動記録 (`forge_launches.jsonl`) に実効値。
- **既定化の provenance** (codex M7: 現行 `stage_manifest.py:148-157` は同じ config hash に最後の起動値だけを対応させ、旧既定と新既定を誤連結する — codex が実証): 起動記録を (cfg_fnv, 起動順) で段に対応させ、hard キーは実効 `scalarGradient` + バイナリ id (git rev)。旧 manifest のキー無しは `gg (legacy)` とし `inferred` と区別して連結しない。修正は別 plan `tooling-stage-manifest-launch-binding` で行う (`slauWallNormalChi` の既定化にも同じ欠陥が効いている)。設計 DB: `runner_sern.py` の `FLAG_POLICY` を既定化日に更新し、行に実効 `scalarGradient` を持つ。`procedures/recommended-settings.md` §9 に「〜既定化の日付までは gg」。

## 5. 実装ステップ

1. `methods/gradient.md`・`methods/discretization.md` §7.3 を更新 (opt-in の LSQ 経路、既知の制約) — **済 2026-09-26**。
2. Phase 1: `cInt` の公開、汎用 LSQ gather カーネル、3 wrapper の分岐、roK/roOmega の初期ミラー、キー・エコー・provenance、回転周期の起動警告。
3. Phase 1: ハーネス拡張 (#2d) と S0/S1、起点の固定 (#2e) と S2・S3 (AWS)。codex result 1 回目。
4. Phase 2: 別 plan #2g の完了を待って既定を `lsq` に、S4、設計 DB、docs、codex result 2 回目。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 (**完了 2026-09-26**) | §4・§6 の確定 | 判断: 2026-09-26 `diagnostician`・§4 案 1–3 採用、周期は並進のみ合併、既定は opt-in から開始、§6 の合否を固定 | F |
| 2 (**完了 2026-09-26**) | codex plan 段 | GO-with-changes C0/M7/m2 → 判断: 2026-09-26 `diagnostician`・全件採用 (M2 の壁 ω 同値性は既存問題として診断項目に格下げ)、分割せず 2 段ゲート | F |
| 2a | M1 訂正の反映 | §3・§4.2・S2 表 (本改訂で済)、`axisymMethod: 1` の S0 変種 | O |
| 2b | roK/roOmega 初期ミラー | `main.cpp:1257` 直後に `periodicMirrorScalarState`。合格: R3 規則でビット不変 | O |
| 2c | S0-e 同値性診断 | channel と case/39 起点で 1 step、applyBconds 直後の k・ω・wall_y_eff の group 内差。非零なら `wall_y_eff` group-min を F 項目に | O (結論 F) |
| 2d | ハーネス拡張 | 非合併参照 (軸対称×周期)、ξ 定数場 (壁込み)、Y 5 種 (チャンク境界)、負の対照、gather 前配列のダンプ | O |
| 2e | 起点の固定 | §6 の起点表 (所在・sha256・バイナリ・実効設定) を AWS に転送し sha256 照合 | O |
| 2f | S2 収束規則のスクリプト | plateau 許容 + 床比 (末尾平均 ≤ 起点の 1.5 倍、ピーク除外・再進入 step 記録) + STEADY 閾値 | O |
| 2g | **別 plan 起票** `tooling-stage-manifest-launch-binding` | 起動順対応・バイナリ id・legacy 区別・S4 試験。本 plan #6 の前提 | F (§4) / O (実装) |
| 2h | 設計 DB | `runner_sern.py` FLAG_POLICY 更新と実効 `scalarGradient` 列 (#6 と同時) | O |
| 2i | S3 の測定手順 | native・同一 GPU・同一 BLOCKSIZE・ウォームアップ 500 後 2000 step × 3 の中央値、REG/spill | O |
| 3 | 実装 (Phase 1) | §4、§5 の 2 | O |
| 4 | S0/S1 | ハーネス (AWS) | O |
| 5 | S2/S3 | §6 の表 (**すべて AWS**、ユーザ指示で 2D も AWS) | O (結論 F) |
| 5r | codex result 1 回目 | Phase 1 | F |
| 6 | 既定の切り替え (Phase 2) | **前提**: 5r が GO、2g 完了、S3 ≤ 5 %。S4、2h、docs、codex result 2 回目 | F |
| 7 (**完了 2026-09-26**: `boundary-node-rotational-periodic` §5.1 #0a) | 回転周期 plan への登録 | | O |

## 6. 検証 (2026-09-26 `diagnostician` 確定、同日 codex plan 対応で改訂。測る前に固定)

- **S0 作用素** (ハーネス `case/09.Taylor-Green/_g0_lsq_seam/g_suite.py` の 7 変種 + 軸対称×周期 1 本 + `axisymMethod: 1` 1 本):
  - **S0-a 純作用素 (BC 前)**: k・ω・ξ・Y (5 種、チャンク境界 N=5 を跨ぐ) の LSQ 勾配が double 参照と ≤ 1e-5·S。参照は並進周期で `lsq_merged_ref`、軸対称×周期は**非合併** (root = identity) の同関数。非線形場 (sin) で識別。
  - **S0-b 定数場**: ξ (channel は入口が無くピン無し) で全節点・壁節点込みで勾配 == 0 (厳密)。k/ω は非ピン節点で == 0、壁隣接節点はピン値 (k=0, ω_w) からの差分形の解析値と ≤ 4ε。
  - **S0-c NS との一致**: 同じ場を ro に入れた NS 勾配と、**周期 gather 前の局所配列がビット同一**。gather 後は 2 member group がビット同一、3 member 以上 (多重周期の辺・角、atomicAdd の順序が非決定) は |差| ≤ 4ε·Σ|部分和|。
  - **S0-d 検出力 (負の対照)**: jitter32 で GG 参照との差が 1e-5·S を超える (閾値が GG と LSQ を識別する証拠)。
  - **S0-e BC 後の同値性 (診断)**: channel と case/39 起点で 1 step、applyBconds 直後の k・ω・wall_y_eff の group 内差 == 0 (記録。非零なら F 項目)。
- **S1 非干渉**: `scalarGradient: gg` で実装前バイナリと全配列が R3 規則で一致 (旧経路の保存、roK/roOmega 初期ミラー追加も含めて)。`lsq` で NS の勾配・リミタ配列が `gg` と同じ step でビット一致 (1 step、初回ダンプ)。
- **S2 物理 A/B — 共通規則**: 起点の最終場から `restart_field.py` で gg/lsq の 2 本、同じバイナリ・同じ step 数。収束は (i) `check_convergence` が DIVERGED でなく RISING 列が無いこと (plateau は既知床として可、先例 case/48 `run_0025`)、(ii) 各列の末尾平均が起点の末尾平均の ≤ 1.5 倍 (ピークは除外し、再進入 step を記録)、(iii) 比較量が `check_quasisteady` で STEADY (閾値は表)。~~`--from-floor`~~ は使わない (codex M4: 全期間ピークを見るので作用素切替直後の跳ねで落ち、case/48・case/40 の起点は plateau で参照側が REFUSED)。ただし case/16 `run_0476` は起点が通常判定 PASS なので `--from-floor` を使える。未達なら同一設定で延長 (初回予算は表の step 数)。**物理ゲートは独立の必須条件**。差が上限を超えたら自動不合格にも自動合格にもせず、格子対診断 (両方向・両経路) を回して F 判断 (codex M5: 縮むかだけでは欠陥と離散化差を識別できない)。

  | ケース | 起点 (所在・sha256 先頭) | step | 必須ゲート (既存) | 差 lsq−gg の上限 (超えたら格子対 → F) | STEADY 閾値 |
  | --- | --- | --- | --- | --- | --- |
  | case/48 冷却平板 | main ワークツリー `/home/sano/work/forge/case/48.flat_plate_cooled_m4/run_0005_B_tw300_y3/res_48000.h5` `9a9022d1…`、mesh `10a13dd0…` | 24000 | Cf/VD-II ∈ [0.96, 1.00] @x 0.3/0.6/0.9、2St/Cf ∈ [1.13, 1.19]、エネルギー閉合 ∈ [1.00, 1.03]、NaN 0。壁ノード F1 = 0 (両経路、記録)、第一内層 F1 の L∞ 差 (記録) | Cf・q_w・δ*・θ 各 ≤ 1 % (3 station)。格子対: 流れ方向 `fp_y1_3um_nx1500` (1.5×、期待比 ≤ 0.75)・壁方向 y1 6 µm (2×、期待比 ≤ 0.6)。差が反復床 (1.8e-6 相対) の 10 倍以上のときだけ縮小率を評価 | `--drift 0.002 --osc 0.005 --tail 0.4` |
  | case/40 軸対称ノズル (axisymMethod 0) | main ワークツリー `/home/sano/work/forge/case/40.nozzle_design_tool/run_0045_node_yp1_outletfix_cont/res_12000.h5` `15e47af0…`、nozzle.h5 `849660af…` | 12000 | η_CF ∈ 0.978 ± 0.003 (README:172、抽出は `design/forge_design/metrics/extract.py` の `thrust_metrics` → `eta_cf`)、ṁ ∈ 1.29–1.30 kg/s、NaN 0、出口列に負値なし | η_CF ≤ 0.1 %、壁温 L∞ ≤ 15 K、ṁ ≤ 0.2 % | `--drift 0.0005 --osc 0.002` (η_CF, ṁ) |
  | case/39 周期丘 | sern ワークツリー `case/39.periodic_hills/run_0039_r1_gradfix_new_ext` の最終 res (AWS、`--from-floor` PASS 済) | 200k | 継ぎ目比 r_gradk・r_gradw ∈ [0.9, 1.1]、NaN 0、S0-e で同値 | Cf 3 点・x_r は記録のみ。dF1_inf: lsq で STEADY なら「**case/39 のこの構成で**の GG 経路の性質」と限定して前提 plan #7(1) を閉じる、DRIFTING なら継続 | `--drift 0.002 --osc 0.005 --tail 0.4` |
  | case/16 `run_0476` 双子 (Y 5 種 + ξ、SFR 2、Euler) | `run_0476/res_24000` | 24000 + 継続 24000 | 継続区間で `--from-floor run_0476` PASS、0 ≤ Y、\|ΣY−1\| ≤ 1e-6、0 ≤ ξ ≤ 1、floor 補正累積 ≤ 1e-6 | vs gg 固定点: Y0 max\|Δ\| ≤ 3.7e-5、Xi ≤ 9.4e-4 (S3 信号の 10 %。反復床は 2.0e-7 / 1.5e-5) | `--quantity machmax,pmax` STEADY |
  | case/16 `run_0482` 双子 (凝縮 S3、SST TP、モーメント 4) | `run_0470/lim1e` 入力 + IC `run_0213` (0482 と同一) | 48000 | `check_convergence` PASS (0482 と同基準)、ALL STEADY、モーメント ≥ 0、NaN 0 | vs 0482: onset (中心線 g = 1e-3 の線形補間交差、格子間隔 0.426 mm) ≤ 0.3 mm (S3 効果 0.72 mm の 4 割)、壁 p/p0 L∞ ≤ 0.3 % (`extract_wall_pp0.py` の CSV、S3 効果 0.81 % 比)。「実験に近づく」は判定に使わない | `compare_condfix --series-csv` ALL STEADY |
  | dual-time FCT smoke (`run_0486` プロトコル + `run_0476` 入力のトレーサ版) | `run_0471/res_12000`、`run_0476/res_24000` | 200 | NaN 0、\|ΣY−1\| ≤ 1e-7、0 ≤ ξ ≤ 1 | vs gg 双子 (同一物理時刻): ro ≤ 3e-6、roY ≤ 5e-6、roUy ≤ 1e-4 (反復差 0486 vs 0487 の 2 倍) | — |

- **S3 性能**: case/48・case/39 で lsq/gg の step 時間比 ≤ 1.05 (AWS native、同一 GPU・BLOCKSIZE、ウォームアップ 500 step 後 2000 step × 3 反復の中央値)。**超過は #6 を保留**。REG・spill・実測を記録。
- **S4 provenance (#6 の前提、別 plan #2g で実施し結果を写す)**: 同一 YAML を旧バイナリ (既定 gg) → 新バイナリ (既定 lsq) で 2 段起動し `stage_manifest.py --segments` が **2 区間**を返す。旧省略 / 新省略 / 明示 gg / 明示 lsq / cell (実効 gg) の 5 通りで実効値と由来 (launch / explicit / legacy / inferred) が正しい。設計 DB で FLAG_POLICY 旧行が選別から外れる。
- **前提 plan #7(4) (ジッタ格子の LSQ 次数)**: ~~本 plan で「欠陥ではなく作用素の性質」として閉じる~~ (codex M3 で訂正) → 本 plan では閉じない。スカラーは NS と同じ LSQ の次数の性質を継承する (S0-c のビット同一が根拠)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-unification-plan.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-unification-plan.md) | **GO-with-changes**, C0/M7/m2 (差分形 LSQ の係数共有は支持) | **全件採用** (2026-09-26 `diagnostician`)。M1 → 壁 F1 = 0 と axisymMethod 1・σ 補間の訂正 (§3・§4.2)。M2 → roK/roOmega 初期ミラー、スカラー gather の述語を `periodicSeamMergeActive` に統一 (§4.3)、壁 ω の同値性は既存問題として S0-e の診断に格下げ。M3 → S0 を a–e に分割。M4 → `--from-floor` を捨て plateau 許容 + 床比 + STEADY、起点を main ワークツリーで確認し sha256 固定。M5 → 物理ゲートは独立の必須条件、格子対は診断で F 判断、STEADY 閾値を締める。M6 → case/16 の 2 双子 + FCT smoke を追加。M7 → 既定化の provenance を別 plan (#2g) に切り出し #6 の前提に。m8 → アクセサの契約 (§4.1)。m9 → S3 を #6 の保留条件に。分割はしない (2 段ゲート) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu`、`ransTransport_d.cu`、`speciesTransport_d.cu`、`periodicNode_d.cu`、`input/solverConfig.{hpp,cpp}`。
- **既存の node SST・化学種再構成・受動種 run の結果が変わる** (勾配作用素の変更)。

## 8. 完了条件

- [ ] `methods/` 更新
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — codex plan (GO-with-changes C0/M7/m2) を全件採用して §3–§6 を改訂 (`diagnostician`)。前回の「壁 F1 は 1 に張り付く」は誤りで訂正 (k=0 ピンで F1=0)。2 段ゲート (Phase 1 opt-in → Phase 2 既定化) に。

- `2026-09-26` — §4・§6 を `diagnostician` 判断で確定 (opt-in から開始、S0–S2 通過後に本 plan 内で既定化、既存ゲート + 格子対の合否規則)。
- `2026-09-26` — 起票 (ユーザ決定「スカラー勾配も LSQ に揃える」、前提 plan の accepted を受けて)。
