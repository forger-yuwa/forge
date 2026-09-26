# node のスカラー勾配を LSQ に統一する (k/ω・化学種・受動種・凝縮モーメント)

## メタ

- **area**: `gradient`
- **status**: `in_progress` (Phase 1 合格・Phase 2 既定化の実装と確認済み 2026-09-27、codex result (Phase 2) 待ち)
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
- 境界の扱いは NS と完全に同一: 内部隣接のみ、境界 incidence は係数 0 を掛けるのでなく `ip >= nNormalPlanes` で **skip** (NS `calcGradient_d.cu:869` と同じ)、疑似点なし。ghost 出力のゼロ初期化の契約は維持 (`speciesTransport_d.cu:1265-1268`)。化学種 dY の ghost は両経路とも書かない (GG も `ic1 < nCells` のみ、`:662`) ので未定義のまま。根拠: accepted の node 境界勾配の原則 (DOF-only) と同じで、変数ごとに境界の扱いが違う現状の方が不整合。GG の内部面値が使う幾何 fx (高 AR 曲面壁で 0.07–0.96 に振れる既知欠陥) も使わない。
- 軸対称: 係数は planar LSQ のまま (NS と同じ)。GG の `A_planar` 除算は不要。
- 壁節点の $k,\omega$: 低 Re + `sstNodeWallKPin` (既定 ON) では壁ノードは k=0・ω=ω_w にピン (`ransBoundary_d.cu:63-79`) され、残差・対角も 0 (`ransSource_d.cu:303-307`)。**壁ノードの F1 は k=0 により arg1_a = arg1_c = 0 で厳密に 0** (`ransSource_d.cu:343-350`、勾配作用素に依らない。~~1 に張り付く~~ は codex plan M1 で訂正)。作用素の差が入るのは (i) 第一内層ノードの $CD_{k\omega}$・$F_1$、(ii) それを通じた W–I 面の σ(F1) 補間、(iii) `axisymMethod: 1` の半径方向ソース (既定 0、使用 run 0 件)。S2 case/48 で第一内層 F1 の L∞ 差を記録する。
- **起動時に roK/roOmega も root→member ミラーする** (`periodicMirrorScalarState` を `main.cpp:1257` の直後に追加。現状の初期ミラーは NS・化学種・roXi のみで、roK/roOmega は更新後 `main.cpp:1722,1947` だけ)。同期入力ではビット不変 (S1 で確認)。

### 4.3 周期

- lsq 経路のスカラー gather (和 → broadcast) は係数合併と同じ述語 `periodicSeamMergeActive` (`periodicNode_d.cu:215-224`: node ∧ 並進 type 0 ∧ 非軸対称) を条件にする。k/ω の専用 gather (`ransTransport_d.cu:246`) は既にこの述語。回転周期・軸対称×周期は片側 LSQ (NS と同じ)。GG 経路の登録 (`periodicGradientGather`) は現行どおり。**lsq では `periodicGradientGather` に dY・dξ・モーメント勾配を登録しない** (`periodicNode_d.cu:186-188`、2026-09-26 実装時に追加・`diagnostician` 採用: 登録すると wrapper 内の合併と二重になり継ぎ目が 2 倍、また同 gather は回転周期でも合併してしまう)。NS gather の述語統一は回転 plan #0a。
- 壁∩継ぎ目の ω_w (`wall_y_eff` は部分 stencil の最短距離、`ransBoundary_d.cu:153-203`) の group 同値性は **本 plan の必要条件ではない既存問題** (BC 値の問題で作用素と無関係、GG も同じ露出。整合格子では一致する見込み)。S0-e の診断で確認し、非零なら §5.1 に F 項目 (`wall_y_eff` を group-min に)。
- node + 回転周期 (type 1) は起動時に警告 (`[config] node 回転周期は未対応 (plan boundary-node-rotational-periodic)`)。エラー化は回転周期 plan の #0a。

### 4.4 切り替えと provenance

- **2 段ゲート**: Phase 1 = opt-in (`mesh.scalarGradient: gg|lsq`、既定 `gg`) の実装と S0–S3 → codex result 1 回目 → Phase 2 = 既定を `lsq` に切り替え (S4、設計 DB、docs) → codex result 2 回目。**既定化 (#6) の前提**: Phase 1 の result が GO、`tooling-stage-manifest-launch-binding` (別 plan、#2g) の完了、S3 ≤ 5 %。cell は GG 固定。
  **未合格項目の移管は Phase 1 GO を意味しない** (2026-09-26 codex (diagnose) 5): 床移動の真因調査を調査メモへ移しても、S2 の収束ゲート未達と FCT smoke の未合格は Phase 1 の未達として残る。
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
| 2a (**完了 2026-09-26**: §3・§4.2・S2 表を改訂済、`axisymMethod: 1` の S0 変種は `S0_axi_m1.txt`) | M1 訂正の反映 | §3・§4.2・S2 表 (本改訂で済)、`axisymMethod: 1` の S0 変種 | O |
| 2b (**実装済 2026-09-26**、合格判定は S1) | roK/roOmega 初期ミラー | `main.cpp:1257` 直後に `periodicMirrorScalarState`。合格: R3 規則でビット不変 | O |
| 2c (**完了 2026-09-26**: channel・case/39 とも group 内差 0、`S0_channel.txt`・`S0Y_channel.txt`・`S0e_case39.txt`) | S0-e 同値性診断 | channel と case/39 起点で 1 step、applyBconds 直後の k・ω・wall_y_eff の group 内差。非零なら `wall_y_eff` group-min を F 項目に | O (結論 F) |
| 2d (**完了 2026-09-26**: `g_suite.py`・`s0y_species.py`・`pregather_check.py`) | ハーネス拡張 | 非合併参照 (軸対称×周期)、ξ 定数場 (壁込み)、Y 5 種 (チャンク境界)、負の対照、gather 前配列のダンプ | O |
| 2e (**完了 2026-09-26**: 起点 5 件の sha256 を AWS で照合、各 run の `IC_FROM.txt`) | 起点の固定 | §6 の起点表 (所在・sha256・バイナリ・実効設定) を AWS に転送し sha256 照合 | O |
| 2f (**実装済 2026-09-26**: `solver_density_cuda/tools/check_floor_ratio.py --start <起点 run> <run>...`、`--selftest` 6 通り PASS。比較量の STEADY は従来どおり `check_quasisteady.py`) | S2 収束規則のスクリプト | plateau 許容 + 床比 (末尾平均 ≤ 起点の 1.5 倍、ピーク除外・再進入 step 記録) + STEADY 閾値 | O |
| 2g (**起票済 2026-09-27**: [`tooling-stage-manifest-launch-binding.md`](tooling-stage-manifest-launch-binding.md)、draft) | **別 plan 起票** `tooling-stage-manifest-launch-binding` | 起動順対応・バイナリ id・legacy 区別・S4 試験。本 plan #6 の前提 | F (§4) / O (実装) |
| 2h (**完了 2026-09-27**: `runner_sern.py` の `FLAG_POLICY` を 2026-09-27、行に実効 `scalar_gradient_effective`、`driver_sern.py` の学習は全作動点 lsq の行だけ (codex diagnose 2026-09-27 の Major を採用)、`design/tests/run_sern_scalar_gradient_gate_tests.py` 5 例 PASS = 判別 A) | 設計 DB | `runner_sern.py` FLAG_POLICY 更新と実効 `scalarGradient` 列 (#6 と同時) | O |
| 2i (**完了 2026-09-26: S3 PASS** — 現行の値は PID 照合版の取り直し (#5i): case/48 0.997・case/39 1.018。以下は旧監視版の記録: lsq/gg 比 case/48 1.010・case/39 1.017 (上限 1.05)。AWS g5・`f99f236d`・block 128・2500 step の step 501–2499 平均・交互 3 反復の中央値 (`s3_perf.py`、`S3_case{48,39}.txt`)。case/48 の初回は別 run と重なったので無効 (`S3_case48_contaminated.txt`) にして走行中監視つきで取り直した。「native」は AWS の native ビルド) | S3 の測定手順 | native・同一 GPU・同一 BLOCKSIZE・ウォームアップ 500 後 2000 step × 3 の中央値、REG/spill | O |
| 3 (**完了 2026-09-26**、判断: 2026-09-26 `diagnostician`・diff レビューで欠陥なし、periodicGradientGather の登録変更を採用) | 実装 (Phase 1) | §4、§5 の 2。REG 38/39/40/48 (NV 1–4)、spill 0。AWS 最小確認: lsq の線形場誤差 ≤ 丸め床、NS 配列は lsq/gg でビット一致、gg の面寄与ダンプは HEAD と不一致 0 | O |
| 4 (**完了 2026-09-26: S0/S1 PASS**。S0 は全変種 PASS (軸対称 2 変種の S0-c は §6 の適用除外)。S1 (1) 4 ケース PASS。S1 (2) は §6 の訂正を適用: tgv の NS 18 勾配は gather 前で gg/lsq 12 本ビット一致、group 外節点は厳密一致、2 member 以上は全点が順列和の集合内 (#4d・#4f)。初回 FAIL の記録は §6.2 相当として `S1.txt` に残す) | S0/S1 | ハーネス (AWS)。結果 `case/09.Taylor-Green/_g0_lsq_seam/{S0_*.txt,S0e_case39.txt,S1.txt}` | O |
| 4a (**完了 2026-09-26**: Minor 2 件を修正して commit (v2 バイナリ `f1703bf8…` で pre-gather ダンプ 4 tag が v1 とバイト一致、E/P = A・N PASS、`~/sglsq/DUMP_RECHECK_tgv_lsq_v2.txt`)。**commit 可 — 判断: 2026-09-26 codex (diagnose) `notes/reviews/2026-09-26-gradient-scalar-lsq-4a-closeout-diagnose.md`**。Minor 2 件 (null 成分をゼロで埋めず診断失敗として tag を捨てる、env 無しでは呼び出し側で vector を組まない) を直してから commit。S1(1) 再実行 4 ケース PASS。「ダンプ非干渉検証済み」とは記録しない — 判定は #4d) | 出力 hook (判断: 2026-09-26 `diagnostician`、既定 off・出力専用・数値不変) | (a) `FORGE_DUMP_PREGATHER=<path>`: `periodicGradientGather` の直前 (`main.cpp:1304`・`:1477`) に NS 18 本 + dY + dξ/モーメントの局所配列 [nVar][nCells][3] を、`ransGradient` の lsq 分岐では gather 直前の dK/dΩ を非 atomic に書く。gg 経路でも同じ位置で書ける。既存 `FORGE_DUMP_SCALARGRAD` (面寄与) とは別名。(b) `dY{s}d{x,y,z}` と `wall_y_eff` を output の extraFields に登録 (post-gather 値)。**入れた commit で S1(1) (gg 旧 vs 新、4 ケース × 3 本) を再実行して不変を示し、ダンプ有効時に res が変わらないことを 1 ケースで見る** | O |
| 4b (**完了 2026-09-26: A**。gather 前 NS 18 配列は 6 本ビット同一、gather 後の差 (res_0 124・res_1 232 節点、全て 4 member) は全点が部分和の順列和の集合内。`PREGATHER_4b.txt`) | tgv S1(2) FAIL の判定 A/B (測る前に固定) | gg 3 本・lsq 3 本で pre-gather ダンプを取り、(i) NS 18 配列の gather 前がビット同一か、(ii) 不一致 36 節点の gg・lsq の gather 後の値が member 部分和の float32 順列和の集合に含まれるか。**A**: (i) 同一かつ (ii) 全点で含まれる → 「atomicAdd の順序差」として閉じ、§6 S1 (2) を下記に訂正。**B**: (i) 不一致 → lsq 経路が NS の入力を壊している → 実装を止めて調べる。**B′**: (i) 同一だが (ii) が外れる → gather 以外の非決定源があるので順序差の説明は採らない | O (結論 F) |
| 4c (**完了 2026-09-26: 4 変種 PASS**。`S0Y_*.txt`、`S0e_case39.txt`) | Y の未測定項目 (S0-a の Y 5 種、S0-c の dY と dξ のビット一致、チャンク境界 4+1、S0-e の wall_y_eff) | 5 種は case/16 `run_0471` の組 [H2O, N2, O2, AR, CO2] (TP、README:30 で PASS 済み)。Y0 = Y4 = q (ξ と同じ量子化場、q < 0.5)、Y1..Y3 は正の sin 場で Σ = 1 − 2q を double で作ってから量子化。Y4 = ξ でチャンク境界の孤立変数を、Y0 で先頭チャンクを見る。S0-a の参照入力は res_0 の読み戻し (`VALUE/Y{s}`・`Xi`)。追加: Σ_s dY_s = 0 (≤ 4ε Σ\|dY_s\|) を全節点で | O |
| 4d (**完了 2026-09-26: E/P = A、N = FAIL 1 配列 (tgv lsq res_1 `res_ro`: 群内 11・14、交差 40、最大差 3 組とも 2.98e-8)。ダンプ全体の非干渉は未確認**。case48 は E/P/N すべて通過。`DUMP_RECHECK_{tgv,case48}.txt`。N の扱いは #4f) | ダンプ非干渉と S1 (2) の再判定 (**事後規則**、判断: 2026-09-26 codex (diagnose)。§6「ダンプ非干渉」) | 追加計算なし、既存の 1 step 出力を §6 の 3 区分 (E 厳密 / P 順列集合 / N 両側ノイズ) で再判定。tgv: ダンプ無し `~/sglsq/s1h/s1_tgv_{gg,lsq}_{a,b,c}` とダンプ有り `~/sglsq/pg/pg_tgv_{gg,lsq}_{a,b,c}` の全 12 本の res_0・res_1 の NS 勾配を、ダンプ有り gg_a の pre-gather 部分和に照合 (非周期節点は部分和そのもの、2 member は唯一の和、3 member 以上は順列集合)。case48: 無し `s1h/s1_case48_gg_*` vs 有り `pg/dc_case48_*`。N 区分は不一致数と最大差を 3 組 (無し同士・有り同士・無し–有り) とも出す。**A**: E が全て厳密一致かつ P が全点集合内 → S1 (2) 訂正を適用して tgv PASS、#4a の非干渉を記録。**B**: E の不一致または P の集合外 → 順序差だけの説明を棄却し S2 を保留 (F)。結果 `DUMP_RECHECK.txt` | O (結論 F) |
| 4f (**完了 2026-09-26: A** — E/P 16 本厳密・集合内、N 全配列 2 倍以内 (`res_ro`: 無し同士 49・有り同士 42・交差 47)。記録は「固定 8 本の規則では超過しない」まで。初回 FAIL は #4d に保存。`DUMP_RECHECK_tgv_lsq8.txt`、run `~/sglsq/pg8/lsq_{off,on}_{d..h}` (AWS スクラッチ、破棄可)。測る前に固定、判断: codex (diagnose) `notes/reviews/2026-09-26-gradient-scalar-lsq-4d-nfail-diagnose.md`) | #4d の N FAIL (tgv lsq res_1 `res_ro`) の一度限りの追加試験 | 変更点は `FORGE_DUMP_PREGATHER` の無効/有効だけ。同一バイナリ (`~/sglsq/forge_hook4a` sha `a398afdc…`、既存 12 本と同じ)・同一入力 (`pg_tgv_lsq_a` の入力を複製、mesh.h5 の sha256 を記録)・同一 GPU・BLOCKSIZE 128/128、1 step。lsq の無し・有りを各 +5 本 (d–h)、実行順は ABBA で固定: d無 d有 e有 e無 f無 f有 g有 g無 h無 h有。計 8 本ずつで、**res_0・res_1 の全 N 配列**について群内 28 対 + 28 対の最大と交差 64 対の最大を比べ、最大差・不一致数が**両方とも** 2 倍以下。E/P も 16 本で再確認し、各 run の sha256・実効 scalarGradient・mesh sha256 を記録。**A**: E/P 通過・全 N 通過 → 「固定 8 本の規則では超過しない」、初回 FAIL を持続的干渉の証拠とする解釈を退ける (初回 FAIL は保存)。**B**: E/P 通過・N FAIL 継続 → ダンプ非干渉は未確認のまま (#4a は診断専用)。**本数をさらに増やして合格を探さない**。どちらでも S2 は進める (ダンプ無し)。E/P 失敗 → S2 保留。gg の不一致数を借りて合格にしない。結果 `DUMP_RECHECK_tgv_lsq8.txt` | O (結論 F) |
| 4e | S0-c 軸対称 2 変種の適用除外 (判断: 2026-09-26 codex (diagnose) の指摘を受けて理由を記録) | §6 S0-c に除外理由を書く (済)。測定は追加しない | O |
| 5 (**S2 = PASS (ユーザ決定の基準差し替え後、2026-09-26)**: 現行 gg 双子を基準に lsq の末尾残差 5 ケースとも ≤ 1.5 倍 (最大 case/48 `rms_roe` 1.42、`S2_GGREF.txt`)、gg 双子は plateau で NaN/RISING なし、FCT smoke の収支は lsq・gg 同値 (remainder 8.45e-2)。物理ゲートと lsq−gg 上限は FCT smoke の精度差 (審査しない) 以外 PASS、S3 PASS。旧判定: S2 = 条件付き・未合格 2026-09-26、判断: codex (diagnose) `notes/reviews/2026-09-26-gradient-scalar-lsq-s2-result-diagnose.md`。物理ゲートと lsq−gg の上限は FCT smoke 以外 PASS、**起点床に対する収束ゲートは case/48・40・0476・0482 で gg・lsq とも不成立 (case/39 は両方 PASS)**、FCT smoke は旧上限 (反復差 × 2) 超過。旧判定は書き換えない。結果の表は各 case README の `run_095x_sglsq_*` 行) | S2/S3 | §6 の表 (**すべて AWS**、ユーザ指示で 2D も AWS)。**着手条件** (2026-09-26 codex (diagnose) 2 回): #4d・#4f の **E/P** が A、#2e (起点の sha256 照合) と #2f (収束規則スクリプト) の完了。N 区分の結果 (#4f の A/B) は S2 を保留する根拠にしない (S2 はダンプ無しで回す) | O (結論 F) |
| 5a (**完了 2026-09-26: B** — `run_0956_sglsq_s2_gg_chi0` (48000 step、実効 chi 0 を起動エコーで確認) も起点床に戻らない: rms_ro 1.82×・roUy 9.08×・roe 1.85× (chi 1 の gg は 1.85×/9.1×/1.88×)。物理量は STEADY、chi 0/1 の差は Cf・q_w・θ ≤ 0.007 %・δ* ≤ 0.085 %。**chi 単独で床移動を説明する仮説を棄却**、自動延長はしない。残る候補は起点 (09-12) 以降の他のバイナリ変更と restart そのもの。測る前に固定、判断: codex (diagnose) S2) | case/48 の chi 一因子試験 | `run_0956_sglsq_s2_gg_chi0`: `run_0950_sglsq_s2_gg` の開始時入力 (同じ IC = run_0005/res_48000 のビット一致コピー) を複製し `space.slauWallNormalChi: 0` だけ変更、同一バイナリ `f99f236d`・BLOCKSIZE 128、48000 step (対照 0950+0952 の合計)。判定: `check_floor_ratio --start run_0005 --factor 1.5 --tail 0.2` の全残差列、`check_convergence` の NaN/RISING、既存の物理ゲート (Cf/VD-II・2St/Cf・閉合) と系列 STEADY。**A**: 全条件を満たし旧床へ戻る → case/48 では床移動を chi 変更で説明する (他ケースへ一般化しない)。**B**: 物理量が STEADY でも床を外れる → chi 単独の説明を棄却。まだ過渡なら判定不能とし自動延長はしない | O (結論 F) |
| 5b (**完了 2026-09-26**: 継続末尾で Y0 1.63e-5・Xi 3.39e-4 (上限 3.7e-5/9.4e-4)、ΣY・ξ・floor 補正 (rel 1.3e-27) とも範囲内、machmax/pmax STEADY。`--from-floor run_0476` は両方不成立のまま) | 0476 双子の継続 24000 | §6 表どおり `run_0950/0951_sglsq_s2_0476_*` の res_24000 から同一設定で +24000 (`run_0956/0957_sglsq_s2_0476_*_cont`)。比較量・上限は §6 表のまま | O |
| 5c (**完了 2026-09-26: 上限内**。onset (中心線 g = 1e-3) 0482・gg・lsq とも 23.23 mm、壁 p/p0 L∞ gg vs 0482 0.165 %・lsq vs 0482 0.163 % (上限 0.3 %、最大は出口 x = 95 mm)、lsq−gg 0.003 %) | 0482 双子の vs 0482 比較 | §6 表の比較対象は **vs run_0482** (lsq−gg ではない): onset と壁 p/p0 を gg・lsq それぞれ run_0482 と比べ、lsq−gg も併記する。旧バイナリ差 (gg vs 0482) と分けて書く | O |
| 5d (**測定済 2026-09-26、判定は #5r**: floor 補正累積 rel 7e-19・limCorr 0・上限違反 0 セル・E_ρ 最大 5.5e-6 (gg) / 4.4e-6 (lsq)、トレーサ総量 200 step 後 4.1705e-5 で gg/lsq の差 rel 9e-8。流入流出の収支照合は未実施) | FCT smoke の保存・有界性 | 旧上限 (反復差 × 2) 超過の記録は残す。**0476 の定常差を上限に流用しない** (codex 却下)。受動種の保存収支 (floorCorr・limCorr・total の時系列) と有界性を独立に判定し、差は記録として #5r で審査する。`nSub 15` の十分性は未確認と明記 | O (結論 F) |
| 5e | `slauWallNormalChi` の文書不整合 | `procedures/recommended-settings.md` §1.0a は既定 0 と書くが `f99f236d` のコードは node+nodeWallDirichlet+SLAU で auto = 1 (`solverConfig.cpp:700`)。担当 plan (`convection-slau-wall-normal-chi-default`) 側で同期してもらう。本 plan の run は実効値を `IC_FROM.txt`・起動エコーで固定済み | O |
| 5r | codex result 1 回目 | Phase 1。**S2 は「条件付き・未合格」のまま審査に回してよい** (codex (diagnose) S2)。旧ゲート不成立・FCT 旧上限超過・#5a–#5d の結果を併記する | F |
| 5f (M1、採用。**取り直し済 2026-09-26**: 修正後ツールで S2 の全床比判定を再実行、判定は不変 (case/39 のみ両方 PASS、判定不能なし)、`FLOOR_RECHECK.txt`) | `check_floor_ratio.py` の誤合格 | 列単位の判定不能・起点で活動していた列の全ゼロ・起点系列の NaN/Inf/必須列欠落/床 ≤ 0 を不合格へ伝播 (**修正済 2026-09-26**、`--selftest` 9 例 PASS)。**修正後に S2 の全床比判定を取り直す** | O |
| 5g (**完了 2026-09-26: A**。床比 f_old/f_new 0.989–1.092 (全活動列)、物理量差 ≤ 0.023 % (x 0.6 の δ*、区間変動 ≤ 0.09 % で判別可)、旧版も系列 STEADY・RISING なし・起点床を外れる (rms_ro 1.83×・roUy 9.10×) → **この case・手順・閾値では本 plan の commit 群による大幅な追加変化を検出せず、旧床逸脱は実装直前版でも再現**。`M2_case48.txt`。規則は測る前に固定、判断: codex (diagnose) 4) | 実装直前版との比較 (case/48) | 旧版 `36d8ba03` (AWS `~/sglsq/forge_36d8ba03`、sha256 `29e8f590…`、`f99f236d` と同じ CXXFLAGS・arch 86) を、現行対照 (`run_0950`+`run_0952`) と**同じ 24000 + restart + 24000 step** で (`run_0957_sglsq_s2_base36d8` → `run_0958_sglsq_s2_base36d8_ext`)。IC・メッシュ・BC・実効設定・GPU・BLOCKSIZE を揃え、`scalarGradient` キーは旧版に無いので書かない。chi は両バイナリの起動エコーで実効値を確認。判定: 延長段の全活動残差列の末尾 20 % 平均 f で**各列 2/3 ≤ f_old/f_new ≤ 1.5**、Cf・q_w・δ*・θ (3 station) の固定末尾区間 (**延長段 step 12000–24000 のスナップショット 7 枚の平均**、2026-09-26 結果を見る前に固定、`m2_compare.py`) 平均で**旧版を分母に相対差 ≤ 0.1 %** (同区間の (max−min)/|mean| が 0.1 % を超える量は判定不能)、両者の系列 STEADY・既存物理ゲート・NaN/Inf/RISING なし。**A**: 全条件成立かつ旧版も起点床を外れる → 「この case・手順・閾値では大幅な追加変化を検出せず、旧床逸脱は実装直前版でも再現」(「床を変えていない」「原因は過去の commit」とは書かない)。**B**: 床比か物理差が上限超過 → commit 群の調査へ戻る。過渡・設定不一致・証拠不足は判定不能 (自動延長・閾値緩和はしない)。case/48 の結果だけで他ケースの収束ゲート未達は閉じない | O (結論 F) |
| 5h (**試験完了・ゲート未合格 2026-09-26: B**。(a) `check_passive_budget --mode fct` は gg・lsq とも FAIL (remainder 8.45e-2、閉合 3e-10・総量照合 4e-9 は通過)、nSub 30 でも同値で 4 本とも FAIL。(b) E_q/D_q ≈ 2 (roY0 のみ 0.2) で感度条件不成立 → **30 の双子を確定値にせず、FCT smoke の精度差の審査は保留**。remainder が作用素・nSub に依らず同値である原因は未特定。`M3_fct.txt`。規則は測る前に固定) | FCT smoke の収支と nSub 感度 | (a) `check_passive_budget.py --mode fct` を `run_0954/0955` に適用 (収支閉合・独立総量照合・低次/HO 残差)。(b) gg・lsq それぞれ `nSubIterDualTime` 15→30 だけを変え、同じ開始状態から同じ 200 物理 step。同じ物理時刻・節点・固定正規化で D_q = ‖q_lsq15 − q_gg15‖、E_q = ‖q_gg30 − q_gg15‖ + ‖q_lsq30 − q_lsq15‖ (q は ro・roY・roUy を含む)。**A**: 全 q で E_q < 0.1 D_q かつ 4 本とも保存・有界性ゲート通過 → 「15→30 の反復感度では差の大部分を説明できない」と記録 (作用素差を支持するが正しい離散化差とは確定しない)。**B**: 感度条件不成立 → 30 の双子を確定値にせず精度差の審査は保留。D_q = 0・ノイズ以下・記録不足・保存ゲート不成立は判定不能。旧上限超過は残し、A でも旧精度ゲートを自動合格にしない | O (結論 F) |
| 5i (M4、採用。**完了 2026-09-26**: PID 照合版で競合検出試験 (別名バイナリの run を測定中に割り込ませ検出・取り直し) を通し、S3 を取り直して case/48 0.997・case/39 1.018 で PASS (`S3_case{48,39}.txt`、`S3_v2_*.txt`)。0.5 s より短い競合は排除できない) | S3 の競合監視 | GPU 計算プロセスの PID を実際の solver 子プロセス (run_case.sh の bash ではない) と照合し、対象の観測と正常終了を必須に。競合検出試験を通してから S3 を 2 case とも取り直す。ポーリングで短時間の競合まで排除したとは書かない | O |
| 5j (M5、採用・条件付き暫定。**修正済 2026-09-26**: `YAML_HARD_PATHS` に `mesh.scalarGradient`、YAML 解析不能時は本文ハッシュを key に入れて連結しない、`test_stage_manifest_scalar_gradient.py` 4 例 PASS・既存 chi テスト PASS。本 plan の S2 run は単一段で manifest を使っていない) | stage_manifest の区間 | `YAML_HARD_PATHS` に `mesh.scalarGradient` を追加 (明示 gg↔lsq を別区間に)。既存 manifest の保存済み key は直らないので、既存記録は設定原本から再生成するか判定区間を明示。YAML 解析不能時に黙って連結しない。**#2g の完了が既定化の前提**のまま | O |
| 5k (m6、採用。**済 2026-09-26**: `methods/gradient.md`・`methods/discretization.md` §7.3・`procedures/solver-settings.md` に `mesh.scalarGradient` 節) | methods / procedures の記述 | `methods/discretization.md`・`methods/gradient.md` を「既定 gg / opt-in lsq (実装済み・検証未完了) / cell は GG 固定」、非合併条件の lsq は「片側 LSQ」に。`procedures/solver-settings.md` に `mesh.scalarGradient` | O |
| 5l (m7、採用。README 行は済 2026-09-26。**床移動の真因は調査メモ [`notes/investigations/2026-09-26-node-steady-floor-shift.md`](../../notes/investigations/2026-09-26-node-steady-floor-shift.md) へ移管** (codex (diagnose) 5: 移管は S2 の免除ではない、#5・§6 の未合格条件は残す)) | 索引と未解決項目 | case/48・case/16 README の追加 run 行と結果。**床移動の真因** (設定差のある 4 case で gg・lsq とも起点床に戻らない、chi 単独は棄却) を独立 F 項目に (別項目化で S2 未合格を解消した扱いにはしない) | F |
| 5m (**完了 2026-09-26: A** — 旧版 `36d8ba03` gg (`run_0960`、sha `29e8f590…`) も現行 gg キー省略 (`run_0961`、実効 gg (default) を確認) も `check_passive_budget --mode fct` が同じ remainder 8.45e-2 で FAIL (閉合 5e-10/−4e-10・総量照合 3.6e-9・有界性 [0, 0.9998] は同水準)。→ 本 plan の commit 群がこの FAIL を初めて生んだ仮説を棄却。lsq の精度ゲートは未合格のまま。測る前に固定、判断: codex (diagnose) 5 `notes/reviews/2026-09-26-gradient-scalar-lsq-phase1-status-diagnose.md`) | FCT smoke の実装直前版 / 現行版 gg | 変更因子はバイナリだけ (`36d8ba03` / `f99f236d`)。`run_0954_sglsq_s2_fct_gg` の開始入力を複製し、**両方とも `scalarGradient` キーを省略** (現行の実効 gg を起動エコーで確認)、nSub 15・dt・IC・BC・GPU・BLOCKSIZE 固定で各 200 物理 step (`run_0960_sglsq_s2_fct_gg_base36d8` / `run_0961_sglsq_s2_fct_gg_nokey`)。同じ checker (`check_passive_budget --mode fct`)・同じ閾値で収支閉合・独立総量照合・remainder・有界性を比較。**A**: 旧版も同じ remainder 項で FAIL → 「本 plan の commit 群が初めてこの FAIL を生んだ」仮説を棄却 (lsq の精度ゲートは未合格のまま)。**B**: 旧版 PASS・現行 FAIL → 既存問題として切り離す判断を棄却し commit 群の共通経路へ戻る。診断列の意味・実効設定が揃わなければ判定不能 (数値の近さで代用しない) | O (結論 F) |
| 5n (codex result-2、採用。**完了 2026-09-26**) | result-2 の指摘対応 | M1: `stage_manifest.py` は PyYAML が無ければ止まる (**済**、回帰試験 5 例)。M2: 一次記録を `case/09.Taylor-Green/_g0_lsq_seam/s2_evidence/` に回収 (**済**: `JUDGEMENTS.txt` = 元ファイルに対する判定ツールの全出力、`INDEX.md` = 元の残差 CSV の sha256。残差の gzip (59 MB) はリポジトリに入れずワークツリーと AWS に保存)。M3: 本表と §4.4 の 2 段ゲートどおり **plan 全体は `in_progress`**、Phase 1 の判定は §6.1 に記録。m4: 文書・索引・本表の状態を同期 | O |
| 5o (**決着 2026-09-27 ユーザ決定「a」: 合格**。lsq−gg の remainder 差 +2.8e-12 は gg 同士の揺れ 1.9e-12 と同程度なので「lsq が追加の不合格を作らない」とみなす。許容幅を後から決めた判断であることを記録) | FCT smoke の収支 remainder の比較規則 | 差し替え規則は「lsq の remainder が gg 以下」と書いたが許容幅を決めていなかった。丸め前 remAbs: gg `run_0954` 3.523596956e-06 / lsq `run_0955` 3.523599715e-06 (lsq が +2.8e-12、相対 8e-7)。同じ現行 gg の 2 本 (`run_0954` と `run_0961`) の差は 1.9e-12、nSub 30 の gg/lsq は 3.523593235e-06 / 3.523597981e-06。**文言どおりなら不成立、差は gg 同士の揺れと同程度**。許容幅を後から決めるので、ユーザ判断とする | F |
| 6 (**実装・確認済 2026-09-27** (commit `24365d72`、AWS バイナリ `fcfdc826…`): node 省略 → `lsq (default)` + 警告、cell 省略 → `gg (default)`、明示 gg → `gg (explicit)`。省略の既定 lsq と明示 lsq の 1 step 比較は tgv・case/48 とも E/P = A・N PASS (`DEFSW_{tgv,case48}.txt`、`defsw.py`)。S4 は下の差し替えにより「警告 + 運用ルール」で代替。**前提の差し替え 2026-09-27 ユーザ決定「B」**: #2g (起動記録の結び付け) の完了を前提にしない。代わりに (i) node で `scalarGradient` 省略時に既定変更の警告を出す、(ii) **既定切り替え日をまたぐ run は途中から再開せず最初から回し直す**を運用ルールにする (`procedures/solver-settings.md`・`recommended-settings.md` §9)。#2g の本格修正は `tooling-stage-manifest-launch-binding` に残して後回し) | 既定の切り替え (Phase 2) | **前提**: 5r が GO、2g 完了、S3 ≤ 5 %。S4、2h、docs、codex result 2 回目 | F |
| 7 (**完了 2026-09-26**: `boundary-node-rotational-periodic` §5.1 #0a) | 回転周期 plan への登録 | | O |

## 6. 検証 (2026-09-26 `diagnostician` 確定、同日 codex plan 対応で改訂。測る前に固定)

- **S0 作用素** (ハーネス `case/09.Taylor-Green/_g0_lsq_seam/g_suite.py` の 7 変種 + 軸対称×周期 1 本 + `axisymMethod: 1` 1 本):
  - **S0-a 純作用素 (BC 前)**: k・ω・ξ・Y (5 種、チャンク境界 N=5 を跨ぐ) の LSQ 勾配が double 参照と ≤ 1e-5·S。参照は並進周期で `lsq_merged_ref`、軸対称×周期は**非合併** (root = identity) の同関数。非線形場 (sin) で識別。
  - **S0-b 定数場**: ξ (channel は入口が無くピン無し) で全節点・壁節点込みで勾配 == 0 (厳密)。k/ω は非ピン節点で == 0、壁隣接節点はピン値 (k=0, ω_w) からの差分形の解析値と ≤ 4ε。
  - **S0-c NS との一致**: 同じ場を ro に入れた NS 勾配と、**周期 gather 前の局所配列がビット同一** (周期の無い変種 — box の全面 slip — では gather が no-op なので `res_1` の最終配列同士を直接比較。化学種は Y を ξ と同じ場にして dξ と dY のビット一致、5 種でチャンク境界 4+1 を通す)。gather 後は 2 member group がビット同一、3 member 以上 (多重周期の辺・角、atomicAdd の順序が非決定) は |差| ≤ 4ε·Σ|部分和|。
  - **S0-c の軸対称 2 変種 (`axi`・`axi_m1`) は適用除外** (2026-09-26 記録、codex (diagnose) の指摘で明文化): ハーネスは TP 混合気の状態を ρ = q で作り直せない (`g_suite.py:452,527`)。除外の根拠: (1) 軸対称では周期 gather も継ぎ目合併も走らない (`periodicGradientGather_d_wrapper` の `isAxisymmetric` 早期 return、`periodicSeamMergeActive` は非軸対称限定) ので、比べる対象はカーネルと係数だけ; (2) カーネルは非軸対称 8 変種の S0-c で NS とビット同一; (3) 係数は NS と同じ `cInt` 配列 (軸対称でも planar LSQ); (4) 軸対称固有の取り違え (A_planar 除算の残り等) は O(1) の誤差になるが、S0-a で double 参照との差/S ≤ 7.6e-8 (`S0_axi.txt`・`S0_axi_m1.txt`)。
  - **S0-d 検出力 (負の対照)**: jitter32 で GG 参照との差が 1e-5·S を超える (閾値が GG と LSQ を識別する証拠)。
  - **S0-e BC 後の同値性 (診断)**: channel と case/39 起点で 1 step、applyBconds 直後の k・ω・wall_y_eff の group 内差 == 0 (記録。非零なら F 項目)。
- **S1 非干渉**: `scalarGradient: gg` で実装前バイナリと全配列が一致 (旧経路の保存、roK/roOmega 初期ミラー追加も含めて)。判定規則 (**測る前の訂正 2026-09-26**、`diagnostician`): 旧同士がビット一致する配列は旧新もビット一致。旧同士でも一致しない配列 (周期 gather の atomicAdd、3 member 以上の group) は、旧 3 本・新 3 本の同一設定反復からノイズ対 (旧旧 3 対・新新 3 対) を取り、旧新 (9 対の最大) が (a) 最大差 ≤ ノイズ対の最大差の 2 倍、(b) 不一致数 ≤ ノイズ対の不一致数の最大の 2 倍。不一致の位置が 3 member 以上の group に限られるかは記録 (判定外)。~~R3 の「不一致数の桁が同じ」規則~~は廃止 (桁境界 8↔12・93↔108 で非決定的に反転する。前提 plan #6a と本 plan の Phase 1 最小確認で再現)。前提 plan の #6a の FAIL 記録は書き換えない。`lsq` で NS の勾配・リミタ配列が `gg` と同じ step でビット一致 (1 step、初回ダンプ)。
- **S1 (2) の判定手順の訂正 (#4b が A の場合に適用。2026-09-26 `diagnostician`、測定後だが決定的な部分には厳しくする方向)**: ~~NS の勾配・リミタは~~ **NS の 18 勾配配列は gather 前の局所配列でビット一致** (決定的段、必須。リミタは gather の**後**に評価される (`main.cpp:1479,1485`) ので gather 前の比較対象にしない — 2026-09-26 codex (diagnose) で訂正。リミタは評価後の配列を下の E/N 区分で比べる)。gather 後は 2 member 以下の group はビット一致、3 member 以上は部分和の順列和の集合に含まれる (代替: \|差\| ≤ 4ε Σ\|p\|)。「旧同士がビット一致なら旧新もビット一致」は gather を含む配列には適用しない (一致は scheduling の偶然で配列の性質でない)。tgv の初回 FAIL (dUxdy 36・dUxdz 32、3 member 以上のみ、最大 1.19e-7) は §6.2 に記録として残す。
- **ダンプ非干渉 (#4a) と S1 の 1 step 比較の区分 (事後規則。2026-09-26 codex (diagnose)、`notes/reviews/2026-09-26-gradient-scalar-lsq-4a-closeout-diagnose.md`)**: 規則が事前に無く、S1 規則の流用 (片側条項「旧同士ビット一致なら旧新もビット一致」) で tgv (res_0 dUxdz・res_1 dUxdy) と case48 (res_1 roe) が FAIL した (`PREGATHER_4b.txt` 末尾・`DUMPCHECK_case48.txt`。**この初回 FAIL は書き換えない**)。観測値で決まる「両側ともビット一致なら厳密」ではなく、**コード上の依存関係で区分を先に決める**:
  - **E (厳密一致)**: res_0 の状態量 (保存量・原始量。初期値 + BC で atomicAdd を通らない)。NS 18 勾配配列 (+ divU) の周期 group に属さない節点 (LSQ gather は atomic なし)。周期の無いケースのリミタ配列 (入力が E、`limiter_d.cu` の atomicAdd は診断 `g_limDiag` 既定 off のみ)。
  - **P (順列集合)**: NS 18 勾配配列の 2 member 以上の節点。ダンプ有りの pre-gather 部分和から作った float32 順列和の集合に含まれること (2 member は唯一の和)。
  - **N (両側ノイズ)**: 上以外 (面ループ atomicAdd の残差・GG スカラー勾配・res_1 の状態量、周期ケースのリミタ等)。S1 の訂正規則 (両側 3 本のノイズ対、最大差 ≤ 2 倍かつ不一致数 ≤ 2 倍) を適用し、不一致数と最大差を 3 組とも記録する。片側条項は N に適用しない。
  - やらないこと: FAIL の配列だけの例外登録、両側ノイズの 2 倍以内だけを根拠に全配列を「数値不変」と書くこと、判定を通すために gather を決定的実装へ変えること (§7 の後続項目は別途)。
- **S2 物理 A/B — 共通規則**: 起点の最終場から `restart_field.py` で gg/lsq の 2 本、同じバイナリ・同じ step 数。収束は (i) `check_convergence` が DIVERGED でなく RISING 列が無いこと (plateau は既知床として可、先例 case/48 `run_0025`)、(ii) 各列の末尾平均が起点の末尾平均の ≤ 1.5 倍 (ピークは除外し、再進入 step を記録)、(iii) 比較量が `check_quasisteady` で STEADY (閾値は表)。~~`--from-floor`~~ は使わない (codex M4: 全期間ピークを見るので作用素切替直後の跳ねで落ち、case/48・case/40 の起点は plateau で参照側が REFUSED)。ただし case/16 `run_0476` は起点が通常判定 PASS なので `--from-floor` を使える。未達なら同一設定で延長 (初回予算は表の step 数)。**物理ゲートは独立の必須条件**。差が上限を超えたら自動不合格にも自動合格にもせず、格子対診断 (両方向・両経路) を回して F 判断 (codex M5: 縮むかだけでは欠陥と離散化差を識別できない)。

  | ケース | 起点 (所在・sha256 先頭) | step | 必須ゲート (既存) | 差 lsq−gg の上限 (超えたら格子対 → F) | STEADY 閾値 |
  | --- | --- | --- | --- | --- | --- |
  | case/48 冷却平板 | main ワークツリー `/home/sano/work/forge/case/48.flat_plate_cooled_m4/run_0005_B_tw300_y3/res_48000.h5` `9a9022d1…`、mesh `10a13dd0…` | 24000 | Cf/VD-II ∈ [0.96, 1.00] @x 0.3/0.6/0.9、2St/Cf ∈ [1.13, 1.19]、エネルギー閉合 ∈ [1.00, 1.03]、NaN 0。壁ノード F1 = 0 (両経路、記録)、第一内層 F1 の L∞ 差 (記録) | Cf・q_w・δ*・θ 各 ≤ 1 % (3 station)。格子対: 流れ方向 `fp_y1_3um_nx1500` (1.5×、期待比 ≤ 0.75)・壁方向 y1 6 µm (2×、期待比 ≤ 0.6)。差が反復床 (1.8e-6 相対) の 10 倍以上のときだけ縮小率を評価 | `--drift 0.002 --osc 0.005 --tail 0.4` |
  | case/40 軸対称ノズル (axisymMethod 0) | main ワークツリー `/home/sano/work/forge/case/40.nozzle_design_tool/run_0045_node_yp1_outletfix_cont/res_12000.h5` `15e47af0…`、nozzle.h5 `849660af…` | 12000 | η_CF ∈ 0.978 ± 0.003 (README:172、抽出は `design/forge_design/metrics/extract.py` の `thrust_metrics` → `eta_cf`)、ṁ ∈ 1.29–1.30 kg/s、NaN 0、出口列に負値なし | η_CF ≤ 0.1 %、壁温 L∞ ≤ 15 K、ṁ ≤ 0.2 % | `--drift 0.0005 --osc 0.002` (η_CF, ṁ) |
  | case/39 周期丘 | sern ワークツリー `case/39.periodic_hills/run_0039_r1_gradfix_new_ext` の最終 res (AWS、`--from-floor` PASS 済) | 200k | 継ぎ目比 r_gradk・r_gradw ∈ [0.9, 1.1]、NaN 0、S0-e で同値 | Cf 3 点・x_r は記録のみ。dF1_inf: lsq で STEADY なら「**case/39 のこの構成で**の GG 経路の性質」と限定して前提 plan #7(1) を閉じる、DRIFTING なら継続 | `--drift 0.002 --osc 0.005 --tail 0.4` |
  | case/16 `run_0476` 双子 (Y 5 種 + ξ、SFR 2、Euler) | `run_0476/res_24000` | 24000 + 継続 24000 | 継続区間で `--from-floor run_0476` PASS、0 ≤ Y、\|ΣY−1\| ≤ 1e-6、0 ≤ ξ ≤ 1、floor 補正累積 ≤ 1e-6 | vs gg 固定点: Y0 max\|Δ\| ≤ 3.7e-5、Xi ≤ 9.4e-4 (S3 信号の 10 %。反復床は 2.0e-7 / 1.5e-5) | `--quantity machmax,pmax` STEADY |
  | case/16 `run_0482` 双子 (凝縮 S3、SST TP、モーメント 4) | `run_0470/lim1e` 入力 + IC `run_0213` (0482 と同一) | 48000 | `check_convergence` PASS (0482 と同基準)、ALL STEADY、モーメント ≥ 0、NaN 0 | vs 0482: onset (中心線 g = 1e-3 の線形補間交差、格子間隔 0.426 mm) ≤ 0.3 mm (S3 効果 0.72 mm の 4 割)、壁 p/p0 L∞ ≤ 0.3 % (`extract_wall_pp0.py` の CSV、S3 効果 0.81 % 比)。「実験に近づく」は判定に使わない | `compare_condfix --series-csv` ALL STEADY |
  | dual-time FCT smoke (`run_0486` プロトコル + `run_0476` 入力のトレーサ版) | `run_0471/res_12000`、`run_0476/res_24000` | 200 | NaN 0、\|ΣY−1\| ≤ 1e-7、0 ≤ ξ ≤ 1 | vs gg 双子 (同一物理時刻): ro ≤ 3e-6、roY ≤ 5e-6、roUy ≤ 1e-4 (反復差 0486 vs 0487 の 2 倍) | — |

- **S2 収束ゲートの基準の差し替え (2026-09-26 ユーザ決定「古い基準なんて気にしなくていい」、判定を見る前に固定)**: 起点 run (旧バイナリ・旧既定) の残差床は使わない。**基準は同じ起点から同じ現行バイナリで回した gg 双子**とし、(i) gg 双子に NaN/Inf・RISING が無い、(ii) lsq 双子の各残差列の末尾 20 % 平均が、**同じ step 区間**の gg 双子の末尾 20 % 平均の 1.5 倍以内 (`check_floor_ratio --start <gg 双子> <lsq 双子>`)。区間は各ケースの最終段 (case/48・case/40 は延長段 `run_0952`/`0953`、0476 は継続段 `run_0956`/`0957`、0482 は `run_0952`/`0953`、case/39 は `run_0950`/`0951`)。物理ゲートと lsq−gg の上限は従来どおり。**FCT smoke**: 保存収支は同じ checker で lsq が gg と同じ判定・同じ remainder 以下 (lsq が追加の不合格を作らない) を合格条件とし、nSub 感度 (#5h = B) のため lsq−gg 差の精度審査は行わない (記録のみ)。旧判定 (起点床・旧上限) の不成立は記録として残す。
- **S3 性能**: case/48・case/39 で lsq/gg の step 時間比 ≤ 1.05 (AWS native、同一 GPU・BLOCKSIZE、ウォームアップ 500 step 後 2000 step × 3 反復の中央値)。**超過は #6 を保留**。REG・spill・実測を記録。
- **S4 provenance (#6 の前提、別 plan #2g で実施し結果を写す)**: 同一 YAML を旧バイナリ (既定 gg) → 新バイナリ (既定 lsq) で 2 段起動し `stage_manifest.py --segments` が **2 区間**を返す。旧省略 / 新省略 / 明示 gg / 明示 lsq / cell (実効 gg) の 5 通りで実効値と由来 (launch / explicit / legacy / inferred) が正しい。設計 DB で FLAG_POLICY 旧行が選別から外れる。
- **前提 plan #7(4) (ジッタ格子の LSQ 次数)**: ~~本 plan で「欠陥ではなく作用素の性質」として閉じる~~ (codex M3 で訂正) → 本 plan では閉じない。スカラーは NS と同じ LSQ の次数の性質を継承する (S0-c のビット同一が根拠)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-unification-result.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-unification-result.md) | **NO-GO**, C0/M5/m2 (実装の中核に欠陥なし。M1 床比ツールの誤合格・M2 収束ゲート未達の根拠不足・M3 FCT 収支未確認・M4 S3 監視漏れ・M5 manifest に scalarGradient なし) | **全件採用** (codex (diagnose) 4 で採否、M2・M3 の規則を修正して固定) → §5.1 #5f–#5l。既定 gg・`in_progress` を維持、Phase 2 へ進まない |
| result | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-unification-result-2.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-unification-result-2.md) | **NO-GO**, C0/M3/m1 (新基準の S2 残差比 PASS は記録と整合・中核実装に欠陥なし。M1 manifest の PyYAML 欠如経路・M2 一次記録の回収・M3 Phase 1 完了と plan 完了の混同、m4 文書同期) | **全件採用** (いずれも数値判断を変えない手続き・ツール修正のため親判断で採用、上位諮問は省略) → §5.1 #5n。**Phase 1 判定 (2026-09-26)**: 現行 gg 基準の S2・S0/S1・S3 は合格、FCT smoke の収支比較規則のみ #5o で未決。plan 全体は `in_progress` |
| diagnose (既定切り替えの実装方針) | `2026-09-27` | [2026-09-27-scalar-lsq-default-switch-diagnose.md](../../notes/reviews/2026-09-27-scalar-lsq-default-switch-diagnose.md) | Major 1 (`FLAG_POLICY` の日付だけでは gg 評価が学習に混ざる) | **採用**: 実効 scalarGradient の学習ゲート (#2h)、判別試験 A |
| diagnose (諮問 5、Phase 1 の現状) | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-phase1-status-diagnose.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-phase1-status-diagnose.md) | Major 4 (合格を求める再審査・FCT の完了条件外し・「物理量を変えない」主張を却下、床移動の調査移管は採用だが S2 免除は却下) | **全件採用**。`in_progress`・既定 gg 維持、#5m (FCT smoke の旧版/現行 gg A/B) を先に。result 再審査の提出文は「測定した 2 ケースで S3 を通過したが、S2 は収束・FCT ゲート未達で Phase 1 は未合格」 |
| diagnose (諮問 4、result 1 の採否) | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-result1-disposition-diagnose.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-result1-disposition-diagnose.md) | C0、7 件採用 (M1 の修正に残る穴・M2 は双方向床比と restart 手順を揃える・M3 は作用素差と確定しない・M4 は PID 照合) | 全件採用 |
| diagnose (諮問 3、S2 解釈) | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-s2-result-diagnose.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-s2-result-diagnose.md) | C0、Major 4 (現行 gg 床への差し替え却下・床移動は要再検証・FCT は 0476 差の流用却下・case/16 の不足区間と比較対象)、Minor 1 (chi 文書不整合) | **全件採用**。S2 は条件付き・未合格、#5a–#5e を追加 |
| diagnose (諮問 2) | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-4d-nfail-diagnose.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-4d-nfail-diagnose.md) | C0、Major 2 (#4d の A をダンプ非干渉と記録しない・N FAIL 継続を「ダンプが分布を変えた」と確定しない) | **全件採用**。gg/lsq のノイズ対プールは却下 (親の案 (a) も採らない)、#4f の 8 本試験を測る前に固定、S2 は E/P と #2e・#2f が条件 |
| diagnose (諮問、レビュー 2 回とは別) | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-4a-closeout-diagnose.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-4a-closeout-diagnose.md) | C0、Major 2 (ダンプ判定規則の原案却下・S0/S1 は再判定が要る)、Minor 2 (hook)、Minor 1 (S1 (2) のリミタは gather 後) | **全件採用** (Fable 上限中のため codex で諮問)。#4a commit 可 (Minor 修正後)、§6 に E/P/N 区分、#4d・#4e、#5 の着手条件 |
| plan | `2026-09-26` | [2026-09-26-gradient-scalar-lsq-unification-plan.md](../../notes/reviews/2026-09-26-gradient-scalar-lsq-unification-plan.md) | **GO-with-changes**, C0/M7/m2 (差分形 LSQ の係数共有は支持) | **全件採用** (2026-09-26 `diagnostician`)。M1 → 壁 F1 = 0 と axisymMethod 1・σ 補間の訂正 (§3・§4.2)。M2 → roK/roOmega 初期ミラー、スカラー gather の述語を `periodicSeamMergeActive` に統一 (§4.3)、壁 ω の同値性は既存問題として S0-e の診断に格下げ。M3 → S0 を a–e に分割。M4 → `--from-floor` を捨て plateau 許容 + 床比 + STEADY、起点を main ワークツリーで確認し sha256 固定。M5 → 物理ゲートは独立の必須条件、格子対は診断で F 判断、STEADY 閾値を締める。M6 → case/16 の 2 双子 + FCT smoke を追加。M7 → 既定化の provenance を別 plan (#2g) に切り出し #6 の前提に。m8 → アクセサの契約 (§4.1)。m9 → S3 を #6 の保留条件に。分割はしない (2 段ゲート) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu`、`ransTransport_d.cu`、`speciesTransport_d.cu`、`periodicNode_d.cu`、`input/solverConfig.{hpp,cpp}`。
- 後続項目 (本 plan では入れない): `periodicGather1ToRoot_d` (`periodicNode_d.cu:57`、atomicAdd) を root スレッドが member を index 順に足す決定的 gather に置き換える。全 node 周期 run が ulp で動くので、別項目として accepted plan の R3 型で検証する。
- **既存の node SST・化学種再構成・受動種 run の結果が変わる** (勾配作用素の変更)。

## 8. 完了条件

- [ ] `methods/` 更新
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-27` — **Phase 2: node の既定を lsq に切り替え** (`24365d72`)。省略時の警告、設計 DB の実効値ゲート、文書 (`procedures/solver-settings.md`・`recommended-settings.md` §9・`methods/`)。既定 lsq ≡ 明示 lsq を tgv・case/48 で確認。

- `2026-09-27` — ユーザ決定「B」: 既定化の前提を #2g 完了から「省略時の警告 + 切り替え日をまたぐ run は最初から回し直す運用ルール」に差し替え。

- `2026-09-27` — #5o をユーザ決定で合格 (FCT smoke の remainder 差は gg 同士の揺れと同程度)。**Phase 1 は全項目合格**、次は Phase 2 の前提 #2g (別 plan 起票)。

- `2026-09-26` — codex result 2 回目 NO-GO (C0/M3/m1、新基準の S2 PASS は記録と整合) を全件採用: manifest の PyYAML 必須化、一次記録の回収 (`s2_evidence/`)、Phase 1 判定を §6.1 に記録し plan は `in_progress`、FCT の remainder 比較規則を #5o (ユーザ判断) に。

- `2026-09-26` — 差し替え後の基準で S2 PASS (lsq/gg 末尾残差 ≤ 1.42 倍、5 ケース)。

- `2026-09-26` — ユーザ決定: S2 の収束ゲートの基準を起点 run から**現行バイナリの gg 双子**へ差し替え (§6 S2 の差し替え規則、判定前に固定)。

- `2026-09-26` — #5m = A (FCT smoke の収支 FAIL は実装直前版 gg でも同値で再現)。

- `2026-09-26` — S3 を PID 照合版で取り直し (競合検出試験込み) case/48 0.997・case/39 1.018 で PASS。codex (diagnose) 5 で Phase 1 は `in_progress`・未合格のまま、床移動の真因は `notes/investigations/2026-09-26-node-steady-floor-shift.md` へ移管 (S2 免除ではない)、#5m を追加。

- `2026-09-26` — #5g = A (実装直前版 36d8ba03 と現行 gg の case/48 床比 0.99–1.09・物理差 ≤ 0.023 %)、#5h = B (FCT smoke は nSub 15→30 の感度が lsq−gg 差の約 2 倍、収支 remainder は 4 本とも FAIL で作用素・nSub に依らない) → FCT smoke の精度差は審査保留。#5f 取り直しで判定不変。S3 は M4 修正版で取り直し待ち (別セッションが GPU 使用中)。

- `2026-09-26` — codex result 1 回目 **NO-GO** (C0/M5/m2、実装の中核に欠陥なし)。codex (diagnose) 4 で 7 件とも採用、#5f–#5l を追加 (M2・M3 の判定規則は測る前に固定)。M1 は修正済。

- `2026-09-26` — #5a = B (chi 0 でも case/48 は起点床に戻らない → chi 単独説明を棄却)、#5b 0476 継続は上限内、#5c vs 0482 は上限内、#5d FCT 保存・有界性を測定、S3 PASS (1.010 / 1.017)。

- `2026-09-26` — S2 を AWS で実施 (4 case README の `run_095x_sglsq_*`)。物理ゲートと lsq−gg の上限は FCT smoke 以外 PASS (case/48 \|Δ\| ≤ 0.016 %、case/40 Δη 0.000 %・壁温 0.62 K、0476 Y0 1.2e-5・Xi 3.1e-4、0482 onset 差 0・壁 p/p0 0.003 %、case/39 全量 STEADY で dF1_inf も STEADY)。起点床の収束ゲートは設定差のある 4 case で gg・lsq とも不成立、case/39 (設定差なし) は両方 PASS。codex (diagnose) S2 で「条件付き・未合格」、#5a–#5e を追加。

- `2026-09-26` — #4f = A (8 本ずつ、E/P 厳密・N 全配列 2 倍以内)。#4a の hook を Minor 2 件修正して commit。S0/S1 を PASS で閉じた (#4)。#2f の `check_floor_ratio.py` を実装。S2 は #2e (起点の AWS 転送と sha256 照合) から。

- `2026-09-26` — #4d: E/P = A (tgv 12 本・case48 6 本)、N は tgv lsq res_1 `res_ro` の 1 配列が 2 倍規則 FAIL。codex (diagnose) 諮問 2 (`notes/reviews/2026-09-26-gradient-scalar-lsq-4d-nfail-diagnose.md`、全件採用) で「ダンプ非干渉は未確認」と記録し、一度限りの 8 本試験 (#4f) を測る前に固定。

- `2026-09-26` — 引き継ぎの未決 3 点を codex (diagnose) に諮問 (`notes/reviews/2026-09-26-gradient-scalar-lsq-4a-closeout-diagnose.md`、全件採用): (1) #4a hook は commit 可 (Minor 2 件修正後)、(2) ダンプ非干渉は事後規則 E/P/N 区分で再判定 (#4d)、(3) S1 (2) 訂正の適用は支持するがダンプ無し出力の集合検査を先に (#4d)、S0-c 軸対称は適用除外を明記 (#4e)、S2 は #4d = A と #2e・#2f の後。

- `2026-09-26` — #4b = A、#4c PASS、#4a は実装済み・未 commit。hook の diff レビューと「ダンプ有効で res 不変」の判定規則、§6 S1 (2) 訂正の適用可否は上位モデルの上限到達で未決 (引き継ぎ `notes/sessions/2026-09-26-gradient-scalar-lsq-handoff.md`)。

- `2026-09-26` — S0 (測れた範囲 PASS) と S1 (tgv の NS 配列で FAIL) の結果を受け、`diagnostician` 判断で出力 hook (#4a)、tgv の決定的 A/B (#4b)、Y の未測定項目 (#4c) を追加。

- `2026-09-26` — Phase 1 実装 (#3・#2b)。diff を `diagnostician` がレビューし欠陥なし。S1 の判定規則を測る前に訂正 (「不一致数の桁」廃止、両側 3 本のノイズ対)。

- `2026-09-26` — codex plan (GO-with-changes C0/M7/m2) を全件採用して §3–§6 を改訂 (`diagnostician`)。前回の「壁 F1 は 1 に張り付く」は誤りで訂正 (k=0 ピンで F1=0)。2 段ゲート (Phase 1 opt-in → Phase 2 既定化) に。

- `2026-09-26` — §4・§6 を `diagnostician` 判断で確定 (opt-in から開始、S0–S2 通過後に本 plan 内で既定化、既存ゲート + 格子対の合否規則)。
- `2026-09-26` — 起票 (ユーザ決定「スカラー勾配も LSQ に揃える」、前提 plan の accepted を受けて)。
