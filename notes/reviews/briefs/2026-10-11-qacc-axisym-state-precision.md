# 諮問ブリーフ: float のずれは「状態の保存・commit の丸め」か「残差の評価」か — qAccumulatorFP64 を軸対称に対応させて軌跡の腕を足してよいか (2026-10-11)

- 依頼者: 主セッション (AGENTS.md のエスカレーション条件 1: §6 を足す、4: 承認済みの手順に無い run の前、6: `cuda_forge/` の数値の振る舞いを変える編集の前)
- plan: `plans/active/architecture-float-state-double-geometry.md` (§5.1 #15、§6.26〜§6.33)。関連: `plans/active/time_integration-fp64-accumulator-rollout.md` (§4.1 の S1b-① 軸対称、§5.1 #2。**この plan の担当は別のセッション (すきま加熱)。いまは動いていない**)
- ユーザの承認: 「この A/B を諮ってから implementer に任せる」流れに「よい」(2026-10-11)。実装・run の前に本諮問を通す。

## 0. 問い
1. 下の §2 の実装 (軸ピンを FP64 正本にも当てる) は正しいか。抜けている writer・順序の穴はないか。
2. 下の §3 の腕と判定は、「float のずれは平均流の状態の commit の丸め (1 ULP 未満の増分の吸収) が通り道か、残差の評価 (流束・勾配・リミタ・EOS) の側か」を分けるのに足りるか。判定の式の倍率 (§3.4) をどうするか。
3. §6.26 の A (`run_0487`、FP64)・B/B′ (`run_0488`・`0489`、float) を比べる相手として使い回してよいか (バイナリが変わる。§3.3 の対照 F で確かめる案)。

## 1. 観測事実 (float plan の §6.27〜§6.33、run はすべて case/45、インスタンス B)
- §6.27 (確定): 共通の Q32 の起点 (`_tr/q32_init.h5`、sha ba859e02: 保存量 9 つを float32 に丸めて double に戻した場) から 30,000 step。FP64 の A は θ_r が ≤ 0.12 % しか動かず、float の B・B′ は +1.0〜1.9 % 動く。
- §6.29: SST の輸送の更新を止めても float は凍結の FP64 から離れる (判別不能)。
- §6.31: 最初の 1 回の更新の段ごとの差は再実行の揺れの 1.5〜3 倍で、段を分けられない (判別不能)。
- §6.33: block の系の組み立てを double にした C・C′ (`implicitSolvePrecision: 1`) は B とほぼ同じに離れた (mC は mB の 94〜97 %)。判定は判別不能 (結果 B の 6 条件のうち 5 が成立、W2 θ_r(40) だけ |mC − ΔA| 1.860 が 10·N_C 1.87 をわずかに下回った)。再実行の差 N_B は 0.07〜0.14 ポイント、N_C は 0.02〜0.19 ポイント。D = (+1.892, +1.647, +1.431) %。
- 設定 (`run_0483_hp7_k1/solverConfig.yaml` と同じ): node、軸対称 (`isAxisymmetric: 1`, `axisSegmentRWeight: 1`, `axisCentroidShift: 1`, `nodeWallDirichlet: 1`)、SLAU、TP 多成分 (MIXDRY + H2O、`thermalMethod: 2`)、SST (`sstEnergyIncludesK` は既定 0)、定常陰解法 (`timeIntegration: 11`, `blockDPLUR: 1`, `lineImplicit: 1`, `nStepInner: 5`, `implicitRelax: 0.7`, CFL 4)、`limiter: 2`。
- 既存の切り替え `time.deltaT.qAccumulatorFP64: 1` (v1a、accepted): 平均流の 5 保存量 (ρ, ρu, ρv, ρw, ρE) の正本を内点 CV の FP64 配列 Qacc に置き、commit を `Qacc += dq; Q = (float)Qacc` で行う。Q (float) は演算に使うミラーのまま。reconcile (`update_d.cu:51-59`、step の末尾の `updateVariablesOuter` の頭) は「Q ≠ (float)Qacc のセルだけ Qacc = (double)Q」。残余ゼロなら OFF とビット同一 (単体試験 `tests/unit/test_qacc_commit.cu` (d))。乱流 k/ω・化学種 ρY は float のまま。
- **v1a は軸対称を拒否している** (`main.cpp:3404-3407`)。理由: 軸ピン `enforceAxisSymmetry_d` (`axisymmetricSource_d.cu:303-325`) は step の頭 (`assembleResidualPre`、`main.cpp:1895`) で軸の節点の Q と commit の基準 roN・roUyN・roeN を射影する (roUy = 0、roe から 0.5 roUy²/ro を引く)。ON では commit の基準が Qacc なので、Q 側の射影は commit で `(float)(Qacc + d)` に上書きされて消える。reconcile は commit の後に走るので拾えない。rollout plan の反例: ρ=1, ρu_r=0.125, ρe=2, dq=0 → 正しくは Q・Qacc とも ρu_r=0, ρe=1.9921875。
- 1 step の中の順序 (定常、`advanceImplicitSteady`、`main.cpp:2518-2552`): `assembleResidualPre` (updateVariablesInner → 壁 no-slip `enforceWallNoSlip_d` (Q だけ、Q_N は触らない) → 軸ピン (Q と Q_N) → 等温壁ピン → EOS `dependentVariables` (TP は `dependentVariables_d.cu:209` で roe を `(float)(ro·(e_mix(T)+ek))` に書き戻す) → 物性 → 境界) → 残差 → block DPLUR + ライン Thomas → commit `applyBlockImplicitCorrection_d` (軸の分岐は nullptr で無効、`update_d.cu:335`) → SST・化学種の更新 → 等温壁ピン (もう一度) → `updateVariablesOuter` (reconcile → Q_N = Q)。
- 上の順序から: commit の前に Q だけを書く writer (壁 no-slip、EOS の roe の書き戻し) は OFF でも ON でも commit で消える (OFF は基準 Q_N が書かれていないので `Q_N + d`、ON は `Qacc + d`)。commit の前に **Q_N も**書くのは軸ピンだけ。commit の後・reconcile の前に Q を書くのは等温壁ピン (reconcile が採用する)。化学種の roe の補正 `species_energy_correction_kernel` は呼び出しが無い (死んだコード)。`sstEnergyIncludesK` は 0。周期なし。

## 2. 実装の案 (呼び出し側の案、未検証)
- `enforceAxisSymmetry_d` に `double* qacc_ro, qacc_roUy, qacc_roe` を足す (nullptr = OFF、OFF の演算列は今のまま)。ON の軸の節点:
  1. KE を今と同じ float の演算で作る: `ke = 0.5f * roUy[ic] * roUy[ic] / ro[ic]` (この時点で Q = (float)Qacc が成り立つ。commit の後の writer は等温壁ピンだけで、reconcile が揃える)。
  2. `qacc_roe[ic] -= (double)ke; qacc_roUy[ic] = 0.0;`、そのあと `roe[ic] = (float)qacc_roe[ic]; roUy[ic] = 0; Uy[ic] = 0`。
  3. roN・roUyN・roeN の射影は今のまま (ON では commit に使われないが、診断と今の作り (Q_N = Q) を保つ)。
  - 残余ゼロのとき (Qacc = (double)Q)、`(float)((double)roe − (double)ke)` は float の `roe − ke` と一致する (最近接丸めで f32(a+b) = f32(f64(a)+f64(b)))。よって OFF とビット同一。残余があるときは 1 ULP 未満の残余を保つ。
  - KE を Qacc から double で作る案もあるが、ビット同一を失う。軸の KE は d2 の 2 乗の大きさで小さいので float で十分と見た。
- `main.cpp` の軸対称の拒否を外す (軸ピンが Qacc を扱うようになったので)。ほかの拒否 (周期・`sstEnergyIncludesK`・dual-time・陽解法) は残す。
- 単体試験 (`tests/unit/test_qacc_axis.cu` を新設、または `test_qacc_commit.cu` に追加): 1 セルごとの処理を `__device__ inline` に切り出して駆動する。
  (a) rollout plan の反例: ρ=1, ρu_r=0.125, ρe=2, dq=0 → ピン → commit で Q・Qacc とも ρu_r=0, ρe=1.9921875。**旧い作り (Qacc に当てない) では ρu_r=0.125, ρe=2 に戻る**ことも同じ試験で示す。
  (b) 残余の保持: Qacc_roe = 2 + 1e-9 (1 ULP 未満) → ピンの後も残余 1e-9 が Qacc に残る。
  (c) ビット同一: 残余ゼロの乱数のセル 20 万で、ON のピンの後の Q が OFF のピンとビット一致。
- 採用の計器: 既存の `[qAccumulatorFP64] step N: reconcile 採用 X` (monitorInterval ごとの合計) に加え、**どの節点で採用したか**を見る診断を足す案。`FORGE_DIAG_QACC_ADOPT=<path>` と `<step>`: その step の `updateVariablesOuter` の直前で Q と Qacc を host に写し、量ごとに Q ≠ (float)Qacc の mask と節点の座標を h5 に書く (host だけ、カーネルは変えない)。期待: 採用は等温壁の節点の roe だけ (軸・内点は 0)。

## 3. 腕と判定の案
### 3.1 腕 (インスタンス B、3 本同時、各 30,000 step・500 ごと)
- E = float + `qAccumulatorFP64: 1`、E′ = E の再実行。§6.26 と同じ共通の入力 `_tr/q32_init.h5`・格子・キー 1・環境変数 (`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500`)。B との設定の差は `qAccumulatorFP64: 1` の 1 語だけ。
- F = 新しいバイナリの float、`qAccumulatorFP64: 0` (バイナリの変更の対照)。
- バイナリ: 新しい commit を B でビルドした float (`~/forge-fgeom9-f32`)。

### 3.2 前提 (総合の判定の前に確かめる。欠ければ判定不能)
- 単体試験 (a)(b)(c) が PASS。
- E・E′・F の全期間の残差・全出力が有限で 61 出力、`check_convergence` が実行できる。初期化後の保存量 (res_0) が B と E・E′・F でビット一致。
- E の採用の mask (step 1 と step 30,000 の 2 回): 軸の節点・内点の採用が 0 で、採用は等温壁の節点の roe だけ。monitorInterval ごとの採用数が (等温壁の節点数 × monitorInterval) 以下。
- 新しいバイナリと §6.26 の float のバイナリで、`enforceAxisSymmetry_d` 以外のカーネルの SASS が一致 (`cuobjdump -sass` の比較)。
- F が B・B′ の範囲に入る: 両窓・3 断面で min(ΔB, ΔB′) − N_B ≤ ΔF ≤ max(ΔB, ΔB′) + N_B。入らなければ、新しいバイナリの float の比べる相手を F に替えるか判定不能にする (事前にどちらかを決める)。

### 3.3 判定 (θ_r の 3 断面、W1 20,000〜25,000・W2 25,000〜30,000、Δ は §6.26 と同じ)
- mB = mean(ΔB, ΔB′)、mE = mean(ΔE, ΔE′)、N_B = |ΔB − ΔB′|、N_E = |ΔE − ΔE′|。
- 結果 S (状態の精度で抑制): max(|ΔE − ΔA|, |ΔE′ − ΔA|) + k·N_E ≤ 0.1|D|、かつ sign(D)(mB − mE) ≥ 0.5|D|、かつ |mB − mE| > k·max(N_B, N_E)。
- 結果 R (状態の精度では抑えられない): E・E′ がそれぞれ A から D の向きへ 0.5|D| 以上、かつ |mE − ΔA| > k·N_E。
- 両窓・3 断面で揃わなければ判別不能。Q_w は副判定。
- k は §6.32 では 10。§6.33 では N_C の 10 倍 (1.87) が D の大きさ (1.9) と同じ桁になり、ほぼ同じ挙動でも結果 B が 1 条件で外れた。**k を 10 のままにするか、再実行の差を 3 本目 (E″) で見積もって k を下げるか**を決めてほしい。結果を見る前に決める。

### 3.4 解釈の範囲
- 結果 S: 「平均流の 5 保存量の commit の丸め (1 ULP 未満の増分の吸収) が、登録の窓の偏りの通り道に要る」。k/ω・ρY は float のままなので、それらの丸めが無関係とは言わない。
- 結果 R: 「平均流の状態を FP64 で積んでも偏りが育つ」→ 残る候補は残差の評価 (流束・勾配・リミタ・粘性・EOS の float の演算) と、float のままの k/ω・ρY の状態。ON でも残差は float のミラー Q から評価されるので、「ミラーの丸めが残差を通して効く」経路はこの試験では消えない。
- どちらでも演算の特定ではない。

## 4. 読んでよいファイル
- float plan §5.1 #15・#17、§6.26〜§6.33
- rollout plan §4.1 (表の「node 軸対称」の行)・§4.1a・§5.1 #2
- `solver_density_cuda/cuda_forge/axisymmetricSource_d.cu` (300〜370)、`update_d.cu` (1〜100、160〜360)、`qAccumulator_d.cuh`、`qAccumulator_d.cu`、`nodeWallDirichlet_d.cu` (1〜140)、`dependentVariables_d.cu` (100〜300)、`speciesTransport_d.cu` (355〜390、845〜862)
- `solver_density_cuda/main.cpp` (1880〜1915、2350〜2420、2518〜2560、3380〜3450)
- `solver_density_cuda/tests/unit/test_qacc_commit.cu`
- `case/45.isobutane_m6_d155/tr.sh`・`ds.sh`・`ds_an.py`
