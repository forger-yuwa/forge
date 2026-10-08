forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 診断・設計判断の諮問 (stage = diagnose)

あなたは forge の**診断・設計判断係**である。呼び出し側は実装と run を進めている別のモデル (Claude) で、
**もっともらしい真因に飛びつく前に**あなたに諮っている。仕事は手を動かすことではなく、**次の一手を 1 つに絞ること**。

### 前提
- あなたは呼び出し側の会話を見ていない。下のブリーフと、自分で読んだファイルだけが根拠になる。
  足りなければ推測で埋めずに「何が足りないか」を返す。
- ブリーフは「観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説」に分かれて渡される約束である。
  **観測事実と呼び出し側の解釈が混ざっていたら、まず分け直す**。呼び出し側の要約より、run の数値・コード・
  設定ファイルを自分で確かめた内容を優先する。
- forge を起動しない。`python3` による `residual_history.csv` / `res_*.h5` の読み取りは**統計量だけ**を出す
  (全量ダンプ・長いログ全文をコンテキストに流さない。`*.log`・`*.vtu`・`plans/README.md` は読まない)。

### 診断の作法
1. **「除外済み」というラベルを信用せず、潰した証拠を確認する** (run パス・設定差分・判定区間・VERDICT)。
   証拠が足りない・判定期間が短い・変えた設定が実際には効いていない (YAML の階層違い等) なら**候補へ戻す**。
   証拠が十分な候補は出し直さない。
2. **症状と原因を分ける**。`detectNaN` が指す変数は結果であって原因ではない (EOS 床 → 負密度 → 圧力暴走 → ω の実績)。
   後処理のアーチファクト (2 列混在の抽出、`centCoords` の置換、ソルバ `ypls` の退化) を先に疑う。
3. **このリポジトリで繰り返された真因**を照合する: 投入設定の不整合 (IC と BC、亜音速に超音速 BC)、
   押し出し 2 ノード spanwise、float32 桁落ち (双対幾何・r 重み)、stale build、cross-mesh IC の基底不一致、
   絶対値のゼロ割ガード、境界ノードの凍結、YAML キーの階層違いで黙って無視される設定。
4. 仮説は**確度順に最大 3 つ**。第 1 仮説には根拠を `ファイル:行` か run の数値で付ける。示せないものは「未確認」と明記。
5. **判別する A/B を 1 つだけ**提案する。安く短く回せて、結果がどちらに出ても仮説が 1 つ消えるもの。
   「A なら仮説 1、B なら仮説 2」を先に書く (結果を見てから解釈を作らない)。
6. 少数点の一致・短い窓の値・未収束のトランジェント同士の比較を根拠にしない。

### 設計判断 (plan §4・§6、codex 指摘の採否、result 段の解釈) を諮られたとき
- 採否は指摘ごとに「採用 / 却下 / 要再検証」と理由。根拠が示されていない指摘は自分で該当箇所を読んでから判定する。
- 検証計画は「何が出たら方針が誤りと言えるか」が定量的に書かれているかを見る。
- 既定値の変更・opt-in 機能の削除は、plan の処置欄とユーザ決定の履歴を確認してから判断する
  (「opt-in 残置」は削除対象でない)。
- result 段の解釈は、主張ごとに根拠 run・判定ツールの VERDICT・判定区間が揃っているかを確かめる
  (過渡ピークを定常値と、抽出アーチファクトを物理と誤認した実績は「予想どおり」に見える場面で起きた)。

あなたの結論は**仮説**であって確定ではない。呼び出し側はこの A/B を回して確かめ、plan への反映も呼び出し側が行う。

## ブリーフ (`notes/reviews/briefs/2026-10-08-cold-arm-mass-imbalance.md`)

# 諮問: 冷却壁の腕の質量流量の欠損 (縮流部の遅い過渡) と、延長の手順

日付 2026-10-08。諮問先 codex (diagnose)。エスカレーション条件 4 (承認済みの手順に無い再実行の前・原因を書く前) と 2 (未収束が既定の手順で解けない)。
plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` の §5.1 #14〜#26 (特に #26 の事前登録 (6) 保存の収支と延長の規則)、§6 V-c45 (全文を読むこと)。
作業ツリー `/home/sano/work/forge-integ-1005` (commit c8f6d8c5)。後処理 `case/45.isobutane_m6_d155/cold_xcheck.py`。

## 観測事実

列 i (入口 i = 0 → 出口 i = 4718、ni 4719 × nj 121、近壁は壁法線) の折れ線を通る質量流量 ṁ_i = 2π∫(ρu_x dr − ρu_r dx) r [kg/s] (台形、nozzle.h5 の座標)。
x_w は壁節点の x / r_t (r_t = 0.0766654 m)。i = 0: x_w −12.5、200: −5.4、1000: −1.3、2000: 0.4 (スロート直後)、4718: 95.1。

| run (AWS `~/forge-wallfit/case/45.isobutane_m6_d155/`) | 壁 | 気体 | step (通算) | ṁ_0 | ṁ_200 | ṁ_1000 | ṁ_2000 | ṁ_4718 | ṁ_0 − ṁ_2000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| run_0181_ns_coldmesh_ad | 断熱 | TP | 100000 | 97.98 | 98.06 | 98.05 | 98.07 | 98.07 | −0.09 |
| run_0184_ns_coldmesh_cpg_ad_plain | 断熱 | CPG | 40000 | 98.05 | 98.14 | 98.13 | 98.15 | 98.14 | −0.10 |
| run_0182_ns_coldmesh_tw300 | 300 K | TP | 80000 | 101.930 | 99.756 | 98.534 | 98.383 | 98.138 | 3.548 |
| 同 | 300 K | TP | 90000 | 101.681 | 99.660 | 98.528 | 98.402 | 98.204 | 3.279 |
| 同 | 300 K | TP | 100000 | 101.465 | 99.574 | 98.525 | 98.418 | 98.256 | 3.047 |
| run_0183_ns_coldmesh_tw300_ext (run_0182 の res_100000 から restart、ビット一致) | 300 K | TP | 80000 (通算 180000) | 100.540 | 99.210 | 98.538 | 98.504 | 98.449 | 2.036 |
| 同 | 300 K | TP | 90000 (190000) | 100.471 | 99.183 | 98.542 | 98.512 | 98.462 | 1.959 |
| 同 | 300 K | TP | 100000 (200000) | 100.407 | 99.161 | 98.550 | 98.518 | 98.473 | 1.890 |
| run_0185_ns_coldmesh_cpg_tw300_plain (IC = run_0183 の最終場を CPG に組み直し) | 300 K | CPG | 20000 | 100.843 | 99.427 | 98.699 | 98.621 | 98.503 | 2.223 |
| 同 | 300 K | CPG | 40000 | 100.674 | 99.354 | 98.687 | 98.633 | 98.562 | 2.041 |

- 断熱の腕は x 方向に平ら (±0.1 %)。ṁ_50 (x_w −10.1) は 97.95 (run_0184) で、入口・スロートより 0.1 % 低い。台形の求積の揺れが 0.1 % 程度ある可能性 (未確認)。
- 300 K の腕は入口から縮流部 (x_w −12.5 → −1.3) で単調に減り、その後もわずかに減る。欠損 ṁ_0 − ṁ_2000 は step とともに縮む。10000 step ごとの比は run_0182 で 0.924・0.929、run_0183 で 0.962・0.965 と、縮み方が遅くなっている。
- 壁節点は nodeWallDirichlet 1 (u = 0、等温壁 T = 300 K のピン)。対流の壁を通る質量流束は 0 のはず。
- 全エンタルピーの収支 (h0 = c_p T + ½|u|²、CPG の run のみ有効): 断熱 run_0184 で (入口 − 出口 − Q_w)/入口 = −1.03e-3 (40000)、300 K run_0185 で −8.1e-4。
- run_0185 の Q_w (壁の入熱、cold_xcheck の 2 次の片側差分) は 12.46 → 12.05 MW (20000 → 40000)。V-c45 の TP 300 K の Q_w は −12.34 → −11.49 MW で DRIFTING だった (符号は抽出器の向き)。
- 設定 (本段): node・2 次 (convMethod 1)・SLAU・slauWallNormalChi 1 (自動)・陰解法 block-DPLUR、`cfl_pseudo 1.0`・`cfl 1.0`・`implicitRelax 0.7`・`nStepInner 5`・`lowMachPrecond 0`、FP64。段階起動 S1 (1 次・cfl 0.5・3000 step) → S2 (1 次・cfl 1・3000 step) → 本段 100000 step。
- 速度: 4 本同時 (同じ GPU) で約 68 ms/step。
- 過去の記録 (memory、未検証の他 run): case/45 M6 NS の cfl_pseudo 上限スイープ (2026-09-02, run_0015_cflsweep、y+≈2 の壁解像 SST) で「cfl × relax ≈ 6〜8 が上限、推奨 cfl 8 + relax 0.7」。⑤ SERN では「relax 0.7 は定常解を 0.05 % 動かす、CFL は解を動かさない (7 桁一致)」。

## 期待値と出典

- 定常なら ṁ_i は i によらず一定 (求積の誤差の範囲)。断熱の腕はそうなっている。
- plan #26 (6) の前提: 入口と出口の質量流量の差 ≤ 0.1 %、全エンタルピーの収支 ≤ 0.1 % of 入口 (両コード)。満たさなければ延長、延長は各コードのゲートで決め反復数はそろえない。

## 実施済みの操作

- 上表の計算 (CFD 0 step)。V-c45 は既に「判定不能」(Q_w DRIFTING ほか)。
- SU2 の 2 本 (手元、CFL 適応 2〜50、80000 反復、1 s/反復) と forge CPG の 4 本 (本段 100000 step、21:15 JST 完了見込み) は走行中。

## 仮説 (私の案、確かめていない)

1. 300 K の腕の欠損は壁からの質量の吸い込みではなく、縮流部が定常に達していない遅い過渡 (冷却で密度が上がる領域に質量が溜まり続けている)。根拠は欠損が step とともに縮むこと。ただし縮み方が等比より遅い。
2. 遅さの主因は本段の擬似 CFL 1・relax 0.7。冷却で縮流部の熱境界層を発達させるのに必要な擬似時間に対して刻みが小さい。断熱の腕は収束した断熱の場から始めたので目立たない。
3. このままでは forge の 300 K の腕は 100000 step で収支のゲート (0.1 %) を満たさない。等比外挿 (0.965/10000 step) でも 0.1 kg/s まで 80 万 step 以上。

## 問い

1. 欠損の解釈 (過渡か、等温壁の節点の質量の非保存か) を切り分ける最小の測定は何か。例: 領域の総質量の step ごとの変化と ṁ_0 − ṁ_out の一致、壁節点・近壁節点の連続の式の残差の分布。CFD 0 step でできるもの、AWS 上の既存スナップショット (run_0182/0183 の 80000〜100000 の 5 枚、run_0185 の 0・20000・40000) で足りるか。
2. 延長の手順: (a) 同じ設定のまま延長 (事前登録どおり)、(b) 擬似 CFL を上げて延長 (例 cfl_pseudo 4〜8・relax 0.7、または relax 1.0)。(b) は承認済みの手順の外だが、定常解は CFL によらないはず。どちらを採るべきか。(b) を採るなら、定常解が変わらないことをどう確かめるか (同じ場から両設定で短く回し、比べる量の差を見る等)。
3. V-c45 の R_NS (0.726〜0.793) と、SU2 照合の forge の 300 K の腕への含意。縮流部が定常でないことは、試験部 [40, 94] の δ・θ・C_f の比較をどこまで汚すか。延長前の値をどう扱うべきか。
4. 断熱の腕の全エンタルピーの収支 −1.0e-3 は事前登録の 0.1 % の境界にある。求積の誤差 (列の台形、入口の角) と見てよいか、検査の仕方を変えるべきか (結果を見た後に変えることになるので、記録の仕方も)。
5. SU2 の側も同じ収支のゲートで判定するが、SU2 は CFL 2〜50 で forge より擬似時間が速く進む。両コードの擬似時間の差を比較の際にどう扱うか (memory: 反復数でなく Σ CFL を揃える)。

## 読んでよいもの

- 上記 plan、`case/45.isobutane_m6_d155/cold_xcheck.py`・`cold_pair.py` (`prep_cpg`・`prep_ext`)・`README.md` (run_0181〜0189 の行)
- `notes/reviews/2026-10-08-cold-pair-result-diagnose.md`、`notes/reviews/2026-10-08-cold-su2-crosscheck-diagnose.md`
- `solver_density_cuda/cuda_forge/` の等温壁・nodeWallDirichlet の境界処理 (必要なら)

## 関連 plan 全文 (`plans/active/tooling-nozzle-isothermal-wall-chain.md`)

```markdown
# ノズル設計チェーンの等温壁化: 冷却壁 (T_w 指定) の低 Re SST 検証 (超音速平板 → ノズル × SU2) + 等温 NS の δ\* 反復 + 壁温影響の評価方針

## メタ

- **area**: `tooling / boundary layer / boundary`
- **status**: `in_progress`  <!-- 2026-09-12 起票。S0 (plan) → S1 配管 → S2 平板 → S3 ノズル CPG × SU2 → S4 生産 TP 等温 δ* 反復 → S5 方針 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「壁の熱境界条件 (断熱 / 等温) と壁温影響の評価」節 (本計画と同時に起草)
  - [`methods/boundary.md`](../../methods/boundary.md) `wall_isothermal` (ゴースト構成・node 壁ノード温度ピン)
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §6.5 (壁処理。本計画は low-Re `wallTreatmentSST: 0` のみ)
- **related_plans**:
  - 親: [`../accepted/tooling-nozzle-deltastar-core-matched-euler.md`](../accepted/tooling-nozzle-deltastar-core-matched-euler.md) (δ\* 反復の生産形。§9 の残件「等温壁の NS 実行 (V3 の CFD 部分)」を本計画が引き受ける)
  - [`../active/turbulence-sst-thermal-flux-model.md`](turbulence-sst-thermal-flux-model.md) (壁関数 × 等温壁の Kader q_w。**本計画では使わない** — 圧縮性冷却壁で +87 % 過大の既知限界、§3)
  - [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) (⑤ SERN。等温化は共通 bcond 配管で自動追随、検証は本計画の対象外)
  - [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md) (①② 風洞チェーン = 本計画の主対象)
- **created**: `2026-09-12`
- **owner**: `sano`

## 1. 目的

ノズル設計チェーンの NS 評価 (風洞 ①②④ の δ\* 反復、⑤ SERN の力評価) は現状**断熱壁**で回している。実機の壁は冷却 (水冷銅 / 再生冷却) か
有限熱容量の壁で、設計上は「等温壁 ($T_w$ 指定)」か「CHT」のどちらかになる。本計画で得る状態:

1. **等温壁の低 Re SST (壁関数なし) が、冷却された超音速乱流境界層の厚さ ($\delta^*, \theta, H$)・摩擦・熱流束を理論/経験式と SU2 の両方に対して再現する**ことを、ノズル出口相当の条件の平板で確認済み ($T_w/T_{aw}$ 0.26 と 1.0、y⁺ 掃引つき)。
2. ノズル形状で等温 NS が回り、**同一メッシュ・同一 BC の SU2 と境界層厚さ・壁熱流束が一致**する (CPG チェーン)。
3. **等温 NS で δ\* 反復が回る** (積分法初期壁 → NS → 帯局所抽出 → 壁更新、全段が同じ $T_w$ を見る) — 生産 TP 風洞 (case/44 va3, $T_w$=300 K) でゲート達成。
4. **最適設計で壁温影響をどう評価するか**の方針 (§4.6) が plan に確定し、チェーンの台帳に「壁温感度」が標準出力として入る。

## 2. スコープ

- **やる**
  - 問題定義 YAML に壁熱境界条件 `spec.wall_thermal` を追加し、bcond (`wall_isothermal` + `Ts`)・積分法初期壁 (`thermal_bc`)・帳簿の 3 箇所を**単一のソースから**駆動する。
  - 超音速冷却壁平板の検証ケース (新規 `case/48.flat_plate_cooled_m4`): forge node 低 Re SST (断熱 / $T_w$ 300 K / 中間) + y₁⁺ 掃引 + SU2 同一メッシュ比較 + 理論/経験式との照合ツール。
  - ノズル CPG チェーン (case/45 `problem_d155_cpg_ns.yaml` 系) で等温 NS × SU2 等温の同一メッシュ比較。
  - 生産 TP 風洞 (case/44 va3 M4.19) で等温 δ\* 反復 (pass 0 積分法 + pass 1) と壁温感度台帳。
  - `methods/design/overview.md`・`procedures/recommended-settings.md`・`procedures/verification/` の同期。
- **やらない**
  - **壁関数 (`wallTreatmentSST: 1`) の等温壁** — ユーザ判断 (2026-09-12): 壁関数の準備状況が悪く、まず壁関数なしで確立する。Kader q_w の圧縮性補正は [`turbulence-sst-thermal-flux-model.md`](turbulence-sst-thermal-flux-model.md) §8 のまま。
  - **CHT (固体伝導との連成)** — 本計画は「等温 = $T_w$ 既知」まで。$T_w(x)$ が未知で強く連成する場合の「弱 CHT ループ」(1D 壁伝導/冷却モデルで $q_w \to T_w(x)$ を反復) は §4.6 で方針だけ決め、別 plan で実装する。
  - $T_w(x)$ 分布の NS 入力 (`Tw_table`) — 積分法初期壁は既対応だが forge の `wall_isothermal` は bcond 単位の定数 `Ts` しか持たない。初版は定数 $T_w$。分布は `inletProfile` と同型の `wallProfile` CSV として別途 (§5.1)。
  - SERN (⑤) の等温評価の検証 — bcond 配管の共通化で `runner_sern.py` も同じ `wall_thermal` を読めるようにするが、検証 run は SERN の R1–R7 の後。
  - 化学非平衡・凝縮との組合せ (等温壁 × 凝縮 ON は case/44 の凝縮 restart 手順で後続)。

## 3. 関連 docs と前提 (既存資産の照合)

| 部品 | 既存資産 | 状態・ギャップ |
| --- | --- | --- |
| 等温壁 BC | `wall_isothermal` (floats `Ts`, `Ux/Uy/Uz`)。cell: ゴースト $T_R = 2T_w - T_L$ (2026-07-20 修正、純伝導厳密解 +0.02 %)。node: 壁ノード温度ピン + `res_roe` 0 化 + DPLUR 行 decouple ([`methods/boundary.md`](../../methods/boundary.md) 「node 等温壁の壁ノード温度ピン」) | **再利用**。case/40 `run_0048` (node y⁺1, $T_s$ 1000 K) で運用実績。y₁⁺≲1 の step1 発散は cross-mesh IC の P 段差が真因で解決済 ([[node-isothermal-wall-thin-cell-mass-source]]) — **同一メッシュ index コピー warm start なら問題なし** |
| 低 Re SST | `wallTreatmentSST: 0` (`_config_sst_node`, ノズル NS の生産設定) | 再利用。平板の既往検証は M 0.2・$T_w$ 320 K (case/26 `run_0023`/`run_0027`: q_w が Colburn と 1.5 %) のみ。**超音速・強冷却は未検証** — 本計画 S2 |
| 壁関数 × 等温 | `sstEnergyWallFunction: 1` (Kader q_w)。case/40 ベル部 (M≈4, $T_w/T_{aw}$≈0.4) で **+87 % 過大** | **使わない** (§2) |
| 積分法初期壁 | `feedback/deltastar_integral.py` (`thermal_bc.mode: adiabatic | prescribed_temperature`, `Tw` / `Tw_table`)。V3 単体試験のみ | 再利用。**NS 側と同じ $T_w$ を渡す配管が無い** (現状 `deltastar_loop` の既定は断熱固定) — S1 |
| δ\* 抽出 | `metrics/deltastar.py::deltastar_from_run` (ρu 質量収支、帯局所参照) | 熱境界条件に依らない定義なので**そのまま使える**。冷却で $\rho u$ 欠損が減り δ\* が小さく (強冷却では負にも) なるのは物理 |
| bcond 生成 | `runner_wt._bcond(p, euler)`: `wall_kind = "slip" if euler else "wall"` 固定。`runner_sern.py` L169 も同形 | S1 で `p.spec.wall_thermal` を読み `wall_isothermal` + `Ts` を書く |
| メッシュ y⁺ | `wall_first_frac` (CPG 6.5e-5 = y⁺≈2 断熱、TP 4.5e-5 = y⁺≈1.4)。AR ≤ 1000 ゲート | **冷却で y⁺ は上がる** (§4.2)。第一セルを詰めると AR が上がるので、平板の y₁⁺ 掃引で許容 y₁⁺ を決めてからノズルの `wall_first_frac`/`ni` を選ぶ |
| SU2 | `.external/su2/bin/SU2_CFD` v8.5、平板 (case/26 `run_0049`) とノズル CPG (case/45 `run_0012`) の SST cfg・msh→su2 変換 ([`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md)) | 再利用。等温は `MARKER_ISOTHERMAL= ( wall, 300.0 )`。**素 SST (`dilatationCorrection: 0, katoLaunder: 0`) で比較する** (case/45 `run_0013`: 素 SST で δ99 ≤3 %・θ ≤0.4 %・δ\* ≤1.3 % 一致、`dilatationCorrection: 2` は BL −16 % のモデル形式差) |
| 後処理 | `tools/flatplate_bl.py` (θ/δ\*/K–S $C_f$/壁法則)、case/26 `tools/cf_node.py`・`cf_retheta_analysis.py`、case/45 `compare_bl_su2.py` | 平板は**圧縮性版** (van Driest II, 回復温度, Crocco–Busemann, van Driest 変換) を `case/48/tools/` に新設し `flatplate_bl.py` の抽出を流用 |

## 4. 設計方針

### 4.1 問題定義: `spec.wall_thermal` を単一ソースにする

```yaml
spec:
  wall_thermal: {mode: isothermal, Tw: 300.0}      # 既定 (省略時) = {mode: adiabatic}
```

- `runner_wt._bcond(p, euler)`: `mode == isothermal` かつ NS のとき `wall: {kind: wall_isothermal, floats: {Ux: 0, Uy: 0, Uz: 0, Ts: <Tw>}}`。Euler は従来どおり `slip`。`runner_sern.py` の壁行も同じヘルパを通す。
- `feedback/deltastar_loop.run_pass0_integral`: `thermal_bc` は**常に** `spec.wall_thermal` から作る (`isothermal` → `{"mode": "prescribed_temperature", "Tw": Tw}`)。`--init-thermal` の上書きは**廃止** (codex m1: NS と積分法が別の壁温を読む状態を作らない。不一致の初期化実験が要るなら明示の例外として run に記録する)。
- `prepare_ns` の `prepare_info.json` に `wall_thermal` を記録し、`collect` の metrics に **壁熱流束の積分 $Q_w=\int q_w\,2\pi r\,ds$ ($ds=\sqrt{1+(dr/dx)^2}\,dx$; 平面は単位幅 $\int q_w\,ds$) と $q_w(x)$ のピーク位置** を追加 (等温のときのみ。断熱は 0)。
- **低 Re の $q_w$ は解像勾配から後処理で取る** (codex M2 採用): 既存の bvar `qwall` は壁関数経路専用で低 Re では 0 のまま (case/40 `run_0048` の `res_wall_3_12000.h5` で全点 0 を確認)。node では壁ノード $T_w$ と壁法線方向の第 1・第 2 内点から 2 次片側差分で $\partial T/\partial n$ を作り $q_w = -\lambda_w \partial T/\partial n$ ($\lambda_w = \mu(T_w) c_p/Pr$)。符号は壁へ入る向きを正。**閉合検証**: 平板で $\int q_w\,ds$ と入口・出口の全エンタルピー流束差 (node の壁 Dirichlet が落とす壁エネルギー残差込み) を突き合わせ 5 % 以内。
- **符号付き δ\*** (codex M1 採用): 生産抽出器 (`metrics/deltastar.py` L379 `negative_deficit` hard 不合格、L443 平滑化 `positive=True` が負値を 0 に丸める) は冷却壁で破綻する。$T_w < T_e$ の区間 (ノズルのチャンバ〜スロート: $T_e$≈1000 K に対し $T_w$ 300 K) では $\rho_w/\rho_e$≈3 で**質量欠損は負** (δ\* < 0 = 壁が実効的に外へ動く) が物理。S1 で抽出ゲート・等価半径変換・平滑化 (`positive=False`)・壁更新を符号付きに拡張し、正/零/負を横断する回帰 (合成プロファイル) を追加する。
- 変換器 (`convertGmshToForge`) の wall_dist は `wall`/`wall_isothermal` 両方の bcond を壁とみなすので変更不要 (`runner_wt.py` L210 の注記どおり)。

### 4.2 冷却壁と y⁺ (メッシュ要件)

壁単位の $y^+ = y_1\sqrt{\rho_w\tau_w}/\mu_w$ は、同じ第一セル高さ $y_1$ に対して**冷却で上がる** ($\rho_w \propto 1/T_w$ で増え、$\mu_w \propto T_w^{0.7}$ で減る。$\tau_w$ も冷却で増える)。
ノズル出口相当 (M 4.19, $T_e$ 283 K, $T_{aw}$≈1170 K) で $T_w$ 300 K なら $\rho_w$ ×3.9、$\mu_w$ ×0.39、$\tau_w$ ×1.2〜1.4 の見積りで **$y^+$ は $\sqrt{3.9\times(1.2\text{〜}1.4)}/0.39$ = ×5.5〜6** (codex m2 で算術を訂正。実測は case/48 run A の y₁⁺ 0.10 [3 µm, 断熱] と run B で確定する) — 断熱で y⁺1 のメッシュは冷却壁で y⁺5〜6 になり low-Re SST の前提を外れる。
対策は第一セルを詰めることだが AR ゲート (≤1000) と競合するので、**S2 の平板で冷却壁の y₁⁺ 掃引 (0.5 / 1 / 2 / 4) を取り、δ\*・q_w が y₁⁺ に依らなくなる上限を実測してからノズルの `wall_first_frac` と `ni` を決める**。
run ごとに実測 y₁⁺ (壁 $\tau_w, \rho_w, \mu_w$ から) を台帳に出し、上限超えは `SUSPECT` にする。

### 4.3 検証 1: 超音速冷却壁平板 (`case/48.flat_plate_cooled_m4`)

**条件** = case/44 va3 風洞の試験部壁を模す: 空気 CPG (γ 1.4, R 287, Sutherland, Pr 0.72, Pr_t 0.9)、$M_e$ 4.19, $P_e$ 5037 Pa, $T_e$ 283 K ($T_{aw} = T_e(1 + r\frac{\gamma-1}{2}M_e^2)$, $r = Pr^{1/3}$ → ≈1170 K)、$Re/m ≈ 5\times10^6$、平板長 1 m ($Re_L$ 5e6, $Re_\theta$ 帯 ≈ 2000–6000)。
メッシュ = case/26 `flat_plate_planar.geo` を派生 (平面 2D、node。押し出し 2 ノードは MUSCL 散逸消滅で不可 [[node-2node-spanwise-muscl-zero-dissipation]])、上流 slip 助走 0.1 m、上面 slip (前縁波の反射は $x$≈1.6 m で板外)、入口 = 超音速一様 Dirichlet、出口 = `outlet_statPress` ($P_s = P_e$ 一致)。
壁法線は冷却壁で y₁⁺ 0.5 になる第一セル (≈1 µm) を基準に、y₁⁺ 1 / 2 / 4 の粗化メッシュを同じ `.geo` のパラメータで生成。

**run 行列** (すべて node・低 Re SST・SLAU・MUSCL 2 次・陰解法; 段階起動は [`procedures/divergence-and-startup.md`](../../procedures/divergence-and-startup.md)):

| run | 壁 | 目的 |
| --- | --- | --- |
| A | `wall` (断熱) | 基準。$T_w = T_{aw}$ の実測 (回復係数) |
| B | `wall_isothermal` $T_w$ 300 K ($T_w/T_{aw}$ 0.26) | 主対象 (強冷却) |
| C | `wall_isothermal` $T_w$ 700 K (≈0.6) | 中間点 — $T_w$ 依存の傾きを 3 点で取る |
| B-y⁺ | B を y₁⁺ 1 / 2 / 4 メッシュで | §4.2 の許容 y₁⁺ 決定 |
| A-plain / B-plain | A / B で素 SST (下の対応表) | SU2-A / SU2-B の**基準対** (codex M7: 素 SST は A/B 両方に要る)。生産設定 (dilatation 2) との差は別途 |
| SU2-A / SU2-B | 同一 `.su2`、SST-2003m、`MARKER_HEATFLUX 0` / `MARKER_ISOTHERMAL 300` | コード間比較 |

**SST 条件の対応表 (codex M7 採用: 2 キーだけでは「同じ SST」にならない)**: forge 素 SST = `dilatationCorrection: 0, katoLaunder: 0, sstOmegaProdFromPk: 0, sstSigmaBlend: 0, sstEnergyIncludesK: 0, sstNodeWallKPin: 1` (2026-09-08 の既定変更 2 件を明示的に旧値へ)、`wallTreatmentSST: 0`、乱流輸送 1 次風上 ↔ SU2 `KIND_TURB_MODEL= SST` + `SST_OPTIONS= V2003m`, `MUSCL_TURB= NO`。物性: Sutherland ($\mu_0$ 1.716e-5, $T_0$ 273.0, $S$ 111.0 — forge の定数に SU2 を合わせる)、$Pr$ 0.72 (constant-Pr 伝導)、$Pr_t$ 0.9。入口 k/ω は同じ $(k_\infty, \omega_\infty)$ = (75 m²/s², 26000 1/s) (SU2 は TI 0.5 % / $\mu_t/\mu$ 10.1 で同値)。**流れ方向の格子感度** (nx 1000 vs 1500) を 1〜3 % 比較の前提として 1 回取る。

**理論・経験式 (何を再現できれば合格か)**:

1. **摩擦 $C_f(Re_\theta; M_e, T_w/T_{aw})$ = van Driest II** (Hopkins & Inouye 1971 が冷却壁データで最良と評価した標準形、散布 ±10 %)。非圧縮基準は case/26 で採用済みの Kármán–Schoenherr、変換は
   $F_c = (T_{aw}/T_e - 1)/(\sin^{-1}\alpha + \sin^{-1}\beta)^2$, $F_\theta = \mu_e/\mu_w$, $C_f = C_{f,i}(F_\theta Re_\theta)/F_c$
   ($\alpha, \beta$ は $T_w/T_e, T_{aw}/T_e$ の標準式)。判定は 2 段: (a) 絶対値 ±10 % (VD-II 自身の散布)、(b) **比 $C_f^{B}/C_f^{A}$ が VD-II の比と ±5 %** (forge 固有の −6 % 級バイアス [[reichardt-5pct-gap-not-forge]] は比で相殺。壁温影響そのものの検証)。
2. **熱流束**: $St = q_w/[\rho_e u_e c_p (T_{aw} - T_w)]$、Reynolds アナロジー係数 $2St/C_f$ が 1.0〜1.2 (Chi–Spalding 1.16、Colburn $Pr^{-2/3}$ 1.24 を上限側) に入る。$T_{aw}$ は run A の実測壁温を使う (回復係数の検証を兼ねる: $r$ = 0.88–0.90)。
3. **積分厚さ $\delta^*, \theta, H$**: (a) **CONTUR 積分法** (`deltastar_integral.py` の平面極限 $r_w \to \infty$、同じ $T_w$) の $\delta^*(x), \theta(x)$ と ±15 % (積分法の精度)。**比 $\delta^{*B}/\delta^{*A}$ は ±10 %** — δ\* チェーンの初期壁がそのまま冷却壁でも使える証拠になる。(b) **温度–速度関係は診断** (codex M5 採用: 古典 Walz/CB 形は CONTUR の閉包と同じ式なので独立参照にならず、冷却壁では Duan–Martín 型 [線形項係数 $C_T$=0.8259, Chen–Gan–Fu JFM 2025 式 1.1a] と $u/u_e$=0.5 で約 8 % 違う)。forge の $T(u)$ を Walz 形と Duan–Martín 形の両方と重ね、**どちらに近いか**を記録する (合否にしない)。CONTUR 平面版は巨大 $r_w$ の流用ではなく、一定外縁条件で $d\theta/dx = C_f/2$ を積分する平板入口 (`flat_plate_integral`) を用意する。(c) $H$ の圧縮性関係 (Walz): $H = H_i T_w/T_e + (T_{aw}/T_e - 1)\cdot(\ldots)$ を CB 求積で作った値と比較 (診断)。
4. **速度分布**: van Driest 変換 $u^+_{VD}$ が対数則 (κ 0.41, B 5.0) に乗るか (診断。強冷却では Trettel–Larsson 変換のほうが良いことが知られており、両方を図示するが合否にしない)。
5. **SU2 同一メッシュ**: 素 SST 同士で $C_f(x)$・$q_w(x)$・$\delta^*(x)$・$\theta(x)$ が **3 % 以内** (x = 0.3–0.9 m 平均。case/45 断熱の実績 δ\* ≤1.3 %)。生産設定 (dilatation 2) との差は「モデル形式差」として台帳に別掲。
6. **ゲート** (codex M3 採用: 「全列 falling/flat」だけでは高い残差の停滞も通るので定量化する): (i) `check_convergence.py` PASS、または warm 床のときは**全列の最終残差が断熱基準 run A の床以下**かつ NaN 0 (元の `NOT CONVERGED` 表示は残す)。(ii) 報告量の時系列判定は既存の `check_quasisteady.py --quantity theta,cf_retheta` (平板専用) に加え、本 case の `tools/cooled_plate_eval.py --series` で $q_w$・$Q_w$・δ\*・θ・$C_f$ の全スナップショット時系列を出し、**末尾 50 % の drift が比較公差の 1/3 以下** (δ\*/θ: 1 %、$q_w$: 1.5 %、$C_f$: 1.5 %) を STEADY とする。(iii) `check_mesh_quality.py` PASS、(iv) 実測 y₁⁺ ≤ 掃引で決めた上限 (超過は生産ゲートで不合格)、(v) **起動ゲート**: step 0/1 の壁 P・ρ・T・ω が有限で、壁 P が自由流の 0.5〜2 倍 (codex M6)。

### 4.4 検証 2: ノズル形状 (CPG) × SU2 等温

case/45 CPG チェーン (`problem_d155_cpg_ns.yaml`, forge `run_0013_cpg_ns_plainsst` ↔ SU2 `run_0012_su2_sst` が断熱で一致済) に `wall_thermal: isothermal 300 K` を足し、
forge (素 SST, `run_0013` から index コピー warm start) と SU2 (`MARKER_ISOTHERMAL`, `run_0012` から restart) を同じメッシュで回す。
メッシュは §4.2 の結果で決める: 現行 y⁺≈2 (断熱) は冷却で y⁺≈10 になるので、**平板の掃引で許容 y₁⁺ を超えるなら `wall_first_frac` を詰めた新メッシュを両者に使う** (AR 超過時は `ni` を増やす)。
**引き継ぎ手順は 3 種を分ける** (codex M6 採用: 「同一 index なら安全」は温度変更と再メッシュに当てはまらない): (a) **同一メッシュで壁温変更** — 温度ピンは ρ を保って $P=\rho R T_w$ を再設定するので 1170→300 K の切替直後に壁圧が 0.26 倍に落ちる過渡が出る。切替後は必ず soft 段 (1 次, cfl 0.5, 2000 step) を挟み、起動ゲート (step 0/1 の壁 P/ρ/T/ω) を通す。必要なら中間温度 (700 K) を経由する。(b) **同一トポロジで形状変更** (δ\* 反復の pass 間) — index コピー。(c) **解像度変更** (`wall_first_frac`/`ni`) — `interp_field.py` (原始変数補間) + soft 段。SU2 側は `SU2_SOL`/内蔵補間でなく新メッシュで cold start (段階 CFL) し、両者の収束を独立に確認する。
比較量 = `compare_bl_su2.py` の 4 ステーション (x/r_t 40/60/80/94) の δ99/δ\*/θ に $q_w(x)$ と $Q_w$ を追加。判定: δ\* ≤3 %・θ ≤1 %・$q_w$ ≤5 % (等温は温度場が新たに効くので断熱より緩める)。

### 4.5 生産: TP 風洞 (case/44 va3) の等温 δ\* 反復

2 段に分ける (codex M4 採用: 等温で初期壁を作り直すと形状補正が壁温影響を打ち消し、「壁温だけの感度」にならない):

- **S4a 固定形状の壁温感度**: 断熱の生産 run `run_0107` (case/44 V4 pass 0) と**同じ物理壁・同じメッシュ**で `wall_isothermal` 300 K を回す (手順 §4.4 (a): 同一メッシュ + soft 段)。台帳 = δ\*_exit / ṁ 比 / 出口面コア M / 軸 M 波 / **出口の全温・全圧分布** (壁熱移動で全エンタルピー・全圧が変わる分。「風洞では δ\* 経由でしか効かない」は撤回) / $Q_w$ / $q_w$ ピーク位置 / 実測 y₁⁺。これが §4.6 の「感度ブラケット」の実測。
- **S4b 300 K での再設計**: `problem_va3_M4.19_Lc8_iso300.yaml` (`wall_thermal` 300 K) で生産レシピ [[deltastar-production-recipe]] を回す: pass 0 = 積分法初期壁 (`thermal_bc` = 300 K、符号付き δ\*) + NS (Euler 基準 `run_0091`) → pass 1 (ω=1.0)。ゲートは同じ (|ṁ_NS/ṁ_E − 1| ≤ 0.3 %、出口面コア M ±0.1 %)。台帳には **S4a の感度と S4b の回復量を別項目**で記す。
- 積分法 (`Tw` 指定) は候補選別の近似にとどめ、最終採否はエネルギー込みの NS 評価に置く。

### 4.6 最適設計で壁温影響をどう評価するか (方針の提案)

壁温は目的関数に **(i) δ\* (実効輪郭 → 出口 M・一様性・推力係数)、(ii) $C_f$ (摩擦損失 → SERN の $C_T$)、(iii) $q_w$ (熱負荷 → 冷却設計の制約)** の 3 経路で効く。
一方 $T_w$ は設計者が選ぶ量ではなく**冷却方式と運転で決まる環境量** (断熱 ↔ 300 K の間のどこか、しかも $x$ 分布)。したがって:

1. **$T_w$ は dv にしない。作動点と同じ「環境シナリオ」として扱う** (SERN の多作動点重み付けと同じ枠組み)。既定シナリオ = {断熱, 等温 $T_w^{\rm nom}$ (冷却設計の公称値。水冷銅なら 300–400 K)}。
2. **まず感度ブラケットを取る**: 公称形状 1 点で断熱と $T_w^{\rm nom}$ の 2 run を回し、目的量の差 $\Delta f = f(T_w^{\rm nom}) - f({\rm ad})$ を台帳に出す (§4.5 が風洞の実測)。
   - $|\Delta f|$ が設計公差より小さい (風洞: 出口コア M ±0.1 %・ṁ 比 0.3 %; SERN: $\Delta C_T$ < 0.002 [SERN plan R6(c)]) → **断熱で設計し、台帳に壁温ロバスト性を記す**だけでよい。
   - 大きい → **公称 $T_w$ で設計する** (本計画の等温 δ\* 反復)。オフノミナル ($T_w$ の不確かさ ±ΔT) は同じブラケットで再評価し、公差内なら終了。
   - それでも公差外 (壁温不確かさが目的量を支配する) → **ロバスト設計**: 既存 MOO の目的ベクトルに $T_w$ シナリオを作動点として束ね (期待値 or 最悪値)、パレートを取る。
3. **δ\* 経路の近似評価が使える場面**: 風洞の一様性目的は主に δ\* 経由で壁温を見る (熱負荷は制約側) が、壁熱移動は出口の全温・全圧分布も変える (S4a で実測)。δ\* の $T_w$ 依存は積分法 (CONTUR, `Tw` 指定) で NS なしに見積れるので、**MOO の内側では積分法で $T_w$ 感度を先に篩い、NS は採用点だけ**にする (S2 で積分法の比 δ\*_cold/δ\*_ad が ±10 % で当たることを確認するのはこのため)。最終採否は NS。
4. **CHT が要るのはいつか**: $T_w(x)$ が未知で、しかも δ\*・$q_w$ が $T_w(x)$ の分布形に敏感なとき (薄肉・再生冷却・局所ホットスポット)。その場合も**フル CHT (固体伝導ソルバの連成) の前に「弱 CHT ループ」** — NS の $q_w(x)$ → 1D 壁伝導 + 冷却剤熱伝達モデル → $T_w(x)$ → `wallProfile` で NS 再実行 — を推奨する。等温壁機構と積分法の `Tw_table` がそのまま使え、固体側は解析式なので実装コストが小さい。フル CHT へ進む判断は**反復の収束/不収束では決めない** (反復が収束しないことからモデル不足は判定できない — codex 2026-09-19)。判断基準は [boundary-conjugate-heat-transfer](boundary-conjugate-heat-transfer.md) §4.8 の**モデル感度と厚さ方向近似の評価** (`local1d` vs `shell2d` の差、$\mathrm{Bi}_t=ht/k_s>0.1$、角部の熱橋) に置く。
5. **台帳の標準項目** (全チェーン共通): `wall_thermal` / 実測 y₁⁺ (max) / $Q_w$ / $q_w$ ピーク位置 / δ\*_exit / 目的量の断熱比。

### 4.7 CONTUR の熱力学・輸送の整合 (第 1 層、2026-10-08 起案)

背景と諮問: §5.1 #10 (壁温が頻繁に変わる前提、§8-4)。codex diagnose 2026-10-08 ([記録](../../notes/reviews/2026-10-08-contur-property-temperature-diagnose.md)) の 2 層分離に従う。
第 1 層は「NS と同じ気体・同じ輸送物性で書く」整合で、選択の余地を持たせない。第 2 層 (圧縮性変換の参照温度、h(v) の 2 次分布・r・a・N の閉包) は今の van Driest II 型を基準に据え置く。
§5.1 #10a の A/B で、熱閉包の書き方だけで冷却壁の δ_r が最大 2.8 % 動くことを確かめた (第 1 仮説を支持)。

1. **熱閉包をエンタルピーで書く**: $h(v) = h_w + a(h_{aw}-h_w)v + [h_e - a(h_{aw}-h_w) - h_w]v^2$、$T = h^{-1}(h)$、$\rho/\rho_e = T_e/T$ (境界層内で圧力・組成が一定)。
   $h(T)$ は縁の状態と同じ気体 (semi-perfect なら NASA-9 の $c_p(T)$ の積分)。CPG では今の温度形と式として一致する。
2. **断熱壁の回復をエンタルピーで**: $h_{aw} = h_e + r(h_0 - h_e)$、$r = 0.72^{1/3}$。r は「NS の断熱壁温に合う回復の閉包」として記録し、分子 Pr とは呼ばない (#10d)。
3. **粘性を NS と同じにする** (#10c): `gas.transport` があれば NS の viscMethod 2 と同じ CEA の種ごとの粘性と混合則 (データは `solver_density_cuda/data/species/forge_transport_v1.yaml`)。
   TP で `gas.transport` が無ければ止める。CPG と明示した旧方式だけ Sutherland。効く場所は μ_e・μ_w (F_Rδ)・N の R_δ・入口 θ₀・R_θi の床。
4. **運動量式の加速の項を気体に整合させる** (#10b): Eq. 61 の $(2 - M^2 + H)/(M(1 + (\gamma-1)M^2/2))\,dM/dx$ を
   $(2 + H - M^2)\,d\ln u_e/dx$ に置き換える (等エントロピー・組成一定なら $d\ln\rho_e = -M^2 d\ln u_e$)。
   上流の亜音速枝は面積から $d\ln u_e/dx = -(2 r_w'/r_w)/(1 - M^2)$ (γ を使わない)、下流は縁の $u_e(M)$ を気体から作って微分する。CPG では今と一致する。
5. **版で切り替える**: `deltastar_initializer.closure_version` (`contur_v1` = 今 / `contur_v2` = 1〜4)。既定は当面 `contur_v1` のまま
   (今の生産の壁をビット同一で再現できる)。既定を `contur_v2` に替えるのは §6 V-c45 で熱閉包の形が NS で裏付けられた後に、ユーザの決定で行う
   (そのとき生産の YAML は `contur_v1` を明示して再現性を残す)。版は `delta_r_initial.json` の settings と prepare_info に残す。
6. **較正の流用の制限** (#11) は別項目。`contur_v2` で較正し直した k_f は `contur_v1` の k_f と混ぜない。

## 5. 実装ステップ

1. **S1 配管** — `design/forge_design/probdef.py` (`spec.wall_thermal` 既定・検証)、`evaluate/runner_wt.py::_bcond` (+ 共通ヘルパ `wall_bcond_line`)、`evaluate/runner_sern.py` L169/L181、`feedback/deltastar_loop.py` (initializer 既定を spec から)、`evaluate/runner_axismach.py::prepare_ns/collect` (帳簿: `wall_thermal`, y₁⁺, $Q_w$)。
2. **S2 平板** — `case/48.flat_plate_cooled_m4/`: `mesh/flat_plate_cooled.geo` (パラメータ化 y₁)、`gen_runs.py` (config/bcond/IC 生成 + 段階起動)、`tools/cooled_plate_eval.py` (VD-II / RAF / CONTUR 平面 / CB / VD 変換 / SU2 読み込み)、SU2 cfg。README に run 一覧。
3. **S3 ノズル CPG × SU2** — case/45 に `problem_d155_cpg_ns_iso300.yaml`、forge run + SU2 run、`compare_bl_su2.py` に $q_w$ 追加。
4. **S4 生産 TP 等温 δ\* 反復** — case/44 `problem_va3_M4.19_Lc8_iso300.yaml`、`deltastar_loop` pass 0/1、壁温感度台帳 (`case/44/README.md`)。
5. **S5 docs** — `methods/design/overview.md` 節の実装同期、`procedures/recommended-settings.md` の壁行に等温レシピ (y₁⁺ 要件)、`procedures/verification/48-flat-plate-cooled.md`、`design/CAPABILITIES.md`。

### 5.1 残作業 (優先順)

**残作業の正本はこの表**。`notes/sessions/` の引き継ぎ文書には写しとポインタだけを置く。

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~S1 配管~~ **実装済 (2026-09-12)**: `Problem.wall_thermal` / `wall_thermal_bc_integral` / `wall_bcond_line`、`runner_wt._bcond`・`runner_sern` の壁行、`prepare_ns` (thermal_bc を spec から強制・`prepare_info.wall_thermal`)、`deltastar_loop` (`--init-thermal` 廃止)。回帰: Euler `run_0037`/`run_0091`・NS `run_0107` の壁行が文字列一致 | §5-1。`wall_thermal` 未指定は断熱で**ビット同一** (回帰: Euler `run_0037`/`run_0091` に加え **NS** の `run_0107` の bcond 文字列・積分法初期壁が一致すること [codex m1]) |
| 1a | ~~**符号付き δ\*** (codex M1)~~ **実装済 (2026-09-12)**: `_negative_delta_r` (参照 ρu を壁値で外へ延長し $r_{eff}^2 = r_w^2 - D/(\pi q_w)$)、`negative_deficit` を soft 化、P-spline `positive` を `wall_thermal` で切替 (断熱は従来どおり ≥0)。`design/tests/run_deltastar_tests.py` に正/零/負回帰 (負 δ_r 復元 2〜3 %) ALL PASS | `metrics/deltastar.py` の `negative_deficit` ゲート・等価半径変換・平滑化 `positive`・`feedback/deltastar_loop` の壁更新を符号付きに。正/零/負の合成プロファイル回帰 |
| 1b | ~~**低 Re 熱流束の後処理 + 閉合** (codex M2)~~ **平板側は実装済** (`cooled_plate_eval.py`: 2 次片側差分 $q_w$、閉合 B 1.017 / C 1.018)。残 = ノズル軸対称 ($ds$ 重み) の `metrics/extract.py` | `case/48/tools/cooled_plate_eval.py` (平板) と `metrics/extract.py` (ノズル軸対称 $ds$ 重み) に $q_w$ 抽出 (2 次片側差分)。平板で全エンタルピー流束差との閉合 5 % |
| 1c | ~~**定量ゲートの実装** (codex M3)~~ **平板側は実装済** (`cooled_plate_eval.py --series`: 末尾 50 % drift, δ\*/θ 1 %・$q_w$/$C_f$ 1.5 %)。残 = ノズル用 (出口コア M・ṁ 比・$Q_w$ の時系列) と起動ゲート | `cooled_plate_eval.py --series` (末尾 50 % drift 公差)、warm 床の定量条件、起動ゲート (step 0/1 壁 P/ρ/T/ω) |
| 2 | S2 平板 — **A/B/C 完了 (2026-09-12, case/48 run_0004/0005/0006)**: B ($T_w$ 300 K): $C_f$/VD-II 0.97–0.99・$2St/C_f$ 1.16・δ\* CONTUR 比 +3.5〜6 %・θ ±0 %・閉合 1.017・series STEADY。C (700 K): 0.91–0.95・1.15・+8〜16 %・±3 %・1.018。A (断熱): 0.87–0.92 (低 $Re_\theta$ 1400–3500)。**比 B/A の $C_f$ は 1.50 (VD-II 1.35–1.39, +8 %: 目標 ±5 % 外。C/A は 1.20 vs 1.16 で内)** — forge の断熱側の不足 (既知の Reichardt 系ギャップ) が冷却で縮む。SU2 の比で切り分ける。**y₁⁺ 掃引完了 (run_0007/0008/0009)**: 冷却壁 y₁⁺ 0.46→0.9/1.8/3.6 で δ\*/θ は −1/−2/−2.2 %、$C_f$ −1/−2.6/−3.1 %、$q_w$ −1.3/−3.7/−7.0 % → **δ\* チェーンは y₁⁺ ≤ 3.6 で 2 % 内 (許容上限 3.5)、熱負荷は y₁⁺ ≤ 1 (1.5 %) / ≤ 2 (4 %)**。**SU2 平板 (中間, it≈3000)**: B の積分量 CD (∫Cf dx) forge/SU2 = 0.998・HF (∫q_w dx) 0.982 (forge 生産 SST vs SU2 V2003m)。**A-plain/B-plain 完了 (run_0010/0011)**: 生産 SST との差は $C_f$ +1〜1.6 %・δ\* +1 % (平板では dilatation 補正は 1〜2 %)。**B-plain vs SU2-B: CD +0.8 % / HF −0.8 %** (積分量; 3 % 目標内)。**SU2-B ステーション比較 (it 5000)**: $C_f$ 1.000、$q_w$ 1.008、θ 0.998、**δ\* +4 %** (3 % 目標をわずかに超過: H 4.66 vs 4.48 の壁近傍密度分布差。$C_f$/θ/$q_w$ が一致しているので離散化でなく温度–速度関係のモデル差)。**SU2-A ステーション比較**: $C_f$ 0.997–1.000, θ 0.98–1.02, δ\* 0.96–1.04 (SU2-A は $T_w$ 1145 K でまだ上昇中)。**nx 感度 (run_0012, nx 1500)**: $C_f$ +0.2 %, δ\*/θ ±0.1 % → 流れ方向は収束。**S2 の残は無し** | §4.3。合否は §6。**起動レシピ (2026-09-12 実測)**: 一様 IC + SST は前縁で step 114 NaN → 層流暖機 (model none, 1 次, cfl 0.2, 2000) → SST soft (1 次, cfl 0.3) → mid (1 次, cfl 1.0) → 2 次ランプ cfl 0.5/1/2 (各 2000) → 本段 cfl 2 (cfl 4 は前縁 x≈1〜3 mm で P 床→NaN)。IC は壁近傍 tanh ランプ (δ₀ 300 µm) |
| 3 | S3 ノズル CPG × SU2 等温 | §4.4。**メッシュ確定 (2026-09-12, `case/44/mesh_probe_iso.py`)**: 冷却壁の y₁⁺ 倍率は run_0107 の壁データから 5.0–5.4 (スロート y₁⁺ 7→37)。`mesh2d` に x 依存の第一セル (`wall_first_frac_throat` + 上下流ブレンド) を追加し、**ni 4001 / nj 113 / wff 1.8e-5 / wfft 2.5e-6 / throat_refine 30 で AR max 846 PASS、y₁⁺_cold 1.5–4.5 (452k 節点)**。AR ≤ 1000 を守ると ni ≥ 3600 が要る (y₁⁺ 1 にするには ni 20000 級で不可: §8-3)。対象は case/45 でなく **case/44 va (Tt 1060, M4.19)** に変更 (SU2 の規模; 断熱 SU2 も同メッシュで取り直す = codex M7)。**投入 (2026-09-12)**: `run_iso_chain_cpg.py` — forge `run_0112` (CPG 断熱・素 SST, 完了: 出口コア M 4.1813, y₁⁺ スロート 0.88) → SU2 `su2_va_cpg_ad` (forge 場から restart_flow_in.csv を座標対応で生成 [float32 座標ずれ 2.4e-7 m、対応一意]、`READ_BINARY_RESTART= NO`、6000 it ≈ 2 s/it) 実行中 → forge `run_0113` (300 K 同壁, 完了: 出口コア M 4.1982 = 断熱比 +0.40 %, TP 系と同じ傾き) → SU2 `su2_va_cpg_iso300` (実行中)。**SU2 6000 it 完了・暫定比較 (`compare_nozzle_su2.py`, 2026-09-12)**: 300 K 対 (run_0113 / su2_va_cpg_iso300): ṁ 0.999、出口コア M +0.28 %、q_w ピーク +4 %、局所 q_w +7〜10 %、δ\* +6〜8 %、θ 0〜9 %、τ_w +3〜6 %、δ99 一致 (両者 y₁⁺ 2〜4)。断熱対: ṁ +0.2 %、出口コア M +0.28 %、δ\* −2〜+9 %、τ_w −2 %、ただし**両コードとも壁温が未発達 (出口 T_w 739 / 784 K vs T_aw ≈ 970 K)** = 熱場の緩和が遅い。**SU2 は未収束 (rms[RhoE] 300 K +1.15 / 断熱 −0.89)** → 各 10000 it を継続中 (`_su2_cont.sh`, 順次)。§6 の判定 (δ\* ≤3 %・q_w ≤5 %) は継続後に出す。現時点の +7 % は 3 % 目標外で、y₁⁺ 2〜4 の影響 (平板掃引: y₁⁺ 1.8→3.6 で q_w −4〜−7 %) と SU2 未収束が交絡 |
| 4 | S4a 固定形状の壁温感度 (同一壁・同一メッシュ, 断熱 vs 300 K) → S4b 300 K 再設計 | §4.5 (codex M4)。台帳は感度と回復量を別項目。**投入 (2026-09-12)**: `case/44/run_iso_chain.py` → run_0108 (断熱, 細分メッシュ, run_0107 と同じ壁 `delta_r_initial.csv`) → run_0109 (300 K 同壁) → run_0110 (300 K 積分法初期壁 pass 0, 符号付き δ\*)。**run_0108 完了**: 品質 PASS (AR 846)、実測 y₁⁺ max 2.8 (上流収縮部) / スロート 0.61 / 出口 0.39 (断熱)、ṁ_NS/ṁ_E 0.9988 (run_0107 0.9992)、出口コア M 4.1869 (run_0107 4.1900, −0.07 % = メッシュ感度)、出口 T0 1060.0 K / P0 1.139 MPa (h0 逆算)、`wall_thermal_ledger.py --series` STEADY (M_core drift 0.000 %)。台帳ツール `case/44/wall_thermal_ledger.py` (y₁⁺(x)・q_w(x)・Q_w・出口 T0/P0・ṁ 比・series)。**run_0109 (S4a-B, 300 K 同壁) 完了 (24000 step, 残差 still converging, series: M_core drift 0.07 %・Q_w 9.6 % = 熱場が未静定 → run_0114 で +36000 継続)**。**固定形状の壁温感度 (暫定, run_0108 → run_0109)**: 出口コア M 4.1869 → **4.2071 (+0.48 %)**、ṁ_NS/ṁ_E 0.9988 → 1.0011 (+0.23 %)、出口 T0 コア 1060 K 不変・質量平均 1060 → 1048 K (−1.2 %, BL の熱損失)、P0 コア不変、**Q_w 7.8 MW、q_w ピーク 2.83 MW/m² @ x −0.41 r_t (スロート直前)**、実測 y₁⁺ スロート 3.0 / 出口 1.8 / 上流収縮部 max 10.8。抽出 δ_r (同壁) は冷却で 0.6〜0.75 倍 (出口 0.091 vs 0.127 r_t)、スロートで ≈0 (−2e-5)。**壁温は風洞の出口 M を設計公差 (±0.1 %) の 5 倍動かす → §4.6-2 の判定は「公称 T_w で設計」側**。run_0110 (S4b pass 0): 積分法初期壁 (Tw 300) は δ_r 出口 0.115 (断熱 0.122, −6 %) と冷却効果を過小評価 (CFD 抽出は −29 %) → pass 1 (run_0115, 符号付き抽出) で補正する設計どおり。**run_0114 (run_0109 +36000 step) 完了**: 出口コア M 4.2016 (24k 時 4.2071 → まだ −0.13 % 動く、series drift 0.03 % STEADY)、ṁ 比 1.0024、Q_w 6.40 MW (tail drift 6 %: 熱場はなお緩和中)、出口質量平均 T0 1043.8 K。**固定形状の壁温感度 = 出口 M +0.35〜0.48 %、ṁ +0.35 %、質量平均 T0 −1.5 %** (符号・桁は確定、Q_w は ±10 % 幅)。残差は rms_roOmega が末尾で上昇 (スロート極薄セル 0.5 µm の ω が要因の疑い、要確認) |**run_0110 完了**: ṁ_NS/ṁ_E 1.0017 (ゲート 0.3 % 内)、出口コア M 4.1998 (+0.23 %: ゲート ±0.1 % 外 → pass 1 へ)、Q_w 7.36 MW、q_w ピーク 2.74 MW/m²。帯局所抽出 (符号付き) はスロートで δ_r ≈ 0 (入力 0.00115)、下流で入力の 0.79 倍、上流収縮部で負値 251 点を保持 (hard 不合格 0)。**run_0115 (pass 1, ω 1.0, 36000 step) 完了 — 等温 δ\* 反復の生産ゲート達成**: ṁ_NS/ṁ_E **1.0002** (≤0.3 %)、出口面コア M **4.1909 (+0.02 %, ±0.1 % 内)**、series M_core drift 0.02 % STEADY・Q_w 6.8 MW (drift 8.6 %, 熱場は緩和中)、固定点比 (抽出/入力) 0.94〜0.99 (x ≥ 8)、スロート近傍はさらに縮む (x=1 で 0.38、スロート −3e-4)。→ **S4b 完了 (pass 1 で固定点近傍)**。残 = 熱負荷 Q_w の静定 (継続 run) と S3 (SU2 対) |
| 5 | S5 docs・CAPABILITIES | §5-5 |
| 9 | ~~**AR 緩和の A/B**~~ **完了 (2026-09-12)**: AR 4140 メッシュ (run_0116 断熱 / run_0117 300 K, 24000 step) vs AR 846 (run_0108 / run_0109, 同 step): 断熱 ṁ 比 −0.002 %・出口コア M −0.02 %・T0 同一; 300 K ṁ 比 −0.02 %・出口コア M +0.02 %・Q_w +10 % (7.8→8.6 MW)・q_w ピーク +8 % (y₁⁺ スロート 3.0→1.04, 出口 1.8→0.81 の壁解像向上分で、平板掃引の y₁⁺ 依存 [3.6→0.5 で q_w +7 %] と整合)。残差は全列 falling、発散なし。**→ 壁法線構造層の AR ≤ 5000 は平均流に無害 (0.02 %) で熱負荷はむしろ改善** (§8-3 確定) | `problem_va_R2_LU6_Lc8_ns_ar5k{,_iso300}.yaml` (ni 2401 / nj 121 / wff 8e-6 / wfft 8e-7, AR max 4140, スロート y₁⁺_cold ≈ 1.0, 290k 節点) で run_0116 (断熱, 同壁, IC run_0108) → run_0117 (300 K 同壁) → run_0118 (pass 0)。判定: run_0108/0114 との差 (出口 M ≤0.1 %, ṁ ≤0.1 %, δ_r ≤5 %, Q_w ≤10 %) と発散/残差床の有無 | `run_iso_chain_ar5k.py` (2026-09-12 投入) |
| 6 | ~~`wallProfile` CSV ($T_w(x)$ 分布の NS 入力)~~ **[boundary-conjugate-heat-transfer](boundary-conjugate-heat-transfer.md) へ移管 (2026-09-19)** | 同 plan §4.5 / §5.1 #4。`applyInletProfiles` を一般化して `ints: {wallProfile: 1}` で `wall_isothermal` の per-face `Ts` を埋める (Ts は `valueTypes==1` なのでカーネル無改修で効く) |
| 7 | ~~弱 CHT ループ (1D 壁伝導 + 冷却剤モデル)~~ **[boundary-conjugate-heat-transfer](boundary-conjugate-heat-transfer.md) へ移管 (2026-09-19)** | §4.6-4 の方針を同 plan が引き取り、Phase 1 (外部弱連成) → Phase 2 (ソルバ内薄肉シェル) の 2 段に具体化した。**種別は `wall_isothermal` のままで `ints: {conjugate: 1}` 属性**とする (新種別 `wall_conjugate` は codex レビューで撤回: `iso_wall_flag`・T ピン・粘性壁・壁距離・DPLUR エネルギー行切離しの 5 経路から漏れるため) |
| 8 | SERN (⑤) の等温評価 | 配管は #1 で共通化。検証は SERN plan R1–R7 の後 |
| 10 | **CONTUR の物性をどの温度・どのモデルで評価するか (2026-10-08 起票、担当 F)** — 壁温が頻繁に変わる前提 (§8-4) で、積分法の初期壁が壁温に正しく追随するようにする。**判断: 2026-10-08 codex diagnose ([記録](../../notes/reviews/2026-10-08-contur-property-temperature-diagnose.md)) — 「2 層分離は修正付き採用。第 1 層 = 熱力学・輸送の整合 (同じ組成・h(T)・EOS・輸送物性)、第 2 層 = 閉包 (h(u) の 2 次分布・回復係数 r・a・N) と適用範囲。まず熱閉包の表現だけを温度形→エンタルピー形に替える CFD 0 step の A/B (#10a)。#8d の再較正と第 2 層の方式選定は保留。第 2 層は当面 van Driest II 型を基準に維持」**。指摘 10 件 (Major 9・Minor 1) の採否は下の #10a〜#10f と #11・#12 に書いた。§4.7 は #10a の結果を見てから書く | 実測 (case/45、CFD 0 step、[verification-m6-axis-wave-mesh-su2](verification-m6-axis-wave-mesh-su2.md) §9 の 2026-10-08): (i) 断熱壁温 T_aw = T_e(1 + r(γ_e − 1)/2 M²) が局所 γ_e のため全温を超える (試験部 1654〜1663 K、Tt 1600 K、NS 1470〜1482 K)。エンタルピー形 h_aw = h_e + r(h_0 − h_e) (r = 0.72^{1/3}) なら NS と +1.5〜4 K。(ii) T_aw の式の違いで δ_r が 3〜5.5 % 動き (等温でも Eq. 69 の (T_aw − T_w) 項で残る)、現行の式では冷やすと出口 δ_r が厚くなる (+0.95 %、全温基準では −1.1 %)。(iii) μ は空気の Sutherland で、NS (種ごとの CEA) と 250 K で −4.4 %、1470 K で +8.9 % ずれる。案: 第 1 層 = NS と同じ気体での整合 (エンタルピー形の Walz・h_aw・NS と同じ μ(T))、第 2 層 = 圧縮性変換の参照温度を壁温を跨いだ交差検証で選ぶ |
| 10a | ~~**熱閉包の A/B (CFD 0 step、2026-10-08 事前登録、担当 O)** — codex diagnose の判別 A/B~~ **済み (2026-10-08): 第 1 仮説を支持** | A = 今の温度形、B = エンタルピー形: h_aw = h_e + r(h_0 − h_e)、h(v) = h_w + a(h_aw − h_w)v + [h_e − a(h_aw − h_w) − h_w]v²、T(v) = h⁻¹(h(v))。変える因子は熱閉包の表現だけ。r = 0.72^{1/3}・a = 1・k_f = k_N = 1・今の μ・入口 θ₀・縁条件・運動量式・求積は両腕で固定。断熱 / 1000 / 600 / 300 K を入口から出口まで積分。**前提**: CPG 極限で A と B の差 < 0.1 %、積分精度 (rtol 1e-6 と 1e-8) の差 < 0.1 %。**判定 (事前登録)**: 試験部 [40, 94] で、B による 300 K/断熱の δ_r 比の変化 (B の比 / A の比 − 1) の絶対値の最大が 1 % 以上なら第 1 仮説 (今の熱閉包は冷却による δ_r の変化に 1 % 以上の系統差を生む) を支持、1 % 未満なら棄却。NS (断熱) との差が減るかは別欄に記録するだけで、再較正はしない。スクリプト `case/45.isobutane_m6_d155/delta_contur_compare.py hform`。**結果 (2026-10-08、`_band_ab/delta_contur/hform_ab.json`)**: 前提は満たした (CPG 極限の A/B 差 最大 0.020 %、rtol の差 最大 0.0074 %)。**判定: 第 1 仮説を支持** — 300 K/断熱の δ_r 比の変化は最大 2.80 % (x = 40、基準 1 %)。600 K で 1.43 %、1000 K で 0.18 %。温度形 (A) では 300 K に冷やすと試験部で δ_r が +1.0〜+4.8 % 厚くなり、エンタルピー形 (B) では +1.9 % (x = 40) → −1.5 % (x = 94)。どちらが NS に近いかは冷却壁の NS が無いので未確定。記録のみ (NS 断熱、k_f = 1): 出口の δ_E/δ_C は A 1.047 → B 1.032、試験部の振れは 2.93 → 3.17 %。次 (第 1 層の実装) は新しい設計の壁を変えるので、ユーザの了承待ち |
| 10b | **運動量積分式の圧力勾配項の定比熱近似 (codex Major、採用、担当 O、#10a の後)** | Eq. 61 の `(2 − M² + H)/(M(1 + (γ−1)/2 M²)) dM/dx` (`deltastar_integral.py` の `rhs`) は一定 γ の d ln u_e/dM を局所 γ_e で使っている。case/45 の気体で、この式と熱力学的に整合な d ln u_e/dx の比は M 3 で 0.964、M 6 で 0.987 (2026-10-08 手元で再計算、codex の 0.964/0.985 と一致)。`d ln ρ_e/dx + (2 + H) d ln u_e/dx + d ln r_w/dx` で書く。上流の `dMdx` の解析式も同じ近似なので、片側だけ替えない |
| 10c | **μ を NS と同じ輸送モデルに (codex Major、採用、担当 O)** | `gas.transport` があれば混合気の μ(T) (NS の viscMethod 2 と同じ CEA 輸送物性)。TP で `gas.transport` が無ければ止める (`probdef.transport_for_ns` の契約と同じ)。Sutherland は CPG と明示した旧方式だけ。効く場所は μ_e・μ_w・N の Re・入口 θ₀・Re_θi の床。注意: F_Rδ·Re_θc = ρ_e u_e θ_c/μ_w なので、F_Rδ の 11 % 差がそのまま摩擦の誤差にはならない |
| 10d | **r は「NS と同じ Pr」ではない (codex Major、採用・記録)** | viscMethod 2 は `prandtlLam` を使わない (μ と λ を混合輸送から直接計算)。混合気の Pr は 0.749〜0.758。r = 0.72^{1/3} は「NS の断熱壁温に合う回復の閉包」として記録し、分子 Pr とは呼ばない。分子 Pr を NS に近いという理由で 0.72 に置き換えることもしない |
| 10e | **訂正: 「δ の不足と傾きは T_aw の誤差では説明できない」は言い過ぎ (codex Major、採用)** | `tw` は壁温だけを NS 値に替え、温度分布の中の T_aw は今の式のままだった。`taw` も全温基準の近似で、エンタルピー形の分布ではない。言えるのは「壁温だけを替えても傾きは消えない」まで。熱閉包全体は原因の候補に残る ([verification-m6-axis-wave-mesh-su2](verification-m6-axis-wave-mesh-su2.md) §9 の 2026-10-08 に訂正を記録) |
| 10f | **θ の抽出の単位の誤り (codex Major、採用、修正済み 2026-10-08)** | `delta_contur_compare.py extract` が CONTUR の θ (既に θ/r_t) をさらに r_t で割っていた (13.04 倍)。δ 比・C_f 比の集計には使っていない。コードは直した。既存の `extract.npz` の `cont_*_theta_rt` は読むときに r_t を掛けて戻す |
| 11 | **較正の壁温条件の記録とガード (2026-10-08 起票、担当 F)** | `deltastar_initializer` (k_f 等) は較正した壁温条件でしか意味を持たない (case/45 で NS の壁温を与えるだけで出口合わせの k_f が 1.0565 → 1.0682)。較正時の `wall_thermal` を記録し、prepare で `spec.wall_thermal` と照合する案。**判断: 2026-10-08 codex diagnose — 修正付き採用**: 壁温だけの照合は弱く、全 prepare を一律に止めるのは強すぎる。較正係数の**流用と生産採用**を、熱条件・組成と物性 DB・輸送モデル・熱閉包の版・入口 θ₀・検証済みの作動範囲に結びつけて制限する。較正していない初期壁での検証計算はできるようにする。壁温を変えた後の NS 1 回は再評価の開始であって、再較正の完了の保証ではない |
| 12 | **壁温分布の単一ソース (2026-10-08 起票、担当 O、#10 の後)** | `spec.wall_thermal` に分布 (表) を足し、NS の `wallProfile` CSV (ソルバは実装済み、`methods/boundary.md`「壁温分布の入力」) と CONTUR の `Tw_table` の両方をそこから作る。**判断: 2026-10-08 codex diagnose — 採用**: 物理長 [m] (スロート原点) を正本にし、CONTUR に渡すときだけ x/r_t に変換する。同じ条件かの判定は平均温度でなく、座標原点・単位・補間/外挿の規則を含む区分関数の一致で行う。表が範囲を覆わないときに `np.interp` の端値保持で黙って外挿しない。`Tw_table` は局所の T_w を代入するだけで熱境界層の発達の履歴を持たないので、急な冷却の開始や加熱への切り替えでの妥当性は別に確かめる |
| 13 | **冷却壁用の格子の表 (2026-10-08、担当 O)** | `meshing/mesh2d.py` に `wall_first_frac_table` (第一セル厚/局所半径を log 線形の表で) と `x_density_table` (x 方向の相対密度の表) を追加 (どちらも opt-in、表が範囲を覆わなければ例外)。既存の格子はビット同一 (4 構成で確認: 生産の throat/blend・一様・axis_cap・axis_gap)。`runner_axismach._mesh_params_from` が NS の mesh ブロックから渡す。冷却で y1+ は断熱の 8.5〜9.4 倍 (case/45、case/44 の 5〜5.4 倍より大きい) で、smoothstep 1 本の第一セルでは y1+ ≤ 1 と AR ≤ 5000 を同時に満たせなかった (AR 最大の見込み 19000) |
| 14 | **冷却壁の NS の対 (2026-10-08、担当 O、判定は F)** | §6 V-c45。`case/45.isobutane_m6_d155/cold_pair_mesh.py` (格子の設計: 断熱 NS run_0179 の y1+ 分布 × 局所の冷却倍率 → y1+_cold ≤ 0.8・AR ≤ 4500 の見込み、ni 4496 × nj 121)、`cold_pair.py`・`run_cold_pair.sh` (準備・実行)。run: `run_0181_ns_coldmesh_ad`・`run_0182_ns_coldmesh_tw300` |
| 15 | **冷却壁の格子の作り直し (2026-10-08 codex plan M1・M2、担当 O、方式の決定は F)** | 初版 (ni 4496) は品質 FAIL (AR 最大 6014、float32 座標)。倍精度でも傾斜壁の接線長を見ていなかった (縮流部 AR 5488)。接線長を入れた版 (ni 4719) は倍精度で AR 最大 4573 だが、**格子の線が半径方向なので傾斜壁ではスキューがある** (縮流部で最大 0.44)。AGENTS.md の AR ≤ 5000 の例外は「壁法線に沿いスキューの無い層」だけなので、スキューの有る層は AR ≤ 1000。**生産の格子も同じ状態** (スキュー > 0.2 かつ AR > 1000 が 2193 セル、最大 3571、縮流部 x −7.7〜−0.8; 記録)。厳密に守ると x 方向の細分で ni 14653 (177 万節点、スキュー ≤ 0.1 を「無い」とした場合)、縮流部だけなら 106 万節点。代案は近壁の層を壁法線に沿わせる格子 (mesh2d の改修)。どちらにするかは諮問とユーザ判断。y1+ の見積もりは法線距離 (半径方向間隔 × cos θ_w) で行い、投入前に変換後の格子 (FP64) で第一内部点の距離と y1+ の見積もりを検査する |
| 16 | **FP64 のビルド (2026-10-08、諮問中)** | 生産と同じソース (`~/forge-wallfit-bin`、HEAD e2696d8f0) の `flowFormat.hpp` の typedef 4 行だけを double にした `~/forge-wallfit-bin-fp64` (forge sha256 65be5e28…)。float32 の座標では第一セル厚 / 局所半径 3.4e-7 が 3.6 ulp で、第一セル厚が最大 25 % ずれる (倍精度の格子を float32 に丸めて実測)。速度は生産の格子 (19.4 万節点) で FP32 3.44 ms/step → FP64 6.94 ms/step (2.0 倍)。残差の床は rms_ro 3.4e-7 → 7.4e-8 (400 step の比較)。**codex diagnose 2026-10-08 ([記録](../../notes/reviews/2026-10-08-cold-pair-precision-diagnose.md)、全件採用)**: typedef だけでは不十分 — `mesh/gmshReader.hpp` が座標を `stof` で読むので double のビルドでも入力時に float に丸まる → double のビルドだけ `stod` に (float のビルドは従来どおり、c096d66c)。msh の座標は 17 桁 (`mesh.msh_digits: 17`、10 桁だと第一層厚が 0.22 % 動く)。品質は厳密な `VERDICT: PASS` (SOFT-PASS を通さない)。「同じ FP64 なら比の誤差は消える」「生産の格子は 40 ulp なので影響なし」は却下 (生産の第一層厚の量子化 −2.1〜+2.0 %、流れの量への影響は未評価として記録)。**A/B (CFD 0 step、事前登録)**: 同じ FP64 typedef・同じ 17 桁の msh・同じ変換設定で、座標の読み込みだけ A = `stof` (`convertGmshToForge_A_stof`)、B = `stod`。生成時の倍精度座標に対し、B の第一層厚の相対誤差 ≤ 1e-6・非正の層厚なし・品質 `VERDICT: PASS` なら採用。B が超えたら残る切り詰めを追い、冷却 NS は投入しない。スクリプト `case/45.isobutane_m6_d155/fp64_reader_ab.py` |
| 17 | **壁解像の判定を面積で (2026-10-08 codex plan M3、担当 O)** | `check_wall_resolution.py` の超過率・評価率は点数の割合。生産の断熱では点数 3.6 % に対し面積 (2πr ds) では約 0.2 %。境界面の実面積で数え、領域 (縮流部・スロート・試験部) ごとの値と超過位置を出す。V-c45 の壁解像のゲートは面積と領域別で登録し直す |
| 18 | **V-c45 の判定基準の書き直し (2026-10-08 codex plan M4・M5、担当 O、投入前)** | ゲートに両腕の δ_E (x = 40・70・94)・断熱壁温・冷却側の Q_w の定常を加え、`check_quasisteady` の条件 (窓の 5 枚すべて、drift・osc の許容) を固定する。plateau を許す場合も残差床の上限を数値で決める。R_NS = 各腕の 5 枚平均の δ_E の比 (共通の x 点・補間・平滑化を固定)。不確かさは相対で揃え、時間 (5 枚の幅)・抽出 (帯の係数の感度)・CONTUR の積分と表 (≤ 0.1 %)・格子 (生産の格子との断熱の δ_E の差を目安) を合成する。R_NS の区間から e_A・e_B の上下限を作り、**区間が分離したときだけ支持**、1 % 以内の主張も上限で判定 |
| 19 | **熱閉包の A/B と contur_v2 全体を分ける (2026-10-08 codex plan M6、担当 O)** | V-c45 が判定するのは熱閉包の形 (温度形 / エンタルピー形) だけ。k_f は診断用の固定値として両腕に同じ値 (1 と生産の 1.0541) を使い、新版の較正値とは扱わない。`delta_contur_compare.py hform` は `cf_scale` を渡しておらず常に k_f = 1 だった → 生産の k_f でも計算するよう直す。contur_v2 全体 (熱閉包 + 粘性 + 加速の項) の評価と再較正は別に行い、既定を替える条件にする |
| 20 | **contur_v2 の試験 (2026-10-08 codex plan m7・m8、担当 O)** | v1 のビット再現、同じ物性・同じ微分の条件での CPG 極限 (熱閉包の一致)、NS と同じ物性への一致 (μ を独立参照実装と照合) を別の試験にする。Sutherland の定数は今の CONTUR (273.15 / 110.4) と NS (273 / 111) で違うので「CPG + Sutherland で v1 と一致」とは書かない。上流の面積の式は M = 1 で 0/0 になるので音速点の近くの接続を決める。B 腕の積分・表の細分の感度を確かめ、温度・エンタルピーの表の範囲外は例外にする。**実装・試験済み (2026-10-08)**: `feedback/deltastar_integral.py` に `closure_version` (`contur_v1` 既定 = 従来どおり / `contur_v2`)、`gas/transport.py` (NS と同じ CEA の混合気の μ)、`runner_axismach._closure_kw`・`contur_mu_fn` (TP は gas.transport 必須、CPG は NS の Sutherland 1.716e-5/273/111)。`design/tests/run_contur_v2_tests.py` 全件合格: (a) v1 がビット同一、(b1) CPG 極限の式の一致 (同じ点で 1.6e-15)、(b2) CPG 極限の積分の v2 − v1 が v1 自身の縁の格子の誤差以下、(c) 生産の燃焼ガスで v2 の断熱壁温が NS と試験部で最大 0.281 %、(d) 表の範囲外は例外、(e) 混合気の μ が独立参照と 2.8e-10。**(b) の初版の基準「積分した δ_r が相対 1e-4」は不合格だった** (断熱 2.1e-4、300 K 5.2e-4): 式は (b1) で一致し、差は縁の状態の補間と数値微分の離散化から来る (v1 自身の縁の格子 4000 → 64000 の差が 1.0〜2.3e-3、積分の許容差 1e-6 → 1e-9 の差が 2〜3e-4)。1e-4 は積分器の誤差より小さく、試験の設計の誤りとして (b1)・(b2) に分けた。既存の design の試験 (deltastar・mesh_params・mesh_euler・wall_single_bspline・pw_upstream_poly・run_tests) も合格 (mesh_euler の (d) は新しい格子のキーを使う 2 本を「変更前と同一」の対象から外し、キーが入ることを別に確認) |
| 21 | **抽出の座標の A/B (2026-10-08 事前登録、codex diagnose、担当 O、CFD 0 step)** | 判断: 2026-10-08 codex diagnose ([記録](../../notes/reviews/2026-10-08-cold-pair-result-diagnose.md)) — 「判定不能を維持。閉包を変える前に抽出の座標の A/B を」。近壁を壁法線にした格子は列の中で x が変わる (x ≈ 40 で列内の幅 0.153 r_t) のに、抽出器 (`metrics/deltastar.py`) は列全体を `x[i, 0]` として積分していた。A = 今の抽出 (列をそのまま)、B = 各行を x 方向に補間して一定 x の断面にそろえる (`_load_structured(remap_constant_x=True)`、opt-in)。Euler 参照・帯 E の選び方・求積・符号・平滑化は同じ。対象は run_0181・run_0182 の判定窓の 5 枚 (CFD 0 step)。**判定 (事前登録)**: 5 枚・試験部 [40, 94] (0.25 刻み) の max \|R_B/R_A − 1\| が 1 % 以上なら「抽出の座標の不整合の寄与を支持」、1 % 未満なら棄却 (棄却しても摩擦の閉包が原因とは確定しない)。`cold_pair.py extract-b`・`ab-compare` → `_band_ab/cold_pair/extract_ab.json`。**結果 (2026-10-08、`_band_ab/cold_pair/extract_ab.json`)**: max \|R_B/R_A − 1\| = **0.060 %** (5 枚・試験部; 平均 0.019〜0.020 %) で**第 1 仮説を棄却**。各腕の δ_E は A/B で最大 0.64 % (断熱)・0.70 % (300 K)、平均 0.13〜0.15 % 違うが、比ではほぼ打ち消す (R は x = 40 で 0.7250 → 0.7254、x = 94 で 0.7869 → 0.7869)。冷却で δ_E が 21〜27 % 薄くなるのは抽出の座標の不整合ではない |
| 22 | **V-c45 の残りの指摘の採否 (2026-10-08 codex diagnose)** | 採用: judge のゲートの集約 (NaN 検査の記録を読む、収束判定の文言で判定不能・RISING・DIVERGED を不合格に、記録の欠落も不合格) と、U の移動平均の端を有効点数で割る (x = 40 の U 0.85 % → 1.24 %) — `cold_pair.py` を修正済み。ゲート 2 (残差の床): 「格子・精度の違う run の絶対値との比較」は不適切という問題提起は採用、私の代案 (同じ run の中で横ばい + 派生量の STEADY) は却下 (高い残差での停滞を通す) — **同じ格子・精度で、残差の正規化・保存の収支・対象量の精度と結びついた上限が要る** (未設計)。ゲート 4: 縮流部の超過は不成立のまま、緩めるには格子感度の裏付けが要る。「熱閉包の形を選ぶ」問いは置き換えない (相対の優劣と絶対の精度は別の問い; 今の `decision_by_kf` の「支持」は条件付きの計算結果で正式な結論ではない)。**設計の手順: 壁温を変えたら NS で再評価は必須、壁の作り直しは設計公差を外れたとき** (§4.6 の感度の判定と同じ)。contur_v2 は opt-in のまま。基準を改訂するときは旧 V-c45 と「判定不能」を残し、改訂日・観測済みのデータ・理由・新しい条件を別版で書き、同じデータの再解析は事後解析と明記する。延長 (run_0183) の結果だけでは旧基準を通せない (断熱側の残差の床の不成立が残る) |
| 23 | **摩擦の冷却応答 (2026-10-08、記録、CFD 0 step)** | NS の壁面せん断 (接線、`res_wall_3_100000.h5`) を CONTUR と同じ非粘性の縁の状態 ½ρ_e u_e² で割った C_f を両腕で比較 (`_band_ab/cold_pair/cf_compare.json`)。試験部の平均で **C_f の比 (300 K/断熱) は NS 1.487、CONTUR 温度形 1.513・エンタルピー形 1.466・v2 1.436 — 冷却による摩擦の増え方はどの形も ±4 % で合う**。絶対値は NS/CONTUR (生産の k_f) が断熱 0.89〜0.94・300 K 0.91〜0.92。NS の 300 K の C_f は判定窓の中で −0.4 %。CONTUR は冷却で θ が 1.9〜2.0 倍・H が 0.52〜0.53 倍になって δ* がほぼ変わらない。**ずれ (NS の δ_E の比 0.75 vs CONTUR 約 1.0) は摩擦でなく、θ の発達か δ*/θ の関係の側** (どちらかは NS の θ が要る; 原因は未確定) |
| 24 | **NS の θ の抽出 (2026-10-08、記録、担当 O)** | 一定 x の断面 (`remap_constant_x`) で、δ_E と同じ帯の外縁を縁として運動量の欠損から θ を作り、両腕の θ の比と δ_E/θ の比を CONTUR (θ 1.9〜2.0 倍、H 0.52〜0.53 倍) と比べる。δ_E/θ は CONTUR の H と定義が違う (codex diagnose) ので、比どうしだけを比べる。**結果 (2026-10-08、run_0181/0182 の res_100000、`_band_ab/cold_pair/theta_*.npz`)**: 300 K/断熱の比は x = 40 → 94 で、NS の θ_r が 2.66 → 2.40 (CONTUR は温度形 2.05 → 1.94、エンタルピー形 1.94 → 1.85、v2 1.93 → 1.83)、NS の δ_loc/θ_r が 0.25 → 0.33 (CONTUR の H の比 0.51〜0.53)、NS の δ_loc が 0.66 → 0.79 (抽出の δ_r 0.73 → 0.79 と整合)。NS の δ_loc/θ_r の値は断熱 11.6〜16.3 (CONTUR の H 12〜14.8 と同程度)、300 K 2.9〜5.4 (CONTUR の 300 K の H は約 7)。**冷却による形状係数の下がり方が NS (1/3〜1/4) と CONTUR (1/2) で違い、θ の増え方も NS が約 3 割大きい** (摩擦の比は ±4 % で合う; θ の差は運動量式の加速の項に H が入ることと整合)。冷却平板 case/48 (一定 M) では CONTUR の δ*・θ が NS と 6 % 以内だった。原因は未確定 (候補: 強い加速の下での形状の閉包 — 速度分布の N(R_δ) と温度分布の組み合わせ — と縮流部からの履歴)。300 K の run はまだ DRIFTING なので、延長 run_0183 で測り直す |
| 25 | **冷却壁で効く CONTUR の係数の試算 (2026-10-08、記録、CFD 0 step)** | ユーザ「NS の反復をしながら CONTUR の係数を調整する感じで合ってる?」への回答として測った。生産の CONTUR (温度形) の 300 K/断熱の δ_r の比を係数を 1 つずつ動かして NS の R (延長後、x = 40/60/80/94 で 0.726/0.762/0.785/0.793) と比べた: 今 (a 1・k_N 1) 1.05→1.01・ずれ最大 45 %、a 0.4 で 0.84→0.81・16 %、a 0.2 で 0.71→0.69・13 % (C_f の比 2.08)、k_N 1.6 で 1.04→1.00・43 %、k_f 0.75 で 0.78→0.76・8 % (C_f の比 1.13、NS 1.49)。**NS の比は下流へ 1 に近づくのに、CONTUR はどの係数でも平らか下流で下がる — 係数 1 つでは冷却壁の x 方向の分布を再現できない**。k_N はほぼ効かず、a は効くが摩擦も動かす。出口だけ k_f で合わせると x = 40 で 1 割近くずれ、摩擦も NS から離れる。どの形で冷却壁の壁を作るか (NS の δ_E を直接使う / CONTUR に x の傾きの補正を足す / 冷却の閉包を直す) は未決 (ユーザ判断・諮問) |
| 26 | **冷却ノズルの SU2 との照合 (2026-10-08 ユーザ指示「SU2 の計算も同時に進めて。SU2 と比較するために。dilatation 補正は本来は入れるのが正しいと思う」、事前登録、担当 O、判定は F)** | 目的: 「forge が間違っている説」の切り分け (冷却で δ_E が 21〜27 % 薄くなる、が forge の解き方のせいでないか)。SU2 は燃焼ガスの多成分を同じ形で解けないので **CPG (γ 1.27354、R 292.13、Sutherland 1.716e-5/273/111、Pr 0.72、Pr_t 0.9) でそろえる**。壁と格子は生産の TP の対 (run_0181/0182) と同じ (CPG の設計は CFD ピンの凍結源が TP なので拒否される — 2026-10-08 確認) で、**CFD だけを CPG にする** (準備済みの run の設定を書き換え: thermalMethod 0・viscMethod 1・cp 1360・γ 1.27354、化学種を外す)。初期値は収束した TP の解 (run_0181 / run_0183 の最終場) の ρ・U・P・k・ω を CPG の保存量に組み直し、段階起動 full。SU2 は同じ格子 (`nozzle.msh` を 17 桁のまま su2 形式に変換) で同じ節点を解く。**run**: forge (FP64、AWS) = CPG × {断熱, 300 K} × {素の SST (dilatationCorrection 0・katoLaunder 0、SU2 と同じ形), 生産の SST (dilatationCorrection 2・katoLaunder 1)} の 4 本、SU2 v8.5 (手元の PC、6 スレッド × 2) = CPG × {断熱, 300 K} の 2 本 (SST V2003m、ROE + MUSCL、`MARKER_ISOTHERMAL` 300 K、Sutherland 273/111)。**比べる量 (Euler の参照が要らないもの)**: SU2 の解を同じ節点で forge の形式に写し、(a) 帯の外縁 y_b(x) を全 run で共通 (TP の断熱 run_0181 の抽出の band_y_b) にした δ_loc と θ_r (`theta_diag` と同じ式、一定 x の断面)、(b) 壁の C_f と q_w、(c) x = 40・70・94 の断面の ρu・T の分布。**判定 (事前登録)**: forge 素 SST と SU2 で、冷却の比 R_loc = δ_loc(300 K)/δ_loc(断熱) の差が試験部 [40, 94] で max \|R_forge/R_SU2 − 1\| ≤ 2 %、各腕の δ_loc・θ_r の差 ≤ 3 %、C_f の差 ≤ 3 % (前回の断熱 CPG の照合の実績: δ* ≤ 1.3 %・θ ≤ 0.4 %・δ99 ≤ 3 %) なら「forge の解き方は冷却ノズルでも SU2 と一致」。超えたら forge 側の原因を調べる (諮問)。記録のみ: forge 生産 SST − 素 SST の R_loc の差 (dilatation 補正の効き)、forge CPG と TP の R_loc の差 (気体の効き)、q_w。**比べる前提**: 全 run で NaN なし・RISING/DIVERGED なし・比べる量 (各腕の δ_loc の x = 40/70/94 と Q_w) が判定窓 5 枚で STEADY (SU2 は最後の 5 つの出力)。満たさなければ延長してから比べる (同じ回数ずつ)。未収束どうしを「一致」と呼ばない。**codex diagnose 2026-10-08 ([記録](../../notes/reviews/2026-10-08-cold-su2-crosscheck-diagnose.md)) の採否 (全件採用) と改訂 (投入前)**: (1) 物性: TP の設定は `thermCondMethod` を消しているので、CPG では `thermCondMethod: 1`・`prandtlLam: 0.72` を明示し `transport`・`species` を外す (prep-cpg が実効の physProp を検査して記録)。(2) SST の形: SU2 V2003m は ω の生産を制限後の P_k から作り σ を F1 でブレンドする (SU2 v8.5 `turb_sources.hpp`・`turb_diffusion.hpp` で確認) — forge の `sstOmegaProdFromPk: 1`・`sstSigmaBlend: 1` (今の既定) に当たる。**§4.3 の対応表 (両方 0) は誤りとして訂正**し、素の SST = dilatationCorrection 0・katoLaunder 0・他は既定。入口の k/ω (forge は k 1・ω 18000 固定、SU2 は TI 0.077・粘性比 12.1)、壁の ω、エネルギー中の k (forge 0、SU2 は含む — 2026-10-08 にソースで確認) の違いは実効条件の表に記録。SU2 の出口は 2237 Pa (旧 cfg の 1389 Pa は流用しない)。(3) 結論の範囲は「この CPG 条件・格子・比較量で forge と SU2 の差が許容内か」まで (TP の h(T)・組成・混合気の物性と冷却の相互作用は検証の外)。不一致でも forge の誤りと即断しない。(4) 帯の外縁: 300 K の y_b は断熱の 1.21〜1.24 倍なので、**y_b = 各 x で TP の断熱 (run_0181) と 300 K (run_0183) の y_b の大きい方**、感度として ×1.25 も記録。δ_loc・θ_r は固定した観測量であって δ_E の検証の代わりではない。(5) 許容差は暫定の照合目標: R_loc ≤ 2 %、各腕の δ_loc・θ_r ≤ 3 %、C_f ≤ 3 %、**q_w ≤ 5 % (§4.4 の目標) を判定に入れる**。(6) 比べる前提に**保存の収支**を加える: 入口と出口の質量流量の差 ≤ 0.1 %、全エンタルピーの収支 (入口 − 出口 − 壁の入熱) ≤ 0.1 % of 入口 (両コード)。比べる量 (δ_loc・θ_r・C_f・断熱の壁温・Q_w) が判定窓 5 枚で STEADY。延長は各コードのゲートで決め、反復数はそろえない。(7) 生産の SST の対は dilatation と Kato-Launder の 2 因子を同時に変えるので「生産の SST 設定一式の感度」と呼ぶ **(8) 量の測り方 (2026-10-08、どちらのコードの結果も見る前に固定)**: 後処理は `case/45.isobutane_m6_d155/cold_xcheck.py` (forge は res、SU2 は restart の CSV を同じ節点の (ni, nj) 構造に並べ、同じ関数で処理)。δ_loc・θ_r は `theta_diag` と同じ式 (run_0181 res_100000 で TP の θ の npz と全桁一致を確認)、y_b は (4) のとおり。C_f = τ_w / (½ρ_e u_e²)、縁の値は y_b での断面の値を壁の x へ補間。τ_w = μ(T_w)\|∂u_t/∂n\|、q_w = μ(T_w) c_p/Pr ∂T/∂n (壁から流体へ向かう法線、流体から壁へ入る向きを正)、微分は壁節点と壁法線上の第 1・第 2 内部節点の 2 次の片側差分、μ は Sutherland 1.716e-5/273/111。Q_w = ∮ q_w 2πr ds。C_f・q_w・δ_loc・θ_r の差は試験部 [40, 94] の 0.25 刻みの点ごとの \|forge/SU2 − 1\| の最大で判定する (q_w は 300 K の腕だけ)。収支は入口と出口の列 (出口の列は壁法線の層で一定 x でない) を通る流束 2π∫(F_x dr − F_r dx) r。全エンタルピーは h0 = c_p T + ½\|u\|² に、**SU2 だけ k を足す** (SU2 の全エネルギーは k を含む: `CNSVariable::SetPrimVar` で静的エネルギー = E − ½u² − k を確認。forge は sstEnergyIncludesK 0)。準定常は判定窓の各スナップショットで δ_loc・θ_r・C_f の x = 40/70/94、断熱の T_w の x = 40/70/94、300 K の Q_w を `check_quasisteady.classify` (5 枚、drift・osc 0.1 %) にかける。**SU2 の初期値**は forge の run_0184/0185 の段 S1 直後の場 (README 参照)。SU2 の C_f・熱流束の表面出力は compact restart のため無いので、SU2 の履歴の `HF` と自前の Q_w の差を記録する (後処理の式の確認、判定には使わない) |

## 6. 検証

- **単体 / 回帰**: `wall_thermal` 省略で `_bcond` 出力が現行と文字列一致 (Euler: case/45 `run_0037` / case/44 `run_0091`、NS: case/44 `run_0107` の `bcondConfig.yaml` と diff 0、積分法初期壁も一致)。`deltastar_integral` の `prescribed_temperature` 単体 (既存 V3) と新設の平板入口 (`flat_plate_integral`) の単体。符号付き δ\* の正/零/負回帰。
- **平板 (S2)** — 合否ライン (§4.3):
  - $C_f$: VD-II 絶対 ±10 %、比 B/A ±5 %。$2St/C_f$ ∈ [0.95, 1.25]。回復係数 (run A) 0.88–0.90。
  - δ\*/θ: CONTUR (平板入口) 絶対 ±15 %、比 B/A ±10 %。温度–速度関係 (Walz / Duan–Martín) は診断 (どちらに近いかを記録)。
  - $q_w$ 閉合: $\int q_w\,ds$ と入口・出口の全エンタルピー流束差が 5 % 以内 (抽出の検証)。
  - SU2 素 SST 同一メッシュ (A-plain↔SU2-A, B-plain↔SU2-B): $C_f, q_w, \delta^*, \theta$ ≤3 %。流れ方向格子感度 (nx 1000→1500) の差が同じ 3 % 未満であること。
  - y₁⁺ 掃引: δ\*・q_w の y₁⁺ 依存が ≤2 % になる上限 y₁⁺ を決める (期待 1〜2)。
  - 全 run: `check_convergence.py` / `check_quasisteady.py --quantity theta,cf_retheta` / `cooled_plate_eval.py --series` / `check_mesh_quality.py` の VERDICT を README に貼る (§4.3-6 の定量条件)。
- **ノズル CPG (S3)**: δ\* ≤3 %・θ ≤1 %・$q_w$ ≤5 % (4 ステーション)。壁温 = 300 K がノード値で再現 (ピン)。
- **生産 TP (S4)**: 生産ゲート (ṁ 比 ≤0.3 %、出口面コア M ±0.1 %) 達成、`check_quasisteady --series` STEADY、壁温感度台帳。
- **V-c45 冷却壁の NS の対 (2026-10-08 事前登録、§5.1 #14・#15〜#19 で改訂、投入前)**: case/45 の生産の壁 (`_band_ab/prod_confirm/prep/wall_repr.json`、run_0167 の入力をビット同一で再現した準備) を
  冷却壁用の格子 (`cold_pair_mesh.py`、ni 4719 × nj 121、近壁は壁法線 `wall_normal_layer [0.02, 0.8]`、msh 17 桁) に載せ、FP64 のビルドで断熱 (`run_0181_ns_coldmesh_ad`) と 300 K (`run_0182_ns_coldmesh_tw300`) を回す。
  **この判定は「この格子・FP64」の条件での比較**であり、生産の FP32 の格子に対する予測精度の主張ではない (codex diagnose 2026-10-08)。
  - **バイナリ**: ソース = 生産のバイナリの commit e2696d8f0 + `flowFormat.hpp` の typedef 4 行 (double) + `mesh/gmshReader.hpp` の座標読み込み (double のビルドは `stod`、c096d66c と同じ差分)。
    Release・CUDA arch 86・CXX flags は生産と同じ。forge sha256 65be5e28ca1aed9f…、変換器 sha256 ac88861fa04ebdf4…。両 run で同じものを使い、run の記録に全差分と sha256 を残す。
  - **投入の前提 (不成立なら回さない)**: 品質は厳密な `VERDICT: PASS` (`--ar-max 5000`、SOFT-PASS 不可)。スキュー > 0.1 かつ AR > 1000 のセルが 0。
    変換後の座標の第一層厚が生成時の倍精度座標と相対 1e-6 以内 (`fp64_reader_ab.py` で B 腕は誤差 0 を確認済み)。物理壁の差 |Δr| ≤ 1e-8 r_t。壁距離が変換し直しと相対 1e-6 以内。
    壁の bcond が断熱 `wall` / 300 K `wall_isothermal` Ts 300。熱境界条件以外 (格子・物性 DB・実効設定) が 2 本で同じ。
  - **ゲート (結果を使う前提、不成立なら判定不能として延長・見直し)**:
    1. NaN・Inf なし (全段の残差と全 res)。
    2. `check_convergence.py --segment` で RISING・DIVERGED なし。plateau は許すが VERDICT の文言をそのまま残し「収束」と書かない。
       残差の床の上限: 本段の最後の 5000 step の rms_ro・rms_roUx・rms_roUy・rms_roe・rms_roK・rms_roOmega の中央値が、生産の run_0179 の同じ量以下。
    3. 準定常: 判定窓 80000〜100000 の 5 枚すべてで、x = 40, 41, …, 94 (1 r_t ごと) の両腕の δ_E、R_NS、断熱の壁温、300 K の壁の熱流束の積分 Q_w を
       `check_quasisteady.classify` (5 枚全部、drift ≤ 0.1 %・osc ≤ 0.1 %) で STEADY。
    4. 壁解像 (`check_wall_resolution.py --weight area`、300 K の run、断熱も記録): y1+ > 1 の面積割合が全壁で ≤ 5 %、試験部 [40, 94] で ≤ 1 %、x ∈ [−1, 40) で ≤ 5 %、
       縮流部 x < −1 で ≤ 10 % (入口の角 0.05 r_t は別に記録)。超過位置を記録する。
  - **判定の量**: R_NS(x) = (5 枚平均の δ_E,300 K) / (5 枚平均の δ_E,断熱)、δ_E は生産の抽出の経路 (`extract_and_merge`、帯 E、300 K は符号付き) の平滑後の値、x は [40, 94] の 0.25 r_t 刻み。
    CONTUR の予測 R_A (温度形)・R_B (エンタルピー形) は熱閉包だけを替えた `delta_contur_compare.py hform` の経路で、同じ壁・同じ k_f を両腕に使う。k_f は診断用の固定値として 1 と 1.0541 の両方。
  - **不確かさ U(x) (相対、線形和)**: u_t = 5 枚それぞれの比 R_k/R_NS − 1 の最大の絶対値。u_ext = 抽出の感度 (`delta_r_sens`、帯の係数 1.125 / 2.25) による δ_E の相対の振れの最大を、
    2 本で足したもの (2 r_t の移動平均)。u_c = 0.001 (CONTUR の積分・表、A/B の前提で ≤ 0.1 %)。U = u_t + u_ext + u_c。格子・精度を替えた感度
    (この格子の断熱の δ_E と生産の FP32 の δ_E の比) は記録するが U には入れない (判定の範囲を「この格子・FP64」に限るため)。
  - **判定**: 各 x で R_NS が [R_NS(1 − U), R_NS(1 + U)] を動くときの |R_m/R − 1| の下限・上限を取り、e_m の下限 = 下限の最大、上限 = 上限の最大 (m = A, B)。
    **e_B の上限 < e_A の下限ならエンタルピー形を支持**、e_A の上限 < e_B の下限なら温度形を支持、それ以外は判定保留。k_f の 2 値で判定が違えば保留。
    勝った側の e の上限が 1 % 以下なら「この格子で冷却の効果を 1 % 以内で当てる」、超えるなら「冷却の効果は当てきれない (壁温を変えたら NS の δ* 反復が要る)」と記録する。
  - **記録のみ**: θ・H・C_f・q_w の 300 K/断熱の比 (NS と CONTUR A/B)、NS の断熱壁温と CONTUR の T_aw、断熱の run と生産 (run_0179、FP32・別格子) の δ_E の差、
    contur_v2 全体 (熱閉包 + 粘性 + 加速の項) の予測 (§5.1 #19、判定には使わない)。
  - **延長の決め方 (2026-10-08、結果の 1 回目の後・延長の前に追記)**: ゲート 3 が不成立の腕だけを、その run の最後の res から `cold_pair.py prep-ext`
    (restart_field、同一格子・ビット一致・FP64 の型のまま、nStepOuter だけ変える、段なし) で 100000 step 延長する。判定窓は**延長した run の最後の 5 枚**。
    STEADY だった腕 (断熱) は元の判定窓をそのまま使う。延長してもゲート 3 が不成立なら、その時点の値と等比外挿の漸近値を記録して判定は保留し、
    さらに延長するかは諮問で決める。ゲート 2 (残差の床) とゲート 4 (縮流部の壁解像) の不成立の扱いは諮問で決める (結果を見た後に基準を緩めない)。
  - **1 回目の結果 (2026-10-08、`_band_ab/cold_pair/V_c45.json`)**: **判定不能 (ゲート不成立)**。
    ゲート 1 (NaN) 合格。ゲート 2: 両 run とも `NOT CONVERGED (stalled/plateau)`、RISING なし、**残差の床が生産の run_0179 を上回る** (rms_ro 7.0e-6・9.1e-6 vs 3.4e-7 ほか全量、10〜40 倍)。
    ゲート 3: 断熱の腕は δ_E・壁温とも全点 STEADY、**300 K の腕は δ_E が全点 DRIFTING** (5000 step ごとの増分が 0.65〜0.85 倍ずつ減る単調増加)、
    **Q_w も DRIFTING** (−12.34 → −11.49 MW、増分の比 0.91〜0.92)。ゲート 4: 300 K の全壁 1.8 %・試験部 0 %・[−1, 40) 0 % は合格、**縮流部 (入口の角を除く) 24 % で上限 10 % を超える** (最大 y1+ 1.23)。
    記録 (判定には使わない): R_NS = 0.725 (x = 40) → 0.787 (x = 94)、**冷却で δ_E が 21〜27 % 薄くなる**。CONTUR の予測は温度形 1.05 → 1.01、エンタルピー形 1.02 → 0.99、
    contur_v2 1.01 → 0.98 で、どの形でも冷却の効果をほぼ再現しない (e は 0.40〜0.46)。300 K の δ_E の等比外挿の残りは x = 94 で約 +0.17 mm (R にして +0.004)。
    格子と精度を替えた感度 (この格子の断熱の δ_E / 生産の FP32 の δ_E) は試験部の平均 +0.32 %・最大 0.52 %。原因は未確定 (諮問の前に書かない)。
    **codex diagnose (2026-10-08、[記録](../../notes/reviews/2026-10-08-cold-pair-result-diagnose.md))**: 判定不能を維持。「冷却で 21〜27 % 薄くなる」「CONTUR が冷却の効果を外す」は保存した窓の抽出値とモデルの差であって、定常解に対する確定した予測誤差ではない。次は抽出の座標の A/B (§5.1 #21)。採否は §5.1 #22。
  - **2 回目 (延長、2026-10-08、300 K = `run_0183_ns_coldmesh_tw300_ext` の最後の 5 枚、断熱 = run_0181 の元の窓、`_band_ab/cold_pair/V_c45.json`; 1 回目は `V_c45_v1.json`)**: **判定不能 (ゲート不成立)**。
    ゲート 1 合格 (延長も 21 枚 CLEAN)。ゲート 3: STEADY でないのは x = 40 の δ_E・R (TRANSIENT-UNSETTLED) と Q_w (−9.79 → −9.55 MW、DRIFTING、増分の比 0.956、等比外挿の残り約 1.3 MW) の 3 つ (1 回目は 111)。
    ゲート 2: 延長の run で rms_roOmega が RISING (`check_convergence --segment`: init 1.17e3・peak 1.92e3・fin 9.1e2)、残差の床は不成立のまま。ゲート 4: 縮流部 23 % (上限 10 %) で不成立のまま。
    記録: R_NS = 0.726 (x = 40)・0.762 (60)・0.785 (80)・0.793 (94)、1 回目からの変化は +0.1〜+0.8 %。e は温度形 0.44〜0.46・エンタルピー形 0.40〜0.41 (k_f 1 と 1.0541)。U 0.6〜1.4 %。
    延長後の θ の比 (`theta_run_0183_*.npz`) は x = 40 → 94 で θ_r 2.85 → 2.53、δ_loc/θ_r 0.23 → 0.32、δ_loc 0.66 → 0.80 (1 回目より θ の差がさらに大きい)。**延長しても「冷却で δ_E が 21〜27 % 薄くなり、CONTUR はどの形でもそれを再現しない」は変わらない**。
    さらに延長するか、ゲート 2・4 をどう作り直すかは未決 (諮問・ユーザ判断)。

### 6.1 レビュー記録 (codex)

[`AGENTS.md`](../../AGENTS.md) 「codex レビュー」の記録表 (`solver_density_cuda/tools/codex_review.py`)。

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-12 | [2026-09-12-tooling-nozzle-isothermal-wall-chain-plan.md](../../notes/reviews/2026-09-12-tooling-nozzle-isothermal-wall-chain-plan.md) | GO-with-changes, C0/M7/m2 | **全件採用**: M1 符号付き δ\* → §4.1 + §5.1-1a / M2 低 Re $q_w$ 後処理 + $ds$ 重み + 閉合 → §4.1 + §5.1-1b / M3 定量ゲート (`--quantity cf`・`--series` の誤記も訂正) → §4.3-6 + §5.1-1c / M4 固定形状感度と再設計の分離 → §4.5 (S4a/S4b) / M5 温度–速度関係は診断 + Duan–Martín + 平板積分入口 → §4.3-3 / M6 引き継ぎ 3 種 + 起動ゲート → §4.4 / M7 A-plain 追加 + SST 条件表 + nx 感度 → §4.3 / m1 `--init-thermal` 廃止 → §4.1 / m2 y⁺ ×5.5〜6 → §4.2。スポット検証 (codex): `run_0048` の `qwall` 全点 0、`run_0013`/`run_0107`/`run_0048` は既定条件で NOT CONVERGED (stalled) — 本計画の warm 床条件はこれを踏まえて定量化 |
| plan | 2026-10-08 | [2026-10-08-tooling-nozzle-isothermal-wall-chain-plan.md](../../notes/reviews/2026-10-08-tooling-nozzle-isothermal-wall-chain-plan.md) (§4.7・§5.1 #13・§6 V-c45 に集中) | GO-with-changes, C0/M6/m2 | **全件採用** (§5.1 #15〜#20): M1 格子の AR (傾斜壁の接線長・半径方向の線はスキュー有りで AR 5000 の例外外) → #15 / M2 float32 の数 ulp・法線距離と半径方向間隔の区別 → #15・#16 (FP64 ビルドの諮問中) / M3 壁解像は点数でなく面積で・領域別 → #17 / M4 比だけでなく両腕の δ_E・断熱壁温・Q_w の定常・残差床の上限 → #18 / M5 e と u の尺度・不確かさの寄与・区間の分離で判定 → #18 / M6 熱閉包だけの A/B と contur_v2 全体を分ける・hform の k_f が常に 1 → #19 / m7 CPG の一致と v1 再現・Sutherland 定数・M=1 の 0/0 → #20 / m8 B 腕の精度確認・表の範囲外 → #20。V-c45 の判定基準は #18 を反映して投入前に書き直す |

## 7. 影響範囲

- `design/forge_design/{probdef.py, evaluate/runner_wt.py, evaluate/runner_sern.py, evaluate/runner_axismach.py, feedback/deltastar_loop.py}`
- 新規 `case/48.flat_plate_cooled_m4/`、case/45・case/44 の問題 YAML 変種と run
- docs: `methods/design/overview.md` (新節)、`procedures/recommended-settings.md` (壁行)、`procedures/verification/README.md` + 新ファイル、`design/CAPABILITIES.md`
- forge 本体は変更しない (等温壁・低 Re SST は既存)。変更が要ると分かった場合 (例: 冷却壁での EOS 床・ピンの不具合) は本 plan に §4.7 として追記してから触る

## 8. 未確定事項 (ユーザ確認)

1. **公称 $T_w$**: 初版は 300 K (ユーザ発言 2026-09-12「壁面温度 300 K くらいでまずは」)。実機の冷却方式が決まったら $T_w^{\rm nom}$ を差し替える。
2. **CHT の要否**: §4.6-4 の判断基準 (弱 CHT ループが収束するか) で決める。フル CHT の実装は本計画の外。
3. ~~**AR ゲートと y₁⁺**: 冷却壁で y₁⁺ ≤ 上限を守ると AR が 1000 を超える形状があり得る。その場合 `ni` を増やす (計算時間 ×数倍) か、AR 上限を「壁法線方向の構造格子は AR 2000 まで可」に緩めるかは S2/S3 の結果でユーザ判断。~~ **決着 (2026-09-12, ユーザ決定「AR 上限緩めようか」)**: 壁法線の構造格子層に限り **AR ≤ 5000** まで緩和 (AGENTS.md / calculation-workflow / recommended-settings に反映、`check_mesh_quality.py --ar-max`、問題 YAML `mesh.ar_max`)。裏付け = §5.1-9 の A/B (AR 846 メッシュ run_0108/0114 vs AR 4140 メッシュ run_0116/0117): **確定 — 出口 M・ṁ の差 0.02 %、熱負荷は y₁⁺ 改善分 +8〜10 %、発散なし。冷却壁ノズルの生産メッシュは `problem_va_R2_LU6_Lc8_ns_ar5k*.yaml` (AR ≤ 5000, スロート y₁⁺ ≈ 1) とする。**
4. **壁温は頻繁に変わる前提で設計チェーンを作る (2026-10-08 ユーザ発言)**: 「これから壁温設定はガンガン変わり得る」「等温壁にするかもしれないし、壁温分布を与えるかも」「CONTUR で物性をどの温度で評価するかは任意性がある話なので、どうするべきかは考えて」。対応は §5.1 #10〜#12。

## 9. 完了条件

- [ ] `methods/design/overview.md` の新節を実装と同期
- [ ] S1–S5 完了、§6 の判定を満たす
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映
- [ ] `status: done`、§10 変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 10. 変更ログ

- `2026-10-08` — **壁温が頻繁に変わる前提 (§8-4、ユーザ発言) を受けて §5.1 #10〜#12 を起票**。case/45 の実測で、CONTUR の断熱壁温の式が燃焼ガスで全温を超えること (+190 K)・その式の違いで δ_r が 3〜5.5 % 動くこと・μ が NS と最大 9 % ずれること・出口合わせの k_f が壁温条件で変わることを確認した (数値は §5.1 #10、詳細は [verification-m6-axis-wave-mesh-su2](verification-m6-axis-wave-mesh-su2.md) §9 の 2026-10-08)。物性の評価の方針は codex (diagnose) に諮問中で、§4 は未反映。
- `2026-10-08` — **codex (diagnose) に諮った**: [`notes/reviews/2026-10-08-contur-property-temperature-diagnose.md`](../../notes/reviews/2026-10-08-contur-property-temperature-diagnose.md) — 結論「係数と輸送物性を固定し、熱閉包だけを温度形からエンタルピー形へ替える CFD 0 step の A/B を先に」。指摘 10 件は全件採用 (§5.1 #10〜#12 に採否)。運動量式の項の近似 (#10b) と θ の単位 (#10f) は手元で再現してから採用した。第 2 層 (圧縮性変換の参照温度) の選定は保留: case/48 の断熱基準の壁温ドリフトと case/44 run_0114 の熱負荷の変動のため、今のデータでは決めない。後で足す NS の最小候補は case/45 の生産形状の 300 K (メッシュが冷却壁の壁解像を満たせば 1 本、満たさなければ同じメッシュで断熱と 300 K の 2 本)。#10a を事前登録した。
- `2026-10-08` — **#10a の A/B を実施: 第 1 仮説を支持** (300 K/断熱の δ_r 比の変化 最大 2.80 %、基準 1 %)。熱閉包の書き方だけで冷却壁の δ_r が数 % 動く。詳細は §5.1 #10a。説明ページ https://claude.ai/artifact/5dFnxMbYCD4DzbmXLk6QX7 (CONTUR の式と較正の係数)。
- `2026-09-12` — **AR 上限を壁法線構造層で ≤5000 に緩和 (ユーザ決定)**。AGENTS.md・手順書・runner (`mesh.ar_max`) に反映。A/B (§5.1-9) で平均流 0.02 %・熱負荷 +8〜10 % (y₁⁺ 改善分)・発散なしを確認し確定。冷却壁の生産メッシュを AR 5000 版 (スロート y₁⁺ ≈ 1) に切替。
- `2026-09-12` — **S4 完了 (case/44)**: S4a 固定形状の壁温感度 (run_0108 断熱 vs run_0109/0114 300 K, 同壁同メッシュ): 出口コア M +0.35〜0.48 %・ṁ +0.35 %・出口質量平均 T0 −1.5 %・Q_w 6.4〜7.8 MW (q_w ピーク 2.3〜2.8 MW/m² @ スロート直前)。S4b 300 K 再設計 (run_0110 pass 0 → run_0115 pass 1): ṁ 比 1.0002・出口コア M +0.02 % で生産ゲート達成、積分法初期壁は冷却効果を過小評価するが pass 1 で固定点近傍。S3 は run_0112 (CPG 断熱) 完了・SU2 対は実行中。
- `2026-09-12` — S1 配管・符号付き δ\*・平板評価ツール実装、case/48 run A/B/C 完了 (§5.1-2 に数値)。y₁⁺ の冷却倍率は実測 ×5.7 (3 µm: 断熱 0.08 → 300 K 0.46)。
- `2026-09-12` — codex plan 段レビュー (GO-with-changes, C0/M7/m2) を §6.1 に記録、**全件採用**して §4.1/§4.2/§4.3/§4.4/§4.5/§5.1/§6 を改訂 (符号付き δ\*、低 Re $q_w$ 後処理と閉合、定量ゲート、S4a/S4b 分離、SST 条件表、引き継ぎ 3 種)。case/48 の起動レシピ実測 (層流暖機 + 2 次 cfl ランプ、本段 cfl 2) を §5.1-2 に記録。
- `2026-09-12` — 起票 (ユーザ依頼: ノズル設計の壁を断熱から等温/CHT へ。等温 low-Re SST を平板で理論・経験式と SU2 に対して検証 → ノズルで SU2 と比較 → 等温 NS で δ\* 反復 → 壁温影響の評価方針)。ユーザ判断 2 件を反映: 壁関数は使わない (準備状況が悪い)、冷却で y⁺ が変わる点を平板の y₁⁺ 掃引で扱う。
```

## 参考: `case/45.isobutane_m6_d155/cold_xcheck.py`

```
"""冷却ノズルの forge と SU2 の照合 (plan tooling-nozzle-isothermal-wall-chain §5.1 #26) の後処理。両コードの場から同じ式で同じ量を作る。
SU2 は forge と同じ節点 (格子の節点の順番 = msh の順番 = forge の node 番号) を解くので、場を同じ (ni, nj) の構造に並べて扱う。

量 (各スナップショット):
  - δ_loc・θ_r: 一定 x の断面 (各行 j を x 方向に補間) で、帯の外縁 y_b(x) の値 (ρ_e, u_e) を縁として面積で等価な厚さ
    (cold_pair.theta_diag と同じ式)。y_b = 各 x で TP の断熱 (run_0181) と 300 K (run_0183) の抽出の band_y_b の大きい方、感度に ×1.25。
  - 壁の τ_w・q_w: 壁の法線に沿った壁節点・第 1・第 2 内部節点から 2 次の片側差分で ∂u_t/∂n・∂T/∂n、μ は Sutherland
    (1.716e-5/273/111)、λ = μ c_p/Pr (Pr 0.72)。q_w は流体から壁へ向かう向きを正。Q_w = Σ q_w 2π r ds (台形)。
  - 入口・出口の質量流量と全エンタルピー流量 (h0 = c_p T + ½|u|²、k は含めない): 列 i = 0 と i = ni−1 で 2π∫ρu_x (·) r dr。
  - x = 40・70・94 の断面の ρu・T の分布。
CPG の定数は γ 1.27354、c_p 1360、R = c_p(γ−1)/γ。TP の run を読むと h0 は CPG の式になるので、収支は CPG の run だけで見る。

usage:
  python3 cold_xcheck.py reduce-forge <run> [--last 5]      (AWS: 最後の 5 枚の res → _band_ab/cold_pair/xcheck_<run>.npz)
  python3 cold_xcheck.py reduce-su2 <run> [--last 5]        (手元: 最後の 5 つの restart_flow_<iter>.csv → 同上)
  python3 cold_xcheck.py gates-forge <run>                   (AWS: NaN の全件検査 + check_convergence --segment → xcheck_gates_<run>.json)
  python3 cold_xcheck.py gates-su2 <run>                     (手元: history.csv を forge の残差 CSV の列名に変換して check_convergence、
                                                              履歴と restart の CSV の非有限値 → xcheck_gates_<run>.json)
  python3 cold_xcheck.py compare                             (手元: 判定 → _band_ab/cold_pair/xcheck.json)
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUTD = HERE / "_band_ab" / "cold_pair"
GAM, CP, PR = 1.27354, 1360.0, 0.72
RGAS = CP * (GAM - 1.0) / GAM
XS = (40.0, 70.0, 94.0)
XE = np.arange(40.0, 94.0 + 1e-9, 0.25)


def mu_suth(T):
    return 1.716e-5 * (T / 273.0) ** 1.5 * (273.0 + 111.0) / (T + 111.0)


def mesh_info(n_nodes):
    z = np.load(OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz")
    ni = len(z["x"]); S = float(z["scale"])
    if n_nodes % ni:
        raise SystemExit(f"節点数 {n_nodes} が ni {ni} で割り切れない (格子が違う)")
    return ni, n_nodes // ni, S


def common_yb():
    a = np.load(OUTD / "theta_run_0181_ns_coldmesh_ad_100000.npz"); b = np.load(OUTD / "theta_run_0183_ns_coldmesh_tw300_ext_100000.npz")
    x = np.asarray(a["x"]); yb = np.maximum(np.asarray(a["band_y_b"]), np.interp(x, b["x"], b["band_y_b"]))
    return x, yb


def reduce_fields(xy, ro, ux, uy, T, k, h0_with_k, ni, nj, S, yb_x, yb):
    x = xy[:, 0].reshape(ni, nj) / S; r = xy[:, 1].reshape(ni, nj) / S
    RO = ro.reshape(ni, nj); UX = ux.reshape(ni, nj); UY = uy.reshape(ni, nj); TT = T.reshape(ni, nj); KK = k.reshape(ni, nj)
    xt = x[:, 0].copy()
    R = np.empty_like(r); Q = {"ro": np.empty_like(RO), "ux": np.empty_like(UX), "T": np.empty_like(TT)}
    for j in range(nj):
        R[:, j] = np.interp(xt, x[:, j], r[:, j])
        Q["ro"][:, j] = np.interp(xt, x[:, j], RO[:, j]); Q["ux"][:, j] = np.interp(xt, x[:, j], UX[:, j]); Q["T"][:, j] = np.interp(xt, x[:, j], TT[:, j])
    out = {"x": xt}
    yb_i = np.interp(xt, yb_x, yb)
    for fac, tag in ((1.0, ""), (1.25, "_s125")):
        th = np.full(ni, np.nan); dl = np.full(ni, np.nan); qe = np.full(ni, np.nan)
        for i in range(ni):
            rr = R[i]; rw = rr[-1]; rb = rw - fac * yb_i[i]
            if not np.isfinite(rb) or rb <= rr[0]:
                continue
            rf = np.linspace(rb, rw, 4001)
            rho = np.interp(rf, rr, Q["ro"][i]); u = np.interp(rf, rr, Q["ux"][i])
            re_, ue_ = rho[0], u[0]
            qm = np.trapezoid((re_ * ue_ - rho * u) * rf, rf); qq = np.trapezoid(rho * u * (ue_ - u) * rf, rf)
            dl[i] = rw - np.sqrt(max(rw ** 2 - 2.0 * qm / (re_ * ue_), 0.0))
            th[i] = rw - np.sqrt(max(rw ** 2 - 2.0 * qq / (re_ * ue_ ** 2), 0.0))
            qe[i] = 0.5 * re_ * ue_ ** 2
        out["delta_loc" + tag] = dl; out["theta_r" + tag] = th; out["qdyn_e" + tag] = qe
    # 壁の τ_w・q_w (実際の節点、壁法線の線)
    P0 = np.stack([x[:, -1], r[:, -1]], 1) * S; P1 = np.stack([x[:, -2], r[:, -2]], 1) * S; P2 = np.stack([x[:, -3], r[:, -3]], 1) * S
    d1 = np.linalg.norm(P1 - P0, axis=1); d2 = np.linalg.norm(P2 - P0, axis=1)
    tv = np.gradient(P0, axis=0); tv /= np.linalg.norm(tv, axis=1)[:, None]
    ut = lambda k: UX[:, k] * tv[:, 0] + UY[:, k] * tv[:, 1]  # noqa: E731
    c0 = -(d1 + d2) / (d1 * d2); c1 = d2 / (d1 * (d2 - d1)); c2 = -d1 / (d2 * (d2 - d1))
    dudn = c0 * ut(-1) + c1 * ut(-2) + c2 * ut(-3)
    dTdn = c0 * TT[:, -1] + c1 * TT[:, -2] + c2 * TT[:, -3]
    Tw = TT[:, -1]; muw = mu_suth(Tw)
    out["tau_w"] = muw * np.abs(dudn); out["q_w"] = muw * CP / PR * dTdn; out["T_w"] = Tw; out["x_w"] = x[:, -1]
    ok = np.isfinite(out["qdyn_e"])                       # C_f は帯の外縁 (y_b) の ρ_e u_e²/2 で割る (壁の x へ補間)
    out["c_f"] = out["tau_w"] / np.interp(x[:, -1], xt[ok], out["qdyn_e"][ok], left=np.nan, right=np.nan)
    ds = np.linalg.norm(np.diff(P0, axis=0), axis=1); rwm = P0[:, 1]
    fq = 2.0 * np.pi * out["q_w"] * rwm
    out["Q_w"] = np.array(float(np.sum(0.5 * (fq[1:] + fq[:-1]) * ds)))
    # 入口・出口の流量: 列の折れ線を通る流束 2π∫(F_x dr − F_r dx) r (出口の列は壁法線の層で一定 x でない)
    for i, tag in ((0, "in"), (ni - 1, "out")):
        xx = x[i] * S; rr = r[i] * S
        h0 = CP * TT[i] + 0.5 * (UX[i] ** 2 + UY[i] ** 2) + (KK[i] if h0_with_k else 0.0)
        fl = lambda q: 2 * np.pi * float(np.sum(0.5 * ((RO[i] * q * UX[i] * rr)[1:] + (RO[i] * q * UX[i] * rr)[:-1]) * np.diff(rr)  # noqa: E731
                                                - 0.5 * ((RO[i] * q * UY[i] * rr)[1:] + (RO[i] * q * UY[i] * rr)[:-1]) * np.diff(xx)))
        out["mdot_" + tag] = np.array(fl(1.0)); out["Hdot_" + tag] = np.array(fl(h0)); out["Kdot_" + tag] = np.array(fl(KK[i]))
    for xs in XS:
        i = int(np.argmin(np.abs(xt - xs)))
        out[f"prof_r_{int(xs)}"] = R[i]; out[f"prof_rhou_{int(xs)}"] = Q["ro"][i] * Q["ux"][i]; out[f"prof_T_{int(xs)}"] = Q["T"][i]
    return out


def load_forge(run, step):
    import h5py
    with h5py.File(run / "nozzle.h5", "r") as h:
        xy = np.array(h["MESH/COORD"], dtype=float).reshape(-1, 3)[:, :2]
    with h5py.File(run / f"res_{step}.h5", "r") as h:
        return xy, *(np.array(h["VALUE/" + k], dtype=float) for k in ("ro", "Ux", "Uy", "T", "k"))


def load_su2(run, it):
    A = np.loadtxt(run / f"restart_flow_{it:06d}.csv", delimiter=",", skiprows=1)
    xy = A[:, 1:3]; ro = A[:, 3]; ux = A[:, 4] / ro; uy = A[:, 5] / ro
    e = A[:, 6] / ro - 0.5 * (ux ** 2 + uy ** 2) - A[:, 7]            # SU2 の Energy は k を含む
    return xy, ro, ux, uy, e * (GAM - 1.0) / RGAS, A[:, 7]


SU2_COLS = {"rms[Rho]": "rms_ro", "rms[RhoU]": "rms_roUx", "rms[RhoV]": "rms_roUy", "rms[RhoE]": "rms_roe", "rms[k]": "rms_roK", "rms[w]": "rms_roOmega"}


def _conv(target: str, extra=()) -> dict:
    import subprocess
    tools = HERE.parents[1] / "solver_density_cuda/tools"
    cc = subprocess.run([sys.executable, str(tools / "check_convergence.py"), target, *extra], capture_output=True, text=True)
    lines = cc.stdout.splitlines()
    verdict = [l for l in lines if "->" in l]
    return {"cmd": f"check_convergence.py {target} {' '.join(extra)}".strip(), "rc": cc.returncode, "verdict": verdict[-1] if verdict else "(判定行なし)",
            "rising": [l.strip() for l in lines if "RISING" in l], "diverged": any("DIVERGED" in l for l in lines), "tail": lines[-12:]}


def gates_forge(run: Path) -> Path:
    """AWS: 全段の残差と全スナップショットの NaN 検査 (ns_n012.nan_scan) と、本段区間の check_convergence。"""
    import ns_n012 as NS
    nan = NS.nan_scan(run)
    NS.jdump(run / "NAN_SCAN.json", nan)
    out = {"run": run.name, "kind": "forge", "nan_verdict": nan.get("VERDICT"), "nan_first": nan.get("first_nonfinite"),
           "convergence": _conv(str(run), ("--segment",))}
    p = OUTD / f"xcheck_gates_{run.name}.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return p


def gates_su2(run: Path) -> Path:
    """手元: SU2 の history.csv (log10 の rms) を forge の残差 CSV の列名・線形値に直して check_convergence にかける (2D なので
    rms_roUz は 0)。履歴の非有限値と、全 restart の CSV の非有限値も数える。"""
    import csv
    with open(run / "history.csv") as fh:
        rd = csv.reader(fh); head = [h.strip().strip('"') for h in next(rd)]; rows = [r for r in rd if r]
    ic = {h: i for i, h in enumerate(head)}
    it = [int(float(r[ic["Inner_Iter"]])) for r in rows]
    cols = {dst: np.array([float(r[ic[src]]) for r in rows]) for src, dst in SU2_COLS.items()}
    nonfin_hist = int(sum(np.count_nonzero(~np.isfinite(v)) for v in cols.values()))
    conv_csv = run / "residual_history_su2.csv"
    with open(conv_csv, "w") as fh:
        fh.write("step," + ",".join(["rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega"]) + "\n")
        for n in range(len(it)):
            v = [10.0 ** cols[c][n] for c in ("rms_ro", "rms_roUx", "rms_roUy")] + [0.0] + [10.0 ** cols[c][n] for c in ("rms_roe", "rms_roK", "rms_roOmega")]
            fh.write(f"{it[n]}," + ",".join(f"{x:.9e}" for x in v) + "\n")
    nf = {}
    for f in sorted(run.glob("restart_flow_[0-9]*.csv")):
        A = np.loadtxt(f, delimiter=",", skiprows=1)
        nf[f.name] = int(np.count_nonzero(~np.isfinite(A)))
    itmax = None
    for line in (run / "sst.cfg").read_text().splitlines():
        if line.strip().startswith("ITER="):
            itmax = int(line.split("=")[1])
    out = {"run": run.name, "kind": "su2", "last_iter": it[-1], "ITER": itmax, "reached_ITER": itmax is not None and it[-1] >= itmax - 1,
           "nonfinite_history": nonfin_hist, "nonfinite_restart": nf,
           "nan_verdict": "CLEAN" if nonfin_hist == 0 and not any(nf.values()) else "NONFINITE",
           "convergence": _conv(str(conv_csv))}
    p = OUTD / f"xcheck_gates_{run.name}.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return p


def reduce_run(run: Path, kind: str, last: int = 5) -> Path:
    yb_x, yb = common_yb()
    if kind == "forge":
        steps = sorted(int(re.match(r"res_(\d+)\.h5$", p.name).group(1)) for p in run.glob("res_[0-9]*.h5"))[-last:]
        loader = lambda s: load_forge(run, s)  # noqa: E731
    else:
        steps = sorted(int(re.match(r"restart_flow_(\d+)\.csv$", p.name).group(1)) for p in run.glob("restart_flow_[0-9]*.csv"))[-last:]
        loader = lambda s: load_su2(run, s)  # noqa: E731
    rec = {}
    for k, st in enumerate(steps):
        xy, ro, ux, uy, T, tke = loader(st)
        rec[f"s{k}_nonfinite"] = np.array(int(sum(np.count_nonzero(~np.isfinite(a)) for a in (ro, ux, uy, T, tke))))
        ni, nj, S = mesh_info(len(ro))
        o = reduce_fields(xy, ro, ux, uy, T, tke, kind == "su2", ni, nj, S, yb_x, yb)   # SU2 のエネルギーは k を含む
        for key, v in o.items():
            rec[f"s{k}_{key}"] = np.asarray(v)
        print(f"[reduce] {run.name} {st}: Q_w {float(o['Q_w']) / 1e6:.4f} MW、ṁ 入口 {float(o['mdot_in']):.4f} 出口 {float(o['mdot_out']):.4f}", flush=True)
    p = OUTD / f"xcheck_{run.name}.npz"
    np.savez(p, steps=np.array(steps), kind=np.array(kind), **rec)
    return p


def summarize(name):
    """xcheck_<run>.npz → 判定窓の平均と時系列 (x = 40/70/94 の δ_loc・θ_r・τ_w、Q_w、断熱の T_w、収支)。"""
    Z = np.load(OUTD / f"xcheck_{name}.npz"); n = len(Z["steps"])
    g = lambda k, key: np.interp(XE, Z[f"s{k}_x"], Z[f"s{k}_{key}"])  # noqa: E731
    gw = lambda k, key: np.interp(XE, Z[f"s{k}_x_w"], Z[f"s{k}_{key}"])  # noqa: E731
    d = {key: np.array([g(k, key) for k in range(n)]) for key in ("delta_loc", "theta_r", "delta_loc_s125", "theta_r_s125")}
    d.update({key: np.array([gw(k, key) for k in range(n)]) for key in ("tau_w", "c_f", "q_w", "T_w")})
    d["Q_w"] = np.array([float(Z[f"s{k}_Q_w"]) for k in range(n)])
    for key in ("mdot_in", "mdot_out", "Hdot_in", "Hdot_out"):
        d[key] = np.array([float(Z[f"s{k}_{key}"]) for k in range(n)])
    d["steps"] = np.array(Z["steps"], dtype=float)
    d["nonfinite"] = np.array([int(Z[f"s{k}_nonfinite"]) for k in range(n)])
    return d


def compare() -> dict:
    sys.path.insert(0, str(HERE.parents[1] / "solver_density_cuda/tools"))
    from check_quasisteady import classify
    runs = {"forge_plain_ad": "run_0184_ns_coldmesh_cpg_ad_plain", "forge_plain_tw": "run_0185_ns_coldmesh_cpg_tw300_plain",
            "forge_prod_ad": "run_0186_ns_coldmesh_cpg_ad_dilat2", "forge_prod_tw": "run_0187_ns_coldmesh_cpg_tw300_dilat2",
            "su2_ad": "run_0188_su2_coldmesh_cpg_ad", "su2_tw": "run_0189_su2_coldmesh_cpg_tw300"}
    D = {k: summarize(v) for k, v in runs.items() if (OUTD / f"xcheck_{v}.npz").is_file()}
    out = {"runs": runs, "available": sorted(D), "gates": {}, "metrics": {}}
    ix = [int(np.argmin(np.abs(XE - x))) for x in XS]
    # 前提: 準定常 (5 枚、drift・osc 0.1 %) と収支
    for k, d in D.items():
        bad = []
        for key in ("delta_loc", "theta_r", "c_f"):
            for i in ix:
                v = classify(d["steps"], d[key][:, i], 1.0, 0.001, 0.001, 5)[0]
                if v != "STEADY":
                    bad.append(f"{key} x={XE[i]:.0f} {v}")
        if k.endswith("_ad"):
            for i in ix:
                v = classify(d["steps"], d["T_w"][:, i], 1.0, 0.001, 0.001, 5)[0]
                if v != "STEADY":
                    bad.append(f"T_w x={XE[i]:.0f} {v}")
        else:
            v = classify(d["steps"], d["Q_w"], 1.0, 0.001, 0.001, 5)[0]
            if v != "STEADY":
                bad.append(f"Q_w {v}")
        mb = float(np.mean(d["mdot_in"] - d["mdot_out"]) / np.mean(d["mdot_in"]))
        hb = float(np.mean(d["Hdot_in"] - d["Hdot_out"] - d["Q_w"]) / np.mean(d["Hdot_in"]))
        gp = OUTD / f"xcheck_gates_{runs[k]}.json"
        g = json.loads(gp.read_text()) if gp.is_file() else None
        nonfin = int(sum(int(Z) for Z in d["nonfinite"]))
        g_ok = (g is not None and g["nan_verdict"] == "CLEAN" and not g["convergence"]["diverged"] and not g["convergence"]["rising"]
                and nonfin == 0 and (g["kind"] == "forge" or g["reached_ITER"]))
        out["gates"][k] = {"not_steady": bad, "mass_balance": mb, "enthalpy_balance": hb, "nonfinite_window": nonfin,
                           "gates_file": gp.name if g else "(無い → 判定不能)", "nan": g and g["nan_verdict"],
                           "convergence": g and g["convergence"]["verdict"], "rising": g and g["convergence"]["rising"],
                           "ok": g_ok and (not bad) and abs(mb) <= 1e-3 and abs(hb) <= 1e-3}
    m = lambda k, key: D[k][key].mean(0)  # noqa: E731
    t = (XE >= 40) & (XE <= 94)
    if all(k in D for k in ("forge_plain_ad", "forge_plain_tw", "su2_ad", "su2_tw")):
        Rf = m("forge_plain_tw", "delta_loc") / m("forge_plain_ad", "delta_loc"); Rs = m("su2_tw", "delta_loc") / m("su2_ad", "delta_loc")
        met = {"R_loc_max_rel": float(np.max(np.abs(Rf / Rs - 1.0)[t])), "R_loc_forge": Rf[ix].tolist(), "R_loc_su2": Rs[ix].tolist()}
        for arm in ("ad", "tw"):
            for key in ("delta_loc", "theta_r", "c_f", "q_w"):
                if key == "q_w" and arm == "ad":
                    continue
                a = m(f"forge_plain_{arm}", key); b = m(f"su2_{arm}", key)
                good = t & np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 0)
                met[f"{key}_{arm}_max_rel"] = float(np.max(np.abs(a[good] / b[good] - 1.0))) if good.any() else None
        met["Q_w_tw_rel"] = float(np.mean(D["forge_plain_tw"]["Q_w"]) / np.mean(D["su2_tw"]["Q_w"]) - 1.0)
        a = m("forge_plain_ad", "T_w"); b = m("su2_ad", "T_w")              # 記録のみ (判定に使わない)
        met["T_w_ad_max_rel_record"] = float(np.max(np.abs(a / b - 1.0)[t]))
        lim = {"R_loc_max_rel": 0.02, "delta_loc_ad_max_rel": 0.03, "delta_loc_tw_max_rel": 0.03, "theta_r_ad_max_rel": 0.03,
               "theta_r_tw_max_rel": 0.03, "c_f_ad_max_rel": 0.03, "c_f_tw_max_rel": 0.03, "q_w_tw_max_rel": 0.05}
        met["pass"] = {k2: (met[k2] is not None and met[k2] <= v) for k2, v in lim.items()}
        gates_ok = all(out["gates"][k]["ok"] for k in ("forge_plain_ad", "forge_plain_tw", "su2_ad", "su2_tw"))
        out["VERDICT"] = ("判定不能 (前提不成立)" if not gates_ok else
                          ("forge の解き方は冷却ノズルでも SU2 と一致 (この CPG 条件・格子・比較量で)" if all(met["pass"].values())
                           else "不一致 (forge 側の原因を調べる — 諮問)"))
        out["metrics"]["plain_vs_su2"] = met
    for tag, (a, b) in {"prod_minus_plain_R": ("forge_prod", "forge_plain")}.items():
        if all(f"{p_}_{w}" in D for p_ in (a, b) for w in ("ad", "tw")):
            Ra = m(f"{a}_tw", "delta_loc") / m(f"{a}_ad", "delta_loc"); Rb = m(f"{b}_tw", "delta_loc") / m(f"{b}_ad", "delta_loc")
            out["metrics"][tag] = {"R_prod": Ra[ix].tolist(), "R_plain": Rb[ix].tolist(), "max_rel": float(np.max(np.abs(Ra / Rb - 1.0)[t]))}
    (OUTD / "xcheck.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    for c in ("reduce-forge", "reduce-su2"):
        p = sp.add_parser(c); p.add_argument("run"); p.add_argument("--last", type=int, default=5)
    for c in ("gates-forge", "gates-su2"):
        p = sp.add_parser(c); p.add_argument("run")
    sp.add_parser("compare")
    a = ap.parse_args()
    if a.cmd == "reduce-forge":
        print(reduce_run(HERE / a.run, "forge", a.last))
    elif a.cmd == "reduce-su2":
        print(reduce_run(HERE / a.run, "su2", a.last))
    elif a.cmd == "gates-forge":
        sys.path.insert(0, str(HERE))
        print(gates_forge(HERE / a.run))
    elif a.cmd == "gates-su2":
        print(gates_su2(HERE / a.run))
    else:
        compare()
```

## 参考: `case/45.isobutane_m6_d155/README.md`

```
# case/45 — イソブタン燃焼ガス M6 風洞・出口径 1.55 m の最短全長スタディ

Pt 5.5 MPa / Tt 1600 K / φ=0.9 燃焼ガス (semi-perfect NASA-9) / M_design 6 /
**出口径 1.55 m は δ\* 込み物理壁で定義**。全長 (スロート→物理出口) の最短化と
粘性壁 (δ\* 物理壁 + SST NS) までの検証。計画:
[`plans/accepted/design-isobutane-m6-d155.md`](../../plans/accepted/design-isobutane-m6-d155.md)。

- **現行の生産の問題 (2026-10-07 ユーザ決定)**: `problem_d155_ns_prod.yaml` (+ 凝縮 `problem_d155_ns_prod_cond.yaml`、`make_prod_problems.py` で N2 から生成)。
  上流は 5 次多項式 (`pw_upstream: poly`)、MOC は analytic + converge、物理壁は全域 1 本の B スプライン (`physical_wall_repr: single_bspline`)。
  較正値 6.8825162455159465e-06・r_t 0.07666536551630307 m・k_f 1.054129117086371。生産の run は dry `run_0167_ns_n012_N2` + 延長 `run_0179_ns_n012_N2_ext`、
  凝縮 `run_0170_ns_n012_N2_cond`。旧生産 (`problem_d155_ns_finemesh_recal_final_mono.yaml`、run_0147+0149/0148) は置き換えた。
  生産の確認 (準備が run_0167 の入力をビット同一で再現) `_band_ab/prod_confirm/PROD_CONFIRM.json`、STEP `_band_ab/prod_confirm/step/wall_physical.step`、標準報告 `_band_ab/prod_confirm/reports/*/…_report.pptx`。
- 最短探索: `search_shortest.py` → `search_shortest.json`
  (基準 = shortest-robust study: margin≥1°/topo≥0.02/hard gate/単峰、n1200 確認)。
  **勝者 R2 / L_c 39.3 / M_K 2.7 → x_F = 95.104 r_t** (R3: 39.7/2.8 → 95.433)。
- r_t = (0.775 − δ\*_exit)/(r_F/r_t), 初期推定 0.0771 m (δ\*_exit ≈ 0.66 r_t)。
  A/A\*(M6, Tt1600) = 88.21。
- 凝縮事前見積り: 出口 T 233.5 K vs Tsat(H₂O) 263.9 K = **−30 K 過冷却** →
  dry 本体 + 採用点の凝縮 ON 再評価 (方針 a)。

> **2026-09-04 更新**: 粘性トリム (Md 6.0144, run_0005–0007) は撤回。新チェーン (積分法初期壁 + 固定 Euler 基準・帯局所抽出) で
> **Md 6.0 トリムなし**のまま ṁ_NS/ṁ_E 1.0002・出口面コア M 6.0002 (`run_0038_ns_final_rt77p02` = **最終形**: r_t 77.02 mm、出口半径 0.7750 m、全長 7.325 m; 点列 `points_d155_final_ns.csv`)。旧 v3 のスロート δ\* は実効値の 8 倍過大だった
> ([調査ノート](../../notes/investigations/nozzle-deltastar-throat-review.md), [plan](../../plans/accepted/tooling-nozzle-deltastar-core-matched-euler.md))。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_euler_shortest_dry` | 最短点 (R2/L_c39.3/M_K2.7) の node Euler TP dry (`problem_d155_R2_Lc39.3_tp.yaml`, split_h2o, 段階起動 12000 step) | **dM_max 0.036 % M_d・軸出口 M 6.00003・overshoot 0.034 %・ε_M_rms 0.0057 %・ε_θ 0.0054°・コア 64/65 面**。品質 PASS (AR OK/skew 0.44)・NaN 0・quasisteady ALL STEADY・軸 M 凍結 8k→12k 5.6e-5。残差は warm 床 (rms_ro 1.1e-6, `check_convergence` NOT CONVERGED 判定 = M6 系列の既知パターン)。`metrics.json`/`residual_history.png` | active (**最短点 Euler 検証の正本**) |
| `run_0002_ns_coarse` | NS 中継 (y+~50, `problem_d155_ns_coarse.yaml`, 物理壁 δ\* v1 相関, IC=run_0001) | 完走 12000 step・NaN 0。残差 2–3 桁降下 (ro/roe は低レベルでプラトー = 中継品質として想定内)。出口壁半径 0.7771 m (v1 相関 δ\* ≈ 0.71 r_t) | active (中継) |
| `run_0003_ns_v1` | NS 本計算 v1 (y+~1.4, `problem_d155_ns.yaml`, 相関 δ\*, 24000 step, IC=run_0002) | 完走・NaN 0・品質 SOFT-PASS。**軸出口 M 6.00098・dM_max 0.327 % M_d (ゲート内)**・ε_M_rms 0.075 %。残差 2 桁降下 warm 床。δ\* v3 抽出元 (`dstar_v3.csv`: 相関は中流過大 median 比 0.909, x14 で 0.669) | active (v1 記録・v3 抽出元) |
| `run_0004_ns_v3` | NS 本計算 v3 (CFD 抽出 δ\* 全域採用, IC=run_0003) | 完走・NaN 0・品質 SOFT-PASS・16k→24k 軸 M 凍結 1.2e-4。**δ\* 固定点確認 (2巡目比 median 1.001 [0.981,1.005])**。出口面軸 M **5.9856 (−0.24 % M_d)**、x_E ディップ 5.956 (law 側残差)、一様区間 5.975–5.996 のうねり。metrics の `M_axis_exit` は x_E 評価値であることに注意 | active (**dry NS 固定点・トリム根拠**) |
| `run_0005_euler_trim` | **粘性トリム版** Euler (`problem_d155_trim_tp.yaml`: Md 6.0144, L_c 39.7/M_K 2.7 [トリム後最短 x_F 96.00 r_t], r_t 0.076511) | dM_max 0.036 % M_d・x_E で M 6.01445 (=目標)・ε_M_rms 0.0056 %・コア 64/65。NaN 0 | active (トリム Euler) |
| `run_0006_ns_trim_v1` | トリム版 NS v1 (相関 δ\*, `problem_d155_trim_ns.yaml`, IC=run_0004 cross-mesh) | 完走・NaN 0・凍結 2.7e-4。出口面軸 M **5.9810** (= 旧 v1 5.9667 + トリム 0.24 % — トリムは線形に作用)。δ\* v3' 抽出元 | active (トリム v1・抽出元) |
| `run_0007_ns_trim_v3` | トリム版 NS v3 (CFD 抽出 δ\*, IC=run_0006) | 完走・NaN 0・凍結 1.2e-4・**δ\* 固定点 (2巡目比 median 0.999)**。**出口面軸 M = 6.00079 (+0.013 % M_d)** — トリム狙い通り。quasisteady machmax/pmax STEADY (shock 指標は無衝撃のため対象外)。出口壁半径 0.77591 m → 最終 r_t 補正 −0.12 % で 0.775 m 合わせ | active (**dry 最終形の正本**) |
| `run_0017_laminar_coarse` | **モデル梯子用 laminar NS** (`model: "none"` 化 + warm interp, coarse y+~50, 18000 step 段階起動) — cfl 上限スイープの probe 基地 | 完走・NaN 0。probe: cfl4 OK / 8@127 / 16@16 (SST・Euler と同一境界) | active (モデル梯子) |
| `run_0018_ns_trim_cfl6_r07` | **生産推奨設定の段階起動検証** (run_0007 と同一入力: IC=run_0006, dstar_v3+blend(−1,−0.5), 24000 step。差分は本段 `cfl 6 + implicitRelax 0.7` のみ) | 完走・NaN 0 (P min 2188 Pa)。**metrics が run_0007 と 5 桁一致** (出口軸 M 5.970306 vs 5.970316)。quasisteady machmax/pmax STEADY。残差は run_0007 の 24000 step 到達水準に 6.4k–11.3k step で到達 (2〜3.8 倍速)・roK/roOmega はさらに 5.8×/3.7× 深い。本段 192.9 s (cfl1 は 184.1 s = relax の追加コスト +5 %) | active (**生産設定検証の正本**) |
| `run_0019_ns_cm_pass1` | **排除厚さ更新計画 V1 pass 1** ([plan](../../plans/active/tooling-nozzle-deltastar-core-matched-euler.md)): pass 0 = `run_0004` の場から固定 Euler (`run_0001`) 基準で δ_r を抽出 (この pass はコア全体 α 方式)、ω=0.5、半径方向オフセット、**Md 6.0 トリムなし**、IC=run_0004 cross-mesh、`problem_d155_ns.yaml`、24000 step | 完走・NaN 0・品質 SOFT-PASS・quasisteady STEADY・残差 2 桁 warm 床 (NOT CONVERGED 判定 = 系列共通)。**ṁ_NS/ṁ_E 1.020 → 1.010** (スロート補正 0.0114 → 0.0065 r_t)、x≤30 の軸残差 −0.4 → −0.1 %。しかし試験部軸 M に +0.8 % (x 55–75) / −1.6 % (x 85–94) の波 = コア全体 α 抽出が上流壁誤差の波を欠損に取り込んだ (x=20–50 で δ +15 %)。→ 抽出器を帯局所参照に変更 (plan §4.3) | active (pass 1 記録・帯局所化の根拠) |
| **`run_0020_ns_cm_pass2`** | **V1 pass 2**: `run_0019` の場から**帯局所参照** (`band_local_deficit`) で δ_r を全域抽出、ω=1.0、半径方向オフセット、Md 6.0 トリムなし、IC=run_0019 | 完走・NaN 0・品質 SOFT-PASS・quasisteady ALL STEADY・残差 2 桁 warm 床 (系列共通)。**ṁ_NS/ṁ_E = 0.9999** (スロート補正 0.0014 r_t = 質量流量実効値)・**出口面コア M 6.0008〜6.0011 (+0.01 %) — トリムなしで目標達成**。軸 M は x=20–75 で Euler 設計 +0.05〜0.17 %、x≥80 の軸上に ±0.5〜1 % の波 (pass1→2 の壁変化 x≈40 が末端軸へ届いたもの、コア平均には出ない)。固定点差: use/in median 1.003、x=40–60 で +4 % → pass 3 で収縮確認 | active (**新チェーンの到達点 (暫定)**) |
| `run_0021_ns_cm_pass3` | **V1 pass 3**: `run_0020` の場から帯局所抽出、ω=0.5、IC=run_0020 (収縮確認) | 完走・NaN 0・SOFT-PASS・ALL STEADY・残差 warm 床。**ṁ_NS/ṁ_E = 1.0000、出口面コア M 6.005/5.996/5.998 (x=80/90/94, ±0.08 %)**、軸波 ±1 % → ±0.5 % に減衰。固定点差は x=20–60 で ±5〜7 % (帯幅依存の下流揺らぎ、帯パラメータ再検討中) | active (pass 3) |
| **`run_0022_ns_ib_pass0`** | **V2 pass 0: 積分法 (CONTUR 運動量積分, 断熱壁) 初期壁のみ、CFD 帰還なし** (`--init-integral`, 半径方向オフセット, Md 6.0 トリムなし, IC=run_0021 cross-mesh) | 完走・NaN 0・SOFT-PASS・ALL STEADY・残差 warm 床。**ṁ_NS/ṁ_E = 0.9991、出口面コア M 6.0026〜6.0037 (+0.04〜0.06 %)、軸 M は x=20〜94 で Euler +0.09〜0.21 % で平坦 (波なし)**。抽出 δ_r vs 積分法 δ_r: x≤60 で ±2 %、x=90 で −1.7 % (固定点にほぼ乗っている)。初期壁だけで主ゲート達成 | active (**積分法初期壁の到達点**) |
| **`run_0023_ns_ib_pass1`** | **V2 pass 1**: `run_0022` (積分法初期壁) の場から帯局所抽出 (適応帯)、ω=1.0、単調性ガード (λ=1.0 追加平滑化)、IC=run_0022 | 完走・NaN 0・SOFT-PASS・ALL STEADY・残差 warm 床。**ṁ_NS/ṁ_E = 0.9999、出口面コア M 6.0008〜6.0013 (+0.01〜0.02 %)、固定点差 use/in p10–p90 0.997〜1.002 (収束)**。軸 M は x=20–94 で Euler −0.14〜+0.45 % (x=70 に小さな山)。V1 pass 3 (相関初期壁から 3 pass) と δ_r が x=20–60 で 3 % 以内で一致 = **初期推定器に依らない固定点** | active (**新チェーンの正本 (積分法初期壁 + 1 pass)**) |
| `run_0024_ns_cm_pass4_q` | **δ_r 平滑化 (5 次 P-spline) の検証, V1 系 pass 4**: `run_0021` の場から帯局所抽出 → 5 次 P-spline (3 階差分ペナルティ, ノット 2 r_t, λ=1) → ω=0.5 (緩和後も再平滑化) | 完走・NaN 0・SOFT-PASS・ALL STEADY。**壁曲率の高周波ノイズ 7.7e-3 → 3.8e-4 [1/r_t]** (積分法初期壁と同水準)。ṁ 0.9997、出口面コア M ±0.1 %、軸 M の試験部最大偏差 0.49 % (run_0021 の 0.40 % と同等 = 軸の波は壁の凸凹ではない) | active (平滑化検証) |
| **`run_0025_ns_ib_pass2_q`** | **同 V2 系 pass 2**: `run_0023` の場から抽出 → P-spline → ω=1.0 | 完走・NaN 0・SOFT-PASS・ALL STEADY。**ṁ 1.0002、出口面コア M +0.01 %、固定点差 p10–p90 0.999〜1.006 (収束)**。軸 M は x=70 に +0.42 % (run_0023 の +0.45 % と同じ) — 滑らかな壁でも残るので、積分法壁→抽出壁の ~1 % の δ 差 (x≈15–25) への軸の集束応答。コア平均は平坦 | active (**生産形の最終 (平滑化込み)**) |
| `run_0026_ns_ib_pass3_nosoft_cfl5` | **起動レシピ A/B (a)**: run_0025 の固定点から `--stages none` (soft/mid なし) + **cfl 5 + implicitRelax 0.7**、24000 step | 完走・NaN 0・SOFT-PASS・ALL STEADY。metrics は run_0025 と同一 (ṁ 1.0002、出口面コア M 6.0000)、roK は 6e-5 と 3 倍深い。**全残差列が run_0025 (cfl1, 30000 step) の最終水準に 5000〜8900 step で到達、出口 M は 8000 step で凍結** → 12000 step で十分 (NS ≈ 95 s)。NS 壁時計 191 s (24000 step) | active (**warm start の生産レシピ根拠**) |
| `run_0027_ns_ib_pass3_ramp_cfl5` | **起動レシピ A/B (b)**: `--stages ramp --ramp 1,2,3.5 --ramp-steps 1000` → 本段 cfl 5 + relax 0.7 | 完走・NaN 0・SOFT-PASS。metrics は (a) と同一 (差 1e-4 以下)。ランプの利得なし (warm start では不要; cold start 用の保険として残す)。本段 step が outStepInterval の倍数でなく最終 res が 16000 止まり → runner 側で丸めるよう修正済み | active (A/B 記録) |
| `run_0028_euler_cfl6_r07` | **Euler の CFL 引き上げ**: run_0001 と同一設計、`--cfl 6 --implicit-relax 0.7 --stages soft` (soft 3000 + 本段 12000; run_0001 は soft+mid+本段 cfl2 18000 step) | 完走・NaN 0・`check_convergence` **ALL PASS**・ALL STEADY。metrics は run_0001 と同一 (dM 0.0363 % vs 0.0361 %、出口軸 M 6.00007 vs 6.00003、ε_M 0.0058 %)。残差は run_0001 最終水準に 4700〜8900 step で到達 → 本段 8000 step で可 (run_0032 で確認) | active (Euler 起動レシピ) |
| **`run_0029_ns_ib_pass0_norelay`** | **中継なし NS pass 0**: 積分法初期壁、IC = Euler `run_0001` を直接 y+~1 メッシュへ cross-mesh (k/ω 貼り + ω 床)、既定 `full` 起動 (soft/mid + 本段 cfl1 24000) | 完走・NaN 0・SOFT-PASS・ALL STEADY。**ṁ 0.9999、出口面コア M +0.00〜+0.04 %、固定点差 0.985〜1.002 = 中継あり (run_0022) と同等 → 中継 (y+~50 の踏み台 run) は不要**。NS 壁時計 514 s は別セッションの GPU 共有中の値 (通常 ~190 s) | active (**中継不要の根拠**) |
| `run_0030_ns_ib_pass0_norelay_cfl5` | 中継なし + soft/mid なし + 本段 cfl 5 (cold start を一気に) | **step 1 で NaN** → cold start には 1 次・低 CFL の soft/mid が必要。ディレクトリ削除済み (記録のみ、`_logs/run_0030.log`) | 削除済 |
| **`run_0031_ns_ib_pass0_norelay_full_cfl5`** | **cold start の短縮版**: 中継なし (IC = Euler run_0001)、`--stages full` (soft/mid 6000) + **本段 cfl 5 + implicitRelax 0.7、12000 step** | 完走・NaN 0・SOFT-PASS。**ṁ 1.0000、出口面コア M +0.00〜+0.03 %、固定点差 0.982〜0.998 = run_0029 (cfl1 24000) と同等**。outStepInterval 8000 が 12000 を割り切らず最終 res が 8000 止まり (prepare_ns 側で修正済み) → quasisteady はスナップショット不足で判定不能、残差は roK 2.8 桁 falling。NS 壁時計 411 s は GPU 共有中の値 | active (**cold start 生産レシピ候補**) |
| `run_0032_euler_cfl6_r07_8k` | Euler `--cfl 6 --implicit-relax 0.7 --stages soft --steps 8000` (run_0028 の短縮版) | 完走・NaN 0・品質 PASS。metrics は run_0001 と同一 (dM 0.0361 %、出口軸 M 6.00005、ε_M 0.0057 %) だが `check_convergence` は未達判定 (12000 step の run_0028 は ALL PASS) → Euler 本段は 12000 を推奨 | active (Euler 短縮の記録) |
| `run_0033_euler_Lpipe10` | 上流配管の影響テスト用 Euler: `problem_d155_R2_Lc39.3_tp_Lpipe10.yaml` (直管 L_pipe 0.5 → 10 r_t)、cfl 6/relax 0.7/soft | 完走・NaN 0 (run_0034 の固定 Euler 基準) | active (基準) |
| `run_0035_euler_rt77p05` / `run_0036_ns_final_rt77p05` | 最終 r_t 補正の 1 回目: r_t 77.05 mm (run_0025 の**生抽出** δ_r(x_F−0.3)=0.684 で解いた値) の Euler + NS (run_0025 から warm start, none/cfl5/12000) | 完走・NaN 0・SOFT-PASS・ALL STEADY。ṁ 1.0002、出口面コア M 6.0000、固定点 0.996〜1.002。**物理出口半径 0.7753 m と 0.3 mm 残った** = 壁に載るのは平滑化後の δ_r,exit 0.688 なので、solve_rt をその値を使うよう修正して r_t 77.02 mm でやり直し (run_0037/0038) | active (補正 1 回目の記録) |
| `run_0034_ns_ib_pass0_Lpipe10` | **上流配管の影響**: 直管 L_pipe 0.5 → 10 r_t (入口境界層が 6 倍厚い)、他は run_0031 と同じ (中継なし, full, cfl5, 12000; ni 1400) | 完走・NaN 0・SOFT-PASS・ALL STEADY。**スロート以降の δ_r は L_pipe 0.5 (run_0031) と 0.3 % 以内で同一 (x=0: 0.0014/0.0014, x=90: 0.659/0.661)、出口面コア M +0.00〜0.03 % も同一**。ṁ 比 0.9992 (スロート実効 δ 0.0018 vs 0.0014 = 0.03 mm)。→ 収縮部の加速で上流の境界層履歴は消える (積分法の θ0 無記憶と整合) | active (**配管非依存の根拠**) |
| `run_0037_euler_rt77p02` | 最終 r_t 77.02 mm の Euler (固定基準, cfl 6/relax 0.7/soft) | 完走・NaN 0・PASS・**ALL PASS**・ALL STEADY | active (最終形の基準) |
| **`run_0038_ns_final_rt77p02`** | **最終形**: r_t 77.02 mm (solve_rt: run_0025 の壁に載せた δ_r,exit 0.688 で 0.775 m 合わせ)、δ_r = run_0025 抽出 (P-spline)、warm start (none/cfl5/relax0.7/12000) | 完走・NaN 0・SOFT-PASS・ALL STEADY (残差 warm 床)。**物理出口半径 0.7750 m (spec 一致)、ṁ_NS/ṁ_E 1.0002、出口面コア M 6.0002 (+0.00 %)、固定点 0.993〜1.002、全長 (スロート→出口) 7.325 m、物理スロート半径 77.14 mm**。点列 `points_d155_final_{ns,euler}.csv` | active (**最終形の正本**) |
| **`run_0039_ns_final_cond`** | **最終壁 (run_0038) の凝縮 ON restart** (`problem_d155_ns_rt77p02_cond.yaml`: Kw+HK condModel1+Kantrowitz, 蒸発 ON; IC=run_0038 同一メッシュ, none/cfl1, 12000 step; 評価 `eval_cond.py`) | 完走・NaN 0・SOFT-PASS・**ALL STEADY (4k/8k/12k で同一)**。軸 onset x≈57.9 r_t、出口 (x=94.1) g 0.27 % (H₂O の 3 %)、S_max 17.4。**出口軸 M 5.918 (−1.36 %)、出口面コア M 5.989 (−0.18 %)** (dry run_0038: 6.0017 / 6.0002)。旧トリム壁 run_0008 (onset 68.6, 軸 −0.99 %, コア −0.17 %) と同じ結論 = Tt 1600 K では凝縮不可避、影響は軸に集束して見えコア平均は −0.2 % | active (**最終形の凝縮評価**) |
| `run_0042_ns_restart_ctrl` | **軸 M の山の格子感度 A0′** ([plan](../../plans/active/verification-m6-axis-wave-mesh-su2.md)): run_0038 と同一メッシュ (nj 97, 軸側 0.075 r_w) を `restart_field` で継続、新バイナリ (`bc1c84d1`, AWS `~/forge-axiswave`)、12000 step・出力 1000 毎 | NaN 0・SOFT-PASS・NOT CONVERGED (収束場からの継続で全列プラトー, rms_ro 4.0e-7)。軸の山 b(0)=0.2228 %pt @x=69.4、末尾変動 0.0007。run_0038 比で M が全 η で −0.08〜−0.11 % (バイナリ差: limiterScaled 既定変更後) | active (腕 A の基準) |
| `run_0040_ns_axgap025` | 同 A1: `axis_gap_frac` 0.025 (nj 110, η≲0.69 一様)、IC=run_0038 を interp_field | NaN 0・SOFT-PASS・NOT CONVERGED (roUy プラトー 2.9 桁)。b(0)=0.1877 %pt @69.2、末尾変動 0.0006 | active (腕 A) |
| `run_0041_ns_axgap0125` | 同 A2: `axis_gap_frac` 0.0125 (nj 141, η≲0.84 一様)、IC=run_0038 を interp_field | NaN 0・SOFT-PASS・**ALL PASS** (3.1–3.4 桁)。b(0)=0.1963 %pt @69.3、末尾変動 0.0004。A0′ 比: 軸ノードが一様に −0.12 %、x∈[60,80] の局所差 max 0.041 % (η=0)/≤0.021 % (η≥0.05) → 事前登録基準では**判定保留** | active (腕 A) |
| `run_0043_ns_restart_ctrl_rep` / `run_0044_ns_axgap0125_rep` | 再現性 A/B: A0′ (run_0042) と A2 (run_0041) を元の入力から再実行 (同一バイナリ sha256 6d5933e1…, BLOCKSIZE 128) | NaN 0。b(0) 末尾平均 0.22296 / 0.19619 (元 0.22268 / 0.19620)、Δb −0.02677 (元 −0.02648)、D̃ の再実行差 ≤ 0.0011 % → 許容 0.003 %pt 内で**再現**。判定は保留のまま | active (腕 A 再現性) |
| `run_0045_ns_band_adaptive` | 帯修正の CFD A/B **B0** ([plan](../../plans/active/verification-m6-axis-wave-mesh-su2.md) §4.4): run_0042 の場から現行 adaptive 抽出で作った壁、IC=run_0042 を interp_field、12000 step (`bc1c84d1`) | NaN 0・SOFT-PASS・NOT CONVERGED (横ばい)。P(0.1) 0.350 %・b70(0) 0.359 %pt (軸 x≈56 +0.48 %、x≈84 −0.39 %)。延長は run_0048 | active (帯 A/B) |
| `run_0046_ns_band_edge` | 同 **B1**: 縁アンカー・x 平滑帯 (band_select="edge", c 1.25) で作った壁、他は B0 と同一 | NaN 0・SOFT-PASS・NOT CONVERGED (横ばい)。**P(0.1) 0.079 %・b70(0) 0.013 %pt** (局所の山谷なし、x≈48 の小山 +0.035 % と x 56→80 の −0.04 % の傾きのみ)。場から再抽出で δ_r/δ_in 1.0017 (固定点)。事前登録の判定は保留 (P(0.1) ≤ 0.06 % 不成立) | active (帯 A/B) |
| `run_0047_euler_rt77p02_newbin` | run_0037 (Euler 参照) を restart_field で新バイナリ `bc1c84d1` 継続 12000 step | NOT CONVERGED (rms_ro 3.3e-6 横ばい)。旧 Euler との差は試験部で ≤0.01 % | active (参照の同一バイナリ化) |
| `run_0048_ns_band_adaptive_ext` | B0 (run_0045) を restart_field で 6000 step 延長 (事前登録の延長規定) | 末尾 5 枚の変動 ≤0.0008、P(0.1) 0.350 % (不変) | active (帯 A/B) |
| `run_0049_ns_contur_cal_full` | CONTUR (積分法) を B1 の抽出 δ に x∈[8,90] で較正した壁 (k_f 1.050・k_N 0.999・a 1.0、壁入力は 5 次補間/関数渡し)、IC=run_0046、12000 step ([plan](../../plans/active/verification-m6-axis-wave-mesh-su2.md) §5.1 #8f) | NaN 0・SOFT-PASS・NOT CONVERGED (横ばい)。r/r_w=0.1 で波 0.007 %・**オーバーシュート +0.064 %** (不合格)・出口コア M 6.0034 | active (CONTUR 較正) |
| `run_0050_ns_contur_cal_exit` | CONTUR を**出口の δ だけ**に合わせた壁 (ユーザ提案; k_f 1.027 のみ)、他は run_0049 と同一 | NaN 0・SOFT-PASS・NOT CONVERGED (横ばい)。**波 0.007 %・オーバーシュート +0.019 %・出口コア M 6.0001** (事前登録の判定 合格)。出口半径 0.7771 m (r_t 補正前) | active (**CONTUR 出口較正の候補**) |
| `run_0051_ns_final_c2` (旧最終形、Hall 初期線・補間壁) | **最終設計 (C2 方式)**: CONTUR を出口 δ に較正した壁 (k_f 1.0257)、**r_t 76.807 mm** (出口半径 0.774955 m)、全長 7.305 m、IC=run_0050、12000 step (`problem_d155_ns_c2final.yaml`、[plan](../../plans/active/verification-m6-axis-wave-mesh-su2.md) §5.1 #9) | NaN 0・SOFT-PASS・NOT CONVERGED (横ばい)。波 0.003 %・オーバーシュート +0.018 %・**出口コア M 6.00005**・出口 δ_E/δ_C 0.9994・ṁ 比 0.99932 — 事前登録の判定すべて合格。報告 `report/run_0051_ns_final_c2_report.pptx` | active (**最終形の正本**) |
| **`run_0052_ns_final_c2_cond`** | 最終設計の**凝縮 ON** (run_0051 を restart_field、cfl 1、12000 step、run_0039 と同設定) | NaN 0。凝縮の始まり (軸) x=59.4 (run_0039 57.8)、出口 g 軸 0.0032・コア 0.00035、最大 S 16.5、出口コア M 5.9858。報告 `report/run_0052_ns_final_c2_cond_report.pptx` | active (**最終形の凝縮評価**) |
| `run_0053`〜`0055_euler_wallfit_interp_r{1,2,3}` | 壁表現 Euler A/B の A (補間壁、MOC 2400 点)、同一入力の 3 回再実行。問題 `problem_d155_euler_c2final_n2400.yaml`、等エントロピー IC + soft 3000 → 本段 cfl 6/relax 0.7 × 12000 (plan verification-m6 §5.1 #15) | 完走・残差 NOT CONVERGED (plateau)・評価量は全量 STEADY。判定は `_band_ab/wallfit_euler_ab.json` | active |
| `run_0056`〜`0058_euler_wallfit_fit_r{1,2,3}` | 同 B (位置+壁角の同時当てはめ壁 λ=1e-9)、3 回再実行。壁節点の B−A 最大 1.96e-5 r_t (x=0.022) | 同上 | active |
| `run_0059`〜`0064_*_ext6k` | 事前登録の延長 (準定常の前提未達 → 6000 step × 1 回): 0053〜0058 の res_12000 から restart_field (9 量ビット一致) | 残差 NOT CONVERGED (plateau)。**判定: 保留** (M 波の末尾変動 T 3.6e-4 > Δ/4、U > Δ/2 の量あり)。オーバーシュートは B が +0.0029 %pt (η=0)・+0.0023 (η=0.1)、再実行幅 ≤1.4e-4。再集計 `_band_ab/wallfit_euler_ab_v2.json`・刻み診断 `_band_ab/wallfit_euler_ab_diag_fixedcoef.json` でも保留 (残りは M 波の時間変動)。**生産は補間壁を維持** (plan §9)。2026-10-05: AWS のディスク逼迫のため 0053〜0055・0059〜0061 の中間 res_*.h5 を削除 (res_0・最終場・残差・VERDICT・wallfit_series.csv は保持) | active |
| `run_0065`〜`0067_euler_wallfit_v4_r{1,2,3}` + `run_0068`〜`0070_*_ext6k` | 壁表現 V4 (始点 r′(0)・r″(0) 自由の当てはめ + 縮流部 Hermite を (1, tan 0.091°, 0.4736) に接続)、3 回再実行 + 延長 (plan verification-m6 §5.1 #16 段 2) | 残差 NOT CONVERGED (plateau)。**縮流部の形が V0 から最大 0.069 r_t (x=−4.8) 動いていた** (端曲率 0.5→0.474 が 12 r_t の Hermite で増幅) ため対 V0 比較は交絡。試験部 P の傾き η0 +0.154 → −0.231 %pt、オーバーシュート η0 0.035 → 0.022 %、出口コア M −6e-5。判定 `_band_ab/wallfit_euler_ab_{fit,interp}_vs_v4_diag_fixedcoef.json` (保留) | active |
| `run_0071`〜`0073_euler_wallfit_v4b_r{1,2,3}` + `run_0074`〜`0076_*_ext6k` | 壁表現 V4b (x<−1.5 は V0 と同一の縮流部、[−1.5,0] だけ 5 次 Hermite で V4 の始点 (r′ 0.091°, r″ 0.4736) へ C² 接続、x≥0 は V4 と同一)、3 回 + 延長 | 残差 NOT CONVERGED (plateau)。V0 → V4b → V4: オーバーシュート η0 0.0352→0.0285→0.0224 %・η0.1 0.0100→0.0030→0.0023、試験部 P の傾き η0 +0.154→−0.146→−0.231・η0.1 +0.056→−0.160→−0.230 %pt。判定 `_band_ab/wallfit_euler_ab_{v4_vs_v4b,fit_vs_v4b}_diag_fixedcoef.json` (登録の「中間」) | active |
| `run_0077`〜`0079_euler_wallfit_pin_r{1,2,3}` + `run_0080`〜`0082_*_ext6k` | CFD ピン (plan tooling-nozzle-cfd-pinned-initial-line): 初期線・m*・アンカーを V0 場 run_0062 res_6000 から凍結、V0 型当てはめ壁、3 回 + 延長 | 残差 plateau。V0 比: |P 傾き| η0 0.154→0.120・η0.1 0.056→0.042 (改善)、出口 M 規格化オーバーシュート η0 0.0372→0.0356・η0.1 0.0120→0.0075 (改善)、生のオーバーシュート η0 0.0352→0.0425 (悪化)、出口コア M 5.99988→6.00042 (|M−6| 悪化)、M 波 η0 0.0053→0.0069・P 波 η0 0.034→0.043 (悪化方向)。登録判定 V3 不採用・V2 超過 (`_band_ab/wallfit_euler_ab_fit_vs_pin_diag_fixedcoef.json`, `_band_ab/cfdpin_v2.json`) | active |
| `run_0083`〜`0085_euler_wallfit_pincal_r{1,2,3}` + `run_0086`〜`0088_*_ext6k` | CFD ピン + 出口較正 (M_design 5.999584、他は pin と同一; plan tooling-nozzle-cfd-pinned-initial-line §6 V3′)、3 回 + 延長 | 残差 plateau。出口コア M 6.000000。V0 比: |P 傾き| η0 0.154→0.121・η0.1 0.056→0.042、オーバーシュート η0 0.0352→0.0356・η0.1 0.0100→0.0074、出口規格化オーバーシュート η0 0.0372→0.0356・η0.1 0.0120→0.0074、P 波 η0 0.034→0.043 (非劣化を示せず)、M 波 η0.1 0.0060→0.0064 (非劣化を示せず、T 超過)。登録判定 保留 (`_band_ab/wallfit_euler_ab_{fit,pin}_vs_pincal_diag_fixedcoef.json`) | active |
| `run_0089_ns_c2pin_pass1` | P5 ①: CFD ピン設計 (initial_line cfd・joint・Md_moc_offset −4.16e-4・物理壁解析経路 A″)、r_t 76.8075 mm・k_f 1.02573、等エントロピー IC + 段階起動 full (`problem_d155_ns_c2pin.yaml`, plan tooling-nozzle-cfd-pinned-initial-line §5.1 #7) | 残差 plateau。E で出口 δ_E 0.7198 r_t → solve_rt r_t 76.773 mm・k_f 1.03287 (`_band_ab/c2pin_solve.json`)。出口コア M 5.9971 (未発達の可能性) | active |
| **`run_0090_ns_c2pin_final`** | P5 ③: 最終 NS (r_t 76.773 mm・k_f 1.03287、IC=run_0089、`problem_d155_ns_c2pin_final.yaml`) | NaN 0・SOFT-PASS・NOT CONVERGED (plateau)。出口半径 0.7749999 m、δ_E/δ_C(x_F) 1.0036、波 η0.1 0.0034 %・オーバーシュート η0.1 0.0044 %・**出口コア M 5.99848 (登録 6.000 ± 2e-4 を外れる、末尾も +9e-5/4000 step で漸増)**、ṁ 比 0.99939。報告 `report/` | active |
| `run_0091_ns_c2pin_final_cond` | P5 ④: 凝縮 ON (run_0090 から convert_species_field conserve、cfl 1、12000 step) | NaN 0。凝縮の始まり x 59.66 (run_0052 59.42)、出口 g 軸 0.0030・コア 0.00033、最大 S 16.45、出口コア M 5.9855 (run_0052 5.9858) | active |
| `run_0092_ns_c2pin_pass2` | P5 手 1: C2 の 2 pass 目 (run_0090 から E で δ_E 0.7225 → k_f 1.03736・r_t 76.753 mm、IC=run_0090、12000 step; `_band_ab/c2pin_solve_pass2.json`) | NOT CONVERGED (plateau)。出口コア M 5.99926 (漸近 ≈ 5.9993)、δ_E/δ_C(x_F) 1.0019、波 η0.1 0.0034 %・オーバーシュート η0.1 0.0140 % (末尾変動 0.0029)、出口半径 0.775 m | active |
| `run_0093_ns_c2pin_pass2_cond` | 最終設計 (run_0092 の壁) の**凝縮 ON** (convert_species_field で run_0092 から、cfl 1、12000 step) | NaN 0。凝縮の始まり x 59.39、出口 g 軸 0.0032・コア 0.00035、最大 S 16.50、出口コア M 5.98548 (値は登録範囲内だが保存 3 枚で準定常判定未完了・凝縮残差 RISING)。報告 `report/run_0093_ns_c2pin_pass2_cond_report.pptx` | active (**最終形の凝縮評価**) |
| **`run_0094_ns_c2pin_pass2_ext6k`** | **暫定最終 (CFD ピン; 壁解像 FAIL・δ_E 取り違えの訂正あり、plan §9 2026-10-05)**: run_0092 を登録どおり 6000 step 延長 (restart_field 9 量ビット一致、1000 step ごと出力)。壁点列 `points_d155_final_c2pin_ns.csv` (r_t 76.7531 mm、出口半径 0.775 m、全長 7.31 m) | NOT CONVERGED (plateau)。評価量は全量 STEADY: 出口コア M 5.99926 (−0.012 %)、波 η0.1 0.0036 %・オーバーシュート η0.1 0.0137 %、δ_E/δ_C(x_F) 1.0038 (訂正; 1.0019 は緩和後の値の取り違え)、壁解像 FAIL (y1+>1 が 29.4 %、run_0051 は 31 %)。報告 `report/run_0094_ns_c2pin_pass2_ext6k_report.pptx` | active (**最終形**) |
| `run_0095_ns_finemesh_pass` | #11 ① 1 回目: 細分メッシュ (`problem_d155_ns_finemesh_pin.yaml`、ni 2000 × nj 97、スロート第 1 セル 4.5e-6 r 比、AR max 4237 で `--ar-max 5000` PASS) に run_0092 を interp_field、本段 cfl 5 から直行 | **DIVERGED step 20** (x/r_t 55〜70 の壁際の薄セルで T が上限 6000 K → ro NaN; `res_nan_20.h5`)。段階起動で run_0098 へ | 破棄予定 |
| `run_0096`・`run_0097` | (欠番: 0095 の発散で番号を送った) | — | — |
| `run_0098_ns_finemesh_pass_staged` | #11 ① 2 回目: run_0095 と同じ IC に段階起動 (soft 1 次 cfl 0.5 → mid 1 次 cfl 1 → 本段 2 次 cfl 5・relax 0.7) | soft・mid 完走、**本段 step 35 で DIVERGED** (x/r_t 79〜94 の出口側壁際で T → 6000 K・ω 120〜155 倍; `res_nan_35.h5`)。発散 2 回目 → codex 諮問 (`notes/reviews/briefs/2026-10-05-m6-finemesh-divergence.md`) | 破棄予定 |
| `run_0099`・`run_0100` | (欠番: #11a の切り分けで番号を送った。③④ は後続番号で回す) | — | — |
| `run_0101_ns_finemesh_diag_cfl5` | #11a 腕 A: run_0098 の本段開始場 (nozzle.h5) から 2 次・cfl 5、100 step (2 step ごと) | step 35 で再現発散。前駆: 軸の P 過渡 (16 %) → step 28 から x/r_t 88.3 の壁第 3 節点で T 指数成長 | 削除済み (2026-10-05、ディスク) |
| `run_0102_ns_finemesh_diag_cfl1` | #11a 腕 B: 同起点、2 次・cfl 1、1000 step | 停止なし、全残差の末尾 200 step が −0.03〜−0.09 桁 (短期合格) | 削除済み (2026-10-05、ディスク) |
| **`run_0103_ns_finemesh_pass_cfl1`** | #11 ①: 同起点、2 次・cfl 1、60000 step (累積 CFL を run_0092 にそろえる)、5000 step ごと | NOT CONVERGED (plateau 2.1〜2.8 桁)、machmax・pmax STEADY、δ_E(x_F) 0.733096 (5000 step ごとの変化 ≤ 0.003 %)、**壁解像 PASS (y1+>1 3.6 %)**、細分前後 (vs run_0094 0.725282) **+1.077 % (登録 ≤ 1 % FAIL)**; δ_E 時系列 STEADY、末尾 5 枚 0.733114 | active |
| `run_0104_ns_coarse_cfl5` | #11b 腕 A: run_0094 res_6000 を restart_field、粗格子・cfl 5・12000 step | NOT CONVERGED (plateau)、δ_E STEADY 0.725280 | ref |
| `run_0105_ns_coarse_cfl1` | #11b 腕 B: 同起点、粗格子・cfl 1・60000 step | NOT CONVERGED (plateau)、δ_E STEADY 0.725276 (CFL 5/1 差 +0.0006 %) → 細/粗 +1.08 % は格子差 | ref |
| `run_0106_ns_finemesh3_pass_cfl1` | #11c 第三水準: 細分格子の全方向 1/1.5 (ni 3000 × nj 145、第 1 セル 8.667e-6・スロート 3.0e-6、431856 セル、AR max 4252)、IC = run_0103 res_60000 を interp_field → 段階起動 → 2 次 cfl 1・60000 step。バイナリ `~/forge-wallfit-bin` (sha256 6b47811b…、plan §9) | δ_E STEADY 0.734629、run_0103 比 +0.21 % [+0.19, +0.22] % → **格子ゲート合格**、壁解像 PASS (0.0 %)、NOT CONVERGED (plateau) | ref |
| `run_0107_ns_finemesh_final` | #11 ③: ② の k_f 1.055734・r_t 76.6715 mm で最終 NS (細分格子、IC run_0103 → 段階起動 → cfl 1・60000 step) | 出口半径 0.7749995 m・δ_E/δ_C 1.0008・波 η0.1 0.0068 %・オーバーシュート η0.1 −0.0024 %・壁解像 PASS 3.6 %、**出口コア M 5.99831 (−0.028 %、登録 ±0.02 % FAIL)**、NOT CONVERGED (plateau) | active |
| `run_0109_ns_finemesh_final_ext` | #11 ③ の延長 1 回 (restart_field、cfl 1・60000 step) | 出口コア M 5.99833 (変わらず) → codex 諮問 | active |
| `run_0110_euler_pin_G0` | #11f Euler 格子 A/B 腕 G0 (旧較正格子 1100 × 65、生産経路の CFD ピン壁 = 旧 pincal と同一、新バイナリ、IC run_0086、cfl 2 × 12000) | 出口コア M 5.999996、STEADY | ref |
| `run_0111_euler_pin_G1` + `run_0112_euler_pin_G1_ext6k` | #11f Euler 腕 G1 (生産 NS 細分格子の格子パラメータ) + 延長 6000 | 出口コア M 5.999207 (STEADY) → G1 − G0 = −0.000789 [−0.000800, −0.000777] (判定 A) | ref |
| `run_0113_euler_pin_G1_recal` + `run_0114_euler_pin_G1_recal_ext6k` | #11f E2: Md_moc_offset +3.770e-4 で作り直した壁の Euler (G1) + 延長 6000 | 出口コア M 5.999998 (STEADY) → 合格。**run_0114 が新しい固定 Euler 参照** | ref |
| `run_0115_ns_recal_pass` | #11f NS ①: 新しい壁の細分 NS (k_f 1.055734、r_t 76.6715 mm、cfl 1・60000) | 完走 → ② k_f 1.054129・r_t 76.6539 mm (`c2pin_solve_recal.json`) | ref |
| `run_0116_ns_recal_final` + `run_0117_ns_recal_final_ext` | #11f NS ③: 最終 NS (k_f 1.054129・r_t 76.6539 mm、段階起動 → cfl 1・60000) + 延長 60000 | **出口コア M 5.998887 (STEADY) = 目標比 −0.019 %、許容 ±0.02 % の下限境界上 (格子・標本による不確かさ ±1e−4 と同程度; ユーザ決定 B で記録して設計は変えない)**・δ_E/δ_C 0.9998・オーバーシュート η0.1 0.0079 %・出口半径 0.7749995 m・壁解像 PASS 3.6 %; Mach 波 η0.1 0.0065 % は DRIFTING (ユーザ決定 A で未達のまま記録); NOT CONVERGED (plateau)。報告 `run_0117_ns_recal_final_ext/report/run_0117_ns_recal_final_ext_report.pptx` (VERDICT 入り、ローカル) | ref |
| `run_0118_ns_recal_final_cond` | #11 ④: run_0117 の凝縮 ON (convert_species_field conserve、cfl 1・18000) | 凝縮 4 量 STEADY: 開始 x 57.94・S_max 16.86・出口コア g 3.13e-4・出口コア M 5.9864 → 合格。報告 `run_0118_ns_recal_final_cond/report/run_0118_ns_recal_final_cond_report.pptx` (ローカル) | active |
| `run_0119_rerun_ctrl` 〜 `run_0123_rerun_fullpath` | rerun_conditions の検証 (plan tooling-rerun-conditions §6): 0119 無変更・0120 Euler Pt 0.8 (scale)・0121/0122 NS Pt 0.8 (scale あり/なし、stages none・cfl 5)・0123 Tt 1500 + H2O 0.10 の full 経路 | 0119 合格 (参照を再現)、0120 合格 (流量 0.8 倍に相対 7.5e−7)、0121/0122 発散 (出口壁際の角 / 入口)、0123 経路合格・量は DRIFTING | ref (監査: 各 run の `QS_VERDICT.txt`・`quantities_series.csv`・`stage_manifest.json`・`RERUN_CONDITIONS.json`、要約 `_band_ab/rerun_audit.txt`; 原データは主ツリー `/home/sano/work/forge/case/45.isobutane_m6_d155/`) |
| `run_0128_rerun_fullpath_ext`・`run_0135`〜`run_0138_rerun_fullpath_blk1..4`・`run_0134_rerun_euler_tt1500` | (iv)(iv″): 0123 の延長 (cfl 5・6000 × 5) と同条件 Euler 参照 | 出口 M は blk4 で STEADY (6.0327)、流量・δ_E は DRIFTING (ブロック間の増分は減衰中; 流量は簡易予想値 ≈ 17150 との差 約 0.1 %、整定余量は未確定 — plan tooling-rerun-conditions §5.1 #11)。**Euler 参照 run_0134 は 6000 step で DRIFTING** (δ_E の参照が動いていた) → ユーザ決定「1」で生産利用へ持ち越し | ref |
| `run_0129_rerun_pt08_scale_full` / `run_0130_rerun_pt08_noscale_full` | (ii′): Pt 0.8、scale あり/なし、stages full・本段 cfl 5 | 両方とも本段で発散 (step 468 出口壁際の角 / step 2 入口) | ref |
| `run_0131_rerun_pt08_scale_full_cfl1` + `run_0132_..._ext` / `run_0133_rerun_pt08_noscale_full_cfl1` + `run_0139_..._ext` | (ii′) 補足 A3 / (ii″) B3: Pt 0.8、full・本段 cfl 1・60000 + 延長 | A3 STEADY (出口 M 5.99212、δ_E 0.74945); B3 は規定時間 (60000 + 6000 step) 内に準定常に達せず、入口配管壁際の逆流域 (779 節点) が残った → Pt 変更は scale-ic pt を推奨 (確定) | ref |
| `run_0124_ns_axiscap020` / `run_0125_ns_axiscap0133` + 延長 `run_0126_ns_axiscap020_ext` / `run_0127_ns_axiscap0133_ext` | #11h 軸側上限 A/B: run_0117 の壁、nj 257・axis_cap_frac 0.02 / 0.0133、IC run_0117、段階起動 → cfl 1・60000 + 延長 60000 | 軸上オーバーシュート −0.007 / −0.009 % (生産格子 0.231 %)、η0.1 の Mach 波 0.0041 / 0.0040 %・オーバーシュート −0.022 / −0.023 %; A vs B は非軸のオーバーシュート・出口 M (共通標本) で差なし、波は準定常未達で保留、軸は保留。出口コア M 自格子平均 5.998888 (cap 0.02、帯内) / **5.998775 (cap 0.0133、下限外)** — 評価格子・標本で合否が入れ替わる (plan §9 の監査)。壁解像 PASS | ref |
| `run_0108_ns_finemesh_final_cond` | #11 ④: 0107 の凝縮 ON (cfl 1・18000 step) | 保留 (③ の出口 M FAIL) | active |
| `run_0140`〜`run_0142_euler_wallfit_pinG1_r{1,2,3}` | plan [tooling-nozzle-throat-monotone-r2](../../plans/accepted/tooling-nozzle-throat-monotone-r2.md) §6 E1 腕 A (現行壁、生産 Euler 格子 G1)、各 3 回の独立再実行。IC = run_0114 最終場を restart_field (ビット一致)、soft 3000 → 本段 2 次 cfl 2・relax 0.7・18000、バイナリ ~/forge-wallfit-bin (AWS)、commit bc08200f | 完走 2026-10-06 (forge exit 0、res 0〜18000)。本段区間の check_convergence は 3 本とも NOT CONVERGED (plateau、全保存量 0.9〜1.5 桁で横ばい; run_0114 と同じ既知の挙動)。最終場は NaN・Inf なし (ρ ≥ 0.0326、T ≥ 232 K、P ≥ 2208 Pa)。評価量の準定常と A/B の判定は腕 B とそろってから | active |
| `run_0143_euler_wallfit_monoG1_r1` / `run_0146_euler_icab_monoG1_nn` | 同 plan §6 E1 腕 B r1 (単調壁、IC = run_0114 最終場の保存量を検証付き番号写像) と予備 A/B の α (同じ格子、最近傍対応で直接転送)。段は腕 A と同じ。commit b80b14f5 | 完走 2026-10-06 (forge exit 0)。本段区間は 2 本とも NOT CONVERGED (plateau)。予備 A/B は**判別不能**: 11 量すべてが末尾 10 枚で STEADY でなく、U が閾値 Δq/10 の数十〜数百倍。step 6000〜18000 の β−α の平均差はどの量も 2 SE 以内 (例: オーバーシュート η0.1 +1.4e-3、SE 9.2e-4)。判定方法を諮問中 (ブリーフ `notes/reviews/briefs/2026-10-06-throat-mono-noise-limited.md`) | active |
| `run_0144_euler_wallfit_monoG1_r2` / `run_0145_euler_wallfit_monoG1_r3` | 同 plan §6 E1 腕 B r2・r3 (run_0143 と同じ準備から独立再実行)。予備 A/B の関門は B23_OVERRIDE で通した (理由はログ) | 完走 2026-10-06。本段区間 NOT CONVERGED (plateau)、最終場 NaN・Inf なし。**実務判定 §6 E′ (A 0140〜0142 vs B 0143〜0145)**: 6 量すべて (D + 2SE)/Δq ≤ 0.42 で許容幅未満 → 単調壁は候補形状。差を検出したのは \|P 傾き\| η0.1 だけ (0.1926 → 0.2008 %pt、Δq の 37 %)。出口コア M は両腕 5.99999。結果 `_band_ab/throat_mono_practical_eval.json` (AWS) | active |
| `run_0147_ns_mono_final` + 延長 `run_0149_ns_mono_final_ext` / `run_0148_ns_mono_final_cond` | plan tooling-nozzle-throat-monotone-r2 §6 N・K: 単調壁の dry NS (問題 `problem_d155_ns_finemesh_recal_final_mono.yaml`、IC = run_0117 res_60000 の検証付き番号写像、段階起動なし、cfl 1・60000 + 延長 20000) と凝縮 NS (IC = run_0147 res_60000 を convert_species_field、18000)。`run_mono_ns_chain.sh`・`run_mono_ns_ext.sh` | 完走 2026-10-06。dry: 窓 60000〜80000 で全ゲート合格・4 量 STEADY (出口コア M 5.998871、旧壁 run_0117 5.998887)。凝縮: 4 量 STEADY (開始 57.965、旧壁 57.938)。NaN なし、残差 plateau (RISING なし)。報告 `run_0149_…_report.pptx`・`run_0148_…_report.pptx`、要約 `_band_ab/throat_mono_compare/throat_mono_summary_report.pptx` (ローカル) | active |
| `run_0150_euler_wallfit_mocG1_r1` / `run_0151_euler_wallfit_mocG1_r2` / `run_0152_euler_wallfit_mocG1_r3` | plan discretization-moc-axis-limit-and-corrector §6 V5 腕 M: 単調壁 + MOC の軸上の解析極限・収束する修正子 (問題 `problem_d155_euler_pin_G1_recal_mono_moc.yaml` = 腕 B の問題 + `moc_axis_limit: analytic`・`moc_corrector: converge`)。生産 Euler 格子 G1、IC = run_0114 res_6000 の検証付き番号写像 (設計壁の変化最大 6.2 µm のため上限 6.21 µm、検査 C1〜C5 は `_band_ab/moc_v5_ic_inspection.json`)、soft 3000 → 本段 18000。比較相手は腕 B (run_0143〜0145)。`run_moc_v5_euler.sh` | 投入 2026-10-07 (AWS)。完走 (RUN_RC=0、soft 3000 + 本段 18000 step)、NaN・Inf なし (全 `residual_history*.csv` の `rms_*` 全列と `res_0`〜`res_18000` の `VALUE/*`)、最終場 ro・P・T > 0。本段区間 (main、`stage_manifest`) の check_convergence は 3 本とも `NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)` (全列 STALLED、低下 0.9〜1.5 dec)。成果物 (AWS `~/forge-wallfit/case/45.isobutane_m6_d155/<run>/`): `res_18000.h5`・`residual_history.png`・`CONVERGENCE_VERDICT_segment.txt`。V5 の判定は評価器の修正後。**V5 判定 (2026-10-07): 保留 (前提不成立)** — 7 本すべてが窓 6000〜18000 の準定常の前提を満たさない (`_band_ab/moc_v5_euler_eval.json`)。参考値は plan §9 | active |
| `run_0153_euler_icdep_mocG1_isen` | 同 §6 V5 の IC 依存の確認: 腕 M と同じ壁・格子で、等エントロピー IC から段階起動 (soft 3000 → 本段 18000) | 投入 2026-10-07 (AWS)。完走 (RUN_RC=0、soft 3000 + 本段 18000 step)、NaN・Inf なし、最終場 ro・P・T > 0。本段区間 (main、`stage_manifest`) の check_convergence は `NOT CONVERGED (still converging — run more steps)` (全列 falling、低下 2.6〜2.7 dec)。成果物 (AWS `~/forge-wallfit/case/45.isobutane_m6_d155/run_0153_euler_icdep_mocG1_isen/`): `res_18000.h5`・`residual_history.png`・`CONVERGENCE_VERDICT_segment.txt`。V5 の判定は評価器の修正後。**V5 判定 (2026-10-07): 保留 (前提不成立)** — 7 本すべてが窓 6000〜18000 の準定常の前提を満たさない (`_band_ab/moc_v5_euler_eval.json`)。参考値は plan §9 | active |
| `run_0154〜0156_euler_wallfit_monoG1_r{1,2,3}_ext36k` / `run_0157〜0159_euler_wallfit_mocG1_r{1,2,3}_ext36k` / `run_0160_euler_icdep_mocG1_isen_ext36k` | plan discretization-moc-axis-limit-and-corrector §6 V5b (延長の診断): 親 run_0143〜0145・0150〜0153 の `res_18000.h5` から同一メッシュ restart (ビット一致) で追加 36000 step (通算 54000)、同じ設定・同じバイナリ、soft 段なし。窓 A 24000〜36000 / B 42000〜54000 で減衰か持続振動かを判別。`run_moc_v5b_ext.sh` | 完走 2026-10-07 (7 本 RUN_RC 0、`EARLY_STOP.txt` なし)、NaN・Inf なし (全 `residual_history*.csv` の `rms_*` 全列と `res_0`〜`res_36000` の `VALUE/*`・境界パッチ出力)、最終場 ro・P・T > 0・P ≤ P0。本段区間 (main、`stage_manifest`、108000 行 = 36000 step × 3 行) の check_convergence は 7 本とも `NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)` (全列 STALLED、低下 0.1〜0.2 dec)。成果物 (AWS `~/forge-wallfit/case/45.isobutane_m6_d155/<run>/`、1 本 957 MB): `res_36000.h5`・`residual_history.png`・`CONVERGENCE_VERDICT_segment.txt`。**V5b 判定: 保留** — 7 本すべてが減衰側ではない (2 本は exit_M_dev が持続振動側、他は判別不能)。窓 A→B でも量が動き続ける。結果 `_band_ab/moc_v5b_ext_eval.json`、参考値は plan §9 | active |
| `run_0161_euler_t0cluster_g1` / `run_0162_euler_t0cluster_u5em3` | plan verification-case45-euler-total-enthalpy §6 E2: Euler G1 (2000 × 97) で半径方向の配点だけを変える A/B。A = G1 のまま (`wall_first_frac` 1.3e-5・`wall_first_frac_throat` 4.5e-6)、B = 両方 0.005。壁は MOC の V5 の腕 M と同じ、等エントロピー IC、soft 3000 → 本段 54000。スロート付近の全温の超過が配点に依存するかを見る。`run_euler_t0_e2.sh` | 完走 2026-10-07 (RUN_RC 0、NaN なし、本段 NOT CONVERGED stalled/plateau)。**E2 判定: 判別不能** — B は窓 42000〜54000 の全時点で全領域の偏差 ≤ 0.103 K・時間の幅 ≤ 0.1 K、A はスロートで > 182 K だが時間変動が大きく整定の前提を満たさない。結果 `_band_ab/euler_t0_e2_eval.json`、plan §9 | active |
| `run_0163_euler_e4_recal_d0` | plan verification-case45-euler-total-enthalpy §6 E4 (出口較正のやり直し、段 1): 単調壁・legacy MOC・`mesh_euler` (2000 × 97・全域 0.005)・`Md_moc_offset` +3.770e-4 (δ₀)、等エントロピー IC、soft 3000 → 本段 54000。出口コア M を G1 の η の列で評価。`run_e4_recal.sh d0` | 完走 2026-10-07 (RUN_RC 0)。**段 1 は判別不能** (全温の時間の幅が 27 列中 1 列で 0.1065 K > 0.1 K。出口 M は両窓 STEADY・M_common の平均 6.000370)。結果 `_band_ab/e4_recal_eval.json` | active |
| `run_0164_euler_e4_recal_d1` | 同 plan §6 E4V (候補 δcand = 6.8825162455159465e-6 の独立の検証、1 回だけ): run_0163 と同じ条件で `Md_moc_offset` だけを δcand に。`run_e4_recal.sh d1 6.8825162455159465e-06` | 完走 2026-10-07 (RUN_RC 0)。**E4V 合格**: M_common の平均 5.99999968・max\|M − 6\| 3.1e-6、全温 ≤ 0.134 K → 新しい較正値 6.8825162455159465e-06・新しい Euler 参照 (`_band_ab/e4_recal_eval.json`) | ref |
| `run_0171_euler_v5d_B_r1` / `run_0172_euler_v5d_B_r2` (+ run_0164 を腕 B の 3 本目に共用) / `run_0174〜0176_euler_v5d_M_r{1,2,3}` | plan discretization-moc-axis-limit-and-corrector §6 V5d: 新しい Euler の格子 (`mesh_euler` 全域 0.005) で、腕 B = 単調壁・legacy MOC、腕 M = 単調壁・analytic + converge、両腕とも較正値 6.8825162455159465e-06 (E4V)、等エントロピー IC、soft 3000 → 本段 54000。`run_moc_v5d.sh --share-e4v` | 完走 2026-10-07 (RUN_RC 0)。**V5d 判定: 保留 (前提不成立)** — 窓の中の時間変動の条件 (幅 ≤ Δq/10・全温 0.1 K) が一部の run で不成立。参考値: (D + 2SE) は全量で許容幅の内、ただし出口 M の誤差は 4.1e-6 → 4.5e-5 と増える向き (2026-10-07 訂正: 初稿の「悪化の向きの量なし」は誤り)、オーバーシュート 0.0103 → 0.0042・P 傾き 0.057 → 0.0075、腕 M の出口 M は 3 本とも \|M − 6\| ≤ 5.3e-5。結果 `_band_ab/moc_v5d_eval.json` | active |
| `run_0165_ns_n012_N0` / `run_0166_ns_n012_N1` / `run_0167_ns_n012_N2` (+ 凝縮 `run_0168〜0170_ns_n012_N{0,1,2}_cond`) | plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4 (NS の 3 条件): N0 = 今の生産壁 (ramp・legacy MOC)、N1 = 上流を多項式に、N2 = さらに MOC の analytic + converge。3 条件とも較正値 6.88e-6・r_t 0.07666537 m・k_f 1.054129、IC は run_0149 から cross-mesh + 段階起動、本段 80000。`run_ns_n012.sh main` | 完走 2026-10-07 (RUN_RC 0、NaN なし、本段 plateau・RISING なし)。**3 条件とも dry 未達**: 出口コア M 5.99851〜5.99854 (STEADY、旧ゲートの下限 5.9988 を約 3e-4 下回る。**2026-10-07 ユーザ決定でゲートを 6.000 ± 0.05 % [5.997〜6.003] に変更、新しい帯では入る**)、オーバーシュートの準定常 DRIFTING/OSCILLATING (値は 0.002 % 前後)。他のゲート (波・δ_E/δ_C 1.0017・出口半径・壁解像 3.6 %) は合格。凝縮は保留 (作っていない)。N1−N0 は全量が時間の幅以下、N2−N1 はオーバーシュート −0.0037 %。結果 `_band_ab/ns_n012_eval.json` | active |
| `run_0177_ns_n012_N0_ext` / `run_0178_ns_n012_N1_ext` / `run_0179_ns_n012_N2_ext` | 上の 3 本の登録どおりの延長 1 回 (各 res_80000 から `restart_field.py`、20000 step・5000 ごと、設定は同一)。諮問 `notes/reviews/2026-10-07-ns-n012-exitM-deficit-diagnose.md`。判定窓は通算 80000〜100000 | 完走 2026-10-07 (RUN_RC 0、NaN なし、plateau・RISING なし)。新ゲートで出口コア M 5.9985 合格・波 STEADY、残る未達だったオーバーシュートの準定常 (DRIFTING、値 0.0019/0.0021/−0.0015 %、窓の動き約 0.0005 %) は、2026-10-07 ユーザ決定 (窓の絶対の幅 ≤ 0.0035 % で可) で合格 → **3 条件とも dry 全ゲート合格**。結果 `_band_ab/ns_n012_eval.json` | active |
| `run_0168_ns_n012_N0_cond` / `run_0169_ns_n012_N1_cond` / `run_0170_ns_n012_N2_cond` | 凝縮 (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4、monotone plan §6 K の条件)。IC は各延長 run (0177/0178/0179) の res_20000 を `convert_species_field`、18000 step。`run_ns_n012.sh cond` | 完走 2026-10-07 (RUN_RC 0、NaN なし)。**N1・N2 合格、N0 未達** (凝縮のモーメントの残差 `rms_rog_0`・`rms_roQ2_0` が RISING、末尾でゆっくり増加。4 量は STEADY)。出口コア M 5.98647/5.98647/5.98662、出口 g 3.01e-4/3.01e-4/2.98e-4。結果 `_band_ab/ns_n012_eval.json` | active |
| `run_0180_ns_n012_N0_cond_ext` | N0 の凝縮 (run_0168) の延長 1 回 (登録、諮問 `notes/reviews/2026-10-07-u4-v5prime-result-interpretation-diagnose.md`)。res_18000 から 13 量をビット一致で継いで 20000 step。`ns_n012_cond_ext.py`・`_band_ab/cond_ext_N0.sh` | 完走 2026-10-07 (RUN_RC 0、NaN なし)。**K 合格**: 通算 34000〜38000 で 4 量 STEADY、連結した残差は plateau・RISING なし。結果 `_band_ab/ns_n012_cond_ext_run_0180_ns_n012_N0_cond_ext.json` | active |
| `run_0181_ns_coldmesh_ad` / `run_0182_ns_coldmesh_tw300` | **冷却壁の NS の対** (plan tooling-nozzle-isothermal-wall-chain §5.1 #14、§6 V-c45 事前登録): 生産の壁 (`_band_ab/prod_confirm/prep`) を冷却壁用の格子 (`problem_d155_ns_prod_coldmesh{,_tw300}.yaml`: ni 4719 × nj 121、第一層と x 密度の表、近壁は壁法線、msh 17 桁) に載せ、**FP64 のビルド** (`~/forge-wallfit-bin-fp64`: 生産のソース + typedef double + 座標の `stod`、forge 65be5e28…・変換器 ac88861f…) で断熱 / 300 K を回す。IC = run_0179 の最終場の cross-mesh、段階起動 full → 本段 100000 step・5000 ごと (2 次・cfl 1・relax 0.7)。`run_cold_pair.sh` | 投入 2026-10-08。準備の前提は全て成立 (品質 PASS: AR 最大 4573・skew 0.457、第一層の誤差 0、高 AR のスキューセル 0、壁の差 0、壁距離一致)。序盤 NaN なし (各 1200 step)。判定は V-c45。ディスク確保のため、判定済みの run_0165・0166・0168・0169・0180 の中間の全場スナップショット (res_0 と最終以外の `res_<n>.h5`) を AWS で削除 (壁・出口ダンプ・残差・判定ファイル・時系列 CSV は残す) | 実行中 |
| `run_0183_ns_coldmesh_tw300_ext` | run_0182 (300 K) の延長 100000 step (V-c45「延長の決め方」、res_100000 から restart_field ビット一致、FP64) | 完走・NaN なし。V-c45 の 2 回目も判定不能 (x = 40 の δ_E と Q_w が未定常、rms_roOmega RISING、縮流部の壁解像 23 %)。R_NS 0.726 (x = 40)〜0.793 (x = 94) | active |
| `run_0184_ns_coldmesh_cpg_ad_plain` / `run_0185_ns_coldmesh_cpg_tw300_plain` / `run_0186_ns_coldmesh_cpg_ad_dilat2` / `run_0187_ns_coldmesh_cpg_tw300_dilat2` | **SU2 との照合の forge 側** (plan tooling-nozzle-isothermal-wall-chain §5.1 #26): 壁・格子は run_0181/0182 と同じ、CFD だけ CPG (γ 1.27354・cp 1360・Sutherland 1.716e-5/273/111・定 Pr 0.72)。素の SST (dilat 0・KL 0) / 生産の SST (dilat 2・KL 1) × 断熱 / 300 K。初期値は TP の解 (run_0181・run_0183 の res_100000) の ρ・U・P・k・ω を CPG に組み直し、段階起動 full → 本段 100000 step・20000 ごと。FP64。`run_cold_cpg.sh` | 投入 2026-10-08 | 実行中 |
| `run_0188_su2_coldmesh_cpg_ad` / `run_0189_su2_coldmesh_cpg_tw300` | **SU2 v8.5 (手元の PC、6 スレッド × 2)**: 同じ格子 (`_band_ab/cold_pair/su2/nozzle.su2`、`cold_su2_mesh.py` で msh から 17 桁のまま変換・座標の完全一致を確認)、CPG・SST V2003m・ROE + MUSCL・`MARKER_HEATFLUX 0` / `MARKER_ISOTHERMAL 300`、入口乱流 TI 0.077・粘性比 12.1、出口 2237 Pa、80000 反復 (10000 ごとに restart の CSV、CFL 5 適応 (0.5, 1.5, 2, 50))。**初期値は forge の CPG の run_0184 / run_0185 の段 S1 (soft) 直後の `nozzle.h5` の場** (`cold_su2_ic.py`、`_band_ab/cold_pair/su2/ic_run_018{4,5}_*.npz`。Energy = roe + ρk、k・ω は原始変数。ρ は元の TP の場 run_0181 res_100000 とも forge の本段の開始場 res_0 とも一致しない — 2026-10-08 確認)。自由流からの起動 (1.4〜2 s/反復) は遅すぎたので止めて、この初期値から 19:40 JST に再投入した | 投入 2026-10-08。1.1〜1.5 s/反復 (i5-12400F の 6 物理コアに 6 スレッド × 2 本)、80000 反復で 25〜33 時間の見込み。CFL は 2〜50 を行き来する。表面の CSV は restart が compact (既定) なので解の量だけになり C_f・熱流束を含まない (SU2 `COutput.cpp` の SURFACE_CSV で確認) → 壁の量は両コードとも体積の場から `cold_xcheck.py` で同じ式で作り、SU2 の履歴の `HF` (積分した熱流束) で式を確かめる | 実行中 |
| (注) | 2026-10-05: AWS の空き不足のため run_0053〜0105 の中間 `res_<n>.h5` を削除 (res_0・最終場・`delta_E_series.csv` は残る) | — | — |
| `run_0008_ns_trim_cond` | 凝縮 ON restart (`problem_d155_trim_ns_cond.yaml`: Kw+HK condModel1+Kantrowitz, 蒸発 ON, IC=run_0007, 12000 step) | 完走・NaN 0・**STEADY** (4k/8k/12k で M_exit 差 5e-4)。軸 onset x≈69 r_t、出口 g 0.20 % (H₂O の 2 %)、**出口軸 M 5.9273 (−1.2 %)**・試験区間に M 低下勾配 (x60→96 で 6.00→5.93)。dry の軸は x≈24 r_t (M5.5) で飽和線越え S≈14 (`axis_values.csv` の Tsat_post) | active (**凝縮評価の正本**) |

## SU2 クロスチェック (境界層厚さ, 2026-09-01)

δ\* 抽出の軸対称化 (円環重み) と合わせて、**同一メッシュ・同一 BC・同一ガス (CPG γ\*=1.27354)** で forge node SST と SU2 v8.5 axisym SST の境界層厚さを比較する ([procedures/su2-cross-check.md](../../procedures/su2-cross-check.md) 準拠)。

| run | 内容 | 結果 | 状態 |
| --- | --- | --- | --- |
| `run_0009_cpg_euler` | CPG 版 Euler (`problem_d155_cpg_euler.yaml`) | 完走・NaN 0。CPG は A/A*=159 (semi-perfect 88.2) でノズルが伸びる (出口壁半径 1.03 m) — 比較チェーン専用で製品形状ではない | active (比較チェーン) |
| `run_0010_cpg_ns_coarse` | CPG NS 中継 (y+~50) | 完走・NaN 0 | active (中継) |
| `run_0011_cpg_ns` | **CPG NS 本計算** (wall_first 6.5e-5=y+~2 [4.5e-5 は CPG 長尺で AR1403 FAIL→緩和], 相関 δ\* 物理壁, 24000 step) — 比較の forge 側 | 完走・NaN 0・品質 PASS。残差 2 桁 warm 床 (既知パターン) | active (**比較 forge 側**) |
| `run_0016_chamber_cpg` | **チャンバー A/B** (Codex レビュー準拠): 完全形状 (入口管+収縮+ノズル) + 出口後方 4D_e×外径 3D_e チャンバー+フランジ壁の 4 ブロック transfinite (`gen_chamber_mesh.py` → `mesh_chamber/`)。品質 PASS (AR296)。静止雰囲気 1389 Pa IC + 段階起動 18000 step で plume 確立 (ジェット軸 M5.97, 外周静止)。**probe (cfl4 OK/8@181/16@21) は素メッシュと同一挙動 = 出口境界は無罪**。**line-implicit 併用の 4 象限も全て同一** (coarse×line 8@216/16@15, chamber×line 8@194/16@15; チャンバーのライン被覆 81.9%・噴流域 point fallback) — 壁法線×出口境界の直積消去で streamwise 内部モード律速が確定 | active (**A/B の正本**) |
| `run_0012_su2_sst` | **SU2 v8.5 RANS-SST** (同一メッシュ msh→su2 変換, axisym, ROE+MUSCL, SST V2003m, Sutherland/Pr0.72/Prt0.9, 40k+20k iter) | 収束: rms[Rho] −6.5・rms[RhoE] −2.3・出口 massflux/Mach ドリフト 0.003 %/5k (手順書合格実績以上)。`history.csv` は継続分のみ | active (**比較 SU2 側**) |
| `run_0013_cpg_ns_plainsst` | forge 素 SST 化 A/B (`dilatationCorrection: 0, katoLaunder: 0` — prepare 後に solverConfig 直接編集。YAML への独自キー追記は黙って無視されるため不可), IC=run_0011, 12000 step | 完走・NaN 0。**素 SST にすると SU2 と一致: δ99 ≤3 %・θ ≤0.4 %・δ* ≤1.3 %** | active (**帰属 A/B の正本**) |
| `run_0014_cpg_ns_dilat0_kl1` | 帰属分離 (dilatation OFF / KL ON) | run_0013 と ≤0.7 % 差 → **Kato-Launder は無関係、差は dilatationCorrection 単独** | active (帰属分離) |

### SU2 クロスチェックの結論 (2026-09-01)

- **forge 生産設定 SST は SU2 比で境界層が薄い**: δ99 −17〜21 %・θ −16 %・δ* −5 % (4 ステーション一貫)。
- **帰属は `dilatationCorrection: 2` (圧縮性生産項の正確形: deviatoric trace 除去 + 等方項 −⅔ρk∇·u) 単独**。素 SST 化した forge は SU2 と δ99 ≤3 %/θ ≤0.4 %/δ* ≤1.3 % で一致 = **ソルバ・離散化は無罪、意図した乱流モデル差**。Kato-Launder の寄与 ≤0.7 %。
- ノズルは全域 ∇·u>0 (強膨張) のため等方項が k のシンクとして常時働く。急膨張の乱流抑制は物理的に実在する効果で forge 形は正当だが、**SST のモデル形式差として境界層厚に ±16 % 級の不確かさ**があると認識すること。δ* への影響は ±5 % (出口 M ±0.1〜0.15 %) に留まり、Md トリムは forge 自身の NS で較正しているため製品性能の結論は不変。
- 比較図・数値: `compare_bl_su2.{png,json}` (3 者重ね描き)、40k 時点の記録 `compare_bl_su2_iter40k.json`。

## 陰解法 cfl_pseudo 上限スイープ (2026-09-02)

`run_0015_cflsweep/` (probe 22 本, `sweep_cfl_implicit.py` + `summary.json`。中間 res は削除済み・各 probe の config/residual/NaN ダンプは保持)。warm 場 (run_0007 TP / run_0011 CPG) restart 2000 step + coarse IC 収束レース 4000 step。

- **素の上限 (implicitRelax=1): TP は cfl 2 (3 で発散)、CPG は 3 (4 で発散)**。既定挙動としては単桁前半で頭打ち。
- **原因は低マッハ調整ではない** (仮説棄却): ① 現行 config は `lowMachPrecond: 0` で低マッハ処理は不活性、② 発散箇所は亜音速チェンバでなく**下流端の超音速壁 BL 内** (壁から ~0.003 r_t の対数層)、③ `lowMachPrecond: 2` を入れると逆に悪化 (step 97→13)。
- **終端症状 = 陰的更新の P アンダーシュート → EOS 床洗浄 → NaN** (ω 爆発は増幅表示)。**帰属の最終確定 (2026-09-02, チャンバー A/B 完了)**: 発散種の位置はメッシュ依存 (fine=壁∩出口コーナー / coarse=出口手前コア / 高cfl=スロート軸上) で、①背圧値 100/10 Pa・②逆流 Tt 1500K・③出口近傍局所 dt キャップ・④**チャンバー付加 (出口境界を 8 m 後方へ, run_0016)** の全てに**不感** (発散 cfl・step ほぼ不変: chamber cfl4 OK/8@181/16@21 vs 素 cfl4 OK/8@215/16@17)。⑤ **乱流状態更新の凍結 (`FORGE_FREEZE_TURB=1`) でも発散 step が 1 step 単位で一致** (cpg@98/28/18, tp@53/24)。**確定したのは**: k-ω segregated 更新は発散に不要 / ω 先行 NaN は症状 / 流れブロックだけで発散する。**未消去**: μt は凍結 roK,roOmega から毎 step 再計算される (dependentVariables/turbulent_viscosity) ため SST closure・μt を含む粘性剛性の寄与は残る。⑥ **laminar NS (`model: "none"`, run_0017 で定常化後 probe)**: cfl4 OK / 8@127 / 16@16 — SST と同じ上限。⑦ **Euler (run_0009 base, 粘性・乱流・BL なし)**: cfl4 OK / 8@150 / 16@55 / 32@19 — **やはり同じ 4→8 境界** (種は x 4–12 r_t の近軸内部 = 最急膨張域)。⇒ 判定木の「Euler まで同じ」分岐に着地: **律速は対流の streamwise defect-correction 反復不安定** (SST closure・μt・粘性剛性は全て消去済み)。根治は streamwise ADI / LU-SGS 流下順序が妥当。
- **切り分け**: `nStepInner` 10/20 無効 (DPLUR 内部収束でない)・`implicitRelaxSST` 0.5 無効 (SST 方程式でない) — **`implicitRelax` (流れ方程式の緩和) だけが効く**。
- **緩和込みの上限**: relax0.7 → cfl 8 (12 で発散)、relax0.5 → cfl 16+。積 cfl×relax ≈ 6-8 で飽和。
- **収束レース**: 実効速度も cfl×relax でスケールし **cfl8+relax0.7 が最速** (4000 step で rms_roUx 現行 cfl1 比 1/3、rms_roOmega 1/4)。cfl16+relax0.5 は頭打ち。図 `run_0015_cflsweep/race_residuals.png`。
- **生産推奨**: CPG `cfl 8 + implicitRelax 0.7` / TP `cfl 6 + implicitRelax 0.7`。現行 cfl1 比で 2〜4 倍速。**TP cfl6+r0.7 は実段階起動系列で検証済み** (`run_0018_ns_trim_cfl6_r07`: run_0007 と metrics 5 桁一致・machmax/pmax STEADY・NaN 0。素の TP cfl3 は発散するため **implicitRelax 0.7 とのペアで固定**すること)。停止条件 (十分速く安定に答えが出る) を満たしたため **ADI/LU-SGS の実装は棚上げ** — 速度が再びボトルネックになったら着手 (GPU 可否は検討済み: 第2ライン族 ADI は v1 機構流用で可、古典 LU-SGS は 2D 波面並列が細く不適 = DPLUR の存在理由)。
- **line-implicit 追試 (2026-09-02, probes7–10, [plan](../../plans/accepted/time_integration-line-implicit.md))**: 壁法線ライン block-Thomas (`lineImplicit: 1`) を実装して検証。ライン構築は完璧 (1250 本・被覆 100 %) だが **cfl 上限は不変 (積 ≈6–8)** — cfl8-line の発散は出口域全断面で、**律速は壁法線でなく streamwise lag** と判明 (壁法線剛性は粘性対角が既に吸収)。効果は同一設定の収束 −14 %/step のみ。実装過程で lu5 ピボット罠 (LASWP 先行必須)・device printf 引数上限罠を踏んで記録。
- **追試 (2026-09-02, probes5/6)**: commit 時の正値性ガード (`updateGuardAlpha`, 新実装 opt-in) は**上限を上げない** (発散遅延のみ。α=0.5 は毎 step 半減を許すため床への歩行が続く)。**lowMachPrecond=2+ガード併用も不成立** — P が床に触れないまま SST チャネル (ω 2.8e22) で爆発。⇒ 真因は「defect-correction 反復 (近似 Jacobian + SST 隔離更新) の不安定」で、床 NaN は終端症状。`implicitRelax` が効くのは反復スペクトル半径を縮めるから。詳細 [plan](../../plans/accepted/time_integration-update-positivity-guard.md)

## 結論 (2026-08-31)

- **最短全長: スロート→物理出口 7.337 m** (x_F 96.00 r_t, r_t 76.42 mm)、入口配管端→出口 8.292 m。うち E→F 一様化区間 ≈4.36 m は物理の床。
- **出口軸 M (dry): 6.0008 (+0.013 %)** — Md 6.0144 トリム (NS 固定点の −0.24 % 欠損を打ち消し) で達成。
- **凝縮リスク**: Tt 1600 K では dry 不成立 (S≈14)。凝縮 ON で出口 M 5.927・軸勾配残り。凝縮フリーは Tt ≳ 1820 K。
- 結果ページ: https://claude.ai/code/artifact/f1cfbdf8-2415-4bfc-ae2d-156793569bd0

注: `run_0006_ns_trim` (旧設計 run_0004 の δ\* CSV を新設計に流用するショートカット) は物理壁フィルタ不合格 (スロート下流非単調) で prepare 段階に失敗し**削除済み**。トリム版も正規フロー (v1→抽出→v3) で回す。

2026-10-07: run_0150〜0163 (V5・V5b・E2・E4 段 1、判定済み) の中間の全場スナップショット `res_<n>.h5` を AWS で削除した (res_0・最終の res・境界の出力・時系列・判定の記録は残した)。
```

## 参考: `notes/reviews/2026-10-08-cold-su2-crosscheck-diagnose.md`

```
# codex 諮問 (diagnose): cold-su2-crosscheck

- **brief**: [`notes/reviews/briefs/2026-10-08-cold-su2-crosscheck.md`](../../notes/reviews/briefs/2026-10-08-cold-su2-crosscheck.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-08
- **commit**: `5eb46ef0` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.3 min, rc=0
- **結論**: **最初に、TP→CPG 変換後の実効設定で熱伝導則だけを比較する CFD 0 step の A/B を行い、定 Prandtl 数への復元を確認してください。**
- **extra**: `case/45.isobutane_m6_d155/cold_pair.py`, `case/45.isobutane_m6_d155/run_0012_su2_sst/sst.cfg`, `procedures/su2-cross-check.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**CPG での SU2 照合は採用しますが、現状の事前登録は修正が必要です。特に、記載どおりの TP→CPG 書換えでは熱伝導モデルが一致しません。**

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **CPG への設定変換：要修正** | TP 設定の生成時に `thermCondMethod` を削除しています（[runner_axismach.py:152](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:152)）。省略時は `0`（[solverConfig.cpp:957](/home/sano/work/forge-integ-1005/solver_density_cuda/input/solverConfig.cpp:957)）なので、`viscMethod: 1` に戻すだけでは熱伝導率は定数のままです。**`physProp.thermCondMethod: 1`・`prandtlLam: 0.72` を明示し、TP 専用の `transport` も外す**必要があります。Sutherland の粘性を合わせるだけでは足りません。 |
| **Major** | **「素 SST＝SU2 と同じ形」：要再検証** | [plan §4.3:98](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:98) の `sstOmegaProdFromPk: 0`・`sstSigmaBlend: 0` という対応表は採用できません。SU2 v8.5 の V2003 は制限後の P_k から ω 生産を作り、拡散係数もブレンドします（[SU2 公式ソース](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/include/numerics/turbulent/turb_sources.hpp#L880)）。対応する forge の設定はこの部分では **両方 `1`** です（[ransSource_d.cu:227](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransSource_d.cu:227)、[ransTransport_d.cu:110](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransTransport_d.cu:110)）。入口 k/ω、壁条件、エネルギー中の k の扱いまで実効条件表を作ってください。旧 SU2 cfg の出口圧も `1389 Pa`、TP 準備側は `2237 Pa` なので、そのままの流用は不可です。 |
| **Major** | **Q1：壁・格子を維持した CPG 比較：範囲を限定して採用** | 同じ形状の CPG 問題を両コードで解く比較として妥当です。ただし **TP の冷却応答全体の正しさは証明できません**。定比熱化により h(T)、温度反転、組成・種輸送、混合気の μ・λ、およびそれらと冷却・乱流との相互作用が検証対象から外れます（[solver-settings.md:233](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:233)）。合格時の結論は「この CPG 条件・格子・比較量でコード間差が許容内」。不合格時も forge の誤りと即断せず、条件・モデル・離散化の不一致として扱ってください。 |
| **Major** | **Q2：共通 y_b の δ_loc・θ_r：条件付き採用** | 共通区間の比較量としては妥当ですが、**境界層全体の厚さを代表する保証はありません**。保存された `theta_*.npz` では、冷却 `run_0183` が選んだ y_b は断熱 `run_0181` の **1.235／1.211／1.223 倍**（x=40／70／94）。`theta_diag` は y_b の点値を ρ_e・u_e に使うため、帯の変更は積分範囲と正規化を同時に変えます（[cold_pair.py:575](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair.py:575)）。これは切断誤差の証明ではありませんが、断熱 TP の帯を無条件に外縁とは呼べません。固定した観測量として定義し、CPG 各腕でも外縁の位置・帯への感度を確認してください。δ_E の検証を完了した扱いにはできません。 |
| **Major** | **Q3：R 2%・厚さ/C_f 3%：目標として採用、実績による裏付けは却下** | 旧比較の θ は `(ρu/ρ_eu_e)(1−ρu/ρ_eu_e)` を積分した代理量です（[compare_bl_su2.py:82](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/compare_bl_su2.py:82)）。今回の θ_r は `ρu(u_e−u)` であり、**「前回 θ≤0.4%」をその精度根拠にできません**。3% は暫定の照合目標とし、時間・抽出・写像の不確かさを含めて判定してください。各腕3%から比2%は保証されないので、両条件を残します。冷却の解き方を検証するなら q_w を記録だけに落とさず、既存の [§4.4](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:112) の **5%目標**と熱収支確認も維持すべきです。 |
| **Major** | **Q4：TP 場からの初期化：採用／「収束済み」と現ゲート：却下** | 初期推定に未収束場を使うこと自体は可能ですが、`run_0181/0183` は収束済みではありません。下記の判定が根拠です。同一節点対応で保存量を CPG の EOS に整合させ、TP の `roe` はコピーしないこと。ρ・P を維持すると T が変わるため、壁温ピン後の状態も確認が必要です。初期化経路が違っても比較はできますが、**NaN/RISING なし＋5枚 STEADY だけでは高残差停滞を排除できません**。これは [§5.1 #22](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:204) で既に却下されています。残差上限・保存収支・壁解像を登録し、θ_r・C_f・断熱壁温も含む全比較量の定常性を確認してください。延長は各コードのゲートで決め、同じ反復数には揃えません。 |
| **Major** | **Q5：素 SST と生産 SST の両方：採用／差を dilatation 単独に帰属：却下** | 素 SST の対はコード間比較、生産 SST の対は生産条件への接続に必要です。ただし `(dilat, KL)=(0,0)→(2,1)` は2因子変更です。[ransSource_d.cu:127](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransSource_d.cu:127) と [同:134](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransSource_d.cu:134) で別々に作用します。現行4本から言えるのは **「生産 SST 設定一式の感度」**までです。補正を物理的に必要と考えることと、比較用に同じモデルへ揃えることは両立します。 |

観測事実の確認として、保存 NPZ に `check_quasisteady.classify` を実行しました。判定窓は各 run の **80000〜100000、5枚全部、drift・osc 各0.1%**です。

- `case/45.isobutane_m6_d155/run_0181_ns_coldmesh_ad/`：δ_E の x=40／70／94 は **`STEADY`**。保存済み収束判定は **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。
- `case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/`：δ_E は x=40 が **`TRANSIENT-UNSETTLED`**、70／94 が `STEADY`。Q_w は **`DRIFTING`**、窓内変化約2.5%。保存済み収束判定には **`rms_roOmega … RISING (divergent)`** があります。
- 保存された総合判定も **`判定不能 (ゲート不成立)`**です。「冷却効果約25%」は現在の抽出値であり、確定した定常応答ではありません。

根拠成果物は [_band_ab/cold_pair/](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/)、run 索引は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:25) です。

結論: **最初に、TP→CPG 変換後の実効設定で熱伝導則だけを比較する CFD 0 step の A/B を行い、定 Prandtl 数への復元を確認してください。**

第 1 仮説: **記載された変換案をそのまま適用すると、forge は定数熱伝導率、SU2 は定 Prandtl 数となり、冷却比較に物性の不一致が混入する。** 確度: **高**。  
根拠: 上表の設定生成・既定値に加え、[gasProperties_d.cu:97](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/gasProperties_d.cu:97)。同じ式を Python で評価すると、cp=1360・Pr=0.72 の λ は **300 K で0.034886、1500 K で0.099507 W/(m·K)**。テンプレートの定数 `0.0257` を残した場合、それぞれ **26.3%／74.2%低い**値になります。これは物性差であり、流れ場の誤差量ではありません。  
反証条件: 実際の生成済み設定では既に `thermCondMethod: 1` が有効で、対象温度範囲の μ・λ が SU2 指定則と一致していること。

第 2 仮説: **断熱 TP 由来の固定 y_b が、冷却・CPG の外縁を代表せず、δ_loc・θ_r の比較に基準位置依存を残す。** 確度: **中**。帯位置の差は確認済み、積分量への寄与は未確認。

第 3 仮説: **未静定の熱場・乱流場が、現在の TP 冷却比に残っている。** 確度: **中**。上記 VERDICT が根拠で、最終的な比への寄与量は未確認です。

判別 A/B:

- **A**：案どおり変換した実効設定。**B**：そこから **`physProp.thermCondMethod` だけを `1`** に変更。
- **追加 CFD は0 step**。同じ温度点・cp・Sutherland 定数・Pr で、設定が選ぶ μ(T)・λ(T) を評価し、SU2 指定則と照合する。
- **A が不一致、B の最大相対差が1e−6以下なら第1仮説を支持**。両方が一致すれば「今回の変換に熱伝導復元漏れがある」という仮説は棄却。B も不一致なら他の物性設定が未整合であり、CFD 比較は判定不能とする。

やらない方がよいこと: **`viscMethod: 1` だけで物性が揃ったとすること、旧 θ の実績を新 θ_r の保証にすること、不一致を直ちに forge のバグとすること、素／生産 SST の差を dilatation 単独の効果と呼ぶこと。**

呼び出し側の前提への異議: **「収束した TP 解」「2キーだけで同じ SST」「5枚 STEADY なら比較可能」は受け入れません。** 一方、既実施の座標 A/B は保存記録で最大 **0.059756%**、登録した1%未満を確認できました。この座標不整合を主因候補へ戻す理由はありません。

不足情報: 新しい6 run の実効設定・物性変換と初期化の実装、元 run の全場 HDF5・全残差系列、使用バイナリと対応ソースが手元にありません。したがって物性復元漏れは**案を記載どおり実施した場合の指摘**であり、既に新 run で発生したとの認定ではありません。

ファイル変更・forge 起動はしていません。**plan 未反映**。反映先は `plans/active/tooling-nozzle-isothermal-wall-chain.md` の **§4.3、§5.1 #26、§6 の新しい CPG 照合条件**です。
```

## 出力形式 (この形のまま)

```
結論: <次にやる一手を 1 文で>
第 1 仮説: <内容>  確度: <高/中/低>
  根拠: <ファイル:行 / run パスと数値>
  反証条件: <何が観測されたらこの仮説は誤りか>
第 2・第 3 仮説: <あれば 1 行ずつ>
判別 A/B: <変える設定 1 点、回す長さ、見る量>  → A なら … / B なら …
やらない方がよいこと: <呼び出し側が取りそうな誤った一手>
呼び出し側の前提への異議: <ブリーフの枠組み・除外判断・指標の定義で受け入れなかったものと理由。無ければ「無し」>
不足情報: <あれば>
```
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
