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

## ブリーフ (`notes/reviews/briefs/2026-10-15-m55-two-coef-contur.md`)

# 諮問ブリーフ: M5.5 の境界層の厚み δ を CONTUR × 倍率 c + 一定のずらし δ₀ で作る方式 (§4.3) の点検 (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1: plan §4・§6 を新しく書いた、4: 承認済みの手順に無い再実行の前)
- plan: `plans/active/tooling-nozzle-card-local-mesh.md` **§4.3 (新規)・§6 V5 (新規、許容はユーザ確認待ち)・§5.1 #3k〜#3o**
- 直前の諮問: `notes/reviews/2026-10-15-m55-p1-axis-wave-diagnose.md` (pass 1 の軸の波。採否は §6.1)
- 関連: `procedures/nozzle-design-workflow.md` 4. (C2 方式 = M6 の生産: CONTUR の摩擦係数倍率 k_f を出口で合わせる)、`plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #8d・#10〜#12 (CONTUR の熱閉包の危険、4 係数化は不採用見込み)

## 経緯 (観測事実)

1. pass 0 `run_0003_ns_cont` (壁 = 理想の形 + CONTUR の δ、cf_scale 1): 試験部の軸の M は +0.35〜+0.42 % でほぼ平ら、出口コア M 5.52042 (+0.37 %)。
2. pass 1 `run_0005_ns` (壁 = 理想の形 + pass 0 の NS の抽出 δ_E をそのまま、`wt_card --ns-correction`): 出口コア M 5.50145 (+0.03 %)、試験部の軸に M の波 −0.76〜+0.55 % (r/r_w ≈ 0.27 で −0.32〜+0.18 %)。設計の理想の形の Euler は ±0.02 %。前の諮問の結論は「pass 2 は保留、抽出の径方向の標本化の A/B (#3k) を先に」。
3. M6 の生産 (case/45、10-07) は C2 で、NS の δ_E を壁に直接入れていなかった。2026-09-04 の生産レシピの記録にも「反復ごとの壁変化が特性線で末端の軸に届き ±0.5 % の軸の波を作る」とある。
4. **ユーザ決定 (2026-10-15)**: 「摩擦係数でなくて、倍率かけて合わせ込むほうが、割り切っててよい」「コンターの倍率だと、そのマイナスの厚み再現できないよな」→ 案 (1)「倍率 c + 全長で一定のずらし δ₀」を選んだ (「1かええな」)。

## 方式 (plan §4.3)

- δ(x) = c·δ_C(x) + δ₀。δ_C は `runner_axismach.integral_delta_r` の出力 (`integral_bl` → 5 次 P-spline 平滑化、cf_scale 1)。δ₀ は実寸 [m] で全長一定 (壁の傾きを変えない)。スロートより上流は今どおり δ を足さない (`pw_upstream: poly`)。
- 較正: pass 0 の NS (`run_0003_ns_cont` res_200000、r_t 0.0999084914 m) の平滑化後の δ_E のスロート (x = 0) と出口 (最後の列、x 8.000044 m) の 2 点。δ_C は同じ r_t の値 (= pass 0 の物理壁の δ = `delta_r_next.csv` の `delta_in_prev`)。
  - δ_C: スロート 0.09159 mm・出口 49.3301 mm。δ_E: −0.24396・41.7109 mm。→ **c = 0.852073、δ₀ = −0.32201 mm**。
  - 途中の差 (2 係数 − NS、pass 0 の比較): x 0.25 m +0.19、0.5 m +0.50、1 m +1.20、2 m +1.52、3 m +2.04、4 m +1.99、5 m +1.97、6 m +1.66、7 m +1.12 mm (2 係数が厚い)。
  - 参考: pass 1 の NS の δ_E はスロート −0.2433・出口 41.433 mm (壁が違う)。
- 実装 (#3n): `deltastar_initializer.delta_scale` (c、既定 1)・`delta_offset_m` (δ₀、既定 0) を `integral_delta_r` の平滑化の後に当てる。同じ関数を寸法の逆算 (`deltastar_loop._sizing_delta_supplier` → `solve_rt_throat`、カードの物理スロート径 200 mm の固定点)・`prepare_ns`・`_mesh/hand_mesh.py walls` が使う。カードから `overrides` で渡す (`_OVERRIDE_SECTIONS` に `deltastar_initializer` がある)。
- 流れ (#3o): カード → 格子 (`hand_mesh.py build`、#3e と同じ係数) → Euler (理想の形、NS の判定のものさし) → NS `deltastar_loop --init-integral` (初期場は run_0005_ns の最後の場から補間、`--stages full`、精度は既定 float、出力 2,000、停止規則 窓 2 万・保持 4 万・上限 20 万) → §6 V5 → 成果物。δ* の反復はしない。出口コア M が V5b を外れたときだけ、その NS で同じ 2 点の較正をやり直して 1 回だけ作り直す。
- スロートの δ が NS に合うので、r_t は pass 1 の 0.10024473616 m にほぼ戻る見込み (pass 1 の r_t は δ_E のスロート −0.244 mm で解いた)。

## §6 V5 (事前登録の案、許容はユーザ確認待ち)

- V5a 寸法: 成果物の寸法の関門のカードの 5 項目 PASS。
- V5b 出口コア M (停止規則の量): 5.5 ± 0.05 % (M6 の生産のゲートと同じ幅)。
- V5c 試験部の一様性 (前の諮問の提案): x 4.0〜7.95 m の M の最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %。軸の節点の既知の癖 (静圧 0.5〜0.7 % 低) を避けて r/r_w = 0.1 の線で測り、軸の節点の値は併記。
- V5d 流れ: 停止規則・`check_convergence --segment`・系列の `check_quasisteady`・壁解像 (y₁⁺ > 1 の面積 ≤ 2 %)・非有限 (#3e (6) と同じ)。

## 問い

1. §4.3 の方式・較正 (2 点、pass 0 を元にする、δ₀ を実寸で一定、平滑化の後に当てる、上流は足さない) に穴はあるか。特に: (a) 較正の元を pass 0 にする根拠 (同じ CONTUR の形の壁の上の境界層) は妥当か、pass 1 や両方の平均にすべきか。(b) δ₀ を当てると上流の 5 次多項式 (`pw_upstream: poly`) とスロートでの接続 (値・傾き・曲率) に問題が出ないか。(c) 寸法の逆算 (CFD 前の許容差 1e−7 m) と組み合わせて物理スロート 200 mm が保たれるか。
2. 途中の +1〜+2 mm の厚すぎ (x 1〜7 m) が試験部の一様性 (V5c) に効く見込み。M6 の C2 はこの種の途中の差を持ったまま合格した。事前に見積もる方法はあるか (なければ「NS で判定」でよいか)。
3. V5 の基準 (特に V5c の線と幅、V5b の幅) は妥当か。前の諮問の提案 (最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %) を r/r_w = 0.1 の線に移すのは緩め・きつめのどちらか。pass 0 (r/r_w ≈ 0.27 で −0.16〜+0.40 %) と pass 1 (同 −0.32〜+0.18 %) はどちらも V5c を満たさない。
4. Euler: 新しい r_t が pass 1 の r_t (Euler `run_0004_euler` の r_t) と何桁一致したら使い回してよいか、それとも常に回し直すべきか (約 4 分)。抽出・停止規則は NS と Euler をそれぞれの r_t で割った単位で比べ、理想の形は実寸の長さ (拡大部 8 m など) を r_t 単位に直して作るので、r_t の差は拡大部の長さの r_t 単位の差になる。
5. 前の諮問で決めた #3k (抽出の径方向の標本化の A/B) は、2 係数の方式では δ_E を 2 点 (スロート・出口) でしか使わない。#3k を後回し・取りやめにしてよいか。代わりにスロートと出口の δ_E の値の標本化の感度だけを見るべきか。
6. 初期場を run_0005_ns (波のある場) から補間するのと、run_0003_ns_cont (pass 0) からにするのとどちらがよいか。

## 関連 plan 全文 (`plans/active/tooling-nozzle-card-local-mesh.md`)

```markdown
# 要求カードの格子を場所ごとの第一層と流れ方向の密度で作る

## メタ

- **area**: `tooling` (`design/forge_design/intake/wt_card.py` の `estimate_mesh`・雛形 `intake/recipes/axismach_*`、`meshing/mesh2d.py` は使うだけ)
- **status**: `draft` (2026-10-15 起票)
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「要求カード」「格子」の節
  - skill `nozzle-design`・`design-intake`
- **related_plans**:
  - [`tooling-generator-defaults-and-launch-gate.md`](tooling-generator-defaults-and-launch-gate.md) — 段取り 2 の「生成器の既定 = 最新の生産レシピ」と `mesh2d` の壁法線の近壁層の既定化 (§6.4)。この plan はカードの格子の作り方だけを先に揃える
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) §5.1 #36r — M6 の H の上流の第一層 (NS の壁量から表を作り直した例)
  - [`tooling-nozzle-mknot-fast.md`](tooling-nozzle-mknot-fast.md) — この枝の元 (M_knot の新しい規則)
- **created**: `2026-10-15`
- **owner**: `M6ノズルの設計` セッション (作業ツリー `forge-cardmesh` / `feature/card-local-mesh`、feature/mknot-fast 7713f84c から)

## 1. 目的

要求カードの雛形 (`axismach_explore_v1`・`axismach_final_v1`) の NS 格子は、第一層を 2 値 (スロート / それ以外) + ブレンドで持ち、流れ方向はスロートの細分だけ、壁に直交する近壁層も無い。
M6 の生産の格子 (case/45 の H: `wall_first_frac_table`・`x_density_table`・`wall_normal_layer`) より古い作り方で、M5.5 (case/71) では最も厳しい場所の薄い第一層が拡大部 8 m 全体に掛かり、
壁際の AR ≤ 5000 のために ni 10,177・約 99 万点になった (ユーザ 2026-10-15「あまりにメッシュ多すぎる、古いやり方でメッシュ作ってない？」「ほかのケースでちゃんと流れ方向に分布つけてたし、壁面に直交するようにメッシュ作ってた」)。
**目標**: カードから作る格子を、場所ごとの第一層の表・流れ方向の密度の表・壁に直交する近壁層で自動で作る。M5.5 は ユーザ決定「a」(2026-10-15) で、見積もりの式の偏りを M6 の実績で補正した約 45 万点の格子で Euler → NS まで回す。

## 2. スコープ

- **やる**: `estimate_mesh` が `wall_first_frac_table`・`x_density_table`・`ni` を作る (作り方 `local_v1`)、`wall_normal_layer` を雛形の既定に、見積もりの補正の係数、記録、M5 デモ (case/70) と M5.5 での確認。
- **やらない**: 見積もりの式そのもの (Eckert の参照温度 + 平板の摩擦係数) の置き換え、`nj` と径方向の配点 (core-grid の範囲)、Euler の格子 (`mesh_euler`、全域 0.005 のまま)。

## 3. 関連 docs と前提 (観測事実)

- **見積もり** (2026-10-15、AWS A `~/m55_mesh_est/`、手元 scratchpad `m55_mesh_est/`; コードは変えずに `estimate_mesh` (7713f84c) の式で計算): M5.5 (r_t 0.0999088 m・M_knot 2.80・nj 97)、y₁⁺ 目標 0.8、AR ≤ 0.97 × 5000 = 4850。

  | 作り方 | ni | 節点 |
  | --- | --- | --- |
  | (A) 今: 第一層 2 値 + ブレンド | 10,177 | 987,169 |
  | (B) 第一層を表、流れ方向は今の分布 | 8,749 | 848,653 |
  | (C) 第一層を表、流れ方向の密度も表 (dx ≤ 4850·y₁·cosθ_w、上限は雛形の ni 2000 相当の間隔) | 3,163 | 306,811 |

  (C) の格子を `mesh2d.generate_axisym_mesh` で実際に作り、AR 最大 4850.0・隣り合う dx の比の最大 1.031 (径方向 1.134)。
  前提を変えると: y₁⁺ 目標 1.0 で 271,503 点、第一層を全域 ×0.5 で 490,335 点。入口の端は y₁,req → 0 (y₁ ∝ s^0.1) なので、断面の 1 点目 (s = 0.011 r_t) より手前はその値で置き換える (生成器の y₁⁺ の検査も入口の節点を除いている)。
- **見積もりの式の偏り** (M6 の H の表は NS の壁量から作ったもの): H の表 ÷ 見積もりの要求 = 配管 2.2〜3.6、縮流部 0.61〜2.2、スロート付近 0.51〜0.61、拡大部 0.50〜0.69。見積もりはスロート・拡大部で第一層を約 2 倍厚く見積もる (壁のせん断を低く見積もる)。配管では逆に薄く見積もる。
  H の表は断熱の NS の y₁⁺ に冷却の倍率 8.5〜9.4 を掛けて y₁⁺_cold ≤ 0.8 にした 300 K 用 (YAML のコメント)。M5.5 は同じ燃焼ガスの組成・壁 350 K。
- M6 の H は今の配点で 46 % の列が AR 4000 超 (既に限界近く)。(C) にしても ni 4,719 → 4,294 (1 割減)。
- `mesh2d` は `wall_first_frac_table` と `wall_first_frac_throat`・`wall_first_up_*` の併用を拒む (`mesh2d.py:100`)。
- `wall_normal_layer` は PR #8 で main に入った (c096d66c)。M6 の H は `[0.02, 0.8]`。カードの雛形には無い。

## 4. 設計方針 (2026-10-15 起案、上位に諮る前)

### 4.1 作り方 `local_v1` (雛形の `mesh.method: local_v1` で選ぶ; 新しい雛形の既定)

1. 設計壁に沿って、今の `estimate_mesh` と同じ式で要求の第一層 y₁,req(s) を出す (y₁⁺ 目標 `Y1PLUS_TARGET` 0.8、入口の端は s = 0.011 r_t の値で頭打ち)。
2. **補正**: y₁,req に区間ごとの係数 c(x) を掛ける。係数は M6 の H の表 ÷ 見積もりから、保守側に取る (表が見積もりより厚い区間は 1.0 = 補正しない): 配管・縮流部 1.0、スロート付近と拡大部 0.5 (H の比の最小 0.50 に合わせる)。区間の境 (縮流部の終わり x = −2 r_t 付近) は滑らかにつなぐ。係数は雛形の `mesh.first_layer_calibration` に表で持ち、出所 (M6 の H の表と見積もりの比、2026-10-15) を記録する。
3. 第一層の表 `wall_first_frac_table`: 0.125 r_t (x < 10) / 1 r_t (x ≥ 10) の点で、前後の区間の c·y₁,req の最小 (下包絡) を局所半径で割る (log 線形補間で要求を超えない)。
4. 流れ方向の密度の表 `x_density_table`: 目標の間隔 h(x) = min(4850·y₁,表·cosθ_w, h_cap(x)) の逆数の上包絡。h_cap は雛形の分布で ni 2000 にしたときの各位置の間隔 (主流の解像の上限)。AR ≤ 4850 と dx ≤ h_cap を満たす最小の ni を二分法で求め、10 の倍数に切り上げる。`throat_refine` は 1。
5. `wall_normal_layer: [0.02, 0.8]` (M6 の H と同じ) を雛形の既定にする。
6. 節点数は `check_node_budget` で判定 (上限を超えたらユーザの了承と `max_nodes`・`max_nodes_reason`; 今のまま)。
7. 記録: card_log の `mesh_estimate` に y₁,req・係数・表・ni・節点数・AR の最大・dx の比の最大。

### 4.2 後で直す仕組み

NS の後に `check_wall_resolution.py` (壁ごと・面積の重み・停止の step) で y₁⁺ を確かめ、超えていれば NS の壁量から表を作り直す (M6 の #36r の台本 `h370_wall_envelope.py` の `design` の手順)。カードの流れでは「NS の後の壁解像の確認 → 必要なら表を直して NS をもう 1 回」を標準の手順として skill に書く。

### 4.3 M5.5 の境界層の厚み δ の作り方: CONTUR の積分式 × 倍率 c + 全長で一定のずらし δ₀ (2026-10-15 ユーザ決定「1かええな」、上位に諮る前)

- **経緯**: pass 1 (#3e) はカードの `--ns-correction` で NS の抽出 δ_E をそのまま壁に足し、試験部の軸に M の波 (−0.76〜+0.55 %) が出た。M6 の生産は C2 (CONTUR の摩擦係数倍率 k_f を出口で合わせる) で、NS の δ_E を直接入れていなかった (#3m)。ユーザ: 「摩擦係数でなくて、倍率かけて合わせ込むほうが、割り切っててよい」「コンターの倍率だと、そのマイナスの厚み再現できないよな」→ 案 (1) 倍率 + 一定のずらし。
- **式**: δ(x) = c·δ_C(x) + δ₀。δ_C は今の CONTUR の積分法 (`integral_delta_r`: `integral_bl` → 5 次 P-spline、cf_scale 1) の出力、c は無次元、δ₀ は実寸 [m] で全長一定 (壁の傾きを変えないので波を作らない)。スロートより上流は今どおり δ を足さない (`pw_upstream: poly`)。
- **較正**: 同じ形の壁 (CONTUR の厚み、cf_scale 1) で回した NS の δ_E (平滑化後、未緩和) のスロート (x = 0) と出口 (最後の列) の 2 点に合わせる: c = (δ_E,e − δ_E,t)/(δ_C,e − δ_C,t)、δ₀ = δ_E,t − c·δ_C,t。δ_C はその NS と同じ r_t で評価した値 (= その NS の物理壁の δ)。pass 0 (`run_0003_ns_cont` res_200000、r_t 0.0999084914 m) から: δ_C スロート 0.09159・出口 49.3301 mm、δ_E −0.24396・41.7109 mm → **c = 0.852073、δ₀ = −0.32201 mm**。新しい r_t では CONTUR を新しい r_t で計算し直してから c・δ₀ を当てる (CONTUR の r_t の依存は式が持つ)。
- **割り切り**: 係数 2 つなので途中は合わない (pass 0 の比較で x 1 m +1.2・3 m +2.0・5 m +2.0・7 m +1.1 mm、NS より厚い)。壁は CONTUR のなめらかさのまま。
- **実装**: `deltastar_initializer` に `delta_scale` (c、既定 1) と `delta_offset_m` (δ₀、既定 0) を足し、`runner_axismach.integral_delta_r` で平滑化の後に当てる (既定のままの YAML はビット同一)。同じ関数を寸法の逆算 (`_sizing_delta_supplier`)・`prepare_ns`・`hand_mesh.walls` が使うので、スロート径 200 mm に合わせた r_t の解き直しも同じ δ で行う。カードからは `overrides: {deltastar_initializer.delta_scale: c, deltastar_initializer.delta_offset_m: δ₀}`。
- **NS の位置づけ**: 壁を決めた後の確認の 1 本 (δ* の反復はしない)。出口コア M が §6 V5b を外れたときだけ、その NS で同じ 2 点の較正をやり直して 1 回だけ作り直す。
- **Euler**: 新しい壁の理想の形 (MOC) で 1 本回す (NS の判定のものさし)。スロートの δ が NS に合うので r_t は pass 1 (0.10024473616 m) にほぼ戻る見込み — 使い回すかは r_t の差を見て諮問で決める。

## 5. 実装ステップ

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~plan 段レビュー~~ | 判断: 2026-10-15 NO-GO (C0/M7/m1) を全件採用 (§6.1)。§4・§6 は書き直して再レビュー | F |
| 3a | **M5.5 の手組みの格子 (先行)** | ユーザ決定「a」(2026-10-15): case/71 の問題 YAML の mesh 節を、正規化した要求 f_req = c(x)·y₁,req/r_w の下包絡の表 (補間の区間の中も検査・物理壁で再検査)・流れ方向の密度の表・`throat_refine 1`・`wall_normal_layer [0.02, 0.8]`・`msh_digits 17`・`max_nodes` (ユーザ了承「a」、約 45 万) で手で組む。係数 c: 配管 1.0、縮流部 0.6、スロート・拡大部 0.5 (M6 の H の表 ÷ 見積もりの比の最小、初期格子の経験の係数)。品質: `check_mesh_quality --ar-max 5000`・`mesh_ar_relax_check`・変換後の第一層厚の照合。生成器の既定は変えない。**実装 (2026-10-15、implementer)**: `case/71.m55_combustion_wt/problem_m55_ns.yaml` (`_mesh/hand_mesh.py build`; 生成は `card_m55.yaml` に `mesh.ni: 3000` を足した `problem_m55_card.yaml`、M_knot 自動 2.80)。ni 4,800 × nj 97 = 465,600 (NS・Euler)。係数の拡大部は比の最小 0.49874 < 0.5 なので 0.49 (0.5 なら ni 4,780)。CFD 前の物理壁の格子 (forge 不起動、A `~/m55_mesh_hand/check/`): 品質 PASS (AR 最大 4849.4)・AR > 1000 の 47,607 セルでスキュー最大 0.0015 (違反 0; 出口面の非平面 0.62 mm は記録だけ)・第一層の実寸 / (c·y₁,req) 最大 0.999999 (超過 0)・見積もりの y₁⁺ 最大 0.797・dx の比の最大 1.033。Euler PASS (AR 最大 6.1)。再現の入力と記録は `case/71.m55_combustion_wt/_mesh/` **判断 (2026-10-15、主セッション)**: (1) 拡大部の係数 0.49 (M6 の比の最小 0.49874 を 2 桁に切り下げ) — 採用。(2) 正規化の半径を max(設計壁, 物理壁 (CONTUR の初期壁)) — 採用 (NS の pass で壁が動けば半径の比だけ第一層がずれる: NS の後の V4 で見る)。(3) FP64 のビルドで回す (第一層/半径が 3.7e−7 まで下がる): AWS A の `~/forge-integ-dflow2-fp64` (2fd36b0c、solver の源は main と同じ) を使う。(4) 探索の雛形 (`axismach_explore_v1`) の物理壁では停止規則が使えない → カードの雛形を `axismach_final_v1` に替えて問題 YAML と `_euler.yaml` を作り直し、手組みの mesh を当て直す (#3b) | O |
| 3b | M5.5 を最終の雛形で作り直し、Euler → NS | カードの `recipe: axismach_final_v1` で問題 YAML・`_euler.yaml` を作り、`hand_mesh.py build` で mesh・mesh_euler を当て直して品質を確かめる。Euler `run_0001_euler` → NS `run_0002_ns` (停止規則、上限 20 万) を AWS A の FP64 のビルドで。**途中 (2026-10-15、implementer、YAML と Euler の準備)**: カード → `axismach_final_v1` で生成 (A `~/forge-cardmesh` 94790833): `case/71.m55_combustion_wt/problem_m55_final_card.yaml`・`_euler.yaml` (M_knot 自動 2.80 で不変、r_t 0.0999084914252905 m、寸法 PASS)。`hand_mesh.py build --plan-item "§5.1 #3b" --axis-cap-frac 0.03` → `problem_m55_final.yaml` (NS)・`problem_m55_final_euler.yaml` (Euler; `wt_card.write_euler_companion` の写し、生成器の写しと mesh 以外が同じことを照合)。係数・正規化・ni 4,800 は #3a と同じ。NS の径方向はユーザ決定 (2026-10-15「上限3%」) で `axis_cap_frac 0.03`、nj 92 (全断面の隣接比の最大 1.1974 ≤ 1.2、nj 91 は 1.2013) → NS 441,600 点、Euler は nj 97・465,600 点のまま。検査 (FP64 の変換器、forge 不起動): NS 品質 PASS (AR 最大 4849.0)・AR 緩和 違反 0・変換後の第一層 実寸/(c·y₁,req) 最大 0.9999997 (超過 0)・見積もりの y₁⁺ 最大 0.797・軸側のセル / r_w 0.030; Euler PASS (AR 最大 6.1)。NS の準備の試行で停止規則 use (上限 200,000)。`run_0001_euler` は準備のみ (A、未起動)。記録 `case/71.m55_combustion_wt/_mesh/hand/final/` **径方向の決定 (ユーザ 2026-10-15「上限3%にしてみますか」)**: 雛形の等比 (壁から軸まで比 1.10〜1.14、nj 97、軸側のセルが局所半径の 11〜12 %) を、`axis_cap_frac 0.03` + 近壁の比 ≤ 1.2 (nj は全断面で比 ≤ 1.2 の最小 = 92) に替えた。軸側のセルは全区間で 3.0 %、節点 441,600。Euler の格子は壁面垂直方向がほぼ一様 (第一層 0.5 %・比 1.014) なので変えない。**投入 (主セッション)**: `run_0001_euler` を AWS A で起動 (FP64 `~/forge-integ-dflow2-fp64`、soft → mid → 本段 cfl 1.0・12,000 step、約 12 ms/step) ユーザが nj 92・上限 3 % の NS 格子を ParaView で確認 (2026-10-15「オッケーそうでした」)。`run_0001_euler` 完走 (227 s、出口の軸の M 5.5002、dM 最大 0.06 % Md、NaN なし; `check_convergence --segment` は NOT CONVERGED = 残差 2.1〜2.5 桁で停滞、M6 の Euler 参照と同じ型 → δ* の参照としては未適格と記録)。`run_0002_ns` を起動 (停止規則 use・上限 20 万、物理壁 single_bspline) **NS の未到達の診断 (2026-10-15、implementer、forge 0 step)**: codex 諮問 `notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.md` (incomplete: 証拠不足) の証拠を A から回収して `case/71.m55_combustion_wt/_runs/` に (h5 なし: 全残差履歴の gz・stage_manifest・solverConfig・bcond・prepare_info・MAIN_STAGE_STARTED・forge_launches・stop_rule.json の全体・帯の抽出・ドライバの log・forge の log の先頭と末尾、Euler の同類)。バイナリは FP64 (`flowFormat.hpp` double)、sha256 3a9e347b…、`~/forge-integ-dflow2-fp64` (HEAD 2fd36b0c、ビルドより新しいソルバの源 0 件)。判定 (原文 `_runs/run_0002_ns/judge/`): (1) 記録の stage_manifest に本段が無い (見張りの `save_record` が OSError で落ち、`run_staged_ns` の `_record(cfg_main, "main")` に届かなかった) ので `--segment` をそのまま回すと S1 → S2 を判定する。本段を足した写しで区間 = 本段 0〜183,998 → `check_convergence --segment` NOT CONVERGED (still converging; ρ 1.8・ρE 2.2・運動量 2.4・k 2.9・ω 2.7 桁)。各列の最小は 157,795〜172,251 step で、直近 2 万 step の log10 の平均はその前の 2 万 step より 0.07〜0.44 桁高い。(2) 停止規則の系列 (1,000〜183,000) の `check_quasisteady --series-csv` DRIFTING (r DRIFTING、θ_r 0.4・0.7 TRANSIENT-UNSETTLED、θ_r 0.95・ṁ_in・ṁ_out・出口コア M STEADY)。(3) res_183000: 非有限 0、ρ・P・T > 0、P > P0 は壁の行 (U = 0) の 298 点で相対 ≤ 5.8e−8、P0 > Pt (1 + 1e−6) が壁から 27 点以上離れた主流〜軸の 50,441 点 (最大 +0.24 %)。(4) 壁解像 (`wall`、area、step 183000) PASS: y1+ > 1 の面積 0.0 %・最大 1.004 (停止の step ではないので V4 の判定ではない)。(5) 実効の設定: 問題 YAML から生成器で作り直した S1・S2・本段の solverConfig と bcond が run の記録と一致 (sha256・cfg_fnv・bcond_fnv)。(6) θ_r の測定帯の A/B (同じ保存場 100k〜183k の 12 枚・同じ Euler res_12000・同じ 3 断面、帯だけ変える; `_runs/run_0002_ns/theta_band_ab.json`): D_θ [%/2 万 step] A 0.97・0.91・0.54、B 2.33・1.96・1.75 → 事前に固定した判定で **「帯の選び方だけで未到達を説明できる」を棄却** (3 断面とも両腕で > 0.05 %)。別に、帯 A (res_1000 から) は後期の縁 y* を使用域の全列で覆っていない (帯 B/A 1.40〜5.52 倍)。容量: A の空き 16 GB、全場 1 枚 100,605,412 B、壁の出力 1 枚 1,170,416 B。記録 `case/71.m55_combustion_wt/_runs/README.md` **判断 (2026-10-15、codex diagnose `notes/reviews/2026-10-15-m55-ns-unsettled-2-diagnose.md`、全件採用)**: 長い延長は保留し、先に既存の保存場で全圧の A/B (A = 出力の原始量・h0 から、B = 保存量 ro・roU・roe・roY からソルバと同じ規約で再構成した状態から、forge 0 step、100k〜183k の 12 枚) を行う (#3c)。manifest の写しに本段を足した判定は採用 (原本を残し「事後復元・本段 0〜183,998」と明記、正常終了の証明ではない)。追加 約 16.8 万 step は到達の予測に使わない。延長に進むときは: 新しい run で同じメッシュの保存量を index restart (前段は回さない)・出力 2,000 step・`evaluate.stop_rule: true`・保持 4 万・**追加 20 万 step を予算の上限** (到達の予測ではない、2 万 step ごとに再評価)・新しい帯と Euler の res_12000 (sha256) を固定して θ_r は別の系列・残差は同じ設定の区間だけ連結 (183k 以後を重ねない)・容量は約 10.5 GB で投入の時に測り直す | O |
| 3c | 全圧の超過の A/B (forge 0 step、**済み: 複合要因、超過の大部分は保存量の状態にある; 3 回目の諮問は thermoFloat の 1 更新の A/B を提案 → ユーザ決定で延長を全部 double にして代える (#3d)**) | 上の判断のとおり。節点ごとの ΔP0/Pt、超過の点数・分位点・位置 (入口・軸・内部)、エネルギーの逆算の残差。判定 (事前に固定): 両経路の差が全時点で ≤ 1e−6 Pt かつ超過が残る → 出力の状態の不整合を主因とする仮説を棄却; B で超過が消え A/B の差が超過を説明 → 保存場の側という仮説を棄却し出力・後処理の経路を直す; 一部だけ → 複合要因として記録。**結果 (2026-10-15、implementer、forge 0 step)**: 台本 `case/71.m55_combustion_wt/_runs/run_0002_ns/p0_ab.py` (A の `~/m55_diag_0002/p0ab/` で実行、ツールは A の `~/forge-cardmesh` b0d528bf = 手元と solver_density_cuda/ の差 0; 原文 `_runs/run_0002_ns/p0_ab.json`・`p0_ab_series.csv`・`p0_ab_series_quasisteady.txt`・`p0_ab.stdout.txt`)。同じ 12 枚 (100k〜180k の 1 万 step ごと + 181k〜183k)・同じ物性の記録 (res の属性 → `resolved_species_fc8ba1a8828147e3_70cc58aa292c9504.yaml`、sha256 70cc58aa…)・datum 298.15 K。B は `dependentVariables_d.cu` の thermoFloat 1・凝縮 off の経路 (ρ の床・float32 の組成の非負化と正規化・e = roe/ρ − \|U\|²/2・T は e_mix の Newton・P = ρ R_mix T・h0 = e_mix + P/ρ + \|U\|²/2、h0_includes_k 0) を写し、全温・全圧は A と同じ式。A は既報 (183k: 50,441 点・最大 +0.23580 %) を再現。(1) \|P0_A − P0_B\|/Pt の全節点の最大は 8.59e−6 (160k)〜1.03e−4 (100k)、1e−6 を超える節点は 7,995〜116,254 (内部が 97〜99.5 %、壁の行 0; x 0〜1 m の割合は時点で 10〜75 %)、中央値 2.5e−8〜2.4e−7。(2) B でも超過が残る: 49,862〜79,964 点 (A 50,161〜80,210、A∩B は A の 97.7〜99.5 %)、B の最大 P0/Pt − 1 は +0.2232〜+0.6166 % (A の最大との差 ≤ 1.4e−5)。位置は A と同じで、入口の列 0・壁の行 0・壁から 27〜91 点・x −1〜8 m (軸の行は B で 100k〜160k に 1,947〜3,782 点、180k〜183k に 51〜107 点; A は 1,947〜3,803 点と 64〜150 点)。(3) A の超過の点で A − B が超過量を説明する割合 (Σ) は 0.11〜0.65 %。(4) エネルギーの逆算の残差 \|e_mix(T_B) − e\| ≤ 2.6e−9 J/kg (T の床・P の床・ρ の床・Y の負値の切り捨て 0 点)。出力の T を保存量の組成と e に当てた差は −35〜+32 J/kg (中央値 0.05〜0.12 J/kg)。T0 > Tt + 1 K は A・B とも入口の列の 9 点 (100k だけ A 59・B 57 点)、B の静温も同じ 9 点で Tt を最大 2.62 K 超える。組成の float32 の丸めの寄与 (B と、double で正規化した B_dbl の差) ≤ 3.4e−7 Pt。**判定 (事前に固定した規則どおり): 複合要因** — 条件 1 (全時点で差 ≤ 1e−6 Pt) も条件 2 (B で超過が消える) も満たさない。差は超過量の 1 % 未満で、超過の大部分は保存量から作り直した状態にも残る (過渡の原因の確定には進まない)。時系列 (最後の時点の A の超過の集合 S、50,441 点を固定、12 枚) の `check_quasisteady --series-csv` (既定の閾値): **DRIFTING** — 超過の点数は A・B とも STEADY、S の平均の ΔP0/Pt は DRIFTING (−5.06e−5 → +1.46e−4、182k で 1.0e−7 下がった 1 回を除いて増加)、S の平均の h0 − h0_in は DRIFTING (−55 → −2,152 J/kg で単調に減少)、軸の近く (S ∩ j ≤ 2) の h0 は OSCILLATING。A と B の列の判定は全て同じ | O |
| 3d | **NS の延長を全部 double で** (ユーザ決定 2026-10-15「どうせやるなら全部doubleかな」) | 3 回目の諮問 (`notes/reviews/2026-10-15-m55-ns-unsettled-3-diagnose.md`) は `thermoFloat` 1/0 の 1 更新の A/B を先に、延長は保留としたが、ユーザの決定 (「どうせやるなら全部doubleかな」「thermoだけとは限らないよね」「化学種の拡散とかはきにしなくていい」) で**延長を切り替えのある double の経路を全部入れて回す**: `physProp.thermoFloat: 0`・`time.deltaT.implicitSolvePrecision: 1`・`FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1` (FP64 のビルド)。切り替えの無い float (化学種の拡散の D_s・h_s・X、リミッタの長さの `sqrtf`、XMF の型) は残す (ユーザ了承)。混合精度の経路を延長から外すことで、H2 のうち `thermoFloat` の寄与を延長の結果に混ぜない。条件は 2 回目の諮問のとおり: 新しい run `run_0003_ns_cont` に run_0002_ns の res_183000 を同じメッシュで index restart (`restart_field.py`)、前段は回さない、出力 2,000 step・`evaluate.stop_rule: true`・保持 4 万・追加 20 万 step の予算 (2 万 step ごとに再評価)、帯と Euler の参照 (run_0001_euler/res_12000、sha256) を固定し θ_r は別の系列、判定区間は延長の本段だけ (設定が違うので run_0002_ns とは連結しない)、容量を投入の時に測る。監視: 残差・r・θ_r、S の h0・P0、壁の積分熱流、入口・出口のエネルギー流束、入口の 9 点の静温の超過。入口の面の実効の流入とエネルギー流束は延長の序盤に確かめる (違反なら止めて直す) **投入 (2026-10-15、implementer)**: AWS A `case/71.m55_combustion_wt/run_0003_ns_cont` (台本・記録 `case/71.m55_combustion_wt/_runs/run_0003_ns_cont/`、`cont_0003.py`、A の木 f85e4f1b)。res_183000 (sha256 c29044a8…) を `restart_field.py --keep-src-dtype` で同じ nozzle.h5 に index restart (VERDICT: OK、9 量、SRC とビット一致)、`run_staged_ns(stages="none", stop_rule=True)`。run_0002_ns の本段から変えたのは `physProp.thermoFloat: 0`・`time.deltaT.implicitSolvePrecision: 1` (timeIntegration 11・blockDPLUR 1・lowMachPrecond 0 を確認)・`outStepInterval` 2,000 と、環境変数 `FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1` (`run_case.sh` の記録の対象に足して `RUN_PROVENANCE.txt` に記録; 4 つとも forge の起動の行で効きを確認)。停止規則は明示・窓 2 万・保持 4 万・上限 200,000 (= 追加分)、Euler の参照 run_0001_euler/res_12000.h5 (sha256 93cb88b0…) に固定、帯は res_2000 から。FP64 のビルド (sha256 3a9e347b…)、26.9 ms/step (run_0002_ns の約 1.9 倍、上限まで約 90 分)、序盤 (step 0〜14,478) の残差・出力に非有限 0。容量: 投入の前の空き 16.36 GB、見込み約 10.4 GB。**残る float**: 化学種の拡散の D_s・h_s・X の float の写し (`speciesTransport_d.cu:240-356`)、リミッタの長さの `sqrtf` (`limiter_d.cu:398`)、XMF の型の宣言。入口の面の確認 (res_183000、forge 0 step、`inlet_face_check.py`): 92 面のうち 91 が流入・逆流 0 (壁との角の 1 点は U = 0 で ṁ = 0)、流入の面の Ht_R − h0_in は相対 1e−13 以下・合計のエネルギー流束 = ṁ·h0_in (相対 −9.5e−14)・mach ≤ 0.066、T0 > Tt + 1 K の 9 点の面も全て流入で境界の状態は Ts_b 1,249.25 K — 事前に置いた違反の目安 (逆流・Ht_R ≠ h0_in・合計 ≠ ṁh0_in・mach ≥ 1) に当たらない。壁から 1〜5 点の面は節点が遅く冷たいので Ps_b が節点の静圧より最大 +0.27 % (面積の割合 1.4e−5) **結果 (2026-10-15、判定区間は延長の本段 0〜199,999; 原文 `case/71.m55_combustion_wt/_runs/run_0003_ns_cont/judge/`)**: 停止規則 **CENSORED** (200,000、`STOP_RULE_STATE.json` 03:15:22 UTC; θ_r 0.4・0.7 の D 6.32e−4・5.33e−4 > 5e−4、0.95 は 4.84e−4、r の平均 1.07e−4; D は 1 万 step ごとに約 0.82 倍で単調に減少)。`check_convergence` NOT CONVERGED (ρ 1.7・k 1.8・Y 1.7 桁 falling、ρUx・ρUy 2.0・ρE 1.6 桁 flat、ω 3.1 桁)。停止規則の系列の `check_quasisteady --series-csv` DRIFTING (θ_r 3 断面 STEADY・単調増加、漸近値 0.0536868・0.0798884・0.0902438 = 最終の +0.135・+0.116・+0.106 %; ṁ_in・ṁ_out・出口コア M STEADY 121.358・121.355 kg/s・5.52042; r DRIFTING)。壁解像 (`wall`、area、200,000) PASS (y1+ 平均 0.479・p99 0.780・最大 1.002 @ x −1 m、超過 0.0 %)。最後の場: 非有限 0、P0 > Pt (1 + 1e−6) 46,936 点・最大 +0.462 % (切り替え直後の res_2000 で 50,069 点・+0.237 % = double にしても消えない)、入口の列の T0 > Tt + 1 K 9 点 (最大 +3.39 K)。監視 (`monitor_series.csv`): 最大 ΔP0/Pt STEADY (漸近 +0.463 %)、S の平均 ΔP0/Pt・h0 − h0_in DRIFTING、E_in・E_out・Q_wall STEADY、(E_in − E_out + Q_wall)/E_in DRIFTING (−2.5e−4 → +7.37e−4、Q_wall は温度拘束の反力を含まない qwall の台形積分なので非保存の証拠にしない; #3f)。δ の 2 窓 (`two_window_check`、160k〜180k・180k〜200k、許容 0.01 mm) **FAIL (measured)**: δ_target の窓の中の振れ 0.022・0.031 mm・平均の差 0.021 mm、δ_E 0.051・0.068・0.045 mm、最大は x 14.9〜21.7 r_t (δ_target 5.3〜9.2 mm の 0.3〜0.5 %、上がって下がる揺れ)、x ≥ 3 m は 0.002 mm 以下 (`judge/two_window_check.txt`)。**float (run_0002_ns) との比較** (設定が違うので連結した判定ではない; 精度の効果と追加の 20 万 step の効果は分けられない): 残差の末尾 ρ 1.56e−6 → 3.58e−8・ρUx 6.6e−4 → 7.4e−6・ρE 0.98 → 0.0275・k 6.6e−4 → 9.0e−6・ω 816 → 0.48 (1.6〜3.2 桁下)。float は各列が 158k〜172k で最小の後に上がり、θ_r 0.4 の D は 100k 4.2e−3 → 140k〜160k 1.1e−2 → 182k 7.5e−3 と下げ止まった; double は 40k 1.05e−2 → 200k 6.3e−4 で単調に減少。r 3.3e−3 → 9.1e−5、出口コア M 5.52214 → 5.52045、θ_r (延長の帯) は延長の間に +4.9・+4.2・+3.9 %、26.9 ms/step (約 1.9 倍)。codex 諮問 4 回目 `notes/reviews/2026-10-15-m55-ns-unsettled-4-diagnose.md` (採否は §6.1)。**判断 (ユーザ 2026-10-15「aかな」; 前に「微々たるものじゃね？」)**: 延長しない。受け渡しの関門の明示の上書き 3 つ (`allow_not_reached`・`allow_unconverged`・`allow_unsettled`、いずれも実測できた超過) で δ* の補正 (#3e) に進む。基準は緩めない (`HANDOFF_GATE.json` に理由・実測の値・`settled: false` が残る)。壁の熱流の A/B (#3f) は後回し | O |
| 3e | **δ* の補正の 1 回目 (pass 1)** (ユーザ決定 2026-10-15「aかな」) | 事前に固定する手順と判定: (1) 関門: `deltastar_loop.handoff_gate(run_0003_ns_cont, problem=problem_m55_final.yaml, euler_ref=run_0001_euler/res_12000.h5 (sha256 93cb88b0…), ω 0.5)` に上書き 3 つ (理由はユーザ決定と実測の値) → run の直下に `HANDOFF_GATE.json`・`delta_r_next.csv` (δ_E は最後の出力 res_200000 の 1 枚)。上書きで通らない条件 (局所 y₁⁺・参照・ω・被覆・held・次の壁) は通常どおり (y₁⁺ は #3g の修正の後)。(2) カード: `wt_card card_m55.yaml --out problem_m55_p1_card.yaml --ns-correction run_0003_ns_cont --delta-next run_0003_ns_cont/delta_r_next.csv` (物理スロート径 200 mm で r_t を解き直す; M_knot は自動で選び直す; 補正の表は YAML の横)。(3) 格子: `_mesh/hand_mesh.py build` を #3b と同じ係数・ni 4,800・`axis_cap_frac 0.03` で → `problem_m55_p1.yaml`・`problem_m55_p1_euler.yaml`、品質 (`--ar-max 5000`・AR 緩和の検査)・変換後の第一層の照合。(4) Euler `run_0004_euler` (#3b と同じ段・本段 12,000 step、FP64 のビルド)。(5) NS `run_0005_ns`: `deltastar_loop --delta-r-csv <表> --ic-from run_0003_ns_cont --euler-ref run_0004_euler --euler-ref-res res_12000.h5 --euler-ref-sha256 …`、~~全部 double~~ → **float** (ユーザ決定 2026-10-15「あんまりマッハ数変わらんようだし、いまのとこはfloatでいいんじゃない？20万ステップ上限のルールで、とりあえず作ってもらえたら。フロー回して。」: 精度の切り替えは生成器の既定のまま = run_0002_ns と同じ (thermoFloat 1・implicitSolvePrecision 0・輸送物性の表・面の h は float)、FP64 のビルド (格子の第一層のため); 出口コア M は float 5.52214 と double 5.52045 で −0.03 %。台本 `_runs/run_0005_ns/p1_ns.py --precision float` (既定)、起動の環境変数 `FORGE_TRANSPORT_TABLE`・`FORGE_DIAG_FACE_H_DOUBLE` は無いことを検査)、停止規則 (窓 2 万・保持 4 万・上限 20 万)、出力 2,000。(6) 判定 (回す前に固定): 停止規則の状態・`check_convergence --segment`・停止規則の系列の `check_quasisteady`・壁解像 (`wall`、局所 y₁⁺ ≤ 1、超過の面積 ≤ 2 %)・最後の場の非有限と P ≤ P0・δ の 2 窓 (0.01 mm) → 関門。届かなければ閾値を動かさず、延長か上書きかをユーザに聞く。寸法 (物理スロート径 200 mm・長さ) の最終の合格は成果物の出荷ゲート。(7) 容量: A の空き 5.1 GB → 投入の前に run_0003_ns_cont の途中の全場を消す (res_0・res_140000・160000・180000・200000、壁・出口の出力、履歴、判定の記録は残す: #3f の再生用)。(8) 成果物 (ユーザ「とりあえず作ってもらえたら」): 停止規則の到達・上限 (CENSORED) のどちらでも run_0005_ns の最後の場から標準の報告 (`nozzle_report`、`_after_ns` が自動) と成果物一式 (`python -m forge_design.export.deliverables run_0005_ns --card card_m55.yaml --problem problem_m55_p1.yaml`) を作り、停止の状態・(6) の判定・関門の結果を添えて報告する。壁法線の近壁層の格子の成果物は #5 が未対応なので、作れない部分 (格子の書き出しなど) はそのまま報告する。次の δ* の pass に進むかは δ_E と入れた δ の差を見てユーザが決める **結果 (2026-10-15、判定区間は本段 0〜199,999; 原文 `case/71.m55_combustion_wt/_runs/run_0005_ns/judge/`)**: (1)〜(4) は済み (関門は上書き 3 つで通過・`settled: false`、r_t 0.10024473616136853 m、格子 442,520 点 品質 PASS、Euler `run_0004_euler` 参照 res_12000 sha256 78294389…)。(5) `run_0005_ns` 7,627 s。(6) 停止規則 **CENSORED** (200,000; θ_r の D 6.55e−4・6.02e−4・5.31e−4 > 5e−4、2 万 step ごとに約 0.65 倍で単調に減少、r の平均 1.36e−4)。`check_convergence --segment` NOT CONVERGED (ρ 2.7・ρUx 3.2 falling、ρUy 2.8・ρE 2.6・ω 3.3 flat)。停止規則の系列 DRIFTING (r だけ; θ_r 3 断面 STEADY・単調増加・漸近値は最終の +0.04〜+0.06 %、ṁ_in・ṁ_out・出口コア M STEADY)。壁解像 (`wall`、area、200,000) PASS (y1+ 平均 0.478・最大 0.999、超過 0.0 %)。最後の場: 非有限 0、P0 > Pt (1 + 1e−6) 46,609 点・最大 +0.458 %、入口の列の T0 > Tt + 1 K 9 点。δ の 2 窓 FAIL (measured、δ_E の振れ 0.060 mm @ x 23.4 r_t)。**関門 (上書きなし) は止めた** (事前の規則どおり延長か上書きかはユーザに聞く)。**物理量**: 出口コア M 5.50145 (目標 +0.03 %; pass 0 は 5.52042)、出口面 r ≤ 0.8 r_w の M 平均 5.5006・範囲 5.4914〜5.5098。**試験部の軸に M の波** (M/5.5 − 1 が −0.76〜+0.55 %、極値 x 4.76・5.90・7.02 m、res_100000 と同じ; pass 0 の NS は +0.35〜+0.42 % で平ら、設計壁の Euler は −0.007〜+0.024 %)。pass 1 の物理壁は pass 0 より x 2.0〜3.5 m で局所に約 −0.08° 余分に曲がる (pass 0 の生の抽出の階段が平滑化の後に残った傾きの変化、仮説)。δ_E − δ_in (診断の抽出、run の外 `judge/diag_extract/`) は x ≥ 0 で最大 0.52 mm (出口 −0.39 mm = −0.93 %、x 1〜8 m に波長約 0.8 m・振幅 0.2〜0.25 mm の振動)、pass 0 は 8.1 mm (出口 −15 %)。設計を終える条件 (固定点 0.05 mm) には未達。**解釈と次の pass は codex 諮問中** (`notes/reviews/briefs/2026-10-15-m55-p1-axis-wave.md`)。(8) 成果物 (手元、`case/71.m55_combustion_wt/run_0005_ns/deliverables/`): cad (STEP の転送 PASS)・mesh (構造 v2、STEP からの再構築の照合)・fluent・excel はエラーなし。寸法の関門 FAIL — カードの 5 項目は PASS (スロート径 199.99999984 mm・L_div 8,000.42 mm 等)、整合の `straight_end_in_mesh_interval` だけが落ちた (#3j)。A では FreeCAD・gmsh の venv が無く cad・mesh・fluent が失敗 (想定どおり)。ṁ の帳簿の誤り (#3i) | O |
| 3f | 壁の熱流の定義の A/B (codex 4 回目の提案、**後回し**: ユーザ決定 2026-10-15「aかな」) | `output.interfaceDiag` 0/1 だけを変え、run_0003_ns_cont の res 140k・160k・180k・200k から同じメッシュの index restart で各 1 更新。D_A = E_in − E_out + Q_wall (qwall の台形積分)、D_B = E_in − E_out − 2πΣ iface_Qf_eff (`conjugateWall.cpp:277`、積分済みなので面積や r を掛けない)、あわせて 2πΣ ifaceFw。4 入力の ‖D_B‖₂ ≤ 0.5‖D_A‖₂ なら「収支の残り 1.9 % の主成分は壁熱量の定義の差」を支持、それ以外は棄却。どちらでもソルバの保存性の証明とは扱わない。全圧の超過 (+0.46 %) は未解明のまま | F |
| 3g | ~~`check_wall_resolution.py` が壁でない境界のダンプを壁として読む~~ (2026-10-15 済み) | 生成器は出口に `outputHDFflg: 1` を書く (`runner.py:363`・`runner_wt.py:284`) ので `res_outlet_2_*.h5` があり、既定 (`--groups` 無し) の判定が「判定不能: outlet」になって受け渡しの関門 (d) が上書きできない形で止まった (run_0003_ns_cont)。修正: `--groups` が無いとき bcondConfig.yaml で kind が `wall*` でない境界のダンプを「対象外」と表示して飛ばす (bcond を読めない・名前が無いときは従来どおり評価)。確かめたこと: run_0003_ns_cont で PASS (wall の値は `--groups wall` と同じ)、負例 3 つ (bcond を外す・出口のダンプだけ・`--groups outlet`) はどれも INDETERMINATE のまま; `test_wall_resolution.py`・`test_gate_bad_input.py` PASS | O |
| 3h | ~~関門の「次の壁」と `prepare_ns` が表の端の丸めで止まる~~ (2026-10-15 済み) | #3e (1) の関門 (上書き 3 つ) は (h) 次の壁で止まった: `SingleBSplinePhysicalWall: δ_r の表の範囲 [-10.0091592, 80.073719] が設計スロート 0 から出口 80.073719 までを覆わない`。表の終点 80.07371897756356 r_t が 設計の x_e 80.07371897756357 r_t より 1 ulp 短い (格子の x [m] / r_t の丸め)。NS 後の補正の経路 (`wt_card --ns-correction` の `snap_table_ends`、#18) は端を ε_f32 以内で揃えるが、同じ設計に表を当てる `next_wall_check` と `prepare_ns(delta_r_csv=...)` は揃えていなかった。修正: `deltastar_loop.snap_table_ends_same_design` (端ごとに |Δ| ≤ ε_f32 なら設計の入口端・出口に置き換え、超えた端はそのまま = 例外にせず今までどおり壁の構築が判定) を両方で使い、記録 (`next_wall` の `table_ends`・prepare_info の `delta_r_table_ends`) に残す。確かめたこと: `tests/run_deltastar_root_zone_tests.py` に 1 ulp 短い終点・ε を超える端 (揃えない)・一致する端の 4 件を足して ALL PASS、`run_handoff_diag_mark_tests.py` ALL PASS | O |
| 3i | `metrics.json` (と抽出の要約) の ṁ を一定 x の断面で積分する | `runner_axismach.collect` (`deltastar.massflow_ratio(run_dir, euler_ref)`、1232 行) と `deltastar.py:900` (抽出の要約の `massflow`) が `remap_constant_x` 既定 False のまま。壁法線の層の格子では列が一定 x でないので列の 2π∫ρU_x r dr は断面の流量でない (同関数の docstring の既知の型; ③ベルの経路 `bell_deltastar.py:54` だけ True)。run_0005_ns: 列のまま 120.19 kg/s (散らばり 2.1 %)・比 0.98977、一定 x 121.350 kg/s (散らばり 6.4e−5、停止規則の出口の ṁ 121.350 と一致)・比 0.99930。直し方 (案): 両方の呼び手を True にする (列が一定 x の格子では np.interp が節点の値をそのまま返すので値は不変 — 試験で確かめる)。判断: 2026-10-15 codex diagnose (`notes/reviews/2026-10-15-m55-p1-axis-wave-diagnose.md`) — 判定の緩和でなく積分断面の訂正として採用。両方の呼び手で `remap_constant_x=True` を明示し、列が一定 x の既存格子 (Euler の run) で値が変わらないことを試験で確かめる | O |
| 3j | 成果物の寸法の関門の `straight_end_in_mesh_interval` が滑らかな接合で誤って落ちる | `export/dimension_gate.py` `mesh_dimensions` は入口半径との差 ≤ `NODE_EQ_MM` (1e−6 mm) の節点を直管と数える。run_0005_ns の格子は接合 (x −500.000 mm) の 0.08 mm 下流に節点 (x −499.920 mm、r − 300 = −7.3e−9 mm) があり直管と数えたので、区間 [−499.920, −498.591] が接合を含まず FAIL。カードの 5 項目と弦の検査 (最大 0.00076 mm) は PASS。判断: 2026-10-15 codex diagnose — 最後の 1 点だけの一致で置き換える案は却下 (他の節点の不整合を見逃す)。**接合の位置は STEP の曲線から決め、メッシュの壁の節点の x で STEP を評価し、直管側から接合を越える区間全体の半径を照合する** (許容は `NODE_EQ_MM`)。カードの寸法・弦の条件は維持し、接合をずらした負例 (接合を ±数節点動かした曲線) を拒否することを試験に入れる | O |
| 3k | **抽出の径方向の標本化の A/B (forge 0 step、事前登録 2026-10-15、codex diagnose の判別)** | 変えるのは `band_local_deficit` の `refine` だけ: A = 4 (今)、B = 32。対象は pass 1 の壁の元 `run_0003_ns_cont/res_200000.h5`、参照 `run_0001_euler/res_12000.h5` (sha256 93cb88b0…)。固定: 一定 x への補間 (`remap_constant_x`、`remap_edge` missing)、各列の保存済みの `band_y_b` (`y_b_fixed`、帯の梯子は選び直さない)、Euler の補間、平滑化のノット 2 r_t・λ 1・重み。手順: まず A が保存済みの raw (`_runs/run_0005_ns/judge/diag_wave/pass0_run_0003_delta_r_equiv.csv`) を再現することを確かめる。判定 (x 1〜3.5 m): raw の段差、平滑化後の δ_E、補正の角 atan(r′_design + δ′_E) の差。**許容 max|Δδ_E| = 0.01 mm、max|Δθ| = 0.01°**。結果 A (どちらかを超える) → 「今の抽出は径方向の標本化に十分独立」を棄却し抽出を直す。結果 B (両方以内) → 固定帯の標本化を壁角の窪みの主要因から外す (帯の選び方の感度は除外しない)。どちらでも軸の波の原因の確定にはしない。AWS A (CPU、保存元)。やらない: 物理壁の Euler 1 本で形状/境界層を二分する・抽出の検査の前に λ を上げる・ω だけ変えて pass 2 を投入する | O |
| 3l | 軸の波の時系列の準定常 | codex diagnose の Major (「過渡ではない」は 2 枚の比較では未立証)。`run_0005_ns` の軸の M の最大絶対偏差・半振幅・極値の位置を各全場から作り `check_quasisteady --series-csv` に渡す。**A の途中の全場は判定の後に消した**ので、残る 5 枚 (100k・140k・160k・180k・200k) で行い、点数が足りなければ「判定不能」と書く | O |
| 3m | 次の壁の作り方の選択 (ユーザ判断待ち) | (a) 直接の δ* (今の経路) を #3k の結果で抽出を直してから pass 2、(b) CONTUR の積分式を NS の δ* で較正した壁 (ユーザ 2026-10-15 の問い「conturの積分式を補正して壁作ってるか？」; チェーンには `deltastar_initializer` の k_f の較正がある、M6 は 10-05 に出口だけの較正の壁が滑らかさ・波・オーバーシュートの基準を通った)、(c) pass 0 または pass 1 を暫定の成果物で止める。試験部の合格条件 (codex 提案: x 4.0〜7.95 m の最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %) もユーザに聞いてから事前登録する **ユーザ 2026-10-15**: 「いつもこうしてなかったか？ 摩擦係数でなくて、倍率かけて合わせ込むほうが、割り切っててよい気はしてます」。事実: M6 の生産 (case/45、10-07) は C2 方式 (CONTUR の摩擦係数倍率 k_f = `deltastar_initializer.cf_scale` を出口の δ_C = δ_E に合わせ、r_t を解き直す; `procedures/nozzle-design-workflow.md` 4.) で、NS の δ_E を直接壁に入れていない。case/71 の pass 1 はカードの `--ns-correction` (直接の δ_E の表、10-13 に足した経路) を使った = M6 の作り方から外れた。2026-09-04 の生産レシピの記録にも「反復ごとの壁変化が特性線で末端の軸に届き ±0.5 % の軸の波を作る」とあった。pass 0 の比 δ_E/δ_CONTUR (cf_scale 1): x 1 m 0.47・2 m 0.67・3 m 0.71・4 m 0.76・6 m 0.80・8 m 0.85、スロート CONTUR +0.09 mm 対 NS −0.24 mm (符号が逆)。案 (b′): CONTUR の δ(x) に一様な倍率 c を掛けて出口で合わせる (コードに出力の倍率は無い — `deltastar_initializer` にキーを足す)。§4 の方式の変更なので、決まったら §4 を書いて codex に諮ってから実装 **ユーザ 2026-10-15「コンターの倍率だと、そのマイナスの厚み再現できないよな」**: そのとおり (正の倍率では CONTUR のスロートの +0.09 mm を NS の −0.24 mm にできない)。出口の δ で倍率を決めると、スロートの差 0.32 mm (半径) が実効のスロート面積で約 0.64 % → 出口 M で約 −0.1 % (面積比の感度 ΔA/A ≈ 5〜6 ΔM/M、γ 1.3 前後の概算) になる。案: (b″) 倍率 c + 全長で一定のずらし δ₀ (2 係数、ずらしは壁の傾きを変えないので波を作らない)、pass 0 の NS のスロートと出口に合わせると c ≈ 0.852・δ₀ ≈ −0.32 mm (概算); (b‴) 倍率 1 つを出口の δ でなく出口 M に合わせる (スロートの差を出口で吸収)。どれを採るかはユーザ判断 **判断: 2026-10-15 ユーザ「1かええな」 — 案 (1) 倍率 c + 一定のずらし δ₀ (§4.3)**。§4.3 を codex に諮ってから #3n・#3o へ | F |
| 3n | §4.3 の実装 (`deltastar_initializer.delta_scale`・`delta_offset_m`) | `runner_axismach.integral_delta_r` で平滑化の後に δ = c·δ_C + δ₀/r_t (r_t 単位) を当て、`init_info` に c・δ₀・当てる前のスロート・出口の値を記録。既定 (1, 0) はビット同一 (試験)。カードの `overrides` で渡せること、寸法の逆算・`prepare_ns`・`hand_mesh.walls` が同じ値を使うこと (スロート径 200 mm に解き直した r_t で物理スロートが 200 mm) を試験で確かめる。**§4.3 の codex 諮問の後** | O |
| 3o | §4.3 の壁で M5.5 を作り直す (Euler → NS → 成果物) | カード (`overrides` に c・δ₀) → `hand_mesh.py build` (#3e と同じ係数) → Euler (使い回しは諮問で決める) → NS `--init-integral`、初期場は run_0005_ns の最後の場から補間、`--stages full`、精度は既定 (float)、出力 2,000、停止規則 (窓 2 万・保持 4 万・上限 20 万) → §6 V5 の判定 → 成果物。**V5 の許容のユーザ確認と §4.3 の諮問の後** | O |
| 5 | ~~**成果物の対応 (M2)**~~ (2026-10-15 済み) (ユーザ 2026-10-15「成果物ツールを対応させて」) | `export/mesh_structure.py`・`deliverables.py` を壁法線の近壁層の格子でも作れるように (既定化の前提)。**通らないのは格子の構造ファイルと STEP からの作り直しだけ** (`mesh_structure.extract_structure` の縦線の検査で StructureError → `part_mesh` が「作らない」と書く; Fluent・CGNS・Excel・STEP は msh・壁ファイルを直接読むので影響なし)。**仕様**: 節点は `mesh2d.py:360-387` の式どおり (i,j) = P_i + (1 − β_ij)·d_ij·(0, −1) + β_ij·d_ij·n_i (P_i = (x_i, r_w,i)、n_i = (r′_i, −1)/√(1 + r′_i²)、d_ij = r_w,i (1 − s_ij))。(1) `extract_structure` に壁の傾き r′_i (station ごと、無次元) を渡せるようにし、縦線の格子は今の v1 の出力をビット同一に保つ。縦線でなく傾きがあれば `forge.mesh_structure/v2` (節点の規則 `wall_normal_blend_v1`、v1 の項目 + `wall_slope`・`s` (ni × nj)・`beta` (ni × nj)・傾きの出典と sha256) を作る。x_i は壁の節点の x (軸の節点の x と一致を検査)、β・d は dX = X − x_i・dR = r_w − R から解く (|r′| ≤ 1e−12 の station は β = 0)。作った後に座標を作り直して元の節点との差を自己検査 (最大 ≤ 1e−7 mm かつ station の第一層厚の 1e−3 以下、外れたら StructureError)。傾きが無い非縦線の格子は今までどおり StructureError (理由に「壁の傾きを渡す」)。(2) `deliverables.part_mesh` は run の壁ファイル (`load_wall_file(run)` の physical の r(x, 1)、x は r_t 単位 = x_m / scale_m) から r′_i を作って渡す (読めなければ今までどおり「作らない」と理由)。(3) `rebuild_from_step.py` (numpy・gmsh だけ) は v1・v2 を読み、v2 では各 station の u_i で STEP の曲線の導関数から r′ = (dr/du)/(dx/du) を出して同じ式で節点を置く。入口・出口は OCC の直線のまま (端の列の節点が直線から外れる量を比較の JSON に記録)。比較の基準は V3 のまま (距離 ≤ 1e−5 mm かつ第一層厚の 0.01 以下・接続・物理グループ・forge の形)。**合格条件**: `tests/run_mesh_structure_tests.py`・`run_deliverables_tests.py` が通る (v1 はビット同一)、mesh2d の `wall_normal_layer` の小さな格子 (端の断面の戻し `wall_normal_end_stations` の有り無し) で v2 の作り直しが 1e−12 L 以内、傾きを 1e−6 ずらすと自己検査が落ちる、実データ (case/71 run_0003_ns_cont の nozzle.msh・wall_repr.json) で v2 の抽出と STEP からの作り直しの比較が PASS **実装 (2026-10-15、implementer、未 commit)**: v2 (`wall_normal_blend_v1`) の抽出・自己検査・`rebuild_from_step.py` の v2・`part_mesh` の傾き。v1 は 6d4334e3 とバイト一致 (`run_mesh_structure_tests.py` ok 66・`run_deliverables_tests.py` ok 8)。実データ (run_0003_ns_cont の格子を prepare_info と壁ファイルから再生成、msh の sha256 が記録と一致) で 自己検査 2.7e−12 mm、STEP からの作り直しの比較 PASS (最大 1.07e−9 mm・第一層厚比 2.2e−6・出口の端の列の直線からの外れ 0.62 mm = #3a の出口面の非平面)。**判断 (2026-10-15、主セッション)**: 自己検査は渡した傾きで (β, s) を解くので、傾きが r′ と同じ符号にずれる誤り (定数倍・単位) を見つけられない (×(1 + 1e−6) が通る; 検出は β を [0, 1] に収める片側だけ)。→ mesh2d の規則のパラメータ (`wall_normal_layer` の d_n・d_b、`wall_normal_end_stations`; run の prepare_info の格子の記録) を v2 に記録し、規則の β(d) で節点を置き直して元と照合する両側の検査を自己検査に足す (パラメータが無ければ「規則の照合は未実施」と記録し STEP からの作り直しの比較だけに頼る)。[0, 1] の制約は残す。v2 は 17 桁の msh が前提 (10 桁は丸めで StructureError; `wall_normal_layer` を雛形の既定にするなら `msh_digits 17` も既定に) **規則の照合の追加 (2026-10-15、implementer)**: v2 に `rule_params` (run の prepare_info.json の `mesh.params` の `wall_normal_layer`・`wall_normal_end_stations`、出典の sha256 と記録の msh の sha256 の一致) を書き、規則の β(d) と渡した傾きで置き直した節点を元と照合する (`checks.rule_check`、許容は自己検査と同じ)。パラメータが無ければ「未実施」と書いて通す。確かめたこと (主セッションで回し直し): `run_mesh_structure_tests.py` FAIL 0 (ok 78、v1 は 6d4334e3 とバイト一致、合成格子で傾き ×(1 + 1e−6)・d_n 1 %・端の断面の数・d_n > d_b が StructureError)、`run_deliverables_tests.py` FAIL 0 (ok 13、FP64 の変換器; 実データ run_0003_ns_cont の格子で規則の照合 6.4e−12 mm、傾き ×(1 + 1e−6) は 2.2e−5 mm・d_n 1 % は 1.1e−2 mm で StructureError、STEP からの作り直し PASS)、implementer で `run_fluent_msh_tests.py`・`run_deliverables_walltemp_tests.py` FAIL 0。残り: 実データの fixture は scratch にある (`make_deliverables_fixture.py` で作り直せる)、`deliverables.GMSH_PYTHON` と FreeCAD は手元の固定パスなので成果物は手元で作る | O |
| 2 | 実装 | `estimate_mesh` の `local_v1`、雛形 (`axismach_explore_v1`・`axismach_final_v1`) の `mesh.method`・`first_layer_calibration`・`wall_normal_layer`、記録、試験 | O |
| 3 | M5.5 のカードの生成 (AWS A) と Euler → NS | ユーザ了承の節点数で | O |
| 4 | skill の手順 | `nozzle-design`・`design-intake` に格子の作り方と NS の後の壁解像の確認 | O |

## 6. 検証

事前登録 (2026-10-15、上位に諮る前の案):

- **V1 見積もりの再現**: 係数 1.0 で M5.5 の ni が 3,163 (±10)、節点 306,811 前後になる (見積もりの台本と同じ入力)。
- **V2 格子の品質**: M5.5 と M5 デモ (case/70 の c2b のカード) の `local_v1` の格子で `check_mesh_quality.py --bcond --ar-max 5000` が PASS、隣り合う dx の比 ≤ 1.1、壁直交性 (`wall_normal_layer`) の判定 PASS。
- **V3 第一層**: 表の値 ≤ c(x)·y₁,req(x) (全断面)、見積もりの y₁⁺ の最大 ≤ 0.8。
- **V4 NS の後** (M5.5 の NS、停止規則まで): 壁ごとの y₁⁺ > 1 の面積 ≤ 2 % (`--groups <全 no-slip 壁> --weight area --step <停止の step>`)。超えたら §4.2 の手順で表を直す (この plan の不合格ではなく、係数の見直しの材料として記録)。
- **V5 M5.5 の 2 係数の壁 (§4.3) の NS** (事前登録 2026-10-15; 数値の許容は**ユーザの確認待ち**、確認の前に結果を見ない): 判定区間は本段 (`--segment`)。
  - V5a 寸法: 成果物の寸法の関門のカードの 5 項目 PASS (物理スロート径 200 ± 0.02 mm ほか)。
  - V5b 出口コア M (停止規則の量): 5.5 ± 0.05 % (案: M6 の生産のゲートと同じ幅)。外れたら §4.3 の再較正を 1 回。
  - V5c 試験部の一様性 (案: codex diagnose 2026-10-15 の提案): x 4.0〜7.95 m の M の最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %。軸の節点は既知の値の癖 (静圧 0.5〜0.7 % 低) があるので r/r_w = 0.1 の線で測り、軸の節点の値は併記する。
  - V5d 流れの判定: 停止規則・`check_convergence --segment`・停止規則の系列の `check_quasisteady`・壁解像 (y₁⁺ > 1 の面積 ≤ 2 %)・非有限 — #3e (6) と同じ。届かなければ閾値を動かさずユーザに聞く。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-15 | [2026-10-15-tooling-nozzle-card-local-mesh-plan.md](../../notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan.md) | NO-GO (分類 incomplete: 見積もりの台本・出力と M6 の比の元の記録がこの作業ツリーに無い), C0/M7/m1 | **全件採用、進め方を 2 段に分ける**: M5.5 は先にケースだけの手組みの格子 (M6 の H と同じやり方; #3a) で進め、生成器の既定の変更 (#2) は指摘を直して再レビューの後。M1 (実寸の下包絡を半径で割ると補間で要求を超える: 反例 1.0607 倍) → f_req = c·y₁,req/r_w の下包絡を表にし、補間の区間の中も検査、物理壁でも実寸を再検査。M2 (`wall_normal_layer` で成果物 `mesh_structure.json`・STEP の再構築が作られない) → 既定化の前に成果物の対応 (#5)。M3 (補正の係数の根拠: H の表は冷却倍率を掛けた見込み、縮流部の比 0.61〜2.2 に係数 1.0 は保守側でない) → 「初期格子の経験の係数」と位置付け、元の run・step・比の分布を固定して commit、縮流部も比の最小 (0.6) に、補間後の c(x) が比を超えない検査。M4 (V2 の `--bcond` がこの枝の検査器に無い・AR 緩和の条件が未定量) → `runner.mesh_ar_relax_check` (AR > 1000 の全セルで skew ≤ 0.1) と変換後の第一層厚の照合・`msh_digits 17` を V2 に。M5 (V4 を記録だけで閉じられる) → 最終の完了は再調整の後の V4 PASS、全 no-slip 壁の同じ step の出力の存在、`check_convergence --segment`・壁量の系列の `check_quasisteady` を事前登録。M6 (見積もりと実生成の設定の契約・旧い overrides) → 見積もり・節点数・品質・実生成で同じパラメータの組み立て、雛形の旧キーの除去、明示の overrides の衝突は拒否、`mesh.ar_max` に従う。M7 (Euler の ni も変わる) → Euler の ni が NS と結合していることを明記し、M5 で出口 M・流量・δ_r の基準格子との差と準定常を受入に。m8 (参照の解決) → 再現の入力を case/71 に commit、関連 plan は統合の後に参照 |
| diagnose | 2026-10-15 | [2026-10-15-m55-ns-unsettled-4-diagnose.md](../../notes/reviews/2026-10-15-m55-ns-unsettled-4-diagnose.md) (ブリーフ [briefs/2026-10-15-m55-ns-unsettled-4.md](../../notes/reviews/briefs/2026-10-15-m55-ns-unsettled-4.md)) | 結論: 長い延長と δ* の補正を保留し、壁の熱流の A/B を先に。Major 4・Minor 1 | (M1) すぐ約 6 万 step 延長する → **採用** (延長しない; ユーザ決定「aかな」)。(M2) θ_r・出口 M の STEADY を理由に δ* へ進む → **採用** (STEADY は根拠にしない): 2 窓の判定を完了させた (FAIL、measured; #3d) うえで、ユーザ決定で関門の明示の上書き 3 つで進む (#3e、基準は緩めない)。(M3) 収支の残り 1.9 % を非保存と解釈する → **採用** (非保存の証拠と書かない); A/B は #3f に入れ、ユーザ決定で後回し。(M4) double の延長で thermoFloat の寄与を除外する → **採用** (「double でも全圧の超過は残る」だけを書き、旧経路の寄与は除外しない)。(m) 終了状態が未確定 → 03:15:22 UTC に CENSORED で確定 (対応済み)。報告の言い方 (予算の上限で終了・停止条件に未到達) は採用 |
| diagnose | 2026-10-15 | [2026-10-15-m55-p1-axis-wave-diagnose.md](../../notes/reviews/2026-10-15-m55-p1-axis-wave-diagnose.md) (ブリーフ [briefs/2026-10-15-m55-p1-axis-wave.md](../../notes/reviews/briefs/2026-10-15-m55-p1-axis-wave.md)) | 結論: pass 2 は保留し、pass 0 の保存場・Euler 参照・帯を固定して抽出の径方向細分 `refine` 4/32 だけを変える forge 0 step の A/B を先に。Major 8・Minor 1 | (H1 要再検証) → **採用**: #3k (refine A/B)。(H2 却下: 生産の抽出は `method="band"` で core_rms は診断値) → **採用** (core_rms で重み 0 にする案 (e) は採らない)。(pass 2 の補正方法: (a)〜(c) 単独は推奨しない、カード経路は delta_target でなく delta_E を読むので ω は効かない) → **採用** (pass 2 は #3k の後、緩和を入れるなら寸法の逆算と壁の生成が同じ表を使う設計が要る)。(成果物: pass 1 は途中の成果、合格品ではない; 試験部の一様性を優先するなら暫定の候補は pass 0 = 軸の最大偏差 0.418 % 対 0.764 %) → **採用** (ユーザに両方を示して選んでもらう)。(E(i) ṁ の 2 つの呼び手) → **採用**: #3i。(E(ii) 直管終端の単一点の例外) → **却下を採用**: #3j を STEP の曲線から接合を決めて区間全体を照合する方式に変更。(「過渡ではない」は未立証) → **採用**: #3l (軸の波の系列の `check_quasisteady`; ただし A の途中の全場は判定の後に消したので残る 5 枚 100k・140k・160k・180k・200k で行う)。(H3 の因果は未分離) → **採用** (出口コア M の改善を全体の補正量だけの効果と書かない)。(m: 「±0.6 %」は不正確) → **採用** (最大絶対偏差 0.764 %・半振幅 0.656 % と分けて書く)。(試験部の合格条件の提案: x 4.0〜7.95 m の最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %) → **保留** (ユーザ指定値ではないので、事前登録の前にユーザに聞く) |

## 7. 影響範囲

- `design/forge_design/intake/wt_card.py`・雛形 2 本・試験。既存のカード (`mesh.method` を書いていないもの) は新しい既定で格子が変わる (ユーザの方針「標準と決めた方式はコードの既定にする」)。
- case/71 の M5.5 の生成。

## 8. 完了条件

- [ ] §6 V1〜V3 を満たし、V4 を記録
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status: done`、§9 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-15` — ユーザ決定で M5.5 の δ を CONTUR × 倍率 c + 一定のずらし δ₀ に (§4.3、#3m〜#3o)。pass 0 から c 0.852073・δ₀ −0.32201 mm。§6 V5 を事前登録 (許容はユーザ確認待ち)。codex diagnose (軸の波) の採否を §6.1 に。
- `2026-10-15` — δ* の補正の 1 回目 (#3e) の NS `run_0005_ns` を判定: CENSORED・NOT CONVERGED・δ の 2 窓 FAIL → 関門は止めた。出口コア M 5.50145 (pass 0 は 5.52042)、試験部の軸に M の波 −0.76〜+0.55 % (新しく出た)。δ_E − δ_in 最大 0.52 mm (pass 0 は 8.1 mm)。成果物一式を手元で作成 (寸法の関門は整合の 1 項目で FAIL、#3j)。ṁ の帳簿の誤りを #3i に。軸の波の解釈と次の pass を codex に諮問中。
- `2026-10-15` — 成果物の構造ファイル v2 (`wall_normal_blend_v1`) と STEP からの作り直しを壁法線の近壁層の格子に対応 (#5)。規則のパラメータで節点を置き直す照合で傾きの誤りを両側で検出。v1 は 6d4334e3 とバイト一致。
- `2026-10-15` — 関門の次の壁・`prepare_ns` で表の端を ε_f32 以内なら設計の端に揃える (#3h; 1 ulp 短い表で 1 本の B-spline の壁が止まった)。
- `2026-10-15` — run_0003_ns_cont (全部 double の延長) を判定: CENSORED・NOT CONVERGED・δ の 2 窓 FAIL (measured、δ の 0.3〜0.5 %)。ユーザ決定「a」で関門の明示の上書きで δ* の補正 (#3e) へ。codex 4 回目の採否を §6.1 に。`check_wall_resolution.py` の壁でない境界のダンプの扱いを修正 (#3g)。
- `2026-10-15` — plan 段レビュー NO-GO (C0/M7/m1) を全件採用。M5.5 は手組みの格子 (#3a) で先行、生成器の既定化は書き直して再レビュー。
- `2026-10-15` — 起票 (ユーザ「古いやり方でメッシュ作ってない？」、決定「a」)。見積もり (3 通りの作り方の ni・節点) と見積もりの式の偏り (M6 の H の表との比) を §3 に。
```

## 参考: `procedures/nozzle-design-workflow.md`

```
# 軸対称風洞ノズル (axismach) の設計・計算手順

問題タイプ `wind_tunnel_axisym_axismach` (軸 Mach 則 → 逆 MOC で壁を作る風洞ノズル) を対象に、
よくある 3 つのユースケースの作業手順をまとめる。最新の実例は `case/45.isobutane_m6_d155`
(M6、出口直径 1.55 m、イソブタン燃焼ガス、2026-10 時点の生産レシピ)。

- 関連: 要件インタビュー [`/design-intake`](../.claude/skills/design-intake/SKILL.md)、
  標準出力と pptx [`nozzle-design-outputs.md`](nozzle-design-outputs.md)、
  計算の一般手順 [`calculation-workflow.md`](calculation-workflow.md)、
  発散対処 [`divergence-and-startup.md`](divergence-and-startup.md)、AWS 投入 [`/forge-aws-run`](../.claude/skills/forge-aws-run/SKILL.md)。
- 設計判断の正本: [`plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md`](../plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md)
  (CFD ピン・当てはめ壁・細分格子・Euler 較正、§5.1 #11〜#11f)、
  [`plans/active/tooling-design-problem-campaign-recipe.md`](../plans/active/tooling-design-problem-campaign-recipe.md)
  (再現可能なパイプライン化、§4.9 C2 方式 — **未実装部分が多い**)。
- コマンドは `design/` 直下で `.venv-opt/bin/python` を使う。run は検証 run でも AWS で回す。
- **計算環境 (2026-10-07)**: 積分法の δ_r (scipy の RK45) は numpy・scipy の版で 1e-6 r_t の桁 (case/45 の出口付近で約 0.5 µm) 動く。
  設計として問題になる大きさではないが、壁・メッシュの**ビット同一を要する照合** (保存した壁の作り直し・W3/W4 型の入力の同一性・生産の
  問題の準備が生産の run の入力を再現するかの確認) は、**照合する相手と同じ環境で**行う。
  - 生産の準備と照合は、生産の run を回す AWS の環境で行う (例: `case/45.isobutane_m6_d155/prod_confirm.py`)。
  - 手元では `design/.venv-opt` の Python を使う。システムの `python3` は numpy・scipy の版が違い (numpy 1.26 には `numpy.trapezoid` も無い)、
    手元の過去の成果物とも一致しない。
  - `prepare`・`prepare_ns` は `prepare_info.json` の `environment` (Python・numpy・scipy の版、実行ファイル、ホスト) に環境を残す。
    照合が合わないときは、まず両方の `environment` を比べる。

## 0. 共通: problem YAML (契約書)

すべての出発点は `case/NN/problem_*.yaml`。読み込みと検査は
`design/forge_design/probdef.py` の `load_problem()` (未知の type・欠けた spec・組成の不整合で例外)。

| ブロック | 主なキー | 注意 |
| --- | --- | --- |
| `spec` | `Pt`, `Tt`, `p_ambient`, `r_throat` [m], `M_design`、任意で `sizing` | **出口径は直接書けない** (`D_e` は未実装)。出口半径は `r_throat` × MOC の面積比 + δ で決まる。物理出口径に合わせるには `solve_rt`、物理スロート径から決めるには `solve_rt_throat` (§1 ④)。どちらで決めたかを `sizing: {method: exit \| throat, target_m: …}` に書くと `prepare_info.json` の `sizing` に目標との差が出る |
| `gas` | `model: semiperfect`, `species` (組成), `transport` | Tt > 600 K・燃焼ガスは semiperfect。条件を変えたら `geometry.r_inlet` も付け替える |
| `dv` | `L_c.value` | 軸 law の長さ。runner は `min/max` を使わない |
| `geometry` | `R` (スロート曲率半径/r_t)、`r_inlet` (入口配管半径/**r_t 単位**)、`L_U` (縮流部長)、`L_pipe`、`axis_law: knot`、`M_knot`、`Lc_mode`、`start_line: throat_char`、`wall_mode: cplus`、`n_axis_inv` | `Lc_mode`: `explicit` (L_c を使う、許容窓の外は例外) / `max` / `from_length` (`dv.L_total` から L_c を解く)。r_t を変えても r_inlet・L_total は r_t 単位のまま (実寸がずれる) **要求カードでは M_knot は既定で自動** (最大壁角が最小・成立の境目から 0.05 以上、2026-10-14; plan `tooling-nozzle-mknot-auto`)。 |
| `geometry` (生産レシピ) | `initial_line: cfd` + `initial_line_run`/`initial_line_res`、`wall_repr: joint`、`Md_moc_offset`、`pw_upstream: poly` | §2 ⑤ 参照。新規設計の初回は `initial_line: hall`・`wall_repr: interp` (既定) でよい。`pw_ramp` は旧設定 (`pw_upstream: ramp` と組で使う。`poly` と併記すると例外) |
| `mesh` | `ni`, `nj`, `wall_first_frac` (+ `_throat`・ブレンド)、`throat_refine`、`throat_width`、`axis_cap_frac`、`ar_max` | **NS (`prepare_ns`) の格子**。壁解像 PASS の格子を使う (§2 ⑥)。Euler には効かない (2026-10-07 から) |
| `mesh_euler` | `mesh` と同じキー (未知のキーは例外) | **Euler (`prepare`) の格子** (2026-10-07、plan `verification-case45-euler-total-enthalpy` §4)。`mesh` からは補完しない。既定は全断面で `wall_first_frac` 0.005 の等比 (スロートの別指定なし・軸側の cap なし、ni × nj は 321 × 65)。**`mesh_euler` の無い問題は Euler の `prepare` が移行先を示して止まる** (旧格子の記録を再現するなら `mesh` の全キーを写す)。case/45 は ni 2000 × nj 97・軸方向は G1 と同じ (`throat_refine 4`・`throat_width 3`)・全域 0.005 |
| `evaluate` | `nStepOuter` (**Euler の本段だけ**)、`ns_steps` (NS の本段、既定 60,000)、`outStepInterval`, `cfl_main`, `tp_species`, `condensation` | 凝縮は `tp_species: split_h2o` + `condensation`。**NS の本段の step 数は `--steps` > `evaluate.ns_steps` > 既定 60,000** (2026-10-14 ユーザ決定。以前は雛形の `nStepOuter` 12,000 が NS にも効いて短すぎた)。`prepare_info.json` の `nStepOuter_source` に出所。**NS の本段は既定で停止規則つき・上限 20 万 step** (2026-10-14 ユーザ決定; `evaluate.stop_rule`・`stop_rule_cap`、記録は run の `stop_rule.json`、plan `tooling-nozzle-cooled-wall-ns-settling` §4.6) |
| `mesh`・`mesh_euler` の上限 | `max_nodes`, `max_nodes_reason` (組で書く) | **2D の節点数の上限 30 万点** (2026-10-14、plan `tooling-nozzle-mesh-size-guard`)。生成器・`prepare`・`prepare_ns` が格子を作る前に実効の節点数 (`axis_gap_frac` の導出 nj を含む) で判定して止まる。引き上げはユーザの了承を取ってから理由と組で書く (AGENTS.md「メッシュの規模は作る前に示す」)。`prepare_info.json` の mesh 欄に `n_nodes` |

新しい問題は既存の YAML を複製して作る。複製元のコメントに残る古い設定 (`M_design` の試作較正値、
`mesh.bndFirstOrder`、廃止キー) を持ち込まない ([`recommended-settings.md`](recommended-settings.md) §9)。
検査: `cd design && .venv-opt/bin/python -c "from forge_design.probdef import load_problem; load_problem('../case/NN/problem.yaml')"`。

## 1. 新規の非粘性ノズルを設計する (軸長・出口径・入口径・入口条件を指定)

1. **要件を決める** — skill `design-intake` (要件インタビュー → 検証済み YAML・方針メモ。本書は以降の実行手順を担う) で、Pt・Tt・組成・設計マッハ数 `M_design`・入口配管半径・軸長
   (または全長)・出口径を確定する。凝縮の可能性 (H₂O を含み膨張で低温) をここで確認する。
   **実寸で与えられた要求** (物理スロート径・入口半径・直管長・縮流部長・設計スロート→出口の長さ・壁温) は要求カード
   (`forge.wt_card/v1`) に書き、`-m forge_design.intake.wt_card CARD --out problem.yaml` で ② の YAML を作る
   (r_t を解いた後の長さの割り直しを生成器が行う。④ で r_t を解き直したら `--r-throat` で作り直す)。
2. **YAML を作る** — 近い既存ケースの `problem_*.yaml` を複製し、`spec`・`gas`・`geometry.r_inlet`
   (= 入口半径 / r_t)・`dv.L_c` (または `Lc_mode: from_length` + `dv.L_total`) を書き換える。
   `r_throat` は仮の値でよい (④ で出口径から解く)。
3. **MOC だけで壁を作る** — `design_chain(p)` (`design/forge_design/evaluate/runner_axismach.py`) が
   Hall のスロート初期線 → 軸 Mach law → 逆 MOC → 壁の QA を決定的に実行する (CFD なし)。
   ファイルに落とすには Euler run の準備だけを行う:
   ```
   .venv-opt/bin/python -m forge_design.evaluate.runner_axismach ../case/NN/problem.yaml ../case/NN/run_XXXX_moc --prepare-only
   ```
   → `wall_design.csv` (x, r, 壁角, 壁 M)・`target_axis_M.csv`・`nozzle.msh`/`nozzle.h5`・`MESH_QUALITY.txt`。
   設計が成立しないと例外になる (`law.gates()`・`wall_qa`: 半径の単調性・出口角 ≤ 0.2°・r_F の 1D 比 ±3 %・
   スロート曲率、`wall.validate()`: スプラインのうねり・C¹/C² 接続、縮流部 μ ≤ 20)。
   **fold・特性線トポロジ・壁マージン (μ_w − θ_w ≥ 1°) は `design_chain` に入っていない** —
   `design/forge_design/geometry/moc_diagnostics.py` で別に確認する (§2 ①)。
4. **寸法 (r_t) を決める** — 設計は r_t 単位なので、寸法は r_t の 1 変数。**出口径から**決めるか、**物理スロート径から**決めるかを選ぶ
   (同時に両方を固定条件にしない)。解いた r_t で `spec.r_throat` を手で書き換え、`spec.sizing` とコメントにどちらで決めたかを残す
   (plan `tooling-nozzle-upstream-poly-and-throat-sizing` §4.2):
   ```
   # 出口径から (物理出口半径 0.775 m)
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ../case/NN/problem.yaml --euler-ref X --run-dir Y --solve-rt 0.775
   # 物理スロート径から (物理スロート半径 [m])
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ../case/NN/problem.yaml --euler-ref X --run-dir Y --solve-rt-throat 0.07676
   ```
   ```yaml
   spec:
     r_throat: 0.0766539     # [m] solve_rt_throat (物理スロート半径 0.07676 m) で決めた (YYYY-MM-DD)
     sizing: {method: throat, target_m: 0.07676}
   ```
   (`--euler-ref`/`--run-dir` は必須引数だがこのモードでは使われない。) CFD 前の δ は `prepare_ns` と同じ経路 (YAML の
   `deltastar_initializer` の k_f・熱条件・平滑化、2026-10-07 から。以前の `solve_rt` は未較正・未平滑化で、生産の壁と出口半径が 2.4 mm
   ずれた) で、物理壁も `prepare_ns` と同じ構築 (`pw_upstream` に従う)。**k_f を較正し直したら寸法も解き直す** (k_f はスロートの δ も変える)。
   **CFD 前の寸法は初期見積もり**: 残差はスロート 1e-7 m・出口 1e-5 m (`SIZING_TOL_PRE_CFD_THROAT_M`・`SIZING_TOL_PRE_CFD_EXIT_M`、
   引数 `tol_R_m` で明示も可) で止め、30 回で収まらなければ不合格 (例外)。CFD 前の δ (積分法) は r_t に対して滑らかでなく
   (r_t を 1e-12 m 変えるだけで出口の δ が ~1e-5 r_t 揺れる)、局所の散らばりはスロート 7.4e-9 m・出口 2.6e-6 m (2026-10-07 U2c) で、
   精密な寸法・感度の評価には未検証 (plan `tooling-nozzle-upstream-poly-and-throat-sizing` §6 U2d・既知の制約)。NS 後の `solve_rt_throat` は
   `SIZING_TOL_M` (1e-9 m) で止め、NS 後の `solve_rt` (`--prev`) は固定回数の反復で許容差を判定しない。
   NS 後に物理スロート径から決め直すときは `--prev RUN --delta-next RUN/_extract_edge/delta_r_next.csv --delta-r-out FILE`
   (δ_E の全分布で寸法を決め、同じ補正の表を FILE に書く。次の壁は `prepare_ns(delta_r_csv=FILE)` でその表から作る)。
   **実寸の長さを保って作り直すとき (要求カードから; 2026-10-13、plan `tooling-nozzle-deltastar-band-root-zone`)**: r_t が動くと r_t 単位の設計を
   据え置く上の経路では実寸の長さが公差を外れる (case/70: r_t ×1.0041 で拡大部長 +4.4 mm)。カードから作り直す:
   ```
   .venv-opt/bin/python -m forge_design.intake.wt_card CARD --out P.yaml --ns-correction PREV_RUN --delta-next PREV_RUN/<抽出>/delta_r_next.csv
       → P.yaml・P_delta_r.csv (最終の r_t での補正の表)・P.yaml.card_log.json
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem P.yaml --euler-ref EULER_RUN --delta-r-csv P_delta_r.csv \
       --ic-from PREV_RUN --run-dir RUN [--stages none --cfl 1 --steps 60000]
   ```
   `--prev` / `--ns-correction` の前の run は `wall_repr.json` (前の設計の区切りの点) が必須。表の端は前の設計の入口端・出口に float32 の丸めの範囲で揃え、
   区間で正規化した座標 s で今の設計へ写す。使用域 [0, x_e] に hard 不合格 (held) の点があれば止まる。`--delta-r-csv` は表の sha256 をカードの記録と照合する。
   カードの寸法ゲートは**設計壁から測った値**で、最終の合格は成果物 (STEP・メッシュ) の出荷寸法ゲート。
   **冷却壁の注意**: case/70 (T_w/T_aw ≈ 0.3) の NS は、残差が床に着いた後も出口 M・流量が 10 万 step 単位で動き続けた (plan §5.1 #11)。
   作り直しの入力の δ_E と最終の出口 M は、長い窓の系列で動きが止まったか (出口 M は末尾の変化と推定の残りの変化 ≤ 2.5e-4) を確かめてから使う。
   **作り直しの入力の δ_E は、整定した場の 5 枚の窓 (5k ごと) で各枚の抽出 (帯は各枚で選ぶ、既定) の点ごとの平均にする** (2026-10-13、plan `tooling-nozzle-cooled-wall-ns-settling` §4.4)。1 枚だけでは、梯子が平らにならない区間 (冷却壁のスロート付近・直管・縮流部) で帯位置が切り替わり、局所で δ_E が 5e-4 r_t 動く。5 枚の平均でも 1 枚の切り替わりで 1e-4 r_t 程度は残る (既知の限界)。
5. **形を見る・出す** — 最終 NS の run から成果物一式を 1 コマンドで出す (plan
   [`tooling-nozzle-wt-deliverables`](../plans/accepted/tooling-nozzle-wt-deliverables.md)、2026-10-13):
   ```
   .venv-opt/bin/python -m forge_design.export.deliverables ../case/NN/run_XXXX --card ../case/NN/card_<tag>.yaml --problem ../case/NN/problem_<tag>.yaml
   ```
   → `run_XXXX/deliverables/`: STEP (2D の壁曲線・3D の内面)・STL・解析メッシュ (msh4.1/2.2)・STEP から同じ構造のメッシュを
   作り直す gmsh スクリプトと構造ファイル・Fluent 2D 軸対称メッシュ (`.msh`、境界の名前と型つき)・CGNS・Fluent の境界条件表
   (動作圧 0 Pa)・Excel 向けの形状表と壁面量表 (CSV、`xlsxwriter` があれば `.xlsx`)。README の冒頭に**寸法の出荷ゲート**
   (STEP とメッシュから測り直した寸法と要求カードの公差) と、壁面の熱量を熱設計に使えるか (`thermal_use`) が出る。
   物理壁の表現が `single_bspline` の run (最終 NS 用の雛形 `axismach_final_v1`) でないと STEP は出ない。
   壁の 2 階微分 (r″) の連続性は `nozzle_report` の「壁の形と微分」で見る (ユーザが重視する観点)。
   旧来の case 側のスクリプト (`case/42.isobutane_wt/export_wall_m42.py`、`case/44.vitiated_air_wt/export_wall_va.py`、
   `case/45.isobutane_m6_d155/cad_deliverables.py`・`viz_*.py`) は絶対パスを含むので新しい設計には使わない。

## 2. 出口・入口・条件を固定し、軸長を MOC でパラスタ → Euler → NS → 形状調整 → 凝縮

M6 (case/45) で実際に通した順番。各段の「何で判定するか」を先に決めてから回す (事前登録、AGENTS.md)。

1. **MOC パラスタ (CFD なし)** — 汎用ツールは無い。case のスクリプトを複製して使う:
   `case/45.isobutane_m6_d155/search_shortest.py` (R × L_c × M_knot の格子で、ハードゲート +
   `moc_diagnostics` の壁マージン ≥ 1°・特性線トポロジ ≥ 0.02 を満たす最短 x_F を探す)、
   `case/42.isobutane_wt/sweep_axislaw_A_MK.py`、`case/44.vitiated_air_wt/study_design_va.py`。
   定数 (条件・パス) がスクリプト内に書かれているので、問題に合わせて書き換える。
   有望形状の基準: 全長 x_F、壁マージン、特性線の健全さ、壁の r″ のなめらかさ。
2. **Euler (slip) で候補を評価** — 候補ごとに problem YAML を作り
   ```
   .venv-opt/bin/python -m forge_design.evaluate.runner_axismach ../case/NN/problem_X.yaml ../case/NN/run_XXXX_euler_X [--stages soft --cfl 2 --implicit-relax 0.7]
   ```
   → `metrics.json` (`collect()`: 軸 M の偏差・出口 M・オーバーシュート・出口一様性・流量比)。
   形状比較は `nozzle_report` の評価量 (r/r_w = 0.1 のMach 波・オーバーシュート・出口コア M) で行う。
   複数候補は `case/44.../run_study_va.py` のようにバッチで回す。
3. **NS + 境界層補正 (δ* ループ)** — `design/forge_design/feedback/deltastar_loop.py`:
   ```
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ../case/NN/problem_ns.yaml --euler-ref ../case/NN/run_EULER --init-integral --ic-from ../case/NN/run_EULER --run-dir ../case/NN/run_XXXX_ns_pass0
   .venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem ... --euler-ref ... --prev ../case/NN/run_XXXX_ns_pass0 --run-dir ../case/NN/run_YYYY_ns_pass1 [--stages none --cfl 5 --implicit-relax 0.7 --steps 12000]
   ```
   物理壁 = 設計壁 + δ_r。δ は CONTUR 積分法で初期化し、NS から E 法 (`band_select="edge"`) で δ_E を測る。
   **注意 (2026-10-13)**: `--ic-from <Euler run>` は、`interp_field.py` が化学種の照合で拒否する (Euler は輸送の記録が無く、NS は
   `transport` を持つので species_hash が違う)。新しい設計の最初の NS は `--ic-from` を付けず、prepare_ns の 1 次元等エントロピー場から
   段階起動 (`--stages full`) で始める (case/70 run_0005)。2 回目以降の NS は前の NS の場から (`--prev`)。
4. **C2 方式 (出口 δ で CONTUR を較正し r_t を解く)** — CONTUR の摩擦係数倍率 k_f を「出口での δ_C = 測った δ_E」
   に合わせ、出口径 0.775 m (例) から r_t を解き直す。**k_f は YAML キーが無く case スクリプト経由**:
   `case/45.isobutane_m6_d155/prep_c2pin.py` (run 準備、k_f を引数で渡す)、`c2pin_solve.py PASS1_RUN EULER_REF
   [--base YAML] [--out NAME] [--cond-steps N] [--cond-out N]` (未緩和 δ_E から k_f・r_t を一発で解き
   `<NAME>.yaml`・`<NAME>_cond.yaml` を書く)、バッチ例 `run_recal_chain.sh`。
   δ は **`delta_r_next.csv` の `delta_E` 列 (未緩和)** を使う (第 2 列 `delta_target` は ω 緩和後 — 取り違えた前例あり)。
5. **形状を整える (生産レシピ、2026-10)** — 探索段階では不要、形状を詰める段階で入れる:
   - `wall_repr: joint` — MOC 点を位置 + 壁角で同時に当てはめた 5 次 B-spline (壁の r″ が連続でなめらか)。
   - `initial_line: cfd` — スロート初期線と軸アンカーを Hall 解でなく Euler の遷音速場から取る
     (`initial_line_run`/`initial_line_res` で凍結源を指定; 凍結源は Hall 初期線の当てはめ壁を Euler で解いた場)。
   - `Md_moc_offset` — Euler の出口コア M を 6 に合わせる 1 係数の較正。**Euler は `mesh_euler` の壁に寄せない配点で決める**
     (2026-10-07 から。case/45 は 2000 × 97・軸方向は NS の G1 と同じ・全断面で壁の第 1 間隔の比 0.005)。NS と同じ壁に寄せた
     G1 の配点 (1.3e-5・スロート 4.5e-6) の Euler は、スロート付近の全温が Tt を数百 K 超えたまま整定しなかった
     (plan `verification-case45-euler-total-enthalpy` §9 E2。全域 0.005 の配点は同じ窓で全領域 |T₀ − Tt| ≤ 0.103 K)。出口コア M は最終断面の
     η ∈ [0.05, 0.7] で、NS の基準格子 G1 の固定の η の列に線形補間した平均 (M_common) で較正し、自格子の平均も残す (同 plan §6 E4)。
     Euler の較正の合格を NS に移せるとは限らない (配点も粘性も違う) ので、**NS の出口 M は NS で判定する**。
     旧手順「生産 NS と同じ格子パラメータの Euler で決める」(plan `tooling-nozzle-cfd-pinned-initial-line` §5.1 #11f。粗い Euler 格子
     1100 × 65 で決めると細分格子の NS で出口 M が約 −0.0008 ずれた経緯) は、G1 の Euler の全温が整定しない (E2) ため取り下げた。
     case/45 の旧値 +3.770e-4 は G1 の Euler (run_0113+0114) 由来で、新しい配点でやり直した (同 plan §6 E4・E4V)。現行値は 6.8825162455159465e-06 (run_0164)。
   - `pw_upstream: poly` (既定。新しい YAML には明記する) — スロートより上流に δ_r を足さず、配管〜設計スロートを 5 次多項式 1 本にする
     (plan `tooling-nozzle-upstream-poly-and-throat-sizing`、2026-10-07)。`pw_ramp` は書かない (併記すると例外)。
     旧設定: `pw_upstream: ramp` + `pw_ramp: [−11, −6]` (縮流部側で δ_r をなめらかに入れる区間)。case/45 の生産 YAML は NS・凝縮の再評価
     (同 plan §6 U4) で採否を決めるまでこの旧設定を明示している。
   - `wall_fit_mono_r2: [0.0, 1.5]` — joint 壁のスロート直後の r″ の山を消す単調拘束 (case/45 の生産、2026-10-07)。
     Euler は実務判定 (全量で 差 + 2SE ≤ 許容幅) に合格、dry NS は登録したゲートに合格、凝縮 NS は 4 量の準定常に合格。新旧の場の差は参考値 (plan `tooling-nozzle-throat-monotone-r2`)。
   - `deltastar_initializer: {model: contur, a_crocco: 1.0, cf_scale: <k_f>, n_scale: 1.0}` — C2 で較正した k_f を YAML に書く。
     `deltastar_loop --init-integral` も読む (書かないと k_f = 1 の物理壁になる)。最終 NS の数値レシピ (CFL 1・60000 step など) は
     YAML ではなく投入スクリプトが持つ (case/45 の手本 `run_mono_ns_chain.sh`、延長 `run_mono_ns_ext.sh`)。
6. **NS の格子** — 壁解像 (y1+ ≤ 1、`check_wall_resolution.py --over-frac 5`) を満たす細分格子を使い、
   格子ゲート (細分前後で δ_E(x_F) の差 ≤ 1 %、三水準目で確認) を通す。case/45 の生産格子 (NS の `mesh`、G1):
   ni 2000 × nj 97、第 1 セル 1.3e-5 (スロート 4.5e-6)、`throat_refine 4`・`throat_width 3`、AR ≤ 5000 (境界層層の例外)。
   この壁に寄せた配点は NS 用で、Euler には使わない (Euler は `mesh_euler`、§0)。
   **細分格子の NS は本段 CFL 5 で発散する** — 段階起動 (soft → mid → 本段) + 本段 cfl 1・60000 step
   (累積 CFL を cfl 5 × 12000 にそろえる)。
7. **最終判定** — 事前登録したゲート (case/45 の例: 出口半径 0.775 m ± 0.1 mm、δ_E/δ_C 1 ± 0.5 %、出口コア M 6.000 ± 0.05 % [2026-10-07 ユーザ決定、旧 ± 0.02 %]、
   r/r_w = 0.1 の波 ≤ 0.01 %・オーバーシュート ≤ 0.035 %、壁解像 PASS) を、各量の時系列に
   `check_quasisteady.py --series-csv` をかけて判定する (`nozzle_report` の末尾幅は VERDICT の代わりにならない)。
   時系列は `case/45.isobutane_m6_d155/exitM_sampling_ab.py` (環境変数 `EULER_REF`・`SOLVE_JSON`・`FINAL_PROBLEM`)。
8. **凝縮 ON** — 最終 NS の場を `convert_species_field.py --mode conserve` で液相付きの種構成に変換して引き継ぐ
   (restart_field は種構成が違うと拒否する):
   ```
   python3 solver_density_cuda/tools/convert_species_field.py LAST_res.h5 RUN_COND/nozzle.h5 --meta RUN_COND/species_meta.yaml --src-run RUN_NS --dst-run RUN_COND
   ```
   problem は `<NAME>_cond.yaml` (`tp_species: split_h2o` + `condensation`)、cfl 1。凝縮 4 量
   (始まり位置・S_max・出口 g・出口コア M) が末尾 5 枚で STEADY になるまで回す。
9. **報告** — `nozzle_report` → pptx ([`nozzle-design-outputs.md`](nozzle-design-outputs.md))。
   case README の run 一覧と plan §9 に記録する。

## 3. 作成済みのノズルで Euler / NS を回す

**まず「形状を固定するか」を決める。** problem YAML の `spec` (Pt・Tt)・`gas` を変えて `runner_axismach` や
`deltastar_loop` に渡すと、`design_chain` が MOC と δ 補正をやり直し、**壁そのものが変わる** (= 作り直し、§1・§2)。
既にあるノズルで条件だけ変えたいときは、YAML を通さず run ディレクトリの config を直接書き換える (3a)。

### 3a. 形状は固定、入口条件・背圧・入口分率だけ変えて回す (よくある用途)

ツール **`solver_density_cuda/tools/rerun_conditions.py`** を使う (plan [`tooling-rerun-conditions.md`](../plans/accepted/tooling-rerun-conditions.md)、仕様は `methods/design/overview.md`「既存 run の条件変更」)。
既存 run の `nozzle.h5` (物理壁の格子) を再利用し、bcond の floats だけを書き換え、参照 run の最終場を初期場に入れ、変更点を `RERUN_CONDITIONS.json` に残す。**forge は起動しない。**

```
python3 solver_density_cuda/tools/rerun_conditions.py REF_RUN NEW_RUN [--Pt P --Ps P | --keep-Ps] [--Tt T] [--Y H2O=0.09] \
    [--Tw T | --keep-Tw] [--scale-ic pt] [--steps N --out-interval N --cfl C] [--dry-run]
```

- 対応する入力: 単一の `inlet_Pressure` (Y{s} 形式、inletProfile なし)・`outlet_statPress`・`wall`/`wall_isothermal`/`slip`・`axis`、`nozzle.h5` を mesh/value に使う run。外れれば作る前に止まる。
  NS / Euler は実効 config で判定する (粘着壁があれば NS; 全壁 `slip` は粘性・熱伝導・種拡散・乱流が無効な config だけ Euler として受理し、輸送が有効なら止まる)。重複した YAML キーを含む config は止まる (ツールとソルバで解釈が違うため)。
- **Pt を変えるときは背圧の指定 (`--Ps` か `--keep-Ps`) が必須**。node の出口は壁際の列が常に背圧を見るので、Pt だけ下げると出口の壁際から崩れる。
  ツールは参照の出口断面の静圧に対する `Ps/(f·P_exit_ref)` を表示する (相似にするなら Ps も同じ比で)。
- 乾き成分 (lump) の組成変更は v1 では止まる (化学種の定義が変わり、どの引き継ぎ経路も拒否する)。H2O 分率 (`--Y H2O=...`) は可。
- 入口の k・ω は変えない (指定時のみ)。等温壁で Tt を変えるときは `--Tw` か `--keep-Tw` を明示する。
- **回し方** (`RERUN_CONDITIONS.json` の `recommended_stages` に従う): 条件の変更なし → `run_staged_ns(run, stages="none")`; 条件を変えた → `stages="full"` (段階起動);
  **Pt だけを変えた → full + 本段 cfl 1 (`--cfl 1.0 --steps 60000`)、`--scale-ic pt` (保存量を Pt 比で一様スケールする初期場変換) を推奨** (2026-10-06 確定: scale なしは規定時間内に準定常に達せず、入口配管の壁際に逆流域が残った)。**Pt と Tt・組成・壁温を同時に変える、または凝縮 run で Pt を変えるときは `--scale-ic pt` が使えないか推奨外で、起動・整定は未検証** (full・本段 cfl 1・60000 step と警告が記録される)。
  粗格子の Pt 0.8 倍は、段階起動でも本段 cfl 5 では発散した (plan §6 (ii′))。Tt・H2O を変えたときは本段の量の整定に時間がかかる (粗格子で本段 36000 step でも流量・δ_E は DRIFTING; 整定長は未確定、plan §5.1 #11)。
- **Euler 参照**: δ_E の抽出や流量比には、同じ条件にした Euler の rerun を対で作る (旧条件の Euler 参照だと δ が数 % 動く)。
- 判定: `check_convergence.py` (`--segment` で本段だけ)・`check_quasisteady.py --series-csv` (報告する量)・NS は `check_wall_resolution.py` (Pt・Tt が変わると Re と y1+ が変わる)。
- 注意: 壁は元の条件で設計したもの。条件を変えると境界層の厚さが変わり、出口 M・試験部の一様性は設計値からずれる (それを見るのがこの計算の目的のことが多い)。ずれを壁に返したいなら §2 の作り直し。

### 3b. 既存の problem YAML から同じ設計を回し直す (作り直しを含む)

- Euler: `.venv-opt/bin/python -m forge_design.evaluate.runner_axismach PROBLEM RUN_DIR [--stages full|soft|none] [--cfl C] [--implicit-relax 0.7] [--steps N] [--ic-from RUN]`
- NS: runner に NS フラグは無い。`deltastar_loop` (`--init-integral` で CONTUR の δ から、`--prev RUN` で前の NS の δ_E から)、
  または C2 で解いた k_f で `prep_c2pin.py RUN K_F --problem YAML --ic RUN --stages full` → `run_staged_ns`
  (`case/45.isobutane_m6_d155/run_recal_chain.sh` が手本)。
- YAML の Pt・Tt・ガスを変えると壁が変わる (上記)。同じ壁を再現したいだけなら YAML を変えない。
- **コードの版も壁を変える**: 2026-10-13 の出口の終点の修正 (de35cb0b) から、終端の壁点を丸めで落としていた設計は同じ YAML でも壁・メッシュ・STEP が 1 点間隔長くなる (case/70 c2b +0.846 mm、case/45 生産 +5.353 mm)。古い run と同じ壁を再現するには、その run の `wall_repr.json`・`wall_physical.csv` (保存された壁) と `RUN_PROVENANCE.txt` のコードの版を確かめる (YAML だけでは再現しない)。

### 3c. 外部の壁座標 (CSV・CAD) しか無い場合

axismach パイプラインに取り込むキーは**未実装** (`wall_csv` 相当なし)。Gmsh `.geo` → `gmsh` → 変換器 → config を手で組む
一般経路 ([`calculation-workflow.md`](calculation-workflow.md)) を使い、config は skill `forge-config` で組む。

### 共通

- 同じ格子の続きは `restart_field.py` (ビット一致を検査)、格子が違えば `interp_field.py`。一様初期場からは始めない。
- 判定: `check_mesh_quality.py` (投入前)、`check_convergence.py` (同一設定の区間)、`check_quasisteady.py` (報告する量そのもの)、
  NS は `check_wall_resolution.py`。残差が plateau (`NOT CONVERGED`) になるのはこの系の常態 — そのまま記録し、量の準定常で判断する。

## 4. 未実装・注意 (2026-10-06 時点)

- MOC パラスタ・作図・凝縮の引き継ぎは **case スクリプト頼み** (汎用 CLI なし)。壁の出力は `forge_design.export.deliverables` (§1 ⑤)、
  C2 の k_f は `deltastar_initializer.cf_scale`。
  統一は campaign-recipe plan の残作業 (#6〜#8・#14・#15)。
- 出口径・スロート径を spec に書けない (r_throat + M_design → `solve_rt` / `solve_rt_throat` → 手で書き換え。`spec.sizing` は記録だけ)。
- `design_chain` には fold・トポロジ・壁マージンの検査が入っていない。
- AWS の共有バイナリは他セッションが作り直すことがある — 自分の worktree でビルドしたバイナリを使い、
  `RUN_PROVENANCE.txt` の sha256 を記録する (case/45 は `~/forge-wallfit-bin`)。
```

## 参考: `case/71.m55_combustion_wt/README.md`

```
# case/71 — M5.5 の燃焼ガスの風洞ノズル (要求カードから NS まで)

2026-10-14 のユーザ依頼:「mach 5.5、入口配管長 0.5 m、収縮部 0.5 m、スロート曲率 5 mm、スロート直径 Φ200、スロート〜出口 8 m のノズルを作って、
凝縮計算はせずとも、NS 計算までやりきる」。うまくいったら、何も知らない低級なモデル (Haiku) に同じことをさせ、つまずいた所を流れの側で直す。

## 要求 (カード `card_m55.yaml`)

| 項目 | 値 | 出所 |
| --- | --- | --- |
| 出口マッハ数 | 5.5 | ユーザ指示 |
| ガス | イソブタンの燃焼ガス (当量比 0.9、case/45 と同じ組成: 質量分率 CO₂ 0.1677・H₂O 0.0858・O₂ 0.0220・N₂ 0.7245) | ユーザ選択 (提案値) |
| 全圧・全温 | 3.5 MPa・1,250 K | ユーザ指示 |
| スロート | 物理の直径 0.200 m、曲率半径 0.500 m (R ≈ 5) | 直径はユーザ指示 (物理)、曲率は**暫定** (ユーザ「500 mm かも、自信なし」; 5 mm は R = 0.05 で拒否される) |
| 入口 | 半径 0.300 m (提案値、ユーザ了承)、直管 0.50 m、縮流部 0.50 m | 長さはユーザ指示 |
| スロート〜出口 | 8.0 m | ユーザ指示 |
| 壁 | 等温 350 K | ユーザ指示 |
| 許容差 | 雛形の仮の値 (長さ 1 mm・径 0.02 mm) | ユーザ「よくわからん」 |
| 凝縮 | 計算しない (dry の NS まで) | ユーザ指示 |
| 格子の節点数 | NS・Euler とも ni 4,800 × nj 97 = 465,600 (手組みの格子 `problem_m55_ns.yaml`)。#3b (最終の雛形) から NS は ni 4,800 × nj 92 = 441,600 (径方向 `axis_cap_frac 0.03`・比 ≤ 1.2、ユーザ 2026-10-15「上限3%」; `problem_m55_final.yaml`)、Euler は 465,600 のまま (`problem_m55_final_euler.yaml`) | ユーザ決定「a」(2026-10-15): 見積もりの偏りを M6 の実績で補正した約 45 万点で Euler → NS (上限 30 万を超えることは了承済み)。それ以前に ni 10,200・約 99 万点を一度了承 → 同日取り下げ (「あまりにメッシュ多すぎる、古いやり方で作ってない？」) |

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_euler` | Euler (slip) の参照 (AWS A)。`problem_m55_final_euler.yaml` (最終の雛形の Euler の写し)、格子 ni 4,800 × nj 97 = 465,600 (全域 0.005)、`runner_axismach` の `--stages full` (soft 3,000 → mid 3,000 → 本段 cfl 1・12,000)、FP64 のビルド (`~/forge-integ-dflow2-fp64`、2fd36b0c)・`FORGE_CUDA_BLOCKSIZE=128` | 完走 (227 s): 出口の軸の M 5.5002、dM 最大 0.06 % Md、NaN なし。`check_convergence --segment` NOT CONVERGED (残差 2.1〜2.5 桁で停滞; M6 の Euler の参照と同じ型 → δ* の参照としては未適格と記録、plan §5.1 #3b)。δ* の参照は `res_12000.h5` (sha256 93cb88b0…) に固定。記録 [`_runs/run_0001_euler/`](_runs/README.md) | ref (run_0002_ns・run_0003_ns_cont の Euler の参照) |
| `run_0002_ns` | NS pass 0 (AWS A): `deltastar_loop --init-integral --stages full` (CONTUR の積分法の初期壁、停止規則・上限 200,000)、`problem_m55_final.yaml` (手組みの格子 ni 4,800 × nj 92 = 441,600、`axis_cap_frac 0.03`)、Euler の参照 `run_0001_euler`、FP64 のビルド | 本段 0〜183,998 step でディスクが満杯 (OSError) になり終了、停止規則は未到達。本段の `check_convergence --segment` (manifest の写しに本段を足した事後復元) NOT CONVERGED (各列 158k〜172k で最小の後に上昇)、停止規則の系列 DRIFTING (θ_r 0.4・0.7 TRANSIENT-UNSETTLED)、壁解像 PASS (最大 1.004)、res_183000 で P0 > Pt が 50,441 点・最大 +0.24 %。全圧の A/B (#3c) と診断は plan §5.1 #3b・#3c。記録 [`_runs/run_0002_ns/`](_runs/README.md) | ref (run_0003_ns_cont の初期場 res_183000) |
| `run_0003_ns_cont` | NS の延長を全部 double で (AWS A、plan §5.1 #3d): `run_0002_ns/res_183000.h5` を同じ格子に index restart (`restart_field.py --keep-src-dtype`、9 量ビット一致)、前段なし (`run_staged_ns` stages none)。本段の設定から `physProp.thermoFloat: 0`・`time.deltaT.implicitSolvePrecision: 1`・`outStepInterval` 2,000 だけを変え、環境変数 `FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1` (`RUN_PROVENANCE.txt` に記録)。停止規則は明示 (窓 2 万・保持 4 万・上限 = 追加 200,000 step)、Euler の参照は `run_0001_euler/res_12000.h5` (sha256 93cb88b0…) に固定、帯はこの run の最初の出力 res_2000 から (θ_r は run_0002_ns と別の系列)。残る float: 化学種の拡散の D_s・h_s・X の float の写し (`speciesTransport_d.cu:240-356`)・リミッタの長さの `sqrtf`・XMF の型の宣言 | 投入 2026-10-15 00:24 UTC、26.9 ms/step。**判定 (本段 0〜199,999)**: 停止規則 CENSORED (200,000; θ_r 0.4・0.7 の D 6.3e−4・5.3e−4 > 5e−4、単調に減少)、`check_convergence` NOT CONVERGED (ρ・k・Y falling、運動量・ρE flat; float の末尾より 1.6〜3.2 桁下)、θ_r 3 断面・ṁ・出口コア M (5.52045) STEADY、壁解像 PASS (y1+ 最大 1.002、超過 0 %)、P0 > Pt 46,936 点・最大 +0.462 %、δ の 2 窓 FAIL (measured、δ_target の振れ 0.031 mm = δ の 0.3〜0.5 %、x 1.5〜2.2 m)。ユーザ決定「a」で関門の明示の上書きで δ* の補正 (plan §5.1 #3e) へ。判定の原文 `_runs/run_0003_ns_cont/judge/`、記録 [`_runs/run_0003_ns_cont/`](_runs/run_0003_ns_cont/README.md) | ref (δ* の補正の 1 回目の入力) |
| `run_0004_euler` | δ* の補正の 1 回目の Euler の参照 (AWS A、plan §5.1 #3e (4)): `problem_m55_p1_euler.yaml` (`wt_card --ns-correction run_0003_ns_cont` で物理スロート径 200 mm に r_t を解き直したカード `problem_m55_p1_card.yaml` → `_mesh/hand_mesh.py build`、r_t 0.10024473616136853 m = 前の +0.34 %、M_knot 2.80)、格子 466,570 点 (品質 PASS、AR 最大 ≤ 1000)、`runner_axismach` の `--stages full` (本段 cfl 1・12,000)、FP64 のビルド。台本 `_runs/run_0004_euler/p1_euler.py`・`_runs/run_0005_ns/p1_chain.sh` | 完走 (228 s、非有限なし): 軸の出口 M の最大 5.5013、dM 最大 0.058 % Md。`check_convergence` NOT CONVERGED (1.2〜1.8 桁で横ばい; run_0001_euler と同じ型)。参照は `res_12000.h5` (sha256 78294389…) | ref (run_0005_ns の Euler の参照) |
| `run_0005_ns` | δ* の補正の 1 回目の NS (AWS A、plan §5.1 #3e (5)): run_0003_ns_cont の δ_E から作った補正の表 `problem_m55_p1_card_delta_r.csv` の物理壁、`problem_m55_p1.yaml` (手組みの格子 442,520 点、品質 PASS・AR 緩和 違反 0)、初期場は run_0003_ns_cont の最後の場から補間、`--stages full`、**精度は既定 (float; ユーザ決定 2026-10-15)**・FP64 のビルド、出力 2,000、停止規則 (窓 2 万・保持 4 万・上限 200,000)、Euler の参照 `run_0004_euler/res_12000.h5`。台本 `_runs/run_0005_ns/p1_ns.py`・`p1_follow.sh` (終わったら成果物を作る) | 2026-10-15 07:37 UTC に終了 (7,627 s)。**判定 (本段 0〜199,999)**: 停止規則 CENSORED (200,000; θ_r の D 6.6e−4・6.0e−4・5.3e−4 > 5e−4、2 万 step ごとに約 0.65 倍で減少)、`check_convergence --segment` NOT CONVERGED (ρUy・ρE・ω flat)、停止規則の系列 DRIFTING (r だけ; θ_r・ṁ・出口コア M STEADY)、壁解像 PASS (y1+ 最大 0.999、超過 0 %)、非有限 0・P0 > Pt 46,609 点・最大 +0.458 %、δ の 2 窓 FAIL (δ_E の振れ 0.060 mm @ x 2.34 m) → 関門は止めた (上書きなし)。**出口コア M 5.50145** (pass 0 は 5.52042)。**試験部の軸に M の波 −0.76〜+0.55 %** (x 4.8・5.9・7.0 m、res_100000 と同じ; pass 0・Euler は平ら) — 原因は諮問中 (plan §5.1 #3e)。δ_E − δ_in (診断の抽出) は x ≥ 0 で最大 0.52 mm (出口 −0.39 mm、pass 0 は 8.1 mm)。成果物は手元で作成 (`run_0005_ns/deliverables/`: cad・mesh・fluent・excel はエラーなし、寸法の関門はカードの 5 項目 PASS・整合の 1 項目で FAIL、plan #3j)。`metrics.json` の ṁ は列の積分で −0.96 % (plan #3i)。判定の原文 `_runs/run_0005_ns/judge/`。A の途中の全場 95 枚は削除 (res_0・100000・140000・160000・180000・200000 を残す)、手元の写しは res_200000 だけ | ref (δ* の補正の 1 回目の結果) |

- **AWS A の途中の全場の削除 (2026-10-15、判定済み)**: `run_0003_ns_cont` の 96 枚 (res_0・140000・160000・180000・200000 は #3f の再生用に残す)、`run_0002_ns` の 21 枚 (res_0・res_183000 = run_0003 の初期場を残す)、`run_0001_euler`・`run_0004_euler` の各 11 枚 (res_0・res_12000 = Euler の参照を残す)。壁・出口の出力・残差履歴・判定の記録は残した。

## 経緯

- 2026-10-14〜15: 旧の M_knot の規則 (`min_theta_max_v1`) でカードを生成 (AWS A `~/forge-m55`、約 98 分) → 最終の組の 4800 の網の長さ合わせが全点未収束で失敗 (plan `tooling-nozzle-mknot-fast` §3)。
- 2026-10-15: 新の規則 `coarse_min_theta_v2` (feature/mknot-fast) で M_knot 2.80 (144 s) → NS の格子の節点数の上限 (30 万) で停止 → ユーザ了承で ni・上限を明示。
- 2026-10-15: ユーザが 99 万点の了承を取り下げ。今のカードの格子の見積もりは第一層を 2 値で持つ古いやり方で、拡大部全体が最も厳しい場所の薄い第一層になる。場所ごとの第一層 (`wall_first_frac_table`) と流れ方向の密度 (`x_density_table`) で見積もり直してから格子を決める。新の規則での生成は M_knot の最終の組の途中で止めた (A `~/m55_card/`、結果なし)。
- 2026-10-15 (plan `tooling-nozzle-card-local-mesh` §5.1 #3a): カードに `mesh.ni: 3000` (生成のためだけ) を足して AWS A `~/forge-cardmesh` (479723ab) で生成 → `problem_m55_card.yaml` (M_knot 自動 2.80、r_t 0.09990877093936801 m、177 s)。その `mesh`・`mesh_euler` を `_mesh/hand_mesh.py build` で手組みの格子 (場所ごとの第一層の表・流れ方向の密度の表・壁法線の近壁層・17 桁) に置き換えて `problem_m55_ns.yaml` (Euler・NS の両方に使う)。CFD 前の物理壁で作った格子の検査 (forge は起動しない、A `~/m55_mesh_hand/check/`): NS `check_mesh_quality --ar-max 5000` PASS (AR 最大 4849.4)・AR > 1000 の 47,607 セルのスキュー最大 0.0015 (違反 0; 出口面は非平面 0.62 mm で記録だけ)・変換後の第一層厚 ≤ c·y₁,req (全 4,800 断面)・見積もりの y₁⁺ 最大 0.797。Euler PASS (AR 最大 6.1)。再現の入力と記録は `_mesh/` ([`_mesh/README.md`](_mesh/README.md))。
- 2026-10-15 (plan §5.1 #3b): カードの雛形を `axismach_final_v1` に替えて (探索の雛形の物理壁では NS の停止規則が使えない) AWS A `~/forge-cardmesh` (94790833) で生成 → `problem_m55_final_card.yaml`・`_euler.yaml` (M_knot 自動 2.80 で不変、r_t 0.0999084914252905 m、寸法 PASS、194 s)。`_mesh/hand_mesh.py build` で手組みの格子に置き換えて `problem_m55_final.yaml` (NS)・`problem_m55_final_euler.yaml` (Euler)。NS の径方向はユーザ決定 (「上限3%」) で `axis_cap_frac 0.03`・nj 92 (隣接比の最大 1.1974)。検査 (FP64 の変換器、forge は起動しない): NS 品質 PASS (AR 最大 4849.0)・AR 緩和 違反 0 (AR > 1000 の 36,115 セル、スキュー最大 0.0015)・第一層 実寸 / (c·y₁,req) 最大 0.9999997 (超過 0)・見積もりの y₁⁺ 最大 0.797・軸側のセル / r_w 0.030 (全断面); Euler PASS (AR 最大 6.1)。`run_0001_euler` を準備 (未起動)。記録は `_mesh/hand/final/` ([`_mesh/README.md`](_mesh/README.md))。
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
証拠: 足りる | 証拠不足: <要るもの>
```
最後の行 (証拠の欄) は必須で、どちらか一方だけを書く: 結論を出すのに要る証拠がそろっていれば `証拠: 足りる`、
足りなければ代わりに `証拠不足: <要るもの>` を足りないものごとに 1 行 (この行があると諮問は不成立として扱われる)。
`不足情報:` は結論を妨げない追加の希望の欄で、それだけでは不成立にしない。
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
