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

## ブリーフ (`notes/reviews/briefs/2026-10-15-m55-settle-and-monitor.md`)

# 諮問ブリーフ: NS の本段・判定・出力の新しい既定 (plan §4.4) の設計の点検、と試験部の残りの傾き (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1: plan §4 を新しく書いた; 7: 結果の解釈 (残りの傾き))
- plan: `plans/active/tooling-nozzle-card-local-mesh.md` **§4.4 (新規)**・§5.1 #3q・#3r・#3t・#3v・#3x・#3y・§6 V5・V7
- 前の諮問: `notes/reviews/2026-10-15-m55-three-coef-diagnose.md` (3 係数・リミッタの基準値 (b)・B0 の既定化の範囲)

## ユーザの決定 (2026-10-15、原文)

- 「もうええやろ、line implicitで十分だよ」「line implicitを規定化してほしい」
- 「ns計算で、生の排除厚さ入れるやり方、もう取る意味ないよね」
- 「落ちた桁数…近々で、これで語るのやめにしようか。残さが準定常になり、流量バランスが取れ、目的とする物理量が準定常担ったら収束と判定しないか」
- 「入口出口の流量差や、断面平均のマッハ数なりを定期的にprobeとして出力する機能欲しい」「出口境界と入り口境界と軸中心のh5を出し続けてくれてもいいな…そこから後処理的にマッハ数・流れ角出したりとか」「結果丸ごと出し続けて容量圧迫されるのがな」
- 「そのプローブから見る方式だと、境界層厚さの時系列チェックはできなくなるからご認識あれ。それは私は別にいいと思ってます。最終だけあればね」

## 観測事実

1. 3 係数の壁の NS `run_0010_ns_b0` (AWS B、B0、初期場 = 2 係数の NS run_0008 の最後の場から補間): 停止規則は本段 46,000 で条件、86,000 で到達 (θ_r の D は 7 万 step で 1e−8 台)、停止規則の系列・V5c の系列 ALL STEADY、壁解像 PASS、非有限 0。`check_convergence` (本段の履歴) は NOT CONVERGED: 本段の始めの残差が前の解から始めたので小さく、低下 2.3〜3.0 桁 (末尾は run_0008 の床の 0.62〜1.27 倍)。`--from-floor run_0008` も過渡の山で NOT CONVERGED (`_runs/run_0010_ns_b0/judge/check_convergence_main*.txt`)。
2. 停止規則の見張りは forge の終了の後に δ の 2 窓 (受け渡しの関門用) で全場 21 枚から δ* を抽出し、1 時間強かかった (run_0007・0008・0010 はユーザ了承で止めた)。抽出は高速化した (113.8 s → 15.9 s 手元 12 コア、ビット同一、AWS 3 並列で 71 s)。
3. 全場 1 枚 約 100 MB、2,000 ごとに出すと 20 万 step で約 10 GB。境界の h5 は出口 29 KB・壁 1.2 MB (中身: Ps・Pt・Ts・Tt・Ux・Uy・Uz・ro・roU・roe・k・ω、壁は qwall・twall・utau・ypls なども)。境界ごとの `outputHDFflg` (`bcondConfig.yaml`) で選び、間隔は全場と同じ `outStepInterval` (`output.cpp:315`・`:333`)。M5.5 は出口・壁が 1、入口・軸が 0。
4. forge の point probe (`probe/point_probes.cu`、`probe.yaml`) は T・P・Ux・Uy・Uz を点ごとの文字の file に書く (音速・密度は無い)。
5. 試験部の M (r/r_w = 0.1、x 4.0〜7.95 m): run_0010 は平均 −0.002 %、x 4.0 m の −0.094 % から出口 +0.058 % へ単調に上がる傾き (波なし); 2 係数の run_0008 は平均 +0.147 %・x 4.4 m の +0.224 % から +0.067 % へ下がる。δ の誤差 (入れた δ − NS の δ_E): run_0010 は x 0.5 m +0.30・1 m +0.62・2 m +0.30・3 m +0.17・4 m −0.13・5 m −0.15・6 m −0.09・7 m +0.15・7.9 m +0.14 mm。同じ規則でもう 1 回較正し直しても最大 0.624 → 0.612 mm・x ≥ 3 m 0.257 → 0.214 mm にしかならない。Euler (設計壁) の試験部は ±0.02 % で平ら。図 `run_0010_ns_b0/report/fig_axis_deviation.png` (NS の軸と r/r_w = 0.1 の線は x ≈ 40 r_t で Euler より遅れて立ち上がり、出口へ上がる)。
6. 停止規則の量 (`stop_rule.recorded_quantities`): 出口コア M = 最後の列の r/r_w 0.05〜0.7 の節点の M の単純平均 (軸の節点は入らない)。r = (ṁ_in − ṁ_out)/ṁ_in は全場の入口・出口の列の折れ線で 2π∫(ρu dr − ρv dx) r。

## 設計 (plan §4.4 の 1〜7)

本文は plan §4.4 を読むこと。要点: forge に `time.fieldOutStepInterval` (出力だけ)、雛形の既定で入口・出口・壁・軸の境界の h5 を 2,000 ごと・全場を 20,000 ごと + 最後、本段 B0 の既定化 (リミッタの基準値は最初の本段の forge の自動の値を記録して継承)、判定の道具 `check_settled.py` (残差の末尾の準定常・境界の h5 の流量の収支・目的の量の準定常の AND)、停止規則も同じ 3 条件で境界の h5 から、θ_r と δ の時系列はやめ最後の全場から δ_E を 1 回抽出、直接の δ_E の経路は既定から外す、AGENTS.md は道具を試してから直す。

## 問い

1. §4.4 の設計に穴はあるか。特に: (a) 残差の「準定常」の判定 (log10 の末尾の窓の傾き・振幅) で、発散の兆し・ゆっくりした上昇・リミットサイクルを取りこぼさないか。窓・閾値をどう置くか。(b) 流量の収支を境界の h5 で測るときの注意 (出口の面が非平面 (x の広がり 5.3e−4 m)、軸の点、node の境界の値と流束の値の違い)。(c) 目的の量に軸の h5 の M を使うと軸の節点の癖 (M +0.1〜0.17 %) が入る — 判定は時系列の準定常なので癖そのものは効かないと考えたが、r/r_w = 0.1 の線を point probe で取る方がよいか。(d) θ_r を判定から外したことで、境界層のゆっくりした発達 (M6 で問題になった遅いモード) を見逃す危険をどう評価するか (ユーザは最終だけでよいと言っている)。
2. 停止規則の「保持 4 万 step」を B0 で短くしてよいか (run_0008・0010 は条件を満たしてから D が 1e−7 台まで 2 万 step で落ちた)。事前に決めるならどの値・根拠か。
3. リミッタの基準値を「最初の本段の開始時の forge の自動の値を記録して継承」にする設計で、次の壁の pass (形が少し変わる) に継承してよいか、毎回の新しい問題では自動に戻すか。
4. 試験部の残りの傾き (観測 5) の原因の第 1 仮説と、安い切り分けの方法。δ の誤差の分布 (x 1 m で +0.6 mm 厚い、x 4〜6 m で −0.13〜−0.15 mm 薄い、x 7〜8 m で +0.14 mm 厚い) と M の傾き (x 4 m で低く出口で高い) の関係をどう読むか。3 係数の形の限界 (もう 1 回の較正でほぼ変わらない) を踏まえ、4 係数目 (例: x の 2 次) や局所の補正を足すべきか、別の要因 (B0 の解・格子・出口の境界条件) を疑うべきか。

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
- **式**: δ(x) = c·δ_C(x) + δ₀。δ_C は今の CONTUR の積分法 (`integral_delta_r`: `integral_bl` → 5 次 P-spline、cf_scale 1) の出力、c は無次元、δ₀ は実寸 [m] で全長一定 (拡大部の壁の傾きを変えない)。スロートより上流は今どおり δ を直接足さない (`pw_upstream: poly`) が、**poly は下流の値・傾き・曲率を端条件に作り直すので δ₀ で上流の形は変わる** (「上流が変わらない」ではない; codex 2026-10-15 ①)。寸法は今どおり物理壁の大域の最小半径で解く (`deltastar_loop` 443 行)。
- **較正**: 同じ形の壁 (CONTUR の厚み、cf_scale 1) で回した NS の δ_E (平滑化後、未緩和) のスロート (x = 0) と出口 (最後の列) の 2 点に合わせる: c = (δ_E,e − δ_E,t)/(δ_C,e − δ_C,t)、δ₀ = δ_E,t − c·δ_C,t。δ_C はその NS と同じ r_t で評価した値 (= その NS の物理壁の δ)。pass 0 (`run_0003_ns_cont` res_200000、r_t 0.0999084914 m) から: δ_C スロート 0.09159・出口 49.3301 mm、δ_E −0.24396・41.7109 mm → **c = 0.852073、δ₀ = −0.32201 mm**。新しい r_t では CONTUR を新しい r_t で計算し直してから c・δ₀ を当てる (CONTUR の r_t の依存は式が持つ)。
- **割り切り**: 係数 2 つなので途中は合わない (pass 0 の比較で x 1 m +1.2・3 m +2.0・5 m +2.0・7 m +1.1 mm、NS より厚い)。壁は CONTUR のなめらかさのまま。
- **実装**: `deltastar_initializer` に `delta_scale` (c、既定 1) と `delta_offset_m` (δ₀、既定 0) を足し、`runner_axismach.integral_delta_r` で平滑化の後に**実効の scale (試しの r_t を含む) で** c·δ_C + δ₀/scale を当てる。値・導関数・`.spline`・出力の表を同じ変換の結果から作る。係数を当てる前の CONTUR (δ_C) も記録に残し、再較正はその δ_C を基底にする (新しい run の `delta_in_prev` は補正済みなのでそのまま基底にしない)。既定のままの YAML はビット同一。カードの `overrides` には **`deltastar_initializer.model: contur` も**書く (係数の 2 キーだけでは `resolve_integral_initializer` が ValueError; codex ①)。同じ関数を寸法の逆算 (`_sizing_delta_supplier`)・`prepare_ns`・`hand_mesh.walls` が使うので、スロート径 200 mm に合わせた r_t の解き直しも同じ δ で行う。カードからは `overrides: {deltastar_initializer.delta_scale: c, deltastar_initializer.delta_offset_m: δ₀}`。
- **NS の位置づけ**: 壁を決めた後の確認の 1 本 (δ* の反復はしない)。出口コア M が §6 V5b を外れたときだけ、その NS で同じ 2 点の較正をやり直して 1 回だけ作り直す。
- **Euler**: 新しい理想の形で回し直す (codex ④: r_t の一致の桁で安全な閾値は無い; カードは R・入口半径・各長さも r_t で割り直すので設計・格子・参照の補間の同一性は r_t の近さで保証できない)。参照の出力と sha256 を固定する。
- **初期場**: run_0005_ns の最後の場から補間 (codex ⑥: 新しい候補の補正に近い; 速い整定の保証ではなく、波が初期場由来かを短い窓で判定しない)。
- **前提の確認 (#3k、codex ⑤)**: 係数は全分布の平滑化の端点から決まるので、抽出の径方向の標本化の感度を全分布で確かめてから係数を確定する。

- **3 係数への拡張 (2026-10-15 ユーザ決定「bがよさそう」、§5.1 #3s・#3t、上位に諮る前)**: δ(x) = (c₀ + c₁·x)·δ_C(x) + δ₀ (x [m] は設計スロートから、c₁ [1/m])。2 係数 (c₁ = 0) の NS (`run_0008_ns_b0`) で比 δ_E/δ_C が x 1 m 0.47 → 8 m 0.84 と一様でなく、途中で入れた δ が NS より 1.2〜2.1 mm 厚かったため。較正: スロート (x = 0) と出口 (最後の列) の δ_E に**等式で**合わせ、残る 1 自由度 c₁ を x ≥ 0 の最小二乗で決める (モデルは c₁ に線形なので閉形式)。元は**今の壁 (2 係数) の上の NS** (`run_0008_ns_b0` res_90000、平滑化後の δ_E、係数前の δ_C = (入れた δ − δ₀)/c): ~~等重み c₀ 0.6743・c₁ 0.02149~~ → **Δx 重み (台形) の最小二乗 (codex diagnose 2026-10-15 m55-three-coef M1 採用): c₀ = 0.6779482479245829、c₁ = 0.021036335310539667 /m、δ₀ = −3.057205936466481e−4 m**。この NS の δ_E に対する誤差 (旧 NS への射影の誤差で、新しい壁の NS の誤差ではない) は最大 0.612 mm (x 1.07 m)、x ≥ 3 m は 0.198 mm (2 係数は最大 2.14 mm)。比 δ_E/δ_C が一定でないことだけでは c₁ の必要は示せない (2 係数でも c + δ₀/δ_C は変わる; 根拠は 2 点で較正し直しても残る途中の誤差、codex m2)。今の壁からの壁角の変化は最大 0.06°。較正し直しは V5b・V5c の両方を判定し、前の更新で誤差が縮み、残る形の誤差を同じ 3 係数でさらに減らせるときだけ 1 回まで (codex M2; 出口 M だけで決めない、V5c だけ落ちて形の残差が小さいときは回さない)。リミッタの基準値は run_0008 の固定値を継承する (codex M3; 3 係数化と同時に方式を変えない)。実装は `deltastar_initializer.delta_scale_slope_per_m` (c₁、既定 0) を足すだけ (2 係数と同じ場所、既定はビット同一)。本段は B0 (§5.1 #3r の既定化の前なので `main_line` を明示)。

### 4.4 NS の本段・判定・出力の新しい既定 (2026-10-15 ユーザ決定の束、上位に諮る前; §5.1 #3q・#3r・#3v・#3x・#3y)

決まったこと (ユーザ 2026-10-15): 本段はライン陰解法 B0 で十分 (#3r)、生の排除厚さを壁に直接入れる経路は生産で使わない (#3v)、収束は「残差の準定常・流量の収支・目的の物理量の準定常」で判定する (#3x)、全場を出し続けず境界の h5 で監視する (#3y)、境界層の厚みの時系列は判定から外し最後の全場 1 枚から抽出する (#3y)。以下はその実装の設計 (案)。

1. **forge の出力 (出力だけの変更、数値に触らない)**: `time.fieldOutStepInterval` を足す (既定 = `outStepInterval` で今とビット同一の出力)。全場 `res_<step>.h5` (`output.cpp:315` `outputH5_XDMF`) はこの間隔と最後の step、境界の h5 (`output.cpp:333` `outputBconds_H5_XDMF`、`outputHDFflg: 1` の境界) は今どおり `outStepInterval`。
2. **生成器の既定 (要求カードの雛形 `axismach_final_v1`・`axismach_explore_v1` の NS)**: (a) 入口・出口・壁・軸の `outputHDFflg` を 1 (今は入口・軸が 0)、(b) `outStepInterval` 2,000 (境界)・`fieldOutStepInterval` 20,000 (全場は 20 万 step で 10 枚 + 最後)、(c) 本段は B0 (`evaluate.main_line` 相当を既定に: dir・cap 50・itj 5・cfl 4・緩和 0.7)、前の段 (S1・S2) は point のまま、Euler の cfl は別に保持 (codex m55-three-coef M4)、(d) リミッタの基準値は最初の本段の開始時に forge の自動の値で決め (今の自動と同じ式: 前の段の終わりの状態の全域の体積平均)、forge の log から読んで run の記録に残し、継続・再開・次の壁の pass はその値を継承して固定する (codex M3 (b))。`main_line.limiter_ref` は省略を「自動で決めて記録」の意味にする。
3. **判定の道具 (新規 `check_settled.py`、`check_convergence.py` の後継の既定)**: 3 条件の AND。(i) 残差の全列の末尾が準定常 — log10(rms) の末尾の窓で `check_quasisteady` と同じ傾き・振幅の検定 (上昇・ドリフトは不合格、平らな床は合格; 低下桁数は見ない)。(ii) 流量の収支 — 入口と出口の境界の h5 から 2π∫(ρu dr − ρv dx) r (停止規則の `mass_flows` と同じ式) で ṁ_in・ṁ_out を測り、|ṁ_in − ṁ_out|/ṁ_in の窓の平均 ≤ 5e−4・振幅 ≤ 5e−3 (今の停止規則の r と同じ値)、両方の系列が準定常。(iii) 目的の量の準定常 — 問題ごとに登録: ノズルの既定は出口コア M (出口の h5、r/r_w 0.05〜0.7)・出口の最大の流れ角・試験部の軸の M の最大偏差と片振幅 (軸の h5; 軸の節点の癖は量の定義に書いて判定は差の時系列で見る)・壁の熱流の積分 (壁の h5)。判定は `check_quasisteady` の STEADY。どれかが判定不能なら合格にしない。
4. **停止規則 (実行中の見張り)**: 同じ 3 条件を境界の h5 と残差履歴で毎出力 (2,000 step) 判定し、保持 (今の 4 万 step) の後に止める。θ_r と δ の時系列 (全場を読む・δ の抽出を 21 回) はやめる。止めた後に最後の全場 1 枚から δ_E を抽出する (高速化済み、数十秒)。
5. **直接の δ_E の経路 (#3v)**: 要求カードの流れの既定は較正した CONTUR (2・3 係数)。`--ns-correction` は既定から外して警告。受け渡しの関門の 2 窓は直接の経路だけ。コードの削除は他の plan (isothermal-wall-chain §4.9・③ベル) と相談してから。
6. **AGENTS.md の「収束確認」の節**: 道具 (3) を作って M5.5 (run_0010) と過去の run (到達した run・未到達の run) で試してから、判定の正本を `check_settled.py` に替える文に直す。他のセッションに効くので調整板に書く。
7. **試験部の残りの傾き** (run_0010、x 4 → 8 m で −0.094 → +0.058 %) の原因は別に諮る (4.3 の 3 係数の残り 0.6 mm (x ≈ 1 m) か、別の要因か)。

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
| 3i | ~~`metrics.json` (と抽出の要約) の ṁ を一定 x の断面で積分する~~ (2026-10-15 済み) | `runner_axismach.collect` (`deltastar.massflow_ratio(run_dir, euler_ref)`、1232 行) と `deltastar.py:900` (抽出の要約の `massflow`) が `remap_constant_x` 既定 False のまま。壁法線の層の格子では列が一定 x でないので列の 2π∫ρU_x r dr は断面の流量でない (同関数の docstring の既知の型; ③ベルの経路 `bell_deltastar.py:54` だけ True)。run_0005_ns: 列のまま 120.19 kg/s (散らばり 2.1 %)・比 0.98977、一定 x 121.350 kg/s (散らばり 6.4e−5、停止規則の出口の ṁ 121.350 と一致)・比 0.99930。直し方 (案): 両方の呼び手を True にする (列が一定 x の格子では np.interp が節点の値をそのまま返すので値は不変 — 試験で確かめる)。判断: 2026-10-15 codex diagnose (`notes/reviews/2026-10-15-m55-p1-axis-wave-diagnose.md`) — 判定の緩和でなく積分断面の訂正として採用。両方の呼び手で `remap_constant_x=True` を明示し、列が一定 x の既存格子 (Euler の run) で値が変わらないことを試験で確かめる **実装 (2026-10-15)**: `runner_axismach.collect` は `massflow_ratio(..., remap_constant_x=True)`、`deltastar_from_core_matched_euler` の帳簿は抽出と同じ `remap_constant_x` を渡す (旧い呼び方 = remap なしは今までどおり)。③の `extract_bell` の `massflow_column` は列のままを明示して測る。確かめたこと: `tests/run_massflow_constant_x_tests.py` ALL PASS (列が一定 x の格子でビット同一・傾いた格子で一定 x の比が列より真値に近い・抽出の帳簿の断面)、`run_bell_deltastar_tests.py` FAIL 0 (#3i の意図した変更は帳簿だけ — 帳簿以外はビット同一と明示に直した)、`run_deltastar_tests.py`・`run_deltastar_next_columns_tests.py`・`run_deltastar_root_zone_tests.py`・`run_handoff_diag_mark_tests.py` PASS、case/45 run_0038/0037 (列が一定 x) で一定 x の帳簿がビット同一 (ṁ 98.97772 kg/s・比 1.0002256) | O |
| 3j | ~~成果物の寸法の関門の `straight_end_in_mesh_interval` が滑らかな接合で誤って落ちる~~ (2026-10-15 済み) | `export/dimension_gate.py` `mesh_dimensions` は入口半径との差 ≤ `NODE_EQ_MM` (1e−6 mm) の節点を直管と数える。run_0005_ns の格子は接合 (x −500.000 mm) の 0.08 mm 下流に節点 (x −499.920 mm、r − 300 = −7.3e−9 mm) があり直管と数えたので、区間 [−499.920, −498.591] が接合を含まず FAIL。カードの 5 項目と弦の検査 (最大 0.00076 mm) は PASS。判断: 2026-10-15 codex diagnose — 最後の 1 点だけの一致で置き換える案は却下 (他の節点の不整合を見逃す)。**接合の位置は STEP の曲線から決め、メッシュの壁の節点の x で STEP を評価し、直管側から接合を越える区間全体の半径を照合する** (許容は `NODE_EQ_MM`)。カードの寸法・弦の条件は維持し、接合をずらした負例 (接合を ±数節点動かした曲線) を拒否することを試験に入れる **実装 (2026-10-15、implementer)**: `dimension_gate.junction_from_mesh` (壁節点の x で STEP の曲線を Newton の逆引きで評価、入口端〜設計スロート x = 0 の全節点を照合) と `evaluate_gate(junction=...)` の `straight_junction_pass` (最大差 ≤ `NODE_EQ_MM`・照合した接合 = STEP の直管の終端・x = 0 まで照合・接合の上流と下流に節点・定義域の外 ≤ `END_TOL_MM`) を必須に、旧の `straight_end_in_mesh_interval` は `reference_only` の記録に。`deliverables.part_dimensions` が照合を渡す。確かめたこと: case/71 run_0005_ns の実物でゲート全体 PASS (照合 1,653 節点・最大 7.9e−12 mm)、M6 の fixture (msh 10 桁) 最大 8.5e−8 mm で PASS、負例 (接合を ±1 節点の間隔ずらす → 0.51〜1.59 mm で FAIL、直管・接合の前後の節点の半径を 1e−5 mm 動かす → FAIL、照合が欠ける → 未判定、別の位置・x = 0 の手前までの照合 → FAIL)、`run_deliverables_tests.py` ok 68・FAIL 0 (`FORGE_DELIV_FIXTURE` で fixture を作り直して渡した; 既定の fixture のパスは今は無い)、`run_mesh_structure_tests.py` ok 94・FAIL 0 | O |
| 3k | ~~**抽出の径方向の標本化の A/B (forge 0 step、事前登録 2026-10-15、codex diagnose の判別)**~~ (2026-10-15 済み: 結果 B・2 係数はすべて以内) | 変えるのは `band_local_deficit` の `refine` だけ: A = 4 (今)、B = 32。対象は pass 1 の壁の元 `run_0003_ns_cont/res_200000.h5`、参照 `run_0001_euler/res_12000.h5` (sha256 93cb88b0…)。固定: 一定 x への補間 (`remap_constant_x`、`remap_edge` missing)、各列の保存済みの `band_y_b` (`y_b_fixed`、帯の梯子は選び直さない)、Euler の補間、平滑化のノット 2 r_t・λ 1・重み。手順: まず A が保存済みの raw (`_runs/run_0005_ns/judge/diag_wave/pass0_run_0003_delta_r_equiv.csv`) を再現することを確かめる。判定 (x 1〜3.5 m): raw の段差、平滑化後の δ_E、補正の角 atan(r′_design + δ′_E) の差。**許容 max|Δδ_E| = 0.01 mm、max|Δθ| = 0.01°**。結果 A (どちらかを超える) → 「今の抽出は径方向の標本化に十分独立」を棄却し抽出を直す。結果 B (両方以内) → 固定帯の標本化を壁角の窪みの主要因から外す (帯の選び方の感度は除外しない)。どちらでも軸の波の原因の確定にはしない。AWS A (CPU、保存元)。やらない: 物理壁の Euler 1 本で形状/境界層を二分する・抽出の検査の前に λ を上げる・ω だけ変えて pass 2 を投入する **拡張 (codex diagnose 2026-10-15 ⑤、2 係数の較正の前提)**: 同じ A/B で、各腕の全分布を平滑化し、同じ x = 0・出口で §4.3 の較正をして c・δ₀ を出し、δ_C (pass 0 の `delta_in_prev`) と設計壁を固定して 2 係数の補正 c·δ_C + δ₀ と壁角を作る。追加の許容: **端点の max|Δδ_E| ≤ 0.01 mm、下流全域 (x ≥ 0) の max|Δ(c·δ_C + δ₀)| ≤ 0.01 mm、2 係数の壁の max|Δ壁角| ≤ 0.01°**。いずれか超過 → 「2 点の較正なら抽出の感度は無視できる」を棄却し抽出を先に直す。すべて以内 → 標本化の依存を新方式の阻害要因から外して #3n へ (元の直接 δ_E の更新の原因の診断は未決着のまま)。hard の判定が腕で変わったら単純な比較を合格にしない。A が保存の raw を再現することを必須にする **結果 (2026-10-15、AWS A、台本 `_runs/run_0005_ns/refine_ab_3k.py`、原文 `_runs/run_0005_ns/judge/refine_ab/`)**: A が保存の raw・平滑化後の δ_E・hard をビット同一で再現、hard は両腕で同じ。元の判定 (x 1〜3.5 m): max|Δδ_E| 0.0087 mm (x 3.13 m)・max|Δθ| 0.0018° → **結果 B** (固定帯の径方向の標本化を壁角の窪みの主要因から外す; 帯の選び方の感度は除外しない、窪みの原因は未決着)。raw の隣接の段の最大は 0.263 → 0.055 mm (refine 32 で小さくなるが平滑化の後は効かない)。2 係数: c 0.852073 / 0.852110、δ₀ −0.32201 / −0.32209 mm、端点の |Δδ_E| スロート 8e−5・出口 0.0017 mm、x ≥ 0 の max|Δ(c·δ_C + δ₀)| 0.0017 mm、壁角 1.6e−5° → **すべて以内** (標本化の依存を 2 係数の阻害要因から外して #3n へ) | O |
| 3l | 軸の波の時系列の準定常 | codex diagnose の Major (「過渡ではない」は 2 枚の比較では未立証)。`run_0005_ns` の軸の M の最大絶対偏差・半振幅・極値の位置を各全場から作り `check_quasisteady --series-csv` に渡す。**A の途中の全場は判定の後に消した**ので、残る 5 枚 (100k・140k・160k・180k・200k) で行い、点数が足りなければ「判定不能」と書く | O |
| 3m | 次の壁の作り方の選択 (ユーザ判断待ち) | (a) 直接の δ* (今の経路) を #3k の結果で抽出を直してから pass 2、(b) CONTUR の積分式を NS の δ* で較正した壁 (ユーザ 2026-10-15 の問い「conturの積分式を補正して壁作ってるか？」; チェーンには `deltastar_initializer` の k_f の較正がある、M6 は 10-05 に出口だけの較正の壁が滑らかさ・波・オーバーシュートの基準を通った)、(c) pass 0 または pass 1 を暫定の成果物で止める。試験部の合格条件 (codex 提案: x 4.0〜7.95 m の最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %) もユーザに聞いてから事前登録する **ユーザ 2026-10-15**: 「いつもこうしてなかったか？ 摩擦係数でなくて、倍率かけて合わせ込むほうが、割り切っててよい気はしてます」。事実: M6 の生産 (case/45、10-07) は C2 方式 (CONTUR の摩擦係数倍率 k_f = `deltastar_initializer.cf_scale` を出口の δ_C = δ_E に合わせ、r_t を解き直す; `procedures/nozzle-design-workflow.md` 4.) で、NS の δ_E を直接壁に入れていない。case/71 の pass 1 はカードの `--ns-correction` (直接の δ_E の表、10-13 に足した経路) を使った = M6 の作り方から外れた。2026-09-04 の生産レシピの記録にも「反復ごとの壁変化が特性線で末端の軸に届き ±0.5 % の軸の波を作る」とあった。pass 0 の比 δ_E/δ_CONTUR (cf_scale 1): x 1 m 0.47・2 m 0.67・3 m 0.71・4 m 0.76・6 m 0.80・8 m 0.85、スロート CONTUR +0.09 mm 対 NS −0.24 mm (符号が逆)。案 (b′): CONTUR の δ(x) に一様な倍率 c を掛けて出口で合わせる (コードに出力の倍率は無い — `deltastar_initializer` にキーを足す)。§4 の方式の変更なので、決まったら §4 を書いて codex に諮ってから実装 **ユーザ 2026-10-15「コンターの倍率だと、そのマイナスの厚み再現できないよな」**: そのとおり (正の倍率では CONTUR のスロートの +0.09 mm を NS の −0.24 mm にできない)。出口の δ で倍率を決めると、スロートの差 0.32 mm (半径) が実効のスロート面積で約 0.64 % → 出口 M で約 −0.1 % (面積比の感度 ΔA/A ≈ 5〜6 ΔM/M、γ 1.3 前後の概算) になる。案: (b″) 倍率 c + 全長で一定のずらし δ₀ (2 係数、ずらしは壁の傾きを変えないので波を作らない)、pass 0 の NS のスロートと出口に合わせると c ≈ 0.852・δ₀ ≈ −0.32 mm (概算); (b‴) 倍率 1 つを出口の δ でなく出口 M に合わせる (スロートの差を出口で吸収)。どれを採るかはユーザ判断 **判断: 2026-10-15 ユーザ「1かええな」 — 案 (1) 倍率 c + 一定のずらし δ₀ (§4.3)**。§4.3 を codex に諮ってから #3n・#3o へ | F |
| 3n | ~~§4.3 の実装 (`deltastar_initializer.delta_scale`・`delta_offset_m`)~~ (2026-10-15 済み) | `runner_axismach.integral_delta_r` で平滑化の後に δ = c·δ_C + δ₀/r_t (r_t 単位) を当て、`init_info` に c・δ₀・当てる前のスロート・出口の値を記録。既定 (1, 0) はビット同一 (試験)。カードの `overrides` で渡せること、寸法の逆算・`prepare_ns`・`hand_mesh.walls` が同じ値を使うこと (スロート径 200 mm に解き直した r_t で物理スロートが 200 mm) を試験で確かめる。カードの overrides に `model: contur` も入れる。物理壁の大域最小で寸法を解くこと・上流 poly の C² 接続を確かめる。**#3k (拡張) の結果の後** (codex diagnose 2026-10-15) **実装 (2026-10-15)**: `runner_axismach.integral_delta_r` で平滑化の後に `_two_coef` (正の有限の c・有限の δ₀、既定 (1, 0) は None) を当て、係数前の δ を `delta_r_contur` に、`init_info.two_coef` に係数・scale・係数前のスロートと出口の値を残す。寸法の逆算の記録 (`_sizing_delta_supplier`) にも係数を出す (使うときだけ)。係数は CONTUR の関数そのものの値で取り直した (codex ①「係数前の CONTUR に固定」): pass 0 の r_t で CONTUR スロート 0.091630・出口 49.32933 mm、δ_E −0.243965・41.71089 mm → **c = 0.852087987522164、δ₀ = −3.220418618084591e−4 m** (格子から測った delta_in_prev で取った値との差 1.5e−5・3e−5 mm)。カード `card_m55_2coef.yaml` (card_m55.yaml + overrides の 3 キー)。確かめたこと: `tests/run_two_coef_initializer_tests.py` ALL PASS (既定の明示はビット同一・表と関数の値の式・導関数 = c·δ_C′・記録・寸法の逆算の δ_r 関数 = prepare_ns の経路・不正な値の ValueError・既定の CONTUR と pass 0 の delta_in_prev が 8.4e−4 mm で一致)。関連の試験: `run_pw_upstream_poly_tests`・`run_wall_single_bspline_tests`・`run_mknot_v2_card_tests` FAIL 0、`run_physical_upstream_tests`・`run_euler_h_ref_tests`・`run_h370_prelaunch_tests` ALL PASS、`run_wall_temperature_tests` は 4 FAIL (V1d 2・V4 2: `--ic-from` の場の化学種の記録の検査 FORGE_ALLOW_UNVERIFIED_SPECIES で止まる) だが**変更前の 330edb5a でも同じ 4 FAIL** (今回の変更と無関係、環境・fixture 側) | O |
| 3o | §4.3 の壁で M5.5 を作り直す (Euler → NS → 成果物) | カード (`overrides` に c・δ₀) → `hand_mesh.py build` (#3e と同じ係数) → Euler (新しい理想の形で回し直す、参照の出力と sha256 を固定) → NS `--init-integral`、初期場は run_0005_ns の最後の場から補間、`--stages full`、精度は既定 (float)、出力 2,000、停止規則 (窓 2 万・保持 4 万・上限 20 万) → §6 V5 の判定 → 成果物。V5 の許容はユーザ了承済み (2026-10-15)。#3n の後 **結果 (2026-10-15、原文 `_runs/run_0007_ns/judge/`・`_runs/run_0008_ns_b0/judge/`)**: カード r_t 0.100243850935207 m (寸法 ok)、格子 NS 442,520 点 (PASS・AR 緩和 違反 0)・Euler 466,570 点、`run_0006_euler` 完走 (軸の出口 M 5.49977、参照 res_12000 sha256 cbd236ec…)。`run_0007_ns` (point cfl 1.0): 上限 200,000 で停止規則に届かず (D 7.2〜8.6e−4; 後処理をユーザ了承で止めたので STOP_RULE_STATE.json は無い = CENSORED と記録から判断)。本段の残差 `check_convergence` (本段の履歴だけ) PASS (3.8〜4.5 桁)。系列 DRIFTING (r だけ、θ_r は単調減少・漸近 −0.09〜−0.11 %)。壁解像 PASS (最大 0.999)。非有限 0、P0 > Pt 47,111 点・最大 +0.464 %。**V5b FAIL (出口コア M 5.503168、許容 5.49725〜5.50275)、V5c FAIL (r/r_w = 0.1: 最大偏差 0.224 %・片振幅 0.079 %; 軸 0.224 %・0.091 %)**。試験部の M は波でなく、x 4.4 m で +0.224 % の山から出口 +0.067 % へ単調に下がるなめらかな分布 (平均 +0.147 %)。pass 1 の波 (片振幅 0.465 %) は消えた。`--segment` は本段の記録が書かれる前に回したので前の段 (0〜5,999) を判定した (無効、本段の履歴で回し直した)。V5a は未判定 (成果物を作っていない) | O |
| 3p | **本段の陰解法の比較 point / B0 (§6 V6)** (ユーザ 2026-10-15「line implicitを使う方法でどうなるか、も同時に確かめてほしい」) | `run_0008_ns_b0` を AWS B で: 問題 `problem_m55_p2_b0.yaml` (= `problem_m55_p2.yaml` の evaluate の cfl_main を外して `main_line` {dir, cap 50, itj 5, cfl 4, relax 0.7, limiter_ref = run_0007_ns の本段の自動の値 ro 0.5810042584・p 186732.713・a 341.8295597} を足した写し)、台本 `_runs/run_0008_ns_b0/p3_ns.py` (p2_ns.py の写し、前の段の設定 `solverConfig_prestage.yaml` にも出力 2,000 を当て、本段に B0 のキーがあることを検査)。A のバイナリ (sha256 3a9e347b…) を B に写して同じものを使う。run_0005_ns の最後の場・run_0006_euler を B に写す。判定は §6 V6 と V5 (judge_0007.sh と同じ一式) **結果 (2026-10-15)**: `run_0008_ns_b0` (B0、AWS B) は本段 50,000 で停止規則の条件を満たし 90,000 で到達 (保持 4 万)、D は 7〜8 万 step で 1e−7 台。本段の残差 PASS (3.0〜5.3 桁)、系列 ALL STEADY、V5c の系列 ALL STEADY、壁解像 PASS。V5b FAIL (5.503217)・V5c FAIL (0.225 %・0.079 %) は point と同じ。**V6 (止まった場の差)**: 出口コア M +0.0009 % (許容 0.01 %、内)、V5c の 2 量の差 0.0008・0.0003 % (内)、θ_r (同じ帯で測り直し) −0.12〜−0.15 % (許容 0.1 %、**外**; point は CENSORED で θ_r が単調減少中、その漸近値なら 0.05 % 以内)、δ_E は未測定 (ユーザ「もうええやろ、line implicitで十分だよ」で比較を打ち切った)。残差は B0 が約 3 万 step・約 12 分で床、point は 14 万 step・34 分 (1 step B0 24.9 ms・point 14.3 ms、A・B は別の機械)。**判断 (ユーザ 2026-10-15)**: line implicit で十分 → B0 を採用 (#3r) | O |
| 3q | **δ の時系列の抽出を要るときだけにし、速く・道具にする** (ユーザ 2026-10-15「ここらへんほんまにいるんかいな」「抽出するとしても高速化・ツール化をきっちりやってほしい」) | 事実: 停止規則の見張りは forge の終了の後に 2 窓 (δ_E・δ_target、受け渡しの関門 §4.9 用) の全場 21 枚から δ* を 1 コアで抽出し、run_0007_ns・run_0008_ns_b0 で 30 分以上かかった (その間ドライバが終わらず stage_manifest の本段の記録も書かれない → `check_convergence --segment` が前の段を判定した)。2 係数の流れ (§4.3) では停止に δ の時系列は要らず、再較正にスロートと出口の δ_E (最後の出力の 1 枚) だけが要る。方針 (案、§4 を書いて諮る): (1) δ の時系列は直接の δ_E の更新 (受け渡しの関門) を使う流れだけで有効にし、2 係数の流れでは既定で切る、(2) 抽出そのものを測ってから速くする (列ごとの並列・edge の梯子のベクトル化など)、(3) run の外から 1 コマンドで回せる道具にする (最後の場・指定の出力・時系列) **ユーザ 2026-10-15「抽出1回は約3.7分で ささっと早くしてくれんか」**: 先に (2) 抽出の高速化をやる (プロファイル → 断面ごとに Euler の参照の補間を 1 回にする・断面を並列に; 出力は今とビット同一であることを試験で確かめる)。(1) の切り替えは §4 を書いてから **(2) 実装 (2026-10-15)**: プロファイル (B、run_0008_ns_b0 res_90000、cProfile 331 s): `band_local_deficit` 57.2 万回 (edge の梯子が 1 断面あたり約 120 の帯の位置)、うち Euler の補間 54 s・lstsq 48 s・台形 40 s・細分の補間 19 s が帯に依らない量の作り直しを含む。直し: `band_local_pre` (細分の点・NS と Euler の値・内側コアの倍率を断面ごとに 1 回、`_pre` で渡す; 帯の位置を固定する経路だけ) と `_map_columns` (梯子を列ごとに fork の子で並列、`FORGE_DSTAR_WORKERS`、既定 min(8, CPU − 1)、他のスレッドがある・デーモンの中なら逐次)。**結果 (手元 12 コア、同じ入力)**: 元 113.8 s → 前計算 66.5 s → 前計算 + 8 並列 15.9 s、出力 4 ファイル (equiv・next・edge・summary) の sha256 は 3 つとも同一。AWS B (混んでいる 3 並列) 124 s でも 5 ファイルが元とビット同一。試験: `tests/run_dstar_parallel_tests.py` ALL PASS (前計算あり = なし、並列 4 = 逐次 1 の全配列、スレッドがあるとき逐次)、`run_bell_deltastar_tests.py` の (4) 基準コミットとのビット同一 (合成の edge・実データ) OK、`run_deltastar_tests`・`run_deltastar_root_zone_tests`・`run_deltastar_next_columns_tests`・`run_handoff_diag_mark_tests`・`run_massflow_constant_x_tests` PASS。さらに: 帯の候補をまとめて配列で計算すればもう 1 桁縮むがビット同一でなくなる (許容を決めて確かめる別の項目) | O |
| 3r | **ライン陰解法 B0 を要求カードの雛形の既定にする** (ユーザ 2026-10-15「line implicitを規定化してほしい」) | 根拠 (§6 V6、同じ条件の 2 腕): B0 は本段 5 万 step で停止規則の条件を満たし 9 万 step で到達、point (cfl 1.0) は 20 万 step で上限 (D 7〜9e−4)。残差は B0 が約 3 万 step・約 12 分で床、point は 14 万 step・34 分 (1 step は B0 24.9 ms・point 14.3 ms)。止まった場の差 (同じ帯): 出口コア M 0.0009 %、θ_r −0.12〜−0.15 % (point は CENSORED で θ_r が単調に減少中)。方針 (案、§4 を書いて諮る): `axismach_final_v1`・`axismach_explore_v1` の雛形に `evaluate.main_line` (dir・cap 50・itj 5・cfl 4・relax 0.7) を既定で入れる。リミッタの基準値 (`limiter_ref`) は今は問題ごとの値を明示する必要があるので、自動で決める方法 (前の段の終わりの場から forge と同じ式で出す等) を設計する。SERN (⑤) は cfl 1.0 が生産 (CFL 4 で発散) なので対象外。生成器の既定の全体は plan `tooling-generator-defaults-and-launch-gate` (作業ツリー `forge-gate`) の段取り 2 と整合させる **リミッタの基準値の自動 (ユーザ「これは大事かも」)**: forge の自動 (`main.cpp:1775-1825`) は開始した場の全域の体積平均 |ρ|・|P|・|a| で、開始場に依存する作用素になる (再開のたびに変わる; そのため B0 は値を固定して書く)。案: 問題そのものから決まる値 (設計した形の 1 次元等エントロピー解に同じ体積平均の式を当てる等) を準備のときに計算して書く — 初期場によらず同じ問題なら同じ値。§4 で比べて諮る **codex (m55-three-coef) M3・M4**: 基準値は初回の本段の開始時に一度決めて以後継承 (値・由来・入力の場の hash を記録、`stages=none` は既存値を継承; forge と同じ値にするには境界条件・壁温拘束の後の状態と `CELLS/volume` で平均)、既定化は NS 本段だけ (Euler の cfl は別に保持) | F |
| 3s | **2 係数の壁の残りの偏り (V5b・V5c の FAIL) の扱い** | 事実 (`run_0008_ns_b0` の res_90000 から抽出、run の外 `~/m55_diag_0008/extract` on B): スロートの δ_E は入れた δ と +0.0003 mm、出口は −0.27 mm、**途中は入れた δ が NS より 1.2〜2.1 mm 厚い** (x 1 m −1.20・3 m −1.92・4 m −2.01・5 m −2.01・7 m −1.29 mm)。比 δ_E/δ_C は x 1 m 0.47 → 3 m 0.72 → 8 m 0.84 と一様でない。試験部の M の過大 (平均 +0.147 %) は途中の厚すぎ (実効の流路が開きすぎ) と向きが合う (原因の確定ではない)。**同じ 2 点で機械的に再較正しても** c 0.84656 (−0.65 %)・δ₀ はほぼ同じで出口の 0.27 mm しか変わらず、途中の 2.0 mm は残る (codex diagnose 2026-10-15 ③ の「端点差と確認せずにもう一本回さない」に当たる) → 回さない。候補: (a) 今の結果を暫定の成果物にする、(b) 3 係数 (c₀ + c₁·x)·δ_C + δ₀ (なめらかさは保つ; この NS の δ_E に最小二乗で途中の誤差 最大 0.47 mm・rms 0.19 mm、スロート −0.19・出口 −0.02 mm)、(c) 最小二乗の 2 係数 (途中は良くなるが出口 −1.76 mm)。**判断: 2026-10-15 ユーザ「bがよさそう」 — (b) 3 係数 (§4.3 の拡張)**。抽出 1 回は約 3.7 分 (B0 の res_90000、1 コア; 停止規則の後処理は 21 回で 1 時間強) | F |
| 3t | 3 係数の実装と M5.5 の作り直し (Euler → NS B0 → 判定 → 成果物) | `runner_axismach.integral_delta_r` の `_two_coef` を c₁ (`delta_scale_slope_per_m`) に広げ、δ = (c₀ + c₁·x_m)·δ_C + δ₀ (x_m = x·scale)。試験: 既定ビット同一・式・導関数 (c₁ の項の導関数 c₁·δ_C + (c₀ + c₁x)·δ_C′)・寸法の逆算に同じ値・不正な値の ValueError。カード `card_m55_3coef.yaml` (overrides に model・delta_scale・delta_scale_slope_per_m・delta_offset_m)、`main_line` (B0) を明示、初期場は `run_0008_ns_b0` の最後の場。判定は §6 V5 (と同じ台本)。諮問済み (m55-three-coef、全件採用)。**実装 (2026-10-15)**: `_slope_coef` (c₁)、`tests/run_two_coef_initializer_tests.py` (5) ALL PASS (式・導関数 6.6e−13・c₁ = 0 はビット同一・寸法の逆算・不正値)。カード `card_m55_3coef.yaml` (Δx 重みの係数)、台本 `_runs/run_0010_ns_b0/p4_all.sh`・`p4_ns.py` (AWS B; Euler `run_0009_euler` → NS `run_0010_ns_b0`) **結果 (2026-10-15、AWS B、原文 `_runs/run_0010_ns_b0/judge/`)**: カード r_t 0.1002434996788185 m (949 s、うち M_knot 682 s = 取り込みの前のコード)、格子 PASS、`run_0009_euler` 完走 (参照 res_12000 sha256 2b2d02ce…)。`run_0010_ns_b0` (B0、初期場 run_0008 の最後の場) は本段 46,000 で停止規則の条件・86,000 で到達、系列 ALL STEADY、V5c の系列 ALL STEADY、壁解像 PASS、非有限 0。**V5b FAIL (出口コア M 5.502841、許容 5.50275 を 0.0016 % 超え)、V5c: 最大偏差 0.0940 % PASS・片振幅 0.0761 % FAIL**。試験部の M (r/r_w = 0.1) の平均 −0.002 % (2 係数は +0.147 %)、x 4.0 m の −0.094 % から出口の +0.058 % へ単調に上がる傾き (波なし)。δ の誤差 (新しい壁の NS の δ_E に対し) 最大 0.624 mm @ x 1.06 m・x ≥ 3 m 最大 0.257 mm (2 係数 2.14 mm)。**V7: δ の誤差は半分以下・最大偏差も半分以下 (0.2247 → 0.0940 %)・片振幅は 0.0788 → 0.0761 % → 事前の条件 (両量とも半分以下) を満たさないので厳密には結果 A でない** (平均の過大は途中の厚すぎが主因だったが、傾きは別の要因で未特定)。同じ規則の較正し直しの効果は小さい (最大 0.624 → 0.612 mm・x ≥ 3 m 0.257 → 0.214 mm) → 規則 (形の誤差をさらに減らせるときだけ) に当たらず回さない。`check_convergence` (本段) は NOT CONVERGED (収束した前の解から始めたので本段の始めの残差が小さく 2.3〜3.0 桁; 末尾は run_0008 の床の 0.62〜1.27 倍)、`--from-floor run_0008` も過渡の山で NOT CONVERGED (→ #3x の判定の変更へ)。後処理 (δ の時系列) は止めた。抽出は新しいコード (3 並列) で 71 s。**成果物 (手元、`run_0010_ns_b0/deliverables/`) PASS**: cad・mesh・fluent・wall_compare・excel・寸法 (物理スロート径 200.000072 mm @ x −0.231 mm、L_div 8,000.40 mm) すべて合格。設計報告 `run_0010_ns_b0/report/` (図 9 枚・pptx 14 枚)。成果物の格子の品質の確認が B0 の本段の設定 (node 専用の `implicitThermalJacobian`) で変換器に拒まれたのを、B0 を当てる前の設定 `solverConfig_prestage.yaml` を使うように直した (`deliverables._mesh_quality`) | O |
| 3u | ~~**成果物に「理想の壁 (Euler・設計壁) と物理壁の座標の比較」の表と図を足す**~~ (2026-10-15 済み) (ユーザ 2026-10-15「オイラーで作った壁と、物理壁の座標をくらべたデータとプロットも、最終的な成果物としてほしい」) | `export/deliverables.py` の成果物一式に: 同じ x [mm] の点で設計壁 (wall_repr.json の design) と物理壁 (physical) の半径 [mm]・差 (= 物理 − 設計、境界層の分) の表 (CSV と Excel のシート) と図 (両方の壁の形・差の分布・スロート付近の拡大) を足し、MANIFEST に載せる。x の点は `wall_shape_mm.csv` と同じ刻み (粗 1 mm・スロート付近 0.1 mm)。設計壁と物理壁の定義域の違い (物理壁は上流の多項式・設計壁の範囲) を表に書く。試験: 既存の成果物の試験が通ること、差が prepare の δ の記録 (スロート・出口) と一致すること **実装 (2026-10-15、implementer)**: `export/wall_compare.py` (`compare_table`・`delta_record_check`・`plot_compare`・`export_compare`) と成果物の part `wall_compare/` (CSV `wall_design_vs_physical_mm.csv`: x・r_design・r_physical・Δr・両方の壁角、片方の定義域の外は空欄、注記に区間と Δr の意味 (x ≥ 0 = 境界層の補正、縮流部 = 上流の 2 曲線の形の差で境界層でない、直管 = r_U の差); 図 3 枚 (全長・Δr・スロート付近 −300〜500 mm); JSON; xlsx のシート `wall_compare`; README の行)。δ の記録との照合: 出口は 1e−6 mm、スロートは表 (`delta_r_csv`・`delta_r_initial.csv`) があれば壁と同じ 5 次補間で照合、表が無い積分法の run は未照合 (記録の `delta_r_throat` が 1500 点の直線補間なので 1e−5 mm 級で合わない; 許容は広げない)。**判断 (2026-10-15 主セッション)**: 未照合の扱い・置き場所 `wall_compare/` は implementer の既定どおり採用。確かめたこと: `run_deliverables_tests.py` ok 89・FAIL 0、`run_mesh_structure_tests.py` FAIL 0、`run_wall_tables_tests.py` ok 62・FAIL 0 (形状表はバイト同一、`part_excel` 単体の出力は変更前と同じ)。run_0005_ns の実物: スロートの Δr −0.24462 mm・出口 41.823 mm、照合 PASS、一括書き出し PASS | O |
| 3v | **NS の生の排除厚さを壁に直接入れる経路を生産から外す** (ユーザ 2026-10-15「ns計算で、生の排除厚さ入れるやり方、もう取る意味ないよね」) | 根拠: pass 1 (直接の δ_E、`wt_card --ns-correction`) で試験部の軸に ±0.6 % の波 (`run_0005_ns`)、2 係数・3 係数の較正した CONTUR では波が出ない (`run_0007_ns`・`run_0008_ns_b0`)、2026-09-04 の生産レシピの記録 (反復で ±0.5 % の軸の波)。抽出 (δ_E を測ること) は較正の材料として残す。方針 (案、§4 を書いて諮る): (1) 要求カードの流れの既定を「較正した CONTUR (係数を NS の δ_E から決める)」にし、`--ns-correction` (直接の表) は既定から外して使うと警告、(2) 停止規則の δ の 2 窓の時系列は直接の経路だけで有効 (#3q (1))、(3) コードの削除は、直接の経路を使う他の plan (`tooling-nozzle-isothermal-wall-chain` §4.9 の受け渡しの関門・③ベル など) と範囲を確かめてから (他のセッションの plan なので相談する) | F |
| 3w | ~~(低い優先度) 物理壁の最小径を設計スロート x = 0 に置く選択肢~~ (2026-10-15 取りやめ) (ユーザ 2026-10-15「まあそれなら曲率ゼロにしてもいい気がいたしましたが。つよい要望でもありませんが」; 「x = 0 で傾き 0」の意味と受け取った、要確認) | 事実: 今は x = 0 で物理壁の傾き = δ′(0) = +6.8e−4 (run_0005_ns)・曲率 2.0e−3 /mm なので最小径が x ≈ −0.34 mm (半径の差 1.15e−4 mm)。NS の δ_E も x = 0 で傾き +7.3e−4・最小は x ≈ −105 mm (run_0008_ns_b0) なので、今の形は NS に合っている。上流の多項式 (`_init_poly`) は x = 0 の物理壁の (r, r′, r″) を引き継ぐだけで、x < 0 の差は境界層でなく 2 曲線の形の差 (−300 mm で NS と 0.19 mm 違うが亜音速なので流れにほぼ効かない)。案: δ の x = 0 の傾きを 0 にする (スロート近くの短い区間でなめらかに寝かせる、C² は保てる) — 実効のスロートが約 0.35 mm 下流にずれる。3 係数の run の後に必要なら **決着 (2026-10-15 ユーザ「ごめん、やっぱやめて。あなたの言ったとおりにしましょう」)**: やらない。今の形 (上流の多項式が x = 0 の物理壁の値・傾き・曲率を引き継ぎ、最小径は NS の境界層の成長どおり x ≈ −0.3 mm、実効のスロートは設計の x = 0) のまま | O |
| 3x | **収束の判定を「残差の準定常・流量の収支・目的の物理量の準定常」に変える** (ユーザ 2026-10-15「落ちた桁数…近々で、これで語るのやめにしようか。残さが準定常になり、流量バランスが取れ、目的とする物理量が準定常担ったら収束と判定しないか」) | 事実: `check_convergence.py` は系列の最大からの低下桁数 (既定 3 桁) で判定するので、収束した前の解から始めた run (`run_0010_ns_b0`: 本段の始めの残差が小さく 2.3〜3.0 桁、末尾は run_0008 の床の 0.6〜1.3 倍) が NOT CONVERGED になる。`--from-floor` も過渡の山で落ちる。方針 (案、§4 を書いて諮る): 判定を 3 条件の AND にする — (1) 残差の全列の末尾が準定常 (`check_quasisteady` と同じ傾き・振幅の検定、上昇・ドリフトは不合格)、(2) 入口と出口の流量の差 |ṁ_in − ṁ_out|/ṁ ≤ 許容 で準定常、(3) 目的の量 (出口コア M・試験部の一様性・流れ角・軸の M など、問題ごとに登録) が準定常。AGENTS.md の「収束確認」の節と `check_convergence.py` の既定を替えるのは全セッションに効くので、道具を作って試してから AGENTS.md を直す (plan `tooling-convergence-and-wall-resolution-gates` と整合) | F |
| 3y | **全場を出し続けずに、境界の h5 で監視・判定する** (ユーザ 2026-10-15「入口出口の流量差や、断面平均のマッハ数なりを定期的にprobeとして出力する機能欲しい」「出口境界と入り口境界と軸中心のh5を出し続けてくれてもいいな…そこから後処理的にマッハ数・流れ角出したりとか」「結果丸ごと出し続けて容量圧迫されるのがな」) | 事実: 境界の出力は `bcondConfig.yaml` の `outputHDFflg` で境界ごとに選べる (`output/output.cpp:333` `outputBconds_H5_XDMF`、`res_<名前>_<physID>_<step>.h5`)。M5.5 は出口 (2)・壁 (3) が 1、入口 (1)・軸 (4) が 0。ただし境界の出力は全場と同じ `outStepInterval` (全場 1 枚 約 100 MB で 20 万 step・2,000 ごとに約 10 GB)。停止規則の θ_r は全場から測っている。方針 (案、§4 を書いて諮る): (1) 入口・軸の `outputHDFflg` を 1 に (生成器の既定に)、(2) forge に全場の出力の別の間隔 (例 `time.fieldOutStepInterval`、既定 = outStepInterval で今と同じ) を足す — 出力だけの変更、(3) 停止規則・判定 (#3x) の量を境界の h5 から取る (入口・出口の流量、出口の M の分布・流れ角、軸の M、壁の量)、(4) 全場は粗い間隔と最後だけ **ユーザ 2026-10-15「そのプローブから見る方式だと、境界層厚さの時系列チェックはできなくなるからご認識あれ。それは私は別にいいと思ってます。最終だけあればね」**: 了承 — 境界層の厚み (θ_r・δ) の時系列は判定から外し、最後の全場 1 枚から抽出する (#3q の (1) とも一致)。全場は最後 (と粗い間隔の数枚) だけ | F |
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
- **V5 M5.5 の 2 係数の壁 (§4.3) の NS** (事前登録 2026-10-15; V5b・V5c の数値は**ユーザ了承 2026-10-15「よい」**、結果を見る前に固定): 判定区間は本段 (`--segment`)。
  - V5a 寸法: 成果物の寸法の関門の全体 (カードの 5 項目・メッシュと STEP の整合・弦) が PASS (#3j の直しの後に判定; カードの 5 項目だけで合格にしない、codex 2026-10-15 ③)。
  - V5b 出口コア M (停止規則の量): 5.5 ± 0.05 % = 5.49725〜5.50275 (M6 の生産のゲートと同じ幅)。外れたら §4.3 の再較正を 1 回。
  - V5c 試験部の一様性 (codex diagnose 2026-10-15 の提案をユーザ了承): x 4.0〜7.95 m の M の最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %。軸の節点は既知の値の癖 (静圧 0.5〜0.7 % 低) があるので r/r_w = 0.1 の線で測り、軸の節点の値は併記する。窓の中で直接 max|100(M/5.5 − 1)| と 100(max M − min M)/(2·5.5) を測る (`nozzle_report` の `wave_pct` = 平滑曲線からの最大差とは別の量)。**両量の時系列 (全場の出力ごと) を `check_quasisteady --series-csv` に渡し STEADY を要る** (停止規則の系列は波高を含まないので、その STEADY では一様性の定常を保証しない; codex ③)。
  - V5d 流れの判定: 停止規則・`check_convergence --segment`・停止規則の系列の `check_quasisteady`・壁解像 (y₁⁺ > 1 の面積 ≤ 2 %)・非有限 — #3e (6) と同じ。届かなければ閾値を動かさずユーザに聞く。
- **V6 本段の陰解法の比較: point (`run_0007_ns`、AWS A) と B0 = 壁法線のライン陰解法 (`run_0008_ns_b0`、AWS B)** (事前登録 2026-10-15、ユーザ「line implicitを使う方法でどうなるか、も同時に確かめてほしい」): 壁・格子・初期場 (run_0005_ns の最後の場から補間)・Euler の参照 (run_0006_euler res_12000、sha256 cbd236ec…)・バイナリ (FP64、sha256 3a9e347b…)・精度 (float)・前の段 (S1・S2)・停止規則 (窓 2 万・保持 4 万・上限 20 万)・出力 2,000 は同じ。違いは本段の陰解法だけ: point は cfl 1.0、B0 は case/45 の生産と同じ (lineImplicit 1・lineDtDirectional 1・lineDtDirectionalCap 50・implicitThermalJacobian 5・cfl 4・implicitRelax 0.7・リミッタの基準値は run_0007_ns の本段の自動の値に固定)。
  - 測るもの: (1) 停止規則の状態と到達の step (または上限)、本段の壁時計の時間と ms/step (A・B は同じ g5.xlarge で各 1 本だけ、別の機械なので参考)。(2) 各腕の V5a〜V5d。(3) 各腕の止まった場 (到達なら到達の出力、上限なら 200,000) の差: 出口コア M、θ_r の 3 断面、V5c の 2 量 (r/r_w = 0.1)、δ_E (スロート・出口、run の外の診断の抽出)。
  - 「同じ解」とみなす許容: |Δ出口コア M| ≤ 0.01 %、|Δθ_r| ≤ 0.1 % (core-grid の反復の不確かさ ±0.1 %)、V5c の 2 量の差 ≤ 0.02 % (絶対)、|Δδ_E| ≤ 0.05 mm (設計を終える条件と同じ)。すべて以内なら「同じ解」として速さ (step・時間) で比べる。超えたら差を報告し、どちらが正しいかは言わない (真値の参照が無い)。
  - 既知の注意: case/45 では B0 が境界層の縁の乱れを増幅して下流へ運んだ記録がある (上限 1 で消えた)。
- **V7 3 係数の壁の判別 (事前登録 2026-10-15、codex diagnose m55-three-coef)**: 対照は 2 係数の B0 (`run_0008_ns_b0`、回し直さない)、介入は壁の補正の規則だけ (3 係数)。B0・精度・リミッタの基準値・格子の規則は固定、新しい理想の形の Euler を作る、新しい NS は run_0008 の最後の場から補間。停止規則 (窓 2 万・保持 4 万・上限 20 万)。新しい NS の max|δ_in − δ_E| (x ≥ 0) が旧 (2.14 mm) の半分以下 **かつ** V5c の 2 量 (0.2247 %・0.0788 %) がともに半分以下 → 結果 A (途中の過大補正が主因の仮説を支持)。δ の誤差が半分以下でも V5c が同程度に良くならない → 結果 B (主因説を退ける)。δ の誤差が減らない → 壁の形で境界層が変わる側 (第 2 仮説)。未収束・未定常なら判別不能。製品の合格は別に V5a〜V5d。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-15 | [2026-10-15-tooling-nozzle-card-local-mesh-plan.md](../../notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan.md) | NO-GO (分類 incomplete: 見積もりの台本・出力と M6 の比の元の記録がこの作業ツリーに無い), C0/M7/m1 | **全件採用、進め方を 2 段に分ける**: M5.5 は先にケースだけの手組みの格子 (M6 の H と同じやり方; #3a) で進め、生成器の既定の変更 (#2) は指摘を直して再レビューの後。M1 (実寸の下包絡を半径で割ると補間で要求を超える: 反例 1.0607 倍) → f_req = c·y₁,req/r_w の下包絡を表にし、補間の区間の中も検査、物理壁でも実寸を再検査。M2 (`wall_normal_layer` で成果物 `mesh_structure.json`・STEP の再構築が作られない) → 既定化の前に成果物の対応 (#5)。M3 (補正の係数の根拠: H の表は冷却倍率を掛けた見込み、縮流部の比 0.61〜2.2 に係数 1.0 は保守側でない) → 「初期格子の経験の係数」と位置付け、元の run・step・比の分布を固定して commit、縮流部も比の最小 (0.6) に、補間後の c(x) が比を超えない検査。M4 (V2 の `--bcond` がこの枝の検査器に無い・AR 緩和の条件が未定量) → `runner.mesh_ar_relax_check` (AR > 1000 の全セルで skew ≤ 0.1) と変換後の第一層厚の照合・`msh_digits 17` を V2 に。M5 (V4 を記録だけで閉じられる) → 最終の完了は再調整の後の V4 PASS、全 no-slip 壁の同じ step の出力の存在、`check_convergence --segment`・壁量の系列の `check_quasisteady` を事前登録。M6 (見積もりと実生成の設定の契約・旧い overrides) → 見積もり・節点数・品質・実生成で同じパラメータの組み立て、雛形の旧キーの除去、明示の overrides の衝突は拒否、`mesh.ar_max` に従う。M7 (Euler の ni も変わる) → Euler の ni が NS と結合していることを明記し、M5 で出口 M・流量・δ_r の基準格子との差と準定常を受入に。m8 (参照の解決) → 再現の入力を case/71 に commit、関連 plan は統合の後に参照 |
| diagnose | 2026-10-15 | [2026-10-15-m55-ns-unsettled-4-diagnose.md](../../notes/reviews/2026-10-15-m55-ns-unsettled-4-diagnose.md) (ブリーフ [briefs/2026-10-15-m55-ns-unsettled-4.md](../../notes/reviews/briefs/2026-10-15-m55-ns-unsettled-4.md)) | 結論: 長い延長と δ* の補正を保留し、壁の熱流の A/B を先に。Major 4・Minor 1 | (M1) すぐ約 6 万 step 延長する → **採用** (延長しない; ユーザ決定「aかな」)。(M2) θ_r・出口 M の STEADY を理由に δ* へ進む → **採用** (STEADY は根拠にしない): 2 窓の判定を完了させた (FAIL、measured; #3d) うえで、ユーザ決定で関門の明示の上書き 3 つで進む (#3e、基準は緩めない)。(M3) 収支の残り 1.9 % を非保存と解釈する → **採用** (非保存の証拠と書かない); A/B は #3f に入れ、ユーザ決定で後回し。(M4) double の延長で thermoFloat の寄与を除外する → **採用** (「double でも全圧の超過は残る」だけを書き、旧経路の寄与は除外しない)。(m) 終了状態が未確定 → 03:15:22 UTC に CENSORED で確定 (対応済み)。報告の言い方 (予算の上限で終了・停止条件に未到達) は採用 |
| diagnose | 2026-10-15 | [2026-10-15-m55-p1-axis-wave-diagnose.md](../../notes/reviews/2026-10-15-m55-p1-axis-wave-diagnose.md) (ブリーフ [briefs/2026-10-15-m55-p1-axis-wave.md](../../notes/reviews/briefs/2026-10-15-m55-p1-axis-wave.md)) | 結論: pass 2 は保留し、pass 0 の保存場・Euler 参照・帯を固定して抽出の径方向細分 `refine` 4/32 だけを変える forge 0 step の A/B を先に。Major 8・Minor 1 | (H1 要再検証) → **採用**: #3k (refine A/B)。(H2 却下: 生産の抽出は `method="band"` で core_rms は診断値) → **採用** (core_rms で重み 0 にする案 (e) は採らない)。(pass 2 の補正方法: (a)〜(c) 単独は推奨しない、カード経路は delta_target でなく delta_E を読むので ω は効かない) → **採用** (pass 2 は #3k の後、緩和を入れるなら寸法の逆算と壁の生成が同じ表を使う設計が要る)。(成果物: pass 1 は途中の成果、合格品ではない; 試験部の一様性を優先するなら暫定の候補は pass 0 = 軸の最大偏差 0.418 % 対 0.764 %) → **採用** (ユーザに両方を示して選んでもらう)。(E(i) ṁ の 2 つの呼び手) → **採用**: #3i。(E(ii) 直管終端の単一点の例外) → **却下を採用**: #3j を STEP の曲線から接合を決めて区間全体を照合する方式に変更。(「過渡ではない」は未立証) → **採用**: #3l (軸の波の系列の `check_quasisteady`; ただし A の途中の全場は判定の後に消したので残る 5 枚 100k・140k・160k・180k・200k で行う)。(H3 の因果は未分離) → **採用** (出口コア M の改善を全体の補正量だけの効果と書かない)。(m: 「±0.6 %」は不正確) → **採用** (最大絶対偏差 0.764 %・半振幅 0.656 % と分けて書く)。(試験部の合格条件の提案: x 4.0〜7.95 m の最大絶対偏差 ≤ 0.1 %・半振幅 ≤ 0.05 %) → **保留** (ユーザ指定値ではないので、事前登録の前にユーザに聞く) |
| diagnose | 2026-10-15 | [2026-10-15-m55-two-coef-contur-diagnose.md](../../notes/reviews/2026-10-15-m55-two-coef-contur-diagnose.md) (ブリーフ [briefs/2026-10-15-m55-two-coef-contur.md](../../notes/reviews/briefs/2026-10-15-m55-two-coef-contur.md)) | 結論: 2 係数方式は採用し、実装・新規 CFD の前に #3k の固定帯 refine 4/32 の A/B を全分布で行い、平滑化後の端点から得る係数と 2 係数の壁の感度を確かめる。Major 6 | (① 方式・pass 0 の較正) → **採用** (c 0.852073279755・δ₀ −0.322005099 mm を codex が再現)。(① 「波を作らない」) → **採用** (§4.3 を「拡大部の傾きを変えない、上流 poly は δ₀ で変わる」に直した、寸法は大域最小のまま)。(① 実装: `model: contur` を overrides に、実効 scale で変換、係数前の δ_C を基底に) → **採用** (§4.3・#3n)。(② 途中の厚すぎは NS で判定) → **採用**。(③ V5b の範囲 5.49725〜5.50275、V5c は窓で 2 量を直接測り両量の時系列を準定常の判定に入れる) → **採用** (§6 V5c)。(③ V5a はゲート全体、停止規則の系列の STEADY は一様性の定常でない) → **採用** (§6 V5a・V5c)。(④ Euler は回し直す) → **採用**。(⑤ #3k の取りやめは却下、全分布で端点・2 係数の壁・壁角を見る) → **採用** (#3k を拡張)。(⑥ 初期場は run_0005_ns) → **採用** |
| diagnose | 2026-10-15 | [2026-10-15-m55-three-coef-diagnose.md](../../notes/reviews/2026-10-15-m55-three-coef-diagnose.md) (ブリーフ [briefs/2026-10-15-m55-three-coef.md](../../notes/reviews/briefs/2026-10-15-m55-three-coef.md)) | 結論: Δx 重み・端点 2 拘束の 3 係数の壁を採用し、B0 の設定とリミッタの固定値を継承した新しい NS を 1 本、2 係数 B0 を対照に。Major 4・Minor 2 | (M1 Δx 重み) → **採用** (係数を c₀ 0.677948・c₁ 0.0210363・δ₀ −0.305721 mm に; 自分で再計算して一致)。(M2 較正し直しは V5b・V5c の両方・形の残差が減らせるときだけ) → **採用** (§4.3)。(M3 リミッタは (b) 初回の本段の開始時に一度決めて以後継承、今回は今の固定値を継承、forge の平均は境界条件・壁温拘束の後の状態と `CELLS/volume`) → **採用** (#3r の設計に入れる)。(M4 B0 の既定化は NS 本段だけ、`cfl_main` を消すと Euler が既定 cfl 4 になるので Euler の設定は別に持つ) → **採用** (#3r)。(m1 x ≥ 3 m の誤差は 0.2422 mm) → **採用** (訂正)。(m2 比の非一定だけでは c₁ の必要は示せない) → **採用** (§4.3 の根拠を直した)。(判別 A/B と結果 A/B の条件) → **採用** (§6 V7)。(B0 採用と V6 合格は別、δ_E 比較は未実施) → **採用** (V6 は「θ_r 外・δ_E 未測定」のまま) |

## 7. 影響範囲

- `design/forge_design/intake/wt_card.py`・雛形 2 本・試験。既存のカード (`mesh.method` を書いていないもの) は新しい既定で格子が変わる (ユーザの方針「標準と決めた方式はコードの既定にする」)。
- case/71 の M5.5 の生成。

## 8. 完了条件

- [ ] §6 V1〜V3 を満たし、V4 を記録
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録
- [ ] `status: done`、§9 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-15` — 2 係数の壁の 2 本 (point・B0) を判定 (V5b・V5c FAIL、なめらかな過大)。ユーザ決定で B0 採用 (#3r) と 3 係数 (#3s・#3t、§4.3 の拡張)。
- `2026-10-15` — ユーザ指示で本段の陰解法の比較 (point / B0、§6 V6・#3p) を事前登録し、B0 の腕 `run_0008_ns_b0` を AWS B で起動。
- `2026-10-15` — #3k (結果 B・2 係数はすべて以内)、#3n (2 係数の実装)、#3j (寸法の関門の接合の照合) 済み。§4.3 の壁で Euler → NS (#3o) を AWS A で起動。
- `2026-10-15` — ユーザ決定で M5.5 の δ を CONTUR × 倍率 c + 一定のずらし δ₀ に (§4.3、#3m〜#3o)。pass 0 から c 0.852073・δ₀ −0.32201 mm。§6 V5 を事前登録 (許容はユーザ確認待ち)。codex diagnose (軸の波) の採否を §6.1 に。
- `2026-10-15` — δ* の補正の 1 回目 (#3e) の NS `run_0005_ns` を判定: CENSORED・NOT CONVERGED・δ の 2 窓 FAIL → 関門は止めた。出口コア M 5.50145 (pass 0 は 5.52042)、試験部の軸に M の波 −0.76〜+0.55 % (新しく出た)。δ_E − δ_in 最大 0.52 mm (pass 0 は 8.1 mm)。成果物一式を手元で作成 (寸法の関門は整合の 1 項目で FAIL、#3j)。ṁ の帳簿の誤りを #3i に。軸の波の解釈と次の pass を codex に諮問中。
- `2026-10-15` — 成果物の構造ファイル v2 (`wall_normal_blend_v1`) と STEP からの作り直しを壁法線の近壁層の格子に対応 (#5)。規則のパラメータで節点を置き直す照合で傾きの誤りを両側で検出。v1 は 6d4334e3 とバイト一致。
- `2026-10-15` — 関門の次の壁・`prepare_ns` で表の端を ε_f32 以内なら設計の端に揃える (#3h; 1 ulp 短い表で 1 本の B-spline の壁が止まった)。
- `2026-10-15` — run_0003_ns_cont (全部 double の延長) を判定: CENSORED・NOT CONVERGED・δ の 2 窓 FAIL (measured、δ の 0.3〜0.5 %)。ユーザ決定「a」で関門の明示の上書きで δ* の補正 (#3e) へ。codex 4 回目の採否を §6.1 に。`check_wall_resolution.py` の壁でない境界のダンプの扱いを修正 (#3g)。
- `2026-10-15` — plan 段レビュー NO-GO (C0/M7/m1) を全件採用。M5.5 は手組みの格子 (#3a) で先行、生成器の既定化は書き直して再レビュー。
- `2026-10-15` — 起票 (ユーザ「古いやり方でメッシュ作ってない？」、決定「a」)。見積もり (3 通りの作り方の ni・節点) と見積もりの式の偏り (M6 の H の表との比) を §3 に。
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
| `run_0006_euler` | 2 係数の壁 (plan §4.3、c 0.852088・δ₀ −0.32204 mm、`card_m55_2coef.yaml` → `problem_m55_p2_euler.yaml`、r_t 0.100243850935207 m) の Euler の参照 (AWS A)、格子 466,570 点、本段 cfl 1・12,000 step、FP64 | 完走 (軸の出口 M 5.49977、dM 最大 0.058 % Md、非有限なし)。参照 res_12000 (sha256 cbd236ec…)、他の全場は削除 | ref (run_0007・0008 の Euler の参照) |
| `run_0007_ns` | 2 係数の壁の NS (AWS A、plan §5.1 #3o): `problem_m55_p2.yaml` (格子 442,520 点、品質 PASS)、初期場 run_0005_ns の最後の場から補間、`--stages full`、float、本段 point cfl 1.0 (雛形の既定)、停止規則 (窓 2 万・保持 4 万・上限 20 万)、台本 `_runs/run_0007_ns/p2_*.sh`・`p2_ns.py` | 上限 200,000 で停止規則に届かず (D 7.2〜8.6e−4)。本段の残差 PASS、壁解像 PASS、非有限 0。**V5b FAIL (出口コア M 5.503168)・V5c FAIL (最大偏差 0.224 %・片振幅 0.079 %、波は消えなめらかな過大)**。後処理 (δ の時系列) はユーザ了承で止めた。判定の原文 `_runs/run_0007_ns/judge/` | ref (point の腕) |
| `run_0008_ns_b0` | run_0007_ns と同じ条件で本段だけ B0 (ライン陰解法 dir・cap 50・itj 5・cfl 4・緩和 0.7、リミッタの基準値は run_0007 の本段の自動の値に固定; `problem_m55_p2_b0.yaml`)、**AWS B** (A のバイナリの写し sha256 3a9e347b…)、plan §5.1 #3p・§6 V6 | **本段 90,000 で停止規則に到達** (5 万で条件、保持 4 万)。残差は約 3 万 step・12 分で床 (point は 14 万 step・34 分)。本段の残差 PASS、系列 ALL STEADY、壁解像 PASS。V5b FAIL (5.503217)・V5c FAIL (0.225 %・0.079 %)。point との差: 出口コア M 0.0009 %、θ_r (同じ帯) −0.12〜−0.15 %。ユーザ決定で B0 を採用。判定の原文 `_runs/run_0008_ns_b0/judge/` | ref (B0 の腕、2 係数の壁の今の最良) |
| `run_0009_euler` | 3 係数の壁 (plan §4.3 の末尾・#3t、`card_m55_3coef.yaml` c₀ 0.677948・c₁ 0.0210363 /m・δ₀ −0.305721 mm → `problem_m55_p3_euler.yaml`、r_t 0.1002434996788185 m) の Euler の参照 (**AWS B**)、格子 466,570 点 | 完走、参照 res_12000 (sha256 2b2d02ce…) | ref (run_0010 の Euler の参照) |
| `run_0010_ns_b0` | 3 係数の壁の NS (**AWS B**、plan §5.1 #3t・§6 V5・V7): `problem_m55_p3_b0.yaml` (本段 B0、リミッタの基準値は run_0008 と同じ固定)、初期場 run_0008_ns_b0 の最後の場から補間、float、停止規則 (窓 2 万・保持 4 万・上限 20 万)、台本 `_runs/run_0010_ns_b0/p4_all.sh`・`p4_ns.py` | **本段 86,000 で停止規則に到達**、系列 ALL STEADY、壁解像 PASS、非有限 0。**V5b FAIL (出口コア M 5.502841、+0.052 %)・V5c 最大偏差 0.094 % PASS・片振幅 0.076 % FAIL** (試験部の平均 −0.002 %、x 4→8 m で −0.094 → +0.058 % の傾き、波なし)。δ の誤差 最大 0.62 mm (2 係数 2.14 mm)。V7 は厳密には結果 A でない。`check_convergence` は NOT CONVERGED (前の解からの開始で桁数が出ない; 判定の変更 #3x へ)。**成果物 PASS** (`run_0010_ns_b0/deliverables/`、理想の壁と物理壁の比較を含む)、報告 `run_0010_ns_b0/report/` (pptx)。判定の原文 `_runs/run_0010_ns_b0/judge/` | ref (今の最良の壁) |

- **AWS A の途中の全場の削除 (2026-10-15、判定済み)**: `run_0003_ns_cont` の 96 枚 (res_0・140000・160000・180000・200000 は #3f の再生用に残す)、`run_0002_ns` の 21 枚 (res_0・res_183000 = run_0003 の初期場を残す)、`run_0001_euler`・`run_0004_euler` の各 11 枚 (res_0・res_12000 = Euler の参照を残す)。壁・出口の出力・残差履歴・判定の記録は残した。
- **AWS A の全場の削除 (2026-10-15 19:45 JST、ユーザ指示「関係ないスナップショットを丸ごと消しといて」)**: `run_0001_euler`・`run_0002_ns`・`run_0003_ns_cont`・`run_0004_euler`・`run_0005_ns` の残っていた全場 (res_0 を含む 17 枚) と `run_0006_euler` の res_12000 以外 (12 枚) を消した (計 29 枚・約 2.9 GB)。残したのは `run_0006_euler/res_12000.h5` (run_0007_ns の Euler の参照、sha256 cbd236ec…) と走行中の `run_0007_ns`。壁・出口の出力・履歴・判定の記録は残した。これで run_0003_ns_cont の場からの再生 (#3f の壁の熱流の A/B) と run_0005_ns の抽出のやり直し (Euler の参照 run_0004 の res_12000) はできない。run_0005_ns の res_200000 は手元と AWS B に写しがある

## 経緯

- 2026-10-14〜15: 旧の M_knot の規則 (`min_theta_max_v1`) でカードを生成 (AWS A `~/forge-m55`、約 98 分) → 最終の組の 4800 の網の長さ合わせが全点未収束で失敗 (plan `tooling-nozzle-mknot-fast` §3)。
- 2026-10-15: 新の規則 `coarse_min_theta_v2` (feature/mknot-fast) で M_knot 2.80 (144 s) → NS の格子の節点数の上限 (30 万) で停止 → ユーザ了承で ni・上限を明示。
- 2026-10-15: ユーザが 99 万点の了承を取り下げ。今のカードの格子の見積もりは第一層を 2 値で持つ古いやり方で、拡大部全体が最も厳しい場所の薄い第一層になる。場所ごとの第一層 (`wall_first_frac_table`) と流れ方向の密度 (`x_density_table`) で見積もり直してから格子を決める。新の規則での生成は M_knot の最終の組の途中で止めた (A `~/m55_card/`、結果なし)。
- 2026-10-15 (plan `tooling-nozzle-card-local-mesh` §5.1 #3a): カードに `mesh.ni: 3000` (生成のためだけ) を足して AWS A `~/forge-cardmesh` (479723ab) で生成 → `problem_m55_card.yaml` (M_knot 自動 2.80、r_t 0.09990877093936801 m、177 s)。その `mesh`・`mesh_euler` を `_mesh/hand_mesh.py build` で手組みの格子 (場所ごとの第一層の表・流れ方向の密度の表・壁法線の近壁層・17 桁) に置き換えて `problem_m55_ns.yaml` (Euler・NS の両方に使う)。CFD 前の物理壁で作った格子の検査 (forge は起動しない、A `~/m55_mesh_hand/check/`): NS `check_mesh_quality --ar-max 5000` PASS (AR 最大 4849.4)・AR > 1000 の 47,607 セルのスキュー最大 0.0015 (違反 0; 出口面は非平面 0.62 mm で記録だけ)・変換後の第一層厚 ≤ c·y₁,req (全 4,800 断面)・見積もりの y₁⁺ 最大 0.797。Euler PASS (AR 最大 6.1)。再現の入力と記録は `_mesh/` ([`_mesh/README.md`](_mesh/README.md))。
- 2026-10-15 (plan §5.1 #3b): カードの雛形を `axismach_final_v1` に替えて (探索の雛形の物理壁では NS の停止規則が使えない) AWS A `~/forge-cardmesh` (94790833) で生成 → `problem_m55_final_card.yaml`・`_euler.yaml` (M_knot 自動 2.80 で不変、r_t 0.0999084914252905 m、寸法 PASS、194 s)。`_mesh/hand_mesh.py build` で手組みの格子に置き換えて `problem_m55_final.yaml` (NS)・`problem_m55_final_euler.yaml` (Euler)。NS の径方向はユーザ決定 (「上限3%」) で `axis_cap_frac 0.03`・nj 92 (隣接比の最大 1.1974)。検査 (FP64 の変換器、forge は起動しない): NS 品質 PASS (AR 最大 4849.0)・AR 緩和 違反 0 (AR > 1000 の 36,115 セル、スキュー最大 0.0015)・第一層 実寸 / (c·y₁,req) 最大 0.9999997 (超過 0)・見積もりの y₁⁺ 最大 0.797・軸側のセル / r_w 0.030 (全断面); Euler PASS (AR 最大 6.1)。`run_0001_euler` を準備 (未起動)。記録は `_mesh/hand/final/` ([`_mesh/README.md`](_mesh/README.md))。
```

## 参考: `solver_density_cuda/tools/check_convergence.py`

```
#!/usr/bin/env python3
"""
収束判定ツール (AGENTS.md「収束確認 (必須)」の実体化)。

forge の run ディレクトリの residual_history.csv を読み、**全保存量の残差列**
(rms_ro, rms_roUx, rms_roUy, rms_roUz, rms_roe, RANS時 rms_roK/rms_roOmega, 化学種 rms_roY*,
凝縮 rms_rog_*/rms_roQ*_*) について
初期値・最終値・低下桁数・末尾トレンド (falling/flat/rising) を出し、明確な VERDICT を返す。

目的: 「rms_ro と NaN だけ見て収束と判断する」ことを防ぐ (AGENTS.md 違反の常習を防止)。
結果を「収束した」「一致した」と報告する前に必ず本ツールを通すこと。

使い方:
  python3 tools/check_convergence.py <run_dir> [<run_dir2> ...]
  python3 tools/check_convergence.py --drop 4 --tail 0.2 run_a run_b   # 判定閾値を調整
終了コード: 全 run が PASS なら 0、1つでも未収束/NaN があれば 1 (CI/スクリプトで使える)。

判定基準 (既定):
  - NaN/Inf が無い。
  - 各非ゼロ保存量残差が**系列のピークから** >= --drop 桁 (既定 3) 低下している
    (step 0 基準ではない: IC の作り方で step 0 が過渡ピークより桁違いに小さい
     成分があり [node の準1D IC は Uy≈0]、初期比だと誤って停滞判定になる)。
  - 末尾 (--tail, 既定 20%) が rising でない (flat/falling)。
"all-zero" 列 (例: 2D の rms_roUz) は判定から除外する。

収束場からの restart (交差 restart など) 用: --from-floor REF_RUN (別名 --reference-floor)
  ピークからの低下桁数は収束場から再開した run では原理的に 0 なので、代わりに **参照 run の末尾床**
  (REF_RUN の residual_history.csv の末尾 --tail 平均) を基準に、各残差列が run の全期間で床の
  --floor-factor 倍 (既定 1.5) 以内に留まり (ピーク ≤ factor×床)、末尾平均も factor×床 以内なら PASS。
  列ごとの床比 (末尾平均/床, ピーク/床) を表示する。NaN/Inf は DIVERGED、超過は NOT CONVERGED。
  例: python3 tools/check_convergence.py --from-floor run_0476 run_0478 run_0479
"""
import csv, math, sys, argparse, os

CONSERVED = ['rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe',
             'rms_roK', 'rms_roOmega',
             'rms_roGamma', 'rms_roReth']   # 遷移モデル (turbulence.transition: lm2009)


def load_series(path):
    rows = [r for r in csv.DictReader(open(path)) if r.get('phase') == 'outer_end']
    if not rows:
        # outer_end が無い構成 (純 explicit) では最終 inner / 全行を使う
        rows = list(csv.DictReader(open(path)))
    cols = {}
    for c in CONSERVED:
        if rows and c in rows[0]:
            cols[c] = [float(r[c]) for r in rows]
    # 化学種 (rms_roY*)・受動トレーサ (rms_roXi) と凝縮 (rms_rog_*, rms_roQ{0,1,2}_*) の保存量残差も検査する (存在時)。
    # 凝縮 run で NS 5 本 + SST 2 本だけ見て「収束」と判定していた穴 (codex 指摘 2026-09-10) を塞ぐ。
    if rows:
        for c in rows[0].keys():
            if c in cols or not c.startswith('rms_') or c.startswith('rms_dq_'):
                continue
            if c.startswith('rms_roY') or c.startswith('rms_rog_') or c.startswith('rms_roQ') or c == 'rms_roXi':
                cols[c] = [float(r[c]) for r in rows]
    return rows, cols


def _tail_rise(ser, tail_frac):
    """末尾 2*tail_frac 窓を 2 通りに評価し `(spike, slow)` を返す。

    - `spike`: **末端の急増**。窓末尾の数点の最大が窓中央値の 10 倍を超えるか。
      緩やかなトレンドの検定では捕まらない「最後の数 step での爆発」を独立に見る
      (実例: `case/37.pintle_nozzle/run_0009` の `rms_roUx` は末尾 5 点で 2.77e-4 → 7.58e22 と
      26 桁跳ねるのに、窓 444 点の回帰では上昇 0.13 桁・散らばり 1.58 桁に埋もれる)。
    - `slow`: **持続的な緩い上昇**。log10|値| の線形回帰の傾きが有意 (t > 3) で、窓全体の
      上昇が 10 % を超えるか。散らばりが大きいことは上昇が無い証明にならないので、
      「上昇 > 散らばり」ではなく回帰の有意性で見る (codex plan レビュー Major 1)。

    点数が足りない・分散が無い場合は **None** を返す (判定は既存の 2 窓平均比だけに委ねる。
    数 step で落ちた run の系列がこれに当たる)。"""
    w = [abs(x) for x in ser[int(len(ser) * (1 - 2 * tail_frac)):] if x > 0.0]
    m = len(w)
    if m < 8:
        return None
    y = [math.log10(v) for v in w]
    xb = (m - 1) / 2.0
    yb = sum(y) / m
    sxx = sum((i - xb) ** 2 for i in range(m))
    if sxx <= 0.0:
        return None
    slope = sum((i - xb) * (yi - yb) for i, yi in enumerate(y)) / sxx
    resid = [yi - (yb + slope * (i - xb)) for i, yi in enumerate(y)]
    var = sum(r * r for r in resid) / (m - 2) if m > 2 else 0.0
    se = math.sqrt(var / sxx) if var > 0.0 and sxx > 0.0 else 0.0
    rise = slope * (m - 1)
    tval = (rise / (m - 1)) / se if se > 0.0 else (math.inf if slope > 0.0 else 0.0)
    slow = tval > 3.0 and rise > math.log10(1.10)

    ntip = max(3, int(math.ceil(0.02 * m)))
    srt = sorted(w)
    med = srt[m // 2] if m % 2 else 0.5 * (srt[m // 2 - 1] + srt[m // 2])
    tip = max(w[-ntip:])
    spike = med > 0.0 and tip > 10.0 * med
    return spike, slow


# **必須の保存量残差列**。これが 1 つも無い / 欠けている CSV を合格にしてはいけない
# (2026-09-19 codex: `step,phase` だけの CSV や `rms_ro` だけの CSV が ok=True になっていた)。
REQUIRED_COLS = ('rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe')
# 遷移モデル (turbulence.transition: lm2009) が有効な run では、この 2 列も必須 (codex plan M5, 2026-09-22)。
# 有効かどうかは残差 CSV の隣の solverConfig.yaml から読む。
TRANSITION_COLS = ('rms_roGamma', 'rms_roReth')


def transition_active(csv_path):
    import re
    cfgp = os.path.join(os.path.dirname(os.path.abspath(csv_path)), 'solverConfig.yaml')
    try:
        txt = open(cfgp, encoding='utf-8').read()
    except OSError:
        return False
    txt = re.sub(r'#.*', '', txt)
    m = re.search(r'\btransition\s*:\s*["\']?([A-Za-z0-9_]+)', txt)
    return bool(m) and m.group(1).lower() not in ('none', '0')


def analyze(path, min_drop, tail_frac):
    rows, cols = load_series(path)
    if not rows:
        return None
    laststep = rows[-1].get('step', '?')
    report = {}
    ok = True
    # --- 入力の完全性検査 (合格の前提。欠けていたら「判定不能」= ok=False) ---
    if not cols:
        report['(入力)'] = ('残差列が 1 つも無い  <-- 判定不能', False)
        return laststep, report, False, False, False, False
    need = REQUIRED_COLS + (TRANSITION_COLS if transition_active(path) else ())
    missing = [c for c in need if c not in cols]
    if missing:
        report['(入力)'] = ('必須の保存量残差列が無い: %s  <-- 判定不能' % ', '.join(missing), False)
        ok = False
    any_nan = False
    any_stalled = False
    any_converging = False
    for c, ser in cols.items():
        nz = [v for v in ser if v != 0.0]
        if not nz:
            report[c] = ('all-zero (inactive, skip)', True)
            continue
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            report[c] = ('NaN/Inf present  <-- DIVERGED', False)
            ok = False; any_nan = True
            continue
        init, fin = ser[0], ser[-1]
        n = len(ser)
        a = ser[int(n * (1 - tail_frac)):]
        b = ser[int(n * (1 - 2 * tail_frac)):int(n * (1 - tail_frac))] or a
        ma = sum(abs(x) for x in a) / len(a)
        mb = sum(abs(x) for x in b) / len(b)
        # rising = 窓平均が 5% 増 **かつ** プラトーを実際に離脱している (系列最小の
        # 2 倍超)。後者のガードが無いと、深く収束したプラトー (例: 3.4 桁低下後の
        # 1e-7 台) のリミットサイクル呼吸 ±5% を「発散傾向」と誤判定する
        # (run 毎の非決定性で同一設定が PASS/NOT CONVERGED に割れる実害が出た
        # 2026-08-15)。本物のリバウンド発散は最小値の 2 倍を速やかに超えるので
        # 検出力は保たれる。
        smin = min(abs(x) for x in ser if x != 0.0) if any(x != 0.0 for x in ser) else 0.0
        # rising の 3 条件目: **末端の急増か、有意な緩い上昇があること**。これが無いと、
        # プラトー自体が数倍の幅で揺れている列では 2 窓平均の大小がジッタの位相だけで決まり、
        # 判定が再現しない。実証 (2026-09-19): case/46 の run_0193 と run_0195 は設定重複で
        # 同一形状・同一レシピになっており rms_roY1 の分布も同一 (後半中央値 1.37e-6 / 1.33e-6、
        # p5-p95 一致、帯 8.5 倍) だったのに、ma/mb が 0.961 と 1.053 に割れて flat / rising に
        # 分かれ、ゲートが PASS / FAIL に反転した。判定は `_tail_rise` に分離してある。
        rise = _tail_rise(ser, tail_frac)
        over_jitter = True if rise is None else (rise[0] or rise[1])
        trend = ('rising' if (ma > mb * 1.05 and ma > 2.0 * smin and over_jitter)
                 else ('flat' if ma > mb * 0.9 else 'falling'))

        # init==0 (例: アライン格子で Uy が初期厳密 0) も下のピーク基準で判定する。旧特例 (rising でなければ
        # 合格) は [0,1,1,...] のように立ち上がって落ちない列を合格にしてしまった (codex 指摘 2026-09-10)。
        # 全期間ゼロの列だけを上の all-zero で除外する。

        # 低下桁数は **step 0 ではなく系列のピーク**から測る。IC の作り方によっては
        # ある成分の step 0 残差が過渡ピークより桁違いに小さいことがあり (node の
        # 準 1D IC は半径方向速度がほぼ厳密 0 なので rms_roUy の init が 4.5e-6、
        # 一方 step 2 のピークは 4.8e-2 = 4 桁上)、初期比だと実際に 5.0 桁落ちて
        # いる収束列を「0.9 桁で停滞」と誤判定する (2026-08-17 実測)。ピークは
        # 「その成分が実際にどこから落ちたか」なので物理的にも正しい尺度。
        # step 0 がピークの通常ケースでは値は変わらない (後方互換)。
        peak = max(abs(x) for x in ser)
        # 低下桁数は **末尾窓の代表値**で測る。最終 1 点だけを使うと、[1]*19+[0] のように
        # 末尾が 1 点だけ 0 の系列で drop=inf になり合格してしまう (2026-09-19 codex)。
        # 代表値は末尾窓の |値| の中央値。窓に 0 が混じっても中央値は 0 になりにくく、
        # 本当に全体が 0 まで落ちた列は下の all-zero か「窓中央値 0」で判定不能にする。
        tail_abs = sorted(abs(x) for x in a)
        fin_rep = tail_abs[len(tail_abs) // 2]
        if fin_rep == 0.0:
            report[c] = ('末尾窓の代表値が 0 (残差が数値的にゼロ)  <-- 判定不能', False)
            ok = False
            continue
        drop = math.log10(peak / fin_rep)
        col_ok = drop >= min_drop and trend != 'rising'
        ok = ok and col_ok
        # 未達の理由を区別: falling=収束途中(あと steps)、flat=停滞、rising=発散傾向
        status = ''
        if not col_ok:
            if trend == 'rising':   status = '  <-- RISING (divergent)'
            elif trend == 'flat':   status = '  <-- STALLED (plateau)'; any_stalled = True
            else:                   status = '  <-- still converging'; any_converging = True
        pk = "" if abs(peak - init) <= 1e-30 * max(abs(init), 1.0) else f" peak={peak:.2e}"
        report[c] = (f"init={init:.2e}{pk} fin={fin:.2e} drop={drop:4.1f}dec {trend:7s}{status}", col_ok)
    return laststep, report, ok, any_nan, any_stalled, any_converging


def build_segment_csv(run_dir):
    """`stage_manifest.json` の**最後の区間**を step オフセット付きで連結した CSV を書いて返す。

    段名でなく**実効設定 (hard キー)** が同じ連続区間だけを繋ぐので、別の方程式・BC の過渡を
    本段の低下桁数の基準にしてしまう事故が起きない (2026-09-19 codex Major 1)。
    列が段で違う場合は**共通列に落とさず**、区間内で列集合が一致することを要求する
    (共通列への縮退は検出したい誤合格を再導入する: 同 Major 2)。
    """
    import csv as _csv
    import json as _json
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from stage_manifest import segments, load_launches
    p = os.path.join(run_dir, 'stage_manifest.json')
    if not os.path.exists(p):
        return None
    segs = segments(_json.load(open(p)), load_launches(os.path.dirname(p)))
    if not segs:
        return None
    seg = segs[-1]
    if any(st.get('chi_source') == 'ambiguous' for st in seg):
        # 同じ設定で chi の実効値が違う起動が混ざる段 (stage_manifest.load_launches)。どの履歴がどの起動か分からない
        # ので連結しない (判定不能。codex result 2026-09-27 chi-default M3)
        print(f"[{run_dir}] 最後の区間に chi の実効値が食い違う起動の段がある -> --segment は使えない "
              f"(判定区間を人が明示すること)  <-- 判定不能")
        return None
    rows, hdr, off = [], None, 0
    used = []
    for st in seg:
        f = os.path.join(run_dir, st['history'])
        if not os.path.exists(f):
            # 段が 1 つだけの run (継続 run の本段のみ) は段別履歴を分けないので、
            # 素の residual_history.csv がその段の履歴そのものである (2026-09-20)。
            # ここで黙って飛ばすと区間が空になり「stage_manifest が無い」と同じ扱いになっていた。
            alt = os.path.join(run_dir, 'residual_history.csv')
            if len(seg) == 1 and os.path.exists(alt):
                f = alt
            else:
                continue
        r = list(_csv.reader(open(f)))
        if len(r) < 2:
            continue
        h = [c.strip() for c in r[0]]
        if hdr is None:
            hdr = h
        elif h != hdr:
            print(f"  [segment] 段 {st['tag']} の列が区間内で一致しない -> 連結しない "
                  f"(共通列に落とすと誤合格を再導入する)")
            return None
        body = [x for x in r[1:] if x and x[0].strip().lstrip('-').isdigit()]
        if not body:
            continue
        n1 = int(body[-1][0])
        for x in body:
            x = list(x); x[0] = str(int(x[0]) + off); rows.append(x)
        used.append(st['tag'])
        off += n1 + 1
    if not rows:
        return None
    out = os.path.join(run_dir, 'residual_history_segment.csv')
    with open(out, 'w', newline='') as f:
        w = _csv.writer(f); w.writerow(hdr); w.writerows(rows)
    print(f"  [segment] 判定区間 = {' -> '.join(used)}  ({len(rows)} 行) -> {out}")
    return out


def reference_floor(path, tail_frac):
    """参照 run の各残差列の末尾床 (末尾 tail_frac の |値| 平均)。
    戻り値 (floor, zero_cols): floor は非ゼロ列だけ、zero_cols は all-zero (非活性) 列の集合。"""
    rows, cols = load_series(path)
    floor = {}
    zero_cols = set()
    for c, ser in cols.items():
        if not any(v != 0.0 for v in ser):
            zero_cols.add(c)
            continue
        a = ser[int(len(ser) * (1 - tail_frac)):]
        floor[c] = sum(abs(x) for x in a) / len(a)
    return floor, zero_cols


def analyze_from_floor(path, floor, zero_cols, factor, tail_frac, ref_path=None):
    """収束場からの restart 判定: 全期間ピークと末尾平均が参照床の factor 倍以内なら列 PASS。
    列対応 (codex result-2 m1): 参照に無い列・参照で all-zero だった列が対象で非ゼロなら FAIL (前提不成立を PASS にしない)。
    参照の非ゼロ列が対象に無い場合も FAIL。"""
    rows, cols = load_series(path)
    if not rows:
        return None
    laststep = rows[-1].get('step', '?')
    report = {}
    ok = True
    any_nan = False
    # 必須列の検査は通常判定と共通 (codex result M2: --from-floor が遷移残差の欠落を通していた)。
    # 参照と対象で方程式系 (遷移モデルの有無) が違えば、床の比較そのものが成り立たない。
    if ref_path is not None and transition_active(ref_path) != transition_active(path):
        report['(参照)'] = ('参照 run と遷移モデルの有無が違う (別の方程式系)  <-- 判定不能', False)
        ok = False
    need = REQUIRED_COLS + (TRANSITION_COLS if transition_active(path) else ())
    missing = [c for c in need if c not in cols]
    if missing:
        report['(入力)'] = ('必須の保存量残差列が無い: %s  <-- 判定不能' % ', '.join(missing), False)
        ok = False
    for c in floor:
        if c not in cols:
            report[c] = ('missing in target (reference has a nonzero floor)  <-- COLUMN MISMATCH', False)
            ok = False
    for c, ser in cols.items():
        if not any(v != 0.0 for v in ser):
            report[c] = ('all-zero (inactive, skip)', True)
            continue
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            report[c] = ('NaN/Inf present  <-- DIVERGED', False)
            ok = False; any_nan = True
            continue
        if c not in floor:
            if c in zero_cols:
                report[c] = ('reference column was all-zero but target is nonzero  <-- COLUMN MISMATCH', False)
            else:
                report[c] = ('no such column in reference  <-- COLUMN MISMATCH', False)
            ok = False
            continue
        peak = max(abs(x) for x in ser)
        a = ser[int(len(ser) * (1 - tail_frac)):]
        tail_mean = sum(abs(x) for x in a) / len(a)
        r_tail = tail_mean / floor[c] if floor[c] > 0 else float('inf')
        r_peak = peak / floor[c] if floor[c] > 0 else float('inf')
        col_ok = (r_tail <= factor) and (r_peak <= factor)
        ok = ok and col_ok
        status = '' if col_ok else ('  <-- ABOVE FLOOR (peak)' if r_peak > factor else '  <-- ABOVE FLOOR (tail)')
        report[c] = (f"ref_floor={floor[c]:.2e} tail={tail_mean:.2e} (x{r_tail:.2f}) peak={peak:.2e} (x{r_peak:.2f}) "
                     f"fin={ser[-1]:.2e}{status}", col_ok)
    return laststep, report, ok, any_nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dirs', nargs='+')
    ap.add_argument('--drop', type=float, default=3.0, help='required orders-of-magnitude drop')
    ap.add_argument('--tail', type=float, default=0.2, help='tail fraction for trend check')
    ap.add_argument('--from-floor', '--reference-floor', dest='from_floor', default=None, metavar='REF_RUN',
                    help='収束場からの restart 判定: REF_RUN の末尾床 (tail 平均) の --floor-factor 倍以内に全期間留まれば PASS')
    ap.add_argument('--floor-factor', type=float, default=1.5, help='--from-floor の許容倍率 (既定 1.5)')
    ap.add_argument('--segment', action='store_true',
                    help='段階起動の run で **stage_manifest.json の最後の区間** '
                         '(方程式・BC・空間離散化が同一の連続区間) を連結して判定する。'
                         '段名でなく実効設定で区間を決めるので、別 BC の過渡を基準にしない')
    args = ap.parse_args()

    floor = None
    if args.from_floor:
        ref_path = args.from_floor if args.from_floor.endswith('.csv') else os.path.join(args.from_floor, 'residual_history.csv')
        if not os.path.exists(ref_path):
            print(f"[{args.from_floor}] NO residual_history.csv (reference)"); sys.exit(1)
        # 参照 run 自身が通常判定 (--drop/--tail) を PASS していることを要求する (codex result-2 m1: 未収束 run を参照にすると
        # 「床に留まった」だけで PASS になる)。不成立なら判定を拒否 (exit 2)。
        ref_res = analyze(ref_path, args.drop, args.tail)
        if ref_res is None:
            print(f"[{args.from_floor}] empty reference residual file"); sys.exit(2)
        ref_laststep, ref_report, ref_ok, ref_nan, _, _ = ref_res
        if not ref_ok:
            print(f"=== reference {args.from_floor}  [last step {ref_laststep}]  -> "
                  f"{'DIVERGED (NaN/Inf)' if ref_nan else 'NOT CONVERGED'} ===")
            for c, (msg, _) in ref_report.items():
                print(f"  {c:12s}: {msg}")
            print("\nREFUSED: --from-floor requires a reference run that PASSES the normal criterion "
                  f"(--drop {args.drop}, --tail {args.tail}); its residual floor is not a converged floor.")
            sys.exit(2)
        floor, zero_cols = reference_floor(ref_path, args.tail)
        print(f"reference floor from {args.from_floor} (normal criterion PASS at step {ref_laststep}; tail {args.tail:.0%} mean), "
              f"factor {args.floor_factor}; columns {sorted(floor)}; inactive {sorted(zero_cols)}")

    all_pass = True
    for rd in args.run_dirs:
        if args.segment and not rd.endswith('.csv'):
            path = build_segment_csv(rd)
            if path is None:
                # **判定不能は合格ではない** (AGENTS.md)。以前は未定義の `worst` を触って
                # UnboundLocalError で落ち、呼び出し側から見ると「判定なし」で素通りしていた。
                print(f"[{rd}] stage_manifest.json が無い -> --segment は使えない "
                      f"(判定区間を人が明示すること)  <-- 判定不能")
                all_pass = False
                continue
        else:
            path = rd if rd.endswith('.csv') else os.path.join(rd, 'residual_history.csv')
        if not os.path.exists(path):
            print(f"[{rd}] NO residual_history.csv"); all_pass = False; continue
        if floor is not None:
            res = analyze_from_floor(path, floor, zero_cols, args.floor_factor, args.tail, ref_path)
            if res is None:
                print(f"[{rd}] empty residual file"); all_pass = False; continue
            laststep, report, ok, any_nan = res
            verdict = ('PASS (within %.1fx of reference floor %s)' % (args.floor_factor, args.from_floor) if ok
                       else 'DIVERGED (NaN/Inf)' if any_nan
                       else 'NOT CONVERGED (residual left the reference floor, or column mismatch)')
            print(f"\n=== {rd}  [last step {laststep}]  -> {verdict} ===")
            for c, (msg, _) in report.items():
                print(f"  {c:12s}: {msg}")
            all_pass = all_pass and ok
            continue
        res = analyze(path, args.drop, args.tail)
        if res is None:
            print(f"[{rd}] empty residual file"); all_pass = False; continue
        laststep, report, ok, any_nan, any_stalled, any_converging = res
        if ok:
            verdict = 'PASS (converged)'
        elif any_nan:
            verdict = 'DIVERGED (NaN/Inf)'
        elif any_stalled:
            verdict = 'NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)'
        elif any_converging:
            verdict = 'NOT CONVERGED (still converging — run more steps)'
        else:
            verdict = 'NOT CONVERGED'
        print(f"\n=== {rd}  [last step {laststep}]  -> {verdict} ===")
        for c, (msg, _) in report.items():
            print(f"  {c:12s}: {msg}")
        all_pass = all_pass and ok
    print(f"\nOVERALL: {'ALL PASS' if all_pass else 'CHECK FAILURES ABOVE'}")
    sys.exit(0 if all_pass else 1)


if __name__ == '__main__':
    main()
```

## 参考: `design/forge_design/evaluate/stop_rule.py`

```
r"""NS の本段の停止規則 (plan tooling-nozzle-cooled-wall-ns-settling §4.6・plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a・#36e、2026-10-14)。

**由来と統合** (2026-10-14、main への統合 1 回目で 2 つの版を 1 つにした):
forge-f5 の feature/nozzle-cooled-wall-settling (66364fff に θ_r の断面の基準を設計壁の出口に直した 60de8125 を重ねた版、plan
tooling-nozzle-cooled-wall-ns-settling §4.6、ユーザ決定・諮問 stop-rule を全件採用) を、M6 の feature/nozzle-wall-fit-and-pipeline が
写して変えた版 (plan tooling-nozzle-isothermal-wall-chain §5.1 #36a・#36e) を土台にし、使うかどうかの既定を forge-f5 の版に戻した。
  - **既定で使う** (forge-f5、ユーザ決定): 問題 YAML に `evaluate.stop_rule` が無ければ使う。使えない (euler_ref が無い・見つからない・
    出力の間隔・上限が出力の間隔の倍数でない・設計壁の出口が無い) ときは forge-f5 の版の動き = UNAVAILABLE を記録し、本段を `fallback_steps`
    (`--steps` > `evaluate.ns_steps` > 60,000) の長さで回す。
  - **`evaluate.stop_rule: true` を明示** (または関数の引数 stop_rule=True) したときだけ、使えない設定を `prepare_ns` が run dir を作る前に・
    `run_staged_ns` が forge を起動する前に例外にする (M6 の版の動き; 記録 `explicit`)。`evaluate.stop_rule: false`・CLI `--no-stop-rule` は使わない。
  - 到達の後の継続 `evaluate.stop_rule_hold_steps` (M6 の版で足した) の**既定は `HOLD_STEPS_DEFAULT` = 40,000 step** (ユーザ決定、forge-f5 経由;
    値はこの 1 か所の定数)。最初に条件を満たした出力 s₀ の後も回し、s₀ + H までの全ての出力で条件を保ったら REACHED (止める出力は s₀ + H
    以後の最初の出力)。途中で外れたら到達を取り消して `hold.breaks` に記録し、次の到達から数え直す。上限 (`cap`) は到達 + 保持を含む本段の
    step 数で、保持の途中で上限に届いたら CENSORED。0 を書けば最初の到達で止める (forge-f5 の版の動き)。
  - M6 の版で足したものは残す: 終了の状態の記録 `STOP_RULE_STATE.json` (`write_state`・`read_state`; `stop_rule.json` は見張りの途中も
    書き換わる詳しい記録、こちらは本段が終わった後に 1 度だけ書く状態の印)、問題 YAML の `evaluate.stop_rule_thresholds` (`thresholds` と
    同じキー; 例 case/45 の試験部の断面 x = 40・70・94 r_t は `x_fracs` に x / x_w を書く。既定は forge-f5 の版のまま)、上限の倍数の検査、
    (§5.1 #36e) 帯の位置の抽出を受け渡しの抽出と同じ条件 (`deltastar_loop.handoff_extract_settings`: 一定 x の断面への補間
    `remap_constant_x=True`・端の補間範囲の不足は欠測 `remap_edge="missing"`) で行うこと、Euler の参照の固定 {run, res, sha256}
    (`deltastar.resolve_euler_ref`、sha が合わなければ使わない)、(§5.1 #36e (7)) 見張りの中で出力ごとに δ_E・δ_target (抽出器の量、使用域
    [0, x_w] の全列) を受け渡しと同じ条件で抽出して `stop_rule.json` の行の `delta` に残すこと (**判定には使わない記録の量**。受け渡しの
    関門の 2 窓の判定がこの記録を読み、計算の後に抽出し直さない)。抽出は別スレッド 1 本で行い (`DeltaRecorder`)、停止の判定を待たせない。
    追いつかないときは最新の出力から 2 窓 (2W) より古い出力を飛ばす (`skipped`; 停止の時刻 s ≥ 最新の出力なので、[s − 2W, s] に入らない
    出力だけ)。要求の断面の位置 (x_fracs × x_w) は `sections` に残す。
  - θ_r の断面の基準は設計壁の出口 `physical_wall.domain[1]` (60de8125 の修正; 軸 Mach 則の終点 x_E ではない)。判定・量・見張り・forge の
    止め方は forge-f5 の版のまま。
  - (§5.1 #36g) 帯の生成と δ の記録は抽出器 (remap_edge="missing") を通るので、Euler の参照の補間できる範囲は共通の
    `deltastar.euler_cover_range` (端の外側だけに幾何の許容 = 端の列の間隔 × 1e−4、許容の内側で端を越える列は端の列の値) になる。
    2 窓に入るべき出力は `required_window_steps` (出力の間隔から列挙; 記録にある行だけから作らない) で決める。

**新しい規則であり、旧 `m9_watch.py` (§6.11、`case/45.isobutane_m6_d155/m9_watch.py`) と同等とは言わない**: m9_watch は窓の両端の差で傾向と振動を
混ぜて判定し、断面・量・連続の到達の条件もここと違う (codex diagnose 2026-10-14、`notes/reviews/2026-10-14-half-inlet-370k-design-diagnose.md`)。
ここは M6 の考え方を移したうえで、**平均の変化と振幅を分けて**判定する別の規則として登録する。m9_watch の到達の step・結果をこの規則の結果として読み替えない。
REACHED は「この規則の量が停止の条件を満たした」の意味に限り、収束・設計の量の定常・壁の更新に進む条件 (§4.9) の合格ではない。

判定量 (出力 res_<N>.h5 ごと、N > 0):
  - 質量の収支 $r = (\dot m_{in} - \dot m_{out})/\dot m_{in}$: 入口 (列 i = 0)・出口 (列 i = ni − 1) の列の折れ線を通る流束
    $2\pi\int(\rho u_x\,dr - \rho u_r\,dx)\,r$ (M6 の `cold_xcheck.reduce_fields` と同じ式)。
  - 境界層の運動量厚さ $\theta_r$ を $x/x_w$ = 0.4・0.7・0.95 の断面で ($x_w$ は設計壁 (物理壁) の出口 `prepare_info.physical_wall.domain[1]`、r_t 単位;
    `prepare_info.x_E` は軸 Mach 則の終点で別物 — 2026-10-14 まで誤って x_E を使っていた、M6ノズルの設計・forge-25 の指摘)。各行 j を x 方向に補間した
    一定 x の断面で、帯の外縁 $y_b$ の $(\rho_e, u_e)$ を縁にした面積で等価な厚さ
    $\theta_r = r_w - \sqrt{r_w^2 - 2\int_{r_b}^{r_w}\rho u(u_e-u)\,r\,dr/(\rho_e u_e^2)}$ (`reduce_fields` と同じ式、断面は最も近い列)。
  - $y_b(x)$ は**本段の最初の出力**を固定 Euler 参照 (`prepare_info.euler_ref`) で `extract_and_merge(.., band_select="edge")` と同じ抽出にかけた
    帯の位置 (縁 y* の 1.25 倍、さらに掛けない; 諮問 stop-rule (a)) を本段の間は固定する (`fix_band`)。
  - 記録だけ (判定に使わない; 諮問 (d)): 出口コア M (最後の列の η = r/r_w ∈ [0.05, 0.7] の M の平均) と無次元の流量
    (2π∫ρU_x r dr の x ∈ (−2.5, x_E) の中央値) — `report/nozzle_report.py` の `metrics` と `case/70.m5_air_wt_demo/qs_series.py` と同じ定義。
    非有限の数 (VALUE の全ての浮動小数点の量)。

判定 (`judge_rows`、閾値・窓・断面は引数、既定は §4.6):
  - 窓 = 直近 W = 2 万 step の出力 (step ∈ [s − W, s]、出力の間隔に依らない)。本段の履歴が実際に W を覆い (最初の出力 ≤ s − W)、窓に 9 点以上。
  - r: |窓の平均| ≤ 0.05 %・振幅 (最大 − 最小)/2 ≤ 0.5 %。
  - θ_r (各断面): $D_\theta$ = |後半の平均 − 前半の平均| / |窓の平均| × W / (後半の平均の step − 前半の平均の step) ≤ 0.05 % (2 万 step あたりに換算;
    諮問の Major)、振幅 (最大 − 最小)/2 / |窓の平均| ≤ 0.5 %。前半・後半は窓の中点 s − W/2 で分け、**中点ちょうどの出力は後半に入れる**
    (点数が奇数で等間隔のとき。諮問の計算と同じ割り方。換算するので線形の系列では割り方に依らない)。換算しない値も `D_raw` として記録する。
  - 満たした最初の出力で REACHED (連続は要らない; 諮問 (b))。ただし hold_steps (既定 `HOLD_STEPS_DEFAULT`) > 0 のときは上の「到達の後の継続」。
    非有限で DIVERGED。REACHED は「指定の量の停止の条件を満たした」の意味に限り、収束・設計の量の定常の VERDICT とは別。

状態: REACHED / CENSORED (上限まで正常に回って届かない) / DIVERGED (非有限、または forge の detectNaN) / EXEC_ERROR (forge が上限より前に異常終了・
見張りの例外) / DATA_ERROR (出力が読めない・判定量が有限でない) / UNAVAILABLE (停止規則を使えない: euler_ref が無い・抽出できない・出力の間隔が
窓に 9 点を入れられない — 本段は ns_steps の長さで回す)。記録は `<run>/stop_rule.json` (一時ファイル → rename で原子的に書く)。
明示 (`evaluate.stop_rule: true`) の run は準備・起動の時点で使えないと分かったら例外で止めるので、UNAVAILABLE は見張りの途中で帯の位置を
決められないとき (本段は prepare の nStepOuter = fallback_steps の出力で止める) だけ出る。既定 (キー無し) の run は準備・起動の時点で使えなくても
UNAVAILABLE を記録して fallback_steps の長さで回す。本段が終わったら `run_staged_ns` が終了の状態を `<run>/STOP_RULE_STATE.json` に
1 度だけ書く (`write_state`)。`run_staged_ns` の終了コード (REACHED・CENSORED・一部の UNAVAILABLE で 0) で合否を決めず、この状態を直接見る
(`feedback.deltastar_loop.handoff_gate`)。

見張り (`watch_main_stage`): `run_case.sh` を子プロセス (新しいプロセスグループ) で持ったまま、書き終わった res を判定し、REACHED・DIVERGED で forge に
SIGTERM → 待って SIGKILL (止まったことを /proc で確かめる)、止めた出力より後の res を消す。forge の終了コードは `run_case.sh` の
`[run_case] forge exit=N` の行から読む (`FORGE_STRICT_CONVERGENCE` の終了コードと取り違えない)。forge だけを止めるので `run_case.sh` の後処理
(収束判定・残差の png) はそのまま走る; 走らなかったら同じ後処理を呼ぶ。
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

# §4.6 の既定 (閾値は M6 の実績を移しただけ; 風洞の形状での妥当性は §6 S-stop (c) で見る)
WINDOW = 20000
MIN_POINTS = 9
R_MEAN_MAX = 5e-4          # |r の窓の平均| ≤ 0.05 %
R_AMP_MAX = 5e-3           # r の振幅 ≤ 0.5 %
THETA_DRIFT_MAX = 5e-4     # D_θ ≤ 0.05 % / 2 万 step
THETA_AMP_MAX = 5e-3       # θ_r の振幅 / |平均| ≤ 0.5 %
X_FRACS = (0.4, 0.7, 0.95)
CAP_DEFAULT = 200_000      # 停止規則を使うときの本段の上限 (諮問 (c)、M6 と同じ)。到達 + 保持を含む本段の step 数
HOLD_STEPS_DEFAULT = 40_000   # 到達の後の継続の既定 (2026-10-14 ユーザ決定、forge-f5 経由)。evaluate.stop_rule_hold_steps で上書き (0 = 最初の到達で止める)
STATES = ("REACHED", "CENSORED", "DIVERGED", "EXEC_ERROR", "DATA_ERROR", "UNAVAILABLE")
STATE = "STOP_RULE_STATE.json"      # 本段の終了の状態 (run_staged_ns が本段の後に 1 度だけ書く; M6 の版で足した)
PLAN = ("plans/active/tooling-nozzle-cooled-wall-ns-settling.md §4.6・plans/active/tooling-nozzle-isothermal-wall-chain.md §4.9・§5.1 #36a "
        "(2026-10-14 統合: 既定で使う・明示のときだけ使えない設定を例外・保持の既定 40,000)")
RECORD = "stop_rule.json"
BAND_DIR = "stop_rule_band"
_RES_RE = re.compile(r"res_(\d+)\.h5")
_OUT_RE = re.compile(r"res_(?:(?!nan_).+_)?(\d+)\.(?:h5|xmf)")   # res_<N>.h5/.xmf と境界の res_<名前>_<id>_<N>.h5/.xmf (res_nan_* は除く)
_FORGE_EXIT_RE = re.compile(r"\[run_case\] forge exit=(-?\d+)")


class Unavailable(Exception):
    """停止規則を使えない (例外として外へ出さず、UNAVAILABLE として記録する)。"""


def thresholds(**kw) -> dict:
    """判定の閾値・窓・断面 (既定は §4.6)。未知のキーは拒否する。値の検査 (M6 の版で足した、問題 YAML の `evaluate.stop_rule_thresholds` も
    ここを通る): window・min_points は正の整数 (min_points ≥ 2)、4 つの閾値は有限の正の数、x_fracs は (0, 1] の有限の数の空でない
    重複の無いリスト (設計壁の出口 x_w に対する比; 断面は x_fracs × x_w に最も近い列)。"""
    d = dict(window=WINDOW, min_points=MIN_POINTS, r_mean_max=R_MEAN_MAX, r_amp_max=R_AMP_MAX,
             theta_drift_max=THETA_DRIFT_MAX, theta_amp_max=THETA_AMP_MAX, x_fracs=list(X_FRACS))
    bad = set(kw) - set(d)
    if bad:
        raise ValueError(f"停止規則の未知の設定 {sorted(bad)}")
    d.update(kw)
    for k in ("window", "min_points"):
        v = d[k]
        if isinstance(v, bool) or not isinstance(v, int) or v <= 0 or (k == "min_points" and v < 2):
            raise ValueError(f"停止規則の {k} は正の整数{' (≥ 2)' if k == 'min_points' else ''} ({v!r})")
    for k in ("r_mean_max", "r_amp_max", "theta_drift_max", "theta_amp_max"):
        v = d[k]
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) or v <= 0:
            raise ValueError(f"停止規則の {k} は有限の正の数 ({v!r})")
    xf = d["x_fracs"]
    if (not isinstance(xf, (list, tuple)) or not xf
            or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) or not 0.0 < v <= 1.0 for v in xf)
            or len({float(v) for v in xf}) != len(xf) or len({theta_key(v) for v in xf}) != len(xf)):
        raise ValueError(f"停止規則の x_fracs は (0, 1] の数の空でない重複の無いリスト ({xf!r})")
    d["x_fracs"] = [float(v) for v in xf]
    return d


def theta_key(frac: float) -> str:
    return f"theta_r_{float(frac):g}"


def required_window_steps(step: int, window: int, out_int: int, n_windows: int = 2) -> list:
    """停止の出力 s で終わる n_windows 窓 [s − n W, s] (閉区間) に入るべき出力の step を、出力の間隔から列挙する (昇順)。
    plan tooling-nozzle-isothermal-wall-chain §5.1 #36g (codex diagnose 2026-10-14 の Major「判定窓の完全性」): 記録 (stop_rule.json の行) に
    ある出力だけから窓を作ると、行そのものの欠けを検出できない。forge は outStepInterval の倍数で res を書き、s もその倍数なので、
    必要な出力は s − k·Δ (k = 0..⌊n W/Δ⌋)。V-h370 (W 20,000・Δ 2,500・2 窓) では 17 出力。引数が整数でなければ ValueError。"""
    vals = {"step": step, "window": window, "out_int": out_int, "n_windows": n_windows}
    for k, v in vals.items():
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)) or (v < 0 if k == "step" else v <= 0):
            raise ValueError(f"required_window_steps: {k} は {'0 以上' if k == 'step' else '正'}の整数 ({v!r})")
    s, n, d = int(step), int(n_windows) * int(window), int(out_int)
    return sorted(s - k * d for k in range(n // d + 1))


# ------------------------------------------------------------------ 設定 (prepare_ns が prepare_info.json の stop_rule に書く)
EVALUATE_KEYS = ("stop_rule", "stop_rule_cap", "stop_rule_thresholds", "stop_rule_hold_steps")


def check_evaluate_keys(evaluate: dict) -> tuple:
    """問題 YAML の `evaluate.stop_rule` (true / false)・`evaluate.stop_rule_cap` (正の整数)・`evaluate.stop_rule_thresholds`
    (`thresholds` のキーの mapping)・`evaluate.stop_rule_hold_steps` (0 以上の整数) を検査して (stop_rule, cap, 閾値の上書き, hold) を返す
    (無ければ None・None・{}・`HOLD_STEPS_DEFAULT`; hold は null も既定)。不正値は ValueError (prepare_ns は run dir を作る前に呼ぶ)。"""
    req = cap = None
    if "stop_rule" in evaluate:
        req = evaluate["stop_rule"]
        if not isinstance(req, bool):
            raise ValueError(f"evaluate.stop_rule は true / false ({req!r})")
    if evaluate.get("stop_rule_cap") is not None:
        cap = evaluate["stop_rule_cap"]
        if isinstance(cap, bool) or not isinstance(cap, int) or cap <= 0:
            raise ValueError(f"evaluate.stop_rule_cap は正の整数 ({cap!r})")
    over = evaluate.get("stop_rule_thresholds") or {}
    if not isinstance(over, dict):
        raise ValueError(f"evaluate.stop_rule_thresholds は mapping ({over!r})")
    thresholds(**over)                          # 未知のキー・不正値はここで ValueError
    hold = evaluate.get("stop_rule_hold_steps")
    if hold is None:
        hold = HOLD_STEPS_DEFAULT
    if isinstance(hold, bool) or not isinstance(hold, int) or hold < 0:
        raise ValueError(f"evaluate.stop_rule_hold_steps は 0 以上の整数 ({hold!r})")
    return req, cap, dict(over), int(hold)


def is_explicit(sr: dict | None) -> bool:
    """prepare_info の stop_rule が**明示の要求** (`evaluate.stop_rule: true` か引数 stop_rule=True) か。明示なら使えない設定は例外、
    既定 (キー無し) なら UNAVAILABLE を記録して fallback_steps で回す。統合 (2026-10-14) より前の記録 (`explicit` が無い) は
    requested_source で決める ("evaluate.stop_rule"・"引数" は明示、"既定 (使う)" は既定)。"""
    if not isinstance(sr, dict) or not sr.get("requested"):
        return False
    if "explicit" in sr:
        return bool(sr["explicit"])
    return sr.get("requested_source") in ("evaluate.stop_rule", "引数")


def resolve_settings(evaluate: dict, *, nsteps_arg, n_main: int, n_main_source: str, out_int: int, euler_ref,
                     requested: bool | None = None) -> dict:
    """停止規則の設定を決める。requested: None なら `evaluate.stop_rule` (**無ければ使う** = 既定)、True / False は関数の引数・
    CLI の --no-stop-rule。`explicit` は明示の要求 (引数 True か `evaluate.stop_rule: true`) — 明示なら呼ぶ側 (`prepare_ns`・`run_staged_ns`)
    は使えない設定を例外にし、既定なら `use` = False と理由を記録して本段を fallback_steps で回す。
    上限: 引数 nsteps_arg (`--steps`) があればそれ、無ければ `evaluate.stop_rule_cap` (既定 200,000)。上限は到達 + 保持を含む本段の step 数。
    使わない・使えないときの本段の長さは n_main (prepare_ns の `--steps` > `evaluate.ns_steps` > 60,000) — `fallback_steps` に記録する。
    使えない条件 (`use` = False と理由): euler_ref が無い、出力の間隔が窓に MIN_POINTS 点を入れられない (間隔 > W/(MIN_POINTS − 1))、
    上限が出力の間隔の倍数でない (forge は outStepInterval の倍数でしか res を書かないので、上限の出力が無く CENSORED を見分けられない)。"""
    from ..feedback.deltastar_loop import handoff_extract_settings
    from ..metrics.deltastar import resolve_euler_ref
    yaml_req, yaml_cap, th_over, hold = check_evaluate_keys(evaluate)
    th = thresholds(**th_over)
    if requested is None:
        if yaml_req is not None:
            requested, src, explicit = yaml_req, "evaluate.stop_rule", bool(yaml_req)
        else:
            requested, src, explicit = True, "既定 (使う)", False
    else:
        requested, src = bool(requested), "引数"
        explicit = requested
    # euler_ref: run dir (str) か固定のスナップショット {run, res, sha256} (§5.1 #36e (6); ここでは形だけ検査、sha は prepare_ns・起動前に照合)
    eref = resolve_euler_ref(euler_ref, verify=False)["record"] if euler_ref else None
    rec = {"requested": requested, "requested_source": src, "explicit": explicit,
           "fallback_steps": int(n_main), "fallback_steps_source": n_main_source,
           "out_step_interval": int(out_int), "euler_ref": eref, **th,
           "thresholds_source": ({"evaluate.stop_rule_thresholds": sorted(th_over)} if th_over else "既定"),
           "hold_steps": hold,
           "hold_steps_source": ("evaluate.stop_rule_hold_steps" if evaluate.get("stop_rule_hold_steps") is not None
                                 else f"既定 {HOLD_STEPS_DEFAULT}"),
           "extract": {**handoff_extract_settings(), "source": "feedback.deltastar_loop.handoff_extract_settings (§5.1 #36e)"}}
    if not requested:
        rec.update(use=False, reason="使わない指定", cap=None, cap_source=None)
        return rec
    if nsteps_arg:
        cap, cap_src = int(nsteps_arg), "引数 (--steps)"
    elif yaml_cap is not None:
        cap, cap_src = yaml_cap, "evaluate.stop_rule_cap"
    else:
        cap, cap_src = CAP_DEFAULT, f"既定 {CAP_DEFAULT}"
    why = []
    if not euler_ref:
        why.append("euler_ref が無い (y_b を決められない)")
    if int(out_int) * (th["min_points"] - 1) > th["window"]:
        why.append(f"出力の間隔 {int(out_int)} step では窓 {th['window']} step に {th['min_points']} 点が入らない "
                   f"(evaluate.outStepInterval ≤ {th['window'] // (th['min_points'] - 1)})")
    if int(cap) % int(out_int):
        why.append(f"上限 {int(cap)} step が出力の間隔 {int(out_int)} step の倍数でない (上限の出力が書かれず CENSORED を見分けられない)")
    rec.update(use=not why, reason=("; ".join(why) or None), cap=cap, cap_source=cap_src)
    if cap < th["window"] + int(out_int) + hold:
        rec["note"] = (f"上限 {cap} では最初の出力から窓 {th['window']} step を覆って保持 {hold} step を終えられない (必ず CENSORED)")
    return rec


# ------------------------------------------------------------------ 量
def _read_coords(run_dir: Path):
    with h5py.File(run_dir / "nozzle.h5", "r") as h:
        return np.array(h["/MESH/COORD"][:]).reshape(-1, 3)


def read_values(path) -> dict:
    """res の VALUE の全データセット (保存の型のまま)。"""
    with h5py.File(path, "r") as h:
        return {k: np.array(d[:]) for k, d in h["VALUE"].items()}


def nonfinite_count(V: dict) -> int:
    return int(sum(np.count_nonzero(~np.isfinite(a)) for a in V.values() if a.dtype.kind == "f"))


def mass_flows(X, R, RO, UX, UY, S: float) -> tuple:
    """入口・出口の列の折れ線を通る質量流量 [kg/s] (`cold_xcheck.reduce_fields` の式; X, R は r_t 単位の (ni, nj))。"""
    out = []
    for i in (0, X.shape[0] - 1):
        xx = X[i] * S; rr = R[i] * S
        a = RO[i] * UX[i] * rr; b = RO[i] * UY[i] * rr
        out.append(2 * np.pi * float(np.sum(0.5 * (a[1:] + a[:-1]) * np.diff(rr) - 0.5 * (b[1:] + b[:-1]) * np.diff(xx))))
    return out[0], out[1]


def theta_r_sections(X, R, RO, UX, yb_x, yb, xs_list) -> list:
    """一定 x の断面 (各行 j を x 方向に補間) で、帯の外縁 y_b の (ρ_e, u_e) を縁にした運動量厚さ θ_r (`cold_xcheck.reduce_fields` の式、
    断面は xs に最も近い列)。戻り [(x_断面, θ_r)]。y_b が軸を越える断面は NaN。"""
    xt = X[:, 0]
    nj = X.shape[1]
    out = []
    for xs in xs_list:
        i = int(np.argmin(np.abs(xt - xs)))
        x0 = float(xt[i])
        rr = np.array([np.interp(x0, X[:, j], R[:, j]) for j in range(nj)])
        ro = np.array([np.interp(x0, X[:, j], RO[:, j]) for j in range(nj)])
        ux = np.array([np.interp(x0, X[:, j], UX[:, j]) for j in range(nj)])
        rw = rr[-1]; rb = rw - float(np.interp(x0, yb_x, yb))
        if not np.isfinite(rb) or rb <= rr[0]:
            out.append((x0, float("nan")))
            continue
        rf = np.linspace(rb, rw, 4001)
        rho = np.interp(rf, rr, ro); u = np.interp(rf, rr, ux)
        re_, ue_ = rho[0], u[0]
        qq = np.trapezoid(rho * u * (ue_ - u) * rf, rf)
        out.append((x0, float(rw - np.sqrt(max(rw ** 2 - 2.0 * qq / (re_ * ue_ ** 2), 0.0)))))
    return out


def recorded_quantities(Xf, Rf, V: dict, x_E: float) -> dict:
    """記録だけの量: 出口コア M と無次元の流量 (`nozzle_report.metrics`・`qs_series.py` と同じ定義; 座標は nozzle_report.load_field と同じ
    float32 の COORD / S、VALUE は float64)。"""
    ux = V["Ux"].astype(float); uy = V["Uy"].astype(float)
    M = np.hypot(ux, uy) / V["sonic"].astype(float)
    eta = Rf[-1] / Rf[-1, -1]; core = (eta >= 0.05) & (eta <= 0.7)
    q = V["ro"].astype(float) * ux
    m = 2 * np.pi * np.trapezoid(q * Rf, Rf, axis=1)
    xs = Xf[:, 0]; sel = (xs > -2.5) & (xs < x_E)
    return {"exit_core_M": float(M[-1][core].mean()), "mdot_rt": float(np.median(m[sel]))}


class RunGeometry:
    """run の格子 (prepare_info・nozzle.h5) を一度だけ読む。"""

    def __init__(self, run_dir):
        self.run_dir = Path(run_dir)
        info = json.loads((self.run_dir / "prepare_info.json").read_text())
        self.info = info
        self.ni, self.nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
        self.S = float(info["scale_m"]); self.x_E = float(info["x_E"])
        pw = info.get("physical_wall") or {}
        dom = pw.get("domain") if isinstance(pw, dict) else None
        if not (isinstance(dom, (list, tuple)) and len(dom) == 2):
            raise ValueError(f"{self.run_dir}: prepare_info.physical_wall.domain が無い (θ_r の断面の基準 = 設計壁の出口を決められない)")
        self.x_wall_exit = float(dom[1])     # θ_r の断面の基準 (設計壁の出口)。x_E (軸 Mach 則の終点) ではない
        nc = _read_coords(self.run_dir)
        if nc.shape[0] != self.ni * self.nj:
            raise ValueError(f"{self.run_dir}: 節点数 {nc.shape[0]} が ni × nj = {self.ni * self.nj} でない (構造格子の node run でない)")
        # θ_r・r は M6 と同じ float64 の座標、記録の量は nozzle_report.load_field と同じ float32 の COORD / S
        self.X = (nc[:, 0].astype(float) / self.S).reshape(self.ni, self.nj); self.R = (nc[:, 1].astype(float) / self.S).reshape(self.ni, self.nj)
        self.Xf = (nc[:, 0] / self.S).reshape(self.ni, self.nj); self.Rf = (nc[:, 1] / self.S).reshape(self.ni, self.nj)


def snapshot(geo: RunGeometry, step: int, band: dict | None, x_fracs=X_FRACS, values: dict | None = None) -> dict:
    """出力 res_<step>.h5 の 1 行。band が None (停止規則を使えない) なら θ_r は出さない。非有限があれば量は出さない (DIVERGED の行)。"""
    V = values if values is not None else read_values(geo.run_dir / f"res_{step}.h5")
    for k in ("ro", "Ux", "Uy", "sonic"):
        if k not in V or V[k].shape[0] != geo.ni * geo.nj:
            raise ValueError(f"res_{step}.h5 の VALUE/{k} が無い・長さが違う")
    row = {"step": int(step), "nonfinite_all": nonfinite_count(V)}
    if row["nonfinite_all"]:
        return row
    sh = (geo.ni, geo.nj)
    RO = V["ro"].astype(float).reshape(sh); UX = V["Ux"].astype(float).reshape(sh); UY = V["Uy"].astype(float).reshape(sh)
    mi, mo = mass_flows(geo.X, geo.R, RO, UX, UY, geo.S)
    row.update(mdot_in=mi, mdot_out=mo, r=(mi - mo) / mi if mi != 0 else float("nan"))
    if band is not None:
        xs = [f * geo.x_wall_exit for f in x_fracs]
        for f, (x0, th) in zip(x_fracs, theta_r_sections(geo.X, geo.R, RO, UX, np.asarray(band["x"]), np.asarray(band["y_b"]), xs)):
            row[theta_key(f)] = th; row[f"x_{theta_key(f)}"] = x0
    row.update(recorded_quantities(geo.Xf, geo.Rf, {k: V[k].reshape(sh) for k in ("ro", "Ux", "Uy", "sonic")}, geo.x_E))
    return row


# ------------------------------------------------------------------ 帯の位置 y_b (本段の最初の出力 × Euler 参照)
def resolve_path(p) -> Path | None:
    """prepare_info の euler_ref (相対なら cwd、無ければリポジトリ直下からの相対) を解決する。無ければ None。固定のスナップショットの指定
    {run, res, sha256} (§5.1 #36e) は run のディレクトリを解決する (sha の照合は `check_euler_ref`)。"""
    if isinstance(p, dict):
        p = p.get("run")
    if not p:
        return None
    cand = [Path(p)]
    if not Path(p).is_absolute():
        cand.append(Path(__file__).resolve().parents[3] / p)
    for c in cand:
        if c.is_dir():
            return c.resolve()
    return None


def euler_ref_spec(euler_ref):
    """prepare_info の euler_ref を抽出器に渡す指定にする: run dir (str、`resolve_path` で解決) か固定のスナップショット
    {run (解決した絶対パス), res, sha256}。見つからなければ None。"""
    eu = resolve_path(euler_ref)
    if eu is None:
        return None
    if isinstance(euler_ref, dict):
        return {"run": str(eu), "res": euler_ref.get("res"), "sha256": euler_ref.get("sha256")}
    return eu


def check_euler_ref(euler_ref) -> dict:
    """euler_ref を解決して固定のスナップショットなら sha256 を照合する (§5.1 #36e (6))。見つからない・sha が違えば ValueError。
    戻り: `deltastar.resolve_euler_ref` の戻り。"""
    from ..metrics.deltastar import resolve_euler_ref
    spec = euler_ref_spec(euler_ref)
    if spec is None:
        raise ValueError(f"euler_ref が無い・見つからない ({euler_ref!r})")
    return resolve_euler_ref(spec, verify=True)


def band_extract_kw(sr: dict) -> dict:
    """帯の位置の抽出 (`fix_band`) に渡す条件: 受け渡しの抽出と同じ一定 x の断面への補間と端の扱い (prepare_info の stop_rule.extract、
    無ければ今の既定 `handoff_extract_settings`)。帯の位置は抽出器の縁の推定 (δ_r の平滑化・緩和の前) で決まるので、ω・平滑化・昇格の上限は
    渡さない (昇格の上限で帯の位置が使えなくならないように)。"""
    from ..feedback.deltastar_loop import handoff_extract_settings
    ex = sr.get("extract") or handoff_extract_settings()
    return {k: ex[k] for k in ("remap_constant_x", "remap_edge") if k in ex}


def fix_band(run_dir, step: int, euler_ref, out_dir, **extract_kw) -> dict:
    """本段の最初の出力 res_<step>.h5 を Euler 参照で `extract_and_merge(.., band_select="edge")` と同じ抽出にかけ、帯の位置 y_b(x) [r_t] を返す
    (`delta_r_edge.csv` の y_b = 縁 y* の edge_c (1.25) 倍; さらに掛けない)。抽出は run の他の出力を読まないよう、prepare_info.json・
    nozzle.h5・その 1 つの res だけを一時ディレクトリに symlink して行う (`_load_structured` は最新の res を読むため)。
    使えない (euler_ref が無い・抽出が失敗・y_b が有限でない) ときは Unavailable。euler_ref が固定のスナップショット {run, res, sha256} なら
    その res を読み、sha が合わなければ Unavailable (§5.1 #36e (6))。extract_kw は `band_extract_kw` (一定 x の断面への補間・端の扱い)。"""
    run_dir = Path(run_dir).resolve(); od = Path(out_dir)
    eu = euler_ref_spec(euler_ref)
    if eu is None:
        raise Unavailable(f"euler_ref が無い・見つからない ({euler_ref!r})")
    res = run_dir / f"res_{step}.h5"
    t0 = time.time()
    try:
        from ..feedback.deltastar_loop import extract_and_merge
        with tempfile.TemporaryDirectory(prefix="stop_rule_band_") as td:
            src = Path(td) / run_dir.name
            src.mkdir()
            for name in ("prepare_info.json", "nozzle.h5"):
                (src / name).symlink_to(run_dir / name)
            (src / res.name).symlink_to(res)
            extract_and_merge(src, eu, band_select="edge", out_dir=od, **extract_kw)
        edge = np.loadtxt(od / "delta_r_edge.csv", delimiter=",", skiprows=1, ndmin=2)
    except Exception as e:  # noqa: BLE001 — 抽出の失敗は「使えない」として記録する (例外として外へ出さない)
        raise Unavailable(f"抽出が失敗: {type(e).__name__}: {e}") from e
    x, yb = edge[:, 0], edge[:, 5]
    if len(x) < 2 or not np.all(np.isfinite(yb)) or not np.all(yb > 0):
        raise Unavailable(f"帯の位置 y_b が有限・正でない (点 {len(x)}、非有限 {int(np.count_nonzero(~np.isfinite(yb)))})")
    return {"x": x.tolist(), "y_b": yb.tolist(), "step": int(step), "res_sha256": _sha256(res),
            "euler_ref": (eu if isinstance(eu, dict) else str(eu)),
            "extract": {"function": "feedback.deltastar_loop.extract_and_merge", "band_select": "edge", **{k: repr(v) for k, v in extract_kw.items()}},
            "out_dir": str(od), "seconds": round(time.time() - t0, 1)}


# ------------------------------------------------------------------ 出力ごとの δ_E・δ_target の記録 (§5.1 #36e (7); 判定には使わない)
DELTA_WAIT_S = 7200.0      # 本段の後に残りの抽出を待つ上限 [s] (抽出 1 本 ~2 分 × 追いつかない分; 飛ばす規則で 2W/出力の間隔 + 1 本まで)


def use_columns(run_dir) -> dict:
    """δ の記録の列: NS の列の x (抽出器の `_load_structured` と同じ式 COORD / S) のうち使用域 [0, x_w] (x_w = 設計壁の出口
    physical_wall.domain[1]) の列。戻り {"x" (r_t、列の昇順), "i" (列の番号), "x_w", "scale_m"}。"""
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S = float(info["scale_m"])
    xw = float((info.get("physical_wall") or {})["domain"][1])
    with h5py.File(run_dir / "nozzle.h5", "r") as h:
        nc = h["/MESH/COORD"][:].reshape(-1, 3)
    x = (nc[:, 0].reshape(ni, nj) / S)[:, 0]
    sel = np.where((x >= 0.0) & (x <= xw + 1e-9))[0]
    return {"x": x[sel], "i": sel, "x_w": xw, "scale_m": S}


def delta_snapshot(run_dir, step: int, euler_ref, extract: dict, x_use) -> dict:
    """出力 res_<step>.h5 1 つを受け渡しの抽出と同じ条件 (`extract` = `handoff_extract_settings`: edge・一定 x の断面への補間・端は欠測・
    ω・平滑化の λ・昇格の上限、固定した Euler の参照) で抽出し、使用域の列 x_use の δ_E (未緩和の平滑化抽出)・δ_target (緩和・再平滑化の
    後の次の壁の入力)・held を返す (r_t 単位、抽出しなかった列は NaN と `missing`)。抽出は run の他の出力を読まないよう一時ディレクトリに
    prepare_info.json・nozzle.h5・その res だけを symlink して行い、出力の表は捨てる (行の値だけを stop_rule.json に残す)。"""
    from ..feedback.deltastar_loop import extract_and_merge, read_delta_r_next
    run_dir = Path(run_dir).resolve()
    res = run_dir / f"res_{step}.h5"
    t0 = time.time()
    kw = {k: v for k, v in extract.items() if k != "source"}
    with tempfile.TemporaryDirectory(prefix="stop_rule_delta_") as td:
        src = Path(td) / run_dir.name
        src.mkdir()
        for name in ("prepare_info.json", "nozzle.h5"):
            (src / name).symlink_to(run_dir / name)
        (src / res.name).symlink_to(res)
        summ = extract_and_merge(src, euler_ref, out_dir=Path(td) / "out", **kw)
        nx = read_delta_r_next(Path(td) / "out" / "delta_r_next.csv")
    x_use = np.asarray(x_use, dtype=float)
    xe = np.atleast_1d(np.asarray(nx["x_rt"], dtype=float))
    if len(xe) == 0:
        raise ValueError(f"res_{step}.h5 の抽出の表が空")
    # 列の対応は x の一致で取る (抽出器の x と use_columns の x は同じ式 COORD / S。表の x は '%.18e' で往復しても同じ値)
    k = np.clip(np.searchsorted(xe, x_use), 0, len(xe) - 1)
    hit = xe[k] == x_use
    dE = np.where(hit, np.atleast_1d(nx["delta_E"])[k], np.nan)
    dT = np.where(hit, np.atleast_1d(nx["delta_target"])[k], np.nan)
    held = np.where(hit, np.atleast_1d(nx["held"])[k] != 0, False)
    return {"delta_E": dE.tolist(), "delta_target": dT.tolist(), "held": [int(v) for v in held],
            "missing": [int(j) for j in np.where(~hit)[0]], "n_missing": int((~hit).sum()), "n_held": int(held.sum()),
            "lam": summ["smooth"]["lam"], "res_sha256": _sha256(res), "seconds": round(time.time() - t0, 2)}


class DeltaRecorder:
    """出力ごとの δ_E・δ_target の抽出 (`delta_snapshot`) を別スレッド 1 本で回す (停止の判定を抽出で待たせない)。submit した順に処理し、
    最新の submit から span (= 2 窓 2W) より古い出力は抽出せず `skipped` にする (停止の時刻 s は最新の出力以後なので [s − 2W, s] に入らない)。
    例外は行の `error` に残す (記録の量なので見張りは止めない)。"""

    def __init__(self, run_dir, euler_ref, extract: dict, x_use, span: int) -> None:
        import queue
        import threading
        self.run_dir, self.euler_ref, self.extract, self.x_use, self.span = Path(run_dir), euler_ref, dict(extract), list(x_use), int(span)
        self._q = queue.Queue()
        self._lock = threading.Lock()
        self._done: dict = {}
        self._latest = None
        self._stop = False
        self._t = threading.Thread(target=self._work, name="stop_rule_delta", daemon=True)
        self._t.start()

    def submit(self, step: int) -> None:
        with self._lock:
            self._latest = int(step) if self._latest is None else max(self._latest, int(step))
        self._q.put(int(step))

    def _work(self) -> None:
        while True:
            step = self._q.get()
            try:
                if step is None:
                    return
                with self._lock:
                    latest, stop = self._latest, self._stop
                if stop:
                    out = {"error": "見張りが終わったので抽出しない"}
                elif latest is not None and step < latest - self.span:
                    out = {"skipped": f"最新の出力 {latest} から 2 窓 {self.span} step より古い (受け渡しの 2 窓に入らない)"}
                else:
                    try:
                        out = delta_snapshot(self.run_dir, step, self.euler_ref, self.extract, self.x_use)
                    except Exception as e:  # noqa: BLE001 — 記録の量の抽出の失敗は行に残す (見張りは止めない)
                        out = {"error": f"{type(e).__name__}: {e}"}
                with self._lock:
                    self._done[step] = out
            finally:
                self._q.task_done()

    def collect(self) -> dict:
        """終わった抽出 {step: 結果} を取り出す (取り出した分は消す)。"""
        with self._lock:
            out, self._done = self._done, {}
        return out

    def finish(self, wait_s: float, abandon: bool = False) -> dict:
        """残りの抽出を wait_s まで待って止める (abandon なら待たずに未着手を捨てる)。戻り: 取り出していない結果。"""
        if abandon:
            with self._lock:
                self._stop = True
        self._q.put(None)
        self._t.join(timeout=max(float(wait_s), 0.0))
        out = self.collect()
        out["_unfinished"] = bool(self._t.is_alive())
        return out


def _attach_deltas(rows: list, got: dict) -> int:
    """`DeltaRecorder` の結果を行の `delta` に付ける。付けた数を返す。"""
    by = {r["step"]: r for r in rows}
    n = 0
    for st, d in got.items():
        if isinstance(st, int) and st in by:
            by[st]["delta"] = d
            n += 1
    return n


# ------------------------------------------------------------------ 判定
def judge_rows(rows: list, k: int | None = None, *, window=WINDOW, min_points=MIN_POINTS, r_mean_max=R_MEAN_MAX, r_amp_max=R_AMP_MAX,
               theta_drift_max=THETA_DRIFT_MAX, theta_amp_max=THETA_AMP_MAX, x_fracs=X_FRACS) -> dict:
    """rows[k] (step 昇順、本段の出力、step > 0) の時点の判定と内訳。ok なら停止の条件を満たす。"""
    k = len(rows) - 1 if k is None else int(k)
    s_end = rows[k]["step"]
    win = [r for r in rows[:k + 1] if r["step"] >= s_end - window]
    out = {"step": s_end, "n_window": len(win), "history_span": s_end - rows[0]["step"], "covers": s_end - rows[0]["step"] >= window,
           "ok": False}
    if not out["covers"]:
        out["why"] = f"本段の履歴 {out['history_span']} step が窓 {window} step を覆わない"
        return out
    if len(win) < min_points:
        out["why"] = f"窓の点数 {len(win)} < {min_points}"
        return out
    s = np.array([w["step"] for w in win], float)
    second = s >= s_end - window / 2.0           # 中点ちょうどの出力は後半
    first = ~second
    if not first.any() or not second.any():
        out["why"] = "窓の前半か後半に出力が無い"
        return out
    rv = np.array([w.get("r", np.nan) for w in win], float)
    out["r_mean"] = float(rv.mean()); out["r_amp"] = float((rv.max() - rv.min()) / 2.0)
    out["r_ok"] = bool(np.isfinite(rv).all() and abs(out["r_mean"]) <= r_mean_max and out["r_amp"] <= r_amp_max)
    ds = float(s[second].mean() - s[first].mean())
    th = {}
    for f in x_fracs:
        key = theta_key(f)
        v = np.array([w.get(key, np.nan) for w in win], float)
        mean = float(v.mean())
        d_raw = float(abs(v[second].mean() - v[first].mean()) / abs(mean)) if mean else float("nan")
        amp = float((v.max() - v.min()) / 2.0 / abs(mean)) if mean else float("nan")
        D = d_raw * window / ds
        th[key] = {"mean": mean, "D": D, "D_raw": d_raw, "amp": amp,
                   "ok": bool(np.isfinite(v).all() and D <= theta_drift_max and amp <= theta_amp_max)}
    out["theta"] = th; out["half_step_gap"] = ds
    out["ok"] = bool(out["r_ok"] and all(t["ok"] for t in th.values()))
    if not out["ok"]:
        bad = (["r"] if not out["r_ok"] else []) + [k_ for k_, t in th.items() if not t["ok"]]
        out["why"] = "条件を満たさない: " + ", ".join(bad)
    return out


def row_data_problem(row: dict, x_fracs=X_FRACS) -> str | None:
    """場が有限なのに判定量が有限でない (y_b が軸を越える・流量 0 など) → 理由。問題なければ None。"""
    keys = ["r"] + [theta_key(f) for f in x_fracs]
    bad = [k for k in keys if not np.isfinite(row.get(k, np.nan))]
    return f"判定量が有限でない: {', '.join(bad)}" if bad else None


def _hold_step(hs: dict, ok: bool, step: int, hold: int, why=None) -> bool:
    """到達の後の継続 (`hold_steps`; M6 の版で足した、既定 `HOLD_STEPS_DEFAULT`) の状態 hs = {start, reached, breaks} を出力 1 つぶん
    進め、REACHED にするかを返す。
    条件を満たした最初の出力 s₀ から数え、step − s₀ ≥ hold になった出力で True。途中で条件が外れたら到達を取り消して breaks に記録し、
    次に満たした出力から数え直す。hold 0 は条件を満たした出力で True (forge-f5 の版の「最初の到達で止める」と同じ)。"""
    if not ok:
        if hs.get("start") is not None:
            hs["breaks"].append({"start": hs["start"], "broken_at": int(step), "why": why})
            hs["start"] = None
        return False
    if hs.get("start") is None:
        hs["start"] = int(step)
        hs["reached"].append(int(step))
    return int(step) - hs["start"] >= int(hold)


def first_stop(rows: list, hold_steps: int = 0, **th) -> dict | None:
    """保存済みの系列の最初の終わり {status: REACHED / DIVERGED / DATA_ERROR, k, step, judge}。無ければ None。
    hold_steps > 0 は到達の後 hold_steps の間すべての出力で条件を保った出力を REACHED にし、`hold` に経過を付ける。**関数の既定 0 は
    最初の到達** (系列の再生用); 見張り (`watch_main_stage`) と同じ結果にするには prepare_info の stop_rule.hold_steps (既定
    `HOLD_STEPS_DEFAULT`) を渡す。"""
    fr = th.get("x_fracs", X_FRACS)
    hs = {"start": None, "reached": [], "breaks": []}
    for k, row in enumerate(rows):
        if row.get("nonfinite_all"):
            return {"status": "DIVERGED", "k": k, "step": row["step"], "judge": None}
        p = row_data_problem(row, fr)
        if p:
            return {"status": "DATA_ERROR", "k": k, "step": row["step"], "judge": None, "note": p}
        j = judge_rows(rows, k, **th)
        if _hold_step(hs, j["ok"], row["step"], hold_steps, j.get("why")):
            out = {"status": "REACHED", "k": k, "step": row["step"], "judge": j}
            if hold_steps:
                out["hold"] = {"hold_steps": int(hold_steps), **hs}
            return out
    return None


# ------------------------------------------------------------------ 記録・後処理
def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def save_record(run_dir, rec: dict) -> Path:
    """`<run>/stop_rule.json` を原子的に書く (一時ファイル → rename)。"""
    p = Path(run_dir) / RECORD
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=_json_default))
    os.replace(tmp, p)
    return p


def write_state(run_dir, state: dict) -> Path:
    """本段の終了の状態 `<run>/STOP_RULE_STATE.json` を原子的に書く (一時ファイル → rename; M6 の版で足した)。state は少なくとも status を持つ。
    `run_staged_ns` が本段の後に 1 度だけ書く。見張りの途中の記録 (stop_rule.json、status は RUNNING から変わる) と違い、このファイルが
    あれば本段は終わっている (無いのに prepare_info が停止規則を要求していれば、本段が終わっていないか記録が失われた)。"""
    if state.get("status") not in STATES + ("NOT_USED",):
        raise ValueError(f"停止の状態 {state.get('status')!r} は {STATES + ('NOT_USED',)} のどれか")
    p = Path(run_dir) / STATE
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps({"version": 1, "plan": PLAN, **state}, indent=1, ensure_ascii=False, default=_json_default))
    os.replace(tmp, p)
    return p


def read_state(run_dir) -> dict | None:
    """`<run>/STOP_RULE_STATE.json` を読む (無ければ None)。"""
    p = Path(run_dir) / STATE
    return json.loads(p.read_text()) if p.exists() else None


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def res_steps(run_dir) -> list:
    return sorted(int(m.group(1)) for f in Path(run_dir).iterdir() if (m := _RES_RE.fullmatch(f.name)))


def delete_after(run_dir, step: int) -> list:
    """止めた出力 step より後の出力 (res_<N>.h5/.xmf と境界の res_<名前>_<id>_<N>.h5/.xmf、N > step) を消す。消した名前を返す。"""
    gone = []
    for f in sorted(Path(run_dir).iterdir()):
        m = _OUT_RE.fullmatch(f.name)
        if m and int(m.group(1)) > step:
            f.unlink(); gone.append(f.name)
    return gone


def postprocess(run_dir) -> dict:
    """`run_case.sh` と同じ後処理 (check_convergence.py → CONVERGENCE_VERDICT.txt、plot_implicit_residuals.py → residual_history.png)。"""
    from .runner import FORGE_TOOLS, _ENV
    run_dir = Path(run_dir)
    with open(run_dir / "CONVERGENCE_VERDICT.txt", "w") as fh:
        r = subprocess.run(["python3", str(FORGE_TOOLS / "check_convergence.py"), "."], cwd=run_dir, env=_ENV, stdout=fh, stderr=subprocess.STDOUT)
    with open(run_dir / "CONVERGENCE_VERDICT.txt", "a") as fh:
        fh.write(f"CONVERGENCE_EXIT: {r.returncode}\n")
    subprocess.run(["python3", str(FORGE_TOOLS / "plot_implicit_residuals.py"), "--input", "residual_history.csv", "--output", "residual_history.png"],
                   cwd=run_dir, env=_ENV, capture_output=True)
    return {"by": "watcher", "convergence_exit": r.returncode}


# ------------------------------------------------------------------ forge のプロセス
def _proc_stat(pid: int):
    try:
        s = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return None
    rest = s[s.rfind(")") + 2:].split()
    return rest[0], int(rest[1])          # 状態, ppid


def _alive(pid: int) -> bool:
    st = _proc_stat(pid)
    return st is not None and st[0] not in ("Z", "X")


def _children(pid: int) -> list:
    out = []
    for d in Path("/proc").iterdir():
        if d.name.isdigit():
            st = _proc_stat(int(d.name))
            if st is not None and st[1] == pid:
                out.append(int(d.name))
    return out


def forge_exe() -> Path:
    """`run_case.sh` が起動する forge (FORGE_BIN、無ければ <FORGE_ROOT>/solver_density_cuda/build/forge)。"""
    from .runner import FORGE_ROOT, _ENV
    return Path(_ENV.get("FORGE_BIN") or (FORGE_ROOT / "solver_density_cuda" / "build" / "forge"))


def _is_forge(pid: int, exe: Path) -> bool:
    try:
        argv = Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
        cwd = os.readlink(f"/proc/{pid}/cwd")
    except OSError:
        return False
    want = os.path.realpath(exe)
    for a in argv[:2]:
        if a and os.path.realpath(os.path.join(cwd, a.decode(errors="replace"))) == want:
            return True
    try:
        return os.path.realpath(os.readlink(f"/proc/{pid}/exe")) == want
    except OSError:
        return False


def _forge_exit(run_dir) -> int | None:
    p = Path(run_dir) / "run_case_stdout.log"
    if not p.exists():
        return None
    m = _FORGE_EXIT_RE.findall(p.read_text(errors="replace"))
    return int(m[-1]) if m else None


def stop_forge(bash_pid: int, run_dir, term_wait_s: float = 180.0, kill_wait_s: float = 60.0, find_wait_s: float = 10.0) -> dict:
    """`run_case.sh` (bash_pid) の子の forge だけに SIGTERM → term_wait_s 待って残れば SIGKILL → kill_wait_s 待ち、止まったことを確かめる。
    forge だけを止めるので run_case.sh の後処理はそのまま走る。forge が見つからず、まだ終わった記録も無いときは (find_wait_s 待った後)
    プロセスグループごと止める (後処理は呼ぶ側が補う)。"""
    exe = forge_exe()
    t0 = time.time()
    pids = [p for p in _children(bash_pid) if _is_forge(p, exe)]
    while not pids and _forge_exit(run_dir) is None and _alive(bash_pid) and time.time() - t0 < find_wait_s:
        time.sleep(0.2)
        pids = [p for p in _children(bash_pid) if _is_forge(p, exe)]
    rec = {"forge_pids": pids, "signals": [], "group_killed": False}
    if not pids:
        if _forge_exit(run_dir) is None and _alive(bash_pid):
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(bash_pid, sig)
                except OSError:
                    pass
                rec["signals"].append(f"killpg {sig.name}")
                t1 = time.time()
                while _alive(bash_pid) and time.time() - t1 < kill_wait_s:
                    time.sleep(0.2)
                if not _alive(bash_pid):
                    break
            rec["group_killed"] = True
        rec["stopped"] = not _alive(bash_pid) or _forge_exit(run_dir) is not None
        return rec
    for sig, wait in ((signal.SIGTERM, term_wait_s), (signal.SIGKILL, kill_wait_s)):
        live = [p for p in pids if _alive(p)]
        if not live:
            break
        for p in live:
            try:
                os.kill(p, sig)
            except OSError:
                pass
        rec["signals"].append(sig.name)
        t1 = time.time()
        while any(_alive(p) for p in pids) and time.time() - t1 < wait:
            time.sleep(0.2)
    rec["stopped"] = not any(_alive(p) for p in pids)
    return rec


# ------------------------------------------------------------------ 見張り (run_staged_ns の本段)
def watch_main_stage(run_dir, sr: dict, *, poll_s: float = 10.0, settle_s: float = 60.0, term_wait_s: float = 180.0,
                     kill_wait_s: float = 60.0, post_wait_s: float = 1800.0, read_tries: int = 3, delta_wait_s: float = DELTA_WAIT_S) -> dict:
    """本段の forge を `run_case.sh` の子プロセスで起動し、停止規則で見張る。sr は prepare_info.json の stop_rule (use = True、cap が
    solverConfig の nStepOuter に書いてあること)。戻り: {status, stop_step, forge_rc, rc}。rc は呼ぶ側に返す終了コード
    (REACHED・CENSORED・UNAVAILABLE で ns_steps まで回った = 0、それ以外は非 0)。記録は `<run>/stop_rule.json`。

    書き終わった出力: それより後の出力がある、forge が終わった、または最後に書かれてから settle_s 以上。読めない出力は read_tries 回まで
    次の見回りで読み直す (それでも読めなければ DATA_ERROR)。"""
    from .runner import FORGE_TOOLS, _ENV
    run_dir = Path(run_dir).resolve()
    th = {k: sr[k] for k in ("window", "min_points", "r_mean_max", "r_amp_max", "theta_drift_max", "theta_amp_max", "x_fracs")}
    cap = int(sr["cap"]); fallback = int(sr["fallback_steps"])
    # 到達の後の継続 (0 = 最初の到達で止める)。prepare_info の stop_rule.hold_steps (resolve_settings の既定 HOLD_STEPS_DEFAULT)。
    # 統合 (2026-10-14) より前に準備した run の記録にキーが無ければ 0 (準備した時の動きのまま)
    hold = int(sr.get("hold_steps") or 0)
    pre = [n for n in res_steps(run_dir) if n > 0]
    if pre:
        raise RuntimeError(f"{run_dir}: 本段の起動前に出力 res_{pre[0]}.h5 ほか {len(pre)} 個がある — 停止規則の系列に混ざるので起動しない")
    geo = RunGeometry(run_dir)
    rec = {"version": 1, "plan": PLAN, "settings": sr, "status": "RUNNING",
           "started": _now(), "band": None, "unavailable": None, "rows": [], "judge_last": None, "tries": {},
           "timing": {"poll_s": poll_s, "settle_s": settle_s, "term_wait_s": term_wait_s, "kill_wait_s": kill_wait_s},
           "hold": ({"hold_steps": hold, "start": None, "reached": [], "breaks": []} if hold else None),
           # 要求の断面の位置 (x_fracs × x_w)。実際に選ばれた列の x は各行の x_theta_r_<比> (§6 V-h370)
           "sections": [{"key": theta_key(f), "x_frac": float(f), "x_requested": float(f) * geo.x_wall_exit} for f in th["x_fracs"]]}
    # 出力ごとの δ_E・δ_target の記録 (§5.1 #36e (7); 判定には使わない)。列は使用域 [0, x_w] の NS の列、条件は受け渡しの抽出と同じ
    from ..feedback.deltastar_loop import handoff_extract_settings
    ex = sr.get("extract") or {**handoff_extract_settings(), "source": "既定 (prepare_info の stop_rule に extract が無い)"}
    uc = use_columns(run_dir)
    rec["delta_series"] = {"x": uc["x"].tolist(), "use_domain": [0.0, uc["x_w"]], "scale_m": uc["scale_m"], "extract": ex,
                           "euler_ref": sr.get("euler_ref"), "span_steps": 2 * int(th["window"]),
                           "note": "判定には使わない記録の量。行の delta に δ_E (未緩和の平滑化抽出)・δ_target (緩和・再平滑化の後)・held [r_t]"}
    recorder = DeltaRecorder(run_dir, euler_ref_spec(sr.get("euler_ref")), ex, uc["x"], 2 * int(th["window"]))
    save_record(run_dir, rec)
    t_launch = time.time()
    with open(run_dir / "run_case_stdout.log", "w") as log:
        proc = subprocess.Popen([str(FORGE_TOOLS / "run_case.sh"), str(run_dir)], env=_ENV, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
    rec["run_case_pid"] = proc.pid
    stop = None                      # {"status", "step", ...}: 見張りが止めると決めた出力
    rows = rec["rows"]
    seen = set()

    def complete(n, steps_now, done):
        if done or any(m > n for m in steps_now) or _forge_exit(run_dir) is not None:
            return True
        try:
            return time.time() - (run_dir / f"res_{n}.h5").stat().st_mtime >= settle_s
        except OSError:
            return False

    try:
        while True:
            if _attach_deltas(rows, recorder.collect()):
                save_record(run_dir, rec)
            done = proc.poll() is not None
            steps_now = res_steps(run_dir)
            for n in [m for m in steps_now if m > 0 and m not in seen]:
                if not complete(n, steps_now, done):
                    break
                try:
                    V = read_values(run_dir / f"res_{n}.h5")
                    if rec["band"] is None and rec["unavailable"] is None:
                        try:
                            rec["band"] = fix_band(run_dir, n, sr.get("euler_ref"), run_dir / BAND_DIR, **band_extract_kw(sr))
                        except Unavailable as e:
                            rec["unavailable"] = {"reason": str(e), "step": n,
                                                  "action": f"停止規則を使わず fallback_steps {fallback} step で止める"}
                        save_record(run_dir, rec)
                    row = snapshot(geo, n, rec["band"], th["x_fracs"], values=V)
                except Exception as e:  # noqa: BLE001
                    rec["tries"][str(n)] = rec["tries"].get(str(n), 0) + 1
                    rec["last_read_error"] = f"res_{n}.h5: {type(e).__name__}: {e}"
                    if rec["tries"][str(n)] >= read_tries:
                        stop = {"status": "DATA_ERROR", "step": n, "note": f"出力が読めない ({rec['last_read_error']})",
                                "read_error": True, "forge_done": done or _forge_exit(run_dir) is not None}
                    save_record(run_dir, rec)
                    break
                seen.add(n)
                if row["nonfinite_all"]:
                    stop = {"status": "DIVERGED", "step": n, "note": f"非有限 {row['nonfinite_all']} 個"}
                elif rec["unavailable"] is None:
                    p = row_data_problem(row, th["x_fracs"])
                    if p:
                        stop = {"status": "DATA_ERROR", "step": n, "note": p}
                rows.append(row)
                if not row["nonfinite_all"]:
                    recorder.submit(n)
                if stop is None and rec["unavailable"] is None:
                    j = judge_rows(rows, **th)
                    row["judge"] = j
                    rec["judge_last"] = j
                    if hold:
                        if _hold_step(rec["hold"], j["ok"], n, hold, j.get("why")):
                            stop = {"status": "REACHED", "step": n,
                                    "note": f"最初の到達 {rec['hold']['start']} の後 {hold} step のすべての出力で条件を保った"}
                    elif j["ok"]:
                        stop = {"status": "REACHED", "step": n}
                if stop is None and rec["unavailable"] is not None and n >= fallback:
                    stop = {"status": "UNAVAILABLE", "step": n, "note": f"停止規則を使えないので fallback_steps {fallback} の出力で止めた"}
                save_record(run_dir, rec)
                if stop is not None:
                    break
            if stop is not None:
                break
            if done and all(m in seen for m in res_steps(run_dir) if m > 0):
                break
            time.sleep(poll_s)
    except BaseException as e:
        # 見張りの例外でも forge を止めて状態を残す
        _attach_deltas(rows, recorder.finish(0.0, abandon=True))
        rec["stop_forge"] = stop_forge(proc.pid, run_dir, term_wait_s, kill_wait_s) if proc.poll() is None else None
        rec.update(status="EXEC_ERROR", note=f"見張りの例外: {type(e).__name__}: {e}", finished=_now())
        save_record(run_dir, rec)
        if not isinstance(e, Exception):
            raise
        return {"status": "EXEC_ERROR", "stop_step": None, "forge_rc": _forge_exit(run_dir), "rc": 1}

    t_stop = time.time()
    if stop is not None and proc.poll() is None:
        rec["stop_forge"] = stop_forge(proc.pid, run_dir, term_wait_s, kill_wait_s)
        if not rec["stop_forge"]["stopped"]:
            stop = {"status": "EXEC_ERROR", "step": stop["step"], "note": f"forge を止められない (止める理由: {stop['status']})"}
    try:
        proc.wait(timeout=post_wait_s)
    except subprocess.TimeoutExpired:
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(proc.pid, sig)
            except OSError:
                pass
            try:
                proc.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                continue
        rec["run_case_killed"] = True
    rec["run_case_rc"] = proc.returncode
    forge_rc = _forge_exit(run_dir)
    rec["forge_rc"] = forge_rc
    # run_case.sh の後処理が止めた後にも走ったか (CONVERGENCE_VERDICT.txt が止めた後に書かれたか)。走っていなければ同じ後処理を呼ぶ
    vp = run_dir / "CONVERGENCE_VERDICT.txt"
    ran = (vp.exists() and vp.stat().st_mtime >= (t_stop if stop is not None else t_launch) - 1.0
           and "CONVERGENCE_EXIT" in vp.read_text(errors="replace"))
    rec["postprocess"] = {"by": "run_case.sh"} if ran else postprocess(run_dir)

    last = max([n for n in res_steps(run_dir) if n > 0], default=0)
    if stop is not None and stop.get("read_error") and stop.get("forge_done") and forge_rc not in (None, 0):
        # forge が異常終了した後に読めない出力 (書きかけ) が残った: 異常終了を DATA_ERROR と取り違えない
        status, step = "EXEC_ERROR", last
        rec["note"] = f"forge が上限 {cap} より前に異常終了 (forge rc {forge_rc}); 書きかけの出力が読めない ({stop['note']})"
    elif stop is not None:
        status, step = stop["status"], stop["step"]
        if stop.get("note"):
            rec["note"] = stop["note"]
    elif forge_rc == 0 and last >= cap:
        status, step = ("UNAVAILABLE" if rec["unavailable"] is not None else "CENSORED"), last
        if status == "CENSORED" and hold and rec["hold"]["start"] is not None:
            rec["note"] = f"到達 {rec['hold']['start']} の後の {hold} step の保持が上限 {cap} までに終わらない (保持の途中で上限)"
    elif forge_rc not in (None, 0) and sorted(run_dir.glob("res_nan_*.h5")):
        status, step = "DIVERGED", None
        nan_files = sorted(run_dir.glob("res_nan_*.h5"))
        rec["note"] = "forge の detectNaN が止めた (" + ", ".join(f.name for f in nan_files) + ")"
        st_nan = [int(m.group(1)) for f in nan_files if (m := re.fullmatch(r"res_nan_(\d+)\.h5", f.name))]
        rec["fail_step"] = min(st_nan) if st_nan else None
    else:
        status, step = "EXEC_ERROR", last
        rec["note"] = (f"forge が上限 {cap} より前に終わった (forge rc {forge_rc}、最後の出力 {last})"
                       + (f"; 未読の出力あり ({rec.get('last_read_error')})" if rec.get("last_read_error") else ""))
    if status in ("REACHED", "DIVERGED", "UNAVAILABLE") and step is not None and stop is not None:
        rec["deleted_after_stop"] = delete_after(run_dir, step)
    if step is not None and (run_dir / f"res_{step}.h5").exists() and status in ("REACHED", "DIVERGED", "UNAVAILABLE", "CENSORED"):
        rec["stop_sha256"] = _sha256(run_dir / f"res_{step}.h5")
    # 残りの δ の抽出を待って行に付ける (止めた出力より後の出力は行に無い。抽出の時間は delta_series に)
    left = recorder.finish(delta_wait_s)
    rec["delta_series"]["unfinished_at_end"] = left.pop("_unfinished")
    _attach_deltas(rows, left)
    secs = [r["delta"]["seconds"] for r in rows if isinstance(r.get("delta"), dict) and "seconds" in r["delta"]]
    rec["delta_series"]["seconds"] = ({"n": len(secs), "mean": float(np.mean(secs)), "max": float(np.max(secs))} if secs else {"n": 0})
    rec.update(status=status, stop_step=step, last_output=last, finished=_now())
    if status == "REACHED":
        rec["recorded_at_stop"] = window_summary(rows, step, th["window"])
    save_record(run_dir, rec)
    ok = status in ("REACHED", "CENSORED") or (status == "UNAVAILABLE" and (stop is not None or forge_rc == 0))
    return {"status": status, "stop_step": step, "forge_rc": forge_rc, "rc": 0 if ok else (forge_rc if forge_rc not in (None, 0, 143, 137) else 1)}


def window_summary(rows: list, step: int, window: int = WINDOW) -> dict:
    """記録だけの量 (出口コア M・流量) の、止めた出力までの窓の平均 ± 振幅 (判定には使わない; §4.6 の記録)。"""
    win = [r for r in rows if step - window <= r["step"] <= step and "exit_core_M" in r]
    out = {}
    for k in ("exit_core_M", "mdot_rt"):
        v = np.array([r[k] for r in win], float)
        if len(v):
            out[k] = {"mean": float(v.mean()), "amp": float((v.max() - v.min()) / 2.0), "n": int(len(v))}
    return out
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
