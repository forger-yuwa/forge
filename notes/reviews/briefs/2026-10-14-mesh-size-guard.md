# 諮問ブリーフ: メッシュの規模を抑える (plan tooling-nozzle-mesh-size-guard §4・§6)

AGENTS.md 条件 1 (plan §4・§6 を新規に書く)。判断役のスイッチが codex。

## 観測事実

- ユーザ (2026-10-14): 別の AI が skill で風洞ノズルを設計させたところ、黙って 50 万点を超える Euler の格子を作った。
- 原因の見立て (コードで確認): 要求カードの生成器 (`design/forge_design/intake/wt_card.py:562-563`) は `mesh_euler.ni` を None にして、
  NS の `mesh.ni` (壁際のセルの AR ≤ 0.97 × 5000 から見積もる。上限 `NI_CAP` 6000、`wt_card.py:67`) を写す。nj は両方 97。
  NS の第 1 層が薄い (Re が高い) と ni は上限 6000 に張り付き、Euler も 6000 × 97 ≈ 58 万点。Euler の格子は壁に寄せない全域 0.005 の等比で AR の制約は無い。
- 写した経緯: 雛形の注記「E4 は生産の G1 と同じ 2000 点」(`design/forge_design/intake/recipes/axismach_explore_v1.yaml:48-49`)。
  case/45 の Euler の評価 (plan verification-case45-euler-total-enthalpy) で NS の G1 に合わせただけで、その plan は「2000 × 97 は今回の評価条件で、汎用の既定にする根拠はない」と書いている。
- 節点数の上限はどこにも無い。
- δ* の抽出 (`design/forge_design/metrics/deltastar.py:465-473` `euler_q_at`) は Euler の ρU と壁の半径を NS の断面の x へ列の間で線形に補間する。
- **V1a (実測 2026-10-14)**: Euler `case/70.m5_air_wt_demo/run_0012b_euler_c2b_cont12k` (ni 1800) の列を 2・3・4 本おきに間引いて (Euler の解は同じで、抽出が使う列だけ)
  NS `run_0010b_ns_c2b_cont60k` (ni 1800、冷却壁 SST) を抽出 (`extract_and_merge(.., band_select="edge")`): 使用域 [0, x_e = 22.23] の max|Δδ_E| は
  1.7e-5・2.2e-5・2.7e-5 r_t (x ≈ 6〜11)、流量の比 1.0076108 の変化 ≤ 1.5e-7。台本 `case/70.m5_air_wt_demo/euler_ni_interp_ab.py`。
- 既知の別の感度: Euler の反復数だけを変えた A/B (36k vs 48k) で δ_E は縮流部 (x ≈ −8) で 7.2e-3 r_t 動く (使用域の外)。
- Euler の所要: 1790 × 97 で 3.56 ms/step、48,000 step で約 3 分 (手元の GPU)。

## 方針 (plan §4) と検証 (§6)

- §4.1: 生成器で Euler の ni = min(NS の ni, 2000)。V1b で Euler の解そのものの ni 依存 (1800・900・600、同じ段取り 48,000 step) を測り、
  比 3 まで合格なら上限 2000、比 2 だけなら max(min(ni,2000), ⌈ni/2⌉)、どちらも不合格なら切り離さない。判定 δ_E ≤ 1e-4 r_t・出口コア M ≤ 2.5e-4・流量の比 ≤ 1e-4。
- §4.2: prepare (Euler・NS) で run dir を作る前に節点数 = ni × nj を数え、既定の上限 (Euler 25 万・NS 50 万、2D) を超えたら止める。
  問題 YAML の `max_nodes` と `max_nodes_reason` があればその値まで許す。
- §4.3: AGENTS.md と skill に「メッシュを作る前に節点数と所要の見込みをユーザに示し、雛形の規模 (Euler 15.5 万・NS 20 万) を超えるなら了承を取る」。

## 問い

**この方針のうち、実装後に失敗する可能性が最も高い点を 1 つ挙げ、防ぐ変更を 1 つ推奨してほしい。** 失敗の形は
「Euler の ni を減らしたことで δ* の帰還 (δ_E・流量の比・r_t の固定点) が黙って変わる」か「上限が素通りされる/正当な作業を止める」。
候補 (縛られなくてよい): (a) V1b の判定量が足りない (スロート付近の Euler の解・r_t の固定点の入力 = 流量の比・`snap_table_ends` の端の値)、
(b) 比 2・3 の試験が ni 1800 の 1 ケースだけで、NS 6000 のケース (高 Re) に外挿できない、(c) `max_nodes_reason` は AI がいくらでも書けるので安全網にならない
(ユーザの了承を機械的に確かめる方法が要るか)、(d) 上限の値 (25 万・50 万) の根拠、(e) スロートの x 方向の細分 (throat_refine 4・width 3) が ni とともに粗くなる。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,114p' plans/active/tooling-nozzle-mesh-size-guard.md`
- `sed -n '545,575p' design/forge_design/intake/wt_card.py` / `sed -n '60,70p' design/forge_design/intake/wt_card.py`
- `sed -n '40,62p' design/forge_design/intake/recipes/axismach_explore_v1.yaml`
- `sed -n '700,800p' design/forge_design/evaluate/runner_axismach.py` / `sed -n '840,870p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '455,520p' design/forge_design/metrics/deltastar.py` / `sed -n '90,130p' design/forge_design/metrics/deltastar.py`
- `cat case/70.m5_air_wt_demo/euler_ni_interp_ab.py`
- `grep -n "2000 × 97\|mesh_euler" plans/accepted/verification-case45-euler-total-enthalpy.md | head -20`
