# 諮問ブリーフ: 長さ合わせののこぎり歯の直し方 (終端 C⁺ の補間) と #6d の扱い (2026-10-16)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件:
  - 1: plan §4 の設計方針 (§4.5) と §6 の検証計画 (V9) を新規に書いた。
  - 4: plan §5.1 #6e は「直し方は上位に諮る」としていた。
- plan: `plans/active/tooling-nozzle-card-local-mesh.md`
  - §4.5 (今回の案)、§6 V9 (事前登録)、§5.1 #6d・#6e
- commit: `a36a87ff` (feature/card-local-mesh)。今回の案はまだ実装していない (試作だけ)。

## 観測事実

### 原因 (forge-7d の掃引、記録 `case/71.m55_combustion_wt/_runs/lcscan_forge7d/README.md`)

- `wall_mode: cplus`・`exit_mode: characteristic` の出口 F は、x_E に**最も近い**軸の節点 m_E の網の C⁺ と壁の交点
  (`design/forge_design/geometry/moc_inverse.py` `_design_cplus`、773 行付近)。
- カード B で L_c を固定して掃くと x_F(L_c) はのこぎり歯 (歯の中の傾き 0.29、段の高さは軸の網の間隔 1 つ、段の位置は x_E が節点の中点をまたぐ所)。
- from_length (`design/forge_design/evaluate/runner_axismach.py` 619〜676 行) は傾き 1 を仮定した固定点反復で、改善が止まると止まる。
  許容 (カードの公差 1 mm → L_total_tol = 0.5 mm / r_t ≈ 0.0048 r_t) に届くかは位相しだい (2400 で約 42 %)。

### 試作 (記録 `case/71.m55_combustion_wt/_runs/lcscan_e6/`)

- 変更は `prototype.diff` (51 行)。x_E を挟む 2 節点 (下流 m_hi・上流 m_hi + 1) の C⁺ と壁の交点を求め、x_E の位置の重みで線形補間する。網・壁は変えない。
- `INTERP=0` (今の最寄り) は forge-7d の基準 x_F 76.24832719272557 をビット再現した。
- `INTERP=1`・n_axis_inv 2400・L_c 37.17〜37.26 の 181 点:

| | 局所の傾き | 直線の当てはめからの残差の最大 |
| --- | --- | --- |
| 今 (最寄り、`../lcscan_forge7d/fine2400.json`) | 中央値 0.2920 (段で 105.7) | 2.84e−2 r_t |
| 補間 (`interp2400.json`) | 0.9999〜1.0001 | 1.34e−6 r_t |

- 4800 の掃引・wt_card の通し (V9 (b) の後半・(c)〜(e)) はまだ。

### #6d (今の自動調整)

- M_knot の選択の中で長さ合わせが未収束だと、`n_axis_inv` 2400 → 4800、それでも駄目なら内部の許容差を 0.05 r_t に緩めて粗い掃引から選び直す
  (`design/forge_design/intake/wt_card.py` の `mk_stall`・`design/forge_design/evaluate/mknot_select.py`)。
- 穴: 最初の r_t で緩めた許容差が、解き直した最終の r_t まで持ち越される。カード B の寸法 FAIL (拡大部の長さ −1.366 mm > 公差 1 mm) はこれが原因
  (最終の r_t なら公差のまま 2400 で 3 回で収束)。
- 4800 の網は 1 評価に約 4.5 GB 使い、M_knot に数分掛かる。

## 期待値と出典

- ユーザの了承 (2026-10-16): 直し方で「全部の壁がビット単位で変わります」に「ぜんぜんいいよ」。設計を早く進めるのが第一。
- 新しい標準は既定にする (旧既定を残さない; ユーザの方針)。

## 仮説

- 終端 C⁺ (E から出る) は、隣の 2 節点から出る網の C⁺ の間にある。x_F を x_E で線形補間すれば、最寄り (誤差 O(Δx)) より良い近似 (O(Δx²) の見込み) で、
  x_F(L_c) が連続になり、from_length が 1〜2 回で収束する。
- 補間で x_F は同じ L_c で最大 (1 − 0.29)·Δx/2 ≈ 0.019 r_t (2400) 動く。

## 問い

1. 案 b (補間) は案 a (網の節点を x_E に乗せる) より良いか。物理的・数値的に見落としはないか
   (例: 一様域の始まりを決める C⁺ が網の折れ線であること、x_E の上流側の C⁺ がまだ一様でない域を通ること、壁の切り方 `_cut_wall_at_exit` との整合)。
2. 既定にしてよいか (キーを作らない)。影響する他の経路 (`terminal_exit` の streamline の経路は x_E から直接追うので対象外のはず) の確認で足りないものはあるか。
3. #6d の自動調整は、#6e の後に保険として残すか、外すか。残すなら穴 (r_t をまたぐ許容差の持ち越し) を直す (r_t ごとに元の許容差から選び直す) でよいか。
4. V9 (plan §6) の事前登録で足りるか。閾値 (傾き 0.98〜1.02、残差 ≤ 1e−4 r_t、4 枚とも自動調整なし・寸法 PASS) に問題はないか。

## 読んでよいファイル

- `plans/active/tooling-nozzle-card-local-mesh.md` (§4.5・§5.1 #6d・#6e・§6 V9)
- `design/forge_design/geometry/moc_inverse.py` (`_design_cplus`・`_axis_grid`・`cplus_lines`・`inverse_design`・`terminal_exit`)
- `design/forge_design/evaluate/runner_axismach.py` (480〜700 行: 軸則・from_length)
- `design/forge_design/intake/wt_card.py` (`mk_stall` の周り)、`design/forge_design/evaluate/mknot_select.py`
- `case/71.m55_combustion_wt/_runs/lcscan_forge7d/` (README・scan.py・analyze.py・grid_check.py・JSON)
- `case/71.m55_combustion_wt/_runs/lcscan_e6/` (README・prototype.diff・scan_interp.py・cmp.py・interp2400.json)
- `procedures/design-chain-quickstart.md` §2 (試走の表)
