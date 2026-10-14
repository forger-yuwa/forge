# 諮問ブリーフ: 拡大部長 (設計スロート→物理出口) を実寸の公差で合わせる方法 (plan tooling-nozzle-wt-deliverables §4.1・§4.2a)

AGENTS.md「モデル分担とエスカレーション」条件 1 (§4 の方針を変える) と 4 (承認済みの手順に無い修正の前)。判断役のスイッチが codex。

## 背景

ユーザは風洞ノズルを実寸で指定する: 出口長さ = **設計スロート (x=0) → 物理出口**の長さ (2026-10-12 ユーザ回答)。要求カードの生成器
(`design/forge_design/intake/wt_card.py`) が `Lc_mode: from_length`・`dv.L_total = L_div/r_t` に写し、`geometry.L_total_tol` を
**カードの長さ公差 / 2 / r_t** にしている (codex 諮問 2026-10-13 `notes/reviews/2026-10-13-wt-deliverables-diagnose.md` の採用:
既定 0.05 r_t は r_t 25 mm で 1.25 mm で寸法を保証しない)。最終の合否は、配布する STEP とメッシュから測り直した寸法をカードの公差で
判定する出荷ゲート (`design/forge_design/export/dimension_gate.py`) が決める。

## 観測事実 (2026-10-13、implementer の測定)

1. `from_length` の解法 (`design/forge_design/evaluate/runner_axismach.py:562-620`) は「ステップ = −残差」の固定点反復で、**最初の非改善で打ち切る** (`:608`)。
   届かなければ `:611-616` で例外。
2. case/45 の生産の実寸 (r_t ≈ 0.0767 m、長さ公差 1 mm → L_total_tol 0.00652 r_t) で、同じ形状の最良の残差 |x_F − L| [r_t] は r_t と n_axis_inv で
   単調でなくばらつく:

   | r_t [m] | n_axis_inv | 最良の \|x_F − L\| [r_t] | 所要 |
   | --- | --- | --- | --- |
   | 0.0767718 | 2400 | 0.015 | |
   | 0.0766654 | 2400 | 0.00728 | |
   | 0.0766700 | 2400 | 0.00141 | |
   | 0.0766700 | 4800 | 0.0115 | |
   | 0.0766700 | 9600 | 0.00058 | 492 s |

   生成器は plan §4.1 1. どおり n_axis_inv を 1 回 2 倍にし、届かなければ例外にしている (`wt_card.py:942` 付近)。この場合 case/45 の実寸では YAML ができない。
3. **from_length が合わせるのは MOC の x_F だが、設計壁 (= メッシュ・STEP) の終点 x_e は壁の点列を x_F の手前で切った点** (`design/forge_design/geometry/moc_inverse.py:903-904`:
   `i_cut = searchsorted(wall_x, x_F); wall = wall[:i_cut]`)。設計スロート径の M5 デモで x_F − x_e ≈ 0.85 mm (壁の点の間隔 0.034 r_t, r_t 25 mm)。
   from_length の残差 +0.47 mm と合わせて x_e は −0.37 mm。最悪は L_total_tol + 1 間隔で、長さ公差 1 mm を超えうる。
4. 物理スロート径の M5 デモでは x_F = x_e だったが、これは x_F が MOC の場の終点に張り付いていたため (雛形 x_end_margin 2.3、case/45 由来)。
   x_end_margin 3.0 にすると x_F は +0.18 mm 動き、x_e は x_F より 1.02 mm 短くなった (張り付きの検出と自動の引き上げは生成器で実装中)。
5. M5 デモ (物理スロート径、r_t ≈ 0.02497 m、公差 1 mm) では、測った拡大部長の差は +0.066 mm で合格した (張り付きの状態での値)。

## 制約

- 本 plan (`plans/active/tooling-nozzle-wt-deliverables.md`) のスコープは「設計チェーン本体 (MOC・δ\* ループ) は変えない」。チェーンを変えるなら別 plan
  (候補: `plans/active/tooling-design-problem-campaign-recipe.md` §5.1 #15「寸法固定 (全長・入口径・出口径を実寸で固定、長さの起点を problem キーで定義)」)。
- 生産の壁 (case/45) は `Lc_mode: explicit` で、この問題の影響を受けない。

## 問い

**実寸の拡大部長をカードの公差 (例 1 mm) で満たすために、今回の plan で取るべき方針を 1 つ推奨してほしい。** 候補 (縛られなくてよい):
(a) 生成器は L_total_tol を「到達できる値」(例: 既定 0.05 r_t か、公差に合わせて試し、届かなければ緩めて記録) にして YAML は必ず作り、合否は出荷ゲートだけが決める
(ゲート FAIL ならユーザに出して公差か形状パラメータを選ばせる)。
(b) 今回の plan の中で from_length の解法 (打ち切り規則・合わせる量を x_F から x_e に・壁の点列の終点を x_F に揃える) を直す (チェーン本体の変更になる)。
(c) 別 plan に切り出し、今回は (a) で出荷ゲートの FAIL を正直に出す。
あわせて、観測 3 (x_F と x_e の差) を「壁の点列の終点を x_F に置く」で消すことが設計 (MOC の出口の定義) として正しいか、ゲートで何を「物理出口」と測るべきかを示してほしい。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.msh.h5`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,140p;250,330p' plans/active/tooling-nozzle-wt-deliverables.md`
- `sed -n '455,625p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '860,930p' design/forge_design/geometry/moc_inverse.py`
- `sed -n '1,30p;380,440p;900,990p' design/forge_design/intake/wt_card.py`
- `sed -n '1,40p;60,110p' design/forge_design/export/dimension_gate.py`
- `grep -n "terminal_exit" -A40 design/forge_design/geometry/moc_inverse.py | head -80`
