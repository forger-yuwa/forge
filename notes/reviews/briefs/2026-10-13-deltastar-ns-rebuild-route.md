# 諮問ブリーフ: NS 後の δ_E 全分布でカードから作り直す経路 (plan tooling-nozzle-deltastar-band-root-zone §4.2・§4.3 の改定)

AGENTS.md 条件 1 (plan §4 の方針を変える)。判断役のスイッチが codex。前回の諮問 (`notes/reviews/2026-10-13-deltastar-band-root-zone-diagnose.md`) の推奨した寸法の A/B を回した結果、
経路を変えたので、その改定を諮る。

## 観測事実

- **寸法の A/B** (`case/70.m5_air_wt_demo/sizing_ab.py`・`_sizing_ab/sizing_ab.json`、forge 0 step): 同じ `problem_m5_air_c2b.yaml`・物理スロート半径 25 mm・同じ壁の構築。

| 腕 | 補正 | r_t [m] | ×run_0010 | 物理スロート x / r [r_t] | 長さの誤差 直管 / 縮流部 / 拡大部 [mm] | (r_t/S_Euler run_0003)² |
| --- | --- | --- | --- | --- | --- | --- |
| A | run_0010 の壁 (CONTUR k_f 0.8246) | 0.024973561 | 1.000000 | −0.00175 / 1.001058 | 0 / 0 / −0.485 | 1.000367 |
| B | run_0010 の修正後 (zone) の δ_E 全分布、終点を x_e に揃えた写し (差 1.54e-6 ≤ ε_f32 8.6e-6) | 0.025075725 | 1.004091 | −0.0011 / 0.99698 | +0.818 / +1.227 / **+4.422** | 1.008568 |

  事前の判定「A が公差内・B が公差外 → r_t 単位の設計を据え置く経路を棄却」に当たった。前回の codex の近似予測 (×1.00403・+4.8 mm) と一致。
- 前回の諮問の Major 2 件 (x = 10 で δ_E −9.3 %・寸法の変化の過小見積もり) は採用済み (plan §6.1)。抽出の修正 (§4.1) は未実装 (測定は `ab_root_zone.py` の差し替え)。
- 壁温分布の配管 (別 plan) で、区間で正規化した座標 s の写し (`design/forge_design/wall_temperature.py` の `anchors_from_design`・`to_x_m`) を実装・検証済み。

## 改定した経路 (plan §4.2・§4.3 を読むこと)

1. Euler (前の設計 = run_0010 と同じ c2b) を回す → 2. run_0010 を修正後の抽出に通す (V3) → 3. 要求カード (実寸の長さ) の r_t の固定点で、内側の寸法を
   δ_E を s で写した補正で解き、最終の表を YAML の横に書く → 4. 新しい設計の Euler → 5. 新しい設計の NS (IC = run_0010 の cross-mesh、none・cfl 1・60000 step)
   → 6. V5 (出口コア M 4.9975〜5.0025・各 VERDICT)、V5b (1 回の作り直しで固定点に近いかの記録)。
- #18 は「前の設計の終点に float32 の 2 回の丸めの範囲で揃える」→ s で写す、で扱う。

## 問い

**改定した経路 (§4.2・§4.3・§6 V4・V5・V5b) のうち、実装後に「出口コア M が 4.9975〜5.0025 に入らない、または入っても理由を取り違える」形で失敗する可能性が最も高い点を 1 つ挙げ、
防ぐ変更を 1 つ推奨してほしい。** 候補 (縛られなくてよい): (a) 補正 1 回 (固定点の反復をしない) で足りるか — 壁の変化で境界層が変わる、
(b) 補正の元の NS (run_0010) の壁が CONTUR の物理壁で、新しい壁とスロートの δ が逆符号 (+0.001 → −0.003 r_t)、(c) s の写しと Re^−0.2 の換算、
(d) 新しい設計の Euler を参照にする (前の設計の Euler でなく) ことと流量の比の解釈、(e) IC を run_0010 から cross-mesh で移すこと。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.cgns`, `*.stl`, `plans/README.md`, `delta_r_equiv_diag.json` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,200p' plans/active/tooling-nozzle-deltastar-band-root-zone.md`
- `cat case/70.m5_air_wt_demo/_sizing_ab/sizing_ab.json`
- `sed -n '288,357p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '1253,1300p;1360,1420p' design/forge_design/intake/wt_card.py`
- `sed -n '100,165p' design/forge_design/wall_temperature.py`
- `sed -n '1,40p' notes/reviews/2026-10-13-deltastar-band-root-zone-diagnose.md`
