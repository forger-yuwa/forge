# 諮問ブリーフ: 壁温分布の配管の result 段の解釈 (plan tooling-nozzle-wall-temperature-profile §6)

AGENTS.md 条件 7 (result 段の解釈を確定する前、`--stage result` に回す前)。判断役のスイッチが codex。

## 観測事実 (commit 4ea3626d、feature/nozzle-deliverables)

| 検証 | 結果 | 出典 |
| --- | --- | --- |
| V1 単体 | `run_wall_temperature_tests.py` 100 ok・`run_wt_card_profile_tests.py` 87 ok | design/tests/ |
| V1b 実効の r_t | 積分法の `res["Tw"]` と原表を `scale × res["x"]` で評価した値の最大差 5.7e-14 K (scale = r_t と 1.2 r_t、s と x_m)。旧半径で換算する誤実装なら 22.4 K。x_m の表で 1.2 r_t だけ範囲外 → 例外 | run_wall_temperature_tests.py |
| V1c 見積もり | 局所の T_w と丸め前の第 1 層厚が期待値と相対差 0 (8 ケース)。誤実装 (a) 代表温度・(b) 無次元の x・(c) カードの長さで区切り を検出 (y1 のずれ 5.6e-5〜4.5e-1、(c) の一部は例外) | run_wt_card_profile_tests.py |
| V1d 正本 | CSV 1 値の書き換え・別 run の JSON・スロート同温で分布だけ違う JSON の差し替えで、`run_staged_ns` は forge の起動前に、成果物は出力を作る前に止まる | 両試験 |
| V3 成果物 | Fluent profile 1799 点 = 壁面数、最近傍の割り当て値と run の実寸の表の相対差 0 (元の s の表を独立に写した値とは 6.4e-11)。壁ノードで書く誤実装なら 585 面が 0.1 K 超 | run_deliverables_walltemp_tests.py 36 ok |
| V4 後方互換 | 一様 (case/70 c2b) と断熱 (case/45 生産) の prepare の生成物が 138e83b2 とバイト同一。成果物は 5a6dd76f と同一 (xlsx・CGNS の時刻を除く)。カードは時刻の値を除き同一 | 各試験 |
| **V2 デモ** | `case/70.m5_air_wt_demo/run_0011_ns_walltemp_demo`: s の表 (−2/−1/−0.3/0/0.2/0.5/1 → 320/330/420/500/450/380/330 K)。壁 1800 節点の場の `VALUE/T` (res_60000.h5、`BCONDS/3/iCells` の節点、`MESH/COORD` の x) と run の JSON の実寸の表の線形補間の**最大差 1.5e-5 K** (許容 0.1 K)、両端一致。`delta_r_initial.json` の `thermal_bc.Tw_table` の x 終点 48.0615 = 1.20009 m / r_t 0.0249699 m。NaN 0、品質 PASS。収束 (合否と別): `check_convergence --segment` main NOT CONVERGED (ρ 1.2 桁、ρUy 横ばい)、`check_quasisteady machmax,pmax` ALL STEADY、壁解像 PASS。成果物 PASS (`fluent/wall_temperature.prof`) | case/70 README |

- 実装の途中で plan の穴が 1 つ出た: 見積もりの s の区切りの点をカードの長さで取ると、公差の内で設計壁の終点がカードより長い (+0.089 mm) ときに表が壁を覆わず止まる → 設計壁から取るように §4.3 を訂正 (plan §9)。

## 当方の解釈 (確定前)

- §6 V1〜V4 は合格。「分布を与えたつもりが NS と積分法・成果物で別の壁温になる」形の失敗は、実効の r_t (V1b)・正本の照合 (V1d)・面中心 (V3)・場の壁ノード (V2) で塞いだ。
- V2 は node の壁が温度をピンするので配管の検査であって、流れの妥当性 (出口 M 4.9855 など) は対象外。

## 問い

**上の結果のうち、「壁温分布が単一ソースで NS・積分法・成果物に入る」という完了の主張を支えていない (あるいは過大に言っている) 点を 1 つ挙げ、accepted に移す前に足すべき確認を 1 つ推奨してほしい。**
候補 (縛られなくてよい): (a) V2 は 1 run・1 つの表で、x_m の表や長さを変えた表を NS で通していない、(b) V1b は積分法の内部の `Tw` を見ているが、その δ_r が物理壁に反映されたこと (壁の形の差) は見ていない、
(c) restart・`rerun_conditions.py`・SERN は範囲外 (止まる側)、(d) Fluent での読み込みが未確認。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.cgns`, `*.stl`, `*.xlsx`, `*.prof`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,170p' plans/active/tooling-nozzle-wall-temperature-profile.md`
- `sed -n '1,313p' design/forge_design/wall_temperature.py`
- `sed -n '1160,1180p;1330,1390p;1526,1560p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '220,260p;310,320p' design/forge_design/export/deliverables.py`
- `grep -n "def \|check(" design/tests/run_wall_temperature_tests.py | head -80`
- `grep -n "run_0011" case/70.m5_air_wt_demo/README.md`
