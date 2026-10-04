forge (CFD ソルバ, 有限体積, node 中心 median-dual) の ⑤ SERN ノズル設計チェーンで、
3D メッシュの作り方を **(a) 現行の自作メッシャを直し続ける** と **(b) 3D CAD を起こして Salome でメッシュを切る**
のどちらにするかを決めたい。推奨を 1 つに絞って理由を述べてほしい。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `plans/active/tooling-nozzle-sern-3d.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの
- `sed -n '1,60p' design/forge_design/meshing/mesh_sern3d.py`   (モジュール docstring + パラメータ定義)
- `sed -n '130,180p' design/forge_design/meshing/mesh_sern3d.py` (station 列と z 分布の構築)
- `sed -n '230,300p' design/forge_design/meshing/mesh_sern3d.py` (dup1/dup2, coords, node(), hex 生成)
- `sed -n '300,350p' design/forge_design/meshing/mesh_sern3d.py` (境界 quad の生成)
- `sed -n '425,450p' design/forge_design/meshing/mesh_sern3d.py` (_add_ext_top3d)
- `sed -n '495,520p' design/forge_design/meshing/mesh_sern3d.py` (_add_vehicle_side3d)
- `sed -n '110,175p' design/forge_design/evaluate/runner_sern3d.py` (prepare: メッシュ生成→変換→初期場→帳簿)
- `sed -n '280,300p' design/forge_design/evaluate/runner_sern3d.py` (forces3d の入口項)
- `notes/reviews/2026-09-21-tooling-nozzle-sern-3d-plan.md` (今日の設計レビュー全文, 10 KB)
- `sed -n '290,330p' procedures/calculation-workflow.md` (Fluent CFF HDF5 取り込みの節)
- grep は自由 (例: forge が読めるメッシュ形式を探す)

## 背景 — 現行 (a) の実装
`mesh_sern3d.py` (596 行) は **Python 手書きの構造化ヘキサメッシャ**。x-station × y-band × z の
テンソル積で節点を張り、スリット壁 (カウル板・側壁) は **同一座標の双子ノード** (`dup1`/`dup2`) で表現する。
出力は msh4.1 → `convertGmshToForge` で forge の HDF5 に変換。
これで 2D/3D とも生産 run を回してきた実績がある (壁解像 y⁺≈1、力係数の格子収束まで到達済み)。

## 指標の定義 (今日の数値)
- **メッシュ品質**: `check_mesh_quality.py` の per-face equiangle skew と edge aspect ratio。
  合格は AR ≤ 5000 (壁法線構造層に限る緩和) かつ skew ≤ 0.90。>0.1 % 超過で FAIL、<0.1 % で SOFT-PASS。
- **格子収束**: 形状・BC を固定し壁第一層厚を 4 倍ずつ細かくした 3〜4 点で、力係数の隣接差を許容と比べる。
  許容 |ΔC_T| ≤ 0.002, |ΔC_L| ≤ 0.002, |ΔC_M| ≤ 0.05。

## 今日ここまでで到達したこと (現行 (a) のまま)
| 項目 | 結果 |
| --- | --- |
| メッシュ品質 | SOFT-PASS (AR>5000 が 0.02 %、skew max 0.701、>0.90 が 0) |
| 壁解像 y⁺ | ramp 0.549 / cowl_in 1.234 / cowl_out 0.185 (側壁は 188 で未解像) |
| 格子収束 (g3→g4) | C_T −0.00013 / C_T_with_shear −0.00024 / C_L +0.00145 / C_M −0.04173 → **4 係数とも許容内** |
| 摩擦 | C_T_friction が g2→g3 +0.00407 → g3→g4 −0.00011 で **収束** |
| 規模 | 3,222,028 節点 / 3,122,688 セル |

## 潰した候補 (既に棄却済み — 再提案しないでほしい)
- リミッタのチャタリング説: `venkatK` を 100 倍振って残差プラトーは 4 % しか動かない → 棄却
- 単純な dt 比例の限界周期: CFL 半分で 0.74 倍 (比例なら 0.5) → 棄却
- 丸め誤差のノイズ床: 同一入力 5 反復の全ペア最大に対し場の変化は 3〜62 倍 → 棄却
- 断熱壁が原因説: 生産の等温壁 1000 K にしても残差プラトーは消えず → 棄却
- 残差プラトー自体: **ユーザ判断で「変動が小さいだけ」として追跡打ち切り** (力は 7 桁で定常、派生量 ALL STEADY)
- 外側 z の細分で AR を下げる案: `nz_out` 15→40 で max AR は 29712 のまま不変 → 棄却 (効いたのは x ブレンド)

## 残っている欠陥と、その修正設計が受けた指摘
残る欠陥は 1 つ: **厚さ 0 の面が 2 つ交わる線**で ρ が床に張り付く。
- 側壁は厚さ 0 のスリット (`sidewall_in`/`out` の z 差が min=中央値=max=**0.000e+00**、双子 5893 対)
- カウル板は側壁ちょうどで厚さ 0 に潰れる (z=0.0980 で 5.000e-04、**z=0.1000 で 0.000e+00**)
- その交線の双子は **内側 ρ 2.498e-01 / P 85018.6 Pa、外側 ρ 1.000e-04 (床) / P 28.7 Pa** (外部静圧 2851 Pa の 1/100)
- **格子細分で悪化**: ρ min が 1.95e-03 → 1.61e-03 → 9.51e-04 → **1.00e-04 (床)**、床ノード 0/0/0/**12**

これを「側壁に物理厚み + z 分布の station 依存化 + カウル側端テーパ廃止」で直す設計を書いたところ、
外部レビューで **NO-GO / Major 7** を受けた (全文は上の notes/reviews のファイル)。要旨:
1. `L_sw`(0.8) < `L_cowl`(1.2) なので、側壁が終わった後の**露出カウル側端を閉じる面が無い** (未分類の単独所有面 12)
2. `dup2` は `j > jm` のみなので **側壁下端 (j=jm) に厚みを与えられない**
3. `_add_ext_top3d` / `_add_vehicle_side3d` も `zs[k]` を使うので**接合が格子依存**になる
4. `i_sw = argmin(|xs − L_sw|)` で **`L_sw` が station に無い** (共有開始 0.7735 H vs 厚み則の要求 0.00332 H)
5. 厚み 0 でも `sz` をステップに替えると **1108 節点の座標が動く** (「既定でビット一致」が成立しない)
6. **入口幅が変わるのに `half_W_m = W·H/2` が力の帳簿に固定**で、格子列を取り直しても系統誤差が残る
7. 既存試験はカウル厚 0・`L_sw` 未指定で**新モードを覆っていない**

## 問い (1 つだけ)
**(a) 自作メッシャを直し続ける / (b) 3D CAD + Salome に移る のどちらを推すか、1 つに絞って答えてほしい。**
判断にあたって次を具体的に評価してほしい:
- 上の Major 7 件のうち、(b) に移ると**構造的に消えるもの**と**形を変えて残るもの**の切り分け
  (特に (1)(2)(4) のような「スリット壁と共有ノードの設計」が、非構造/CAD ベースで本当に消えるのか)
- **forge 側の入力制約**: node 中心 median-dual・スリット壁 (同一座標の双子ノード) を要求する境界・
  壁距離・physID による BC 対応。grep で確かめてよい。既存の外部メッシュ経路は
  gmsh msh4.1 (`convertGmshToForge`) と Fluent CFF HDF5 (`fluent_h5_to_forge.py`)。
  **Salome の出力 (MED/UNV) から forge へ入れる経路は現状ない**。その新設コストも見積もりに入れてほしい
- **壁解像 y⁺≈1 の境界層層数**を (b) で作れるか (現行は構造化なので自然に積めている)
- **設計チェーン特有の要件**: dv (M_c, f, θ_r0, θ_c0, L_cowl) を振って**形状が毎回変わる**ので、
  メッシュ生成は**完全自動・無人**である必要がある。MOO で数十〜数百ケース回す
- 移行するなら**何を捨てることになるか** (今日到達した格子収束・壁解像・起動レシピの再取得コスト)

(a)(b) 以外に筋の良い第三案があるなら挙げてよいが、**推奨は 1 つに絞ること**。
