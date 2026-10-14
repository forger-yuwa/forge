# 諮問ブリーフ: 上流の第一層を薄くするには流れ方向の細分も要る — 案 X08 で登録してよいか (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1・4: 方針の変更、plan に無い修正の前)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #36r (解析の結果を追記済み)、§6「#36o・#36p の解釈」
- 前の諮問: `notes/reviews/2026-10-15-h370-wallres-timeband-result-diagnose.md` (直すのは上流の表の値、既存の H 系列の局所 τ_t・y₁⁺ の包絡から表を設計)
- commit: 31e5795e (台本 `case/45.isobutane_m6_d155/h370_wall_envelope.py`)、7e3e26c9 (記録 `_band_ab/half_inlet/h370_wall_envelope.json`・`h370_wall_table_design.json`)、b60c650b (plan)

## 観測事実

- 包絡 (run_0951、112,500〜152,500 の 17 出力): 時間でほぼ動かない (y₁⁺ の相対の幅 ≤ 2.7e−5、1 超えの面積は全出力 2.932 %)。最大 x −10.53 r_t (y₁ 1.152 µm、τ_t 50.4 Pa、y₁⁺ 2.766)。入口の節点だけ τ_t 146.8 Pa (x −12 で 63 Pa)。
- 対照: 元の表で作り直した格子は run_0951 と生成器の sha・変換後の座標の sha・品質の出力 (AR 4574.0、skew 0.320、PASS) が一致。
- 表の上流 (x < −3.5) の値だけを y₁⁺ 目標 0.8 / 0.9 から薄くする案 (T08・T09、点の置き方 2 通り): **AR 最大 14,433 / 12,829 で `check_mesh_quality --ar-max 5000` FAIL** (5000 超え 813〜1,023 セル、x −12.5〜−6.19)。流れ方向の間隔は変わらないので AR = dx/d₁ が約 3 倍になる。表だけを AR 5000 の上限まで薄くしても 1 超えの面積 2.164 % (> 2 %)。
- 案 X08 (表の T08 + 上流 x ≤ −5.5 の `x_density_table` を最大 3.19 倍、ni 4719 → 4909、節点 755,040 → 785,440): 品質 PASS (AR 4,572.3、AR > 1000 のセルの skew 最大 0.00067)、見積もり y₁⁺ 最大 0.8 (τ_t を包絡の最大に固定した一次の見積もり)、x ≥ 0 の列は元と 1.1e−12 m 以内 (ビット同一ではない)。壁の直交性: 第一内部節点は壁法線上 (角度 約 8e−6 rad)。
- その他: 今のツリーの prepare_ns は節点数の上限 30 万で止まる (今の H の YAML 755,040 も。run_0951 は上限の無い版で準備)。`max_nodes` の上書きにはユーザの了承が要る。run_0951 の格子は出口面が平面でない (x の幅 0.879 mm) ので `mesh_ar_relax_check` は対照も `ok=False` (prepare_ns はこの検査を呼ばない)。

## 期待値と出典

- 壁解像: 局所 y₁⁺ ≤ 1、超過の面積 ≤ 2 % (壁ごと)。メッシュ品質: AR ≤ 5000 (壁法線の構造層の緩和)、skew ≤ 0.9。
- 前の諮問: 「`wall_first_up_*` の追加はしない」「表の上流の値を直す」。

## 呼び出し側の案

1. X08 (流れ方向の細分を上流だけに足す) で新しい NS 1 本を事前登録する。差は上流 (x ≤ −5.5) の流れ方向の間隔と第一層だけ、下流の列は元と 1e−12 m で一致。
2. 節点数の上限: H の格子 (今の 755,040 も X08 の 785,440 も) は上限 30 万を超える。`max_nodes` の上書きはユーザに了承を求める。
3. 入口の節点の τ_t の飛び (146.8 Pa) は前縁に似た格子依存の値。X08 の倍率の最大と見積もりの最大は入口で決まるが、入口の 1 節点のために全体を細かくするのは過剰かもしれない。

## 問い

1. 「表だけ」から「表 + 上流の流れ方向の細分」へ方針を変えてよいか。他に、より小さい変更 (例: 第一層の目標を入口の 1〜2 節点だけ緩める、入口の配管を短くして上流のせん断の大きい区間を減らす、など) で 2 % 以内に入る道があるか。
2. X08 と X09 (目標 0.9、節点 778,720) のどちらで登録するか。
3. 新しい NS の合格条件 (前の諮問の案: 全保存量の `check_convergence --segment --drop 3`、NaN/Inf なし、停止の step と保持の窓の壁の解像、局所 τ_t・y₁⁺ の系列の準定常、壁ごとの超過の面積 2 % 以下) に足すものはあるか。下流の列が 1e−12 m で一致するが restart はどうするか (cross-mesh の `interp_field` を使うか)。
4. 出口面の非平面 (0.879 mm) は今回の判定に入れるべきか。

## 読んでよいファイル

- plan §5.1 #36q〜#36s、§6 の V-h370 系
- `case/45.isobutane_m6_d155/h370_wall_envelope.py`、`_band_ab/half_inlet/h370_wall_envelope.json`・`h370_wall_table_design.json`
- `design/forge_design/meshing/mesh2d.py` (`_x_stations`・`wall_first_frac_table`・`x_density_table`)、`design/forge_design/evaluate/runner_axismach.py` (`check_node_budget`・`MAX_NODES_2D`)、`design/forge_design/evaluate/runner.py` (`mesh_ar_relax_check`)
- `notes/reviews/2026-10-15-h370-wallres-timeband-result-diagnose.md`
