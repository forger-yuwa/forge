# 諮問ブリーフ (2 回目): M5.5 の NS (case/71 run_0002_ns) の未到達 — 証拠をそろえた上での次の一手 (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 2・4)
- 前の諮問: `notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.md` (分類 incomplete: 証拠不足 → 回収して判別 A/B を実施)
- 前のブリーフ: `notes/reviews/briefs/2026-10-15-m55-ns-unsettled.md` (条件・観測の抜粋)
- 回収した記録と判定の原文: `case/71.m55_combustion_wt/_runs/` (索引 `_runs/README.md`、commit c4dfdedb)。plan `plans/active/tooling-nozzle-card-local-mesh.md` §5.1 #3b の行に要約。

## 観測事実 (今回そろえたもの)

- 実効の設定: 問題 YAML (sha256 f1ae7ee1…) から生成器で作り直した S1・S2・本段の solverConfig・bcond・species_meta が run の記録と一致 (文字列・sha256・cfg_fnv × 3・bcond_fnv)。バイナリは FP64 (flow_float・geom_float = double、sha256 3a9e347b…、2fd36b0c、ビルドより新しい源なし)。停止規則は既定の要求 (explicit false)。`thermoFloat 1` は既定で有効。
- `stage_manifest.json` に本段が無い (見張りの save_record が OSError で落ち、`_record(cfg_main,"main")` に届かなかった)。本段を足した写しで `check_convergence --segment` (区間 = 本段 0〜183,998): NOT CONVERGED (still converging)、低下 ρ 1.8・ρE 2.2・運動量 2.4・k 2.9・ω 2.7 桁、全列 falling。各列の最小は 157,795〜172,251 step にあり、直近 2 万 step の log10 の平均はその前の 2 万 step より 0.07〜0.44 桁高い。
- `check_quasisteady --series-csv` (停止規則の系列 1,000〜183,000): DRIFTING (r DRIFTING、θ_r 0.4・0.7 TRANSIENT-UNSETTLED、θ_r 0.95・ṁ_in・ṁ_out・出口コア M は STEADY (既定の drift 5 %))。
- res_183000: 非有限 0、ρ・P・T > 0。静圧 ≤ 全圧は壁の行の 298 点で破れ (相対 ≤ 5.8e−8、U = 0)。**主流〜軸の 50,441 点で P0 > Pt (1 + 1e−6)、最大 +0.24 %**。入口の列の 9 点で T0 > Tt + 1 K (最大 1,253.4 K)。
- 壁解像 (step 183,000、全 no-slip 壁、面積の重み): PASS、y₁⁺ 平均 0.485・p99 0.821・最大 1.004、超過の面積 0.0 % (停止の step ではない)。
- θ_r の測定帯の A/B (forge 0 step、同じ 12 枚の場 100k〜183k、同じ Euler の参照 res_12000、同じ 3 断面、帯だけ変える): D_θ [%/2 万 step] 断面 0.4/0.7/0.95 = 腕 A (res_1000 の固定帯) 0.969/0.912/0.538、腕 B (res_183000 の帯) 2.332/1.960/1.748。登録どおり「帯の選び方だけで未到達を説明できる」を棄却。帯 A は後期の境界層の縁 y* を使用域の全列で覆っていない (帯 B は A の 1.40〜5.52 倍の幅)。
- 停止規則の系列 (前のブリーフの表): 質量の収支 r が 3 万 step で約 0.62 倍 (時定数 約 6 万 step) で単調に減少、ṁ_in が下がり ṁ_out が上がる。
- A の空き 16.2 GB、全場の出力 1 枚 100.6 MB、壁の出力 1.17 MB。

## 問い

1. 続けて伸ばすか、先に何かを切り分けるか。特に主流・軸の P0 > Pt (+0.24 %) は物理的にありえない (等エントロピーの核で全圧は入口を超えない) が、過渡 (1 次元の初期場からの未収束)・数値 (SLAU・`thermoFloat 1`・軸)・後処理 (`total_quantities.py` の h0 と物性) のどれを疑うべきか、判別の最小の手。
2. 伸ばす場合の条件: 出力の間隔 2,000 step・`evaluate.stop_rule: true` の明示 (前の諮問)・上限・帯の扱い (続きの run は帯を取り直す; 元の系列との連結の可否)・ディスクの管理 (1 枚 100 MB、2,000 step ごと)。前の諮問の見込み「質量の条件だけで追加 約 16.8 万 step」をどう扱うか。
3. 本段の manifest の欠落の扱い (写しで足して判定してよいか)。

## 読んでよいファイル

- `case/71.m55_combustion_wt/_runs/` 一式 (README・judge/・theta_band_ab.json・stop_rule.json・residual_history.csv.gz・設定)
- `design/forge_design/evaluate/stop_rule.py`、`runner_axismach.py` (`run_staged_ns`)、`solver_density_cuda/tools/total_quantities.py`
- `plans/active/tooling-nozzle-cooled-wall-ns-settling.md`、`plans/active/tooling-nozzle-card-local-mesh.md`
- 前の諮問の記録
