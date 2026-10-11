# 諮問ブリーフ: ③ベル (D) の MOO — warm start の発散と第一セル厚の余裕

- 日付: 2026-10-11
- 依頼元: 主セッション (Opus)。エスカレーション条件 4 (plan §4.1・§5.1 #8a の手順に無い変更と、キャンペーンの回し直しの前)。このセッションの諮問は diagnostician (Fable) — ユーザ指示
- 作業ツリー: `/home/sano/work/forge-thruster` (HEAD `aed11dc5`)。run は AWS B の `~/forge-thruster/case/68.bell_two_mode/` (手元に結果ファイルは無い。数値は下に写した)
- 対象 plan: `plans/active/tooling-nozzle-bell-two-mode.md` §4.1 (「最初の TP の seed は段階起動、形状が変わる warm start は `interp_field.py` の cross-mesh restart」)・§5.1 #8a

## 1. 前提 (決定済み、覆さない)

- ③は (D) 直接モードと (M) 特性曲線モードの両方を提供する。(D) の MOO は `problem_bell_tp_d_moo.yaml` (dv = θ_max・θ_e・L/r_t、目的 −η_exit と L/r_t、LHS 40 + EHVI 10 回 × 2)。1 評価 = 段階起動 soft 2000 → mid 2000 → 本段 60000 (陰解法 cfl 4・2 次・SST 低 Re)。ゲートは §4.1 の組み替え後 (i)〜(v)。
- 段の長さは warm start の有無で変わらない (どちらも同じ 3 段を回す。warm start は初期場だけを差し替える)。

## 2. 観測事実

キャンペーン `run_0010_moo_d` (B、GPU は単独、1 評価 約 2〜3 分):

| 評価 | x (θ_max°, θ_e°, L/r_t) | 起動 | 結果 |
| --- | --- | --- | --- |
| doe_000 | 24.03, 16.81, 18.20 | 冷間 (段階起動、準 1D の等エントロピー IC) | PASS、η_exit 0.968948、y1+ 最大 0.846 |
| doe_001 | 31.96, 12.39, 16.20 | warm (doe_000 の最終場から `interp_field.py`、最近傍。最近傍距離の最大 0.138 m ≈ 1 r_t) | 本段で発散 (RETRYABLE)。soft・mid は下がった (mid 末 rms_ro 7e-5)。本段の rms_ro は step 2000 で 9.0e-6 → 4000 で 1.6e-5 → 6000 で 2.8e-4 → 10000 で 2.4e-2 と上昇し、高いまま暴れて step 41558 で NaN |
| doe_002 | 37.10, 2.63, 20.09 | warm (同上) | 本段は完走、η_exit 0.965720、収束系のゲートは通過、**壁解像で FAIL**: y1+ 最大 1.096、1 を超えた点 7.3 % (許容 2 %)、位置 x 1.822 m・r 0.923 m (拡大部の下流) |

doe_001 の本段のスナップショット (3000 step ごと) で、拡大部 (x > 0.01 m) の Ux < −1 m/s の節点数: step 0 で 0 → 3000 で 32 (出口と壁の角、x 2.254〜2.268・r 0.987〜0.990) → 6000 で 457 → 9000 で 1566 (x 2.09〜2.27) → 12000 で 7998 (x 0.67〜2.27、r 0.30〜0.97)。step 39000 では軸上 x 1.72 m で T 6000 K (上限)・P 1.17e7 Pa (燃焼室 1e7 を超える)・Ux −1884 m/s。本段の開始時点 (res_0) で出口の角の帯 (出口から 2 cm・壁から 1 cm) は Ux 最大 1568 m/s・T 最大 3552 K (doe_000 の同じ帯は 3147〜3854 m/s・3053〜3267 K)。出口 BC は `outflow` (全量外挿)。

**解釈 (主セッション、仮説)**: θ_max が大きい doe_001 の壁は、同じ x で doe_000 の壁より外にある。最近傍の補間では、その間の帯 (最大 ≈ 1 r_t) に doe_000 の壁の値 (速度 0・回復温度) が入り、出口の壁際が遅く熱い状態で本段が始まった。そこから出口の亜音速の層を通る逆流が育った。

**A/B (変えたのは warm start の有無だけ)**: `run_0011_tp_d_doe001_cold` = doe_001 と同じ problem YAML を冷間の段階起動 (`runner --staged`)。結果: **全ゲート PASS** — 本段の rms_ro は 2.2e-6 (step 6000) → 6.5e-7 (12000) → 3.4e-7 (60000) で単調に下がって横ばい、ω は 6.1e3 → 1.2e2。η_exit 0.970260、η_wall 0.970269、CV 残差 8.8e-6、ṁ/ṁ_1D 0.99073、y1+ 最大 0.853、推力系 4 列と ṁ が STEADY (1.5e-4 / 1e-3)、残差は `NOT CONVERGED (stalled)` (併記)。

**第一セル厚の試し** (forge は回さず格子だけ): doe_002 の形で下流の第一セル厚 `mesh.wall_first_frac` を 2.5e-5 → 1.75e-5 (×0.7、局所半径比) にすると、AR 最大 3770.1 (変わらない — 最大はスロート側のセル、`wall_first_frac_throat` 3.5e-6)、AR p99 1539 → 1549、平均 142 → 150、AR > 1000 の全セルでスキュー最大 0.0020 (AR 緩和の条件 OK)、品質 VERDICT PASS。y1+ は第一セル厚にほぼ比例するので、doe_002 の 1.096 → 約 0.77 の見込み (未測定)。

過去の経緯: CPG の MOO (case/40) でも cross-geometry warm start を陰解法 cfl 4 に直接入れると発散・不始動に落ちた (`opt/driver.py` 冒頭の docstring、run_0019・run_0073)。そのため CPG 経路は soft 段を常用している。TP 経路の冷間の段階起動は run_0008 (TOP)・run_0009 ((D) 1 点)・doe_000・run_0011 の 4 本とも PASS。

## 3. 主セッションの案 (諮りたいこと)

1. **TP の評価は全点を冷間の段階起動にする** (driver の TP 経路で warm start を使わない)。段の長さは同じなので費用は増えない。warm start の方を直す (例: x/L・r/r_w(x) の正規化座標での補間) のは、利得が無いので今はやらない。plan §4.1 の該当文を差し替える。
2. **下流の第一セル厚を 1.75e-5 (×0.7) にする** (スロート側は据え置き)。キャンペーンの問題 YAML の変更。
3. **キャンペーンは新しい run ディレクトリで最初から回し直す** (`run_0012_moo_d`)。run_0010 (doe_000〜002) は記録として残し、混ぜない (格子の規則が変わるため)。

## 4. 問い (推奨を 1 つずつ)

1. 案 1〜3 で進めてよいか。warm start を TP で残す理由 (例: 不始動の偽解の回避は CPG の文脈だった) があれば。
2. 冷間で全点を回すとき、別の失敗の仕方 (例: θ_e が小さい・θ_max が大きい形で冷間の準 1D IC が悪い) を警戒すべきか。警戒するなら、キャンペーン中に何を見るか (例: RETRYABLE が何点続いたら止める)。
3. 第一セル厚は ×0.7 で足りるか (y1+ の余裕をどこまで取るか)。それとも壁解像のゲートの方を見直すべきか (許容 2 % を変えない前提で)。
4. RETRYABLE (本段の発散) の点を MOO はどう扱うべきか。今の driver は FAIL として台帳に入れて先へ進む (再投入しない)。

## 5. 読んでよいもの

- `plans/active/tooling-nozzle-bell-two-mode.md` (§4.1・§5.1 #8a・#3f)
- `design/forge_design/opt/driver.py` (冒頭 1〜60 行の docstring、`_warm_for` 409 行付近、`_evaluate_tp` 529 行付近、`run` 605 行付近)
- `design/forge_design/evaluate/runner.py` (`run_staged_tp`・`_prepare_tp` の warm_from の扱い 660〜680 行付近)
- `case/68.bell_two_mode/problem_bell_tp_d_moo.yaml`、`case/68.bell_two_mode/README.md`
- `procedures/divergence-and-startup.md` §0・§1 (A)
- 禁止: `*.log`, `res_*.h5`, `*.pdf`, `plans/README.md`。ファイルは変更しない。AWS には入らない
