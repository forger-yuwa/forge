# 諮問: 単調壁 Euler A/B の投入前の判断 4 件 (腕 B の IC 写像、近零量の準定常判定、評価器の前提条件規則、δ_r の環境依存)

日付 2026-10-06。諮問先 codex (diagnose)。エスカレーション条件 4 (plan §6 の手順から外れる修正の前)。
plan: `plans/active/tooling-nozzle-throat-monotone-r2.md` (§6 E1・E4・N、§6.1)。作業ツリー `/home/sano/work/forge-integ-1005`、commit d4364794 (実装・形状ゲート PASS)。

## 観測事実 (実装担当の報告、主セッションで再確認したのはテスト 5 本の FAIL 0 と形状ゲート JSON の VERDICT PASS)

1. **腕 B の IC 写像** (`case/45.isobutane_m6_d155/throat_mono_ab.py prep`、ローカルの dry prep で実測。メッシュは生産 Euler 格子 G1 = ni 2000 × nj 97):
   - 腕 B の変換後メッシュは、腕 A (= run_0114 と座標一致) と節点数・接続が同じ。違いは壁と近傍の節点座標だけで、壁節点の差は最大 0.4992 µm (解析値 0.5009 µm、x 0.0634 r_t)。
   - interp_field (最近傍) で run_0114 の最終場を移すと、194,000 節点中 37 節点 (一致率 0.99981) が j−1 (1 層内側) の値を取った。場所は x 0.032〜0.094 r_t、j = 86〜96 (壁から 11 層)。原因は、スロート近傍の第 1 セル厚 0.34〜0.35 µm が壁の移動 0.5 µm より小さいこと。
   - その結果、密度が局所で最大 23 % ずれる (max|Δρ| / max|ρ| = 0.11、roUx 0.22)。原始量からの組み直しのため、roUx は 130,588 節点で丸め差が出る。
   - plan §6 E1 は「B は interp_field (規則どおり。両腕への interp_field や座標検査の緩和はしない)」「IC 写像を記録し、IC の影響は限界として扱う」(前回諮問 ④ の採用)。
2. **近零量の準定常判定**: run_0114 のデータで評価器を通すと、`exit_M_dev` = |出口コア M − 6| ≈ 3e-6 は check_quasisteady の相対基準で DRIFTING になる。E4 では「末尾 5 枚 STEADY でない量は判定に入れず保留」なので、この量はほぼ常に保留に落ちる見込み。Δq は 1.8e-4。
3. **評価器の追加規則**: 実装担当は、壁の証拠照合の不一致・腕 B の mono_r2 欠落・出口の既定帯の不一致があれば、総合判定を「保留 (前提不成立)」にする規則を足した。E3 には書いていない。
4. **δ_r の環境依存**: 生産経路で δ_r を再計算すると、run_0117 の `delta_r_initial.csv` と 1.38e-5 r_t (venv、numpy ≥ 2) / 8.7e-6 r_t (system numpy 1.26) 違う。差は `integral_bl` の生出力 (dstar_n) の時点で出ている。§6 N の NS は run_0117 の δ_r を約 1 µm の精度でしか再現しない (設計壁の変更 0.5 µm と同程度)。形状ゲート S6 の「単調版 vs 現行の物理壁」最大 Δr 7.78e-6 r_t (x 90.27) もこの再計算差による。

## 選択肢と当方の推奨

1. IC 写像:
   - (i) plan どおり interp_field で、記録を限界とする。
   - (ii) **推奨**: 「節点数・接続が同一で、座標の移動が ≤ 1 µm」を検査したうえで、保存量を番号でコピーする (case/60 `tools/stack_init.py` と同じ考え方)。腕 A の restart_field と同じ「各節点が自分の層の値を持つ」IC になり、腕間の IC の扱いがそろう。実装は case 内のスクリプト (restart_field 本体は変えない) とし、移動量の分布 (最大・位置・節点数) を記録する。前回諮問で「座標検査の緩和は推奨しない」とされたのは 1e-6 r_t 案 (0.08 µm で 0.5 µm を許容できない) についてで、接続の同一性検査は含んでいなかった。
   - (iii) 両腕とも等エントロピー IC から段階起動 (対称だが過渡が長い)。
2. 近零量: |出口 M − 6| のように末尾平均が Δq より十分小さい量は、相対基準ではなく絶対基準で判定する。「末尾 5 枚の幅と直前 5 枚との平均差が、どちらも ≤ Δq/10」なら定常扱いとする。他の量は従来どおり check_quasisteady。
3. 追加規則: 採用する (保守側)。E3 に明記する。
4. δ_r: N は AWS の同じ環境 (numpy ≥ 2 の venv) で回し、run_0117 の δ_r との差を記録する。N の比較は「設計壁の変更 + δ_r の再計算差」の合算と明記する。もしくは N だけは run_0117 の δ_r 表を固定で使い (`delta_r_csv`)、設計壁の変更だけを分離する。

## 問い

1. IC 写像は (i)〜(iii) のどれにするか。(ii) の検査条件 (接続同一・移動 ≤ 1 µm) は十分か。AGENTS.md の「同一メッシュは restart_field、別メッシュは interp_field」とどう整合させるか。
2. 近零量の絶対基準 (Δq/10) は妥当か。
3. 評価器の追加規則を採用してよいか。
4. N の δ_r は、生産経路で再計算する (= 生産そのもの) のと、run_0117 の表に固定する (= 設計壁の変更だけを分離) のと、どちらを主にするか。両方回すべきか。

## 読んでよいもの

- 上記 plan
- `case/45.isobutane_m6_d155/{throat_mono_ab.py, eval_wallfit_euler.py, throat_mono_judge.py, run_throat_mono_ab.sh, throat_mono_shape_gate.py}`
- `_band_ab/throat_mono_shape_gate.json`
- `solver_density_cuda/tools/{restart_field.py, interp_field.py}`
- `design/forge_design/evaluate/runner_axismach.py` (`integral_delta_r`、`prepare_ns`)
- 前回の諮問記録 `notes/reviews/2026-10-06-throat-monotone-r2-diagnose.md`

編集は禁止。
