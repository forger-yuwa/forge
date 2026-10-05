# 諮問: rerun_conditions.py の設計 (plan §4 草案の要判断 4 点と §6 検証)

日付 2026-10-06。diagnostician (ユーザ指定)。エスカレーション条件 1 (plan §4・§6 の新規作成)。
plan: `plans/active/tooling-rerun-conditions.md` (全文が草案)。動機はユーザの「単純に入口条件変えて計算回したい」(形状を固定)。
関連: `procedures/nozzle-design-workflow.md` §3a、`solver_density_cuda/tools/{restart_field.py,convert_species_field.py,interp_field.py}`、
`design/forge_design/evaluate/runner_axismach.py` (`prepare_ns`・`run_staged_ns`・`run_staged`)、`design/forge_design/report/nozzle_report.py`、
参照例 `case/45.isobutane_m6_d155/run_0094_ns_c2pin_pass2_ext6k/{bcondConfig.yaml,solverConfig.yaml,prepare_info.json}` (ローカルにある)。

## 観測事実

plan §3 のとおり。bcond は flow 形式の 1 行 1 境界、入口は `inlet_Pressure` (Y0 MIXDRY, Y1 H2O, Pt, Tt, k, omega)、出口 `outlet_statPress` (Ps, Pt, Tt)。
restart_field は化学種の互換ハッシュが違えば書き込まずに停止。細分格子 NS は本段 cfl 5 で発散した (段階起動 + cfl 1 が実績)。

## 問い

Q1. 置き場所 (`solver_density_cuda/tools/` の汎用ツール) と CLI の形は妥当か。
Q2. H2O 分率を変えたとき乾き成分 Y0 を 1 − ΣY で自動に合わせる既定は妥当か (エラーにすべきか)。
Q3. lump (乾き成分の組成) を変えたときの初期場: convert_species_field (conserve/reinit) に自動で切り替えるか、停止してユーザに選ばせるか。
Q4. 条件の大きな変更に対する初期場のスケーリング (Pt 比で保存量を一様スケール、Tt 比で速度・温度) を入れるべきか — NS (SST、等温/断熱壁) で何が壊れうるか。入れないなら、起動は段階起動 (full) を既定にすべきか。
Q5. 入口の k・ω を Pt・Tt に合わせて自動調整すべきか (乱流強度・粘性比の保持)。
Q6. §6 の検証 (単体 4 項目 + case/45 の無変更再実行・Pt 0.8 倍) の合否の数値。無変更再実行の再現幅は forge の再実行揺れ (run_0104 で δ_E 0.003 %) をどう扱うか。
Q7. 見落としている罠 (例: `prepare_info.json` の `ic_from`・`Md`、`nozzle_report` の Euler 参照、壁温 BC、probe、凝縮 ON の run を参照にした場合)。
