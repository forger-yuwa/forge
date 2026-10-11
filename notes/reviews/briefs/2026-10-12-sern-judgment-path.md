# 諮問ブリーフ: SERN 生産の判定区間の経路 (plan tooling-sern-te-wake-grid §5.1 R2) と受入れへの進め方

- 日付: 2026-10-12、ブランチ `feature/sern-design` (作業ツリー `/home/sano/work/forge-sern-design`)
- 関連 plan: [`plans/accepted/tooling-sern-te-wake-grid.md`](../../../plans/accepted/tooling-sern-te-wake-grid.md) (§4・§5.1 R1〜R3・#7・§6・§6.1)
- 前回の result レビュー: [`notes/reviews/2026-10-11-tooling-sern-te-wake-grid-result.md`](../2026-10-11-tooling-sern-te-wake-grid-result.md) (NO-GO、M1〜M3 を全件採用 → R1〜R3)
- 諮る理由: R2 の担当は「O (経路の設計は F)」。実装担当が経路の設計で 8 点を選んだので、その採否を決める (AGENTS.md エスカレーション条件 1・7)。あわせて、この plan を result 段の再レビューに回してよいかを決める。

## 1. 観測事実 (実装と試験)

### R1 (学習の採否に床事象の要件)
- `design/forge_design/opt/driver_sern.py:167 floor_requirement_of`: 基準 YAML の `evaluate.floor_events` / `floor_events_window_steps` から要件 {必須, 区間長の下限, 検査更新数の下限} を作る (0 なら要件なし)。
- `:182 _floor_op_ok` / `:195 _floor_row_ok` / `:203 _learnable(..., floor=None)`: 各作動点で verdict PASS・区間 b−a ≥ 要件・`n_updates_checked` ≥ 要件。欠測・短窓・FAIL・判定不能は除外 (台帳には残す)。
- 要件は評価方式 ID ではなく**キャンペーンの要求**に置いた (より長い窓の行を採れる。旧評価の持ち越し承認は要件を緩めない)。
- 試験 `design/tests/run_sern_floor_learn_tests.py` 51 PASS (前回レビューの再現 = ゲート任意で FAIL の行も拒否)。

### R3 (評価方式 ID から設計変数を外す)
- dv の 5 変数は格子パラメータに MOC 経由の 2 つの角度でしか入らず、角度は既に "design" 置換済み。残る設計変数は `mesh3d.L_sw` だけ。
- `design/forge_design/evaluate/runner_sern.py:452 mesh_recipe`: 値つきの L_sw を "design" に置換 (None = 側壁がカウル後縁まで、は方式として残す)。値は `eval_method.design_values` に記録。`EVAL_METHOD_VERSION` は据え置き。
- 試験 `run_sern_eval_method_dv_tests.py` 33 PASS (run_1084 / run_1096 の YAML で ID 一致、L_sw_exact・層数・band/te_wake 0・変換器 sha では変わる、R3 前の定義の ID の行は除外)。

### R2 (判定区間の段)
- `runner_sern.py:1630〜` `judgment_options` / `judgment_config`、`:1796 run_judgment`、`sern_gates.py:304〜` `judgment_window_gate` (#3)・`judgment_temperature_gate` (#2)・`judgment_gate`・`judgment_extendable`、`runner_sern3d.py` の te_monitor の呼び出しと main。
- 生産 YAML (`case/46.sern_design/problem_3d_prod_3op_wallres_lswx08.yaml` / `lswx10.yaml`) に `judgment_steps: 20000`・`judgment_extend_once: 1`。
- 試験 `run_sern_judgment_tests.py` 134 PASS: **保存済みの検証 run 12 本**に判定を通すと、a・d・STEADY が手作業の判定 (`V3_EVAL` / `QUASISTEADY`) と一致。run_1087 は低温 FAIL (節点 517159、102.53 K)、run_1085・1088・1097 は未定常、1084・1086・1091・1096・1098・1099 は PASS。
- GPU 接続試験 `run_sern_judgment_connection_test.py` 16 PASS (2.4 万節点の粗い 3D、窓 40・N 80): 床事象の判定 PASS・区間 (40, 80]・更新 40 (前回レビューが再現した判定不能は解消)、config の差は nStepOuter だけ、最後の 1 枚だけ残る。延長しない経路: 1 回目が rms_roK の上昇 (RESIDUAL_RISING) で FAIL → 規則どおり延長せず。延長する経路: 延長の判断だけを強制 → 2 回目が 1 回目の res_80 から継続、1 回目の記録が `*_judgment1.*` に退避されることを確認 (2 回目自体は粗い流れの発散で落ちた)。
- **延長を経て PASS する経路は偽の forge でしか確かめていない**。
- 既存の試験 10 本 (zte_wiring 150・floor_events 26・te_wake 122・te_wake_local 68・mesh_provenance 149・mesh3d 64・twin_ic・floor_learn・eval_method_dv・judgment) と `test_gate_bad_input.py` は手元で再実行して PASS (2026-10-12)。範囲外の `run_sern_gates_tests.py` の FAIL は HEAD でも同じ (親 plan R13)。

## 2. 期待値と出典

- plan §5.1 R2: 「段階起動 → 同一設定 20000 step の判定区間 → §6 #1〜#3 の判定 → 窓が未達なら +20000 を 1 回」を標準経路に。**窓を短くして通さない**。判定区間の run を同じディレクトリで続けるか別にするかを決め、床事象の記録の session・区間と整合させる。生産 YAML で標準経路が床判定まで通ることを小さい格子で接続試験。
- plan §6 (事前登録): 各 run 20000 step、窓未達なら +20000 を 1 回、判定は末尾 10000 step。窓 a = max|q − mean|、d = 前後の窓の平均差、許容 5e-5 (C_T・C_T_ws・C_L)・5e-4 (C_M)。#2 は後縁まわりの監視領域で T < 150 K の節点 0 と温度の STEADY。
- 検証での手作業の手順: `case/46.sern_design/diag/tewake_ext.sh` (延長は前の run の最終場から `restart_field.py` で新しい run へ)、`te_monitor.py` (R_TE/R_SE/RET/ALL の Tmin CSV)。

## 3. 実装担当が選んだ経路の設計 8 点 (採否を求める)

1. **ディレクトリ**: 同じ run ディレクトリで forge の別 session として続ける (`layout: same_run_dir`)。判定区間が最後の session になり、`floor_events.csv` の session と collect の判定が一致する。検証では延長を**別ディレクトリ**で回していた。
2. **段の移行**: `restart_by_index` + 保存量のビット一致の検査 (`_restart_mismatch`)。`res_*` を消してから回す (検証では `restart_field.py`)。
3. **窓**: W = `floor_events_window_steps` を #1〜#3 で共通に使い、必須。N ≥ 2W、outStepInterval は W と N を割り切り、窓あたり 20 標本以上。満たさない設定は起動前に止める (窓は短くしない)。
4. **metrics の正本**: 最後の判定区間の `metrics.json`。延長前の記録は `*_judgment1.*` に退避、経過は `JUDGMENT.json`。
5. **延長の条件**: 未達が「未定常」(UNSTEADY と JUDGMENT_UNSETTLED) だけのとき。低温・床事象・判定不能・rc≠0・残差の上昇があれば延長しない。延長の 2 回目は 1 回目の最終場から継続し、#3 の d は 2 回目の中の前後の窓で取る (検証の +20000 と同じ)。
6. **延長後も未達**: `JUDGMENT.json` は INDETERMINATE (§6 #3 の判定不能)、ゲートは FAIL/JUDGMENT_UNSETTLED で学習に入らない。
7. **te_monitor の適用範囲**: te_monitor は g3 生産形状の座標 (後縁、x=0.08 のカウル y、W/2、L_sw ≥ 0.08 m) を定数で持つ。run の形状が合わなければ段を回さずに INDETERMINATE。RET の長さは L_b (0 の格子は検証と同じ 1.0)。
8. **2D**: #2 の監視領域が無いので、段を回さずに INDETERMINATE (2D plan の「必要な領域が空なら判定不能」に合わせた)。2D の生産 YAML には判定区間のキーを入れていないので、現行の 2D 経路は変わらない。

## 4. 残っている穴 (実装担当の報告)

- (a) te_monitor の領域を形状から作る処理が無い → 生産 YAML 以外の設計点 (L_cowl・H・W が違う) は判定不能。
- (b) `driver_sern._eval_op` (2D キャンペーン) は判定区間の段を通らない。2D の #2 も未実装。
- (c) L_sw は driver の x (`DV_ORDER`) に入っていない。3D で L_sw を変えて学習するなら x に入れる必要がある (台帳には `design_values` で残る)。
- (d) 既存の run_1084 / run_1096 などの prepare_info は旧定義の ID のまま。取り込むには承認による持ち越しか prepare のやり直しが要る。
- (e) `FORGE_CUDA_BLOCKSIZE=128`・`FORGE_ALLOW_UNVERIFIED_SPECIES=1` は runner が設定せず環境に依存 (検証の手順と同じ)。

## 5. 仮説・問い

- Q1. 設計 8 点のうち、変えるべきものはあるか。特に 1 (検証は別ディレクトリ、生産は同じディレクトリ)・2 (`restart_field.py` でなく index コピー + ビット一致検査)・5 (延長の条件)・7・8 (判定不能に倒す) の妥当性。
- Q2. 「延長を経て PASS」の経路を偽の forge だけで受け入れてよいか。それとも本物の GPU で通す確認 (例: 生産 YAML を g3 で AWS-B、数時間) を plan の受入れ条件にすべきか。粗い格子では 2 回目が発散した (流れの問題で経路の問題ではない、と実装担当は見ている)。
- Q3. 残る穴 (a)〜(e) のうち、この plan の受入れ (status done) の前に塞ぐべきものはどれか。後回しにするなら行き先 (親 plan R7b・2D plan・新規 plan) をどう分けるか。とくに (d) は R7b の設計判断 (L_sw 0.8 対 1.0) の材料に直結する。
- Q4. 上を反映した後、この plan を codex の `--stage result` の再レビューに回してよいか (AGENTS.md 条件 7)。

読んでよいファイル: plan 本体、`design/forge_design/evaluate/runner_sern.py`・`runner_sern3d.py`・`design/forge_design/metrics/sern_gates.py`・`design/forge_design/opt/driver_sern.py`、`design/tests/run_sern_{floor_learn,eval_method_dv,judgment}_tests.py`・`run_sern_judgment_connection_test.py`、`case/46.sern_design/diag/{te_monitor.py,tewake_ext.sh}`、`case/46.sern_design/README.md`、前回の result レビュー。
