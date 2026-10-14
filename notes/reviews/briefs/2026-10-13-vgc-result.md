# 諮問ブリーフ: V-g0c-chain の結果の解釈と、壁温を変える計算への進み方 (2026-10-13)

- 依頼者: 主セッション (エスカレーション条件 7: result の解釈の確定の前、条件 1: 次の方針)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.8・§5.1 #29〜#31・§6 V-g0c-chain (事前登録と結果)、§8 (未確定事項、壁温)。関連 `plans/active/architecture-float-state-double-geometry.md` §6.58・§6.61・§6.62 の結果、§5.1 #25・#28
- ユーザ: 「この計算をもう少し早くしないと、先に進みたくない」(10-12)。「ちょいとマッハ数高くなってるけど、壁の調整をしてないからだよね。壁温上げてから。」(10-12、冷却壁 300 K の NS の出口マッハ数が設計の 6 を少し超えていることについて)。G0c は「問題がおきなかったならこれでもよい」。低マッハの前処理の調査は別セッションの担当 (10-13)。

## 観測事実

- V-g0c-chain (`case/45.isobutane_m6_d155/run_0810_vgc_chain`、判定 `_band_ab/cold_pair/vgc_judge.json`): 結果 A。チェーンの `prepare_ns` が G0c を作り (`mesh.x_respace`、.msh が既存の G0c とバイト一致)、`run_0179` res_20000 から格子をまたいで移し、段階起動 S1・S2 (1 次、各 3,000 step) の後、本段だけ B0 (`evaluate.main_line`) で 120,000 step で停止規則 (§6.11) に到達。T2a・T2b と M1〜M7 が許容内 (θ_r +0.076/+0.068/+0.062 %、代理 Q_w −0.040 %、q_w の最大の相対差 0.69 %・超過 0、出口 M の最大 |ΔM| 0.0005)。対照 `run_0804_sp_g0c` (同じ G0c・B0 を Q32 の起点から) との差は θ_r 0.003 %、Q_w −0.001 %。前提の検査 (1)〜(6) 合格。`check_convergence --segment` (本段) は ALL PASS だが起点の残差が大きいため (`run_0804` は収束に近い起点から FAIL)、準定常は θ_r DRIFTING・Q_w TRANSIENT-UNSETTLED。本段 8,884 s (GPU を共有、参考)。
- 準備で見つけたこと: (i) チェーンの初期線の凍結源 `run_0062` の `nozzle.h5` が 10-10 の掃除で消えていた → B の写しに `res_6000.h5` の `MESH` だけを写して補った (正しさは .msh のバイト一致で確認)。(ii) 変換器の終了時の既知の GPU エラーでチェーンが止まる → `runner.run_converter` で、書き出しの行・その GPU エラー・出力の HDF5 の 3 つがそろうときだけ通す (018d993d)。
- §6.61 (占有、登録の Q32 の起点): G0c・B0 は約 40.8 分で T2 との比較の要件を満たし、G0・B0 は止まるが θ_r が +0.95〜+1.23 % 離れる。§6.62 (試験 P): point + `lowMachPrecond 3` は資格成立のうえ持続する悪化 → 長い試験から外した。
- 未決: main への統合は保留 (float の移植の回帰が判別不能、§5.1 #27)。`thermoFloat 1` が FP64 のビルドでも float の熱物性の演算を残す (forge-25 の知らせ、こちらでは未確認)。

## 問い

1. V-g0c-chain の結果 A から書いてよい範囲。§4.8 の G0c・B0 のチェーンへの組み込み (#29〜#31) を「完了」としてよいか。残すべき但し書き (準定常が DRIFTING、時間は共有の参考、run_0062 の補い、AR 緩和の範囲)。
2. ユーザの次の段 (壁温を上げて壁を調整する = 壁温を変えた問題で δ* ループを回す) に G0c・B0 を使ってよいか。壁温を変えると旧 T2 は参照にならない (前の諮問の注意)。壁温を変えた問題で G0c・B0 を使う前に、何を事前登録して確かめるべきか (例: 新しい壁温で G0 の上限 1 の参照を 1 本作る、または G0c の上限 1 と上限 50 の対だけで足りるか)。費用を抑える最小の確かめ方を提案してほしい。
3. δ* ループの各 pass (壁が変わる) で G0c・B0 を使うときの扱い。`x_respace` の区間 (x/r_t −8〜−2) は r_t で正規化されているが、壁が変わると縮流部の形も変わりうる。pass ごとに何を確かめれば足りるか。
4. いま決めずに残すべきこと (float のビルド、main への統合、thermoFloat の基準)。

## 読んでよいファイル

- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.8・§5.1 #29〜#31・§6 V-g0c-chain と結果・§8
- `plans/active/architecture-float-state-double-geometry.md` §6.58・§6.61・§6.62 の結果、§5.1 #25〜#28
- `case/45.isobutane_m6_d155/vgc_run.py`・`vgc_judge.py`・`_band_ab/cold_pair/vgc_judge.json`・`problem_d155_ns_prod_coldmesh_tw300_g0c_chain.yaml`
- `design/forge_design/meshing/mesh2d.py` (`_x_stations_respaced`)、`design/forge_design/feedback/deltastar_loop.py` (`run_pass`)
- `notes/reviews/2026-10-13-g0c-chain-design-diagnose.md`・`2026-10-13-speed-s-result-diagnose.md`
