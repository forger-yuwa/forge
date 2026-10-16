# 諮問ブリーフ: M5.5 の NS の全部 double の延長 (case/71 run_0003_ns_cont) が予算の上限 20 万 step で停止規則に僅かに届かず — 次の一手と解釈 (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 3・4・7: 事前の規則に届かない、plan に無い再実行の前、結果の解釈の前)
- plan: `plans/active/tooling-nozzle-card-local-mesh.md` §5.1 #3b・#3c・#3d (#3d が今回の run の条件と監視の項目)
- 前の諮問: `notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.md`・`-2-diagnose.md`・`-3-diagnose.md` (3 回目は thermoFloat 1/0 の 1 更新の A/B を先に、延長は保留を提案 → ユーザの決定で延長を全部 double にして代えた)
- 記録 (手元): `case/71.m55_combustion_wt/_runs/run_0003_ns_cont/` (台本 `cont_0003.py`・`judge_0003.sh`・`monitor_0003.py`、判定の原文 `judge/`、`stop_rule.json` の全体、`solverConfig.yaml`・`bcondConfig.yaml`・`RUN_PROVENANCE.txt`・`CONVERGENCE_VERDICT.txt`)

## 再現条件

- run: AWS A `case/71.m55_combustion_wt/run_0003_ns_cont`。run_0002_ns の res_183000 を同じ nozzle.h5 に `restart_field.py --keep-src-dtype` で index restart (9 量、ビット一致)、前段なし、`run_staged_ns(stages="none", stop_rule=True)`。
- run_0002_ns の本段からの差: `physProp.thermoFloat: 0`・`time.deltaT.implicitSolvePrecision: 1`・`outStepInterval` 2,000、環境変数 `FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1`。FP64 のビルド (flowFormat.hpp double、sha256 3a9e347b…)。float のまま残るのは化学種の拡散の D_s・h_s・X、リミッタの長さの sqrtf、XMF の型だけ。
- 停止規則: 窓 2 万・保持 4 万・上限 200,000 (= 追加の予算、plan #3d)。r ≤ 5e−4 (平均)・振幅 ≤ 5e−3、θ_r の D ≤ 5e−4・振幅 ≤ 5e−3。帯は延長の res_2000 から、Euler 参照 run_0001_euler/res_12000 (sha256 93cb88b0…)。θ_r は run_0002_ns とは別の系列。
- 26.9 ms/step、200,000 step で forge は rc 0 で終了 (01:54 UTC)。見張りの台本はその後も記録用の δ の抽出を待っていて、stop_rule.json の status はまだ RUNNING (終了の状態は CENSORED になるはず)。

## 観測事実 (判定区間は延長の本段だけ: step 0〜199,999、設定が違うので run_0002_ns とは連結しない)

1. 停止規則 (最後の出力 200,000 の判定): **条件を満たさない (theta_r_0.4, theta_r_0.7)**。r の平均 1.07e−4 (ok)・振幅 1.76e−5。D: 0.4 L 6.32e−4、0.7 L 5.33e−4、0.95 L 4.84e−4 (ok)。D は 1 万 step ごとに約 0.815 倍で減っている (30k 1.21e−2 → 100k 4.14e−3 → 200k 6.32e−4)。r も単調に減る (10k 2.76e−3 → 100k 5.51e−4 → 200k 9.07e−5)。
2. `check_quasisteady --series-csv` (停止規則の系列 100 行、既定の閾値) DRIFTING: r DRIFTING (0 へ向かう単調減少)、θ_r 3 断面 STEADY (単調増加、漸近値は最終の +0.135 %・+0.116 %・+0.106 %)、ṁ_in・ṁ_out STEADY (漸近 121.358・121.355 kg/s)、出口コア M STEADY (5.52045、漸近 5.52042)。
3. `check_convergence` (履歴そのまま = 本段だけ; stages=none なので stage_manifest.json が無く `--segment` は判定区間を決められない) NOT CONVERGED: ρ 1.7 桁 falling、ρUx・ρUy 2.0 桁 flat、ρE 1.6 桁 flat、k 1.8 桁 falling、ω 3.1 桁、Y 1.7 桁 falling。基準は restart 直後の値 (延長の最初の値 1.70e−6 ≈ 最大 1.74e−6)。
4. 壁解像 (`check_wall_resolution --groups wall --weight area --step 200000`) PASS: y1+ 平均 0.479・p99 0.780・最大 1.002 (x −1 m、入口の角)、超過の面積 0.0 %。
5. 最後の場 res_200000 (`final_field_check.py`): 非有限 0、ρ・P・T > 0、P > P0 0 点。**P0 > Pt (1 + 1e−6) が 46,936 点・最大 P0/Pt − 1 = +0.462 %** (run_0002_ns の res_183000 は 50,441 点・+0.236 %)。位置は入口の列以外・壁から 27〜90 点・軸の行 0 点、x −1.0〜8.0 m。T0 > Tt + 1 K は入口の列の 9 点 (最大 +3.39 K)。
6. 監視の時系列 (`monitor_0003.py`、11 時点 2,000・2 万・…・20 万、`judge/monitor_series.csv`; `check_quasisteady` の判定を括弧に):

   | step | 超過の点数 | 最大 ΔP0/Pt [%] | S の平均 ΔP0/Pt | S の平均 h0 − h0_in [J/kg] | E_in [W] | E_out [W] | Q_wall [W] | (E_in − E_out + Q_wall)/E_in |
   | --- | --- | --- | --- | --- | --- | --- | --- | --- |
   | 2,000 | 50,069 | 0.237 | 1.44e−4 | −2,208 | 1.395028e8 | 1.338610e8 | −5.677e6 | −2.53e−4 |
   | 40,000 | 49,102 | 0.308 | 1.04e−4 | −3,088 | 1.393835e8 | 1.338092e8 | −5.541e6 | +2.37e−4 |
   | 100,000 | 47,359 | 0.453 | 6.15e−5 | −3,821 | 1.393002e8 | 1.337687e8 | −5.449e6 | +5.90e−4 |
   | 160,000 | 46,964 | 0.460 | 4.57e−5 | −4,089 | 1.392735e8 | 1.337545e8 | −5.420e6 | +7.09e−4 |
   | 200,000 | 46,936 | 0.462 | 4.16e−5 | −4,155 | 1.392675e8 | 1.337510e8 | −5.414e6 | +7.37e−4 |

   判定: 超過の点数 STEADY (漸近 46,924)、最大 ΔP0/Pt STEADY (漸近 +0.463 %)、S の平均 ΔP0/Pt DRIFTING (減少中)、S の平均 h0 − h0_in DRIFTING (−2,208 → −4,155 で減少中、増分は減衰)、E_in・E_out・Q_wall STEADY、(E_in − E_out + Q_wall)/E_in DRIFTING (−2.5e−4 → +7.4e−4、増分は減衰)、入口の列の静温の超過 STEADY (+2.64 K)。
   - S = run_0002_ns の res_183000 で P0 > Pt (1 + 1e−6) の節点 (50,441 点、#3c と同じ定義)。h0_in = 1,147,573.50 J/kg。
   - E は入口・出口の列の折れ線を通る ∫ρ h0 (U·n) 2πr ds (`stop_rule.mass_flows` の式に h0 を掛けた)、Q_wall は `res_wall_3_<n>.h5` の qwall を壁の節点の x の順に 2πr ds で台形積分 (符号はソルバのまま、負 = 気体から壁へ)。**この Q_wall がソルバの壁の CV の熱流 (node の半割面) と同じ量かは確かめていない**。
   - (E_in − E_out + Q_wall) = +1.03e5 W (200k) = 気体が入口と出口の間で失うエンタルピー 5.517 MW のうち 1.9 % が壁の熱流で説明されない。r (9.1e−5) × E_in ≈ 1.3e4 W。
7. 容量: A の空き 5.3 GB、run_0003 は 10.6 GB (全場 101 枚 × 100.6 MB)。

## 期待値と出典

- 停止規則の条件 (plan #3d、`stop_rule.py`)。予算は追加 20 万 step (plan #3d: 「到達の予測ではない、2 万 step ごとに再評価」)。
- 全圧: 冷却壁のノズルで主流の P0 は Pt を超えないのが期待 (等エントロピー以下)。ただし冷却 (熱を奪う) はエントロピーを下げうる。超過の点は壁から 27 点以上離れた主流〜軸の側。
- エネルギー: 定常なら E_in − E_out = −Q_wall (壁の熱流の向きの規約による)。

## 呼び出し側の仮説・案 (結論ではない)

1. 停止規則: D は約 0.815 倍 / 1 万 step で減るので、同じ規則で 0.4 L が 5e−4 を切るのは約 1.2 万 step 先。新しい run は窓を数え直すので、到達は 2 万 step 前後、保持 4 万で計 約 6 万 step (約 27 分)。ただし予算の外。
2. 代わりに、θ_r・出口 M・ṁ が STEADY (漸近値と最終の差 ≤ 0.14 %) であることを使い、停止規則は CENSORED のまま、最初の pass として δ* の補正に渡す道。
3. 全圧の超過は全部 double でも残り、最大は 0.24 % → 0.46 % に増えて STEADY。thermoFloat (混合精度) は主因ではない見込み。S では h0 が不足 (−0.36 %) しつつ P0 が超過している。エネルギーの収支の残り 1.9 % (壁の熱流の積分の定義の差か、数値の非保存か) は切り分けていない。

## 問い

1. 停止規則に届かせるために延長を足すか (案 1、約 6 万 step、予算の追加はユーザの了承が要る)、案 2 で先へ進むか。案 1 なら上限と判定の置き方 (窓の数え直し、再開の跳ね) をどうするか。
2. 全圧の超過 (+0.46 %) と、S の h0 不足の同居をどう扱うか。「ソルバに欠陥は無い」とは言えないとして、δ* の補正 (θ_r と出口 M を使う) の判断に影響するか。切り分けに使える最小の A/B (forge 0 step で既存の場から、または短い run) は何か。
3. エネルギーの収支の残り (E_in − E_out + Q_wall)/E_in = +7.4e−4 (DRIFTING、増分減衰) の評価。qwall の積分の定義を確かめる方法 (ソルバの node の壁の熱流と一致させる手順)。
4. 結果の書き方: 残差 NOT CONVERGED (運動量・ρE が 2.0・1.6 桁で flat) を、θ_r などの量の STEADY とどう並べて報告するか。

## 読んでよいファイル

- `plans/active/tooling-nozzle-card-local-mesh.md` (§5.1 #3a〜#3d、§6)
- `case/71.m55_combustion_wt/_runs/run_0003_ns_cont/` 全部 (`judge/` の原文、`stop_rule.json`、`monitor_0003.py`、`judge_0003.sh`)、`case/71.m55_combustion_wt/_runs/run_0002_ns/` (`p0_ab.py`・`p0_ab.json`・`final_field_check.py`・`judge/`)
- `design/forge_design/evaluate/stop_rule.py`、`design/forge_design/evaluate/runner_axismach.py` (`run_staged_ns`)
- `solver_density_cuda/tools/total_quantities.py`、`solver_density_cuda/tools/check_quasisteady.py`
- `solver_density_cuda/cuda_forge/` の壁の熱流 (`qwall` を書く箇所、`viscousFlux_d.cu` など) と等温壁の境界 (node)
- 前の諮問の記録 3 本
