# 諮問ブリーフ: 移植の回帰の回し直し (§6.60) の結果と、main への統合に進んでよいか (2026-10-13)

- 依頼者: 主セッション (エスカレーション条件 3: 事前登録の比較が判別不能、条件 7: 解釈の確定前)
- plan: `plans/active/architecture-float-state-double-geometry.md` §5.1 #26、§6.59・§6.60 と結果。統合の記録 `notes/investigations/2026-10-13-main-integration/README.md`
- ユーザ: main へ統合してよい (2026-10-13)。前回の諮問 (`2026-10-13-port-reg-redo-diagnose.md`) は「統合は回し直しの結果を待つ」。

## 観測事実 (§6.60、判定 `_band_ab/cold_pair/port_reg2_judge.json`)
- 各組で移植の前 3 本・後 3 本、同じ入力から 10 step、A で GPU を占有して前・後を交互に。D = max|群の平均の差|、E = 各群の中の (最大 − 最小) の最大。τ: 状態・原始量 1e-4 × max|平均|、`res_*`・`omg_*`・`dt_local` 1e-1 × max|平均|、rms 列 step ごとに 1e-2 × |平均| の最小。
- FP64 の 3 組 (R1 point・R1 B0・R3 平板): 結果 A ((D + E_A + E_B)/τ の最大 0.07・0.22・1.2e-4)。
- float の 3 組: 判別不能。許容を超えた量 6/48・7/49・6/31。どの量も D < E_A + E_B (D/(E_A+E_B) の最大 0.50・0.40・0.35)、E_A ≈ E_B。結果 B (D > τ + E) の量は無い。例: R1 B0 float の `vis_turb` D 1.17e-5・E_A 2.34e-5・E_B 2.70e-5・τ 2.19e-6、`res_roOmega` D 7.0e3・E 1.74e4・1.77e4・τ 3.6e3、`roUy` D 7.3e-2・E 0.143・0.130・τ 0.165 (この量は許容内)。R3 float の `rms_roUy` D 1.4e-6・E 8.8e-6・7.9e-6・τ 9.9e-7。
- 補助の証拠: SASS の比較で、`lowMachPrecond < 2` の経路のカーネルは移植の前と同一 (前処理のカーネルと新しいカーネルを除く)。R2 (§6.59) は結果 A。
- 実行の記録: A のディスクを一杯にした (入力の写しと res_0 を全部残した)。自分の run の不要なものを消して回復、ほかの forge は走っていなかった。

## 問い
1. §6.60 の結果から書いてよい範囲。
2. 統合 (main) に進んでよいか。進むなら、float の判別不能をどう記録するか。進まないなら、何を足せば判断できるか (例: float 用の許容を別に事前登録して新しい試験として回す、本数を増やす、など。結果を見た後で許容を決める問題に注意)。
3. 統合の記録の確認 4 (生産の問題の準備の再現) が入力の欠損で実行不能だった件の扱い。

## 読んでよいファイル
- float plan §5.1 #26、§6.59・§6.60 と結果、`case/45.isobutane_m6_d155/port_reg2_an.py`・`_band_ab/cold_pair/port_reg2_judge.json`
- `notes/investigations/2026-10-13-main-integration/README.md`
