# 諮問ブリーフ: ③ベル (M) の V2 (Euler の推力と MOC の推力の整合) が登録の線の外

- 日付: 2026-10-12
- 依頼元: 主セッション (Opus)。エスカレーション条件 3 (事前に書いた参照値・許容差との比較が FAIL)・7。このセッションの諮問は diagnostician (Fable) — ユーザ指示
- 作業ツリー: `/home/sano/work/forge-thruster` (HEAD `7e89630d` + README の行の追記のみ)。run は AWS B の `~/forge-thruster/case/68.bell_two_mode/` (手元に結果ファイルは無い、数値は下に写した)
- 対象: plan `plans/active/tooling-nozzle-bell-two-mode.md` §6 V2 (「(M) の非粘性壁を §4.1 の Euler 設定で回し、推力が MOC の値と相対 0.2 % 以内 (粗い整合確認)」)・§5.1 #6 (W2・W3)・§4.4 (初期線は CFD ピンが既定、Hall は検算)

## 観測事実

- run: `run_0015_tp_m_euler_l100` (壁 `rao_m/rao_t6_parts/L100_n321_mq2_wall.npz`) と `run_0016_tp_m_euler_l70` (`L70_n321_mq2_wall.npz`)。どちらも (M) の Hall 初期線 (n 321) の壁、`evaluate.model: euler` (`viscMethod 0`・`visc 0`・`thermCond 0`・`turbulence none`・物理壁 slip)、格子 ni 1201・nj 121・全断面一様の `wall_first_frac 0.005` (AR 最大 2.75/3.32)、soft 2000 → mid 2000 → 本段 30000、出力 1500 ごと。コード `7e89630d`。
- ゲート (`driver --gate-only`): 両方 PASS — 残差 `NOT CONVERGED (stalled)` (併記)、推力系 4 列と ṁ ALL STEADY、CV 残差 +8.42e-7 / −9.19e-7、メッシュ品質 PASS、壁解像は Euler なので対象外。NaN なし。
- 最終値 (step 30000):

| | Euler η_exit | Euler η_wall | Euler ṁ/ṁ_1D | MOC η_cs | MOC ṁ/ṁ_1D | F_Euler/F_MOC − 1 |
| --- | --- | --- | --- | --- | --- | --- |
| L100 | 0.990362 | 0.990363 | 0.992630 | 0.987530 | 0.989706 | +0.287 % |
| L70 | 0.977316 | 0.977315 | 0.992627 | 0.974509 | 0.989706 | +0.288 % |

- 流量で割った比 (η / (ṁ/ṁ_1D)、流量 1 単位あたりの推力の比): L100 で Euler 0.997715・MOC 0.997801 (差 −8.6e-5)、L70 で 0.984573・0.984645 (差 −7.3e-5)。つまり +0.29 % の差はほぼ全部が流量比の差 (+0.295 %) で、流量あたりの推力は 1e-4 以内で一致。
- MOC の流量 = Hall 初期線 (R_t 1、n 321) を横切る流量。Hall の級数は R_t 1 で C⁻ の適合式と 1.1e-2 rad 合わない (plan §5.1 #3b)。Euler の流量は 2 本で 2.6e-6 しか違わない (収縮部・スロートが同じなので当然)。
- MOC の参照値の換算: F[N] = F_nd·ρ*·a*²·r_t² (真空)、ρ* a*² は同じ凍結 TP の 1D スロート状態 (`ideal_1d_tp`)、F_1D 1,159,112.9 N。

## 主セッションの案

V2 を「登録の線では FAIL (+0.29 %)」と記録し、差の内訳 (スロートの流量 = Hall 初期線の流量の誤差 −0.29 %、流量あたりの推力は 1e-4 以内) を併記する。予定どおり W3 (この Euler 場を凍結源にした CFD ピンの初期線で (M) を作り直す) に進み、CFD ピン版の (M) 壁で V2 をもう一度回す (同じ線 0.2 %)。

## 問い (推奨を 1 つずつ)

1. この V2 の結果をどう判定・記録するか (登録の線は推力そのもの)。流量あたりの推力の一致を V2 の判定に使ってよいか、それとも線は変えず「FAIL・原因は初期線の流量」として W3 の後に回し直すか。
2. W3 の後の V2 の回し直しの事前登録: CFD ピン版の (M) 壁で、推力の線 0.2 % をそのまま使うか、流量あたりの推力の線を足すか (値も)。
3. Hall の流量の誤差 −0.29 % は、(D) の MOO と (M) の比較 (V3) にどう効くか — (D) は CFD で評価しているので Hall に依らない。(M) の L_k の形 (θ_B・θ_E) は初期線の流量で変わるはず (W3 で測る)。V3 の前に必要なことはあるか。
4. この 2 本の Euler 場を CFD ピンの凍結源に使ってよいか (`feedback/cfd_initial_line.py` の「凍結源: 同じ形・ガスで Hall 初期線の V0 型壁を Euler で解いた場 (ピン壁自身の Euler は循環定義になるので使わない)」— この 2 本は Hall 版の (M) 壁で、ピン版の壁とは別物)。格子 (全断面一様 0.005、ni 1201) はスロート付近の初期線の抽出に十分か。

## 読んでよいもの

- `plans/active/tooling-nozzle-bell-two-mode.md` (§4.4 の「初期線」・§5.1 #3b・#6・§6 V2・V3)
- `design/forge_design/feedback/cfd_initial_line.py` (冒頭 1〜30 行、`CFDPinnedThroat`・`check_source_matches`)
- `design/forge_design/geometry/moc_forward.py` (`hall_initial_line` 119〜141)
- `case/68.bell_two_mode/problem_bell_tp_m_euler_l100.yaml`
- 禁止: `*.log`, `res_*.h5`, `*.pdf`, `plans/README.md`。ファイルは変更しない。AWS には入らない
