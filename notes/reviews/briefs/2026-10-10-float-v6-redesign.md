# 諮問ブリーフ: V6 (標準ケースの回帰) を「float と FP64 の比較」に作り替え、case/48 から始めてよいか (2026-10-10)

- 依頼者: 主セッション (AGENTS.md のエスカレーション条件 1: plan §6 の検証計画を変える)
- plan: `plans/active/architecture-float-state-double-geometry.md` (§5.1 #10・#15、§6 の V6、§6.11・§6.16・§6.18〜§6.20)
- 関連: `plans/active/axisymmetric-freestream-hoop-gauge.md` §4.10 (hoop の修正の確認)、`case/48.flat_plate_cooled_m4/README.md`

## 0. 問い

1. ユーザの目標は「なんとか double をやめたい」(2026-10-10、FP64 のビルドをやめて float で生産を回したい)。今の V6 は「新しい float の経路が古い float の結果を壊していないか」の回帰で、この目標に答えない。V6 を「各標準ケースで、double の幾何の格子の上で float が FP64 と同じ所に止まるか」に作り替えるのは妥当か。
2. 最初の 1 本を case/48 (冷却平板 M4、平面の 2D・SST・等温壁 300 K) にするのは、case/45 の float のずれ (θ_r +1.3〜1.8 %) を「軸対称・軸の近くに固有」と「冷却壁の境界層」に分ける切り分けとして成り立つか。交絡 (§3) をどう扱うべきか。
3. case/48 の組み方と事前登録すべき判定 (§4 の呼び出し側の案) に穴はないか。閾値・反復の数・長さ・物差し。
4. 順番: float plan §6.20 (hoop の修正の格子での V4 のやり直し) の結果が 2〜3 時間後に出る。case/48 はその結果を待たずに回してよいか。

## 1. 観測事実

- **case/45 (軸対称・TP・ライン陰解法・冷却壁 300 K、57 万節点)**: float plan §6.16 (V4): float と FP64 が同じ停止規則の水準に到達したとき、float は θ_r(40/70/94) が +1.82 / +1.43 / +1.29 %、Q_w が +1.43 % 違う所で止まった。差の大きい場所は収縮部 (x/r_t −5〜−1、j 20〜60: |U| 最大 45 %、μt 最大 30 倍の差)。
- §6.19 (精度の切り替えの A/B): float の到達の場から FP64 に切り替えると、3 万 step で θ_r は FP64 の値へ戻った (−1.4〜−2.2 %)。float で続けると ±0.1 % 以内に留まった。Q_w は再実行の揺れ (A と A2 の差 0.225 %) で事前登録の条件を外れ、判定は「判別不能」。
- **hoop の欠陥** (hoop plan §4.10、修正済み・commit `0b4dff4e`): 旧来の格子 (キー 0) では、FP64 でも一様な静止場で軸の近く j 2〜8 に |res_roUy|/(P·A) = 3.3e-3 の偽の力が出ていた (float でも同じ大きさ)。修正後 (キー 1) は FP64 で 5e-14、10 step 後の最大 |u| は 0.257 → 6.4e-5 m/s。V4 は修正前の格子だった。
- **実行中**: float plan §6.20 (事前登録済み) で、修正後の格子 (キー 1) の float (`run_0484_fl1_f32_k1`) と FP64 (`run_0483_hp7_k1`) を V4 と同じ停止規則で回している。結果の区分 (a) 差が消える / (b) V4 と同符号で半分以上残る / (c) その間、を事前に決めてある。
- **case/48 の既存の基準** (README): `run_0025_B_tw300_y3_fx05` (float の旧バイナリ・float32 の格子、手順書どおりの冷間起動 48000 step): Cf/VD-II 0.969–0.994、2St/Cf 1.15–1.16、エネルギー閉合 1.016、`SERIES VERDICT: STEADY`。`check_convergence` は `NOT CONVERGED (stalled/plateau)` (前縁特異点の既知の残差床。基準 run も同じ)。
- **case/48 の再実行の揺れ** (README): 同一設定の 5000 step の反復で場の差が `ro` 5.9e-4・`roUy` 3.4e-3 (`run_0028`/`0029`、残差プラトーの有界振動で atomicAdd の非決定性が育つ)。200 step では場の差 1.8e-6 (`run_0013`/`0016`)。1 step でもビット再現しない (`run_0030`/`0031`)。S2 の双子 (gg/lsq) の物理量の差は Cf・q_w・δ*・θ で ≤ 0.016 %。
- **double の幾何の格子** (float plan §6.14・§6.15): 新しい変換器で作った case/48 の格子 (`_conv614/c48/new64/m.h5`、hoop の修正後の変換器とビット一致) を、新しい float / FP64 のソルバが読めることは確かめた。今日、この格子に `run_0025` の res_48000 を `restart_field.py` で移して (7 量、ビット一致) 新旧の float のバイナリで 20 step 回し、新旧の差 / 再実行の差 1.95 だった (hoop plan §4.6 の 6)。
- 1 step の時間: case/48 は float で約 1.5 ms/step (8.9 万節点)。FP64 は未計測 (2〜3 倍と見込む)。

## 2. 期待値と出典

- V4 の閾値 (§6.11): θ_r 0.05 %・Q_w 0.1 % (到達時)。case/48 には θ_r の断面の代わりに `cooled_plate_eval.py` の Cf・q_w・δ*・θ (x 0.3/0.6/0.9)、積分の CD・HF がある。
- 回帰の規則の罠 (メモリ): 「新旧の差 ≤ 旧の再実行の範囲」の形は同じ分布でも 8 割 FAIL する。閾値は再実行の揺れの何倍か、または絶対値で決める。

## 3. 交絡 (case/48 と case/45 の違い)

case/48 は平面 (軸対称でない)・CPG (TP の多成分でない)・ブロック DPLUR の点陰解法 (ライン陰解法でない)・cfl 2・境界層の外は一様な超音速流 (収縮部・スロート・軸がない)・格子は 1001 × 約 90 の直交。共通は冷却壁 300 K・SST (`dilatationCorrection 2`・`katoLaunder 1`)・低 Re 壁・node・SLAU・2 次。

## 4. 呼び出し側の案 (検証していない)

- **V6 の作り替え**: 各ケースで、新しい変換器 (double の幾何) で作った同じ格子・同じ初期場から、float (`~/forge-fgeom7-f32`) と FP64 (`~/forge-fgeom7-fp64`) を同じ設定で回し、物差しの量の差を比べる。今の V6 の「旧 float との回帰」は、§4.6 の 6 と段 ①〜④ の V0〜V3 で代える。
- **case/48 の組み方**:
  - 格子 `_conv614/c48/new64/m.h5`、初期場 `run_0025` の res_48000 (`restart_field.py`、7 量)。設定は `run_0048` (`run_0954` の系列、B 300 K、生産 SST、`blockDPLUR 1`・cfl 2)。
  - 腕: float ×2、FP64 ×2 (再実行の揺れを両方で測る)。各 48000 step (再開直後の残差の跳ねが数千 step で戻る既知の挙動を超える長さ。起点が float で収束した場なので、FP64 が離れていくかを見る)。
  - 物差し: `cooled_plate_eval.py --json --closure --series` の Cf・q_w・δ*・θ (x 0.3/0.6/0.9)、CD・HF、閉合。系列を `check_quasisteady` に渡す。`check_convergence` は記録 (既知の床で PASS しない)。
  - 判定案: 各量の float と FP64 の差 (2 本ずつの平均の差) が、(i) 0.05 % 以下 (θ・δ*)・0.1 % 以下 (Cf・q_w・HF) かつ再実行の揺れの 10 倍以下なら「一致」。(ii) 0.5 % 以上なら「ずれる」。その間は「判別不能」。両腕とも系列が STEADY でなければ判定不能。
- **解釈の案**: case/48 が「一致」なら、case/45 のずれは冷却壁の境界層そのものではなく、軸対称・ライン陰解法・TP・収縮部のいずれかに絞る。「ずれる」なら、冷却壁の境界層 (SST・低 Re 壁) が float で別の所に止まる。
- 残りの 4 ケース (05・09・44・16) は、case/48 と §6.20 の結果を見てから組む。

## 5. 読んでよいファイル

- `plans/active/architecture-float-state-double-geometry.md` (§5.1、§6 の V4〜V6、§6.11、§6.16〜§6.20)
- `plans/active/axisymmetric-freestream-hoop-gauge.md` (§4.10)
- `case/48.flat_plate_cooled_m4/README.md`
- `case/45.isobutane_m6_d155/fl1_judge.py`・`v4_judge.py`・`v6ab_judge.py`
