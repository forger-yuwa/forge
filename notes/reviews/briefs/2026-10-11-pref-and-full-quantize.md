# 諮問ブリーフ: ρE 単独の量子化では出なかった — 次は pRef (圧力の基準値) と 9 保存量の量子化の短い試験 (2026-10-11)

- 依頼者: 主セッション (エスカレーション条件 1・4・7)
- plan: `plans/active/architecture-float-state-double-geometry.md` (§6.38・§6.40・§6.41)
- ユーザの見立て (2026-10-11): 「ふつうに残差計算に誤差が入るのではないか」。短い run で進める方針は継続。

## 0. 問い
1. §6.41 から書いてよい範囲 (下の §2) に穴はないか。
2. §3 の 2 つの腕 (pRef・9 保存量の量子化) を 1 回の事前登録で回してよいか。判定の式をどう置くか。
3. pRef の効き方について、下の §2 の見積もり (p̃ の表現の誤差は残る) は妥当か。

## 1. 観測事実
- §6.40/§6.41 (AWS A、各 2,000 step、x/r_t < −2 の内点 83,310 節点、参照の型への射影 a、既存の float − FP64 の平均の応答が 1):
  - C (FP64 ×3): a 0.2164〜0.2166。step 2,000 の温度の C の平均との差 RMS 0.002〜0.003 K・最大 0.05〜0.08 K。
  - E (FP64 + ρE だけ step 末尾に float32 に量子化 ×3): a 0.218〜0.225、e_E = +0.006 [−0.007, +0.018]。RMS 0.16〜0.19 K・最大 4.4〜5.5 K。
  - F (float ×3、§6.26 と同じバイナリ): a 1.04〜1.23、e_F = +0.91 [+0.57, +1.25]。RMS 3.7〜4.5 K・最大 119〜134 K。
  - 判定は結果 B (ρE 単独の量子化では初期の温度差の主な成分は出ない)。
- 入口〜収縮部の静圧 (`run_0487` の step 0): 入口 x/r_t < −10 で 5.4989e6〜5.4995e6 Pa (幅 約 550 Pa)、収縮部で 5.317e6〜5.4995e6 Pa。float の 5.5e6 の刻み (ulp) は 0.5 Pa。出口は約 2.2e3 Pa。
- ρE の float への丸めを圧力に換算すると 0.03〜0.14 Pa (§6.38 の訂正、場の値)。
- `space.pRef` (既定 0 = オフ、case/45 は未設定): SLAU の内部の面 (`convectiveFlux_slau_d.inc.cuh:618` `p_tilde_r = p_tilde − d_pRef`)、境界の半割面 (`convectiveFlux_boundary_d.inc.cuh:197`)、軸対称の hoop ソース (`axisymmetricSource_d.cu:54` `(P + pk − pRef)·A_closure`) に入る。もともと非直交格子の一様流の float の桁落ちの対策 (hoop と同じゲージにそろえる修正済み)。pRef は `flow_float` (float のビルドでは float32)。
- §6.35: 平均流の 5 保存量を FP64 で積んでも (残差は float のミラーから float で計算) 大幅な抑制は観測されなかった (判別不能)。

## 2. 呼び出し側の読み (未検証)
- E (状態の丸めだけ、演算は double) では差がほとんど出ず、F (演算も float) では出る。よって主な成分は float の演算 (残差の組み立て: 流束の積・足し込み・勾配・リミタ・粘性・ソース・EOS) の側にある、と読みたい。ただし E は ρE だけで、ほかの保存量 (ρ・ρu・k・ω・ρY) の丸めは試していない。
- 入口の静圧 5.5 MPa に対しセルごとの正味の圧力差は数 Pa。float で p̃·S (5.5e6·S) を作って足すと、面ごとに約 0.3 Pa·S の丸めが入り得る。pRef = 5.499e6 Pa にすると (p̃ − pRef) は数百 Pa になり、積と足し込みの桁落ちは小さくなる。ただし p̃ 自体は float の P (刻み 0.5 Pa) から作るので、その表現の誤差 (約 0.25 Pa) は残る。効果は 2〜数倍程度かもしれず、桁で効くとは限らない。

## 3. 次の試験の案 (未検証)
- 各 2,000 step、§6.40 と同じ入力・設定・環境変数・指標 (参照の型 v_s と a、温度の K の表)、AWS A で 1 本ずつ。
- **pRef の腕**: Cp = FP64 (`~/forge-fgeom11-fp64`) + `space.pRef: 5.499e6` ×3、Fp = float (`~/forge-fgeom7-f32`) + 同じ pRef ×3。比較は Fp 対 Cp (同じ作用素どうし)。Cp と既存の C の差 (FP64 で pRef が作用素を変えないこと) も記録。
- **9 保存量の量子化の腕**: G9 = FP64 + `FORGE_DIAG_QUANTIZE_STATE=1` (全 9 本) ×3。比較は G9 対 C (既存)。
- 判定の案: e_Fp = mean(a_Fp) − mean(a_Cp)、e_G9 = mean(a_G9) − mean(a_C)、Welch の区間、3 つの比較 (Fp・G9・既存 F の再掲は無し) で全体 95 %。
  - pRef: 「抑制」= e_Fp の区間の上限 ≤ 0.5 (既存の e_F 0.91 の約半分以下)、「抑制不十分」= 下限 ≥ 0.5。
  - G9: §6.40 と同じ (結果 A: 下限 ≥ 0.5、結果 B: 区間 ⊂ [−0.1, 0.1])。
- 解釈の範囲の案: pRef で抑制なら「圧力の積・足し込みの桁落ちが、入口の初期の温度差の主な通り道」(下流の θ_r は別に 30,000 step で確かめる)。G9 が結果 B なら「状態の丸め (9 保存量、step 末尾) だけでは出ない → float の演算の側」。

## 4. 読んでよいファイル
- float plan §6.38・§6.40・§6.41
- `solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh` (600〜640)、`convectiveFlux_boundary_d.inc.cuh` (190〜225)、`convectiveFlux_common_d.cuh`、`solver_density_cuda/cuda_forge/axisymmetricSource_d.cu` (25〜60)
- `case/45.isobutane_m6_d155/rq.sh`・`rq_an.py`・`rq_dT.py`
