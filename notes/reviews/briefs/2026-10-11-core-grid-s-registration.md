# 諮問ブリーフ: Grefx と Gref の反復条件の感度 s の測り方の登録

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.4 の「反復条件の感度」・§4.21〜§4.25、§5.1 #13)
- **date**: 2026-10-11
- **諮問の条件**: AGENTS.md「エスカレーション条件」1 (§6 の検証計画の登録)。§5.1 #13 (f) は担当 F。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
sed -n '145,186p' plans/active/tooling-nozzle-core-grid.md     (§4.4 反復条件の感度・§4.5 量と判定 (J1 に s が入る))
sed -n '494,596p' plans/active/tooling-nozzle-core-grid.md     (§4.21 事前登録・§4.22〜§4.25 結果とユーザ決定)
grep -n 'def cmd_judge2' -A 75 case/45.isobutane_m6_d155/cg_gref_judge.py
```

## 1. 観測事実

- §4.24・§4.25: Grefx (`run_0610_cg_grefx_ext_k1`、3086 × 160、通算 120k step) と Gref の継続 (`run_0608_cg_gref_cont_k1`、4719 × 160) は、両方とも資格を満たした (全量 STEADY、d ≤ τ/8、U ≤ τ/4)。ユーザ決定で τ は θ_r・δ_loc・Q_w 1.5 %、M 0.05、判定の量から x∈[40,50] の窓を外した。結果は A (全量で |E| + U_A + U_B + p_A + p_B ≤ τ)。
- 判定の量ごとの |E| + U + p (s 抜き) と、τ までの残り:
  - θ_r [65,75]: 0.816 (残り 0.684)。θ_r [84,94]: 0.732 (残り 0.768)。
  - δ_loc [65,75]・[84,94]: 0.271・0.239 (残り 1.23 以上)。
  - Q_w: 0.712 (残り 0.788)。
  - 出口の主流 M: 0.008 (残り 0.042)。軸の M: 0.009 %pt (残り 0.041)。η 0.1 の M: 0.010 %pt (残り 0.040)。
- Grefx の壁解像は PASS (y₁⁺ > 1 の面積 0.0 %、最大 1.21)。
- 両腕の設定は B0 (ライン陰解法 dir・lineDtDirectional 1・implicitThermalJacobian 5・lineDtDirectionalCap 50・cfl 4・implicitRelax 0.7・`mesh.axisSegmentRWeight 1`・FP64 `~/forge-fgeom7-fp64`)。
- 既存の s の実例 (G0、旧バイナリ): `run_0262` は B0 系の場から point (ライン 0・implicitThermalJacobian 0・cfl 4・緩和 0.7) で 4 万 step。θ_r(70) −0.081 %、Q_w +0.009 %。
- 両腕の最後の場: `run_0608` res_60000、`run_0610` res_60000 (AWS A、全出力 20k〜60k が残っている)。
- 費用: 単独なら Gref 41.8 ms/step・Grefx 約 28 ms/step。今は別セッションの 4〜5 本と GPU を共有していて 2〜3 倍遅い。

## 2. 期待値と出典

- §4.4: 各腕の窓の最後の場から point で 4 万 step 続け、判定の量の変化 s を測る。運用上の頑健性の確認であり、真の反復誤差の上限ではない。
- §4.5 の J1: B に s_G + s_ref を足して τ 以下。s が無ければ候補にしない。
- 前回の諮問 (§4.18): point 切替の前後の差には自然な時間変化が混ざる。基準の窓が準定常でなければ反復方式の効果と帰属しない。

## 3. 判断してほしいこと

推奨を 1 つに絞る。候補を挙げるが、縛られなくてよい。

- (a) **s の定義と窓**: 候補は、両腕の res_60000 から `restart_field.py` (ビット一致) で新しい run を作り、point (lineImplicit 0・lineDtDirectional 0・implicitThermalJacobian 0・cfl 4・緩和 0.7、キーとバイナリは固定) で 4 万 step、2500 ごと。s = |μ_point − μ_base| とし、μ_point は point の 20k〜40k の時間窓平均、μ_base は §4.24 で使った窓平均 (Gref は `run_0608` の 40k〜60k、Grefx は `run_0610` の 40k〜60k)。point 側にも資格 (STEADY・d ≤ τ/8・U ≤ τ/4) を課す。
- (b) **合否**: 候補は、全量で |E| + U_A + U_B + p_A + p_B + s_A + s_B ≤ τ なら「暫定の許容差での生産の候補」としてユーザに諮る。point 側の資格が無ければ判定不能。
- (c) point の 4 万 step が足りないとき (資格なし) に延長するか、判定不能で止めるか。
