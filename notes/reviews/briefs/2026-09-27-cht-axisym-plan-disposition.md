# 諮問ブリーフ: 軸対称 fem2d の plan レビュー (Major 5) を全件採用して実装に進んでよいか (2026-09-27)

AGENTS.md 条件 **5** (codex の Major の採否)。

## 読んでよいファイル
- レビュー: `notes/reviews/2026-09-27-boundary-cht-axisymmetric-fem2d-plan.md`
- plan: `plans/active/boundary-cht-axisymmetric-fem2d.md` (全文。§4.4〜§4.4d、§5.1、§6 が反映先)
- 必要なら `case/52.conjugate_slab/README.md` (一様 IC 単一区間で PASS した `run_0002` の記録)

## 採否 (全件採用・却下 0)
- M1: 丸め床による例外を撤去。V-ax2 は**一様 IC (325 K) からの単一区間**で `check_convergence` PASS を要求。
  根拠: case/52 の `run_0002` (一様 IC から) は PASS、継続 run は下げ幅が取れず NOT CONVERGED だった。
- M2: V-ax2 に全節点 `iface_q_eff`/`_raw` の対解析値 ≤0.5 %、V-ax2b に恒等式 $q_{\rm eff}A^r=Q_f$ と分母を外した模擬出力の負例。
- M3: 再開状態に `geometry`・`load_unit`・`state_contract`、軸対称で属性の無い旧状態は拒否、平面の旧状態は受理 (§4.4c、V-ax4)。
- M4: `check_cht_balance.py` の r 重み化と状態照合、`cht_loop.py` の軸対称 fem2d 拒否 (§4.4d、V-ax4)。
- M5: V-ax2b の格子を事前固定 ($N_r$ 8/16/32 一様 + 32 非一様、対角方向固定)、判定は $N_r$=32 で ≤0.5 % かつ倍増ごとに 1/2.5 以下。
- m6: 保証範囲 (定常・node・method 0・fem2d・FP64 で検証済み) を明記、FP32 は受理するがログで未検証と出す。

## 諮りたいこと (推奨を 1 つに)
1. 採否に問題はあるか (特に M1 の「一様 IC 単一区間で PASS」が正当な代替か、M5 の判定 1/2.5 が codex の予測列 0.94/0.30/0.088 % と整合するか)。
2. `status: in_progress` にして §5.1 #3 (固体の r 重み) から実装に入ってよいか。
