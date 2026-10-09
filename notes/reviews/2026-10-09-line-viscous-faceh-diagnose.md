# codex 諮問 (diagnose): line-viscous-faceh

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-faceh.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-faceh.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `00cc0831` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.6 min, rc=0
- **結論**: **§6.15を限定付きの支持として記録し、値2・3の探索を終了して、既定の面エンタルピー精度を維持した本線の「方向別dt＋上限＋point仕上げ」の総壁時計評価へ戻る。**
- **extra**: `case/45.isobutane_m6_d155/_jprobe/fh_judge.json`, `case/45.isobutane_m6_d155/_jprobe/s0p7hb_compare.txt`, `plans/active/time_integration-line-implicit-speed.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0、Major 3）

| 論点 | 採否 | 根拠と対案 |
|---|---|---|
| ①「精度依存を支持、H-c の説明としても支持」 | **採用** | 保存された集計では、両側とも差分再現・再評価ノイズのゲートを通過。J_t の変化は **0.96576**、近似作用素との相対差は **0.97476 → 0.38926**、比は **0.39934**。§6.15 の分岐を満たす。[判定値](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_jprobe/fh_judge.json:5)。対象は **S0・p7・ライン2183・自由なエネルギー行・ε = 1e-6** に限定する。 |
| 「h の凍結で ṁ∂h が抜けた。近似作用素の欠陥ではなかった」 | **Major・要再検証** | 切替は温度入力だけでなく、組成・係数・演算・戻り値の精度も変える。[float係数生成](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/thermo_d.cuh:133)、[入力変換](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/thermo_d.cuh:1079)。面ごとの Δh と流束微分の内訳は未測定で、**凍結と相殺欠落という機構の確定には足りない**。Bでも相対差 **38.9%** が残る。「約36倍の大きさの隔たりには面エンタルピー評価精度が強く寄与した」と記録する。 |
| ②「収束した S0 で残差が 0.55% 変化し、残差床の可能性」 | **Major・前提を訂正して課題登録を採用** | S0 の出所 `case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/` は、保存された判定で **NOT CONVERGED (stalled/plateau)**。区間は `main`、連結300000行で、`rms_roOmega` は RISING。[判定記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/gates_aws.json:32)。また **0.005456 は ‖R_B−R_A‖/‖R_A‖** であり、残差ノルムが0.55%増減したという意味ではない。[集計式](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jprobe.py:253)。現時点では「固定状態の残差評価が変わる」まで。 |
| B の有限振幅誤差 0.232、実補正の予測差 0.406 | **Major・却下** | 両側の `pp` に旧 float 経路の `run_0317_jp_s0p7_pp` を指定している。[指定箇所](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jprobe_fh.sh:42)。Bでは **R_A(Q+0.7p)−R_B(Q)** を比較してしまう。[計算箇所](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/jprobe.py:175)。Bの `finite_vs_linear` と `op*_secant_rel` は無効として掲載を外す。**§6.15 の正式判定はこれらを参照しないため維持できる**。今回は追加評価せず区切る。 |
| ③ 値2・3を保留し、本線へ戻る | **採用** | [§6.15](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:381)の停止条件どおり。面エンタルピーのdouble化による発散回避・長期収束・速度改善は測っていない。既定値と本線の残差評価精度を維持する。 |
| ④ plan を閉じ、値2・3を診断として残す | **条件付き採用** | コード残置は[既存の処置](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:241)と整合する。ただし[§5.1](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:94)の **U0・FP32未実施**を完了扱いにしない。resultレビューで各項目の処置を確定した後、現役の「本線不採用・診断として残置」という判断文書として **accepted** に移すことを推奨する。 |

結論: **§6.15を限定付きの支持として記録し、値2・3の探索を終了して、既定の面エンタルピー精度を維持した本線の「方向別dt＋上限＋point仕上げ」の総壁時計評価へ戻る。**

第 1 仮説: **面エンタルピーのfloat評価に含まれる量子化が、今回のエネルギー行の大きな方向微分と近似作用素との差に強く寄与した。** 確度: **高。ただし詳細機構は中。**

根拠: [切替コード](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:399)は面エンタルピー評価を変更し、固定した方向・作用素との比較で ‖J_t p‖ が **148.862 → 5.27835（約1/28.2）**、相対差が **0.97476 → 0.38926** になった。差分幅の再現性は、量子化された評価経路の局所応答を測れていることを示しても、連続な熱力学モデルの微分を保証しない。

反証条件: 元データの監査で、入力状態・方向・実効設定が共通でない、または保存集計が再現せず、有効な差分で J_t の変化が10%未満だった場合。**「ṁ∂h が消えた」という狭い説明は、面ごとの応答を測っていないため未確定のまま残す。**

第 2 仮説: **残る38.9%の差には、近似作用素と再構成を含む実残差とのモデル差や、他の混合精度経路が寄与する。** 確度: 中。今回スカラー対角を引いても **0.38926 → 0.38874** なので、この対象ではそのスカラーだけを主要因とする説明は支持されない。

第 3 仮説: **面エンタルピー精度がpoint仕上げの残差停滞にも寄与する。** 確度: 低。固定状態の残差差だけでは、反復の到達限界との因果を示せない。

判別 A/B: **今回追加では回さず、別課題を再開するときの一組として登録する。**

- 登録先は `plans/active/time_integration-implicit-thermal-jacobian.md` **§5.1** の「point仕上げにおける残差評価精度の監査」。速度planから参照し、今回の粘性Jacobian planには移管先を残す。
- 同じpoint仕上げの保存状態・同じバイナリ・同じ設定から、A＝現行float、B＝`FORGE_DIAG_FACE_H_DOUBLE=1` のみ変更。**各2000 step、全残差を毎step**、開始・1000・1500・2000 stepの状態を保存する。
- 末尾500 stepの全残差の水準・傾きに加え、保存状態を**共通の残差評価精度でも再評価**する。表示する残差の定義が違うだけの改善を分離する。
- **Bだけで共通評価のエネルギー残差が10%以上低下し、再評価ノイズを十分上回るなら**、この期間の停滞への寄与を支持する。**自方式の表示残差だけが変わり、共通評価では10%未満なら**、この期間の主要因説を支持しない。まだ減衰中なら残差床の判定は不能とし、自動延長しない。収束・定常性は別途各判定ツールで扱う。

やらない方がよいこと:

- 今回の支持を、§6.13の正式判定「S0・S1とも(d)判別不能」の上書きに使う。
- 面エンタルピーを直ちにdoubleへ既定化する、値2・3を本線へ戻す、熱伝導Kの列削除を続ける。
- 「datum込みだからhが大きい」を実測なしに原因へ採用する。係数にはsensible datumが反映されており、対象面でのhと微分項の大きさは未測定。
- 不合格・未実施項目を消してplanを成功完了にする。

呼び出し側の前提への異議: **S0は収束済みではない。精度依存の支持は、近似作用素の健全性や発散原因の確定ではない。Bの有限振幅指標には評価経路の混在がある。**

planを閉じる際は、§6.15の結果、上記の留保、§6.13の(d)、各試験の不合格・判定不能・未実施の処置、診断コードの残置範囲、再開条件を残す。U0は既定経路の回帰確認として扱い、FP32など延期する検証は適用範囲の制限と移管先を明記する。§5.1に残る「実装経路は正しい」「1 ulp」等の古い要約も、本文の訂正に合わせて更新してからresultレビューへ渡す。

不足情報: ローカルには `fh_judge.json` と比較TXTがある一方、元のNPZ、ライン行列、対象run、実効設定、今回の比較JSONはない。**集計の閾値判定とコードは確認したが、元配列からの独立再計算・AWSバイナリの同一性確認はできていない。** `forge`は起動せず、ファイルも変更していない。**plan未反映。**
