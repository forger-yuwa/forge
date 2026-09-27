# codex 諮問 (diagnose): graetz-plan-disposition

- **brief**: [`notes/reviews/briefs/2026-09-27-graetz-plan-disposition.md`](../../notes/reviews/briefs/2026-09-27-graetz-plan-disposition.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-graetz.md`](../../plans/active/boundary-cht-axisymmetric-graetz.md)
- **date**: 2026-09-27
- **commit**: `6d32c645` (feature/cht-axisym-graetz)
- **codex**: effort `high`, 4.4 min, rc=0
- **結論**: **`cfl_pseudo=2`の本番固定を保留し、FP64・`N_r=64`・共役ΔT=10 Kで、擬似CFLだけを2対1に変える起動A/Bを先に行ってください。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**元レビューの6指摘は全件採用。ただし、M2・M3・M4の処置案は修正が必要です。Critical はありません。**

| 指摘・重大度 | 採否と根拠・対案 |
|---|---|
| **M1・Major：入口温度検査** | **採用。** 加熱入口面の検査と、加熱開始断面の参考比較を分ける案でよい。ただし「対照の `x=0` は入口BCの非一様だけを見る」は言い過ぎです。対照にも散逸・膨張の効果が残ります。現行 [eval_graetz.py:168](/home/sano/work/forge-cht/case/63.graetz_cht/eval_graetz.py:168) は対照断面の温度幅しか検査していません。**実入口の温度幅を追加し、300 Kからの平均偏差も記録**してください。一様な温度ずれは温度幅では検出できません。 |
| **M2・Major：格子間差と固体層数** | **指摘は採用、処置案は要再検証。** モデル差を含む解析解誤差から格子間差へ変更するのは妥当。ただし固体16層を加熱runだけ追加すると、[差し引き式:122](/home/sano/work/forge-cht/case/63.graetz_cht/eval_graetz.py:122) の対照は8層のままになります。これは「16層の差し引き `Nu`」ではありません。**16層のΔT=0/10 Kを対で追加**し、8層の対と比較してください。0.05%は事前登録する感度許容として採用できますが、反復不確かさを含めて判定する必要があります。格子間差がその不確かさ以下なら、厳密な大小関係を認定せず判定保留とします。 |
| **M3・Major：欠損・準定常** | **指摘は採用、処置案は要修正。** 全節点の欠損・重複拒否は必要です。現行 [eval_graetz.py:85](/home/sano/work/forge-cht/case/63.graetz_cht/eval_graetz.py:85) はダンプ内座標を読むだけで、メッシュとの完全性照合がありません。**対照の `q0`・`Tw0−Tb0` も全対象節点で検査**してください。また、下記の反例から、末尾変動0.03%を反復誤差0.03%と扱うことは却下します。 |
| **M4・Major：G-if/G-cons** | **指摘は採用、`q_floor`案は却下して修正。** [check_cht_balance.py:6](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:6) は `ε ≤ min(tol_abs, tol_rel × max(Σ|Qf|, q_floor))`。提案値では小熱量の対照に対する実効許容は最小 **8.4e−8 W/rad**となり、意図した加熱尺度の許容 **8.4e−5 W/rad**より1000倍厳しくなります。**`q_floor=0.084 W/rad`、`tol_rel=1e−3`、`tol_abs=8.4e−5 W/rad`**を推奨します。G-ifの提示値は界面残差の登録閾値として採用可能ですが、`Nu`精度の保証値ではありません。 |
| **m5・Minor：固体外温** | **採用。** 現行 [make_run.py:119](/home/sano/work/forge-cht/case/63.graetz_cht/make_run.py:119) には `ROBIN/TC` の一致検査と実値記録が既にあります。実装面の指摘は対応済みです。 |
| **m6・Minor：主張の範囲** | **採用。** [plan:89](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:89) の「連成が入れた誤差」は未訂正です。提案どおり、**登録したnode・FP64条件での差し引き局所 `Nu`**に限定し、UWT近似を含む差と明記してください。 |

**M3の追加根拠：小さい変動幅は、残った過渡の上限ではありません。**

`Nu(s)=4[1+0.01 exp(−s/10⁶)]` を0～60000 step、2000 step間隔で与え、現行 `check_quasisteady.py` の `classify_series` を登録引数で実行しました。

- 判定：**`STEADY`**
- 末尾半分の変動：**0.028409%**
- 既知の漸近値4に対する最終値の誤差：**0.941765%**
- ツール自身も漸近値4、最終値との差 **−0.933%** を表示

根拠は [check_quasisteady.py:280](/home/sano/work/forge-cht/solver_density_cuda/tools/check_quasisteady.py:280) と [同:299](/home/sano/work/forge-cht/solver_density_cuda/tools/check_quasisteady.py:299)。**0.03%は変動ゲートとして残し、反復不確かさとは分ける**べきです。単調系列では漸近差も判定対象とし、残存過渡を評価できなければ0.05%の感度比較を確定しないでください。

さらに、[CHTテンプレート:24](/home/sano/work/forge-cht/case/63.graetz_cht/template/solverConfig_cht.yaml:24) の保存間隔2000は、[同:41](/home/sano/work/forge-cht/case/63.graetz_cht/template/solverConfig_cht.yaml:41) の連成間隔50の整数倍です。同位相の出力だけでは更新周期内の振動を見逃します。末尾の診断区間では更新間の位相も採取してください。

結論: **`cfl_pseudo=2`の本番固定を保留し、FP64・`N_r=64`・共役ΔT=10 Kで、擬似CFLだけを2対1に変える起動A/Bを先に行ってください。**

第 1 仮説: **主な起動障害は等温壁を含む離散系の擬似時間反復の不安定性で、CFL低減が有効である。**　確度: **中**

  根拠: 次のrunで、メッシュ座標・初期保存量・入口プロファイルのハッシュ一致を確認しました。

| run（`case/63.graetz_cht/`配下） | 確認した保存場の数値 |
|---|---|
| `run_0003_smoke_iso0_r16` | step 100：`min P=−650349 Pa`、`min ro=−7.5534`、`max|Uy|=161.003 m/s` |
| `run_0004_smoke_walls_iso_r16` | step 14：`max|Uy|=0.665 m/s`、step 16：`min P=−59488 Pa`。全壁等温でも破綻 |
| `run_0004_smoke_walls_adiab_r16` | step 200：`min P=101325 Pa`、`max|Uy|=0.001260 m/s` |
| `run_0005_smoke_iso0_cfl2_r16` | step 400：`min P=101325 Pa`、`max|Uy|=0.001258 m/s` |
| `run_0005_smoke_iso0_cfl1_r16` | step 400：`min P=101325 Pa`、`max|Uy|=0.000987 m/s` |

  これは保存時点の起動診断値であり、定常値ではありません。`check_convergence.py` を全履歴に実行した結果、**上記すべて `NOT CONVERGED`**。加熱した `run_0006_smoke_iso10_cfl2_r16` もstep 399で **`NOT CONVERGED`**でした。

  反証条件: 本番構成でCFL 1でも同様に破綻し、2→1で増幅が抑制されない場合、「CFL低減だけで起動できる」という仮説を棄却します。なお、現在の証拠だけではエネルギー行の処理まで真因を特定できません。

第 2 仮説: **FP32の丸め・幾何精度が不安定性を誘発または増幅している。未確認。** 全スモークの保存場はfloat32で、FP64対照がありません。

第 3 仮説: **細格子の共役更新開始後に別の不安定性が現れる。未確認。** 現スモークは非連成で、テンプレートの`warmup=5000`後を検証していません。

判別 A/B: **A=`cfl_pseudo: 2`、B=`1`**。同一FP64バイナリ、`N_r=64`、共役ΔT=10 K、同一IC、`implicitRelax=0.7`、`Df_scale=5`、`interval=50`を固定。**各10000 step**とし、連成開始後も約100更新を観測します。全残差列、`min ro/P/T`、壁近傍の摂動増幅、毎更新の界面残差を確認します。

- **Aだけ破綻、Bは破綻しない** → CFL 2の共通採用を却下。CFL低減策は支持。
- **両方とも破綻しない** → この起動区間でのCFL 2不適合を除外し、**3格子共通の起動候補として2を採用**。収束・精度は登録ゲートで別に判定。
- **両方破綻、またはBだけ破綻** → 単純なCFL上限説明では足りず、固定案を保留。

同じCFLは格子比較の必須条件ではありません。局所刻みは格子ごとに変わります（[solver-settings.md:43](/home/sano/work/forge-cht/procedures/solver-settings.md:43)）。今回は比較条件を揃えるため共通値を狙い、未収束場同士の差から解のCFL非依存性を主張しないでください。

やらない方がよいこと: **400 stepの無破綻を「解決」として全本番runへ展開すること、弱形式壁へ切り替えること、0.03%の変動幅を誤差上限として0.05%の感度差を認定すること。**

呼び出し側の前提への異議: **熱的な断熱／等温の継ぎ目が必須原因という説明は後退しましたが、「強制形そのものが真因」は未確定です。** 全壁等温でも境界IDの分割と格子集中は残ります。また、同じ絶対閾値を指定しても、G-ifの相対条件は各runの荷重で規格化されるため、対照が加熱尺度だけで判定されるわけではありません（[check_cht_interface.py:14](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:14)）。

不足情報: 各スモークの実行バイナリの識別・ビルド精度・鮮度記録、本番FP64での細格子・連成開始後の履歴。`run_0003`の保存場は0/100/200 stepなので、「step 14から指定内部節点で発生」は当該runから独立確認できません。**ファイル変更なし、plan未反映。反映先は対象planの§4.6・§5.1・§6・§6.1です。**
