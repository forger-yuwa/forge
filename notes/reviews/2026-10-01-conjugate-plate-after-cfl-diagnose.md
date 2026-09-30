# codex 諮問 (diagnose): conjugate-plate-after-cfl

- **brief**: [`notes/reviews/briefs/2026-10-01-conjugate-plate-after-cfl.md`](../../notes/reviews/briefs/2026-10-01-conjugate-plate-after-cfl.md)
- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **date**: 2026-10-01
- **commit**: `1f6e7ed9` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 3.6 min, rc=0
- **結論**: **本番 6 本を保留し、CFL 0.5 の同一起点から `space.limiter: 2 / 0` だけを変える短い A/B を行い、残差の空間分布も同時採取する。**
- **extra**: `case/65.conjugate_flat_plate/ab_judge_cfl.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 判断 | 根拠と対案 |
|---|---|---|
| Major | **「後縁に残差が局在し、刻みに依らない」は要再検証** | 保存判定では下流 `P` の振幅比が 0.841、全体残差比が 0.428〜0.510。しかし圧力変動の最大位置は残差の発生位置を保証しない。CFL 2 点だけで刻み非依存とも言えない。[判定記録:10](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/AB_JUDGE_CFL_run0015_0016.txt:10)。対案は、残差二乗和への領域別寄与を測ること。 |
| Major | **「前縁の変動減衰」は限定採用、「G-if ① も解消」は要再検証** | 保存判定の前縁 `q` 標準偏差は `0.1139 → 2.285e-6`。書けるのは「固定壁温・当該観測区間で大幅に減衰」。G-if は固体作用素・荷重と流体熱量の**差**であり、時間変動が小さくても平均的な不釣合いは残り得る。[conjugateWall.cpp:946](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:946)。対案は同一節点・評価時点の界面残差を再構成すること。 |
| Major | **本番 6 本の CFL 0.5 再走は現時点では却下** | `run_0016` の保存された判定は **`NOT CONVERGED (stalled/plateau)`**、区間は再開後 step 0〜47999、能動全列の低下は 0.3〜0.4 桁。[CONVERGENCE_CHECK.txt:1](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/run_0016_cfl_c1_n64_cfl05/CONVERGENCE_CHECK.txt:1)。対案は C1 n64 の診断を先行し、6 本全合格の条件を維持すること。 |
| Minor | **`limiter_T` の変動ゼロを物理的な証拠にすることは却下** | `reconT` の既定は 0。この場合 `limiter_T` は計算されず 0 のまま。[variables.hpp:247](/home/sano/work/forge-cht/solver_density_cuda/variables.hpp:247)。対案は実際に使う `ro/Ux/Uy/Uz/P` の係数を評価すること。 |

結論: **本番 6 本を保留し、CFL 0.5 の同一起点から `space.limiter: 2 / 0` だけを変える短い A/B を行い、残差の空間分布も同時採取する。**

第 1 仮説: **リミッタを含む再構成作用素が、CFL 0.5 に残った停滞の維持に関与している。** 確度: **中**  
  根拠: `case/65.conjugate_flat_plate/run_0016_cfl_c1_n64_cfl05/` の保存判定では、下流 `P` の変動が残り、登録節点で `limiter_Uy/P` の時間標準偏差最大が `0.4095/0.2838`。[判定記録:19](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/AB_JUDGE_CFL_run0015_0016.txt:19)。ただし、係数変動と残差の位置的対応は未確認で、**「係数の切り替わりが真因」とまでは言えない**。  
  反証条件: 無制限再構成でも、同じ場所の残差と下流 `P` 変動が対照と同程度に持続するなら、「リミッタの作動がこの停滞に必要」という仮説を棄却する。

第 2 仮説: 後縁の no-slip／slip 接合部の境界処理が停滞を維持し、リミッタ変動はその応答である。確度: **低、未確認**。  
第 3 仮説: CFL 0.5・DPLUR 4 sweep の更新作用素にも反復振動が残っている。確度: **低、未除外**。`nStepInner` は固定残差に対する線形 sweep 数である。[main.cpp:1509](/home/sano/work/forge-cht/solver_density_cuda/main.cpp:1509)。

判別 A/B:

- **起点・変更点**: `run_0016` の最終 FP64 場と固定壁温を共通にして、新規 run を 2 本作る。A は `space.limiter: 2`、B は **`space.limiter: 0`**。`convMethod: 1`、`cfl_pseudo: 0.5`、他の数値設定・BC は共通。
- **意味**: B は無制限の 2 次再構成。係数凍結でも 1 次化でもない。既存コードは係数を 1 に充填して戻る。[limiter_d.cu:450](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/limiter_d.cu:450)。係数凍結は未実装で、既存 plan でも導入を見送っている。[limiter-config-simplify.md:67](/home/sano/work/forge-cht/plans/active/limiter-config-simplify.md:67)。
- **長さ**: 各 **12000 step、4000 step × 3 区間**。これは診断用で、収束保証の長さではない。
- **採取**: 両側とも毎 step の既存系列に、全域の `res_ro/res_roUx/res_roUy/res_roe` と使用中の `limiter_*` の空間統計を追加する。ゴーストを除き、各保存量について「後縁帯・その外側」が占める **時間積算残差二乗和の割合**を出す。全域 RMS は CSV と照合する。ソルバの RMS は `sqrt(sum(res²)/nCells)` で、体積で割った残差とは別物。[residualMonitor_d.cu:137](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/residualMonitor_d.cu:137)。
- **事前判定**: A が従来の残差・下流 `P` 振幅を各 ±10% 以内で再現することを前提にする。  
  → **B の最後の 2 区間で、能動全残差と下流 `P` 振幅が A の 1/10 以下へ減衰**すれば、リミッタを含む再構成への依存を支持する。ただし時間変動だけの因果証明ではない。  
  → **B でも全指標が A の 1/2 以上残り、3 区間の max/min ≤1.1**なら、リミッタ単独原因説を棄却する。境界欠陥の確定にはしない。  
  中間・発散・対照再現失敗は判定不能。収束と準定常の VERDICT は別途併記する。

やらない方がよいこと: `limiter: 0` で改善しただけで本番設定に採用すること、後縁の境界を先に変更すること、G-if の窓除外や閾値緩和で完了させること。

呼び出し側の前提への異議: **「圧力が最も揺れる場所＝全体残差の発生源」「係数が揺れる＝係数が原因」は受け入れない。** 案 (a) の局在確認は採用するが、上記 A/B の共通観測として行う。一方、B1 の限定結論は維持できる。手元の `run_0011/0012` の生系列を再集計し、B の壁温時間変化は 0、下流 `P` は両側とも `check_quasisteady.classify_series` で 20000 step × 3 区間すべて **`OSCILLATING`** だった。壁温の動的更新だけが維持原因という候補は出し直さない。

不足情報: **B2′ の `residual_history.csv`、`ab_series.npz`、保存 HDF5、準定常判定が手元にない。** `check_convergence.py` の今回の実行結果は両 run とも `NO residual_history.csv`。上記 B2′ 数値は保存された判定テキストの照合までで、独立再集計・NaN/Inf 検査は未実施。これらを揃えてから A/B の対照基準と起点を確定する。run 索引は [case README](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:9)。**plan 未反映**。反映先は呼び出し側の `boundary-cht-conjugate-flat-plate.md` §4.3・§5.1 とする。
