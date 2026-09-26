# codex 諮問 (diagnose): gradient-scalar-lsq-s2-result

- **brief**: [`notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-s2-result.md`](../../notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-s2-result.md)
- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **date**: 2026-09-26
- **commit**: `f99f236d` (feature/sern-design)
- **codex**: effort `high`, 2.4 min, rc=0
- **結論**: **S2 は「条件付き（未合格・#5rでの審査は可）」とし、次は case/48 の gg で `space.slauWallNormalChi: 0` だけを変更した追加1本を回す。**
- **extra**: `case/48.flat_plate_cooled_m4/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0。実測値は README・ブリーフの記録であり、AWS の CSV・VERDICT 原本は未確認）

| 重大度 | 対象 | 採否・根拠と対案 |
|---|---|---|
| **Major** | 「両双子で外れたので、現行バイナリの正常な床」 | **要再検証**。gg 側にも原因があることは支持するが、共通の回帰不具合も同じ観測を作る。`case/40` の床約56倍を正常とする根拠にはならない（[README:112](/home/sano/work/forge-sern-design/case/40.nozzle_design_tool/README.md:112)）。対案は下記の `slauWallNormalChi` 一因子試験。 |
| **Major** | 現行 gg 初段を床にして既存結果を PASS に変更 | **現時点では却下**。現行規則は起点床の1.5倍以内（[plan:123](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:123)）。変更するなら「現行 gg に対する追加悪化の検査」として別判定にし、旧ゲート不成立を残す。参照 run・末尾窓・倍率を固定し、gg 自身の妥当性を独立に確認する。通常の収束 PASS の代用にはしない。 |
| **Major** | FCT のノイズ上限超過を離散化差として閉じる | **要再検証**。「同一設定反復のノイズ上限が作用素変更の許容差として不適切」は**採用**。元の0486は SFR 0、今回の0954/0955は SFR 2（[README:393](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:393)）。ただし「0476の定常差を上限にする」は**却下**。別の物理時刻・時間積分の差を拘束しない。旧上限超過を保存し、FCT の保存収支・有界性を独立に検証する。`ΣY≈1` は受動種の保存証明ではない（[solver-settings:174](/home/sano/work/forge-sern-design/procedures/solver-settings.md:174)）。 |
| **Major** | case/16 の既存ゲートを満たしたという整理 | **却下**。0476系は計画の「24000＋継続24000」に対して記録は24000のみ。0482系は計画が **vs 0482** なのに主表は **lsq−gg**（[plan:130](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:130)、[README:393](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:393)）。対案は不足区間・比較対象を明記し、未達のまま審査する。 |
| **Minor** | `slauWallNormalChi` の設定根拠 | **文書不整合の指摘を採用**。指定文書は既定0・case/48も0と記載（[recommended-settings:70](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:70)）。一方、`f99f236d` のコードでも node＋壁Dirichlet＋SLAUなら自動1を確認した（[solverConfig.cpp:700](/home/sano/work/forge-sern-design/solver_density_cuda/input/solverConfig.cpp:700)）。対案は実効値を provenance で固定し、文書を同期する。 |

結論: **S2 は「条件付き（未合格・#5rでの審査は可）」とし、次は case/48 の gg で `space.slauWallNormalChi: 0` だけを変更した追加1本を回す。**

第 1 仮説: **case/48 の gg が旧床に戻らない主因は `slauWallNormalChi` の実効0→1であり、0に戻せば元の床ゲートを満たす。** 確度: **中**  
  根拠: [case/48 README:55](/home/sano/work/forge-sern-design/case/48.flat_plate_cooled_m4/README.md:55) の記録では、`run_0952_sglsq_s2_gg_ext` は旧床比 `rms_ro=1.85`、`rms_roUy=9.1`、判定は `NOT BACK ON FLOOR`。このフラグは壁隣接面の質量流束に使う `chi` を実際に変更する（[convectiveFlux_slau_d.inc.cuh:544](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:544)）。単なる設定表示の差ではない。  
  反証条件: **同じ初期場・現行バイナリで0に戻し、量が静定しても旧床の1.5倍を超える列が残る。** この場合、「この変更だけで床移動を説明できる」は棄却する。

第 2 仮説: FCT smoke の差には、SFR 2での意図した面再構成変更が寄与する。確度: **中**。受動種勾配が `Pface_out` に入る経路は確認済み（[同ファイル:592](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:592)）。ただし差の全量が正常な離散化差か、`nSub=15`で十分かは未確認。

第 3 仮説: case/48 の lsq 側の追加床上昇は SST の勾配依存項による。確度: **低・未確認**。F1の1 step診断と物理量の小差だけでは、減衰機構も「結果に効かない」も確定できない。

判別 A/B: **既存の chi=1 の gg 系列を対照に、chi=0を追加1本。**

- 出力先案: `case/48.flat_plate_cooled_m4/run_0956_sglsq_s2_gg_chi0/`（AWS側で番号未使用を確認）。
- `run_0950_sglsq_s2_gg` の**開始時入力**を複製し、同じバイナリ・GPU・BLOCKSIZE・BC・ICで実施。変更は `space.slauWallNormalChi: 0` のみ。長さは対照の合計と同じ **48000 step**。
- 事前固定する判定: `check_floor_ratio --start run_0005_B_tw300_y3 --factor 1.5 --tail 0.2` の全残差列、`check_convergence` の NaN/Inf・RISING、既存の物理ゲートと系列の STEADY。
- **A: 全条件を満たして旧床へ戻る** → case/48では「chi変更が無関係」を退け、床移動の説明として採用する。他ケースへは一般化しない。
- **B: 物理量がSTEADYでも床ゲートを外れる** → 第1仮説の単独説明を棄却する。まだ過渡なら判定不能とし、合格を探す自動延長はしない。

やらない方がよいこと: **合格する参照床への差し替え、FCT差への定常差上限の流用、今回のAだけでS2全体をPASSにすること。** FCTの保存・有界性だけを満たしても、精度差の審査まで閉じない。

呼び出し側の前提への異議:  
「4ケース共通」は不正確で、case/39は両方とも床ゲートPASS。また、case/48の「約1.4倍」は参照窓依存であり、ブリーフの延長末尾値同士では `rms_ro≈1.49倍`、`rms_roe≈1.54倍`になる。これらが判定用平均と同じ定義か確認が必要。`nodeAxisDirichlet` のキー削除も、旧バイナリでの実効動作を確認せず床移動の原因とは扱わない。

不足情報: AWSの対象 run の `residual_history.csv`、判定区間付きの `CONVERGENCE_VERDICT.txt`・準定常判定原本、実効設定記録、FCTの保存収支判定、0482との直接onset比較。ローカルには対象 `run_095*_sglsq*` がなく、今回は実測を再計算していない。

**plan 未反映**。呼び出し側で §6 に旧判定と追加試験規則、§5.1 #5/#5r に未達項目を残すこと。
