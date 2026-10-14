# codex 諮問 (diagnose): test-p-design

- **brief**: [`notes/reviews/briefs/2026-10-13-test-p-design.md`](../../notes/reviews/briefs/2026-10-13-test-p-design.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-13
- **commit**: `edc2dce7` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.2 min, rc=0
- **結論**: **G0c・共通Q32起点・移植後FP64・`cfl_pseudo=4`・`precondEps=0.15`で、序盤の資格確認を組み込んだP0/P3の二腕試験を事前登録し、まず各20,000 stepの応答を比較する。**
- **extra**: `notes/reviews/2026-10-12-speedup-tests-diagnose.md`, `procedures/solver-settings.md`, `case/45.isobutane_m6_d155/sp_judge.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | G0c・point・P0/P3 の比較：**条件付き採用** | 主対照は **P0**。`run_0804_sp_g0c` の B0 は、生産設定を置き換える価値を見る別の基準とする。P3 が P0 より速くても B0 より遅ければ、生産の高速化にはならない。比較は「登録した停止・出力要件までの時間」に限定する。 |
| **Major** | 「20,000 stepの残差が2倍なら、速度改善なし」：**却下** | 異なる反復経路の残差比は、到達時間の大小を保証しない。既存の停止規則も収束解の同等性を保証しないと明記されている〔[plan §6.11](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:494)〕。**持続悪化を理由とする試験予算の打ち切り**として登録し、「この条件では継続対象から外した」と報告する。具体的な基準は下記。 |
| **Major** | 移植後 FP64 バイナリを両腕に使用：**採用**。β<1・多成分の資格確認を省略：**却下** | 保存済み `port_reg2_judge.json` では FP64 の3組が `A_ok: true`。一方、β=1では化学種の刻みを分離する効果を検査できない。`speciesPrecondDt=1` は拡大前の刻みを保存するが、適用対象は化学種の scalar-DPLUR に限る〔[setDT_d.cu:417](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:417)、[solver-settings.md:168](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:168)〕。**同じP0/P3試験の序盤に実効刻み・拘束行・ρY更新の確認を組み込む**。 |
| **Major** | `implicitSolvePrecision=1` の追加腕で前処理単独に帰属：**却下**。今回は追加しない案：**採用** | 精度以外にも作用素が違う。通常経路の軸対称ソースにはない `chi_hoop` が前処理経路にはある〔[timeIntegration_d.cu:1085](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1085)、[同:1398](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1398)〕。追加腕で分けられるのは通常経路の解法精度への感度まで。今回は**精度・Jacobian・擬似時間・化学種の刻みを含む設定一式の効果**を測る。 |
| **Major** | forge-25 の結果をこちらの不合格判定に流用：**却下** | ブリーフ自身が判定出力を未読としており、対象の `run_0688_pp_c`・`run_0689_pp_p3` はこの環境にも存在しない。格子・起点・バイナリ・`cfl_pseudo` が違う。**悪化を監視する動機には使うが、G0c の試験を省略する根拠や閾値の校正には使わない**。 |
| **Minor** | CFL・`precondEps` の探索：**今回は行わない** | 最初は両腕とも **`time.deltaT.cfl_pseudo: 4`、`implicitRelax: 0.7`、`precondEps: 0.15` を明記**する。定常の刻みを決めるのは `cfl_pseudo` である〔[solver-settings.md:35](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:35)〕。CFLは1水準とし、失敗後に片腕だけ1へ下げて同じ試験として扱わない。結論もこの設定に限定する。 |

結論: **G0c・共通Q32起点・移植後FP64・`cfl_pseudo=4`・`precondEps=0.15`で、序盤の資格確認を組み込んだP0/P3の二腕試験を事前登録し、まず各20,000 stepの応答を比較する。**

第 1 仮説: **P3はこの設定では縮流部の残差を持続的に悪化させ、長時間試験の継続条件を満たさない。** 確度: **低**。

  根拠: P3は擬似時間項と刻みを変える一方、物理FVSの対角・近傍結合は残す〔[timeIntegration_d.cu:1345](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1345)、[setDT_d.cu:250](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:250)〕。刻みの拡大だけでは有限回の内部反復の改善を保証しない。forge-25の悪化報告はこの予想と整合するが、独立検証できていないため補助情報にとどめる。**G0cでの悪化は未確認**。

  反証条件: 資格確認に合格し、下記の持続悪化条件に該当しないこと。これは短期悪化の予想を退けるだけで、高速化の証明ではない。

第 2 仮説: **初期過渡では残差が高くても、P3は登録した出力要件へP0より早く到達する。** 確度: **低・未確認**。段2で判別する。

判別 A/B:

- **変更する設定は `time.deltaT.lowMachPrecond: 0 / 3` だけ。** 両腕に同じバイナリ `75ee4675b0c90f41`、G0c入力、point設定、内部反復数、出力・診断条件を使う。`implicitSolvePrecision: 0` は明記して固定する。P3側では別経路へ分岐するため、これで解く精度まで同じになるとは解釈しない。
- **序盤の資格確認を同じ二腕に含める。** 初期化後の9保存量を照合し、最初の更新で β と `dt_local` の拡大がコードの式に従うことを確認する。`speciesImplicitCoupling` の実効値が1以上であること、`speciesPrecondDt: 1` により `dt_local_sp` が拡大前の刻みになることを確認する。BC再適用前の拘束成分が0であること、ρYの更新・再正規化前後の補正量、組成・ρ・P・Tの有限性と物理範囲も記録する。**確認不能・異常なら、速度試験の結果にせず資格不成立として止める。**
- **20,000 stepでの打ち切り条件を固定する。** 以下は実測で校正済みの許容値ではなく、今回の継続判断用の提案値。
  - 窓は W1＝15,000以上17,500未満、W2＝17,500以上20,000以下。
  - 全域は保存された全 `rms_*` 列。局所は固定領域 **x/r_t ∈ [−8, −2)** の全内部節点について、各保存量の `Σ|res_q|` を使う。領域・節点集合は起動前に固定し、悪化を見てから絞らない。
  - 局所統計も500 step以下の間隔で採取し、各窓5点以上を確保する。停止規則の評価間隔は既存どおり2,500 stepとする。
  - **同じ指標で**両窓とも `中央値(P3) / 中央値(P0) ≥ 2`、かつ `中央値(P3,W2) / 中央値(P3,W1) ≥ 0.9` なら、持続悪化として段2へ進まない。P0の中央値が0なら比を作らず別記し、それだけで打ち切らない。欠損・非有限は通常の比判定から分離する。
  - **A：持続悪化あり → 第1仮説を支持し、この設定の長時間試験を打ち切る。B：持続悪化なし → 第1仮説を退け、段2へ進む。** Bを「P3が優れる」とは呼ばない。
- **段2は同じ二腕を通算20万stepまで。** §6.11の停止規則、T2の両参照とのM1〜M7を適用する。両腕が出力要件を満たした場合だけ `T_reach(P3)/T_reach(P0)` を評価し、≤0.9を改善、≥1を改善なし、その間を未確定と事前登録する。片腕だけ合格なら「今回の要件を満たしたのはその腕」と報告し、速度倍率を出さない。
- B0との比較は別に示す。基準は `case/45.isobutane_m6_d155/run_0804_sp_g0c/` の **2,448.5秒**。起動・出力・診断・再開を含む時間を数え、GPU専有と計時条件を確認する。旧バイナリによる既存値なので、僅差を精密な速度差として断定しない。

やらない方がよいこと: **残差2倍から最終的な速度不改善を断定すること、θ_rがT2の値を通過した時刻を到達時間にすること、失敗を見てCFL・ε・判定窓を変更すること、`=2`をこちらでも追加すること。** β=1で確認した境界修復を、証拠なく再び主要因に戻す必要もない。

呼び出し側の前提への異議: **今回の設計で答えられるのは「このP3設定一式で、登録した出力要件までの所要が短くなるか」であり、前処理単独の因果効果や収束した定常解への加速ではない。** 保存済み [sp_judge.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/sp_judge.json) の `run_0804_sp_g0c` は、停止step 80,000について次の判定である。

```text
check_convergence: OVERALL: CHECK FAILURES ABOVE
check_quasisteady: θ_r(40/70/94), Q_w はすべて DRIFTING
判定区間: 準定常は末尾20,000 step
```

この既存結果を「収束済みの基準」と呼ばない。新しい二腕でも、停止・出力要件の合否と収束・準定常のVERDICTを分けて残す。

不足情報: 実際の `run_0713_pt_d1` のYAML、G0c入力、対象runの生CSV/HDF5、実バイナリはこの環境にない。特に `speciesImplicitCoupling`、初期組成分布、β<1での化学種更新は未確認。既存数値は保存済み判定JSONを確認した範囲であり、forgeは実行していない。**ファイル変更なし・plan未反映**。採用する内容は呼び出し側が `architecture-float-state-double-geometry.md` §5.1 #25・#26と試験Pの事前登録節へ反映する。
