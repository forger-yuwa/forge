# codex 諮問 (diagnose): lvc-faceh-result

- **brief**: [`notes/reviews/briefs/2026-10-10-lvc-faceh-result.md`](../../notes/reviews/briefs/2026-10-10-lvc-faceh-result.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `5e2ac3c5` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 2.8 min, rc=0
- **結論**: **棄却を切替単独の十分性に限定して記録し、次は同じ薄層モデルに対するU-Jの列別・壁拘束照合を行う。**
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcfh_judge.json`, `notes/reviews/2026-10-10-lvc-faceh-design-diagnose.md`, `plans/accepted/time_integration-line-viscous-jacobian.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 論点 | 採否 | 根拠と対案 |
|---|---|---|
| 登録の「棄却」 | **採用** | [保存判定 JSON](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcfh_judge.json) は45ゲート通過、4本とも `DIVERGED`。代表 step は A1/A2 = 122/122、B1/B2 = 121/122。[§6 分岐2](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:79)どおり、**「この条件・期間では切替だけで非有限化を回避できない」**と記録してよい。「遅らせる」の登録条件にも該当しない。生データの独立監査は未実施。 |
| 「4桁一致」「増幅率不変」から精度依存を除外する | **却下・Major 1** | [ブリーフ:18](/home/sano/work/forge-faceh/notes/reviews/briefs/2026-10-10-lvc-faceh-result.md:18)自身が step 50 の `rms_ro` を A1 = 3.900e−3、B1 = 3.898e−3 としており、4桁一致ではない。全残差の最大/開始や閾値通過 step は、増幅率の推定でもない。**「両側で早期増幅と非有限化を観測」**までに限定する。旧 `run_0311` との少数点の近さも、LAYOUT2・バイナリ差を除外する証拠にはしない。 |
| 「局所応答と、増幅を支配する有限振幅応答は別だった」 | **要再検証・Major 2** | [親 §6.16](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:400)の148.9→5.28は、S0・p7・ライン2183・自由エネルギー行・ε = 1e−6に限る測定。しかも[B側の有限振幅評価は無効](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:416)。**「局所応答の精度依存は確認されたが、今回の破綻回避には十分でなかった。両者の因果関係は未同定」**と書く。 |
| `ro` が最初の原因、ρv・ωは面エンタルピーと無関係 | **却下・Major 3** | [`main.cpp:2723`](/home/sano/work/forge-faceh/solver_density_cuda/main.cpp:2723)は検査順を `ro` から固定し、[GPU検査](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/residualMonitor_d.cu:86)も最初に見つかった変数で戻る。`ro` 表示は時間的な発生順ではない。また、エネルギーの変化は状態更新を通じて他成分へ伝わる。**「検知時に報告された変数は `ro`」**と訂正し、原因の順序は局所状態・床・補正の記録で確かめる。ログとCSVの1 step差には[表示が `iStep+1` である事情](/home/sano/work/forge-faceh/solver_density_cuda/main.cpp:2771)もある。 |
| 次はU-Jの列ごとの照合 | **採用** | [単体試験:69](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:69)は現在もD・K全体の最大値で正規化しており、小さい列の相対誤差を保証しない。まず登録済みの列別判定を完了する。ただし、未検証であること自体を実装バグの証拠にはしない。 |
| 証拠ゲート通過＝帳簿の全記録がそろう | **要再検証・Minor 1** | [`ledger_check:143`](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:143)は節点集合と呼び出し集合を別々に検査するだけ。メモリ内の人工入力で、**各呼び出し1節点・1フィールド、計123行**でも通過した。112節点×123呼び出しの網羅性は保証されない。実帳簿の `(call, tag, node, field)` の欠損・重複・値の解釈可能性を監査する。今回の帳簿が欠けていると断定するものではない。 |
| 証拠の保存範囲を縮小する | **却下** | [保存規定:96](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:96)を維持する。今回は `diagnose` であり、result段の監査は未完了。4本すべての初期場・`res_100`・`res_nan`・格子・帳簿に加え、設定・残差CSV・実行来歴を保持する。 |

結論: **棄却を切替単独の十分性に限定して記録し、次は同じ薄層モデルに対するU-Jの列別・壁拘束照合を行う。**

第 1 仮説: 薄層モデルの微分自体は整合していても、実残差・境界・分離更新を含む反復との不整合が残り、面エンタルピーの精度変更後も増幅する。確度: **中**。特定の項への帰属は未確認。

- 根拠: 今回の4本は精度切替の両側で `DIVERGED`。また、[共通関数の定義](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94)は薄層の前処理近似であり、実残差の厳密なJacobianではない。旧U-Jの合格も全行列正規化に限られる。
- 反証条件: 同じ薄層モデルとの列別照合で、差分幅を変えても残る系統的な誤り、または壁拘束の不整合が見つかれば、「微分実装は整合している」という前提を撤回する。

第 2 仮説: 小さい列、面の向き、壁拘束の処理に実装上の不整合があり、従来の判定で隠れている。確度: **低・未確認**。現時点の根拠は試験の保証範囲の不足であり、誤りを検出したわけではない。

判別 A/B: **同一状態・同一薄層流束モデルで、微分の評価方法だけを変える。A = 製品の共通関数によるD/K、B = 独立な流束実装の中心差分。時間発展は0 step、host上の200組と短いラインで完結させる。**

- 列ごとに正規化し、doubleで相対誤差 ≤ 1e−6。差分幅 h と h/2 の再現性を確認し、零列は事前固定した無次元の絶対誤差で扱う。差分自体が解像できない列は判別不能とする。
- 同じ面の両向き、異なる密度、高速流、速度固定・温度固定を含める。壁拘束は自由度を消去した系でも照合し、登録済みの零空間 ≤ 1e−12、壁温拘束 ≤ 1e−12、短いラインの解の差 ≤ 1e−10も維持する。
- **AとBが有効な全列・拘束で整合** → 第2仮説を試験範囲内で退け、第1仮説の検証へ進む。実残差との整合や安定性を証明したとは扱わない。
- **再現可能な不一致** → 第1仮説の前提を撤回し、第2仮説を優先する。ただし、その不一致が今回の破綻原因かは別に確認する。

U-J合格後は、製品経路の[物性入力・面重み・K配置](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:970)と[壁の行置換](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1144)まで照合する。そのうえで、共通の床到達前状態・実補正方向について、実残差の微小応答と有限補正応答を同じ評価経路内で測る。時間項・値3の追加スカラー・薄層項を分け、**近似作用素の差が実際の一回の更新で増幅につながるか**を確かめる。長期runや速度比較は、この確認の代用にならない。

やらない方がよいこと: 値2の追加runを先に回す、熱伝導Kの列を削る、CFL・緩和・LHS精度を同時に変える、double評価を既定化する、診断機能を消す、今回の結果を理由に証拠を削除する。

呼び出し側の前提への異議: **観測されたのは切替後にも破綻したことであり、精度依存の不在、同じ破綻機構、増幅率の同一性ではない。** 親§6.16の限定付き支持と今回の棄却は両立する。

不足情報: 対象の `case/45.isobutane_m6_d155/run_0540_lvcfh_v3_a1/`、`run_0541_lvcfh_v3_b1/`、`run_0542_lvcfh_v3_a2/`、`run_0543_lvcfh_v3_b2/` と出発 `run_0183` はローカルにない。したがって、生CSV・HDF5・帳簿・実効設定の再照合、および `check_convergence` の独立実行はしていない。確認したのは保存JSON、コード、文書、判定関数の人工入力である。run索引は [case README:191](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/README.md:191)。

ファイル変更・`forge` 起動はしていない。**plan未反映**。呼び出し側で対象planの§6.2・§6.1・§5.1、およびcase READMEの「4桁一致」を更新すること。
