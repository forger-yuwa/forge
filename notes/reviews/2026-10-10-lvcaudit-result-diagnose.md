# codex 諮問 (diagnose): lvcaudit-result

- **brief**: [`notes/reviews/briefs/2026-10-10-lvcaudit-result.md`](../../notes/reviews/briefs/2026-10-10-lvcaudit-result.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `fe878c4b` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.7 min, rc=0
- **結論**: ?
- **extra**: `case/45.isobutane_m6_d155/lvcaudit_judge.py`, `solver_density_cuda/tools/line_audit_helper.cpp`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。コードは直接確認しましたが、`run_0540`〜`0543`・`run_0550`〜`0552` の実データはこの作業ツリーにありません。以下の run 数値は plan §6.8 の記録に基づき、再計測した値ではありません。

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 登録の `FAIL` | **採用、原文保存**。C の不一致だけでなく、D の判別不能 3 件と T の不合格も残るため、C の尺度を説明できても総合を `PASS` に読み替えられない。[plan:268](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:268)。追加調査は別欄に記録する。 |
| Major | H1「C は相殺によるもので組立の誤りではない」 | **要再検証**。特に「V と音速の相殺」という説明は、`w = nz = 0` ならコードと合わない。音響固有ベクトルの z 成分が零となり、対流 `K[3][3] = S·max(−V, 0)` に簡約される。[共通関数:62](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:62)、[同:82](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:82)。まず既存記録で `w,nz` と `u·nx,v·ny` を確認する。相殺候補は **V の内積評価**であり、両者が double から外れるだけでは CUDA–host 間の差を説明したことにならない。 |
| Major | H3「D の 1 ulp 差は FMA で、数値の性質を変えない」 | **コード生成差は仮説として採用、無害という結論は却下**。記録上、通常ビルドの再実行では D が一致し、監査用だけ異なる。これは計装による差を疑う根拠だが、FMA の特定や補正への影響の小ささは未証明。さらに rhs も異なるので、dq の差を D だけに帰属できない。[判定器:369](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:369)。T の保留を維持する。 |
| Major | A の意味と次の一手 | **幾何精度の不整合は採用、破綻への因果は要再検証。候補 (a) を推奨**。LHS は座標を ST に落としてから引き、残差は FP64 ビルドで double の座標差を使う。[LHS:1021](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1021)、[残差:137](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/viscousFlux_d.cu:137)。候補 (b) は対流 Jacobian・粘性 Jacobian・sweep の演算精度まで変えるため、この原因の切り分けには広すぎる。[精度分岐:1902](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1902)。 |
| Minor | 元 plan の対象範囲 | **「LHS を対象へ追加してもらう」は不要**。元 plan §4.2a は ISP 0 の LHS が FP64 ビルドでも変わることを明記し、段③に `timeIntegration_d.cu` を含めている。[元 plan:71](/home/sano/work/forge-faceh/plans/active/architecture-float-state-double-geometry.md:71)、[同:141](/home/sano/work/forge-faceh/plans/active/architecture-float-state-double-geometry.md:141)。申し送りは対象追加ではなく、観測と検証条件の共有にする。 |
| Minor | S の「sweep の有限性も合格」 | **記述補正と検査追加を採用**。`rhs/dqnew/dqold` は存在確認されるが、S の有限性検査の対象に入っていない。[判定器:137](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:137)、[同:156](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:156)。同じ NaN は T のビット比較を通ることも確認した。各 sweep の形・有限性を明示的に検査する。今回の記録が非有限だったという指摘ではない。 |

**結論:** 候補 **(a)**、すなわち **LHS の座標差だけを double で引いてから ST に丸める診断切替**を別登録し、値 3・マスク 7 の早期非有限化を回避できるかを確かめる。

**第 1 仮説:** 座標差の量子化による LHS 係数のずれが今回の早期非有限化を支配し、その修復だけで登録期間内の非有限化を回避できる。  
**確度: 低**。精度不整合の存在はコードで確認できるが、破綻との因果は未確認。

- **根拠:** `timeIntegration_d.cu:1021` の変換順序と、`case/45.isobutane_m6_d155/run_0550_lvcaudit/` の記録。β・κ の差が 1e-3 超の面寄与は 294 件、最大 10.8 %。ただし係数は `δ/dcc` に依存し、ガード非作動時には `S²/|e·S|` となるため、**dcc の長さだけでなく e·S の誤差**として扱う。
- **反証条件:** 変更後に対象係数の幾何由来の誤差が解消したことを確認しても、変更腕が再実行を含め非有限化するなら、「この修復だけで回避できる」は棄却する。精度不整合の寄与自体をゼロとは結論しない。

**第 2 仮説:** 幾何を修復しても、薄層近似と実残差・境界・分離更新との不整合が残り、破綻する。確度: 中。共通関数自身も実残差の厳密な微分ではなく前処理近似と定義されている。[共通関数:94](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94)。

**第 3 仮説:** C・T の外れには演算順序やコード生成の違いが寄与する。確度: 中。ただし、原因の特定と「無害」の証明は別であり、現時点ではどちらも未完了。

**判別 A/B:** 同じ FP64 バイナリの診断切替で、次の **座標差の評価方法だけ**を変える。

- **A:** `ST(cc[o]) − ST(cc[i])`
- **B:** `ST(double(cc[o]) − double(cc[i]))`

x・y・z の3成分を一組として切り替え、その後の ST 演算・ガード・面の対象・物性・加算順序は固定する。`β,κ` だけでなく、同じ幾何を使うスカラー対角と熱伝導対角にも一貫して効かせる。ISP 0、値 3、マスク 7、キー 5、方向別 dt、上限なし、CFL 4、緩和 0.7、5 sweep、面エンタルピーの既定評価を固定する。

同じ `run_0183` の保存場から、**各腕2本、最大2000 step、非有限検出で終了**とする。最初の factor で、対象係数の double 参照に対する差が **1e-5 以下**へ減ることを介入の成立条件として登録する。監査ビルドの時間発展を通常経路の代用にせず、変更したバイナリの A 腕で既知の破綻を再現する。

見る量は全残差列の有限性、最初の非有限 step、保存場の有限性・ρ/P/T の正値。収束判定の VERDICT も区間付きで残すが、2000 step 有限を「収束」と呼ばない。

→ **A が2本とも非有限、B が2本とも2000 step 有限:** この条件・期間では座標差修復だけによる回避を支持。  
→ **A・B がともに2本非有限:** 修復だけで十分という第1仮説を棄却。  
→ A が再現しない、B の結果が分かれる、係数修復の成立条件を満たさない場合は判別不能。事後延長や閾値変更はしない。

**やらない方がよいこと:** C の許容を広げて元の `FAIL` を消すこと、音速の大きさを根拠に `K[3][3]` の誤差を許容すること、ISP 1 の成功から座標差を原因と断定すること、マスク 5 への退避や既定値変更を同時に行うこと。

**呼び出し側の前提への異議:** B の合格が示すのは「記録された入力から float の式を再現できた」ことであり、その式が実残差に十分整合することではない。また、本線 B0 への影響はスカラー対角だけではない。**キー 5 の熱伝導対角も同じ `δ/dcc` を使う**。[実装:1120](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1120)。

元セッションへの申し送り案:

> FP64 ビルドでも ISP 0 の block DPLUR は、座標を float にしてから差を取ります。case/45 の監査記録では壁法線ライン面の β・κ に最大10.8 %の差がありました。同じ幾何は B0 のスカラー対角とキー5の熱伝導対角にも使われます。元 plan §4.2a・段③の対象には既に含まれています。今回の発散への因果と B0 の収束挙動への影響は未確定です。LHS の座標差だけを変える診断結果を共有します。残差の固定状態での不変性と、LHS・更新履歴の変化を分けて検証してください。

**不足情報:** `run_0550`〜`0552` の `linedump`・判定 JSON・ビルド条件、および登録外の `convk_probe.cpp` が必要。特に不一致125件の `w,nz,V` と内積の各項がないため、C の相殺機構は確定できない。

ファイル変更・forge 起動は行っていません。**plan 未反映**。提案と採否の反映は、依頼どおり呼び出し側が行う前提です。
