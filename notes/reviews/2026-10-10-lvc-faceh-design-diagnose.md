# codex 諮問 (diagnose): lvc-faceh-design

- **brief**: [`notes/reviews/briefs/2026-10-10-lvc-faceh-design.md`](../../notes/reviews/briefs/2026-10-10-lvc-faceh-design.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `0be79325` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **判定器・入力ゲート・証拠保存を修正し、値3・マスク7で面エンタルピー精度だけを変える各2本・最大2000 stepのA/Bを先に行う。**
- **extra**: `plans/accepted/time_integration-line-viscous-jacobian.md`, `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| このA/BでH1とH2を判別する | **要再検証・Major** | [plan:72](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:72)の分岐で測れるのは、切替による**2000 step以内の非有限化回避**まで。Bが有限でも、近似Jacobianの問題と精度依存は併存し得る。両側が発散しても、[plan:22](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:22)の「熱伝導K側に原因が残る」には絞れない。支持は「この期間の回避を支持」、棄却は「この切替だけでは回避できない」と記録する。 |
| `FINITE`判定の実装 | **要再検証・Major** | [判定器:61](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:61)は存在する残差列だけを検査する。メモリ内の模擬入力で、**残差列なし／負の`rms_ro`だけ**のCSVが、ともに`FINITE`になった。[判定器:87](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:87)。必須の平均流・乱流・化学種残差列、全行の有限性・非負性、欠損・破損ファイルを検査し、異常は理由付き`INVALID`にする。さらに[191行](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:191)はrun分類が`INVALID`でも入力ゲートが通れば終了コード0になるため、これも直す。 |
| 200 stepごとの場と判定後の削除 | **却下・Major** | 既往の破綻は29、122〜150 stepなので、[現在の出力間隔](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:66)では初期場と非有限化後の場しか残らない可能性が高い。[後片付け規定](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:85)はそれらも消す。既存の[序盤ledger採取](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/e1b.sh:24)を両側同条件で追加し、初期場・床到達前の状態・`res_nan_*`・対応する格子をresultレビューまで保持する。選択節点のledgerだけで全域の原因を確定しない。 |
| 初期場・実効設定のゲート | **採用、補強・Minor** | [準備コード:111](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/cold_cfl.py:111)は親runの最新出力を選ぶ。現判定器には[指定の`res_100000`との事後照合](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:155)があるが、実行前にもファイル名・事前固定したハッシュを照合する。またLAYOUT2には[自動フォールバック](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2965)があるため、実効の並びと`implicitSolvePrecision=0`を確認する。FP64ビルドはLHS全体の倍精度を意味しない。 |
| 値3を主、各2本、最大2000 step | **採用** | 値3は[スカラー対角を保持](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:969)し、既存の方向微分診断ともつながる。2000 stepは早期破綻の回避を調べる期間として妥当。ただし各2本は探索的な再現確認で、破綻確率や長期安定性は保証しない。副の値2・1/2本は探索扱いに限り、今回は主の4本を先に判定する。 |

結論: **判定器・入力ゲート・証拠保存を修正し、値3・マスク7で面エンタルピー精度だけを変える各2本・最大2000 stepのA/Bを先に行う。**

第1仮説: 面エンタルピー評価の混合精度が、値3・マスク7の早期増幅に寄与している。確度: **中**。非有限化を回避できるかは未確認。

- 根拠: [切替実装:401](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:401)はTP多成分の面エンタルピー評価を変更する。親planの[§6.16記録](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:398)では、S0・p7・ライン2183の自由エネルギー行で‖Jₜp‖が148.9→5.28、近似作用素との相対差が0.975→0.389。ただし、これは局所応答の記録であり、破綻回避の測定ではない。
- 反証条件: ゲートを満たしたBが2本とも2000 step以内に非有限化すれば、**「この切替だけで非有限化を回避できる」という強い予測**を棄却する。精度依存そのものや、元の破綻機構への寄与までは否定しない。

第2仮説: 薄層近似・壁拘束・ライン外結合・SST更新などとの不整合が、精度変更後も増幅を残す。確度: **中**。熱伝導K単独への帰属は未確認。実装自身も[厳密Jacobianではなく前処理近似](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94)としている。

第3仮説: 精度変更は非有限化を2000 stepより後へ遅らせるだけである。確度: **低・未確認**。Bが2000 step有限でも除外できない。

判別A/B:

- 対象は`case/45.isobutane_m6_d155/run_0540_lvcfh_v3_a1/`〜`run_0543_lvcfh_v3_b2/`の4本。Aは切替なし、Bは`FORGE_DIAG_FACE_H_DOUBLE=1`。他の実効設定・入力・バイナリを固定する。
- 全残差を毎step検査し、最初の増幅・非有限化のstep、保存場の有限性・正値性、序盤の局所状態を記録する。`FINITE`側には区間0〜1999の`check_convergence --segment`のVERDICTを添える。
- **Aが2本とも非有限、Bが2本とも有限** → この条件・期間の非有限化回避を支持。H2の否定、H3の除外、量子化による増幅機構の確定には使わない。
- **両側が2本とも非有限** → 切替単独による回避を棄却。破綻stepの延長は観測として記録し、同じ機構かは別途確認する。
- Aに有限が混じれば帰属不能、Bが分かれれば判別不能、入力・証拠に不備があれば`INVALID`。欠けた反復を旧runで補わない。

支持された場合の次の確認は、**保存した増幅前の共通状態・共通方向で、double経路の実残差応答と近似作用素を照合すること**。差分幅の再現性と再評価ノイズを先に検査し、過去に無効だった有限振幅評価も各評価経路内で閉じる。

棄却された場合の次の確認は、**未完了のU-J列ごとの照合**。[既存の合格は行列全体の最大値による正規化](/home/sano/work/forge-faceh/methods/time_integration/implementation.md:543)に限られる。まず同じ薄層モデルに対するD/K・壁拘束の整合を確認し、それが通ってから実残差との差を調べる。どちらも長期runより先に行う。

やらない方がよいこと:

- Bの2000 step完走を「安定」「収束」「H1の機構が証明された」とする。
- 両側の発散から、熱伝導Kの式・符号が誤りだと決める。
- マスク・CFL・緩和・LHS精度まで同時に変える。
- 判定後すぐに破綻の証拠を削除する、結果を見て期間や閾値を変更する。
- この結果だけでdouble評価を既定化する、診断機能を削除する。

呼び出し側の前提への異議: **H1とH2は排他的ではない。** また、マスク5の比較は親planで探索的観測に格下げされており、「熱伝導Kが必要」は確定していない。point仕上げの保存判定`NOT_SUPPORT`は確認できたが、今回のライン破綻について面エンタルピー精度を除外する証拠にはならない。「LHSは変えない」は固定状態での式についてであり、反復で状態が変わった後のLHS数値まで同一という意味ではない。

不足情報: 対象の旧run・初期HDF5・方向微分の元配列・AWSバイナリはローカルにない。旧run数値はplanの記録として扱い、独立再計算はしていない。コード、保存されたpoint判定集計、判定器の模擬入力を確認した。`forge`は起動せず、ファイルも変更していない。**plan未反映。反映先は対象planの§4・§6・§5.1。**
