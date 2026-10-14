# codex 諮問 (diagnose): port-reg-redo

- **brief**: [`notes/reviews/briefs/2026-10-13-port-reg-redo.md`](../../notes/reviews/briefs/2026-10-13-port-reg-redo.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-13
- **commit**: `96515d46` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.2 min, rc=0
- **結論**: **判定器の入力検査を修正し、量別の固定許容を事前登録した Q2 の群比較で R1・R3 を再検証することを推奨する。main 統合は保留する。**
- **extra**: `case/45.isobutane_m6_d155/port_reg_an.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/port_reg_judge.json`, `notes/investigations/2026-10-13-main-integration/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| R1・R3 の旧判定 | **採用：不合格を維持** | 登録式は [plan §6.59](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1513)。ただし、この不合格だけで実装の回帰が確定したとは言えない。新しい試験は別の判定として記録する。 |
| Q1：D ≤ E_A + E_B | **却下・Major** | [依頼案](/home/sano/work/forge-integ-1005/notes/reviews/2026-10-13-port-reg-redo-diagnose.prompt.md:70) は、変更後のばらつきが増えるほど通りやすい。例えば旧群 `[0,0,0]`、新群 `[-100,10,120]` でも D = 10 ≤ 220 で合格する。**ばらつきは許容枠を広げる側ではなく、固定許容を消費する側に置く**。 |
| Q2：事前に固定した許容 | **方式を採用・数値は要再検証（Major）** | [提示された相対許容](/home/sano/work/forge-integ-1005/notes/reviews/2026-10-13-port-reg-redo-diagnose.prompt.md:71) の 1e−6／1e−4 には、量ごとの根拠がない。ゼロ近傍、乱流量、局所残差、`dt_local` も未定義。固定許容 τ_q を量ごとに定義し、**D_q + E_A,q + E_B,q ≤ τ_q** を用いる。これは観測した標本に対する保守的な判定であり、信頼区間ではない。 |
| Q3：SASS 同一性で回帰を代替 | **代替は却下・補助証拠として採用（Major）** | [plan #26](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:188) の SASS 比較は有力だが、同じ命令でも入力・呼び出し条件・配列の初期化が違えば結果は変わる。`test_check_solver_config.py` の PASS だけで構造体配置の整合を保証したとは扱えない。[CUDA 単体試験](/home/sano/work/forge-integ-1005/solver_density_cuda/tests/unit/test_precond_node_rows.cu:4) も対象カーネルの直接起動である。実効設定・入力・呼び出し側の確認を補い、run 比較は残す。 |
| 現行判定器の再利用 | **却下・Major** | [判定処理](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/port_reg_an.py:29) を抽出して確認したところ、後群の NaN、空の比較対象、後群だけの追加量、前bの broadcast 可能な形状不一致がいずれも `pass=True` になった。まず非有限値・欠損・空配列・キー集合／形状／step 列の不一致を拒否する。実 run にこれらがあったと確認したわけではない。 |
| main 統合を先行 | **却下・Major** | [統合条件6と停止条件](/home/sano/work/forge-integ-1005/notes/investigations/2026-10-13-main-integration/README.md:26) が未達。再設計した回帰の結果を待つ。また、同文書の確認4は入力欠損で未検証であり、R1・R3 の解決だけで統合条件全体を満たしたことにはならない。 |

結論: **判定器の入力検査を修正し、量別の固定許容を事前登録した Q2 の群比較で R1・R3 を再検証することを推奨する。main 統合は保留する。**

第 1 仮説: **旧基準が再実行のばらつきを許容限界として扱ったことが、不合格の主要因である。** 確度: **中**

  根拠: [判定 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/port_reg_judge.json:1) は全6組で `bitident_pre=false`、`pass=false`。判定器は各量について、旧版のたった1組の差を上限にしている。同じ分布からの交換可能な標本でも、後−前aの差が前b−前aを超えることはある。これは同等性の許容限界ではない。ただし、**差の原因が atomicAdd だけという説明は未確認**。

  反証条件: 入力・実効設定・比較対象の同一性を確認した再試験で、ある量の群間差が **D_q > τ_q + E_A,q + E_B,q** となること。この場合、「判定基準だけの問題で、許容外の変化はない」という説明を退ける。

第 2・第 3 仮説:
- 第2：ホスト側の設定受け渡し・初期化・呼び出し条件に、SASS 比較で検出できない変更がある。確度: **低、未確認**。
- 第3：実装変更が加算順などを介して再実行のばらつきを増やした。確度: **低、未確認**。平均差だけでは除外できない。

判別 A/B: **旧版／移植版のバイナリだけを変える、各3本・10 step の比較を、既存の6組で行う。**

- R1：case/45 の point／B0 × float／FP64。R3：case/48 × float／FP64。合計36本、各回は同じ初期入力から始め、継続計算にはしない。
- GPU は占有し、旧・新を交互に実行する。既存の台本は [「逐次だが試験GとGPU共有」](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/port_reg.sh:9) なので、外部負荷の交絡を外す。
- 入力・BC・実効設定・実行環境を照合し、`VALUE` の全比較対象と全 `rms_*` 列を同じ節点・step で比較する。量別の固定許容には絶対許容も含め、結果を見る前に確定する。
- **結果A：全量で D + E_A + E_B ≤ τ** → この入力・10 step・登録許容の範囲で回帰合格。第1仮説を支持する。
- **結果B：いずれかで D > τ + E_A + E_B** → 第1仮説の「基準だけで説明できる」を退け、実装側の差を調べる。
- 中間は**判別不能**。ばらつきが大きいだけの結果を合格や実装回帰に二分しない。結果を見ながら本数や許容を増やして通さない。

やらない方がよいこと: Q1 で合格へ読み替えること、SASS 一致だけで閉じること、許容を今回の最大差に合わせること、10 step の比較を「収束解の一致」と呼ぶこと。

呼び出し側の前提への異議:

- **Major：「すべて1〜2倍」は誤り。** R1 point FP64 の `res_roUy` は 2.135e−4 対 2.015e−5、**約10.6倍**。[JSON:40](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/port_reg_judge.json:40)。絶対的な重要性は別途評価が必要だが、要約から落としてはいけない。
- **Minor：量名の取り違え。** B0 float の 5.138e−3 対 2.641e−3 は `Uy`。`roUy` は 0.1110 対 0.1080。[JSON:55](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/port_reg_judge.json:55)。判定資料と plan の表記を訂正する。

不足情報: 元の `case/45.isobutane_m6_d155/run_0763〜0779` と `case/48.flat_plate_cooled_m4/run_0065〜0070` はこの作業ツリーに存在せず、生の HDF5・残差履歴・実効設定を再確認できない。量別許容の根拠と、実際に使用したバイナリに結び付く SASS 比較資料も必要。**本回答は設計提案であり、plan 未反映。依頼どおりファイル変更・forge 起動は行っていない。**
