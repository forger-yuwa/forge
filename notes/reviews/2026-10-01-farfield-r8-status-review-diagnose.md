# codex 諮問 (diagnose): farfield-r8-status-review

- **brief**: [`notes/reviews/briefs/2026-10-01-farfield-r8-status-review.md`](../../notes/reviews/briefs/2026-10-01-farfield-r8-status-review.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-10-01
- **commit**: `2c7cd880` (feature/sern-design)
- **codex**: effort `high`, 4.9 min, rc=0
- **結論**: **同じ新バイナリ・同じ g3 初期場で旧輸送の対照を追加し、既存の新輸送 run と比較して、3D の変化を輸送モデルに帰属できるか確かめる。**
- **extra**: `plans/active/tooling-nozzle-sern-chain.md`, `case/46.sern_design/README.md`, `case/58.farfield_verification/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

数値は対象 run の一次データが手元にないため、台帳・plan の報告値として扱った。コードと判定経路は確認した。`forge` は起動していない。

| 判断 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| V3 の必要幅は 2.50 H | **限定付き採用** | 「当該形状・m6_on・g3・旧輸送・固定リミッタ基準値で、試験した幅のうち最小の適合幅」とする。新輸送・別格子・別作動点への適用は未検証。[farfield plan:205](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:205) |
| 格子差 G の改善は farfield 化による | **限定付き採用** | 同じ新バイナリで slip 対照を取った点は妥当。台帳値では C_M の G が slip −0.04071、farfield +0.01590。**今回の格子間差の縮小**までは帰属できるが、真値への接近・収束次数改善・反射除去の証明ではない。[case README:348](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:348) |
| R8 段 (i) の合格 | **限定付き採用** | 演算経路による保存係数の再現を追加した処置は、前回指摘に沿う。(0) の FAIL を保持し、丸め差の例外受理と (1) の準定常回帰 PASS を分ける。報告上、残差は6本とも `NOT CONVERGED` 相当のプラトーであり、収束解の一致とは呼ばない。[R8 行:432](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:432) |
| 3D の輸送感度は 2D より一桁大きい | **解釈を却下・Major** | 変化率の分母が違う。報告値から、ΔC_L は 2D −8.0e−5／3D +1.410e−4、ΔC_M は +1.65e−3／−3.4931e−3。絶対変化は約 **1.8倍／2.1倍で、符号も逆**。「約10倍」の大部分は基準係数が約4.4〜4.5倍違うため。絶対差と面別寄与で比較する。[R8 行:432](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:432)、[case README:354](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:354) |
| 生産への二つの切替で精度確認も完了 | **要再検証・Major** | ユーザ決定による限定運用は維持してよい。ただし C_M の D=0.0041 が「εの82%」という比較は、**輸送モデル差を領域感度の許容で測っている**。これだけで合否や残り18%の余裕を定義できない。新輸送での格子差・幅感度を確認し、同じモデルの G+D で評価する。[R4d:1942](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1942) |
| V2c はソルバの角の非一意性と確定 | **却下・Major** | 絶対残差が凸角近傍に 4.6e−3／7.3e−3 残り、plan 自身が「二つの離散定常解」を棄却している。観測は**初期場に依存した停滞状態**であり、定常解の非一意性や farfield 無関係までは確定しない。[farfield plan:148](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:148) |
| V2c・V2d-2 を切り出せば accepted にできる | **却下・Major** | 汎用境界として検証するというユーザ決定が残る。診断作業を別 plan に移すことは可能だが、元 plan の未達条件・依存関係は残す。生産切替の承認は検証範囲縮小の承認ではない。[farfield plan:22](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:22) |
| `--force-species` は一回だけで済む | **却下・Major** | Python は `species_input_unverified=1` を記録照合前に未検証と分類し、既定で停止する。許可すると宛先属性を消すため、ソルバ起動にも許可が必要。**標準の restart 経由では繰り返し必要になる**。対してソルバ直接読込はハッシュ一致なら印を継承して通す。既知の未決事項 #3d である。[forge_species.py:810](/home/sano/work/forge-sern-design/solver_density_cuda/tools/forge_species.py:810)、[speciesDB.cpp:1340](/home/sano/work/forge-sern-design/solver_density_cuda/input/speciesDB.cpp:1340)、[種DB plan:227](/home/sano/work/forge-sern-design/plans/active/thermophysics-solver-owned-species-db.md:227) |
| EXH 構成種がソルバ内蔵という説明 | **訂正を採用・Minor** | 読込側は `legacy_builtin: solver` の種だけを採用する。現行の外部生係数ファイルは必要。依頼元へ訂正を伝え、不要として削除しない。[speciesDB.cpp:118](/home/sano/work/forge-sern-design/solver_density_cuda/input/speciesDB.cpp:118) |

結論: **同じ新バイナリ・同じ g3 初期場で旧輸送の対照を追加し、既存の新輸送 run と比較して、3D の変化を輸送モデルに帰属できるか確かめる。**

第 1 仮説: 輸送モデル変更による壁圧分布の変化が、今回の ΔC_L・ΔC_M の主成分である。 **確度: 中**

根拠: 台帳上の変化は上表のとおり。ただし変更対象は μ だけでなく熱伝導率 λ も含む。[runner_sern.py:313](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:313) また、報告する C_L・C_M はノズル面の**圧力積分**で、壁摩擦の直接寄与ではない。輸送による変化なら壁圧への間接効果として説明する必要がある。[runner_sern3d.py:303](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:303)

反証条件: 同じ新バイナリでは輸送を変えても差が再現せず、旧輸送の対照自体が旧 `run_0997` から今回の変化に近い量だけ動くこと。

第 2 仮説: マージ後のコード・生成設定の変更が混在している。 **確度: 中、未確認。** `run_0997` と `run_1011` はバイナリ世代が違い、準備スクリプトも `physProp`・`turbulence` 全体と BC を生成し直す。「輸送だけの差」は生成後の照合が必要。[prod_restart_setup.py:28](/home/sano/work/forge-sern-design/case/46.sern_design/prod_restart_setup.py:28)

第 3 仮説: 輸送変更の応答に格子・側方幅との相互作用がある。 **確度: 低、未確認。** 新輸送は g3・2.50 H の一点なので、旧輸送の G・幅感度を移せる証拠がない。

判別 A/B: **変える要因は輸送モデルだけ。**

- A＝現行バイナリ・AMB/lump・farfield のまま旧 Sutherland 輸送。
- B＝同じ条件で種ごとの輸送。既存 `case/46.sern_design/run_1011_prod3d_ff_tr` → `run_1012_prod3d_ff_tr_cont20k` は、実バイナリ・初期保存量・他設定の一致を確認できれば再利用する。
- A は新しい `run_*` に同一格子 restart で作り、まず20000 step、500 step間隔で保存。既存規則の前後10000 step窓を満たさなければ20000 step延長し、なお未達なら判定不能とする。
- 見る量は4係数の平均・振幅・前窓差、およびノズル面別の圧力力・モーメント。`check_convergence` と `check_quasisteady` の原本、NaN・置換回数も残す。
- 診断許容は事前に τ=0.2ε と固定する。**A が旧 `run_0997` を全量 D≤τ で再現し、B−A が既報の符号・変化量を τ 内で再現するなら**輸送主因を支持する。**B−A が τ 内に収まり、A 自体が既報変化を再現するなら**輸送主因を棄却する。中間結果は複合要因として保留する。

やらない方がよいこと: 「0.22%だから異常」「ε以内だから新設定も認定済み」と判断すること、V2c の出口端を除いて合格にすること、時間刻み細分で誤差が増えた事実だけから float32 を真因と断定すること。未検証履歴を消すために `species_input_unverified=0` を書くことも避ける。

呼び出し側の前提への異議:

- **`GATES PASS` は残差収束ではない。** `require_residual_pass` が無効なら残差 PASS を必須にしない。[sern_gates.py:292](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:292) 今回の限定的な生産運用と、収束解・汎用境界の認定は分ける。
- 未検証の印を保持すること自体は正しい。問題は「一回だけ移行」という説明と標準ツールの挙動の不一致である。メモリ上の属性を使った読み取り判定でも、既定拒否・許可時に継承属性なしとなることを再現した。対案は #3d で、**記録完全性と宛先互換性を確認できる印付き場は、印を1のまま継承する**仕様を検討すること。現行仕様を変更する判断として明記する。
- 優先順位は、上記の輸送対照 → 新輸送での格子・幅感度 → V2 未達と独立参照 → result レビュー → main 統合。#3d の運用訂正と内蔵種の説明訂正は、次の restart 前に行う。
- accepted への最小経路は、V2c の必要条件と全線判定、V2d-2 の時間精度を満たし、独立参照の精度確認・比較記録、残る実装項目と docs を完了して result レビューを受けること。独立参照は境界合否とは別だが、計画上の未完了作業である。別 plan への移管だけでは完了にならない。

不足情報: 対象の `case/58.farfield_verification/run_*` と `case/46.sern_design/run_0991`〜`run_1012` は手元にない。実 config 差分、バイナリ識別、初期保存量の一致記録、残差・準定常性 VERDICT、力時系列、R8 演算経路照合の出力原本が必要。したがって、数値の PASS は独立再認定していない。

**plan 未反映。** 依頼どおりファイルは変更していない。反映先は farfield plan §5.1・§6.1、SERN chain §5.1 R8、種DB plan §5.1 #3c残（#3d）。
