# codex 諮問 (diagnose): qacc-axisym-state-precision

- **brief**: [`notes/reviews/briefs/2026-10-11-qacc-axisym-state-precision.md`](../../notes/reviews/briefs/2026-10-11-qacc-axisym-state-precision.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `5a6d3e8f` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **現案のままの実装・投入は見送り、作業用 `Q` と更新基準 `Q_N/Qacc` の射影を分けて修正し、単体試験を通した後に `k=10` 固定の E・E′・F 比較を行う。**
- **extra**: `plans/active/time_integration-fp64-accumulator-rollout.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **現案のままの実装・投入は見送り、作業用 `Q` と更新基準 `Q_N/Qacc` の射影を分けて修正し、単体試験を通した後に `k=10` 固定の E・E′・F 比較を行う。**

指定ファイル全文を読み、コードと照合した。ファイル変更・`forge` 起動は行っていない。参照 run はローカルに存在せず、過去の実測値・VERDICT は計画文書の記録として扱った。

採否表:

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| **Major 1** | 現在の `Q` から KE を求め、`roe=float(Qacc_roe)` とする軸ピン | **却下。境界の重複で成立しない。** 壁 no-slip が先に `Q` だけを書き換えるため、軸ピン直前の `Q=float(Qacc)` は一般には偽。[main.cpp:1892](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1892)、[nodeWallDirichlet_d.cu:24](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:24)。`Qacc` には **射影前の `Q_N` から求めた KE** を反映し、作業用 `Q` の既存の射影は別に保つ。 |
| **Major 2** | S を「sub-ULP 増分の吸収が偏りに必要」、R を「残差評価側」と解釈 | **却下。感度試験として採用。** `Qacc` は増分の完全吸収だけでなく、切上げ・切捨てを含む累積丸め全体を変える。[qAccumulator_d.cuh:25](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/qAccumulator_d.cuh:25)。S は「この変更で登録窓の偏りを抑制」、R は「この変更だけでは抑制不十分」とする。 |
| **Major 3** | 旧 A・B・B′ の流用 | **条件付き採用。** F の６条件のいずれかが外れたら**判定不能**とし、結果を見て F 単独へ対照を差し替えない。SASS に加え、ホスト側の初期化・引数・実効設定の差分も確認する。初期化後のエネルギーは精度間で一致しない記録があるため、A は「保存済みの参照軌跡」である。[float plan:799](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:799)。 |
| **Major 4** | 今回の試験で軸対称の拒否を全面解除 | **保留。** rollout の回帰・改善ゲートを今回の単体試験と F だけで置き換えない。[rollout plan:170](/home/sano/work/forge-integ-1005/plans/active/time_integration-fp64-accumulator-rollout.md:170)、[同:200](/home/sano/work/forge-integ-1005/plans/active/time_integration-fp64-accumulator-rollout.md:200)。診断用変更と対応完了を分け、境界重複も検証対象にする。 |
| **Minor 1** | E″ を足して `k` を下げる | **却下。今回は `k=10` を維持。** ３本のばらつきから小さい倍率の妥当性は決まらない。§6.33 の僅差の未達は、その判定を変更する根拠にならない。[float plan:820](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:820)。倍率は信頼係数でなく暫定の診断基準と明記する。 |
| **Minor 2** | reconcile の２時点の mask と累計上限 | **採用するが、保証範囲を限定。** 累計は保存量別・位置別ではないため、途中の想定外採用を否定できない。[update_d.cu:51](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:51)。両 ON 腕で確認し、２時点の観測を「全期間で軸・内部の採用ゼロ」と表現しない。 |

**Major 1 の反例は CPU の float32 演算で再現した。**

`axis_flag=wall_flag=1`、ρ=1、ρuₓ=1、ρuᵣ=0.125、ρE=2、残余ゼロ、dq=0 とする。

| 時点・量 | 現行 OFF | 提案 ON |
|---|---:|---:|
| 壁 no-slip 後の `roe` | 1.4921875 | 1.4921875 |
| 軸ピン後、残差評価へ渡す `roe` | 1.4921875 | **2.0** |
| commit 後の `roe` | 1.9921875 | **2.0** |

壁ピン後の `roUy` はゼロなので、提案の KE はゼロになる。一方、現行は別に保持した `roUyN=0.125` から更新基準の KE を引く。[axisymmetricSource_d.cu:315](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:315)、[同:319](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:319)。

さらに `roe=float(Qacc_roe)` が、直前の壁による KE 除去まで取り消す。**step 末尾の reconcile では、この残差評価前の誤りを検出できない。** 今回の `case/45` に重複ノードがあるかは未確認だが、一般の軸対称対応としては反例が成立する。

修正は次の契約を推奨する。

- 作業用 `Q` と `Q_N` への既存の射影を維持する。
- `Q_N` を射影する**前**に、既存と同じ `flow_float` の式・密度ガードで `ke_N` を計算する。
- ON では `Qacc_roe -= double(ke_N)`、`Qacc_roUy=0` を追加する。
- ここで作業用 `roe` を `Qacc` から上書きしない。ミラーへの書戻しは既存の commit が担う。

これなら変更対象を更新基準の累積に限定できる。単体試験には提示済みの (a)〜(c) に加え、上の**壁→軸→commit の順序**、非軸ノードの不変性、非ゼロ dq、残余を持つ複数 step を含める。「ビット一致」は既存試験の数値比較 `on != off` ではなく、生のビットで検査する。[test_qacc_commit.cu:56](/home/sano/work/forge-integ-1005/solver_density_cuda/tests/unit/test_qacc_commit.cu:56)。

第 1 仮説: **平均流の累積丸めを除くだけでは、登録窓の偏りを十分には抑えられない。** 確度: **低**  
　根拠: commit 前に EOS の保存量再構成があり、要求増分も `flow_float` のまま入力されるため、ON 後も精度依存の経路が残る。[dependentVariables_d.cu:209](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/dependentVariables_d.cu:209)、[update_d.cu:291](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:291)。ただし、それらが θ_r の偏りを支配する証拠は未確認。  
　反証条件: 修正済み実装と全前提を通した E・E′ が、両窓・３断面すべてで結果 S を満たすこと。

第 2 仮説: **平均流の累積丸めの変更だけで偏りを抑制できる。** 確度: **低**。既存の累積機構はこの介入を可能にするが、case/45 の帰属は未測定。  
第 3 仮説: **float のまま残る乱流・化学種の更新と平均流の相互作用が効く。** 確度: **低**。SST 凍結試験は判別不能の記録であり、除外しない。

判別 A/B: **変更は `time.deltaT.qAccumulatorFP64: 0 → 1` の１点。** 修正と単体ゲート通過後、共通入力から E・E′・F を各30,000 step、500 step ごとに出力する。W1・W2、θ_r の３断面、D、S/R の式は事前固定し、`k=10` とする。

- **S なら**第１仮説を棄却し、「平均流の累積方法の変更で、登録窓の参照からの偏りを抑制できた」。
- **R なら**第２仮説の「この変更だけで十分」を棄却する。残差評価、増分生成、スカラー更新のどれかは決まらない。
- 中間なら判別不能。延長・窓の選び直し・倍率の変更はしない。

旧対照の流用は、F の登録条件と初期化後の９保存量の照合が通る場合だけ認める。S の効果対ばらつき条件では、`N_B` に代えて **B・B′・F の Δ の最大値−最小値**を使い、F で観測した差も含める。A の再実行差は未測定なので、A への近さを FP64 一般への同等性とは呼ばない。

やらない方がよいこと: 軸の commit 分岐の `nullptr` を有効化すること、軸ピン直前に全 `Qacc` を `Q` へ reconcile して残余を捨てること、R を根拠に左辺・増分生成を除外すること。現在の commit 分岐は明示的に無効であり、rollout の「軸上 `roUy` 厳密ゼロ」は**測定時点を軸ピン直後と区別して定義し直す必要がある**。[update_d.cu:345](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:345)。

呼び出し側の前提への異議: 「状態保存か残差評価か」の二択は受け入れない。今回測れるのは平均流の累積方法への感度である。また、§6.31 の正式な判定不能理由は初期化後の Q₀ 不一致であり、「差が再実行の揺れ程度」は登録外の記録である。過去の除外判断を強める根拠にはしない。

不足情報: AWS-B の `case/45.isobutane_m6_d155/run_0487_tr_fp64/`、`run_0488_tr_f32/`、`run_0489_tr_f32b/`、`run_0590_ds_f32/`、`run_0591_ds_f32b/` の実データ、正式 VERDICT、実効設定・バイナリ照合、および今回の格子の軸・壁マスクの重複数。所在の索引は [case README:212](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:212)。文書上の `check_convergence` は **CHECK FAILURES** であり、本レビューは収束を認定しない。

**plan 未反映**。依頼どおり読み取り専用の諮問結果として出力した。
