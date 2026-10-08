# codex 諮問 (diagnose): te-wake-ab-result

- **brief**: [`notes/reviews/briefs/2026-10-08-te-wake-ab-result.md`](../../notes/reviews/briefs/2026-10-08-te-wake-ab-result.md)
- **plan**: [`plans/active/convection-zero-thickness-edge-reconstruction.md`](../../plans/active/convection-zero-thickness-edge-reconstruction.md)
- **date**: 2026-10-08
- **commit**: `bdefe043` (feature/sern-design)
- **codex**: effort `xhigh`, 10.2 min, rc=0
- **結論**: **§6.0 の結果 A を採用し、次は `te_wake_blend_H=1.0`・`w` 無効を固定した m10_on の g3→g4 格子感度確認を行う。**
- **extra**: `case/46.sern_design/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 2）

| 重大度 | 採否 | 根拠と対案 |
|---|---|---|
| **Major** | **§6.0 の「結果 A」は採用。生産採用完了の認定は却下** | 床カウンタの代用は[事前登録で許容済み](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:195)。B は保存40標本すべてで床近傍0節点、全域で T < 200 K が0節点だった。ただし、これを「全更新の床補正0件」と読み替えてはいけない。恒久採用では[§6.2 #1](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:223)の毎更新カウンタが必要。 |
| **Major** | **3D メッシャの既定を直ちに1.0へ変更する案は却下** | 現在、格子署名の区間識別への追加は [`w` 有効時だけ](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:458)。[学習データの採否](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:91)にも今回の格子変更を区別する条件がない。既定0を維持し、生産候補の YAML に1.0を明示する。採用前に、曲線版・実格子署名・変換器の来歴を評価系列へ結び付ける。座標・接続の署名だけでは双対幾何の変更を検出できない点にも注意する。 |
| **Major** | **格子修正の受入れは別 plan に分ける。`w` は予備として保持し、§6.1 は保留** | 現在の[§6.2](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:211)は、2リング・マスク幅・発火試験を含む `w` 用の受入れ表である。遡及的に格子修正へ読み替えない。格子側には床・温度・力の窓・格子感度・旧評価の判定を引き継ぎ、マスク固有項目を明示的に除く。`w` は[予備として残す既決方針](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:117)を維持する。小さい2Dでの加熱懸念が未解決のため、今回の成功後に追加投入する理由はない。 |
| **Major** | **力の差を「冷点除去による物理効果」と認定すること、格子感度許容による旧評価の持ち越しは却下** | 物理形状・BCを固定して格子を変えた比較なので、観測したのは離散化変更への応答である。採用判断に原因別の一意な分解は不要。m6_on・m4_off の旧評価は[§6.2 #6・#7](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:233)で別途判定する。旧 m10_on は床到達により不成立のまま。旧データは保存し、修正後の学習系列への混入を止める。 |
| **Minor** | **「A′は全列プラトー」を訂正** | [A′の保存VERDICT](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-te-wake-ab/run_1078_tewake_A0_m10/CONVERGENCE_VERDICT.txt:9)では `rms_roK` が **0.5桁低下・falling**。両 run の総合判定が `NOT CONVERGED` であることは変わらない。対案は「A′は多くの列がプラトー、`rms_roK` は低下中」と書くこと。 |
| **Minor** | **力の比較を末尾窓平均へ統一し、`a` の定義を訂正** | 末尾20標本の平均差 B−A′は、C_T **+6.544×10⁻⁷**、C_T_with_shear **−7.900×10⁻⁵**、C_L **+2.084×10⁻⁴**、C_M **−4.526×10⁻³**。ブリーフは最終標本の差である。また [`v3_farfield_eval.py:25`](/home/sano/work/forge-sern-design/case/46.sern_design/v3_farfield_eval.py:25) の `a` は半レンジではなく **max｜q−平均｜**。今回はこちらの厳しい定義でも窓条件を満たすため、合否は変わらない。 |

結論: **§6.0 の結果 A を採用し、次は `te_wake_blend_H=1.0`・`w` 無効を固定した m10_on の g3→g4 格子感度確認を行う。**

第 1 仮説: 局所格子変形と現行の再構成の組合せで冷点を解消でき、その効果は g4 でも維持される。確度: **中**。g3での十分性は支持されたが、g4への適用は未確認である。

根拠:

- A′：`case/46.sern_design/run_1078_tewake_A0_m10/` は、step 500～20000 の全40標本で節点517160が50 K。
- B：`case/46.sern_design/run_1079_tewake_B10_m10/` は、同じ40標本で床近傍0節点。判定区間 **10000 < step ≤ 20000** の最低温度は、R_TE **227.364 K**、R_SE **227.326 K**、RET **224.811 K**、全域 **210.606 K**。
- 保存CSVから正式ツールを再実行し、両 run の温度4量・力4量について **`OVERALL: ALL STEADY`** を確認した。Bの窓条件で最も許容に近い C_M も、`a=1.731×10⁻⁴`、`d=4.185×10⁻⁴ < 5×10⁻⁴`。
- 保存された残差判定は両方 **`NOT CONVERGED (stalled/plateau)`**。メッシュ品質は両方 **`PASS (AR≤5000, skew≤0.90)`**、双対閉性も **`PASS`**、Bの投入条件は **`ADMISSION VERDICT: PASS`**。

数値の根拠は [A′の記録](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-te-wake-ab/run_1078_tewake_A0_m10/)・[Bの記録](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-te-wake-ab/run_1079_tewake_B10_m10/)。run索引は [case README](/home/sano/work/forge-sern-design/case/46.sern_design/README.md)。

反証条件: 品質・初期場・判定窓が成立したg4で、冷点が再発・移転する、床補正が残る、または力の格子感度が登録許容を超える。その場合、**この修正をg3の生産評価へ適用する十分性**を棄却する。g3で観測した冷点消失自体は撤回しない。

第 2 仮説: g3での改善には、せん断層に対する節点配置・数値拡散の変化が強く寄与し、細分化すると改善が失われる。確度: **中、未確認**。

第 3 仮説: 下流77,181節点の補間が別の状態への移行を促し、初期場履歴が結果に残っている。確度: **低、未除外**。ただし冷点そのものを補間で消した説明は、下記のstep 0記録と整合しない。

判別 A/B: **A＝既存Bの修正済みg3、B＝同じ形状・作動点・1.0H曲線のg4**。変更対象は格子仕様だけとし、BC・数値設定・物性・介入の有無を固定する。

追加g4は `run_1079` の最終場からcross-mesh restartし、**20000 step、未定常なら＋20000を一度**。採否は末尾10000 step、500 step間隔の20標本で行う。既存g4は複数の格子パラメータを変える仕様なので、差分を列挙し、**複合的な格子感度**と呼ぶ。単一方向の収束次数は推定しない。

見る量は床補正、固定領域と復帰区間の最低温度・低温体積、4力量。力の平均差の許容は推力2量・C_Lが **0.002**、C_Mが **0.05**。両側で温度・窓・STEADY条件を満たすことを前提にする。生産受入れ用には毎更新床カウンタを備え、計測追加でバイナリを変えた場合はg3側も共通バイナリへ揃える。

→ **冷点なし・感度許容内なら**「改善はg3だけの現象」という仮説を、この比較範囲で退ける。  
→ **g4で冷点再発または感度超過なら**生産候補としての十分性を退ける。  
→ 品質不合格・転送不備・延長後も未定常なら判定不能。

やらない方がよいこと: 既定値の一括変更、結果を見て遷移長を調整すること、`w`との同時変更、3Dの成功による2Dの自動受入れ。g4確認後に、2Dへの同方式の実装と独立した検証、m6_on・m4_offおよび設計差の再評価を進める。旧評価の持ち越しは、その結果から判定する。

呼び出し側の前提への異議:

- **短い過渡だけで判断した試験ではない。** [Bのstep 0記録](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-te-wake-ab/run_1079_tewake_B10_m10/TE_MONITOR_step0.csv:2)にも、同じ座標・同じ節点517160の50 K冷点が残る。その後、床消失を20000 stepまで確認している。「500 step以内に消えた」ことだけを理由に保留する必要はない。
- **R_M未定義は、今回の低温移転判定を覆さない。** 全域でT < 200 Kが0節点なので、未抽出の領域へ150 K未満の冷点が移った可能性も保存時点では排除できる。ただし、R_M固有の最低温度をSTEADYと判定した事実はない。格子修正用planでは監視領域を明示し直す。
- **同一GPUでの同時実行だけでは比較を無効にしない。** 性能比較やビット再現性の根拠には使えないが、今回の温度・力の登録条件を取り消す理由は見つからない。
- **力の差の原因別分解は、採用の前提にしなくてよい。** 説明を進めるなら、同じ壁面上で圧力・摩擦・タグ別モーメント差を積分し、共通物理断面のT・P・U・組成・h0とエネルギー収支を照合する。これは変化の所在を測る方法であり、「冷点除去の物理寄与」と「格子寄与」を一意に分離する方法ではない。

不足情報: ローカルの写しには入力・最終HDF5、残差CSV、実config、最初100 callの再構成・エネルギー帳簿がない。温度・力はCSVから再検証したが、残差と最終場の健全性は保存VERDICT・`metrics.json`を根拠にした。また、`RUN_PROVENANCE.txt`だけではバイナリと`aa712de0`の対応、単面介入OFFを独立に確認できない。これらは機序の確定・生産受入れ前に一次記録を補う対象である。

**plan未反映。ファイル変更・`forge`起動は行っていない。** 呼び出し側で既存planの§4.2・§6.0・§5.1 #2eに結果を記録し、格子修正の受入れを別planへ分けること。
