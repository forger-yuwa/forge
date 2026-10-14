# codex 諮問 (diagnose): cold-pair-result

- **brief**: [`notes/reviews/briefs/2026-10-08-cold-pair-result.md`](../../notes/reviews/briefs/2026-10-08-cold-pair-result.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-08
- **commit**: `3ea01d72` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.8 min, rc=0
- **結論**: **既存の両腕5枚を使い、格子列の直接積分から一定 x 断面への補間にだけ切り替える、CFD 0 step の抽出 A/B を先に行ってください。**
- **extra**: `case/45.isobutane_m6_d155/cold_pair.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/V_c45.json`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/gates_aws.json`, `design/forge_design/feedback/deltastar_integral.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**現行の「判定不能」を維持します。次は閉包の変更より先に、抽出座標の A/B を行ってください。**

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **NS の抽出値をそのまま閉包比較に使う：要再検証** | 今回の壁法線格子は列内で x が変わります（[mesh2d.py:241](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:241)）。しかし抽出器は `N["x"][i,0]` を列全体の x とし、異なる x の `ρUx` を半径方向に積分しています（[deltastar.py:437](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:437)）。保存された壁・格子設定から座標を再生成すると、x≈40 の列内の x 幅は **0.1533 r_t＝11.75 mm**。**一定 x 断面への補間を先に検査**すべきです。この幅だけから δ_E の誤差量は決められません。 |
| **Major** | **ゲート2の見直し：問題提起は採用、提示された代案は却下** | 別格子・別精度の絶対残差床を普遍的な上限とする根拠は不足しています。しかし「横ばい／低下、RISING なし＋派生量 STEADY」だけでは、高残差の停滞を通します。`run_0181` は既に δ_E が STEADY でも、`rms_roe=10.16` に対し基準は `0.545`。**同じ格子・精度で、残差の正規化、保存収支、対象量の精度と結び付いた上限**を定める必要があります。その根拠は今回の資料にはありません。 |
| **Major** | **ゲート4の縮流部超過を無害とする：要再検証** | CSV の再集計でも冷却側の超過面積は **24.253%＞10%**。全壁の `PASS` は領域別ゲートの合格ではありません。最大 y₁⁺≈1.23 が小さいことや、試験部の超過が0%であることから、上流履歴の影響を除外できません。**現基準の不成立を残し、緩和には当該条件の格子感度の裏付けを要求**してください。 |
| **Major** | **「熱閉包の形を選ぶ」問いを置き換える：却下** | 相対的な優劣と絶対的な予測精度は別の問いで、[plan §6:238](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:238) は既に両方を判定しています。ゲート成立後なら「Bが相対的に近いが、両者とも精度不足」は有効な結論です。現在の `decision_by_kf` 内の支持は**条件付き計算結果**で、正式結論ではありません。 |
| **Major** | **`judge` によるゲート集約：修正が必要** | 最終判定は準定常・残差・壁解像だけを論理積にし、**ゲート1の全場 NaN/Inf 検査を読み込んでいません**（[cold_pair.py:395](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair.py:395)）。また収束判定の終了コードを保存しても採否に使わず、「判定不能」も明示的に拒否していません（[同:278](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair.py:278)、[同:336](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair.py:336)）。**検査結果の欠落・判定不能を不合格にし、許容する plateau と区別**してください。今回の結果が誤って合格したという指摘ではありません。 |
| **Minor** | **不確かさ U の端点処理：修正を採用** | [cold_pair.py:368](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair.py:368) の `convolve(..., mode="same")` は外側をゼロ埋めします。有効点数で正規化すると、x=40 の U は **0.845%→1.236%**。生産 k_f の誤差区間は A **[0.42964, 0.46542]**、B **[0.38956, 0.42435]**となり、区間分離は残りますが余裕は縮みます。帯外の実データを使うか、有効点数で正規化してください。 |
| **Major** | **壁温変更時の NS 反復必須化・v2 の既定化：修正付き採用／現時点の既定変更は却下** | **壁温変更後の NS 再評価は必要、壁更新の反復は設計公差を外れた場合に必要**としてください。「温度を変えたら必ず形状を反復」は [plan §4.6:133](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:133) の感度判定より強すぎます。v2 の熱力学・輸送整合の価値は残りますが、熱閉包だけの比較は v2 全体の検証ではありません。**v2 は opt-in を維持**し、全体の検証・較正条件の記録・ユーザ決定を既定変更の条件にします（[plan:156](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:156)、[同:201](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:201)）。 |

基準の改訂自体は可能です。ただし **旧 V-c45 の登録条件と「判定不能」を保存し、改訂日・観測済みデータ・変更理由・新しい条件を別版として記録**してください。同じデータでの再解析は事後解析と明記します。改訂を延長前の事前登録として扱ってはいけません。

`case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/` の延長結果だけでは、旧基準を通せません。断熱側 `run_0181` の残差床不成立が残るためです。断熱の元窓を保持する運用は登録どおりですが、冷却側は δ_E だけでなく **Q_w・比の準定常性、残差、壁解像も再判定**する必要があります。

結論: **既存の両腕5枚を使い、格子列の直接積分から一定 x 断面への補間にだけ切り替える、CFD 0 step の抽出 A/B を先に行ってください。**

第1仮説: **抽出座標の不整合が、冷却比 R_NS の1%精度判定に無視できない寄与を持つ。** 確度: **中**。  
根拠: 上表の格子生成と抽出コードの不整合。構造的な不整合は確認済みですが、R_NS への寄与は未測定です。**e≈0.4の主因と断定していません。**  
反証条件: 下記 A/B で、全5枚・試験部全域の `|R_B/R_A−1|` が1%未満となること。「1%以上の寄与」という仮説を棄却します。

第2仮説: **抽出を整合させても残る差は、CONTUR の分布閉包・発達履歴にある。** 確度: **中、原因箇所は未確認**。θの約2倍化だけでは摩擦過大を証明できません。運動量式には H と加速項も入ります（[deltastar_integral.py:273](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:273)）。

第3仮説: **冷却 NS の未発達または解像誤差が比較値に残っている。** 確度: **中、寄与量は未確認**。冷却側の `DRIFTING` と縮流部ゲート不成立が根拠です。

判別 A/B:

- **対象**：`case/45.isobutane_m6_d155/run_0181_ns_coldmesh_ad/` と `case/45.isobutane_m6_d155/run_0182_ns_coldmesh_tw300/` の **80000〜100000 step、各5枚**。追加 CFD は **0 step**。
- **A**：現行どおり格子列を直接使用。**B**：実際の `nozzle.h5` 座標を使い、隣接列から `ρUx` と半径位置を同じ物理 x へ補間して使用。
- Euler 参照、帯 E の選択アルゴリズム、求積、符号処理、平滑化条件は固定。各腕の δ_E と冷却比を、x/r_t=40〜94 の0.25刻みで比較します。
- **最大差≥1%なら抽出不整合の寄与を支持／最大差<1%なら第1仮説を棄却**。後者でも摩擦閉包が原因と確定するわけではありません。

この比較は同じ保存場に対する後処理の感度試験なので、過渡場でも実施できます。定常的な冷却効果の認定には使いません。θ・H・C_f の抽出は、その後の閉包診断に必要です。その際、**円環の根から得る δ_E と、CONTUR の積分 δ* を区別し、δ_E/θをそのまま同じ H と呼ばない**でください。

やらない方がよいこと: **残差上限や縮流部上限を今回の値に合わせて緩めること、δ_Eだけの静定で延長を合格にすること、C_f・N・回復係数を同時に較正すること。**

呼び出し側の前提への異議: **「冷却で21〜27%薄くなる」「CONTURが冷却効果を外す」は、現段階では保存窓における抽出結果とモデルとの差です。定常解に対する確定した予測誤差ではありません。**

保存 NPZ を `check_quasisteady.classify` で再判定した結果は、全5枚、drift・osc許容0.1%、試験部217点で以下でした。

- 断熱 δ_E・壁温：**`STEADY`**
- 冷却 δ_E・冷却比：**`DRIFTING`**
- 冷却 Q_w：**`DRIFTING`**（窓内 drift 約7.1%）

保存済み収束判定は両腕とも、本段区間で **`NOT CONVERGED (stalled/plateau)`**。総合 **`判定不能 (ゲート不成立)`** を支持します。

不足情報: 元 run の HDF5・全残差系列・実効設定・NaN 検査記録・延長結果はこの作業ツリーにありません。したがって、元場の健全性と上述の抽出 A/B は未確認です。

ファイル変更・forge 起動はしていません。**plan 未反映**。採用時の反映先は `plans/active/tooling-nozzle-isothermal-wall-chain.md` の **§4.6・§4.7、§5.1 #18〜#20、§6 V-c45** です。
