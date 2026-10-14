# codex 諮問 (diagnose): cold-arm-mass-imbalance

- **brief**: [`notes/reviews/briefs/2026-10-08-cold-arm-mass-imbalance.md`](../../notes/reviews/briefs/2026-10-08-cold-arm-mass-imbalance.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-08
- **commit**: `c8f6d8c5` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.7 min, rc=0
- **結論**: **CFL を変える前に、保存場の台形積分流量と、同じ状態におけるソルバの数値質量流束・連続残差を照合する収支監査を行ってください。**
- **extra**: `case/45.isobutane_m6_d155/cold_xcheck.py`, `case/45.isobutane_m6_d155/README.md`, `notes/reviews/2026-10-08-cold-su2-crosscheck-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **「欠損＝冷却による質量蓄積」：要再検証** | 壁温ピンは `ρ` を変更せず、壁の残差射影も `res_ro` を消していません（[nodeWallDirichlet_d.cu:57](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:57)、[同:84](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:84)）。ただし、現在の ṁ は節点値の台形積分です（[cold_xcheck.py:94](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:94)）。**ソルバの数値流束と符号付き連続残差で照合するまで、物理的な蓄積とも非保存とも確定できません。** |
| **Major** | **総質量の step 差を入口−出口と直接比較／両コードの ΣCFL を揃える：却下** | `dt_local` は節点ごとに異なり、DPLUR は `V/dt_local` を含む連成系を緩和付きで解きます（[setDT_d.cu:215](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:215)、[timeIntegration_d.cu:827](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:827)、[同:1108](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1108)）。**共通の物理時間はなく、ΔM/Δstep は kg/s の収支ではありません。** ΣCFL もコード間の時間尺度にはなりません。各コードで収支・残差・比較量の静定を独立に満たし、効率はゲート到達までの壁時計時間で比べてください。 |
| **Major** | **直ちに CFL 4〜8／relax 1 で延長：却下** | 古い粗格子の成功記録だけでは不足です。同じ case の細分格子には、`cfl 5・relax 0.7` で本段 step 35 発散、同一起点の CFL 1 は短期安定という後続記録があります（[README.md:85](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:85)）。今回の上限を示す実測ではありませんが、旧上限の流用を正当化できません。**延長設定は (a) を維持し、次の作業は収支診断に限定**します。追加100000 stepの盲目的な投入も推奨しません。 |
| **Major** | **下流の R_NS・δ・θ・C_f を定常比較に使う：却下** | 再判定でも冷却側の x=40 の δ_E は `TRANSIENT-UNSETTLED`、Q_w は `DRIFTING`。保存済み収束判定には `rms_roOmega … RISING` があります。**現在の値は未収束場の参考値として保持**してください。上流の流量欠損率から下流の厚さ誤差を換算することはできず、一部の δ_E が `STEADY` でも影響上限は定まりません。 |
| **Major** | **断熱の熱収支 −1.03e−3 を求積誤差として通す：却下** | 報告値どおりなら絶対値は登録上限 1e−3 を超えています。評価式は対流エンタルピー流束と壁勾配熱流束であり、離散エネルギー収支そのものではありません。等温壁ではエネルギー残差の射影も入ります（[nodeWallDirichlet_d.cu:94](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:94)、[methods/boundary.md:680](/home/sano/work/forge-integ-1005/methods/boundary.md:680)）。**入口・出口の伝導／粘性仕事、壁の拘束寄与を含む離散収支と、現在の後処理指標を分けて監査**してください。forge の h0 は保存された `VALUE/h0` を使います。 |
| **Major** | **現在の比較ゲート実装：要修正** | [cold_xcheck.py:250](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:250) は `diverged`・`rising` の不在だけを確認し、収束判定の実行失敗・判定行欠落・判定不能を明示的に拒否しません。**判定失敗を不合格にし、plateau の受理条件を別に定義**してください。今回すでに誤合格したという指摘ではありません。 |

結論: **CFL を変える前に、保存場の台形積分流量と、同じ状態におけるソルバの数値質量流束・連続残差を照合する収支監査を行ってください。**

第 1 仮説: **冷却側には縮流部の連続残差が残っており、壁温ピンによる直接の密度消去より、擬似時間反復の未収束が欠損の主要因である。**　確度: **中**

- 根拠：温度ピンは `ρ` を保持し、DPLUR も壁の連続行を残します（[timeIntegration_d.cu:1029](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1029)）。これらの関連ソースには、申告されたビルド元 `e2696d8f0` と現 HEAD の差分がありません。ブリーフの欠損減少も整合しますが、その質量流量は手元で再計算できていません。
- 反証条件：**数値流束で評価すると欠損が消える**、または欠損を説明する密度の非保存変更・壁境界流出が検出されること。
- 「未収束」と「さらに回せばゼロまで減衰する」は別です。現在の観測から後者や必要反復数は保証できません。

第 2 仮説: **節点台形積分と数値流束の差が収支判定を汚している。**　確度: **中、寄与量は未確認**。境界数値流束は境界状態から `mdot = |S|ρ_R U_n,R` として作られます（[convectiveFlux_boundary_d.inc.cuh:192](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:192)）。内部断面にも再構成・数値散逸の寄与があります。

第 3 仮説: **EOS 床など、温度ピン以外の状態修正が非保存寄与を持つ。**　確度: **低、未確認**。密度床の処理は存在します（[dependentVariables_d.cu:76](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/dependentVariables_d.cu:76)）。今回作動した証拠はありません。

判別 A/B: **変更点は流量の評価方法だけ**にします。

- **対象**：`case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/res_100000.h5`。断熱対照は `case/45.isobutane_m6_d155/run_0181_ns_coldmesh_ad/res_100000.h5`。
- **A**：現在の節点値の台形積分。
- **B**：同じ状態・実効設定で評価した面数値流束。全領域と縮流部を覆う双対 CV 集合について、入口・出口・壁・内部切断面の流束と **符号付き `res_ro` の総和**を照合します。RMS は収支には使いません。断面位置と双対面集合の違いを明示し、軸対称の半径重み・2π・体積の規約を揃えます。
- **長さ**：必要なのは各状態の残差評価1回です。既存ダンプがなければ、呼び出し側で別の診断 run に同一メッシュ・ビット一致 restart を作り、**1 step の最初の評価だけ**を使えます。既存の `FORGE_DUMP_MASSFLUX` は最初の面流束を保存します（[convectiveFlux_d.cu:565](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:565)）。FP64 ビルドでは実装どおり `flow_float` の幅で読み、コメントの「float32」を信用しないでください。
- 同時に更新前後・BC 適用前後・EOS 適用前後の密度変更を集計し、**陰的補正で説明される変更と、その後の上書きを分離**します。

**事前に固定する解釈**：

- **A だけに約2%の欠損があり、B は0.1%以内** → 観測欠損の主因を未収束とする第1仮説を棄却し、流量評価の問題を支持。
- **B にも同程度の欠損があり、連続残差と閉合し、壁流出・密度上書きが無視できる** → 「壁から消えた質量」という説明を棄却し、未収束の離散収支を支持。ただし減衰完了は未証明。
- **閉合しない／状態上書きが説明量を持つ** → 単なる遅い過渡という説明を採用せず、その位相を原因候補にする。

閉合の監査許容は、例えば **入口流量の1e−6以下**と先に固定します。瞬時の離散恒等式の検査なので、未収束の保存場でも有効です。

**CFD 0 step で既存スナップショットからできること**は、総質量・縮流部質量・壁／近壁密度の時系列、現在の流量指標の再計算です。提示された5枚はその傾向確認に使えます。しかし、通常の場出力だけでは、各更新位相の上書き量や数値面流束は復元済みではありません。`run_0185` の3枚も、登録した5枚の準定常判定には不足します。

やらない方がよいこと: **短い CFL A/B の場の一致を定常解不変の証明にすること、relax も同時に変えること、等比外挿を延長完了の保証にすること、0.103%を丸めて0.1%合格とすること。** 将来 CFL 変更を評価する場合も、同じ場からの両系列が独立にゲートを満たした後で比較する必要があります。

呼び出し側の前提への異議: **「断熱側は収束済み」「下流はほぼ止まったので比較可能」「ΣCFL がコード間の擬似時間を表す」は受け入れません。**

保存 NPZ に `check_quasisteady.classify` を実行した結果は、各 run の **80000〜100000、全5枚、drift・osc 各0.1%**で次のとおりです。

| 根拠 run | 判定 |
|---|---|
| `case/45.isobutane_m6_d155/run_0181_ns_coldmesh_ad/` | δ_E の x=40／70／94：`STEADY`。保存済み本段収束判定：`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)` |
| `case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/` | δ_E：x=40 `TRANSIENT-UNSETTLED`、70／94 `STEADY`。Q_w：`DRIFTING`、窓内約2.5%。保存済み本段収束判定：`NOT CONVERGED`、`rms_roOmega … RISING (divergent)` |

R_NS の保存窓平均は x=40で **0.72578**、94で **0.79293**ですが、総合 **`判定不能 (ゲート不成立)`** を維持します。成果物は [V_c45.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/V_c45.json)、[gates_aws.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/gates_aws.json)、run 索引は [README.md:122](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:122) です。

既実施の座標 A/B は最大 **0.059756%**を確認できました。これを R_NS 差の主因候補へ戻しません。ただし、この試験は質量流束の台形積分と離散流束の同等性を検証したものではありません。

不足情報: **AWS の元 HDF5・実効設定・全残差系列・数値面流束・更新位相別の密度変更量**が手元にありません。したがって欠損の帰属と延長所要量は未確定です。

ファイル変更・forge 起動はしていません。**plan 未反映**。採用時は `plans/active/tooling-nozzle-isothermal-wall-chain.md` の **§5.1 #26、§6 V-c45** に反映してください。収支指標や判定方法を改訂する場合は、旧値・旧判定を残し、観測済みデータと変更理由を記録した別版にして、今回の再解析を事後解析と明記してください。
