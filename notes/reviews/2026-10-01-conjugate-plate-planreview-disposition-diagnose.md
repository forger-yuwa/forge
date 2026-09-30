# codex 諮問 (diagnose): conjugate-plate-planreview-disposition

- **brief**: [`notes/reviews/briefs/2026-10-01-conjugate-plate-planreview-disposition.md`](../../notes/reviews/briefs/2026-10-01-conjugate-plate-planreview-disposition.md)
- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **date**: 2026-10-01
- **commit**: `b1573788` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 2.8 min, rc=0
- **結論**: **B1の採取仕様とrestart受入条件を上記に修正し、同一起点からの200 step A/B予備確認を次の一手にしてください。**
- **extra**: `notes/reviews/2026-10-01-boundary-cht-conjugate-flat-plate-plan.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 指摘 | 採否 | 根拠と必要な処置 |
|---|---|---|
| M1〔Major〕観測・判別仕様 | **採用。ただし今回の具体案は修正が必要** | 1010 step 間隔は `gcd(1010,50)=10` なので **5 位相だけ**を採ります。各20000 step区間も19～20点です。単位振幅・周期1010 stepの正弦波をこの間隔で抽出すると、標準偏差は数値実験で **1.53e−14**になりました。実際の周期が1010とは主張しませんが、「50の倍数でない」だけでは振幅判定の根拠になりません。下記の密な採取確認を先に行ってください。設定の根拠：[ブリーフ:7](/home/sano/work/forge-cht/notes/reviews/briefs/2026-10-01-conjugate-plate-planreview-disposition.md:7)。 |
| M2〔Major〕6本全合格 | **採用** | C1/C2 × n16/32/64すべてに前提ゲートと主判定を要求してください。メッシュ品質・参照解の自己検査・格子間差の減少・C2の固体効果と不確かさも残します。[発注元:150](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:150)。Aの限定閉鎖はCの不合格を免除する根拠になりません。限定結果への変更は別の判断です。 |
| m3〔Minor〕restart成立確認 | **採用。ただし更新番号の継続条件は却下** | 対象は `flux_avg: 1`。現行コードでは累積stepは復元しますが、`nUpdate` の復元は `SOLID/QBUF` がある分岐だけです。[conjugateWall.cpp:655](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:655)、[同:914](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:914)。**番号の継続ではなく、累積更新step・更新回数・壁温の変化**で確認してください。この診断のためにバイナリを修正する必要はありません。 |
| m4〔Minor〕丸めの分離 | **採用** | 診断残差だけでなく更新右辺も相殺を含みます。[conjugateWall.cpp:953](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:953)、[同:975](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:975)。実際の非一様保存温度を固定し、同じ荷重・作用素・分解条件で残差と `Δu` を比較してください。 |

M1 の具体値は、次のように修正してください。

- **採取間隔**：200 stepの予備確認では両側とも毎step出力。本段は **11 step間隔を候補**にします。連成50位相をすべて巡回し、60000 stepで定期出力5454組です。予備確認の毎step系列と11 stepに間引いた系列で、判別対象の標準偏差が10%以内に一致することを最低条件とし、満たさなければ本段へ進めません。短い予備確認だけで長周期の採取精度まで保証したとは扱わないでください。
- **節点集合・尺度**：提案の帯は採用できます。`run_0007_c1_n64/mesh.h5` の `MESH/COORD` と `BCONDS/*/iCells` から確認した期待数は、板の前縁・後縁が**各10点**、上流・下流slipが**各21点**、直上内部点もそれぞれ同数、第一内部点の高さは **2.5e−6 m**です。`CELLS/centCoords` ではy=0が **0点**なので使わないこと。`Uy` は「+y方向成分」と定義すればよく、漏れ流束とは呼ばないでください。圧力尺度 **709.275 Pa**、温度10 K、速度 `U∞` は妥当です。熱流束尺度は元runの固定した評価時刻・空間窓・積分方法を記載し、平均の**絶対値**を使用します。
- **振幅・区間**：20000 stepずつの算術平均、横ばい `max/min ≤ 1.1` は採用。振幅は**量ごと・節点群ごと**に標準偏差の節点最大値を求め、異なる量を一つの最大値にまとめないでください。「単調減少」は非増加か厳密減少かを明記し、既に測定床に達した系列を比だけで評価しないこと。
- **Aの再現条件〔Major〕**：「熱流束の振幅最大が前縁」は**要再検証**です。元の根拠は**界面残差**の最大位置であり、熱流束の時間標準偏差ではありません。[対象plan:30](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:30)。両者は別量で、界面残差には固体側の項も入ります。[conjugateWall.cpp:946](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:946)。対案は、残差水準の1/2～2倍条件に加え、**同じ界面残差の定義・区間で位置を比較**すること。熱流束は全界面で系列を保持してください。
- **正式ツールとの関係**：`--abs-scale 1 --drift 0.002 --osc 0.002` は補助診断として採用。ただし温度では0.02 K相当なので、報告された0.003 K程度の変動を `STEADY` とする可能性があります。B/A振幅比の判別を置き換えません。`--series-cols` を指定し、各20000 step区間を別CSVにして `--tail 1` とすれば判定区間も揃います。[check_quasisteady.py:278](/home/sano/work/forge-cht/solver_density_cuda/tools/check_quasisteady.py:278)。

結論: **B1の採取仕様とrestart受入条件を上記に修正し、同一起点からの200 step A/B予備確認を次の一手にしてください。**

第1仮説: 壁温更新を止めても、流体側に同程度の停滞が残る。 **確度: 中**  
　根拠: `case/65.conjugate_flat_plate/run_0006_c1_n32/` の `conjugate_history.csv` を再集計すると、累積step **596000～599950**の `dTw_max` 最大は **7.0047e−5 K**。一方、保存された `CONVERGENCE_CHECK.txt` は **`NOT CONVERGED (stalled/plateau)`**、末尾値は `rms_ro=5.04e−8`、`rms_roUy=6.69e−6`、`rms_roe=1.54e−2`です。小さい壁温更新でも停滞していることは支持材料ですが、連成から独立である証明ではありません。索引：[case README](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:9)。  
　反証条件: Aが元の停滞を再現し、Bの末尾2区間で3残差すべてがAの1/10以下になり、登録した局所変動も減衰すること。

第2仮説: 動的な連成が停滞を維持している。**確度: 中、未確認**。n64の大きな界面残差だけでは原因と結果を分離できません。  
第3仮説: restartまたは抽出の不整合がA/Bの差を作る。**確度: 低、未確認**。一次データと予備確認で先に除外します。

判別A/B: **変更は `conjugate.warmup` だけ**。Aは5000、Bは延長後も更新しない値、例えば900000。予定先は `case/65.conjugate_flat_plate/run_0011_ab_c1_n64_coupled/` と `case/65.conjugate_flat_plate/run_0012_ab_c1_n64_fixed_tw/`。予備確認後、本段は同じ元checkpointから各60000 step、必要なら120000 stepまで延長します。  
→ **Bが減衰する場合**：壁温固定の流体系だけで同程度の停滞を維持する仮説を棄却。  
→ **Bでも3残差がAの1/2以上かつ横ばいの場合**：動的連成だけが停滞を維持する仮説を棄却。slip原因の確定にはなりません。  
予備確認では初期壁温差≤1e−9 K、Bの壁温固定、Aの累積更新stepを確認します。出力名のstepは残差のstepより1大きく、出力は連成更新前です。[main.cpp:1808](/home/sano/work/forge-cht/solver_density_cuda/main.cpp:1808)、[同:2086](/home/sano/work/forge-cht/solver_density_cuda/main.cpp:2086)。

やらない方がよいこと: B1前のslip修正、`Df_scale`変更、④の閾値緩和、および熱流束振幅の位置を未確認のままrestart合否条件にすること。

呼び出し側の前提への異議: **観測事実は「界面残差が前縁に集中」であり、「熱流束の変動最大も前縁」は未検証の解釈**です。また、更新番号の継続は現行 `flux_avg: 1` の再開契約と一致しません。

不足情報: `run_0007_c1_n64/` の元残差CSV、`res_600000.h5`、`conjugate_state_5.h5`、全界面の時系列が手元にありません。今回の正式収束ツール再実行は **`NO residual_history.csv`**。保存済み評価窓内77節点の温度・熱流束系列は正式ツールで **`OVERALL: ALL STEADY`**を再確認しましたが、前後縁を含む全界面の証拠にはなりません。

ファイル変更なし。**plan未反映**です。
