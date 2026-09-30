# codex 諮問 (diagnose): conjugate-plate-sampling-redesign

- **brief**: [`notes/reviews/briefs/2026-10-01-conjugate-plate-sampling-redesign.md`](../../notes/reviews/briefs/2026-10-01-conjugate-plate-sampling-redesign.md)
- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **date**: 2026-10-01
- **commit**: `e2e3af06` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **旧 FAIL と本段停止を維持し、既存200 stepの生系列で「FP64値／現行probe相当の6桁出力」の比較を先に行う。**
- **extra**: `case/65.conjugate_flat_plate/ab_series.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 判定・重大度 | 根拠と対案 |
|---|---|---|
| 旧 FAIL を維持し、本段を止める | **採用** | [`AB_PRE_CHECK_run0013_0014.txt:62`](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/AB_PRE_CHECK_run0013_0014.txt:62) は `VERDICT: FAIL`。採取方式を改訂しても旧判定は残す。これは**観測計画の不合格**であり、連成原因の判定ではない。 |
| 案1: 現行 probe をそのまま使用 | **却下・Major** | [`point_probes.cu:207`](/home/sano/work/forge-cht/solver_density_cuda/probe/point_probes.cu:207) の新規 `ofstream` は精度指定がなく、既定の有効6桁。約10万 Paでは刻み約1 Pa、300 Kでは約0.001 Kとなる。記録された A の後縁変動は **P 約0.00627 Pa、T 約0.0000205 K**で、これより小さい（集計ファイル11–12行）。対案は **FP64のまま保存する小容量系列**。毎 step 採取自体は妥当。 |
| 案2: 101 step 間引きを2020 stepのスペクトル等で承認 | **却下・Major** | 101≡1 mod 50なので、2020 stepの20標本は50位相中20位相しか覆わない。「互いに素」は短窓での代表性を保証しない。更新周期は [`conjugateWall.cpp:1082`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1082)。対案は界面量も毎 step保存し、101 stepは全場の確認用出力にだけ使う。 |
| 案3: 標本数による95%区間へ変更 | **却下・Major** | 実装は単一系列の標本標準偏差ではなく、**固定位相で間引いた系列の、節点別標準偏差の最大**を比較している（[`ab_series.py:127`](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/ab_series.py:127)、175行）。独立正規標本を前提とする約17%の式を、そのまま使えない。区間が広くなって包含することも「10%以内の精度」の証明ではない。 |
| 200 stepの界面変動が近いことを原因除外に使う | **却下・Major** | 集計32・59行の0.1245／0.1206は短窓の振幅。Bには再開時の流体変動が残るため、壁温更新の長期的な寄与は除外できない。本段の登録済み3区間の減衰判定に使う。 |

結論: **旧 FAIL と本段停止を維持し、既存200 stepの生系列で「FP64値／現行probe相当の6桁出力」の比較を先に行う。**

第 1 仮説: 現行probeへの置換は、採取間隔の問題を出力の量子化誤差に置き換え、微小振幅の判別を壊す。 **確度: 高**
  
  根拠: 上表の `point_probes.cu:207–219` と後縁の実測振幅。これは**提案された採取方式の欠陥**であり、既存HDF5から得た FAIL の原因ではない。
  
  反証条件: 同じ生系列を有効6桁に丸めても、測定床を超える全登録系列で振幅差が10%以内に収まり、非ゼロ変動の消失も起きないこと。

第 2 仮説: 既存 FAIL は、短窓・固定開始位相による標準偏差推定の不安定さを含む。 **確度: 中**。ただし、生系列が手元にないため、標本数不足と周期的な採取偏りの寄与は未確認。

判別 A/B: **変更するのは数値の保存精度だけ**。既存の `run_0013_pre_c1_n64_coupled`／`run_0014_pre_c1_n64_fixed_tw` の同一200 stepについて、Aは元のFP64値、Bは現行probe相当の有効6桁に変換して読み戻した値とし、同一節点・同一stepの振幅、最大値を取る節点、変動消失を比較する。**追加計算は0 step**。
  
  → 差が10%を超える、または変動が消失すれば第1仮説を支持し、現行probe案を棄却する。全対象が反証条件を満たせば、この窓について第1仮説を棄却する。連成原因については、どちらの結果からも結論しない。

やらない方がよいこと: FAILを「標本数で説明できる」として合格へ読み替えること、基準を緩めて101 step採取へ進むこと。

  **本段の推奨は、全判別量の毎 step・FP64小容量保存に一本化する。** 登録集合は124流体節点、手元のn64メッシュの界面は161節点だった。`P/T/Uy` と全界面の `q/Tw` なら、60000 stepの数値本体は非圧縮でも **約0.333 GB/run**。座標・節点IDは一度だけ保存すればよい。
  
  界面採取は `warmup` に依存しない出力経路に置く必要がある。連成更新処理の中に置くと、Bでは [`conjugateWall.cpp:1081`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1081) の早期returnで採れない。新経路は既存壁ダンプと同じ診断式・出力時点を使い、短い同時出力で値を照合してから本段へ進む。

呼び出し側の前提への異議:

- **Minor:** 「164点」は登録内容と不一致。`2×(10+10+21+21)=124`点。
- **Minor:** 0.1245／0.1206の分母は各200 step窓の平均ではない。元runのstep 600000、空間評価窓の平均熱流束 **1379.36 W/m²**を共通尺度にしている（[`ab_series.py:112`](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/ab_series.py:112)）。
- probeの採取頻度と位置については除外できた。毎 step呼ばれ、nodeでは探索前に `centCoords` が節点座標へ置換される。KD-treeという説明は現行CUDA実装と違い、実際は全探索だが、それ自体は今回の障害ではない。
- 200 stepの観測は本段の注目量を維持する根拠にはなるが、判定基準を変える根拠にはならない。

不足情報: 手元には予備確認2本のディレクトリがなく、集計記録だけがある。生系列、実際の設定差分、AWSバイナリとのソース対応は独立確認できていない。収束・準定常性を新たに認定する材料もない。

**plan未反映**。読み取り専用の指示に従い変更していない。呼び出し側で採用する場合は、§4.1.1・§4.1.2に旧FAIL、事後改訂理由、採取精度の照合条件を残す。
