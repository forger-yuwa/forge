# codex 諮問 (diagnose): gradient-scalar-lsq-result1-disposition

- **brief**: [`notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-result1-disposition.md`](../../notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-result1-disposition.md)
- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **date**: 2026-09-26
- **commit**: `1cfd5623` (feature/sern-design)
- **codex**: effort `high`, 3.5 min, rc=0
- **結論**: **まず M1 の残る誤合格を塞ぎ、その後は restart 手順を揃えた case/48 の実装直前版／現行 `gg` 比較を一度だけ行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**Critical 0／Major 5／Minor 2。7 件とも採用。ただし、親案の M2・M3 の結論づけは修正が必要です。**

確認中に共有ファイルへ修正が入りました。以下は修正後も確認した判断です。AWS の実測・VERDICT 原本は今回未確認であり、plan の数値を再計算した結果とは区別します。

| 指摘 | 採否・根拠・対案 |
|---|---|
| **M1・Major** | **採用。修正は未完了。** 修正版のメモリ上試験では「末尾ゼロ」「参照 Inf」は拒否されました。しかし、必須5列すべてについて**参照が定数 `1e-6`、対象が全点ゼロでも PASS**。`all-zero` を列判定より先に除外しています（[check_floor_ratio.py:55](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_floor_ratio.py:55)）。**参照で活動していた列が対象で全ゼロになった場合も判定不能として拒否**し、負の試験へ追加する。両側で非活動の列とは区別する。修正後、S2 の全床比判定を取り直す。 |
| **M2・Major** | **採用。ただし片方向の床比と「実装は床を変えていない」という結論は却下。** 親案は旧版床／現行床の上限しか見ないため、旧版が現行の `0.1` 倍でも通ります。実際、合成系列 `1e-6→1e-7` は PASS（[check_floor_ratio.py:67](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_floor_ratio.py:67)）。**双方向の床比と restart 手順を揃えた比較**に変更する。詳細は下記。case/48 の結果だけで、他ケースの収束ゲート未達は閉じない（[plan:102](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:102)）。 |
| **M3・Major** | **採用。ただし「作用素差と確定」「不成立なら nSub 30 で審査」は却下。** 保存収支未確認・反復十分性未確認は [plan:106](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:106) のとおり。`check_passive_budget --mode fct` は採用し、収支閉合・独立総量照合・低次／HO 残差まで通す（[docstring:10](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_passive_budget.py:10)）。倍増感度が大きければ **30 も十分とは言えず、精度差の判定を保留**する。 |
| **M4・Major** | **採用。現在の修正版も親案を満たしていない。** `len(now)>1` と `len(seen)<=1` だけで、測定対象 PID と照合していません（[s3_perf.py:47](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/s3_perf.py:47)）。**実際の solver 子プロセスを特定して GPU PID と照合**する。`Popen` の対象は `bash` なので、その PID を無条件に solver PID としない。対象の観測・正常終了も必須とし、競合検出試験後に両 case を取り直す。ポーリングだけで短時間の競合まで排除したとは記録しない。 |
| **M5・Major** | **採用。暫定修正は条件付きで可。** 明示 `gg→lsq` が `segments()` で **1 区間**になることをメモリ上で再確認しました（[stage_manifest.py:225](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:225)）。`YAML_HARD_PATHS` 追加で新規 manifest を分離する方針は可。ただし既存 manifest の保存済み key は自動修復されません（[stage_manifest.py:275](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:275)）。**既存記録は設定原本から再生成するか、判定区間を明示**する。YAML 解析不能時も黙って連結させない。#2g は「起票」だけでなく**完了が既定化の前提**。 |
| **m6・Minor** | **採用。親案でよい。** 不整合の指摘箇所は [レビュー:69](/home/sano/work/forge-sern-design/notes/reviews/2026-09-26-gradient-scalar-lsq-unification-result.md:69)。既定 `gg`／opt-in `lsq`／cell は GG 固定を明記する。非合併条件の `lsq` は「片側 LSQ」とする。文書修正の完了確認は今回の限定読取範囲外。 |
| **m7・Minor** | **採用。親案でよい。** 索引不足の記録は [レビュー:75](/home/sano/work/forge-sern-design/notes/reviews/2026-09-26-gradient-scalar-lsq-unification-result.md:75)、床移動の未解決は [plan:103](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:103)。run 索引・結果・成果物・状態を同期し、真因調査を独立した F 項目にする。ただし別項目化によって S2 未合格を解消した扱いにはしない。 |

M3 の事前規則は次で固定するのが妥当です。**これは後続試験の規則で、次に投入する試験は下記 M2 に絞ります。**

- `gg`・`lsq` それぞれで変更は `nSubIterDualTime: 15→30` のみ。**元の開始状態から同じ200物理 step**を回す。物理刻み・BDF/FCT 履歴・実効設定を揃える。
- 各比較量 \(q\) について、同じ物理時刻・同じ節点・固定した正規化で  
  \(D_q=\|q_{\mathrm{lsq},15}-q_{\mathrm{gg},15}\|\)、  
  \(E_q=\|q_{\mathrm{gg},30}-q_{\mathrm{gg},15}\|+\|q_{\mathrm{lsq},30}-q_{\mathrm{lsq},15}\|\)  
  を評価する。旧比較対象の `ro`・`roY`・`roUy` はすべて含め、保存された共通時刻の系列で確認する。
- **A：全対象で \(E_q<0.1D_q\)、かつ4本とも保存・有界性ゲートを通過** → 「15→30 の反復感度では差の大部分を説明できない」と記録する。作用素差を支持するが、正しい離散化差とは確定しない。
- **B：感度条件を満たさない** → 反復数への依存が無視できない。**30 の双子を確定値として採用せず、精度差の審査は保留**。
- `D_q=0`、ノイズ以下、記録不足、保存ゲート不成立は判定不能。旧上限超過は残す。A でも旧精度ゲートを自動合格にしない。

結論: **まず M1 の残る誤合格を塞ぎ、その後は restart 手順を揃えた case/48 の実装直前版／現行 `gg` 比較を一度だけ行う。**

第 1 仮説: **case/48 の旧床逸脱は実装直前版でも再現し、本 plan の commit 群による大幅な追加変化ではない。** 確度: **低・未確認**  
  根拠: S1 の短時間比較は PASS と記録されている一方、chi=0 でも `rms_roUy` の旧床比は `9.08`（[plan:103](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:103)）。これは比較する理由にはなるが、長時間の無回帰を証明しない。  
  反証条件: **入力・実効設定・restart 手順を揃え、静定を確認した比較で、下記の床比または物理差の上限を超える。**

第 2 仮説: **現行版の commit 群が `gg` の長時間挙動を変えた。** 確度: **低・未確認**。短時間の S1 だけでは除外できない。  
第 3 仮説: **連続48000 step と24000＋restart＋24000の違いが比較を交絡する。** 確度: **低・未確認**。親案にはこの手順差がある（[ブリーフ:25](/home/sano/work/forge-sern-design/notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-result1-disposition.md:25)）。

判別 A/B: **変更因子をバイナリだけにする。**

- 旧版 `36d8ba03` を、現行対照と同じ **24000＋restart＋24000 step** で実行する。出力先候補は `case/48.flat_plate_cooled_m4/run_0957_sglsq_s2_base36d8/`。番号の未使用を確認する。
- IC・メッシュ・BC・実効設定・ビルド条件・GPU・BLOCKSIZEを揃える。chi は「既定だから同じ」とせず、**両バイナリの実効値で確認**する。
- 全活動残差列の末尾20%平均を \(f_{\rm old},f_{\rm new}\) とし、**各列で \(2/3\le f_{\rm old}/f_{\rm new}\le1.5\)**。非活動列の不一致・入力不正は拒否する。
- `Cf`・`q_w`・`δ*`・`θ` の3 stationは、固定した末尾区間の平均について **旧版を分母とする相対差 ≤0.1%**。両者の系列 STEADY、既存物理ゲート、全残差の NaN/Inf・RISING なしも必須。時間変動が0.1%の判別を妨げる場合は判定不能とする。
- **A：全条件成立し、旧版も元の起点床ゲートを外れる** → 「この case・手順・閾値では大幅な追加変化を検出せず、旧床逸脱は実装直前版でも再現」とする。**「床を変えていない」「原因は過去の commit」とは断定しない。**
- **B：有効な比較で床比または物理差が上限超過** → 追加変化なしの仮説を棄却し、commit 群の調査へ戻す。
- 過渡・設定不一致・証拠不足は **判定不能**。自動延長や閾値緩和はしない。

やらない方がよいこと: **現行床への参照差し替えで S2 を PASS にすること、M3 の A を正常な作用素差の証明にすること、今回の診断だけで既定化へ進むこと。**

呼び出し側の前提への異議: **`check_floor_ratio --start run_0950+0952` は現行 CLI の連結指定ではありません。** `--start` は単一パスで、複数 target も個別判定です（[check_floor_ratio.py:110](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_floor_ratio.py:110)）。同一設定の区間を明示し、step の重複を処理した連結 CSV を別途用意してください。

不足情報: AWS の残差 CSV、判定区間付き `CONVERGENCE_VERDICT.txt`、準定常 VERDICT・比較系列、両バイナリの実効設定・ビルド識別、restart 手順、FCT 収支判定原本。

**plan 未反映。** 呼び出し側で §5.1・§6・§6.1 に採否と事前規則を記録してください。ファイル変更・`forge` 起動は行っていません。
