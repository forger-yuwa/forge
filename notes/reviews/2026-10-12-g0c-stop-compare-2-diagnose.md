# codex 諮問 (diagnose): g0c-stop-compare-2

- **brief**: [`notes/reviews/briefs/2026-10-12-g0c-stop-compare-2.md`](../../notes/reviews/briefs/2026-10-12-g0c-stop-compare-2.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `894e5a43` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.1 min, rc=0
- **結論**: **試験 G を「登録した停止規則での出力互換性試験」に限定して組み直し、G0c の上限 1 対 50 を比較する。ただし採用候補と G0 の直接比較を追加する。**
- **extra**: `notes/reviews/2026-10-12-g0c-stop-compare-diagnose.md`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/rg0c_judge.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **試験 G を「登録した停止規則での出力互換性試験」に限定して組み直し、G0c の上限 1 対 50 を比較する。ただし採用候補と G0 の直接比較を追加する。**

採否表（Critical 0 / Major 5 / Minor 1）:

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| Major | **停止規則での比較への変更：採用** | [plan §6.11:494](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:494) は、当初から定常解の同等性を保証しない運用試験を定義している。§5.1 #25 のユーザ決定も停止値との比較である。**この目的に対する前回の「両腕の収束合格を必須」は撤回する。** 収束・準定常の VERDICT は独立して残し、量の窓判定で残差判定を置き換えない。 |
| Major | **「G0 は構成上合格できない」：要再検証** | `STALLED` は低下桁数不足と末尾トレンドから付く分類であり、原因の診断ではない。[check_convergence.py:193](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_convergence.py:193)。**「確認した履歴では未合格」と訂正する。** 同一実効設定の継続履歴を適切に連結して判定するが、未合格参照の `--from-floor` 拒否は回避しない。 |
| Major | **G1 対 G0 を純粋な格子誤差と呼ぶ案：却下** | G0c の初期場は移植を含み、区間内の 44 節点は別の壁法線行から値を受け取っている。[plan:1212](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1212)。未収束では反復経路・停止時刻の影響も残る。**「格子変更と登録した初期化・停止手順を含む出力差」と呼ぶ。** 定常解の格子非依存性は判定しない。 |
| Major | **二段の相対比較だけで G2 を採用する案：却下** | [ブリーフ:18](/home/sano/work/forge-integ-1005/notes/reviews/briefs/2026-10-12-g0c-stop-compare-2.md:18) の G1/G0 と G2/G1 が各 +0.1 % なら、G2/G0 は +0.2001 %。**G2a・G2b を G0 の FP64 二本それぞれと直接比較する。** 出口 M・収支・せん断の許容値、比較座標、局所量のゼロ近傍の絶対許容も事前登録する。未定義の指標は合格扱いにしない。 |
| Major | **代理物性による比較：限定採用。エネルギー収支への流用：却下** | q_w は Sutherland・固定 CP/Pr、h0 未指定時は定数比熱の式になる。[cold_xcheck.py:95](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:95)、[同:104](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:104)。**代理指標として明記して比較することは可。** 実熱流束との相対差が保存されるには係数比の共通性の確認が必要。TP の収支には `VALUE/h0` と実効物性・離散化に整合した熱流束を用いる。符号付きせん断の追加は採用。 |
| Minor | **数値収集は三本並行、速度は占有：採用** | [plan §6.11:512](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:512) と整合する。総時間は共通起点から数え、起動・出力費用も含める。G2 二本の幅は観測した再実行差であり、信頼区間ではない。G1 一本では上限 1 側の再現幅は未評価と残す。 |

第 1 仮説: **G0c でも、上限による反復経路の違いが停止時の出力に残り、0.1 % の互換性判定に影響する。** 確度: 中。ただし G0c の停止時については未確認。  
根拠: 上限は平均流・SST の擬似時間刻みを変える。[solver-settings.md:75](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:75)。停止判定は θ_r 三断面の両端変化と質量欠損だけで、残る反復誤差を制限していない。[t2_judge.py:17](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/t2_judge.py:17)。実際に `case/45.isobutane_m6_d155/run_0754_t2_d2/` は停止時にも θ_r が +0.038〜+0.045 %/2 万 step 動いている。保存された判定は **`OVERALL: CHECK FAILURES ABOVE`**、再開後 122,500〜142,500 step の θ_r 全断面・Q_w は **`DRIFTING`**。[t2_judge.json:82](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/t2_judge.json:82)。  
反証条件: 同じ G0c・起点・停止規則で、上限 1 と 50 の停止出力が、事前登録した全指標の許容内に収まること。この場合に退けるのは「今回の有限停止出力に許容超過が残る」という予測だけである。

判別 A/B: **変更する設定は `time.deltaT.lineDtDirectionalCap: 1 / 50` の一点。**

- 提案の G1・G2a・G2b を使う。全腕 FP64、同じ `_g0c/q32_init.h5`、同じバイナリ・環境変数・その他の実効設定。入力と初期化後の九保存量を照合する。
- 各腕は §6.11 の最初の到達で停止し、共通起点から最大 200,000 step。2,500 step 間隔、末尾九出力を保存する。比較対象に近づけるための延長はしない。
- 主比較は G2a/G1・G2b/G1。θ_r・δ_loc・代理 Q_w は 0.1 %、代理 q_w は全区間 1 %を暫定許容とする。出口分布・収支などは上表の定義を埋める。全残差列と各量の VERDICT、判定区間を残す。
- **結果 A**：全腕が到達し、全指標が許容内 → この起点・停止手順での「上限変更による許容超過」を支持しない。
- **結果 B**：到達した腕で許容超過、または予算内に未到達 → 上限 50 が今回の運用要件を満たすという案を退ける。再実行間で採否が割れれば再現性未確認、欠損・非有限なら判別不能／発散として分ける。
- G0c の採用判断には、さらに **G1・G2a・G2b と `run_0752_t2_d1/`・`run_0754_t2_d2/` の直接比較**を使う。A だけでは G0 との互換性を認めない。

やらない方がよいこと: `--drop` を事後的に緩めること、`STEADY` を残差合格の代用にすること、G0c の密度残差の低下桁数だけで「停滞しない」と結論すること。T2 の結果 B も保持する。

呼び出し側の前提への異議: **運用上の停止出力比較は成立し得るが、それを収束解の比較に読み替えることはできない。** 今回は比較対象の変更であり、収束の定義変更ではない。

不足情報: 対象五 run の実ディレクトリ・残差 CSV・HDF5 はこのワークスペースに存在しないため、今回確認できたのはコードと保存済み判定 JSON。G0c の不合格列、最新の全列 `STALLED`、`--from-floor` 拒否の実出力は独立再検証できていない。**ファイル変更なし・plan 未反映。反映先は呼び出し側の plan §5.1 #25 と試験 G の事前登録節。**
