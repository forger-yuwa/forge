# codex 諮問 (diagnose): conjugate-plate-fluid-stall-next

- **brief**: [`notes/reviews/briefs/2026-10-01-conjugate-plate-fluid-stall-next.md`](../../notes/reviews/briefs/2026-10-01-conjugate-plate-fluid-stall-next.md)
- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **date**: 2026-10-01
- **commit**: `f9579878` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **壁温固定の同一起点から、`time.deltaT.cfl_pseudo` だけを 2 → 0.5 に下げる A/B を行い、有限の擬似時間刻みが停滞を維持しているかを先に判別する。**
- **extra**: `case/65.conjugate_flat_plate/ab_judge.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 「連成は停滞の原因ではない」 | **却下**。B1 が支持するのは「この再開状態では、壁温の動的更新を止めても停滞が持続する」。形成過程への連成の寄与や、固定した非一様壁温の影響までは除外していない。[plan:102](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:102) の最後の断定を、この限定文へ置き換える。 |
| Major | `limiter: 0` を「1 次」とする案 | **却下**。実装はリミッタ配列を **1 に充填して戻る＝無制限再構成**。[limiter_d.cu:450](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/limiter_d.cu:450)。1 次化は `space.convMethod: 0`。ただし散逸も変わるため、改善してもチャタリングの特定にはならない。今回は採らない。 |
| Major | 前縁の熱流束変動が G-if ①の主成分 | **要再検証**。毎 step の標準偏差と、固体更新時の界面不釣合い最大値は別の量。[ab_judge.py:59](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/ab_judge.py:59)。さらに `q_eff` は質量残差による蓄積補正を含む半離散式の推定値。[methods/boundary.md:640](/home/sano/work/forge-cht/methods/boundary.md:640)。同一節点・同一更新時刻で固体側熱量との差を再構成するまで「主成分」としない。 |
| Minor | 「前縁の熱流束は 11 % 揺れる」 | **尺度を明記して採用**。`run_0012` の末尾 20000 step を再計算すると、最大標準偏差は **158.287 W/m²**、同じ節点の平均は **3033.749 W/m²**。1379 W/m² 基準なら約11.5 %、局所平均基準なら約5.22 %。「評価窓基準の規格化標準偏差」と書く。 |

結論: **壁温固定の同一起点から、`time.deltaT.cfl_pseudo` だけを 2 → 0.5 に下げる A/B を行い、有限の擬似時間刻みが停滞を維持しているかを先に判別する。**

第 1 仮説: 有限の擬似時間刻みと現在の更新作用素の組合せが、減衰しない反復振動を維持している。**確度: 中**。ブリーフの H2 を、内部反復不足とはまだ特定せず採る。  
  根拠: `case/65.conjugate_flat_plate/run_0012_ab_c1_n64_fixed_tw/ab_series.npz` の60000点を再計算し、壁温の最大時間変化 **0 K**、上流 `P` の卓越周期 **7.2446 step**・単一周波数のパワー比 **91.83 %**を確認した。これは物理時間の周期ではない。`nStepInner` は固定残差に対する Jacobi sweep 数であるため、この周期だけから「内部反復不足」とは言えない。[main.cpp:1509](/home/sano/work/forge-cht/solver_density_cuda/main.cpp:1509)  
  反証条件: 下記の時間幅を合わせた比較で、CFL 0.5 でも残差・主要変動が対照の半分以上残り、3区間で横ばいなら、**「CFL 2 から0.5への縮小で主要な停滞を解消できる」仮説を棄却**する。時間積分全般の無罪までは意味しない。

第 2 仮説: 前後縁付近の再構成・リミッタの非線形変動。**確度: 低、未確認**。現在の系列にリミッタ値がなく、チャタリングの直接証拠がない。ソース検索でも既存のリミッタ凍結機能は見つからなかった。

第 3 仮説: 既知の node slip 欠陥。**確度: 低、今回との同一性は未確認**。既知症状は接線速度の空間的な符号交替であり、現在の毎 step 抽出は `P/T/Uy` のみで `Ux` がない。[調査記録:9](/home/sano/work/forge-cht/notes/investigations/node-slip-tangential-density-spurious-flow.md:9)、[ab_extract.py:74](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/ab_extract.py:74)。低マッハ性は第1仮説を増幅する背景として残すが、現段階で独立した真因には格上げしない。

判別 A/B:

- **対照A**: B1 と同じ累積600000 step の再開状態・固定壁温、`cfl_pseudo: 2`。**12000 step＝4000×3区間**。
- **試験B**: 同じ状態から `cfl_pseudo: 0.5` のみ変更。**48000 step＝16000×3区間**。他の数値設定・バイナリは共通。CFL×step を揃え、刻み縮小による単なる進行遅延を避ける。ただし局所刻みなので、代表節点の累積 `dt_local` でも時間幅を確認する。
- 全保存量残差、B1 の毎 step 観測量を測る。両側に同じ観測追加として、slip 上と直上の **`Ux`・`ro`・リミッタ値**を採り、接線速度の交互成分も見る。
- **Aだけで主要変動が維持され、Bの最後の2区間で能動残差全列と主要変動がAの1/10以下、かつ減衰傾向** → 有限刻み依存を支持。「空間的欠陥だけで同程度の振動が不可避」を棄却する。H1・H3の潜在的欠陥までは除外しない。
- **Bでも半分以上残り、3区間平均・振幅のmax/min≤1.1** → 上記のCFL縮小仮説を棄却する。**H3確定とはしない**。
- 中間・発散・対照の再現失敗は判定不能。短い診断から収束を宣言しない。準定常判定を併記し、収束は別途 `check_convergence.py` で判定する。

`ab_extract.py` の方式は使えるが、`ab_judge.py` は60000点・20000 step区間・連成Aの再現条件を固定しているため、そのまま流用できない。判定条件の変更を先に登録する。

やらない方がよいこと: `implicitRelax` とCFLの同時変更、slip→no-slipで原因を確定すること、1次化で改善しただけでH1を採用すること。**既存キー一つの変更でH3だけを無効化する、検証済みの対照は確認できない**。CFLのA/Bをその代用品とは呼ばない。

呼び出し側の前提への異議: **「どちらの結果でもH1/H3のどちらかが必ず消える」という二分法は成立しない**。境界閉包・リミッタ・時間更新は相互作用する。今回排除できるのは上記の限定仮説である。

また、`run_0012` の前縁最大変動節点について、`check_quasisteady.py` の `classify_series` を3区間へ適用した結果はすべて **`OSCILLATING`**だった（登録尺度1379、`tail=1`、`drift=osc=0.002`）。これは「固定壁温でも診断熱流束が振動する」根拠であり、物理的な熱流束変動やG-ifの原因確定ではない。

Cの完了条件は **6本全合格を維持する判断を採用**する。[plan:109](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:109)。停滞が解けなければ「前提ゲート不合格、Cは判定不能・未完了」とする。事後の窓除外・時間平均化・閾値緩和で完了させず、次の診断を別途登録する。

不足情報: 手元には対象3本の `residual_history.csv`、B1の全場HDF5、元runの尺度算出用壁HDF5がない。収束ツールの再実行は **`NO residual_history.csv`**、B1判定器は **`REFUSED`**となった。保存済み `CONVERGENCE_CHECK.txt` はA/Bとも **`NOT CONVERGED (stalled/plateau)`**だが、今回独立に再現できたのはNPZの系列統計とYAML差分である。H3の照合には `Ux/ro` の空間分布、G-ifの内訳には同期した界面ログが必要。

**plan未反映**。ファイルは変更していない。採用時は呼び出し側が§4・§5.1・§6へ反映する。
