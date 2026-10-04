# codex 諮問 (diagnose): twophase-onoff-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-03-twophase-onoff-interpretation.md`](../../notes/reviews/briefs/2026-10-03-twophase-onoff-interpretation.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-03
- **commit**: `f8e761ba` (feature/species-transport)
- **codex**: effort `high`, 5.5 min, rc=0
- **結論**: **次の一手は、既存保存場を使った 0 step の旧／新拡散作用素 A/B と同じ領域の液・総水分収支の後処理に絞り、物理機構の確定はその結果まで保留する。**
- **extra**: `case/16.nozzle_wys/README.md`, `methods/condensation.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（実測値は依頼文・plan の記録に基づく。対象 run はローカルに無く、VERDICT・場の独立再検証は未実施）

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | 差を「物理的効果として確定」と記載 | **要再検証**。記録された `OVERALL: ALL STEADY` は指定した報告量の準定常性を支持する。一方、全残差は `STALLED`、独立監査は `NOT CONVERGED`。[plan:130](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:130) のユーザ了承済み受入条件は尊重するが、これは数値誤差や機構の検証完了を意味しない。**「現設定・現格子で受入条件を満たした ON−OFF 差」**として記録し、原因は仮説と明記する。 |
| **Major** | H3「総水分拡散は本来 ON/OFF で同じ」 | **却下**。OFF は ∇Y_w、ON は分子成分に ∇z_v、乱流成分に ∇Y_w を使う。[speciesTransport_d.cu:303](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:303)、[twoPhaseDiffusion_d.cuh:86](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:86)。壁側の総水分低下はこの変更と整合する方向であり、それ自体はバグの兆候ではない。対案は壁法線方向の分子・乱流流束の分離測定。 |
| **Major** | 「液移流流束 +3.24 %＝凝縮総量 +3.24 %」 | **要再検証**。[wflux.py:17](/home/sano/work/forge-species/notes/investigations/2026-10-03-twophase-onoff/wflux.py:17) は移流成分しか積分していない。拡散を含む境界流束、相変化ソース、離散残差の収支が必要。収支が閉じても言えるのはまず**正味凝縮＝凝縮−蒸発**の増加であり、正の凝縮ソース総量の増加とは区別する。 |
| **Major** | 再正規化は ON/OFF 同程度なので無害 | **却下**。OFF は再正規化で液・Q を変更しないが、ON は変更する。[speciesTransport_d.cu:194](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:194)、[同:2157](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2157)。最大係数偏差が同程度でも、作用する成分・場所・符号が違う。補正が主因という証拠も無いので、対案は液・水の局所収支不整合を観測差と比較すること。 |
| **Minor** | 出口分布・固定マスクの扱い | **要確認**。[bands.py:34](/home/sano/work/forge-species/notes/investigations/2026-10-03-twophase-onoff/bands.py:34) は最後の 0.3 mm を抽出し、壁距離だけで並べるため、単一断面・片側壁の分布を保証しない。また系列抽出の M0 は既定で各 run の restart 場から作るが、比較表は OFF の M0 を両者に使う。[twophase_ab_series.py:213](/home/sano/work/forge-species/case/16.nozzle_wys/twophase_ab_series.py:213)、[twophase_ab_compare.py:122](/home/sano/work/forge-species/case/16.nozzle_wys/twophase_ab_compare.py:122)。対案は断面を一意に選び、ON 系列にも run_0482 由来の同じ M0 を明示すること。 |

結論: **次の一手は、既存保存場を使った 0 step の旧／新拡散作用素 A/B と同じ領域の液・総水分収支の後処理に絞り、物理機構の確定はその結果まで保留する。**

第 1 仮説: **気相基準の分子拡散が蒸気を壁側から外側へ運び、液・モーメントの乱流拡散が液を壁側へ運ぶことで、境界層の組成・相変化分布が変わった。** 確度: **中**。流束モデルの差は確認済みだが、温度差・液流束差の定量的な因果は未確認。

根拠: OFF は総水分が一様なら、その分子・乱流拡散ともほぼゼロになる。ON は

`z_v = (Y_w − g)/(1 − g)`

を使うので、同じ一様 Y_w でも

`∇z_v = −(1 − Y_w)/(1 − g)² · ∇g`

となる。液が多い外側ほど蒸気組成が低く、分子拡散は蒸気を外側へ、液の乱流拡散は逆向きへ運び得る。これに対抗する総水分の乱流流束が生じるには、壁側の Y_w が低くなる方向が整合する。[twoPhaseDiffusion_d.cuh:97](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:97)、[同:114](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:114)。

したがって H3 は、単に「H2O の D が大きいから」ではなく、**拡散の駆動変数が総水分から気相蒸気組成へ変わったこと**で説明するのが適切。−20 %という大きさの妥当性までは、式だけでは保証できない。

反証条件: 壁側の総水分減少域で、新作用素による蒸気分子流束が外向きにならない、またはその寄与が局所の収支不整合に埋もれるなら、この説明を主因とする仮説は棄却する。

第 2 仮説: **液・モーメントの混合が、従来液滴の乏しかった過飽和層へ凝縮核・表面積を供給し、そこでの成長を増やす。** 確度: 中〜低。ブリーフの OFF で S>1・液ほぼゼロの帯と整合する。ソースは Q2 と成長速度に依存するため、g の増減だけでは判断できない。[methods/condensation.md:695](/home/sano/work/forge-species/methods/condensation.md:695)。最内層での蒸発と、その外側での追加凝縮は両立する。

第 3 仮説: **再正規化・モーメント射影と未解消の離散残差が、ON−OFF 差の一部を作っている。** 確度: 低〜中。主因の証拠は無いが、既存のゲート記録では除外できない。DPLUR の右辺は全残差であり、読んだ箇所にはソースを別前処理で分割する旧懸念の再発は見当たらない。[condensationTransport_d.cu:870](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:870)。

判別 A/B: **OFF の `res_16000.h5` を共通入力に固定し、診断上の拡散作用素だけ A＝OFF、B＝ON に変える。追加計算は 0 step、状態更新・再正規化は行わない。** 係数・幾何・境界処理は本番と揃え、壁距離帯ごとの総水分の分子／乱流流束、液流束、補正流束を比較する。

- **B で明瞭な蒸気外向き・液壁向き流束が現れ、A の総水分流束がほぼゼロ**なら、「OFF と ON は同じ総水分拡散のはず」を棄却し、第 1 仮説の輸送方向を支持する。
- **B でもその方向が現れない、または差が評価誤差以下**なら、第 1 仮説を棄却する。後から符号の説明を変えない。

これは輸送方向の判別であり、−11 K や +3.24 %全体の原因を証明する試験ではない。同じ後処理で追加する測定は次の三つまででよい。

1. **ON 自身の相変化の符号と大きさ**  
   固定した壁距離帯で `S_ON`、`condDrdt_0`、液ソース S_g の正負別体積積分を出す。利用可能なら既存診断量を使い、不足分は本番ソース関数で再評価する。`S_OFF<1` は ON の蒸発の証拠ではなく、`g_ON>0` もその場で凝縮した証拠ではない。`S_ON>1` だけでも、臨界半径などの成長条件を満たした証拠にはならない。[methods/condensation.md:700](/home/sano/work/forge-species/methods/condensation.md:700)。

2. **液・総水分の閉じた収支**  
   ON/OFF それぞれについて、移流＋拡散の境界流束、正負別の相変化ソース、残差を同じ制御体積で比較する。単なる ON/OFF 同一断面での流束一致を「保存」と呼ばない。事前基準案として、収支不整合の大きさが観測された液流束差の 10 %未満なら正味相変化量の議論に進み、超えれば定量解釈を保留する。

3. **報告対象と同じ定義の系列確認**  
   M0 を統一した系列、局所温度差、p95|ΔT|、片側液存在体積について保存済み系列と VERDICT を確認する。case README では中間 `res` は削除済みと記録されているため、新しい局所量の時系列を最終場一枚から保証しない。[case README:390](/home/sano/work/forge-species/case/16.nozzle_wys/README.md:390)。

やらない方がよいこと: **既定値変更、再正規化の撤去、機構を見ない長時間延長を先行させない。** また、「液が移動したから LΔg 相当だけ冷えた」と説明しない。液輸送の潜熱エネルギーと EOS の液量依存には相殺があり、蒸発による冷却とは区別する必要がある。[methods/condensation.md:1110](/home/sano/work/forge-species/methods/condensation.md:1110)。外縁の +2〜+4 K や壁圧上昇を、現時点で潜熱放出・熱的閉塞に一意に帰属させるのも早い。

呼び出し側の前提への異議:

- **「ON のみ凝縮の節点」**は、現状では「ON のみ液が閾値を超える節点」。輸送された液の存在と、局所の正の凝縮ソースを分ける。
- `g_exit_mw` と液移流流束の増加は独立した二つの証拠ではない。前者は後者を質量流束で割った量である。[twophase_ab_series.py:177](/home/sano/work/forge-species/case/16.nozzle_wys/twophase_ab_series.py:177)。単なる分母変化という H4 は提示数値と合わないが、相変化収支まで確定したことにはならない。
- 報告可能なのは、条件を付けた**数値比較値**。記録どおりなら、7 報告量の指定区間に対する `ALL STEADY` を添え、出口質量重み g 約 +3.2 %、平均壁温約 −2.4 K、壁圧差などを記載できる。`T_cond0_vmean_K` は M0 同一性の確認が条件。
- p95|ΔT|＝4.08 K、9576 節点、局所最大差、壁の −20 %、追加の断面流束は、今回提示された7量の定常判定に含まれないため、**最終保存場での値**と明記する。
- onset −0.0012 mm、pw21 +0.003 %を、有意な物理変化として扱う根拠は不足。各系列の変動幅との比較が必要。全列 `STALLED` を「収束した」、`STEADY` を「誤差が無害」と言い換えない。

不足情報: AWS の run 設定、判定出力・系列CSV、M0 指定、ON のソース診断・収支。根拠となる保存場は `case/16.nozzle_wys/run_0520_twophase_off_16k/res_16000.h5` と `case/16.nozzle_wys/run_0521_twophase_dplur_nn0_cont/res_16000.h5`、ON の前段は `case/16.nozzle_wys/run_0519_twophase_dplur_nn0_16k/`。本レビューでは forge を起動せず、ファイルも変更していない。**plan 未反映**。
