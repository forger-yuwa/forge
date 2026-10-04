# codex 諮問 (diagnose): condensation-diffusion-error

- **brief**: [`notes/reviews/briefs/2026-09-27-condensation-diffusion-error.md`](../../notes/reviews/briefs/2026-09-27-condensation-diffusion-error.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-09-27
- **commit**: `65af7812` (feature/gap-heating-precision)
- **codex**: effort `xhigh`, 8.2 min, rc=0
- **結論**: **`run_0482` の同一保存量から、拡散作用素だけを「現行／保存的な二相輸送」に切り替える A/B を行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 流束比から「数 K 動く」とする解釈：**却下** | 二相 EOS は `e=e_gas+g(R_wT−L)`（[condensationEOS_d.cuh:381](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:381)）。液拡散による保存量変化率を `r_g` とすると、局所的に L 一定の乱流輸送では追加エネルギー項は `r_E=−Lr_g`。温度への作用は `ρc_v,eff δṪ=r_E+(L−R_wT)r_g=−R_wT r_g` となり、**潜熱部分が相殺される**。したがって `L·g/cp≈26 K` は今回の温度変化予測にならない。対案は、追加する全保存量残差を EOS に投影した温度変化率を測ること。 |
| **Major** | 分子拡散の絶対量：**要再検証** | 後処理の [analyze_liquid_diffusion_error.py:70](/home/sano/work/forge/case/16.nozzle_wys/analyze_liquid_diffusion_error.py:70) は係数 `1.8583e−7` にさらに `1e−4` を掛けている。ソルバは `1.8583e−3×1e−4`（[thermo_d.cuh:563](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:563)）。同じ入力で後処理の D は **10⁴ 倍小さい**。対案はソルバと同じ係数評価を使うこと。なお、この誤りは `E1/F1` では相殺され、乱流流束比にも影響しない。 |
| **Major** | §4.2 の補正速度・保存性：**要再検証** | [plan:63](/home/sano/work/forge/plans/active/condensation-two-phase-transport.md:63) の「気相の流束で ΣJ=0」を乱流にも適用すると、液流束を加えた全相の和はゼロにならない。対案は、**分子流束の補正と、全相で相殺する共通 `Sc_t` の乱流流束を分離**し、面ごとに `ΣJ_gas+J_l=0`、`J_w=J_v+J_l` を満たすこと。エネルギーは `Σh_sJ_s(gas)+h_lJ_l`。 |
| **Major** | 有限 `Sc_l` でも蒸気の非負性が構造的に保たれるとの説明：**却下** | [plan:61](/home/sano/work/forge/plans/active/condensation-two-phase-transport.md:61) の式では、液の分子流束が総水分に入らず、差 `J_v=J_w−J_l` に `+(μ/Sc_l)∇g` が残る。独立した蒸気拡散方程式ではない。対案は液流束を総水分・全相収支にも整合させること。**opt-in の残置方針は維持し、今回の A/B では液分子拡散をゼロに固定する。** |

結論: **`run_0482` の同一保存量から、拡散作用素だけを「現行／保存的な二相輸送」に切り替える A/B を行う。**

第 1 仮説: 保存的に二相拡散を入れた場合、凝縮域の定常温度差は小さい――判別用に、共通凝縮領域での体積重み `p95|ΔT| < 1 K` を仮置きする。確度: **低**。  
  根拠: 上記 EOS と流束の組合せでは、潜熱による直接温度作用が相殺される。局所定係数近似なら残る係数は L ではなく `R_wT` で、T≈210 K では提示された L の約4%。ただし、これは**局所的な式の見積もりであり、定常場の実測ではない**。1 K は今回の事前判別閾値である。  
  反証条件: 両 run が収束・定常性の判定を通り、末尾変動を十分上回る `p95|ΔT| ≥ 1 K` が残ること。

第 2・第 3 仮説: 第2＝液・モーメントの再分布と蒸気の分子拡散が凝縮ソースを変え、上記相殺があっても定常温度差を 1 K 以上にする。確度: **低・未確認**。第3は置かない。

判別 A/B: **両方が判定を通った後、温度差が1 K未満なら第2仮説を棄却、1 K以上なら第1仮説を棄却する。**

- 共通初期場は `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/` の最終保存量。同一メッシュなので `restart_field.py` を使い、新しい2つの run に複製する。同じバイナリ・S3・BC・CFL・物性設定を使う。
- **変更点は拡散モデルの切替1点**。A＝現行。B＝蒸気の分子＋乱流拡散、液および `Q0/Q1/Q2` の同じ `Sc_t` による乱流輸送、整合したエンタルピー流束を一括で有効化する。`h_l=h_v−L(T)` は既存 EOS と同じ基準を使う。これは試験用実装が必要で、既存設定だけでは成立しない。
- 最初は **4000 step、200 step ごとに保存**。未達なら4000 stepずつ延長する。比較は各モデルを固定した区間で `check_convergence.py` が PASS、報告量を含む `check_quasisteady.py` が STEADY になってから行う。A の床を B の収束基準として流用しない。
- 測り直す量は **`δṪ = (∂T/∂U)·(R_B−R_A)/V` の1つ**。共通初期場で、実際の双対面流束・補正・境界処理を使って求める。エネルギーだけでなく、総水分と液の残差も EOS に通す。これを「熱流束の発散との比」に置き換えない。
- 定常比較では温度差に加え、同じ抽出方法による onset、出口 g、壁圧を記録する。壁熱流束を報告するなら、その時系列も別途定常判定に含める。温度差が小さくても onset や g が小変化とは限らない。

  初期場の選択根拠は、既存記録の `run_0482: PASS`、`run_0483: NOT CONVERGED`。後者は `rms_roe` が2.9桁低下で停滞している（各 [0482 判定](/home/sano/work/forge/case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/CONVERGENCE_VERDICT.txt:2)、[0483 判定](/home/sano/work/forge/case/16.nozzle_wys/run_0483_passive_wys_s0_sfr0/CONVERGENCE_VERDICT.txt:2)）。

やらない方がよいこと: **`rog` を拡散させず、エネルギー式だけに潜熱流束を足すこと。** 上記の相殺を壊し、測りたいモデルとは異なる温度応答を作る。`Sc_t` 全体の変更や SFR の変更も今回の対照にならない。

呼び出し側の前提への異議: 「現行では拡散で動く水は蒸気として h_v を運ぶ」「総水分の乱流流束を維持したときの追加流束は `LρD_t∇g`」は採用する。ただし、**大きな流束比から温度影響を読む H1 と、発散比に替えればよいという H2 の枠組みは採用しない**。後処理は勾配の大きさしか残さず（[解析スクリプト:63](/home/sano/work/forge/case/16.nozzle_wys/analyze_liquid_diffusion_error.py:63)）、実装の面差分流束とも異なる（[speciesTransport_d.cu:263](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:263)）。`E1/F1=1` 単独では、補正後の現行水流束がゼロである証明にもならない。

  **Minor:** plan の「S3 面正規化は未確認」は更新すべき。正規化は [convectiveFlux_slau_d.inc.cuh:364](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:364) に存在する。

不足情報: EOS に投影した追加残差と、二相拡散を有効にした応答の実測。既存の `ALL STEADY` は4000～48000 stepの12点による onset・壁圧偏差・出口量の判定であり、今回の流束比や局所温度場の定常性を保証しない（[判定対象](/home/sano/work/forge/case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/QUASISTEADY_SERIES_VERDICT.txt:2)）。禁止された run データは読んでおらず、forge 実行・ファイル変更もしていない。**plan 未反映。呼び出し側で §4.2・§6 に反映する内容である。**
