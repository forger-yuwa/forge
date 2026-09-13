# codex レビュー: condensation-air (plan)

- **plan**: [`plans/active/condensation-air.md`](../../plans/active/condensation-air.md)
- **stage**: `plan`
- **date**: 2026-09-12
- **commit**: `d01a4118` (feature/condensation-air)
- **codex**: effort `high`, 7.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M8/m1
- **focus**: §4.1 n2_latent の低温線形化 (c_l=2 kJ/kg/K) の妥当性と Arthur 検証への影響、§4.2 pure CPG への二相音速と境界 (inlet_uniformVelocity / outlet_statPress 超音速) の扱い、§4.3 空気擬似種の過飽和度定義 2 案 (N2 分圧基準 vs 理想溶液露点線; 露点線は N2 線より 4.5 K 高温側で Daum & Gyarmathy の air≈N2 と合わない可能性) と N2 核生成・Iland 補正の流用、§6 の Daum & Gyarmathy onset 比較 (Ṗ を揃えた理論線 ±3 K) の合否条件が検証として適切か
- **extra**: `notes/investigations/condensation-carrier-kantrowitz-air-survey.md`, `case/34.arthur_n2_nozzle/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

目的は妥当で、既存の accepted plan とも重複しません。しかし、空気モデルの物質収支、CPG の境界・陰解法との整合、検証基準に実装前の修正が必要です。
特に「Arthur 実験との一致」「理論 onset 線との差」は、リポジトリ内の数値が計画の説明を支持していません。ファイルは変更していません。

1. **Critical — N2 選択凝縮と「pure 空気」の EOS が両立していない。**

   **根拠:** [plan:67](/home/sano/work/forge-cond/plans/active/condensation-air.md:67) は \(R_{\mathrm{eff}}=(1-g)R_{\mathrm{air}}\) とし、[plan:83](/home/sano/work/forge-cond/plans/active/condensation-air.md:83) は O2 をキャリア、凝縮相を N2 としています。後者なら必要なのは、
   \[
   p_{N_2}=\rho(Y_{N_2}-g)R_{N_2}T,\qquad
   R_{\mathrm{eff}}=R_{\mathrm{mix}}-gR_{N_2}
   \]
   です。凝縮が進めば気相の \(y_{N_2}\) も変わります。現在の [cond_vapor_state:19](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cu:19) の pure 経路は、全気相を凝縮可能として数えます。

   また、核生成率は \(S\) だけでなく `R`、分子質量、蒸気密度に依存します（[cond_nucleation:75](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cuh:75)）。これらを空気の値に替えて Iland 係数だけ残しても「N2 核生成」にはなりません。§5 には、SLAU の物性選択が `condModel != 1` を N2 に落とす箇所も抜けています（[流束:347](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:347)）。

   **対案:** 初版を **CPG の N2 選択凝縮＋非凝縮 O2 キャリア**に固定してください。`g` は総混合物に対する液体 N2 質量分率とし、EOS・核生成・成長・蒸発・枯渇上限・流束で同じ定義を使う。これは `nCondSpecies: 2` を必要としません。露点線モデルは別の混合液モデルとして後続に分けるべきです。

2. **Major — 潜熱修正は飽和圧も変える。現在の R2 は「潜熱だけ」の A/B にならない。**

   **根拠:** [n2_psat:108](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:108) は低温 C–C 外挿の潜熱に `n2_latent(50)` を直接使います。コードの係数で再計算すると、旧値は **203.94 kJ/kg**、提案値は **228.59 kJ/kg**。ここを連動変更すると、38 K の飽和圧は旧値の **0.592 倍**になります。核生成の駆動力も大きく変わります。

   また、30 K の潜熱増加は **32.2%**で、[plan:62](/home/sano/work/forge-cond/plans/active/condensation-air.md:62) の「+5〜+25%」を超えます。\(c_l=2000\) は低温側の正の熱容量を保証する有用な閉包ですが、63–77 K の測定値を30 Kまで延長する物性精度は、壁圧への適合だけでは検証できません。

   **対案:** 新モデルでは同じ \(L(T)\) を使って
   \[
   \ln\frac{p_{\mathrm{sat}}(T)}{p_{\mathrm{sat}}(T_s)}
   =\int_{T_s}^{T}\frac{L(\tau)}{R\tau^2}\,d\tau
   \]
   を低温外挿に適用し、接続点と微分の整合を検査する。旧物性一式を回帰用に保持し、潜熱のみ変更する診断 run と、飽和圧まで整合させる最終モデルを区別してください。\(c_l\) の感度評価も必要です。

3. **Major — 「CPG では `sonic` のみ変更」は block-DPLUR の既存構造と整合しない。**

   **根拠:** [timeIntegration_d.cu:773](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/timeIntegration_d.cu:773) は CPG でも `gamma_arr` を読みます。ただし [gasProperties_d.cu:47](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/gasProperties_d.cu:47) が config γ を設定し、[呼出箇所:1392](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1392) は一般 EOS 分岐を TP のみに限定しています。CPG 固有系は実際の `Ht` を使わず、\(\chi=0\) を仮定します（[Jacobian:49](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:49)）。

   提案物性の \(T=45\) K、\(g=0.1\) では、固定 \(g\) の EOS から **\(\gamma_{2\phi}=1.30781\)、\(\kappa=0.30781\)、\(\chi=8516\ \mathrm{J/kg}\)**。config γ=1.4、\(\chi=0\) とは異なります。これは直ちに定常解の誤りを意味しませんが、整合した frozen Jacobian ではありません。

   **対案:** CPG 二相でも局所 \(\gamma_{2\phi}\)・実 `Ht`・一般 EOS 固有系を使い、`gasProperties` による上書きも解消する。既存の実 Jacobian と流束有限差分の比較を N2/空気へ拡張し、double と float32 の両方で確認してください。音速 on/off の場差を「観測」するだけでは不足です。

4. **Major — `slip` が既に二相非整合。出口の事後確認だけで自動 ON にできない。**

   **根拠:** CPG `slip` は ghost を
   `roe=P/(γ−1)+ρek`、`T=P/(ρR)`、`sonic=sqrt(γP/ρ)`
   で再構成します（[boundaryCond_d.cu:101](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/boundaryCond_d.cu:101)）。内部が \(p=(1-g)\rho RT\) なら、ghost 温度は内部温度の \(1-g\) 倍になります。`g` は別途コピーされるため、二相 EOS と矛盾します。Arthur は側壁・対称面・面外両面すべてがこの境界です。

   `outlet_statPress` の超音速分岐が内部値をコピーする点は正しい。しかし、その判定は **法線 Mach ではなく速度絶対値の Mach**です（[境界:616](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/boundaryCond_d.cu:616)）。最終時刻だけの亜音速セル数では、起動中の不整合も検出できません。

   **対案:** `slip` は内部の熱力学状態を保持し、速度反射に伴う運動エネルギーだけ更新する。出口は全実行期間の \(u_n/c\) を監視し、未対応状態への遷移を検出する。これらの境界試験が通るまで自動 ON を広げないでください。

5. **Major — 理論 onset の参照値が本文と CSV で約4.5 K違い、±3 Kの合否が逆転する。**

   **根拠:** [理論線 CSV:5](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/daum_gyarmathy_theory_onset_Pdot20000_n2.csv:5) は100 Paで36.5 K、1000 Paで46 Kです。対数圧力補間すると700 Paで **44.53 K**となり、[plan:132](/home/sano/work/forge-cond/plans/active/condensation-air.md:132) の「約40 K」と一致しません。

   `case/34.arthur_n2_nozzle/run_0008_ref_n2/res_8000.h5` に既存の onset 抽出を再実行すると、**678.3 Pa、37.85 K、\(\dot P=1.73\times10^4\ \mathrm{s^{-1}}\)**。CSV の理論値44.40 Kとの差は **−6.55 K**です。これは後述のとおり定常性未確認の抽出値ですが、少なくとも「参照 run は理論線より約1.5 K低い」という記述は成立しません。

   \(\dot P\) を考慮する方針自体は妥当です。一方、Grossir は最小 onset 曲線を最大過冷却の限界として説明しており、「低膨張率の風洞を含む包絡だから」という説明も修正が必要です。[Grossir & Rambaud, §III.D](https://dipot.ulb.ac.be/dspace/bitstream/2013/208925/3/AIAA_2014_1153_Grossir.pdf)

   **対案:** 原図の線種・単位・読み取り点を再確認して参照データを一本化する。原図そのものは今回取得できず、どちらの読みが正しいかまでは確認できません。理論線±3 Kはモデル間比較の基準とし、空気モデルの実験検証合格とは分ける。同じ一点でモデルを選んで合否まで判定せず、独立した複数条件を用意してください。

6. **Major — Arthur の「実験との1–2%一致」は、保存された比較表では Lin の計算曲線との比較。**

   **根拠:** [arthur_fig2_digitized.csv:2](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/arthur_fig2_digitized.csv:2) は `cond` を計算曲線、`exp` を Arthur 実験記号と区別しています。3/4/5 in の比は、

   | 比較対象／dry | 3 in | 4 in | 5 in |
   |---|---:|---:|---:|
   | Lin 計算曲線 | 1.200 | 1.308 | 1.447 |
   | Arthur 実験記号 | 1.133 | 1.250 | 1.500 |

   計画が引用する1.20/1.31/1.45は前者です。したがって、これを根拠に「実験±5%を維持」とは言えません。

   **対案:** 実験の `P/P0` と直接比較し、Lin 計算曲線との差は別に報告する。読み取り誤差、dry 側の誤差、onset 周辺を含めた評価点を先に固定してください。潜熱モデルを旧計算曲線への適合で選び直すべきではありません。

7. **Major — 現在の保存頻度と後処理では、§6の準定常合格を実行できない。**

   **根拠:** ツールを再実行した結果は次のとおりです。

   | run（`case/34.arthur_n2_nozzle/` 配下） | `check_convergence.py` | `check_quasisteady.py` |
   |---|---|---|
   | `run_0008_ref_n2/` | `NOT CONVERGED (stalled/plateau)` | `TRANSIENT-UNSETTLED` |
   | `run_0013_air_dry/` | `NOT CONVERGED (stalled/plateau)` | `TRANSIENT-UNSETTLED` |

   準定常判定は `machmax,pmax` に対するもので、両 run とも保存場が0/4000/8000の **3枚しかない**ためです。壁圧比・onset・`g_exit` の判定機能や `--series` は現行ツールにありません（[check_quasisteady.py:280](/home/sano/work/forge-cond/solver_density_cuda/tools/check_quasisteady.py:280)）。保存された全 HDF5 の数値配列には NaN/Inf を認めませんでしたが、それで定常性は示せません。

   **対案:** 報告量を抽出して同ツールの判定へ渡す処理を、検証 run より先に実装する。初期過渡後に十分な枚数を保存し、onset 温度の drift は K、壁圧比は相対値で許容幅を事前定義する。`NOT CONVERGED` を準定常判定で「収束」に読み替えないことも明記してください。

8. **Major — 空気 dry は同じ \(P_0,T_0\) になっていない。密度だけの変換は誤り。**

   **根拠:** [空気入口設定:7](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/run_0013_air_dry/bcondConfig.yaml:7) は密度を6.161に変更していますが、速度はN2と同じ329.94 m/sです。設定値から CPG 関係で逆算すると、

   | 設定 | \(M_{\mathrm{in}}\) | \(T_0\) | \(P_0\) |
   |---|---:|---:|---:|
   | N2参照 | 1.05004 | 290.010 K | 844.128 kPa |
   | 空気dry | 1.06768 | 291.798 K | 862.362 kPa |

   これは流れの定常性と独立した、入力境界条件の不一致です。

   **対案:** 同じ \(P_0,T_0,M_{\mathrm{in}}\) から入口全量を再生成する。N2場を初期推定に変換するなら、密度に加えて速度を \(\sqrt{R_{\mathrm{air}}/R_{N_2}}\) 倍し、保存量を再構成する必要があります。対応する入口速度は約 **324.48 m/s**です。修正した E2 を先に検証してから E1 を開始してください。

9. **Major — cellだけの検証と広い自動適用が釣り合っていない。**

   **根拠:** §6はcellのArthurのみですが、変更対象のEOS・流束・境界・陰解法はnode/cell共有です。[verification/README.md:39](/home/sano/work/forge-cond/procedures/verification/README.md:39) は両離散化の検証を必須としています。また、単体 sweep の30–120 K、\(g\le0.3\) は、現在のEOSが許す \(g\le0.99\) や参照dryの最低温27.33 Kを覆いません。

   **対案:** node/cell双方のdry・N2・空気回帰、EOS反転往復、音速有限差分、エネルギー流束、境界状態の試験を追加する。周期・軸対称・他流束まで検証しない場合は、その構成への自動適用を制限する。既存メッシュは今回 **`VERDICT: PASS`、最大AR=21.9、最大skewness=0.381**でしたが、新しいnodeメッシュも投入前に確認してください。

10. **Minor — 単体試験の期待値と空気の組成定義が揃っていない。**

    **根拠:** [plan:97](/home/sano/work/forge-cond/plans/active/condensation-air.md:97) は1 atmの露点を78.8 Kとしていますが、指定されたJacobsen＋Antoine＋理想溶液式を再計算すると **82.236 K**です。NISTのAntoine係数自体は記載どおりでした。[NIST Oxygen](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=4)

    また、N2/O2=0.79/0.21の二成分混合なら \(M=28.8503\) g/mol、\(R=288.193\) J/kg/Kで、計画の28.9647、287.05とは異なります。

    **対案:** 二成分空気か実在乾燥空気の擬似物性かを統一し、定数を同じ組成から生成する。§1・§2・§4.3・§5と `methods/condensation.md` の期待値・既定モデルも同期してください。

推奨は、**N2 の熱力学・境界・Jacobian を整合させたうえで、CPG の「N2選択凝縮＋O2キャリア」だけを初版として実装すること**です。露点線とN2分圧線を一つの適合試験で選ぶ方針は取り下げ、修正した実験基準とnode/cell検証が通ってから既定化してください。本レビューの提案は **plan未反映**です。

指摘数: Critical 1 / Major 8 / Minor 1
