# codex レビュー: condensation-air (plan)

- **plan**: [`plans/active/condensation-air.md`](../../plans/active/condensation-air.md)
- **stage**: `plan`
- **date**: 2026-09-12
- **commit**: `d01a4118` (feature/condensation-air)
- **codex**: effort `high`, 6.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **focus**: 1 回目 NO-GO (C1/M8/m1) を全件採用した v2: CPG carrier 形 (N2 選択凝縮 + O2 キャリア, Y_w 定数, 全経路で同じ g) の EOS/流束/実現可能性の整合、潜熱と飽和圧 C–C 外挿の連動 (閉形式) と 2 段 A/B、slip の状態保持、CPG 二相音速を範囲外にした判断、理論線の再読み値 (1 kPa 38.5 K) と ±3 K をモデル間比較基準に限定した合否、node/cell と --series 判定。残る Critical/Major があれば根拠つきで
- **extra**: `notes/investigations/condensation-carrier-kantrowitz-air-survey.md`, `notes/reviews/2026-09-12-condensation-air-plan.md`, `case/34.arthur_n2_nozzle/daum_gyarmathy_theory_onset_Pdot20000_n2.csv`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

前回の Critical「N2 選択凝縮と EOS の矛盾」は、v2 の CPG carrier 形で解消しています。目的は妥当で、accepted plan に同等の実装はありません。  
ただし、物性接続・温度反転・流束・検証ツールに Major が残り、このまま実装へ進むことは推奨しません。ファイルは変更していません。

1. **Major — 提示式は、要求している「値・微分連続」を満たさない。**

   **根拠:** [plan §4.2:86](/home/sano/work/forge-cond/plans/active/condensation-air.md:86) の式を、[既存物性の係数:69](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:69)で再計算しました。

   | 接続条件 | 低温側 | 高温側 |
   |---|---:|---:|
   | \(L'(70\,K)\) [J/kg/K] | −961.20 | −1071.64 |
   | \((\ln p_{\rm sat})'(50\,K)\) [1/K] | 0.308073 | 0.278926 |

   値は連続ですが、微分は不連続です。[§5.7:113](/home/sano/work/forge-cond/plans/active/condensation-air.md:113) の微分連続試験は、この設計では通りません。

   また、提示された**線形潜熱を積分する閉形式**では、38 K の飽和圧は旧値の **0.5176 倍**です。本文の **0.59 倍**は新潜熱を一定値として外挿した場合に相当し、v2 の式の数値ではありません。

   **対案:** 初版では提示式を **C0 接続**と明記し、微分連続という合格条件を撤回してください。片側微分・正の熱容量・飽和圧の単調性を検査し、厳密な C–C 整合を主張する範囲を \(T<50\,K\) に限定する。数値例も閉形式で更新するのが最小の修正です。

2. **Major — `cond_T_from_e_cpg` の流用だけでは、指定した掃引範囲を反転できない。**

   **根拠:** [既存 Newton:75](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:75) は、30 回で収束しなくても温度を返します。新潜熱式と、この反復・初期値・物性上限クランプを Python で忠実に再現すると、次の反例が得られます。

   | 入力 | 値 |
   |---|---:|
   | \(g\) | 0.75 |
   | 正しい \(T\) | 122 K |
   | 入力 \(e\) | 61,173.786 J/kg |
   | 30 回後に返す \(T\) | **99.226 K** |
   | 返した温度から再計算した \(e\) | **−28,444.678 J/kg** |

   これは計画の \(g\le0.99Y_w=0.759429\)、25–125 K の範囲内です。反復が臨界側の物性クランプをまたいで往復します。さらに [dependentVariables:226](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/dependentVariables_d.cu:226) は反転後のエネルギーで `roe` を再設定するため、反転失敗を保存量の変更に変えてしまいます。

   **対案:** 単調な \(e(T,g)\) に対する**括弧付き Newton＋二分法退避**へ変更し、エネルギー残差と成功状態を返してください。失敗した温度で保存量を上書きしないことも設計に含める。単体試験は悪い初期推定を含め、double の反転精度と float32 保存量からの復元精度を別々に定義してください。

3. **Major — SLAU の補正式は正しいが、既存コードでは「同じ面状態の温度」で評価されない。**

   **根拠:** [SLAU:340](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:340) は単相項を再構成した `P_L/ro_L` で計算し、補正項には `T_cell` を使います。係数だけを carrier 形へ変える [plan §5.5:111](/home/sano/work/forge-cond/plans/active/condensation-air.md:111) では、この混在が残ります。

   例えば \(g=0.1\)、セル温度40 K、再構成した \(p_f/\rho_f\) が示す二相温度45 Kでは、提案係数を既存処理へ入れた静エンタルピーは **21,051.87 J/kg**、同じ面状態の EOS は **22,051.87 J/kg**。差は **−1,000 J/kg**です。単一温度を使った補正式の代数試験では検出できません。

   **対案:** \(g_f\) を一次のセル値として保持する場合も、
   \[
   T_f=\frac{p_f}{\rho_f(R_{\rm air}-g_fR_w)},\qquad
   h_f=c_{p,\rm air}T_f-g_fL(T_f)
   \]
   として**面状態を一貫して構成**してください。実際の面再構成を通す試験と、node/cell の質量・全エネルギー流束収支、`h0` の検査を §6 に加えるべきです。

4. **Major — 新 carrier モードを受け付ける構成と、実装・検証範囲が一致していない。**

   **根拠:** [plan §5.2:108](/home/sano/work/forge-cond/plans/active/condensation-air.md:108) の拒否条件には流束の制限がありません。しかし [Roe:169](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/convection/convectiveFlux_roe_d.inc.cuh:169) は既存の pure 補正 \(g(c_pT-L)\) を使い、[KEEP:97](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/convection/convectiveFlux_keep_d.inc.cuh:97) は CPG 内部エネルギーを \(p/[\rho(\gamma-1)]\) から作ります。SLAU だけの変更では対応しません。

   隣接 plan との接点も残っています。[source:101](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cu:101) の Feder 補正は、CPG では `carrierSum=0`。新しい空気 carrier で `condKantrowitz: 2/3` を指定しても、O2 の衝突冷却は入りません。

   **対案:** 初版の新 carrier モードを **SLAU・CPG N2・単一凝縮種・検証する境界構成**に限定し、未対応の組合せは拒否してください。Feder 2/3 も今回は拒否し、O2 寄与の実装を後続へ明記する。加えて、`Yw` の有限性・範囲と \(R_{\rm air}-Y_wR_w>0\) を検証する必要があります。

   CPG 二相音速の延期自体は、この限定の下では妥当です。ただし既存 block-DPLUR は**近似 Jacobian**として扱い、熱力学的に整合した frozen Jacobian と呼ばないでください。

5. **Major — `--series` が、計画で合否を判定する量を正しく評価していない。**

   **根拠:** 作業ツリーの [onset_analysis.py:45](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/onset_analysis.py:45) は「壁圧比」を中心線の `P[C]/pd` から計算しています。[出口処理:48](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/onset_analysis.py:48) は出口近傍の `Ux/sonic` であり、境界面法線による評価ではありません。

   また、[判定対象:57](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/onset_analysis.py:57) には出口 Mach 条件・`p_on`・`Pdot_on` が含まれません。`g_exit=0` から dry と推定すると、凝縮 ON でも出口まで未発達な run を dry 扱いできます。

   再実行した既存 run の判定は以下です。保存済み `CONVERGENCE_VERDICT.txt` とも整合しています。

   | run（リポジトリルートからの相対パス） | `check_convergence.py` | `check_quasisteady.py`¹ | `--series` |
   |---|---|---|---|
   | `case/34.arthur_n2_nozzle/run_0008_ref_n2/` | `NOT CONVERGED (stalled/plateau)` | `TRANSIENT-UNSETTLED` | `TRANSIENT-UNSETTLED` |
   | `case/34.arthur_n2_nozzle/run_0013_air_dry/` | `NOT CONVERGED (stalled/plateau)` | `TRANSIENT-UNSETTLED` | `TRANSIENT-UNSETTLED` |

   ¹ `machmax,pmax`。これらの判定から、onset や壁圧の定常性は主張できません。

   **対案:** 壁・中心線・出口境界を別々に抽出し、dry/cond は config から決める。報告する量を `check_quasisteady.py` の判定へ接続し、絶対・相対閾値を plan に固定してください。出口 \(u_n/c>1\) は全保存時刻の独立した失敗条件にする。1000 step 間隔の保存だけでは「実行全期間」の保証にはならないため、その表現も限定すべきです。

6. **Minor — `slip` の ghost と `bvar` では、保持すべき全エネルギーが異なる。**

   **根拠:** [plan §4.3:96](/home/sano/work/forge-cond/plans/active/condensation-air.md:96) は `roeb` も内部値としています。しかし [slip_d:70](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/boundaryCond_d.cu:70) では、ghost は速度反射、`bvar` は法線速度を除去した接線速度です。後者の運動エネルギーは内部より \(U_n^2/2\) 小さくなります。

   **対案:** ghost の `roe` コピーは維持し、境界値は
   \[
   \rho E_b=\rho E_i-\tfrac12\rho U_n^2
   \]
   としてください。理想的な slip の質量流束はゼロなので主要流れへの影響は限定的ですが、状態整合の修正としては必要です。CPG/TP の非零法線速度を使った境界単体試験を加えてください。

7. **Minor — R2 と R4 を再現する設定経路が、実装ステップにない。**

   **根拠:** [plan §4.2:91](/home/sano/work/forge-cond/plans/active/condensation-air.md:91) は一つの `latentLowT` で潜熱・飽和圧を同時に切り替えます。この二値だけでは R2「潜熱だけ新」を選べません。`c_l=1500/2500` を選ぶ経路も §5 の追加キーにありません。

   **対案:** 診断専用の飽和圧切替と液比熱パラメータを定義し、各 run に解決済み設定を記録してください。物性設定は EOS・流束・成長・蒸発・実現可能性処理へ共通に渡す。`air_kgas` も追加だけで終わらせず、[成長・蒸発が呼ぶ箇所](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cuh:232)への伝播を実装表へ書く必要があります。

**推奨は、CPG carrier＋SLAU の限定実装を維持し、上の 1 → 2 → 3・4 → 5 の順で plan を修正してから着手することです。** 6・7 も同時に実装仕様へ反映してください。node/cell の追加、空気入口速度の修正、Arthur 実験記号と Lin 計算曲線の分離、±3 K をモデル間比較に限定した判断は適切です。

理論線の新 CSV と本文の内部整合は改善しています。ただし、今回 [Grossir 原図のリンク](https://dipot.ulb.ac.be/dspace/bitstream/2013/208925/3/AIAA_2014_1153_Grossir.pdf)から図を取得できず、**1 kPa＝38.5 K の再読み値そのものは独立確認できていません**。本レビューは `HEAD d01a4118` と確認時の作業ツリーを対象とし、提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
