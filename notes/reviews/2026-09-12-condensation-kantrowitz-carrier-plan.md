# codex レビュー: condensation-kantrowitz-carrier (plan)

- **plan**: [`plans/active/condensation-kantrowitz-carrier.md`](../../plans/active/condensation-kantrowitz-carrier.md)
- **stage**: `plan`
- **date**: 2026-09-12
- **commit**: `0512823d` (feature/condensation-air)
- **codex**: effort `high`, 7.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m3
- **focus**: §4.1 の Feder carrier 形 (q, b² の定義とキャリア量のセル値からの逆算) が Feder 1966 / Wedekind 2008 と整合するか、mode 2/3 の分け方、§4.2 σ 倍率による感度試験の妥当性と Tolman 補正を入れない判断、§6 の検証 (Wysłouzil 2D 6 run と onset 序列の合否条件) の十分性
- **extra**: `notes/investigations/condensation-carrier-kantrowitz-air-survey.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

carrier 冷却の追加は未解決の課題であり、目的は妥当です。  
ただし、キャリア物性の逆算、陰解法への組み込み、onset の合否条件は実装前に修正が必要です。表面張力の文献調査には結論を逆に読んだ箇所があります。

1. **Major — `cp_cell` からのキャリア比熱逆算は、純蒸気極限で破綻する。**

   **根拠:** [plan:64–68](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:64) は、混合物の比熱から定数 `CondSpeciesProps.cv=1393.5` を引いています。しかし、[dependentVariables_d.cu:142](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/dependentVariables_d.cu:142) の `cp_cell` は、各種の温度依存 NASA-9 比熱を合成した値です。差し引く蒸気比熱の定義が一致しません。

   `run_0335_condfix_new/species_db.yaml` の係数による再計算では、230 K の真の N2 比熱は `742.209 J/kg/K`。提案式は `Yw=0.99` で `649.076`、`Yw=0.999999` で **−940000 J/kg/K** を返します。これは倍精度でも発生する定義の不整合で、float32 の差し引き誤差がさらに加わります。`Yc=0` では明示的なゼロ除算にもなります。

   また、複数キャリアの衝突項は種別の和です。平均分子量と平均比熱だけでは、`sqrt(Mv/Mi)` を含む和を一般には復元できません。[Horsch らの一次論文、式14](https://mb.uni-paderborn.de/fileadmin-mb/tdy/Publikationen/Veroeffentlichungen/co2air.pdf) とも、この点で異なります。

   **対案:** mode 2/3 は種 DB と種組成から衝突項を直接集計してください。蒸気枯渇時の無限大を避けるには、
   \[
   a_v=\frac{Y_w-g}{M_v},\quad
   \theta=
   \frac{a_v\hat q^2}
   {a_v(\tilde c_{v,v}+\tfrac12)+
   \sum_{i\ne v}\frac{Y_i}{M_i}\sqrt{\frac{M_v}{M_i}}
   (\tilde c_{v,i}+\tfrac12)}
   \]
   と評価できます。ここで \(\hat q=q/(k_BT)\)。内部計算は既存どおり double とし、mode 0/1 は従来経路を維持します。今回の `MIXDRY` は実際に N2 ですが、一般の混合物をまとめた擬似種まで厳密に扱えるとは主張しないでください。

2. **Major — 新しい補正を陰解法の差分評価へ伝える実装ステップが欠けている。**

   **根拠:** 本体の核生成評価は [condensationSource_d.cu:169](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cu:169) ですが、`src_jac` は別途、[同:207](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cu:207) と [同:218](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cu:218) から `cond_source_vector` を呼びます。その内部にも [condensationSource_d.cuh:151](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cuh:151) の核生成呼び出しがあります。

   §5 は `cond_nucleation` の引数追加だけを明記し、この経路を扱っていません。新引数を既定値で補う実装にすると、本体は carrier 補正、差分側は純蒸気補正となり、**異なるソースモデル同士を引いた値が Jacobian に入ります**。

   **対案:** `cond_source_vector` と全呼び出し元を変更対象に明記し、本体・温度摂動・モーメント摂動で同一モデルを使ってください。温度摂動時に比熱を再評価するか凍結するかも定義し、選んだ定義に対してソースと差分係数を検証します。block-DPLUR 全体の再設計は不要ですが、この伝播は必須です。

3. **Major — 「IAPWS 外挿は241.8 Kまで実測で支持」という文献要約が誤っている。**

   **根拠:** [plan:51](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:51)、[調査ノート:49](/home/sano/work/forge-cond/notes/investigations/condensation-carrier-kantrowitz-air-survey.md:49) は Vinš 2020 を「外挿と一致、異常なし」に含めています。しかし当該論文は、**−20 ℃未満で外挿 IAPWS 式から有意な偏差を検出し、深い過冷却域で異常の余地がある**と報告しています。[Vinš et al. 2020、著者抄録](https://pubmed.ncbi.nlm.nih.gov/32419467/)

   Tolman 長の引用も不正確です。Wilhelmsen 2015 の中心的手法は **square-gradient theory＋CPA EOS** であり、単に「SPC/E・TIP4P/2005 の MD 結果」ではありません。[Wilhelmsen et al. 2015](https://diposit.ub.edu/dspace/bitstream/2445/67221/1/653251.pdf)

   **対案:** methods とノートを訂正し、平面界面の温度外挿と曲率依存を分けてください。**本段階で Tolman 補正を実装しない判断には賛成**です。ただし根拠は「異常が実測で否定されている」ではなく、対象温度・臨界核サイズで採用する曲率モデルと係数の検証が不足していることです。半径50 nmでの曲率補正が小さくても、低温の平面 σ 自体の精度は保証されません。

4. **Major — 「onset 序列が逆転したら実装誤り」は成立しない。**

   **根拠:** [plan:117](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:117) は局所の \(\theta\) の大小を、結合流れ場の onset 序列に直結させています。しかし報告する onset は [compare_condfix.py:71](/home/sano/work/forge-cond/case/16.nozzle_wys/compare_condfix.py:71) の **中心線 `g=1e-3` 到達位置**です。核生成率そのものではありません。

   実装では核生成・成長が合算され、潜熱を介して温度・圧力・過飽和度が変わり、さらに Δg・ΔT・蒸気量による制限が掛かります。[condensationSource_d.cu:181](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationSource_d.cu:181)  
   σ 変更も核生成障壁だけでなく臨界半径と成長率を変えます。この系に onset の比較原理は示されていません。

   **対案:** 序列の合否判定は、**同一の `T,pv,rhov,Y,g` を固定した局所評価**へ移してください。想定する適用域で `J0 ≥ J3 ≥ J2 ≥ J1` を検査します。結合 run の onset と σ 応答の符号は観測項目とし、逆転時は局所状態・成長・制限係数を調べる条件にします。時間変動や空間分解能以下の差を、有意な序列として扱わないことも必要です。

5. **Major — mode 1 の σ 感度だけでは、新しい carrier モデルの不確かさを評価できない。**

   **根拠:** [plan:112](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:112) は σ±3% を旧純蒸気形だけで試します。新モデルでは核生成域と成長履歴が変わるため、その感度を mode 2/3 に転用できません。また、一定倍率は温度依存・半径依存の誤差を表さず、[plan:41](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:41) の「Tolman 補正の代替」にはなりません。±3%を対象温度域の不確かさ幅とする根拠もありません。

   **対案:** **6 run の予算は維持し、σ±3%の2本を mode 3 基準へ変更**してください。結果は「一定倍率による局所感度」と呼び、信頼区間や曲率モデルの妥当性確認とは区別します。核生成・成長・有効時の蒸発 Kelvin 項に同じ倍率を適用する設計自体は妥当です。

6. **Major — node の6 runと単一点の単体試験だけでは、共有実装の検証が不足する。**

   **根拠:** [plan:103](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:103) 以下には cell 検証も、実 CUDA 経路の係数照合もありません。[verification/README.md:39](/home/sano/work/forge-cond/procedures/verification/README.md:39) は共有コードを node/cell 双方で検証するよう要求しています。

   230 K の手計算は再現できました。コード物性で計算すると `θ1=167.404、θ2=4.046、θ3=2.982` です。ただし参照 run の `res_48000.h5` では、中心線で最初に `g>1e-3` となる節点は **T=213.84 K、ln S=5.032**。230 K・ln S=3.4 の一点だけでは、実際に通る状態域を覆いません。

   **対案:** 実装前に次を検証表へ追加してください。

   - 温度・過飽和度・組成の掃引、`Yc→0`、`Yv→0`、pure N2、mode 0/1 回帰。
   - CPU 参照値と CUDA ソース・差分係数の照合。CPU の double 式の許容誤差と、float32 経路の許容誤差を分ける。
   - cell の mode 1/3 対照、凝縮 OFF 回帰、σ=1 の全経路回帰。
   - 核生成域でソース制限が結果を支配していない確認。支配する場合は CFL を下げた対照。
   - メッシュ品質、IC の同一メッシュ・同一 index 対応、NaN、収束、`check_quasisteady.py` と報告量用 `--series` の判定条件を明記。

   `cond_sigma` と構造体の実在場所は [condensationProperties_d.cuh:252](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:252) です。これも変更対象と full rebuild の対象に含めます。局所ソース変更だけを理由に周期・軸対称の大規模検証を全追加する必要はありませんが、それらまで検証済みとは扱えません。

7. **Minor — mode の位置付けと pure 極限の記述が矛盾している。**

   **根拠:** [plan:35](/home/sano/work/forge-cond/plans/active/condensation-kantrowitz-carrier.md:35) と43行は「pure で mode 1 と同値」、47行と67行は「約2%差」としています。後者が正しい説明です。mode 3 はキャリアがなくても表面仕事項が残ります。

   [Wedekind 2008、式8–10](https://arxiv.org/pdf/0804.1516) に対応するのは表面項を含む mode 3 です。mode 2 はその項を落とした比較モデルです。なお `q<0` を `fmax(q,0)` で潰す操作も、提示した \(q^2\) の式とは異なります。

   **対案:** **mode 3＝物理モデル、mode 2＝表面項を除く比較用、mode 1＝旧結果再現用**と明記してください。pure は各式の解析値を個別に検査します。`q<0` はモデルの適用域を定義して扱い、黙って等温補正に変更しないでください。

8. **Minor — 「既定値は1のまま」はコードと異なる。**

   **根拠:** [solverConfig.hpp:394](/home/sano/work/forge-cond/solver_density_cuda/input/solverConfig.hpp:394)、[solverConfig.cpp:689](/home/sano/work/forge-cond/solver_density_cuda/input/solverConfig.cpp:689) の既定値は **0** です。1 は Wysłouzil 参照 config の明示指定です。

   **対案:** 「グローバル既定値0を維持、参照 run は1を明示」と訂正してください。既定不変を検証するため、キー省略時のテストも加えます。

9. **Minor — 小半径の見積りに算術上の不整合がある。**

   **根拠:** [調査ノート:58](/home/sano/work/forge-cond/notes/investigations/condensation-carrier-kantrowitz-air-survey.md:58) の「半径1 nmで30–60分子」は、コードの液密度約994 kg/m³では **約139分子**です。また障壁 \(A=\Delta G^*/k_BT=50–70\) に σ+10%を適用すると、障壁増分は
   \[
   [(1.1)^3-1]A=16.55–23.17
   \]
   であり、「+5–6 \(k_BT\) と整合」とはいえません。

   **対案:** 半径・分子数・障壁変化を同じ条件で再計算してください。引用論文の曲率補正と、平面 σ の一定倍率変更を同じ計算として説明しないことが必要です。

実測確認には、隣接 checkout `/home/sano/work/forge` の [`case/16.nozzle_wys/run_0335_condfix_new/`](/home/sano/work/forge/case/16.nozzle_wys/run_0335_condfix_new) を使用しました。全保存時刻の `VALUE/*` は NaN/Inf 0。再実行した判定は以下です。

| 確認 | 結果 |
|---|---|
| `check_convergence.py` | `PASS (converged)`、凝縮残差を含む |
| `check_quasisteady.py --quantity pmax,machmax` | `OVERALL: ALL STEADY` |
| `compare_condfix.py --series` | `VERDICT(series): STEADY` |
| 保存済み `MESH_QUALITY.txt` | `PASS`、AR最大724.2、skewness最大0.203 |

onset は22.513 mm、末尾30000–48000 stepの変動幅は `9.07e-5 mm`。参照 run の記述は実データで裏付けられています。[run 索引](/home/sano/work/forge-cond/case/16.nozzle_wys/README.md:334)  
`--series` は描画初期化が read-only 環境で失敗したため、描画部分のみメモリ上で除外し、計量・判定処理はそのまま実行しました。ファイル変更はありません。

**推奨は、mode 3 を中心とする carrier 補正の計画へ修正して進めることです。** 実装前の優先順は上記 **1→2→3→4→5→6** とし、7–9も同時に訂正してください。目的は既存の γ 修正や accepted の凝縮計画と重複していません。Tolman 補正・既定変更・分圧スイープを後続に置く範囲設定も妥当です。

Feder 1966 本文は取得できず、直接照合は未完了です。式の確認は Wedekind 2008 と、多原子・複数キャリアを記述する Horsch らの一次論文で行いました。レビューのみの依頼のため、以上は **plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 3
