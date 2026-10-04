# codex レビュー: boundary-weak-isothermal-wall (plan)

- **plan**: [`plans/active/boundary-weak-isothermal-wall.md`](../../plans/active/boundary-weak-isothermal-wall.md)
- **stage**: `plan`
- **date**: 2026-09-20
- **commit**: `6178dd7a` (feature/sern-design)
- **codex**: effort `high`, 4.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m0

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

熱的壁閉包の A/B は実施価値があります。ただし、現稿にはヤコビアンの定義、距離、熱流束診断の不整合があり、そのまま実装へ進めません。既存 run は「他仮説を棄却済み」「どちらの結果でも原因が確定する」という主張を支持していません。

1. **Major — §4.3 は、§4.1 の残差の厳密ヤコビアンではない。**

   **根拠:** 物性を凍結して \(g=k_{\rm eff}A/d_1\) と置くと、提案する壁寄与は
   \[
   R_W^{wall}=-g(T_I-T_w).
   \]
   したがって、この寄与の \((\rho e)_W\) 微分は **0** で、温度を介する微分は **内部点 \(I\) の列**にあります。§4.3 の非零対角項とは一致しません。

   [SU2 の実装](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp:715)も、内部点温度による残差に対して壁点の対角項を加えています。これは移植候補となる**近似線形化**ですが、厳密微分として説明するのは誤りです。また forge は `res_roe` を右辺に取り、行列へは残差微分の符号を反転して組みます。[timeIntegration_d.cu:819](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:819)  
   計画の負符号をそのまま `diag_block` に加えてはいけません。

   **対案:** §4.3 を「SU2 型の近似対角項」と明記し、残差・厳密微分・実際の近似行列を分けて定義してください。CPG ではエネルギー対角への追加は \(+g/(\rho c_v)\)。壁点と内部点を別々に摂動する差分検査と、小規模伝導問題で符号を確認します。

   初版は **`thermalMethod: 0` に限定**するのが妥当です。TP では `e(T)` を使っており、密度微分の CPG 式をそのまま一般化できません。[thermo_d.cuh:919](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:919)

2. **Major — 第一内部点と距離の仕様が不足し、記載どおりでは SU2 と異なる。**

   **根拠:** 計画の \(d_1=|\boldsymbol{x}_I-\boldsymbol{x}_W|\) に対し、SU2 は `GetNormal_Neighbor()` と **法線投影距離** `NormalDistance()` を使います。[CNSSolver.cpp:679](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp:679)  
   forge の既存 `firstInterior()` も `d1` に法線投影距離を保存しています。[conjugateWall.cpp:118](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:118)

   点間距離を使うと、法線から30°傾いた辺では、同じ温度差に対する熱流束が投影距離版より約13.4%小さくなります。C3X の直交した第一層だけでは、この実装差を検出できません。

   **対案:** 壁ノード法線、内部点の選択、投影距離、半割面積、同一ノードへの加算方法を先に確定してください。既存探索は隣接点から最大 alignment を選ぶため、**別の壁点を選ばない検査、内部点なし・退化距離の拒否**も必要です。

   幾何と接続は初期化時に確定し、残差と陰解法が同じ値を参照する構成にします。初版は平面・静止壁に限定し、軸対称と壁ノードが周期同一視される構成は拒否するのが妥当です。翼列の遠方周期境界まで禁止する必要はありません。

3. **Major — `Qw_Wall` の流用説明が「内部双対面は不変」と衝突する。**

   **根拠:** 現行 `Qw_Wall` は、片端が壁の **W–I 内部双対面**の伝導を置換する機構です。[viscousFlux_d.cu:261](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:261)  
   壁半割面の置換は、別の `qwall_b`／`wallTreatment` 経路です。[viscousFlux_d.cu:530](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:530)

   さらに `Qw_Wall > -0.5` が有効判定なので、符号付き熱流束とマーカを兼用する構造です。冷却壁の負の流束へ無条件に流用できません。

   **対案:** 新経路は `viscousFlux_wall_d` の伝導項だけを置換し、内部カーネルには新しい `Qw_Wall` ポインタを渡さない設計にしてください。選択フラグと符号付き流束は分離し、加熱・冷却の両方向を検証します。`nodeWallDirichlet==1` も注意書きだけでなく起動時に必須化してください。

4. **Major — `ifaceRraw=0` だけでは、CHT と A/B の評価量が整合しない。**

   **根拠:** 現行 `iface_q_compact` は
   ```cpp
   keff * (T[j] - T[ic]) / fi.d1[ib]
   ```
   であり、指定壁温ではなく**壁ノード保存温度**を使います。`iface_q_2nd` も同様です。[conjugateWall.cpp:219](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:219)

   強制側では `T[ic]=Tw` なので問題が隠れますが、弱形式側では新しい境界流束と異なる量になります。外部 CHT ループの既定入力も `q_compact` です。[cht_loop.py:157](/home/sano/work/forge/solver_density_cuda/tools/cht_loop.py:157)

   **対案:** `Tw_bc` と `T_W` を別々に記録し、比較用のコンパクト熱流束は両枝とも指定 `Tw_bc` で定義してください。保存温度からの勾配は別名の診断量にします。CHT 検証は `--flux q_eff` を明示し、`ifaceFw` が実際に残差へ加えた新流束を保存することも確認します。

   **熱的拘束を完全に外した場合の \(C=0\) 自体は妥当**です。ただし、それだけでは流体側が定常になったことも、固体側との収支が閉じたことも保証しません。

5. **Major — 目的と前提の因果主張が、実測と前回レビューを越えている。**

   **根拠:** 今回、`check_convergence.py` を再実行しました。

   | run | VERDICT と主要結果 |
   |---|---|
   | `case/53.c3x_vane_cht/run_0009_uniformTw/` | **`NOT CONVERGED (stalled/plateau)`**。`rms_roUy` 2.9桁、`rms_roK` 2.4桁低下で停滞 |
   | `case/53.c3x_vane_cht/run_0029_roe/` | **`NOT CONVERGED (stalled/plateau)`**。`rms_roe` 最終 \(3.84\times10^4\)、低下0.0桁 |

   また `run_0009_uniformTw/wall_series.csv` に対する `check_quasisteady.py` は **`q_total: STEADY`** ですが、このCSVには交番振幅の列がありません。積分熱量の定常性を交番振幅へ転用できません。run の索引は [case/53 README](/home/sano/work/forge/case/53.c3x_vane_cht/README.md:293)です。

   前回レビューも「温度の単純な量子化では説明不足」と述べ、**残差組立ての桁落ちは未検証**としています。[レビュー記録:13](/home/sano/work/forge/notes/reviews/2026-09-20-codex-c3x-checkerboard-triage.md:13)

   **対案:** 「棄却済み」を訂正し、目的を**熱的壁閉包変更の効果測定**に限定してください。結果の扱いは、改善・不変に加えて **帰属不能** を設けます。Aが再現しない、片枝が未収束、振幅が非定常なら帰属不能です。不変でも「この変更だけでは解消しない」までしか言えません。

6. **Major — §6 はケース選択は適切だが、実装の正しさと切り分け成立を判定する条件が足りない。**

   **根拠:** [計画§6](/home/sano/work/forge/plans/active/boundary-weak-isothermal-wall.md:133)の「−24%より良い」「有意に下がる」だけでは、整合性・精度次数・有意差が定義されません。G-cons の `q_floor`・絶対許容値、CFL掃引の固定条件も未指定です。

   `run_isoT_cond*` には cell の run と `cfl_pseudo: 100` の旧入力が混在しています。実行対象をワイルドカードで指定できません。case/48 の既存手順には、段階起動と **3%以上の変化を回帰とする基準**があります。[48-flat-plate-cooled.md:9](/home/sano/work/forge/procedures/verification/48-flat-plate-cooled.md:9)

   **対案:** 実装前に次を検証表へ追加してください。

   - **局所演算:** 定数温度場の熱流束ゼロ、加熱／冷却の符号、壁寄与の一回加算、近似行列の符号。
   - **純伝導:** 使用する node run を特定し、直交・傾斜格子と格子細分化で境界整合性を確認。厳密解シードに加え摂動からの緩和も確認。
   - **C3X:** 同一ソース・同一初期場のA/B。交番振幅の抽出式、正規化、保存間隔、末尾窓、有意差閾値を固定し、その系列を `check_quasisteady.py --series-csv` へ渡す。
   - **CFL・CHT:** CFL以外を固定し、同一解への到達と安定限界を分けて判定。G-cons の単位・絶対床・許容値と、解析解誤差の基準を登録。
   - **投入条件:** メッシュ品質、適用可能なIC・段階起動、NaN検査、判定区間、新規runとREADME索引を明記。

   case/24・52・53・48という選択は対象機能に合っています。cell は起動拒否とコード経路の確認に留める方針で、現行の node 主体の検証規則と整合します。

**推奨は、対象を絞った SU2 型弱閉包の診断用 opt-in を、上記修正後に実装することです。** `accepted/` の DOF-only 化や境界スカラー拡散の修正は、このエネルギー閉包を解決しておらず、直接の重複ではありません。

実装前の優先順は **①残差・近似行列と対応範囲 → ②幾何・壁専用流束 → ③診断／CHT契約 → ④定量ゲートと因果主張の訂正**です。その後、局所検査→case/24→case/52→C3X→case/48の順で進めてください。

ファイル変更なし。提案は **plan 未反映（read-only）**です。

指摘数: Critical 0 / Major 6 / Minor 0
