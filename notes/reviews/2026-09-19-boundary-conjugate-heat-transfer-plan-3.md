# codex レビュー: boundary-conjugate-heat-transfer (plan)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `66b6ff36` (feature/sern-design)
- **codex**: effort `high`, 8.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m1
- **focus**: これは 1 巡目 NO-GO (notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan.md, C1/M9/m2) を全件採用して全面改訂した 2 巡目である。なお 2 巡目の初回実行は利用上限で中断し、そのとき codex 自身が検算で見つけた 3 件 (界面熱量の符号・組合せ sum(F)-C、成分セカント D_f の過小評価とスペクトル半径 1.818、放射 1500 K = 287 kW/m2) は既に反映済み (notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-2.md)。まず 1 巡目の C1/M1-M9/m1-m2 が実際に解消されたかを 1 件ずつ検証してほしい (特に C1 実効界面熱量の定義と G-cons 収支ゲート、M2 の固定点保存更新式と D_f の上界既定、M3 の R_tot=t/k_s+R_back と断熱孤立系の零空間、M4 の wall_isothermal+conjugate 属性で 5 経路に乗るか、M7 の V1)。次に §4.9 の一次適用先 (NASA CR-168015 の Mark II 超音速出口 case 5411 / C3X 亜音速) と §4.4c の fem2d バックエンドが検証計画として現実的か、記述が一次資料に照らして妥当かを見てほしい。最後に実装に進んでよいかを判定すること。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
固定点保存の更新式、界面熱量の符号、固体抵抗、等温壁の配管は改善されています。  
しかし、`D_f` の上界という主張は反例で崩れ、必須検証と依存診断の対応範囲も矛盾しています。現稿での実装開始は推奨しません。

ファイルは変更していません。コード照合、CPU 検算、保存済み run の収束再判定を実施しました。新規 CFD 計算は実行していません。

前回原文の番号は **C1、M2–M10、m11–m12** です。以下ではその番号を使います。「解消」は計画上の修正を指し、未実装機能の動作保証ではありません。

| 前回指摘 | 再評価 |
|---|---|
| C1 実効界面熱量 | **式と符号は解消**。`Q_f=ΣF−C` は正しい。G-cons の必須化も適切。ただし診断の対応範囲、接合・転送、判定閾値は未完結。今回 #1・#3・#5。 |
| M2 更新式 | **部分解消**。固定点は保存する。しかし `k_eff A/d1` の上界保証は成立しない。今回 #2。 |
| M3 固体抵抗・零空間 | **主要問題は解消**。`R_tot=t/k_s+R_back` は正しい。断熱孤立系の説明は、指定荷重で解く固体単独問題として妥当。流体との連成が定数モードを拘束する場合まで「解なし」と一般化しないこと。 |
| M4 等温壁の配管 | **解消**。`wall_isothermal` を維持すれば、壁フラグ・温度ピン・粘性壁・壁距離・DPLUR エネルギー行切離しに入る。 |
| M5 プロファイル座標 | **解消**。node 座標での評価と `VALUE/T`・EOS の検証は適切。 |
| M6 幾何・接合 | **部分解消**。primal facet と共有温度 DOF の分離は正しいが、実効節点荷重から固体への転送契約が不足。今回 #3。 |
| M7 V1 | **解析解は解消**。CPU 検算で `Tw=800 K`、両側 `q=100 kW/m²`。ただし dual-time と診断依存が衝突。今回 #1。 |
| M8 SU2 | **計画上は解消**。1D multizone 先行とシェル近似誤差の分離は適切。ローカル SU2 の引数逆転も残っており、`DIRECT_*` 選択は合理的。 |
| M9 旧 run | **解消**。予備資料への格下げは妥当。今回も両 run は `NOT CONVERGED`、壁ダンプの `qwall` は全値ゼロ。 |
| M10 界面ゲート | **部分解消**。検査対象は追加されたが、合否を再現できる数式・数値閾値がない。今回 #5。 |
| m11 費用・性能主張 | **主要問題は解消**。A→B の選択と侵襲性を理由とする判断を支持。性能は測定で評価すべき。 |
| m12 docs 先行 | **順序は解消**。S0 の実施を実装開始条件として維持すること。 |

M4 は [mesh.cpp:939](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:939)、[nodeWallDirichlet_d.cu:130](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:130)、[viscousFlux_d.cu:954](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:954)、[calcWallDistance_kdtree.cpp:129](/home/sano/work/forge/solver_density_cuda/input/calcWallDistance_kdtree.cpp:129)、[timeIntegration_d.cu:1025](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1025) で確認しました。

M9 の再実行結果は、各 run の本段 `residual_history.csv` に対する次の判定です。段階起動全体を連結した判定ではありません。

| run | 今回の VERDICT・実測 |
|---|---|
| `case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/` | `NOT CONVERGED (stalled/plateau)`。壁9枚、`qwall` 9,009値すべてゼロ。 |
| `case/44.vitiated_air_wt/run_0114_va_ns_fine_iso300_samewall_cont/` | `NOT CONVERGED (stalled/plateau)`。`rms_roUy`・`rms_roOmega` は `RISING`。壁9枚、36,009値すべてゼロ。 |

**今回の指摘一覧**

1. **Major — 必須検証を、依存先の診断が受理しません。**

   **根拠:** [依存診断 plan:32](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:32) は **dual-time・周期・軸対称を初版で明示拒否**しています。一方、本計画は [V1:299](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:299) に dual-time、V5 に周期翼列を要求し、全 run に G-cons を課しています。依存計画を現在の仕様どおり完成させても、V1 と V5 を合格判定できません。

   また、dual-time の実装は空間残差に BDF 項を追加します（[main.cpp:1804](/home/sano/work/forge/solver_density_cuda/main.cpp:1804)）。過渡での拘束反力は一般に
   \[
   C_i=D_t(V_iE_i)-R_i^{raw}
   \]
   であり、定常式 `C=−Rraw` をそのまま物理的な瞬時入熱にはできません。

   **対案:** 依存計画と共同で、**V1 前に dual-time、V5 前に周期の診断を解除するマイルストーン**を追加してください。空間残差・蓄積・拘束反力の採取位相と単位を固定し、軸対称・壁モデルは対応する解除試験まで明示的に対象外とします。「依存先の完了」だけでは不十分です。

2. **Major — `D_f=k_eff A/d1` は一般の上界ではなく、既定反復が発散する反例があります。**

   **根拠:** [plan:91](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:91) は、この係数を上界と呼んでいます。しかし実効界面熱量には、壁 CV の接線方向輸送と内部場の応答も入ります。実装も内部面ごとの伝導を組み立てています（[viscousFlux_d.cu:228](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:228)）。

   **CPU で検算した反例:** `x={0,2}`、`y={0,4,8}` の直交 median-dual 格子、単位スパン、`k=1`、側面断熱、上端温度固定、下端2節点を界面とします。AR は2です。内部温度を消去した応答は
   \[
   H=-\frac{\partial Q_f}{\partial T_w}
   =\begin{pmatrix}
   1.180556&-1.055556\\
   -1.055556&1.180556
   \end{pmatrix}.
   \]
   計画の係数は `D_f=0.25 I`。`A_s=0.1 I` とすると、
   \[
   \operatorname{eig}\!\left[(A_s+D_f)^{-1}(D_f-H)\right]
   =\{0.357143,\,-5.674603\}.
   \]
   **スペクトル半径は5.675で発散**します。これは伝導離散の反例であり、forge の実行結果ではありません。

   [plan:103](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:103) の「未収束流体では応答が小さい」も一般には逆です。第一内部点を凍結すると応答は `kA/d1`、再収束すると流体層全体の抵抗が効いて小さくなる場合があります。

   **対案:** 固定点保存式は維持し、`k_eff A/d1` を**初期推定値**へ格下げしてください。未緩和残差を評価した受理判定と、`D_f` 増大・流体状態を含む退避を既定経路に入れます。V1b に接線方向の交互温度モードを追加し、Aitken は基礎反復の安定性確認後に限定してください。

3. **Major — `fem2d` の採用は妥当ですが、連成契約と単体検証が未設計です。**

   **根拠:** [plan:176](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:176) は一般2D固体を追加していますが、§4.2 の未知数は表面温度だけです。FEM では内部節点・冷却孔節点が増え、界面への制限作用素、荷重転送、Robin 境界の組立てが必要です。

   また、[plan:162](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:162) の facet 荷重と、CV に一つしかない拘束反力の関係が未定義です。依存診断は共有角の反力を**接合部の別勘定**にします（[依存 plan:70](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:70)）。それをどの固体 DOF が受け取るかまで必要です。

   **対案:** 初版は流体・固体の外周節点を一致させ、固体全温度 `u` と界面抽出 `E` を使って
   \[
   (K_s+E^\mathsf TD_fE)u^{k+1}
   =b_s+E^\mathsf T[Q_f(Eu^k)+D_fEu^k]
   \]
   と定義してください。共有反力は一度だけ組み込み、積分済み荷重に面積を再乗算しない契約にします。平面2Dの単位は `W/m` です。

   単体検証には、**孔の Robin 条件を持つ円環伝導解析解、非一様界面荷重、共有角の保存試験**を追加してください。現状の直列抵抗・フィン試験だけでは新バックエンドを検証できません。

4. **Major — Phase 2 の指定フックは、block-DPLUR の実行経路にありません。**

   **根拠:** [plan:198](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:198) が指定する `main.cpp:1684` は、実際には **`advanceExplicitRK` のステージループ**です。定常陰解法は [advanceImplicitSteady:1740](/home/sano/work/forge/solver_density_cuda/main.cpp:1740) → `implicitNonlinearUpdate` を通ります。指定どおり実装すると、主対象の陰解法で連成が作動しません。

   **対案:** 行番号ではなく関数・評価位相で設計してください。定常陰解法の完了ステップ数で `K` を数え、最後に採取した**同一状態の `Q_f` と `Tw`**を対で固体へ渡し、次の残差組立て前に新しい `Ts` を適用します。dual-time は別契約とし、subiteration 中の壁温変更を混入させないでください。

5. **Major — G-cons／G-if は、まだ機械的な合否仕様になっていません。**

   **根拠:** [plan:306](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:306) の「熱量の0.5%」は、符号相殺時の分母と絶対許容値が未定義です。G-if には局所ノルム、絶対温度許容値、必要連続反復数がありません。

   また、`check_quasisteady.py` の既定値は **drift 5%、変動幅10%**（[コード:543](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:543)）。これだけの `STEADY` を、0.2～0.5%の比較精度の根拠にはできません。

   **対案:** `r=A_sT−b_s−Q_f` に対する局所最大ノルム、固体内部残差、温度更新量を別々に定義してください。収支の規格化には正味熱量だけでなく絶対値和と事前登録した絶対床を使います。準定常許容は各比較許容より十分小さく設定し、`Tw` と実効 `Q_f` の系列を明示的に検査します。

   さらに **G-if・区間ハッシュの実装を S6 から Phase 1 の合格判定前へ移す**必要があります。現状の順序では、V1–V3 を通した後に判定機能を作ることになります。

6. **Major — V5 は連成系の検証として有効ですが、連成誤差だけを分離する試験にはなっていません。**

   **根拠:** [plan:229](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:229) は、平板・すきまで流体側を検証し、翼で連成側を検証すれば誤差を分離できる構成です。しかし翼の圧力勾配・遷移・衝撃による流体側誤差は、平板検証では拘束できません。壁温が合っても、流体 HTC と内部冷却条件の誤差相殺を排除できません。

   一次資料の選択自体は支持します。NASA 表VIII・IXでは、Mark II `5411` は **run 42、`M2=1.04`**、C3X `4411` は **run 108、`M2=0.90`**。報告のデータ処理も2D伝導 FEM を使っており、`fem2d` は現実的です。[NASA CR-168015](https://ntrs.nasa.gov/api/citations/19830020105/downloads/19830020105.pdf)

   ただし同報告 p.21 では、内部 HTC は入口助走補正付き相関、計測断面の冷却剤温度は入口・出口測定値から推定しています。**HTC 感度帯だけでは入力不確かさを覆いません**。また、表VIIの温度比不確かさは±2%で、HTC の不確かさとは別です。[NASA 本文](https://ntrs.nasa.gov/api/citations/19830020105/downloads/19830020105.pdf?attachment=true)

   **対案:** 選定を維持し、同じ翼・同じ試験条件で、**実測壁温を与えた流体計算 → 公開条件による固体単独検証 → 壁温を未知とした CHT**の順にしてください。冷却孔ごとの HTC・計測断面温度・材料物性と不確かさを計算前に固定し、合うまで帯を広げる運用を禁止します。case/49・51 の誤差幅は、その用途のモデル感度から別途評価してください。

7. **Minor — 計画索引と親 plan に撤回済みの方針が残っています。**

   **根拠:** [plans/README.md:28](/home/sano/work/forge/plans/README.md:28) は依然として「コンパクト差分形を一次」「Bi 非依存の安定化」を記載しています。[親 plan:165](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:165) にも `wall_conjugate` が残っています。

   **対案:** S0 で両方を同期してください。今回訂正した契約が索引から再び旧仕様へ戻る状態を残さないことです。

NASA の数値・手順は公開検索索引に収録された本文で確認しました。PDF の直接取得・原画像照合は403等でできておらず、翼型座標全量や物性表の転記精度まで確認済みとはしていません。

**推奨は一つです。A→B と Mark II 主・C3X 先行を維持し、実装前に界面契約を完成させて再レビューしてください。**

優先順は、①診断の対応範囲、②`D_f` と受理・退避、③`fem2d` の荷重転送、④陰解法フック、⑤定量ゲート、⑥翼の同一条件での検証分解です。目的の重複実装は見当たらず、親 plan からの移管も確認できました。止める理由は目的ではなく、現在の解除条件と数値契約です。レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 1
