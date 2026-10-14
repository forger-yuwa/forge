# codex 諮問 (diagnose): lowmach-origin-next-test

- **brief**: [`notes/reviews/briefs/2026-10-11-lowmach-origin-next-test.md`](../../notes/reviews/briefs/2026-10-11-lowmach-origin-next-test.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `a5491d23` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.9 min, rc=0
- **結論**: **次は、同じ FP64 バイナリで「step 末尾の ρE だけの量子化 OFF／ON」を各 2,000 step 比較し、入口〜収縮部の初期応答を再現する十分条件かを調べる。**
- **extra**: `case/45.isobutane_m6_d155/c45_diffmap.py`, `case/45.isobutane_m6_d155/c45_diffmap_loc.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **次は、同じ FP64 バイナリで「step 末尾の ρE だけの量子化 OFF／ON」を各 2,000 step 比較し、入口〜収縮部の初期応答を再現する十分条件かを調べる。**

第 1 仮説: ρE の step 末尾の量子化だけで、入口〜収縮部に観測された初期の温度差の主要成分を生じさせられる。確度: **低**
  
  根拠: [plan §6.38](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:933) は、`case/45.isobutane_m6_d155/run_0488_tr_f32/` と `run_0487_tr_fp64/` の差が step 1,000 で上流側の閾値を超えると記録している。コード上も、ρE は内部エネルギー・温度反転の入力になる（[dependentVariables_d.cu:83](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/dependentVariables_d.cu:83)）。ただし、発生源を ρE に特定する証拠はない。対象 run の実データは今回のワークスペースに無く、位置・数値は独立再検証できていない。

  反証条件: 下記の陽性対照が初期応答を再現し、ρE の量子化も実効性検査を通ったうえで、介入による応答の区間が基準効果の ±10 % 内に収まること。棄却できるのは**「この時点での ρE 単独量子化で十分」**という限定仮説であり、途中の丸めを含むブリーフの H1 全体ではない。

第 2 仮説: 他の保存量、原始量・物性・残差などの途中の丸めとの組み合わせが必要。未確認。

第 3 仮説: 上流域では微小な摂動を増幅し、観測の一部は初期化差や再実行の揺れを表している。未確認。float 同士でも局所差が育つという記録と整合する。

判別 A/B:

- **対照 C**: FP64、量子化 OFF、独立に 3 本。
- **介入 E**: 同じ FP64 バイナリで、step 末尾に `roe ← double(float(roe))` だけを適用、独立に 3 本。
- **陽性対照 F**: float、独立に 3 本。短い期間・選んだ指標で既知の症状が再現するかを確認するために使う。因果を判別する比較は C 対 E の一つとする。
- 全格子・共通入力 `_tr/q32_init.h5`・実効設定を維持し、各 2,000 step 固定、step 0 と 500 ごとの出力を保存する。提案の R2 は今回は省く。

介入位置は現在のフック、すなわち等温壁のピン後・`updateVariablesOuter` 前でよい（[main.cpp:2538](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2538)）。これは**全 CV の末尾量子化の試験**であり、上流域だけが発生源だという試験ではない。

判定は次のように事前登録する。

1. 共通の double 座標と節点対応を検査し、評価域を `x/r_t < −2` の内部節点に固定する。入口・壁・軸の境界節点は別集計とする。
2. 既存の A・B・B′から、各時点の符号つき温度差を固定した参照パターン `v_s` とする。分母は共通の A の初期温度。新しい run の結果で評価域・符号・パターンを選び直さない。
3. 各 run の温度変化 `z_r,s = (T_r,s − T_r,0)/T_A,0` を用い、指標を  
   `a_r = mean_s〈v_s, z_r,s〉/〈v_s, v_s〉`  
   とする。時点は事前に `s = 1,000、1,500、2,000` と固定する。この定義では、既存の float−FP64 の平均応答が 1 になる。絶対値の分布・速度・圧力は副記録とする。
4. 独立標本は **run ごとの a_r**。節点数や保存時点を反復数に数えない。`e_E = mean(a_E) − mean(a_C)`、`e_F = mean(a_F) − mean(a_C)` に腕別の標本分散から Welch 区間を付け、二つの比較について全体 95 % に補正する。

→ **結果 A**: 陽性対照の区間下限が 0.5 以上、かつ介入の区間下限も 0.5 以上なら、「ρE の末尾量子化だけで、登録した初期温度差の成分を少なくとも半分の規模で生じさせられる」。第 2 仮説の「その成分には他の丸めも必須」を退ける。

→ **結果 B**: 陽性対照が同条件を満たし、介入の区間全体が `[−0.1, +0.1]` 内なら、第 1 仮説を退ける。

中間・逆符号・陽性対照の不成立は判別不能。3 反復で検出力が足りる保証はなく、有意差が出ないことを B と扱わない。0.5／0.1 は診断用の暫定閾値である。この判定から下流の θ_r、定常解、全場の再現までは主張しない。

前提として C・E の**初期化後**の保存量一致、量子化直前・直後の変更対象、`roeN/roeM` への引き継ぎを確認する。全残差と保存場の有限性、`check_convergence`・`check_quasisteady` の結果も記録する。有限期間の応答試験なので、収束・準定常の合格とは分ける。

やらない方がよいこと:

- R0〜R3 の四腕を一度に投入し、「R1・R2 が離れなければ H1 を棄却」と読むこと。
- 2,000 step の上流応答を、30,000 step の下流 θ_r の原因確定に使うこと。
- 現段階で修正候補に順位を付けること。**順位は保留**する。ρE 単独量子化の十分性が成立しても、double 化の必要範囲や基準値差分化の有効性は別に検証する必要がある。

呼び出し側の前提への異議:

1. **Major／要再検証 — 「入口で生まれて下流へ運ばれる」は観測を超えている。**  
   [plan:937](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:937) 自身、step 1,000 の下流に `3e−6` の速度相対差を記録している。「差がない」ではなく「指定の閾値未満」である。また、[c45_diffmap.py:38](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/c45_diffmap.py:38) は局所初期値で割るため、低速域の速度相対差は大きく見えやすい。定常の局所擬似時間を使う計算の step 上の先後だけでは、物理的な移流経路も確定しない。  
   **対案:** §6.38 は「指定の保存時点・正規化・閾値では、上流域で先に閾値を超え、後に下流域でも超えた」と記録する。原因・伝播経路は仮説として分離する。

2. **Major／要再検証 — ρE の大きさと生成エンタルピーの前提が未確認。**  
   [thermo_d.cu:91](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/thermo_d.cu:91) は、`thermoHrefTemp > 0` なら種ごとの基準エンタルピーを除去する。したがって「TP だから生成エンタルピー込み」「ρE ≈ p/(γ−1)」は、この run の量子化幅の根拠にならない。  
   **対案:** 実効設定・物性記録・局所 `roe` で確認する。単相・組成固定の局所線形化なら、  
   `δT ≈ δ(ρE)/(ρc_v)`、`δp ≈ (R_mix/c_v)δ(ρE)`。  
   実際の丸め誤差と ULP 幅を区別し、局所動圧 `ρ|U|²/2` と比較する。ブリーフの代表値でも M = 0.005〜0.007 の動圧は約 6.6〜13 Pa であり、M = 0.0005 の約 0.066 Pa と同一視できない。

3. **Major／却下 — R1・R2 の結果から「量子化一般」や H1 全体を判定できない。**  
   現行診断は平均流・SST・化学種を含む保存量を対象とする（[quantizeStateDiag.cpp:162](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/quantizeStateDiag.cpp:162)）。R1・R2 はその全体を覆わず、両者の同時作用も測らない。どちらも陰性でも、組み合わせや途中の丸めが残る。  
   **対案:** 今回は上記の ρE 単独の十分性に問いを限定する。R1 が陽性でも、float の実際の原因や必要性とは呼ばない。

4. **Major／要再検証 — 絶対差の割合・中央値だけでは再現を判定できない。**  
   [c45_diffmap_loc.py:16](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/c45_diffmap_loc.py:16) の集計は符号を捨てる。同じ割合でも逆向き・別の場所の差になり得る。さらに、共通入力から初期化後の ρE が精度別に変わった実績がある（[plan §6.31:801](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:801)）。自分の step 0 を差し引くだけでは、その後の増幅への初期化差の影響は消えない。  
   **対案:** C・E の初期化後一致を必須とし、符号と場所を保持した指標を使う。F は陽性対照に限り、初期化差も記録する。

5. **Major／却下 — 「H1 なら低マッハ前処理では直らない」とはまだ言えない。**  
   前処理は失われたビットを復元しないが、圧力差に対する流束応答や増幅を変え得る。実装も `lowMachPrecond` 1/2 で質量流束の圧力差項を変更する（[methods/convection/implementation.md:155](/home/sano/work/forge-integ-1005/methods/convection/implementation.md:155)）。  
   **対案:** 効果は未検証と書く。ただし今回は変更しない。`lineImplicit` と `lowMachPrecond >= 2` は併用不可であり、単純な一キー比較にもならない（[solver-settings.md:67](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:67)）。

6. **Minor／採用条件 — 丸め対象の選択には自己点検も追従させる。**  
   現行コードは対象と履歴のリストを別々に持つ（[quantizeStateDiag.cpp:162](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/quantizeStateDiag.cpp:162)）。  
   **対案:** `roe` 選択時は対象と履歴検査を整合させ、非選択量はフック直前・直後で不変と確認する。他の量まで float32 表現可能になることを要求しない。

不足情報: 対象の `run_0487`〜`run_0489`、`_tr/q32_init.h5`、`c45_diffmap.npz`、実効 config・物性記録・元の VERDICT がローカルに無く、初期マップの再計算と局所 EOS 感度の評価は未実施。

保存済みの [case/48 判定 JSON](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/c48_li_judge.json) では、θ の `max(|e| + E) = 9.23538e−5` を再計算でき、記録の判定は「結果 B（再現せず）」だった。収束判定は四本とも **`OVERALL: CHECK FAILURES ABOVE`**。これは限定された有限期間の再現試験の結果として採用する。[qAccumulator 判定 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/qa_judge.json) も **「判別不能」**であり、平均流の丸めを除外する材料にはしない。

**plan 未反映。** 読み取り専用のレビュー依頼に従い、ファイル変更・forge 起動は行っていない。採用時は呼び出し側で §6.38 の解釈と新試験の事前登録、§5.1 #19 の対象選択仕様へ反映すること。
