# codex レビュー: tooling-nozzle-isothermal-wall-chain (plan)

- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **stage**: `plan`
- **date**: 2026-10-08
- **commit**: `656bd37f` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m2
- **focus**: 2026-10-08 に追加した §4.7 (CONTUR の熱力学・輸送の整合 = 第 1 層、closure_version contur_v2) と §6 V-c45 (冷却壁の NS の対の事前登録: 投入の前提・ゲート・熱閉包の判定 e_A/e_B と 2u の基準) と §5.1 #13 (mesh2d の表の格子) に集中。判定基準が結果を見る前の登録として十分か、e の定義・u の取り方・k_f の扱い、格子の設計 (y1+ の見積もり方・AR)、§4.7 の 4 項目の数式と CPG での一致の主張、版の切り替えの方針を点検してほしい
- **extra**: `case/45.isobutane_m6_d155/cold_pair_mesh.py`, `case/45.isobutane_m6_d155/cold_pair.py`, `design/forge_design/meshing/mesh2d.py`, `design/forge_design/feedback/deltastar_integral.py`, `notes/reviews/2026-10-08-contur-property-temperature-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
§4.7 の整合化は進めてよいですが、現行の格子と V-c45 の判定条件は修正が必要です。  
特に格子の実形状、壁解像の面積重み、不確かさ、`contur_v2` 全体の検証が不足しています。

目的・スコープは妥当です。親計画の等温 NS 残件を引き継いでおり、既存の `wallProfile`・CHT とは責務を分けられています。§4.7 のエンタルピー分布は壁・外縁の端点条件を満たし、加速項の変換も等エントロピー・組成一定という前提で成立します。上流の面積–速度関係も整合しています。[NASA の導出](https://www.grc.nasa.gov/www/k-12/airplane/nozzled.html)

以下は重大度内でも優先順です。

1. **Major — 登録した格子は、AR ≤ 5000 を満たしていません。**

   根拠: [cold_pair_mesh.py:57](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair_mesh.py:57) の設計は軸方向間隔を使い、同:70 の見積もりも `Δx/(fr·r)` です。傾斜壁の接線長が入りません。

   保存済みの生産壁 `_band_ab/prod_confirm/prep/wall_repr.json` と `mesh_cold_tables.json` から、現行 `generate_axisym_mesh` で座標をメモリ上に再生成し、品質ツールの計算関数で確認しました。

   | 座標 | AR 最大 | AR > 5000 のセル数 |
   |---|---:|---:|
   | float64 | 5488.48 | 230 |
   | float32 に丸めた座標 | 6013.67 | 703 |

   さらに、倍精度で AR 最大のセルは skew ≈ 0.437 です。[mesh2d.py:220](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:220) は半径方向に格子を伸ばすため、収縮部は壁法線に沿いません。「壁法線・スキューなし」に限る AR 緩和を全域へ適用できません。

   **対案:** 傾斜域は既定の AR 制限を守るか、壁法線に沿う格子へ変更し、実座標の辺長・角度で設計を反復してください。また [ns_n012.py:217](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012.py:217) は部分文字列で `SOFT-PASS` も受理します。V-c45 が要求する厳密な `PASS` と照合してください。上表は変換済み HDF5 に対する正式 VERDICT ではありません。

2. **Major — 第一層が float32 の数 ULP になり、設計した層厚を保証できません。**

   根拠: 同じ再生成格子の最小半径方向間隔は **2.67×10⁻⁸ m**、局所座標の **3.58 ULP** です。float32 化による第一層間隔の相対差は最大 **24.7%** でした。幾何型は [flowFormat.hpp:7](/home/sano/work/forge-integ-1005/solver_density_cuda/flowFormat.hpp:7) で `float` です。

   また [cold_pair_mesh.py:44](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair_mesh.py:44) は、壁法線距離 `y1` を半径方向の `wall_first_frac` として扱っています。傾斜域では両者が異なり、今回の扱いは余分な細分化にもなります。

   **対案:** 法線距離と半径方向間隔を分け、変換後・実効精度での第一内部点距離を投入条件にしてください。座標重複、層厚の丸め差、局所 y₁⁺ 見積もりを検査し、1% の閉包比較に対する格子感度へ含める必要があります。24.7% は幾何間隔の差であり、そのまま流束誤差を意味するものではありません。

3. **Major — 「y₁⁺ 超過面積 ≤ 5%」を、現行ツールでは判定できません。**

   根拠: [check_wall_resolution.py:377](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_wall_resolution.py:377) の評価率・超過率は、どちらも**点数割合**です。

   保存された `y1p_ad_wall.csv` を再集計すると、y₁⁺ > 1 は点数で **3.60%**、2πr ds による面積重みでは約 **0.20%** でした。今回追加する `x_density_table` では、この差を無視できません。

   **対案:** 境界面の実面積で超過率と評価可能率を計算してください。軸対称の周長重みを二重計上しないことも確認が必要です。全壁の5%だけでは試験部の局所不足を隠せるので、収縮部・スロート・試験部ごとの値と超過位置も事前登録してください。

4. **Major — 比の STEADY だけでは、比較する両方の熱場が発達したと判断できません。**

   根拠: [plan:214](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:214) は残差の上昇・発散を除外し、δ_E の比を3地点だけで判定します。

   現行 `check_quasisteady.classify` に5点の合成系列を与えると、断熱・冷却の δ_E がともに約4%増加する場合でも、比例していれば次の判定になります。

   ```text
   断熱 δ_E: DRIFTING
   冷却 δ_E: DRIFTING
   比 R_NS: STEADY
   ```

   **対案:** 両 run の δ_E 自体、断熱壁温、冷却側 Q_w の定常性をゲートに追加してください。最大誤差を区間全体で取るなら、定常性も同じ評価点集合で確認します。5枚全部を使う `--tail 1 --min-snaps 5 --drift 0.001 --osc 0.001` 等の実効条件、延長時の窓移動規則も固定すべきです。plateau を許容する場合も、`NOT CONVERGED` 表示を残し、残差床の定量上限と判定不能の拒否を明記してください。

5. **Major — e と u の正規化が異なり、時間幅だけでは「2u」を判定の根拠にできません。**

   根拠: [plan:218](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:218) の e は R_NS に対する**相対誤差**、u は R_NS の**絶対幅**です。例えば R_NS ≈ 0.7、幅 0.005 なら、相対幅は約0.00714です。現在の `2u=0.010` と尺度が合いません。

   また、5枚の時間幅には格子差・帯 E の抽出感度・CONTUR の積分誤差が入りません。同一格子の対でも、それらが壁温間で相殺する保証はありません。

   **対案:** R_NS の代表値を「各腕の5枚平均の比」などに固定し、共通 x 点・補間・平滑化・分母の有効条件を登録してください。不確かさは同じ相対尺度に統一し、時間・格子・抽出・積分の寄与を評価します。保守的には R_NS の許容区間から各 e の上下限を作り、区間が分離した場合だけ支持とします。**1%精度の主張にも誤差上限を使う**べきで、識別できなければ判定保留です。

6. **Major — V-c45 の B 腕は、`contur_v2` 全体の検証になっていません。**

   根拠: [plan:155](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:155) の `contur_v2` は熱閉包・輸送・加速項をまとめて変更します。一方、[hform_patches:403](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/delta_contur_compare.py:403) が変更するのは熱閉包だけです。そこで B が勝っても、完全な新版やその較正係数の採用根拠にはなりません。

   加えて [同:452](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/delta_contur_compare.py:452) は `cf_scale` を渡さず、現在の `hform` は常に k_f=1 です。計画の「生産 k_f でも比較」は未実装です。

   **対案:** 熱閉包だけの A/B と、全変更を含む v2 の評価を分離してください。前者では旧 k_f を**診断用の固定係数**として両腕に使うことは妥当ですが、新版の較正済み係数とは扱わないこと。k_f=1 と生産値の両方で同じ支持が成立する条件を登録します。その後に v2 全体の物性照合・加速項検証・再較正を行い、既定変更の条件にしてください。

7. **Minor — CPG の「式として一致」と、旧実装の再現を分ける必要があります。**

   根拠: 熱閉包は同じ一定 c_p なら代数的に一致します。しかし [deltastar.py:18](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar.py:18) の Sutherland 定数は 273.15 / 110.4、NS 側は [gasProperties_d.cu:98](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/gasProperties_d.cu:98) の 273 / 111 です。[methods/design/overview.md:1243](/home/sano/work/forge-integ-1005/methods/design/overview.md:1243) の無条件の一致表現は強すぎます。

   また、新しい上流式は M=1 で 0/0 になり、現行コードの `M < 0.98` に相当する接続処理が必要です。

   **対案:** v1 のビット再現、同一物性・同一微分条件での CPG 極限、新しい NS 物性への一致を別試験にしてください。音速点の有限極限と枝の接続も明文化すべきです。

8. **Minor — 既存 A/B の数値精度確認は B 腕を検査していません。**

   根拠: [delta_contur_compare.py:467](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/delta_contur_compare.py:467) の厳しい rtol 比較は A 腕だけです。B のエンタルピー表・反転の分解能は未評価で、[同:407](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/delta_contur_compare.py:407) の 100–1800 K 外は `np.interp` が端値に固定します。

   **対案:** 両腕について積分・求積・熱物性表の細分化感度を確認し、比較量への誤差を0.1%未満に抑えてください。新版では温度・エンタルピーの範囲外を明示的に拒否します。

**推奨は、v1 を固定保存したまま、上記1〜6を先に計画と評価コードへ反映し、その後に固定形状の冷却対で検証することです。** 熱閉包 A/B の支持と、v2 全体の生産採用を別の判定にしてください。

対象の `run_0181_ns_coldmesh_ad`・`run_0182_ns_coldmesh_tw300` と元 run の HDF5 は、この作業ツリーにはありません。したがって実 run の収束・準定常 VERDICT は再確認していません。ファイル変更・forge 起動はなく、**plan 未反映**です。反映先は対象 plan の §4.7・§5.1 #10〜#14・§6 V-c45 です。

指摘数: Critical 0 / Major 6 / Minor 2
