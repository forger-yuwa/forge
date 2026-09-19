# codex レビュー: boundary-conjugate-heat-transfer (plan)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `21796754` (feature/sern-design)
- **codex**: effort `high`, 7.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M9/m2
- **focus**: §4.1 の作戦選択 (A 外部弱連成 → B ソルバ内薄肉シェルの 2 段、monolithic 不採用) が forge の node median-dual 構造に対して妥当か。§4.2 の抵抗加重 (Robin) 更新の安定性論、§4.3 の界面 q_w 定義 (コンパクト差分形を一次、再構成勾配形と拘束反力を audit) が viscousFlux_d.cu / nodeWallDirichlet_d.cu の実装と整合するか。Ts が per-face bvar なので CUDA 改修なしで壁温分布を与えられるという前提 (mesh/mesh.cpp:64, cuda_forge/boundaryCond_d.cu:308, cuda_forge/nodeWallDirichlet_d.cu:75, boundaryCond.cpp の applyInletProfiles) を実機で検証してほしい。§6 の V1 (解析解) と V2 (SU2 CHT 同一メッシュ) が実行可能か。親 plan は plans/active/tooling-nozzle-isothermal-wall-chain.md §4.6-4 と plans/active/case-hypersonic-gap-heating-validation.md で、リポジトリ内にあるので必要なら直接読むこと。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
A→B の段階導入は妥当ですが、界面の保存性、固体抵抗、V1 の解析解に修正必須の問題があります。  
特に「コンパクト差分を一次、離散収支は任意 audit」の契約では、熱量が閉じない解を合格させます。

ファイルは変更していません。コード追跡と保存済み run の再評価を実施しました。非一様 `Ts` の GPU 実行確認は、`nvidia-smi` が `GPU access blocked by the operating system` を返したため未実施です。

1. **Critical — 界面熱量の定義が、forge の離散エネルギー保存と整合していません。**

   **根拠:** [viscousFlux_d.cu:529](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:529) は再構成勾配の熱流束を壁 CV に加算し、[nodeWallDirichlet_d.cu:90](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:90) がその後のエネルギー残差を消します。内部 W–I 面の伝導は、面物性と非直交補正を含む別の離散です（[viscousFlux_d.cu:223](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:223)）。

   したがって、**実効壁熱量は拘束反力だけでも、任意の片側差分だけでもありません**。依存先 plan 自身が「物理境界流束＋拘束反力」と明記しています（[tooling-energy-balance-diagnostics.md:58](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:58)）。本 plan §4.3 はこの訂正を取り込んでいません。

   また、固体向きを正とするなら、コンパクト形は `k_eff(T1−Tw)/d1` です。記載の `k_eff(Tw−T1)/d1` は逆符号です。

   **対案:** 保存的な実効界面熱量を連成の正本とし、コンパクト差分は Robin 係数の近似・精度診断に使ってください。拘束診断の完成を**連成結果の合格条件**にし、同一状態での流体・固体・接合部を含む熱収支に定量ゲートを設けるべきです。

2. **Major — Robin 更新の最適性を、異なる熱抵抗モデルへ拡張しています。**

   **根拠:** [plan:90](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:90) の D–N の解析は正しいですが、そこでの `h` と、次式の `kf/d1` は同じではありません。流体を再収束させると `T1` も `Tw` に応答します。

   例えば `h=100`、`gf=kf/d1=1000`、`gs=1/Rs=100` なら、局所差分と整合する流体応答に対して、記載の更新の誤差増幅率は

   \[
   \frac{g_f-h}{g_f+g_s}=0.81818
   \]

   です。`ω*=1/(1+Bi)=0.5` による一回収束とは違います。分布・温度依存物性・未収束流体を `K` step ごとに更新する Phase 2 への無条件安定性も導けません。

   さらに、面内伝導を解いた後に局所抵抗平均を掛ける手順では、**シェル方程式の固定点を保存する更新式が未定義**です。

   **対案:** 固体の離散作用素を `As`、荷重を `bs`、積分界面熱量を `Qf` として、例えば

   \[
   (A_s+D_f)T^{k+1}=b_s+Q_f(T^k)+D_fT^k
   \]

   のように、元の収支を固定点として保つ式を明記してください。Aitken には係数制限・分母退化処理・残差増加時の退避を設け、V1b の「全 Bi で単調収束」は、この限定モデルの試験と分布問題の試験に分けるべきです。

3. **Major — 固体モデルが厚さ方向の抵抗を落とし、CG の成立条件も満たしていません。**

   **根拠:** [plan:154](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:154) は `R_back=0` を背面等温、`1/hc` を冷却剤条件としています。しかし未知数がガス側表面温度 `Tw` なら、背面等温でも表面から背面までの抵抗 `t/ks` が必要です。記載どおりでは `Tw=Tb` となり、共役スラブが通常の等温壁に退化します。

   また、背面断熱・端部断熱の独立固体を、指定熱流束だけで先に解けば、Laplacian は定数零空間を持ちます。正味入熱が非零なら定常解自体がありません。「CG、対称正定値」は全入力に対して成立しません。

   **対案:** `R_back` をガス側表面から背面環境までの**全抵抗**として定義し、背面等温なら `t/ks`、冷却剤なら `t/ks+1/hc` としてください。温度依存物性・多層はその定義に沿って扱い、断熱孤立系は Robin 項込みで解くか、適合条件と零空間処理を明示してください。

4. **Major — `wall_conjugate` の登録だけでは、既存の等温壁経路に入りません。**

   **根拠:** 種別名による判定が各所に直書きされています。

   - 壁・等温壁フラグ: [mesh.cpp:929](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:929)
   - 温度ピン・残差射影: [nodeWallDirichlet_d.cu:130](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:130)
   - 粘性壁処理: [viscousFlux_d.cu:954](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:954)
   - 壁距離: [calcWallDistance_kdtree.cpp:129](/home/sano/work/forge/solver_density_cuda/input/calcWallDistance_kdtree.cpp:129)

   特に `iso_wall_flag` が立たなければ、[timeIntegration_d.cu:1025](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1025) の DPLUR エネルギー行切離しが効きません。§7 の変更一覧には、この配管が欠けています。

   **対案:** 初版では内部の物理種別を `wall_isothermal` に維持し、連成を別属性にしてください。これなら既存カーネルの再利用方針に合います。また、cell では `vizBfaceNodes` が空なので、シェルとの一対一対応を両離散化共通とは主張できません。初版の実装保証範囲も node に限定するのが妥当です。

5. **Major — `wallProfile` の機構は既にほぼ汎用ですが、そのまま流用すると温度を違う座標へ課します。**

   **根拠:** [boundaryCond.cpp:244](/home/sano/work/forge/solver_density_cuda/boundaryCond.cpp:244) は `inletProfile` フラグだけを見ており、入口種別に限定していません。`valueTypes==1` の `Ts` は既存の補間・H2D 対象になります。これは新規一般化の範囲を小さくできます。

   しかし補間位置は [boundaryCond.cpp:312](/home/sano/work/forge/solver_density_cuda/boundaryCond.cpp:312) の**面重心**で、温度を課す位置は壁ノードです。実際に `case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/mesh.h5` では、壁ノード `x=1 m` に対応する境界面重心は `x=0.999320626 m`、差は **0.679 mm** でした。

   壁ダンプの `Ts` だけを見る verify も不十分です。出力は `bvar` の再出力なので、入力が入ったことと、流体場のピンが効いたことを区別できません。

   **対案:** node はノード座標、cell は面座標で評価し、CHT 内部転送は安定したノード ID を正本にしてください。非一様プロファイル試験では、`VALUE/T` と EOS 整合まで確認する必要があります。

6. **Major — 境界面接続の存在だけでは、シェルの幾何・接合条件は決まりません。**

   **根拠:** [gmshReader.hpp:1618](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:1618) は半割面ベクトルをノードへ合算し、[同:2250](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:2250) はその合成ベクトルの大きさを面積にしています。角・曲面では一般に

   \[
   \left|\sum_f\mathbf S_f\right|\ne\sum_f|\mathbf S_f|
   \]

   です。単一ノード温度の共有はできても、各面の法線・面積・第一内部点・熱荷重の帰属まで一意にはなりません。周期接合と、軸対称シェルの伝導作用素にも規定がありません。軸対称は最終積分だけに `r` を掛ければ済む問題ではありません。

   **対案:** primal facet ごとの幾何と荷重を保持し、温度 DOF の共有とは分離してください。第一内部点の選択、非直交補正、端部 BC、周期同一視、軸対称の面内伝導を定義し、未対応構成は起動時に拒否するべきです。角ノード試験は #10 から Phase 0 へ前倒ししてください。

7. **Major — V1 の解析解が逆で、試験問題も未確定です。**

   **根拠:** [plan:279](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:279) の式は抵抗の重みが逆です。正しくは

   \[
   T_w=\frac{T_{aw}/R_f+T_b/R_s}{1/R_f+1/R_s}.
   \]

   CPU 上で `Rf=0.002`、`Rs=0.005 m²K/W`、`Taw=1000`、`Tb=300 K` を代入すると、計画式は **Tw=500 K、qf=250、qs=40 kW/m²**。正しい式は **Tw=800 K、両側とも100 kW/m²** です。

   また「既知 h の層流平板」は非一様壁温で局所独立の抵抗になるとは限らず、Couette は粘性発熱の扱いが必要です。

   **対案:** V1 は定物性・静止流体層＋固体層の一次元純伝導に固定し、厚さから `Rf` を厳密に定義してください。フィン問題は格子収束で連続解析解への接近を調べ、離散残差の機械精度検査と分けるべきです。

8. **Major — V2 は実現可能ですが、現在の記述では同一問題の比較になりません。**

   **根拠:** SU2 CHT は固体ゾーン・界面マーカー・ゾーン別設定を要します。現行 `procedures/su2-cross-check.md` はその手順を提供しておらず、case/48 の既存 SU2 ケースも等温流体ケースです。[SU2 公式 CHT 手順](https://su2code.github.io/tutorials/Static_CHT/)も複数ゾーンを前提としています。

   さらに、forge の縮約シェルと SU2 の体積固体では、同じ物性・流体メッシュでも同じ離散問題ではありません。シェル近似誤差をソルバ差へ混ぜることになります。

   SU2 の先例にも注意が必要です。ローカルの [CNSSolver.cpp:593](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp:593) の引数順は `thermal_conductivity, dist_ij` ですが、[呼出し:718](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp:718) は逆です。これは引用した `AVERAGED_*` の係数に影響します。V2 指定の `DIRECT_*` 分岐には、この係数誤りは直接作用しません。

   **対案:** V2 をまず V1 と同じ一次元スラブで成立させ、SU2 の版・バイナリ、固体メッシュ、界面対応、物性、両ゾーンの残差を固定してください。case/48 への拡張では固体格子収束とシェル近似誤差を別に測るべきです。

9. **Major — 指定された既存 run は、予定している界面量監査の材料が揃っていません。**

   **根拠:** 保存済み壁ダンプを全時刻確認しました。

   | run | 壁ダンプ | `qwall` |
   |---|---:|---:|
   | `case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/` | 9枚・9,009値 | 全値0 |
   | `case/44.vitiated_air_wt/run_0114_va_ns_fine_iso300_samewall_cont/` | 9枚・36,009値 | 全値0 |

   `T1,d1,k_eff` と拘束反力も保存されていません。現在のソースに診断書込みがあることは、旧 run にその値がある証拠になりません。

   `check_convergence.py` を各継続 run の本段 CSV に実行した結果は、**両方 `NOT CONVERGED`**。run_0114 は `rms_roUy`、`rms_roOmega` が `RISING` でした。収束済み参照による `--from-floor` 判定は別途必要です。

   一方、run_0011 の `x≈0.5 m`、9時刻に対して、コンパクト差分・二次片側差分・`Tw` は `check_quasisteady.py --series-csv` で **`STEADY`**。最終値はそれぞれ **96.184 / 98.820 kW/m²** で、約 **2.67%** 違います。二次片側差分はカーネルの再構成勾配とは別物ですが、差分形式の選択だけでも V2 の3%許容に近い差が出ます。

   **対案:** 診断実装後、既存 run を入力参照にした新規 run で監査してください。旧 run は準備用資料に留め、収束と界面精度の合格根拠にはしないでください。元の索引は [case/48 README](/home/sano/work/forge/case/48.flat_plate_cooled_m4/README.md:37)、[case/44 README](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:695) です。

10. **Major — 既存ゲートへ列を追加するだけでは、界面収束を検査できません。**

    **根拠:** [check_convergence.py:45](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:45) は既知の保存量列だけを取り込みます。`max_delta_Tw` や任意の界面残差列を出しても検査されません。

    また、小さい緩和係数や float32 の丸めで `ΔTw=0` になっても、熱収支は閉じていない可能性があります。**float32 の温度を host double に変換しても、失った温度差は戻りません**。1000 K 付近の1 ULPは約 `6.10e-5 K` です。正負の熱流が相殺する場合、`|ΔQ|/Q` も破綻します。

    [stage_manifest.py:75](/home/sano/work/forge/solver_density_cuda/tools/stage_manifest.py:75) は設定テキストしか識別せず、外部 `wallProfile` や `solid.json` の内容変更を捉えません。

    **対案:** 既存ツールを拡張し、未緩和の局所界面不釣合い・固体方程式残差・絶対温度更新量を明示的に判定してください。絶対・相対許容値を併用し、必要熱流束を温度量子化で解像できるか確認するゲートも必要です。区間識別には外部入力の内容ハッシュと実効連成状態を含めてください。

11. **Minor — 作戦選択の結論は支持しますが、不採用理由と性能効果は言い過ぎです。**

    **根拠:** [plan:74](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:74) の「restart を挟まないので N 倍のコストが消える」は成立しません。消えるのは起動・I/Oなどで、流体の熱場緩和は残ります。

    また、界面に共有温度未知数を置くこと自体は、伝導の温度連続と保存的流束の表現を禁止しません。forge の保存変数・EOS・輸送・陰解法に固体を混ぜる費用が大きいことが、不採用の十分な理由です。

    **対案:** A→B と monolithic 不採用は維持し、理由を既存実装への侵襲性と検証費用に修正してください。Phase 2 の速度向上と host 求解費用は測定項目に変更すべきです。

12. **Minor — docs 更新の順序がリポジトリの開発フローと逆です。**

    **根拠:** [plan:244](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:244) は仕様文書を実装・検証後へ置いていますが、[AGENTS.md:155](/home/sano/work/forge/AGENTS.md:155) は仕様更新後の実装着手を求めています。

    **対案:** 界面契約・固体式・対応範囲を実装前に `methods/` へ反映し、完了時は実装との整合確認を行ってください。

**推奨は、A→B を維持し、Phase 0 を「保存的な界面契約の確立」に組み直してから再レビューすることです。** 優先順は、①実効壁熱量と符号、②固体抵抗と固定点を保存する Robin 式、③等温壁経路・座標・接合の配管、④V1/V2とゲート、です。

目的の切出しには賛成です。`plans/accepted/` に既存 CHT 実装は見当たらず、親 plan の弱連成項目を引き取る関係も確認できました。node のみの検証も現行規則と整合しています。ただし、現稿のまま実装へ進むことは推奨しません。レビュー提案は **plan 未反映**です。

指摘数: Critical 1 / Major 9 / Minor 2
