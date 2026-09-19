# codex レビュー: case-plate-annular-cavity-m5 (result)

- **plan**: [`plans/active/case-plate-annular-cavity-m5.md`](../../plans/active/case-plate-annular-cavity-m5.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-19
- **commit**: `fde8f3fc` (feature/sern-design)
- **codex**: effort `high`, 9.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M10/m1
- **focus**: 重点 2 点。(1) **熱伝達率の基準温度に使う「すきま中央面」が意図どおりか**: ユーザ指定は『すきま中央断面の最近傍点の全温を駆動温度、断面の底面はすきま幅の半分だけ切り取り』。実装は tools/cavity_eval.py の mid_surface() と wall_href()、面の形状は tools/geom_common.py の gap_center_radius/inner_radius_at (偏心対応、θ=0 が上流)。偏心 x_off>0 のとき中央面が本当に内外壁の中間に来ているか、底の切り取り量が すきま幅の半分 (θ 依存) になっているか、壁ノードから最近傍点を引く cKDTree の使い方が妥当か、T0 を VALUE/h0 から total_quantities.total_state で逆算している点が AGENTS.md の規約に合っているかを、コードを読んで検証してほしい。h_ref がノードごとに q''/(T_ref-T_w) を取って面積平均する定義であることの妥当性 (面積平均 q'' / 面積平均 dT との違い) も評価してほしい。(2) **収束確認がちゃんとできているか**: 段階起動のため本段だけの check_convergence.py は必ず plateau になる。tools/stack_residuals.py で全段を連結して判定し (run_0103_hex_s14 は起動全体で 3.3〜8.0 桁低下)、結論量は tools/check_cavity_steady.py で準定常判定している (最細 run で h_ref drift 0.0%/tail STEADY、底面 dT のみ TRANSIENT-UNSETTLED)。この 2 段構えで AGENTS.md の『収束確認』『準定常確認』の要求を満たしているか、連結履歴で判定することの落とし穴 (段ごとに数値設定が違う) をどう扱うべきかを見てほしい。格子収束は tools/compare_runs.py --richardson (ASME/Celik) で 3 点、開口リップの特異点のため総 Q と q''max は収束しない (GCI 32%/224%) と結論し、リップ帯を除いた Q (GCI 0.31%) を収束量としている。この扱いが妥当かも見てほしい。
- **extra**: `case/49.plate_annular_cavity_m5/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

偏心時の中央面と開口積分に座標の符号誤りがあり、`h_ref` と偏心比較の再評価が必要です。収束・準定常・格子収束の判定にも合格を裏付けない経路があり、現状では `accepted` に移せません。

指定された `git diff main...HEAD -- solver_density_cuda design methods procedures` と履歴を取得し、case 側の実装も確認しました。ただし、この作業ツリーにある実測は `run_0001`〜`run_0006` と前駆計算です。`run_010x`・`run_011x`・`run_020x` は存在せず、それらの数値は原データで再検証できていません。以下では、実測の再確認と文書記載値の検算を区別します。ファイルは変更していません。

1. **Critical — 偏心時の中央面が逆向きに配置され、固体内へ入る**

   **根拠:** [geom_common.py:69](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/geom_common.py:69) の `inner_radius_at()` は上流基準で正しく、
   \[
   r_i=-x_{\rm off}\cos\theta+\sqrt{R_i^2-x_{\rm off}^2\sin^2\theta}
   \]
   を返します。しかし [cavity_eval.py:239](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:239) は、それを `x=+R*cos(theta)` で座標化しています。必要なのは `gc.ray_dir()` の **負の cos** です。[開口積分の199行](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:199) にも同じ誤りがあります。

   関数を実行した結果、`x_off=1.5 mm` では次の誤配置になりました。

   | 物理的位置 | 正しい中央面下端 `(x,z)` [mm] | 実装 |
   |---|---:|---:|
   | 下流側 | `(24.5, −49.5)` | `(23.0, −48.0)` |
   | 上流側 | `(−23.0, −48.0)` | `(−24.5, −49.5)` |

   中央面は **181 方位中56方位で固体内**に入ります。床からの切り取りも、下流で必要な **0.5 mm が2.0 mm**、上流で必要な **2.0 mm が0.5 mm** になります。開口積分の面積重みのうち固体内に置かれる割合は、偏心0.5/1.0/1.5 mmで **12.09/24.14/36.25%** でした。

   `cKDTree.query()` 自体は点群へのユークリッド最近傍検索として妥当ですが、検索対象が誤っています。さらに `Field` が固体内の点にもガス値を返すため、異常が隠れます。

   **対案:** 中央面・開口積分とも `gc.ray_dir()` に統一し、点群が流体内にあることを検査する。偏心runの `Tref`・`h_ref`・開口流束・収支を再生成してください。半径方向の中点という定義は妥当ですが、偏心時の「両壁への最短距離が等しい面」とは区別して記述すべきです。

2. **Major — `h_ref` の分母を1 Kへ切り上げており、指定された定義ではない**

   **根拠:** [cavity_eval.py:259](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:259) は
   `dT = np.maximum(Tref - Tw, 1.0)` としています。小さい正値だけでなく、負の温度差も正の1 Kに置換します。

   `case/49.plate_annular_cavity_m5/run_0003_stageB_cpg/res_24000.h5` を再評価すると、この補正の対象面積は外筒 **21.22%**、内円柱側面 **21.12%**、床 **15.31%**。無視できるゼロ割防止ではありません。

   また、[plan:662](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:662) の「基準温度は2%しか変わらないので分母の効果ではない」は誤りです。分母は絶対温度ではなく、**524.1−500=24.1 K → 513.6−500=13.6 K、43.6%減少**です。

   面積平均した局所係数
   \[
   \langle h\rangle_A=\left\langle q''/\Delta T\right\rangle_A
   \]
   はユーザー指定の定義として成立します。ただし、総入熱を再現する係数
   \[
   h_{\rm eff}=Q/\!\int_A\Delta T\,dA
   \]
   とは異なります。[compare_offsets.py:54](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/compare_offsets.py:54) の `h_ref@ref` は壁面平均温度差で割っており、局所分母を固定した比較になっていません。

   **対案:** 符号を保存し、温度差が評価誤差以下の点は「係数を定義できない」として面積割合を報告する。局所係数平均と `h_eff` を別名で出し、分母寄与は対応する壁位置の `ΔT_ref` を固定して比較する。

   なお、[total_state の呼び出し](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:429) は `VALUE/h0` を使っており、この点は AGENTS.md に適合しています。中央面の有限点群への近似にも感度があり、Stage Bで分点数を181×200から721×800へ増やすと床の実装値は **44.541→44.903 W/m²K** に変わりました。

3. **Major — 異なる方程式・境界条件の履歴を連結しても、本段の収束証明にはならない**

   **根拠:** [stack_residuals.py:25](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/stack_residuals.py:25) は slip・層流・等温化・SST・次数変更を連結します。一方、[check_convergence.py:105](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:105) は系列全体の最大値を低下桁数の基準にします。別の境界条件や物理モデルの大きな過渡を、本段の基準値にできます。

   ローカルデータをメモリ上で連結し、正本ツールの `analyze()` で再判定しました。

   | run | S6単独 | 全段・共通保存量 | 最終空間設定のS5以降 |
   |---|---|---|---|
   | `case/49.plate_annular_cavity_m5/run_0003_stageB_cpg/` | `NOT CONVERGED` | `NOT CONVERGED` | `NOT CONVERGED` |
   | `case/49.plate_annular_cavity_m5/run_0004_stageC_cpg/` | `NOT CONVERGED` | `NOT CONVERGED` | `NOT CONVERGED` |
   | `case/49.plate_annular_cavity_m5/run_0005_stageD_cpg/` | `NOT CONVERGED` | `NOT CONVERGED` | `NOT CONVERGED` |

   全段でも Stage C の `rms_roUy` は **1.7桁**、Stage D は **2.5桁**低下にとどまります。「本段だけだから不合格」という説明では解消しません。

   **対案:** 全段履歴は起動診断用に残し、収束判定は同じEOS・BC・乱流・空間離散化の区間で行う。本件なら設定の同一性を確認したS5以降が候補です。CFL変更だけの連結は可能ですが、残差の定義と正規化が同じことを記録してください。既収束場からの継続は、検証済み参照に対する `--from-floor` を使うべきです。

4. **Major — 準定常ツールが欠損・NaN・必須量を見落とす**

   **根拠:** [check_cavity_steady.py:25](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_cavity_steady.py:25) の既定対象には壁別入熱・`mdot_out`・`mdot_imbalance`・`Q_nolip` がありません。[60行](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_cavity_steady.py:60) は欠けた列を無言で飛ばし、63行でNaNを除去する `classify()` を呼びます。

   再実行結果は以下です。

   - Stage B：`q_wall_sum` と `h_ref` の列がなくても **`VERDICT: STEADY (全量)`**。
   - Stage D：`h_ref` は12点すべてNaN。case側は **`TRANSIENT-UNSETTLED`**、正本のCSV判定は **`NONFINITE`**。
   - `[1,1,1,1,1,NaN]` という最小例も、使用中の `classify()` は **`STEADY`** を返します。

   計画の絶対許容値も読み込んでいません。さらに、単に系列を626 Kなどで割っても、平均値で正規化する判定では係数が相殺され、ゼロ近傍の問題は解消しません。

   **対案:** `classify_series()` または正本の `--series-csv` 経路を使い、必須列の欠損・非有限値を失敗にする。固定基準尺度によるdrift、2 K・0.5 mm・入熱1%・質量不釣合い0.1%を明示的に判定する。現状の二段構えは AGENTS.md の要求を満たしていません。

5. **Major — 「局所熱流束が特異だから総入熱も原理的に収束しない」は誤り**

   **根拠:** [plan §4.4.2](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:333) と [cavity_eval.py:461](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:461) は両者を同一視しています。しかし、主張されている特異性そのものでも
   \[
   q''\sim s^{-1/2}
   \quad\Rightarrow\quad
   \int_0^\varepsilon q''\,ds\sim2\sqrt{\varepsilon}
   \]
   で、積分は有限です。

   **局所最大値の発散と、総入熱の有限値への遅い収束は両立します。** 記載された総Qの観測次数0.24・GCI32%は、総Qの発散の証明ではありません。

   **対案:** 結論を「総Qの格子不確かさが大きく、現時点では精度を確定できない」に訂正する。リップ除外Qは限定領域の別指標として併記し、目的にある総入熱の代替として完了扱いしないでください。

6. **Major — GCIの適用不可判定を無視し、`0.00%` を収束の証拠にしている**

   **根拠:** [compare_runs.py:83](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/compare_runs.py:83) は振動系列や異常次数でもGCIを計算・表示します。README記載の中央温度 **511.911→511.915→511.740 K** を同じ関数に入力すると、**p=12.11、判定「振動」、GCI表示0.00%** になります。これは原データの再検証ではなく、掲載値を使った再現ですが、「中央TのGCI0.00%=収束」という解釈が成立しないことは確認できます。

   また、[cavity_eval.py:324](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:324) の `Q_nolip` はノードの `z<-lip` だけで面積重み全体を採否判定します。境界を横切る面を切断しないため、格子ごとに実効的な積分範囲がずれます。

   **対案:** 適用条件を満たさない系列はGCIを「推定不可」とする。固定した物理領域を面単位で切断積分し、反復・準定常・後処理の誤差を抑えた上で格子列を再評価してください。3格子は入力条件であり、漸近域の確認も必要です。[NASAの格子収束手順](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html)

7. **Major — CV閉性から、収支誤差の原因を後処理だけに断定できない**

   **根拠:** [plan:492](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:492) の断定に対し、[cavity_eval.py:223](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:223) は離散保存の独立検算が未実装と明記しています。幾何学的閉性は、収束した離散エネルギー収支の証明ではありません。

   実装の開口エネルギー流束は [217行](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:217) の `ΣρUz h0 dA` だけで、内部評価面を通る熱伝導・粘性仕事を含みません。さらに、構造格子では `zpick` の値を使う一方、壁積分の切断位置には指定値 `−zf` を使っています。

   README自身の記載でも、最細の収支 **6.7% > 5%**、質量不釣合い **0.26% > 0.1%** でゲート未達です。ローカルStage Dの保存CSVでも、末尾の `mdot_imbalance` 平均は **−0.02446** でした。

   **対案:** 同一の検査面・非拘束CV群について離散流束と蓄積項を収支計算する。物理面積分では伝導・粘性仕事を含めるか、その寄与が小さいことを測って示す。原因の断定はそれからです。

8. **Major — 旋回平均の終点一致だけで「半割で可」を返している**

   **根拠:** [check_symmetry.py:205](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_symmetry.py:205) は、旋回平均の最初と最後の基準値との差だけで成功を返します。途中の成長判定、末尾振幅、温度・熱流束の非対称性は、この成功経路を止めません。

   したがって **0.0077** は終点で基準値に近いことを表すだけで、減衰包絡線や他の反対称モードの安定性を保証しません。平均が相殺される非対称モードも検出できません。定常解の非ゼロ旋回を「メッシュ床」として引くことも、その原因を独立に検証しない限り循環論法です。

   加えて [plan:750](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:750) が必須とするdt半減・`nSub`倍増・細分・観測窓倍増は、§5.1 #9の残作業から抜けています。

   **対案:** 末尾窓の反対称成分のノルム・振幅・包絡線で判定し、計画済みの時間・空間感度を完了させる。それまでは「この擾乱・離散条件で旋回平均が基準へ近づいた」に結論を限定してください。

9. **Major — メッシュ品質と壁解像の合格条件が実装と一致しない**

   **根拠:** 計画のCV閉性閾値は [283行](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:283) の `1e-6`、実装既定は [check_mesh_extra.py:60](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_mesh_extra.py:60) の `5e-5` です。

   `case/49.plate_annular_cavity_m5/mesh/stageD.h5` の再実行結果：

   - `check_mesh_quality.py`：**`VERDICT: PASS`**、AR127.6・skew0.849。
   - `check_mesh_extra.py --closure-max 1e-6`：**`VERDICT: FAIL`**、最大 **5.114e−6**、超過 **22,019 CV**。

   また、prism検査は頂点Jacobianや隣接層厚比の検査になっておらず、上下面積比の超過数さえ失敗条件に入りません。[check_mesh_extra.py:103](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_mesh_extra.py:103)

   壁解像も、評価しているのは壁ダンプの `ypls` です。その生成式は **ghost間距離 `dcc` と `mu_total`** を使っており、計画が要求する第一内部ノードの壁距離による `y₁⁺` ではありません。[viscousFlux_d.cu:573](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:573)

   **対案:** 合格条件を計画と実装で統一し、閾値緩和は解への影響で裏付ける。未実装の品質検査を補い、`y₁⁺` は実際の第一内部ノードから算出してください。

10. **Major — run固有の条件を保存しても、評価時に共有条件へ戻る**

    **根拠:** [cavity_eval.py:368](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:368) はrunのmanifestを読みますが、371行の `load_conditions()` は [setup.py:92](/home/sano/work/forge/case/49.plate_annular_cavity_m5/setup.py:92) 経由で共有 `case.json` を読み直します。Mや壁温を変更すると、過去runの評価が変わります。生成側も [gen_runs.py:45](/home/sano/work/forge/case/49.plate_annular_cavity_m5/gen_runs.py:45) でmanifestと共有条件を別々に読みます。

    さらに [cavity_eval.py:301](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:301) は、CPG runでも存在する `Taw_tp` を優先します。今回の温度差では、CPG基準の `h_aw` より約 **9.8%大きく**なります。

    **対案:** 条件・壁温・EOSをrunに固定されたmanifest/configから解決する。`h_aw` は実際のEOSに対応する回復温度を選択し、共有入力の変更で過去runの評価が変わらないようにしてください。

11. **Major — 必須検証が未完了で、残作業表の「完了」が実態を表していない**

    **根拠:** [§5.1](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:701) 自身が、塞ぎ3DのBL比較、TP本計算、領域・入口BL・緩和感度を未実施としています。これは§8の完了条件と両立しません。

    一方、評価ツールを完了扱いにしたため、ソースに残る離散保存検算、解析解による評価検証が残作業表から見えなくなっています。全周URANSの必須感度も同様です。段間移行も [gen_runs.py:309](/home/sano/work/forge/case/49.plate_annular_cavity_m5/gen_runs.py:309) は主にstep到達だけを見ており、計画の残差・物理範囲・キャビティ蓄積量による移行判定ではありません。

    **対案:** §5.1の評価・URANS・格子収束を再開し、抜けた検証を優先順付きで戻す。TPとBL移送を含む必須ゲートを完了し、対象runの原データ・VERDICT・評価条件をレビュー可能な場所に揃えてください。

12. **Minor — ソルバ変更の影響範囲と現在仕様の記載が未同期**

    **根拠:** [plan §7](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:762) は「ソルバコード変更なし」としていますが、実際には [viscousFlux_d.cu:565](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:565) に診断出力が追加されています。`methods/`・`procedures/` には対応する差分がなく、[plans/README.md:26](/home/sano/work/forge/plans/README.md:26) も生産ヘキサへの変更を反映していません。

    追加された `heatflux/sss` の単位・符号と、壁関数経路を避ける条件はコード上整合しています。ただし、nodeの壁温ピン後の実効授受と等しいことや「解はビット不変」は、今回の確認だけでは実測保証できません。

    **対案:** §7と現在仕様に出力の意味・適用経路を追記し、対象となるnode壁処理の回帰結果を残す。「ビット不変」は比較証拠を付けるか、コード上の意図として表現してください。

**推奨は、計画を `active` に留め、評価系を修正してから必須検証を完了し、resultレビューを再度受けることです。** 優先順は、①偏心座標と `h_ref` の修正、②欠損を合格にしない判定系と再後処理、③同一方程式での収束・離散収支確認、④TP・BL・URANS・感度ゲートの完了です。本レビューの指摘・採否は、ファイル変更禁止に従い **plan未反映**です。

指摘数: Critical 1 / Major 10 / Minor 1
