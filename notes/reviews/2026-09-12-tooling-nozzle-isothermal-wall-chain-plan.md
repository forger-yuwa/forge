# codex レビュー: tooling-nozzle-isothermal-wall-chain (plan)

- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **stage**: `plan`
- **date**: 2026-09-12
- **commit**: `0512823d` (feature/sern-design)
- **codex**: effort `high`, 5.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m2
- **extra**: `methods/design/overview.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
課題の同定と「低Re SSTの平板検証→ノズル比較→設計チェーン」の順序は妥当です。  
ただし、符号付きδ*、熱流束の取得、合否判定、壁温感度の実験設計に欠落があり、現状のまま実装へ進むことは推奨しません。

`plans/README.md` と accepted plan を確認しました。親計画に残る「等温壁のCFD検証」を引き受ける位置づけは正しく、壁関数の別計画とも対象が分かれています。以下はファイルを変更せず、コード読解・既存runの再判定・純関数の実行で確認した指摘です。

1. **Major — 「δ*抽出はそのまま使える」は現行実装と矛盾する**

   **根拠:** [plan:52](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:52) は負のδ*を物理として認めています。しかし、生産抽出器は [deltastar.py:379](/home/sano/work/forge/design/forge_design/metrics/deltastar.py:379) で `negative_deficit` を hard不合格にします。また、[同:443](/home/sano/work/forge/design/forge_design/metrics/deltastar.py:443) の平滑化は既定 `positive=True` で、[同:475](/home/sano/work/forge/design/forge_design/metrics/deltastar.py:475) が負値を0にします。実行確認でも、全点 `δ_r=-0.01` の入力が全点0になりました。初期壁・反復更新の双方がこの既定を使います。

   これは既存の[親plan:124](/home/sano/work/forge/plans/accepted/tooling-nozzle-deltastar-core-matched-euler.md:124)に沿った動作であり、単なる配管不足ではありません。なお、引用された `deltastar_from_run` は生産の帯局所抽出器ではありません。

   **対案:** S1に符号付きδ*の設計変更を含め、抽出の合否、負欠損の等価半径変換、平滑化、壁更新まで一貫して扱ってください。正・零・負を横断する質量欠損の回帰を追加します。今回の300 K条件で負値が発生することを確認した、という指摘ではありません。

2. **Major — 低Reの熱流束取得方法が未設計で、`Q_w` の積分式も不正確**

   **根拠:** [plan:69](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:69) は熱流束の出所を定義していません。既存の `qwall` は壁モデル用で、低Re経路では解像熱流束を書き込みません。[ransWallFunction_d.cu:107](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransWallFunction_d.cu:107)、[viscousFlux_d.cu:526](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:526)で確認できます。

   実ファイル `case/40.nozzle_design_tool/run_0048_node_yp1_isoT_ref/res_wall_3_12000.h5` でも、`Ts=1000 K` に対して **`qwall` は全点0**でした。これを積分しても低Reの熱負荷にはなりません。

   また、壁面積当たりの熱流束なら軸対称積分は
   \[
   Q_w=\int q_w\,2\pi r\,ds,\qquad ds=\sqrt{1+(dr/dx)^2}\,dx
   \]
   です。計画の `2πr dx` は傾斜壁の面積を過小評価します。nodeでは壁エネルギー残差が[強制的に0化](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:85)される点も、熱収支評価に関係します。

   **対案:** S1で低Re専用の熱流束抽出、必要な出力場、法線・符号・面積重みを定義してください。局所勾配による熱流束と、Dirichlet拘束が除去するエネルギー収支を照合し、入口出口の全エンタルピー流束との差との閉合も検証します。平面形状では単位幅の積分に分岐させます。

3. **Major — 現行の収束・準定常ゲートでは、計画の合格を判定できない**

   **根拠:** `check_convergence.py` を既定条件で再実行した結果は次のとおりです。

   | 基準run | VERDICT | 代表的な最終残差 |
   |---|---|---|
   | `case/45.isobutane_m6_d155/run_0013_cpg_ns_plainsst/` | `NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)` | `rms_roe=1.03`、低下2.1桁 |
   | `case/44.vitiated_air_wt/run_0107_va_R2_LU6_Lc8_ns_ib_pass0/` | 同上 | `rms_roe=0.341`、低下2.5桁 |
   | `case/40.nozzle_design_tool/run_0048_node_yp1_isoT_ref/` | 同上 | `rms_roK=0.0124`、ピークから低下0.8桁 |

   これだけで場が誤りとは断定しません。しかし、[plan:105](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:105) の「全列 falling/flatなら可」には残差の絶対水準がなく、高い残差の停滞も通ります。

   準定常コマンドも実行すると、`--quantity cf,theta` は **`unknown quantity ['cf']`**、`--series` は **`unrecognized arguments`** でした。[check_quasisteady.py:282](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:282)、[同:419](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:419)。既定のドリフト許容5%も、出口M±0.1%やθ比較1%には粗すぎます。

   **対案:** 検証ツールの拡張をS1の前提に追加し、`q_w`、`Q_w`、δ*、θ、出口コアM、流量比を実際の報告定義で時系列判定してください。時間変動の許容は各比較公差より十分小さく設定します。warm床を認める場合も全列の定量的な床と物理量の定常性を事前規定し、元の `NOT CONVERGED` 表示は残してください。

4. **Major — §4.5の比較は、§4.6が必要とする「壁温だけの感度」にならない**

   **根拠:** [plan:117](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:117) は冷却条件で積分法初期壁を作り直してδ*反復します。一方、[同:127](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:127) は公称形状を固定した感度を要求しています。実装も[runner_axismach.py:689](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:689)で壁温に応じた壁補正を生成します。

   断熱設計と等温再設計が両方とも出口Mゲートを満たせば、形状補正が壁温影響を打ち消します。その差を「壁温感度が小さい」と読むことはできません。また、[同:131](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:131) の「風洞ではδ*経由でしか効かない」は、壁熱移動による全エンタルピー・全圧の変化を除外する根拠がありません。

   **対案:** S4を「同一の物理壁・メッシュで断熱／300 Kを比較」→「300 Kで再設計」の順に分けてください。台帳では固定形状感度と再設計による回復量を別項目にします。積分法は候補選別に限定し、最終採否はエネルギーを含むNS評価に置きます。

5. **Major — 経験的な温度関係を3%の厳密ゲートにする根拠が足りない**

   **根拠:** [plan:102](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:102) の温度式は端点整合しています。しかし、強冷却でも3%以内という保証は別問題です。[既存積分法:49](/home/sano/work/forge/design/forge_design/feedback/deltastar_integral.py:49)も同じ形を使うため、CONTURとの比較は独立した温度モデル検証にはなりません。

   冷却壁向けのDuan–Martín型関係は、線形項に `C_T=0.8259` を使います。[Chen・Gan・Fuの一次論文、式1.1a](https://doi.org/10.1017/jfm.2025.312)。計画の条件を代入すると `u/u_e=0.5` で、計画式は **514.15 K**、同関係は **476.13 K**。近似関係同士ですでに約8%違います。これはCFD誤差ではありません。

   **対案:** 古典式との3%一致は診断に変更し、冷却壁データに基づく関係を独立参照として追加してください。評価位置・誤差の分母・適用範囲を明記し、CONTURの精度評価とソルバの合否を分離します。CONTUR平面版は単に巨大な `r_w` を渡さず、一定外縁条件で `dθ/dx=C_f/2` を積分する入口を用意する方が明確です。

6. **Major — 「同一indexなら起動に問題なし」は、壁温変更と再メッシュを考慮していない**

   **根拠:** [plan:48](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:48) は安全性を一般化していますが、[同:111](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:111) は第一層と `ni` の変更を予定しています。`ni` を変えれば旧runとのindex対応は成立しません。現在の `prepare_ns(ic_from=...)` は[runner_axismach.py:784](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:784)で `interp_field.py` を実行します。

   さらに温度ピンは密度を保持して `P=ρRT_w` を再設定します。[nodeWallDirichlet_d.cu:75](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:75)。同一メッシュでも、例えば1000→300 Kの変更直後には壁圧が元の0.3倍になります。indexコピーだけではこの過渡は消えません。

   **対案:** 同一メッシュでの温度変更、同一トポロジの形状変更、解像度変更の引き継ぎ手順を分けてください。温度変更は段階起動または壁温continuationを規定し、再メッシュは両ソルバの移植方法を明示します。step 0/1の壁圧・密度・温度・ωを起動ゲートに追加してください。

7. **Major — SU2との比較行列が「同じSST」を保証していない**

   **根拠:** [plan:89](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:89)以降では、素SST版は `B-plain` だけで、`SU2-A` に対応する `A-plain` がありません。また、`dilatationCorrection` と `katoLaunder` の2キーだけではモデルの全条件を固定できません。

   現行コードでは `sstOmegaProdFromPk` と `sstSigmaBlend` が2026-09-08に既定1へ変更されています。[solverConfig.hpp:206](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:206)。旧基準runの設定はこれらを明示していません。エネルギーへのkの取り込みも[別の選択肢](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:210)です。旧runの比較誤差を、そのまま現行バイナリの実績にはできません。

   **対案:** `A-plain/B-plain` と `SU2-A/SU2-B` を基準対にし、生産設定との差を別途測ってください。SST変種、エネルギー定式化、入口k/ω、物性、乱流輸送の離散化を対応表で固定し、現行バイナリの断熱基準を取り直します。1～3%比較には、壁法線だけでなく流れ方向の格子感度確認も必要です。

8. **Minor — 「単一ソース」と上書き仕様、回帰対象が整合していない**

   **根拠:** [plan:68](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:68) は `--init-thermal` による初期壁だけの上書きを残すため、NSと積分法が異なる壁温を読む状態を作れます。また、[同:160](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:160)の回帰対象 `run_0037` と `run_0091` はEuler基準で、NSの `wall` 分岐の互換性を確認できません。

   **対案:** 上書きは有効壁温を全段へ適用し、不一致の初期化実験だけ明示的な例外として記録してください。省略／断熱／等温×Euler／NSの生成結果と、NS基準のbcond・初期壁を回帰対象にします。

9. **Minor — y⁺倍率の算術が合っていない**

   **根拠:** [plan:75](/home/sano/work/forge/plans/active/tooling-nozzle-isothermal-wall-chain.md:75)の係数なら、
   \[
   \frac{y^+_{\rm cold}}{y^+_{\rm ad}}
   =\frac{\sqrt{3.9(1.2\text{～}1.4)}}{0.39}
   =5.55\text{～}5.99
   \]
   で、記載された4～5倍ではありません。

   **対案:** 見積りを修正し、第一内部nodeまでの法線距離と実測壁物性を用いてメッシュを決めてください。許容上限超過は生産ゲートでは不合格とし、前縁・角点の扱いを別途定義します。

**推奨は、低Re方針を維持して計画を修正することです。** 実装前に、優先順として①符号付きδ*、②熱流束と熱収支、③実行可能な定量ゲート、④固定形状の壁温比較、⑤独立参照・SST条件・引き継ぎ手順を§4・§5.1・§6へ反映してください。その後にS1→S2→S3→S4へ進めます。ソルバ共有部の変更が必要になった場合は、node/cell双方の検証を追加する条件も明記してください。

本レビューはread-onlyのため、**plan未反映**です。

指摘数: Critical 0 / Major 7 / Minor 2
