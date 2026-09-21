# codex レビュー: discretization-node-face-weight-midpoint (plan)

- **plan**: [`plans/active/discretization-node-face-weight-midpoint.md`](../../plans/active/discretization-node-face-weight-midpoint.md)
- **stage**: `plan`
- **date**: 2026-09-21
- **commit**: `18916ed9` (feature/sern-design)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m0
- **focus**: 問いは 2 つ。(1) §4.2: node の内部双対面で fx=0.5 に固定してよいか。fx を読むのは粘性流束の面値と面勾配補間・リミッタ・Ducros・interpVelocity_c2p。中点固定で精度や保存性を損なう箇所が無いか、伸長格子 (壁法線方向の成長率 1.1–1.2) で fx=0.5 が 2 次精度を保つかを見てほしい。(2) 同梱の判断 (plans/active/boundary-conjugate-heat-transfer.md §5.1 #61。grep -n '^| 61 ' で該当行のみ読むこと): 等温壁 Dirichlet 節点の保存的壁熱流束を q_eff=(Rraw-Fw-H_w*R_rho)/A に改め、未収束の擬似時間で壁 CV の質量残差 R_rho が運ぶ全エンタルピー H_w=(roe+P)/ro を蓄積項として引く。実測は対流エネルギー残差と R_rho の相関 1.000・傾き 571 kJ/kg (c_p T_w=568.5)。引く係数は H_w でよいか、内部エネルギー e_w=c_v T_w にすべきか (差 R T_w R_rho は等温圧縮の流動仕事)。定常推定量としてどちらが妥当か、推奨を 1 つに絞ってほしい。実装は solver_density_cuda/conjugateWall.cpp の sed -n '183,225p' のみ読むこと。禁止: *.log, residual_history.csv, res_*.h5, *.vtu, plans/README.md, plans/active/boundary-conjugate-heat-transfer.md の全文を読まない。
- **extra**: `solver_density_cuda/cuda_forge/calcStructualVariables_d.cu`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
`fx=0.5` の恒久化は推奨します。ただし、伸長格子での二次精度の根拠と影響範囲が不足しています。  
同梱の `H_w R_rho` を「保存的な蓄積補正」とする判断は誤りです。補正係数には条件付きで **`e_w`** を推奨します。

指定された禁止ファイルは読まず、ファイル変更もしていません。`plans/README.md` は個別の禁止指定を優先して未読です。以下の実測確認は既存 VERDICT に限定し、熱流束の数値・相関は提示値として評価しました。

1. **Major — `H_w R_rho` はエネルギーの蓄積項ではなく、圧力仕事まで除去する。**

   **根拠:** [CHT plan:445](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:445)、指定範囲の [conjugateWall.cpp:188](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:188) は、対流残差との相関を根拠に `H_w R_rho` を蓄積と解釈しています。しかし、実際の壁状態ピンは [nodeWallDirichlet_d.cu:83](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:83) の `roe = ro*(g.e + ek)` です。

   固定体積・静止壁・固定壁温・固定組成なら、壁 CV のエネルギーは \(V\rho e_w\)。半離散連続式 \(V\dot\rho=R_\rho\) に沿った蓄積は
   \[
   \frac{d(V\rho e_w)}{dt}=e_wR_\rho
   \]
   です。[既存の残差符号規約](/home/sano/work/forge/methods/boundary.md:537) に従えば、
   \[
   q_{\mathrm{eff},e}
   =\frac{R^{raw}-F_w-e_wR_\rho}{A}.
   \]
   `H_w` を使った値との差は
   \[
   q_{\mathrm{eff},H}-q_{\mathrm{eff},e}
   =-\frac{P_w}{\rho_w}\frac{R_\rho}{A}.
   \]
   この項は消してよい蓄積ではありません。固定体積への質量流入では、流入エンタルピーと内部エネルギー蓄積の差が、等温を保つための放熱になります。保存変数が \(\rho E\)、対流エネルギー流束が \((\rho E+p)\mathbf u\) であることとも整合します。[SU2 の支配方程式](https://su2code.github.io/docs_v7/Theory/)

   提示された相関 `1.000`・傾き `571 kJ/kg` は、**対流がエンタルピーを運ぶこと**の裏付けです。蓄積係数がエンタルピーである証拠にはなりません。実際、[run_0122 の設定:5](/home/sano/work/forge/case/53.c3x_vane_cht/run_0122_fxhalf/solverConfig.yaml:5) の \(c_p=1004.5,\gamma=1.4\) と提示壁温 `566 K` なら、\(h_w=568.55\)、\(e_w=406.11\)、差は **162.44 kJ/kg**。無視できる差ではありません。

   **対案:** この CPG ケースの残差補正係数は `e_w` にする。ただし、block-DPLUR は [5×5 系を解いて更新する](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1101) ため、`R_rho` は実際の反復更新 \(V\Delta\rho/\Delta\tau\) と一般には一致しません。`e_w R_rho` も「半離散式に基づく定常推定」と明記し、実際の反復収支には状態更新から蓄積を測ること。組成変化・壁温変化・`roe` に含む追加エネルギーがある場合は、対応する蓄積も必要です。

   `H_w` 補正も \(R_\rho\to0\) では同じ極限になりますが、未収束時の推定誤差が小さい保証はありません。スパイクが消えることを保存性の証明にしてはいけません。

2. **Major — 「辺中点の補間が二次精度」と「伸長格子上の演算子・解が二次精度」を混同している。**

   **根拠:** [対象 plan:67](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:67) の説明では不十分です。滑らかな関数について
   \[
   \frac{\phi_i+\phi_j}{2}
   =\phi\!\left(\frac{x_i+x_j}{2}\right)+O(|x_j-x_i|^2)
   \]
   は成立します。成長率 `1.1–1.2` でも、**各辺中点の値補間**として問題ありません。

   一方、[viscousFlux_d.cu:282](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:282) の直交・一定熱伝導率への縮約は二点差分です。節点を \((-a,0,b)\)、\(T=x^3\) とすると、中央節点の体積正規化演算子は
   \[
   L_hT=
   \frac{(T(b)-T(0))/b-(T(0)-T(-a))/a}{(a+b)/2}
   =2(b-a).
   \]
   \(b=ra\)、`r=1.2` のまま局所尺度を半減すると、実行した代数チェックで誤差は **0.04 → 0.02 → 0.01**。点値に対する局所打切り誤差は一次です。解の収束次数は別途格子系列で測る必要があります。

   また、双対面のパッチ重心と辺中点は一般に異なります。[幾何定義](/home/sano/work/forge/methods/discretization.md:59) に対して、辺中点の値が面積分まで無条件に二次精度になるわけではありません。

   **対案:** §4.2 を「辺中点で線形厳密・滑らかな場の値補間は二次」に限定する。§6 に、成長率 `1.1/1.2` の伸長格子、曲面高 AR 格子、回転した同一格子での製造解検証を追加する。少なくとも三格子で温度・速度・壁熱流束の誤差次数を測り、二次を主張する範囲では例えば \(p\ge1.8\) を要求する。細分化時に成長率をどう変えるかも明記すること。

   **保存性については変更を支持します。** 内部面流束は [viscousFlux_d.cu:380](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:380) で両 CV に同じ値を逆符号で加算しており、`fx=0.5` はこの保存構造を壊しません。float32 の集積丸めは別途残ります。

3. **Major — `fx` の参照箇所の調査が不完全で、回帰範囲が狭すぎる。**

   **根拠:** [plan §3:30](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:30) にない実使用があります。

   | 経路 | 根拠と影響 |
   |---|---|
   | SST 勾配 | [ransTransport_d.cu:39](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransTransport_d.cu:39)：`k/omega` の面補間。さらに [同:218](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransTransport_d.cu:218) は GG を実行しており、流れの LSQ 固定だけでは保護されない |
   | SST 拡散 | [scalarTransport_d.cu:218](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:218)：面係数と陰的対角が変わる |
   | 多成分拡散 | [speciesTransport_d.cu:229](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:229)：密度・温度・圧力・組成・輸送係数が変わる |
   | 受動スカラー／FCT | [passiveKernels_d.cuh:158](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:158)、[passiveFct_d.cuh:26](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveFct_d.cuh:26)：拡散係数が変わる |
   | CFL 評価 | [setDT_d.cu:76](/home/sano/work/forge/solver_density_cuda/cuda_forge/setDT_d.cu:76)：面速度・密度・渦粘性が変わり、反復経路も変わる |

   逆に `limiter_d.cu` と `ducrosSensor_d.cu` の `fx` は引数として渡されていますが、カーネル内で参照されていません。「引数にある」と「値を読む」を区別する必要があります。

   **対案:** §3・§7 に実使用一覧を反映し、[検証手順の node 一覧](/home/sano/work/forge/procedures/verification/README.md:10) に沿って、周期・受動スカラーの `case/09` と軸対称・多成分 TP・凝縮の `case/44` を追加する。平板は `case/48` など対象 run と壁処理を特定する。

   周期 node は [main.cpp:1463](/home/sano/work/forge/solver_density_cuda/main.cpp:1463) のとおり両側内部面の残差を gather するため、内部面限定の変更は構造と整合しています。ただし、周期保存性・軸近傍・SST の回帰を C3X と一次元 slab だけで代表させることはできません。

4. **Major — 合格条件が未収束同士の比較を許し、熱流束指標の定常性と定義変更を管理していない。**

   **根拠:** [plan:104](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:104) の「VERDICT が対照より悪化しない」では、両方が未収束なら通ります。確認した記録は次のとおりです。

   | run（`case/53.c3x_vane_cht/` 配下） | 既存 VERDICT | `rms_roe` 最終代表値 |
   |---|---|---:|
   | [run_0122_fxhalf](/home/sano/work/forge/case/53.c3x_vane_cht/run_0122_fxhalf/CONVERGENCE_VERDICT.txt:2) | `NOT CONVERGED (stalled/plateau)` | 4.22 |
   | [run_0123_fxgeom_ctrl](/home/sano/work/forge/case/53.c3x_vane_cht/run_0123_fxgeom_ctrl/CONVERGENCE_VERDICT.txt:2) | `NOT CONVERGED (stalled/plateau)` | 4.27 |
   | [run_0124_fxhalf_rro](/home/sano/work/forge/case/53.c3x_vane_cht/run_0124_fxhalf_rro/CONVERGENCE_VERDICT.txt:2) | `NOT CONVERGED (stalled/plateau)` | 4.19 |

   これらは「後縁の物理的渦放出」が原因である証拠にはなりません。対象設定は `unsteady: 0` です。また、確認対象の三 run では準定常 VERDICT ファイルを確認できませんでした。

   さらに、`q_eff` の定義を同時に変更すると、「うねりが半減した」の原因が `fx` 変更か残差補正か分離できません。うねり RMS のトレンド除去・空間重み・時間窓も未定義です。

   **対案:** 実装前に次を §6 の合格条件として固定すること。

   - `fx` の A/B は**同じ熱流束定義**で比較し、`raw` と補正後を両方保存する。
   - 定常回帰には `check_convergence.py` の `PASS` を要求する。収束床参照は、参照 run 自身の妥当性を確認して使う。
   - 壁熱流束 RMS・領域別熱伝達率・積分熱量の時系列を作り、[check_quasisteady.py:519](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:519) の `--series-csv` で判定する。既定の量だけでは壁熱流束を判定できません。
   - C3X が未収束のままなら機構診断として扱う。定常値や物理的な時間平均として生産報告を更新しない。
   - メッシュ品質、restart 元、判定区間、使用バイナリ・設定を残作業表に追加し、生産値の更新は結果レビュー後に置く。

**推奨は、`fx=0.5` を固定する方針を維持し、上記を直してから実装することです。** 回転非不変の欠陥はコードから確認でき、同じ幾何を 45° 回転させる代数チェックでも現行重みが `0.8 → 0.636608` に変わりました。辺中点評価は [SU2 の離散化方針](https://su2code.github.io/docs_v7/Theory/) とも整合します。

既に [accepted の median-dual plan:271](/home/sano/work/forge/plans/accepted/discretization-median-dual.md:271) に同じ対処の履歴があり、今回は撤去時の誤った前提を訂正する位置付けです。cell の変更を別扱いにするスコープも妥当です。

実装前の優先順は **① `H_w` 蓄積説の撤回と `e_w` 補正の適用条件明記 → ② 参照箇所・回帰範囲の補完 → ③ 精度主張と格子検証の修正 → ④ 収束・準定常・比較量の合格条件固定**。レビューは read-only のため **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 0
