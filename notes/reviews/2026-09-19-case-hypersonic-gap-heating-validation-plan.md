# codex レビュー: case-hypersonic-gap-heating-validation (plan)

- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `3d4b0d25` (feature/sern-design)
- **codex**: effort `high`, 7.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M8/m1
- **focus**: §4.2 条件再構成 (燃焼生成物試験気体・3点同時一致ゲート)、§4.3 分母の3系統、§4.6 すきま内定常化とIC独立性、§6 の合否ライン (帯の妥当性)、そして T1 層流→case/49 の SST 3D への外挿の論理に重点を置いてほしい
- **extra**: `notes/investigations/hypersonic-gap-cavity-heating-survey.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

T0→T1 の層流検証には価値がありますが、現計画では比較量の定義と合否基準が不十分です。  
特に「T1/T1b の差＝case/49 のモデル形式不確かさ」は成立せず、T2 の分母・相関選択にも修正が必要です。  
ファイルは変更していません。以下の提案は **plan 未反映**です。

`plans/README.md` と `plans/accepted/` を確認した範囲では、Wieting の深キャビティ検証と重複する完了計画はありません。一方、流入 BL 前駆は case/49 に既存実装があります。平面 node・SLAU・block-DPLUR・段階起動という基本構成は現行方針と整合します。今回の主な問題は検証の組み立て方です。

1. **Critical — 層流→SST の差を、case/49 のモデル形式不確かさとして移植できない**

   **根拠:** [plan:87](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:87)、[plan:213](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:213) はこの差を「実測値」としています。しかし T1b は、層流実験に乱流モデルを適用した感度試験です。実際の乱流キャビティに対する SST の誤差を測っていません。

   同一メッシュでも、SST を入れると前縁から開口までの BL 自体が変わります。差には流入 BL の変化とキャビティ内の乱流輸送が混在します。さらに case/49 には環状流路、周方向輸送、偏心、異なる壁温履歴があります。T1/T1b の差が小さくても、これらの誤差を拘束できません。

   **対案:** T1b は「層流条件での SST 適用感度」としてのみ報告する。case/49 の定量的不確かさは、適用条件を合わせた乱流実験検証と case/49 自身の数値・モデル感度から評価する。それまでは **モデル誤差は未定量**と明記する。

2. **Major — キャビティ位置の再構成で、圧力孔位置と開口位置を混同する余地がある**

   **根拠:** [plan:132](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:132) は前縁から 3.5 in の圧力孔を位置決定の根拠に挙げています。[W70 の Fig.1、Table III](/home/sano/work/forge/papers/gap_heating/Wieting_1970_NASA-TN-D-5908_deep_cavities_hypersonic.pdf) を画像で確認すると、3.5 in は模型面圧の測定位置です。後壁位置は前縁から **6.24 in**、したがって分離位置は \(L=6.24\mathrm{\,in}-w\)、中点は \(6.24\mathrm{\,in}-w/2\) です。

   最狭幅でも中点は約158 mmです。これを88.9 mmに置けば、他条件一定の平板則 \(q_{fp}\propto x^{-1/2}\) だけで約33%変わり、G1/G2 の許容幅を超えます。Fig.1 の公称前縁半径も **0.03 in＝0.762 mm**で、現在の感度点に含まれていません。

   **対案:** 前縁・前壁・後壁・中点・圧力孔を別々の座標として固定する。前縁は図面の公称形状を基準とし、鋭前縁と試験後の鈍化を感度条件にする。条件再構成より先に幾何表を確定する。

3. **Major — 「3点同時一致」は条件の一意性や試験気体の再現を保証しない**

   **根拠:** [plan:99](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:99) では \(\gamma\) を未知としていますが、TP では組成と温度から決まります。[`GasSemiPerfect`:194](/home/sano/work/forge/design/forge_design/gas/semiperfect.py:194) も、エンタルピー差と温度依存比熱から速度・音速を求めています。

   組成 \(Y\) と輸送則を固定すれば、再構成は例えば

   \[
   h(T_\infty,Y)+\frac{M^2\gamma(T_\infty,Y)R(Y)T_\infty}{2}
   =h(T_t,Y),\qquad
   \rho_\infty=\frac{Re'_\infty\mu(T_\infty,Y)}{U_\infty}
   \]

   です。ここで合わせた \(M,Re'_\infty\) は入力の再現確認であり、独立した検証点ではありません。文献の \(q_{fp}\) も解析値なので、G1 は主に条件・物性・基準式の整合確認です。

   また、CEA に必要な燃料空気比、反応物温度、凍結条件の決め方がありません。\(p_t\) を拘束から外すこと自体は可能ですが、組成誤差を未特定の全圧損失へ吸収してはいけません。

   **対案:** W70 が物性評価に用いた文献18を基に燃焼生成物モデルを再構成し、組成・輸送則・低温域の扱いを先に固定する。\(M,Re',T_t\) の再現と \(q_{fp}\) の照合を別ゲートにし、再構成した全圧と実測全圧の差も診断量として残す。空気 TP は代替モデルの感度として扱う。

4. **Major — 分母3系統を列挙しても、比較時の共通分母と平均面積が未定義**

   **根拠:** [plan:117](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:117) は実験と CFD 内部比に異なる分母を割り当てています。異なる分母の比を直接比較すると、smooth 側の誤差が相殺されます。

   さらに [plan:183](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:183) は平均量を「面積分」とするだけです。W70 本文 p.10 は明確に

   \[
   \bar q_c=\frac{Q_c}{w\ell},\qquad
   \frac{\bar q_c}{q_{fp}}=\frac{Q_c}{Q_{fp}}
   \]

   と定義しています。**平均の分母は開口面積です。** 濡れ面積 \((2d+w)\ell\) で割ると、文献の \(1.1\to0.5\) は約 \(0.0336\to0.1038\) となり、幅に対する傾向まで逆になります。G6 と G7 も独立した検証量ではありません。

   **対案:** 一次比較は実験・CFD とも Table IV の同じ \(q_{fp}\) を使い、次元付き \(q_s,Q_c\) も併記する。再計算した解析値と CFD smooth 値は分母診断として別表示する。開口面積平均と濡れ面積平均を別名で定義し、前壁理論補完版にはそのモデル依存性も付記する。

5. **Major — T2 の \(h_{smooth}\) が定義できず、比較相関も対象と一致していない**

   **根拠:** [plan:140](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:140) の平板は断熱です。実装も [`viscousFlux_d.cu`:523](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:523) で断熱熱流束を厳密にゼロにしています。この smooth run から通常の \(q/(T_{aw}-T_w)\) で有限の熱伝達係数は得られません。

   [CK75](/home/sano/work/forge/papers/gap_heating/Christensen_Kipp_1975_NASA-CR-134345_RSI_gap_heating_analysis_vol1.pdf) の引用式 \(0.366(Z/W)^{-0.96}\) は**層流・in-line gap**用です。2D 横すきまの乱流比較に対応するのは §5.6／Fig.191 です。引用式は \(Z/W=4\) で9.67%、10で4.01%となり、TH75 の2%を同時に普遍的上限として課す設計とも整合しません。

   [TH75 本文 p.15](/home/sano/work/forge/papers/gap_heating/Throckmorton_1975_NASA-TN-D-7939_RSI_tile_array_gaps_M10.3.pdf) の2%は、特定の厚い BL・配列・中心線外の観測です。本文には低熱流束点の精度への留保もあります。M5 の単一すきまへの一般的合格条件ではありません。

   **対案:** 熱的履歴を明記した有限温度差の smooth 基準を別途用意する。乱流横すきま用相関について、式番号、分母、必要な BL 量、適用範囲、残差帯を固定する。条件を合わせない TH75 比較は参考比較に限定し、G12 の合格条件から外す。

6. **Major — TP・燃焼生成物の検証なのに、熱伝導モデルと抽出器が閉じていない**

   **根拠:** [plan:148](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:148) は粘性則のみ指定しています。[`gasProperties_d.cu`:58](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:58) では `viscMethod: 1` の熱伝導率は別設定で決まり、[`solver-settings.md`:353](/home/sano/work/forge/procedures/solver-settings.md:353) によれば `thermCondMethod` の既定は定数です。組成 DB を替えるだけでは燃焼生成物の輸送物性を替えたことになりません。

   土台とする [`cooled_plate_eval.py`:25](/home/sano/work/forge/case/48.flat_plate_cooled_m4/tools/cooled_plate_eval.py:25) は `CP=1004.5, PR=0.72` 固定で、[105行](/home/sano/work/forge/case/48.flat_plate_cooled_m4/tools/cooled_plate_eval.py:105) の熱流束評価もそれを使います。

   **対案:** `conditions.json` から EOS、\(\mu(T,Y)\)、\(\lambda(T,Y)\)、Pr、エンタルピー基準を一貫して解決する。抽出器も同じ物性を読む。空気／燃焼生成物比較では「熱力学のみ変更」と「輸送物性も変更」を区別し、T0 で両経路を検査してから T1 に進む。

7. **Major — G4/G6 の±25%は、§4.9 の不確かさ予算から導けない**

   **根拠:** [plan:209](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:209) は実験＋digitize の帯を合否基準としていますが、[plan:260](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:260) は±25%です。

   代表 \(q_{fp}=30.4\) kW/m²、digitize を局所値の3%と解釈すると、単純加算でも比の許容幅は
   \[
   \Delta R=0.68/30.4+0.03R
   \]
   です。\(R=0.6\) で相対6.7%、\(R=1.5\) で4.5%。±25%とは一致しません。一方 \(R=0.02\) では実験誤差だけで局所値の約112%なので、深部を高精度の相対誤差で評価できません。

   **対案:** 工学的許容誤差と実験との整合判定を分ける。測定点ごとの絶対誤差、読取誤差、共有分母の相関、前壁理論補完の不確かさを伝播する。深部は検出限界を含む絶対量で判定する。±25%を残すなら「実験不確かさ帯」とは呼ばず、用途から根拠を示す。

8. **Major — G9/G10 は、深部の定常化と IC 独立性を保証するには不足**

   **根拠:** [plan:167](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:167) の指定量は、現行 [`check_quasisteady.py`:377](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:377) の組込み量にありません。`cavity_eval.py` の時系列 CSV と `--series-csv` を接続する実装が必要です。`--segment` に必要な `stage_manifest.json` 生成も実装ステップにありません。

   また、現行判定器は漸近値を合否に使いません。[258行以降](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:258) を実行して確認した結果、既定閾値では合成系列 \(v(s)=2-e^{-s/1000}\)、\(s=0,\ldots,11\) に対して、

   ```text
   STEADY
   drift=0.4%/tail
   漸近値 2 (最終比 +97.836%)
   ```

   となりました。**漸近値を併記するだけでは、未緩和の場を合格にできます。**

   既存データも確認しました。`case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/` の本段単独判定は `NOT CONVERGED`、`rms_roe` は最終0.142・低下0.1桁です。`case/49.plate_annular_cavity_m5/run_0003_stageB_cpg/cavity_series.csv` の壁別熱量は各1/13点が非有限で `NONFINITE`。後者は系列欠損の証拠であり、ソルバ発散の証拠ではありません。前者も同一設定区間全体の不合格を断定するものではありませんが、無条件に合格済み基準として継承できません。

   **対案:** 判定区間、必須CSV列、出力間隔、窓延長試験を明記する。IC 差には相対＋絶対許容を設定し、数値誤差をその許容より小さくする。二つの IC は共通の暖機で差が消えないよう、同じ最終方程式の開始時点で与える。代表ケースでは非定常摂動と時間刻み・内部反復感度も確認し、定常反復への収束だけで物理的安定性を認定しない。

9. **Major — 総入熱の GCI だけでは、局所分布・深部熱流束・抽出精度を検証できない**

   **根拠:** [plan:226](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:226) と G8 は代表1幅の総入熱だけです。開口付近が支配する積分量が収束しても、深部熱流束の2%精度は保証されません。

   W70 の最上部熱電対は Table I では **\(x=0.01\) in＝0.254 mm**です。G4 の \(x/d=0\) は実測点そのものではなく、本文・近似曲線の上端値です。これを CFD の鋭角ノード最大値と比較すると、観測位置と格子依存性が混ざります。

   加えて [`nodeWallDirichlet_d.cu`:85](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:85) は等温壁ノードのエネルギー残差をゼロ化します。場から求めた熱流束と離散的な熱授受の一致は、残差ゲートでは確認できません。

   **対案:** 実測座標での局所値、固定した深部領域の熱量、総入熱それぞれに格子誤差を出す。少なくとも最狭幅と広い幅を対象にする。物理熱流束と独立した離散エネルギー収支を照合し、開口の対流・伝導・粘性仕事、非定常時の蓄積を含める。低マッハ A/B も方式番号を固定する。`lowMachPrecond: 2` は [実装上 RHS 散逸も変更](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1433)するため、単なる収束加速の比較ではありません。

10. **Minor — 「上端3wに深さ方向節点の半分」は2幅で成立しない**

    **根拠:** [plan:136](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:136)。幅7.77／10.64 mmでは \(3w=23.31／31.92\) mmとなり、全深さ20.3 mmを超えます。

    **対案:** 開口近傍の細分化高さを深さで制限し、残りの領域、最大成長率、底隅の解像を含む幅別のメッシュ規則にする。

**推奨は、計画を「W70 の T0/T1 一次検証」に絞って再構成することです。**

実装前の優先順は、①幾何・試験条件・平均面積・共通分母の確定、②物性と抽出器の整合、③不確かさと局所格子・定常性・収支ゲートの確定です。その後に T0→T1 を実施する。T1b/T2 は探索的感度として分離し、乱流横すきまの検証が成立するまで case/49 の誤差幅として受け渡さないことを勧めます。

指摘数: Critical 1 / Major 8 / Minor 1
