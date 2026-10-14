# codex レビュー: tooling-nozzle-bell-two-mode (plan)

- **plan**: [`plans/active/tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md)
- **stage**: `plan`
- **date**: 2026-10-11
- **commit**: `75d7eb8e` (docs/thruster-nozzle-design)
- **codex**: effort `high`, 5.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M8/m1
- **focus**: §4.1 推力・評価器の TP 化、§4.3 (D) の拘束の成立、§4.4 (M) の TP Rao 導出と C2 接続 g(κ_B) の根の存在、§4.5/§6 V3 の比較の公平さと U_k、§8-1 の寸法の仮値。前回 diagnose (notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) の指摘が反映されたかも確認。*.pdf・*.log・plans/README.md は読まない

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**
目的と「TP 評価器 → D → M」の順序は妥当です。ただし、現案には推力・収束・不確かさを誤判定する穴があり、そのまま比較キャンペーンへ進めません。
推奨は、共通評価器の受理を先行させ、M は導出・C2 成立性の独立ゲートを通してから組み込むことです。

指定どおり、ファイル変更・forge 起動・PDF・ログ・`plans/README.md` の閲覧はしていません。確認した `accepted/` では、[Phase 2](/home/sano/work/forge-thruster/plans/accepted/tooling-nozzle-moo-loop.md:17) は TOP の直接最適化であり、今回の TP Rao 生成器は未解決の課題です。一方、該当する `run_0074`〜`run_0076` の場・履歴はこの環境に存在せず、過去の「0.133 %」と収束は再認証できません。以下の数値検算は、明記した合成例・幾何計算です。

1. **Major — TP 化の変更範囲が、実際の評価経路を覆っていない。**

   **根拠:** [plan:139](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:139) の変更対象には、推力抽出器と最適化ドライバがありません。しかし、[extract.py:38](/home/sano/work/forge-thruster/design/forge_design/metrics/extract.py:38) は圧力を CPG の `P=(γ−1)(ρE−ρu²/2)` で復元し、同ファイルの `cf_ideal` も CPG です。[driver.py:84](/home/sano/work/forge-thruster/design/forge_design/opt/driver.py:84) の流量ゲート分母も定比熱式です。

   さらに、[runner.py:197](/home/sano/work/forge-thruster/design/forge_design/evaluate/runner.py:197) は既存 IC 関数に TP 引数を渡さず、[同:243](/home/sano/work/forge-thruster/design/forge_design/evaluate/runner.py:243) の壁平滑化も CPG の保存量を再構成します。`FrozenGas` はそのままでは MOC が要求する `nu`・`mach_of_nu` インターフェースも持ちません。

   **対案:** §5・§7 に、推力抽出、固定 ε の TP 理想推力、TP チョーク流量、IC・壁平滑化、物性記録、最適化ドライバを明記してください。同一の解決済み組成・物性 DB・datum を `GasSemiPerfect` と CFD に渡す契約にします。CEA 照合は出口温度だけでなく、ṁ₁D・出口速度・圧力・F₁D まで対象にし、反応物の温度・状態と燃焼室での凍結条件も固定すべきです。

2. **Major — 「閉じた CV」の推力式に出口の粘性応力が欠けている。**

   **根拠:** [plan:52](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:52) は入口に粘性応力を含め、出口には含めていません。SST の粘性流れで、これを無条件に完全な運動量収支とは呼べません。また、[sern_forces.py:73](/home/sano/work/forge-thruster/design/forge_design/metrics/sern_forces.py:73) は `twall` の扱いに node/cell の符号差を明記しています。既存平面 SERN の会計をそのまま軸対称へコピーするのも危険です。

   **対案:** 外向き法線を用いて全境界を

   `Πₓ = ∫[ρuₓ(u·n) + (p−pₐ)nₓ − (τ_eff n)ₓ] dA`

   に統一し、定常時の総和ゼロから両会計を導いてください。軸対称の `dA=2πr ds`、力の向き、CV 残差の正規化分母を固定します。出口応力を省略する場合は、その寄与が誤差予算内であることを別途確認する必要があります。前回 **C1 は入口圧力の追加まで反映済みですが、完全解決ではありません。**

3. **Major — 現行のゲートを継承すると、「未収束」と「推力未確認」を受理できる。**

   **根拠:** [driver.py:173](/home/sano/work/forge-thruster/design/forge_design/opt/driver.py:173) は `conv != DIVERGED` を合格条件にしており、`STALLED` も通ります。準定常判定は引数なしで呼ばれますが、[check_quasisteady.py:377](/home/sano/work/forge-thruster/solver_density_cuda/tools/check_quasisteady.py:377) の既定対象は `shock, asym, machmax, pmax` で、推力・η・流量ではありません。

   [polish.py:41](/home/sano/work/forge-thruster/design/forge_design/opt/polish.py:41) は `runner.prepare` 経由で再補間・壁平滑化を行い、η の前後差だけで終了します。これは同一メッシュの純粋な継続計算ではなく、収束・準定常の再判定もありません。

   **対案:** 受理条件を、同一設定区間の `check_convergence` の **PASS** と、両推力会計・η・ṁ・CV 残差の時系列に対する **STEADY** にしてください。後者は既存の `--series-csv` を使えます。polish は `restart_field.py` による同一メッシュ継続へ変更します。

   §5 には、最初の TP seed の段階起動、形状変更時の cross-mesh restart、`StageManifest`、メッシュ品質・局所 y₁⁺ の受理を追加してください。特に真空評価では、[runner.py:126](/home/sano/work/forge-thruster/design/forge_design/evaluate/runner.py:126) の `outlet_statPress: Ps=p_ambient` を残さず、計画の超音速流出を実際の `outflow` として配線する必要があります。

4. **Major — V3 の U_k は誤差幅として成立せず、「同等」の根拠にできない。**

   **根拠:** [plan:125](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:125) の三つの項に問題があります。

   - **格子:** 各モードの格子差の最大値は、差 Δη の誤差幅ではありません。合成例で D の格子変化が −4×10⁻⁴、M が +4×10⁻⁴ なら、採用値は 4×10⁻⁴ ですが、Δη の変化は **8×10⁻⁴** です。二水準だけでは残存格子誤差も評価できません。
   - **反復:** 連続二点の差が 3×10⁻⁴ 未満でも、漸近値までの残差がその値以下とは限りません。
   - **探索:** 最良と二番目が近いことは、探索されていない良い領域がないことを保証しません。同じ局所解付近を密に評価するだけで、この項を任意に小さくできます。

   **対案:** 各 L_k の採用ペアを系統的な三水準で評価し、**Δη 自体の格子収束**を確認してください。反復幅は延長履歴と漸近推定から、探索の再現性は独立初期点・seed と予算追加から評価します。探索の上限を保証できない場合、結論は「登録した予算で得た設計同士の比較」に限定します。1.5×10⁻³ の工学的同等幅は維持できますが、現在の U_k で分類してはいけません。前回 **B4 は未解決**です。

5. **Major — TP Rao の停留性、最大性、C2 根の存在・一意性が混同されている。**

   **根拠:** [plan:77](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:77) は質量流量・長さの拘束から導出するとし、検証では ε も固定します。固定出口半径が端点変分をどう制限するかが未定義です。[同:80](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:80) の「長さごとに最大推力輪郭が1本」は、未導出の TP・可変 kernel・C2 条件については根拠がありません。

   符号付き g の追加は前回 B1 の適切な修正ですが、B・D の解の枝が切り替われば、符号を挟んでも連続な根探索になりません。また、一次変分ゼロは極小・鞍点でも成立します。

   **対案:** S1 を独立した受理ゲートにしてください。固定 ε を含む変分の許容範囲、kernel の族、κ_B の探索域、B・D の枝の追跡方法を定義し、各 L_k で根と分解能依存性を確認します。複数根の選択規則も必要です。停留性試験には摂動幅の掃引・MOC 細分化と二次変分の符号確認を追加してください。証明範囲が限定される間は「指定 kernel・拘束下の停留輪郭」と表現します。

   CFD ピンの流量にも注意が必要です。[現在仕様](/home/sano/work/forge-thruster/methods/design/overview.md:337) の m* は CFD 保存量の積分ではなく、等エントロピー写像後の流量です。V1 の「スロート流量」の定義を固定し、元の CFD 流量との差を別に検査してください。

6. **Major — δ* 補正後の拘束に、共通の物理スロート半径・位置が抜けている。**

   **根拠:** [plan:79](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:79) と [同:110](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:110) は ε・L・曲率を拘束・検査しますが、物理壁の `r(0)=r_t`、`r′(0)=0` と最小半径を明示していません。曲率を固定しても半径は固定されません。

   この問題は既存の [物理スロート寸法 plan:91](/home/sano/work/forge-thruster/plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md:91) が、実際の物理壁の最小半径を使って扱っている問題と共通です。

   **対案:** 最終壁に `r(0)=r_t`、`r′(0)=0`、`r″(0)=1/R_t`、`r(L)=r_t√ε` と全域の最小半径条件を課してください。共通収縮部との C2 接続も検査対象です。さらに、拘束への射影で δ* 補正をどれだけ戻したかを、位置・壁角・曲率の残差として保存し、受理値を定めます。前回 **B3 は修正方向は正しいものの、寸法保証が未完**です。

7. **Major — δ* 積分法は μ(T) の置換だけでは TP 化されない。**

   **根拠:** [deltastar_integral.py:208](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:208) の運動量厚さ方程式は、速度勾配を `1/[M(1+(γ−1)M²/2)]` で換算しています。[同:116](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:116) の上流勾配も同じ定比熱関係です。局所 γ(T) を代入しても TP の関係にはなりません。

   既存 `FrozenGas` を使った**単一 H2O の合成検算**では、Tt=3500 K、T=1000 K における `d ln u/dM` はエンタルピーからの数値微分で **0.0510884**、現行換算式で **0.0480304**、差は **−5.99 %** でした。これは基準ケースの CFD 結果ではありませんが、μ 以外にも誤差源がある直接の反例です。

   **対案:** 方程式を TP の縁速度から得た `d ln u_e/dx` で評価してください。回復温度もエンタルピーから求める形に揃え、残す閉包近似を明記します。k_f は較正用形状と検証用形状を分けて受理し、単一倍率で形状依存の熱力学誤差まで吸収したと扱わないことです。

8. **Major — D の数学的構成は可能だが、基底と寸法 spec が未確定で設計族を再現できない。**

   **根拠:** [plan:65](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:65) の拘束は七つです。固定した spline 空間で拘束行列が満行ランクなら、目的関数の零空間である二次多項式はスロートの三拘束で除かれるため、一意解を構成できます。したがって、この方式自体を「過拘束」と退ける理由はありません。ただし、次数だけでは係数数・ノット・重複度が決まりません。

   また、[§8-1](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:144) には収縮部長さと入口直管長がありません。[UpstreamThroatPoly:99](/home/sano/work/forge-thruster/design/forge_design/geometry/wall_walldriven.py:99) の単調条件は `L_U²/[R_t(r_U−r_t)]≤20` です。収縮比3では、許容 L_U/r_t は R_t/r_t=1 で **3.826以下**、0.4で **2.420以下**。現行 runner の既定 `L_U/r_t=3` を使うと、後者は実際に `validate()` が拒否しました。

   **対案:** ノット規約・係数数・設計変数範囲・曲率上限・区間内の極値検査を固定し、無次元座標を明示してください。L_U・L_pipe も両モード共通 spec に追加します。

   r_t=0.14 m 自体に桁の破綻はありません。計算上、出口半径は **0.98995 m**、登録した五つの長さは **1.903／2.220／2.538／2.855／3.172 m**、1000 kN に必要な C_F は **1.624** です。ただし、これは推力達成の検証ではなく仮寸法です。

9. **Minor — 旧 η との 3×10⁻⁴ 回帰試験が、配線変更とソルバ変更を分離していない。**

   **根拠:** [plan:105](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:105) は CPG・`viscMethod:1` を維持すれば旧結果を再現できる構成ですが、[recommended-settings.md:348](/home/sano/work/forge-thruster/procedures/recommended-settings.md:348) には `slauWallNormalChi`、`scalarGradient`、SST 既定値の変更が明記されています。

   **対案:** 配線回帰は同一バイナリ・同一入力場・同一実効設定で新旧評価コードを比較してください。過去 run の再現試験は別に扱い、必要な成果物がない間は未実施とします。

推奨は一つです。**共通評価器を先に受理し、D を実装、その後に独立ゲートを通した M を接続する段階方式**にしてください。実装前の修正優先順は、① TP 評価経路・推力会計、②収束と準定常ゲート・U_k、③ D と寸法 spec、④ M の導出・根の成立条件、⑤ δ* と最終物理壁の保証です。M の成立前に大規模 MOO を投入する費用対効果はありません。

本レビューは **plan 未反映**です。反映先は `tooling-nozzle-bell-two-mode.md` §4.1〜4.5、§5・§5.1、§6、§7、§8-1 です。

指摘数: Critical 0 / Major 8 / Minor 1
