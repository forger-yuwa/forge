# codex レビュー: condensation-two-phase-transport (plan)

- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `9909ab7c` (feature/species-transport)
- **codex**: effort `high`, 5.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M6/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
課題の同定は正しいものの、§4.2 の分子拡散モデルには多成分で成立しない点があり、離散更新後の蒸気非負性も保証できません。  
流束・陰的更新・FCT・検証条件を確定してから再レビューすべきです。

レビュー対象は `feature/species-transport`、HEAD `9909ab7c`。ファイルは変更していません。指定された `run_0482`、`run_0483`、case/42 の `run_0071*`、case/45 の `run_0039*` はこのワークツリーに存在せず、`check_convergence.py` も前二者に `NO residual_history.csv` を返しました。したがって、計画記載の実測値は今回再検証できていません。以下の数値反例は、計画式・コードの更新式を Python で評価したもので、forge の実行結果とは区別します。

目的については、既存コードが総水分を拡散させ、液モーメントを拡散させないことを確認しました。[既存の accepted plan](/home/sano/work/forge-species/plans/accepted/species-passive-scalar-unification.md:31) も凝縮モーメントを「拡散なし」と明記しており、本件は解決済み機能との重複ではありません。`h_l = h_v − L` によるエネルギー輸送と、局所定係数での潜熱相殺も [EOS](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:351) と整合します。問題はその先の作用素設計です。

1. **Critical — 分子拡散の駆動勾配が、三成分以上では偽の拡散を作る。**

   **根拠:** [plan:72](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:72) は、混合物基準の気相成分分率に `−ρD_k∇Y_k` を適用し、気相組成で総流束を補正します。しかし、気相組成 z_k が一様で液量 g だけが変化する場合、気相成分の分率は `(1−g)z_k`。計画式の補正後流束は

   `J_k = ρ z_k (D_k − Σ_j z_jD_j) ∇g`

   となり、種別 D が異なるとゼロになりません。**気相組成・T・P が一様でも、液の濃度勾配だけで気相を分離させます。**

   `[0.2, 0.3, 0.5]` の等分子量三成分、二元 D を `[1, 2, 3]×10⁻⁵ m²/s` として、[現行の混合平均式](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/thermo_d.cuh:566) を評価すると、`∇g=0.1 m⁻¹` に対する `J/ρ` は `[-1.3182e−7, −6.1364e−8, 1.9318e−7] m/s`。総和はゼロでも、各成分には偽流束が残ります。二成分・同一 D の試験では検出できません。

   **対案:** 分子拡散は気相基準の組成勾配で定義してください。現行の質量分率 Fick 近似を継承するなら、まず `ρ_g=ρ(1−g)`、`z_v=(Y_w−g)/(1−g)` として `j_k⁰=−ρ_gD_k∇z_k` を作り、`j_k=j_k⁰−z_kΣj_j⁰` とする。係数と駆動変数の対応も明記する必要があります。気相組成を駆動変数とする点は [Cantera の公式拡散式](https://www.cantera.org/stable/reference/onedim/governing-equations.html#diffusive-fluxes) と整合します。上記三成分の「気相一様・液量非一様で分子流束ゼロ」を必須試験にしてください。

2. **Major — 流束差が拡散形でも、別々の陰的更新後には蒸気が負になり得る。**

   **根拠:** [speciesTransport_d.cu:265](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:265) は種別の分子＋乱流拡散を対角に加え、[同:370](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:370) で成分ごとに増分を割ります。液側には分子拡散がないため、同じ乱流残差でも総水分と液の増分は一致しません。

   例えば、相変化・移流・制限を除いた二セル試験で、両セルの蒸気をゼロ、`Y_w=g` を `0.1/0.2`、分子・乱流の面係数と `V/Δτ` を各 1 とします。既存の対角構成を継承すると、流入側は `Y_w'=0.13333334`、`g'=0.15`、**蒸気は float32 で `−0.016666666`**。種の再正規化後も負のままです。

   さらに [main.cpp:2094](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2094) は化学種を再正規化してから液を更新します。連続式の説明だけでは、この更新写像を保証できません。

   **対案:** 保存変数の格納形式は維持しつつ、蒸気・液の連成した更新を設計してください。共通の保存的な制限を総水分・液・対応するエネルギー流束に適用し、`ρY_w−ρg≥0` を更新後に保証する必要があります。陰的対角、種の再正規化、モーメント制限まで検証対象に含めてください。

3. **Major — §4.1 の変更対象では、既定の表引き経路が直らない。**

   **根拠:** [plan:67](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:67) は変更を `gas_transport_cell_Y` と壁組成に限定しています。しかし既定の表引きは [gasProperties_d.cu:100](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/gasProperties_d.cu:100) から **`gas_transport_cell_rY`** を呼びます。[同:27](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/gasProperties_d.cu:27) は `roY` をそのままコピーしており、`gas_transport_cell_Y` の修正は通りません。

   **対案:** 表引き・直接評価・壁・診断 probe が同じ気相組成構築を使う設計に変更してください。試験も表引き ON/OFF、湿潤セル、湿潤壁、表範囲外への退避を含める。直接評価だけの合格では不十分です。

4. **Major — 残差の組立順序と dual-time FCT が、統一拡散の設計から抜けている。**

   **根拠:** [main.cpp:1833](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1833) は化学種輸送の後に凝縮輸送を呼び、[passiveAdvection_d_wrapper:1269](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1269) は液の残差・対角をゼロ初期化します。化学種拡散内から液残差へ加算する実装では、その寄与が消えます。

   また [speciesTransport_d.cu:1719](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1719) の FCT は拡散対象 `qDiff` をトレーサ一つに限定し、[passiveFct_d.cuh:218](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/passiveFct_d.cuh:218) もその成分だけを面流束に含めます。凝縮残差に拡散を足しても、FCT の低次作用素・面流束履歴には反映されません。

   **対案:** 全残差の初期化後、境界ピン・周期集約前に二相拡散を組み立てる順序を確定してください。同じ面流束を液・総水分・エネルギー・FCT 履歴で共有する設計が必要です。CPG／単成分は [speciesEnabled:52](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:52) を通らないため、対応経路または明示的な対象外条件も定義してください。

5. **Major — 有限 `Sc_l` は、質量収支だけ直してもモーメント分布を壊す。**

   **根拠:** [plan:76](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:76) は全相の収支への整合を約束していますが、具体式と `Q0/Q1/Q2` の分子輸送がありません。g だけを拡散すると、単分散で等号が成立していた `Q2²≤Q1Q3` に対し、流出側の Q3 だけが減って直ちに実現可能性を破ります。現行コードは [condensationRealizability_d.cuh:129](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:129) で g から Q3 を作り、違反時に Q1/Q2 を射影します。

   また `μ/Sc_l` なので、**小さい `Sc_l` は強い拡散**です。「微小な液 Sc」を「微小な安定化」と扱えません。

   **対案:** 初回実装は液の分子拡散ゼロに固定してください。有限 `Sc_l` は、気相の補償流束・全液モーメントの輸送・エネルギー・実現可能性を一体で定義する後続項目に分離することを推奨します。

6. **Major — g のクランプ総量だけでは、数値的不整合の監視にならない。**

   **根拠:** [condensationRealizability_d.cuh:148](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:148) の予算は、上限・下限制限に加えて [同:186](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:186) の液滴消滅もまとめて記録します。一方、先行する [増分制限・floor](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:388) や化学種再正規化は別処理です。

   したがって、物理モデル上の消滅で WARN が出る一方、先行制限が違反を隠して最終クランプだけゼロになる可能性があります。

   **対案:** 蒸気上限違反、負値 floor、更新制限、モーメント射影、液滴消滅、化学種再正規化を理由別に計上してください。累積値に加え区間値、制限前の最小蒸気量・違反数を出し、液量ゼロ時の正規化と restart 時の累積の扱いも定義する必要があります。

7. **Major — §6 は実装を受け入れる検証仕様になっていない。**

   **根拠:** [plan:101](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:101) は保存誤差・非負性・物性比較の許容差を未確定にしています。`p95|ΔT|=1 K` は影響の大小を分類する値で、正しい実装の合格条件ではありません。

   また、[準定常ツールの標準量](/home/sano/work/forge-species/solver_density_cuda/tools/check_quasisteady.py:377) には onset・出口 g・壁熱流束がありません。[`--series-csv`](/home/sano/work/forge-species/solver_density_cuda/tools/check_quasisteady.py:519) 用の抽出定義が必要です。旧 `viscMethod:2` は [現行設定では置換済み](/home/sano/work/forge-species/procedures/solver-settings.md:202) なので、旧 run の物性誤差を現行モデルの基準にもできません。

   **対案:** 実装前に、少なくとも次を合否付きで確定してください。

   - 三成分の気相一様試験、蒸気ゼロ近傍試験、相変化を止めた拡散単独の保存試験。
   - 面ごとの `ΣJ=0`、`J_w=J_v+J_l`、エンタルピー輸送と EOS 応答の照合。誤差の分母・絶対許容差・float32 許容差を明記する。
   - node 周期と軸対称の試験。[標準検証表](/home/sano/work/forge-species/procedures/verification/README.md:15) に沿って case/09、case/44 を使う。
   - A/B の物性モデルと初期保存量を固定し、報告量の時系列 CSV、`check_convergence` と `check_quasisteady` の判定条件を明記する。

   cell は現行規則上、実行回帰の対象外です。共有コードの経路確認と「cell 未検証」の明記を行い、node の結果で両方式を保証しないでください。

8. **Minor — 実装順序が本文・残作業表・リポジトリ規則で一致していない。**

   **根拠:** [§5](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:84) は監視を先に入れる一方、§5.1 は試験用拡散実装を要する A/B を先に置き、docs を最後にしています。[AGENTS.md:164](/home/sano/work/forge-species/AGENTS.md:164) は仕様文書と plan を実装前に揃える規則です。

   **対案:** 「仕様・合否確定 → 監視 → 現行バイナリで基準取得 → 物性 → 二相拡散の単体検証 → A/B → 回帰」に統一してください。

**推奨は、§4.2 を“気相内の分子拡散＋全相共通の乱流混合”として再設計し、有限 `Sc_l` を後続へ分離して再レビューすることです。** 優先順は、①分子流束の定義、②保存的で非負性を保つ更新、③組立順序と FCT、④監視・定量ゲート、⑤物性の全経路対応です。現状のまま試験用実装へ進むことも推奨しません。本提案はレビューのみで、**plan 未反映**です。

指摘数: Critical 1 / Major 6 / Minor 1
