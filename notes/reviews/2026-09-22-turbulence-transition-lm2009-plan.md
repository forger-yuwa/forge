# codex レビュー: turbulence-transition-lm2009 (plan)

- **plan**: [`plans/active/turbulence-transition-lm2009.md`](../../plans/active/turbulence-transition-lm2009.md)
- **stage**: `plan`
- **date**: 2026-09-22
- **commit**: `ded72e17` (feature/sern-design)
- **codex**: effort `high`, 7.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M8/m1
- **focus**: 冷却翼 CHT 検証の最大の誤差 (低 Re SST が前縁から乱流) に対して γ–Reθt (Langtry–Menter 2009) を SST に足す計画。問い: (1) モデル選定 (§3 の比較表) は妥当か。1 方程式 γ (Menter 2015) の方が良い理由があれば指摘してほしい。(2) §4 の設計に、node-centred median-dual + float32 + segregated point-implicit という forge の条件で破綻しそうな箇所が無いか (壁での扱い §4.4、√(γ F_onset) の陰的対角、Re_v の y^2、1 次風上)。(3) §6 の検証計画と合格ラインは妥当か。T3A を M≈0.2 に上げて Re を保つやり方の落とし穴、自由流 Tu の減衰の合わせ方。式の正本は methods/turbulence/theory.md の §11 (grep -n '^## 11' で位置を見て、そこから末尾まで sed で読むこと)。参照実装は .external/su2-src/SU2_CFD/include/numerics/turbulent/transition/trans_sources.hpp (sed -n '75,250p') と SU2_CFD/src/solvers/CTransLMSolver.cpp (sed -n '196,262p')。forge の土台は methods/turbulence/implementation.md の §3.3, §4, §6 (sed -n '58,84p' '371,400p' '489,504p')。禁止: *.log, residual_history.csv, res_*.h5, *.vtu, plans/README.md を読まない。推奨は 1 つに絞ること。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

LM2009 の採用は妥当です。ただし、現計画には拡散係数・SST 結合・周期処理・収束ゲートの具体的な抜けがあります。  
モデルを変更するより、これらを実装前に修正し、T3A の条件と検証基準を固定することを推奨します。

ファイルは変更していません。禁止ファイルは読まず、コード、保存済み `CONVERGENCE_VERDICT.txt`、README、一次資料を確認しました。数値検査はメモリ内で実施しています。`plans/README.md` は指定どおり未確認ですが、`plans/accepted/` に同等の遷移モデル実装は見当たりません。CHT plan の残作業 #29 を具体化する計画として位置づけるのが適切です。

1. **Major — descriptor を二本追加するだけでは、`Reθt` の拡散係数が間違う。**

   **根拠:** 計画 §4.1 の `2(μ + μt)` は正しい式です。しかし、現行 [scalarTransport_d.cuh:28](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cuh:28) の descriptor が表せるのは `μ + σ μt` だけです。[scalarTransport_d.cu:216](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:216) でも分子粘性の係数は固定で 1。`sigma=2` を指定すると `μ + 2μt` になり、層流域では必要な拡散の半分になります。

   **対案:** 分子粘性にも独立の係数を設け、`γ` は `(1,1)`、`Reθt` は `(2,2)` にすることを §4・§5 に明記してください。単一版・融合版の流束と `transport_diag` の両方に反映し、`μt=0` の拡散試験で係数 2 を検証する必要があります。

2. **Major — `F1` と `Pω` の結合箇所を指定しないと、LM2009 と異なる方程式になる。**

   **根拠:** [ransSource_d.cu:143](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:143) は `F1` を再計算し、`sstF1` を上書きしたうえで `α`・`β`・交差拡散に使います。計画どおり `ransBlendF1` だけを変更すると、拡散係数には補正後 `F1`、ソースには補正前 `F1` が使われます。

   また、[ransSource_d.cu:205](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:205) の `Pw` は可変の `Pk` を直接参照します。先に `Pk *= gammaEff` を挿入すると、禁止したはずの `ω` 生産抑制も入ります。SU2 は [turb_sources.hpp:971](/home/sano/work/forge/.external/su2-src/SU2_CFD/include/numerics/turbulent/turb_sources.hpp:971) で `pw` を確定した**後**、1003 行で `pk` を補正しています。

   **対案:** `Pk_base` と `Pk_transition` を分け、`Pw` は前者だけを参照してください。`F1` は共通関数にして、拡散・ソースの全利用箇所を同じ値にそろえるべきです。さらに SU2 の版・相関・SST オプション・Kato–Launder・圧縮性補正・一次風上の対象を固定し、「同じモデル」を設定表で定義してください。既定設定の継承では比較条件を保証できません。

3. **Major — 周期境界は残差 gather だけでは足りない。**

   **根拠:** 計画 §4.1・§5 は残差 gather しか明記していません。しかし既存実装には、別途、[periodicNode_d.cu:94](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:94) の輸送対角集約と、[同:137](/home/sano/work/forge/solver_density_cuda/cuda_forge/periodicNode_d.cu:137) の更新後状態ミラーがあります。ソースには [ransSource_d.cu:357](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:357) のとおり**部分体積**を渡さないと、gather 後に seam で二重計上します。

   **対案:** 二変数について、部分体積によるソース組立て、残差・輸送対角の gather、合併体積による更新、状態ミラー、原始量の再生成までを仕様化してください。体積当たりのソース Jacobian は残差と同じように単純加算しないこと。翼へ進む前に、周期 seam を横切る移流・拡散と一様ソースの試験を追加すべきです。

   node 限定は [検証手順:63](/home/sano/work/forge/procedures/verification/README.md:63) と整合します。cell・壁関数・DES・未対応の軸対称などは、未検証のまま動作させず、受付条件を設定検査で明示してください。

4. **Major — 相関の「最大 10 回」と数値安全性の仕様が不足している。**

   **根拠:** SU2 の [trans_sources.hpp:220](/home/sano/work/forge/.external/su2-src/SU2_CFD/include/numerics/turbulent/transition/trans_sources.hpp:220) は最大 **100 回**です。同じ反復式をメモリ内で評価すると、許容される低 `Tu` の状態で次の反例が得られました。これは CFD run の結果ではありません。

   | 入力・結果 | 値 |
   |---|---:|
   | `Tu` | 0.1% |
   | `ν(dU/ds)/U²` | −7.28618×10⁻⁸ |
   | 10 回後の相関値 | 743.2174 |
   | 相対変化 10⁻⁷ まで反復した値 | 745.5042 |
   | 必要反復数／10 回打切り誤差 | 33 回／約 0.307% |

   したがって、float32 の停止閾値を `1e−5` にすることは、10 回で十分という根拠になりません。また、[methods の式:1141](/home/sano/work/forge/methods/turbulence/theory.md:1141) は `U` による除算を含みます。壁だけを除外しても、内部の停滞点でのゼロ割は残ります。

   **対案:** 上限はまず参照実装と同じ 100 回とし、安全な相対誤差判定、未達件数の診断、未達時の処置を定義してください。`U→0`、`Ω→0`、`γ→0`、相関分岐点、正負の圧力勾配を単体試験に含めます。

   負のソース微分だけを陰化する方針自体は妥当です。ただし更新分母は、現コードの [update_d.cu:330](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:330) に合わせて
   `V/Δτ + V max(−Jsource,0) + transport_diag`
   と明記してください。`γ` 下限は原始量基準で定義し、保存量への換算とクリップ量を検査する必要があります。SU2 には既に [CTransLMSolver.cpp:102](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CTransLMSolver.cpp:102) の `γ≥1e−4` があります。

5. **Major — 現行ゲートでは、新しい二方程式が未収束でも検出できない。**

   **根拠:** [check_convergence.py:35](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:35) と 54 行の読込対象に `rms_roGamma`・`rms_roReth` がありません。メモリ内の模擬入力に `rms_roGamma=NaN`、`rms_roReth=1e30` を入れても、両列は検査対象から落ちました。

   また、[stage_manifest.py:81](/home/sano/work/forge/solver_density_cuda/tools/stage_manifest.py:81) に `turbulence.transition` がありません。実関数で確認すると、`transition: none` と `lm2009` の段キーは同一でした。SST 起動段の大きな残差を、遷移モデルの収束判定に混入できます。

   **対案:** 残差出力・NaN 検査・収束対象列・モデル有効時の必須列検査・段キーを実装範囲に追加してください。「新変数が NaN」「新列欠落」「SST→遷移の段切替」を先に回帰試験にします。§6 の「ツールを通す」だけでは、この穴は塞がりません。

6. **Major — `M≈0.2` 化は妥当だが、`Re/m` と前縁 `Tu` 一点だけでは相似条件を固定できない。**

   **根拠:** 計画 [§4.5:75](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:75) は入口粘性比約 12・`Tu_in=3.3–3.5%` と、高速化を別々に記述しています。既存の圧縮性 T3A ベンチマークは `M=0.2`、`Re/m=2.0×10⁵`、入口距離 0.25 m、粘性比 11.9、**入口 `Tu=5.855%`、前縁 3.3%**です。これらを元の実験・論文条件と混用できません。[TMR の条件表](https://tmbwg.github.io/turbmodels/t3_transition_mainpage.html)

   **対案:** まずこの公開条件一式を採用して実装照合を行い、ERCOFTAC 実験再現を別の条件表で定義してください。独自に速度を倍率 `a` で変更する場合、同じ長さでの相似条件には `ν→aν`、`k→a²k`、`ω→aω` が必要です。熱境界条件・`Pr`・粘性則も固定します。

   自由流は前縁だけでなく、遷移域までの複数位置で `Tu(x)` と粘性比を比較してください。入口条件の調整はその減衰データに対して行い、`Cf` に合わせて再調整しないこと。`M=0.2` は近似的な低マッハ条件なので、実験との ±5% 判定には Mach 感度も含めるべきです。

7. **Major — §6 の合格ラインは、格子誤差とモデル誤差を分離できない。**

   **根拠:** [plan:104](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:104) 以降には、一次風上の格子収束、局所壁解像、遷移長さの基準がありません。同一メッシュで SU2 と開始位置が ±10% でも、両者が同じ程度に数値拡散している可能性があります。実際、TMR の LM2009 検証は一次風上を使いながら **8 格子**で収束傾向を確認しています。[TMR 検証結果](https://tmbwg.github.io/turbmodels/t3_transition_LM.html)

   また、「遷移後の `Cf` が前縁から乱流の SST 解に ±5%」は、同じ `Re_x` で境界層の発達履歴が異なるため、直ちには成立しません。`Cf` 最小位置も、実験側と同じ定義で抽出しなければ開始位置の比較になりません。

   **対案:** 一次風上は初版として維持し、少なくとも三段階の格子と局所 `y₁⁺≤1` を要求してください。提案する数値誤差ゲートは、細二格子間で開始位置変化 5% 未満、指定位置の `Cf` 変化 3% 未満です。その後に SU2 ±10%、実験 ±20% を適用します。±20% は出典を示せない限り「文献上の再現度」でなく、暫定的な工学基準と記すべきです。

   `Cf` 分布・開始位置・遷移長さの時系列をそれぞれ判定し、十分下流の比較区間を事前指定してください。翼への移行前には圧力勾配のある検証も一つ必要です。T3A/B だけでは `λθ` と剥離補正の空間的な動作を十分に検証できません。

8. **Major — 動機は有力だが、「最大誤差の原因を確認済み」という断定と翼の評価設計は強すぎる。**

   **根拠:** 保存済みゲートは次のとおりです。

   | 比較基準 run | 保存済み VERDICT | 末尾残差の例 |
   |---|---|---|
   | `case/53.c3x_vane_cht/run_0125_prod_fxhalf/` | `NOT CONVERGED (stalled/plateau)` | `rms_roe=2.76`, `rms_roOmega=11.8` |
   | `case/54.markii_vane_cht/run_0024_tecut_fxhalf/` | `NOT CONVERGED (stalled/plateau)` | `rms_roe=15.1`, `rms_roOmega=2.90e3` |

   出典は各 [C3X VERDICT](/home/sano/work/forge/case/53.c3x_vane_cht/run_0125_prod_fxhalf/CONVERGENCE_VERDICT.txt:2)・[Mark II VERDICT](/home/sano/work/forge/case/54.markii_vane_cht/run_0024_tecut_fxhalf/CONVERGENCE_VERDICT.txt:2) です。README には両者の熱伝達系列について `ALL STEADY` の記録がありますが、**場の収束とは別**です。[C3X README:760](/home/sano/work/forge/case/53.c3x_vane_cht/README.md:760)、[Mark II README:75](/home/sano/work/forge/case/54.markii_vane_cht/README.md:75)

   さらに Mark II の層流対照は `DRIFTING`、同ケースの現行格子系列は壁解像ゲート `FAIL` と記録されています。[README:76](/home/sano/work/forge/case/54.markii_vane_cht/README.md:76)  
   層流と SST が実験を挟むことは遷移導入の有力な動機ですが、誤差の主因を一意に証明しません。

   **対案:** §1 と methods の断定を「最優先で検証する仮説」に変更してください。翼はまず同じ実測壁温・同じ適格メッシュで遷移 ON/OFF を比較し、その後に CHT を比較する順序を明記します。未収束場は機構診断として扱い、必要なら物理時間計算と統計評価へ進めます。

   翼の実験誤差に合格ラインを置かないことは許容できますが、数値品質には必要です。壁解像、派生量の VERDICT、CHT 界面残差・熱収支、内部冷却条件を再同定しないことを必須条件にしてください。

9. **Minor — `y²` と `√γ` のリスク説明が不正確。**

   **根拠:** [plan:77](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:77) の「µm の `y²` は桁落ちしやすい」は、その乗算自体には当たりません。メモリ内 float32 検査で `ρ=1.2, y=10⁻⁶, S=10⁵, μ=1.8×10⁻⁵` とすると、`Re_v=0.00666666683`、相対誤差は約 `2.4×10⁻⁸` でした。

   また、`γ→0` で `Pγ` の微分が発散するのは、`F_onset>0` なら**正側**です。負の部分だけを採る設計では、これがそのまま巨大な陰的減衰対角になるわけではありません。危険なのはゼロでの評価、`0×∞`、陽的に残る成長、更新後の有界性です。

   **対案:** 主なリスクを壁距離の幾何精度、`F_onset2−F_onset3` の差、停滞点の除算、累乗の overflow、有界性に書き換えてください。`Re_v` に恣意的な下限を入れる理由にはしないこと。

壁の扱いについては、**§4.4 の Neumann 方針を支持します**。壁節点をピンせず壁面拡散流束をゼロにすることは、[SU2 の壁処理](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CTransLMSolver.cpp:371) と整合します。ただし、消えるのは物理壁面を通る流束であり、壁 CV と内部を結ぶ双対面の輸送は残します。壁でソースをスキップするときも、ソース対角をゼロ、`gammaEff` を有効な値に毎回設定し、未初期化値を残さない仕様が必要です。

**推奨は、LM2009 を維持して計画を修正することです。** Menter 2015 の一方程式モデルには、輸送式一本の削減、相関反復の除去、Galilean 不変性という明確な利点があります。[原論文](https://link.springer.com/article/10.1007/s10494-015-9622-4)  
それでも、今回の静止翼・SST 資産・手元の SU2 参照実装を考えると、検証可能性を優先する LM2009 が適切です。「一方程式モデルの精度が劣る」ことを不採用理由にはしません。

実装前の修正優先順は次のとおりです。

1. §4 に拡散係数、`F1` 共通化、`Pk_base` 分離、周期の体積・対角・状態処理を確定する。
2. ソースの停止判定・有界性・停滞点処理と、対応する単体試験を定義する。
3. §5 に収束ゲート・段キー・新変数の残差出力を追加する。
4. §6 に固定した T3A 条件、乱流減衰、格子・壁解像・時系列の判定基準を追加する。
5. §1 の因果断定を修正し、翼の固定壁温比較から CHT へ進む条件を定める。

読み取り専用の依頼に従い、以上は **plan 未反映**です。

指摘数: Critical 0 / Major 8 / Minor 1
