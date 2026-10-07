# codex 諮問 (diagnose): throat-monotone-r2

- **brief**: [`notes/reviews/briefs/2026-10-06-throat-monotone-r2.md`](../../notes/reviews/briefs/2026-10-06-throat-monotone-r2.md)
- **plan**: [`plans/active/tooling-nozzle-throat-monotone-r2.md`](../../plans/active/tooling-nozzle-throat-monotone-r2.md)
- **date**: 2026-10-06
- **commit**: `c821b71e` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.3 min, rc=0
- **結論**: **単調拘束の方針は採用するが、Euler 投入の前に、平滑化項の積分方法だけを変える CFD 0 step の形状 A/B を行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **単調拘束の方針は採用するが、Euler 投入の前に、平滑化項の積分方法だけを変える CFD 0 step の形状 A/B を行う。**

採否表（Critical なし。生産化判断には Major の未解決事項あり）:

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| ① 始点 r″ を固定した単調拘束 | **採用** | 始点条件は [wall_axismach.py:739](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:739)。r′(0)=0、r″≤1/R なら積分から r′(x₁)≤x₁/R。したがって第1点の角度誤差は **少なくとも約 0.0219°**。これは当てはめ失敗ではなく、選んだ形状要求の代償。点の重みを下げてもこの下限は消えず、局所 λ 強化では単調性を保証できない。直接の不等式拘束を推奨する。[0,1.5] は今回の固定区間として維持し、境界付近の有効制約と区間外への影響を記録する。 |
| ② r‴ の急変と滑らかさ | **要再検証・Major** | [wall_axismach.py:736](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:736) は全長を8000点で積分するため、スロートでは正則化項の積分が粗い。下記の演算子試験でも大きな重み誤差を確認した。内部ノットが単純な5次スプラインなので、係数拘束だけで内部の C⁴ 連続性は失われない。しかし **C⁴ と「変化が緩やか」は別**。各ノット区間で厳密に求めた ∫(r‴)²dx と max|r⁗| を現行壁と比較する条件を追加し、悪化した場合は滑らかさ改善と呼ばない。 |
| ③ S4・S5 の事前登録 | **独立した検証という扱いは却下・Major** | [対象 plan:112](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:112) の閾値は試し C を見た後のもの。開発用の受入条件にはできるが、C の妥当性を独立に裏付けない。S4 の第1点は上記の解析的下限を根拠にし、その他の誤差は「要求に伴う忠実度の損失」として明示する。旧点上ゲートの不合格を参考値として残した [前 plan:119](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md:119) と同じ扱いが適切。S5 の第1点だけでは、拘束で飽和する量の確認に留まるため、全域の形状・導関数の感度も記録する。 |
| ④ IC の非対称性 | **「18000 step なら無視できる」は却下・Major** | [restart_field.py:6](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/restart_field.py:6) に保存量再構成で差が入った実績がある。A=`restart_field`、B=`interp_field` という使い分け自体は規則どおり。まず B の対応節点・保存量差を記録し、IC 影響は未消去として扱う。各腕3回はこの系統差を消さない。両腕への `interp_field` 適用や座標検査の緩和は推奨しない。なお提案の **1e−6 r_t は 0.0767 µm** で、想定最大変位0.5 µmを許容できない。 |
| ⑤ Δq・U・各腕3回 | **Δq は採用、検出力は要再検証・Major** | 目的は微小な差の検出ではなく、登録した悪化幅を排除することなので、形状差に合わせて Δq を縮める必要はない。3回は再実行ばらつきの確認として妥当だが、`3R` は保証された信頼限界ではない。さらに **出口 M の定義が異なる**。較正側は節点の算術平均 ([euler_grid_ab.py:39](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/euler_grid_ab.py:39))、提案評価器は η 重み付き面積平均 ([eval_wallfit_euler.py:68](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/eval_wallfit_euler.py:68))。両者を同じ「較正済み M=6」と扱わないこと。また後者の出口 M は軸方向評価刻み h に依存せず、現在の E は出口の標本感度を検査できない。比較量・標本・重みを固定し、その評価誤差を U に含める。 |
| ⑥ Euler 非劣化だけで NS 生産化 | **却下・Major** | [対象 plan:139](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:139) の根拠では不足。提示された NS 出口 M=5.998887 の下限5.9988までの余裕は **8.7e−5** しかなく、Euler の許容悪化幅1.8e−4より小さい。δ_r に対する形状変位の比から性能感度は上限評価できない。**候補形状として保持する間は NS 不要、新形状に既存の NS 性能・凝縮結果を引き継いで生産化するなら再評価が必要**。旧壁に対するユーザ決定「下限境界上のまま記録」は、新壁の検証免除ではない ([前 plan:227](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md:227))。 |
| S6 の物理壁検査 | **要修正・Major** | 固定 δ_r なら r_phys″=r_design″+δ_r″ なので、S6 の上限は S1 からほぼ自動的に成立し、物理壁の山の消失を保証しない ([wall_axismach.py:497](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:497))。さらに生産経路は新しい設計壁から `integral_bl` を再計算する ([runner_axismach.py:946](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:946))。固定した run_0117 の δ_r だけでは、生産形状を検査したことにならない。生産経路の物理壁について、r″ の極値・増加区間・全域の形状差を直接検査する。 |
| QP の停止条件と S1 | **要修正・Minor** | 試作は G を行正規化し、`max(1,|Gc|)` を使う ([throat_r2_explainer.py:61](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_r2_explainer.py:61))。plan の停止条件はこれと異なり、S1 の未正規化係数≤1e−12とも対応していない。正規化前後の残差、等式残差、双対実行可能性・相補性を記録し、許容誤差を統一する。 |

第 1 仮説: **試し C の r‴ の急変には、平滑化項の粗い積分が無視できない影響を与えている。** 確度: **中**

  根拠: `joint_fit_wall` と試作の両方が、全長8000点の `Σ(r‴)²·gradient(x)` を使う。評価器記載の x_F=95.22667765 と実装のノット生成式を使い、今回、積分演算子だけを数値確認した。積分刻みは **0.0119048**、最初のノット間隔は **0.0125**。始点の r・r′・r″を変えない基底方向 i=3、4、5で、現行積分／区間ごとの厳密積分はそれぞれ **4.506、0.638、0.846**。これは実際の C の誤差を測った値ではないが、「同じ λ だから十分に平滑化されている」という前提は成立しない。

  反証条件: 同じ MOC 点群・ノット・λ・拘束で積分だけを厳密化し、下記の導関数差が10%以内なら、**急変の主要因が積分誤差である**という仮説を退ける。

第 2 仮説: 急変は r‴≤0 の有効制約から離れる位置で生じる、拘束付き最適解そのものの特徴。確度: 中。試し C の有効制約4本という提示値と整合するが、係数・乗数は未確認。

第 3 仮説: 固定 δ_r で見積もった局所的な設計壁差が、生産の δ_r 再計算を通じて下流の物理壁にも広がる。確度: 低。経路はコードで確認したが、影響量は未確認。

判別 A/B: **単調拘束付き C の平滑化項の積分方法だけを変える。CFD は両腕0 step。**

- A: 現行の全長8000点積分。
- B: 各非零ノット区間で3点 Gauss 積分。r‴ は区間ごとに2次式なので、その二乗の積分を厳密に評価できる。
- MOC 点群、λ=1e−9、ノット、始終点拘束、単調区間、QP 解法は固定する。
- [0,0.3] の Dₖ=max|r_B⁽ᵏ⁾−r_A⁽ᵏ⁾| / max|r_A⁽ᵏ⁾|（k=3,4）、有効制約、全域の Δr・Δθ を測る。極値は区間多項式から求める。

→ **D₃またはD₄>0.10なら積分依存あり**。現在の急変を拘束固有とみなす前提を棄却し、厳密積分で形状ゲートを組み直す。**両方≤0.10なら積分誤差主因説を棄却**し、拘束に伴う形状上の代償として評価する。10%はこの診断の事前閾値であり、性能許容値ではない。

やらない方がよいこと: λ・点重み・拘束区間を同時に調整すること、先にEulerを6本回すこと、粗い壁節点2〜3個から求めた r″だけで山の消失を認定すること、既存 NS 報告を新しい生産 YAML の性能証拠にすること。

呼び出し側の前提への異議: **「山を消す」はユーザ要求、「流れへの影響が小さい」「18000 step でIC差が消える」「NS不要」は未検証の解釈**。また atan(x/R) は円弧の接線角ではなく、r″=1/R一定の放物線の接線角である。ただし今回の単調拘束から導く角度上限としては正しい。

不足情報: 実際の MOC 点列・試し C の係数と KKT 残差、run_0114/0117 の時系列・判定区間・VERDICT、生産経路で再計算した δ_r。run 原データがないため、既存 CFD の収束・準定常判定は独立検証していない。**ファイル変更なし、plan 未反映**。呼び出し側で対象 plan §4.1・§6に反映すること。
