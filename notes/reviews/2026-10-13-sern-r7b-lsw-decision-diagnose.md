# codex 諮問 (diagnose): sern-r7b-lsw-decision

- **brief**: [`notes/reviews/briefs/2026-10-13-sern-r7b-lsw-decision.md`](../../notes/reviews/briefs/2026-10-13-sern-r7b-lsw-decision.md)
- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **date**: 2026-10-13
- **commit**: `69ce3191` (feature/sern-design)
- **codex**: effort `high`, 4.5 min, rc=0
- **結論**: **まず既存の g3 六本と g4 参照一本を、元データ・元設定に結び付けて現行ゲートで再判定し、その成立前には追加の長時間 run を投入しない。**
- **extra**: `plans/active/convection-zero-thickness-edge-reconstruction.md`, `notes/investigations/2026-10-11-lsw10/DESIGN_DIFF.txt`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 重大度 | 論点 | 採否・根拠・対案 |
|---|---|---|
| Minor | Q1：主張できる改善量 | **採用。ただし順位も未認定。** 保存集計は ΔJ = 1.0659×10⁻³、時間変動幅だけで 7.35×10⁻⁶。非負の E_J を引いて 0.002 を超えることはない。順位には **E_J < ΔJ** が必要で、格子・領域感度に残る予算は約 **1.0586×10⁻³**。[DESIGN_DIFF.txt:16](/home/sano/work/forge-sern-design/notes/investigations/2026-10-11-lsw10/DESIGN_DIFF.txt:16) |
| Major | Q1：m10・m4 を省略する案 | **却下。** 小さい設計差は、小さい感度幅の証拠ではない。m10 は圧力推力の差 +1.3123×10⁻³ と、そこから算出される摩擦寄与の差 −1.2913×10⁻³ が相殺している。既存の m10・0.8 の格子感度では、1.0 側も別作動点も上から押さえられない。**同じ作動点の二水準を対にして測る。** [DESIGN_DIFF.txt:2](/home/sano/work/forge-sern-design/notes/investigations/2026-10-11-lsw10/DESIGN_DIFF.txt:2) |
| Major | Q1：同時領域拡大だけで代用 | **却下。** 方向別の影響が相殺しうる。親 plan の一方向ずつ＋同時拡大を省ける定量的根拠は提示されていない。R4 の記載も `run_0102` の Euler 加速点であり、今回の三作動点・二形状の設計差を保証しない。**計算を逐次投入し、判別不能が確定した段階で打ち切る。** [親 plan:428](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:428)、[同:439](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:439) |
| Major | Q2：格子感度の定義・起動履歴 | **設計差の感度を採用。純粋な格子誤差という解釈は却下。** 対象は G_B − G_A。各形状の絶対感度も別途記録する。起動履歴の違いは既に受理 plan が明示しており、STEADY だけでは同じ解への到達を証明しない。**当面は起動履歴を含む複合感度として扱い、g3 全本の取り直しは要求しない。** 純粋な格子効果・初期化非依存を主張するなら、同一格子で起動経路の比較が別途必要。[受理 plan:157](/home/sano/work/forge-sern-design/plans/accepted/tooling-sern-te-wake-grid.md:157) |
| Major | Q3：現行ゲートでの再判定 | **採用。追加 CFD より先。** 確認した保存 `metrics.json` は `acceptance_version`・`judgment` を持たない。現行検査は残差の完全性に加え、壁面部品・摩擦出力の欠落も拒否する。設計判断でもこれらは必要な証拠であり、「学習には使わない」ことは免除理由にならない。**元記録を保持して再判定結果を別記録にする。** [受理 plan:131](/home/sano/work/forge-sern-design/plans/accepted/tooling-sern-te-wake-grid.md:131)、[同:136](/home/sano/work/forge-sern-design/plans/accepted/tooling-sern-te-wake-grid.md:136) |
| Major | Q4：C_M の扱い | **暫定制約として維持。実機採用判断は保留。** 作動点別窓はユーザが根拠なしとして保留している。現行のノズル面・圧力のみ・基準点 −20 H の定義で各候補の絶対値を制約と照合し、設計差は併記する。**ΔC_M を格子感度の許容 0.05 と比較したり、結果を見て窓や罰則を作ったりしない。** 推力順位の認定と実機トリム成立を分ける。[親 plan:438](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:438) |

結論: **まず既存の g3 六本と g4 参照一本を、元データ・元設定に結び付けて現行ゲートで再判定し、その成立前には追加の長時間 run を投入しない。**

第 1 仮説: 既存の判定区間は、追加 CFD なしの再判定で比較の基準列として再利用できる。確度: **中**。  
根拠: `case/46.sern_design/` 配下の `run_1084_tewake_g3_L10_m10`、`run_1086_tewake_g4_L10_cont40k`、`run_1096_tewake_lsw10_m10_cont20k`、`run_1098_tewake_lsw10_m4_cont20k`、`run_1099_tewake_lsw10_m6_cont40k` は、保存 CSV の四力量を再計算して平均・a・d が集計表と整合した。各履歴は 500～20000 step の40標本で、力量に非有限値・予定 step の欠落はなかった。保存された判定は次のとおり。

- `check_quasisteady`: **OVERALL: ALL STEADY**（力量・最低温度）。
- `check_convergence`: **NOT CONVERGED (stalled/plateau)**。
- 力の平均・a は `(10000, 20000]`、d は直前の10000 step 窓との比較。延長 run もその run 内の区間である。

これは**保存された準定常判定の確認**であり、現行ゲートの全面的な再合格や収束解の認定ではない。索引は [case README:401](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:401)。  
反証条件: 必須の残差列・step・壁面部品・traction・場・床記録・来歴のいずれかが欠ける、または現行ゲートが不合格を返す。その run は基準列に戻せない。

第 2・第 3 仮説: 追加しない。現在の証拠では、順位を壊す原因を格子・領域・起動履歴のどれかへ絞れない。

判別 A/B: **同じ保存入力に対し、受理検査だけを旧版／現行版に変える。forge は 0 step。** 対象は上記五本に `run_1090_tewake_m4_L10`・`run_1091_tewake_m6_L10_cont40k` を加えた七本。力量・温度・床・全残差・出力部品・判定区間を照合する。  
→ **結果 A：全必要条件を満たす**なら、既存基準列の再利用を支持し、感度評価へ進む。  
→ **結果 B：欠落・不合格がある**なら、「旧 PASS のまま再利用できる」を棄却し、欠けた証拠の回復を先にする。欠測は物理的不成立と区別し、版の付け替えだけで救済しない。

呼び出し側の前提への異議:

**E_J は現状では「保証された誤差上限」ではなく、登録した変更への総合感度幅である。** 二格子の差だけでは真の誤差上限にならない。以下を今後の事前登録案とする。

- q = `C_T_with_shear`、Δᵢ,g = qᵢ,1.0,g − qᵢ,0.8,g。
- 各 run の時間幅を u = max(a, d, 10⁻⁶)、対の時間幅を Uᵢ,g = uᵢ,1.0,g + uᵢ,0.8,g とする。
- 格子感度幅を **Gᵢ = |Δᵢ,g4 − Δᵢ,g3| + Uᵢ,g4 + Uᵢ,g3** とする。
- 領域感度幅を、事前固定した各拡大条件 s に対する **Dᵢ = maxₛ(|Δᵢ,s − Δᵢ,基準| + Uᵢ,s + Uᵢ,基準)** とする。
- **Eᵢ = Uᵢ,g3 + Gᵢ + Dᵢ** とする。時間幅を重複して数える保守的な式だが、差の差が測定幅より小さいときにゼロ同然と扱う誤りを避ける。各形状の絶対感度に対する既存の四力量の検査も残す。

**最小集合は、採択と打切りで異なる。** 現在の証拠から、順位を採択できる「m6 の g4 二本＋領域二本」という縮約は認められない。

- 格子感度を三作動点で揃えるには、`run_1086` の再利用が成立する条件で、**追加 g4 は五本**。m6 の二本を先にし、m10 の1.0一本、m4 の二本を続ける。
- 領域感度は、今回へ適用できる既存証拠がない条件について、**拡大条件一つ当たり六本**が必要。方向を省くには、その方向・各作動点・両形状について省略可能な感度幅の根拠が要る。
- **一括投入は不要。** 固定した基準 ΔJ に対し、測定済みの非負の感度寄与だけで E_J ≥ ΔJ となれば、その登録規則では順位を支持できないため、残りを打ち切れる。逆に、未測定項をゼロとして早期採択してはいけない。
- 各新規評価の長さは登録済みの判定区間20000 step、窓未達なら＋20000を一回まで。短窓への変更では節約しない。

やらない方がよいこと: 旧格子との設計差の近さを格子独立性の代用にすること、m10 の片側感度を全作動点へ流用すること、同時領域拡大だけで方向依存を除外すること、E_J ≥ ΔJ を「0.8 の優位」と解釈すること。

不足情報: 指定されたローカルの写しには生の残差 CSV・場／壁面 HDF5・完全な実効 config がなく、七本の現行ゲート再判定は未実施。`run_1090`・`run_1091` の個別時系列も今回指定された三ディレクトリにはない。AWS-B に必要入力が残っているかの確認が必要である。

**plan 未反映**。依頼どおりファイル変更・forge 起動は行っていない。採用する場合の反映先は、親 plan §5.1 R7b ④と、そこで参照する E_J・逐次打切り条件である。
