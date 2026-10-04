# codex 諮問 (diagnose): sern3d-r5h-ledger-result

- **brief**: [`../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r5h-ledger-result.md`](../../../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r5h-ledger-result.md)
- **plan**: [`../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md`](../../../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md)
- **date**: 2026-09-27
- **commit**: `fc5483f3` (feature/gap-heating-precision)
- **codex**: effort `high`, 4.0 min, rc=0
- **結論**: **次は同一場・同一格子で、面1535830の低温点側の速度再構成だけをセル値へ戻す1更新のA/Bを行い、H1′の速度寄与を分離する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（実測値はブリーフ記載。対象 run は両 checkout に存在せず、数値の独立再検算はできていない）:

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | 分岐 A 成立、EOS 単独原因の除外 | **対象の初回応答に限定して採用**。`case/46.sern_design/run_0978_r5h_ledger_A/` と `run_0979_r5h_ledger_B/` の報告では、面和と対流残差の差が ρe で ≤2.3e−7、段別差は対流に集中し、EOS・BC 処理の正味保存量変更は 0。これは再構成依存を支持する。ただし「低温点が形成された履歴でも EOS は無関係」まで広げない。また、**面和→残差の閉合と、陰解法による更新の収支は別**である。対案は、この二つを区別して記録すること。 |
| **Major** | H1′を「低温スポットの確定原因」「粘性加熱と釣り合う105 Kの定常」とする | **要再検証**。全域の `convMethod: 1→0` は両側の ρ・P・速度等を同時に変える。+x 面が差の約91%を担うことは、その面の**低温点側速度だけ**が原因である証明ではない。また `res_roe` の粘性寄与は内部エネルギーの加熱量そのものではない。対案は「後流側面の再構成が局所エネルギー残差差を支配し、速度外挿が有力な寄与」と記録し、下記の一変数介入で検証する。根拠: [再構成](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:186)、[運動量・エネルギー流束](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:609)。 |
| **Major** | 「相手値を超えないので Venkat は正常」としてリミッタ実装を除外 | **却下〔除外根拠として〕**。実装は**全接続近傍の成分別 min/max**を使い、その面の相手値との二点間を制限しているわけではない。速度の大きさだけでも判定できない。さらに ε² を持つ Venkat は厳密な最大値原理を保証しない。対案は、実効 `limiterScaled`、各速度成分の近傍範囲、制限前増分、確定 ψ を照合すること。根拠: [近傍範囲](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiter_d.cu:323)、[Venkat の式](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiterFunctions_d.cuh:41)。**現時点でバグの証拠もない**。 |
| **Major** | 面 h₀ を両セル h₀ の範囲に制限すれば解決する | **要再検証**。示されたのは面 h₀ が**自セル**を上回ることだけで、両セルの範囲から逸脱する証拠はない。対案は、まず相手セルの h₀ も出し、候補制限が実際に発火するか確認すること。エネルギー流束だけを切る案は、質量・運動量・組成との整合を別途設計する必要がある。根拠: [面 h₀ の構成](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:394)、[流束への使用](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:614)。 |

結論: **次は同一場・同一格子で、面1535830の低温点側の速度再構成だけをセル値へ戻す1更新のA/Bを行い、H1′の速度寄与を分離する。**

第 1 仮説: **低温点側の大きな速度外挿が、後流側面の移流全エンタルピーと質量流束を変え、観測された対流エネルギー残差差の大部分を生んでいる。** 確度: **中**。

  根拠: `case/46.sern_design/run_0978_r5h_ledger_A/` と `run_0979_r5h_ledger_B/` の報告では、全差6.358のうち面1535830が5.76を担う。ブリーフの丸め値から計算すると、低温点側の運動エネルギーはセルで約1.91e4 J/kg、面で約9.42e5 J/kg。その増加約9.23e5 J/kgは、面とセルの h₀ 差約9.76e5 J/kgの**約95%**に相当する。これは速度外挿を優先して切り分ける根拠になる。

  実装も、再構成速度から運動エネルギーを作り、熱力学的エンタルピーへ加えて流束に使っている（[コード:279](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:279)、[コード:396](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:396)）。ただし、**95%は h₀ の増分の内訳であって、冷却原因の寄与率ではない**。

  反証条件: 他の再構成を維持して低温点側速度だけを戻しても、当該面のエネルギー残差差の大部分が残る、または質量・運動量を含めた内部エネルギー応答が加熱側へ動かないこと。

第 2・第 3 仮説:

- **第2仮説〔中、未確認〕**: 相手側状態・圧力・質量流束の連成、または一次化による拡散増加が主要因であり、低温点側速度だけでは説明できない。
- **第3仮説〔低、未確認〕**: 後縁近傍の勾配・再構成点・リミッタ評価の局所的な問題が、大きな速度外挿を作っている。現資料では実効設定と制限前増分が不足している。

判別 A/B: **変更は面1535830の低温点側速度外挿の係数 α だけ。**

- **A:** α=1、現行の再構成速度。
- **B:** α=0、同面・同側の3速度成分をセル値に戻す。ρ・P・組成、相手側の再構成、他面、粘性、BC は変更しない。
- 同一 restart から各1更新、`nStepInner: 1`。更新後の EOS を評価するところまで記録する。これは**診断用介入**であり、既存の設定キーがあるという意味ではない。
- 速度を変更した後、χ・質量・全運動量・エネルギー・それに従属する種流束を整合して再計算する。同じ内部面流束を両 CV に逆符号で加える。
- 観測対象に相手節点522112を加え、面流束、全保存量残差、実際の保存量増分、EOS前後を記録する。

事前の判定は次に固定する。

**→ 速度支配の結果なら:** 共通状態で、この介入が既報の当該面差5.76の**80%以上**を同方向に再現し、全残差から評価した内部エネルギー応答と実更新も加熱側へ動く。H1′の「局所残差を維持する速度再構成の機序」を採用する。ただし、正しい温度や恒久対策までは認定しない。

**→ それ以外なら:** 「低温点側速度が差の大部分を説明する」という第1仮説を棄却する。全域一次化の結果だけで速度制限を恒久対策にしない。収支が所定の1%基準で閉じなければ判定不能とする。

内部エネルギーの符号確認には、共通セル状態で保存量残差を

`ρV·de_int/dt = R_ρE − u·R_ρu + (|u|²/2 − e_int)R_ρ`

へ射影する。これは半離散残差の診断であり、**陰解法の実増分とは分けて照合する**。温度応答には組成変化も含める。`R_ρE` 単独の符号を加熱・冷却と呼ばない。

やらない方がよいこと:

- **有限厚後縁を先に入れて、低温点の原因が証明されたとすること。** 有限厚化は物理形状の設計課題として扱う。no-slip の有限厚ベースにも低速域は存在するため、「遅い節点の塊がなくなる」は保証されない。
- **h₀ のクリップや共通リミッタを直ちに全ケースへ適用すること。** 次の優先は **(b) の診断**。恒久実装は、診断結果を踏まえた別 plan と回帰検証の対象とする。
- **局所的だから (c) 監視のみへ格下げすること。** 力への影響はまだ測れていない。

呼び出し側の前提への異議: **「105 Kの定常」「粘性加熱との釣合い」は、現在の帳簿から言える範囲を超える。** 言えるのは、評価時点で対流と粘性の**全エネルギー残差寄与**がほぼ相殺していること。起点の plan 記録は `check_convergence: NOT CONVERGED`（[plan:1878](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1878)）で、低温点自身の定常性の VERDICT は提示されていない。

力係数について、最小の追加評価は、機序確認後に**同一形状・同一格子の局所介入を維持した対照**を準定常まで継続し、圧力・摩擦・モーメントを全対象面で比較すること。ただし得られるのは、その介入に対する感度であり、真の誤差上限ではない。両者の `check_quasisteady` の VERDICT と残差判定区間を示し、平均差に双方の時間変動幅を加えた幅で評価する。**1更新の力差、g3/g4差、新旧形状差から影響上限は出せない。**

不足情報: 対象 run の実効 config、局所帳簿の原本、運動量を含む閉合表、陰解法増分の照合、相手セルの h₀ と近傍成分別範囲。現ダンプの `entry→after_eos_bc` は複数処理をまとめた正味比較なので、「EOSとBCがそれぞれ変更ゼロ」とは区別する必要がある（[記録開始](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1427)、[記録終了](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1469)）。

ファイル変更・forge 起動はしていない。**plan 未反映**。呼び出し側の反映先は `plans/active/tooling-nozzle-sern-3d.md` §5.1 R5h。§4.35とR5h旧説明に残る確定的な原因記述も、今回の限定された判断と整合させること。
