# codex 諮問 (diagnose): double-solve-trajectory

- **brief**: [`notes/reviews/briefs/2026-10-11-double-solve-trajectory.md`](../../notes/reviews/briefs/2026-10-11-double-solve-trajectory.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `d2906fff` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.3 min, rc=0
- **結論**: **追加比較は採用するが、「線形 solve 全体の double 化」ではなく「block 系の組み立て精度の感度試験」に訂正し、再現性を含む判定条件を事前登録してから C・C′ を回す。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **追加比較は採用するが、「線形 solve 全体の double 化」ではなく「block 系の組み立て精度の感度試験」に訂正し、再現性を含む判定条件を事前登録してから C・C′ を回す。**

採否表:

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major 1** | 「ライン Thomas もこの切り替えで double になる」— **却下** | `implicitSolvePrecision` は `implicit_defect_correction_block_d<ST>` を切り替えるが、ライン上では点解せず、`diag`・`rhs`・`K` を `flow_float` に保存する。[timeIntegration_d.cu:1184](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1184)、同ファイル:952。Thomas の分解・代入は別カーネルで、既定で **double**。[同:2339](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2339)、同:2417。対案は「行列・ライン外近傍寄与を含む右辺の組み立てを double にする試験」と定義すること。ライン非被覆点があれば、その点解の精度も変わる。 |
| **Major 2** | 結果から必要性・固定点を断定する解釈 — **却下** | C でも係数・右辺・補正の保存は float のまま。差が残っても、線形更新側を除外できない。差が抑えられても、組み立ての変更が他所の誤差への応答や過渡の速度を変えた可能性が残る。[同:1195](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1195)、同:2438。さらに既存 3 本は plan 記録上 `CHECK FAILURES`。[plan:742](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:742)。結論は「登録期間の偏りを、この設定変更で抑制できた／できなかった」に限る。 |
| **Major 3** | C・C′ の平均だけによる「ずれが消える」判定 — **却下** | 提案式には再実行差の制限がない。例えば W1・θ_r(40) で ΔC = ΔA + D、ΔC′ = ΔA − D とすると、両腕とも A から 1.892 ポイント離れるのに、提案の抑制条件は成立することを計算確認した。参照値は [plan:730](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:730)。各腕の距離と再実行差を判定に含め、B→C の効果には B・B′ の再実行差も使う。 |
| **Major 4** | 「共通入力なので精度間で実際の初期状態も同じ」— **要再検証** | §6.31 は初期化後の ρE 不一致を記録している。[plan:799](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:799)。小さい丸めだから長期軌跡への影響も小さい、とは未確認。§6.27 の精度依存の観測を捨てる必要はないが、初期化と反復中の差は分離されていない。今回の主対照は同じ float バイナリの B 対 C とし、初期化後の保存量が揃うことを確認する。A は参照軌跡として扱う。 |

第 1 仮説: **`ST` で行う block 系の組み立ての丸めが偏りを増幅しており、この部分の double 化だけで登録窓の偏りを 90 % 以上抑えられる。** 確度: **低**  
　根拠: 切り替えは Jacobian・近傍寄与の演算に実際に届く一方、Thomas 自体は既に double。§6.31 の更新差の分解は判別不能であり、この仮説を直接支持する実測はまだない。[timeIntegration_d.cu:946](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:946)、[plan:800](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:800)。  
　反証条件: 有効性・再現性の確認を通った C・C′ が、両窓・全断面で A から D の向きに 0.5|D| 以上離れること。この場合、「この変更だけで十分抑制できる」を退ける。

第 2 仮説: 残差評価・SST・状態更新、または変更されない係数・補正の float 保存などが偏りを維持し、組み立て精度だけでは抑制できない。確度: 中、発生箇所は未確認。

第 3 仮説: 初期化の丸めや非決定的な揺れが遅い過渡を選び、今回の設定変更はその発達速度を変える。確度: 低、未確認。固定点の比較とは分ける。

判別 A/B: **float の `time.deltaT.implicitSolvePrecision` だけを 0→1 にする比較を 1 つ行う。**

- 既存 B・B′を対照、C・C′を追加。共通入力・SST 更新あり・`fgeom7-f32`・30,000 step 固定・500 step ごと。W1/W2、θ_r の 3 断面、D と Δ の定義は維持する。Q_w は副判定。
- 各窓・断面で mB = mean(ΔB, ΔB′)、mC = mean(ΔC, ΔC′)、N_B = |ΔB − ΔB′|、N_C = |ΔC − ΔC′| とする。
- **結果 A「抑制」**:  
  max(|ΔC − ΔA|, |ΔC′ − ΔA|) + 10N_C ≤ 0.1|D|、  
  かつ sign(D)·(mB − mC) ≥ 0.5|D|、  
  かつ |mB − mC| > 10max(N_B, N_C)。  
  → 第 1 仮説を支持する。ただし、誤差の発生源や必要性の証明にはしない。
- **結果 B「抑制不十分」**:  
  C・C′がそれぞれ A から D の向きへ 0.5|D| 以上離れ、|mC − ΔA| > 10N_C。  
  → 第 1 仮説の「この変更だけで十分」を棄却する。残差・commit だけが原因とは結論しない。
- 両窓・全断面で揃わなければ判別不能。これらは診断用の暫定基準であり、信頼区間ではない。

**既存 `run_0487`〜`0489` の再利用は条件付きで採用する。** 実行日時が違うことだけを理由に回し直す必要はない。ただし、次を実物で確認する必要がある。

- バイナリの完全なハッシュ、入力 HDF5・BC・実効設定が同じで、B の実効値が `implicitSolvePrecision: 0`、C が `1`。
- `FORGE_CUDA_BLOCKSIZE=128` を含む環境変数が同じ。`FORGE_FREEZE_TURB`・`FORGE_LINE_F32` など余分な設定が混入していない。元の台本は環境を清掃している。[tr.sh:9](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tr.sh:9)。
- B・B′・C・C′ の初期化後の保存量、出力点数、抽出座標・手順を照合する。step 0 欠損を別時刻で代用しない。
- 全期間の有限性と正式ツールの実行を**総合判定の前**に確認する。`tr_an.py` は判定を先に確定し、ツールを後から呼ぶため、そのまま複製しない。[tr_an.py:62](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tr_an.py:62)、同:104。

やらない方がよいこと: Thomas の double 化を新規実装すること、差が残っただけで LHS を除外すること、差が縮んだだけで根治とすること、判定窓の延長・選び直し、未収束の軌跡を固定点の一致と呼ぶこと。

呼び出し側の前提への異議: 観測として残せるのは「共通入力からの精度依存の軌跡差」と「1 step 診断では発生段を分離できなかった」まで。「反復中の特定演算が真因」「定常の固定点の偏り」は未確定。また、§6.27 の `ALL STEADY` は既定の緩い閾値による記録であり、今回の等価性の根拠にはならない。

不足情報: このワークスペースには対象 run、`q32_init.h5`、判定 JSON がなく、実測・VERDICT・実効 YAML・バイナリ同一性を独立確認できなかった。上記の run 数値は plan の記録として扱った。一方、確認した 2 つの C++/CUDA ファイルは `0b4dff4e` との差分がなく、切り替え範囲の指摘は対象ソースでも成立する。**ファイル変更・forge 起動は実施していない。plan 未反映であり、呼び出し側で §6 の事前登録と §5.1 に反映すること。**
