# codex 諮問 (diagnose): twophase-diffusion-4e

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-diffusion-4e.md`](../../notes/reviews/briefs/2026-10-02-twophase-diffusion-4e.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `ad106bfb` (feature/species-transport)
- **codex**: effort `high`, 5.2 min, rc=0
- **結論**: **#4e の受入を保留し、まず液枯渇時の θ の丸めを、同じ1セル入力に対する GPU 更新の A/B で判別する。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 論点 | 重大度・採否 | 根拠と対案 |
|---|---|---|
| Q1：輸送単独の累積保存試験も免除する | **Major／却下** | [plan:121](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:121) の免除は人工ソースを含む S9 に限定され、輸送単独の条件は [137行](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:137) に残っている。**試験側の BE 解法を直す方針を推奨**する。既に液量で改善した `V/Δτ=0` を正式候補とし、BE 方程式・6ε 床・保存許容 1e-6 は維持して、独立残差と全対象量を再判定する。元の 3.14e-6 FAIL は記録に残す。 |
| Q2：`rms_roYv` の低下だけで独立評価を代替する | **Major／却下** | [残差実装:708](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:708) は丸め済み残差の減算。初期値から何桁低下したかでは、局所の演算誤差床との距離を保証できない。最小の対案は、通常監視を維持し、**収束候補を受け入れる段階で格納状態から EOS・流束・ソースを独立再評価する監査**を追加すること。蒸気を含む全成分を個別尺度で判定し、既存 float 残差の double 化で代用しない。 |
| Q3：点対角、φ_N δρ の省略、Σz 正規化、定 Sc | **Minor／採用。ただし説明を限定** | [更新:143](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:143)、[共通再正規化:2120](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2120)、[定 Sc:2079](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2079) は定常初版の選択として妥当。φ_N δρ の省略は共通再正規化を採る更新設計として支持するが、`ρ/ΣρY` は種の更新にも依存するため、「両者は常に同じ密度補正」という代数的同一視はしない。点対角の到達性は CFD では未検証。 |
| Q3：θ_b 撤去・θ≥0 による非負保証 | **Major／要再検証** | 旧 θ_b の復活は不要だが、[162行](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:162) の `double → float` 変換で θ が切り上がり、[173行](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:173) の液量が負になる反例がある。**共通の実効 θ を、安全側へ丸める修正**を推奨する。後段の `passiveBounds` に吸収させて非負保証と呼ばない。 |
| 受入1：面恒等式の PASS | **Major／要再検証** | [試験:315](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_kernel.cu:315) は許容の分母を `8ε(Σ|j⁰|+ta)` に広げているが、表示は `8εΣ|j⁰|` のまま。[事前条件:638](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:638) の達成を示していない。分子恒等式は同じ状態で `ct=0` として直接検査し、乱流の差し引き誤差を混入させない。 |
| Q4：#4e を受け入れて AWS の #1b へ進む | **Major／現時点では却下** | 独立残差契約と上記受入条件が未完了。さらに [更新後処理:819](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:819) には再正規化と bounds がある。前処理だけなら原方程式の零点は変わらないが、**PASS・STEADY だけで補正による固定点の変化まで除外できない**。受入後の場の比較には、両者の独立残差合格、制限解除、再正規化・floor・射影の補正が無視できることも必要。収束速度は比較しない。 |

結論: **#4e の受入を保留し、まず液枯渇時の θ の丸めを、同じ1セル入力に対する GPU 更新の A/B で判別する。**

第 1 仮説: **共通 θ の float 化が非負境界を越え、後段の bounds が液の負値を補正している可能性がある。** 確度: **高。ただし GPU と本番での発生は未確認。**

根拠: [制限・commit:158](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:158) は double で求めた θ を float に変換して更新し、液の丸め負値には補正を設けていない。次の入力で、ホスト上の IEEE float32 演算を検算した。

- `M=V=ρ=ω=1`、輸送対角・ソース Jacobian は0。
- `rYw=0.01f`、`rg=3e-6f`、`Rw=0`、`Rg=-1.45e-4f`。
- `dg_max=0.005`、`dT_max=1`、`L=2.5e6`、`cveff=1000`。Q の残差は0。

このとき θ=0.0206896561673 に対し、float 化後は 0.0206896569580 と大きくなる。

| θ の変換 | 更新後の液量：積と和を別々に丸める | 更新後の液量：FMA |
|---|---:|---:|
| 現行の最近接丸め | −2.2737368e-13 | −1.1464300e-13 |
| 切り上がった場合だけ隣の小さい float へ戻す | +2.2737368e-13 | +1.5544055e-13 |

これは小さい誤差だが、「補正前に液が非負」という契約への反例である。**3セル累積保存 FAIL の原因と同一だとは主張しない。**

反証条件: 本番と同じコンパイル条件の `k_update` にこの入力を渡し、後処理前の現行出力が非負になること。その場合は演算順序との差を特定し、上の検算を GPU の反例としては採用しない。

第 2 仮説: **3セルの累積誤差は、解き残しと有限精度 commit の偏りが、前処理によって変わる現象である。** 確度: **中**。[設計メモ:670](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:670) の 3.14e-6 → 2.83e-7 は支持材料だが、変更後の診断は [試験:423](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_kernel.cu:423) で液量しか集計しておらず、全成分の保存を証明していない。

第 3 仮説: 追加しない。

判別 A/B: **変更点は θ の float への変換方法だけ**。上記1セルを既存の `k_update` で1回更新する。A＝現行の最近接丸め、B＝切り上がった場合だけ安全側へ戻す。状態・残差・前処理・閾値は同一とし、bounds 前の液・蒸気・Q、`qcut`、`vround` を見る。

- **Aだけ液が負、Bは全量非負・補正0なら**、θ の丸めによる穴と判定し、Bを採用する。
- **A・Bとも液が負なら**、θ の変換だけで直る仮説を棄却し、commit を調べる。
- **Aも非負なら**、GPUでの反例は不成立。コンパイルされた演算との差を確認する。

やらない方がよいこと: **累積保存の許容緩和、未承認の残差持ち越し、θ_b の復活、`rms_roYv` を double にキャストしただけの独立監査扱い。**

呼び出し側の前提への異議: **「前処理を変えると改善したので更新写像は原因でない」は採用しない。** 有限精度の commit も写像の一部である。また、T3 は [487行](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_kernel.cu:487) の独立再評価で停止しているため、その合格は本番の float 残差監視の妥当性を保証しない。

不足情報: 本番条件での独立残差と誤差尺度、補正が消える末尾区間、`V/Δτ=0` の全成分監査、上記 GPU 演算反例の結果。今回は指定コード・文書とホスト演算を確認した。forge・GPU 試験は実行せず、禁止された run 出力は未読、ファイル変更なし。**plan 未反映**。採否は呼び出し側で [plan](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md) §4.2・§5.1 #4e/#1b・§6 に反映すること。
