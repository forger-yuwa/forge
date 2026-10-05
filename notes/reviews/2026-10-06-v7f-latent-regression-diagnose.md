# codex 諮問 (diagnose): v7f-latent-regression

- **brief**: [`notes/reviews/briefs/2026-10-06-v7f-latent-regression.md`](../../notes/reviews/briefs/2026-10-06-v7f-latent-regression.md)
- **plan**: [`plans/active/thermophysics-solver-owned-species-db.md`](../../plans/active/thermophysics-solver-owned-species-db.md)
- **date**: 2026-10-06
- **commit**: `9582f0b1` (feature/species-transport)
- **codex**: effort `high`, 4.4 min, rc=0
- **結論**: **V7(f) は現状では閉じず、MW・datum・入力を固定して潜熱実装だけを替える 0-step A/B を次の一手とする。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 重大度 | 提案・前提 | 採否・根拠・対案 |
|---|---|---|
| **Major** | 液のある節点が 200 K 以上なら、#10 は影響しない | **却下**。核生成の Kantrowitz 補正は液がまだない状態でも L を使う。[condensationSourceF_d.cuh:23](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:23)、同 `:41–60`。`g > 1e-8` の抽出は onset 前の過飽和セルを除外する。**対案:** 液相域に加え、核生成対象・中間状態・EOS 反復中の評価温度を確認する。 |
| **Major** | double の ΔL が相対 3e-15 なら、float 経路も不変 | **却下**。表は L を単に float に丸めるのではなく、刻み 1e-4 K の差分から Hermite 係数を作る。[condensationTables_d.cuh:88](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTables_d.cuh:88)。小さい L の差でも微分係数に残り得る。**対案:** 旧新の係数全成分と、実際の評価関数による L・L′を比較する。 |
| **Major** | 現行バイナリと旧実装の差は、200 K 以上で丸めだけ | **要再検証**。#13-3 で MW が変更されている。現行単体試験の「丸め程度」は、**新 L に MW 比を掛けた比較**である。[test_cond_latent_pair.cu:161](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_cond_latent_pair.cu:161)。**対案:** #10 単独は同じ MW・係数・datum で比較し、旧場から現行モデルへの移行差を別記する。 |
| **Major** | ΔL 表と温度分布で V7(f) を完了にする | **却下**。登録条件は onset・g の変化量と時系列判定であり、物性差の評価とは別。[plan:383](/home/sano/work/forge-species/plans/active/thermophysics-solver-owned-species-db.md:383)。**対案:** 影響評価を補助証拠として追加し、V7(f) は未完了のままにする。代替検証へ変更するなら、免除する範囲と未測定の量を明記する。 |

結論: **V7(f) は現状では閉じず、MW・datum・入力を固定して潜熱実装だけを替える 0-step A/B を次の一手とする。**

第 1 仮説: **200 K 以上でも、#10 の演算順序変更による丸め差が float 表の微分に残る。** 確度: **中**（実バイナリ未確認）。

  根拠: 表の構築・評価は `condensationTables_d.cuh:85–98,44–54`。これを Python の scalar double と、FMA を使わない float32 演算に移植して検算した。旧 MW `0.0180153`、新方式の datum `298.15 K` に固定した結果:

  | 比較対象 | 検算結果 |
  |---|---:|
  | 207.5–400 K を覆う 771 区間の係数不一致数（c₀, c₁, c₂, c₃） | 0, 56, 762, 754 |
  | 同範囲の 100,001 点で L の評価値不一致 | 0 点 |
  | 同じ点で L′の評価値不一致 | 12,812 点 |
  | 最大 |ΔL′| | 4.8828125×10⁻⁴ J/(kg·K) |

  **これは式の移植による検算であり、生産ビルドのビット比較でも CFD の影響測定でもない。** 差は小さいが、「L が同じなら L′も同じ」という推論は成立しない。

  反証条件: 同一ビルド条件・同一 MW・同一 datum で、実際の表構築と device 評価を比較し、対象温度域の L′が全点ビット一致すること。

第 2 仮説: **抽出対象外の核生成前セルや EOS の試行温度が 200 K 未満を通る。** 確度: **低・未確認**。核生成は `g` を入力条件に持たず、EOS は反復温度で L を評価する（`condensationSourceF_d.cuh:41`、`condensationEOS_d.cuh:22–35`）。

第 3 仮説: **現行と旧モデルの比較に #13-3 の MW 差が混入する。** 確度: **高・比較対象次第**。MW だけで L は約 +1.11017×10⁻⁶ 相対変化する。今回の式の検算では 300 K で約 +2.70692 J/kg。これは #10 単独の影響ではない。

判別 A/B: **A＝#10 以前の潜熱関数、B＝気液ペア方式。変更点は潜熱の供給関数だけ。時間積分は 0 step。**

- 両腕の MW は `0.0180153` に固定し、気相係数・datum・表格子・コンパイラ設定も揃える。**生成後の表への MW 比掛けでは代用しない**。
- 同じ `cond_tables_fill` で全係数を構築し、成分別のビット不一致数を出す。実際の device 評価で L・L′を比較する。200 K と両隣の float、200 K をまたぐ区間 **199.90–200.15 K**、各区間内部、退避境界を含める。
- 同じ入力状態を固定して、核生成・成長ソース、二相 EOS 残差と傾き、音速まで差を追う。150 K を正の対照に含め、既知の ΔL ≈ −465.920 J/kg が出ることも確認する。
- **対象温度域の L′まで完全一致なら第 1 仮説を棄却。差が残れば「float 化ですべて消える」という閉鎖根拠を棄却。** 係数だけに差があって利用結果に届かない場合も、その区別を記録する。

この比較は旧 CFD 場の変換や `forge` 起動を必要としない。ただし、**0-step の一致だけでは onset・g の湿潤回帰を完了したことにはならない**。

やらない方がよいこと: `condFloat=0/1` を旧新比較の代わりにすること、係数差の存在だけで「CFD に有意な影響」と断定すること、最終湿潤セルの最低温度だけで低温経路を除外すること。

呼び出し側の前提への異議:

- 「200 K 以上」は評価温度そのものについて確認が必要。double の L′は T±0.1 K を使うため、200–200.1 K でも低温側を読む（`condensationEOS_d.cuh:28`）。表の境界横断も同様。
- 二相熱容量・音速には `cp₂ = cp_allvap − gL′` として差が入る。ハイブリッド反転は表の傾きと double の L を併用するため、float の L だけの一致では足りない（[condensationEOS_d.cuh:184](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:184)、同 `:352–359`）。
- Python 変換は潜熱モデルの変更自体でエネルギー再構成を発火させる。気相・g・Rw が同じなら補正は **Δroe/ρ = −gΔL**。200 K 以上でも丸め差や MW 差があれば非ゼロになり得る（[convert_species_field.py:203](/home/sano/work/forge-species/solver_density_cuda/tools/convert_species_field.py:203)、同 `:742–750,847–853`）。補正式は妥当であり、このレビューでは欠陥とは判定しない。

不足情報: AWS の対象 run・HDF5・実効設定・バイナリ識別情報・判定区間と VERDICT。本環境には `run_0567`、`run_0561`、`run_0482`、`run_0524` がなく、提示された最低温度は独立確認できなかった。「現行の全 H2O 凝縮ケース」を網羅した証拠も不足している。

ファイル変更・`forge` 起動は行っていない。**plan 未反映**。呼び出し側で `plans/active/thermophysics-solver-owned-species-db.md` §5.1 #10・§6 V7(f) に採否を記録すること。
