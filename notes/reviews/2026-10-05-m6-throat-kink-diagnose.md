# codex 諮問 (diagnose): m6-throat-kink

- **brief**: [`notes/reviews/briefs/2026-10-05-m6-throat-kink.md`](../../notes/reviews/briefs/2026-10-05-m6-throat-kink.md)
- **plan**: [`plans/active/verification-m6-axis-wave-mesh-su2.md`](../../plans/active/verification-m6-axis-wave-mesh-su2.md)
- **date**: 2026-10-05
- **commit**: `22096e09` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 10.4 min, rc=0
- **結論**: **MOC 網を固定し、壁抽出だけを「累積流束閉包」から「同じ網の θ による流線接続」へ替える形状 A/B を次の一手とする。**
- **extra**: `case/45.isobutane_m6_d155/throat_kink_probe.py`, `case/45.isobutane_m6_d155/_band_ab/throat_kink_probe.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| (a) 始点を除外して曲率を自由化 | **却下・Major** | 原因の切り分けにならない。始点の位置・接線・曲率は現行の明示的な接続条件である（[wall_axismach.py:76](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:76)）。まず入力場と壁抽出を分けて診断する。 |
| (b) 原因側を調べる | **採用。ただし修正案は要再検証・Major** | `m*` は既に `_flux_along` で計算され、初期累積流束にも同じ関数を使う。C⁺ 上も同じ中点式である（[moc_inverse.py:543](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:543)、[同:413](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:413)、[同:308](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:308)）。**「求積を網と同じにする」は変更にならない。** 対案は下記の固定場 A/B。 |
| (c) `a/2` を理論誤差としてゲートを緩和 | **却下・Major** | `a` は始点を除いた有限区間の回帰値であり、独立に導いた理論値ではない（[throat_kink_probe.py:17](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_kink_probe.py:17)）。既知の不整合を、その測定値で許容する循環になる。現行ゲートを維持する。 |
| H2「物理壁への影響は局所だけ」 | **却下・Major** | 推定した κ_t を入口からスロートまでの Hermite に渡している（[wall_axismach.py:339](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:339)）。対案は、物理壁と上流 Hermite の全域を評価対象にすること。CFDへの影響量は未確認。 |

結論: **MOC 網を固定し、壁抽出だけを「累積流束閉包」から「同じ網の θ による流線接続」へ替える形状 A/B を次の一手とする。**

第 1 仮説: **Hall 初期線と MOC 適合条件の不整合が網の角度場へ入り、壁抽出だけでは消えない角度オフセットを生んでいる。** 確度: **中**

  根拠:

- **始点の θ は強制上書きされていない。** `throat_char` は Hall の M・θ をそのまま使い、上書きは `vertical` の場合だけ（[moc_inverse.py:527](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:527)、[同:539](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:539)）。独立再計算で **θ_Hall(0,1)=−1.48e−16°**。局所の流線曲率は **約0.47409**で、Hall 場自体はゼロ角から出発する。**曲率の5%差と有限の折れ角は別問題である。**
- 初期線全体について、実装の C⁻ 適合式（[moc_kernel.py:442](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_kernel.py:442)）に対応する  
  `E = Σ[Δ(θ+ν) − S₋]`  
  を読み取り専用で再計算した。軸端の `sinθ/r` は隣接点から極限を近似すると、次の値だった。

  | `n_start` | TP | Hall と同じ γ の CPG |
  |---|---:|---:|
  | 41 | −0.13508° | −0.11903° |
  | 161 | −0.13231° | −0.11626° |

  最後の一区間の **約−0.009°だけで影響の上限を判断できない**。ただし、この積算値をそのまま折れ角へ換算する根拠もない。
- 生産設定の壁を再計算し、**a_r=0.09098°、a_θ=0.09623°、最初の割線−平均角度差=0.04410°**を再現した。

  反証条件: **網の状態を一切変えず、壁抽出だけの変更で a_θ が80%以上減少するなら、「抽出方法によらず残る角度場の不整合が支配する」という第1仮説を撤回する。**

第 2 仮説: **C⁺ 上の流束積分・交点抽出が、角度オフセットの大部分を生む。** 確度: 低〜中。中点積分と線形交点抽出は確認できるが、寄与量は未確認（[moc_inverse.py:308](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:308)）。

判別 A/B: **変えるのは壁の閉包式だけ。CFD 0 step、同じ網から各1回。**

- 共通: `n_axis_inv=2400`、`n_start=161`、`axis_dx0=0.03`、既存の `ds`・ガス・初期線を固定。161点は近傍回帰に12点を確保するためで、解消策として扱わない。`levels`、`init_cum`、`m*` を両腕で共有する。
- **A:** 現行 `cplus_flux_wall`。
- **B:** 始点 `(0,1,θ_Hall)` から、各 C⁺ 線との交点を順に、  
  `rᵢ−rᵢ₋₁ = ½(tanθᵢ₋₁ + tanθᵢ)(xᵢ−xᵢ₋₁)`  
  で決める。交点の状態は同じ網の線分から補間する。Delaunay・壁当てはめ・曲率拘束変更は入れない。連続する交点が得られなければ診断不能として止める。
- 主指標は同じ `0<x<0.08` の **a_θ**。併せて a_r、κ、位置差、および **B の交点を現行流束式へ戻した `Fᵢ−m*`**を記録する。Bが構成上小さくする割線残差だけで成功判定しない。

  → **Bで |a_θ| が80%以上減れば第1仮説を撤回し、流束閉包側を優先する。減らなければ「壁抽出の変更だけで80%以上解消できる」という第2仮説を棄却する。** 後者でも第1仮説の確定とはしない。これは診断基準であり、生産の形状ゲートは変更しない。

やらない方がよいこと: **κ=0.4735への置換、`m*` の単独調整、始点区間のゲート緩和、AWSでのEuler先行。** 今回の A/B は流束閉包と角度場の寄与を分ける診断であり、Bをそのまま生産方式として採用する試験ではない。

呼び出し側の前提への異議:

- **「内部では整合した流線」「始点だけが折れている」はまだ解釈。** 位置・角度の有限区間回帰が近いことは観測だが、右極限の存在や各区間の整合性までは証明しない。
- **「解像度に依存しない」は広すぎる。** `n_axis` を変えても近軸の `axis_dx0` は固定される（[moc_inverse.py:460](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:460)）。追加再計算では `axis_dx0: 0.03→0.003` により a_θ が **0.09623→0.08146°**へ動いた。主因を解消する変化ではないが、完全な不感ではない。
- **κ_t≈0.494 は局所微分の直接測定ではない。** 保存された δ_r の写しと現行コードから再構成すると **κ_t=0.494227、±0.15窓内8点**。窓±0.08／±0.30では **0.507162／0.480847**だった。実装は位置と半径を別途推定した後の1係数最小二乗である（[wall_axismach.py:379](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:379)）。
  同じ推定位置・半径を固定して分解すると、**設計壁0.491361 + δ_r項0.003748 − 基準半径のずれ0.000881 = 0.494227**。これは再構成形状の診断値であり、`run_0051` 本体の再検証値ではない。上流 Hermite への伝播もあるため、CFD影響を「局所だけ」とは限定できない。

不足情報: **`run_0051_ns_final_c2` 本体がこのワークツリーに無く、CFD影響・収束・準定常性は再検証できない。** Euler比較は、形状候補がゲートを通過してから必要になる。ファイル変更・`forge` 起動は行っていない。**plan 未反映**。呼び出し側で §5.1 #13・§9へ採否と診断基準を記録すること。
