# codex 諮問 (diagnose): sern3d-r4d-domain-design

- **brief**: [`../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r4d-domain-design.md`](../../../forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r4d-domain-design.md)
- **plan**: [`../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md`](../../../forge-sern-design/plans/active/tooling-nozzle-sern-3d.md)
- **date**: 2026-09-27
- **commit**: `d8c11f52` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.0 min, rc=0
- **結論**: ?

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。参照先はブリーフが指定する `/home/sano/work/forge-sern-design/`。**現案のまま6変種を投入することは推奨しません。**

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **g1だけで領域独立性を認定：却下** | [plan:1898](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1898) の g1→g3 差は、圧力のみの `C_T` でも −0.002659、`C_L` は +0.000876。提案の領域許容 0.0005 を超える。「摩擦が主」は、領域感度が格子に依存しない証拠ではない。**g1は予備調査に限り、認定はg3で行う。** |
| **Major** | **V1/V2をYAML指定だけで実施：却下** | [runner:60](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:60) は `top_out_kind == "outlet"` の場合だけ出口処理を呼び、それ以外は `slip`。`side_far` は63行で `slip` 固定。関数を独立実行し、**`top_out_kind: outflow` でも生成結果が `slip`** と確認した。V2の現行指定は `top_out_kind: outlet` と `outlet_kind: outflow`。V1には生成器側の対応が必要。投入前に生成済みBCの差分を検査する。 |
| **Major** | **slip→outflow＝反射の有無：却下** | `outflow` は全量外挿であり、非反射境界の保証ではない。[運用手順:30](/home/sano/work/forge-sern-design/procedures/divergence-and-startup.md:30)。判断に必要なのは全Machでなく **Mₙ = U·n/a**。自由流がx方向なら `side_far` のMₙは0で、流入音響特性が残る。超音速流入がなくても問題になる。**BC変更感度と非反射性の検証を区別し、流入特性に外気状態を与える遠方境界の実装・適用範囲を確認する。** |
| **Major** | **寸法だけ変えれば1因子試験：却下** | [mesher:108](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:108) は出口位置からランプ上・後流の点数配分を再計算する。V6はランプ上の解像度まで変える。V3も[168行](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:168)、V4も[438行](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:438)、V5も[223行](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:223)で分布を再生成する。**既存領域の格子を保持し、外側へセルを追加する。点数の比例増加だけでは不十分。** |
| **Major** | **全係数で§8の1/4を確保：一部却下** | [plan:1865](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1865) の `C_M` 格子差は0.0404948。提案の0.0125を加えると **0.0529948 > 0.05**。振幅を入れる前から予算超過。**`C_M` の領域許容は0.005に下げる。** 他3量の0.0005は試験上の許容として採用する。 |
| **Major** | **各変種が個別合格ならR4d完了：却下** | 個別変化の累積・相互作用を扱えていない。また[§8:1972](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1972)は生産3作動点が対象。**1因子試験は原因分離に採用するが、最終BCを固定した同時拡大領域の確認を追加し、m6だけでR4dを閉じない。** |
| **Minor** | **`top_out` はランプ延長壁という説明：却下** | 現行 `ext_top` では、[mesher:440](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:440)の `b(x)+top_depth` が外側上面となり、460行で `top_out` に分類される。物理的な `vehicle_top` と別面。旧コメントを現トポロジの根拠に使わない。 |

許容と評価方法は、次のように事前固定することを推奨します。

| 量 | 領域感度の許容 ε |
|---|---:|
| `C_T` | 0.0005 |
| `C_T_with_shear` | 0.0005 |
| `C_L` | 0.0005 |
| `C_M` | **0.005** |

比較量は **D = |平均B − 平均A| + aA + aB**、振幅は **a = max|各保存時点の値 − 平均|** と定義します。標準偏差との混用は避けます。

最終認定では、同じ定義で求めた格子差の幅Gについても **G + D ≤ §8の総許容**を要求します。これは観測された感度の予算管理であり、真の誤差上限ではありません。planの掲載値では `C_M` は 0.0404948 + 0.00090 + 0.00034 + 0.005 = **0.0467348**ですが、最終値ではなく比較窓の平均で再計算が必要です。

**結論:** g3で共通領域の格子を保持した **`Z_ext: 1.5 → 2.25` の領域拡大A/Bを1組だけ**行い、現側方境界位置への感度を先に判定する。

**第1仮説:** この側方拡大によるg3の力係数変化が、上記の領域許容を超える。  
**確度:** 低。遮蔽する境界の存在は確認できるが、力への影響は未測定。  
**根拠:** [runner:63](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:63)で `side_far` は `slip` 固定、[mesher:169](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:169)で現在の位置は z/H = W/2 + Z_ext = 2.5。  
**反証条件:** 共通領域の格子を保持し、両者の準定常性を確認した比較で、4量すべてD ≤ εとなること。ただし反証できるのは「この拡大への許容超の感度」であり、無限遠への独立性ではない。

**第2仮説:** 現行メッシャをそのまま用いた寸法比較では、格子再配分が領域感度に混入する。コード上の再配分は確認済み、力への寄与は未確認。  
**第3仮説:** 設定しない。

**判別A/B:**

- A：g3、`Z_ext=1.5`。B：g3、`Z_ext=2.25`。**BCは両方とも現行のまま固定**する。
- 共通領域の座標・接続を保持し、Bの追加領域だけにセルを足す。現行の分布再生成をそのまま使う比較は不可。
- 基準は `case/46.sern_design/run_0971_3d_g3_chidef_cont16k/`。新しいrunを2本作り、Aは同一格子restart、Bはcross-mesh restart。共通部の補間差を記録する。
- 各20000 stepを初回枠、500 step間隔で保存。末尾10000 stepで評価し、その前の10000 stepとの平均差も各量 **0.1ε以下**を要求する。未達なら延長し、短窓で判定しない。
- 両者の4量に対する `check_quasisteady` の **STEADY**、床・NaN検査、メッシュ品質判定を必要条件とする。`check_convergence` のVERDICTも併記する。残差プラトーは[§4.39のユーザ決定](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1593)に従い、それだけで却下しない。
- **分岐A：** 1量でもD > εなら第1仮説を支持し、現側方領域は不合格。  
  **分岐B：** 全量D ≤ εなら第1仮説を棄却し、この拡大に対する感度は許容内とする。  
  準定常性・格子保持条件を満たさなければ判定不能。

**やらない方がよいこと:** V1/V2を「非反射化」として直ちに生産採用すること、g1の無感度をg3へ転用すること、全変種に独立に誤差予算を配って合算を省くこと。V3～V6は省略せず、まずV3だけを上記条件で実施する。

**呼び出し側の前提への異議:** `bottom` のBC適合性も未解決です。[plan:141](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:141)には旧runでMₙ≈0だった記録があり、現runnerでも `outlet` と `bottom` は同じ設定で切り替わります。最大圧力偏差だけでなく、各遠方境界のMₙ分布・流入面積率・組成を確認対象へ追加してください。V1/V2で適切なBCへ変更して係数が大きく動いた場合、その差の小ささを理由に旧BCを残してはいけません。

**不足情報:** 現runの境界Mₙ分布、比較窓の平均・振幅、判定区間付きVERDICT原票。既存runの `GATES: PASS`・`STEADY`・`NOT CONVERGED` は今回planの記載を参照しており、原票は独立確認していません。ファイル変更・forge実行は行っていません。**plan未反映**であり、採用時の反映先は当該planの§5.1 R4d・§6・§8です。
