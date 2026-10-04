# codex 諮問 (diagnose): d7275-thermofloat-result

- **brief**: [`notes/reviews/briefs/2026-09-27-d7275-thermofloat-result.md`](../../notes/reviews/briefs/2026-09-27-d7275-thermofloat-result.md)
- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **date**: 2026-09-27
- **commit**: `ac727692` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.1 min, rc=0
- **結論**: **両腕を `thermoFloat=0` に固定し、SLAU の面エンタルピー評価だけを float／double に切り替える継続 A/B を次の一手とする。**
- **extra**: `case/60.flatplate_d7275_m7/README.md`, `case/60.flatplate_d7275_m7/acceptance.json`, `case/60.flatplate_d7275_m7/tools/tf_ab.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 判断 | 採否 | 理由・対案 |
|---|---|---|
| 自由流域では float 熱力学経路が主因でない | **限定して採用** | **Major:** 棄却できるのは「`thermoFloat` で切り替わるセル側経路が主因」という仮説まで。SLAU の面エンタルピーは両腕とも float のまま。[実装:396](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:396) |
| 第一内部列の流れには寄与が大きい | **観測として限定採用** | 記録上の `roUy` 比 0.028、`roe` 比 0.042 は大きな改善。ただし領域別残差は最終一枚の集計なので、持続する改善として確定するには保存済み時系列で再集計する。[tf_ab.py:55](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/tf_ab.py:55) |
| ω は保留 | **採用** | B の窓間変化 −7.7% は事前登録の 5% を超える。A の残差上位点で改善しても、第一内部列全体の解消にはならない。[acceptance.json:216](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:216) |
| 次は幾何を double にする | **現時点では却下** | ローカルの対象メッシュでは自由流域の閉包誤差が全点ゼロだった。AWS 側との同一性を確認し、先に未検証の面熱力学経路を切り分ける。 |

結論: **両腕を `thermoFloat=0` に固定し、SLAU の面エンタルピー評価だけを float／double に切り替える継続 A/B を次の一手とする。**

第 1 仮説: **`thermoFloat=0` でも残る面エンタルピーの量子化が、エネルギー・圧力の連成を通じて自由流域の残差床を維持している。** 確度: **中**。経路の存在は確認済み、床への因果関係は未確認。

根拠:

- SLAU は `thermalMethod==2` で `thermo_h_mix_f(spf, …)` を呼ぶ。この呼出しに `thermoFloat` の分岐はない。[convectiveFlux_slau_d.inc.cuh:285](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:285)、[同:396](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:396)
- double ビルド用オーバーロードも **Y と T を float に落として float 本体を呼ぶ**。戻り値も float。ブリーフ記載のコミット `4687c3c` にも同じ実装がある。[thermo_d.cuh:959](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:959)
- この面評価は、セル温度反転とは別の高速化として導入されている。[performance-3d-node-sst-speedup.md:74](/home/sano/work/forge/plans/accepted/performance-3d-node-sst-speedup.md:74)
- 約 227.7 K で float32 の温度刻みは **1.526e−5 K、相対 6.70e−8**。記録された変動と近い桁だが、これは原因の証明ではない。

反証条件: **実バイナリで当該経路の切替を確認し、過渡が落ち着いた後も自由流域の流れ・化学種残差と局所振幅が対照の 0.8–1.2 倍に残るなら、「面エンタルピー量子化が主因」を棄却する。**

第 2・第 3 仮説: 現段階では追加しない。幾何の閉包誤差説は、下記の実測を先に照合する。

判別 A/B:

- **変更する一点:** A は現行 `thermo_h_mix_f`、B は同じ面状態 Y・T を既存の double 係数・`thermo_h_mix` で評価する。戻り値だけの double 化では不可。面状態の作り方、再構成、化学種拡散、陰解法は固定する。double 関数は既存。[thermo_d.cuh:241](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:241)
- **共通条件:** `case/60.flatplate_d7275_m7/run_0018_t26_tf0_ext/` の同じ最終場から、新規の二つの run へ `restart_field.py --keep-src-dtype` で保存量をビット一致コピー。メッシュ・BC・DB・CFL・緩和・内反復を固定する。これは診断用変更であり、既定値変更の提案ではない。
- **長さ:** 各 5,000 step。3–4k／4–5k の**領域別指標**の中央値が 5% 超動けば、各 10,000 step まで一度延長。それでも動けば保留。
- **見る量:** 全残差列、自由流域と第一内部列の成分別 √Σres² の時系列、固定節点の高精度 ρ・P・k・ω、比較点 St・位置 II・パネル平均熱流束。局所振幅は末尾で連続反復を採り、100 step 間隔の標本だけで振動消失を判定しない。
- **事前の読み:** **B が A の 1/10 以下となり局所振幅も減る → 第1仮説を支持。0.8–1.2 倍に残る → 主因説を棄却。** 成分限定の改善はその成分に限定し、中間値・過渡継続は保留。G15 の合格判定とは分ける。

やらない方がよいこと:

**Major — 閉包を測らずにメッシュを作り直すこと。** 今回、次を実行した。

```text
python3 solver_density_cuda/tools/check_dual_closure.py \
  case/60.flatplate_d7275_m7/mesh/fp_d7275_y3_wall_top8.h5

CV 394464
閉性 |ΣS|/Σ|S|: median 0、p99 0、max 2.23e−10
VERDICT: PASS
```

さらにブリーフと同じ **y=0.03–0.79 m の 131,488 CV** を抽出すると、**ΣS≠0 の CV は 0、最大 |ΣS| も 0**。既定閾値で PASS だったというだけではない。

ソルバは保存された `PLANES/surfVect` を直接読む。[mesh.cpp:258](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:258)  
したがって、座標・体積の dtype だけでは閉包由来の残差を説明できない。将来精度を検証する場合も、単なる double への型変換では失われた情報は戻らない。変換器の出力型は `geom_float` に従うが、座標読込には `stof` が残る。[gmshReader.hpp:504](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:504)、[同:2394](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:2394)

呼び出し側の前提への異議:

- **Major — 「float 熱力学経路を除外済み」は広すぎる。** セル側の切替だけを試している。対案は採否表の限定表現と上記 A/B。
- **Major — 最終一枚の領域別比と全域履歴の横ばいを合わせても、領域別の定常性は保証されない。** `tf_ab.py` は履歴には窓を使うが、領域別集計には最終場だけを使う。対案は既存スナップショットの領域別時系列化。[tf_ab.py:38](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/tf_ab.py:38)、[同:55](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/tf_ab.py:55)
- **Major — 8 回の変更で熱流束がほぼ動かなかったことを、反復誤差の上限や「既知の数値床」の証明には使えない。** 記録上の判定は両腕とも **`NOT CONVERGED (stalled/plateau)`**。[acceptance.json:219](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:219)

  **追跡はまず上記 A/B 一つで区切る。** 支持された場合は、その経路への感度と残る k/ω の過渡を分けて整理する。棄却・保留の場合は、追加の手当てを連鎖させず、**「G15 未達の探索結果として他項目へ進む」提案をユーザに出す**のが妥当。「比較量に影響しない既知の床」として条件付き合格にする段階ではない。

  G15 の変更を提案するなら、対象量ごとの数値誤差予算、St 等の準定常 VERDICT、壁解像・格子感度、試験 28 の扱いを明記する必要がある。現行 G15 は試験 26・28 と G14 同等の前提ゲートを要求する。[plan:1030](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1030)

不足情報: **AWS の run 原本、実効 YAML、バイナリ差分、残差・場・VERDICT 原本はローカルにない。** A/B の数値は記録に基づく評価であり、独立に再集計できていない。閉包を実測したローカルメッシュの SHA-256 は以下。AWS 入力との一致を確認すること。

```text
4cfe615a7da97b45883fd46db6646cf4060297c127db6ced4dff481726a48eb0
```

ファイル変更・`forge` 起動なし。**plan 未反映**。反映先は呼び出し側の plan §4.12・§5.1 #61・§6 G15。
