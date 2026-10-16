# codex 諮問 (diagnose): m55-p1-axis-wave

- **brief**: [`notes/reviews/briefs/2026-10-15-m55-p1-axis-wave.md`](../../notes/reviews/briefs/2026-10-15-m55-p1-axis-wave.md)
- **plan**: [`plans/active/tooling-nozzle-card-local-mesh.md`](../../plans/active/tooling-nozzle-card-local-mesh.md)
- **date**: 2026-10-15
- **commit**: `36af1655` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 9.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-m55-p1-axis-wave-diagnose.codex.log` (codex: ok)
- **結論**: pass 2 は保留し、pass 0 の保存場・Euler 参照・抽出帯を固定して、抽出の径方向細分 `refine=4/32` だけを変える、forge 0 step の A/B を先に行う。
- **extra**: `case/71.m55_combustion_wt/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **pass 2 は保留し、pass 0 の保存場・Euler 参照・抽出帯を固定して、抽出の径方向細分 `refine=4/32` だけを変える、forge 0 step の A/B を先に行う。**

採否表

| 対象 | 判定・重大度 | 根拠と対案 |
|---|---|---|
| H1「補正壁の波が軸へ伝わった」 | **要再検証・Major** | 壁形状起因は有力だが、「梯子 → 壁角の窪み → 軸の波」の因果は未実証。後述の実測では、抽出の階段と径方向の採用点の切り替わりが対応している。まず抽出の数値感度を調べる。 |
| H2「コア不一致が core-matched 抽出へ直接混入」 | **この説明は却下・Major** | 生産経路は `method="band"`。参照は壁際の帯から作り、内側コアの RMS は診断値である。[deltastar.py:352](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:352)、[同:814](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:814)。帯まで波が及ぶ可能性は残るが、`core_rms` だけで区間の重みを 0 にする案 (e) は採らない。 |
| pass 2 の補正方法 | **(d) の前提診断を採用・Major** | 現段階では (a) 全量更新、(b) 緩和だけ、(c) 平滑化だけ、のいずれも推奨しない。特にカード経路は **`delta_target` でなく `delta_E` を読む**ため、関門の ω を 0.5 にしても壁更新は半減しない。[deltastar_loop.py:351](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:351)、[同:381](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:381)。緩和を導入するなら、寸法逆算と実際の壁生成が同じ緩和表を使う設計が必要。 |
| pass 1 の成果物を渡す | **途中成果として採用、合格品扱いは却下・Major** | `run_0005_ns/HANDOFF_GATE.json` は `proceed:false`、寸法も `deliverables/DIMENSIONS.json` が `FAIL`。限界を添えた比較資料として渡すことは妥当。**風洞の試験部一様性を優先する暫定候補は pass 0** を推す。保存された軸分布では最大絶対偏差が pass 0 の 0.418 % に対し pass 1 は 0.764 %。出口平均だけで pass 1 を上位としない。どちらも検証済み設計とは呼ばない。 |
| E(i) 流量の二つの呼び手 | **採用・Major** | 自分で再計算して、列のままでは **120.192639 kg/s・散らばり 2.1149 %**、一定 x では **121.350003 kg/s・0.00641 %** を再現した。これは判定緩和ではなく積分断面の訂正。[runner_axismach.py:1232](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:1232)、[deltastar.py:900](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:900) で設定を明示し、一定 x の既存格子の非劣化を確認する。 |
| E(ii) 直管終端の関門 | **単一点の例外案は却下・Major** | 半径の許容差から接合位置を逆算する現行処理には問題があるが、最後の一点だけの一致で置き換えると、他の節点の不整合を見逃す。[dimension_gate.py:123](/home/sano/work/forge-cardmesh/design/forge_design/export/dimension_gate.py:123)、[同:238](/home/sano/work/forge-cardmesh/design/forge_design/export/dimension_gate.py:238)。**接合位置は STEP の曲線から決め、同じメッシュ節点位置で STEP を評価し、直管側から接合を越える区間全体の半径を照合する方式**を推す。カード寸法・弦誤差の条件は維持し、接合をずらした負例も拒否すること。 |

第 1 仮説: **帯のフィット点・積分点の離散的な出入りが δ 抽出の階段を作り、その成分が pass 1 の補正壁へ転写された。** 確度: **中**。階段の数値的起源は強く疑うが、軸の波までの因果は未確定。

  根拠:

- [deltastar.py:381](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:381) は帯内に入った点だけで等重みの最小二乗を解く。[同:418](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:418) は `y <= y_b` の点だけで欠損を積分し、積分端 `y_b` の点を追加していない。滑らかに帯が動いても採用点は離散的に変わる。
- `case/71.m55_combustion_wt/run_0005_ns/res_200000.h5` と `_runs/run_0005_ns/judge/diag_wave/pass1_run_0005_delta_r_equiv.csv` を照合した。x＝1〜3.5 m の **527 区間中 86 区間**でフィット点集合が変わり、`δ_raw−δ_in` の隣接差の絶対値上位 **10 箇所すべて**でフィット点・積分点が切り替わった。
- 採用点が不変の 441 区間では隣接差の最大 **0.00604 mm**、切り替わる区間では最大 **0.26365 mm**。pass 0 の保存 CSV にも x＝2.64858 m で **−0.24737 mm** の隣接段差がある。
- `edge` は縁を連続推定して x 方向にも平滑化する実装である。[deltastar.py:492](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:492)。したがって、階段を直ちに「帯の梯子そのもの」と断定するのは早い。

  反証条件: 固定帯で径方向細分を増やしても、平滑化後の補正量・補正角が後述の許容内で不変なら、**固定帯の径方向標本化が壁角の窪みを作る主要因**という仮説は退ける。梯子による帯選択の感度まで除外したことにはしない。

第 2 仮説: 補正量の長波長成分によって NS の有効壁が設計壁からずれ、軸へ波を送った。確度: 中。壁角差はあるが、到達位置・振幅の対応は未確認。

第 3 仮説: 精度・径方向解像・残存過渡が波に寄与した。確度: 低、未確認。pass 0→1 は壁以外にも精度・格子・起動を変えており、対照実験ではない。

判別 A/B: **変えるものは `band_local_deficit` の `refine` だけ。A＝4、B＝32、forge は両腕とも 0 step。**

- 対象は壁更新の元になった `case/71.m55_combustion_wt/run_0003_ns_cont/res_200000.h5`、固定参照は `run_0001_euler/res_12000.h5`。一定 x への補間、各列の保存済み `band_y_b`、Euler の補間、平滑化のノット・λ・重みを固定する。**帯の梯子は再選択しない。**
- 既存の `_load_structured`・`band_local_deficit(y_b_fixed=…)`・`smooth_delta_quintic` を再利用する。まず A が保存済み raw を再現することを確認する。
- x＝1〜3.5 m で、raw の段差、平滑化後の δ_E、補正角 `atan(r′_design＋δ′_E)` の差を見る。判別用の事前許容を **max|Δδ_E|＝0.01 mm、max|Δθ|＝0.01°** とする。
- → **結果 A：どちらかを超える**なら「現抽出は径方向標本化に十分独立」という前提を棄却し、抽出側を直す。  
  → **結果 B：両方以内**なら、今回疑った固定帯の標本化を壁角の窪みの主要因から外す。  
  これは抽出の検査であり、どちらでも軸の波の原因確定とはしない。

やらない方がよいこと: **物理壁 Euler 一本の波の有無だけで「形状／境界層」を二分すること。** それは排除厚を補う物理壁から境界層を除去する試験で、NS の有効壁と異なる。波が出ても今回の局所窪みが原因とは限らず、出なくても相互作用だけが原因とは言えない。また、抽出の検査前に λ を上げて波を消したり、ω だけを変えて pass 2 を投入したりしない。

呼び出し側の前提への異議:

- **Major：「過渡ではない」は未立証。** 再実行した判定は本段 0〜199,999 が `NOT CONVERGED`、停止規則系列 2,000〜200,000 が `DRIFTING`。出口コア M は `STEADY` だが、軸の波高・極値位置の系列は判定対象にない。100k と 200k の二枚比較で過渡を除外しない。軸の波そのものを `check_quasisteady --series-csv` に渡す。
- **Major：H3 の因果は未分離。** 出口コア M が目標へ近づいた観測は採用するが、全体の補正量だけの効果と確定しない。
- **Minor：「軸の波 ±0.6 %」は不正確。** 最終 HDF5 から再確認した範囲は **−0.76357〜＋0.54743 %**、半振幅は **0.65550 %**。目標からの最大絶対偏差と半振幅を分けて報告する。
- **Major：合格条件は出口平均と別に必要。** 新しい提案値として、試験部 x＝4.0〜7.95 m の最大絶対偏差 ≤0.1 %、半振幅 ≤0.05 %を事前登録する。これはユーザ指定値ではない。既存の **未緩和 max|δ_E−δ_in| ≤0.05 mm**、壁の pass 間変化 ≤0.05 mm、二窓 ≤0.01 mm は維持する。現 pass 1 の固定点差は **0.51666 mm＠x＝7.62967 m**。固定点条件は反復の終了条件であり、次の更新の開始条件と混同しない。

不足情報: 軸の波の末尾時系列とその VERDICT。pass 0・Euler の全場はこの手元コピーにはないため、提案 A/B は保存元で行う必要がある。これらは原因確定に必要だが、今回の「次は抽出の感度を調べる」という判断は妨げない。結果の索引は [case README の計算 run 一覧](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:21)。

証拠: 足りる
