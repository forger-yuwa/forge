# codex 諮問 (diagnose): graded-ic-next

- **brief**: [`notes/reviews/briefs/2026-09-27-graded-ic-next.md`](../../notes/reviews/briefs/2026-09-27-graded-ic-next.md)
- **plan**: [`plans/active/axisymmetric-graded-grid-static-gas.md`](../../plans/active/axisymmetric-graded-grid-static-gas.md)
- **date**: 2026-09-27
- **commit**: `a1dbec00` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 5.0 min, rc=0
- **結論**: **純伝導 IC を連成計算の回避策として採用する前に、初期場だけを変える非連成 A/B を各 20000 step まで再登録し、静止保持を判定する。**
- **extra**: `case/62.conjugate_disk/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 重大度 | 諮問事項 | 採否・根拠・対案 |
|---|---|---|
| Major | 純伝導 IC で V-ax2b を直ちに回し直す | **現時点では却下**。`run_0016_ic_{A_uniform,B_conduct}` は再実行した残差判定でも双方 **`NOT CONVERGED`**。B の局所保存場から得た速度・圧力差の準定常判定も **`TRANSIENT-UNSETTLED`（2 snapshots）**。先に非連成の静止保持を判定する。 |
| Major | IC 変更は登録変更に当たらない | **却下**。親 plan の V-ax2b 本文には明記がないが、[case README:22](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md:22) と `run_0016_ic_B_conduct/solverConfig.yaml:4` は一様 IC を明記している。純伝導 IC を採用する場合は**変更理由・保証範囲を追記して再登録**し、既存の失敗を残す。 |
| Major | 非一様格子が圧力の低い偽定常状態に居座った | **要再検証**。`run_0014_disk_r32g_w20k/res_20000.h5` の速度 **0.231109 m/s**、圧力 **728.336–742.394 Pa** は確認したが、単一時点では定着を証明できない。さらに局所擬似時間を使うため、初期質量から計算した **1051.736 Pa** への不一致だけでは保存欠陥を診断できない（[setDT_d.cu:362](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/setDT_d.cu:362)）。本 plan で未解決として追跡し、まず速度・熱流束・圧力不均一の時系列で判断する。 |
| Minor | B は厳密な平衡初期場、step 500 の全速度は `3.3e-4` | **訂正を採用**。B の `res_0.h5` は ΔP **5.26e-5 Pa**、入力 `VALUE/ro,roe` は **float32**。step 500 は **max\|U\| = 4.04865e-4 m/s**、**max\|Uy\| = 3.32140e-4 m/s**。初期場の微小誤差と速度成分を区別して記録する。 |

結論: **純伝導 IC を連成計算の回避策として採用する前に、初期場だけを変える非連成 A/B を各 20000 step まで再登録し、静止保持を判定する。**

第 1 仮説: 大きな起動非平衡が大振幅化を支配し、純伝導 IC では微小な離散誤差への応答に留まる。確度: **中**。後半の長時間保持は未確認。
  
  根拠: `case/62.conjugate_disk/run_0016_ic_A_uniform/` と `run_0016_ic_B_conduct/` は設定が同一、バイナリ SHA256 も同一。入力 HDF5 の差分は `VALUE/ro,roe` のみだった。step 500 の全速度最大値は **39.3098 → 4.04865e-4 m/s**、冷却壁の市松振幅は **3.55916e5 → 3.84905e-3 W/m²**。初期壁ピンが密度を保持して圧力を変更する実装も確認した（[nodeWallDirichlet_d.cu:84](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:84)）。これは起動依存性の証拠であり、定常解同士の比較ではない。

  反証条件: 純伝導側でも、非連成のまま 20000 step 以内に速度・熱流束誤差・圧力不均一が登録閾値を超える。「平衡初期化で静止保持できる」という部分を棄却する。

第 2 仮説: 平衡近傍にも弱い駆動源または遅い不安定モードがあり、500 step は発現前である。確度: **低〜中、未確認**。B の増速だけでは判別できず、運動量残差は実際には低下している。

第 3 仮説: **機構としての slip 欠陥と幾何・離散化誤差の順位は未確定**。現データから EOS 床、リミッタ、陰解法対角のどれかを選ぶ根拠はない。

判別 A/B: **変更点は初期保存量場だけ**。A＝一様 325 K、B＝今回と同じ純伝導 IC。既存入力を複製した新規 run を使い、同じバイナリ・格子・数値設定で各 **20000 step**。両方とも連成を開始させない。`warmup` は終了 step より大きく固定する。コードは `stepAbs == warmup` で連成更新を許すためである（[conjugateWall.cpp:1070](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1070)）。

- **長さ・窓**: 20000 は既存の失敗観測と揃える診断期間であり、収束保証ではない。全期間を保存・監視し、最終比較窓を **10000–20000** に事前固定する。場・両壁の診断量は例えば 100 step ごとに保存する。
- **判定**: 残差は各 run の同一設定区間全体で `check_convergence.py`。準定常は全速度最大値、半径速度最大値、圧力差、市松振幅、**全壁節点の熱流束**の系列を `check_quasisteady.py --series-csv` に渡す。末尾半分、`--drift 0.001 --osc 0.001` を事前登録する。
- **合格条件**: plan §6 の **max\|U\| ≤1e-3 m/s、全壁熱流束誤差 ≤0.5%、相対圧力差 ≤1e-4** と、残差 **PASS**・対象量 **STEADY** を併用する。
- **A が異常を再現し、B が全条件を満たすなら**、この期間での起動依存性を支持し、純伝導 IC に限定した V-ax2b 再登録へ進める。**B も閾値を超えるなら**、初期化だけで解決する案を棄却する。残差の低下桁不足だけ、または閾値内でも `DRIFTING` なら**判別未了**として止める。

やらない方がよいこと: B の小ささを理由に CHT を開始すること、`--drop` を下げて合格させること、未収束の B を `--from-floor` の参照にすること、圧力低下だけから EOS 床や質量漏れを追い始めること。静止保持を通っても、**一様 IC からの失敗は本 plan に未解決として残す**。

呼び出し側の前提への異議: **H2 支持は採用するが、根本原因の特定とは認めない**。既存の閉包補正 A/B が弱めたのは「大振幅の起動異常を閉包補正だけで解消できる」という仮説であり、B の微小な流れに対する幾何誤差の寄与まで除外していない。また、ツールの `needs scheme change` という文言自体は原因診断ではない。

不足情報: AWS にある step 250–500 の集計元系列と、両壁全節点の熱流束系列。ローカルでは登録窓の最大値を全面的には再検証できなかった。**plan 未反映**。呼び出し側で本 plan §4・§5.1・§6 に上記診断を登録し、親 plan の V-ax2b は依存障害として保留を維持する。
