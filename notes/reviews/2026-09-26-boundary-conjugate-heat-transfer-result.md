# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `6b7c6f4a`)
- **date**: 2026-09-26
- **commit**: `55d93fa7` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 8.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m2
- **focus**: 8 巡目。7 巡目以降の変更: V2 case/48 の判定と (ii) 1 点の例外受け入れ (ユーザ決定)、V6′ (case/58) の実装・判定 (§5.1 #86–#103、§6 V6′ と閉じる条件)。閉じる条件 3 つが揃ったとして status: done → accepted/ に移してよいかを判定すること。V6′ は範囲限定の PASS (§6 閉じる条件 2 の 1 文)。solver 変更は conjugate.node_log (診断出力のみ) と check_cht_interface.py --band-y、評価器 eval_v6p.py / sens_v6p.py
- **extra**: `case/58.conjugate_slot/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

V2 の数値と、V6′ の固定深部帯に限定した PASS は再現できた。V2 は「FAIL の1点をユーザ判断で例外受け入れ」とする記録を維持する。`accepted/` への移動前に、以下の評価器と文書の修正が必要。

1. **Major — V6′ 評価器が後壁の熱流束の符号誤りを見逃す。**

   [eval_v6p.py:182](/home/sano/work/forge-cht/case/58.conjugate_slot/eval_v6p.py:182) は両壁の `iface_q_eff` を絶対値にしてから、(b)〜(d) を判定する。基準 `run_0012_v6p_df5_i50_nlog` の後壁流束だけをメモリ上で反転したところ、**(c) 0.2072％、(d) 0.1992％、`VERDICT: PASS` のまま**だった。

   **対案:** 流体から壁へ正という規約を保ち、前壁は `q_eff`、後壁は `-q_eff` と方向を明示して比較する。流体柱の収支は符号付きの和で検査し、後壁符号反転の負例を追加する。実データでは前壁が正、後壁が負であり、今回の結果を否定する指摘ではない。

2. **Major — 帯内 G-if が座標欠損・節点間の時刻不一致を誤合格させる。**

   [check_cht_interface.py:71](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:71) は節点表の有限性を検査せず、[同:94](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:94) は節点IDを検証せず配列へ代入し、時刻はブロック先頭行だけを採る。

   読み取りデータをメモリ上で変更して再現した：

   - 帯上端の節点212の座標を `NaN`、残差を `1000 W/m` にすると、その節点を除いた**210節点で PASS**。
   - 最終更新の節点212だけ `step=0` にしても **PASS**。

   **対案:** 節点表の有限性・正の面積、IDの整数性・一意性・完全性、更新内の全行の時刻一致を必須にする。固定帯の節点集合も照合し、不備は `REFUSED` とする。準定常系列を作る `eval_v6p.py` にも同じ検査を共用する。今回の4 run の実データは、これらの整合を独立確認できた。

3. **Major — 受け渡し用の誤差保証が検証範囲を超えている。**

   [plan §4.11(a):358](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:358) の「SU2 との差 0.0001％」は、SU2―解析解の値との混同。§5.1 #81 と [case/52 README:80](/home/sano/work/forge-cht/case/52.conjugate_slab/README.md:80) にある forge―SU2 差は、ソルバ内 **−0.02471 K ≈ −0.152％ of rise**、外部 **＋0.00409 K ≈ ＋0.025％**である。

   また、[plan:364](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:364) の「連成由来の不確かさは1 K未満、caseを変えても保持」は、有限個の検証からは保証できない。V6′ 自身も深部帯限定で、全域の収束を保証していない。

   **対案:** 比較対象を訂正し、誤差は検証したケース・設定・領域での実測値として記載する。他ケースへの一律の1 K保証を削除し、受け渡し文と #85 も同期する。

4. **Minor — SU2 手順書に撤回済みの収支判定が残る。**

   [su2-cross-check.md:189](/home/sano/work/forge-cht/procedures/su2-cross-check.md:189) は `HF[0]` と `HF[1]` を収支として比較するよう指示している。一方、今回の #91/#93 は両者を異なる定義とし、forge の適用荷重との比較対象を `HF[1]` に訂正している。

   実際、`su2_cht_cont` 最終値は **59725.64 / 60684.57 W/m**で約1.58％違う。この差をそのまま連成不釣合いと扱うと、撤回した診断を繰り返す。

   **対案:** 診断流束と適用荷重を区別し、固体側 `CFL_NUMBER=1000` の検証結果、温度を固体 `restart_1.csv` から採る手順まで反映する。

5. **Minor — 完了時の索引・現在仕様の同期が不足している。**

   [plans/README.md:30](/home/sano/work/forge-cht/plans/README.md:30) は V2・V3・V6 を未完了とし、旧 case/51 適用を記載している。[methods/boundary.md:490](/home/sano/work/forge-cht/methods/boundary.md:490) には `node_log` と帯内判定の説明がない。[case/48 README:45](/home/sano/work/forge-cht/case/48.flat_plate_cooled_m4/README.md:45) にも最終SU2比較と例外受け入れが未反映。

   **対案:** 現行仕様・限定付き承認・残作業の移管先を同期する。SU2 成果物は今回のワークツリーにはなく、実際には `/home/sano/work/forge/case/48.flat_plate_cooled_m4/su2_cht_cont/` にあったため、保存先も台帳に明示する。

**再検証の根拠**

V6′ は、`case/58.conjugate_slot/` 配下の次の4 run を固定211節点で再評価した。

| run | (e) 最大誤差 | 帯内 G-if | 準定常 |
|---|---:|---|---|
| `run_0012_v6p_df5_i50_nlog` | 0.0011％ | PASS | ALL STEADY |
| `run_0013_v6p_df20_i50_ext394k` | 0.0025％ | PASS | ALL STEADY |
| `run_0013_v6p_df5_i200_ext394k` | 0.0015％ | PASS | ALL STEADY |
| `run_0014_v6p_df20_i200_ext688k` | 0.0018％ | PASS | ALL STEADY |

全4条件で (a)〜(f) と感度比較は `PASS`。準定常は登録条件の `tail=0.5, drift=0.001, osc=0.001` で各424系列が `STEADY`。**流体残差は再実行でも `NOT CONVERGED (stalled/plateau)`、全域 G-if も `NOT CONVERGED`**。成果物と恒久索引は [case/58 README](/home/sano/work/forge-cht/case/58.conjugate_slot/README.md:43) を参照。

V2 は forge `run_0036_v2_fem2d_nl16` と上記SU2成果物から、平均 **0.500901％**、局所最大 **2.658843％・超過1点**、積分熱量 **0.572809％**を再現した。

**推奨:** 指摘1→2の修正と負例試験、3→5の文書同期、既存データの再評価を済ませてから、限定付きの成果として `accepted/` へ移す。今回の証拠から V6′ の追加CFD計算を閉じる条件に増やす必要はない。ファイル変更なし・plan 未反映。

指摘数: Critical 0 / Major 3 / Minor 2
