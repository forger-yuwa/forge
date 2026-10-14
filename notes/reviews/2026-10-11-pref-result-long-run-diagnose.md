# codex 諮問 (diagnose): pref-result-long-run

- **brief**: [`notes/reviews/briefs/2026-10-11-pref-result-long-run.md`](../../notes/reviews/briefs/2026-10-11-pref-result-long-run.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `47a6bec6` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.4 min, rc=0
- **結論**: **短期試験の再試行は挟まず、無介入 float の陽性対照1本を加えた計4本で、`space.pRef` だけを変える30,000 stepの有限期間比較を行う。**
- **extra**: `case/45.isobutane_m6_d155/tr_an.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/rpref_judge.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| §6.45 を「判別不能」として残す | **採用** | 保存 JSON の `a` から e_Fp＝0.0920105、d＝0.8205338、SE_d＝0.0550075 を再計算し一致した。L_d＝0.4939798 を丸めて合格にしない。[rpref_judge.json:46](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rpref_judge.json:46) |
| 「型と直交する差が RMS 約1.5 K 残る」 | **却下・Major** | 計算しているのは、正規化した温度変化の射影係数を3時点で平均した値。step 2,000 の絶対温度差の RMS を直交成分と同一視できない。[rpref_an.py:65](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rpref_an.py:65)。**「温度差の RMS は1.5〜1.6 K残る。直交成分の大きさは未測定」**へ直す。 |
| 新規3本と旧 A・B・B′だけで長期効果を判定する | **要再検証・Major** | 旧 A は `fgeom7-fp64`、短期 Cp は `fgeom11-fp64`。短期の同等性は異なるバイナリの30,000 stepの同等性を保証しない。また、別試験の無介入 F は W1・θ_r(40) で旧 B・B′平均より約0.183ポイント低かった。[plan:716](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:716)、[plan:976](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:976)、[plan:877](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:877)。**旧バイナリそのものを使い、同じ実行環境で float・pRefなしの陽性対照を1本追加する。** |
| 「2本とも満たす」を再実行誤差への保証にする | **却下・Major** | 2本の一致は不確かさの上限ではない。旧試験でも再実行差は介入によって変わった。[plan:873](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:873)。10N条件を外すなら、**新しい有限期間の診断基準として事前登録し、統計的同等性・一般的な再現性の証明とは呼ばない。**過去の判定は変更しない。 |
| 短期試験を5本へ増やす操作を先に挟まない | **採用** | 上流の射影に関する精度を上げても、下流 θ_r・Q_w への効果は未測定のまま。目的量を直接測る方を優先する。 |

結論: **短期試験の再試行は挟まず、無介入 float の陽性対照1本を加えた計4本で、`space.pRef` だけを変える30,000 stepの有限期間比較を行う。**

第 1 仮説: `pRef` の変更だけで、登録した下流 θ_r の有限期間の差を参照尺度 D の0.1倍以内へ抑えられる。 **確度: 低**

　根拠: 短期試験では射影の点推定が0.913から0.092へ低下し、FP64側の関門も保存 JSON 上は通過している。ただし正式判定は「判別不能」で、下流への効果は測っていない。`case/45.isobutane_m6_d155/run_0700_rp_cp1`〜`run_0705_rp_fp3` の保存された収束判定は、全6本とも **`OVERALL: CHECK FAILURES ABOVE`**。定常解の改善を示す結果ではない。

　反証条件: 以下の前提を満たした上で、pRefありの float 2本とも、両窓・全3断面で Ap から D の向きへ0.5|D|以上離れる。

第 2・第 3 仮説: 第2は、上流の射影の減少だけでは下流の目標を満たさず、他の演算・初期化との差の相互作用が残る（未確認）。第3は、旧対照との実行条件差が見かけの改善を作る（未確認。今回の陽性対照で点検する）。

判別 A/B: **変更は `space.pRef: 0 → 5499000.0` の1点。**

- 新規は **F0＝float・pRefなし×1、Bp＝float・pRefあり×2、Ap＝FP64・pRefあり×1**。全て共通の `_tr/q32_init.h5` から開始し、30,000 step固定、step 0と500ごとの出力を残す。同じ機械で実行順・同時実行条件を事前固定する。
- バイナリは旧 A・B・B′で使った `fgeom7` の実ファイルを使用し、SHA256を照合する。各精度内で初期化後の9保存量、入力・幾何・物性・BC・環境変数を照合する。短期 Cp のバイナリをそのまま流用しない。
- Δ＝100×〔各runの窓平均−そのrunのstep 0〕／|旧Aのstep 0|。W1＝20,000〜25,000、W2＝25,000〜30,000、D＝(1.892, 1.647, 1.431) %を固定する。独立な反復はrunであり、窓内の11点や6組を反復数に数えない。
- **前提**：全6組で |ΔAp−ΔA|≤0.1|D|。さらに F0 が旧 B・B′の範囲を N_B だけ両側へ広げた区間に入り、sign(D)×(ΔF0−ΔAp)≥0.5|D|。欠ければ今回の比較は判別不能。
- **結果A**：全6組で、Bpの2本がそれぞれ |ΔBp−ΔAp|≤0.1|D|。→「今回の2本・登録窓では目標を満たした」とし、第2仮説の「この条件でも不十分」を退ける。
- **結果B**：全6組で、Bpの2本がそれぞれ sign(D)×(ΔBp−ΔAp)≥0.5|D|。→第1仮説の「この変更だけで十分」を退ける。
- その他は判別不能。各組の残存差と、無介入との差も数値で残す。**0.1|D|以内は「今回観測した差を90%削減」と同義ではない**。Dは別の停止時点に由来する尺度である。
- Q_w は独立した副判定とする。θ_rがAでもQ_wの改善は主張しない。全期間の有限性、`check_convergence`、対象系列の `check_quasisteady` のVERDICTを残し、未収束なら有限期間の応答として報告する。

やらない方がよいこと: L_dの僅差を理由に短期試験を合格扱いすること、結果を見ながら反復を追加して従来の区間を再判定すること、長期試験の陽性を「残差の加算が真因」「doubleを廃止できる」へ広げること。

呼び出し側の前提への異議: **「前提を全て検証済み」はコードだけからは確認できない。** `rpref_an.py:50` は存在する保存量だけを返し、`:55` もその集合だけを比較するため、9量の欠損を確実には拒否しない。入力ハッシュ・座標対応の証跡も保存 JSON にはない。対案は、期待する9キーの完全一致を先に検査し、入力・座標・設定の照合結果を判定記録に残すこと。**Major**。[rpref_an.py:49](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rpref_an.py:49)

不足情報: 元の `run_*`、HDF5、残差履歴、長期試験の `tr_judge.json` はこの作業ツリーにない。今回直接確認できたのは判定コード・短期試験の保存JSON・plan記録であり、元データの有限性や長期VERDICTは独立再検証していない。ファイル変更・forge起動は行っていない。**plan未反映**。呼び出し側で§6.45の表現修正と次試験の事前登録へ反映する。
