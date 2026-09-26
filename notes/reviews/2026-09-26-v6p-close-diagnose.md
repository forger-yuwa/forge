# codex 諮問 (diagnose): v6p-close

- **brief**: [`notes/reviews/briefs/2026-09-26-v6p-close.md`](../../notes/reviews/briefs/2026-09-26-v6p-close.md)
- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **date**: 2026-09-26
- **commit**: `3550f7d5` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **現行説明を同期して、登録した深部帯に限定した V6′ `PASS` を記録し、codex result 8巡目へ進む。**
- **extra**: `case/58.conjugate_slot/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 判断・指摘 | 採否 | 理由・対案 |
|---|---|---|
| V6′ を登録条件どおり `PASS` とする | **採用** | 固定211節点について、4条件とも帯の適格性・(a)〜(f)・局所準定常・帯内 G-if・設定感度を満たす。評価を再実行して確認した。 |
| result 8巡目へ進む | **採用** | V2 の例外受け入れは [plan:744](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:744)、放射の保留は同:754に記録済み。**`accepted/` への移動は result レビューの承認後**。 |
| **Minor：完了前に現行説明を同期する** | **採用** | [plan:706](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:706) は旧 `(e)`、同:708は診断出力を未実装と記載。[README:63](/home/sano/work/forge-cht/case/58.conjugate_slot/README.md:63) も完了済み作業を「次」としている。履歴は残し、現行式 `$Q_{\rm sol}=(K_su-b_s)_{\rm iface}$` と完了結果への参照を付ける。 |

結論: **現行説明を同期して、登録した深部帯に限定した V6′ `PASS` を記録し、codex result 8巡目へ進む。**

記録する主張の範囲は次の1文とする：

> V6′ は、全域 FP64・本ケースの固定深部帯211節点において、登録した伝導漸近、連成保存性・温度連続、局所準定常、帯内 G-if、および `Df_scale` 5/20 × `interval` 50/200 の感度条件を満たして PASS とするが、流体残差と全域 G-if は NOT CONVERGED であり、全域 CHT 解の収束や固定点の設定非依存を実証したものではない。

第 1 仮説: **今回の `20/200` 延長結果は登録条件 A に該当し、V6′ の未完だった感度検証を満たした。** 確度: **高**

根拠：

以下のパスはすべて `case/58.conjugate_slot/` 配下。固定帯は **−19.779487〜−4.166299 mm、211節点、15.613188 W**。

| run | (e) 最大誤差 [%] | 準定常 VERDICT | 帯内 G-if VERDICT | 基準との差 T / q [%] |
|---|---:|---|---|---:|
| `run_0012_v6p_df5_i50_nlog/` | 0.0011 | `ALL STEADY` | `PASS` | 基準 |
| `run_0013_v6p_df20_i50_ext394k/` | 0.0025 | `ALL STEADY` | `PASS` | 0.0000 / 0.0030 |
| `run_0013_v6p_df5_i200_ext394k/` | 0.0015 | `ALL STEADY` | `PASS` | 0.0000 / 0.0017 |
| `run_0014_v6p_df20_i200_ext688k/` | 0.0018 | `ALL STEADY` | `PASS` | 0.0004 / 0.0017 |

- `eval_v6p.py` の評価部分を副産物の書き込みなしで再実行し、全4条件で **`VERDICT: PASS`**。`sens_v6p.py` も **`VERDICT: PASS`**。
- `check_quasisteady.py` を登録どおり `--tail 0.5 --drift 0.001 --osc 0.001` で再実行し、各runの**帯平均2量＋局所422量**すべて **`OVERALL: ALL STEADY`**。`run_0014` の判定区間は累積 **541,000–687,800 step**。
- `run_0014` の帯内 G-if は累積 **672,000–687,800 step の80更新**で **`VERDICT: PASS`**。① **0.19605 W/m²**、② **2.4386e−5**、③ **5.4845e−6 K**、④ **1.3970e−9 W/m**。節点ログの更新ごとの321節点、非有限値なし、②③の値を独立に確認した。
- 全4条件の流体残差は再実行でも **`NOT CONVERGED (stalled/plateau)`**。保存済みの全域 G-if も **`VERDICT: NOT CONVERGED`**。この制限は残る。
- 成果物は各runの `v6p_band.json`、`QUASISTEADY_band.txt`、`CHT_INTERFACE_BAND_VERDICT.txt`、`CONVERGENCE_VERDICT.txt`。恒久索引は [case README の計算run一覧](/home/sano/work/forge-cht/case/58.conjugate_slot/README.md:43)。

反証条件: 固定帯・登録閾値・指定区間で、いずれかの局所量が非 `STEADY`、帯内 G-if が不合格、または感度差が0.5 %を超えること。今回の再評価では該当しない。

第 2・第 3 仮説: 追加なし。緩和遅れとは整合するが、更新数/`Df_scale` の相似則や固定点の一意性まで確定しない。

判別 A/B: **#103で登録した「計算長だけ＋294,000 step」の試験を今回の判別として採用する。** A＝静定・釣合いを満たし基準差≤0.5 %なら登録範囲の感度合格、B＝静定・釣合いを満たして基準差>0.5 %なら「過渡だけ」を棄却。**観測はAであり、追加計算は閉じる条件として要求しない。**

やらない方がよいこと: 全域未収束を省略すること、V2 の既存 `FAIL` を `PASS` に書き換えること、放射保留に伴う受け渡し注記を落とすこと。

呼び出し側の前提への異議: 数値判定への異議は無し。ただし今回の診断承認は result 8巡目の承認を代替しない。「0.0000 %」も丸め表示であり厳密一致ではない。

不足情報: 今回の限定された V6′ 判定を妨げる不足は無し。**ファイル変更なし・plan 未反映**。反映先は呼び出し側の plan §5.1 #103、§6 V6′・レビュー記録、および case README。
