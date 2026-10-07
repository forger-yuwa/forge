# codex 諮問 (diagnose): euler-t0-result-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-07-euler-t0-result-interpretation.md`](../../notes/reviews/briefs/2026-10-07-euler-t0-result-interpretation.md)
- **plan**: [`plans/active/verification-case45-euler-total-enthalpy.md`](../../plans/active/verification-case45-euler-total-enthalpy.md)
- **date**: 2026-10-07
- **commit**: `628a4c91` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.9 min, rc=0
- **結論**: 追加 CFD は行わず、上記の限定と採用根拠を §9 に反映した版を正式な result レビューに回す。
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical なし）

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| ①「G1 格子の soft 段の計算の中で発生」 | **要再検証・Major** | **新配点 B でも soft 終了時に −329.136〜+346.184 K**。G1 の A は −354.816〜+315.382 Kで、両者とも極値は下流軸上。一方、A の本段 54000 ではスロート壁際に −234.926〜+292.609 K がある。[E2 評価記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/euler_t0_e2_eval.json:25897)。**起動時の下流軸上偏差と、後期の G1 壁際偏差を別の観測として書くこと。** B→C の最大 1 ULP は本段起動が数百 K 偏差の発生源でないことを支持するが、soft 起動直後の初期化と反復計算は区別していない。 |
| ②「run_0062 は影響なし」 | **却下・Major** | 確認できたのは、凍結源 `res_6000` の偏差が −0.058〜+0.073 Kであることと、初期線の指定量の変動幅。全温・質量流束の時系列と離散化誤差は未評価と明記されている。[plan:191](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:191)。**「使用する凍結スナップショットには数百 K の異常を認めず、凍結源は維持する。履歴・離散化誤差の影響は未評価」**に限定する。 |
| ②「較正差 −3.70e−4 は格子の違い全体の結果」 | **要再検証・Major** | 全温異常だけに帰属させない判断は採用。ただし旧較正 `run_0113/0114` と E4V は、単調壁への変更・初期化・評価期間も異なる。[run 索引:92](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:92)、[plan:107](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:107)。**「新しい評価条件で再較正した結果の差。格子単独・全温異常単独の寄与は未同定」**とする。 |
| ③「二つのゲートを分ける」「run_0164・0174 は両方満たす」 | **分離は採用、0174 の断定は要再検証・Major** | `e4_recal_eval.json` の対象は `run_0163/0164` のみ。E5 は `run_0174` を含め PASS だが、出口 M・全温・残差の条件を検証するものではない。**0174 自身の条件別記録を添えるまで「両方合格」としない。** V5d 全体の保留から、逆に 0174 単独の不合格とも断定できない。[E4V 採用先](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/e4_recal_eval.json:51954)、[V5d の記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:112)。 |
| ①の「B は整定」、②の「V5/V5b を V5d で置き換え」 | **表現修正・Minor** | E2-B は本段区間 `NOT CONVERGED (stalled/plateau)`、13 枚の準定常判定は `DRIFTING 25 列・OSCILLATING 2 列`。**「登録窓の絶対変動幅条件を満たした」**と書く。V5d は `保留（前提不成立）` なので、**「再評価先を V5d に移したが採否は保留」**とする。 |

確認できた採用根拠は次のとおりです。以下の run はすべて `case/45.isobutane_m6_d155/` 配下です。

- **E4V：採用を支持。** `run_0164_euler_e4_recal_d1` は本段の収束判定が `NOT CONVERGED (stalled/plateau)`。出口 M は 42000〜54000 と 50000〜54000 の両窓で `STEADY`、最大目標偏差 3.0574e−6。全温偏差は最大 0.13413 K、27 列の時間幅の最大はそれぞれ 0.06893 K・0.04810 K。登録された実務較正条件を満たしている。
- **E5：限定された感度試験として PASS を支持。** `run_0165/0166/0167` の NS 最終場を固定したとき、参照13枚による幅は 3.7518e−5・1.5851e−5・4.0780e−5。JSON 内の系列から再計算し、全条件で閾値 5e−4 以下と最終値の照合を確認した。ただし、これは参照選択による変動幅の検査であり、NS の収束や旧参照からの偏りを保証しない。[E5 評価記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/euler_ref_de_sensitivity.json)

**現状のまま `done`・accepted 移行は不可です。** 上記の解釈修正、`run_0174` の証拠提示または未確認への訂正、正式な result レビューの記録が必要です。原因候補と未確定範囲を残して完了すること自体は §8 に適合し、倍精度ビルド等による根因確定を追加の完了条件にする必要はありません。今回の diagnose は result レビューの代替ではありません。

結論: 追加 CFD は行わず、上記の限定と採用根拠を §9 に反映した版を正式な result レビューに回す。

第 1 仮説: soft 終了時の数百 K 偏差には両配点で生じる起動過渡が含まれ、G1 の後期壁際異常とは分けて扱う必要がある。確度: 中  
　根拠: E2 の両腕で soft 終了時に下流軸上の数百 K 偏差があり、B は本段 42000〜54000 の全保存時点で最大絶対偏差 ≤0.103227 K。一方、A は同窓でもスロートの最大絶対偏差が全時点で >182 K。  
　反証条件: soft の反復前から数百 K 偏差が存在すれば、「反復中に形成する起動過渡」という発生時点の解釈は誤り。

第 2 仮説: G1 の後期異常には配点による離散化・幾何精度・反復挙動のいずれかが関与する。確度: 中。ただし E2 の正式判定は `判別不能`、各機構の分離は未確認。

判別 A/B: ①の「計算の中」という断定を残す場合に限り、同一 run の **A＝soft 起動直後 `res_0`、B＝soft 最終 `res_3000`** を同じ保存量復元で比較する。変更点は取得時点のみ、追加 CFD 0 step。A が全節点 ±1 K 以内で B が数百 Kなら起動直後の書換え説を退ける／A に既に >100 Kあれば反復中に初発した説を退ける。中間値・保存物欠落なら限定できない。

やらない方がよいこと: E2 を事後的に合格へ変えること、soft の偏差から G1 の後期異常の原因を断定すること、E5 PASS を `run_0174` の出口較正条件の証明に使うこと。

呼び出し側の前提への異議: 「整定しない収束」は観測された症状であり、離散化・精度問題と排他的な原因候補ではない。単調壁 E′ の維持と NS 許容帯 ±0.05 % はユーザ決定として維持し、その判断を数値的な原因特定とは分ける。

不足情報: soft 起動直後の状態、`run_0174` の本段全残差・両評価窓の出口 M・全温条件の個別記録。ファイル変更・forge 起動は行っていない。plan 未反映。
