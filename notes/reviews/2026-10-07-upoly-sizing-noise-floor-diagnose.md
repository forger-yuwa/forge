# codex 諮問 (diagnose): upoly-sizing-noise-floor

- **brief**: [`notes/reviews/briefs/2026-10-07-upoly-sizing-noise-floor.md`](../../notes/reviews/briefs/2026-10-07-upoly-sizing-noise-floor.md)
- **plan**: [`plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md`](../../plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md)
- **date**: 2026-10-07
- **commit**: `733d8b1d` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.5 min, rc=0
- **結論**: **許容差・壁モデル・反復法を固定し、`integral_bl` の実効 `rtol` だけを変える CFD 0 step の A/B で、微小な寸法変更に対する不安定さが解消するかを先に判別してください。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 論点 | 判定・重大度 | 根拠と対案 |
|---|---|---|
| CFD 前の許容差を直ちに緩める | **要再検証／Major** | 1e−9 m は逆算と壁生成の閉合条件です（[plan:140](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:140)）。加工公差・float32 の精度から変更する根拠にはなりません。**先に積分精度への依存を測る**べきです。変更する場合、旧 U2b の FAIL は保存し、測定点・誤差指標・許容差の算定規則を先に登録した新試験にしてください。観測した散らばりは、真値からの誤差上限ではありません。 |
| 生産と共通の δ_r 経路 | **採用／Major** | 逆算と壁生成が同じ関数を使う構造は妥当です（[deltastar_loop.py:126](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:126)）。対処候補はまず **`rtol` の厳格化**。固定刻みへの変更は保留。Re^−0.2 を厳密な分離として導入する案は却下します。CONTUR の閉包には N(log Re) と対数型摩擦則があり、単純なべき乗則ではありません（[deltastar_integral.py:149](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:149)、同:165）。 |
| `solve_rt` の例外を警告へ変更／旧経路へ戻す | **却下／Major** | [U2b.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/upoly/U2b.json) は3目標とも30回で FAIL、最終残差は +7.210e−7／−5.488e−7／+9.529e−7 m。機能回復は必要ですが、成功扱いや経路不一致の復活では直りません。**例外と共通経路を保持して下記 A/B を先行**してください。 |
| U2 が共通経路の優位を実証した | **その因果解釈は却下／Major** | A は逆算前の壁ゲートで停止しています。さらに A は `cf_scale` を渡さず既定1、B は1.054129…で、plan の「k_f 固定」と一致しません（[upoly_verify.py:303](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/upoly_verify.py:303)）。**B の往復検証と、A の入力不適合を別に記録**してください。共通経路の採用理由は維持できますが、「A の逆算寸法だけが生産壁で外れる」は未検証です。 |

結論: **許容差・壁モデル・反復法を固定し、`integral_bl` の実効 `rtol` だけを変える CFD 0 step の A/B で、微小な寸法変更に対する不安定さが解消するかを先に判別してください。**

第 1 仮説: **積分の適応刻みと閉包・補間の組合せによる数値誤差が、r_t に対する微小変動を増幅している。** 確度: **中**

  根拠: `integral_bl` は RK45、`rtol=1e−6`、`atol=1e−14` を使います（[deltastar_integral.py:210](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:210)）。今回の読み取り専用再計算では、r_t=0.07665396532029795 m の同一入力2回は同値。一方、r_t を +1e−12 m 動かすと、**平滑化前の出口補正が実寸換算 +2.595e−7 m、平滑化後が +2.639e−7 m**変わりました。変動は P-spline より前から存在します。保存 JSON の局所的な散らばりも再集計で確認しました。

  反証条件: `rtol` が実際に変更されたことを確認しても、平滑化前後の変動がほぼ残る場合、**「`rtol` の厳格化だけで解消できる」という説明を棄却**します。RK45 自体が真因と確定したわけではありません。

第 2 仮説: 閉包の打切り・床・補間の折れ目が支配し、外側の積分許容差だけでは不足する。確度: **低・未確認**。候補箇所は [deltastar_integral.py:153](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:153) の打切りと同:165 の Re 床ですが、今回それらが切り替わった証拠はありません。

判別 A/B: **A=`rtol=1e−6`、B=`rtol=1e−10`。変更はこの1点だけ。**

- CFD は0 step。MOC、k_f、熱条件、平滑化、`atol`、`max_step`、逆算許容差は固定する。
- 上記 r_t を中心に、既存 `noise_floor` と同じ11点を両腕で評価する。同一入力の再評価も行う。測るのは生の δ_r、平滑化後の δ_r、物理スロート・出口半径の散らばり。
- 出口目標0.775 mについて、両腕とも固定点反復は最大30回。返った寸法から、**同じ腕の精度設定で**生産壁を作り直す。
- **実装上の注意:** 現在の [runner_axismach.py:1028](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1028) は `rtol` を転送しません。YAML に追加するだけでは無効です。診断用の引数注入などで、実効値を確認・記録してください。

  → **B で両半径の局所散らばり ≤1e−10 m、かつ出口往復誤差 ≤1e−9 mになれば**、現行積分精度が障害だったという仮説を支持します。  
  → **B でもいずれかを満たさなければ**、この精度変更だけで復旧できる案を棄却します。改善途中の結果を PASS にしません。

この基準は新しい診断用の事前登録です。B を本採用するには、±5%目標での往復、さらに厳しい精度での値の安定性、全壁の半径・傾き・曲率の差と既存形状ゲートを確認する必要があります。局所散らばりが小さくても、積分の系統誤差や NS 性能の不変は保証されません。

やらない方がよいこと: **1e−7／1e−5 mへの即時緩和、反復回数を増やして偶然の閾値通過を待つこと、ゲートを外して A を完走させること。** また、`rtol` 変更と Re スケーリング導入を同時に行わないでください。

呼び出し側の前提への異議:

- **「U2 は判定として成り立たない」は広すぎます。** [U2.json:13](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/upoly/U2.json:13) の基準目標は初期値から作られ、1回・残差0で通っています。±5%は6回・7回です。登録済みの往復条件を満たした事実は残し、**微小摂動への頑健性は保証しない**と限定すべきです。
- **「常に例外」は未立証。** 確認できた事実は、この問題の登録3目標で失敗したことです。
- **「雑音の床」は確率的な誤差下限ではありません。** 現状の測定は11個の異なる入力に対する局所変動幅です。同一入力の再現性とは区別してください。
- 今回は幾何・逆算の診断です。保存された「NS 後経路 PASS」を、CFD の収束や派生量の定常性の根拠には用いていません。

不足情報: 実効積分精度を変えた比較、精度を上げた参照解との差、新しい寸法許容差に割り当てられる誤差予算。加工公差は未確認です。確認時の HEAD は `733d8b1d` で、ブリーフの `2c09d681` から対象コード4ファイルに差分はありません。

**plan 未反映**。編集禁止のため、呼び出し側で対象 plan §4.2・§6・§9 に診断基準と採否を記録してください。
