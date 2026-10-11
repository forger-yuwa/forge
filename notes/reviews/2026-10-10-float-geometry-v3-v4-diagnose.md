# codex 諮問 (diagnose): float-geometry-v3-v4

- **brief**: [`notes/reviews/briefs/2026-10-10-float-geometry-v3-v4.md`](../../notes/reviews/briefs/2026-10-10-float-geometry-v3-v4.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-10
- **commit**: `cc0a19b4` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.7 min, rc=0
- **結論**: **V4 の前に、同じ Q32 から FP64 の流束直前入力を「自生した値／float 側から移した値」で切り替える1 step の A/B を行い、粘性・ω拡散を除外した解釈を修正する。**
- **extra**: `plans/active/time_integration-line-implicit-speed.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 諮問事項 | 採否・理由・対案 |
|---|---|---|
| Major | V3 の残差誤差は「粘性の面流束以外」にある | **却下**。`GEOMAB_REF` は原始量・勾配・物性・面ベクトルを float 側の入力で上書きする。測っているのは、その入力を固定した演算差であり、FP64 本番との粘性流束差全体ではない。[plan:107](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:107)。さらに集計対象は運動量とエネルギーだけで、k・ω の拡散は集計していない。[v3_an.py:48](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v3_an.py:48)。**入力生成の丸めを候補へ戻し、下記の短い A/B で確認する。** |
| Major | V3 の PASS と原因の切り分けを同一視する | **相対改善の判定は採用、原因の結論は要再検証**。記載された E_new は旧値の半分以下だが、集計器は参照ノルムゼロを除外し、面流束の非有限をゼロに置換し、内訳の集計失敗でも総合 PASS を出せる。[v3_an.py:28](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v3_an.py:28)、[同:49](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v3_an.py:49)。合成データでは**評価件数ゼロで PASS** を再現した。既知の未評価境界面だけを除外し、必要な量・領域の欠損は判定不能にする。残差改善と内訳の判定を分ける。 |
| Major | E_new = 0.37 が float の ω の床を決める | **要再検証**。これは固定した Q32 における作用素の差であり、反復後の残差床でも、その下限でもない。[plan:396](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:396)、[同:443](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:443)。**V4 の末尾系列で確認し、0.37 から床の値を予測しない。** |
| Major | V4 を「水準到達＋到達時のフラグ」に変更する | **目的を限定して採用**。「水準で止める」はユーザ決定と整合する。[speed plan:420](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:420)。ただし「VERDICT は合否に使わない」は広すぎる。既存の判定器も DIVERGED・記録なし・判定不能を先に除外している。[tt_judge.py:75](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/tt_judge.py:75)。**NOT CONVERGED／非定常を許容することと、発散・不完全データを許容することを分ける。** |
| Minor | FP64 の B0 自身が旧 V4 を満たせない | **根拠不足**。ドリフト ≤0.05 % は ≤0.01 % と両立する。実際、既存 B0 の最初の到達では +0.007/−0.004/+0.009 % と記録されている。[speed plan:431](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:431)。これは窓全体の合格証拠でもないが、「満たせない」の証拠にはならない。**変更理由は不可能性ではなく、評価目的を定常同等性から停止規則下の比較へ変更すること、と書く。** |
| Minor | r0・r1 のメモリを直ちに削減する | **現状の受け入れを推奨**。記載の面数から、float の r0・r1 は約27.5 MB、e を含め約41.3 MBずつ、ホストとデバイスに増える見積もり。[plan:367](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:367)、[同:416](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:416)。**今は配列構成を固定し、見積もりを訂正してピーク使用量と V5 を測る。** ホストの写しの解放ではデバイス容量は減らない。 |

V4 は、次の形なら事前登録できる。

1. **停止**：段④の同じコード系列の float／FP64、同じメッシュ・起点・B0 の実効設定で比較する。2500 step ごとに記録し、水準を2出力連続で満たした**2回目**で停止。上限は両腕とも追加20万 step。参照との値の差や準定常判定を、停止後の延長理由にしない。

2. **比較成立条件**：両腕の到達を保存系列から再計算して確認する。末尾2万 step の9点と全必須残差列がそろい、非有限・発散がないことを要求する。現状の `m9_watch.py` は出力間隔・点数を確認しておらず、合成系列 `[2500, 22500, 25000]` の**3点だけでも REACH** を返した。[m9_watch.py:110](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:110)、[同:130](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:130)。欠損系列を到達扱いしない検査を加える。

3. **フラグ**：各腕自身の到達時を比較し、θ_r 各断面は0.05 %、符号付き Q_w は0.1 %、末尾2万 step の各 `rms_*` 中央値は FP64 の2倍を上限とする。相対差の分母は FP64 到達値の絶対値、と明記する。参照中央値ゼロを無条件に無視しない。R3 を踏襲するなら、元の **2πΣ|res_ro| の2倍条件も残す**。[speed plan:503](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:503)。

4. **判定名**：「登録した水準へ到達し、停止時比較のフラグなし」。これは定常解の同等性や物理精度0.05 %の保証ではない。`check_convergence` の NOT CONVERGED、`check_quasisteady` の非定常判定は、そのまま併記する。旧 V4 の厳しい窓条件は撤回した履歴を残す。

5. **記録・速度**：末尾9点の系列と必要な場を保存する。現状の見張りは最新・到達以外を削除するため、そのまま流用しない。[m9_watch.py:123](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:123)。共有 GPU の実経過時間は参考記録とし、採用に使う総時間は V5 の専有単価・出力費用・起動費用から比較する。

結論: **V4 の前に、同じ Q32 から FP64 の流束直前入力を「自生した値／float 側から移した値」で切り替える1 step の A/B を行い、粘性・ω拡散を除外した解釈を修正する。**

第 1 仮説: **ω拡散へ入る原始量・勾配・物性・幾何入力の丸め差が、残った ω 残差差の過半を説明する。** 確度: **低・未確認**  
　根拠: 現在の診断はこの差を上書きで消している。最大誤差が j=119 にあるという記録は、この候補と整合するが原因の証明ではない。[plan:107](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:107)、[同:447](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:447)。  
　反証条件: 下記で入力差の寄与を符号付きで差し引いても、ω残差差のノルムが元の90 %以上残り、再実行差で説明できない場合。

第 2 仮説: **流束・ソースの演算と残差への加算の丸めが主因。** 未確認。再実行差 N が小さくても、再現性のある丸め誤差は除外できない。  
第 3 仮説: **commit の更新消失が V4 の到達を妨げる。** 未確認。V3 は commit 前なので、今回の E_new の説明には使えない。

判別 A/B: **切り替えるのは FP64 の流束直前入力の由来だけ。**

- 同じ Q32・同じ FP64 バイナリで、A は自生した FP64 入力の `GEOMAB_DUMP`、B は float 入力を移す既存の `GEOMAB_REF`。各1 step、各2回。物理設定は固定する。
- 粘性に加え、**ω拡散の面流束**を必ず集計する。入力差による寄与 C = Σ±(F_B − F_A) を、本番と同じ符号で FP64 集約する。
- 主判定を壁際 j≥117 の ω に固定し、既存の ΔR = R32 − R64 に対して `r = ‖ΔR − C‖₂ / ‖ΔR‖₂` を測る。必要な内部面の欠損・非有限は判定不能。
- **r ≤ 0.5**、かつ再実行差が ‖ΔR‖ の10 %以下 → 第1仮説を支持し、「粘性・拡散以外」という除外を撤回。  
  **r ≥ 0.9**、同じ再現性条件 → 第1仮説の「過半を説明」を棄却。第2仮説が残るが、加算原因とまでは確定しない。  
  中間は判別不能として残す。いずれも長期の残差床の判定にはしない。

やらない方がよいこと: **0.37 を理由に `qAccumulatorFP64` の軸対称拒否を外すこと、原子加算を真因と決めること、メモリ削減を同時に入れること、水準到達後に旧 V4 の条件を満たすまで延長すること。**

呼び出し側の前提への異議: **観測された相対改善と、誤差の発生箇所の推定が混ざっている。** また、「停止時にフラグなし」と「収束した解が同等」は別の判定である。

不足情報: **自生した FP64 入力での面流束、ω拡散の集計、float 実行時の commit 診断、実測ピークメモリ。** 指示に従い巨大ファイルは読んでおらず、run 数値は plan 記載値として扱った。独立確認したのは判定コードと上記の合成データでの挙動である。ファイル変更・forge 起動は行っていない。**plan 未反映**。反映先は float plan §4.4・§5.1 #7〜#12・§6 V4・§6.10。
