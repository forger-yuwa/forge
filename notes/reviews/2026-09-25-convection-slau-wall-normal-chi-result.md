# codex レビュー: convection-slau-wall-normal-chi (result)

- **plan**: [`plans/accepted/convection-slau-wall-normal-chi.md`](../../plans/accepted/convection-slau-wall-normal-chi.md)
- **stage**: `result` (diff base `39526328`)
- **date**: 2026-09-25
- **commit**: `e9f2f821` (feature/sern-design)
- **codex**: effort `high`, 7.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: 3 回目の result 段。§6.2 末尾の『codex result 段の指摘への対応表』の各行 (前 2 回の Major 10 件) が本当に閉じているか、『受入の範囲』の結論文が測った範囲を超えていないか、§2 制限事項と §5.1 #11 の前提が実測と整合するかに集中してほしい。前 2 回の型は『測れた範囲を超えて書いた』。新規: case/16 V2・V3 (2026-09-25, §6 V3 case/16 条件の結果)。前提ゲートの参照値を撤回・差し替えた経緯も同節にある (合否条件は不変と主張している)。
- **extra**: `notes/reviews/2026-09-23-convection-slau-wall-normal-chi-result.md`, `notes/reviews/2026-09-23-convection-slau-wall-normal-chi-result-2.md`, `case/16.nozzle_wys/v3_case16/B_10b_RULE.txt`, `case/16.nozzle_wys/v3_case16/C_CONVERGENCE_RATIO.txt`, `case/16.nozzle_wys/v3_case16/QUASISTEADY_VERDICT.txt`, `case/16.nozzle_wys/_v2/V2_VERDICT.txt`, `case/16.nozzle_wys/_v2/SIDE_new0_vs_new1.txt`, `case/16.nozzle_wys/_v2/PROVENANCE.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
`chi_mass` の実装、V5 の保存ダンプ照合、case/16・case/48 の回帰には裏付けがあります。  
ただし、V1-b・固定点・衝撃足・格子感度には未完了の受入条件があり、対応表の「Major 全件閉」は支持できません。

HEAD `e9f2f821` で指定 diff を取得しました。質量流束のみの変更、再構成後の面速度、壁 mask、境界半割面の除外、既定0、config の拒否条件は設計と整合しています。追加流束式に明白な符号・単位・ゼロ割の誤りは見つかりませんでした。

独立に確認できた範囲は次のとおりです。

| 対象 | 再確認結果 |
|---|---|
| Python 単体試験／case/48 V5 保存ダンプ | `ALL PASS`／P0–P5 `PASS` |
| case/16・case/48 V2 | 両方 `VERDICT: PASS` |
| case/16 `run_0506`／`run_0507` | 残差は両方 `NOT CONVERGED`。指定閾値で輪郭壁197点すべて `STEADY`。末尾平均輪郭の最大相対差 **0.0000464 %** |
| case/48 `_v3/v3_flag{0,1}` | 残差は両方 `NOT CONVERGED`。各11 dump・5地点の `Cf`／`qw` は指定閾値で全て `STEADY` |
| 上記 case/16・case/48 V3 | 全スナップショットの `VALUE/*` に非有限値なし |

case/48 は `/home/sano/work/forge/case/48.flat_plate_cooled_m4/` の実データを使用しました。case/46 の主要な AWS run、SERN 延長 run、V7、V6 の実測ファイルは手元で確認できず、以下ではコード上の事実と plan 記載値に基づく判断を区別します。

1. **Major — V1-b の未測定部分が残り、R2-M4 を閉じられません。**

   **根拠:** [plan:233](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:233) は、起点直後の正味流入と、局所時間刻みを含む正規化収支を要求しています。しかし「完全な記録」は [plan:1153](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1153) で **`Δt_l` 未出力**と明記し、単位 `[1/s]` の値だけを示しています。これでは「1 %/dump 未満」を判定できません。起点直後の flag 1 の負の収支も提示されていません。

   また、[受入結論:1199](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1199) の「V1-b 観測/予測1.000」は **V1-c の検証量**です。V5 の成功は診断式の裏付けになりますが、V1-b の不足を代替しません。

   **対案:** 起点直後の全接続面収支と、末尾3 dump の `|Σṁ|Δt_l/(ρV)` を保存して判定する。床処理を省略するなら、監視壁だけでなく接続内点を含む前処理の不活性も確認し、R2-M4 をそれまで未完了へ戻してください。

2. **Major — V7 は有限区間の小さな差を示していますが、「CFL に依存する固定点」までは示していません。**

   **根拠:** [plan:490](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:490) は線形ドリフトを過渡として不合格にする設計です。一方、[結果:536](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:536) は、窓長に比例するドリフト、2点比1.42の未定義帯、ドリフト側という解釈を明記しています。それでも #10c と R2-M2 を閉じています。

   `--from-floor: PASS` と場の距離が ε 以下であることは有用ですが、ドリフトの減衰や漸近先を証明しません。累積 CFL を揃えた差が持続することも、異なる過渡の差を排除しません。

   **対案:** 現結果を「測定区間の場の変化・CFL間差は規定 ε 内」と限定し、固定点の断定を撤回する。独立した後続窓でドリフトの減衰・飽和を判定し、組成の遅いモードも観察対象として残してください。固定点確認を受入条件として維持する限り、#10c は未完了です。

3. **Major — SERN の衝撃足は、各 run の準定常確認がまだ差系列で代用されています。**

   **根拠:** [plan:950](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:950) は各 run の `C_T` を判定していますが、衝撃足について示すのは **`foot_L2_pct = 0.8385 ± 0.0091 %` と位置差0**です。[生成コード:61](/home/sano/work/forge-sern-design/case/46.sern_design/v3sern_series.py:61) も、両 run の壁圧差ノルムと位置差を出しています。両側が共通して動けば、差は小さく安定していても各場は未定常です。

   微小な差系列に相対閾値を掛けた誤りを直す判断は妥当です。しかし、修正後も**各 run の壁圧分布・衝撃足位置**の判定が必要です。また同節で `C_L` は flag 0 が `DRIFTING` なのに、[結論:974](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:974) は限界サイクル振幅による「実効果」の確定を維持しています。

   **対案:** 各側の壁圧分布と絶対衝撃足位置を指定閾値で判定し、その後に差ノルムを評価する。`C_L`／`C_M` は現区間の差として記録し、§5.1 #11 の独立精度評価には準定常確認を明記してください。`C_T` の改善は認めますが、R2-M1 全体は閉じません。

4. **Major — 3格子の起動成功を、格子感度ゲートの完了として扱っています。**

   **根拠:** [plan:324](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:324) は収束悪化・格子感度を受入対象にしています。しかし [結果:1049](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1049) は格子間の設定一致を保証できず、[同:1059](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1059) は定常解・壁圧・収束率の格子比較が未実施と明記しています。追加したのは起動成否の証拠であり、残っていた懸念の検証ではありません。

   さらに [同:1062](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1062) と1076行には、撤回したはずの「細かいほど早い」が現役の説明として残っています。

   **対案:** #9 を「起動成否3水準は済／壁圧・収束性の格子評価は未」に分ける。設定を照合可能な比較で残りを評価し、R1-M3・R2-M3 と #11 の前提を同期してください。

5. **Major — 新フラグが収束判定区間の識別から漏れています。**

   **根拠:** [stage_manifest.py:80](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:80) の hard キーに `space.slauWallNormalChi` がありません。case/16 の実 config を同じ YAML 書式に揃え、**この値だけ0→1**にしたところ、`stage_key` が同一となり、[segments:214](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:214) は **1区間**を返しました。

   段途中で flag を切り替えると、異なる空間離散化の残差を連結し、変更前の過渡を低下桁数の基準にできます。今回の独立した2本の比較が直ちに無効になるという指摘ではなく、新機能と既存の収束判定機構の接続漏れです。

   **対案:** 新フラグを hard キーへ追加し、省略と明示0は同一、0と1は別区間になることを確認してください。

6. **Minor — 撤回済みの数値・結論と台帳が同期していません。**

   **根拠:** [plans/README:62](/home/sano/work/forge-sern-design/plans/README.md:62) は「本命は未確定」のままです。[case/46 README:344](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:344) には「5桁一致＝独立」、続く行には旧固定点 PASS・格子非依存の説明が残ります。step 800 の記載値 `2.6535e−5` と `2.6542e−5` の差は **0.0264 %**で、5有効桁の一致でもありません。

   case/48 の保存評価を再計算した最大差は **`Cf` 0.001804 %／`qw` 0.002388 %**で、plan の0.000 %／0.011 %と異なります。許容内という結論は変わりません。[methods の仕様節](/home/sano/work/forge-sern-design/methods/convection/theory.md:377) にも、生産検証範囲と周期・軸対称の小規模試験の区別がありません。

   **対案:** 最新結果を正本にし、旧主張には撤回・後継節を付ける。run 索引、計画一覧、methods、対応表を同時に同期してください。

7. **Minor — case/16 の合格から、ノイズや局所機序を言い過ぎています。**

   **根拠:** [受入結論:1200](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1200) は床比 **0.80–1.47**を「反復ノイズ以下」としています。これは正確には「設定したノイズ床2倍以内」です。また場の差が小さいことだけから、[同:1204](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:1204) の原因を `χ_n≈χ` と確定できません。流束差には面圧力差も掛かり、実装が読むのは再構成後の速度です。

   **対案:** 「この500 step試験では反復ノイズと同程度、V3許容内」とする。機序を確定するなら、当該ケースの再構成面状態による `Δχ`・`Δp` を測定してください。

case/16 の前提参照値を変更した経緯は明示されており、**2D層流の稼働確認として使う限り、変更自体を不合格理由にはしません**。新しい帯は実験精度の検証ではなく、V3 の相対差・準定常条件とは別物、という限定を維持してください。§2 の周期・軸対称を小規模試験までに限定する方針も適切です。

**推奨は、既定0・`active / in_progress` を維持し、上記Majorを解消して再レビューすることです。** R1-M1 の診断式修正とR1-M5／R2-M5の適用範囲修正には実質的な裏付けがありますが、対応表全体を「閉」とする段階には達していません。ファイル変更・plan反映はしていません。

指摘数: Critical 0 / Major 5 / Minor 2
