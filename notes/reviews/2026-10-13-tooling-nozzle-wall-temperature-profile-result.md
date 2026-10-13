# codex レビュー: tooling-nozzle-wall-temperature-profile (result)

- **plan**: [`plans/active/tooling-nozzle-wall-temperature-profile.md`](../../plans/active/tooling-nozzle-wall-temperature-profile.md)
- **stage**: `result` (diff base `138e83b2`)
- **date**: 2026-10-13
- **commit**: `cf35b7c1` (feature/nozzle-deliverables)
- **codex**: effort `high`, 4.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m4
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

主要な配管とデモの合格主張は、コード・保存済み数値で裏付けられました。Critical / Major は認めません。
ただし、代表温度の計算と検証記録に不整合があります。以下を直してから `accepted` へ移すことを推奨します。

1. **Minor — `s = 0` を含まない表では、代表温度が実際のスロート温度と異なる。**

   根拠: [wall_temperature.py:102](/home/sano/work/forge-deliv/design/forge_design/wall_temperature.py:102) は元の `s` 座標で代表値を補間します。一方、実際の壁温は同ファイル:139 で実寸へ写してから補間します。

   読み取り専用の再現計算で、`s=[-2,-1,1]`、`Tw_K=[300,400,600]`、直管 0.2 m・縮流部 0.3 m・拡大部 1.2 m の場合、**代表値は 500 K、実際の x=0 の温度は 440 K** でした。入力検査はこの表を受理します。CSV による温度拘束は正しいものの、JSON・bcond・Excel の「設計スロートの値」という表示が誤ります。

   **対案:** 寸法解決後の実寸表から x=0 の温度を計算し、bcond と各帳簿の代表値を統一する。この再現例を回帰試験に追加する。

2. **Minor — V1e の「出口 0.155 mm」は、出口より約 25 mm 上流の値。**

   根拠: [walltemp_wall_ab.py:65](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/walltemp_wall_ab.py:65) の評価終点は `x_e−1`。同ファイル:73 はその値を `dr_ref_at_exit_m` と記録し、[plan:138](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:138) も「出口」と記載しています。

   保存された A/B の表から再計算すると、次の値です。

   | 評価位置 | x [m] | Δr_ref [mm] |
   |---|---:|---:|
   | 登録された検査区間の終点 `x_e−1` | 1.1751193 | 0.1552999 |
   | 実際の出口 `x_e` | 1.2000892 | 0.1590835 |

   **対案:** 既存値を「検査区間終点」と改名し、出口値は別項目で記録する。**登録区間での判定 A 自体は有効**で、伝達誤差 7.0690e−17 m も再現できました。

3. **Minor — plan 内で、確定した API・正本の扱いと検証契約が矛盾している。**

   根拠:
   - [plan:76](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:76) は分布時の `wall_thermal_bc_integral` を例外としますが、:79–80 には旧来の「`spec.r_throat` で換算する互換」が残っています。実装は例外です。
   - 同 plan §4.4 は異なる `--problem` を渡しても run の表を採用するとしていますが、§6 V1d（:141）は別の問題 YAML で成果物出力が停止すると要求しています。実装・試験は前者です。
   - §4.2 に記載された較正スクリプトの分布未対応と、§4.4 の Fluent 読み込み未確認は、§5.1 に明示的な追跡先がありません。

   **対案:** 旧案を履歴として明確に区別し、V1d を「run 内の正本破損は拒否、外部 YAML の壁温では上書きしない」に同期する。較正スクリプトと Fluent 実機確認は、範囲外・後続先を §5.1 に残す。`methods` の説明にも対応範囲を明記する。

4. **Minor — `s` 分布では直管長ゼロを扱えない制限が、入力契約に出ていない。**

   根拠: [wt_card.py:336](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:336) は `pipe_length_m=0` を受理しますが、[wall_temperature.py:151](/home/sano/work/forge-deliv/design/forge_design/wall_temperature.py:151) は入口端と縮流部始点の一致を拒否します。`anchors_m(0, 0.3, 1.2)` と正常な `s` 表で例外を再現しました。

   **対案:** 現行対応範囲を「`s` 分布は直管長 > 0」と入力契約・インテークに明記し、カード検査で早期拒否する。直管なしでは `x_m` を案内する。

検証結果の裏付けは以下のとおりです。対象は `case/70.m5_air_wt_demo/run_0011_ns_walltemp_demo/`、恒久索引は [case README](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/README.md:13) です。

| 項目 | 確認結果 |
|---|---|
| V2、`res_60000.h5` の壁 1800 節点 | 最大温度誤差 **1.5244293e−5 K** |
| V3、Fluent profile の最近傍割り当て | **1799 面、座標差 0、温度相対差 0** |
| 全 61 snapshot の `VALUE/*`・各段の全残差列 | **NaN/Inf 0** |
| 本段の収束判定を再実行 | **NOT CONVERGED (stalled/plateau)** |
| `check_quasisteady --quantity machmax,pmax` を再実行 | **OVERALL: ALL STEADY** |
| 保存された熱量の準定常判定 | **DRIFTING**、総熱量は **TRANSIENT-UNSETTLED** |
| 保存された品質・壁解像判定 | **VERDICT: PASS (AR≤5000, skew≤0.90)**／**PASS (壁解像)** |

未収束を流れの「一致」とする主張は、本計画には認めません。V2 は温度拘束の配管検査として成立しています。`machmax,pmax` の定常性から熱量の定常性は導けませんが、成果物の `thermal_use: FAIL` はこの区別を維持しています。

指定 diff にソルバ本体の変更はなく、今回の変更による境界半割面・周期 seam・軸の演算変更もありません。一様・断熱については、保存された変更前後の生成物を比較し、case/45・case/70 の config・積分法出力・壁表現のバイト一致を確認しました。書き込みを伴う試験一式は再実行していません。

**推奨は、指摘 1 → 2 → 3・4 の順に修正し、対応試験と result レビュー記録を追加してから `accepted` へ移すことです。** この配管検証を完了するための NS 延長は不要です。ファイル変更は行っておらず、今回の指摘は plan 未反映です。

指摘数: Critical 0 / Major 0 / Minor 4
