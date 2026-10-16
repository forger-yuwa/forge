# codex 諮問 (diagnose): m55-ns-unsettled-2

- **brief**: [`notes/reviews/briefs/2026-10-15-m55-ns-unsettled-2.md`](../../notes/reviews/briefs/2026-10-15-m55-ns-unsettled-2.md)
- **date**: 2026-10-15
- **commit**: `9c06989a` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-m55-ns-unsettled-2-diagnose.codex.log` (codex: ok)
- **結論**: 長時間延長を保留し、既存保存場について「出力原始量・h0」と「保存量から再構成した熱力学状態」の全圧 A/B を forge 0 step で行う。
- **extra**: `case/71.m55_combustion_wt/_runs/README.md`, `plans/active/tooling-nozzle-card-local-mesh.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（記録の基準パスは `case/71.m55_combustion_wt/_runs/run_0002_ns/`）:

| 対象 | 判断 | 重大度・根拠・対案 |
|---|---|---|
| 直ちに長時間延長する | **要再検証** | **Major**：本段 0〜183,998 の VERDICT は **NOT CONVERGED**、系列 1,000〜183,000 は **DRIFTING**。両判定を gz の読み取りで再現した。さらに直近20kの残差の平均 log10 は前20k比で運動量が +0.32／+0.44 桁。「全列 falling」は長い判定窓の結果で、最近も改善中という意味ではない。対案は下記の保存場 A/B。 |
| 「測定帯だけが未到達の原因」を棄却する | **採用** | `theta_band_ab.json` の後期固定帯でも Dθ = 2.332／1.960／1.748 %/20k。既存帯より大きく、帯だけでは未到達を説明できない。 |
| manifest の写しに本段を補って判定する | **採用** | **Major**：元の manifest では本段を判定できない。補完した本段の `cfg_fnv=d903391a098b4612` は起動記録と一致し、本段単独判定とも一致する。原本を保持し、「事後復元・本段0〜183,998」と明記する。これは正常終了の証明にはならない。欠落経路は [runner_axismach.py:2316](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2316) → 同 `:2341`。 |
| 追加約168kを到達予測として採用する | **却下** | **Major**：単一指数・漸近値ゼロを仮定した質量条件の試算にすぎない。後期帯の θ_r、最近の残差増加を説明しない。対案は予算上限と到達予測を分け、実測の停止条件で判断すること。 |

結論: **長時間延長を保留し、既存保存場について「出力原始量・h0」と「保存量から再構成した熱力学状態」の全圧 A/B を forge 0 step で行う。**

第 1 仮説: 全圧超過の少なくとも主要部分は保存場側に残る偏差であり、`total_quantities.py` の逆算だけが作ったものではない。確度: 中  
  根拠: `judge/final_field_check.json` は、場に紐づく解決済み物性記録・`h0_includes_k=0` を使用している。さらに入口9点では、逆算した T0 だけでなく**保存された静温自体が最大1,252.62 K > Tt=1,250 K**。ただし、この9点と全圧最大点は同じと確認されていない。後処理の経路は [total_quantities.py:132](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/total_quantities.py:132)、`h0` の出力は [output.cpp:261](/home/sano/work/forge-cardmesh/solver_density_cuda/output/output.cpp:261)。  
  反証条件: 同一保存量を同一物性で再構成すると、元の全圧超過域の超過が消え、A/B差が報告された超過を説明する場合。「保存場側にある」という仮説を撤回する。**現時点では過渡と離散化誤差を区別できない。**

第 2・第 3 仮説:
- 第2：出力 `h0`・T・P・Y と保存量の不整合が全圧超過を作っている。確度: 低、未確認。現検査は有限性・正値を調べるが相互整合を検査しない（`final_field_check.py:36`〜`:64`）。
- 第3：`thermoFloat=1` の混合精度演算が反復中の偏差に寄与している。確度: 低、未確認。FP64ビルドでも組成は float で正規化され、γ・音速にも float 経路が残る（[dependentVariables_d.cu:95](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/dependentVariables_d.cu:95)、同 `:186`・`:225`）。ただし、これだけで +0.24 %の原因とは言えない。

判別 A/B: **変更点は後処理へ渡す状態の取得元だけ。追加計算は0 step。**

- 対象は既に帯 A/B に使った100k〜183kの12枚。同じ節点・物性記録・datum・組成規約で比較する。
- **A**：現行の `total_state` による T0・P0。
- **B**：保存量 `ro, roU*, roe, roY*` から、ソルバと同じ非負化・正規化を適用した組成と内部エネルギーを作る。既存の `forge_species.thermo_gas` で TP の温度を逆算し、P と h0 を再構成して全圧を求める。定数 cp の近似は使わない。
- 記録するのは節点ごとの ΔP0/Pt、超過点数・分位点・位置、入口／軸／内部の別、およびエネルギー逆算残差。最大値同士だけを比較しない。
- **両経路の差が全比較時点で ≤10⁻⁶ Pt、かつ超過が残るなら**、出力状態の不整合を +0.24 %の主因とする第2仮説を棄却する。**Bで超過が消え、A/B差が超過を説明するなら**第1仮説を棄却し、出力・後処理経路を修正する。差が一部しか説明しなければ複合要因として記録し、過渡原因の確定には進まない。

既存の `total_state` は保存量からの独立再構成を提供しないため、この比較部分だけ薄い診断処理が必要。物性計算・定常判定は既存ツールを再利用する。これは収束した解同士の比較ではなく、**同一保存場内の整合検査**である。

やらない方がよいこと:
- `thermoFloat: 0`、CFL、SLAU、軸設定を先に変更すること。保存場の監査前では原因を混ぜる。
- 追加168kで届くと約束すること、一般準定常判定の `STEADY` を停止規則の合格に読み替えること。
- 継続で取り直した帯の θ_r を元の帯の系列へそのまま連結すること。

呼び出し側の前提への異議:
- **Major：「等エントロピーの核だから入口全圧を超えない」は、今回の保存場について未確認の前提。** 観測事実は「183kで計算した P0/Pt 最大1.002358」。定常性も、核のエントロピー保存も未確認なので、これだけからソルバの特定箇所を犯人にしない。対案は上記A/B後、同じ領域の h0・全圧超過の時系列を `check_quasisteady --series-csv` に渡す。
- **Major：`evaluate.stop_rule: true` だけでは記録失敗を防げない。** `save_record` の書込み例外はそのまま伝播する（[stop_rule.py:660](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:660)）。継続前に OSError の原因と子プロセス終了を確認する。
- 延長へ進む場合は、新しい run で同一メッシュの保存量を index restart し、前段を再実行しない。出力間隔2,000、停止規則明示、保持40kを維持。上限は**追加200kの予算上限**とし、20kごとに再評価する。到達予測ではない。
- 新しい帯と Euler の `res_12000`・sha256を固定し、θ_r は別系列にする。旧系列とつなぐなら旧保存場も同じ帯で再計算する。残差は実効設定が同じ区間だけ連結し、183k restartなら旧183k以後を重複させない。
- 容量記録から、追加200k・2k出力は全場・壁・出口・残差で**約10.3 GB**、格子等を含め約10.5 GB。記録時の空き16.23 GBには収まるが余裕は約5.7 GB。投入時に再測定し、停止判定の直近40k、restart元、固定Euler参照を確保する。削除より検証済み退避を優先する。

不足情報: HDF5は手元になく、全圧の場統計は回収済み検査出力に基づく。次のA/BにはAWS側の既存保存場が必要。OSErrorの errno と発生時の容量・quotaは継続前に確認したい。いずれも、今回の「まず保存場のA/B」という判断は妨げない。

証拠: 足りる
