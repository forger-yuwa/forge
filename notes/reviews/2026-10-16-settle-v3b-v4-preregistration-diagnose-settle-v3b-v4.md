# codex 諮問 (diagnose): settle-v3b-v4-preregistration

- **brief**: [`notes/reviews/briefs/2026-10-16-settle-v3b-v4-preregistration.md`](../../notes/reviews/briefs/2026-10-16-settle-v3b-v4-preregistration.md)
- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **date**: 2026-10-16
- **commit**: `6151511f` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.8 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-settle-v3b-v4-preregistration-diagnose-settle-v3b-v4.codex.log` (codex: ok)
- **結論**: 閾値を現行の暫定値に固定し、`case/71.m55_combustion_wt/run_0012_ns_free/` の単独実行で「仮想停止後の同一軌道」を検証する事前登録へ改める。
- **extra**: `case/71.m55_combustion_wt/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（現案のままの閾値確定は **NO-GO**。Critical 0／Major 4／Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| — | 目的量の許容を固定：**採用** | 観測した変動に合わせて設計公差を緩めない方針は妥当。ただし、現在の判定は値の設計適合ではなく**窓内の変化量**を見ている（[settle.py:1311](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1311)、[同:1370](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1370)）。設計適合の検査は別に残す。 |
| **Major 1** | T_q を整定の正解にする：**却下** | 「終端まで各窓が合格」は、累積ドリフトが小さいことを保証しない。既存の `_eval_exit_core` をメモリ上の人工系列に適用すると、**2万 step 当たり ΔM＝0.0002、許容0.000275で、42k〜200kの80窓すべてが SETTLED**。しかし仮想停止82kから200kまで ΔM＝0.00118、許容の4.29倍になる。[実装:1311](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1311)。対案：T_q は有限期間の記述値に限定し、**停止候補時点を基準にした、その後の累積変化**で検証する。 |
| **Major 2** | 残差・流量を一括して最大×1.5：**却下** | [plan:125](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:125) は、流量の保存誤差まで実測値に追随させる。既存の `_eval_mass` で、一定の**流量不釣合い1%**は現行 `mean_abs=5e−4` なら NOT_SETTLED、提案式による0.015なら SETTLEDになった。実装は [settle.py:1148](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1148)。対案：流量の許容は固定する。残差の×1.5は安全率ではなく**候補値の生成法**に限定し、未使用の系列で追加の判別能力を示すまで確定しない。 |
| **Major 3** | q_tol・ε_relも同じ較正：**却下** | 判定は `q_tol·A_i + eps_rel·abs(Q_ref,i)` の和であり、出力される最悪値はこの許容で割った無次元比（[settle.py:1441](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1441)、[同:1454](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1454)）。その最大×1.5から、単位も異なる二つの係数は一意に決まらない。対案：今回は両方を現行値に固定して暫定検証する。確定には熱設計上の絶対・相対許容を別々に与える。 |
| **Major 4** | 同条件の再実行で独立検証と小型ケースを代替：**却下** | [plan:123](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:123) は暫定閾値で停止し、[同:125](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:125) は後で固定した閾値を再生する。必要区間が残る保証がなく、SETTLEDが一度も出なくても「T_qより前には出ない」を満たせてしまう。同条件の反復は再現性の確認には使えるが、別条件への適用性の確認にはならない。対案：自由継続runを先行させ、閾値固定後にV4を実行する。**停止到達の正例が実在すること**を必須にし、小型ケース未実施ならプリセット全体の `provisional: false` は保留する。 |
| **Minor 1** | 同一GPUで同時実行：**却下** | V4は判定時間・停止遅延を測る計画だが、監視にはポーリングと判定処理が入る（[stop_rule.py:1549](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1549)、[同:1580](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1580)）。対案：**1本ずつ実行**し、判定開始〜終了、要求設置〜停止通知を秒とstepで分けて記録する。現在の `judged_at` だけでは判定所要時間は測れない。 |

結論: **閾値を現行の暫定値に固定し、`case/71.m55_combustion_wt/run_0012_ns_free/` の単独実行で「仮想停止後の同一軌道」を検証する事前登録へ改める。**

第 1 仮説: 現案では、短窓では小さい持続的ドリフトを「目的量が整定した正例」に取り込み、その区間に合わせた較正で誤停止を見逃す。確度: **高（判定構造について。case/71での発生は未確認）**  
  根拠: 上記の人工系列で80窓すべてがSETTLEDでも、仮想停止後の累積変化が許容の4.29倍になった。これは既存の判定関数の実行結果であり、CFD実測ではない。  
  反証条件: case/71で、固定閾値による仮想停止後も、事前指定した継続区間全体で目的量の累積変化が許容内に残ること。その場合、このrunで当該誤停止が起きるという仮説を退ける。

第 2・第 3 仮説: 追加しない。現時点で物理・メッシュ・精度を原因に持ち出す根拠はない。

判別 A/B: **変えるのは観測終端だけ**とする。自由継続runを最大200,000 stepまで1本実行し、A＝固定閾値が40,000 step保持された最初の停止候補 t_s、B＝同じ軌道の t_s＋80,000 stepまでを再生する。追加の比較用CFDは不要。  
  出口コアM・線の各位置・流れ角・壁の局所熱流について、**t_sの値からの差の最大**を固定許容で調べ、残差・流量を含む全条件も各出力時点で確認する。  
  → **A判定（全期間で許容内・全条件合格）なら**、この条件・80,000 stepの範囲で早期停止の懸念を退ける。  
  → **B判定（いずれか超過）なら**、現在の停止判定を反証し、閾値を緩めず停止規則を見直す。  
  t_sが120,000を超える、または停止候補が無い場合は**検証未成立**とする。80,000 stepは今回の有限な検証範囲であり、永久整定の保証とは扱わない。

やらない方がよいこと: 暫定閾値で先に止めたrunの欠けた将来を、非決定的な別runで補ったことにすること。×1.5を安全性の根拠にすること。残差条件がANDに入っているだけで「独立の情報を追加した」とすること。追加能力は、目的量だけなら通る不適切な区間を残差が落とせるかで確認する。

呼び出し側の前提への異議: `run_0010` は新方式の正例として確定していない。保存済み判定は、旧系列の **2,000〜86,000 stepで `OVERALL: ALL STEADY`**（[記録:12](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0010_ns_b0/judge/check_quasisteady_series.txt:12)）、一方で残差は **終端88,519で `NOT CONVERGED`**（[記録:2](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/run_0010_ns_b0/CONVERGENCE_VERDICT.txt:2)）。これは旧判定の観測事実であり、新しい数値流束・局所壁熱流を含む全条件の正解ラベルではない。また、`run_0011` の早期停止はV4の動作試験には使えるが、停止後の妥当性検証には使えない。

不足情報: 小型ノズルの具体的なcase名は、指定された閲覧範囲では選定できない。既存候補のREADME・設定・判定記録が必要。残差の確定値と壁熱流の許容値も今回の証拠だけでは決められないが、現案を却下して上記検証を先行させる判断は妨げない。

証拠: 足りる
