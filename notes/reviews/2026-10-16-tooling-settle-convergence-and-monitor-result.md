**未成立 (incomplete)** — incomplete: 証拠不足: V4 `case/71.m55_combustion_wt/run_0011_ns_v2/` の実行・判定記録。手元に実 run がなく、オンライン13件との一致、停止時再評価、容量・所要時間の主張を確認できません。 / V6 `case/71.m55_combustion_wt/run_0014_chain_v2/` の全残差履歴、最終場 `res_46171.h5`、判定区間の節点・境界出力。非有限ゼロと停止時 `SETTLED` を独立に再検証できませ …

# codex レビュー: tooling-settle-convergence-and-monitor (result)

- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-16
- **commit**: `9ec0000a` (feature/card-local-mesh)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 8.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result.codex.log` (codex: incomplete)
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **extra**: `case/71.m55_combustion_wt/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

`accepted` への移行は保留です。本段の照合を回避する経路、出力失敗の見逃し、流量判定のゼロ割が残っています。
保存済みの検証記録には整合する部分がありますが、V4・V6 の一次データが不足しており、結果レビューは未成立です。ファイルは変更していません。

1. **Major — 本段の照合不一致が旧判定器への切替になり、通常の CLI では照合自体も省略される。**  
   [settle.py:1885](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1885) は、起動と本段記録の FNV が異なるだけで `check_convergence.py` を選びます。`run_0014_chain_v2` の記録をメモリ上で差し替え、本段記録の FNV だけを変更すると `not_main / problems=[]` になることを再現しました。起動記録の欠落も `main / problems=[]` になります。さらに [check_settled.py:110](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_settled.py:110) は `--bind` 指定時しか照合せず、生成された基準ファイルに記載する通常コマンドでは照合を省略します。  
   **対案:** 本段の記録がある場合は照合を既定で必須にし、不一致・起動記録欠落は `INDETERMINATE` にする。旧判定への切替は、前段であることを独立に確認できた場合に限定してください。

2. **Major — 途中まで書かれた XMF を出力成功と認定し得る。**  
   [output.cpp:274](/home/sano/work/forge-cardmesh/solver_density_cuda/output/output.cpp:274) の全場 XMF と同ファイルの境界 XMF は、書込み・flush・close の失敗を検査していません。一方、[monitorOutputs.cpp:258](/home/sano/work/forge-cardmesh/solver_density_cuda/output/monitorOutputs.cpp:258) の確認はファイルサイズが正かどうかだけです。非空の不完全な XMF が残る失敗では、`STOPPED_AT` や随時出力の成功通知を出せます。これは「書込み成功時だけ完了通知」という V1 の条件を満たしません。  
   **対案:** 各出力処理で flush・close までの成否を確認し、失敗を完了通知側へ伝播する。非空ファイルが残る途中書込み失敗も試験に追加してください。

3. **Major — 出口流量がゼロだと判定結果を返さず例外になる。**  
   [settle.py:1240](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1240) は、出口流量の平均で無条件に除算します。実 run の基準を使い、入口流量が正・出口流量がゼロの有限な6観測点を与えると、`ZeroDivisionError` を再現しました。停止時観測を加える同様の除算もあります。オンラインでは [stop_rule.py:1630](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1630) の例外処理に入り、計算を強制停止する経路です。  
   **対案:** 分母の成立条件を検査し、流出なし・平均ゼロ・逆流を理由付きの不合格として返す。ゼロ流出と符号反転の負例を追加してください。

4. **Minor — 現在仕様と残作業の整理が追いついていない。**  
   [recommended-settings.md:130](/home/sano/work/forge-cardmesh/procedures/recommended-settings.md:130) は、まとめ前に旧判定器を一律必須としています。[methods/design/overview.md:1224](/home/sano/work/forge-cardmesh/methods/design/overview.md:1224) にも新方式の分岐がなく、更新済みの `calculation-workflow.md` と揃っていません。また、別条件の小さいノズルの確認は [plan:162](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:162) で保留ですが、§5.1 に独立した未完了項目として残っていません。  
   **対案:** 現在仕様・推奨手順を更新し、残る検証は§5.1に明記する。今回受け入れる範囲と、`provisional` 解除に必要な作業を分けてください。

5. **Minor — V3 の参考値「最大の比 0.05」が保存記録と異なる。**  
   [plan:158](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:158) の t_s′＝30,000 からの検証について、保存された `VERDICT.txt` の壁熱流の最大比は run_0012 が **0.06269**、[run_0013 が **0.08724**](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/run_0012_ns_free_judge/v2_replay/run_0013_ns_free260/VERDICT.txt:13) です。約0.05は出口コア Mach の値です。  
   **対案:** 全量の最大比を約0.0873に訂正する。許容内という結論は変わりません。

確認できた保存記録では、V3 の両 run は `VERDICT: A`、t_s＝46,000、検証区間46,000〜126,000でした。V6 の `run_0014_chain_v2/` は `REACHED`、保存済み判定は `SETTLED（閾値は暫定）`、関門は上書きなしで通過しており、`delta_r_next.csv` は停止時の抽出表とバイト同一です。run の索引は [case README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md) にあります。ただし、これらは保存記録の確認であり、以下の一次データの再検証はできていません。

**推奨:** `active` に留め、上記修正と負例検証、証拠の補完を終えてから結果レビューをやり直してください。

証拠不足: V4 `case/71.m55_combustion_wt/run_0011_ns_v2/` の実行・判定記録。手元に実 run がなく、オンライン13件との一致、停止時再評価、容量・所要時間の主張を確認できません。

証拠不足: V6 `case/71.m55_combustion_wt/run_0014_chain_v2/` の全残差履歴、最終場 `res_46171.h5`、判定区間の節点・境界出力。非有限ゼロと停止時 `SETTLED` を独立に再検証できません。

指摘数: Critical 0 / Major 3 / Minor 2
