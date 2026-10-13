# codex レビュー: tooling-nozzle-mesh-size-guard (result)

- **plan**: [`plans/active/tooling-nozzle-mesh-size-guard.md`](../../plans/active/tooling-nozzle-mesh-size-guard.md)
- **stage**: `result` (diff base `88675245`)
- **date**: 2026-10-14
- **commit**: `5a752714` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.4 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-tooling-nozzle-mesh-size-guard-result.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M1/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

通常経路の節点数ガードと既定維持は確認できました。ただし、計数処理自身に巨大配列を確保する経路が残っています。
V1b の数値は再現でき、追加 CFD は不要です。以下を修正してから `accepted` に移すことを推奨します。

1. **Major — 節点数ガード自身が、上限判定前に巨大配列を確保する。**

   根拠: [mesh2d.py:97](/home/sano/work/forge-deliv/design/forge_design/meshing/mesh2d.py:97) の `_capped_layout` は、打ち切り判定より先に `_radial_fracs(nj_base, ...)` と `np.diff` を実行します。[同:76](/home/sano/work/forge-deliv/design/forge_design/meshing/mesh2d.py:76) は指定 `nj` 全体の配列を確保します。

   実行確認では、`ni=2000, nj=1000000000, wall_first_frac=0.005, axis_gap_frac=0.005`、既定上限30万点に対して、**10億要素・8 GB の `np.linspace` 呼び出しへ到達**しました。確保直前で差し替え関数により停止させています。§4.2 の「点数の計算そのものも大きな配列を作らない」を満たしていません。既存試験の `nj=100000`・時間制限だけでは検出できません。

   **対案:** 径方向の計数を、指定 `nj` に比例する配列を確保せずに打ち切れる処理へ変更する。巨大 `nj` で配列確保前に拒否する回帰試験を追加し、§5.1 #3 を修正完了まで再開する。

2. **Minor — V1b を「判定不能で閉じる」という決定が本文・現在仕様へ徹底されていない。**

   根拠: [plan:110](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mesh-size-guard.md:110) は「δ_E・M 合格」、直後の113行は「δ_E も単一時刻なので合格とはしない」と矛盾します。[methods/design/overview.md:362](/home/sano/work/forge-deliv/methods/design/overview.md:362) にも「2次の格子収束に乗る」という断定が残っています。

   再実行した判定は次のとおりです。対象はいずれも `case/70.m5_air_wt_demo/` 配下、本段48,000 stepです。

   | run | `check_convergence` | 出口 M／流量の `check_quasisteady` |
   |---|---|---|
   | `run_0015_euler_c2b_ni1800` | NOT CONVERGED：stalled/plateau | STEADY／OSCILLATING |
   | `run_0016_euler_c2b_ni900` | NOT CONVERGED：still converging | STEADY／OSCILLATING |
   | `run_0017_euler_c2b_ni600` | NOT CONVERGED：stalled/plateau | STEADY／OSCILLATING |

   出口 M の末尾5枚平均は **4.999365354／4.999518066／4.999743543**。plan の数値と整合しますが、これだけで漸近的な2次格子収束を確立したとは言えません。各 run の `CONVERGENCE_VERDICT.txt`、`quantities_series.csv`、[README の run 一覧](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/README.md:33) も照合しました。

   **対案:** 「δ_E は数値上の閾値内だが合否保留」「case/70 の感度評価であり、収束次数は未確定」に統一する。Euler の `ni` を切り離さない決定は維持し、追加計算は要求しません。

3. **Minor — V4 の実 `prepare-only` 検証に成果物への参照がない。**

   根拠: [plan:129](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mesh-size-guard.md:129) は `n_nodes=174600`・品質 PASS を報告していますが、出力先を記載していません。case/70 内では、新キー `max_nodes_source` を含む `prepare_info.json` を確認できず、README にも対応する検証先がありません。

   **対案:** 実行先と `prepare_info.json`・`MESH_QUALITY.txt` を記録へリンクする。一時ディレクトリを破棄済みなら、その事実と保存した検証ログを明記する。現状ではこの実行結果を独立に追跡できません。

回帰については、既存試験の書き込み不要部分を再実行して失敗0。5経路の座標・接続・境界辺は `88675245` とビット一致しました。通常の A/B も194,000点で通過、402,000点で拒否を再現しています。停止対象の既存 YAML 7件と節点数も§7どおりでした。CUDA 流束・境界半割面・周期 seam の計算変更はありません。

V1b の最終 `res_48000.h5` は3本とも `VALUE/*` の NaN/Inf が0、ρ・P・Tは正値。保存された品質判定も全て `VERDICT: PASS` でした。全54項目と既存12本の試験は、書き込みを伴うため今回は再実行していません。既存の未解決事項である他 runner と `nj` 既定値差は、§5.1 #6・#7 に残っています。

**推奨:** Euler の `ni` 維持・30万点上限という方針を採用し、移動前に **①計数時の巨大確保を修正・再試験、②V1b の記述を統一、③V4 の証拠を紐付け**の順で完了させてください。

指摘数: Critical 0 / Major 1 / Minor 2
