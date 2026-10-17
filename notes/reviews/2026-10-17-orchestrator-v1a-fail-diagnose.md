# codex 諮問 (diagnose): orchestrator-v1a-fail

- **brief**: [`notes/reviews/briefs/2026-10-17-orchestrator-v1a-fail.md`](../../notes/reviews/briefs/2026-10-17-orchestrator-v1a-fail.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-17
- **commit**: `28945473` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-orchestrator-v1a-fail-diagnose.codex.log` (codex: ok)
- **結論**: 元の V1a＝FAIL を保持し、比較項目を明確に事前登録したうえで、同一の現行コード・環境による旧呼び出し対 CLI の prepare-only A/B を１組行う。
- **extra**: `case/71.m55_combustion_wt/_runs/v1a_executor/result.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（`28945473` を基準にレビュー。ファイル変更・forge 起動は未実施、plan 未反映）:

| 対象 | 採否 | 根拠・対案 |
|---|---|---|
| 同じ現行コード・同じ環境で旧呼び出しと CLI を比較する | **採用。ただし解釈を限定（Major）** | 両経路は同じ `run_pass0_integral` に到達する（[deltastar_loop.py:1918](/home/sano/work/forge-cardmesh/design/forge_design/feedback/deltastar_loop.py:1918)）。一致すれば、この固定入力について呼び出し経路による差を除外できる。ただし、V6 とはコード以外に環境・変換器も異なるため、**「#6e が唯一の原因」とは確定できない**。また A＝B なら A−V6＝B−V6 は自動的に成立し、独立した追加証拠にはならない。 |
| 元の V1a を FAIL のまま残す | **採用** | [plan:536](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:536) の許容は `source.generator_sha256` のみ。`source.written` は事前登録の漏れとして記録し、元の判定を変更しない。新しい比較で除外する具体的なキーを先に登録し、別の検証結果にする。 |
| 比較台本をそのまま A/B に再利用する | **却下（Major）** | [compare.py:88](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/v1a_executor/compare.py:88) は事前登録の `stop_rule.cap` より広い `stop_rule` 全体を比較する。実際、旧呼び出しの `stop_rule=True` と CLI の既定値で `explicit`・`requested_source` が異なる。さらに部分文字列による除外は、配置パスである `initializer.two_coef.plan` と計測時間 `physical_wall.fit.build_seconds` を除外できていない。**比較対象・除外対象を完全なキーパスで列挙する**こと。(9) の `ok` も forge しか検査していないため、変換器の sha256 の存在確認を追加する（同:103）。 |
| 準備後のファイル変更を起動前に拒否する修正 | **採用** | [plan.py:186](/home/sano/work/forge-cardmesh/design/forge_design/campaign/plan.py:186) を [runner_axismach.py:3027](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:3027) が最初の起動前に呼ぶ。既存の準備と確定計画は通過し、メモリ上で期待ハッシュを変えると対象５ファイルすべての不一致を拒否した。 |
| 初期壁の係数チェックで穴が閉じたとする判断 | **要再検証（Major）** | `28945473:design/forge_design/campaign/recipes/wt_card_chain.py:315–328` は **NS 終了後**に照合する。旧台本は準備時に係数の存在を検査していた（[p5_ns.py:116](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0011_ns_v2/p5_ns.py:116)）。加えて当該比較式を実行すると、期待値 `1.0` に対する実値 `NaN` が不一致リスト `[]` となった。準備後・最初の起動前へ移し、型・有限性・値を検査する。正常値は通過し、欠落・不一致・NaN は起動０回で拒否することを確認する。 |

結論: 元の V1a＝FAIL を保持し、比較項目を明確に事前登録したうえで、同一の現行コード・環境による旧呼び出し対 CLI の prepare-only A/B を１組行う。

第 1 仮説: 座標・節点表・初期壁数値の差の主因は、呼び出し経路ではなく形状生成側の版差であり、`#6e` が最有力候補である。確度: 中  
根拠: `eabbab97` は終端 C⁺ の選び方を最近傍から補間へ変更している（[moc_inverse.py:812](/home/sano/work/forge-cardmesh/design/forge_design/geometry/moc_inverse.py:812)）。実ファイルで、実行役側と `case/71.m55_combustion_wt/run_0014_chain_v2/` の `x_E` はそれぞれ `40.633257251889965`、`40.63079680586963`。`L_c` も同じ `0.0024604460 r_t` だけ異なる。変換前の `nozzle.msh` も異なるため、変換器だけを原因とする説明は成立しない。  
反証条件: コード・環境・参照入力を固定した旧呼び出しと CLI で、形状・格子・初期壁の数値差が残れば、「呼び出し経路に依存しない版差で説明できる」という仮説は棄却する。A/B が一致しても `#6e` 単独原因の確定にはならない。

第 2 仮説: FAIL の一部は比較仕様・実装による差である。確度: 高。`source.written`、計測時間、文書の移動先、停止規則の指定元を、数値差と分けて扱う必要がある。

第 3 仮説: `ic_from` と固定 `ic_file` などの経路差が生成入力に影響している可能性は未確認。今回の A/B で検査する。

判別 A/B: 変える因子は呼び出し経路だけ。A は旧台本の `run_pass0_integral(..., ic_from=…, stop_rule=True, prepare_only=True)`、B は実行役の CLI。コード・Python と依存環境・変換器・数値に効く環境変数・問題 YAML・Euler 参照・IC と継承元の記録を固定する。A が選ぶ最大番号の場も `res_90000.h5`、sha256 `c64a73f1…` であることを確認する。各準備１回、forge は０ step。既存 B は固定条件の証跡が揃う場合だけ再利用する。

見る量は全段設定、BC、基準の実効内容、節点表、座標・トポロジ・msh、初期壁数値、継承値。両方とも準備直後なので、補間後の初期場も比較する。時刻・計測時間・配置パスは完全なキーパスで限定し、`stop_rule.explicit/requested_source` は差を別記して `cap` の一致と混同しない。

→ **数値入力の差０なら**、この入力について呼び出し経路差を除外し、別検証として同等性を記録する。  
→ **数値入力に差が残れば**、版差だけによる説明を棄却し、最初に差が生じた生成段階を調べる。

やらない方がよいこと: 元の許容範囲を後から広げて V1a を PASS にすること、A を「V6 の再現」と呼ぶこと、今回の切り分け前に GPU 計算や AWS のビルド再現へ進むこと。`#6e` の巻き戻し比較は、単独原因まで確定する必要が生じた段階の追加検証でよい。

呼び出し側の前提への異議: **Minor** — 「トポロジも違う」は誤り。既存の `h5_mesh_hashes` で再計算し、両者とも `df45a7db…` と一致した。異なるのは座標と msh。節点表は実行役側 3,830 行、V6 側 3,831 行で、集合差はそれぞれ６個・７個だった。「(5)〜(7) はすべて #6e」と一括して扱わず、数値差・来歴差・比較器の問題を分けること。

不足情報: `#6e` 単独原因の確定には変更を制御した追加証拠が必要。また #6d 全体の完了には、旧台本の検査との対応表と負例の実行結果が必要。いずれも今回の次の一手の判断は妨げない。

証拠: 足りる
