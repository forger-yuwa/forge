# codex 諮問 (diagnose): orchestrator-minimal-campaign

- **brief**: [`notes/reviews/briefs/2026-10-17-orchestrator-minimal-campaign.md`](../../notes/reviews/briefs/2026-10-17-orchestrator-minimal-campaign.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-17
- **commit**: `45c01e6c` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 2.3 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-orchestrator-minimal-campaign-diagnose.codex.log` (codex: ok)
- **結論**: §4.1.1を「準備前の入力固定 → 準備後・各段起動前の実効計画確定 → 計画からの実行と前後照合」に修正し、既存の `ns_stage_configs` を再利用したGPUなしの検証を先に行う。
- **extra**: `plans/active/tooling-design-problem-campaign-recipe.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 指摘 | 採否 | 根拠・理由と対案 |
|---|---|---|
| **Major 1：準備前の入力一覧だけでは、実行計画の契約を満たさない** | **却下**：(b)(d)をそのまま実装する案 | 親 plan は実効設定の確定・計画だけからの実行・各段直前の確定を要求する（[親 plan:92](/home/sano/work/forge-cardmesh/plans/active/tooling-design-problem-campaign-recipe.md:92)）。現在の runner は起動時に設定を読み直し、さらに本段の `nStepOuter` を停止規則の上限へ変更する（[runner_axismach.py:2946](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2946)、[同:3002](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:3002)）。**対案：準備前は入力・生成予定を記録し、準備後・起動前に実効設定と生成物を確定する。各段は確定した設定を消費し、起動直前と終了後に照合する。** |
| **Major 2：「本段との物理の一致」と「予定した全段との一致」は別** | **却下**：(d)の照合範囲と包括的な変更許可 | `check_run_physics` が読むのは現在の `solverConfig.yaml` と `bcondConfig.yaml`（[approval.py:494](/home/sano/work/forge-cardmesh/design/forge_design/intake/approval.py:494)）。`run_0014_chain_v2/stage_manifest.json` の各段の `key` は18項目で、軸対称・輸送物性全体などを含まない。現行 NS の段生成は物理を維持して次数・CFL等を変更しており、乱流・粘性の任意変更まで許す必要はない（[runner_axismach.py:2668](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2668)）。**対案：版付き recipe が段列と許す設定差を所有し、各段の設定全文・境界条件・外部参照を計画と照合する。承認条件との照合も別に残す。** |
| **Major 3：初期場ファイルの固定だけでは入力が閉じない** | **採用**：(c)。ただし付随入力の固定を追加 | 最大番号の snapshot 選択は実在する（[runner_axismach.py:2523](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2523)）。加えて、初期場の元 run の設定・`limiter_ref.json` から基準値を継承する経路がある（[同:1777](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:1777)）。snapshot のハッシュだけでは、この依存を固定できない。**対案：継承値を解決して計画に格納し、出所も固定する。親 plan が要求する形状・メッシュ・輸送種配置・エネルギー基準の適合性を検査する。準備後の初期場も、元 snapshot とは別の生成物として識別する。** |
| **Minor 1：親 #6・#7 を一括して「探索待ち」にすると必要な契約まで延期される** | **採用**：探索機能・汎用 adapter の延期。**却下**：#6・#7全体の延期 | 親 #6・#7には未知キー拒否、入力参照、全段照合も含まれる（[親 plan:181](/home/sano/work/forge-cardmesh/plans/active/tooling-design-problem-campaign-recipe.md:181)）。また親の search 無し campaign は `evaluations` を持つため、今回の簡略形式との対応が必要（[同:75](/home/sano/work/forge-cardmesh/plans/active/tooling-design-problem-campaign-recipe.md:75)）。**対案：#6a/#7aを今回の限定実装、#6b/#7bを探索等の需要待ちに分割する。未対応の `search`・未知 recipe・未知キーは明示的に拒否する。** |

結論: §4.1.1を「準備前の入力固定 → 準備後・各段起動前の実効計画確定 → 計画からの実行と前後照合」に修正し、既存の `ns_stage_configs` を再利用したGPUなしの検証を先に行う。

第 1 仮説: 現案の不足は探索機能ではなく、確定した計画と実際の起動入力を結ぶ契約の欠落である。確度: 高  
根拠: 上記 Major 1・2。`ns_stage_configs` は既に forge なしで各段の実効設定文字列を返すため、runner 全体を分解する前に利用できる（[runner_axismach.py:2676](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:2676)）。  
反証条件: 現案のままで、全段の実効設定・参照入力・実行環境が起動前に確定し、確定後の相違を起動前に拒否する経路と正負例が示されれば、この不足の指摘は撤回する。

第 2 仮説: snapshot 固定後も、元 run からの付随値の再解決によって入力が変わり得る。確度: 高。根拠は Major 3。ただし、`run_0014_chain_v2` で実際に変動したという主張ではない。

判別 A/B: CFD **0 step**、起動を記録するダミーを使う。計画確定後、Aは無変更、Bは起動対象の `S1_soft` 設定のCFLだけを計画値から変更する。見る量は照合結果・不一致パス・ダミー起動回数。**Aが通り、Bが起動回数0で止まれば**計画との結び付けが成立する。**Bが変更後の設定で起動すれば**第1仮説を支持する。Aも止まれば試験不成立。このA/Bは実装後に行う提案であり、今回実施した結果ではない。

やらない方がよいこと: 起動後に実行済み設定を計画へ追記して「一致」とすること、`StageManifest.key` だけを全設定の証拠とすること、汎用 adapter と探索スキーマを先に作ること。数値設定の調整は試行計画を確定する前に行い、確定後の無記録な変更を許さない。

呼び出し側の前提への異議: 「副作用のない `resolve_evaluation` 全体を今は作らない」は採用できるが、親 plan の現行文面とは変更になる。親 §4.3に、今回の限定実装では準備処理を先行させ、生成済み入力から各段の計画を確定する例外と、その照合条件を明記すること。「設定文字列を足す」だけでは足りず、実行側がその確定内容を使う必要がある。将来生成される壁・格子・段間 snapshot は生成元を参照し、消費する段の直前に実体へ結び付ける。1 runに1ファイルでも、各段の確定済み記録を上書きしなければよい。

また、V1aはこの修正だけでは完了しない。[子 plan:439](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:439)が挙げる基準・節点表・`RECIPE.json` 等の比較と、旧台本の検査から実行器の関門への対応表が必要である。合格条件は、許可した正規化以外の入力差0、未対応の検査0、登録した負例の通過0とする。

不足情報: V1a実施時には `p5_ns.py`・`p5_all.sh` の固定版と関門の対応表が必要。今回の設計判断を妨げない。ファイル変更・forge起動・流れ場の合否判定は行っていない。

証拠: 足りる
