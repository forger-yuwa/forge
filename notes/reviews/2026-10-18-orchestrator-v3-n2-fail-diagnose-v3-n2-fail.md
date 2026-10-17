# codex 諮問 (diagnose): orchestrator-v3-n2-fail

- **brief**: [`notes/reviews/briefs/2026-10-18-orchestrator-v3-n2-fail.md`](../../notes/reviews/briefs/2026-10-18-orchestrator-v3-n2-fail.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-18
- **commit**: `b1a02511` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 2.6 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-18-orchestrator-v3-n2-fail-diagnose-v3-n2-fail.codex.log` (codex: ok)
- **結論**: N2 の FAIL と元の上限を維持し、同一環境で旧版・現行版の数値生成経路だけを比較する A/B を次の一手とする。
- **extra**: `case/71.m55_combustion_wt/_mesh/hand_mesh.py`, `case/90.qs_trial_a/README.md`, `design/tests/fixtures/v3_requests/requests.yaml`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 案 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| (a) N2 の FAIL を残す | **採用・Major** | [fixture:24](/home/sano/work/forge-cardmesh/design/tests/fixtures/v3_requests/requests.yaml:24) の期待は NS・Euler の品質 PASS。[結果記録:12](/home/sano/work/forge-cardmesh/notes/investigations/2026-10-18-v3-13a/README.md:12) は生成前の停止なので、N2・13a は未合格のまま。「上限の保護が働いた」は別の評価として記録する。 |
| (b) 旧版と現行版を比較する | **採用。ただし比較条件を具体化・Major** | [試走の一次記録:5](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/qs_trial_a_evidence/README.md:5) では p0 は `a36a87ff`、p1 から `eabbab97`。同じカードでも生成コードは同じではない。Python・依存・外部入力を固定した、下記の A/B を先に行う。 |
| (c) NS 上限を 413,000 にした N2′ | **現案は却下・Major** | **Euler も不足する**。[hand_mesh.py:416](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_mesh/hand_mesh.py:416) は Euler の `ni` を NS とそろえる。`ni=4480`、`nj=97` なら **434,560 > 434,000**。N2′ を設ける場合は、両方の実効節点数を確認し、1,000 単位の切り上げなら `{413000, 435000}` を新しい提案値とする。旧 N2 の合格への読み替えには使わない。 |
| (d) 案内文の重複を直す | **採用・Minor** | [hand_mesh.py:440](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_mesh/hand_mesh.py:440) の `hint` と [runner_axismach.py:830](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:830) の定型文が重複する。了承を求める文は定型側に集約し、呼び出し側は修正対象と手順だけを渡す。数値判定は変えない。 |

結論: **N2 の FAIL と元の上限を維持し、同一環境で旧版・現行版の数値生成経路だけを比較する A/B を次の一手とする。**

第 1 仮説: 生成コードの変更が設計壁・物理壁または配点表を変え、`ni` の切り上げ境界をまたがせた。確度: **中**。  
根拠: 旧 [problem_p0.yaml:69](/home/sano/work/forge-cardmesh/case/90.qs_trial_a/problem_p0.yaml:69) は自動倍化後の `n_axis_inv: 4800`、同ファイルの手組み格子は `ni_min=4470`・`ni=4470`。その後の `eabbab97` は終端 C⁺ を補間する方式へ変更し、該当経路の網の自動倍化も止めている。現行の [ni 決定処理:306](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_mesh/hand_mesh.py:306) は壁際 AR と間隔の条件から探索し、10 の倍数へ切り上げる。ただし、**差の原因を `eabbab97` 単独に特定した証拠はまだない**。  
反証条件: 同じ環境・カードから旧版と現行版が同じ壁・表・`ni_min`・実効節点数を返す場合、コード版の差だけによる説明は棄却する。

第 2 仮説: 旧 AWS 環境と今回の環境の依存ライブラリ、外部入力または実効オプションの差が寄与している。確度: **低、未確認**。旧版の再現対照がないため除外できない。

判別 A/B: 変える因子は**数値生成コードの版だけ**とする。A=`a36a87ff`、B=`b1a02511`。カード原本 `card_p0.yaml`（sha256 `23caf0b7…`）、Python・依存、物性・係数の入力、`axis_cap_frac=0.03`、径方向比の上限 `1.2` を固定する。旧版も読める同一の v1 カードから、カード生成→壁→配点表→節点数決定まで各 1 回実行する。**forge は 0 step、全格子生成・変換・品質検査はこの A/B では不要**。旧版の上限自動拡大・出力処理には進めず、得た両方の節点数を共通の現行の数え方で元の上限と照合する。

見る量は `r_throat`、`M_knot`、`n_axis_inv`、壁の終点、第一層・密度の全表、丸め前の `ni_min`、丸め後の `ni`、実効 `nj`、NS・Euler の節点数。

→ **A が旧値 411,240／433,590、B が現行値 412,160／434,560 を再現すれば**、環境差を必要とせずコード版の変更で差を説明できる。ただし `#6e` 単独への帰属は保留する。  
→ **A と B が一致すれば**、コード版だけの仮説を棄却し、旧実行の環境・入力・オプションを候補へ戻す。  
→ どちらかが対応する記録を再現しなければ、対照不成立として原因を確定しない。

やらない方がよいこと: 上限超過の停止を N2 の PASS に読み替えること、NS 上限だけを増やして再試行すること、`ni=4470` を強制すること、原因を確かめず AR・第一層の条件を緩めること。

呼び出し側の前提への異議: 「160 点不足」は最初に発生した NS の拒否だけを表す。[hand_mesh.py:442](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_mesh/hand_mesh.py:442) は NS→Euler の順に検査するため、Euler の不足はまだ表示されない。読み取り確認では、旧問題をメモリ上だけで `ni=4480` にした共通 `node_count`・`check_node_budget` が、NS **160 点超過**、Euler **560 点超過**をそれぞれ返した。これは生成の再現や品質 PASS の証拠ではない。また、旧 p0 の生成日は YAML 上 **2026-10-16**であり、ブリーフの「10-17 に作った」は修正が必要。

不足情報: 原因の確定には旧実行時の Python・依存ライブラリ・実効コマンドと、今回の丸め前の `ni_min`・生成途中の表が欲しい。今回確認したカードと旧問題 YAML の sha256 は事前登録と一致した。ファイル変更・forge 起動は行っていない。

証拠: 足りる
