# codex 諮問 (diagnose): orchestrator-mesh-limits

- **brief**: [`notes/reviews/briefs/2026-10-18-orchestrator-mesh-limits.md`](../../notes/reviews/briefs/2026-10-18-orchestrator-mesh-limits.md)
- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **date**: 2026-10-18
- **commit**: `0161e51c` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-18-orchestrator-mesh-limits-diagnose-mesh-limits.codex.log` (codex: ok)
- **結論**: 要求条件の一致・承認上限・両格子の実効点数を生成前に照合する共通処理へ設計案を修正し、その接続を検証する A/B を先に事前登録する。
- **extra**: `plans/accepted/tooling-nozzle-mesh-size-guard.md`, `case/71.m55_combustion_wt/_mesh/hand_mesh.py`, `design/forge_design/intake/approval.py`, `design/forge_design/campaign/recipes/wt_card_chain.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論は **修正して採用**です。Q1 は現案のままでは抜け道が残ります。Q2 の変更前後 A/B、Q3 の段の関数としての検証は、保証範囲を限定して採用できます。

確認対象は `HEAD 0161e51c`、提示プロンプト、および読み取り時点の未コミット差分です。ファイル変更・forge 起動・格子生成は行っていません。

| 対象 | 採否 | 根拠と対案 |
|---|---|---|
| (1) 承認から上限を取る | **要再検証 — Major 1** | 承認への参照だけでは、同じ依頼の試行カードとは保証できない。未コミットの [`approval.py:495`](/home/sano/work/forge-cardmesh/design/forge_design/intake/approval.py:495) は参照と `mesh_limits` を照合するが、要求条件を比較しない。**上限を返す前に、試行カードの保護条件も承認と照合する。** |
| (2) 手組み格子の自動拡大を廃止 | **採用。ただし Major 2 の修正が必要** | `hand_mesh.py` は Euler を `ni × nj_e` で数える。`axis_gap_frac` の導出点数を取り逃がす。**出力予定の mesh 節から実効の `Mesh2DParams` を作り、共通の `check_node_budget` を YAML 書き出し前に使う。** |
| (3) 問題 YAML の上限照合 | **採用。ただし Major 3 の修正が必要** | `問題の上限 ≤ 承認の上限` は必要。ただし照合するファイルと時点が不足している。**実際に消費する NS・Euler 各 YAML を、準備の直前にも照合する。** `check_run` にも共通検査を明示的に接続する。 |
| (4) v2 は `mesh_limits` に一本化 | **採用** | 上限の指定箇所を一本化する方向は妥当。値の転記は `estimate_mesh` の前に行う必要がある。同関数は [`wt_card.py:1093`](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:1093) で全格子を生成する。承認済み試行では、その呼び出し前に Major 1 の照合も必要。 |
| (5) 雛形の段の関数 | **採用、recipe v1 への追加は不要** | plan は雛形と手組みの二実装を定めるが、同一 recipe への収容は要求していない（`tooling-design-orchestrator.md:403`）。ただし生成手順のコピーを増やさず、実際の準備経路と処理を共有する。`limits` は検査済みの承認元に結び付ける。 |
| (6) 試験 | **要再検証** | 以下の V2a〜V2c の補強後に採用。特に「最後に例外が出た」だけで生成前停止の合格にしない。 |
| V2a | **採用 — Minor 1** | #8 による退行の検査として変更前後 A/B が適切。`0161e51c` の台本だけを固定し、依存コード・入力・`--axis-cap-frac 0.03`・変換器を共通にする。座標だけでなく接続・境界所属・変換後の不変部分も比較する。**V6 再現、V1a の FAIL 解消とは呼ばない。** |
| V2b | **採用、部品の検証として限定** | NS・Euler の品質 PASS と点数だけでなく、実際の準備経路との格子同一性、承認元、問題 YAML と格子のハッシュの結び付きを確認する。case/70 の campaign 通しの成功は主張しない。 |
| V2c | **要再検証 — Major 1〜3** | 現負例に、①有効な承認参照を付けた要求変更カード、②Euler だけの超過と `axis_gap_frac`、③格子段の後に変更された Euler companion を追加する。正例として較正係数だけを変えた試行も通す。負例は生成関数の呼び出し 0 回・出力 YAML 作成 0 件を確認する。 |

**Major 1：承認の参照と、承認条件の一致を混同している。**

書き込みなしの関数単体確認で、I/O をメモリ上の fixture に置き換え、実際の検証関数を呼びました。試行カードの差は入口全圧だけです。

| 試行 | `conditions.Pt` | `mesh_limits_for_card` | 完全な `approval.check` |
|---|---:|---|---|
| A：承認条件と同じ | 3,500,000 Pa | 443,000／467,000 を返す | `OK` |
| B：全圧だけ変更 | 7,000,000 Pa | **同じ上限を返す** | `MISMATCH`: `requirements.conditions.Pt` |

これは現実のユーザ承認を作った試験ではなく、関数の検査範囲を確認したものです。`approval.check` は保護条件を比較しますが（`approval.py:528`）、上限取得関数は比較しません。

さらに recipe の承認照合は `build`・`check` の後です（[`wt_card_chain.py:290`](/home/sano/work/forge-cardmesh/design/forge_design/campaign/recipes/wt_card_chain.py:290)）。後で拒否しても「格子を作る前に止まる」という V2c は満たせません。`hand_mesh.build` 内の `estimate_mesh` も全格子を作るため、その前に照合する必要があります。

**Major 2：指定点数の積では、実効点数の上限を守れない。**

実際の `mesh_params_euler` と `node_count` を使った、格子配列を作らない確認結果です。

| 条件 | `ni × 指定 nj` | 共通関数が数えた実効点数 | 承認上限 467,000 |
|---|---:|---:|---|
| `axis_gap_frac` なし | 4810 × 97 = 466,570 | 466,570 | 上限内 |
| `axis_gap_frac: 0.005` | 4810 × 97 = 466,570 | **966,810** | 超過 |

未コミットの `hand_mesh.py` の修正も、上限との比較には `ni × nj_e` を使用しています。既存の実効点数検査は YAML を書いた後なので、そこで止まっても「問題 YAML を書かずに停止」は満たしません。[`runner_axismach.py:812`](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:812) の共通処理を、出力予定の両ブロックに対して先に適用してください。

**Major 3：照合済み NS YAML と、実際に使う Euler YAML は別物である。**

recipe の格子段は `sd["problem"]` を照合しますが、Euler 段は `sd["problem_euler"]` を渡します（`wt_card_chain.py:296`・`:302`）。入力ハッシュの変更も、理由付き `note` があれば通る設計です（`:87`・`:100`）。

したがって fixer が Euler companion の上限を変更して `note` を出す経路では、ハッシュ検査だけでは承認上限を保護できません。runner は YAML 自身の上限しか知りません。**`note` の受理と承認条件の照合を分け、実際に使う YAML を準備直前に再検査する**必要があります。また現在の `approval.check_run` は `check` を呼ばないため、`check` への追加が自動的に `check_run` に効くという前提も成立しません。

結論: **要求条件の一致・承認上限・両格子の実効点数を生成前に照合する共通処理へ設計案を修正し、その接続を検証する A/B を先に事前登録する。**

第 1 仮説: **現案では、承認参照を持つ要求変更カードが、完全な承認照合より先に格子生成へ到達する。** 確度: 高  
　根拠: 上記関数単体確認、および `wt_card_chain.py:290–296` の処理順序。  
　反証条件: 要求変更カードが、`estimate_mesh` 内を含む最初の全格子生成より前に拒否され、正当な試行カードだけが生成入口へ到達する。

第 2 仮説: **手組みの点数比較が指定 `nj` に依存すると、Euler の導出点数を取り逃がす。** 確度: 高。466,570 対 966,810 の数えで確認。

第 3 仮説: **格子段後の Euler companion の変更が、理由付き `note` 経由で承認上限の再照合を通らず準備へ渡る。** 確度: 中。コード経路を確認、統合再現は未実施。

判別 A/B: **変更因子は試行カードの `conditions.Pt` だけ。** A は承認値、B はその 2 倍。承認参照・上限・他の入力は固定し、正式な格子段入口を各 1 回呼ぶ。最初の全格子生成関数を到達記録だけの代替にし、forge は 0 step。**A が到達し B が生成呼び出し 0 回・YAML 作成 0 件で拒否されれば、第 1 仮説を棄却。B も到達すれば支持。A も拒否されれば陽性対照不成立。**

やらない方がよいこと: 後段の `MISMATCH` を生成前停止の証拠にすること、V2a の合格で V6 との差を解消済みにすること、V2b の部品試験を campaign 全体の成功と扱うこと。

呼び出し側の前提への異議: 「`approval.check` と同じ規則」は参照照合だけでなく保護条件の一致まで含む。また、上限の転記と後段照合を二重に置くだけでは、生成前停止を保証できない。

不足情報: 実装完了後の統合 A/B と V2a〜V2c の実測結果。今回の設計修正の判断は妨げない。

証拠: 足りる
