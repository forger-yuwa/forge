# 諮問ブリーフ: 承認した要求と判定条件の固定 (§5.1 #2) の設計 — カードの版 2・承認の記録・固定した版での判定 (2026-10-17)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件 1: plan §4 に設計 (§4.7.1) を新しく書いた。実装の前。
- plan: `plans/active/tooling-design-orchestrator.md` §4.0 (4 つと守る仕組み)・§4.2 (判定器と判定条件の固定)・§4.7・**§4.7.1 (今回の案、作業ツリーに未 commit)**・§5.1 #2・§6 V3
- commit: 9c52fdd4 + 作業ツリー (§4.7.1 の追加だけ)

## 観測事実

- `design/forge_design/feedback/wt_qualify.py`:
  - `resolve_tolerances` (165 行〜) は CLI → カードの `qualification` → 問題 YAML の `qualification` → 問題 YAML の `evaluate.settle` → run の `settle_criteria.yaml` → 既定 の順。
  - 54 行の docstring: 要求カードの書式 `forge.wt_card/v1` は `qualification` を受け付けない (生成器に通すカードには書けない)。
  - `deliverables_check` (742 行〜) は `DIMENSIONS.json` の `verdict` と `card_sha256` を読む。`optional_parts` は CLI (`--optional-parts`) から来る。
- `design/forge_design/export/deliverables.py`: `part_dimensions` (701 行〜) が `dimension_gate` の `step_dimensions`・`mesh_dimensions`・`chord_from_mesh`・
  `junction_from_mesh`・`evaluate_gate` で寸法の関門を作る。`MANIFEST.json` の `git_commit` は `git rev-parse --short HEAD` (作業ツリーの変更の有無は記録しない)。
- `design/forge_design/intake/wt_card.py`: `CARD_FORMAT = "forge.wt_card/v1"`、`_TOP_REQUIRED`・`_TOP_OPTIONAL` (111〜112 行)、`_OVERRIDE_FORBIDDEN` (115 行)、
  `validate_card` (278 行) は既定値を埋めた写しを返す。
- 前の諮問 (`notes/reviews/2026-10-17-orchestrator-direction-diagnose.md`): 判定器の sha の記録だけでは保護にならない、依頼の開始で承認した commit・依存環境・判定条件を固定し
  作業側から変更できない場所で判定する、CLI・問題 YAML の相違は拒否する、`DIMENSIONS.json` の判定を作った処理も固定した版に結び付けるか再評価する。

## 案 (§4.7.1 の要約)

- カードの版 2: v1 の欄 + 省略可の節 `qualification`・`budget`・`mesh_limits`・`requirements_without_gate`・`sources`。生成器は v1 と v2 を受け付ける。
- ユーザ欄の葉は `sources` で `user` (ref つき) か `recipe_default` でなければ起動の前に止める。既定で埋めた欄は `recipe_default`。
- 承認の記録 `card_approval.json`: カードの sha256・解決した後のユーザ欄の値の全部・`user_ref`・`gate` (HEAD の完全な hash・作業ツリーがきれい・Python の版・`pip freeze` の sha256)。
  道具が書き、AI はユーザに表を見せて承認を得てから書く (規約。強制は #14)。
- `wt_qualify --approval`: `gate.commit` の別 checkout (`git worktree add --detach`) で判定。HEAD と追跡ファイルの無変更・Python・`pip freeze` を照合。
  公差・窓・`optional_parts` は承認の値だけ、CLI の指定は値が同じでも拒否、問題 YAML の相違も拒否。寸法は固定した版で測り直す。`requirements_without_gate` があれば `incomplete`。
  `--approval` なしの今の使い方は残す。

## 仮説

- 承認の記録にカードの sha256 と解決した値を両方書けば、承認の後のカードの書き換えと、既定の値の黙った変更 (雛形の既定が後で変わる) の両方を検出できる。
- 別 checkout の版の照合 + `pip freeze` の照合で、作業側の作業ツリーの書き換えは判定に効かない。

## 問い

1. §4.7.1 の案に穴はあるか。特に、(a) ユーザ欄の範囲 (`recipe`・`overrides` を外したこと、評価の位置や反復の上限は入れるべきか)、(b) `sources` の粒度 (節と葉)、
   (c) 承認の記録を AI が書く (規約) ことの限界の書き方、(d) 別 checkout の検査 (HEAD・追跡ファイル・未追跡のファイルは見ない) で足りるか、
   (e) 寸法の測り直しで足りるか (MANIFEST の他の部品の判定はどうするか)。
2. `pip freeze` の sha256 は依存環境の固定として妥当か (重すぎ・軽すぎ)。
3. 実装の順と合格の条件 (§5.1 #2) に足すべき負例はあるか (今: 根拠の無いユーザ欄、CLI での公差・窓・必須成果物の変更、判定器なし、codex の判別 A/B: 承認値 0.001 に CLI 0.1 → 拒否)。

## 読んでよいファイル

- `plans/active/tooling-design-orchestrator.md` (作業ツリー; `git diff` で §4.7.1 の追加が見られる)
- `design/forge_design/feedback/wt_qualify.py`、`design/forge_design/intake/wt_card.py`、`design/forge_design/export/deliverables.py`、`design/forge_design/export/dimension_gate.py`
- `notes/reviews/2026-10-17-orchestrator-direction-diagnose.md`、`procedures/design-chain-quickstart.md` (カードの例、277 行)
