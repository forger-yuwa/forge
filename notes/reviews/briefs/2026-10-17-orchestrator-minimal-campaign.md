# 諮問ブリーフ: §5.1 #4 (親 plan の最小の契約) の範囲 — campaign・execution_plan・入力の固定・全段の照合 (2026-10-17)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件 1: plan §4 に設計 (§4.1.1) を新しく書いた。実装の前。
- plan: `plans/active/tooling-design-orchestrator.md` §4.1・**§4.1.1 (今回の案、作業ツリーに未 commit)**・§4.5・§5.1 #4・#6・§6 V1a
- 親 plan: `plans/active/tooling-design-problem-campaign-recipe.md` §4.1〜§4.4・§5.1 #6・#7
- commit: 45c01e6c + 作業ツリー (§4.1.1 の追加だけ)

## 観測事実

- #2・#3 は実装済み (§5.1 の結果の欄): 承認の記録 (`design/forge_design/intake/approval.py`)、判定専用の環境 (`design/forge_design/feedback/gate_env.py`)、
  中間の判定の再検査、run の本段の物理の照合 (`approval.check_run_physics`)。case/71 run_0014 の写しで、正式な入口の判定が accepted
  (`design/tests/run_gate_deliverables_tests.py` (P4))。
- 親 plan の campaign・adapter・`resolve_evaluation`・`execution_plan.json` はコードに無い (`grep` で 0 件)。
- `runner_axismach.prepare` は `ic_from` (ディレクトリ) の最大番号の `res_*.h5` を選ぶ (`runner_axismach.py:1031〜1033`)。準備の記録 (`prepare_info.json`) は
  `ic_from` のディレクトリしか持たず、使ったファイル・step・sha256 を持たない (case/71 run_0014 の `prepare_info.json`)。
- 段の記録 (`stage_manifest.json`) は段ごとに `tag`・`restart_from`・`key` (実効の設定の 18 キー)・`cfg_fnv` を持つ (run_0014: `S1_soft`・`S2_mid`・`main`)。
- 起動の記録 (`RUN_PROVENANCE.txt`) は forge の実行ファイルのパスと sha256・git_head・環境変数 (`FORGE_CUDA_BLOCKSIZE` など) を持つ。
- 改訂した方針 (§4.0 の原則 2): 仕組みは同じことを繰り返すようになってから作る。

## 案 (§4.1.1)

- (a) search の無い campaign: 名前・承認の記録・recipe の名前と版・case のディレクトリ・台帳の置き場所だけ。
- (b) `execution_plan.json` (試行 1 つに 1 つ、段の準備の前に書く): 固定した入力 (カード・問題 YAML・承認の記録・格子の msh・初期場のファイル (パス・step・sha256)・Euler の参照・δ の表・
  forge と変換器の sha256・数値に効く環境変数)、予定の段の列、予定の判定。
- (c) runner に初期場のファイルを渡す引数を足し、準備の記録にパス・step・sha256 を書く。`ic_from` は開発用に残す。
- (d) 走った後の照合 (`verify_plan`): forge の sha256・格子のハッシュ・初期場の sha256・段の列が計画と同じ。全段の設定は、本段 = #3、前の段 = 入口・壁・軸対称・熱力学・輸送が本段と同じ
  (乱流・粘性の有無・数値の設定は段階起動で変わってよい)。
- 作らない: adapter の 3 関数、探索のスキーマ、副作用の無い `resolve_evaluation` (探索の依頼が出てから親 plan で)。

## 問い

1. この範囲で、実行役 (#6) と §6 V1a (決定的な入力の比較・来歴の照合・台本の検査が関門として存在) に足りるか。足りないもの・過剰なものは何か。
2. (b) の計画を「準備の前に書く」形で、副作用の無い `resolve_evaluation` を作らないことは、親 plan の codex の指摘 (計画を各 pass・各段の直前に確定、実行は計画だけを読む) と矛盾しないか。
   矛盾するなら、最小で何を足せばよいか (例: 実効の solverConfig の文字列を準備の後・起動の前に計画に足す)。
3. (d) の前の段の照合 (変わってはいけないもの / 変わってよいもの) の分け方は妥当か。recipe が許す段の列をどこで持つべきか。
4. 親 plan の #6・#7 の残りを「探索の依頼が出たら」に回すことを、親 plan にどう書くべきか (親 plan は別の方針の文書で、所有者はこのセッション)。

## 読んでよいファイル

- `plans/active/tooling-design-orchestrator.md` (作業ツリー)、`plans/active/tooling-design-problem-campaign-recipe.md`
- `design/forge_design/evaluate/runner_axismach.py` (930〜1060 行 `prepare`、`prepare_ns`、2950〜3010 行 段と段の記録)
- `design/forge_design/intake/approval.py`、`design/forge_design/feedback/gate_env.py`
- `case/71.m55_combustion_wt/run_0014_chain_v2/` の `prepare_info.json`・`stage_manifest.json`・`RUN_PROVENANCE.txt`
- `procedures/design-chain-quickstart.md` (今の手作業の流れ)
