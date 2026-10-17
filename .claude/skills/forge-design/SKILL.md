---
name: forge-design
description: 要求カードから承認・recipe・判定・成果物までを実行役 (`forge-design run`) で回す作業の規範 (plan tooling-design-orchestrator)。レシピのある機種 (最初は風洞ノズル、recipe axismach.wt_card_chain/v1) の設計を「承認して回す」「止まった所から再開する」「状態と集計を見る」「直して要求を出す (fixer)」「結果を報告する」ときに使う。要求の聞き取りと問題 YAML の新規の作成は skill design-intake、ノズルの作り方・道具の詳しい説明は skill nozzle-design。
---

# /forge-design — 要求カードから実行役で回す

正本は plan [`tooling-design-orchestrator`](../../../plans/active/tooling-design-orchestrator.md) (§4.0 原則・§4.2 判定・§4.4 直す範囲・§4.5 実行・§4.7 承認)。
この skill は薄く保つ。中身は 5 つだけ。

## 1. 変えてはいけない 4 つ (§4.0)

1. **判定とその条件・公差** (承認の記録の qualification・停止の基準・品質の閾値・`judge_run` の仕様)。
2. **ユーザが決めた欄** (カードの要求・`mesh_limits`・`budget`・判定器の無い要求; 出所 `sources` が user のもの)。変えるならユーザに承認を取り直す。
3. **予算** (承認の GPU 時間・ディスク)。
4. **既存の run と成果物** (上書き・削除しない。直したら新しい試行・新しい出力先)。

## 2. 道具の索引

[`design/CAPABILITIES.md`](../../../design/CAPABILITIES.md) の「8. 設計の実行役」。各道具の `--help` と冒頭の説明に入力・出力・終了コード・所要時間がある。

## 3. カードを作る (出所つき)

- 版 2 のカード (`format: forge.wt_card/v2`) に `sources` (欄ごとの出所: user は依頼の原文の参照、recipe_default は既定) を書く。
  節点数の上限は `mesh_limits` (overrides の `max_nodes` は書けない)。判定器の無い要求は `requirements_without_gate`。
- `python -m forge_design.intake.approval preview <card> <problem>` の表をユーザに見せ、ユーザの言葉の参照つきで `write` (承認の記録)。
  AI が承認の記録を書き換えない (sha256 で検査される)。

## 4. 回す

- `design/bin/forge-design run <campaign.yaml> --role fixed|fixer` (再開は `--resume`)。役割はユーザが決める。
  - `fixed` (決まった流れを回す役): 段が止まったら直さずに止まった理由と記録をユーザに見せる (手順書の「自分で直さない」と同じ)。
  - `fixer` (直してよい役): §4.4 の範囲 (承認済みの要求・判定条件・予算を保ち、既存のデータを変えず、新しい試行で行う設定の調整・格子の作り直し・
    case の中の探索の台本の変更) で直し、`forge-design note <campaign.yaml> --path <直した入力> --reason "<1 行>"` を出してから `--resume`。
    それより外 (ソルバ・判定器・共通の道具・承認の条件の変更、原因や結果の解釈の確定) は AGENTS.md のエスカレーション条件のとおり。
- 状態と集計: `forge-design status <campaign.yaml>` (試行の数・人に返した回数・直した回数・消費した予算)。
- forge を直接起動しない。AWS の実行担当 (`--executor aws`) は plan §5.1 #7 の後。

## 5. 報告する

- 結果は判定の記録で言う: `WT_QUALIFICATION.json` の `result` (accepted / back_to_correction / rejected / incomplete) と、
  解析の成立 (`judge_run` の回し直し; 6 項目と壁解像)。`incomplete` を「達成」と言わない。
- 何を試し、何を直し (note の理由)、どの判定で通った/止まったかを、判定の出力 (VERDICT の行) とファイルの名前を添えて書く。
- case README の「計算 run 一覧」は実行役が台帳から同期する。手で書き足すときは自分の run の行だけ。
