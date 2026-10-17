# 諮問: §6 V3 (依頼の課題集、GPU なし) の事前登録の案と、§5.1 #13 の分け方 (2026-10-18)

plan `plans/active/tooling-design-orchestrator.md` §6 V3・§5.1 #13。#2〜#9 は済んだ (#8 は V2a〜V2c PASS、#9 は CLI・skill・AGENTS.md)。

聞きたいこと: (Q1) 下の 13a の事前登録の案 (機械の側だけ) に穴が無いか、(Q2) モデルを回す側 (13b) を費用の承認の待ちに分けてよいか、
(Q3) 「予算が足りない依頼」の扱い (下の観測事実 3)。

## 1. 観測事実

1. §6 V3 の文面 (2026-10-17 codex で決めた): 正常 3 件 (case/71・試走 A・case/70 M5 の空気) と止まるべき 3 件 (公差が届かない・ユーザ欄の欠け・予算が
   足りない)。GPU なしの到達点はカード → 問題 YAML → 格子 (品質の VERDICT)。流れの公差の未達は固定した既存の不合格のデータの再判定で扱う。安価な負例
   (11 種) は AI を使わずプログラムの検査だけで全部止まる。正例 (承認したカード、格子・係数だけを変えた試行) は止まらない。回し方は賢い AI (`--role fixer`) と
   知識の少ないモデル (Haiku、使えればローカル LLM; `--role fixed`) で各依頼 3 回。合格は機械の誤合格 0・安価な負例が全部止まる。
2. 安価な負例の多くは既存の試験にある: 根拠の無いユーザ欄・出所の未知のパス・既定との不一致・予算の欠落・`requirements_without_gate` の削除・CLI での公差の
   上書き (`run_approval_tests.py`)、外部の Python・import の元の違い (`run_gate_env_tests.py`)、総合 PASS と個別の FAIL の矛盾 (`run_wt_qualify_tests.py` 4c)、
   予算の実行中の超過・再開による二重の起動・既存の出力先 (`run_driver_tests.py`)。
3. recipe v1 (`design/forge_design/campaign/recipes/wt_card_chain.py:203`) は全部の段で `estimate: {}` を返す。driver の起動の前の予算の検査
   (`budget.check_start`) は見込みが空なので止めない。予算の超過は試行の service の監視が実行中に止める (GPU の段が走ってから)。したがって「予算が足りない
   依頼」は GPU なしでは起動の前に止まらない (Euler の段を起動しようとする)。
4. モデルを回す側: 依頼 6 件 × 2 モデル × 3 回 = 36 回のエージェントの実行 (サブエージェントの固定の費用だけで 1 回 約 4.7 万トークン; AGENTS.md の実測)。

## 2. 案

**13a (機械の側、今やる)**: 依頼の原文・期待する結果・固定する条件を `design/tests/fixtures/v3_requests/` に置いて commit し、`design/tests/run_v3_requests_tests.py` で回す。

| 依頼 | 機械で回す範囲 | 期待 |
| --- | --- | --- |
| N1 case/71 (版 2、mesh_limits {443000, 467000}) | 承認 → `wt_card` → `hand_mesh.py build --approval` → `check` | 承認の preview の問題 0、NS・Euler の品質 PASS |
| N2 試走 A (case/90 card_p0、mesh_limits {412000, 434000}) | 同上 | 同上 |
| N3 case/70 M5 の空気 (既定の上限) | 承認 → `wt_card` → `mesh_stage.check_template` | NS・Euler の品質 PASS |
| S1 公差が届かない (試走 A に `qualification.exit_core_M_tol_rel` を厳しく) | 固定の不合格のデータ = 試走 A `run_0003_ns_p0` を `wt_qualify` で再判定 | `back_to_correction` (accepted にならない)、反復の上限なら `rejected` |
| S2 ユーザ欄の欠け (`sources` から `wall` を除く) | 承認の preview | 問題が出て承認を書けない (問題 YAML・格子を作らない) |
| S3 予算が足りない (`budget.gpu_hours` 0.01) | 下の (Q3) | 起動の前に止まる |

安価な負例: 既存の試験の該当の検査を名前で列挙して集め、全部が「止まる・拒否・incomplete」であることを 1 本の集計で確かめる (新しい検査は書かない)。
正例: N1〜N3 と、既存の試験の「格子・係数だけを変えた試行が通る」。

**13b (モデルを回す側、費用の承認の待ち)**: 13a の依頼の原文を使い、Haiku (`--role fixed`) と賢い AI (`--role fixer`) で各 3 回。数えるもの・合格は §6 V3 のとおり。
費用の見積もりを付けてユーザに諮る。

**(Q3) の案**: recipe の各段に GPU 時間・ディスクの見込みを持たせる (実測から: case/71 の Euler 約 4 分・NS 約 1〜3 時間、ディスクは run の大きさ)。
driver の起動の前の検査が見込みで止める。見込みの値は保守的に (大きめに) し、出所を recipe に書く。13a の S3 は「Euler の段の前で予算で止まる (試行を
記録しない)」を期待にする。
