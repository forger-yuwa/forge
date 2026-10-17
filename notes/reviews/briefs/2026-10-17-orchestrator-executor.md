# 諮問ブリーフ: 実行役 (§5.1 #6) の作り — driver と試行の 2 層の service・台帳・予算・recipe の段・役割 (2026-10-17)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件 1: plan §4 に設計 (§4.5.1) を新しく書いた。実装の前。
- plan: `plans/active/tooling-design-orchestrator.md` §4.0・§4.3・§4.4・§4.5・**§4.5.1 (今回の案、作業ツリーに未 commit)**・§4.7.1・§4.1.1・§5.1 #6・#7・#9・§6 V1a・V3・V4
- commit: 7cc17b9f + 作業ツリー (§4.5.1 の追加だけ)
- ユーザ: #5 (forge-gate の統合) より #6 を先に (2026-10-17「じゃあ#6からの流れで」)。

## 観測事実 (済んだ部品)

- 試行ごとの service の契約は手元・AWS B で確かめた (`notes/investigations/2026-10-17-attempt-service-ab/`: 実行担当の SIGKILL の後 0.10 秒で全子孫が消える; Linger=no では ssh を切ると約 12 秒で止まる)。
- 承認の記録・判定専用の環境・中間の判定の再検査・run の物理の照合 (#2・#3): `design/forge_design/intake/approval.py`、`design/forge_design/feedback/gate_env.py`。
  run_0014 の写しで正式な入口の判定が accepted (`design/tests/run_gate_deliverables_tests.py` (P4))。
- 初期場の固定・実行の計画・campaign (#4): `design/forge_design/campaign/ic.py`・`plan.py`・`campaign.py`。run_0014 の写しで確定した計画の段の設定が実際の段の記録と一致、
  S1 の CFL だけを変えると起動 0 回 (`design/tests/run_execution_plan_tests.py`)。
- 今の手順書 (`procedures/design-chain-quickstart.md` §2〜§9) の流れ: `wt_card` → `hand_mesh.py build/check` → `runner_axismach` (Euler) → `deltastar_loop --init-integral`
  (Euler の参照を sha256 で固定、p1 以降は前の NS から `--ic-from`) → `contur_calibrate` → カード `card_p1.yaml` → 2〜5 をやり直す → `wt_qualify` → `deliverables` → `wt_qualify`。
- `deltastar_loop.run_pass0_integral` (`deltastar_loop.py:1708`) が `prepare_ns` → `run_staged_ns` を呼ぶ (計画の確定はこの間に入れる)。

## 案 (§4.5.1)

- driver (campaign 1 つ、user service) が台帳の唯一の書き手で recipe の段を順に回す。試行 (段の試行 1 つ、別の user service) が段の処理と予算の監視を持ち、自分の状態と停止の記録だけを書く。
  CLI は driver を起動するか読むだけ。再開は driver を起動し直し、unit の状態と完了の証拠で照合する。
- recipe の段: problem → mesh → euler → ns → qualify → (calibrate → 次の pass) → deliverables → qualify。合格の印は各道具の出力。
- 予算: GPU 時間 = forge の経過時間、ディスク = run の大きさ。予約 (残り全部を段の上限に)・監視 (超えたら service を止めて停止の記録 budget)・精算。
- 役割: fixed は recipe が作った入力だけ (作り直した sha256 と照合)、fixer は止まった後に直して `--resume`、直した入力と理由を台帳に。
- 分け方: #6a 台帳・記録・予算の勘定、#6b 試行の service と予算の監視 (障害の試験)、#6c recipe の段・driver・CLI・README、#6d §6 V1a。

## 問い

1. 2 層の service の分け方に穴はあるか。特に、driver が死んだ後の試行の終わりの取り込み、driver の二重の起動 (二つの `--resume`)、試行の service の名前の衝突、
   ログの置き場所。
2. 予算の勘定 (GPU 時間を forge の経過時間で数える、段の上限を残り全部にする、超えたら service を止める) は、ユーザ決定 (依頼ごとに見積もって承認、上限は一律にしない) に照らして妥当か。
   見積もりとの差 (例: 見積もりの 1.5 倍で人に聞く) を入れるべきか。
3. 役割 fixed の「recipe が作った入力だけ」の照合 (問題 YAML・格子・試行の入力のカードを前の段の出力から作り直した sha256 と比べる) は実装できる粒度か。より単純で十分な代わりはあるか。
4. #6 の合格の条件 (§5.1 #6 に書いた障害の試験・予算の超過・二重の起動・既存の出力先・停止の記録の欄・README の行) に足すべきもの、#6a〜#6d の順に問題はあるか。
5. GPU を使う通しの確認 (§6 V1b・V4 (a)) は #7 (AWS の idle 自動停止) の後。手元の GPU で V4 (a) を先にやってよいか (手元は 09-27・10-14 に落ちた記録がある; メモリ「Prefer AWS for runs」)。

## 読んでよいファイル

- `plans/active/tooling-design-orchestrator.md` (作業ツリー)
- `design/forge_design/campaign/`・`design/forge_design/intake/approval.py`・`design/forge_design/feedback/gate_env.py`・`design/forge_design/feedback/deltastar_loop.py` (1640〜1810 行)
- `design/forge_design/evaluate/runner_axismach.py` (`run_staged_ns`・`staged_ns_launch_configs`)、`design/forge_design/evaluate/stop_rule.py` (1860〜1990 行)
- `procedures/design-chain-quickstart.md`、`notes/investigations/2026-10-17-attempt-service-ab/README.md`
- `case/71.m55_combustion_wt/_runs/run_0011_ns_v2/p5_ns.py`・`p5_all.sh` (旧の台本)
