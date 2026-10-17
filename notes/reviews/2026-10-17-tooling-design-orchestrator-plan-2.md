# codex レビュー: tooling-design-orchestrator (plan)

- **plan**: [`plans/active/tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md)
- **stage**: `plan`
- **date**: 2026-10-17
- **commit**: `9436b084` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-tooling-design-orchestrator-plan-2.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1
- **extra**: `notes/reviews/2026-10-17-orchestrator-direction-diagnose.md`, `procedures/design-chain-quickstart.md`, `design/forge_design/feedback/wt_qualify.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes

目的と、既存の campaign／recipe・判定器を再利用する方向は妥当です。  
ただし、障害時の実行管理、メッシュ規模の承認、検証の期待結果と完了条件を実装前に修正してください。  
ファイル変更・forge 起動はしていません。以下の対案は plan 未反映です。

1. **Major — 実行器が落ちた場合の監視・予算保護・再接続の契約が不足しています。**

   根拠: [plan:233](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:233) は `run_case.sh` を包み、プロセスが生存していれば監視へ戻る設計です。しかし既存の [stop_rule.py:1901](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1901) は `start_new_session=True` で子を起動し、監視状態と `Popen` を呼び出し元が保持します。監視側が強制終了しても子は残り得ます。一方、同関数の再呼び出しは既存出力を拒否します（[stop_rule.py:1879](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:1879)）。PID の生存確認だけでは停止規則の監視を復元できず、再開されるまでの予算超過も防げません。

   **対案:** attempt ごとの実行担当プロセスに、runner・停止規則・予算監視をまとめて所有させ、CLI の再開はその担当へ接続する方式に絞ってください。担当自体が死亡した場合は子も停止し、新しい attempt に引き継ぐ契約にします。障害試験には「CLI 死亡後も予算監視が続く」「監視担当死亡後に子だけ残らない」「既存 runner を再呼び出しせず再接続する」を追加してください。

2. **Major — `hand_mesh.py` をそのまま包むと、別依頼にも過去のメッシュ拡大承認を流用します。**

   根拠: [plan:263](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:263) はそのまま包む方針ですが、[hand_mesh.py:354](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_mesh/hand_mesh.py:354) は点数から `max_nodes` を自動拡大し、理由に固定文「ユーザ決定『a』2026-10-15」を書きます。既存の [runner_axismach.py:790](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:790) は上限が正整数、理由が非空であることしか確認しません。新しいカードでも、今回の承認と無関係な文字列で上限を通せます。

   **対案:** メッシュ段の契約に、依頼ごとの承認済み節点数上限と承認元を追加してください。生成器が算出する必要点数と、ユーザが許した上限を分離し、格子生成前に照合します。case/71 の承認はその対象に限定し、「別カード＋過去の承認文」「承認上限を超える再生成」を負例にしてください。これは係数一般化の #10 を待たず、#2・#8 に必要です。[既存の承認規則:76](/home/sano/work/forge-cardmesh/plans/accepted/tooling-nozzle-mesh-size-guard.md:76) とも整合します。

3. **Major — V1b の期待結果が参照 run と逆です。**

   根拠: [plan:343](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:343) は V6 の入力を固定して比較しますが、[plan:350](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:350) は「初回壁なので公差未達、補正へ戻す」を正解としています。参照の `case/71.m55_combustion_wt/run_0014_chain_v2/` は**3 係数の補正済み壁**です。保存された [WT_QUALIFICATION.json:50](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/run_0014_chain_v2/wt_qualify_provenance/WT_QUALIFICATION.json:50) は `accepted` で、出口コア M は 5.502841913、V5c 最大偏差は 0.062305 %、片振幅は 0.057285 %。現在の公差を満たしています。

   **対案:** V6 再現の正例は、成果物まで揃えて `accepted` を期待してください。「補正へ戻す」は別の未補正壁の fixture に分離し、カード・公差・参照 snapshot を明示します。同じ参照に相反する期待結果を置かないでください。

4. **Major — 完了条件が、要求カードから成果物まで一度も成功しない実装を許します。**

   根拠: [plan:393](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:393) は V4 を「詰まる所の一覧と、その対処」で完了としています。V1b は既存の前段成果物を使う 1 pass、V3 は格子までなので、カード生成・初回 NS・較正・次の壁への引き継ぎ・最終判定を通す成功条件がありません。正しく停止することだけで、§1 の「要求カードから成果物まで」を達成扱いにできます。

   **対案:** 既知の正常カード 1 件について、決定的な recipe が入口から必須成果物と最終 `accepted` まで到達することを完了条件にしてください。両モデルの全試行成功を要求する必要はありません。モデルの介入回数や停止箇所の測定は、この機能の成功条件と分けて残します。

5. **Minor — campaign 台帳と case の run 索引の同期が実装項目にありません。**

   根拠: [plan:293](/home/sano/work/forge-cardmesh/plans/active/tooling-design-orchestrator.md:293) は台帳からの集計を定めていますが、[AGENTS.md:59](/home/sano/work/forge-cardmesh/AGENTS.md:59) が必須とする case `README.md` の「計算 run 一覧」への反映担当がありません。自動再試行や計時 run が台帳にしか残らなくなります。

   **対案:** #6 に、run 作成・終了時の索引同期を追加してください。台帳から当該 campaign の行を生成し、他の run の行を保持する形で十分です。

推奨は、**現在の小さな campaign 実行器という構成を維持し、上記 1 → 2 → 3 → 4 の順で契約と検証条件を修正してから実装へ進むこと**です。親計画との分担、ソルバ演算を変更しない範囲、node を対象とする検証方針は妥当です。

証拠は commit `9436b084` のコード・計画と既存 run を確認しました。`run_0014_chain_v2/` に対して `check_settled.py --end-step 46171 --no-write` を再実行し、**`VERDICT: SETTLED`** を確認しています。対象は窓 36,000〜46,000 と停止観測 46,171。保存済みメッシュ判定は **`VERDICT: PASS (AR<= 5000, skew<= 0.90)`** です。qualification 全体の再実行は一時ディレクトリ制限で完了せず、保存 JSON の全検査が PASS であることと、現行 `decide` がその入力から `accepted / done` を返すことを確認しました。run の索引は [case README:38](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:38) です。

証拠: 足りる

指摘数: Critical 0 / Major 4 / Minor 1
