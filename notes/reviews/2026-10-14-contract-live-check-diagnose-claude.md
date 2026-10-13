# Claude 代替諮問 (diagnose): contract-live-check

> **Claude の別インスタンスによる代替** (独立性は文脈の分離だけ。依頼側と同じ系統のモデルでありうる; codex・Fable と別系統の判定ではない)

- **brief**: [`notes/reviews/briefs/2026-10-14-contract-live-check.md`](../../notes/reviews/briefs/2026-10-14-contract-live-check.md)
- **date**: 2026-10-14
- **commit**: `d9f554f4` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `claude` — 選んだ理由: 明示 (`--backend claude`)
- **要求したモデル**: `haiku`
- **実効のモデル**: `claude-haiku-5-5` (modelUsage: `claude-haiku-5-5` out 3613)
- **CLI の版**: `2.1.296 (Claude Code)`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `claude-haiku-5-5`) (同じ Claude 系統)
- **permission_denials**: 0 件
- **道具**: `Read,Grep,Glob` (`--safe-mode`、Bash なし) / 管理ポリシー: 無し (`/etc/claude-code/managed-settings.json`, `/etc/claude-code/managed-settings.d`, `/etc/claude-code/managed-mcp.json` が無い)
- **実行**: effort `high`, 0.3 min, rc=0, total_cost_usd 0.0037
- **試行のログ**: `notes/reviews/2026-10-14-contract-live-check-diagnose-claude.claude.log` (claude: ok)
- **証拠一式**: `notes/reviews/2026-10-14-contract-live-check-diagnose-claude.evidence/` (git 追跡外)
- **結論**: `radial_node_count` で `max_nj` を渡し、指定 nj − 1 が `max_nj` を超えるときは `_capped_layout` は呼ばれない。このままでよい。

本文は Claude の別インスタンスの最終メッセージ (JSON の `result`) をそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: `radial_node_count` で `max_nj` を渡し、指定 nj − 1 が `max_nj` を超えるときは `_capped_layout` は呼ばれない。このままでよい。

第 1 仮説: 呼び出し側の想定 (`_capped_layout` が呼ばれる) は誤り。打ち切りが先に効き、`_capped_layout` は呼ばれず `(nj − 1, False)` が返る。 確度: 高
  根拠: `design/forge_design/meshing/mesh2d.py:166` で `if max_nj is not None and int(prm.nj) - 1 > max_nj:` を判定し、真なら `:171` で `return int(prm.nj) - 1, False` と返す。`_capped_layout` の呼び出しは `:172` にあり、これは `:166` の if ブロックの後にあるため、条件が真のときは到達しない。`radial_node_count` 内の `_capped_layout` 呼び出しは `:172` だけ (grep で確認)。
  反証条件: `max_nj` を渡し nj − 1 > max_nj の状態で `_capped_layout` が実行されたら誤り。例えば、`prm.nj` が int 変換で値が変わる、`:166` と `:172` の間に別の分岐が入る、といった変更が入った場合。

第 2 仮説: `:194` の `generate_axisym_mesh` 側の `radial_node_count(prm)` (max_nj なし) では `_capped_layout` が cap 無しで呼ばれる。これは今回の問い (max_nj あり) とは別経路。 確度: 高 (呼び出し箇所のみ確認、中身は未読)
  根拠: `mesh2d.py:194` (grep)。`:172` が `max_nj is None` のとき必ず通る。
  反証条件: `:194` 経路で巨大 nj の配列確保が実際に起きること。これは別件 (M1 の範囲外) として扱う。

第 3 仮説: 無し。

判別 A/B: 設定 1 点 = `prm.nj=50`、`axis_gap_frac` を設定、`wall_first_frac_throat=None`、`max_nj=10` で `radial_node_count` を呼ぶ。`_capped_layout` をモンキーパッチで例外を出すものに差し替えてから呼ぶ。→ A: 例外が出ず `(49, False)` が返れば仮説 1 が成立 (呼ばれていない)。B: 例外が出れば仮説 1 は誤り、`:166` の分岐を読み直す。 (forge の起動は不要。純関数の呼び出しだけなので短く回る。ただし本セッションではシェルが無く、この A/B は私は回していない。)

やらない方がよいこと: `_capped_layout` の中身 (`mesh2d.py:92` 付近) を読んで「cap の計算が誤っている」と判断し、この分岐の修正を考えること。今回の問いには不要で、`:166–171` の打ち切りが先に効く。また `:172` を `max_nj` 付きでも無条件に呼ぶように変えること。配列確保を防ぐ意図の打ち切りが効かなくなる。

呼び出し側の前提への異議: 問いの「診断」としての枠組みは不要。これは静的に読めば答えが出る契約の確認 (判定の行 + 証拠欄の形式の確認) であり、真因の仮説を立てる段ではない。ブリーフの観測事実 (`radial_node_count` が `axis_gap_frac` 時に実効 nj を数える) は `:154–176` と一致する。`:170` のコメントは「nj が大きく第 1 層が薄いと Newton が退化して第 1 層 0 になる既存の欠陥がある」と書いてあるが、これは本問いの範囲外であり、ここでは採らない (未確認)。

不足情報: `_capped_layout` の定義 (`mesh2d.py:92–`) と `:114` の呼び出し元は範囲外のため未読。grep で位置だけ確認した。本問いの結論には不要。

証拠: 足りる (`mesh2d.py:154–176` の全行と `:172` の呼び出し箇所の grep 結果で、問いに答えられる)
