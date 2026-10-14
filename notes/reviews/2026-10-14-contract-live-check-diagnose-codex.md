# codex 諮問 (diagnose): contract-live-check

- **brief**: [`notes/reviews/briefs/2026-10-14-contract-live-check.md`](../../notes/reviews/briefs/2026-10-14-contract-live-check.md)
- **date**: 2026-10-14
- **commit**: `d9f554f4` (feature/nozzle-cooled-wall-settling)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 明示 (`--backend codex`)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 0.5 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-contract-live-check-diagnose-codex.codex.log` (codex: ok)
- **結論**: 指定条件では `_capped_layout` は呼ばれないため、このままでよい。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: 指定条件では `_capped_layout` は呼ばれないため、このままでよい。
第 1 仮説: `axis_gap_frac` 分岐では、`int(prm.nj) − 1 > max_nj` が成立すると配列生成前に返る。確度: 高
  根拠: `design/forge_design/meshing/mesh2d.py:166` で条件判定し、`:171` で `(int(prm.nj) - 1, False)` を返す。`_capped_layout` の呼び出しは後続の `:172`。
  反証条件: この条件を満たす入力で `_capped_layout` の呼び出しが観測されること。
第 2・第 3 仮説: 無し。
判別 A/B: 併用禁止設定を避け、`axis_gap_frac` を有効、`nj=10` に固定し、`max_nj` だけ A=8 / B=9 に変えて関数を各 1 回呼ぶ。呼び出し記録用モックで `_capped_layout` の呼び出し回数を見る。→ A が 0 回、B が 1 回なら早期 return を確認／A が 1 回なら第 1 仮説を棄却。
やらない方がよいこと: 確認のために巨大な `nj` で実配列を確保すること。
呼び出し側の前提への異議: 無し。
不足情報: 無し。A/B は未実施だが、呼び出し順の判断には提示範囲のコードで足りる。
証拠: 足りる
