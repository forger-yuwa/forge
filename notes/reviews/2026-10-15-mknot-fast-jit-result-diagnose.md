# codex 諮問 (diagnose): mknot-fast-jit-result

- **brief**: [`notes/reviews/briefs/2026-10-15-mknot-fast-jit-result.md`](../../notes/reviews/briefs/2026-10-15-mknot-fast-jit-result.md)
- **plan**: [`plans/active/tooling-nozzle-mknot-fast.md`](../../plans/active/tooling-nozzle-mknot-fast.md)
- **date**: 2026-10-15
- **commit**: `a36051fb` (feature/mknot-fast)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 2.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-mknot-fast-jit-result-diagnose.codex.log` (codex: ok)
- **結論**: 既定化の受入を保留し、位置差 0.25 を生じた元の前線を復元して、numpy／numba の A/B を一つ実施する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 判断対象 | 採否 | 理由・対案 |
|---|---|---|
| 現状の V5 で numba の既定化を受け入れる | **要再検証（Major）** | 位置差 0.25 が出た入力を除外しただけで、適用範囲外だと示していない。失敗入力を戻して照合する。根拠: `plans/active/tooling-nozzle-mknot-fast.md:236`、`design/tests/run_moc_jit_tests.py:423`。 |
| 「3 倍」未達を受け入れる | **採用** | 結果後の変更とユーザ決定が明記されている。追加高速化は今回の阻害条件にしない。根拠: `plans/active/tooling-nozzle-mknot-fast.md:214`。 |
| V3(b)③ を未判定のまま plan 完了とする | **却下** | 寸法ゲートは受入条件で、未判定を合格にできない。別 plan への正式な移管が決まるまでは未完了として残す。検証のためだけに節点上限を緩和しない。根拠: 同ファイル `:207`、`:210`、`:227`。 |

結論: **既定化の受入を保留し、位置差 0.25 を生じた元の前線を復元して、numpy／numba の A/B を一つ実施する。**

第 1 仮説: **傾きで交点を求める式が垂直特性線付近で丸め差を増幅し、現行版でも V5 の一致条件を破る入力が残っている。** 確度: 中  
  根拠: `moc_fill_nb.py:124`・`:125` は sin/cos、`:132` は大きな傾きと小さな座標差の積、`:138`・`:139` は再び cos で割る。`:127` の検査は傾きの差であり、**垂直に近い特性線を除外する検査ではない**。位置差 0.25 は plan `:236` の記録値で、受入上限 1e−10 r_t を大幅に超える。ただし、こちらの再構成した小入力では backend 差は 0 であり、元の反例を独立再現したわけではない。  
  反証条件: 元の失敗入力を現行版に戻しても、保存される有限点の位置・θ/ν、欠損位置、対ごとの分類・分岐、`moc_gate` がすべて V5 の条件を満たす。

第 2・第 3 仮説: 追加しない。表の補間は節点・範囲外を含む照合が試験にあり、現時点で最優先に戻す根拠はない（`run_moc_jit_tests.py:372`）。

判別 A/B: **変更する設定は `FORGE_MOC_KERNEL=numpy／numba` だけ。** 元の失敗前線と、そのときの `axis_limit`・`corrector`・反復上限を固定し、`fill_levels` を前線消滅まで各 1 回実行する。CFD は 0 step。既存の `RecDiag`／`run_both` を使い、位置差 ≤ 1e−10 r_t、θ/ν 差 ≤ 1e−10 rad、欠損位置・対の分類・源項分岐・`moc_gate` の一致を判定する。  
→ **A: 全条件を満たす**なら、この既知反例が現在も残るという仮説を棄却し、JIT の result レビューへ進む。  
→ **B: 一つでも破る**なら既定化の受入を止め、当該反例を回帰試験に保持する。入力を削除して合格に戻さない。

やらない方がよいこと: 追加高速化や `L_c` の使い回しを再開すること。現在の判断を分けるのは速度ではなく、既知の数値的不一致の処置である。

呼び出し側の前提への異議: 「照合の入力に向かない」は「実装の適用範囲外」の証明にならない。また `FORGE_MOC_KERNEL=numpy` でも網の走査と壁閉包は numba を通る（`moc_diagnostics.py:160`、`moc_inverse.py:617`）。したがって V5 単独を変更全体の独立照合とは呼べず、共通部分の旧版参照試験を別の根拠として残す必要がある。

不足情報: 元の失敗入力・実行モード・差が最大になった段と対を、再試験の記録に添付してほしい。AWS の所要は今回原票を再測定・閲覧せず、提示された計画内の記録として扱った。

証拠: 足りる
