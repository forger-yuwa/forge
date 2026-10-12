# codex 諮問 (diagnose): bell-m-corrected-wall-prereg

- **brief**: [`notes/reviews/briefs/2026-10-12-bell-m-corrected-wall-prereg.md`](../../notes/reviews/briefs/2026-10-12-bell-m-corrected-wall-prereg.md)
- **plan**: [`plans/active/tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md)
- **date**: 2026-10-12
- **commit**: `6a749467` (feature/thruster-nozzle-design)
- **codex**: effort `high`, 1.8 min, rc=0
- **結論**: **`V-ext/V-zero` の差を悪化許容幅から外し、`judge_corrected_wall.py` の145–146行だけを修正してから判定台本を凍結することを推奨する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **`V-ext/V-zero` の差を悪化許容幅から外し、`judge_corrected_wall.py` の145–146行だけを修正してから判定台本を凍結することを推奨する。**

第 1 仮説: **[Major] 推定区間への感度を許容幅に加えるため、その感度による悪化を隠して「採り方に依らない」と判定できてしまう。** 確度: 高

根拠: [145行](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/judge_corrected_wall.py:145)で `u ≥ |m(Vext) − m(Vzero)|` とし、[151行](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/judge_corrected_wall.py:151)で `m(V) ≤ m(U) + u` を判定している。したがって、ある指標で小さい側が `m(U)` 以下なら、大きい側も**差の大きさにかかわらず必ず合格する**。これでは両方の合格を推定区間への頑健性の根拠にできない。

人工値による再現では、`axM(U)=0.131`、`Vext=0.130`、`Vzero=0.200`、格子差 `0.001` とすると、現行式は `u=0.070` となり両方合格する。壁圧を両方 `1.4 %`、他の指標を改善する値にすると、総合判定も「採り方に依らない」になる。**これは実測結果ではなく、判定式の反例である。**

対案: 145–146行の `u` 定義を次に置き換える。

```python
    u = {k: max(abs(m["Vext_fine"][k] - m["Vext"][k]),
                abs(m["Vzero_fine"][k] - m["Vzero"][k]))
         for k in keys + ("pw_t",)}
```

`Vext/Vzero` の差は推定区間への感度として報告し、合否の許容幅には使わない。この変更では、指標・比較範囲・30 %の改善線・`u_eta` は変更しない。

反証条件: 上記の人工値を現行式に与えて、大きい側の `axM` が不合格になること。読み取った式の再現では両方合格し、反証されなかった。実場での発生有無は未確認。

第 2・第 3 仮説: 挙げない。今回の最重指摘は上記1件。

判別 A/B: **同じ人工値で、許容幅に系列差を含める現行式／除外する修正式を比較する。** CFD は0ステップ、判定計算だけでよい。→ A: 現行は両方合格、修正後は片方だけ合格なら、この欠陥を確認／B: この組合せにならなければ、レビュー対象と実行式の対応を再確認する。再現計算の結果は **A** だった。

やらない方がよいこと: この欠陥への対処として `axM` を除外したり、比較区間を短縮したりすること。同じ問題は `exM`・`edM` にもある。

呼び出し側の前提への異議: `Vext/Vzero` の差を、両者の頑健性を判定する許容幅に兼用する前提は却下する。観測したい感度そのものが合否を緩めている。

不足情報: 修正後の格子差も、`U` 側を含む比較全体の不確かさを保証するものではない。本提案は循環した許容幅を除く最小修正である。補正壁の結果ファイルは未開封、ファイル変更なし。**plan 未反映**：呼び出し側で `plans/active/tooling-nozzle-bell-two-mode.md` §5.1 #7（7b-4）へ採否を記録すること。
