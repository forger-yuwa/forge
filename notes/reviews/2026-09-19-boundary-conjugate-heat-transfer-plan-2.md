# codex レビュー: boundary-conjugate-heat-transfer (plan, 2 巡目) — **中断**

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `plan` (2 巡目。1 巡目 [NO-GO の記録](2026-09-19-boundary-conjugate-heat-transfer-plan.md) を全件採用した改訂稿に対するもの)
- **date**: 2026-09-19
- **結果**: **中断 (判定なし)**。codex が 196 s 時点で
  `ERROR: You've hit your usage limit ... try again at 8:40 PM` を返し、最終メッセージが得られなかった (rc=1)。
- **生ログ**: `notes/reviews/2026-09-19-boundary-conjugate-heat-transfer-plan-2.log` (**git 追跡外**, 464 KB)

## 中断前に codex が実行した検算 (ログから回収)

codex は指摘本文を書く前に Python で数値検算を行っており、その出力がログに残っていた。
**3 件とも plan の記述の誤りを示していたので、判定を待たずに採用して修正した**。
いずれも手元で独立に再計算して一致を確認済み。

| # | codex のログ出力 | 意味 | 反映先 |
| --- | --- | --- | --- |
| 1 | `sign test Fout_wall=80 C_into_fluid=-20: plan -60.0 conservative 100.0` | 界面熱量の**符号と組合せ**が誤り。外向き正の面流束 $\sum F$ と流体への供給が正の拘束反力 $C$ に対し、固体へ入る熱は $\sum F - C$ (=100)。plan の $-(\sum F + C)$ は $-60$ | §4.3 を書き換え。検算例を本文に明記し、V1 で符号と絶対値の再現を解除条件にした |
| 2 | `distributed SPD example: true diag [2. 2.] componentwise secants [1. 1.]` / `spectral radius with secant D 1.818...` | **成分ごとのセカント $D_f$ は真の応答 (非対角結合を持つ) を過小評価**し、反復が発散する (スペクトル半径 1.82) | §4.2: 既定を上界 $k_{\rm eff}A/d_1$ に変更。セカントは $D_f$ を**増やす方向にのみ**使用 |
| 3 | `radiation kW/m2 [3.54, 56.70, 287.06]` | 放射の桁の誤記 (plan は 1500 K を 57 kW/m² と書いていたが 287) | §4.10 |
| — | `V1 analytic 800.0 100000.0 100000.0` | 改訂した V1 の解析解 ($T_w$=800 K, 両側 $q$=100 kW/m²) は**正しい**ことの確認 | 修正不要 |

## 次にやること

- **上限解除後 (20:40 以降) に 2 巡目をやり直す**。上の 3 件は反映済みなので、同じ focus
  (1 巡目指摘の解消確認 + §4.9 一次適用先 + §4.4c `fem2d`) で再実行する。
- 判定が出たら本ファイルを**上書きせず**、正規の 2 巡目記録として別途残す (本ファイルは中断の記録)。
