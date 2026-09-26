# 諮問ブリーフ: result 8 巡目の Major 全件採用と `accepted/` への移動 (2026-09-26)

AGENTS.md エスカレーション条件 **5** (codex の Critical / Major の採否)。

## 読んでよいファイル
- レビュー本文: `notes/reviews/2026-09-26-boundary-conjugate-heat-transfer-result.md`
- 反映: plan `plans/active/boundary-conjugate-heat-transfer.md` の §5.1 **#104 の行だけ** (`grep -n "^| 104 |"`)、§4.11(a) (`grep -n "### 4.11"` から 20 行)、§6「閉じる条件」
- 差分: `git show a73ad11f --stat` と、`git show a73ad11f -- case/58.conjugate_slot/ solver_density_cuda/tools/check_cht_interface.py procedures/su2-cross-check.md methods/boundary.md`
- 負例試験: `python3 case/58.conjugate_slot/test_v6p_negative.py` (実行してよい。一時ディレクトリで動き、元 run は h5 だけ読む)

## 採否 (全件採用・却下 0。すべて修正済み)
- M1 (後壁符号): 符号つき比較に変更。負例 (後壁 q_eff 反転) → FAIL を確認。
- M2 (節点ログの整合): `load_node_log()` で全検査。負例 3 件 (座標 NaN+残差 1000 / 1 節点 step=0 / ID 重複) → REFUSED。正例 2 件 PASS。
- M3 (誤差保証): §4.11(a) と #85 の「SU2 0.0001 %」を forge–SU2 −0.02471 K / +0.00409 K に訂正、「1 K 未満・case 横断」を撤回。
- m4 / m5: SU2 手順・methods・case/48 README を同期。`plans/README.md` は移動時に同期する。
- 修正後に 4 run を再評価: 数値不変・全 PASS。

## 諮りたいこと (推奨を 1 つに)
1. 上の採否 (全件採用) で妥当か。修正が不十分な項目があれば 1 つ挙げる。
2. **`status: done` → `plans/accepted/` へ移してよいか**。result 9 巡目が要るか。
