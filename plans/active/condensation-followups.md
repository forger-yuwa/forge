# 非平衡凝縮モデルの後続課題 (carrier 補正・空気凝縮の引き継ぎ)

## メタ

- **area**: `condensation`
- **status**: `draft`
- **related_docs**:
  - `methods/condensation.md` (現在仕様: 核生成・成長・Kantrowitz/Feder 補正・空気 CPG carrier 形・N2 低温物性)
- **related_plans**: 親 [condensation-kantrowitz-carrier.md](condensation-kantrowitz-carrier.md), [condensation-air.md](condensation-air.md), [condensation-kantrowitz-gamma-twophase-sonic.md](condensation-kantrowitz-gamma-twophase-sonic.md)
- **created**: `2026-09-13`
- **owner**: `sano`

## 1. 目的

2026-09-12/13 の 3 つの凝縮 plan (γ_v・二相音速 / Feder carrier 形 / 空気凝縮) で「後続」に送った課題を 1 か所に集め、
親 plan を accepted へ移した後も残作業の正本が消えないようにする (codex result レビュー 2026-09-13 の要求: 後続を active plan にリンク)。
親 plan の §5.1 で「後続」と書かれた行は本 plan の表を正本とする。

## 2. スコープ

- **やる**: 下表の課題の設計・検証 (各項目は着手時に本 plan を `in_progress` にし、§4 に設計方針を書いてから codex plan レビュー)。
- **やらない**: 親 plan で検証済みの実装の再検証。

## 3. 関連 docs と前提

- Feder carrier 形と σ 感度: [condensation-kantrowitz-carrier.md](condensation-kantrowitz-carrier.md) §9 (Wysłouzil node run_0350–0356)。
- 空気 CPG carrier 形と N2 低温物性: [condensation-air.md](condensation-air.md) §9 (case/34 run_0014–0035)。
- 実験差の現状: Wysłouzil は mode 1 が壁圧 −4.5 % (21 mm)、mode 2/3 は +15〜+18 % (計算上の onset が 8 mm 上流 = モデル間差)。
  Arthur は壁 cond/dry が 3–4 in で実験より ~9 % 過大 (物性修正で解消せず、原因未切り分け)。

## 4. 設計方針

各課題の着手時に記述する (現時点は未着手)。

## 5. 実装ステップ

着手時に記述する。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 (出典) |
| --- | --- | --- |
| 1 | 成長率 dr/dt の α 感度 | Hertz–Knudsen の質量適応係数 α=1 (上限) を `condAccommodation` キーにし、Wysłouzil mode 3 で α=0.5/0.1 の onset・壁圧応答を見る。J と dr/dt の寄与を分離する (carrier plan #5b, ユーザ指摘 2026-09-12) |
| 2 | 分圧スイープ | Wysłouzil の p_v0 0.5 / 0.26 kPa の dry node 場を作り mode 1/2/3 の分圧応答 (実験 Fig.3 との比較) (carrier #5) |
| 3 | `condKantrowitz` 既定値の決定 | 1 と 2/3 のどちらを既定にするかを #1–#2 と実験一致で判断 (carrier #6) |
| 4 | σ 0.97 系の h0 変動 | case/16 run_0354 の中心線 h0 偏差 (平均 1.200 ± 0.015 kJ/kg, 4 点で上昇傾向) が限界サイクルか drift か: 保存間隔を密にして判定し `check_quasisteady.py` に h0err を統合 (carrier #7) |
| 5 | Arthur 3–4 in の壁圧 ~9 % 過大の原因切り分け | 候補: 核生成 J (CNT×Iland), 成長 dr/dt (Goodheart, α), 壁圧の抽出位置 (壁セル列 vs 静圧孔), Arthur 記号の読み取り。J/dr/dt 各 ×0.5/×2 の感度と Fig.2 再デジタイズ (air #10) |
| 6 | CPG 二相音速 | γ_2φ + 実 Ht + 一般 EOS 固有系を CPG 二相 (N2/空気) へ; 流束 FD 照合 (double/float32) (air #7, sonic plan) |
| 7 | 混合液 (露点線) モデル | O2/N2 理想溶液の露点線、2 成分凝縮 (air #8) |
| 8 | 条件の拡張 | Longshot 級 (M 10–14, Ṗ 小) と Daum & Gyarmathy の複数条件で理論線・実験点との比較 (air #9) |
| 9 | 境界試験の実装経由化・流束収支 | slip ghost/bvar を kernel 経由で通す試験、node/cell の質量・エネルギー流束収支 (air #6b) |
| 10 | `check_quasisteady.py` 統合 | onset/壁圧比 (dry 参照場を要する case 固有量) を `--quantity` に (air #6c) |

## 6. 検証

各課題の着手時に記述する。判定の共通ルール: 同一バイナリ反復 ≥3 run の全ペア最大をノイズ床にし `diff_res.py --tolfile --factor 2` の exit code で判定、
凝縮 run は `check_convergence.py` (凝縮列込み) + `compare_condfix.py --series` / `onset_analysis.py --series`、実カーネル照合は `verify_theta.py`。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan 免除 | `2026-09-13` | — | — | 起票のみ (draft; 設計方針未記述)。各課題の着手時に §4 を書いて `--stage plan` を回す |

## 7. 影響範囲

着手時に記述する。

## 8. 完了条件

表の各行が個別 plan に分離されるか、実装・検証済みで decision が書かれること。

## 9. 変更ログ

- `2026-09-13` — 起票 (codex result レビュー ② の要求: carrier / air の後続を active plan にリンク)。

## 10. 未確定事項

- なし (課題ごとに着手時に書く)。
