# 共役平板の検証 (C1/C2) — 前提ゲート FAIL の原因切り分けと再判定

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) — 共役壁 (`fem2d`)、slip 境界 (node)
  - [`notes/investigations/node-slip-tangential-density-spurious-flow.md`](../../notes/investigations/node-slip-tangential-density-spurious-flow.md) — 既知の未修正欠陥 (node slip + 接線方向の密度勾配で偽の流れ)
- **related_plans**:
  - [`boundary-cht-conjugate-benchmarks.md`](boundary-cht-conjugate-benchmarks.md) — **発注元**。C をここから分割 (2026-09-30)。登録条件 (§4.3・§4.5・§4.6・§4.7・§6) はそちらの文をそのまま引き継ぐ
- **created**: `2026-09-30`
- **owner**: `sano`

## 1. 目的

発注元 plan の C (共役平板: C1 厚さ方向の抵抗、C2 板の軸方向伝導) は、主判定 (界面温度・熱流束の独立参照解との比較) が 6 本すべて PASS したが、
**前提ゲートの流体収束と G-if が 6 本とも NOT CONVERGED** で、登録どおり**判定不能**として保留した。
本 plan はその原因を切り分け、登録条件を満たす形で C を再判定する。**完了条件は発注元 plan の §6 (C) のまま**変えない。

## 2. スコープ

- **やる**: 前提ゲート FAIL の原因切り分け (境界実装・連成反復・抽出のどれか)、対処、再判定
- **やらない**: C の窓外を事後にゲートから外して完了扱いにすること、slip を真因と決めて修正に入ること (切り分けの前に)

## 3. 引き継ぐ事実 (発注元 plan §5.1 #6・#6b)

- run: `case/65.conjugate_flat_plate/run_0005`〜`0010` (C1/C2 × n16/32/64、AWS FP64、600000 step)。
- 主判定 PASS (C1 n64: θ 差 6.1e-4、q 0.10 %。C2 n64: θ 1.3e-4、q 0.38 %)。準定常 (評価窓内の全節点) PASS、G-cons PASS。
- **前提ゲート FAIL**: 流体収束 NOT CONVERGED (stalled、rms_roe 1.4 桁で横ばい)。G-if NOT CONVERGED: 界面残差は板の端に集中 —
  n32 は後縁 x/L 1.0 (5.8 W/m²)、**n64 は前縁の直後 x/L 0.004〜0.015 (C1 678 W/m²・C2 607 W/m²)**、評価窓内は 0.14〜0.15 W/m² (許容 1.0 内)。
  細格子ほど悪化。後縁のすぐ下流の slip 境界に最大の時間変動 (P 3 Pa、T 3e-3 K)、前縁上流の slip 節点に法線速度 0.72 m/s。
- 原因は**未確定** (node slip の既知欠陥、前縁の熱流束の特異性と連成反復、抽出のどれか)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 原因の切り分けの設計 | 判別 A/B を諮問で決める (例: 連成を止めた等温壁の同条件 run で流体収束が回復するか、slip を別の境界にした場合の比較) | F |
| 2 | 対処と再判定 | 1 の結論に従う。完了条件は発注元 §6 (C) | O |

## 6. 検証

発注元 plan §6 (C) の登録文をそのまま使う。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 9. 変更ログ

- `2026-09-30` — 発注元 plan から C を分割して起票 (codex diagnose の推奨)。
