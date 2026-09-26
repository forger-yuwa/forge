# case/60 — Deveikis & Hunt TN D-7275 大型校正パネル (M7, 8-ft HTST) = **T4-0b**

計画: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md) §4.12 / §6 G15 / §5.1 #61
精読メモ: [`notes/investigations/tp1187-test-system-reading.md`](../../notes/investigations/tp1187-test-system-reading.md)

TP-1187 (case/56) の乱流分母 q_FP の出所。フェンス付きパネルホルダー・鋭い前縁・0.24 cm 球トリップ (前縁から 12.7 cm)。
T4-0a (case/59、Cary) の比較方法が固まってから計算する。

## 実験データ

- `digitize_d7275_fig20_test26.json`: 試験 26 (Tt 1867 K、Pt 18.06 MPa、α 0°、M∞ 6.64、Re∞ 4.757×10⁶/m) の中心線 St*_l (Fig 20 の □、10 点)、
  Table II/III の試験 26 行、換算式と読み取り幅 (St ±1.3 %、x 対応 ±3 %)。

## 準備 (2026-09-27)

- 自由流の再構成: 燃焼ガスモデル (case/50 の CombustionProducts、φ 0.711) で M∞ 6.64・p∞ 2117 Pa・Tt 1867 K から
  T∞ 227.7 K (表 228.3)、Re∞ 4.71×10⁶/m (表 4.757、−1 %)、(ρVcp)∞ 67.95 kW/m²K (表 68.89、−1.4 %)。
- メッシュ `mesh/fp_d7275_y3` (case/56 の `gen_mesh.py`、前縁から 2.6 m、y1 3 µm、H 0.5 m、37.6 万節点): `VERDICT: PASS` (AR 最大 854)。
- 生成器 `gen_runs.py` (段の引き継ぎは restart_field.py、FP64 アキュムレータ、stage manifest)。
- トリップ (12.7 cm) は 2D で入れられないので、乱流の起点は前縁と、同じ解を後処理で 12.7 cm ずらした場合の 2 通りで示す (有効長の感度、上下界ではない)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| (まだ無い) | | | |
