# case/59 — Cary TN D-5863 冷却平板 (M6, 空気) = **T4-0a (判別)**

計画: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md) §4.12 / §6 G14 / §5.1 #61
精読メモ: [`notes/investigations/tp1187-test-system-reading.md`](../../notes/investigations/tp1187-test-system-reading.md)

## 目的

TP-1187 (case/56) の乱流分母 q_FP の出所 TN D-7275 は、実測が Eckert 参照温度法より 10–30 % 低い。
forge は理論側に乗る。**このずれが冷壁に特有か (モデル側) / 水準だけか** を、壁温比だけを振った実験で判別する。
合格条件は「実測に一致」ではなく、完全乱流域の R = St_forge / St_exp の**壁温比依存**を測ること (plan §6 G14)。

## 実験 (Cary 1970, NASA TN D-5863)

- Langley 20-inch hypersonic tunnel、**空気**、M∞ 6.02 ± 0.02、Tt 533 K、T∞ 65 K
- 鋭い前縁 (厚み 0.0038 cm) の平板、AISI 405、薄板 0.076 cm、液体窒素で内部冷却、両側に端板、α = 0° (Me 6.0)
- Pt 3.55 MPa (Re 2.6–3.1 × 10⁵/cm) と 1.83 MPa (Re 1.4 × 10⁵/cm)
- **Tw/Tt 0.2–0.7 の 9 系列** (Table II、転記済み `conditions.json`)。自然遷移で、加熱ピーク (遷移終端) は x = 14–35 cm
- St_∞ = q / (ρ∞ u∞ cp (Taw − Tw))、乱流の回復係数 0.89

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| (まだ無い) | | | |
