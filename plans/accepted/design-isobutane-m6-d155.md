# M6 出口径 1.55 m イソブタン風洞 — 最短全長スタディ (case/45)

## メタ

- **area**: `design campaign`
- **status**: `done`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (axis-Mach チェーン)
- **related_plans**:
  - [`tooling-nozzle-shortest-robust-axis.md`](../accepted/tooling-nozzle-shortest-robust-axis.md) (最短基準の踏襲元)
  - [`tooling-nozzle-axismach-length-dv.md`](../accepted/tooling-nozzle-axismach-length-dv.md) (全長 dv 化)
  - [`tooling-nozzle-axismach-viscous-deltastar.md`](../accepted/tooling-nozzle-axismach-viscous-deltastar.md) (δ* v1/v3 フロー)
- **created**: `2026-08-31`
- **owner**: ユーザ指示 (2026-08-31)

## 1. 目的

M6・出口径 1.55 m の軸対称風洞ノズルについて「全長をどこまで短くできるか」を、
設計 (逆 MOC) → Euler CFD → 粘性壁 (δ\* 物理壁 + SST NS) まで通して確定する。
出口コア M は可能な限り 6 に近づける (既往 M6 到達点: NS v3 で −0.43 %)。

## 2. インタビュー決定事項 (2026-08-31)

| 項目 | 決定 | 根拠 |
| --- | --- | --- |
| ガス | イソブタン燃焼ガス φ=0.9 (CO2 0.1677 / H2O 0.0858 / O2 0.022 / N2 0.7245), semi-perfect | case/42 M6 系列と同組成・検証資産と直接比較可能 |
| Pt / Tt | 5.5 MPa / 1600 K | ユーザ指定 |
| 出口径 1.55 m の定義 | **δ\* 込み物理壁** | ユーザ選択。r_t = (0.775 − δ\*_exit)/(r_F/r_t)。δ\* は NS v3 実測で 1 回補正 |
| R (スロート曲率) | 2 と 3 で最短探索し短い側 | 最短実績 R3 / 生産推奨 R2 の両にらみ |
| 入口配管 | 物理半径 0.5 m (r_inlet = 0.5/r_t)、L_U 12 | case/42 流用。L_U8 悪化・16 同等の実測から 12 |
| 凝縮 | 後段確認 (dry 本体 + 採用点のみ診断→必要なら凝縮 ON 再評価) | ユーザ選択 (既定 a) |

## 3. 事前見積りの重要所見 (凝縮)

Tt 1600 K の M6 出口は静温 233.5 K、H₂O 分圧 303 Pa の飽和温度は 263.9 K —
**出口で −30 K の過冷却 (S≫1)**。M5 (+29 K で dry 確定) と異なり dry 前提が
出口近傍で破れるため、採用点の凝縮 ON 再評価は必須になる見込み。CFD YAML は
最初から `tp_species: split_h2o` (dry では pseudo と同一 [M4.2 検証済]、凝縮
restart 可能) で組む。飽和到達 M は 5.3–5.5 帯の見込み (dry run の Tsat_post で確定)。

## 4. 設計方針

- 軸則 = knot (A6)。最短判定は shortest-robust study と同一基準:
  hard gate (壁 QA + CFDWall validate + 壁単調 + flip/交差/堆積 0) + 単峰 +
  min(μ_w−θ_w) ≥ 1° + min sin∠ ≥ 0.02。n600 で探索し境界を n1200 で確認。
- 探索は `case/45.isobutane_m6_d155/search_shortest.py` (L_c 0.1 刻みまで。
  0.01 刻みは追わない — 0.1 r_t ≈ 8 mm は工作公差未満)。
- 粘性は δ\* v1 (相関) → v3 (CFD 抽出, `build_dstar_v3_m5.py` 流用)。M5/M6 の
  教訓どおり v1 はゲート外級 (−0.6 %) を想定し v3 を正とする。
- r_t は δ\*_exit 初期推定 0.66 r_t (case/42 run_0105 の 0.783 r_t を Re・全長で
  スケール) で仮置きし、v3 実測後に 1 回補正 (設計は r_t 単位で不変、NS のみ再実行)。

## 5. 実行ステップ

1. 最短探索 (R2/R3) → 勝者確定 ✅ **R2 / L_c 39.3 / M_K 2.7, x_F = 95.104 r_t**
   (R3 は 95.433)。margin 1.002° = ロバスト床の最短境界。
2. Euler CFD (run_0001, dry split_h2o) — dM ≤ 0.5 %・出口一様性・overshoot。
3. NS coarse 中継 → NS v1 (相関 δ\*) → δ\* v3 抽出 → NS v3。
4. r_t 補正 (物理出口半径 0.775 m 合わせ) → 必要なら NS 再実行。
5. 凝縮後段確認: dry の Tsat_post 診断 + 凝縮 ON restart 再評価 (出口 M への影響)。
6. 成果まとめ (最短全長・出口 M・凝縮影響、README run 表・結果ページ)。

## 6. 結果 (2026-08-31 完了)

- **最短全長: スロート→物理出口 7.337 m** (x_F 96.003 r_t, r_t 76.421 mm, 出口径 1.550 m)。
  入口配管端→出口 8.292 m。E→F 一様化区間 ≈55.8 r_t = 4.36 m は物理の床。
- **粘性トリム**: NS 固定点の出口欠損 −0.24 % が線形・再現的 → Md 6.0144 で設計し
  **NS dry 出口軸 M 6.0008 (+0.013 %)** (run_0007, δ* 2 巡固定点 median 0.999)。
- **凝縮**: dry 軸は x≈24 r_t (M5.5) で飽和線越え、S≈14。凝縮 ON (Kw+HK+蒸発) で
  onset x≈69 r_t・出口 g 0.20 %・**出口軸 M 5.927 (−1.2 %)** + 試験区間 M 低下勾配。
  凝縮フリーの目安 Tt ≳ 1820 K。トリムでの打ち消しは正帰還のため不採用。
- run 台帳と VERDICT は case/45 README。結果ページ:
  https://claude.ai/code/artifact/f1cfbdf8-2415-4bfc-ae2d-156793569bd0

## 変更ログ

- 2026-08-31: 起票。最短探索完了 (R2 39.3/2.7 → 95.104 r_t)。
- 2026-08-31: Euler/NS v1/v3 (トリム前後)・凝縮 ON まで完了、status done。
  ショートカット失敗の記録: 旧設計の δ* CSV 流用は物理壁フィルタ不合格 (正規 v1→v3 要)。
- 2026-09-04: **粘性トリム (Md 6.0144) は撤回**。トリムが打ち消していた −0.24 % はスロート δ\* の 8 倍過大 (相関) による有効スロート面積 +2 % が主因と判明
  ([調査ノート](../../notes/investigations/nozzle-deltastar-throat-review.md))。新チェーン (積分法初期壁 + 固定 Euler 基準抽出,
  [plan](tooling-nozzle-deltastar-core-matched-euler.md)) で **Md 6.0 のまま出口面コア M +0.01〜0.02 %** (case/45 `run_0023`)。生産形は run_0023 の壁。
- 2026-09-04: **最終形 = case/45 `run_0038_ns_final_rt77p02`**: r_t 77.02 mm (`deltastar_loop --solve-rt 0.775` で出口半径 0.775 m 合わせ)、
  物理出口半径 0.7750 m、全長 (スロート→出口) 7.325 m、ṁ_NS/ṁ_E 1.0002、出口面コア M 6.0002、Md 6.0 (トリムなし)。点列 `points_d155_final_{ns,euler}.csv`。
  上流直管を 10 r_t にしても δ_r・出口 M は不変 (`run_0034`)。凝縮 ON を最終壁で再評価 (`run_0039`): 軸 onset 57.9 r_t、出口軸 M −1.36 %、出口面コア M −0.18 %、g 0.27 % — 結論不変 (Tt 1600 K は凝縮不可避、凝縮フリーは Tt ≳ 1820 K)。
