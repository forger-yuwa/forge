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
| `run_0001_t26` | 試験 26 (α 0°)。燃焼ガス、壁 300 K、前縁から乱流 (SST、Tu 0.5 %、μt/μ 10)、`fp_d7275_y3`、段階起動 lam→soft→mid→2 次ランプ→本段 60k (cfl 1.5)、**全域 FP64** (`~/forge56-double`)、段の引き継ぎは restart_field | 本段 60k 完走、NaN 0。判定区間 (`--segment`) で `NOT CONVERGED (still converging)` (全列 1.6–3.3 dec、下降中)。壁解像: y1+>1 は前縁の 2 列のみ (最大 3.54、長さ重み 0.03 %)、パネル域 1.07–2.46 m は 0.26–0.28。暫定 R (St*_l、前縁起点) = 1.27–1.43 (平均 1.31)、トリップ起点で平均 1.32、x 対応 ±3 % で 1.307–1.314 (**未収束につき解釈しない**) | 延長元 |
| `run_0002_t26_ext` | run_0001 から同一設定 +60k (restart_field、倍精度のまま 12 量ビット一致)。1 回目は `species_db.yaml` を複製し忘れて forge が起動せず → スクリプトを直して空の run を消して再投入 | 延長区間 `NOT CONVERGED (stalled/plateau)`、run_0001 の判定区間と接続しても 2.2–3.9 dec で `still converging`/プラトー (**全域 FP64 なので丸め床ではない**)。St **ALL STEADY**。**R (St*_l、前縁起点) = 1.266–1.442、平均 1.313** (トリップ起点 1.328、x ±3 % で 1.313–1.319)。forge の q_w: パネル平均 89.0 kW/m²、位置 II 87.8 kW/m² (Tw 300 K) | active (**T4-0b 試験 26 の結果、収束ゲート未達**) |
| `run_0003_t26_probe` | run_0002 の最終場から同一設定 4000 step、1000 毎に出力 (残差の出どころを探す短い継続) | 2000→4000 step の場の相対変化は**下流の slip バッファ (x > 2.55 m、平板の後ろの後流) だけ** (ρ 7e-3、P 8e-3、roK 0.15)。比較域の壁近傍 (wd < 1 mm) は 2e-6 以下、前縁 1e-6 → **残差プラトーは比較域の外に局在** | ref (診断) |
| `run_0004_t26_wdA` / `run_0005_t26_wdB` | **T4-0b-WD** (事前登録 commit fa55f89a): 下流 slip を physID 7 に分けたメッシュで、A = 壁距離そのまま / B = 変換時 `wallDistExtraPhysIDs: [7]` (x ≥ 2.6 m の壁距離だけ変わる)。run_0002 の最終場から restart_field --keep-src-dtype (化学種を含む 12 量ビット一致)、同一設定 20k、1000 毎。**1 回目は新メッシュに roY が無く 7 量しか移らなかったので起動 1 分で止めて作り直した** | 実行中 | active |
