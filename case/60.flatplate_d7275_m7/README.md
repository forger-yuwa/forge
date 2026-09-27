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
| `run_0004_t26_wdA` / `run_0005_t26_wdB` | **T4-0b-WD** (事前登録 commit fa55f89a): 下流 slip を physID 7 に分けたメッシュで、A = 壁距離そのまま / B = 変換時 `wallDistExtraPhysIDs: [7]` (x ≥ 2.6 m の壁距離だけ変わる)。run_0002 の最終場から restart_field --keep-src-dtype (化学種を含む 12 量ビット一致)、同一設定 20k、1000 毎。**1 回目は新メッシュに roY が無く 7 量しか移らなかったので起動 1 分で止めて作り直した** | **1 回目は AWS ディスク満杯で step 1900 に Write failed → 整理して 2 回目を投入** (出力 2000 毎)。2 回目の結果: 両腕 `NOT CONVERGED (stalled/plateau)`。B の末尾残差は A の 0.53–0.73 倍で **1/10 に届かない → 事前登録の読み = 壁距離単独主因説を退ける** (codex 再読: 『壁距離変更だけでは解消しないが寄与はある』)。熱流束の差は比較点で 1.9e-7・パネル域で 2.1e-6 → **比較量は下流の変更に鈍い**。St 両腕 ALL STEADY、R 平均 1.313 | ref (A/B) |

> **AWS ディスク整理 (2026-09-27)**: run_0001 の中間スナップショット、run_0002 の res_0 (run_0001 最終場の複製)、run_0003_t26_probe の場のファイル (結論は上の行に記録済み) を削除。

> **AWS ディスク圧縮 (2026-09-27、ユーザ指示「貴方の関連のディスク容量圧縮の為、ファイルをガンガン削除していい」)**: 判定の根拠 (壁・表面出力 `res_wall_*`/`res_plate_*`/`res_gap_*`、残差履歴、VERDICT、時系列 CSV、設定) と最新の継続用最終場 (Cary: run_0005・0014–0024、D-7275: run_0002・run_0005) だけを残し、それ以外の全場スナップショット・run 内の `mesh.h5`・旧 3D メッシュ `gap3d_prod.h5` を削除した (計約 13 GB)。削除した場が要る再解析は再計算になる。
| `run_0006_t26_resmap` | run_0002 から 2000 step、`FORGE_OUT_RESIDUALS=1` で残差場を狙ったが、`output.level: 1` のフィルタで残差場が落ちて取れなかった | — | 破棄 |
| `run_0007_t26_resmap2` | 同上、`extraFields` に `res_ro/res_roUx/res_roUy/res_roe` を明示 | **残差場の二乗和の 99.99–100 % が平板の後ろ (x ≥ 2.6 m、節点の 6.3 %) にある**。最大は出口∩上端 slip の角 (x≈2.79 m、y≈0.47–0.50 m) と出口∩下流 slip の角 (y=0)。比較域 1.07–2.46 m と境界層 (wd < 3 cm) は 0.00 % → **残差プラトーは出口の角に局在** (状態差の測定 run_0003 と整合、今度は残差そのもの) | ref (診断) |
| `run_0008_t26_outflow` | **T4-0b-OUT** (事前登録 commit 57b94f2e): run_0002 の最終場から出口だけ `outlet_statPress` → `outflow` (現行レシピの超音速 node 出口)。20k、2000 毎、残差場も出力 | `NOT CONVERGED (stalled/plateau)`。末尾残差は同じ長さの対照 run_0004 (statPress) と 1.00 倍 → **出口の種類を替えても改善しなかった** (超音速流出では両者とも内部状態の外挿で実質同じ作用なので、出口閉包そのものは除外していない — codex)。残差は依然 100 % が x ≥ 2.6 m。熱流束の差は比較点 3.9e-6・パネル域 4.9e-6 → **比較量は鈍い**。St ALL STEADY、R 平均 1.313 | ref (A/B) |
| `run_0009_t26_relaxA` / `run_0010_t26_relaxB` | **T4-0b-RELAX** (事前登録 commit cd395f01): run_0002 の最終場から implicitRelax 1.0 / 0.5 だけを変える。20k、残差場 (k/ω 含む) を出力 | 両方 `NOT CONVERGED (stalled/plateau)`。末尾残差 B/A = 0.86–1.19 (1/10 に遠い)、出口近傍 (x ≥ 2.7 m) の 2000 step 状態差は B が A の約 0.5 倍 → **事前登録の読み = 判別保留** (中間的。振れ幅の半減は更新量半減の見かけの可能性)。A の roK は落ち着き条件も未達 | ref (A/B) |
| `run_0011_t26_ext` | **T4-0b-EXT** (ユーザ指示「出口側を slip にして延長したらいいよ」、事前登録 commit 820d20a5): 下流 slip を 0.2 → 1.0 m (出口 x 3.6 m) に延長したメッシュ `fp_d7275_y3_ext` (上流・平板の分割は元と一致、42.5 万節点、`SOFT-PASS` AR 最大 ≈ 1043)。run_0002 の最終場から interp_field (化学種を含む 12 量、ΣY = 1)、全域 FP64、本段設定 40k、5000 毎、残差場も出力 | `NOT CONVERGED (stalled/plateau)`、未達残差は同じ長さの対照の **1.9–14.6 倍に悪化** → 事前登録の読み = 延長だけでは解消しない。**残差の 96–98 % は延長した slip 区間の内部 (2.6–3.5 m)、出口付近は 2–4 %** → 源は出口の角でなく平板の後ろの slip 後流。St ALL STEADY、R 平均 1.313 (延長前と同じ) | ref |
