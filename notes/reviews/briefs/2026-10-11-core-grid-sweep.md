# 諮問ブリーフ: B0 の遅いモードの仕組み — 内部反復の解き残しを直接測る 1 手

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.33〜§4.41)
- **date**: 2026-10-11
- **諮問の条件**: AGENTS.md「エスカレーション条件」4 (登録に無い run の前)
- **ユーザの依頼**: 「なぜ流れ方向に細かいとライン陰解法の反復が遅くなるのか、という仕組みです。つきとめてくれ」

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.f64`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
grep -n '^### 4.3[3-9]\|^### 4.4[01]' plans/active/tooling-nozzle-core-grid.md   (その節を sed で 30 行ずつ)
sed -n '1511,1612p' solver_density_cuda/cuda_forge/timeIntegration_d.cu   (FORGE_LINE_DUMP の書き出し)
sed -n '1900,2000p' solver_density_cuda/cuda_forge/timeIntegration_d.cu   (ライン Thomas、relax を sweep の中で掛ける)
sed -n '3270,3280p;3385,3395p' solver_density_cuda/cuda_forge/timeIntegration_d.cu   (書き出しの呼び出し位置)
sed -n '2052,2100p' solver_density_cuda/main.cpp                          (DPLUR の sweep の輪、毎 step dq=0 から)
sed -n '140,200p' solver_density_cuda/cuda_forge/setDT_d.cu                (方向別 dt と上限)
```

## 1. 観測事実

- B0 = `lineImplicit 1`・`lineDtDirectional 1`・`implicitThermalJacobian 5`・`lineDtDirectionalCap 50`・cfl 4・`implicitRelax 0.7`・`nStepInner 5` (FP64 バイナリ sha 7e95e5bb)。Gref (4719 × 160) で残差プラトー、point の止まった値から θ_r +0.6 %・Q_w +0.5 % 離れて止まる。
- 右辺は不変 (上限・方向別 dt は setDT だけ)。交互更新なし、増分の消失・床の切り詰めなし (§4.35)。
- 遅いモード: 縮流部 x −7〜−2 r_t の冷却壁の境界層の縁 (壁から 60〜100 節点目)、エントロピーが主役 (s の RMS が圧力の 18 倍)、M 0.02〜0.08 (§4.38)。この場所の AR の中央値 2.7、格子に特異なし。
- **効いた**: 上限 50 → 1 (この場所の dt はほぼ不変 0.977、壁際の dt だけ 1/50): 縮流部の残差 1/25〜1/160 (§4.33)。縮流部だけ流れ方向を粗く (Grefc、列 589 → 337): 縮流部の残差 ρ 1/36・ρE 1/62 (§4.37)。
- **効かなかった** (いずれも同一起点 `run_0608` res_60000、対照 `run_0619` の同区間、窓 10k〜15k・15k〜20k):
  - F1 SST の凍結 (`run_0626`): 縮流部の比 ρ 1.00・1.12、ρE 0.94・1.00 (§4.39 結果 B)。
  - M-a ITJ 5 → 0 (`run_0627`): 縮流部の比 ρ 0.97・0.94、ρE 0.84・0.83、集中の場所の比 0.96〜1.01 (§4.41 結果 B)。全域の SST の残差は 1.6〜2.8 倍に増えた。
- **実装 (読んだ)**: `implicitRelax` は各 sweep の中で掛かる (`timeIntegration_d.cu` 1980 付近 `dq_new = relax·dq`)。ライン外の近傍は前の sweep の `dq_old` (relax 込み) を読む。毎 step `dq = 0` から `nStepInner` 回 (main.cpp 2056〜)。無限回の sweep の解は (L/ω + N)⁻¹(−R) (ω = 0.7) で、(L + N)⁻¹ ではない。
- **既存の診断**: `FORGE_LINE_DUMP_DIR`・`_CALL`・`_NODES` で、指定した節点を含むラインの全節点について、factor の直前に D (5×5、時間項込みの対角)・Kprev・Knext・dt_local・体積・状態を、続く各 sweep の後に rhs (ライン外の寄与込み)・dq_new・dq_old を float64 で書く。この FP64 バイナリに入っている (strings で確認)。

## 2. 私の見立て (未確認)

- H2 (優先): ライン外 (流れ方向) の結合を Jacobi 型の遅れで 5 回だけ解くので、境界層の縁で解き残す。方向別 dt の点では V/Δτ が ライン外の面の λ だけで決まり、ライン外の結合に対して対角の余裕が小さい (おおよそ CFL/(1+CFL) の縮小率)。流れ方向に粗いと AR が大きく上限 50 に当たる節点が増え、時間項が大きくなる → 解き残しが減る、という dx 依存の筋。上限 1 も同じ向き。
- H2': sweep の中の relax で、解く行列そのものが (L/ω + N) にずれる。解き残しが小さくても、この行列が遅いモードに対して悪い前処理になる。

## 3. 提案 (M-d、縛られなくてよい)

1 外反復だけ、`nStepInner 40` で回して全 sweep を書き出す (毎 step dq = 0 から始まるので、sweep 0〜4 は本番の 5 回と同じ)。変えるのは nStepOuter・nStepInner・書き出しの環境変数だけ。

| 腕 | 格子・設定 | 状態 |
| --- | --- | --- |
| d1 | Gref、B0 (上限 50) | `run_0619` res_60000 (B0 で停滞した場) |
| d2 | Gref、B0 だが上限 1 | 同じ |
| d3 | Grefc、B0 (上限 50) | `run_0625` res_60000 (Grefc の B0 で 6 万 step) |

書き出すライン: 縮流部 x ∈ [−8, −2] r_t の列を 0.05 r_t ごとに 1 本 (各腕 約 120 本 × 160 節点)。

指標 (ρ・ρE を ρ_ref・ρ_ref a_ref² で無次元化、領域は島 = x −6〜−3.25・壁から 60〜100 節点目、と縮流部の全ライン):
- 解き残し u5 = ‖δ₅ − δ₄₀‖ / ‖δ₄₀‖。
- sweep あたりの縮小率 κ = k = 10〜39 の ‖δ_{k+1} − δ_k‖ / ‖δ_k − δ_{k−1}‖ の中央値。
- ライン外の寄与の大きさ ‖rhs₄₀ − rhs₀‖ / ‖rhs₀‖。
- 時間項の割合 (V/Δτ) / D₀₀ (壁際 10 節点・島)。
- 40 回で収束していること (‖δ₄₀ − δ₃₉‖ / ‖δ₄₀‖ ≤ 1e-3) を前提条件に、満たさなければ判定不能。

読み方の案: d1 の島の u5 ≥ 0.3 かつ d2・d3 の u5 ≤ d1 の 1/3 → H2 の中間量を支持 (因果は次の試験 M-b `nStepInner 40` の 2 万 step で確かめる)。d1 の島の u5 ≤ 0.1 → H2 を退け H2' (解く行列) を優先。

## 4. 判断してほしいこと

推奨を 1 つに絞り (3 腕を 1 組としてよい)、事前の読み方 (閾値) まで示してほしい。M-d が因果を言えないなら、M-b (`nStepInner 5 → 40` で 2 万 step、対照 `run_0619`) を先にすべきか、M-d → M-b の順がよいかも判断してほしい。sweep の中の relax (H2') を切り分ける腕 (relax 1.0) を足すべきかも。
