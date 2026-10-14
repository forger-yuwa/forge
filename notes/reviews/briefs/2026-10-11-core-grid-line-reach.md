# 諮問ブリーフ: ライン陰解法 (B0) のまま point の止まった値に届かせる設定を、軽い計算で探す試験の組み方

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.27〜§4.31、§5.1 #13 (p))
- **date**: 2026-10-11
- **諮問の条件**: AGENTS.md「エスカレーション条件」1 (新しい試験の設計)・4 (登録に無い run の前)
- **ユーザの依頼** (2026-10-11): 「速度は測る必要がなく、line で同じ結果を出すにはどうしたらいいか、を軽い計算で確かめてほしい」

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
sed -n '60,110p' procedures/solver-settings.md                  (lineImplicit・lineDtDirectional・lineDtDirectionalCap・implicitThermalJacobian の説明)
grep -n '^### 4.2[7-9]\|^### 4.3[01]' plans/active/tooling-nozzle-core-grid.md   (その節を sed で 30 行ずつ)
grep -rn 'lineDtDirectionalCap\|lineDtDirectional' solver_density_cuda/cuda_forge/*.cu | head -20
```

## 1. 観測事実

- B0 = `lineImplicit 1`・`lineDtDirectional 1`・`implicitThermalJacobian 5`・`lineDtDirectionalCap 50`・cfl 4・`implicitRelax 0.7`・FP64 `~/forge-fgeom7-fp64`・`mesh.axisSegmentRWeight 1`。point = `lineImplicit 0`・`lineDtDirectional 0`・`implicitThermalJacobian 0`・上限キーなし (ほかは同じ)。
- Gref (4719 × 160、75.5 万節点、単独 41.8 ms/step):
  - B0 (`run_0607` → `run_0608`): 残差プラトー rms_ro 1.6e-5、`check_convergence` NOT CONVERGED、量は STEADY (θ_r の U 0.09〜0.17 %)。
  - point (`run_0608` の res_60000 から `run_0611` → `run_0613` → `run_0614`、通算 20 万 step): 残差 2.4 桁 (7.05e-8 でプラトー)、rms_roOmega は途中で 9e1 → 3.2e3 → 1.1e2。量は資格あり (8 量 STEADY、d ≤ 0.01 %) で、B0 の窓より θ_r [65,75] −0.63 %・[84,94] −0.57 %・Q_w −0.51 %・δ_loc −0.03〜−0.04 %・M ≤ 0.004。最後の場は `run_0614_cg_gref_point_ext2/res_80000.h5` (AWS A)。
- Grefx (3086 × 160): B0 で残差 5 桁 (1.2e-7)、point と 0.002 % 以内。止まった Gref との差は θ_r ≤ 0.05 %・Q_w 0.04 %。
- G0・G1 (ni 4719) も B0 で残差プラトー・揺れが大きかった。G1x (ni 3086 だが半径方向が粗い) も揺れが大きかった (§4.15)。
- 別セッション (forge-8a) は上限 50 で「全部入りは止まる状態」という比較をしている (内容は未確認)。

## 2. 期待値と出典

- 定常解 R(U) = 0 は陰解法の左辺 (反復の方式) に依らないはず。B0 が point と違う位置で STEADY に見えるのは、(i) B0 が収束しない状態 (遅いモード・リミットサイクル) に留まっている、(ii) B0 の設定が右辺 (残差) を変えている (擬似 dt に依る項がある)、のどちらか (私の見立て、未確認)。
- ユーザは軽い計算を望む (速度は測らない)。

## 3. 判断してほしいこと

推奨を 1 つに絞り、判定の条件まで示してほしい。候補 (縛られなくてよい):

- T1 (固定点の試験): point で止まった Gref の場 (`run_0614` res_80000) から B0 で 2 万 step。止まった値に留まるか、B0 の窓の値へ戻るか。留まる → (i) を支持 (B0 は止まった解を保つが、そこへ行けない)。戻る → (ii) を支持 (B0 の右辺か反復が別の状態を作る)。
- T2 (1 つずつ外す): B0 の窓の場 (`run_0608` res_60000) から、`lineDtDirectional 0`・`implicitThermalJacobian 0`・`lineDtDirectionalCap` 外しを 1 つずつ、各 2〜4 万 step。どれで量が止まった値へ動き出すか。
- 費用: Gref は 2 万 step で約 14 分 (単独)。A に同時 2 本まで (メモリ 3.3 GB/本)。
