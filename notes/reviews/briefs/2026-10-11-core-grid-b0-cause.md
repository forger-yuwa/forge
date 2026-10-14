# 諮問ブリーフ: B0 (ライン陰解法) が止まるべき値から離れた所に止まる原因の切り分け

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.27〜§4.33、§5.1 #13 (q))
- **date**: 2026-10-11
- **諮問の条件**: AGENTS.md「エスカレーション条件」1 (新しい試験の設計)・4 (登録に無い run の前)・「原因は○○」と書く前
- **ユーザの依頼** (2026-10-11): 「原因が知りたいんだよな。打つ手無しって感じ？」

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
sed -n '130,200p' solver_density_cuda/cuda_forge/setDT_d.cu          (方向別 dt と上限)
sed -n '360,380p' solver_density_cuda/cuda_forge/setDT_d.cu
sed -n '2025,2060p' solver_density_cuda/main.cpp                     (残差の記録: step 0 の outer_begin が更新前の残差か)
grep -n '^### 4.3[1-3]' plans/active/tooling-nozzle-core-grid.md      (その節を sed で 40 行ずつ)
grep -n 'dt_local\|dtLocal' solver_density_cuda/cuda_forge/*.cu | head -30   (擬似 dt が右辺に入る箇所があるか)
```

## 1. 観測事実

- B0 = `lineImplicit 1`・`lineDtDirectional 1`・`implicitThermalJacobian 5`・`lineDtDirectionalCap 50`・cfl 4・`implicitRelax 0.7`。point = ライン 0・方向別 dt 0・ITJ 0・上限キーなし。FP64 `~/forge-fgeom7-fp64`。
- Gref (4719 × 160): B0 は残差の床 rms_ro 1.6e-5 で量は STEADY。point で通算 20 万 step 解くと rms_ro 7.05e-8 で横ばい、量は資格あり (B0 の窓より θ_r −0.57〜−0.63 %・Q_w −0.51 %)。最後の場 `run_0614_cg_gref_point_ext2/res_80000.h5` (AWS A)。
- Grefx (3086 × 160、流れ方向だけ粗い): B0 で残差 5 桁 (1.2e-7)、point と 0.002 % 以内。
- §4.33: Gref の B0 の場から上限だけ 50 → 1 にすると、壁際の擬似 dt が 1/50 になり、平均流の残差は 2 万 step で上限 50 の 0.6〜4 %、θ_r は point の止まった値へ動き出す (差 0.63 → 0.29 %)。ω の残差は 3〜13 倍、末尾で RISING。
- 他セッション: forge-8a (G0、段 ③): 左辺 (粘性・熱伝導の近傍の Jacobian) で止まる位置が 0.7 % 変わり、残差の床も 4〜5 倍変わる。方向別 dt は壁の第一内部節点で上限なしだと point の約 9000 倍、上限 50 は 605 節点中 231 で効く。M6 (G0): point で θ_r が B0 から −1.1〜−1.5 % 動く (3 万 step、まだ動いている)。B0 では float だけが +0.9〜1.9 % 動いたが、point では float と FP64 が一致。

## 2. 期待値と出典

- 定常解 R(U) = 0 は左辺 (反復の方式) に依らないはず。B0 が別の位置で止まる理由は、(a) B0 は収束しておらず、揺れ (有界の振動・遅いモード) の平均がずれる、または (b) B0 の設定 (方向別 dt・上限・ITJ・ライン) が右辺 R そのものを変える、のどちらか (私の見立て)。
- ユーザは原因を知りたい。計算は軽く。

## 3. 判断してほしいこと

推奨を 1 つに絞り、事前の読み方まで示してほしい。候補 (縛られなくてよい):

- S1 (右辺の比較、数分): `run_0614` res_80000 から同じ格子へ restart し、B0 と point の設定でそれぞれ数 step。step 0 の outer_begin の全 rms (更新前の残差) を比べる。同じ大きさなら (b) を退ける。B0 の方が大きければ (b) を支持し、擬似 dt が右辺に入る箇所を探す。
- S2 (安定性、数十分): 同じ場から B0 と point で各 1 万 step、500 step ごとに `res_ro` の場を出す。B0 が離れていくなら、B0 の反復が解の近くで不安定 (増える揺れがある) ことを支持し、残差が増える場所 (壁からの距離・x・変数) を場で見る。
- S3 (しきい): 同じ場から上限 2・5・10・50 で各 5000 step。どの上限から離れるか。
- 私の仮説 (未確認): ライン陰解法は壁法線の方向だけを陰に解き、流れ方向の結合は陽 (遅れ) のまま。そこに 50 倍まで伸ばした擬似 dt を組み合わせると、流れ方向に細かい格子で増える揺れが出る。Grefx (流れ方向に粗い) で問題が出ないことと整合する。
- 費用: Gref は単独 42 ms/step。A に同時 2 本まで。
