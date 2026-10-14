# 諮問ブリーフ: B0 の遅いモードの仕組み — 平均流の側の候補を切り分ける次の 1 手

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.33〜§4.39)
- **date**: 2026-10-11
- **諮問の条件**: AGENTS.md「エスカレーション条件」4 (登録に無い run の前)
- **ユーザの依頼**: 「なぜ流れ方向に細かいとライン陰解法の反復が遅くなるのか、という仕組みです。つきとめてくれ」。ユーザの提案「低マッハ補正をかけるべきなのかな」(§4.38)。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
grep -n '^### 4.3[3-9]' plans/active/tooling-nozzle-core-grid.md       (その節を sed で 30 行ずつ)
sed -n '130,200p' solver_density_cuda/cuda_forge/setDT_d.cu            (方向別 dt と上限)
grep -n 'implicitThermalJacobian' solver_density_cuda/cuda_forge/timeIntegration_d.cu | head -20
sed -n '1020,1040p;1155,1175p' solver_density_cuda/cuda_forge/timeIntegration_d.cu   (ITJ のビット 1・2)
grep -n 'nStepInner' solver_density_cuda/main.cpp | head
```

## 1. 観測事実 (§4.33〜§4.39)

- B0 = `lineImplicit 1`・`lineDtDirectional 1`・`implicitThermalJacobian 5`・`lineDtDirectionalCap 50`・cfl 4・`implicitRelax 0.7`・`nStepInner 5` (FP64)。Gref (4719 × 160) で残差プラトー (rms_ro 1.6e-5)、point の止まった値から θ_r +0.6 %・Q_w +0.5 % 離れて止まる。
- 右辺は不変 (上限・方向別 dt は `setDT` だけ、`dt_local` は組立に入らない)。交互更新なし、増分の消失・床の切り詰めなし (D1)。
- 遅いモードの場所と形: 縮流部 x −7〜−2 の冷却壁の乱流・温度境界層の縁 (壁から 60〜100 節点目、物理的には壁の近く)。主役はエントロピー (s の RMS が圧力の 18 倍)、速度は縁から主流へ斜めの縞。M 0.02〜0.08、温度の前線 1100〜1300 K → 1600 K、μt/μ 50〜150 → 8〜11。
- **効いたもの**: 上限 50 → 1 (壁際の方向別 dt を 1/50、主流の dt はほぼ不変): 縮流部の残差 1/25〜1/160 (ω は 3〜13 倍)。縮流部 x −8〜−2 だけ流れ方向を粗く (Grefc、列 589 → 337): 縮流部の残差 ρ 1/36・ρE 1/62、量は point の止まった値から 0.05 % 以内。
- **効かなかったもの**: SST の凍結 (F1、`FORGE_FREEZE_TURB=1`、ρk・ρω の更新を止める): 縮流部の残差の比 0.94〜1.12 (結果 B)。
- **他セッション** (forge-8a、G0 の格子・段 ③ のバイナリ): 左辺に熱伝導の近傍の Jacobian K を入れる (ITJ の値 3・マスク 7、上限 50) と止まる位置が θ_r +0.5〜0.7 %・Q_w +0.7〜0.8 % ずれて床が 3.3 倍。上限なしでは熱伝導の近傍 K があると 300 step 前後で発散、外すと有限。
- 低マッハ前処理 (`lowMachPrecond 3`) はライン陰解法と組み合わせ不可、エントロピーの固有値は Uₙ のままで純粋なエントロピーの成分には恒等 (§4.38 の諮問)。

## 2. 期待値と出典

- 私の見立て (未確認): 平均流の側。候補は (a) 熱の Jacobian (ITJ 5) が壁際の大きな方向別 dt と組んで、温度境界層の縁のエントロピーのモードを遅くする、(b) ライン外 (流れ方向) の結合を Jacobi 型の遅れで扱うので、流れ方向に細かい所でライン外の結合が強く (セルが流れ方向に細く)、有限回 (5 回) の内部反復で解き残す、(c) その両方。

## 3. 判断してほしいこと

推奨を 1 つに絞り (2 本を同時に回すなら、その組を 1 組として)、事前の読み方まで示してほしい。候補 (縛られなくてよい):

- M-a: `run_0608` res_60000 から B0 のまま `implicitThermalJacobian 5 → 0` だけを変えて 2 万 step。
- M-b: 同じ起点から B0 のまま `nStepInner 5 → 20` だけを変えて 2 万 step。
- M-c: 同じ起点から B0 のまま `lineDtDirectional 1 → 0` (方向別 dt をやめる、上限キーも外す。ライン被覆 100 % なので全節点が point の dt になる) だけ。上限 1 とほぼ同じになる見込み。
- 対照はいずれも `run_0619` の 0〜20k (領域の残差あり)。判定は縮流部の ρ・ρE の比 (≤ 1/3 / ≥ 2/3)。
