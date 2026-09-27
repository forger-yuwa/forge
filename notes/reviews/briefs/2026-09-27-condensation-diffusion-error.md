# 諮問ブリーフ: 凝縮域の拡散の近似誤差 — 測定の解釈と次の一手 (2026-09-27)

関連 plan: `plans/active/condensation-two-phase-transport.md` (§1–§4)。

## 問い (1 つ)

下の実測 (既存 NS + 凝縮 run の場から後処理で求めた流束比) の**解釈が正しいか**、そして
**「目標モデル (蒸気は蒸気の勾配で分子+乱流拡散、液は乱流拡散のみで h_l を運ぶ) に変えると結果 (onset・g・壁の量) が動く大きさ」を
最小のコストで確かめる A/B を 1 つ**設計してほしい。解釈に誤りがあれば、測り直すべき量を 1 つに絞って示してほしい。

## 観測事実

### コード読み取り (2026-09-27)

1. 総水分 `roY_w` は Fick 拡散 `J_w = −ρ(D_w + D_t)∇Y_w` (D_w: kinetic 混合平均 `speciesDiffusionMethod 1` 既定, D_t = μt/(ρ Sc_t), Sc_t 既定 0.7) と
   ΣJ=0 補正 (面ごと, `speciesTransport_d.cu:254-280`) を受け、エネルギーに `h_s J_s*` を加える。`viscMethod != 0` のとき (`:813`)。
2. **液 `rog` (凝縮モーメント) は拡散しない**: 受動拡散 `passiveDiffusion_d_wrapper` (`speciesTransport_d.cu:1294`) はトレーサ (`tracerTransport_d.cu:176`) からしか呼ばれない。
   case/16 README の run_0482 行も「モーメントに拡散は入っていない」と記す。
3. 実現可能性クランプ `0 ≤ rog ≤ roY_w` が毎 step 入る (`condensationRealizability_d.cuh:92-113`)。
4. S3 (`speciesFaceReconstruction 2`) の面組成は種ごとに `limiter_ro` で再構成 → `[0,1]` クランプ → ΣY=1 に正規化 (double で) → upwind 側を `Yface_out` に保存
   (`convection/convectiveFlux_slau_d.inc.cuh:340-369, 582-586`)。化学種移流は `ṁ·Y_f` (`passiveKernels_d.cuh:10-35`)、node 境界半割面は所有セルの Y。
   → 移流の ΣY は float 丸めの範囲で連続の式と整合 (当方の読み取り)。

### 当方の解釈

- 拡散で動くのは実質「蒸気」(g は不変) で h_v を運ぶので、**エネルギー勘定は蒸気として整合**。問題は駆動勾配が `∇Y_w` であること。
- case/16 Wysłouzil は入口組成一様なので `Y_w` はほぼ一様 → **現行では水の拡散流束がほぼ 0**。一方、目標モデルでは
  - 蒸気の分子+乱流流束 `−ρ(D_w + D_t)∇Y_v` (≠0; `∇Y_v ≈ −∇g`)
  - 液の乱流流束 `−ρD_t∇g`
  - 乱流部分は質量が打ち消し (`∇Y_v + ∇g = ∇Y_w ≈ 0`)、**エンタルピーは `(h_l − h_v)(−ρD_t∇g) = L ρ D_t ∇g` が残る** (乱流による潜熱の輸送)。
  - 分子部分は総水分の正味流束 `−ρD_w∇Y_v` が残る (蒸気が凝縮域へ向かって拡散する)。

### 後処理の実測 (`case/16.nozzle_wys/analyze_liquid_diffusion_error.py`, 最終スナップショット res_48000)

定義: ノード勾配は VIZMESH の四角形を三角形に分割した P1 勾配の面積重み平均。D_w は H2O–N2 二元 Chapman–Enskog (Neufeld Ω11; 希薄 H2O なので混合平均 ≈ 二元)、
D_t = μt/(ρ·0.7)、熱流束 Q = (λ + cp μt/Pr_t)|∇T| (cp 1040, Pr_t 0.9)、L = 2.5 MJ/kg。対象は g > 1e-4 のノード。

| run | 設定 | g>1e-4 ノード (うち μt/μ>10) | E1/F1 (余計な蒸気分子流束 / 蒸気分子流束) | L·ρD_t|∇g| / Q 全体 p50/p95/p99/max | 同 μt/μ>10 のみ p50/p95/max |
|---|---|---|---|---|---|
| `case/16.nozzle_wys/run_0483_passive_wys_s0_sfr0` | SST TP MIXDRY(=N2)+H2O, 非平衡凝縮, 旧経路 (scheme 0 / SFR 0), 48000 step | 9742 (3256) | **1 / 1 / 1 / 1** | 0.255 / 4.17 / 6.78 / 138 | 0.768 / 4.78 / 5.75 |
| `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1` | 同, S3 (scheme 1 / SFR 2 / coupling 1) | 9368 (2992) | **1 / 1 / 1 / 1** | 0.21 / 4.58 / 6.93 / 15.9 | 0.747 / 5.32 / 6.51 |

- E1/F1 = 1 は `|∇g| = |∇Y_v|` (Y_w 一様) を反映 = 現行の蒸気拡散流束は 0 (余計な分が真の流束をちょうど打ち消す)。
- run_0483 の max 138 は x 20.87 mm, y 0.081 mm (wall_dist 2.7 mm, g 5.9e-4, T 213.5, μt/μ 2.1) — ∇T が小さい点の比の発散と思われる。
- 液の乱流混合の相当温度スケール L·g_max/cp ≈ 26 K (g_max 0.0109)。
- 各 run の判定 (既存記録, case/16 README:380): run_0482 `check_convergence` PASS、run_0483 NOT CONVERGED (rms_roe 2.9 桁)、両方 `check_quasisteady --series-csv` ALL STEADY。
- 比は**流束の大きさの比**であり、場への影響 (T・g・onset・壁量の変化) ではない。場への影響は未測定。

## 期待値と出典

- 目標モデルは plan §4.2 (ユーザ決定 2026-09-27)。液滴の懸濁効果は φ ≈ 2e-6 で無視 (plan §3)。
- case/16 は Wysłouzil et al. 2000 の実験比較ケース。forge の onset は実験より ~5 mm 下流 (recommended-settings §3)。

## 実施済みの操作

- 上記の後処理のみ。コード変更・新規 run なし。

## 仮説

- H1: 乱流域で欠けている潜熱の乱流輸送 `L ρ D_t ∇g` が乱流熱流束と同程度 (p50 0.75、p95 ~5) なので、目標モデルに変えると境界層の縁の T・g 分布が数 K・数 % 動く。
- H2: 比が大きいのは ∇T が小さい点で比が発散しているだけで、エネルギー収支上の寄与 (発散量) は小さい。流束比でなく**流束の発散の比**で見るべき。
- どちらか当方では決めきれない。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックス)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。h5 を読む Python の実行もしない。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か本ブリーフの数値で示す。

## 読んでよいもの

- `sed -n '1,80p' case/16.nozzle_wys/analyze_liquid_diffusion_error.py`
- `sed -n '179,285p' solver_density_cuda/cuda_forge/speciesTransport_d.cu`
- `sed -n '795,830p' solver_density_cuda/cuda_forge/speciesTransport_d.cu`
- `sed -n '1285,1310p' solver_density_cuda/cuda_forge/speciesTransport_d.cu`
- `sed -n '165,185p' solver_density_cuda/cuda_forge/tracerTransport_d.cu`
- `sed -n '90,125p' solver_density_cuda/cuda_forge/condensationRealizability_d.cuh`
- `sed -n '340,370p' solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh`
- `sed -n '578,590p' solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh`
- `sed -n '10,40p' solver_density_cuda/cuda_forge/passiveKernels_d.cuh`
- `sed -n '1,120p' plans/active/condensation-two-phase-transport.md`
- `sed -n '380,381p' case/16.nozzle_wys/README.md`
