# G0/G1/G2/G2′: node 並進周期の継ぎ目勾配 (plan boundary-node-periodic-gradient-fix §6)

32³ の三重周期 TGV メッシュ (`case/09.Taylor-Green/mesh`、node 変換) などに線形・二次場を焼き込み、1 step の勾配を
解析値・CPU 参照と比べる。run はすべてスクラッチ (`/tmp/claude-1000/.../scratchpad/g_harness_*`、破棄可)。

## ファイル

| ファイル | 内容 |
| --- | --- |
| `bake_linear.py` / `aggregate.py` | 初版の G0 (Ux=10+y、x 継ぎ目の接線成分のみ) |
| `solverConfig.yaml` / `bcondConfig.yaml` / `solverConfig_sst.yaml` | 初版の設定 (TGV の周期 bcond は以後も共通で使う) |
| `G0_translational.txt` / `G0_translational_m1.txt` | 初版の結果 (修正前 `3de26cba` / 修正後 M1) |
| `gharness.py` | 共通部: h5 読み込み、周期 group の再構成 (面重心 + 並進で相手面を照合、root = 最小 index)、合併 stencil LSQ の double 参照、node GG の float32 再現と double、場の焼き込み、run 作成 (`run_case.sh` を block 128 で起動) |
| `mkmesh.py` | gmsh 構造格子 → 節点ジッタ (±0.2h、index を周期で折り返したハッシュ = 周期像は同じ量) → node 変換 → `check_mesh_quality.py` |
| `g_suite.py <variant>` | G0 (解析値) / G2 (CPU double 参照) / G1-a・G1-b (GG) / 4 量の床 を 1 run で測る。変種は下表 |
| `g0_cpu_cases.py` | G0 の CPU 作用素試験: float32 反例 (root 両順序)、部分 CV の rank 欠損 → 合併で回復、合併後も退化 → 打ち切り参照 |
| `g2p_jitter.py` | G2′: ジッタ格子 16³/32³/64³ の二次場で次数・継ぎ目/内部比・ゼロ成分、一様 32³ の別行 |
| `f1_read_test.py` | F1 の読み出し (codex m3): 1・2 回目の残差組立で輸送が計算値 F1 を読むか |
| `G_<variant>.txt` / `G0_cpu_cases.txt` / `G2p_jitter.txt` / `G_f1_read.txt` | 結果 (各ファイル末尾に VERDICT) |

場は ρ=1 で焼く (channel は ρ も線形)。w は継ぎ目中心の局所座標 (周期軸を中央で折り返し、group 内は root の値)。
折返し不連続に触れる節点は解析値比較から除く。状態は `res_0.h5`、勾配は `res_1.h5` (k/ω の勾配は res_0 では未計算)。

## 結果 (2026-09-26、バイナリ sha256 `bb55e1bf…`、HEAD `1f03c036`)

| 試験 | 変種 / メッシュ | VERDICT | 主要数値 |
| --- | --- | --- | --- |
| G0 | tgv (一様、2/4/8 member) | PASS | 非零成分 最大誤差/\|a\| ≤ 1.34e-6、ゼロ成分 0、ρ=1 の勾配 0 |
| G0 | tgv_mirror (点対称 = root を継ぎ目の反対側へ) | PASS | 同上。(group, 継ぎ目軸) 3072 組で root が下側にあるのは tgv 2080 / mirror 992 = 全組で側が入れ替わった |
| G0 | tgv_shift (原点 (100,100,100)) | PASS | ≤ 8.1e-7 |
| G0 | jitter32 (非対称 stencil) | PASS | 非零成分 ≤ 1.62e-6、ゼロ成分 ≤ 6.2e-8 |
| G0 | channel (x 等比・y 両壁、壁∩継ぎ目) | PASS | ρ 線形で壁∩継ぎ目 (面・辺) も ≤ 1.5e-6 |
| G0 | CPU: float32 反例 | PASS | 実変位 (現行) は root 両順序で 1.000000。代表変位 (旧) は 0.999980 / 1.000020 |
| G0 | CPU: rank 欠損 → 合併で回復 | PASS | 合併 4.7e-8、部分ごとの和 (旧) は y を 2 重計上 (誤差 0.85) |
| G0 | CPU: 合併後も退化 | PASS | 打ち切り参照との差 4.8e-8 |
| root 順序 | tgv_bcswap / tgv_repeat | PASS | bcond の記述順では root は変わらない (最小 index)。run 間差は atomicAdd 順序の床 ≤ 1.75 ε·max\|φ\|/h |
| G2 | 全 7 run (tgv・mirror・shift・bcswap・repeat・jitter32・channel) | PASS | 最大差/S ≤ 1.4e-7 (折返し不連続を除いた S_clean で ≤ 2.4e-7)。G2′ の 3 水準でも ≤ 2.8e-6 |
| G1-a | 全 7 run | **FAIL** | 継ぎ目で 4.7–27 ε·max\|φ\|/h (内部は ≤ 1.7)。下の「G1-a の不合格」 |
| G1-b | 全 7 run | PASS | 継ぎ目 ≤ 26.5 ε (閾値 2·N_max·n_member = 40)。jitter32 の ω 面は 19.75 / 閾値 20.0 |
| G2′ 次数・比 | ジッタ 16/32/64 | PASS | 継ぎ目の次数 0.95, 0.93、継ぎ目/内部 1.00 / 0.85 / 0.83 |
| G2′ ゼロ成分 | ジッタ 16/32/64 | **FAIL** | ∂ψ/∂z (ψ は z に依らない二次場) が 1.5e-2 → 8.4e-3 → 4.3e-3 (O(h) で収束、床 1.7e-6 を超える)。一様 32³ は 0 |
| F1 読み出し | TGV 32³、F1 計算値 0 / 1 | PASS | 1・2 回目とも Δres_ω/Δres_k = 151.893 / 151.895 (期待 151.893) |
| R3 (2026-09-26、旧 `1266aba1` sha `9b45d018…` / 新 sha `88e949fa…` = §4.2a 修正後) | case48 (fp_y1_12um、SST、一様 IC) / case16 (run_0507 の発達場、speciesFaceReconstruction 1) / axi (101×41 軸対称 node、MIXDRY/H2O、焼き込み場。case/44 のメッシュ削除済みの代替) を旧 2・新 2 で 1 step。`r3_prepare.py` → `r3_compare.py` → `R3.txt`、run はスクラッチ `r3_<case>_<old\|new>_<a\|b>` | PASS (3 ケースとも) | res_1 の勾配・リミタ 30 配列: 旧同士ビット一致の配列は旧新もビット一致。case48 の dK/dΩ (GG) は旧同士でも 62–544 点不一致、旧新は 49–570 点・最大差は旧同士と同値 (ノイズ床)。NaN/Inf なし。化学種 ∇Y は出力に無く直接は比べていない (res_1 の roY は旧同士と同程度の不一致)。case48 の res_1 roK/roOmega は旧新で 5.1 万点不一致 (旧同士 0/2 点) |

### 4 量の床 (GPU − CPU double、ε·max|φ|/h 単位、tgv)

| 量 | 内部 最大 | 継ぎ目 最大 |
| --- | --- | --- |
| dUx (LSQ) | 0.4 | 0.4 |
| dK (GG) | 1.1 | 20.8 |
| dΩ (GG) | 1.3 | 26.5 |
| dξ (GG、化学種 Y の代理) | 0.9 | 22.7 |

化学種 Y は `dY{s}d*` が出力変数に無い (`variables.hpp` の output_cellValNames) ので測れない。受動トレーサ ξ は同じカーネル
(`species_gradient_d`、`speciesTransport_d.cu:1229`) と同じ合併 (`periodicGradientGather`) を通るので代理にした。

### 修正後 (2026-09-26、plan §4.2a、バイナリ sha256 `88e949fa…`)

`G_*.txt` 7 本は修正後の再実行で上書きした (修正前は git の `633997bf`)。G0・G2・G1-a・G1-b すべて PASS。
継ぎ目の床 (tgv、ε·max|φ|/h) dK 1.7 / dΩ 2.4 / dξ 1.5 (修正前 20.8 / 26.5 / 22.7)、継ぎ目/内部比 0.83–1.85。上の表の G1-a 行は修正前。

### G1-a の不合格: 周期半割面が GG から除外されていない (修正前、上記で解消)

`ransTransport_d.cu:43` と `speciesTransport_d.cu:649` の除外条件は `ip >= nNormalPlanes && ic1 < nCells` だが、実行時の
`plane_cells` では周期半割面にも ghost が割り当てられ (`mesh.cpp:443-456`、全 bcond の iBPlanes に ghost を作る)、
`ic1 >= nCells` になる。したがって周期半割面は除外されず、node の境界面規則 (面値 φ[ic0]) で積算されている。

- CPU の float32 再現に「周期半割面も φ[ic0] で積算」を足すと、GPU との差は継ぎ目でも ≤ 1.75 ε·max|φ|/h (全変種、`G_*.txt` の診断列)。
- 合併後に残るのは対の半割面の面ベクトル不一致 (S_a + S_b)·φ/V。tgv で max|S_a+S_b|/|S_a| = 3.7e-6 → 定数場の GG 勾配が継ぎ目で
  3.7e-6·φ/h (≈ 31 ε·φ/h)。**φ に比例** (Δφ ではない) し、float32 演算の丸めではない。初版で「継ぎ目法線成分の ~20 ulp·φ/h の丸め」
  と書いたもの (`G0_translational_m1.txt`) はこれ。
- NS の GG 経路 (`calcGradient_b_d`) は host 側で `bcondKind == "periodic"` を飛ばしている (`calcGradient_d.cu:1091`)。

### G2′ ゼロ成分の不合格について (事実のみ)

ψ = w·(a_x, a_y, 0) + ½κ(w_x² + ½w_y²) は z に依らないが、ジッタ格子の非対称 stencil では LSQ の z 成分に二次項の打ち切り誤差が
乗り、他の成分と同じく O(h) で収束する (1.5e-2 → 8.4e-3 → 4.3e-3、継ぎ目・内部とも)。一様格子 (対称 stencil) では 0。
§6 G2′ の「ゼロ成分は絶対誤差 ≤ e_floor」は一様格子・線形場でしか成り立たない条件で、ジッタ格子の二次場には当てはまらない。
次数・継ぎ目/内部比の判定 (PASS) とは別の行にしてある。

## run (スクラッチ、破棄可)

`/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad/` の下:
`g_harness_{tgv,tgv_mirror,tgv_shift,tgv_repeat,tgv_bcswap,jitter32,channel}`、
`g_harness_g2p_k0.2_N{16,32,64}_j0.2`、`g_harness_g2p_k0.2_N32_j0`、`g_harness_f1_{inf_b0,inf_b1,inf_b1_repeat,tiny_b0,tiny_b1}`、
メッシュは `g_harness_mesh/{box16_j0.2,box32_j0.2,box64_j0.2,box32_j0,box32_j (jitter32 用),box16_j (動作確認のみ),channel}` (各 `quality.txt` は PASS)。
