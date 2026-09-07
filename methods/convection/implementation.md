# 対流フラックス — 実装

forge の対流フラックス計算の実装と、ソース上の対応関係をまとめる。
理論的背景は [theory.md](theory.md) を参照。

## ソースファイル

| ファイル | 役割 |
| --- | --- |
| [`solver_density_cuda/convectiveFlux.hpp`](../../solver_density_cuda/convectiveFlux.hpp) | 旧 CPU 経路用ヘッダ (現在ほぼ宣言のみ・コメントアウト) |
| [`solver_density_cuda/convectiveFlux.cpp`](../../solver_density_cuda/convectiveFlux.cpp) | 旧 CPU 実装 (Roe 平均、ヤコビアン構築) |
| [`solver_density_cuda/cuda_forge/convectiveFlux_d.cuh`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cuh) | GPU カーネル / device 関数の宣言 |
| [`solver_density_cuda/cuda_forge/convectiveFlux_d.cu`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu) | スキーム本体・ラッパ (約 2700 行) |
| [`solver_density_cuda/cuda_forge/AUSM_d.cu`](../../solver_density_cuda/cuda_forge/AUSM_d.cu) | AUSM 系の補助関数 (M±, β± など) |

実運用は GPU 経路のみ。CPU の `convectiveFlux.cpp` はリファレンス用途で、
通常のビルドでは呼ばれない。

## エントリポイント

[`convectiveFlux_d_wrapper`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L2518)
が共通入口。`cfg.solver` 文字列でディスパッチする。

| `cfg.solver` | 呼び出すカーネル | 状態 |
| --- | --- | --- |
| `"SLAU"` | `SLAU_d` (`slauVariant=1`) | 有効 |
| `"SLAU2"` | `SLAU_d` (`slauVariant=2`) | 有効 (圧力束のみ SLAU2、低マッハ改良) |
| `"HLLE"` | `HLLE_d` | 有効 |
| `"ROE"`  | `ROE_d`  | 有効 |
| `"KEEP"` | `KEEP_d` | 有効 (純粋 KEEP 中心流束・散逸なし。LES/ILES 向け。cell/node 両対応) |
| `"AUSM+"`, `"AUSM+UP"`, `"KEEP_SLAU"` | (legacy に実装あり) | ラッパに分岐なし・到達不能 (`legacy/`) |

呼び出し直前に `res_ro`, `res_roUx`, `res_roUy`, `res_roUz`, `res_roe` を
`cudaMemset` でゼロクリアする。

## 状態再構成と dispatch

すべてのスキームが共通に [`interp_dispatch`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L142)
を呼び、`cfg.convMethod` (= `scheme`) と `cfg.limiter` (= `limit_scheme`) に応じて
次の device 関数を選ぶ。

| `scheme` | 関数 | 内容 |
| --- | --- | --- |
| `0`, `-1` | `interp_1stUp` | 1 次風上 (`phif = phiC`) |
| `1` | `interp_MUSCL_2nd` | 線形再構成: $\phi_C + \psi\,\nabla\phi_C \cdot \mathbf{cp}$ |
| `2` | `interp_MUSCL_3rd` | $k = 1/3$ で 3 次 MUSCL |
| (上記以外) | `interp_MINMOD` | 上流値 $\phi_U$ を構成して minmod $r$ を計算し制限 |

L 側は `(phiC, phiD) = (ro[ic0], ro[ic1])` ほか、距離ベクトル `dcc`, `dc0p` を渡し、
R 側は左右を反転 (`-dcc`, `dc1p`, `1-f`) して呼ぶ。
`limiter_*[ic]` は別途生成済みのスカラリミッタを各成分独立に渡す。

### Ducros センサの適用

ROE_d など一部のカーネルで [`apply_ducros_limiter`](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L189) を経由して

```cpp
duc = clamp(max(ducros[ic0], ducros[ic1]), 0.0, 1.0);
limiter = (duc > 0.8) ? max(0.0, (1.0 - duc) * limiter) : limiter;
```

を適用する (Ducros 値が高い、すなわち圧縮性ショック近傍で再構成を抑制)。

### 多成分 TP の face 組成整合 (`speciesFaceReconstruction`)

`convectiveFlux_d.cu` の `thermalMethod==2` (TP) 分岐で、face thermo ($R_\mathrm{mix}, T, h$) に使う
組成 $Y_{s,f}$ を切り替える (理論と検証は [theory.md](theory.md) の同名節、設計判断は
[`plans/accepted/convection-multispecies-contact-pressure.md`](../../plans/accepted/convection-multispecies-contact-pressure.md))。

| 値 | 動作 | 位置づけ |
| --- | --- | --- |
| `0` (既定) | face thermo にセル値 $Y_s$ (1 次) — mixed-order | ビット不変の従来挙動 |
| `1` (**S2**) | $\nabla Y_s$ (Green-Gauss) + **$\psi_\rho$ (ρ と同一リミッタ)** で $Y_{s,L/R}$ を再構成し thermo に使用。clamp + $\Sigma Y=1$ 正規化 | **多成分 TP の production 推奨** |
| `2` (S2+S3) | 加えて species 移流も同一 face 組成の upwind (`species_advection_faceY_d`) | **experimental・cfl≤2 限定** (1 次 LHS との defect-correction 不整合で高 CFL 発散) |

関連: `multispeciesRhoYCommonLimiter` (既定 0) は ρ–Y 共通リミッタの診断オプション、
`limiter_Y` ($\psi_Y$) は診断用残置 ($Y$ 再構成には使わない — $\min(\psi_\rho,\psi_Y)$ は発散実績あり)。

## 各カーネルの構造

共通する処理パターン:

1. `ip_orig = blockDim*blockIdx + threadIdx` から `normal_halo_planes_d[ip_orig]` で
   実際の面 ID `ip` を引く (有効面のみで並列化するための間接参照)。
2. `plane_cells[2*ip+0/1]` で両側セル `ic0`, `ic1` を得る。
3. `ic1 >= nCells` (= 非周期境界のゴースト) のとき `conv_scheme = -1` に強制し
   1 次風上に落とす。
4. `interp_dispatch` で L/R 状態を作る。`ro_L = max(ro_L, small_rho)` 等で
   負値 / 0 をクランプ (Roe のみ)。
5. スキーム固有のフラックス計算 (下記)。
6. `atomicAdd(&res_*[ic0], -res_*_temp)` / `atomicAdd(&res_*[ic1], +res_*_temp)`
   で両側に符号反転加算。

### `SLAU_d` ([L218](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L218))

並列単位は `nNormal_ghst_Planes` (`normal_ghst_planes_d[]` で実面 ID へ間接参照)。

1. `interp_dispatch` で `(ro, Ux, Uy, Uz, Ps)` の L/R 状態を構築。
   `velocity2_L/R = |u|²`, `h_p/m = γP/((γ-1)ρ) + 0.5|u|² = H_t` を計算。
2. 法線速度 `Vn_p = u_L·n̂`, `Vn_m = u_R·n̂` (`/sss` で正規化)。
3. **音速平均** `c_hat = 0.5*(sonic[ic0] + sonic[ic1])` (算術平均、Roe 平均ではない)。
4. **Mach 関数** `beta_p, beta_m` を分岐式でその場計算 (`|M|≥1` と `<1` で別式)。
   AUSM ファミリの汎用ヘルパ ([`AUSM_d.cu`](../../solver_density_cuda/cuda_forge/AUSM_d.cu)) は
   呼ばず、SLAU 用に直接ベタ書き。
5. **低マッハ補正係数**:
   ```cpp
   g            = -max(min(M_p, 0), -1) * min(max(M_m, 0), 1);
   Vn_hat_abs   = (ro_L*|Vn_p| + ro_R*|Vn_m|) / (ro_L + ro_R);
   Vn_hat_p_abs = (1-g)*Vn_hat_abs + g*|Vn_p|;
   Vn_hat_m_abs = (1-g)*Vn_hat_abs + g*|Vn_m|;
   M_hat        = min(1, sqrt(0.5*(|u_L|² + |u_R|²)) / c_hat);
   chi          = (1 - M_hat)²;
   ```
6. **圧力束** `p_tilde` と **質量流束** `mdot`。圧力束の**第 3 項のみ** `slauVariant` で
   分岐する (`mdot` は両版共通):
   ```cpp
   // 第 3 項: slauVariant==2 (SLAU2) のときだけ差し替え
   if (slauVariant == 2)               // SLAU2 (低マッハ圧力散逸)
       p_third = (beta_p+beta_m-1)*sqrt(0.5*(velocity2_L+velocity2_R))
                 *0.5*(ro_L+ro_R)*c_hat;
   else                                // SLAU
       p_third = (1-chi)*(beta_p+beta_m-1)*0.5*(P_L+P_R);
   p_tilde = 0.5*(P_L+P_R) + 0.5*(beta_p-beta_m)*(P_L-P_R) + p_third;
   mdot    = sss*0.5*((ro_L*(Vn_p+Vn_hat_p_abs) + ro_R*(Vn_m-Vn_hat_m_abs))
                      - chi/c_hat * (P_R-P_L));
   ```
   SLAU2 の根拠と式は [`theory.md` の SLAU2 節](theory.md#slau2-圧力束の低マッハ改良) を参照。
   `velocity2_L/R`, `c_hat` は既存量、新規は `0.5*(ro_L+ro_R)` のみ。
7. **残差**: `0.5*(mdot ± |mdot|)` で風上選択し、運動量に `p_tilde * S_*` を加え、
   エネルギは `h_p, h_m` (全エンタルピ) を風上にとる。
8. 両側セルに `atomicAdd(±)`。

リミッタは Ducros 補正を通さず生の `limiter_*[ic0/ic1]` をそのまま渡している
(Roe との違い)。

#### 低マッハ前処理 (`lowMachPrecond`)

`cfg.lowMachPrecond == 1`（フラックス散逸前処理）または `== 2`（RHS+LHS 完全前処理）のとき、質量流束の
圧力散逸スケール `c_hat` を**前処理音速** `c_prime` に置き換える (既定 `0` は従来挙動でビット不変)。
`== 3`（LHS のみの完全前処理）は RHS 散逸を触らず `c_hat` のまま (収束解は `0` とビット一致)。共有ヘッダ
[`lowMachPrecond_d.cuh`](../../solver_density_cuda/cuda_forge/lowMachPrecond_d.cuh) の device 関数で

```cpp
// |u|_face は L/R 速度の平均、eps は cfg.precondEps
flow_float Ur  = min(c_hat, max(vel_mag, eps*c_hat));
flow_float bta = (Ur/c_hat)*(Ur/c_hat);
flow_float c_prime = 0.5*sqrt((1-bta)*(1-bta)*Vn*Vn + 4*Ur*Ur);   // → c_hat (M≥1), → Ur (M→0)
...
mdot = sss*0.5*((ro_L*(Vn_p+Vn_hat_p_abs) + ro_R*(Vn_m-Vn_hat_m_abs))
                - chi/c_prime * (P_R-P_L));   // lowMachPrecond=0 のとき c_prime==c_hat
```

を計算する (既定 `precondEps=0.15`)。`SLAU_d` には `int lowMachPrecond, flow_float precondEps` 引数を追加し、
`c_prime` は散逸項 (`chi/·*ΔP`) にのみ用いる。`Vn_hat_p/m_abs` (速度上流項) と `p_tilde` (圧力束) は変更しない。
理論的根拠と検証 (低マッハ自励振動の振幅を ε=0.15 で −32% 減衰) は
[`theory.md` の低マッハ前処理節](theory.md#低マッハ前処理-weisssmith-型散逸スケール是正) を参照。

> 当初は同じ `lowMachPrecond_d.cuh` を `setDT_d` のスペクトル半径・block DPLUR の固有値でも共用する計画
> だったが、いずれも block DPLUR の対角優位性を崩して有害と判明し不採用 (計画 §9)。**現状フラックス散逸前処理のみ**。

#### 低マッハ再構成補正 (`lowMachThornber`)

`cfg.lowMachThornber == 1` のとき、L/R 速度再構成の直後 (まだ `velocity2_*`/`h_*`/`Vn_*` を
組む前) で左右速度ジャンプを局所マッハ `z=min(1, |u|_face/c_hat)` で縮める。`lowMachPrecond`
(圧力散逸側) とは作用する項が異なり**直交**するため独立トグルとし、既定 `0` でビット不変。

```cpp
// L/R 速度再構成 (Ux_L..Uz_R) の直後、velocity2_L/h_p より前に挿入
if (lowMachThornber == 1) {
    flow_float c_hat_th = 0.5*(sonic[ic0] + sonic[ic1]);
    flow_float v2L = Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L;
    flow_float v2R = Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R;
    flow_float z   = min((flow_float)1.0, sqrt(0.5*(v2L+v2R))/c_hat_th);   // M≥1 で z=1 (恒等)
    flow_float um, du;
    um = 0.5*(Ux_L+Ux_R); du = 0.5*(Ux_L-Ux_R); Ux_L = um+z*du; Ux_R = um-z*du;
    um = 0.5*(Uy_L+Uy_R); du = 0.5*(Uy_L-Uy_R); Uy_L = um+z*du; Uy_R = um-z*du;
    um = 0.5*(Uz_L+Uz_R); du = 0.5*(Uz_L-Uz_R); Uz_L = um+z*du; Uz_R = um-z*du;
}
// 以降 velocity2_L/R, h_p/h_m, Vn_p/Vn_m はブレンド後の速度から算出
```

`velocity2_L`/`h_p` は元コードでは L 再構成直後に計算しているが、ブレンドを受けるため
挿入点の**後段へ移動**する (R 側 `velocity2_R`/`h_m` は元から後段)。圧力 `P_{L/R}`・密度
`ro_{L/R}` は不変。`SLAU_d` に `int lowMachThornber` 引数を追加し wrapper の起動に
`cfg.lowMachThornber` を渡す。理論・漸近的根拠は
[`theory.md` の低マッハ再構成補正節](theory.md#低マッハ再構成補正-thornber-型速度ジャンプ縮約) を参照。

> **検証結果 (負)**。`case/23.axi_nozzle` のノズル低マッハ自励振動には**無効〜僅かに悪化**
> (理由・データは theory.md 同節の検証所見および計画 §9 `2026-06-08`)。`lowMachThornber` は
> 正しく opt-in (既定 0・OFF 経路不変) で実装ずみだが、本症状の根治用途では使わない。

#### SST 全エネルギー (`sstEnergyIncludesK`) の p*・H* (2026-09-08)

`CondArgs.energyK != 0` のとき、L/R の面エンタルピーに $\tfrac53 k$ (セル値) を足し、圧力流束 ($\tilde p$ の $P_L, P_R$ と質量流束の $\Delta p$) を
$p^* = p + \tfrac23\rho k$ で評価する。面温度 (TP の $h(T_f)$)・音速・$\chi$/β の Mach スイッチは熱力学圧 $p$ のまま。ghost 側 (ic1 ≥ nCells) の $k$ は
内部値で代用する (node ghostless)。詳細は `methods/turbulence/implementation.md` の同名節。

### `ROE_d` ([L1474](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L1474))

並列単位は `nNormal_halo_Planes`。レジスタ圧が最も高いカーネル。

1. **Ducros 補正リミッタ**: `apply_ducros_limiter(limiter_*[ic0/1], duc)` で
   L/R 別の補正済みリミッタを作る (`duc = clamp(max(ducros[ic0], ducros[ic1]), 0, 1)`)。
2. `interp_dispatch` で `(ro, Ux, Uy, Uz, Ps)` の L/R 状態を構築。
   `ro_*, P_* = max(*, 1e-8)` でクランプ。`roe_*, Ht_*, ca_*` を再計算。
3. **Roe 平均状態** (`sqrt_ro_L/R` から `roa, ua, va, wa, Ha, ca, Ua`)。
4. **保存量差分** `delQ[]` の並びは `(Δρ, Δ(ρE), Δ(ρu), Δ(ρv), Δ(ρw))` で
   一般的教科書順 `(ρ, ρu, ρv, ρw, ρE)` とは異なる:
   ```cpp
   delQ[0] = ro_R - ro_L;                            // Δρ
   delQ[1] = (ro_R*Ht_R - P_R) - (ro_L*Ht_L - P_L);  // Δ(ρE) = Δ(ρe)
   delQ[2] = ro_R*Ux_R - ro_L*Ux_L;                  // Δ(ρu)
   delQ[3] = ro_R*Uy_R - ro_L*Uy_L;
   delQ[4] = ro_R*Uz_R - ro_L*Uz_L;
   ```
5. **圧力の保存変数微分**:
   ```cpp
   P_ro  = 0.5*(γ-1)*(u²+v²+w²);  P_roe = γ-1;
   P_rou = -(γ-1)*u;              P_rov = -(γ-1)*v;   P_row = -(γ-1)*w;
   z   = u²+v²+w² - P_ro/P_roe;   // = 0.5*|u|² (理想気体)
   fai = P_ro - c²;
   ```
6. **右固有ベクトル `R[5][5]`** をその場で構築 (列の順は遅い音波 / せん断 ×3 / 速い音波)。
7. **左固有ベクトル `L[5][5]`** をベタ書きし、最後に二重ループで一括 `L /= c²`
   (実装上のショートカット)。
8. **特性座標** `dW = L * delQ`、**固有値** `lam[] = (|U-c|, |U|, |U|, |U|, |U+c|)`。
9. **Harten 形エントロピーフィックス** (速度比依存型のみ `if (true)`、他はガード `if (false)`):
   ```cpp
   eta_vl = 0.1*(|Ua|/ca + 1.0);
   if (lam[i] <= eta_vl)
       lam[i] = (lam[i]² + eta_vl²) / (2*eta_vl);
   ```
10. **散逸ベクトル** `difQ1 = |lam| * dW`、`difQ2 = R * difQ1`。
11. **中央束** を質量流束 `mdot = 0.5*(ρU_L + ρU_R)*sss` 経由で書き、
    `0.5*(F_L + F_R) - 0.5*difQ2*sss` の形で残差に積算。
12. 両側セルに `atomicAdd(±)`。

### `HLLE_d` ([L1268](../../solver_density_cuda/cuda_forge/convectiveFlux_d.cu#L1268))

並列単位は `nNormal_halo_Planes`。Roe より軽量で堅牢。

1. `interp_dispatch` で L/R 構築 → `max(*, 1e-8)` クランプ。
2. `roe_L/R = P/(γ-1) + 0.5*ρ|u|²`, `Ht_L/R = (roe + P)/ρ`, `U_L/R = u·n̂`,
   `ca_L/R = sqrt(γP/ρ)`。
3. **Roe 平均**で `ua, va, wa, Ha, ca, Ua` を作り、**波速見積もり**:
   ```cpp
   S_L = min(U_L - ca_L, min(U_R - ca_R, Ua - ca));
   S_R = max(U_L + ca_L, max(U_R + ca_R, Ua + ca));
   ```
4. **物理フラックス** `F_L_*, F_R_*` を $(\rho U, \rho u U + P n_x, \dots, \rho H_t U)$ で構築。
5. **三分岐**:
   - `S_L ≥ 0` → `F_L * sss`
   - `S_R ≤ 0` → `F_R * sss`
   - 中間 → `(S_R*F_L - S_L*F_R + S_L*S_R*(Q_R - Q_L)) / max(S_R-S_L, 1e-8) * sss`

   中間状態のエネルギ差は `roe_R - roe_L` (= $\rho e_R - \rho e_L$) を使う。
6. 両側セルに `atomicAdd(±)`。

Ducros 補正は通さず、`limiter_*[ic0/ic1]` を直接 `interp_dispatch` に渡す。

### `KEEP_d` ([convectiveFlux_keep_d.inc.cuh](../../solver_density_cuda/cuda_forge/convection/convectiveFlux_keep_d.inc.cuh))

**純粋 KEEP (Kinetic Energy & Entropy Preserving) 中心流束**。**散逸項を持たない**低散逸対流スキームで、
SGS の散逸は WALE (`turbulence`) が担う構成を想定する。legacy の `KEEP_FVS_d` (中心流束が `if(false)` で
無効・単純中心平均のみ稼働) を modern bundled API (`FaceGeom`/`PrimState`/...) へ移植して KEEP 中心流束を
有効化し、その後 **Roe 行列散逸・MUSCL 再構成・リミタ・Ducros を撤去**して純粋 KEEP に簡素化した。

- **中心流束**: 隣接対 `(ic0,ic1)` の生値で構成。KE は二次分割で保存するが、エントロピーは**厳密保存ではない**:
  算術/交差平均構成のため Tadmor 条件 $\Delta w^{\mathsf T}F^*=\Delta\psi$ を満たさず、エントロピー誤差はジャンプ 3 次
  $O(\Delta^3)$ の**準保存** (滑らかな場では設計次数で整合、大ジャンプでは非保存。
  [tools/verify_keep_tadmor.py](../../solver_density_cuda/tools/verify_keep_tadmor.py) で数値確認済 —
  乱数 O(1) ジャンプで $E=\Delta w^{\mathsf T}F^*-\Delta\psi \ne 0$ (両符号)、微小ジャンプで $|E|=O(\varepsilon^3)$。
  厳密 entropy-conservative には log-mean 平均 (KEPEC 型) が必要)。したがって散逸レイヤ (下記) を足しても
  「流束全体が entropy stable」は主張しない (散逸項単体の ES 性のみ)。
  $\tilde C = \overline{\rho}\,\overline{U}_n S$, 運動量 $\tilde M_i = \tilde C\,\overline{u_i} + \overline{p}\,n_i S$,
  エネルギー $= (\tilde K + \tilde I + \tilde P)S$ ($\tilde K=\tilde C\,\tfrac12\sum u_{i,0}u_{i,1}$,
  $\tilde I=\tilde C\,\tfrac12(p_0/\rho_0+p_1/\rho_1)/(\gamma-1)$, $\tilde P$ は圧力仕事の split)。
- **散逸レイヤ (opt-in・`keepDissType`)**: 既定 (`keepDissType: 0`) は散逸ゼロの純粋 KEEP でビット不変。
  σ=`keepDissCoeff` (既定 0.05)、$c'$ は `lowMachPrecond>=1` で `lowMachCprime` 低マッハスケール。
  - `1` (**scalar ES**): $F \mathrel{-}= \tfrac12\sigma\lambda'\Delta U$, $\lambda'=|U_n|+c'$ を全 5 式に。
    LLF 型は $\Delta w^{\mathsf T}\Delta U\ge0$ で ES。全成分同一 λ' なので渦も食う
    (市松 ~4桁減衰/400step・TGV KE cost 2.7%@σ=0.05)。頑健フォールバック。
  - `2` (**matrix ES, LES 第一候補**): $F \mathrel{-}= \tfrac12\sigma R|\Lambda'|SR^{\mathsf T}\Delta w$
    (entropy-scaled Roe 型)。$w=\partial\eta/\partial U$、$H=RSR^{\mathsf T}$ の $S$ は音響 $\rho/2\gamma$・
    エントロピー $\rho(\gamma-1)/\gamma$・せん断 $p$ ([tools/verify_entropy_scaling.py](../../solver_density_cuda/tools/verify_entropy_scaling.py)
    で数値検証済)。$|\Lambda'|$ は**音響のみ** $|U_n|+c'$、せん断/エントロピーは $|U_n|$ →
    **市松減衰は scalar 同等以上 (7.1e-8/400step) で KE cost 半分 (1.36%)**。
  - **音響波速 $c'$ は `keepDissCprime` (既定 1) で制御し `lowMachPrecond` から独立** (グローバル旗への
    相乗りを廃止)。$c'$ は散逸係数のみに入り伝播は不変。フル $c$ (`keepDissCprime: 0`) は市松掃除最速
    (9.2e-9) だが、渦の $\Delta U_n$ が音響固有ベクトル経由で食われ **KE cost 3.2 倍 (4.35%)** — c' が
    この漏れを塞ぐ (Guillard-Viozat/Thornber の $c\cdot\Delta u$ 汚染の実測)。
  σ は L1 (市松減衰) と L2 (TGV KE) の両ゲートで較正
  ([convection-keep-es-dissipation](../../plans/accepted/convection-keep-es-dissipation.md) 変更ログ参照)。
  `massflux[ip]` は散逸込み総質量流束 (スカラー輸送と整合)。
  - **`keepDissJump` (matrix CPG 枝, 既定 0=生ジャンプ・ビット不変)**: `1` で Δw を**再構成後
    ジャンプ** ($\Delta_{rec}=\Delta_{raw}-\tfrac12(g_i+g_j)\cdot d_{ij}$, κ=0 線形・リミタ無し) から
    組む。純 2Δ モード (市松) は中心勾配が厳密ゼロで再構成が効かずフル減衰のまま、滑らかな場の
    ドレインは桁減 (実測: 市松 1.00 / TGV 層流期 0.006 / 遷移期 0.16 / ピーク 0.45)。
    `2` (**推奨**) はさらに再構成/生の両ジャンプを特性空間へ射影し成分ごと minmod
    (**sign-property クリップ, TeCNO 型**): 各波で $rd_{used}\cdot rd_{raw}\ge0$ が構造保証され
    **エントロピー散逸性が証明付きで成立**。KE コストは jump=1 と生の中間 (64³ TGV 終値 −1.4%)。
    ghost 面 (`ip>=nNormalPlanes`) と再構成正値性 NG の面は生ジャンプへフォールバック。
    設計判断と検証: [convection-keep-diss-recon-jump](../../plans/accepted/convection-keep-diss-recon-jump.md)。
  - **`keepDissPrecond` (matrix CPG 枝, 既定 0=従来・ビット不変)**: `1` で音響対の散逸を
    **Weiss–Smith/Turkel 前処理散逸**に置き換える (**低マッハ壁付き LES の推奨設定**)。
    $c'$ (固有値縮小のみ = 圧力ジャンプ散逸も一緒に弱る) と異なり、前処理は音響 2×2
    サブシステム $D_p=\Gamma|\Gamma^{-1}A_2|$ (閉形 $|M_2|=\varphi_1M_2+\varphi_2 I$,
    固有ベクトル不要) の**行列構造**を変え、**速度ジャンプ散逸を $\rho\,U_r\,\Delta U_n$ へ縮小**
    (低マッハ精度, Guillard–Viozat) しつつ**圧力ジャンプ散逸を $\Delta p/U_r$ へ増強**
    (連続式への圧力 Laplacian = Rhie–Chow 相当の市松キラー)。$M\ge1$ で標準 Roe $|A_2|$ に厳密復帰、
    sym(K) は全域正定値 = **ES 性維持** ([tools/verify_precond_dissipation.py](../../solver_density_cuda/tools/verify_precond_dissipation.py)
    全 6 検証 PASS・完全 Weiss–Smith 5×5 と 6e-4 一致)。実測: **真の市松 (case/35 L1) を
    c' 版の 142 分の 1 へ減衰** (node 4.3桁/400step) しつつ物理コストは c' 版以下
    (backstep Uy rms 6.94 vs 6.89)。有効時 `keepDissCprime` は不使用。せん断/エントロピー波は
    従来どおり $|U_n|$。`keepDissJump` と直交して合成可。なお backstep 型の「ギリギリ解像の
    2-4Δ 物理」による見かけの縞は散逸では消えない (メッシュ解像の仕事)。
    設計・検証: [convection-keep-diss-lowmach-precond](../../plans/active/convection-keep-diss-lowmach-precond.md)。
  - **`keepDissFdBlend` (既定 0=ビット不変, DES 専用)**: `1` で σ を DDES シールド関数 f̃_d で
    face ごとに駆動する: $\sigma_f = \max(\sigma_\mathrm{min},\,(1-\tilde f_{d,f})\,\sigma_\mathrm{max})$
    ($\tilde f_{d,f}$=隣接セル平均、$\sigma_\mathrm{min}$=`keepDissCoeff`、$\sigma_\mathrm{max}$=`keepDissCoeffMax`
    既定 1.0)。RANS 帯 (f̃_d≈0=付着 BL・薄セル) はフル ES 散逸 (他コードの風上相当)、LES 域はフロア
    のみ、という DES 標準の散逸ゾーニング (Travin 型 / SU2 `FD`=max(0.05,1−f_d) / UZEN-Probst の
    KE 保存 flux+f̃_d ブレンドと同型)。`DESmode>0` (model: sst-ddes/sst-iddes) 必須。`fd_shield` は
    `turbulent_viscosity_d` が convectiveFlux より先に毎サブ反復更新するため鮮度整合。
    **`keepDissPrecond=1` との併用 = 演算子ブレンド (opBlend)**: 音響枝 z1/z5 を
    $D_\mathrm{eff}=\sigma_\mathrm{min}D_\mathrm{precond}+(\sigma_f-\sigma_\mathrm{min})D_\mathrm{Roe}^{(c)}$
    と内分する。precond の Δp∝c²/Ur 増幅は較正済み σ_min 固定 (σ_f→1 との積は壁 Ur フロアで
    即発散するため)、RANS/grey 帯の増分は**フル c** の標準 Roe 型 (渦保護不要な帯で音響スケール
    減衰を確保)。precond を外して c' 退避した構成は、乱流発達後に**大波長音響/圧力モードの成長
    不安定** (壁 ΔP が動圧 70 倍まで指数成長; 滑らかモードは jump2 で O(h³)・音響減衰 c'~0.15c
    のためほぼ無散逸の閉箱共鳴) を許すことが case/39 本番で実証されており、fdblend 使用時も
    precond は外さないこと。汚染場減衰試験 (case/39 `run_diag_pdamp_*`): 現行(c'退避)=汚染維持 /
    フル c 一律=22.7 kPa / opBlend=**13.5 kPa 最速**。
  - **`keepDissCbCoeff` / `keepDissCbEps` (既定 0.0 / 0.10, opt-in)**: 高周波圧力欠陥駆動の
    mass-flux 補正 (Rhie–Chow 型市松キラー)。再構成不整合 $\delta p_f^{HF}=p_R^f-p_L^f$
    (keepDissJump≥1 の面再構成を再利用; 線形圧力場で厳密ゼロ、2Δ 市松で生 Δp) を入力に

    $$\delta\dot m_f=-\tfrac12 C_{cb} S_f\,\frac{\delta p_f^{HF}}{U_{r,cb}},\quad
      U_{r,cb}=\min(c,\max(|\boldsymbol u|,\epsilon_{cb}c)),\qquad
      \delta\boldsymbol F_f=\delta\dot m_f\,[1,\bar u,\bar v,\bar w,\bar H_t]^T$$

    を全保存量に整合配布する (スカラー/種/乱流は散逸込み `massflux` 経由で自動整合)。
    **sign gate**: rank-1 構造を利用し、entropy 変数ジャンプとの内積
    $r_f=\Delta\boldsymbol w^T[1,\bar u,\bar v,\bar w,\bar H_t]^T$ に対し
    $\delta p_f^{HF} r_f>0$ の面のみ有効化 → 面ごとに entropy 散逸性を構造的に保証
    (純圧力市松は素通し)。有効面は「内部面かつ両側再構成成功」に限定し、無効面は補正ゼロ
    (生 Δp フォールバック禁止)。`keepDissCbCoeff>0 && keepDissJump<1` は config エラー。
    σ_min×precond の Δp 項 (実効 σ_min/precondEps≈0.333) への**加算**であり、総圧力結合は
    0.333+C_cb/ε_cb。運動量の圧力力 p·n·S には作用しない。既存 σ/f̃_d/opBlend と直交。
    **適用範囲の実測**: 純 2Δ 市松には劇的 (case/35 L1: cell C_cb=0.02 で 5 step 700 倍減衰=基準の
    10 倍速) だが、**2-3Δ の準解像成分には原理的に弱い** (中心勾配が捉えるため δp^HF が生 Δp の
    ~6% に落ちる — これは設計どおりの「解像スケール保護」でもある)。explicit RK3 CFL0.5 の
    安定限界は実効係数 C_cb/ε≲0.2、implicit dual-time (case/39 本番構成) でも C_cb/ε≥1 は発散。
    導入動機だった case/39「丘頂鋸歯」は後に抽出アーチファクトと判明し撤回 (plan 変更ログ)。
    設計・検証: [convection-keep-cb-pressure-correction](../../plans/active/convection-keep-cb-pressure-correction.md)。
  - **`keepDissOpBlendRaw` (既定 0, opBlend 併用時のみ)**: `1` で opBlend の標準 Roe 側 (RANS/grey
    帯の増分) に **jump2 フィルタ非適用の生ジャンプ** (rr1/rr5) を使う。jump2 の再構成フィルタは
    全散逸経路に掛かるため、fdblend で σ_f→1 にしても実効ジャンプが濾され散逸ゾーニングが無効化
    される — 本キーで「RANS 帯=完全風上相当」(DES 定石) を忠実化する。生ジャンプは一様場で厳密 0
    → free-stream 安全。`keepDissFdBlend:1 + keepDissPrecond:1 + keepDissJump:2` 必須 (config 検証)。
    case/39 本番構成で 2000 step 安定を確認済み (run_0019)。導入動機だった「丘頂鋸歯」は
    抽出アーチファクトと判明し撤回されたため採否は保留 (ゾーニング忠実化としての理屈は残る)。
    設計: [turbulence-iddes-sst §4.8](../../plans/active/turbulence-iddes-sst.md)、
    実務相場: [turbulence-des-flux-survey §7](../../notes/investigations/turbulence-des-flux-survey.md)。
- **free-stream 保存 (`space.pRef`)**: 運動量圧力項 Gtilde は `(Ps−pRef)` の面ごと差分で組む
  (SLAU と同処方)。非直交メッシュの float32 桁落ちによる偽運動量源を除去し、
  `case/33` 歪み hex で機械精度 (4e-12) の一様場保存。**非直交メッシュで KEEP を使うときは
  `pRef` に動作静圧を必ず設定** (無指定 0.0 は従来挙動=非直交で発散)。
- **動く一様流の保存 (`space.roRef` + `space.uRef`, advGauge)**: pRef は静止一様場のみ守る。
  U∞≠0 では移流項 $(\rho u,\ \rho u u,\ \rho u H)\cdot r_c$ の桁落ちが効き (エネルギー
  $\rho u H\sim10^7$ が支配的)、M0.1 でも KEEP/SLAU とも数十 step で発散する
  (`case/33` run_0008-0010)。対策として参照一様流 $(\rho_\infty,\mathbf u_\infty,p_\infty{=}\texttt{pRef})$
  の厳密流束 $F_\infty(s)$ (s 線形の定数流束: 保存厳密・セル和 $F_\infty(\Sigma s)=0$ の解析ゲージ)
  を面ごとに引く。要点は**差分形因数分解** ((差分)×(平均)+(参照)×(差分)) で大きな積を作る前に
  引くこと (式と検証は [plan §8](../../plans/active/convection-freestream-preserving-flux.md)、
  数値検証 [tools/verify_advective_gauge.py](../../solver_density_cuda/tools/verify_advective_gauge.py))。
  M0.1 動く一様流で step0 残差 1.3e-2 → 5.6e-9 (機械精度・6.5桁改善)・300step 安定 (run_0011)。
  `uRef` は平均流/フリーストリーム速度。**スコープ: KEEP CPG × cell/node** (node は境界半割面
  `convectiveFlux_boundary_d` にも upwind 差分形ゲージを実装し CV 全面で telescoping を維持
  ([plan §8.5](../../plans/active/convection-freestream-preserving-flux.md)、`case/33` run_0017-0020 で
  cell と同水準の機械精度を検証)。TP は $e_\infty$ の thermo 評価が未対応で off。
  スカラー輸送 `massflux` は物理流束のまま=種の free-stream 桁落ちは残課題)。
  既定 `roRef: 0.0` = off (ビット不変)。
- **cell/node 両対応**: 周回面は `geom.nLoopPlanes` (= `convPlaneBound`)。cell は内部+境界 ghost を
  周回 (境界面は ic1=ghost の生値で中心流束、専用境界カーネルは skip)、node 弱形式は内部双対面のみ
  周回し境界は `convectiveFlux_boundary_d` が担う。`massflux[ip]` に総質量流束を書きスカラー輸送と整合。
- **保存性 (Taylor-Green M0.4, $32^3$)**: cell・node とも運動量 ~1e-7・KE 0.4%・エントロピー ~1e-5 で保存
  (非粘性)。検証 [`case/09.Taylor-Green`](../../case/09.Taylor-Green/README.md)。なお cell 全周期の保存には
  partnerCellID の device 転送修正が前提 ([boundary-cell-periodic-conservation](../../plans/accepted/boundary-cell-periodic-conservation.md))。
- **注意 (低マッハ checkerboard)**: 純粋 KEEP (`keepDissType: 0`) は非散逸ゆえ低マッハ圧力 odd-even
  (2Δ 市松) を**原理的に減衰できない** (中心平均が市松モードを相殺して見えない=離散 null-mode)。
  **`lowMachPrecond` は SLAU 散逸内の c' 置換であり、散逸を持たない純 KEEP には作用しない** (前処理対象が無い)。
  市松が出る場合は `keepDissType: 1` (低マッハスケール scalar ES 散逸) を使う。検証は
  [case/35](../../case/35.uniform_periodic_box/README.md) の市松摂動ケース (L1) 参照。
- **periodic 継ぎ目の P 振動 (node, 解消済)**: node periodic 合併 DOF の保存量 state が master/slave で desync
  していた件は **root→member ミラー** (§4.5.9 / `periodicMirrorNSState`) で解消済。詳細は
  [median-dual-3d plan §4.5.9](../../plans/active/discretization-median-dual-3d.md)。

legacy の `KEEP_SLAU_d` / `AUSMp_d` / `AUSMp_UP_d` は依然 `legacy/` に退避され到達不能。

## 並列化メモ

- 並列単位は面 1。`normal_halo_planes_d` で実際に処理する面のみ列挙し、
  `dimGrid_normal_halo = ceil(nNormal_halo_Planes / blocksize)` をグリッドサイズとする。
- 面ごとに両側セルへ `atomicAdd`。`flow_float = double` の場合 CC 6.0+ が必要。
- 面ループ内で `R[5][5], L[5][5]` などローカル配列を構築するためレジスタ圧が高い。
  Roe カーネルは特に重い。

## 入出力

入力: 保存量 / 原始量 (`ro, roUx, roUy, roUz, roe, Ux, Uy, Uz, P, Ht, sonic`)、
勾配 (`drod*`, `dUxd*`, `dUyd*`, `dUzd*`, `dPd*`)、リミッタ
(`limiter_ro, limiter_Ux, limiter_Uy, limiter_Uz, limiter_P`)、Ducros 値 (`ducros`)、
セル中心・面中心座標、面ベクトル `(sx, sy, sz, ss)`、補間重み `fx`、`plane_cells`。

出力 (累積): セル中心残差 `res_ro, res_roUx, res_roUy, res_roUz, res_roe`。

## 既知の TODO / 注意点

- KEEP (`KEEP_d`) は復活済み (cell/node 両対応)。`KEEP_SLAU` / `AUSM` 系は依然 legacy で無効。
- `c_hat` の中央平均は単純算術平均。Roe 平均の方が低マッハ性に望ましい場面あり。
- リミッタは `interp_dispatch` 経由で各成分 (`ro, Ux, Uy, Uz, P`) 独立。
  限界保証 (TVD) は厳密には未保証 (MUSCL 2/3 次経路は無リミット)。
- 旧 CPU 経路 (`convectiveFlux.cpp`) は使われていない。GPU と同期させる場合は
  ラッパと dispatch 構造を反映する必要がある。
