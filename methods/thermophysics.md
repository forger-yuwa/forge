# 多成分 thermally-perfect gas

本ドキュメントは理論(係数・方程式)と実装(ソース対応)をまとめる。

## 理論

### 1. 目的と適用範囲

本章は forge に導入する**多成分 thermally-perfect gas (熱的完全気体)** の理論的背景をまとめる。対象は**非反応**の多成分混合 (混合・希釈・多成分境界層など) であり、化学反応源項は扱わない (将来拡張点)。

従来の forge は単一成分・**熱量的完全気体 (calorically perfect gas, CPG)** のみを扱っていた。すなわち比熱 $c_p,\,c_v$ と比熱比 $\gamma$ は定数で、温度は内部エネルギーから線形に陽に求まる。本拡張では:

- 比熱 $c_p(T)$ を温度依存とする (**NASA-9 多項式**, CEA 準拠)。
- 複数化学種の質量分率 $Y_s$ を輸送し、混合物性を **ideal-gas mixing** で評価する。
- 粘性・熱伝導率・質量拡散係数を **kinetic theory** (Chapman-Enskog) で評価する。

`thermalMethod==0` (CPG) は完全に保持し、`thermalMethod==2` で本モデルを有効化する。

### 2. 状態方程式と熱力学

#### 2.1 理想気体の状態方程式

各化学種・混合物とも理想気体とする。混合気体定数は質量分率 $Y_s$ から

$$ R_{\mathrm{mix}} = R_u \sum_s \frac{Y_s}{W_s} = \frac{R_u}{W_{\mathrm{mix}}}, \qquad \frac{1}{W_{\mathrm{mix}}} = \sum_s \frac{Y_s}{W_s} $$

ここで $R_u=8.314462618\ \mathrm{J/(mol\,K)}$、$W_s$ は化学種 $s$ の分子量 [kg/mol]。圧力は

$$ P = \rho R_{\mathrm{mix}} T. $$

#### 2.2 NASA-9 多項式による比熱・エンタルピー・エントロピー

CEA (McBride & Gordon 2002) のネイティブ形式である 9 係数多項式を用いる。各化学種・温度域ごとに係数 $a_0\dots a_8$ をもち、モル基準で

$$ \frac{C_{p,s}}{R_u} = a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T + a_4 T^2 + a_5 T^3 + a_6 T^4 $$

$$ \frac{H_s}{R_u T} = -a_0 T^{-2} + a_1 \frac{\ln T}{T} + a_2 + \frac{a_3}{2} T + \frac{a_4}{3} T^2 + \frac{a_5}{4} T^3 + \frac{a_6}{5} T^4 + \frac{a_7}{T} $$

$$ \frac{S^{\circ}_s}{R_u} = -\frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln T + a_3 T + \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3 + \frac{a_6}{4} T^4 + a_8 $$

積分定数 $a_7,a_8$ により、エンタルピーは**標準生成エンタルピーを含む絶対エンタルピー基準**となる。係数は 2 温度域 ($T_{\mathrm{lo}}\le T<T_{\mathrm{mid}}$ と $T_{\mathrm{mid}}\le T\le T_{\mathrm{hi}}$、標準は 200/1000/6000 K) で切り替える。範囲外は端でクランプし、エンタルピーは $h(T)\approx h(T_c)+c_p(T_c)(T-T_c)$ と線形外挿して衝撃波での暴走を防ぐ。

質量基準は $c_{p,s}=C_{p,s}/W_s$、$h_s=H_s/W_s$。

##### エンタルピー基準オフセット (sensible-enthalpy datum, 非反応流)

絶対基準では生成エンタルピーが種ごとに桁違い ($h_{\mathrm{H_2O}}(298\,\mathrm{K})\approx-13.4\,\mathrm{MJ/kg}$ に対し $h_{\mathrm{N_2}},h_{\mathrm{O_2}}\approx0$)。**化学反応が無い流れでは各種のエンタルピー基準を任意の定数 $c_s$ だけ平行移動してよい**。実際、全エネルギー方程式に現れるのは対流・拡散によるエンタルピー輸送 $\rho h_{\mathrm{mix}}\mathbf{u}+\sum_s h_s \mathbf J_s$ と内部エネルギー貯蔵 $\rho e_{\mathrm{mix}}$ であり、$h_s\to h_s-c_s$ とすると流束発散の変化 $+\partial_t(\sum_s c_s\rho Y_s)$ と貯蔵項の変化 $-\partial_t(\sum_s c_s\rho Y_s)$ が (反応が無く種連続式 $\partial_t(\rho Y_s)+\nabla\!\cdot(\rho Y_s\mathbf u+\mathbf J_s)=0$ が成り立つので) **厳密に相殺**する。基準はエントロピー $s^\circ_s$ には影響しない。

そこで $c_s=h_s(T_{\mathrm{ref}})$ を選び全種を $h_s(T_{\mathrm{ref}})=0$ (sensible enthalpy) に揃えると、桁違いの生成エンタルピーが消えて種間のエンタルピーが同程度の大きさになる。これは **多成分 implicit の数値安定化**に効く: $\rho e$ と $\rho Y_s$ を別緩和で解く分離 implicit では緩和ミスマッチが温度反転に伝わるが、その増幅率が $|h_s|$ に比例するため、生成エンタルピーを除くと温度ジャンプが抑えられる (詳細・検証は implementation.md と `case/16.nozzle_wys`)。反応を扱う場合は生成エンタルピー差が物理的に効くためこの平行移動はできない。

#### 2.3 ideal-gas mixing による混合物性

混合比熱・絶対エンタルピーは質量分率重み:

$$ c_{p,\mathrm{mix}}(T) = \sum_s Y_s\,c_{p,s}(T), \qquad h_{\mathrm{mix}}(T) = \sum_s Y_s\,h_s(T). $$

定容比熱・比熱比は

$$ c_{v,\mathrm{mix}} = c_{p,\mathrm{mix}} - R_{\mathrm{mix}}, \qquad \gamma_{\mathrm{mix}} = \frac{c_{p,\mathrm{mix}}}{c_{v,\mathrm{mix}}}. $$

比内部エネルギーと凍結音速は

$$ e_{\mathrm{mix}}(T) = h_{\mathrm{mix}}(T) - R_{\mathrm{mix}} T, \qquad a = \sqrt{\gamma_{\mathrm{mix}} R_{\mathrm{mix}} T}. $$

#### 2.4 温度反転 (Newton 法)

保存変数からは比内部エネルギー $e = \rho e_t/\rho - \tfrac12|\mathbf{u}|^2$ が得られる。CPG では $T=e/c_v$ と陽に解けるが、TP では $e_{\mathrm{mix}}(T)$ が非線形なので Newton 反転で解く:

$$ f(T) = e_{\mathrm{mix}}(T) - e = 0, \qquad f'(T) = c_{v,\mathrm{mix}}(T) = c_{p,\mathrm{mix}}(T) - R_{\mathrm{mix}}. $$

$$ T^{(k+1)} = T^{(k)} - \frac{e_{\mathrm{mix}}(T^{(k)}) - e}{c_{v,\mathrm{mix}}(T^{(k)})}. $$

厳密微分 $c_{v,\mathrm{mix}}$ により 2 次収束する。初期値は前ステップの $T$ をウォームスタートに使い (定常で約 2 反復)、各反復で $T\in[T_{\min},T_{\max}]$ にクランプする。反転の精度は次のとおり (2026-09-12 改訂, plan performance-3d-node-sst-speedup §4.2-3): 既定 `physProp.thermoFloat: 1` では
**ハイブリッド** = float 係数 (`SpeciesThermoF`) の Newton (warm start, 最大 12 反復) で ~1e-6·T まで寄せ、double の Newton を
収束 (|ΔT| < 1e-3 + 1e-6·T) まで最大 3 段当てて研磨する (`thermo_T_from_e_hybrid`; 通常 1 段)。厳密 double 参照との誤差は ≤4e-10·T
(`tools/test_thermo_float.cpp`, 冷間開始・組成端点・50–6000 K 込み)。`thermoHrefTemp>0` が前提で、絶対 datum (H2O の h≈−13 MJ/kg)
では float 段が収束しないので datum 無し config では自動的に従来の **全 double Newton** (`thermoFloat: 0` 相当) に落ちる。
面ごとの h_mix(Y_f, T_f) (SLAU) と化学種拡散の h_s(T_f) も float 係数で評価する (評価点は不変)。

### 3. 化学種輸送方程式

非反応・質量分率 $Y_s$ の保存形 ($\rho Y_s$):

$$ \frac{\partial (\rho Y_s)}{\partial t} + \nabla\!\cdot(\rho \mathbf{u} Y_s) = -\nabla\!\cdot \mathbf{J}_s, \qquad s=1,\dots,N $$

$\mathbf{J}_s$ は化学種拡散流束 (§5)。移流は**連続の式と同じ質量流束** $\dot m$ を再利用して上流化し、$\sum_s$ が連続の式と整合する (質量保存)。実現可能性として $Y_s\ge 0$ をフロアし、各ステージで $\sum_s Y_s = 1$ に再正規化する。

#### 3.1 陰解法での緩和整合 (multi-component implicit coupling)

定常陰解法 (block-DPLUR) では流れ 5 変数 $(\rho,\rho u,\rho v,\rho w,\rho e_t)$ を 5×5 ブロックの近似 Jacobian で擬似時間前進する。化学種 $\rho Y_s$ を**別系**として更新すると、両者の**擬似時間緩和速度がミスマッチ**し、状態方程式 $T=T(e,Y)$ の Newton 反転で温度ジャンプを生む。なぜなら $e=\rho e_t/\rho-\tfrac12|\mathbf u|^2=\sum_s Y_s e_s(T)$ であり、$\rho Y_s$ の増分 $\Delta(\rho Y_s)$ と $\rho e_t$ の増分 $\Delta(\rho e_t)$ が**同じ緩和率で進まない**と、各擬似時間ステップで $(e,Y)$ が相互不整合になり、$T$ がそのズレを吸収して跳ねる。増幅率は $|h_s|$ に比例し、生成エンタルピーの大きい $\mathrm{H_2O}$ ($h\approx-13.4$ MJ/kg) で発散を誘発する ($\mathrm{N_2},\mathrm{O_2}$ は無害)。sensible-enthalpy datum (§2.2) は $h_s$ の桁を下げて増幅を**振幅で**抑えるが、緩和ミスマッチという**構造**は残る。

**緩和整合 (matched relaxation)**: 化学種を流れと**同一の擬似時間緩和**で前進させれば、$\Delta(\rho Y_s)$ と $\Delta(\rho e_t)$ が同じ率で進み、各ステップで $(e,Y)$ が整合したまま $T$ が決まる。$\rho Y_s$ は受動移流 (接触波速 $V=\mathbf u\!\cdot\!\mathbf n$ のみで音響波に結合しない) ので、流れの 5×5 ブロックに同梱せずとも、**スカラ DPLUR** で同一の凍結残差・同一局所擬似時間 $\Delta\tau$・同一緩和係数 $\omega$・同一 sweep 回数で緩和すれば緩和率が一致する。1 次風上移流の凍結 Jacobian は

$$ \Big(\tfrac{V}{\Delta\tau}+\!\sum_f \tfrac{\max(\dot m_f,0)}{\rho}\Big)\,\delta(\rho Y_s)_{\text{(自セル)}}\;-\;\sum_f \tfrac{\max(-\dot m_f,0)}{\rho_{\text{nbr}}}\,\delta(\rho Y_s)_{\text{(隣)}}=\mathcal{R}_{\rho Y_s}, $$

すなわち対角=流出質量流束 (連続式と整合)、非対角=流入質量流束で、流れの DPLUR と同じ M 行列構造を持つ。Jacobi sweep を $n_{\text{StepInner}}$ 回反復し $\omega=\texttt{implicitRelax}$ で緩和、最後に $\rho Y_s=\rho Y_s^{\,n}+\delta(\rho Y_s)$ を commit して $\sum_s\rho Y_s=\rho$ へ再正規化する。$n_{\text{StepInner}}=1,\omega=1,$ 非対角無視で従来の点陰的 (segregated) 更新に一致する。

化学種拡散 ($\mathbf J_s$, §5) の非対角は省いて点陰的のまま残す (拡散は剛性が低く、定常 ($\mathcal R\to0$) では $\delta\to0$ ゆえ収束先は不変)。完全結合 (5+$N$ ブロック) は最も根本的だが、まず緩和整合で「緩和率の統一だけで安定 $\Delta\tau$ 上限が上がるか」を切り分ける。実装は 本ドキュメントの「実装」節。

### 4. kinetic theory による輸送係数

#### 4.1 純成分 (Chapman-Enskog)

換算温度 $T^*=k_B T/\varepsilon_s$ と Lennard-Jones パラメータ ($\sigma_s$ [Å], $\varepsilon_s/k_B$ [K]) を用い、衝突積分 $\Omega^{(2,2)*},\Omega^{(1,1)*}$ は Neufeld ら (1972) の曲線近似で評価する。純成分粘性・熱伝導率 (modified Eucken) は

$$ \mu_s = 2.6693\times10^{-6}\,\frac{\sqrt{W_s[\mathrm{g/mol}]\,T}}{\sigma_s^2\,\Omega^{(2,2)*}(T^*)}\ [\mathrm{Pa\,s}], \qquad k_s = \mu_s\,(1.32\,c_{v,s} + 1.77\,R_s). $$

二成分拡散係数は

$$ D_{ij} = 1.8583\times10^{-7}\,\frac{\sqrt{T^3\,(1/W_i+1/W_j)}}{P\,\sigma_{ij}^2\,\Omega^{(1,1)*}(T^*_{ij})}, \qquad \sigma_{ij}=\tfrac12(\sigma_i+\sigma_j),\ \varepsilon_{ij}=\sqrt{\varepsilon_i\varepsilon_j}. $$

#### 4.2 混合則

粘性は Wilke、熱伝導率は Wassiljewa/Mason-Saxena:

$$ \mu_{\mathrm{mix}} = \sum_i \frac{X_i \mu_i}{\sum_j X_j \phi_{ij}}, \quad k_{\mathrm{mix}} = \sum_i \frac{X_i k_i}{\sum_j X_j \phi_{ij}}, \quad \phi_{ij} = \frac{[1+\sqrt{\mu_i/\mu_j}\,(W_j/W_i)^{1/4}]^2}{\sqrt{8(1+W_i/W_j)}}. $$

ここで $X_i$ はモル分率。質量拡散は混合平均 (Curtiss-Hirschfelder):

$$ D_{i,\mathrm{mix}} = \frac{1-X_i}{\sum_{j\ne i} X_j/D_{ij}}. $$

### 5. 化学種拡散とエネルギー結合

化学種拡散流束は $\mathbf{J}_i = -\rho D_{i,\mathrm{mix}} \nabla Y_i$。混合平均は $\sum_i \mathbf{J}_i\ne 0$ となるため、補正速度で $\sum_i \mathbf{J}_i=0$ を担保する (定数 Schmidt 数フォールバックも用意)。エネルギー方程式には**化学種拡散によるエンタルピー輸送** $\sum_i h_i \mathbf{J}_i$ を熱伝導 $-k\nabla T$ に加える:

$$ \mathbf{q} = -k\,\nabla T + \sum_i h_i(T)\,\mathbf{J}_i. $$

### 6. 数値スキームとの整合 (対流流束)

密度ベース流束 (SLAU, Roe 等) は CPG では被移流全エンタルピーを $\gamma/(\gamma-1)\,P/\rho + \tfrac12|\mathbf u|^2$ と再構成するが、これは $c_p T$ に等しく、TP の絶対エンタルピー $h(T)=\int c_p\,dT'$ とは**一致しない**。TP では面の全エンタルピーを NASA の $h(T_{\mathrm{face}})$ で再構成する (全エネルギー保存を決定する最重要項)。音速は混合整合な $a=\sqrt{\gamma_{\mathrm{mix}}R_{\mathrm{mix}}T}$ を用いる。詳細は [`methods/convection/`](convection) および実装解説を参照。

### 参考文献

- B. J. McBride, M. J. Gordon, *NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species*, NASA/TP-2002-211556 (2002).
- R. A. Svehla, *Estimated Viscosities and Thermal Conductivities of Gases at High Temperatures*, NASA TR R-132 (1962).
- P. D. Neufeld, A. R. Janzen, R. A. Aziz, *J. Chem. Phys.* 57, 1100 (1972).
- C. R. Wilke, *J. Chem. Phys.* 18, 517 (1950).
- M. Vinokur, J.-P. Montagné, *J. Comput. Phys.* 89, 276 (1990) — 一般 EOS の Roe 平均。

## 実装

### 1. 熱力学モジュール `cuda_forge/thermo_d.{cuh,cu}`

#### データ構造
`SpeciesThermo` (POD): 分子量 `MW` [kg/mol]、Lennard-Jones パラメータ `sigma_LJ` [Å]・`eps_kB` [K]、温度域 `Tlo/Tmid/Thi`、NASA-9 係数 `low[9]/high[9]`。

#### device/host 共通関数 (`__host__ __device__`, 内部 double)
- `thermo_cp_molar / thermo_h_molar` — NASA-9 評価。範囲クランプ + エンタルピー線形外挿。
- `thermo_cp_mass / thermo_h_mass / thermo_R_species` — 質量基準・比気体定数。
- 混合則 `thermo_inv_Mmix / thermo_R_mix / thermo_cp_mix / thermo_h_mix / thermo_e_mix / thermo_gamma_mix` — 質量分率 `Y[]` を引数に取る (`__host__ __device__` 共通なので初期条件 host 計算とカーネルで同一)。
- `thermo_T_from_e` — 比内部エネルギーからの Newton 温度反転 (厳密微分 $c_v$、毎反復クランプ、最大 20 反復)。
- `thermo_T_from_h` — 比エンタルピーからの反転 (Roe 平均状態の有効 γ 用)。

#### DB 構築・アップロード (`thermo_d.cu`)
- 内蔵 DB に代表化学種 (N2/O2/Ar/CO2/He/AIR) の NASA-9 + LJ を保持。
- `thermo_init_db(cfg)`: `cfg.speciesNames` 順に host 配列を構築し、`cfg.speciesDBFile` (yaml) があれば上書き、device global memory (`cudaMalloc`) へアップロード。
- アクセサ `thermo_species_device_ptr()` / `thermo_num_species()` / `thermo_species_host()`。化学種データは translation unit をまたぐため `__constant__` ではなく device pointer をカーネル引数で渡す (`-rdc` 非依存)。
- `main.cpp` の `initializeSimulation` で `cfg.read()` 直後に `thermo_init_db(cfg)` を呼ぶ。
- **⚠ TP + 陰解法 (block-DPLUR) では `thermoHrefTemp` を必ず指定すること (2026-08-17 確定)**:
  絶対基準 (既定 0) のままだと陰解法 Jacobian の $\chi_{\rm eos}=c^2-\kappa h$ に生成エンタルピー込みの
  $h$ (燃焼ガスで $-1.8\times10^6$ J/kg) が入り桁違いになって、軸対称ノズルでは軸近傍から
  発散する (case/42 run_0020–0025 で切り分け: 一定 $c_p$ 種・陽解法は完走、実 NASA-9 + 陰解法だけ
  発散、`thermoHrefTemp: 298.15` で完走。低温外挿・軸ホップ Jacobian は無罪)。IC の `roe` も同じ
  datum で作ること (設計側 `paste_isentropic_ic(h_ref_T=...)`)。
- **エンタルピー基準オフセット** (`physProp.thermoHrefTemp`, >0 で有効・既定 0=従来絶対基準でビット不変): host 配列構築後に、各種の NASA-9 係数 `a[7]` (`low[7]`,`high[7]` の両方) に $\Delta a_7=-h_s(T_{\mathrm{ref}})/R_u$ を加算する。NASA-9 では $h_{\mathrm{molar}}=R_u T(\dots)+R_u a_7$ と $a_7$ が定数項 $R_u a_7$ を与えるため、これで全 $T$ で一定オフセット $c_s=-h_s(T_{\mathrm{ref}})$ を与え $T_{\mathrm{mid}}$ 連続性も保たれる ($h_s(T_{\mathrm{ref}})=0$, sensible enthalpy)。係数に焼き込むので `thermo_h/cph/T_from_e/混合則/SLAU 面エンタルピー/BC` が自動的に同一基準になり**一点改変で全経路整合**、エントロピー $a_8$ は不変。非反応流で物理不変・多成分 implicit の安定化目的 (theory.md「sensible-enthalpy datum」参照)。IC 生成器 (`case/16.nozzle_wys/gen_ic_2sp_from_n2.py` 等) も同一 $T_{\mathrm{ref}}$ でオフセットして整合させること (`roe` 不整合は step0 ジャンプを招く)。検証: `case/16.nozzle_wys` run_0101〜0109 で安定 `cfl_pseudo` 上限 1→2、step0 再構成が絶対基準と機械精度一致 (物理不変)。

#### 内蔵 species DB の一覧と出典

内蔵 species と各定数の一次情報は `builtinDB()` ([`thermo_d.cu:53`](../solver_density_cuda/cuda_forge/thermo_d.cu#L53)) にベタ書きされている。各値の出典は下表の通り (NASA-9 係数 / Lennard-Jones パラメータで出典が異なる)。外部 yaml (`physProp.speciesDBFile`, [`thermo_d.cu:142`](../solver_density_cuda/cuda_forge/thermo_d.cu#L142)) で上書き/追加できる。温度域は全種共通で $T_{\mathrm{lo}}/T_{\mathrm{mid}}/T_{\mathrm{hi}}=200/1000/6000$ K。

| species (別名) | MW [kg/mol] | $\sigma_{LJ}$ [Å] | $\varepsilon/k_B$ [K] | NASA-9 出典 | LJ 出典 | 備考 |
| --- | --- | --- | --- | --- | --- | --- |
| N2 | 0.0280134 | 3.621 | 97.53 | CEA [1] | [2] 系 | — |
| O2 | 0.0319988 | 3.458 | 107.4 | CEA [1] | [2] 系 | — |
| Ar (`AR`/`Ar`) | 0.039948 | 3.330 | 136.5 | 単原子理想 [†] | [2] 系 | NASA-9 は $a_2=c_p/R=5/2$ のみ・$a_7,a_8$ は CEA 基準定数 |
| CO2 | 0.0440095 | 3.763 | 244.0 | CEA [1] | [2] 系 | — |
| He (`HE`/`He`) | 0.0040026 | 2.551 | 10.22 | 単原子理想 [†] | [2] | 単原子 (Ar と同形)。$\varepsilon/k_B$ は Svehla [2] に一致 |
| H2O (`h2o`/`WATER`) | 0.0180153 | 2.605 | 572.4 | CEA [1] | [2] 系 | 非平衡凝縮 (Wyslouzil) のキャリア + 凝縮種 |
| AIR (`Air`/`air`) | 0.0289647 | 3.711 | 78.6 | **forge 擬似種 [‡]** | 空気の標準 LJ [2] 系 | 単成分擬似空気。$c_p/R=7/2$ 一定 ($\gamma=1.4$)、CEA 種ではない |

出典記号:

- **[1] CEA**: B. J. McBride, M. J. Gordon, *NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species*, NASA/TP-2002-211556 (2002)。NASA Glenn CEA データベース (`thermo.inp`) の 2 区間 (200–1000–6000 K) 係数。N2/O2/CO2/H2O は published CEA 値そのまま。
- **[†] 単原子理想**: Ar/He は理想単原子気体として $c_p/R=5/2$ を厳密値で与える ($a_0{=}a_1{=}a_3{\dots}{=}0,\ a_2{=}2.5$)。エンタルピー基準定数 $a_7=-745.375$ (生成エンタルピー 0 の単原子基準)・エントロピー定数 $a_8$ は CEA に一致。
- **[‡] forge 擬似種 (AIR)**: CEA の Air 混合に近い熱力学を**単成分**で近似するための簡便種。$c_p/R=3.5$ 一定の擬似多項式 ($a_2{=}3.5$, $a_7=-1043.1$, $a_8=3.0$) で、低温で $c_p\approx1004.5$ J/(kg·K)・$\gamma=1.4$ となる。**CEA 由来ではない** ので、空気を厳密に扱う場合は N2/O2 混合か外部 DB を使う。
- **[2] LJ パラメータ**: R. A. Svehla, *Estimated Viscosities and Thermal Conductivities of Gases at High Temperatures*, NASA TR R-132 (1962) ほか標準の動力学理論輸送データ集に基づく。各 $\sigma_{LJ}$/$\varepsilon/k_B$ は広く使われる輸送データ集 (例: CHEMKIN transport database) とも整合する。He の $\varepsilon/k_B$ は Svehla 値に一致。「[2] 系」は文献単位での整合を意味し、**species ごとの該当表・行レベルの個別照合までは未追跡** — 厳密な監査が必要な場合は一次文献の該当表と突き合わせること。

> 衝突積分の近似式 ($\Omega^{(2,2)*},\Omega^{(1,1)*}$) は Neufeld et al. (1972) 閉形式 (§5c)。理論的背景と全参考文献は 本ドキュメントの「理論」節 を参照。

### 2. 従属変数と温度反転 `cuda_forge/dependentVariables_d.cu`

`dependentVariables_d` に `thermalMethod`・化学種ポインタ・`gam_array/cp_array` を追加。`thermalMethod==2` では:
- 組成 `Y` を構築 (`nSpecies==1` は `Y={1}`、多成分は `roY_s/ρ` をフロア+正規化)。
- `thermo_T_from_e` で温度反転 (ウォームスタートに前ステップ `T[ic]`)。
- `P=ρ R_mix T`、`Ht=h_mix+ek`、`sonic=√(γ_mix R_mix T)`、`gamma[ic]=γ_mix`、`cp[ic]=cp_mix` を更新し、`roe` を (floor 済 ρ, 反転 T) と整合再構成。

`thermalMethod==0` は従来式 (`T=e/(cp/γ)`, `P=(γ-1)(ρe-ρek)`) をビット不変で保持。CPU フォールバック (`dependentVariables.cpp`, `gpu==0`) は CPG のみ対応 (GPU 経路が主)。

物性更新 `gasProperties_d.cu` は `thermalMethod==0` のときだけ `gamma/cp` を定数で埋める。`thermalMethod==2` では `dependentVariables` が毎ステップ `gamma/cp` を設定するため触らない (粘性は viscMethod に従う。kinetic theory は M3)。

#### EOS 正値化フロア (config 化)

`dependentVariables_d` は膨張領域での非物理的な $\rho\to0,\,P\to0$ による速度爆発を防ぐため、密度・圧力・温度に下限フロアを適用する (全 EOS 経路: 単相 CPG / 二相 CPG / TP gas)。フロア値は `physProp.{pMin,roMin,tMin}` で指定でき、既定は従来ハードコード値 (`pMin=1.0` Pa, `roMin=1e-4` kg/m³, `tMin=1e-4` K) と同一のため**未指定なら従来挙動 (ビット不変)**。

- 注意: 既定 `pMin=1.0` Pa は大気スケール想定であり、**無次元・低圧ケースでは場を破壊する**。例えば Taylor-Green は $P_0=1/\gamma\approx0.714$ Pa で圧力場全体が $1.0$ Pa 未満のため、既定フロアが初期場を $P\equiv1.0$ にクランプし保存量も崩れる。ケースの圧力/密度スケールに合わせ `pMin` を下げる (例 `1e-6`)。
- フロア適用後、単相 CPG 経路は `roe = P_temp/(γ-1) + ρ·ek` で保存量 `roe` をフロア後の `P,ρ` と整合再構成するため、フロアが発火するセルでは `roe` (=全エネルギー) が改変される点に注意。
- 設計判断は [`plans/accepted/thermophysics-eos-positivity-floor-config.md`](../plans/accepted/thermophysics-eos-positivity-floor-config.md)。

### 3. 対流流束の TP 整合 `cuda_forge/convectiveFlux_d.cu`

`SLAU_d` / `ROE_d` に `int thermalMethod, const SpeciesThermo* sp, int nSpecies` を追加。ラッパーから `cfg.thermalMethod, thermo_species_device_ptr(), cfg.nSpecies` を渡す。

#### SLAU (最小修正)
被移流全エンタルピー `h_p/h_m` のみ TP 化:面温度 `T_face=P_face/(ρ_face R)` → `h(T_face)`(NASA) + `½|u|²`。音速 `c_hat=½(sonic[i]+sonic[j])` は既に `sonic[]` 配列 (混合整合) のため無改造。圧力束は `P` 直接で EOS 非依存。

**多成分 (M5)**: `h_p/h_m` を L 側 ic0 / R 側 ic1 の**セル組成 `Y`** で再構成 (`thermo_R_mix`/`thermo_h_mix`)。`sp[0]` 固定だと He/N2 のように質量基準 cp が大きく異なる contact でエネルギー束が不整合になり圧力発散する。`SLAU_d` に `flow_float** roY` (= `species_roY_device_ptr()`) を追加。`nSpecies<=1` は従来 `sp[0]` 経路でビット不変。組成は 1 次 (セル中心値) で再構成 (高次は将来)。

#### ROE (段階A)
L/R 状態の `roe_L/Ht_L/ca_L` (および R 側) を NASA で再構成。Roe 平均状態は `thermo_T_from_h` で `h̃=Ha-ek` から `T̃` を反転し、有効比熱比 `ga_eff=cp(T̃)/(cp(T̃)-R)` と音速 `ca=√(ga_eff R T̃)` を用いる。圧力微分ブロック `P_ro,P_roe,...` は `κ=ga_eff-1` を採用 ($\chi=\partial P/\partial\rho|_{\rho e}\ne0$ 項は stage-A で無視。defect-correction なので収束解は不変)。厳密化は Vinokur-Montagné (stage-B, 将来)。

> SLAU は多成分対応済 (M5, セル組成再構成)。ROE/HLLE/AUSM/KEEP は単成分 (`sp[0]`) のまま (多成分化は後続)。

### 4. 陰解法 Jacobian (予定)

`timeIntegration_d.cu` の `build_jacobian_split/build_abs_jacobian` は `chi=(γ-1)/sonic` で EOS をパラメータ化する。TP では `gamma` をスカラ→`gamma[ic]` (=γ_mix) 配列読みに変えれば `κ=γ_mix-1` が厳密一致し、`sonic[ic]` は既に混合整合。frozen-coefficient (毎反復凍結再構成) で収束解は不変。**実装済** (per-cell `gamma[ic]`, M1)。

### 4b. 化学種陰解法の緩和整合 (matched-relaxation scalar-DPLUR)

理論は 本ドキュメントの「理論」節。多成分 TP の陰解法 (block-DPLUR) で、流れ 5 変数 (5×5 block・`nStepInner` sweep・`implicitRelax` 緩和) と化学種 `ρY_s` (従来 commit 後に 1 回だけ点陰的 forward-Euler) の**擬似時間緩和ミスマッチ**を、化学種を流れと同一の緩和で前進させて解消する (要因2 の恒久修正)。

- **設定** (`input/solverConfig.{hpp,cpp}`): `time.deltaT.speciesImplicitCoupling` (既定 0=従来 segregated 点陰的・ビット不変, 1=緩和整合 scalar-DPLUR)。`timeIntegration==11` かつ `nSpecies>=2` でのみ作用。
- **バッファ** (`variables.cpp registerSpecies`): 化学種ごとに `dq_roY{s}` / `dq_roY{s}_old` を追加登録 (Jacobi sweep の new/old)。
- **sweep カーネル** (`speciesTransport_d.cu species_dplur_sweep_d`): `implicit_defect_correction_d` (流れ scalar-DPLUR) を化学種 1 本にミラー。セルが面 (`map_cell_planes_index/_d`) を走査し、凍結 `massflux[ip]` から
  - 対角 `D = V/Δτ + transport_diag_Y{s}`  (`transport_diag` は流出 `max(±ṁ,0)/ρ` + 粘性時 Fick 拡散対角、`speciesTransport_d_wrapper` が確定済)、
  - 非対角 = 流入 `max(∓ṁ,0)/ρ_nbr` × `dq_old[nbr]` (`ρ_nbr` は基準 `roN`=ρⁿ で `transport_diag` の ρ と整合)、
  を集め `dq_new = implicitRelax·(res_roY{s}+Σ非対角)/D` を書く。流れと同じ凍結残差・`dt_local`・`implicitRelax`。
- **ドライバ** (`speciesImplicitDPLURSolve_d_wrapper`): `dq_*` を memset→`nStepInner` 回 sweep (毎回 old↔new swap)→`species_commit_correction_d` で `ρY_s = ρY_s^N + dq` (floor 0) を commit。`main.cpp implicitNonlinearUpdate` は `speciesImplicitCoupling==1` で `speciesUpdateOuter`(ρY_N=ρYⁿ)→本ソルバ→`speciesRenormalize`(ΣρY=ρ)→`speciesPrimitive` を呼ぶ。0 のときは従来経路 (`speciesTimeIntegration(0)`) を維持。
- **退化**: `nStepInner=1`・`implicitRelax=1`・非対角 0 で従来点陰的に一致。化学種拡散の非対角は省く (点陰的のまま。定常で収束先不変)。**完全結合 (5+N block) ではなく**、緩和率統一だけで安定 `cfl_pseudo` 上限が上がるかを切り分ける段階的打ち手。検証は `case/16.nozzle_wys` hot N2+H2O (plan `thermophysics-species-implicit-coupling.md`)。

### 5. 設定 `input/solverConfig.{hpp,cpp}`

`physProp` に追加 (いずれも任意、`thermalMethod==2` で使用):
- `species: [N2, O2, ...]` — 混合構成 (順序が index $s$)。
- `speciesDBFile` — 外部 NASA-9/LJ DB (yaml)。空なら内蔵 DB。
- `speciesDiffusionMethod` (0:定数Sc / 1:kinetic 混合平均)、`Sc`、`Sc_t`。
  乱流シュミット数 `Sc_t` は `physProp.Sc_t` または `turbulence.turbulentSchmidt` で設定する
  (両方あれば後者を優先。`turbulence.turbulentPrandtl` と同じブロックで揃えられる)。既定 0.7。
- `nSpecies` は `species` の要素数。未指定で `thermalMethod==2` なら既定 N2 単成分。

### 5b. 多成分化学種輸送 (M2) `cuda_forge/speciesTransport_d.{cuh,cu}`

化学種は NS の 5 元ブロックには結合させず、RANS k/ω と同じ汎用スカラ輸送コア
`scalarTransport_d` (`ScalarTransportDesc`) を化学種ごとに再利用して **segregated** に解く。
保存質量分率 $\rho Y_s$ の移流方程式 $\partial_t(\rho Y_s)+\nabla\!\cdot(\rho Y_s\mathbf{u})=0$ を、
対流流束カーネルが面で確定した `massflux` (SLAU の $\dot m$) による 1 次風上で組み立てる
(M2 は移流のみ。拡散は M4)。

- **変数登録**: `variables::registerSpecies(nSpecies)` が `cellValNames`/`output_cellValNames`
  (非 const 化) へ 1 化学種あたり `roY{s},Y{s},roY{s}N,roY{s}M,res_roY{s},res_roY{s}_m,`
  `transport_diag_Y{s},src_jac_Y{s}` を追加し、`c`/`c_d` マップにも空エントリを作る。
  `allocVariables` 前に 1 度だけ呼ぶ。**`nSpecies<=1` では一切登録せず M1 と同一経路**
  (回帰がバイト一致)。
- **device roY 配列**: `speciesInit_d` が `flow_float*[nSpecies]` を 1 度構築し、
  `dependentVariables_d` (従来 `nullptr`) へ渡して混合則 thermo ($R_{mix},c_{p,mix},\gamma_{mix}$) を有効化。
- **原始量・BC**: `speciesPrimitive_d` が $Y_s=\rho Y_s/\rho$ を更新。`applySpeciesBoundaries` は
  全境界を Neumann (zero-gradient ghost) で埋める (組成依存エネルギー BC は後続)。
- **時間積分・実現可能性**: `scalarTimeIntegration_d` (floor=0) で更新後、
  `speciesRenormalize_d` が $\rho Y_s\ge0$ にクランプし $\sum_s\rho Y_s=\rho$ へ再スケール
  ($\sum_s Y_s=1$)。`roY{s}N/M` は `speciesUpdateOuter/Inner` が D2D copy で NS の N/M に同期。

### 5c. 輸送係数 — kinetic theory (M3) `cuda_forge/thermo_d.cuh` + `gasProperties_d.cu`

`thermalMethod==2` で `viscMethod==2` を選ぶと、混合粘性 $\mu$ と熱伝導率 $\lambda$ を
Chapman-Enskog + Wilke/Mason-Saxena で per-cell に評価する (LJ パラメータ `sigma_LJ`,`eps_kB` を使用)。

- **単成分** (`thermo_d.cuh`): $\mu_s=2.6693\times10^{-6}\sqrt{M_s[\mathrm{g/mol}]\,T}/(\sigma_s^2\,\Omega^{(2,2)*})$ [Pa·s]、
  $\lambda_s=\mu_s(c_{p,s}+1.25R_s)$ (modified Eucken)。衝突積分 $\Omega^{(2,2)*},\Omega^{(1,1)*}$ は
  Neufeld et al. (1972) 閉形式。
- **混合**: Wilke $\mu_{mix}=\sum_i X_i\mu_i/\sum_j X_j\phi_{ij}$、Mason-Saxena $\lambda_{mix}$ (同 $\phi_{ij}$)。
  モル分率 $X$ は `thermo_X_from_Y` で質量分率から変換。
- **per-cell `thermCond`**: `variables` に `thermCond` を登録し、`gasProperties_d` が全 viscMethod で
  値を書く (0/1 は一定 `cfg.thermCond`、2 は $\lambda_{mix}$)。`viscousFlux_d` の熱流束はスカラ→
  面平均 per-cell 配列に変更。viscMethod 0/1 では一定値のため面平均は同値で**既存挙動不変**。
- **二元/混合平均拡散** `thermo_Dbinary`/`thermo_Dmix_species` も追加 (M4 の化学種拡散で使用)。

### 5d. 化学種拡散 + エネルギー結合 (M4) `cuda_forge/speciesTransport_d.cu`

粘性ケース (`viscMethod!=0`, `nSpecies>=2`) で `species_diffusion_d` が面ループで化学種拡散を加える。

- **Fick 拡散**: $J_s=\rho D_s\,\nabla Y_s$ (over-relaxed 法線, `viscousFlux_d` と同じ幾何)。
  $D_s$ は `speciesDiffusionMethod==1` で混合平均 (`thermo_Dmix_species`)、`==0` で定数 Schmidt
  $D=\mu/(\rho\,\mathrm{Sc})$。乱流寄与 $\mu_t/(\rho\,\mathrm{Sc}_t)$ を全種に加算。
- **質量保存補正**: $J_s^\* = J_s - Y_s\sum_k J_k$ で $\sum_s J_s^\*=0$ を厳密化 (混合平均拡散は
  そのままでは $\sum J\neq0$)。`res_roY_s += J_s^\*`。
- **エネルギー結合**: 化学種拡散が運ぶエンタルピー $\sum_s h_s(T)\,J_s^\*$ を `res_roe` に加える
  ($h_s$ は NASA 絶対基準で `roe` と整合)。
- **汎用拡散との二重計上回避**: 汎用スカラコア `scalarTransport_d` の拡散 (μ ベース) は
  `ScalarTransportDesc.diffusion` フラグで制御し、化学種は `0` (移流のみ) とする。RANS k/ω は `1`。
  これを怠ると層流粘性由来の $\nu$ 拡散が Fick 拡散に重畳する (検証で発見・修正)。

### 5e. 境界条件の TP 整合 `cuda_forge/boundaryCond_d.cu`

`thermalMethod==2` で各 BC の ghost / 境界状態を NASA 熱力学で再構成する (CPG 前提の定数 γ,cp を排除)。
各カーネルは `int thermalMethod, const SpeciesThermo* sp` を受け取り、`thermalMethod==2` 分岐で
`thermo_R_species`/`thermo_cp_mass`/`thermo_h_mass`/`thermo_isentropic_from_total_single`/
`thermo_isentropic_from_total_Ps_single` を用い、`γ_mix(T)`・`R_mix`・`e_NASA(T)` で
`roe`/`sonic`/`T` を整合させる。CPG 分岐は代数的に不変。

- **M1 整合済**: `slip` / `wall`(断熱) / `inlet_Pressure` / `outflow`。
  - **`inlet_Pressure` TP 分岐の閉包 (2026-08-16 修正)**: 旧実装は内部静圧から作る $M_b$ と内部速度から作る
    $M_c$ を rf=0.5 で混ぜていたが、静圧参照は $P_{ic}\leftrightarrow M_b$ の帰還で**大径低速入口 (M≈0.01–0.03)
    が振動する** (case/42 M4.2: 残差床 rms_ro 3e-2・入口 P/Pt 1.03 振動 / M6: 16k step 以降に指数成長し
    P/Pt 6 で発散, run_0040)。CPG 分岐は速度参照のみで同条件を rms_ro 1e-9 まで収束させるので、TP 分岐も
    **rf=0 (速度参照のみ)** に統一した。効果: M4.2 の残差床 3e-2 → 9e-7・入口 P/Pt 0.9986–0.9993 (run_0050)、
    M6 は 36k step で軸 M 変動 3e-5 に凍結 (run_0051)。
- **本節で追加**: `outlet_statPress`(順流 + backflow 等エントロピー), `wall_isothermal`(等温壁 ghost),
  `inlet_Pressure_dir`(全状態と外挿静圧から NASA 等エントロピー反転), `inlet_uniformVelocity`。
- **新規 thermo**: `thermo_isentropic_from_total_Ps_single` (全温・全圧 + 静圧 Ps から
  `s0(Ts)=s0(Tt)+R ln(Ps/Pt)` の Newton 反転で Ts, ρ, |u|)。
- **組成依存超音速入口 (M5)**: `inlet_uniformVelocity` を多成分化。`readBcondConfig` が `inlet_*`
  種別に `Y0..Y{n-1}` を type-1 uniform-float bvar として動的登録 (未指定なら `Y0=1`・他 0 で単成分後方互換)。
  `bcond::Yb_d` (device `flow_float*[nSpecies]`) を `inlet_uniformVelocity_d_wrapper` が遅延構築し、
  カーネル TP 分岐が `thermo_R_mix`/`thermo_e_mix`/`thermo_gamma_mix` で `T=P/(ρR_mix(Y))`,
  `roe=ρ(e_mix+ek)`, `sonic=√(γ_mix P/ρ)` を入口組成で再構成。化学種 ghost は
  `speciesBoundary_d_wrapper` を kind 分岐し**入口は Dirichlet** (`roY_s[ig]=ρ[ig]·Y_s^in`)、他は Neumann。
- **制約**: slip/wall/outlet の ghost roe/sonic はまだ `sp[0]` 固定 (非粘性では SLAU 流束がセル組成で
  再構成するため支配的影響なし)。`inlet_uniformVelocity` 以外の入口種別の組成指定も後続。
- **注意**: TP の `roe` は NASA 絶対基準。CPG res からの warm-start 不可。IC は `roe=ρ(e_NASA(T)+ek)` で生成する。
- **startup 注意**: 静止 IC (u=0) + `outlet_statPress` は出口が backflow 分岐に落ち、`Pt`/`Tt` 未指定で
  0/0 NaN 化する。IC に下流速度を与え出口に `Pt`/`Tt` を指定する (`case/28` 参照)。

### 6. 検証 (M1)

`case/13.nozzle_H/run_0001_tpgas_n2` (TP-N2) と `run_0002_slau_cpg` (CPG 参照) を SLAU・非粘性・warm-start で実行。`validate_tpgas.py` が中心線の M・T・P と全エンタルピー $H_0=h(T)+½|u|^2$ を抽出し、NASA(=CEA) 等エントロピー曲線と比較する。単体の NASA-9 値は NIST と一致 (N2: $R=296.8$, $c_p(300)=1040$, $c_p(2000)=1284$ J/kgK, $h(298.15)=0$、$e\to T$ 反転誤差 $\sim10^{-11}$)。

### 6b. 性能最適化 (consumer GPU の FP64 律速対策, M6)

TP-gas 経路は CPG に比べ大幅に遅い。プロファイル前の分析で**主因はアルゴリズムではなく演算精度**と判明した。

#### 背景: GeForce (consumer Ampere/Ada) は FP64 が FP32 の 1/64

開発機 RTX 3060 (GA106, CC 8.6) は SM あたり FP64 コアが 2 基のみで、**倍精度スループットが単精度の 1/64**。一方 `thermo_d.cuh` の TP コアは桁落ち回避のため**内部計算が全て `double`**。CPG 経路は全て `flow_float` (FP32) の代数演算なので、同一演算でも TP は 1op あたり 32〜64 倍遅い土俵で回る。さらに `log/pow/exp` を毎セル・毎面・毎 Newton 反復で評価するため、ホットパス (`dependentVariables_d` 温度反転 / `SLAU_d` 面エンタルピー / `gasProperties_d` kinetic 輸送) で TP が支配的に重くなる。

#### 施策 (M6)

精度方針を「**桁落ちのある量だけ double、それ以外は FP32**」に改める。NASA 絶対エンタルピーは生成エンタルピーを含み大きさ ~$10^7$ J/kg で増分が小さく FP32 では桁落ちするため double を保持する。一方、輸送係数・cp・混合則・面エンタルピーの増分は FP32 で安全。

1. **輸送係数の FP32 化 (①)**: `thermo_omega22/11`・`thermo_mu_species`・`thermo_lambda_species`・`thermo_wilke_phi`・`thermo_mu_mix`・`thermo_lambda_mix`・`thermo_Dbinary`・`thermo_Dmix_species` を FP32 (`expf/powf/sqrtf`) 化。これらは桁落ちが無く、`gasProperties_d` (viscMethod==2) の $O(n_{sp}^2)$ pow/exp/sqrt と `species_diffusion_d` の拡散係数評価という最重量パスを FP32 へ移す。非粘性回帰 (viscMethod==0) は不変。粘性ケースは μ/λ を NIST 値と再照合 (FP32 で十分な精度)。
2. **Newton ループの cp+h 融合評価 (②)**: `thermo_T_from_e`/`thermo_T_from_h` は毎反復 `thermo_h_mix` と `thermo_cp_mix` を**別々の種ループ**で呼び、係数選択・$T$ 冪・`log(T)` を二重計算していた。両者を 1 スイープで返す `thermo_cph_mass`/`thermo_cph_mix` を追加し、係数・$T$ 冪・`lnT` を共有。反復あたりの多項式仕事がほぼ半減。演算の並べ替えのみで結果はビット同等。`dependentVariables_d` の反転後ブロック (cp/h/gamma 再計算) にも適用。
3. **per-cell R_mix のキャッシュ (③)**: `SLAU_d` の面エンタルピーで `thermo_R_mix(YL/YR)` を**毎面で種ループ再計算**していたが、これは再構成前のセル組成だけの定数。`dependentVariables_d` で `Rmix[ic]` を 1 回計算して `gamma/cp` と並べて格納し (CPG は $R=(\gamma-1)c_p/\gamma$)、`SLAU_d` は読むだけにする。面数 ≫ セル数なので有効。

> 残レバー (M6 後続): NASA エンタルピー反転自体の FP32 化は、生成エンタルピー ≈ 0 の不活性種 (N2/He/air) では安全だが、燃焼種では sensible enthalpy ($h-h(T_{ref})$) への再定式化が必要。最大の追加効果は温度ルックアップテーブル (cp/h/s を温度グリッドで事前計算し線形補間、`log/pow/Newton` を全廃) だが結果が微変するため別マイルストーン。

### 7. マイルストーン状況

- M1 (本実装): 単成分 TP, NASA-9, Newton 反転, SLAU/ROE TP 整合 — 実装済・検証完了 (衝撃管 CEA 照合)。
- M2: 多成分輸送 (`speciesTransport_d`), 実現可能性+再正規化, 化学種 BC/IC/出力 — 実装済・検証完了。
  単成分回帰がバイト一致、[N2,N2] 識別子トレーサで流れが単成分と完全一致、[O2,N2] 実 2 ガスが安定・$\sum Y=1$。
- M3: kinetic theory 輸送係数 (`gasProperties_d` 拡張) — 実装済・検証完了。μ/λ が NIST と一致 (N2: μ(300)=1.81e-5, μ(1000)=4.15e-5, λ(300)=0.0255; air μ=1.86e-5)、device viscMethod=2 が式と FP32 一致、非粘性回帰不変。
- M4: 化学種拡散 + エンタルピー拡散エネルギー結合 — 実装済・検証完了。静止 N2/N2 トレーサの拡散が kinetic theory 自己拡散係数と 1.3% 一致、ΣJ=0 で種質量保存 (7 桁)、O2/N2 粘性で安定。
- M5: 組成依存超音速入口 BC (`inlet_uniformVelocity` 多成分化) + SLAU 対流流束の多成分エンタルピー — 実装済・検証完了。`case/28.cutler_coaxial_jet` の binary [He,N2] 同軸ジェット (Mach1.8) が 1次/2次でともに NaN/発散なし、$\sum Y=1$、中心 He コア (Mach1.8 He, Ux~1280 m/s)・coflow N2 を再現。CPG 単成分回帰不変。定量 Cutler 照合 (粘性 + RANS-SST + 乱流化学種拡散) は後続。
- M6: TP-gas 性能最適化 (§6b, 輸送係数 FP32 化 / cp+h 融合反転 / per-cell `Rmix` キャッシュ) — 実装済・検証完了。① CPG 単成分回帰は M6 前後バイナリの差が CUDA atomicAdd 非決定性 (旧×旧の run 間差) と同等で系統的変化なし。② FP32 輸送が double と相対 ~1e-8 一致 (host 単体テスト `solver_density_cuda/tools/test_m6_transport.cpp`)、NIST と N2 μ 1.5%/λ 1.9%・AIR μ 0.3%/λ 4.3%・He μ 1.7%。③ TP binary [He,N2] スモーク (`case/28` `run_0026_he_n2_m6smoke`, 500 步) が NaN なし・$\sum Y=1$ (1.2e-7)・$0\le Y\le1$・ro/P/T 物理的。速度 (native RTX 3060, 76k セル発達場 per-step): 3 成分粘性 **1.56×** (46.07→29.56 ms)、binary 粘性 1.31×、binary 非粘性 1.02× — 支配項は FP32 輸送 (viscMethod=2 × $O(n_{sp}^2)$)。詳細は plan 変更ログ 2026-06-10。
- M7: zero-gradient/反射系 BC の組成依存エネルギー化 (`slip`/`wall`/`wall_isothermal`/`outlet_statPress` の ghost `R/e/cp/γ` を `sp[0]` 固定→内部セル混合則) — 実装済・検証完了。`bc_cell_Y` ヘルパ + `thermo_isentropic_from_total_mix`/`thermo_s0_mix` 追加。これまで Cutler (species[0]=He) で空気/coflow 領域の壁・出口が He 物性で熱再構成されるバグを解消。単成分 TP 回帰はカオス的非決定性内で不変、Cutler binary [He,N2] スモークが NaN なし・$\sum Y=1$。`inlet_Pressure`/`inlet_Pressure_dir` (指定組成 `Yb` 要・多成分未使用) は後続。詳細は plan 変更ログ 2026-06-10。
