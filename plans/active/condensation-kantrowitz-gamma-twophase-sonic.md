# Kantrowitz 補正の比熱比を蒸気 γ_v に正し、凝縮セルの音速を二相 frozen 音速にする

## メタ

- **area**: `condensation`
- **status**: `in_progress`  <!-- plan 段 codex レビュー 1 回目 NO-GO → 全件採用しスコープ縮小 → 2 回目 GO-with-changes (M5/m2 全件採用) → 実装 -->
- **related_docs**:
  - `methods/condensation.md` (§2 核生成「Kantrowitz 非等温補正」、§5「二相 frozen 音速」、「モデル切替 (config フラグ)」表 — 本 plan 起票時に更新済み)
- **related_plans**:
  - [condensation-nonequilibrium.md](../accepted/condensation-nonequilibrium.md) (親: 4 モーメント凝縮、Kantrowitz と一温度二相 EOS の導入元)
  - [condensation-equilibrium-eos.md](../accepted/condensation-equilibrium-eos.md) (`condEquilibrium: 2`; **本 plan では既定対象外**, §2)
  - [time_integration-general-eos-jacobian.md](../accepted/time_integration-general-eos-jacobian.md) (block-DPLUR の一般 EOS 固有系: per-cell `gamma`/`sonic` を読む)
- **created**: `2026-09-10`
- **owner**: `CFD Dev`

## 1. 目的

非平衡凝縮の 2 つの不整合を正す。

1. **Kantrowitz 非等温補正の比熱比 (純蒸気形の係数修正)**: 現実装 (`condensation_source_d`) は
   $\theta=\frac{2(\gamma-1)}{\gamma+1}b(b-\frac12)$ の $\gamma$ に**セルの気相混合** $c_{p,gas}/c_{v,gas}$
   (carrier=`cp_cell/(cp_cell-Rmix_cell)`, H2O–N2 では N2 支配 ≈1.40) を渡している。Kantrowitz/Feder の**純蒸気形**は
   係数 $\frac{2(\gamma_v-1)}{\gamma_v+1}=\frac{R_v}{c_{v,v}+R_v/2}$ (クラスタが蒸気分子との衝突で捨てられるエネルギー揺らぎ) なので
   **凝縮種 (蒸気) 自身の** $\gamma_v$ (H2O 1.331) を使う式である。純蒸気形の枠内で係数を正す (H2O: $\theta$ −15 %)。
   キャリア衝突を含めた真値との差 (carrier では $\theta$ がさらに小さい) は本 plan の対象外 (§10)。
2. **凝縮セルの音速**: 現実装 (`dependentVariables_d`) は `sonic = √(γ_mix R_mix T)` で液相分率 $g$ を無視した「全蒸気気相」の
   音速を使っている。一温度二相 EOS $p=\rho R_{eff}T$ ($R_{eff}=R_{mix}-gR_w$) と整合する **固定 $g,Y$ の frozen 音速**
   $c^2=\gamma_{2\phi}R_{eff}T$ に改め、per-cell `gamma` も $\gamma_{2\phi}$ にする (block-DPLUR LHS の固定 $g,Y$ frozen 近似)。

完了時: (1) は全凝縮種で既定 (N2 は $\gamma_v$=1.4 で数値不変)。(2) は **TP carrier H2O 非平衡 (`condEquilibrium 0`) かつ検証済み境界構成
(`inlet_Pressure` / `outflow` / `wall` / `slip` / `periodic` のみ)** に限って自動 ON。それ以外 (pure N2 / CPG / `condEquilibrium 1,2` /
`outlet_statPress` や `wall_isothermal` を含む構成) は自動 OFF (明示キー `condSonicModel: 1` で ON 可、ただし未検証で起動ログに警告)。旧挙動は A/B 用キー (`condKantrowitzGammaMode: 1`, `condSonicModel: 0`) で再現できる。
`condensation: 0` はコード経路不変。Wyslouzil 2D (case/16) node/cell で変化量を定量化する。

## 2. スコープ

- **やる**:
  - Kantrowitz $\theta$ の $\gamma$ を `CondSpeciesProps.cp/cv` (蒸気) に変更。H2O の `cp/cv` を NASA-9 200–300 K の値
    (1855 / 1393.5, $\gamma_v$=1.331) に更新 (現状 1880/1418=1.326。`cp/cv` は他で未使用)。N2 は 1.4 のまま。
  - 二相 frozen 音速 $c^2=\gamma_{2\phi}R_{eff}T$ と $\gamma_{2\phi}$ を **TP 分岐** (`thermalMethod 2`) の carrier / pure に実装。
    **CPG 分岐 (`thermalMethod 0`, pure N2 case/34 が使う) は触らない** (旧式のまま; codex C1 の N2 潜熱フィット問題)。
  - `condSonicModel` の既定を構成依存にする (`input/condSonicResolve.hpp` の純関数 `resolveCondSonicModel`; bcond 読込後 `main.cpp` で確定):
    TP carrier H2O・`condEquilibrium 0`・全 bcond が検証済み種別のとき 1、それ以外 0。解決値と理由を起動ログに出す。明示指定は尊重するが
    未検証構成なら警告。選択条件は `tests/unit/test_cond_sonic.cpp` で回帰試験する (codex 2 回目 M1)。
  - 防御: $c_{v,2\phi}\le 0.05\,c_{p,2\phi}$ または $c^2\le0$ のセルは旧式にフォールバック (H2O では $L'=c_{p,v}-c_l<0$ で起きない;
    明示 ON した pure N2 TP のための保険。既定 OFF 構成の代替ではない)。
  - A/B 用キー 2 つ (`condKantrowitzGammaMode`, `condSonicModel`)、単体テスト (`tests/unit/test_cond_sonic.cpp`)、
    `check_convergence.py` に凝縮残差列 (`rms_rog_*`, `rms_roQ*_*`) を追加し、`init==0` 特例 (低下桁数を見ない) を廃してピーク基準に統一
    (全期間ゼロの列だけ除外; codex 2 回目 M2)。報告量の時系列定常判定 (`case/16/compare_condfix.py --series`: 末尾窓 ≥3 枚のトレンド・振幅で
    STEADY / DRIFTING / OSCILLATING / TRANSIENT-UNSETTLED、全報告量に許容値、onset は閾値交差を補間、不合格は非ゼロ終了; M3)。
    cell 後処理 (mesh `CELLS/centCoords` + 壁最近接セル) を run_0197 で検証してから R4/E4 (M4)。
    `tools/test_eos_jacobian.cpp` に固定 $g,Y$ 二相状態 (mode 2) を追加し、実 `accumulate_split_jacobian_cf` (double/float) と流束 FD を照合 (M5)。
  - case/16 Wyslouzil 2D: node 5 run (HEAD 参照 / 反復 / 旧キー回帰 / γ_v のみ / 音速のみ / 両方) + cell 2 run + 凝縮 OFF 回帰 (§6)。
- **やらない**:
  - Feder の carrier 拡張 (キャリア分子との衝突を $b^2$ に数え $\theta$ を下げる形)。§10。
  - 平衡音速 ($S=1$ 拘束下の微分) への置換。
  - **pure N2 (case/34 Arthur, CPG) への二相音速適用**: `n2_latent` の 4 次多項式は低温 (例 45.2 K) で $L'=+4160$ J/kg/K となり
    $c_l=c_{p,v}-L'<0$ (熱力学的に不整合) → $g=0.2$ で $c_{v,2\phi}<0$, $c^2<0$ (codex C1)。潜熱・内部エネルギーの整合した修正と
    有効域全体の $c_{v,2\phi}>0$ 検証が先。別 plan。
  - `condEquilibrium: 2` (case/44) への既定適用: $g$ を状態量として再決定する経路で、frozen 音速の妥当性は別途検証 (§5.1)。
  - 境界 ghost / ピンの二相化 (TP 亜音速 `outlet_statPress` の ghost 再構築、node 等温壁ピンは全蒸気 EOS のまま; §7)。
  - block-DPLUR の $\rho g$ 列を含む厳密 Jacobian。LHS は固定 $g,Y$ の frozen 近似のまま (§4.2)。

## 3. 関連 docs と前提

- [methods/condensation.md](../../methods/condensation.md): §2 核生成 (Kantrowitz 式と $\gamma_v$ の根拠)、§5 一温度二相 EOS と
  「二相 frozen 音速」節 (導出)、config 表 (新キー 2 つ)。本 plan 起票時に更新済み (2 回目レビュー前に本 plan と同期)。
- 一温度二相 EOS の内部エネルギー: carrier $e=e^{全蒸気}_{gas}(Y,T)+g(R_wT-L)$、pure $e=e_v(T)+gR_vT-gL$。
  固定 $g,Y$ で $de/dT=c_p^{全蒸気}-R_{eff}-gL'(T)$ (= 温度反転 Newton の `demix`)。
- 音速の消費側: SLAU `c_hat=0.5(sonic[ic0]+sonic[ic1])` (`convectiveFlux_slau_d.inc.cuh`)、`setDT_d` (CFL)、
  block-DPLUR 固有系 (`block_dplur_jacobian_d.cuh`; TP は `c`,`gamma`,`Ht`、CPG は `c/(γ−1)` で `sonic` を読む)、境界 (超音速判定・TP 出口)、低マッハ前処理。
- **境界の状態再構築** (codex M3): TP 亜音速 `outlet_statPress` は ghost を全蒸気 $R_{mix}$ から再構築 (`boundaryCond_d.cu` ~660, 既存注記「要フォロー」)、
  node 等温壁は `pin_wall_node_temperature_d` で `P/roe/sonic` を全蒸気状態に上書き。**本 plan の検証境界は `inlet_Pressure` (亜音速) /
  `outflow` (超音速: `roe/T/sonic` 内部値コピー) / 断熱 no-slip 壁のみ**。
- H2O 蒸気 $\gamma_v$ (NASA-9, `species_db.yaml`): 200 K 1.332 / 240 K 1.331 / 300 K 1.329。定数 1.331 で十分。
- H2O 潜熱 `h2o_latent`: 273.15 K 未満は $c_l$=4228 J/kg/K 一定で線形外挿 → $L'=c_{p,v}-c_l\approx-2375$ J/kg/K (常に負)。
  273.15–373.15 K は CEA 液相多項式 ($c_l\approx4180$–4220)。よって H2O では $c_{p,2\phi}>c_p^{全蒸気}$, $c_{v,2\phi}>0$ が全域で成り立つ (単体テストで掃引確認)。

## 4. 設計方針

### 4.1 Kantrowitz の $\gamma$ (純蒸気形の係数修正)

$$
\theta=\frac{2(\gamma_v-1)}{\gamma_v+1}\,b\Big(b-\tfrac12\Big)=\frac{R_v}{c_{v,v}+R_v/2}\,b\Big(b-\tfrac12\Big),\qquad b=\frac{L(T)}{R_vT},\qquad \gamma_v=\frac{c_{p,v}}{c_{v,v}}
$$

- `condensation_source_d` で `gamma_kw = (mode==1) ? cpg/cvg : cprops.cp/cprops.cv` を作り `cond_source_vector → cond_nucleation`
  に渡す (関数シグネチャは不変、`gamma_gas` 引数の意味を「Kantrowitz 用 γ (蒸気)」と明記)。
- 既定 `condKantrowitzGammaMode: 0` (蒸気)。`1` で旧 (セル気相混合)。
- 位置づけ: **純蒸気形近似の中での係数修正**。$J_{pure}\le J_{carrier}\le J_{iso}$ は「同一の核生成障壁・前因子を固定し、
  衝突による熱除去だけを追加したモデル間」の関係であって、carrier 中の核生成率の真値の保証範囲ではない (Wedekind et al. は carrier の
  $pV$ 仕事など逆向きの寄与も区別する; codex 2 回目 m2)。carrier 中の真値に対する残差は本修正では確定しない。
- 観測見込み (合否条件ではない): H2O, T=230 K, $b\approx26$ で $\theta$ 228 → 194、$J$ ×1.17。onset の移動は小さい (指数支配) 見込み。

### 4.2 二相 frozen 音速 (固定 $g,Y$)

一温度二相 EOS $p=\rho R_{eff}T$, $e=e^{全蒸気}(T)+g(R_wT-L(T))$ の、**$g,Y$ を固定した**等エントロピー微分 (相変化なし, 液滴は同速度・同温度):

$$
c^2=\Big(\frac{\partial p}{\partial\rho}\Big)_{e,g,Y}+\frac{p}{\rho^2}\Big(\frac{\partial p}{\partial e}\Big)_{\rho,g,Y}
=R_{eff}T\Big(1+\frac{R_{eff}}{c_{v,2\phi}}\Big)=\gamma_{2\phi}R_{eff}T,\qquad
c_{p,2\phi}=c_p^{全蒸気}-gL'(T),\quad c_{v,2\phi}=c_{p,2\phi}-R_{eff}
$$

- carrier: $R_{eff}=R_{mix}-gR_w$、$c_p^{全蒸気}$=`cpmix` (NASA-9 全蒸気混合)。pure TP: $R_{eff}=(1-g)R$。
- $L'(T)$ は温度反転と同じ数値微分 `(L(T+0.1)-L(T-0.1))/0.2` (刻み依存は単体テストで確認)。
- 実装: `condensationEOS_d.cuh` に `cond_twophase_sonic(cp_allvap, R_eff, g, dLdT, T, &gamma2, &c2)` (戻り値: 適用可否)。
  `dependentVariables_d` の **TP 分岐のみ** (`sonic[ic]`, `gam_array[ic]`) で `condSonicModel==1 && g_liq > 1e-12` のとき呼ぶ
  (それ未満は従来式 → dry セルは同一入力状態に対して従来式と bit 同一。**場としての一致は主張しない**: 凝縮帯の変化は圧力・粘性・
  陰解法を介して dry 域にも波及する)。`cp_array` は全蒸気 $c_p$ のまま (粘性・熱伝導の Pr 換算用)。
- 既定解決 (`input/condSonicResolve.hpp`, `main.cpp` の bcond 読込後): `condSonicModel` 未指定 (−1) なら
  `thermalMethod==2 && condGasSpecies>=0 && condModel==1 && condEquilibrium==0 && 全 bcond.kind ∈ {inlet_Pressure, outflow, wall, slip, periodic}` のとき 1、
  それ以外 0。理由つきで起動ログに出す。明示 1 で未検証構成なら警告して従う。
- **陰解法 LHS の扱い (codex M2)**: block-DPLUR の一般 EOS 固有系は per-cell `gamma`/`sonic`/`Ht` から $\kappa=\gamma-1$, $\chi=c^2-\kappa h$ を組む。
  `gamma`=$\gamma_{2\phi}$, `sonic`=$c_{2\phi}$ を渡すことで、**固定 $g,Y$ の frozen 状態に対する** $\kappa=R_{eff}/c_{v,2\phi}=\gamma_{2\phi}-1$ が
  入る。これは NS 5 本を $g$ 固定で線形化した近似ブロックであり、$\rho g$ 列 ($\xi_g=\partial p/\partial(\rho g)|_{\rho,\rho e}=-R_wT+\kappa(L-R_wT)$) を含む
  保存系の厳密 Jacobian ではない (NS 更新後に `rog` を別更新する分離解法のまま)。LHS は収束経路にだけ効き収束解には効かないので、
  合否は `check_convergence.py` PASS (参照・新の双方) で見る。**今回使う固定 $g,Y$ ブロックは `test_eos_jacobian.cpp` mode 2 で検証する**
  (二相 EOS の流束 FD vs 実 `accumulate_split_jacobian_cf`、double と float32 の方向微分)。$\rho g$ 列込みの厳密化のみ §5.1 の後続。
- 見込み (観測項目): Wyslouzil 出口 $g=0.011$, $L'\approx-2.4$ kJ/kg/K で $c_{p,2\phi}$ +2.5 %, $R_{eff}$ −1.7 %, $\gamma_{2\phi}$ 1.401→1.379、$c$ −1.6 %。

### 4.3 変えないもの (確認事項)

- SLAU 面温度 $T=p/(\rho R_{eff})$・面エンタルピー $-gL$ 補正 (実装済み) は不変。
- `condensation: 0` は kernel 内の分岐条件が全て偽で従来式。**CPG 分岐は一切変更しない** (case/34 は sonic も Jacobian も不変)。
- 境界: `outflow` (超音速) は内部 `sonic` をコピーするので二相音速がそのまま出る。TP 亜音速 `outlet_statPress` と node 等温壁ピンは
  全蒸気 EOS で ghost/ピンを作る既存近似のまま (本 plan では未検証構成; §7)。
- 出力 `sonic` (level 2) で実カーネル値を照合できる。

## 5. 実装ステップ

1. `cuda_forge/condensationProperties_d.cuh`: `condProps_H2O` の `cp/cv` を 1855/1393.5 に (コメントに $\gamma_v$=1.331 と NASA-9 根拠)。
2. `input/solverConfig.hpp/.cpp`: `condKantrowitzGammaMode` (既定 0), `condSonicModel` (既定 −1=自動 → 上記規則で解決) を `condensation` に追加。
   **struct 変更なので full rebuild** ([stale-build-struct-layout-trap])。
3. `cuda_forge/condensationSource_d.cu`: `gamma_kw` の選択と kernel 引数 (`kwGammaMode`)。`condensationSource_d.cuh` のコメント更新。
4. `cuda_forge/condensationEOS_d.cuh`: `cond_twophase_sonic` (host/device, 防御付き)。
5. `cuda_forge/dependentVariables_d.cu`: kernel 引数 `condSonicModel`、TP 分岐の `sonic`/`gamma`。
6. `tools/check_convergence.py`: `rms_rog_*` / `rms_roQ*_*` 列を検査対象に追加 (存在時)、`init==0` 特例を廃止 (ピーク基準)。
   `case/16/compare_condfix.py --series` (報告量の時系列定常判定、cell 対応)。cell 後処理を run_0197 で検証 (R4/E4 の前提)。
7. `tests/unit/test_cond_sonic.cpp` (§6) と `tools/test_eos_jacobian.cpp` mode 2 (固定 g,Y 二相)、`main.cpp` の解決ログ、`variables.hpp` の出力に `gamma` を追加 (実カーネル照合用; 無ければ)。
8. 検証 run (§6)、case/16 README の run 一覧、`methods/condensation.md` の検証節、`procedures/recommended-settings.md` §3、本 plan §9。

### 5.1 残作業 (優先順)

**残作業の正本はこの表**。

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー 2 回目~~ | 決着 (2026-09-10, §6.1): GO-with-changes M5/m2 → 全件採用 (自動適用の限定・PASS 必須・定常判定の窓化・cell 後処理・Jacobian 試験復帰) |
| 2 | ~~実装 (§5 1–6)~~ | 決着 (2026-09-10, §9): 実装済み・full rebuild |
| 3 | ~~単体テスト (§6 a–g)~~ | 決着 (2026-09-10, §9): `test_cond_sonic` ALL PASS、`test_eos_jacobian` Level1/2/3/3b PASS |
| 4 | ~~回帰 (§6 R1–R4)~~ | 決着 (2026-09-10, §9): 全て反復ノイズ床以内 |
| 5 | ~~効果 run (§6 E1–E4)~~ | 決着 (2026-09-10, §9): node 48000 step PASS+STEADY、cell 24000 step (既知の床) STEADY |
| 6 | ~~実カーネル照合~~ | 決着 (2026-09-10, §9): `verify_sonic.py` 3e-7 / 5e-7 |
| 7 | **codex result レビュー 2 回目** | 1 回目 NO-GO (2026-09-10, §6.1) は全件対応済み。2 回目は 2026-09-10 23:00 に codex の使用上限で実行できず (log のみ、`ERROR: You've hit your usage limit ... try again at Sep 11th 01:53`) → **2026-09-11 01:53 以降に `codex_review.py --stage result --base main` を再実行し、GO なら `status: done`・accepted へ移動** |
| 8 | (後続) `condEquilibrium: 2` への適用検証 | case/44 で frozen 音速の効果と収束を確認してから既定化 |
| 9 | (後続) 境界の二相整合 | TP 亜音速 `outlet_statPress` ghost / node 等温壁ピンを二相 EOS で再構築 |
| 10 | (後続) pure N2 の潜熱フィット整合 | `n2_latent` の低温 $L'>0$ ($c_l<0$) を直し、CPG 分岐にも二相音速を展開 (case/34 再検証) |
| 11 | (後続) LHS の $\rho g$ 列整合 | 固定 ρg (分離解法の実際) との差 $\xi_g$ を含む Jacobian。固定 g,Y の実ブロック試験は本 plan で実施 (§6 a-f) |
| 12 | (未確定) Feder carrier 拡張 | onset 較正で必要になったら別 plan |

## 6. 検証

- **単体 / ビルド**: `cmake --build solver_density_cuda/build` (native, full rebuild)。`tests/unit/test_cond_sonic.cpp` (nvcc host):
  (a) CPG pure 物性 (H2O; N2 も式レベルで) と TP carrier (NASA-9 N2+H2O, `cond_T_from_e_carrier`) で
      $(\partial p/\partial\rho)_e+(p/\rho^2)(\partial p/\partial e)_\rho$ を中心差分で評価し `cond_twophase_sonic` の $c^2$ と相対誤差 <2e-5;
  (b) 掃引: T∈{200,230,260,273.15±1,300,350} K × $Y_w$∈{0.005,0.011,0.05} × $g$∈{0,0.25,0.5,1}$\times Y_w$ で $c_{v,2\phi}>0$, $c^2>0$,
      $c_{2\phi}<c^{全蒸気}$ (H2O では $L'<0$ なので単調);
  (c) $L'$ の刻み依存: ΔT=0.05/0.1/0.2 で $c^2$ 相対変化 <1e-6 (273.15 K 接続点近傍は <1e-4 を許容し記録);
  (d) $g=0$ で従来式 `gmix*Rmix*T` と diff 0 (等号); $g>0$ の対象域では $c_{2\phi}<c^{全蒸気}$ (厳密不等号); float32 に落とした `gamma−1` の相対誤差 <1e-6;
  (e) Kantrowitz θ: γ_v=1.331 で旧 1.40 比 $(1+\theta_{old})/(1+\theta_v)$ と `cond_nucleation` の J 比が一致;
  (f) `resolveCondSonicModel` の選択条件: 検証済み構成 → 1、`condEquilibrium 1/2`・`outlet_statPress`/`wall_isothermal` 含み・CPG・pure → 0、明示指定は尊重;
  (g) `tools/test_eos_jacobian.cpp` mode 2: TP carrier (N2+H2O) の固定 $g,Y$ 二相状態で ||RΛL−A_FD||/||A_FD|| <1e-5、実 `accumulate_split_jacobian_cf` の
      double 組み立てが閉形式と 1e-12 で一致、float32 組み立ての方向微分 (diag·dq) が double と相対 1e-5 以内。
- **収束・定常判定ツール**: `check_convergence.py` に `rms_rog_0`, `rms_roQ{0,1,2}_0` を追加 (存在時)、`init==0` 特例を廃止。
  `compare_condfix.py --series RUN` は全 `res_*.h5` で onset (中心線 $g=10^{-3}$ の交差を補間)、壁 p/p0 の実験偏差 (x≥10 mm 平均と 21/42/52 mm)、
  出口中心 $g$/$M$/$c$、中心線 $h_0$ の最大偏差 [kJ/kg] を計算し、**末尾窓 4 枚 (最低 3 枚)** の振幅 (max−min) とトレンドで判定する:
  全量が許容値内 (壁偏差 0.1 %pt, onset 0.1 mm, g_exit 0.5 %, M_exit 0.1 %, c_exit 0.1 %, h0 0.01 kJ/kg) → STEADY、単調で超過 → DRIFTING、
  非単調で超過 → OSCILLATING、枚数不足 → TRANSIENT-UNSETTLED。不合格は非ゼロ終了。STEADY と `check_convergence` PASS の両方を満たすまで
  step を延長する (24000 から開始)。`check_quasisteady.py` への統合は後続 (汎用ツールに対象量抽出の仕組みが要る)。cell run は mesh
  `CELLS/centCoords` と壁最近接セルで同じ量を出す (run_0197 で動作確認してから使う)。
- **検証ケース** (`case/16.nozzle_wys/`, 2D Wyslouzil, TP MIXDRY+H2O, `condKantrowitz: 1`, HK 成長, 非平衡, 出力 level 2):
  node は run_0230 のプロトコル (IC=run_0213/res_24000 index コピー, `outflow`, cfl_pseudo 2)、cell は run_0197 のプロトコル。

  | # | run | 離散化 | バイナリ | キー | 目的 |
  | --- | --- | --- | --- | --- | --- |
  | R0 | `run_0331_condfix_head_ref` | node | HEAD (`build/forge.head`) | 既定 | 参照 (run_0230 再現済: onset 23.0 mm, 壁偏差同一) |
  | R1 | `run_0338_condfix_head_repeat` | node | HEAD | 既定 | 同一バイナリの反復ノイズ (回帰許容値の根拠) |
  | R2 | `run_0339_condfix_dry_head` / `run_0340_condfix_dry_new` | node | HEAD / 新 | 凝縮 OFF (run_0213 継続 1000 step) | 凝縮 OFF 回帰 (コード経路不変) |
  | R3 | `run_0332_condfix_legacy_keys` | node | 新 | `condKantrowitzGammaMode: 1, condSonicModel: 0` | 旧キー回帰 (≤ R1 ノイズ) |
  | R4 | `run_0336_condfix_cell_head` / `run_0337_condfix_cell_legacy` | cell | HEAD / 新+旧キー | 旧キー | cell 旧キー回帰 (cell の atomicAdd ノイズ ~6e-4 [cell-atomicadd-nondeterminism] 以内) |
  | E1 | `run_0333_condfix_kwgamma` | node | 新 | `condSonicModel: 0` | γ_v の単独効果 |
  | E2 | `run_0334_condfix_sonic` | node | 新 | `condKantrowitzGammaMode: 1` | 二相音速の単独効果 |
  | E3 | `run_0335_condfix_new` | node | 新 | 既定 (両方) | 新既定 |
  | E4 | `run_0341_condfix_cell_new` | cell | 新 | 既定 | cell 新既定 (node/cell 両方の必須事項) |

- **判定基準 (合否)**:
  1. 単体 (a)–(e) 全 PASS。
  2. R2: 凝縮 OFF で HEAD と新の全 VALUE が一致 (node: 差 0 または R1 と同水準)。
  3. R3 ≤ R1 のノイズ水準 (全 VALUE の max|Δ|/max|ref| と、壁 p/p0・onset・出口 g の差)。R4 ≤ cell ノイズ水準。
  4. 実カーネル照合: E3 の res `sonic` が host 式と相対 ≤1e-5 (float)、dry セルでは従来式と一致。
  5. 有界性: 全 run NaN 0、$c>0$、$g\in[0,Y_w]$、T_min > 150 K (物性有効域)。
  6. 保存性: `total_quantities.py` の $h_0$ 誤差が R0 と同水準以下 (凝縮帯で +0.6 % 型の非保存が無い)。
  7. 収束性: 参照 R0 と効果 run E1–E3 (node) で `check_convergence.py` (凝縮列込み・ピーク基準) が **PASS** かつ `--series` が **STEADY**
     になるまで延長し、その VERDICT を貼る (PASS に達しない場合は「未収束」と明記し比較を過渡比較として扱う)。
     **cell (R4/E4) は生産対象外 (ユーザ決定 2026-09-10: 「cell はもう基本使わない」) なので回帰の参考に留め、PASS を完了条件にしない**:
     dry 場から 48000 step 回して VERDICT をそのまま記録し、PASS 未達なら「未収束の準定常比較」と明記する。
  8. 回帰許容値 (codex result M2): 同一バイナリ反復 (R1: run_0331 vs run_0338, 48000 step) の変数別 max|Δ|/max|ref| を `noise_node_48000.json` に残し、
     `diff_res.py --tolfile … --factor 2` で保存量+原始量+凝縮量 29 変数 (診断量除外) を判定する。cell も同様 (run_0336 vs run_0342)。
- **観測項目 (合否にしない)**: onset x、壁 p/p0 の実験偏差、出口 $g$/$M$/$c$、$\gamma_{2\phi}$。§4 の見込みとの整合は考察に書く
  (合わなければ実装の取り違えを疑って原因を書く)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result | `2026-09-10` | 1 回目 [`notes/reviews/2026-09-10-condensation-kantrowitz-gamma-twophase-sonic-result.md`](../../notes/reviews/2026-09-10-condensation-kantrowitz-gamma-twophase-sonic-result.md) | NO-GO, C0/M5/m2 | **全件採用**。M1 (cell が PASS 未達を「検証完了」) → cell 4 run を dry 場 (run_0196) から 48000 step で取り直し (収束済み場からの継続では落ち幅が測れない)、PASS 未達なら「未収束の準定常比較」と記録して §8 を戻す。M2 (回帰許容値の根拠) → `diff_res.py` を既定キー集合 (保存量+原始量+凝縮量, 診断量除外)・欠落/形状/非有限 FAIL・変数別許容値 (`--tolfile`=同一バイナリ反復の rel × 2) に改め、node 反復も 48000 step で取り直し。M3 (「超音速で全風上」は誤り: 面法線 Mach) → §9・methods・README の因果説明を撤回し「差が小さい (ノイズ床同水準), 原因未特定」に修正。M4 (ツールが NaN を合格) → `compare_condfix --series` は非有限を UNEVALUABLE/非ゼロ終了、`verify_sonic` は gamma も合否に。M5 (methods/plan §7 の既定範囲が resolver と不一致) → `condEquilibrium 0`+検証済み境界に統一、§5 平衡拘束形の記述修正。m1 (h0 0.2669 vs 0.2720) → 数値訂正・許容 0.01 kJ/kg 明記。m2 (carrier 上下限の表現) → methods §2/plan §10 を §4.1 の限定表現に統一 |
| plan | `2026-09-10` | 2 回目 [`notes/reviews/2026-09-10-condensation-kantrowitz-gamma-twophase-sonic-plan-2.md`](../../notes/reviews/2026-09-10-condensation-kantrowitz-gamma-twophase-sonic-plan-2.md) | GO-with-changes, C0/M5/m2 | **全件採用**。M1 (自動適用が境界・緩和形を見ない) → `resolveCondSonicModel` を bcond 読込後に評価し検証済み境界 + `condEquilibrium 0` に限定、選択条件の単体試験 (§6 f)。M2 (PASS 未要求・init==0 特例) → PASS 必須、`check_convergence.py` をピーク基準に統一。M3 (末尾 2 枚判定) → 末尾窓 ≥3 枚・全報告量に許容値・onset 補間・非ゼロ終了 (`check_quasisteady.py` 統合は後続)。M4 (cell 後処理不能) → mesh centCoords 経路を追加し run_0197 で検証を前提に。M5 (Jacobian 試験の後送り) → `test_eos_jacobian.cpp` mode 2 を本 plan に戻す (§6 g)。m1 (g=0 の条件矛盾) → 等号/不等号を分離。m2 (上下限の表現) → §4.1 を限定 |
| plan | `2026-09-10` | 1 回目 [`notes/reviews/2026-09-10-condensation-kantrowitz-gamma-twophase-sonic-plan.md`](../../notes/reviews/2026-09-10-condensation-kantrowitz-gamma-twophase-sonic-plan.md) | NO-GO, C1/M5/m1 | **全件採用**。C1 (N2 で $c^2<0$) → pure N2 / CPG 分岐を対象外 (§2, §5.1 #10)。M2 (Jacobian 厳密整合の主張・CPG も変わる) → §4.2 を固定 g,Y frozen 近似に書き直し、CPG 分岐不変に (§5.1 #11)。M3 (境界 ghost 全蒸気) → §3/§7 に明記し検証境界を限定 (§5.1 #9)。M4 (凝縮残差・定常判定) → `check_convergence.py` 拡張 + `--series` 判定 (§6)。M5 (構成不足) → cell 2 run・凝縮 OFF 回帰・掃引単体・実カーネル照合を追加、未検証構成は既定外。M6 (合否条件) → 反復ノイズ基準・見積りは観測項目へ (§6)。m1 (γ_v の位置づけ) → §1/§4.1 を純蒸気形の係数修正と明記、$J_{pure}\le J_{carrier}\le J_{iso}$ |

## 7. 影響範囲

- 触るファイル: `cuda_forge/condensationProperties_d.cuh`, `condensationSource_d.cu(h)`, `condensationEOS_d.cuh`,
  `dependentVariables_d.cu`, `input/solverConfig.hpp/.cpp`, `tools/check_convergence.py`, `tests/unit/test_cond_sonic.cpp`,
  `case/16.nozzle_wys/compare_condfix.py`。
- 既定変更の影響: **Kantrowitz ON の H2O run** (case/16 非平衡) で $\theta$ −15 %。**TP carrier H2O 非平衡 (`condEquilibrium 0`) かつ検証済み境界のみ**
  で凝縮セルの音速が下がる (≲2 %) と block-DPLUR LHS が変わる。case/44 (`condEquilibrium: 2`)・緩和形 (`condEquilibrium 1`)・case/34 (pure N2 CPG)・
  `outlet_statPress`/`wall_isothermal` を含む構成・凝縮 OFF は**不変** (`resolveCondSonicModel`)。CPG 分岐は明示 1 でも変わらない。
- 境界 (codex M3): TP 亜音速 `outlet_statPress` の ghost 再構築 (`boundaryCond_d.cu` ~660) と node 等温壁ピン (`nodeWallDirichlet_d.cu`
  `pin_wall_node_temperature_d`) は全蒸気 EOS のままなので、内部 (二相) と ghost (全蒸気) の音速が不連続になる。本 plan の検証は
  `inlet_Pressure`/`outflow`/断熱壁のみ。これらの境界で二相音速を使う run は未検証 (§5.1 #9)。
- ドキュメント: `methods/condensation.md` (更新済み; 2 回目レビュー前に同期)、`procedures/recommended-settings.md` §3、case/16 README run 一覧。

## 8. 完了条件

- [x] 関連 `methods/condensation.md` の現在仕様を更新済み
- [x] 実装・検証完了 (§6 合否 1–8: node 全 PASS+STEADY・回帰は変数別ノイズ床 ×2 で PASS; cell は非ゲートの未収束準定常比較として記録)
- [ ] codex レビュー (`plan` 2 回 / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status: done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-10` — 初稿。ユーザ指摘 (Kantrowitz の γ が混合気体、音速が凝縮水を無視) を受け現状確認: `condensation_source_d` は
  `gamma_gas=cpg/cvg` (セル気相混合)、`dependentVariables_d` は `sonic=√(γ_mix R_mix T)` (g 無視) を確認。§4 の設計と §6 の検証計画を起票。
- `2026-09-10` — codex plan レビュー 2 回目 GO-with-changes (M5/m2) を全件採用 (§6.1)。`status: in_progress`、実装着手。
- `2026-09-10` — **実装・検証完了** (case/16 README「Kantrowitz γ_v + 二相 frozen 音速」節、`compare_condfix.{png,txt}`, `compare_condfix_cell.{png,txt}`):
  - 実装: `condProps_H2O` cp/cv 1855/1393.5 (γ_v 1.331); `condKantrowitzGammaMode` (既定 0=蒸気); `condSonicModel` (−1=自動,
    `input/condSonicResolve.hpp` で bcond 読込後に解決, `main.cpp` がログ); `cond_twophase_sonic` (`condensationEOS_d.cuh`);
    `dependentVariables_d` TP 分岐で `sonic`/`gamma` を γ_2φ に; 出力 `gamma` 追加; `check_convergence.py` 凝縮列 + ピーク基準。
  - 単体: `tests/unit/test_cond_sonic.cpp` ALL PASS (CPG/TP の FD 一致 ≤1e-7, 掃引 84/84 で $c_{v,2\phi}>0$, $L'$ 刻み依存 7e-8,
    float32 γ−1 1.7e-7, θ 比 1.172, resolver 12 件)。`tools/test_eos_jacobian.cpp` mode 2 (二相 g=0.009): ||RΛL−A_FD|| ≤2e-8,
    実 `accumulate_split_jacobian_cf` double 3e-16 / float32 方向微分 3e-7 → PASS。
  - 回帰: R1 node 反復ノイズ床 1.6e-5 (Y1; ρ/P/T 1–2e-6, run_0338)。R2 凝縮 OFF (run_0339/0340) core ≤7e-6 (床 2.4e-5) → 経路不変。
    R3 旧キー (run_0332) 1.3e-5 ≤ 床 → PASS。R4 cell (run_0336/0342 床 ro 3.5e-4, Uy 2.0e-3; run_0337 5.2e-4 / 3.0e-3) ≤1.5 倍 → ノイズ水準、報告量同一
    (ただし rog_0 は 3.5 倍で床超え: codex result M2。cell は dry 場から取り直して再判定, 下記 2026-09-10 追記)。
  - 効果 (node 48000 step, 全 5 run `check_convergence` PASS・`--series` STEADY; cell 24000 step は既知の床で plateau, STEADY):
    **γ_v (E1 run_0333): onset 22.96 → 22.51 mm (−0.45 mm 上流)、壁 p/p0 偏差 @21 −4.8 → −4.5 % / @42 +4.8 → +5.1 % / @52 +5.5 → +5.4 % / 平均 +2.28 → +2.46 %**。
    **二相音速 (E2 run_0334): 凝縮帯 c −1.3 % (最大 −1.63 % @ g=0.0109, 見積り −1.6 % と一致)、dt_local +1.0 %。収束場の差は小さい
    (参照との差 P/T ≤2e-6, Uy 1.1e-5, g 1.4e-5 = 同一バイナリ反復ノイズ床と同水準)**。原因は特定していない: SLAU は面法線速度で Mach を作るので
    横向き面は亜音速のままで圧力流束は c に依存する (codex 計測: 凝縮帯内部面 16,726 のうち 8,304 面が両側とも法線亜音速、旧 c に戻すと p̃ が最大 0.024 Pa 変わる)。
    「超音速で全風上だから影響 0」という当初の説明は誤り (codex result M3 で撤回)。M_exit は二相 c 定義で 1.605 → 1.631。
    新既定 (E3 run_0335 / E4 cell run_0341): γ_v の効果 + 二相 c。実カーネル照合 `verify_sonic.py` 2.9e-7 / 4.8e-7 (cell 6.0e-7 / 8.1e-7)。
    h0 中心線最大偏差: 参照 0.2669 → 新既定 0.2720 kJ/kg (Δ0.005 kJ/kg, 許容 0.01 kJ/kg 以内 = 同水準; 「同一」ではない)。NaN 0、T_min 207.5 K。
  - 見積りとの整合: γ_v は「onset ≲1 mm 上流」の範囲内 (−0.45 mm)。音速は c の変化は見積りどおりだが場への影響は「≲1 %」より小さくノイズ床と同水準。
- `2026-09-10` — codex plan レビュー 1 回目 NO-GO (C1/M5/m1) を全件採用: 二相音速を TP carrier H2O 非平衡/緩和形に限定 (CPG/N2/平衡拘束形は対象外)、
  Jacobian は固定 g,Y frozen 近似と明記、境界 ghost の未対応を明記、検証を node/cell・凝縮 OFF・反復ノイズ・掃引単体・実カーネル照合に拡張、
  合否条件を EOS 微分・保存性・有界性・収束性に置換。参照 run_0331 (HEAD, 24000 step) 取得済 (run_0230 を再現)。

- `2026-09-10` — codex result 2 回目は使用上限で未実行 (§5.1 #7)。検証済みマイルストーンとして commit (plan は `in_progress` のまま active/)。
- `2026-09-10` (codex result 1 回目の対応) — node 反復 run_0338 を 48000 step で取り直し、変数別ノイズ床 `case/16.nozzle_wys/noise_node_48000.json`
  (ρ 1.3e-6, P 9.7e-7, T 2.0e-6, U_y 1.1e-5, roQ0_0 1.5e-5, Y1 2.0e-5, h0 2.9e-5, vis_turb 1.0e-3) を作成。`diff_res.py --tolfile noise_node_48000.json --factor 2`
  (保存量+原始量+凝縮量 29 変数、診断量除外、欠落/形状/非有限は FAIL) で **R3 旧キー run_0332 は 29 変数すべて PASS**、
  **E2 音速のみ run_0334 は `sonic` (1.4e-2) 以外の 28 変数が床 ×2 以内** (音速修正は sonic 以外の場を反復ノイズの範囲でしか変えない)。
  cell 4 run (run_0336/0342 HEAD 反復, run_0337 旧キー, run_0341 新既定) を run_0196 の dry 場から 48000 step で取り直した:
  **全 run `NOT CONVERGED (plateau)`** (rms_ro 4e-8 で頭打ち、凝縮列は 5–7 桁低下; cell の既知の atomicAdd 床) → **未収束の準定常比較** として記録
  (`--series` は STEADY)。ノイズ床 `noise_cell_48000.json` (ρ 2.3e-5, P 3.5e-5, U_y 2.6e-4, rog_0 5.0e-5)。旧キー run_0337 は報告量が参照と同一だが
  場の差は床の 2〜5 倍 (ρ 7.3e-5, P 1.7e-4, U_y 1.3e-3) で `--factor 2` は FAIL: プラトー上のリミットサイクルで軌道が分かれるため 1 対の反復では床が
  定まらない。共有カーネルの同一性は node 回帰 (29 変数 PASS) で担保し、cell は非ゲート (ユーザ決定: cell は生産対象外, §6 判定 7)。
  新既定 run_0341: onset 22.94 → 22.49 mm、壁偏差 +2.24/−4.8/+4.8/+5.4 → +2.42/−4.5/+5.0/+5.3 % (node と同じ応答)、`verify_sonic` PASS。

## 10. 未確定事項

- Feder の carrier 拡張 ($b^2=kT^2\sum_i(\beta_i/\beta_v)(c_{v,i}+k/2)$, キャリア衝突で $\theta$ が下がる; Wedekind et al. 式 9–12) を入れるか。
  純蒸気形と等温 CNT の関係は §4.1 の限定 (障壁・前因子固定で熱除去だけを足したモデル間) でしか言えず、carrier 中の真値を挟む上下限ではない。
  現時点では入れない (較正課題が出たら別 plan)。
