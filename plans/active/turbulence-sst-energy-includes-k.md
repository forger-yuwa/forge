# SST: 全エネルギーに乱流運動エネルギー ρk を含める (R3, SU2 形)

## メタ

- **area**: `turbulence / thermophysics / convection`
- **status**: `in_progress` (2026-09-08 実装済・検証中、既定値未決)
- **related_docs**:
  - `methods/turbulence/theory.md` §7 (圧縮性補正・等方項), `methods/turbulence/implementation.md` (整合オプション)
  - `methods/convection/implementation.md` (SLAU 圧力流束・全エンタルピー), `methods/diffusion.md`
- **related_plans**: `turbulence-sst-consistency-options.md` (分離型オプション 4/5 の後継), `accepted/turbulence-sst-node-corner-heating.md` (H2 の出所)
- **created**: `2026-09-08`
- **owner**: `CFD Dev`

## 1. 目的

現行の保存量は E_m = ρ(e + u²/2) で乱流運動エネルギー k を含まない。せん断生産 P_k は μt の粘性仕事として直接熱に入り、
等方応力 −(2/3)ρk δᵢⱼ は入っていないため、平均流と k の間のエネルギー収支は局所平衡 P_k ≈ ε でしか閉じない
(codex 2026-09-08 H2)。分離型の代替 (`sstEnergyKSource`, `sstIsotropicStress`) は E_m 形のまま源・応力で補うもので、
有限刻みで ΣV(ρE_m + ρk) が保存せず、境界閉包も未完備だった。

本計画は保存量を **E_t = ρ(e + u²/2 + k)** (SU2 と同形) に切り替え、k 式の P_k/ε を E_t 内部の配分にする。
完了時: 断熱・定常で h + u²/2 + k が面単位で保存し、全温・全圧の後処理が k 込みで自動的に整合する。
分離型オプション 4/5 は本計画で置換 (廃止)。

## 2. スコープ

- **やる**: E_t 定式化 (config `sstEnergyIncludesK`, 初期既定 0 → 回帰後 1)、p* = p + (2/3)ρk の圧力項一括化、
  エネルギー流束への k 拡散項、T/P の算出経路の全更新、後処理・IC・restart の整合、回帰試験。
- **やらない**: k 式そのものの変更 (生産・散逸・リミッタは現行)。壁関数の壁 k (k_w≠0) の p* 扱いは実装するが物理検証は別。
  LES/WALE (k を解かない) は k=0 で不変。

## 3. 関連 docs と前提

- 理論: Wilcox (2006) §5 / SU2 の `CEulerVariable::SetPrimVar(…, turb_ke)` (P = (γ−1)[E − ½ρ|u|² − ρk] で k を除いて圧力を出す形)。
- 前提: `sstOmegaProdFromPk` / `sstSigmaBlend` 既定 ON (turbulence-sst-consistency-options §2.1) の上に積む。

## 4. 設計方針

### 4.1 収支 (なぜ交換ソースが消えるか)

E_m 形: ∂(ρk)/∂t + … = P_k − ε, ∂E_m/∂t + … = −P_k + ε (これが `sstEnergyKSource`)。
両式の和 ∂E_t/∂t + … = 0 なので、E_t 形では交換ソースが現れない。k 式は引き続き解き、E_t のうち k がいくらかを決める。
温度は残り T = (E_t − ½ρ|u|² − ρk)/(ρ c_v) (TP/多成分は e_mix の反転) から出す。

### 4.2 流束に移る項

| 式 | 追加・変更 | 対になる k 式の項 |
|---|---|---|
| 運動量 | 圧力を p* = p + (2/3)ρk に (SLAU の圧力流束 p̃、境界 bvar/ghost の圧力、node 壁半割面 pressure-only、軸対称 hoop 源 p*/r) | k 生産の等方項 −(2/3)ρk∇·u (dilatation 2) |
| エネルギー対流 | H* = (E_t + p*)/ρ = h + u²/2 + k + (2/3)k | k の対流 |
| エネルギー拡散 | + (μ + σ_k μt)∇k·S (同じ面流束・同じ σ_k ブレンド) | k の拡散 |
| 粘性仕事 | 現行どおり μ_eff τ·u (deviatoric)。内部面の −(2/3)ρk S (`sstIsotropicStress`) は撤去 (p* が担う) | — |

### 4.3 使い分け (p か p* か)

- p* を使う: 力と仕事 (上表)。
- p を使う: EOS (T, ρ)、音速、SLAU の Mach 依存スイッチ (χ, f_p)、壁圧の出力、出口 `Ps` の指定値 (熱力学圧; bvar 側で p* = Ps + (2/3)ρ_b k_b を組む)。
- 面での p* は p と ρk を別々に再構成して面で合成 (k は現行の面平均、新しいリミッタ配列は増やさない)。

### 4.4 陰解法

block-DPLUR の対角ブロックは p = (γ−1)(E_t − ½ρ|u|² − ρk) の ∂p/∂E_t = γ−1 (不変)、∂p/∂(ρk) はラグ (k は分離解法)。
p* の (2/3)ρk も源同様にラグ扱い (定常では問題なし。非定常 dual-time は subiter 内で更新されるので整合)。

### 4.4a 実装形: 分割保持 (2026-09-08 決定, 実装済)

保存量 `roe` は**平均流 E_m のまま格納**し、エネルギー残差だけを E_t の流束 (§4.2) で組む。流れの更新は E_t の更新と見なせる
(残差が E_t のものだから) ので、k 更新後に `roe -= Δ(ρk)` とすると各 CV で E_t^{new} = E_t^{old} + Δt·R(E_t) が**厳密に**
成立する (k 式の点陰化・床置きの如何によらず)。Δ(ρk) の基準は流れの更新形に合わせる: explicit RK は各 stage を N から組み直す
ので Δ = ρk − ρk_N、point-implicit / dual-time subiter は roe += dq の増分更新なので Δ = ρk − (k 更新直前の ρk) (`begin` で退避)。
利点: 熱力学 (`dependentVariables`)・全境界カーネル・IC 貼付・restart・後処理の T は**無改変** (旧 res からの restart も変換不要)。
T は E_m から出るので k を除いた温度そのもの。§4.5 の「旧 res 変換」は不要になった。

### 4.5 後処理・IC・restart

- 全温 T0 = T + |u|²/2c_p + k/c_p、全圧はそこから (T は k 込みの E_t から出た値なので、ツールは T0 に k/c_p を足すだけ)。
- IC 貼付 (`paste_isentropic_ic`)・restart/interp・段間 index コピーは roe に ρk を含めた値を作る/運ぶ。旧 res (E_m 形) からの restart は
  roe += ρk の変換が要る (`restart_field.py` / `interp_field.py` にフラグ)。
- `res_*.h5` の T/P はソルバが書く値なのでツール側は不変。`check_quasisteady` の pmax 等も不変。

## 5. 実装ステップ

1. config: `turbulence.sstEnergyIncludesK` (0/1, 既定 0)。`solverConfig.{hpp,cpp}`, `procedures/solver-settings.md`。
2. 熱力学: `dependentVariables_d.cu` (intE = roe/ρ − ek − k, roe/Ht の組立に +ρk / +k)、TP/多成分 (`thermo_d.cuh` 経由の e_mix 反転) の入口で k を引く。
   `implicitCorrection_d.cu` / `timeIntegration_d.cu` の ∂p/∂q は不変 (確認のみ)。
3. 圧力 p*: `turbulent_viscosity` 後に Pstar 配列を作る新カーネル。`convectiveFlux_d.cu` (SLAU 圧力流束・H*), `convectiveFlux_boundary_d`
   (node 壁 pressure-only), `boundaryCond_d.cu` (入口/出口/slip/対称面/軸の bvar・ghost 圧力), `axisymmetricSource_d.cu` (hoop),
   `nodeWallDirichlet_d.cu` (roe の KE 剥ぎ取りは k を残す)。`sstIsotropicStress` の内部面項を撤去。
4. エネルギー拡散: `viscousFlux_d.cu` 内部面・境界 (node 半割面は k 拡散 skip と同じ扱い) に (μ + σ_k μt)∇k·S。σ_k は `sstF1` ブレンド。
5. 境界・IC・restart: `paste_isentropic_ic` (`design/forge_design/evaluate/ic.py`), `tools/restart_field.py`, `tools/interp_field.py`
   (旧 res 変換フラグ), `periodicNode_d.cu` (res gather に変更なし・確認)。
6. 後処理: `case/*/postproc_centerline.py`, `extract_wall_pp0.py`, `design/forge_design/metrics/extract.py` の T0/P0 に k/c_p。
7. docs: `methods/turbulence/implementation.md` (E_t 定式化の節), `methods/convection/implementation.md` (p*), `procedures/solver-settings.md`。
8. 回帰 (§6) → 既定 1 に切替、`sstEnergyKSource`/`sstIsotropicStress` を削除 (config で指定されたら警告)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~ステップ 1–3 (config・熱力学・p*)~~ | 済 (分割保持形, §4.4a)。キー 0 は前バイナリと 1e-6 で一致 (run_0309) |
| 2 | ~~ステップ 4 (k 拡散のエネルギー流束)~~ | 済 (`viscousFlux_d.cu`, 内部面) |
| 3 | 後処理の全温 T0 + k/c_p (centerline / wall_pp0 / metrics) | 未 (キー 1 を常用にする時点で) |
| 4 | 既定値の判断 (ユーザ懸念: k の跳ねが T に伝播する脆弱性) と分離型 4/5 の撤去 | 未 (§6.1 の結果と対話で決める) |
| 5 | 3D (1.04M, `mesh/nozzle_user_3d_finexy_z25.msh`) で k 込み全温が Tt を超えないこと | 未 (ローカル ~55 分) |

## 6. 検証

- **単体 / ビルド**: フルビルド (hpp 変更)。`sstEnergyIncludesK: 0` で case/16 2D node の継続がビット同一 (ab_run keys0 方式)。
- **偽残差 0**: 一様 ρk (k=const, u=0, p=const) を IC に置いた閉領域/自由流で全残差が丸め誤差 (現行の内部面のみ `sstIsotropicStress` はこれが境界セルで落ちる)。
- **相似試験**: `tools/test_scale_invariance.py _scaletest_tpl1 --alpha 1e-3` PASS。
- **検証ケース**: case/26 平板 node/cell (Cf/Schlichting、現行比 ±0.5 % 以内)、case/16 2D node/cell (壁 p/p0 現行比 ±0.3 %、T0+k/c_p が
  境界層内で Tt に対し ±1 K)、case/23 軸対称 (hoop 整合: 一様 k の自由流で残差 0)、case/16 3D 1.04M (`mesh/nozzle_user_3d_finexy_z25.msh`) で
  k 込み全温が Tt を超えないこと。すべて `check_convergence` / `check_quasisteady` の VERDICT を貼る。
- **判定基準**: 上の許容内で、断熱ケースの max(T0 + k/c_p) − Tt ≤ +1 K (Pr 効果分は別途 2D 参照値と比較)。

### 6.1 結果 (2026-09-08, ローカル RTX 3060)

| 試験 | run | 結果 |
|---|---|---|
| キー 0 のビット同等 | case/16 2D node `run_0309_ek_node_off` vs `run_0305_sstdef_new` (前バイナリ) | P/T/k/ω ≤ 7e-6, 壁 p/p0 同一 |
| 相似試験 | `_scaletest_tpl2` (キー 1) α=1e-3 | PASS (ω 9.2e-4, 他 ≤1e-5)。※最初の FAIL はテンプレートに IC h5 を入れ忘れた不備 (既定 IC のコールドスタート) で、キー 0 でも同じ挙動 |
| 一様乱流減衰 (周期箱 node, RK4 explicit, u=0, k=0.1, ω=1, 400 step) | case/09 `run_0046_sst_hdt_uniform_ek1` / `run_0047_..._ek0` / `run_0048_sst_hdt_sine_ek1` | **ek1: ΣV(E_m+ρk) の変化 ≤1.3e-7、max|u| 1.9e-7 (一様 p* で偽の力なし)、k 0.100→0.0847 の減少分が c_v ΔT に 5 桁一致 (Δe=+0.01534 = −Δk)**。ek0: E_t が −0.81 % 減 (k の散逸が消える)、T 不変。sine (k=0.1(1+0.5 sin x)): p* 勾配で u≈0.02 の流れが立つが ΣV E_t 変化 ≤1.1e-7、質量 4e-8 |
| 一様乱流減衰 **dual-time** (timeIntegration 11 + dualTime 1) | case/09 `run_0049_sst_hdt_uniform_ek1_dual` / `run_0050_sst_hdt_sine_ek1_dual` | エネルギー行の BDF に ρk を含める修正 (`addUnsteadyTimeTerm_d`) 後: E_t 変化 +3.5e-7 / −3.1e-6 (subiter 収束で決まる)、c_p ΔT = +0.02148 (explicit と同一)。**修正前は E_t が −0.81 % で T 不変** (subiter 内の roe 補正は dq→0 で消えるため無効だった) |
| case/16 2D node/cell (陰解法 point-implicit, 増分形補正) | `run_0310_ek_node_on` / `run_0311_ek_cell_on` (キー 0 = run_0305 / run_0308 との差) | 壁 p/p0 @16.4/45.6/85: node 0.3658/0.2575/0.1876 → 0.3659/0.2577/0.1877, cell 0.3644/0.2579/0.1872 → 0.3645/0.2581/0.1874 (**+0.03〜+0.08 %**)。max(T0−Tt) 1.96→1.56 K (node) / 2.00→1.58 K (cell)、境界層内 max(T0+k/c_p−Tt) +1.38 K。k/μt ≤ 数 %。`check_convergence`: node rms_ro 3.0 桁↓・roK 3.6 桁↓・roe 1.9 桁↓ still converging、cell plateau (既知)。切替直後 rms_roK 7e-3 / rms_roe 1e-1 のスパイク→6000 step で床 (rms_roe の床は 1.1e-3 で off の 3e-4 より高い = E_t 残差に k 項が入るため) |
| case/26 平板 node | `run_0027_ek_on` (vs `run_0026_sstdef_new`) | Cf/Schlichting 0.8920/0.9307/0.9562 → 0.8918/0.9305/0.9559 (**−0.03 %**), P ≤2e-4, k 2e-3, μt 9e-3。全列 falling |

## 7. 影響範囲

- 触るモジュール: §5 のファイル。SST 以外は k=0 で不変 (LES/WALE/層流)。
- 既存ケース: 旧 res からの restart は roe 変換が要る (フラグ)。既定 1 に切り替えた後は全 SST ケースの T/T0 が境界層内で最大 k/c_p (1〜4 K) 変わる。
- docs: `methods/turbulence/implementation.md`, `methods/convection/implementation.md`, `methods/index.md` (節追加時)。

## 8. 完了条件

- [ ] 関連 `methods/` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6)
- [ ] `status: done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/` へ移動、`plans/README.md` 同期
- [ ] `turbulence-sst-consistency-options.md` の 4/5 を superseded と記載

## 9. 変更ログ

- `2026-09-08` — 初稿。ユーザ決定「R3 直行」(対話 2026-09-08): 分離型 `sstEnergyKSource`/`sstIsotropicStress` (p* 一括化案を含む) は本計画で置換。
- `2026-09-08` — 分割保持形で実装 (§4.4a; `sstEnergyIncludesK`, 既定 0)。周期箱減衰で E_t 厳密保存 (explicit / dual-time)、相似試験 PASS、2D node/cell/平板の回帰 (§6.1)。経路別の補正: explicit RK = roe −= (ρk − ρk_N)、steady point-implicit = 増分形、dual-time = エネルギー行 BDF に ρk を含める。
