# 二相拡散カーネルの実装設計メモと上位諮問の材料 (2026-10-02)

- plan: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md) §4.2・§6・§5.1 #1/#4
- 仕様: [`methods/condensation.md`](../../methods/condensation.md) 実装 §7c (2026-10-02 追加、未実装)
- 参照実装: [`solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py`](../../solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py)
- 作業ツリー: `feature/species-transport` (HEAD `33520773`)。**cuda_forge は未変更**。本メモはエスカレーション 6 (数値カーネル変更の前) の諮問材料。

本メモは調査メモであり方針の正本ではない。方針は plan §4.2、未確定事項の決定は上位の判断の後に plan へ反映する。

## 1. 現行コードの該当箇所

| 役割 | 箇所 | 現行の振る舞い |
| --- | --- | --- |
| 化学種拡散 (面ループ) | `cuda_forge/speciesTransport_d.cu:208` `species_diffusion_d` | $J_s=\rho_f(D_s+D_t)\Delta Y_s\,\delta/d_{cc}$ (水は総水分)、補正 $J_s^*=J_s-Y_{f,s}\sum J$ (`:307`, $Y_f$ は面の算術平均)、エネルギー $\sum h_s(T_f)J_s^*$ (`:310`)。対角 `transport_diag_Y{s}` に $\rho(D+D_t)$ を各種 (`:300-301`)。液は見ない |
| 化学種残差の初期化と呼び出し | `speciesTransport_d.cu:824` `speciesTransport_d_wrapper` | `res_roY`/`transport_diag_Y`/`src_jac_Y` をゼロ化 (`:828-833`) → 移流 → `viscMethod≠0` で拡散 (`:846-860`) |
| 液・モーメント残差の初期化 | `speciesTransport_d.cu:1334` `passiveAdvection_d_wrapper` (受動種経路) / `condensationTransport_d.cu:400-404` (旧経路) | `res_*`/`h_p_diag`/`h_p_sj` をゼロ化してから移流。**拡散なし** |
| 残差の組立順序 | `main.cpp:1879` speciesTransport → `:1881` speciesPinResidual → `:1884` condensationTransport → `:1885` condensationSource → `:1887` passivePinResidual → … → `:1913` periodicNodeGather | 液残差のゼロ化は化学種のピン除去の**後** |
| 化学種の陰的更新 | `speciesTransport_d.cu:946` `speciesImplicitDPLURSolve_d_wrapper` / `:355` `species_dplur_sweep_d` | coupling 1: 非対角は**移流のみ** (`:385-399`)、対角は `V/Δτ + transport_diag` (`:403-409`)。拡散の非対角なし |
| 再正規化 | `speciesTransport_d.cu:157` `species_renormalize_d` | ρY≥0 化 → Σ=ρ に比例配分。**液・Q には掛けない** |
| 定常の更新順序 | `main.cpp:2127-2145` (化学種更新 → 再正規化 `:2140`) → `:2150-2152` (液・モーメント更新) | 液は「再正規化済みの新しい ρY_w」を見て更新される |
| dual-time の更新順序 | `main.cpp:2355-2373` (化学種) → `:2377-2385` (液) → `:2390-2405` (物理 step 末尾の FCT) → `:2410-2413` (EOS・射影・履歴) | 化学種は FCT を通らない |
| 液の更新クランプ | `condensationUpdateLimiter_d.cuh:22` `cond_moment_update_limited_body` (受動種経路は `:155`、呼び出し `condensationTransport_d.cu:441`) | 増分 (輸送 + ソース) 全体に共通 θ (`dg_max`/`dT_max`/蒸気枯渇 `avail = Y_w − g_old` `:78-80`/蒸発上限)。avail の `roY_w` は**この時点で更新済み** |
| 受動種の FCT | `passiveFct_d.cuh:20` `passive_fct_diff_coef`、`:49` `passive_fct_lowflux`、wrapper `speciesTransport_d.cu:1770` (`qDiff` `:1786`) | 低次作用素の拡散は**トレーサ 1 本だけ** (`q == qDiff`)。液・Q・総水分は拡散を持たない |
| 実現可能性クランプ | `condensationRealizability_d.cuh:103` `cond_realizability_clamp_d` | $0\le\rho g\le\rho Y_w$、Q≥0、射影、消滅 |
| トレーサ拡散 (参考) | `passiveKernels_d.cuh:143` `passive_diffusion_d` | 受動種 1 本の Fick (μ/Sc + μ_t/Sc_t) |
| 対象 run の設定 | 本体 `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/solverConfig.yaml` | node・定常 (`unsteady 0`)・`timeIntegration 11`・`speciesImplicitCoupling 1`・`passiveImplicitCoupling 1`・`passiveScalarScheme 1`・S3 (`speciesFaceReconstruction 2`)・`condLimiterMode 1`・`[MIXDRY, H2O]`・`viscMethod 1` |

## 2. 実装設計 (案; 未確定箇所は §4)

### 2.1 新カーネル `twophase_diffusion_d` (面ループ、float32)

- 入力: `roY[s]`, `rog`, `roQ2/1/0`, `ro`, `T`, `P`, `vis_lam`, `vis_turb`, 幾何 (現行 `species_diffusion_d` と同じ `fx/sx/sy/sz/ss/cc*`)、`GasPhaseLiquid liq`。
- 面ごと: セルの $\rho_g=\rho-\rho g$、$\rho v=\rho Y_w-\rho g$、$z$ ($\Sigma z$ で正規化) → 面の $\rho_{g,f}$・$z_f$ → 気相組成の混合平均 $D_k$ (`species_transport_X_f` の気相組成と同じ) →
  $j_k^0$、$\Sigma j^0$、補正 (面の $z$ は §4-1) → 乱流係数 $c_t=a_f\mu_{t,f}/Sc_t$ で $J_t$ → $J_v=j_v+J_t(v)$、$J_l=J_t(g)$、$J_w=\mathrm{fl}(J_v+J_l)$、$J_Q$ →
  残差 `res_roY[k]`, `res_roY[iw]` ($J_w$), `res_rog`, `res_roQn`, `res_roe` ($\sum h_kJ_k+h_vJ_w-L(T_f)J_l$) へ atomicAdd。
- 対角: 非水 `transport_diag_Y[k]` に $a\rho_{g,f}D_k/\rho_{g,i}+c_t/\rho_i$ (+補正の流出分)、水 `transport_diag_Y[iw]` は**蒸気の対角** $a\rho_{g,f}D_v/\rho_{g,i}+c_t/\rho_i$、液・Q の対角に $c_t/\rho_i$。
- 面係数の共有: dual-time FCT が同じ値を使えるよう、$c_t$ (と必要なら分子の面係数) を面バッファに書く (FCT 側で再計算しない; plan §4.2「面流束を 1 回だけ作り共有」)。
- `species_diffusion_d` 側: TP carrier 凝縮のときは起動しない (水・気相種の分子拡散を二重に足さない)。凝縮 OFF・CPG carrier・pure は現行カーネルのまま (コード経路不変)。
- **g=0 のビット一致**: 面の両セルで `rog == 0` のとき、現行 `species_diffusion_d` の面計算 (乱流を $D$ に足して 1 本の補正に含める式) をそのまま呼ぶ (device 関数に切り出して共有)。
  新式は $g=0$ で代数的に一致するが丸めは一致しない ($\rho(D_s+D_t)$ を 1 つにまとめた現行と、分子と乱流を分けた新式) ので、分岐で担保する。
  注: 残差の atomicAdd 順序はビット再現しない ([res0-is-not-the-flux-state] の既知事項) ので、#5 の「ビット一致」は G0 と同じ 0 step の probe 比較か、
  atomics の影響を受けない量で判定する必要がある (判定法は #5 で決める)。
- CPG carrier: 対象外ログ「CPG carrier: 二相拡散は未対応」を `condensationInit_d` (`condensationTransport_d.cu:152`) で 1 回。cell: 経路は同じ面ループ (ghost で閉じる) だが実行検証しない。

### 2.2 呼び出し位置

- 液・Q の残差ゼロ化 (`passiveAdvection_d_wrapper` `speciesTransport_d.cu:1337-1341`) の**後**でなければならない。現行は `main.cpp:1881` の化学種ピン除去が
  液のゼロ化 (`:1884`) より前にあるので、新カーネルを `condensationTransport_d_wrapper` の末尾 (移流の後) に置き、**化学種のピン除去を新カーネルの後に移す**
  (または `:1887` の passivePinResidual の隣でもう一度呼ぶ)。周期集約 (`:1913`) の前であることは assembleResidual 内なので自動的に満たす。
- FCT 末尾の再評価 (`main.cpp:2394` の `assembleResidual(s, 1)`) でも同じカーネルが走るので、FCT の `r_H` は拡散込みになる。

### 2.3 更新 (蒸気/液の増分)

- 定常 (`main.cpp:2127-2152`) の順序を「非水種と蒸気の増分 δρ_k, δρv を作る → 液・Q の更新 (ソース・θ・床込み) で Δρg を確定 → 水を
  $\rho Y_w=\rho Y_w^N+\delta\rho v+\Delta\rho g$ で commit → 再正規化 (係数を ρg・Q にも) → 実現可能性クランプ」に変える。
  現行は化学種を commit・再正規化してから液を更新するので、**水の commit を液の更新の後へ遅らせる**必要がある (`speciesImplicitDPLURSolve_d_wrapper` の commit
  `speciesTransport_d.cu:981-988` を水だけ外す)。
- 蒸気の残差 $R_v=R_w-R_g$。ここで $R_g$ は相変化ソースを含む (`condensationSource_d_wrapper` が `res_rog` に足す) ので、ソースの扱いは §4-2。
- coupling 1 (scalar-DPLUR) の非対角は移流のみ (`species_dplur_sweep_d:385-399`)。契約 (固定点) には拡散の非対角は不要 (前処理の質だけ) なので、初版は入れない。
  収束が遅ければ後で足す。
- coupling 2 (案C, EOS クロス) は予測 δ(ρY) を流れの RHS に入れる (`speciesTransport_d.cu:1014`) ので、水を液の後に commit する順序と噛み合わない。
  初版は **coupling 2 × 二相拡散を起動時エラー**にする案 (対象 run は coupling 1)。
- dual-time (`main.cpp:2355-2385`) も同じ順序変更。物理 step の保存は「解き切った」ときだけ成り立つ (§3 S3/S5) ので、`nSubIterDualTime` 固定の打ち切りでは
  残差の和がそのまま総量誤差になる。受動種は末尾の FCT が流束形で閉じるが、総水分は FCT を通らない (§4-3)。

## 3. ホスト参照実装の結果 (2026-10-02, `python3 solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py`, 8 分 51 秒, ALL PASS)

判定は float64 の参照実装で行い、float32 は記録 ([INFO])。

| plan §6 の条件 | 実測 | 判定 |
| --- | --- | --- |
| 分子流束の構造 (i) 三成分・気相一様・∇g=0.1、dx=1 m | max\|j\| float64 2.7e-22 / float32 4.3e-13、許容 2.1e-12。旧案 j/ρ = [−1.3182e-7, −6.1364e-8, 1.9318e-7] (上位検算と一致) | PASS (旧案は FAIL で判別可) |
| 同 (i) dx=1 mm・Δg=1e-4/面 (情報) | float32 2.2e-10 (許容の 98 倍) | 許容が ∇g 基準で float32 の z 丸めを考えていない (§4-4) |
| 面恒等式 (ii) | Σj: ≤0.23×許容、J_w−J_v−J_l: ≤0.031×許容 (20 万面 × 気相 2/3 種 × 補正 z 2 通り × 順序 2 通り、すべて ≤1) | PASS |
| 3 セル判別 A | 総液量 +3.3333333 % (1 更新) / +2.8571429 % (1000 更新)、float32 +3.3333346 % / +2.8571441 % (上位検算と一致) | 保存を破る (想定どおり) |
| 3 セル判別 B (float64) | 1 更新 0、1000 更新 総液量 +2.2e-9・各種 ≤2.2e-9、min ρv 0.1、\|ΣρY−ρ\|/ρ ≤2.2e-16 (反復 ≤30) | PASS → 判別規則「A FAIL・B PASS」で契約を支持 |
| 3 セル B (float32, 情報) | 停止の丸め床 16ε: 1000 更新 +9.98e-7 / 4ε: +2.8e-7 | 停止基準しだい (§4-4) |
| 非負 2 セル | 新 (蒸気/液変数, 解き切り): ρv = [0, 0]、クランプ 0。改訂前 (総水分/液の点対角): ρv = −0.016667 (上位検算と一致) | PASS |
| 非負 ランダム 100 セル × 10⁴ | 気相 2/3 種 × 補正 z mean/upwind × float64/32 × 水の戻し increment/rebuild の 16 通りすべてクランプ 0・min ρv 0 | PASS |
| 非負 補正 z の反例 (情報) | D 比 100・組成段差 1 の 2 セル: mean は点対角 1 回で ρv −6.2e-3、upwind は +3.2e-4。どちらも固定点反復は 5000 回で収束せず | §4-1 |
| 保存 1000 更新 (float64) | 各種 ≤6.8e-9、ρE 8.2e-16、初期 0 の量 0、min ρv 4.0e-4 | PASS |
| 保存 (float32, 情報) | 丸め床 16ε: 7.2e-6 (> 1e-6) / 4ε: 5.0e-7 | §4-4 |
| エネルギー 接線投影 | 6.1e-9 (\|ΔT\| 最大 0.103 K)。有限更新 (EOS 反転) 超過 −8.1e-4 K (= 許容内)。判別: 潜熱流束なし 27 倍・エネルギーだけ潜熱 33 倍ずれ | PASS |
| エネルギー (情報) | 面の h・L を T^(n+1) で評価すると 4.4 % (J_l·∇L の 2 次項が r_g→0 の節近くで効く; 物理の項) | 「接線投影」は Tⁿ のまわりの線形化と定義した |

ハーネスで決めた解釈 (plan に書かれていないもの; 上位の確認対象):
- (i) の「D = [1,2,3]e-5」は上位検算に合わせ**二元 D (D_12, D_13, D_23) を混合平均した D_k** と解釈した。
- 「解き切った」は「BE 残差 ≤1e-7 × 反復開始時の残差 (成分ごと)」**または** 丸め床 16ε·max|Mρφ| 以下。一様な蒸気のように開始時の残差自体が丸め程度の成分があるため床が要る。
- 保存試験 (S5) の解き切りは 1D の参照用に「点対角 + 乱流の三重対角」を前処理にした (固定点は同じ)。点 Jacobi では拡散数 ~200 のエネルギー試験で収束しなかった。
- エネルギー試験の有限更新は熱伝導なし (種のエンタルピー流束だけ)。

## 4. 未確定事項 (方針の穴になりうるもの; 実装前に上位の判断が要る)

1. **補正項の面の z** (plan §4.2 は $z_k\Sigma j^0$ の面値を定めていない; `species_diffusion_d:307` の現行は算術平均)。
   算術平均は隣セル分が負の非対角になり M 行列でない。通常の組成 (D 比 ≤3) ではランダム試験で負にならないが、D 比 100・段差 1 の反例で蒸気 0 のセルが負になる。
   補正の質量流束 $-\Sigma j^0$ が出ていく側の z (風上) なら点対角 1 回は非負。H2 を含む燃焼ガス (D 比 ~3–4) で段差の大きい面がどれだけあるかは未調査。
2. **相変化ソースと更新クランプが蒸気の増分にどう入るか** (plan §4.2 の「蒸気/液を変数にした増分」は輸送だけを想定した書き方)。
   現行: ソースは `res_rog` だけ (総水分にソース無し)、液の増分全体に θ (`condensationUpdateLimiter_d.cuh:58-91`)、成分ごとの非負化 (`:92-113`)、$\phi_N\delta\rho$ 項、θ_b、bounds が掛かる。
   $R_v=R_w-R_g$ を蒸気の対角で割ると、ソース分 $S(1/D_g^{\rm eff}-1/D_v)$ が反復の途中で総水分を動かす (固定点では消えるが、θ が固定点でも効いていると消えない)。
   候補: (P1) 蒸気の増分 = 輸送分 $R_v^{\rm tr}/D_v$ − 液の増分のうちソース分 (液と同じ θ・床を通した後の値)、総水分の増分 = 輸送分だけ。
   (P2) 総水分は現行どおり総水分の増分 (ソース無し) で作り、蒸気の非負は液側の θ (蒸気枯渇) とクランプに任せる — ただしこれは codex M2 の 2 セル例で蒸気が負になる元の形。
   (P3) 輸送 (拡散・移流) とソースを分割ステップにする。どれも §4.2 の文言からは決まらない。
3. **総水分の dual-time FCT**: 化学種は FCT を通らない。末尾 FCT で液だけ反拡散補正すると蒸気が動く (移流でも既存の不整合)。
   plan §4.2 は「FCT 履歴を液・Q・総水分へ拡張」と書くが、総水分 (化学種) を受動種 FCT に入れるのか、蒸気を受動種として扱うのか、が決まっていない。
4. **float32 での §6 の許容**: (i) の許容 $10^{-6}\rho_g D|\nabla g|$ は ∇g 基準なので、面あたり Δg が小さい (細かい格子) と float32 の z の丸め
   ($\sim\varepsilon_{32}z/\Delta g$) で満たせない (dx 1 mm・Δg 1e-4 で 98 倍)。保存 1e-6/1000 更新は float32 の停止基準 (残した残差) で決まり、16ε 床では超える。
   カーネル (float32) の受け入れ試験の許容を、ハーネス (float64) の条件と別に事前固定する必要がある。
5. **g=0 のビット一致の判定法** (#5): 残差の atomicAdd 順序は run ごとに変わるので、run 全体のビット一致は分岐を入れても保証できない。0 step probe で判定するか。

## 5. 上位諮問ブリーフの下書き (diagnostician / codex `--stage diagnose` 共通)

> 保存先の想定: `notes/reviews/briefs/2026-10-02-twophase-diffusion-kernel.md` (諮問を回すのは主セッション)。

**問い**: plan `condensation-two-phase-transport` §4.2 の二相拡散カーネル (§5.1 #4) に着手してよいか。§4 の未確定 1–5 をどう決めるか。
特に (2) 相変化ソース・更新クランプと蒸気/液の増分の関係、(3) 総水分の FCT。

**観測事実** (ホスト参照実装のみ。forge の run・カーネルは無い):
- `tests/unit/test_twophase_diffusion_harness.py` (float64 参照, 2026-10-02, HEAD `33520773` + 未 commit のハーネス) が §6 の単体条件を全件 PASS (§3 の表)。
- 3 セル判別: A (点対角 1 回) 総液量 +3.3333333 % / +2.8571429 %、B (解き切った BE) +2.2e-9 (1000 更新)。上位検算 (+3.3333346 % / +2.8571441 %, float32) と一致。
- 2 セル非負: 改訂前の更新で ρv −0.016667 (上位検算と一致)、新 (蒸気/液変数) で 0。
- 補正の面 z: D 比 100・段差 1 の 2 セルで、算術平均は ρv −6.2e-3、風上は +3.2e-4 (点対角 1 回)。通常の組成のランダム試験では両者とも負にならない。
- float32: (i) の許容は dx 1 mm・Δg 1e-4/面で 98 倍超過。保存は停止の丸め床 16ε で 7.2e-6、4ε で 5.0e-7 (1000 更新)。
- エネルギー: 接線投影 6e-9、潜熱流束なし・エネルギーだけ潜熱はそれぞれ 27・33 倍ずれる。面の L を T^(n+1) で評価すると節の近くで 4.4 % (J_l·∇L)。

**期待値と出典**: plan §6 の事前固定の数値 (単体・保存・非負・エネルギー)。判別規則は diagnose 2026-09-27 (`notes/reviews/2026-09-27-condensation-transport-redesign-diagnose.md`)。

**再現条件**: `python3 solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py` (約 9 分, `--only S3` などで個別)。CFD・GPU 不要。

**実施済みの操作と結果**: methods §7c に §4.2 の仕様を記述 (未実装と明記)。ハーネスで決めた解釈は本メモ §3 末尾 (二元 D の混合平均、解き切りの定義と丸め床、1D の前処理、熱伝導なし)。
現行コードの該当箇所は本メモ §1。

**仮説** (検証していない):
- H1: 契約 (保存形 BE の固定点 + 点対角の前処理 + 蒸気/液変数) は、輸送だけなら §6 を満たす (ハーネスが支持)。
- H2: 実装上の主なリスクは輸送ではなく、(a) 相変化ソースと θ・床が蒸気の増分にどう入るか、(b) 水の commit を液の後に遅らせる順序変更 (`main.cpp:2127-2152`, `:2355-2385`) と
  coupling 2 との非互換、(c) dual-time で総水分が FCT を通らないこと。
- H3: 補正の面 z は風上にすると非負が構造で保たれるが、算術平均 (現行) からの変更は g=0 のビット一致を崩さない範囲 (TP carrier 凝縮の新カーネル内だけ) に限られる。

**呼び出し側の前提で疑わしいもの**: 「解き切った」の定義 (相対 1e-7 + 丸め床) と、float32 の受け入れ許容を float64 の条件から分けるべきか。
