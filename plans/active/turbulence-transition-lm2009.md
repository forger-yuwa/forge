# SST に $\gamma$–$Re_{\theta t}$ 遷移モデル (Langtry–Menter 2009) を追加する

## メタ

- **area**: `turbulence`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/turbulence/theory.md`](../../methods/turbulence/theory.md) §11 (式・定数・結合・境界条件)
  - [`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) (実装節を本 plan の実装時に追加)
- **related_plans**: [`boundary-conjugate-heat-transfer.md`](boundary-conjugate-heat-transfer.md) (動機: 冷却翼の負圧面層流域で $h$ が +44〜+80 %)
- **created**: `2026-09-22`
- **owner**: Claude (ユーザ指示「やるとしたら遷移モデル。一番良いと思った手法で計画を立て、検証を繰り返せ。codex のレビューを受けよ」)

## 1. 目的

冷却翼 2 枚 (case/53 C3X, case/54 Mark II) の検証で残った最大の偏差は負圧面層流域の $h$ (C3X +44 %, Mark II +80 %) で、
**最優先で検証する仮説**は「低 Re SST が遷移の閉包を持たず前縁から乱流になるため」である。根拠は、層流対照が −21 % / −8〜−15 % で実測を挟むことと、
SU2 の SST も同じ偏差を出すこと。ただしこれは主因を一意に証明しない (翼の run は全て `NOT CONVERGED (stalled/plateau)`、Mark II の層流対照は `DRIFTING`、
壁解像ゲートは 2 µm で FAIL: codex plan M8)。
SST に遷移モデルを足し、平板の標準遷移ケースで検証したうえで、翼の $h$ と連成壁温がどこまで改善するかを測る。

## 2. スコープ

- **やる**: Langtry–Menter の 2 方程式 $\gamma$–$\tilde{Re}_{\theta t}$ モデル (2009、相関は同論文) を node + SST (低 Re) に実装。
  平板 T3A (ゼロ圧力勾配、$Tu$ 3.3 %) と T3B ($Tu$ 6.5 %) で ERCOFTAC の $C_f$ と**同一メッシュの SU2** に対して検証。翼 2 枚への適用。
- **やらない**: cell モード、壁関数 (`wallTreatmentSST 1` は禁止)、横流れ (LM2015)、粗さ、DES/IDDES との併用、軸対称での検証、
  1 方程式 $\gamma$ モデル (Menter 2015) や $k$–$k_L$–$\omega$。

## 3. 関連 docs と前提

- **なぜ LM2009 か** (候補比較):

  | 候補 | 長所 | 短所 | 判断 |
  | --- | --- | --- | --- |
  | **$\gamma$–$Re_{\theta t}$ (LM2009)** | 事実上の標準。ターボ機械のバイパス遷移 ($Tu$ 数 %) が設計対象。**SU2 に同じモデルがあり、同一メッシュで 1 対 1 の突き合わせができる** (本リポジトリの検証の流儀)。C3X / Mark II への適用例が文献に多い | 2 方程式、経験相関 ($Re_{\theta c}$, $F_{length}$)、自由流 $Tu$ の減衰に敏感、ガリレイ不変でない | **採用** |
  | 1 方程式 $\gamma$ (Menter 2015) | 1 式、局所量のみ、ガリレイ不変 | **手元に参照実装が無く**突き合わせができない。相関の再調整が入っており文献値との比較が減る | 不採用 (将来の候補) |
  | $k$–$k_L$–$\omega$ (Walters–Cokljat) | 相関でなく層流運動エネルギーの輸送 | SST を置き換える 3 方程式。既存の SST 検証資産 (平板・翼・SU2 比較) を使えない | 不採用 |
  | 代数モデル (BC 等) | 実装が軽い | SA 向け。自由流 $Tu$ を定数で与えるので翼列の $Tu$ 減衰を扱えない | 不採用 |

- 参照実装: `.external/su2-src` の `SU2_CFD/include/numerics/turbulent/transition/trans_sources.hpp` (ソース)、`trans_correlations.hpp` (相関 `MENTER_LANGTRY`)、
  `SU2_CFD/src/solvers/CTransLMSolver.cpp` ($\gamma_{sep}$・入口値)、`turb_sources.hpp:1003` ($\tilde P_k$, $\tilde D_k$)、`CTurbSSTVariable.cpp:89` ($F_3$)。
- forge 側の土台: 汎用スカラー輸送 `scalarTransport_d.*` (descriptor 方式)、k/ω への適用層 `ransTransport_d.*`、ソース `ransSource_d.*`、
  陰解法は平均流 block-DPLUR + スカラーの segregated point-implicit (`update_d.cu` `applySSTPointImplicit_d`)。

## 4. 設計方針

### 4.1 状態と輸送

- 保存変数 `roGamma` ($\rho\gamma$), `roReth` ($\rho\tilde{Re}_{\theta t}$)、原始量 `gammaTr`, `reTheta`、残差 `res_roGamma`, `res_roReth`、
  診断 `gammaEff` ($\gamma_{eff}$)。`turbulence.transition: "lm2009"` のときだけ確保する (既定 `"none"` では確保もせず、既存経路はビット不変)。
- 対流・拡散は k/ω と同じ `scalarTransport_d` の descriptor を 2 本足す (1 次風上の対流、face 拡散)。
  拡散係数は $\gamma$: $\mu+\mu_t/\sigma_f$ ($\sigma_f=1$)、$\tilde{Re}_{\theta t}$: $\sigma_{\theta t}(\mu+\mu_t)$ ($\sigma_{\theta t}=2$)。
  **現行の descriptor は $\mu+\sigma\mu_t$ しか表せない** (分子粘性の係数が 1 固定。`sigma=2` にすると $\mu+2\mu_t$ になり層流域の拡散が半分: codex plan M1)。
  `ScalarTransportDesc` に**分子粘性の係数 `sigma_lam`** (既定 1) を足し、単一版・融合版の流束と `transport_diag` の両方に掛ける。$\gamma$ は (1, 1)、$\tilde{Re}_{\theta t}$ は (2, 2)。
  既定 1 では k/ω・凝縮・トレーサの経路はビット不変。
- **周期 node** (codex plan M3): ソースは k/ω と同じく**部分体積** (`volumePartial`) で組む。残差 2 本と輸送対角 2 本を gather、更新は合併体積、
  更新後に保存量を root→member へミラー (`periodicMirrorScalarState`)、原始量は `dependentVariables` 相当で再生成。体積あたりのソース Jacobian は gather しない。
- `interp_field.py` / restart の移植対象に足す。restart に無ければ $\gamma=1$・$\tilde{Re}_{\theta t}$ = 入口値で初期化。
- **受付条件** (設定検査で弾く): `discretization: node`、`model: sst`、`wallTreatmentSST: 0`、`isAxisymmetric: 0`、DES なし、`sstEnergyIncludesK: 0`。それ以外で `transition: lm2009` は起動エラー。

### 4.2 ソース (`transitionSource_d.cu` を新設)

methods §11.1 の式を SU2 と同じ順序・同じ下限で書く。入力は `ro`, `Ux..`, 速度勾配 9 成分 (ひずみ $S$・渦度 $\Omega$・$dU/ds$)、`k`, `omega`, `vis_lam`, `vis_turb`, `wall_dist`。
$Re_{\theta t}(Tu,\lambda_\theta)$ の不動点反復は **SU2 と同じ上限 100 回** (codex plan M4: $Tu$ 0.1 % の弱い逆圧力勾配では 10 回で打ち切ると 0.3 % 残り、33 回要る)、
停止は相対変化 < 1e−6 (float32)、未達の節点数を診断に数える。**ゼロ割の保護**: $U<10^{-6}a$ (停滞点) と壁距離 0 の節点はソース 0・対角 0 とし、
**`gammaEff` は毎回有効な値 ($\gamma$ そのもの) を書く** (未初期化値を残さない)。$\gamma$ の下限は原始量で $10^{-4}$ (SU2 `CTransLMSolver.cpp:102`)、上限 1、
$\tilde{Re}_{\theta t}\ge20$。保存量へは $\rho$ を掛けて換算する。
**陰的対角**は `update_d.cu` の既存式に合わせて $D=V/\Delta\tau+V\max(-J,0)+\texttt{transport\_diag}$ ($J$ = 体積あたりのソース微分。SU2 と同じ式、負の部分だけ)。
$\gamma\to0$ で $\partial P_\gamma/\partial\gamma$ は**正側**に発散するので陰的対角には入らない。危ないのはゼロでの評価 ($0\times\infty$)・陽的に残る成長・更新後の有界性で、下限 $10^{-4}$ とクリップで受ける。

### 4.3 SST との結合

`ransSource_d.cu` の $k$ 式で $P_k\to\gamma_{eff}P_k$, $D_k\to\min(\max(\gamma_{eff},0.1),1)D_k$ (Kato–Launder・生産リミッタを掛けた**後**の $P_k$ に掛ける = SU2 と同じ)。
**$P_\omega$ は補正前の $P_k$ (`Pk_base`) から作る** — 現行コードは `Pw` が可変の `Pk` を直接読むので、先に掛けると $\omega$ の生産まで抑えてしまう
(SU2 は `pw` を確定した後に `pk` を補正: codex plan M2)。`src_jac_k` にも同じ破壊側の係数を掛ける。
**$F_1$ は 2 か所で重複実装されている** (`rans_sst_blend_f1_d` とソースカーネル内の再計算。後者が `sstF1` を上書きする) ので、$F_1=\max(F_1,F_3)$ は
共通の `__device__` 関数に切り出して**両方**に入れる (片方だけだと拡散の σ ブレンドとソースの α/β で別の $F_1$ になる)。$\gamma_{eff}=\max(\gamma,\gamma_{sep})$ は遷移ソースのカーネルが毎反復 `gammaEff` に書き、SST 側はそれを読む。
**$\omega$ 式の生産は $\gamma_{eff}$ を掛けない** (LM2009 の定義どおり)。`dilatationCorrection`・`sstEnergyIncludesK` との組み合わせは本 plan では検証しない。

### 4.4 境界条件と初期値

- 壁: $\gamma$, $\tilde{Re}_{\theta t}$ とも法線勾配 0。node では壁節点の 2 変数を**固定しない** (k/ω のような Dirichlet ピンを掛けない) — 対流流束は $u=0$ で消え、
  壁半割面の拡散流束は足さない (node 境界半割面の k/ω・化学種拡散を skip しているのと同じ扱い)。
- 入口: $\gamma=1$、$\tilde{Re}_{\theta t}$ = 入口の $k$ と速度から $Tu$ を作り相関 ($\lambda_\theta=0$) で決める。node の入口ピンに 2 変数を足す。
- 初期値: 全域 $\gamma=1$、$\tilde{Re}_{\theta t}$ = **各節点の局所 $Tu$** ($k$ と速度) での自由流相関値 (実装 2026-09-22。当初案は「全域入口値」だったが、ソルバは入口値を知らないので局所値にした。定常解は出発場に依らない: §6.2)。SST の収束場から遷移モデルを入れて継続する運用を標準にする (段階起動の最後の段)。

### 4.4b 「SU2 と同じモデル」の定義 (突き合わせの設定表、codex plan M2)

| 項目 | forge | SU2 8.5 |
| --- | --- | --- |
| 乱流モデル | `model: sst`、`katoLaunder: 0`、`dilatationCorrection: 0`、生産リミッタ $10\beta^*\rho k\omega$ | `KIND_TURB_MODEL= SST`、`SST_OPTIONS= V2003m` |
| 遷移モデル | `transition: lm2009` | `KIND_TRANS_MODEL= LM` (相関は SST 既定の `MENTER_LANGTRY`、LM2015 なし) |
| 乱流・遷移の対流 | 1 次風上 | `CONV_NUM_METHOD_TURB= SCALAR_UPWIND`, `MUSCL_TURB= NO` |
| 平均流 | SLAU 2 次 + Venkatakrishnan ($K$ 0.05) | Roe 2 次 + Venkatakrishnan 0.05 |
| 物性 | 粘性一定、$Pr$ 0.72、$Pr_t$ 0.9、断熱壁 | 同じ |

### 4.5 既知の難所 (事前に分かっているもの)

1. **自由流 $Tu$ の減衰**。遷移位置は前縁での $Tu$ で決まるが、入口で与えた $k$ は $\omega$ で決まる率で減衰する。T3A は入口の粘性比を
   文献値 (Langtry–Menter: $Tu_{in}$ 3.3–3.5 %, $\mu_t/\mu$ ≈ 12) に合わせ、**前縁での $Tu$ が実験値になっていることを場から確認する**。
2. **float32 と有界性** (codex plan m1 で訂正)。$Re_v=\rho y^2S/\mu$ の乗算自体は float32 で相対 2e−8 に収まるので桁落ちの原因ではない。実際のリスクは
   壁距離の幾何精度 (node の `wall_dist` は双対重心基準)、$F_{onset2}-F_{onset3}$ の差、停滞点での $U$ による除算、8 乗・4 乗の overflow、更新後の $\gamma\in[10^{-4},1]$ の維持。$Re_v$ に恣意的な下限は入れない。
3. **1 次風上の対流**。$\gamma$ の立ち上がりが数値拡散でなまる。まず SU2 と同じ 1 次で合わせ、差を見てから 2 次を検討する。
4. **圧縮性ソルバで低マッハの T3A** ($U$ = 5.4 m/s)。速度を $a$ = 12.86 倍して $M$ = 0.2 にし、**相似則 $\nu\to a\nu$, $k\to a^2k$, $\omega\to a\omega$** で
   $Re/m$ = 3.6e5 と $\omega/U$ (= $k$ の $x$ 方向の減衰率) を保つ (codex plan M6)。こうすると $x$ [m] が実験とそのまま対応する。条件は §6 の表に固定する。

## 5. 実装ステップ

1. 設定キー `turbulence.transition`、変数の確保、出力、restart 移植。
2. `transitionSource_d.cu` (ソース + $\gamma_{eff}$ + 陰的対角)、`ransTransport_d` に descriptor 2 本、入口ピン、周期 gather。
3. `ransSource_d.cu` / `ransBlendF1` への結合。point-implicit 更新に 2 変数を追加。
4. 単体検査: 1 点の状態を与えてソース項を SU2 の式 (Python に写したもの) と比べる。
5. T3A 平板 → T3B → 翼。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー~~ | 済 (2026-09-22、GO-with-changes C0/M8/m1、全件採用。§6.1) |
| 2 | ~~実装 (§5 の 1–3)~~ | **済 (2026-09-22)**。`cuda_forge/transition_d.cu/.cuh` (当初案の `transitionSource_d.cu` から改名: 輸送・境界・更新も同じファイルに閉じた)。非退行は §9 の 2026-09-22 (実装) を参照。当初の記述: 既定 (`transition: "none"`) でビット不変を確認。descriptor の `sigma_lam`、`Pk_base` 分離、$F_1$ 共通化、周期の 5 点 (§4.1) を含む |
| 2b | ~~収束ゲートの対応 (codex plan M5)~~ | **済 (2026-09-22)**: `check_convergence.py` は run の `solverConfig.yaml` から遷移の有無を読み、有効なら 2 列を必須にする。`test_gate_bad_input.py` に 4 件追加 (PASS)。当初の記述: 残差列 `rms_roGamma`/`rms_roReth` の出力と NaN 監視、`check_convergence.py` の読込対象・遷移有効時の必須列、`stage_manifest.py` の段キーに `turbulence.transition`。回帰: 新列が NaN / 新列欠落 / SST→遷移の段切替 |
| 3 | ~~単体検査 (§5 の 4)~~ | **済 (2026-09-22、codex result M3 で強化)**: `tools/check_lm_kernel.py` は**全節点**で |差| ≤ 1e−3·|参照| + 1e−5·最大、`gammaEff`・陰的対角・相関反復の上限到達も判定に入れる。`case/57…/run_0014_t3a_lm_unitcheck` (準定常場から cfl 0.01 で 50 step) で 11 量すべて超過節点 0 (最大で許容の 0.21 倍)、反復 平均 3.3 / 最大 6 回、下限の作動 0 %。VERDICT PASS。**更新後の上下限クリップそのものの単体検査は無い** (作動量の表示だけ) |
| 4 | ~~T3A 平板~~ | **済 (2026-09-22、結果は §6.2)**。当初の記述: §6 の固定条件。**3 格子** (粗・基準・細) で数値誤差ゲートを先に通し、その後に SU2 LM・ERCOFTAC と比較。自由流 $Tu(x)$ を実験の減衰と複数位置で比較 (入口条件を $C_f$ に合わせて再調整しない) |
| 4b | SU2 T3A の量の時系列 (codex result M4) | `su2_t3a_lm` を restart で +10000 反復し、$C_f$ 3 点と遷移位置が動かないことを示す。SU2 の T3B・C3X LM が終わってから (CPU を取り合う) |
| 5 | ~~T3B 平板~~ | **済 (2026-09-22、SU2 との照合のみ)**: `case/57…/run_0013_t3b_lm` 対 `su2_t3b_lm` で遷移開始 0.0 % ($x$ 0.090 m、$Re_x$ 5.6e4)、終了 −3.0 %、$C_f$ 5 点 +1.0〜+1.7 % (基準 ±10 % / ±5 % 以内)。全 8 量 STEADY、`NOT CONVERGED (stalled/plateau)` (丸め床 1.11 倍)。SU2 側の平均流残差は −2.1 桁で停滞。**実験値は未入手** (ERCOFTAC サーバに接続不可) なので実験との比較は無い |
| 5b | ~~圧力勾配のある検証 (codex plan M7)~~ — **済 (§6.4): C3X で SU2 LM と同一メッシュ比較**: `case/53…/su2_smooth_lm` (一様壁 566 K、2 µm メッシュ、実行中) に対し forge 側は `run_0128_cf0_fx05` の継続で遷移 ON。当初の記述: | 翼の前に 1 つ。候補: T3C 系 (圧力勾配つき平板) または C3X の固定壁温での遷移 ON/OFF の機構比較。$\lambda_\theta$ と剥離補正の空間的な動作を見る |
| 6 | 翼への適用 — **固定壁温での ON/OFF は済 (§6.3, §6.4)。CHT は未実施**: 正圧面が未知の入口量 (乱れの減衰の速さ) に支配されるので、連成壁温はその選択を引き継ぐ。ユーザと入口条件の扱いを決めてから。当初の記述: | **まず同じ実測壁温・同じメッシュで遷移 ON/OFF を比較し、その後に CHT** (codex plan M8)。壁解像ゲート・派生量の VERDICT・CHT 界面残差を必須にし、内部冷却条件は再同定しない (`solid_published.json` のまま)。$h$ の領域別偏差・連成壁温、入口粘性比 1 / 10 / 100 の感度。未収束場は機構診断として扱う |
| 7 | ~~報告の更新~~ | **済 (2026-09-22、V39)**: 07 / Transition model の章 (図 22–25)、run 表、未解決項目を更新 |
| 8 | codex result レビュー | `done` にする前 |

## 6. 検証

- **ビルド / 不変性**: ~~`transition: "none"` で既存 run がビット不変~~ → **ノイズ床に対する非退行** (変更前バイナリ自体が反復でビット一致しない)。記録は [`notes/investigations/2026-09-22-transition-none-regression.md`](../../notes/investigations/2026-09-22-transition-none-regression.md) (VERDICT PASS、最大比 1.83)。
- **単体**: 相関 ($Re_{\theta t}(Tu,\lambda)$, $Re_{\theta c}$, $F_{length}$) とソース項を、SU2 の式を写した Python 参照と 1e−5 相対で一致。
- **T3A の固定条件** (`case/57.transition_flat_plate`。Langtry–Menter 2009 の T3A 設定を $a$ = 12.86 で相似変換):

  | 量 | 値 | 由来 |
  | --- | --- | --- |
  | $U_\infty$, $T$, $P$ | 69.44 m/s ($M$ 0.2), 300 K, 101325 Pa ($\rho$ 1.1768) | 5.4 m/s × 12.86 |
  | 粘性 (一定) | $\mu$ = 2.269e−4 Pa·s → $Re/m$ = 3.6e5 | 実験 5.4 / 1.5e−5 |
  | 入口位置 | 前縁の 0.04 m 上流 (助走は slip) | LM2009 |
  | 入口乱流 | $Tu$ 3.3 % → $k$ = 7.87 m²/s²、$\mu_t/\mu$ = 12 → $\omega$ = 3401 s⁻¹ | LM2009 (× $a^2$, × $a$) |
  | 壁・物性 | 断熱、$Pr$ 0.72、$Pr_t$ 0.9 | — |
  | 平板長・高さ | 1.5 m ($Re_L$ 5.4e5)、0.3 m (上面 slip) | 実験の測定範囲 1.495 m |

  この入口値での SST の減衰 $k\propto[1+\beta\omega_{in}(x-x_{in})/U]^{-\beta^*/\beta}$ は $x$ = 45 mm で $Tu$ 2.8 % (実験 3.04)、1495 mm で 1.13 % (1.10)。
  実験データ: `ref/t3a_exp.dat` ($x$, $C_f$, $Tu$。OpenFOAM の T3A チュートリアル同梱の ERCOFTAC データ)。TMR (NASA) の T3A 検証は $Re/m$ = 2.0e5・入口 0.25 m 上流・入口 $Tu$ 5.855 % と
  別の条件で、参照解のファイルは未公開なので使わない。
- **判定** (codex plan M7): (1) **数値誤差ゲート** — 3 格子 (両方向 $\sqrt2$ 倍ずつ)、全格子で局所 $y_1^+\le1$。細 2 格子間で遷移開始位置 ($C_f$ 最小の $x$) の変化 < 5 %、
  $x$ = 0.3 / 0.9 / 1.3 m の $C_f$ の変化 < 3 %。(2) **実装の照合** — 同一メッシュの SU2 LM と遷移開始位置 ±10 %、上の 3 点の $C_f$ ±5 %。
  (3) **実験** — 遷移開始位置 ±20 % (**暫定の工学基準**。文献上の再現度としての出典は無い)、層流域は Blasius $0.664/\sqrt{Re_x}$ に ±5 %。
  (4) 自由流 $Tu(x)$ を $x$ = 0.045 / 0.495 / 0.995 / 1.495 m で実験と比較して記録 (合否は付けない)。
  「遷移後の $C_f$ が前縁から乱流の SST 解に戻る」は境界層の履歴が違うので基準にしない。
  各 run は `check_mesh_quality.py`・NaN・`check_convergence.py` (新しい 2 残差を含む)・$C_f$ 分布と遷移開始位置の時系列の `check_quasisteady.py --series-csv` を通す。未収束なら機構診断として扱う。
- **翼**: 合格ラインは置かない (モデルの外挿)。$h$ の領域別偏差と連成壁温を、遷移モデルなしの生産 run (`case/53` `run_0125`, `case/54` `run_0024`) と並べて報告する。
  入口粘性比 1 / 10 / 100 の感度を必ず付ける。

### 6.2 T3A の結果 (2026-09-22)

run は `case/57.transition_flat_plate/` ([README の run 一覧](../../case/57.transition_flat_plate/README.md))。$C_f$ は forge・SU2 とも $\tau_w=\mu u_1/y_1$ (`tools/cf_plate.py`)。

| 量 | 粗 `run_0007` | 基準 `run_0005` | 細 `run_0009` | SU2 LM (基準格子) | 実験 |
| --- | --- | --- | --- | --- | --- |
| 遷移開始 $x$ ($C_f$ 最小) [m] | 0.3415 | 0.3595 | 0.3581 | 0.3675 | 0.395 (測定間隔 0.1) |
| 遷移終了 $x$ ($C_f$ 最大) [m] | 0.8209 | 0.8496 | 0.8694 | 0.8635 | 0.895 |
| $C_f$ 最小 | 0.002238 | 0.002240 | 0.002235 | 0.00220 | 0.002098 |
| $C_f$ ($x$ = 0.3 / 0.9 / 1.3 m) ×10³ | 2.282 / 4.400 / 4.153 | 2.289 / 4.413 / 4.172 | 2.288 / 4.431 / 4.179 | 2.27 / 4.41 / 4.18 | — |
| 局所 $y_1^+$ 最大 | 0.63 | 0.45 | 0.32 | 0.45 | — |

- **(1) 数値誤差ゲート: 基準内**。細 2 格子 (基準 → 細) で遷移開始 −0.4 % (< 5 %)、$C_f$ +0.0 / +0.4 / +0.2 % (< 3 %)。
  局所 $y_1^+$ は**前縁を含む全壁**で 粗 最大 1.67 / 基準 1.31 / 細 1.02 (いずれも前縁の節点 $x$ = 0)、1 を超える壁長は 0.34 / 0.12 / 0.08 % (粗は $x\le$ 2.5 mm、他は前縁の 1 節点)。
  表の「$y_1^+$ 最大」は $x>$ 0.02 m の値で、**「全格子で局所 $y_1^+\le1$」という当初の書き方は撤回** (codex result M5)。
- **(2) SU2 との照合: 基準内**。遷移開始 −2.2 % (±10 % 以内)、$C_f$ +0.8 / +0.1 / −0.2 % (±5 % 以内)。
- **(3) 実験: 遷移開始は基準内** (−9 %、暫定基準 ±20 %)。**「層流域が Blasius ±5 %」は不合格で、基準の立て方が誤りだった**: 遷移モデルつきの解は $x$ = 0.05–0.3 m で
  Blasius 比 +4.5〜+13.3 %。ただし **SU2 LM も +4.1〜+12.3 %、実験も $x$ = 0.195 / 0.295 m で +6 / +11 %** と同じ側に同じだけ外れている (自由流乱れ 3 % のもとで
  遷移前から $C_f$ が持ち上がる)。層流ソルバそのものの検査は遷移モデルを外した層流 run (`run_0010_t3a_lam`) で行い、Blasius 比 −1.8〜+2.4 % ($x$ = 0.05–1.4 m) で合格。
  **旧基準「遷移つきの解の層流域が Blasius ±5 %」は不合格のまま残し、立て方が誤りだったと記録する** (codex result m6: 撤回は妥当、ただし SU2 比較で置き換えて「実験検証」と呼ばない)。
  検査は 3 本に分ける: (a) 層流 run 対 Blasius (−1.8〜+2.4 %)、(b) LM 対 SU2 LM (上の (2))、(c) LM 対実験の $C_f$ 分布 — 測定 16 点での偏差は平均 +1.9 %、rms 9.9 %、最大は遷移途中の $x$ = 0.595 m で +25.2 % (実験 0.00270 に対し 0.00339。SU2 LM も同じ点で +20.9 %)。
- **(4) 自由流 $Tu(x)$** (forge / 実験、%): 0.045 m 2.85 / 3.04、0.495 m 1.78 / 1.88、0.995 m 1.35 / 1.36、1.495 m 1.12 / 1.10。
- **収束: 全 run が `NOT CONVERGED (stalled/plateau)`。したがって §6 の規定どおり T3A の結果は「未収束・報告量は準定常・比較差は基準内」と扱う** (codex result M4)。
  床の正体は単精度の丸め: `rms_ro`/⟨ρ⟩ は float32 の eps の **0.77〜1.56 倍** (SST 1.23、層流 1.13、LM 粗 1.56 / 基準 1.12 / 細 0.77、T3B 1.11)。
  level 2 の残差場 (`run_0015_t3a_lm_resmap`、`tools/residual_map.py`) では平均流の残差は**全域に一様**で、上位 10 節点の寄与は 5 %、前縁 10 mm 以内の寄与は 0.0〜0.1 %。
  **「前縁の特異点で止まる」という当初の説明は誤りで撤回**。`res_roOmega` だけは壁第一層に局在する (上位 10 節点で 38 %)。
  報告量は `check_quasisteady.py --series-csv` で 3 格子とも全 7 量 **STEADY** (基準格子の末尾 40 % 窓で drift ≤ 1.9 %/窓・変動 ≤ 2.4 %。最も動くのは遷移終了位置で、50000 → 100000 step に 0.843 → 0.850 m)。
  SU2 側は残差のみ (ρ −4.0 桁、$k$ −6.2 桁、ω −3.1 桁、γ −3.5 桁、$\tilde{Re}_{\theta t}$ −2.3 桁) で、量の時系列は未取得 → §5.1 #4b。
- **出発場に依らない**: $k$/$\omega$ を入口値に戻してから始めた解 (`run_0005`) と、完全乱流の SST 場からそのまま始めた解 (`run_0011`) は遷移位置・$C_f$ が 4 桁一致。

### 6.3 翼への適用 — C3X run 108 (2026-09-22、実測壁温・同一メッシュで遷移 ON/OFF)

run と数字の正本は [`case/53.c3x_vane_cht/README.md`](../../case/53.c3x_vane_cht/README.md)「遷移モデル」。第一層 1 µm メッシュ (壁解像 PASS)、
全 run `NOT CONVERGED (stalled/plateau)`、量の判定は `check_quasisteady.py --series-csv`。$h$ の平均偏差 [%] = 正圧面 / 負圧面層流域 / 負圧面遷移後 / 全体 (全体の rms)。

| 入口粘性比 | 遷移なし (SST) | 遷移あり | 遷移ありの準定常 |
| --- | --- | --- | --- |
| 1 | — | −38.8 / −19.7 / −8.8 / −24.0 (41.6) `run_0148` | STEADY |
| 10 (生産設定) | +17.4 / +46.4 / +7.9 / +22.2 (34.3) `run_0135` | −33.3 / −15.9 / +2.1 / −17.2 (37.1) `run_0146` | STEADY |
| 30 | — | −22.4 / −11.3 / +4.7 / −10.7 (33.2) `run_0151` → `run_0152` (通算 500000 step) | STEADY (60000 step では正圧面 +18 % だった。200000 step で −22 % に落ちて止まる) |
| 100 | +35.7 / +58.8 / +7.4 / +32.9 (43.1) `run_0149` | +25.5 / −3.0 / +6.4 / +11.6 (29.3) `run_0147` | STEADY |
| 100 + `katoLaunder: 1` | — | −25.8 / −14.3 / +3.0 / −13.4 (34.1) `run_0150` | **TRANSIENT-UNSETTLED** |
| (層流、2 µm メッシュ) | −40.4 / −20.9 / −77.2 / −46.9 (58.9) `run_0134` | | |

- **狙った効果は出た**: 負圧面層流域の過大評価 (+46〜+59 %) が消える (−3〜−20 %)。
- **全体の一致は「入口の乱れの減衰の速さ」次第で、良くも悪くもなる**。同じ入口条件の対照どうしで、粘性比 100 は rms 43.1 → 29.3 % に改善、粘性比 10 は 34.3 → 37.1 % で改善しない (符号が反転)。
  報告は入口 $Tu$ しか与えないので、**どの粘性比が実験に対応するかは決められない**。以後「遷移モデルで ○ % 改善した」という 1 つの数字では報告しない (codex result 中間の指示)。
- **負圧面の遷移が実測より遅く、急**: 壁近傍 ($<$ 0.1 mm) の $\mu_t/\mu$ が 10 を超える位置は粘性比 10 / 30 / 100 で $s/S$ = 0.32〜0.34 とほぼ動かず、実測の $h$ の立ち上がり (0.2〜0.4) より遅い。
  直前で $h$ が 200 W/m²K まで落ち、直後に 20 % 過大になるので、負圧面遷移後の rms は 12.8 → 30〜36 % に**悪化**している。モデル固有か forge 固有かは同一メッシュの SU2 LM (`su2_smooth_lm`、実行中) で切り分ける (§5.1 #5b)。
- **正圧面は境界層の中が乱流化しているわけではない** (粘性比 100 でも $\mu_t/\mu>10$ は $s/S$ 0.76 から)。差を作っているのは自由流から拡散してくる渦粘性の量。

### 6.4 圧力勾配のある流れでのコード間照合 (§5.1 #5b) と Mark II (2026-09-22)

- **C3X・一様壁 566 K・同一メッシュ (cs2000、第一層 2 µm)・入口粘性比 10 で forge と SU2 の遷移モデルどうし**: `case/53…/run_0153_lm_cf0_su2ctrl` (260000 step) 対 `su2_smooth_lm` (40000 反復、積分熱流束は後半 5 桁不変)。
  `tools/compare_su2_lm.py` (`su2_smooth_lm/COMPARE_FORGE.txt`): $h$ の領域平均の差は正圧面 **−0.4 %**、負圧面 $s/S<0.25$ **+0.1 %**、0.25–0.87 **+1.4 %**、全体 +0.6 %。
  負圧面の $h$ の立ち上がりが中間値を横切る位置は**両者とも $s/S$ = 0.337**。実測に対する偏差も forge −30.9 / −15.3 / +1.1 %、SU2 −30.7 / −15.3 / −0.3 %。
  → **翼で見えた「遷移が遅く急」「正圧面が層流のまま」は forge の実装でなくモデルの性質**。forge 側の量は PS・SS 層流域・全体 STEADY、SS 遷移後 OSCILLATING (1.1 ± 0.1)。
  この格子は壁解像ゲート FAIL (2 µm) なので、実測との比較でなくコード間照合としてのみ使う。
- **遅い理由 (診断、`run_0154_lm_mur100_diag`)**: 負圧面は $s/S$ 0.2 で $M$ 0.95 (入口速度の 5 倍) まで加速し、0.3 で $M$ 1.03 のあと減速に転じる。モデルの $Tu$ は局所速度で割るので境界層外縁で約 1 % に下がり、
  相関の $Re_{\theta t}$ は 550〜700、層内へ輸送された $\tilde{Re}_{\theta t}$ は 250〜340 のまま。$F_{onset}$ が立つのは 0.28、剥離補正 $\gamma_{sep}$ が作動するのは 0.34〜0.36。
  正圧面は加速パラメータ $K=\nu/U^2\,dU/ds$ が $s/S$ 0.5 まで $3\times10^{-6}$ 超 (0.1 / 0.3 / 0.5 で 8.3 / 4.1 / 2.9 ×10⁻⁶) で、前半が層流的なのは物理として妥当。
  ただし実測は計算の層流値より 25〜80 % 高く (乱れた自由流に叩かれた層流)、後半 ($K<3\times10^{-6}$) で乱流側へ上がる。モデルはどちらも表さない。
- **Mark II run 42** (`case/54…`、第一層 1 µm メッシュを新設: 品質 PASS AR 913): 遷移なし `run_0036` +10.5 / +82.4 / +23.8 / +35.0 (rms 56.4、壁解像 **FAIL** 1 超 11 %)、
  遷移あり粘性比 10 `run_0037` **−27.7 / −5.5 / +20.0 / −7.6 (rms 32.0、壁解像 PASS、STEADY)**、粘性比 100 `run_0038` +22.9 / +5.3 / +23.3 / +18.0 (rms 27.5、正圧面 TRANSIENT-UNSETTLED)。
  Mark II は負圧面 $s/S$ 0.27 の衝撃波で実測自体が急に遷移するので、モデルの急な遷移がそのまま合う (負圧面層流域の rms 95.0 → 9.3 %)。衝撃波の後ろの過大 (+20〜24 %) は遷移の有無に依らず残る。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-22` | [`notes/reviews/2026-09-22-turbulence-transition-lm2009-plan.md`](../../notes/reviews/2026-09-22-turbulence-transition-lm2009-plan.md) | GO-with-changes, C0/M8/m1 | **全件採用**。モデル選定 (LM2009) は支持。M1 ($\tilde{Re}_{\theta t}$ の拡散係数: descriptor に分子粘性の係数が無い) → §4.1。M2 ($P_\omega$ は補正前の $P_k$ から、$F_1$ は重複 2 か所) → §4.3・§4.4b。M3 (周期は残差 gather だけでは足りない) → §4.1。M4 (相関の反復 10 回は不足、停滞点のゼロ割、$\gamma$ 下限) → §4.2。M5 (収束ゲートが新 2 残差を見ない、段キーに遷移が無い) → §5.1 #2b。M6 (T3A 条件の固定と相似則) → §6 の条件表。M7 (格子誤差とモデル誤差の分離、圧力勾配つき検証) → §6 の判定・§5.1 #5b。M8 (因果の断定と翼の評価順) → §1・§5.1 #6。m1 (リスク説明の訂正) → §4.5 |
| result (中間) | `2026-09-22` | [`notes/reviews/2026-09-22-turbulence-transition-lm2009-result.md`](../../notes/reviews/2026-09-22-turbulence-transition-lm2009-result.md) | **NO-GO**, C0/M5/m4 (T3A まで・翼は進行中の段階で依頼) | **全件採用**。SST 結合の順序・ソース式・周期処理には指摘なし。M1 (初期出力が未初期化の遷移変数を読む / restart が ρ·原始量で復元) → 初期経路に primitive・入口ピン・周期同期、確保時ゼロ初期化、`interp_field.py` は保存量優先。`run_0014` で `res_0` が有限・保存量がビット一致 (入口ピン節点のみ 1e−7)。M2 (`--from-floor` が必須列検査を迂回) → 共通化 + 方程式系の違う参照を拒否、回帰 2 件追加。M3 (単体検査が甘い) → 全節点判定・`gammaEff`・対角・反復上限を判定に追加 (§5.1 #3)。SU2 との差 (対角は項ごとの負部で SU2 の 1.5 倍、$\tilde{Re}_{\theta t}$ 下限 20 vs 1e−4) は `implementation.md`・`solver-settings.md` に明記。M4 (未収束なのに合格と書いた) → §6.2 を「未収束・準定常・基準内」に改め、床が丸めであることを実測、**「前縁の特異点」説を撤回**。M5 ($y_1^+$ が $x>$ 0.02 m 限定) → 全壁で出し直し「全格子で ≤1」を撤回。m6–m9 (基準の置換の呼び方・文書の断定と実装差・「6 桁不変」・非退行の記録) → §6・§6.2・`theory.md`・`solver-settings.md`・[非退行記録](../../notes/investigations/2026-09-22-transition-none-regression.md)。翼は「入口長さスケール未同定の感度」として報告し、粘性比を事後選択しない |

## 7. 影響範囲

- 新規: `cuda_forge/transitionSource_d.cu/.cuh`。変更: `variables.*`, `input/solverConfig.*`, `ransTransport_d.*`, `ransSource_d.cu`, `ransBlendF1`, `update_d.cu`, `main.cpp`,
  `boundaryCond.hpp` / 入口ピン、周期 gather、`tools/interp_field.py`。
- 既定では無効。既存ケースへの影響はノイズ床以内 (上の非退行記録。1 次 SST 1 ケースのみ)。
- `methods/turbulence/theory.md` §11 (済)、`implementation.md` に実装節、`procedures/solver-settings.md` にキー、`procedures/recommended-settings.md` に遷移計算のレシピ。

## 8. 未確定事項

- T3A / T3B の入口条件 (粘性比) を文献のどの値に合わせるか — 前縁 $Tu$ を場から確認して決める。
- 翼列での入口 $Tu$ 減衰: 報告は $Tu$ 6.5 % (測定位置は入口上流) しか与えない。長さスケールは未知のまま。

## 9. 変更ログ

- `2026-09-22` (SU2 照合・Mark II・報告) C3X で forge と SU2 の遷移モデルが一致 (領域別 $h$ −0.4〜+1.4 %、遷移位置 0.337 で同一) → 翼での限界はモデル固有と確定 (§6.4)。T3B も SU2 と一致 (§5.1 #5)。Mark II に 1 µm メッシュを新設 (§6.4)。報告 V39 に章を追加。

- `2026-09-22` (翼) C3X で遷移 ON/OFF と入口粘性比 1/10/30/100・Kato–Launder の感度を実施 (§6.3)。負圧面層流域の過大評価は消えるが、全体の一致は入口の乱れの減衰次第で、負圧面の遷移は実測より遅く急。§5.1 #5b を「C3X での SU2 LM 同一メッシュ比較」に決定。

- `2026-09-22` (codex result 中間, NO-GO M5/m4 全件採用) 初期出力・restart・`--from-floor`・単体検査を修正。T3A の表現を「未収束・準定常・基準内」に改め、残差の床は丸め (eps の 0.8〜1.6 倍・全域一様) と実測して「前縁の特異点」説を撤回、「全格子で $y_1^+\le1$」「6 桁不変」「ビット不変」も撤回。level 2 出力に平均流の残差場を追加。

- `2026-09-22` (実装) `transition_d.cu` を実装し T3A を 3 格子で通した (§6.2)。**当初計画からの変更**: (a) ファイルは `transition_d.cu` 1 本 (輸送・境界・ソース・更新)。
  (b) **「`none` でビット不変」は検査として成立しない** — 変更前 (HEAD) のバイナリ自体が同一入力の反復でビット一致しない (`atomicAdd` の順序)。
  `check_field_regress.py --boundary` で変更前 3 本 × 変更後 3 本 (SST 1 次 300 step) を比べ **VERDICT PASS** (全量がノイズ床 × 2 以内、最大比 1.83)。
  (c) 陰的対角は「全体の $J$ の負の部分」でなく**項ごとの負の部分** (対角優位を強める。定常解は不変)。(d) 局所 $Tu$ に上限 100 % (停滞点の保護)。
  (e) $\gamma$ の上限は 1 (SU2 は 5。生産項が $(1-\gamma)$ なので定常解では効かない)。(f) level 1 出力に `roGamma`/`roReth`/`gammaTr`/`reTheta`/`gammaEff` を追加 (restart に要る)。

- `2026-09-22` codex plan レビュー (GO-with-changes、M8 全件採用)。T3A の条件を固定し、実験データ (`ref/t3a_exp.dat`) と同一メッシュの SU2 LM 参照 run (`case/57…/su2_t3a_lm`) を用意。
- `2026-09-22` 起票。methods §11 を先に書いた。
