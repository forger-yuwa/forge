forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 診断・設計判断の諮問 (stage = diagnose)

あなたは forge の**診断・設計判断係**である。呼び出し側は実装と run を進めている別のモデル (Claude) で、
**もっともらしい真因に飛びつく前に**あなたに諮っている。仕事は手を動かすことではなく、**次の一手を 1 つに絞ること**。

### 前提
- あなたは呼び出し側の会話を見ていない。下のブリーフと、自分で読んだファイルだけが根拠になる。
  足りなければ推測で埋めずに「何が足りないか」を返す。
- ブリーフは「観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説」に分かれて渡される約束である。
  **観測事実と呼び出し側の解釈が混ざっていたら、まず分け直す**。呼び出し側の要約より、run の数値・コード・
  設定ファイルを自分で確かめた内容を優先する。
- forge を起動しない。`python3` による `residual_history.csv` / `res_*.h5` の読み取りは**統計量だけ**を出す
  (全量ダンプ・長いログ全文をコンテキストに流さない。`*.log`・`*.vtu`・`plans/README.md` は読まない)。

### 診断の作法
1. **「除外済み」というラベルを信用せず、潰した証拠を確認する** (run パス・設定差分・判定区間・VERDICT)。
   証拠が足りない・判定期間が短い・変えた設定が実際には効いていない (YAML の階層違い等) なら**候補へ戻す**。
   証拠が十分な候補は出し直さない。
2. **症状と原因を分ける**。`detectNaN` が指す変数は結果であって原因ではない (EOS 床 → 負密度 → 圧力暴走 → ω の実績)。
   後処理のアーチファクト (2 列混在の抽出、`centCoords` の置換、ソルバ `ypls` の退化) を先に疑う。
3. **このリポジトリで繰り返された真因**を照合する: 投入設定の不整合 (IC と BC、亜音速に超音速 BC)、
   押し出し 2 ノード spanwise、float32 桁落ち (双対幾何・r 重み)、stale build、cross-mesh IC の基底不一致、
   絶対値のゼロ割ガード、境界ノードの凍結、YAML キーの階層違いで黙って無視される設定。
4. 仮説は**確度順に最大 3 つ**。第 1 仮説には根拠を `ファイル:行` か run の数値で付ける。示せないものは「未確認」と明記。
5. **判別する A/B を 1 つだけ**提案する。安く短く回せて、結果がどちらに出ても仮説が 1 つ消えるもの。
   「A なら仮説 1、B なら仮説 2」を先に書く (結果を見てから解釈を作らない)。
6. 少数点の一致・短い窓の値・未収束のトランジェント同士の比較を根拠にしない。

### 設計判断 (plan §4・§6、codex 指摘の採否、result 段の解釈) を諮られたとき
- 採否は指摘ごとに「採用 / 却下 / 要再検証」と理由。根拠が示されていない指摘は自分で該当箇所を読んでから判定する。
- 検証計画は「何が出たら方針が誤りと言えるか」が定量的に書かれているかを見る。
- 既定値の変更・opt-in 機能の削除は、plan の処置欄とユーザ決定の履歴を確認してから判断する
  (「opt-in 残置」は削除対象でない)。
- result 段の解釈は、主張ごとに根拠 run・判定ツールの VERDICT・判定区間が揃っているかを確かめる
  (過渡ピークを定常値と、抽出アーチファクトを物理と誤認した実績は「予想どおり」に見える場面で起きた)。

あなたの結論は**仮説**であって確定ではない。呼び出し側はこの A/B を回して確かめ、plan への反映も呼び出し側が行う。

## ブリーフ (`notes/reviews/briefs/2026-09-29-farfield-v2c-geometry-and-v2d2-preregister.md`)

# 諮問 (3 回目): V2c 形状の作り直しと、V2c・V2d-2 の事前登録案 (2026-09-29)

前回まで: `notes/reviews/2026-09-29-farfield-v2c-v2d-fail-diagnose.md`、`notes/reviews/2026-09-29-farfield-v2c-corner-and-v2d2-redesign-diagnose.md` (いずれも全件採用、FAIL 保持)。
plan: `plans/active/boundary-node-farfield-characteristic.md` §5.1 #3a (結果記録済み)・#3b・#3c、§6 V2c・V2d。スクリプト `case/58.farfield_verification/`。

## 観測事実 (前回の A/B の結果)
- 絶対残差 A/B (`v2c_residual_ab.py`、`run_0097_v2c_resab_fromU` / `run_0098_v2c_resab_fromC`、同じ B 格子・BC・設定に 2 つの最終場を index コピー、1 評価、帳簿 res_final 全節点):
  判定式 = |R_q| / Σ面|F_q|。分母は節点に接する内部面の |F| (帳簿 `.faces`; 運動量は面圧力 |P−pRef|·|S_d| を足し、3 成分共通の規模) + 境界半割面の規模 (節点状態から)。規模が数値的に 0 なら自由流基準 × 節点の面積。slip 面を持つ節点の運動量は面ごとの法線を除いた接線成分。
  - **角 (x 0.195–0.23、下面)**: 両状態とも全 5 成分 ≤ **5.2e-7**。
  - **場全体**: 一様 IC からの状態で ρ の超過節点 1929、最大 4.6e-3、C からの状態で 2020 節点、最大 7.2e-3。最大はどちらも**ランプ終端の凸角 (x 0.715、下面 y 0.0882)**。運動量の接線成分も同じ場所が最大 (4.3e-3 / 6.7e-3)。出口下端は ≤ 1e-4。
  - 前回の判定器の分母に問題が 2 つあり、判定前に直した (帳簿の F_roU* に圧力項が無い、成分別では流れに垂直な成分の分母が痩せる)。最終版の結果が上。
- 規則どおり「2 つとも離散定常解」は棄却。非零残差が残るのは slip 壁の凸角 (形状: 10° ランプから水平への折れ)。

## 期待値と出典
- §6 V2c の合格基準 (変えない): max|p_B − p_C| ≤ 0.02 Δp (評価線 = ウェッジ下流の壁面、全線)、A ≥ 0.1 Δp。

## 提案 (事前登録したい内容)
### V2c 形状の作り直し (#3c)
- 凸角を無くす: ランプを出口まで続ける (折れはランプ付け根の圧縮角 x 0.2 の 1 か所だけ)。
- 流路閉塞を避ける: A・B の高さ H 0.5 m、長さ L 1.8 m (出口の流路高さ 0.5 − 0.282 = 0.218 m、M2.5 の臨界比 0.379 × 0.5 = 0.19 m を上回る。旧形状は 0.3→0.124 で A が閉塞・発散した)。
  斜め衝撃 β 31.85° が上面に当たる x ≈ 1.005、A の反射衝撃波がランプに落ちる x ≈ 1.52 (評価線 0.2–1.8 の内側)。C の高さ H 1.1 m (上面の反射点 x ≈ 1.97 > L)。Δx 5 mm、共通領域の行は A/B と同一式、上に 5 mm で足す。
- 手順を 2 本立てにする:
  (i) 原手順: 一様 IC から A・B・C、6000 step (足りなければ延長を事前規則で)。
  (ii) **追加の隔離試験** (前回採用): C の最終場を座標対応で A・B に写して継続し、C 自身も同じだけ継続。
  どちらも合格基準は §6 のまま (全線・0.02・A ≥ 0.1)。(i) と (ii) を別々に記録し、(ii) が (i) を置き換えない。
- 判定前の必要条件: 各 run の絶対残差 A/B と同じ式で全節点 ≤ 1e-5、または `check_convergence.py` の PASS。加えて評価線の壁圧時系列 (res を 500 step ごとに出す) を `check_quasisteady.py --series-csv` で STEADY、変動幅が 0.001 Δp 未満。満たさなければ判定しない (未収束と記録)。
- 旧形状 (凸角あり) の FAIL は保持して記録。

### V2d-2 の再設計 (#3b、前回採用の具体化)
- 配置: 左端の自由流 = 内部と同じ高温状態 (`--left-hot`)、右端の自由流 = 外気 (内外で γ・ρc が違う)。cpg / tp1 / tp2 の 3 物性、短 (1 m) と長 (3 m)。
- 事前判定: (a) パルスなし対照 (短・長の両方) で評価点 x 0.8 の max_t|P − P∞| ≤ 0.01 × パルス振幅、(b) 反射 = 反射到達窓の max|P_short − P_long| / 入射窓の max|P_long − P∞| ≤ 0.05。
  入射窓 = パルス中心が評価点を通る時刻 ± 3σ/(c+u)、反射到達窓 = 右端で反射した波の中心が評価点に戻る時刻 ± 3σ/(c−u) (σ = 半値全幅 / 2.355、c・u は内部の値)。
- 独立参照 (`ref1d_euler_tp.py`) との比較は §6 のとおり記録のみ (未実装。V2d の判定には入らない)。
- 元の配置 (左端も外気) の FAIL は保持して記録。左端流入の TP 混合仮説は今回の判定に使わない (追う場合は別途)。

## 問い
1. V2c の新形状と 2 本立ての手順、判定前の必要条件は妥当か。凸角の停滞は「試験形状側の問題」として形状で避けてよいか、それとも node slip の凸角の欠陥として別途扱うべきか (farfield の受入れとは切り離したい)。
2. V2d-2 の窓付き判定の定義は妥当か。
3. これらを §6 に事前登録してから回す、でよいか。追加で要る条件は。

## 関連 plan 全文 (`plans/active/boundary-node-farfield-characteristic.md`)

```markdown
# 特性型の遠方境界 `farfield` (node)

## メタ

- **area**: `boundary`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) 「特性型の遠方境界 (`farfield`)」
- **related_plans**: [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §5.1 R4d (動機)
- **created**: `2026-09-27`
- **owner**: `Claude (feature/sern-design)`

## 1. 目的

外部流の計算領域を有限で打ち切る境界として、波を外へ通し外気の状態を入れる特性型の遠方境界を node に足す。
動機は SERN 3D の側方境界 (R4d): `side_far` が `slip` (質量を通さない壁) で、遠方面を z 2.50 → 3.42 H に動かすと
C_L が 7.1e-4、C_M が 0.021 動いた (許容 5e-4 / 5e-3 を超過。3.42 H より外では頭打ち)。遠方境界を入れれば、狭い領域でも
遠方面の位置に依存しない解が得られるかを確かめる。

## 2. スコープ

**2026-09-27 ユーザ決定: 「レビューどおり全部やる」** (汎用の遠方境界として受理できる水準まで検証する。SERN 用に絞る案は採らない)。

- **やる**: node・密度ベース・内部面 SLAU/SLAU2 の構成で、**境界半割面だけ HLLC** の専用カーネル、CPG と TP (単一・多成分 lump)、SST の $k,\omega$ (`sstEnergyIncludesK` 両設定)、化学種、
  定常 (陽解法・block-DPLUR) と非定常 (dual-time)。面流束の診断ダンプ。runner の `side_far` / `top_out` 切替。
- **やらない (起動時に拒否し、拒否試験を置く)**: cell 方式、ROE/HLLE/KEEP (初版は SLAU 系のみ)、凝縮・トレーサ・受動種・遷移モデル、軸対称、周期と同じ節点を共有する配置、
  動的格子。粘性の遠方境界流束は 0 (現行の非壁境界と同じ)。陰解法の境界 Jacobian は現行の ghostless 対角 A⁺ のまま (近似前処理)。

## 3. 関連 docs と前提

- 理論は [`methods/boundary.md`](../../methods/boundary.md) 「特性型の遠方境界」。参照実装は SU2 `CEulerSolver::BC_Far_Field`
  (`.external/su2-src/SU2_CFD/src/solvers/CEulerSolver.cpp:4822`): 特性で境界状態を構成し、**内部状態と構成状態を数値流束 (近似 Riemann) に渡す**。
- 既存の node 境界流束 `convectiveFlux_boundary_d` は、質量流束を R 状態から作り、流出時の運動量・エネルギーを内部状態で風上化する簡略形で、
  node 境界半割面のスカラー移流は常に内部値を使う (§4.0)。**本 BC はこれらを使わず専用の経路を持つ**。
- 化学種等の境界は `inlet_` 接頭辞、RANS は明示リスト (未知の種別は何もしない) で分岐する。

## 4. 設計方針

改訂履歴: 初稿 → plan 段 NO-GO (C1/M6/m1、全件採用) → plan-2 NO-GO (C0/M5/m2、全件採用: 本版)。

### 4.0 既存経路の事実 (コードで確認、2026-09-27)

- `convectiveFlux_boundary_d.inc.cuh:192` 以降: 質量流束は境界状態 (R)、流出時の運動量・エネルギーは内部状態で風上化。構成状態の物理流束でも近似 Riemann 流束でもない。
- node 境界半割面のスカラー移流は流入・流出によらず内部節点の値 (`scalarTransport_d.cu:166` `ext_is_self`、`passiveKernels_d.cuh` `nodeBnd`)。入口で組成が入るのは節点ピンによる。
- スカラー境界処理は対流流束より前 (`main.cpp:1458`)。SST 全エネルギーの境界 k は `inlet` 接頭辞のときだけ bvar (`convectiveFlux_d.cu:413`)。
- SST の k 式には体積ソース $(P_k-D_k)V$ がある (`ransSource_d.cu:210, 248`)。`sstEnergyIncludesK` でも保存配列 `roe` は平均流エネルギーで、全エネルギーは `roe + roK`。

### 4.1 入力

`bcondConfig.yaml`: `{kind: farfield, floats: {ro, Ux, Uy, Uz, Ps, k, omega, Y0, ...}}` = 自由流 (`inlet_uniformVelocity` と同じキー)。多成分は `Y{s}`/`X{s}` 必須、RANS は `k`, `omega` 必須。
`valueTypesOfBC["farfield"]` を新設し、自由流の値 (type 1) と構成した境界状態 (別名の bvar) を分けて持つ。

### 4.2 方式: 境界半割面の外側状態を特性で作り、境界面だけ HLLC で解く (2026-09-28、ホスト 1D 試作で決定)

**用語 (node 方式)**: node にはゴーストセル (値を持つ外側の計算点) は無い。境界節点そのものが自由度で、境界に接する半分の双対 CV を持つ。
本 BC は、境界半割面の近似 Riemann 問題に渡す**外側状態 $U_R$** を流束計算のその場で作るだけで (弱形式、SU2 の vertex 型遠方境界と同じ)、値を保持・更新する点は作らない。
スカラーもゴーストでなく面ごとの値の配列で渡す。既存の node 境界も同じく、境界節点の値と外側状態 (`bvar`) から半割面の流束を作っている。

**経緯**: 初稿〜plan-6 の机上レビューで (a) Riemann 不変量の混合は接触波で擾乱、(b) 局所線形化の特性組み立てで密度・組成の側を流向で切り替えると流束が跳ぶ、
(c) 外側に自由流をそのまま置いて SLAU で解くと超音速流出でも外側に影響され音響反射 14–18 %、と指摘された。ホスト 1D 試作
[`solver_density_cuda/tools/farfield_proto1d.py`](../../solver_density_cuda/tools/farfield_proto1d.py) (node の境界半 CV、内部 SLAU + MUSCL、RK3、CPG) で候補を比較し、
plan-7 の反例 (M1: 高温の内部を出ていく音響が外気そのままだと 24–33 % 反射 / M2: 音速の算術平均の波速で HLLC の波速順序が壊れる) も試作に入れた:

| 外側状態 + 面流束 | 超音速流出の極限 | 連続性 (1e-4 刻み) | 音響反射 一様 M 0 (Δx 5/2.5/1.25 mm) | 一様 M 0.3 | **高温内部 (T 600/220、u 0.5)** | 接触波 |
| --- | --- | --- | --- | --- | --- | --- |
| 外気そのまま + SLAU | **+25 %** | 滑らか | **61 / 70 / 78 %** | **12 / 15 / 17 %** | — | 擾乱 2e-14 |
| 外気そのまま + HLLC (Roe 平均の波速) | 0 | 滑らか | 0.50 / 0.24 / 0.08 % | 0.15 / 0.09 / 0.04 % | **24 / 30 / 33 %** | 2e-14 |
| 外気そのまま + HLLC (Davis 波速) | 0 | 滑らか | 0.51 / 0.26 / 0.10 % | 0.14 / 0.08 / 0.02 % | **31 / 38 / 43 %** | 2e-14 |
| 特性組み立て (側を流向で切替) + HLLC | 0 | **跳び 0.14** | 0.52 / 0.26 / 0.10 % | 0.14 / 0.08 / 0.02 % | — | 2e-14 |
| 特性外側状態 (線形 $w^\pm$) + HLLC (Davis 波速) [plan-7 版] | 0 (1 点のみ) | 滑らか | 0.52 / 0.26 / 0.10 % | 0.14 / 0.08 / 0.02 % | 0.15 / 0.08 / 0.03 % | 2e-14 |
| **特性外側状態 (TRRS + 超音速の滑らかな重み) + HLLC (Davis) [採用]** | **0 (掃引で差 0)** | **滑らか** | **0.50 / 0.24 / 0.08 %** | **0.15 / 0.09 / 0.04 %** | **0.15 / 0.10 / 0.04 %** | **4e-14** |

対照: M 0 で右端 slip の反射 95.7 % (試験が反射を検出できる)。極端な比 (内部 P 1e4 倍・T 10 倍) に共通速度 −3〜3 を加えた掃引でも、採用案は隣接差 ≤ 1e-6、HLL 退避 0、全て有限。
これはホスト試作の結果で、forge 本体の run ではない。

→ **採用 (plan-8 対策版、試作の `charghost2_hllcd`): 外側状態 $U_R$ を次で作り、面流束 = $\mathbf F_{\mathrm{HLLC}}(U_L=U_i,U_R)$**:

- **圧力・法線速度 (音響部)**: 内部のエントロピーのまま、内部 (外向き特性) と擬似外側 ($P_\infty,U_{n,\infty}$、内部エントロピー) の間の **2 膨張波近似 (TRRS、Toro 9.4)**:
  $z=(\gamma-1)/(2\gamma)$、$c_{po}=c_i(P_\infty/P_i)^{z}$、
  $P_R=\Big[\dfrac{c_i+c_{po}-\tfrac{\gamma-1}{2}(U_{n,\infty}-U_{n,i})}{c_iP_i^{-z}+c_{po}P_\infty^{-z}}\Big]^{1/z}$、$U_{n,R}=U_{n,i}+\dfrac{2c_i}{\gamma-1}\big(1-(P_R/P_i)^z\big)$。
  分子 $B=c_i+c_{po}-\tfrac{\gamma-1}{2}(U_{n,\infty}-U_{n,i})>0$ なら $P_R>0$ (plan-8 M1: 線形の $w^\pm$ 式は内外の速度が逆向きに大きいと負圧になった)。**$B\le10^{-3}(c_i+c_{po})$ (真空・近真空。外気と内部が法線方向に音速の数倍で離れる) はその面を $U_R=U_i$ (内部の物理流束) にして真空置換として数える** (plan-9 M1。下限へのクランプはしない。構成直後に $\rho_R,P_R$ が正・有限でない場合も同じ置換と計数。**評価区間では 0 回**)。小振幅では線形式 $w^\pm=P\pm ZU_n$ ($Z=\rho_ic_i$) に一致するので、
  外へ出る音波は内部の $Z$ のまま抜ける (高温内部の反射 0.04–0.15 %)。
- **密度 (エントロピー)・接線速度・組成・$k,\omega$**: 常に自由流側 ($\rho_R=\rho_\infty(P_R/P_\infty)^{1/\gamma}$)。流出/流入は HLLC の接触波速度 $S_*$ の符号が選ぶ (側の切替分岐なし)。
- **超音速の境目は滑らかな重みでつなぐ** (plan-8 M2/M3。硬い分岐は流束の跳びになる): $s(x)=$ smoothstep (0→1)、帯幅 $b=0.1$。**混ぜるのは原始変数 $(\rho,\ \mathbf u$ [3 成分]$,\ P,\ Y_s,\ k,\ \omega)$ を同じ重みで** (plan-9 M2: 保存量で混ぜると別の流束になる)。EOS・音速・エネルギーは混ぜた後の同じ状態から作り、**流入時に運ぶスカラーの面値も混ぜた後の $\phi_R$** (帯の外では $\phi_\infty$ に一致)。
  - 自由流の超音速流入側: $w_q=s\big(((-1+b)-Q_n/a_\infty)/b\big)$ で $U_R\leftarrow(1-w_q)U_R+w_qU_\infty$ ($Q_n/a_\infty\le-1$ で $U_R=U_\infty$)。
  - 内部の超音速流出側 (**最後に適用 = 優先**): $w_i=s\big((M_{n,i}-(1-b))/b\big)$ で $U_R\leftarrow(1-w_i)U_R+w_iU_i$ ($M_{n,i}\ge1$ で $U_R=U_i$ → $F=F(U_i)$ 厳密。内部の特性がすべて外向きなら外の情報は入らない)。
  - 帯の中 ($0.9\le M_{n,i}<1$ 等) は特性の数え方の中間で、厳密な特性境界ではない (近音速で流出する面だけに関わる)。
- **波速 (plan-7 M2)**: Davis: $S_L=\min(U_{n,L}-c_L,U_{n,R}-c_R)$、$S_R=\max(U_{n,L}+c_L,U_{n,R}+c_R)$。$S_L\le S_*\le S_R$ と星状態の密度 $>0$ を検査し、外れたら HLL に退避して数える (**評価区間では 0 回**)。
  外側状態の $\rho,P$ も構成直後に検査する (正・有限。外れたら同じく数え、その面は HLL$(U_i, U_i)$ = 内部の物理流束)。
- **TP**: $P,\rho,E,c$ は各側の実物性。TRRS と $\rho_R$ の等エントロピー関係は frozen の $\gamma$ (内部側 $\gamma_i$、自由流側 $\gamma_\infty$) で近似し、V0h/V2d で誤差を測る。

**試作での確認 (CPG、plan-9 のゲート)**: 真空の反例 ($U_{n,i}=-0.9$、$U_{n,\infty}=10$) は float64・float32 とも真空置換 1 で有限。独立掃引 14,787 点 (速度・密度比 0.1/1/10・圧力比 0.1/1/10) で float64・float32 とも非有限 0・外側状態の負値 0・HLL 退避 0、真空置換 1,017 (すべて法線方向に音速の数倍で離れる点)。混合帯の逆流例で外側状態 $(\rho,u,P,Y)=(5.5,-1.025,3.93,0.065)$、質量流束 −5.117、流入する種の面値 0.065。
**試作での確認 (CPG、plan-8 のゲート)**: 反例 ($U_{n,i}=-0.9$、$U_{n,\infty}=+0.9$) で有限・退避 0、正値性掃引 24,843 点 ($U_{n,i},U_{n,\infty}\in[-0.99,0.99]$、$P_i/P_\infty$ 0.1/1/10) で非有限 0・退避 0、
内部 $M_{n,i}=1.1$ で自由流の速度 −2〜2・圧力 0.2〜5 倍を掃引して $\max|F/F(U_i)-1|=0$、$Q_n/a_\infty=-1\pm h$ の流束差 7.5e-4 / 7.5e-6 / 7.5e-8 ($h$ 1e-4/1e-6/1e-8、$\propto h$ = 連続)、
内部 $M$ 0.8→1.2 の 1e-5 刻みの隣接差 4.6e-5、音響反射 (一様 M 0: 0.50/0.24/0.08 %、M 0.3: 0.15/0.09/0.04 %、高温内部: 0.15/0.10/0.04 %)、接触波の擾乱 4e-14。

- **既知の限界** (plan-6 M2): TP 多成分では温度・組成の違う気体の保存形の混合だけで圧力が 0.5–0.65 % ずれる (境界と無関係)。境界の試験とは分けて判定する (§6 V0p・V2d)。

### 4.3 面流束 (farfield 専用カーネル `farfield_flux_d`)

- 面流束 = $\mathbf F_{\mathrm{HLLC}}(U_L=U_i, U_R;\hat{\mathbf n})|S|$ ($U_R$ は §4.2、1 次、再構成なし)。3D では法線成分で HLLC を解き、接線速度は星状態で各側の値を保つ (Toro 10.4)。
  HLLC を L/R 状態を引数に取る **`__host__ __device__` の純粋関数**にする (GPU 配列・診断・残差加算は外側)。
- 化学種・$k,\omega$: HLLC の質量流束の符号 (= 接触波速度 $S_*$ の符号) で風上化: $\dot m>0$ なら内部値、$\dot m<0$ なら外側状態の $Y_R,k_R,\omega_R$ (§4.2 の混合後。帯の外では自由流の値)。
  (HLLC の星状態の質量流束は $\rho_*S_*$ なので、スカラー $\dot m\,\phi_{\text{風上}}$ は HLLC の多成分拡張 [星状態で $Y$ を保存] と同じ。)
- 同じカーネルが `massflux[ip]` と、スカラーの面値配列を書く。スカラー移流は node 境界半割面の `ext_is_self` / `nodeBnd` 経路に、**farfield 面だけ面値配列を読む**分岐
  (面フラグで判定、既存面はビット不変)。一次輸送と S3 (化学種) の両経路。呼び出し順は対流流束 → スカラー移流 = 同じ評価時点。
- `sstEnergyIncludesK`: エネルギー・圧力に $k$ を含める経路では、**L = $k_i$、R = 混合・置換後の外側状態の $k_R$** (§4.2。帯の外では $k_\infty$) で $E^*=E+\rho k$、$p^*=p+\tfrac23\rho k$ として HLLC に渡す。移流する $k$ の面値 (流入時) も同じ $k_R$ (plan-10 M1: 右側のエネルギーを $k_\infty$ で作ると混合帯でエネルギー流束が 3.4 % ずれる)。
- ピンしない。非有限の流束はその面の流束を 0 にして面 ID・理由を診断カウンタへ (**検証の評価区間では 0 回**)。
- **面流束の診断ダンプ** (plan-2 M5): env `FORGE_DUMP_FARFIELD=<path>` で、farfield 面ごとに (面 ID、評価回、$\dot m$、運動量 3、エネルギー、化学種、$k,\omega$ の面流束、面積ベクトル $\mathbf S$、圧力基準 `pRef`)。
  運動量は残差に投入した $p-p_{ref}$ 基準のまま書き、SERN 帳簿 (外気圧 $p_a$ 基準) へ渡すときは $(p_{ref}-p_a)\mathbf S$ を加える (plan-3 m7)。
- **帳簿ダンプ (`FORGE_DUMP_LEDGER`) の拡張** (plan-5 m4): 化学種を実際の種数から列挙、全節点指定、一括転送、評価回・段の対応付け。

### 4.4 その他

- ディスパッチ・読込 (`boundaryCond.{hpp,cpp}`): `farfield` を追加、Y/X と k/omega の読込・検査を「入口または farfield」に広げる (`inletCornerWall` は入口専用のまま)。§2 の非対応構成は起動時に拒否。
- 陰解法: 境界半割面は現行の ghostless 対角 A⁺ のまま。
- runner: `evaluate.side_far_kind` / `top_out_kind` に `farfield`、`top_out_kind: outflow` が slip に落ちる不具合も直す。生成 YAML の実効 BC を照合する試験。

## 5. 実装ステップ

1. `boundaryCond.{hpp,cpp}`: `valueTypesOfBC["farfield"]`・読込 (Y/X、k/omega)・ディスパッチ (自由流状態の実物性計算)・起動時の拒否。
2. HLLC の `__host__ __device__` 純粋関数 (CPG/TP、3D 法線、$k$ 込みの $E^*,p^*$ 分岐)。ホスト単体試験から呼ぶ。
3. `farfield_flux_d` (§4.3) と、`convectiveFlux_d.cu` の境界ループでの振り分け (farfield 面は既存の `convectiveFlux_boundary_d` を通らない)。
4. スカラー: `scalarTransport_d.cu`・`passiveKernels_d.cuh`・`speciesTransport_d.cu` の farfield 面値分岐、`ransBoundary_d.cu` に farfield (何もしない分岐を明示)。
5. 診断ダンプ `FORGE_DUMP_FARFIELD`、帳簿ダンプの拡張、収支の積分ツール `solver_density_cuda/tools/farfield_balance.py`。
6. 独立参照解 `solver_density_cuda/tools/ref1d_euler_tp.py` (1D、HLLC、NASA-9 lump、長い領域)。1D 試作 `farfield_proto1d.py` を TP に拡張して流用してよい。
7. runner (`runner_sern3d.py`)・`sern_momentum.py` (farfield 面は診断ダンプを使う)。
8. docs: `methods/boundary.md` (計画中 → 実装済み、ディスパッチ表)、`procedures/recommended-settings.md`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~codex plan 段レビュー~~ **済 (2026-09-27 NO-GO C1/M6/m1、全件採用)** | | F |
| 1b | ~~plan-2~~ **済 (2026-09-27 NO-GO C0/M5/m2、全件採用で本版)** | M1 → §4.2 (判定を $Q_n$ に固定) + §4.3 (近似 Riemann 流束)、M2 → V2d、M3 → V2b、M4 → V2a、M5 → §4.3 ダンプ、m6 → V3 の位置付け、m7 → §2/§3 統一・拒否試験・置換 0 回 | F |
| 1c | ~~plan-3~~ **済 (2026-09-27 NO-GO C0/M6/m1、全件採用)** | | F |
| 1d | ~~plan-4~~ **済 (2026-09-27 NO-GO C0/M3/m2、全件採用)** | | F |
| 1e | ~~plan-5~~ **済 (2026-09-27 NO-GO C0/M2/m3、全件採用)** | | F |
| 1f | ~~plan-6~~ **済 (2026-09-27 NO-GO C0/M3/m1、全件採用)** | | F |
| 1g | ~~ホスト 1D 試作で境界流束の候補を比較~~ **済 (2026-09-28)**。経過: 自由流そのまま + HLLC (plan-7 で高温内部の反射 24–33 % と判明) → 線形の特性外側状態 + HLLC (plan-8 で負圧・超音速流出の独立性・超音速流入の跳び) → **TRRS + 滑らかな重み + HLLC (現採用、§4.2)**。未通過: TP 版 (V0h) (plan-6 M1) | `solver_density_cuda/tools/farfield_proto1d.py`: node 境界半 CV を持つ 1D Euler (CPG、内部 2 次 MUSCL・境界 1 次、RK3)。候補 = (a) 自由流ゴースト + SLAU、(b) 自由流ゴースト + HLLC、(c) 局所線形化の特性振幅で組み立て + HLLC、(d) 同 + SLAU。測る量 = (1) 音響反射率 (M 0 / 0.3、長領域との差、Δx 3 水準)、(2) 超音速流出の極限 (外側状態を変えても流束が内部の物理流束に一致)、(3) 接触波 (P・u 同じで T 違い) の擾乱、(4) 流向反転・音速通過の連続性。合格 = (1) ≤ 5 %、(2) 相対 1e-6、(3) 圧力擾乱 ≤ 1e-4 P∞、(4) 流束の跳びなし。結果で §4 を決め直し、plan-7 に回す | F/O |
| 1h | ~~plan-7~~ **済 (2026-09-28 NO-GO C0/M2/m1、全件採用、試作で対策を確認)** | | F |
| 1i | ~~plan-8~~ **済 (2026-09-28 NO-GO C0/M3/m1、全件採用、試作で確認)** | | F |
| 1k | ~~plan-9~~ **済 (2026-09-28 NO-GO C0/M2/m1、全件採用、試作で確認)** | | F |
| 1l | ~~plan-10~~ **済 (2026-09-28 GO-with-changes C0/M2/m1、全件採用)** → 実装着手可 (V0h ホストゲートから) | | F |
| 1j | ~~V0h ホストゲート (TP 拡張)~~ **済 (2026-09-29、合格)**: `solver_density_cuda/tools/farfield_proto1d_tp.py` (SERN の EXH/AIR 2 擬似種、NASA-9、TP SLAU 内部 + TP HLLC 境界)。超音速流出の極限 (外気の速度 −2〜2 a・圧力 0.2〜5 倍を掃引) で $\max|F/F(U_i)-1|=0$ (外気そのまま + SLAU 16、+ HLLC 40)。連続性 (高温排気の内部、1e-4 a 刻み) 隣接差 ≤ 4e-4、退避 0。正値性 9,477 点で非有限 0・HLL 退避 0・真空置換 242 (極端な離反のみ)。音響反射: 一様な外気 M 0 0.50/0.24 %・M 0.3 0.15/0.09 % (Δx 5/2.5 mm)、**高温排気 (T 600 K・Y 0.13) を出ていく音響 0.16/0.10 %** (外気そのまま + HLLC 32/39 %)、パルスなし対照 3e-11。接触波で境界が加える誤差 4.4e-6/1.3e-6 P (長領域自身 8.0e-6/1.5e-6 P)。V0p (保存形の混合 0.9:0.1) の圧力誤差 +0.654 % (境界と無関係、記録) | O |
| 2 | ~~実装 (§5 の 1–5)、V0・V0u・V0k~~ **済 (2026-09-29)**。V0: farfield を含まない SERN 構成 (run_0971 設定) で初回 massflux・状態ダンプが旧バイナリとビット一致 (case/46 run_0991_ff_v0_old/new)。V0u/V0k: `tests/unit/test_farfield_flux.cu` 全件 PASS — (iv) 一様 1.3e-16、(ii) 接触波 9e-16、(i) 連続性 隣接差 1.1e-5 (≤1e-4)・置換/退避 0、(iii) 243 組合せで帯外の外側値 = 外気・有界、(iii'') 独立掃引 f64/f32 で非有限 0・負値 0 (真空置換 CPG 2144/16875・TP 54/10125、HLL 退避 0)、(iii''') $Y_R$ 0.065・$k_R$ 0.075、V0k 56 面×CPG/TP で相対 2.7e-8。**実装中の不具合**: 外側状態の構造体代入 `R = f` が GPU 上で組成 Y と k・ω を 0 にしていた (ホストでは正常)。明示コピー `ff_copy` と値初期化に替え、V0k で旧版 29/56 面不一致 → 修正版 0 を確認。**V0u (v) 済 (2026-09-29)**: `case/58.farfield_verification/v0u_reject.py` — cell・ROE・凝縮 (CPG N2 担体)・トレーサ・遷移・軸対称・周期共有の 7 構成が拒否メッセージで非 0 終了、対照は完走 (`V0U_REJECT_VERDICT.txt`) | AWS でビルド | O |
| 3 | 独立参照解 (§5 の 6) と V1–V2。**V1 済 (2026-09-29、全 4 本合格、置換・退避 0)**: (a) CPG M0.5 7.5e-7 (`case/58.farfield_verification/run_0001_v1a`)、(b) TP M6 6.7e-6 (`run_0005_v1b_fix`)、(c) 30° 傾き 6.3e-6 (`run_0006_v1c_fix`)、(d) SST 輸送残差の打消し 5.8e-7 (`run_0004_v1d`、帳簿 全節点)。(b)(c) のずれは step 0 で 1e-6 (float32 EOS 往復) → 500 step で 6e-6 に達し以後横ばい。旧 `run_0002_v1b`/`run_0003_v1c` は不具合版 (破棄予定)。**V2a 済 (2026-09-29、全判定合格)**: dual-time、Δx 5 mm・dt 1.25e-6 (音響 CFL 0.11)・nSub 20・cfl_pseudo 12。反射率 M0.3 SLAU **0.099 %** (`run_0040`/`0041`)・SLAU2 0.118 % (`0042`/`0043`)・M0 0.203 % (`0048`/`0047`)、対照 M0 slip 91.3 % (`0046`)。時間精度 nSub 40 との差 0.15 %・dt/2 との差 0.58 % (≤ 1 %)。置換・退避 0。**経緯**: (1) 初回 dt 5e-6 (`run_0010`–`0018`) は dt/2 との差 7.6 % で時間精度不足、(2) dt 1.25e-6 (`run_0020`–`0028`) で差の最大が入射パルス通過時刻 (反射でない) に 1.1 % 出た → 原因は**リミッタ基準値の自動決定** (`limiterRefLength` = 領域の対角長、ρ・P・a = 初期場平均) が短 (L 1.0)・長 (L 3.0) で違い離散化が別物になっていたこと。自由流の値で明示固定して解消。プローブ出力の桁 (既定 6 桁で P が 1 Pa 刻み) も 10 桁にした (`probe/point_probes.cu`、出力専用)。**V3 への含意**: 幅系列も同じ自動決定なので、V3 では `limiterRefLength`・`limiterRoRef`・`limiterPRef`・`limiterARef` を全幅で同一値に固定する (R4d の幅系列は未固定 — 4.35/6.19 H の差 C_L 2.5e-5 から主因ではないと見るが未検証)。<br>**V2b 済 (2026-09-29、6 本 × 3 判定 = 18 件合格)**: `run_0030`–`0035` (+ `_bal` = 最終場から restart した 1 評価)。離散恒等式 (初回評価) 全成分 |差|/規模 ≤ 1e-9、Σ種 = 質量 相対 ≤ 2e-8、定常後の全体収支 |R|/代表量 ≤ 6e-9、Y_EXH 最大 5e-36 (置換)。残差は float 床まで 4.4 桁低下 (`rms_roY0` は恒等的 0 で判定不能 = 完全置換)。**帳簿の読み方の訂正**: 各残差配列は担当段の中でゼロに戻されるので、2 回目以降の評価では「段の前後の差」は寄与にならない (段の前に前回の最終残差が残る)。`eval_v2b.py` は担当段直後の値を寄与とする。<br>**V2f 済 (2026-09-29、合格)**: `run_0036_v2f` — 評価 1 で xmax の 169 面すべてが ṁ<0、その面の外側組成 = 外気、化学種の恒等式で流入面の Y = 外気 (Σ 1.6e-10)・流出面 (評価 2・3) で内部値 (相対 1e-10)、置換・退避 0、定常後は自由流に戻る (Y 1e-36)。<br>**V2e 済 (2026-09-29、合格)**: 音響 + SST dual-time 反射 0.12 % (ek0、`run_0052`/`0053`)・0.33 % (ek1、`0054`/`0055`)、V1b 陽解法 (RK3、局所 dt) ずれ 8.7e-6 (`run_0050`)、V2b(i) 陽解法 恒等式・収支・置換 PASS、残差 5.7 桁低下 (`run_0051` + `_bal`)。<br>**V2c FAIL (2026-09-29、保持)**: `run_0063_v2c_A`/`0064_B`/`0065_C` (ランプは 0.7 m で終わる形。出口まで伸ばした `run_0060`–`0062` は A が流路閉塞で発散、破棄)。B vs C 最大 0.087 Δp (x 0.21、ランプ角直後)、角・出口端を除いても 0.048 Δp (基準 0.02)。A vs C 2.61 Δp (対照 OK)。角の値は A = B (6 桁一致)・C だけ違う、cfl 0.5 (`0066`/`0067`) と再変換 (`0068`) で不変。上面の境界節点で B が C より 0.10 Δp 低く、弱い反射波 (内部 0.06 Δp) が出口寄りの下面に届く。**原因未分離** (codex diagnose 採用)。<br>**V2d**: 接触波 V2d-1 は 3 物性 PASS (`run_0070`/`71`・`75`/`76`・`80`/`81`)。**V2d-2 パルスなし対照 FAIL (保持)**: cpg 2.4e-4 Pa PASS、tp1 1.64 Pa・tp2 2.49 Pa (基準 0.0285 Pa)。長領域でも同じ (`run_0085` 2.4888 Pa)、左端の自由流を内部と同じ高温にすると 0.0000 Pa (`run_0086`/`0087`、隔離対照として別記) → 左端流入で作られる接触面の TP 保存形混合が音源という仮説 (未確定)。V2d-2 の反射 (0.07–0.23 %) は左端擾乱が反射より先に評価点に届くため隔離試験になっていない → 再設計。 | AWS | O |
| 3a | **V2c 角の作用素 A/B** (codex diagnose 2026-09-29): 低領域 (A 格子、上面 slip) と高領域 (C 格子、上面 slip) に、C の最終場を座標対応でそのまま与え、1 step の初回評価だけを比較。対象 = 下面 x 0.195–0.230 と再構成に要る近傍。照合 = 双対体積・面積ベクトル・`after_eos_bc` 状態・面の再構成状態と流束 (帳簿 `.faces`)・`res_after_conv`。**判定 (事前固定)**: 各保存量の残差差 / その節点に接する面流束の絶対和 ≤ 1e-5 (ゼロ規模は自由流基準)。超えたら局所作用素の交絡を採用し最初に異なる幾何・再構成・流束へ絞る、全て許容内なら作用素差を棄却し反復過程へ移る。どちらでも収束・境界性能は判定しない。**結果 (2026-09-29)**: 作用素は同一 (`run_0092`/`0093`、幾何・状態・面の差 0、残差差 ≤ 7.7e-8、`V2C_OPAB_VERDICT.txt`) → 第 1 仮説棄却。C の最終場から低領域 B・A・C 自身を 6000 step 継続 (`run_0094`–`0096`) → 角は 3 本とも C の値を保持。同じ初期場での B vs C は全線 0.0210 Δp (出口端)、除くと 0.0080。codex diagnose 2 回目 ([記録](../../notes/reviews/2026-09-29-farfield-v2c-corner-and-v2d2-redesign-diagnose.md)、全件採用) の絶対残差 A/B (`run_0097`/`0098`、同じ B 格子に 2 つの最終場、事前閾値 \|R\|/Σ\|F\| ≤ 1e-5、slip 節点の運動量は接線成分、運動量は 3 成分共通の規模、規模 0 は自由流基準): **角 (x 0.195–0.23) は 2 状態とも全 5 成分 ≤ 5.2e-7**、しかし**場全体では約 2000 節点が超過し最大はランプ終端の凸角 (x 0.715、下面) で 4.6e-3 / 7.3e-3** → 規則どおり「2 つとも離散定常解」を棄却、非零残差が残る凸角を追う (`V2C_RESAB_VERDICT.txt`)。**次**: 凸角の停滞の扱いを諮る (#3c) | AWS | O |
| 3c | **V2c 凸角 (ランプ終端 x 0.7) の残差停滞**: 形状側 (slip の凸角) の問題として扱うか、境界と無関係に node slip 凸角の既知欠陥か、を切り分ける。案を codex に諮ってから §6 V2c の手順 (共通初期場の追加隔離試験を含む) を事前登録し直す | AWS | F |
| 3b | **V2d-2 の再設計** (codex diagnose 2026-09-29、2 回目で採用確定: `--left-hot` を右端流出の隔離試験として正式化、短・長の両方でパルスなし対照、3 物性、評価器は入射窓で振幅・反射到達窓で短長差を測る形に直す): 右端流出を隔離する対照 (左端 = 内部と同じ高温状態) を正式な隔離試験として別記、元の配置の FAIL は保持。反射の評価は左端擾乱の到達前に反射が評価点へ着く配置に直す (§6 の該当節を改訂してから回す)。左端流入の TP 混合仮説は境界を含まない試験で確かめる | AWS | F |
| 4 | SERN V3 | AWS。R4d へ反映 (restart は `r4d_common_restart.py` 型の index コピー、メッシュ品質、判定区間、case README の run 索引) | O |
| 5 | codex result 段 → accepted | | F |

## 6. 検証

すべて AWS。合否は実行前に固定。判定区間を明記し `check_convergence.py`・`check_quasisteady.py` の出力原本を run に残す。
**非有限流束の置換 (§4.3) は各試験の評価区間で 0 回**を必要条件とする。
V1–V2 は**境界機能の受入れ**、V3 は**SERN での配置 (側方幅) の採否**。

- **V0 既存境界の不変**: farfield を含まない構成 (run_0971 設定) で、同一初期状態からの初回 `massflux` と状態ダンプが変更前バイナリとビット一致。更新後保存量は旧バイナリ 3 回反復の再現性幅以内。
- **V0h ホストゲート (本体実装の前、plan-7 推奨)**: `farfield_proto1d.py` を TP (単成分・多成分 lump、既存 `FrozenGas`) に拡張し、§4.2 の表の全項目 + V2d-2 (TP の高温内部を出ていく音響、反射 ≤ 5 %) + 波速順序・星状態の妥当性・退避 0 を確認してから CUDA に進む。
- **V0p TP の保存形混合 (前提、境界と無関係)** (plan-6 M2): 単一 CV に高温側 (T 600 K、Y 0.13) と外気 (T 220 K、Y 0) の保存量を 0.9:0.1 で混ぜ、EOS から復元した圧力の誤差を記録する (期待 +0.65 %)。これは境界の合否に入れず、V2d の「長領域自身の誤差」の説明に使う。
- **V0u 単体試験** (ホスト、HLLC の `__host__ __device__` 関数): (i) **連続性**: $U_{n,i}$・$\rho_i$・$P_i$・接線速度を、法線 Mach ±1・0 付近で 1e-6 刻みに掃引し、面流束 5 成分・スカラー面流束の隣接点の差 ≤ 1e-4 × 流束規模 (plan-2〜5 の反例入力をすべて含む)、(ii) **接触波**: CPG と TP で P・$U_n$ 同じ・T 600/220 K・Y 0.13/0 の内部/外気、流出と流入の両方で面流束 = 風上側の物理流束 (相対 1e-6)、(iii) **流向の全組合せ**: $Q_n$ の符号 × 実際の $\dot m$ の符号 × ゼロ通過で、$\dot m<0$ のスカラー面値が外側状態の値 (帯の外では外気値)・$\dot m>0$ で内部値、(iii'') **真空と float32** (plan-9 M1): 真空の反例と、速度 ($U_{n,i}$ −3〜3、$U_{n,\infty}$ −3〜10)・密度比・圧力比の独立掃引を float64 と float32 で行い、非有限 0・外側状態の負値 0 を確認、真空置換と HLL 退避の回数を別々に記録 (有限であることと置換が無いことを分けて判定)、(iii''') **混合帯の逆流**: 内部 $(1,0.95,1/\gamma,Y\,0.13)$・自由流 $(10,-3,10/\gamma,Y\,0)$ で、流入の種の面値が外側状態の $Y_R$ (0.065) と一致すること、TP・SST (k 込み) でも同じ。`sstEnergyIncludesK: 1` で $k_i=0.14$、$k_\infty=0.01$ を与え、運動量・全エネルギー・$\rho k$ の流束が同じ $k_R$ (0.075) から作られていること (plan-10 M1)、(iv) 一様 ($U_i=U_\infty$) で $F_{num}=F(U)$ (相対 1e-6)、(v) 非対応構成の起動拒否 (cell・ROE・凝縮・トレーサ・遷移・軸対称・周期共有) がそれぞれエラー終了。
- **V0k CUDA 一致**: 同じ入力を GPU カーネルに与え、ホストの HLLC 関数と面流束が相対 1e-6 で一致。
- **V1 自由流保持** (定常、block-DPLUR): 一様流の 3D hex 直方体 (z 3 層以上)、全境界 farfield、2000 step。(a) CPG M 0.5、(b) TP 2 種 lump (SERN 外気) M 6、(c) (b) を面に 30° 傾ける、
  (d) (b) + SST。合格 (a)–(c): 全節点で $|P/P_\infty-1|,|\rho/\rho_\infty-1|,|\mathbf u-\mathbf u_\infty|/|\mathbf u_\infty|\le10^{-5}$、$|Y-Y_\infty|\le10^{-6}$。(d) は **一様状態での輸送残差の打消し**だけを見る (plan-3 M3: SST は一様流でも消滅項 $-\beta^*\rho k\omega$、$-\beta\rho\omega^2$ があり一様値は保たれない): 初回評価で帳簿ダンプの段別残差のうち対流 + $k,\omega$ 輸送の寄与が全節点で、保存量ごとの規模 (その節点に接する面流束の絶対和) の 1e-5 以下 ($\rho,\rho\mathbf u,\rho E$ は各流束、$k$ は $\dot m k_\infty$、$\omega$ は $\dot m\omega_\infty$。plan-4 m5)。ソース込みの時間発展は V2b で見る。
- **V2a 音響反射 (非定常 dual-time)**: 3D 薄板チャネル (z 3 層以上、長さ 1 m、Δx 5 mm)。
  - 本試験: 一様流 M 0.3 (CPG、空気)、ガウス圧力パルス (振幅 1e-3 P∞、半値幅 40 mm = 8 セル以上)、下流端 farfield。参照 = 同じ背景流で長さ 3 m (反射が評価窓内に戻らない)。
    評価点 = 下流端から 0.2 m 上流、評価窓 = 入射パルス通過後から 0.8 m / (c−u) までの反射到達時刻窓。**反射振幅 (参照との差の最大) ≤ 0.05 × 入射振幅**。
  - 対照 (試験が反射を検出できること): M 0 (静止) で下流端 slip → 反射 ≥ 0.9 × 入射。
  - 時間積分: `unsteady: 1`・dual-time、物理 CFL ≤ 1 (音響)。**時間精度の確認** (plan-3 M6、[`recommended-settings.md`](../../procedures/recommended-settings.md) の dual-time 節): `nSub` を倍にした run と Δt を半分にした run で、評価点の圧力時系列の差が合格許容の 1/5 以下 (= 入射振幅の 1 %)。非定常の精度判定は定常用 VERDICT と別ファイルで保存する。
  - SLAU と SLAU2 の両方。
- **V2b 保存収支** (定常、block-DPLUR): V1 の箱で内部を自由流と異なる状態 (Y 0.13、T 600 K、k・ω 10 倍) から始め、(i) 流出配置、(ii) 流入配置、`sstEnergyIncludesK` 0/1 の 4 本。
  - **配置の定義 (2026-09-29、実行前に確定)**: 代表量が $|\mathbf u_\infty|$ を含むので静止外気は使わない。自由流 = SERN 外気 (TP、220 K、2851 Pa、Y_EXH 0、SST k 479.653・ω 119844)、M 0.5 の通過流。(i) **流出配置** = +x 方向 (xmin だけ流入、xmax 流出、側面 4 面は $Q_n=0$ の接線流 = SERN 側面と同じ状況)、(ii) **流入配置** = 方向 $(\cos30°\cos20°,\sin30°\cos20°,\sin20°)$ (xmin・ymin・zmin の 3 面から斜めに流入、残り 3 面から流出)。初期値は内部全域を P∞・u∞ のまま T 600 K・Y_EXH 0.13・k・ω 10 倍。4 本 = (i)(ii) × `sstEnergyIncludesK` 0/1 (化学種は `speciesFaceReconstruction` 0)、加えて化学種の 2 経路用に (i)(ii) × `speciesFaceReconstruction` 2 (`sstEnergyIncludesK` 0) の 2 本。離散恒等式は初回評価 (非一様の初期場) の帳簿 + 面ダンプで、定常後の全体収支は最終場から同一格子 restart した 1 評価の帳簿で見る。
  - **離散恒等式 (毎評価)**: 全節点の対流残差の和 = −(farfield 面流束の和) (内部面は相殺)。診断ダンプ (§4.3) と帳簿ダンプ (`FORGE_DUMP_LEDGER`、全節点に印) で照合。**判定 (plan-4 m5)**: |差| ≤ 1e-5 × (その保存量の全面流束の絶対和) (正味和で割らない。float32 の atomicAdd 積算の丸め幅を含む。ホスト集計は double)。
  - **定常後の全体収支** (plan-3 M4): 各保存量で残差の正本を $R=-\sum_{faces}F+\sum_{nodes}S$ (外向き流束は残差に $-F$、生成ソースは $+S$) とし、定常で $R\to0$ を見る。$S$ は $k,\omega$ の体積ソース (帳簿の段別残差 `res_after_sources` − `res_after_species`) で、他の保存量は 0。`sstEnergyIncludesK: 1` の全エネルギーは `res_roe` が既に全エネルギーの残差なので **`res_roK` を足さない** ([`turbulence-sst-energy-includes-k.md`](../accepted/turbulence-sst-energy-includes-k.md) §分割保持)。規格化は**事前登録した非零の代表量** (ゼロ成分を落とさない): 質量・化学種 $\rho_\infty|\mathbf u_\infty|A$、運動量 $(\rho_\infty|\mathbf u_\infty|^2+P_\infty)A$、エネルギー $\rho_\infty|\mathbf u_\infty|H_\infty A$、$k$ $\rho_\infty|\mathbf u_\infty|k_\infty A$、$\omega$ $\rho_\infty|\mathbf u_\infty|\omega_\infty A$ ($A$ = farfield 面の総面積)。**$|R|$/代表量 ≤ 1e-4** を全成分に。
  - (ii) で内部が自由流へ置き換わること ($|Y-Y_\infty|\le10^{-4}$)。
  - **化学種の 2 経路** (plan-10 M2): 流入・流出の両配置を `speciesFaceReconstruction: 0` (一次輸送) と `2` (S3、`convMethod ≥ 1`・`speciesImplicitCoupling: 1`) の両方で回し、実効設定を run に保存する。各経路で、流入組成が外側状態の値であること、Σ種流束 = 質量流束 (相対 1e-6)、残差との離散恒等式を確認する (一次経路だけ直した実装を通さないため)。
- **V2f 自由流分類と局所流況が違う面** (plan-5 M2): V1 の箱で、自由流は面から超音速流出 ($Q_n=2a_\infty$) だが、内部の初期状態を局所に逆流 ($U_{n,i}=-0.2a$) させた帯を置く。 **配置 (2026-09-29 確定)**: TP 外気 +x M 2 (xmax で $Q_n=2a_\infty$)、初期値は x ≥ 0.8 m の帯だけ $U_x=-0.95a$・Y_EXH 0.13 (他は自由流)。**逆流の強さはホスト関数で決めた** (CPG、同密度・同圧): $Q_n=2a_\infty$ の自由流は TRRS の J⁻ で外向きに強く引くので、$U_{n,i}=-0.2a$ では $\dot m=+0.29\rho c$ (流出のまま)、−0.5a で +0.13、−0.8a で −0.0065、**−0.95a で −0.018** (置換・退避 0)。帳簿を最初の 20 評価記録し、実際に $\dot m<0$ となった xmax の面を面ダンプから数える。合格 (plan-6 M3): **実際の $\dot m<0$ の面**で化学種面値が外気値・$\dot m>0$ の面で内部値 (内部速度の符号ではなく流束の符号で判定)、NaN・置換 0、定常後は自由流へ戻る。逆流帯の強さは、HLLC で実際に $\dot m<0$ になる面が生じることをホスト関数で事前に確かめて決める。
- **V2e 時間積分経路の網羅** (plan-3 M6): V1(b) と V2b(i) を **定常陽解法** でも回す (合格条件は同じ)。V2a の本試験に SST を足し、dual-time で `sstEnergyIncludesK` 0/1 の両方 (反射 ≤ 0.05、NaN・置換 0)。
- **V2c 超音速の斜め衝撃波** (定常): 3D 薄板、M 2.5、半角 10° ウェッジ、衝撃が上側境界に当たる配置。上側境界 (A) slip / (B) farfield / (C) 上方に、反射衝撃が評価線の下流端より後ろに着く高さまで広げた slip
  (共通領域の格子は同一、`z_append` 型)。評価線 = ウェッジ下流の壁面、全線。合格: $\max|p_B-p_C|\le0.02\,\Delta p_{shock}$ かつ A は 0.1 Δp 以上。
- **V2d TP での境界の独立参照** (非定常、plan-3 M5 で初期値と誤差予算を固定): 1 次元問題を独立コード `ref1d_euler_tp.py` (HLLC、NASA-9 lump) と forge の薄板 (境界 = farfield、対照に同格子で 3 倍長の forge 長領域) で解く。格子 Δx 5 mm、長さ 1 m、自由流 = SERN 外気 (T 220 K、P 2851 Pa、Y_EXH 0)。2 種類に分ける:
  - **V2d-1 接触波 (温度・組成の塊の流出/流入)**: 初期 = P・u 一様 (P∞、u = ±0.5 c∞)、中央にガウス型 (半値幅 0.1 m) の T 600 K・Y_EXH 0.13 の塊 (ρ は P・T・Y から)。連続系では圧力擾乱ゼロの移流。合格: 塊が境界を通過する時間窓で、境界から 5 セル内側の $|P-P_\infty|/P_\infty\le10^{-3}$、T・Y の時系列が参照と一致 (T は塊の振幅の 2 %、Y は 0.002 以内)。
  - **V2d-2 TP 音響**: 同じ高温・高 Y の一様状態 (T 600 K、Y 0.13) の中を進む圧力パルス (振幅 1e-3 P) が、外気 (T 220 K、Y 0) を自由流に持つ farfield に当たる (= 内外の γ が違う)。**パルスなしの対照を必須** (plan-4 M1: 内外の状態差だけで境界が擾乱を出さないこと。評価点の $|P-P_0|$ ≤ 0.01 × パルス振幅) を先に通し、その上で反射振幅 ≤ 0.05 × 入射 (参照との差で)。
  - 誤差予算: 独立参照は Δx を 1/4 にした解との差が合格許容の 1/5 以下であること (参照の収束確認)。**判定を 2 つに分ける** (plan-6 M2): (a) forge 長領域 vs 独立参照 = 既存の内部離散化 (TP 保存形混合を含む) の誤差 — 記録のみ、合格条件に入れない。(b) forge 短領域 vs forge 長領域 = **境界が加える誤差** — 上の合格条件はこちらに適用する。
  - CPG・単成分 TP・多成分 TP の 3 物性で行う。
- **V3 SERN での配置 (側方幅)**: 境界機能の受入れ (V1–V2 合格) の後。g3、同一新バイナリ。
  - V3a 同一格子の対照: 遠方面 2.50 H で `side_far` のみ slip ↔ farfield (他の BC は生成 YAML で照合)。run_0986 最終場から各 20000 step (同一格子 restart)。R4d と同じ ε・D・窓条件で差を**記録** (合否なし)。
  - V3b farfield の幅系列: 2.50 / 3.42 / 4.35 H (`z_append`、restart は `r4d_common_restart.py` 型、メッシュ品質 PASS)。
    **採否** (plan-4 M3・plan-5 m5): 候補幅を**それより広いすべての試験済み幅** (少なくとも 1 つ必須) と比較し、全 4 量で D ≤ ε (C_T・C_T_with_shear・C_L 5e-4、C_M 5e-3) を満たす最小の幅を「farfield での必要幅」とする。最大幅しか残らない (比較相手が無い) ときは「必要幅未確定」とし、系列を広げる。結論は「試験した幅系列で許容内」に限定し、無限遠との一致とは言わない。
  - SERN の帳簿 (C_T 等の境界寄与を使う場合) は診断ダンプの面流束を使う。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan.md) | NO-GO, C1/M6/m1 | **全件採用** (C1・M2 はコードで再確認: `scalarTransport_d.cu:166` の `ext_is_self`、`convectiveFlux_boundary_d.inc.cuh` の流出時内部風上)。§4 を「面流束を自前で組む専用境界」に改訂、§6 に V0u/V2a–c/V3a–b。§5.1 #1 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-2.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-2.md) | NO-GO, C0/M5/m2 | 構造 (専用流束 + 同時刻スカラー面値) は妥当。**全件採用 (2026-09-27、ユーザ決定「レビューどおり全部やる」)**: M1 → §4.2 判定を自由流 $Q_n$ に固定 + §4.3 近似 Riemann (SLAU) 流束、M2 → V2d (独立 1D 参照)、M3 → V2b (離散恒等式 + SST ソース込み全体収支)、M4 → V2a (dual-time、長領域参照、M0 の slip 対照)、M5 → §4.3 面流束ダンプ、m6 → V3 を配置の採否に、m7 → §2/§3 統一・拒否試験・置換 0 回 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-3.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-3.md) | NO-GO, C0/M6/m1 | **全件採用**: M1 → §4.2 判定用 $a_\infty$ を自由流の実物性で固定・V0u(iii')、M2 → V0u(iii) 構成状態と数値流束の極限を分離、M3 → V1(d) を輸送残差の打消しに、M4 → V2b の符号・代表量・res_roK、M5 → V2d-1/2 と誤差予算、M6 → nSub 倍・Δt 半減、V2e、m7 → §4.3 ダンプの圧力基準 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-4.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-4.md) | NO-GO, C0/M3/m2 | **全件採用**: M1 → §4.2 を局所線形化の特性振幅 (Whitfield–Janus/Blazek) に変更・V0u(vi)・V2d-2 のパルスなし対照、M2 → エントロピー部を実際の流向 ($U_{n,b}$、食い違い時は $\dot m$) で選択・V0u(vii)、M3 → V3b の採否を全広幅と比較、m4 → 明示 $Y,T$ の熱物性関数・V0u(viii)、m5 → 規格化の保存量別・絶対和 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-5.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-5.md) | NO-GO, C0/M2/m3 | **全件採用**: M1 (密度側の選び直しが SLAU の $\dot m$ と自己整合しない)・M2 (固定の超音速分類が局所逆流の情報を失う) → §4.2 を「外側に自由流をそのまま置き SLAU で解く」に変更 (組み立て・分類を撤去)、V0u(i)(iii)・V2f。m3 → 流束関数を `__host__ __device__` に・V0k、m4 → 帳簿ダンプ拡張、m5 → V3b の最大幅規則 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-6.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-6.md) | NO-GO, C0/M3/m1 | **全件採用**。M1 (自由流ゴースト + SLAU は超音速流出でも外側速度が流束を 25 % 変え、1D 線形化で音響反射 14–18 % [格子細分で下がらない]) → 境界だけ HLLC を候補にし、**全体実装の前にホストの 1D 試作で候補を比較** (§5.1 #1g)。M2 (TP 多成分の保存形混合で圧力 +0.5〜0.65 %、境界と無関係の既存性質) → V0u に単一 CV 更新 + EOS 復元の前提試験、V2d で長領域自身の誤差と短領域の追加誤差を分離。M3 → V2f の判定対象を「実際の ṁ<0 の面」に。m4 → §5・§6・methods・README の同期 |
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-7.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-7.md) | NO-GO, C0/M2/m1 | **全件採用**。M1 (外気そのままの外側状態は、高温の内部を出ていく音響を 24–33 % 反射) → §4.2 外側状態の音響部を境界節点の $Z$ で線形化した特性量に (試作で 0.03–0.15 %)、V0h ホストゲートに TP 版を追加。M2 (音速の算術平均の波速で $S_*>S_R$・流束の跳び) → Davis 波速 + 順序・星密度の検査 + HLL 退避と計数 (試作の極端な比の掃引で退避 0・跳びなし)。m3 → §7・methods の同期 |
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-8.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-8.md) | NO-GO, C0/M3/m1 | **全件採用、試作で対策を確認**: M1 (線形 $w^\pm$ で負圧・NaN) → TRRS (内部エントロピー) で常に正、M2 (内部超音速流出でも外気で流束が 19 % 変わる) → 内部 $M_{n,i}\ge1$ で $U_R=U_i$ に滑らかに寄せ最後に適用、M3 (自由流の超音速流入判定で跳び 0.215) → 滑らかな重み、m4 → §5.1 #1g と README を同期 |
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-9.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-9.md) | NO-GO, C0/M2/m1 | **全件採用、試作で確認**: M1 (TRRS も真空生成で非正、試作のクランプは float32 で NaN・未計数) → 真空判定 $B\le10^{-3}(c_i+c_{po})$ で $U_R=U_i$ に置換して計数・V0u(iii'') に float32 と独立掃引、M2 (混合変数が未定義で EOS とスカラーの状態が分離) → 原始変数 $(\rho,\mathbf u,P,Y,k,\omega)$ を同じ重みで混合・流入スカラーは $\phi_R$・V0u(iii''')、m3 → methods 一覧表を同期 |
| plan | 2026-09-28 | [2026-09-28-boundary-node-farfield-characteristic-plan-10.md](../../notes/reviews/2026-09-28-boundary-node-farfield-characteristic-plan-10.md) | **GO-with-changes**, C0/M2/m1 | **全件採用**: M1 → §4.3 の SST 右状態を混合後の $k_R$ に統一・V0u(iii''') に k 込みの照合、M2 → V2b に `speciesFaceReconstruction` 0/2 の両経路、m3 → methods・説明ページの非反射の主張を「法線方向の小振幅音波」に限定し、斜入射の線形反射係数 $(1-\cos\theta)/(1+\cos\theta)$ (45° で約 17 %) を明記。codex はホストで反射率・独立掃引を再現 |
| diagnose | 2026-09-29 | [2026-09-29-farfield-v2c-v2d-fail-diagnose.md](../../notes/reviews/2026-09-29-farfield-v2c-v2d-fail-diagnose.md) (brief [2026-09-29-farfield-v2c-v2d-fail.md](../../notes/reviews/briefs/2026-09-29-farfield-v2c-v2d-fail.md)) | V2c・V2d-2 の FAIL を保持、Major 4 | **全件採用**: (1) V2c を「境界の実力 ≈6 % 反射」と確定しない (評価器は反射係数を分離していない)、(2) 角を評価から除外しない (原因を同定・修正して全線で判定し直す)、(3) Z = ρU/√(M²−1) への変更は現段階で却下 (上面 Mₙ≈0.36 で小振幅・接線流の仮定外)、(4) V2d-2 を「右端擾乱ゼロ」で合格にしない (left-hot は隔離対照として別記)。V2d-2 の反射 PASS も左端の擾乱が反射より先に評価点に届くため隔離試験でない → 再設計。次の一手 = V2c 低・高領域に同一状態を与えた 1 評価の作用素 A/B (§5.1 #3a) |

## 7. 影響範囲

- ソルバ: `solver_density_cuda/boundaryCond.{hpp,cpp}`、`cuda_forge/boundaryCond_d.cu` (自由流状態の起動時計算)、`cuda_forge/convection/convectiveFlux_d.cu` (境界ループの振り分け)、
  新規 `cuda_forge/convection/farfieldFlux_d.inc.cuh` (HLLC・外側状態・`farfield_flux_d`)、`cuda_forge/scalarTransport_d.cu`・`cuda_forge/passiveKernels_d.cuh`・`cuda_forge/speciesTransport_d.cu` (farfield 面値)、
  `cuda_forge/ransBoundary_d.cu`、`cuda_forge/convection/convectiveFlux_d.cu` の帳簿ダンプ拡張。
- ツール: `solver_density_cuda/tools/farfield_proto1d.py` (ホスト試作・ゲート)、`farfield_balance.py`、`ref1d_euler_tp.py`。
- 設計側: `design/forge_design/evaluate/runner_sern3d.py`、`design/forge_design/metrics/sern_momentum.py`。
- 既存ケース: 変更なし (新種別を書かなければ既存経路はビット不変、V0)。
- docs: `methods/boundary.md` (一覧表・理論節・ディスパッチ表)、`methods/index.md` (見出しのみ)、`procedures/recommended-settings.md`。

## 8. 完了条件

- [ ] 関連 `methods/` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6 V0–V3)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`
- [ ] [`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-29` — V2e 合格。V2c・V2d-2 が事前基準で FAIL → codex diagnose (全件採用、FAIL 保持、§5.1 #3a・#3b)。R4d 幅系列はリミッタ自動基準値 (p_ref 最大 10 %・L_ref) が幅ごとに違い、別バイナリ (`sglsq/forge_2fa3826c`) だったことを記録 (V3 は新バイナリ・基準値固定で回す)。
- `2026-09-29` — V0u (v)・V2a・V2b・V2f 合格 (§5.1 #2・#3)。V2a でリミッタ基準値の自動決定が領域長に依存し比較を汚すことを発見 → V3 で固定する。
- `2026-09-29` — GPU 実装完了。V0 (ビット一致)・V0u/V0k (単体試験 `tests/unit/test_farfield_flux.cu`)・V1 (a)–(d) 合格。GPU 上の構造体代入で組成・k・ω が落ちる不具合を V1b で検出し修正 (§5.1 #2)。
- `2026-09-29` — V0h (TP ホストゲート) 合格 (§5.1 #1j)。GPU 実装に着手。
- `2026-09-28` — plan-10 **GO-with-changes** (C0/M2/m1) を全件採用。SST の右状態を混合後の $k_R$ に、化学種の一次・S3 両経路を V2b に、非反射の主張を法線方向に限定 (斜入射 45° で線形理論でも約 17 %)。status を in_progress へ。
- `2026-09-28` — plan-9 NO-GO (C0/M2/m1) を全件採用。真空・近真空は $U_R=U_i$ に置換して計数、混合は原始変数を同じ重みで、流入スカラーは混合後の外側状態の値。試作に float32 と独立掃引のゲートを追加し合格。
- `2026-09-28` — plan-8 NO-GO (C0/M3/m1) を全件採用。外側状態の音響部を TRRS (内部エントロピー、常に正) に、超音速の境目を滑らかな重み (内部の超音速流出を優先) に。試作で plan-8 の全ゲート合格。
- `2026-09-28` — plan-7 NO-GO (C0/M2/m1) を全件採用。外側状態の音響部を境界節点の $Z$ で線形化した特性量、密度・組成は常に自由流側、波速は Davis + 検査 + HLL 退避。試作で全項目合格。「ゴースト」の表記を node の弱形式の外側状態に改めた (ユーザ指摘)。
- `2026-09-28` — ホスト 1D 試作で 5 候補を比較し、自由流ゴースト + 境界面だけ HLLC に決定 (§4.2 の表)。§5・§6 を同期 (V0p・V2d の判定分離・V2f の判定対象)。
- `2026-09-28` — plan-6 NO-GO (C0/M3/m1) を全件採用。机上の反復をやめ、ホスト 1D 試作で境界流束の候補を比較してから §4 を決める (§5.1 #1g)。
- `2026-09-27` — plan-5 NO-GO (C0/M2/m3) を全件採用。境界状態の組み立てをやめ、外側に自由流を置いて SLAU で解く方式に変更 (frozen-γ も不要に)。
- `2026-09-27` — plan-4 NO-GO (C0/M3/m2) を全件採用。境界閉包を Riemann 不変量の混合から局所線形化の特性振幅に変更 (接触波で擾乱を出さない)。
- `2026-09-27` — plan-3 NO-GO (C0/M6/m1) を全件採用 (判定用音速の固定、合格条件の修正)。
- `2026-09-27` — plan-2 NO-GO (C0/M5/m2) を全件採用し全面改訂 (ユーザ決定: 検証はレビューどおり全部)。
- `2026-09-27` — codex plan 段 NO-GO (C1/M6/m1) を全件採用し §4/§6 を改訂。
- `2026-09-27` — 初稿 (ユーザ「遠方境界入れたらすっきりかもね。やってみますか」)。`methods/boundary.md` に理論節を追加し、`outflow` の説明 (実装は全量コピー) を訂正。
```

## 参考: `case/58.farfield_verification/setup_v2c.py`

```
#!/usr/bin/env python3
"""V2c 超音速の斜め衝撃波 (plan boundary-node-farfield-characteristic §6 V2c) の run を作る (AWS で実行)。
3D 薄板 (z 3 層、5 mm、slip)。下面 = 平板 (x ≤ 0.2 m) + 10° ランプ (0.2–0.7 m) + 水平 (0.7 m–、slip、x ≥ 0.2 m が評価線)。
ランプを出口まで伸ばすと流路が 0.3 → 0.124 m に絞られ、slip 上面 (A) で閉塞して発散した (2026-09-29) ので 0.7 m で止める、M 2.5 (CPG 空気 300 K、101325 Pa)。
ランプ角の衝撃波 (β ≈ 31.9°) が上面に当たる。上面:
  --case A: H 0.3 m、slip (反射する対照)      --case B: H 0.3 m、farfield
  --case C: H 0.8 m、slip (上面での反射点が x > 1.2 m になる高さ。共通領域 0 ≤ j ≤ 60 の格子は A/B と同一)
格子: x 0–1.2 m を 240 分割 (Δx 5 mm)。y は 0 ≤ j ≤ 60 で下面から H_A = 0.3 m (水平) までを等分、C はその上に 5 mm で 100 行足す。
入口 xmin = farfield (超音速流入 → 自由流)、出口 xmax = outflow。定常 block-DPLUR。
  python3 setup_v2c.py RUN --case A|B|C [--steps 6000] [--cfl 1.0]
"""
import argparse, math, os, subprocess, sys
import h5py, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_box_msh import write_hex_msh  # noqa: E402

XC, XE, TH, LX, HA, NX, JA, DZ, NZ = 0.2, 0.7, math.radians(10.0), 1.2, 0.3, 240, 60, 0.005, 3


def yb(x):
    return (min(max(x, XC), XE) - XC) * math.tan(TH)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--case", choices=("A", "B", "C"), required=True)
    ap.add_argument("--steps", type=int, default=6000); ap.add_argument("--cfl", type=float, default=1.0)
    ap.add_argument("--bin", default=os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff"))
    a = ap.parse_args()
    run = a.run
    os.makedirs(run)
    JC = 100 if a.case == "C" else 0
    ny = JA + JC
    xs = np.linspace(0.0, LX, NX + 1)

    def coord(i, j, k):
        x = xs[i]
        y = yb(x) + (j / JA) * (HA - yb(x)) if j <= JA else HA + (j - JA) * DZ
        return x, y, k * DZ
    write_hex_msh(os.path.join(run, "ramp.msh"), NX, ny, NZ, coord, (LX, HA + JC * DZ, NZ * DZ))
    P, T, gam, cp, M = 101325.0, 300.0, 1.4, 1004.5, 2.5
    R = cp - cp / gam; ro = P / (R * T); c = math.sqrt(gam * R * T); U = M * c
    open(os.path.join(run, "solverConfig.yaml"), "w").write(f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "ramp.h5", valueFileName: "ramp.h5"}}
gpu: 1
solver: "SLAU"
physProp: {{thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: {cp}, gamma: {gam}}}
time:
  unsteady: 0
  dualTime: 0
  last: {{nStepOuter: {a.steps}}}
  deltaT: {{control: 1, dt: 1e-8, cfl: {a.cfl}, cfl_pseudo: {a.cfl}, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {max(a.steps // 6, 1)}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {P}, limiterRefLength: 1.0, limiterRoRef: {ro!r}, limiterPRef: {P!r}, limiterARef: {c!r}}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
""")
    fl = f"ro: {ro!r}, Ux: {U!r}, Uy: 0.0, Uz: 0.0, Ps: {P!r}"
    top = "farfield" if a.case == "B" else "slip"
    kinds = {1: ("xmin", "farfield"), 2: ("xmax", "outflow"), 3: ("ymin", "slip"), 4: ("ymax", top), 5: ("zmin", "slip"), 6: ("zmax", "slip")}
    open(os.path.join(run, "bcondConfig.yaml"), "w").write("".join(
        f"{nm}: {{physID: {p}, kind: {k}, outputHDFflg: 0, ints: , floats: {{{fl if k == 'farfield' else ''}}}}}\n" for p, (nm, k) in kinds.items()))
    open(os.path.join(run, "probe.yaml"), "w").write("outStepInterval: 100000\noutStepStart: 0\npoints:\nsurfaces:\n")
    subprocess.run([os.path.join(a.bin, "convertGmshToForge"), "ramp.msh", "ramp.h5"], cwd=run,
                   stdout=open(os.path.join(run, "convert.log"), "w"), stderr=subprocess.STDOUT)
    if not os.path.exists(os.path.join(run, "ramp.h5")) or "writeInputH5: wrote /VIZMESH" not in open(os.path.join(run, "convert.log")).read():
        raise SystemExit(f"変換に失敗: {run}/convert.log")
    e = P / ((gam - 1.0) * ro)
    with h5py.File(os.path.join(run, "ramp.h5"), "r+") as f:
        V = f["VALUE"]; n = V["ro"].shape[0]
        vals = {"ro": ro, "roUx": ro * U, "roUy": 0.0, "roUz": 0.0, "roe": ro * (e + 0.5 * U * U)}
        for kk, v in vals.items():
            V[kk][...] = np.full(n, v, dtype=V[kk].dtype)
    open(os.path.join(run, "IC_FROM.txt"), "w").write(f"V2c case {a.case}: 一様 M {M} (ρ {ro:.6g}、U {U:.6g})、上面 {top}、H {HA + JC * DZ} m\n")
    print(f"prepared {run}: case {a.case}, nodes ≈ {(NX + 1) * (ny + 1) * (NZ + 1)}")


if __name__ == "__main__":
    main()
```

## 参考: `case/58.farfield_verification/v2c_residual_ab.py`

```
#!/usr/bin/env python3
"""V2c 角の履歴差の切り分け (plan boundary-node-farfield-characteristic §5.1 #3a、codex diagnose 2 回目 2026-09-29)。
同じ低領域 B の格子・BC・設定に、2 つの最終場 (一様 IC から = run_0064_v2c_B/res_6000、C から = run_0094_v2c_B_fromC/res_6000) を
index コピーで与え、1 step の初回評価の**絶対残差** (帳簿 res_final = 壁射影等の後、更新に渡る値) を全節点・全 5 保存量で見る。
診断閾値 (事前固定): |R_q| / Σ面|F_q| ≤ 1e-5。分母 = その節点に接する内部面流束の絶対和 (帳簿 .faces) + 境界半割面の規模
(質量 |ρ u·S|、運動量 |ρ u (u·S)| + |P − pRef| |S|、エネルギー |ρ H u·S|; 節点の状態から)。規模 0 は自由流基準 × 節点の境界面積。
これは収束基準の代用ではない。
  python3 v2c_residual_ab.py setup | compare
"""
import csv, os, shutil, subprocess, sys
from collections import defaultdict
import h5py, numpy as np

SRC = "run_0064_v2c_B"
STATES = {"run_0097_v2c_resab_fromU": "run_0064_v2c_B/res_6000.h5", "run_0098_v2c_resab_fromC": "run_0094_v2c_B_fromC/res_6000.h5"}
BIN = os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff/forge")
CONS = ("ro", "roUx", "roUy", "roUz", "roe")
P_REF, GAM, CP = 101325.0, 1.4, 1004.5
RI = P_REF / ((CP - CP / GAM) * 300.0); UI = 2.5 * (GAM * (CP - CP / GAM) * 300.0) ** 0.5
HI = CP / GAM * 300.0 + P_REF / RI + 0.5 * UI * UI   # e = cv T (CPG、setup_v2c と同じ) + P/ρ + U²/2


def setup():
    with h5py.File(os.path.join(SRC, "ramp.h5")) as f:
        n = f["VALUE/ro"].shape[0]
    for dst, snap in STATES.items():
        os.makedirs(dst)
        for f in ("solverConfig.yaml", "bcondConfig.yaml", "ramp.h5", "probe.yaml"):
            shutil.copy(os.path.join(SRC, f), dst)
        s = open(os.path.join(dst, "solverConfig.yaml")).read().replace("nStepOuter: 6000", "nStepOuter: 1").replace("outStepInterval: 1000", "outStepInterval: 1")
        open(os.path.join(dst, "solverConfig.yaml"), "w").write(s)
        with h5py.File(snap) as s5, h5py.File(os.path.join(dst, "ramp.h5"), "r+") as f:
            assert s5["VALUE/ro"].shape[0] == n, "同一格子でない"
            for q in CONS:
                f["VALUE"][q][...] = s5["VALUE"][q][:]
        open(os.path.join(dst, "IC_FROM.txt"), "w").write(f"{snap} の保存量を index コピー (同一格子 {n} 節点)\n")
        env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_DUMP_LEDGER="ledger.csv", FORGE_DUMP_LEDGER_CALLS="1",
                   FORGE_DUMP_LEDGER_NODES="all", FORGE_DUMP_FARFIELD="ffdump")
        r = subprocess.run([BIN], cwd=dst, env=env, stdin=subprocess.DEVNULL, stdout=open(os.path.join(dst, "forge_run.log"), "w"), stderr=subprocess.STDOUT)
        print(dst, "exit", r.returncode)


def compare():
    with h5py.File(os.path.join(SRC, "ramp.h5")) as f:
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        sv = np.array(f["PLANES/surfVect"]).reshape(-1, 3)
        bpl, slipn = {}, {}
        SLIP = {3, 5, 6} | ({4} if "kind: slip" in open(os.path.join(SRC, "bcondConfig.yaml")).read().split("ymax")[1].split("\n")[0] else set())
        for b in f["BCONDS"].keys():
            ic, ip = f[f"BCONDS/{b}/iCells"][:], f[f"BCONDS/{b}/iPlanes"][:]
            for c, p in zip(ic, ip):
                bpl.setdefault(int(c), []).append(sv[int(p)])
                if int(b) in SLIP:
                    slipn.setdefault(int(c), {}).setdefault(int(b), np.zeros(3))
                    slipn[int(c)][int(b)] += sv[int(p)]
    n = len(xyz)
    ok_all = True
    for run in STATES:
        st, res = defaultdict(dict), defaultdict(dict)
        with open(os.path.join(run, "ledger.csv")) as fh:
            rd = csv.reader(fh); next(rd)
            for row in rd:
                if row[0] != "1":
                    continue
                if row[1] == "after_eos_bc":
                    st[row[3]][int(row[2])] = float(row[4])
                elif row[1] == "res_final":
                    res[row[3]][int(row[2])] = float(row[4])
        sc = defaultdict(lambda: np.zeros(n))
        area = np.zeros(n)
        seen = set()
        with open(os.path.join(run, "ledger.csv.faces")) as fh:
            for r in csv.DictReader(fh):
                if r["call"] != "1" or r["ip"] in seen:
                    continue
                seen.add(r["ip"])
                # 帳簿の F_roU* は圧力項を含まない (接線流の y・z 運動量で分母が 0 に近くなる) ので、面の圧力 × 面積成分を足す
                pf = abs(0.5 * (float(r["P_L"]) + float(r["P_R"])) - P_REF)
                sd = {"roUx": abs(float(r["sx"])), "roUy": abs(float(r["sy"])), "roUz": abs(float(r["sz"]))}
                for c, q in (("F_ro", "ro"), ("F_roUx", "roUx"), ("F_roUy", "roUy"), ("F_roUz", "roUz"), ("F_roe", "roe")):
                    v = abs(float(r[c])) + (pf * sd[q] if q in sd else 0.0)
                    sc[q][int(r["ic0"])] += v; sc[q][int(r["ic1"])] += v
                a = float(r["ss"]); area[int(r["ic0"])] += a; area[int(r["ic1"])] += a
        for c, planes in bpl.items():
            ro = st["ro"][c]; u = np.array([st["Ux"][c], st["Uy"][c], st["Uz"][c]]); P = st["P"][c]
            H = (st["roe"][c] + P) / ro
            for S in planes:
                un = float(u @ S); A = float(np.linalg.norm(S))
                sc["ro"][c] += abs(ro * un); sc["roe"][c] += abs(ro * H * un)
                for d, q in enumerate(("roUx", "roUy", "roUz")):
                    sc[q][c] += abs(ro * u[d] * un) + abs(P - P_REF) * abs(S[d])
        print(f"{run}: 内部面 {len(seen)}、境界半割面を持つ節点 {len(bpl)} (うち slip {len(slipn)}、運動量は接線成分で判定)")
        # slip 面を持つ節点は運動量残差から壁法線成分 (面ごとの法線、Gram–Schmidt) を除く: 更新時に射影で捨てられる成分なので収束の指標でない
        Rm = np.stack([np.array([res["res_" + q].get(i, np.nan) for i in range(n)]) for q in ("roUx", "roUy", "roUz")], 1)
        for c, dct in slipn.items():
            basis = []
            for v in dct.values():
                w = v / np.linalg.norm(v)
                for b in basis:
                    w = w - (w @ b) * b
                if np.linalg.norm(w) > 1e-6:
                    basis.append(w / np.linalg.norm(w))
            for b in basis:
                Rm[c] = Rm[c] - (Rm[c] @ b) * b
        Rproj = {"roUx": Rm[:, 0], "roUy": Rm[:, 1], "roUz": Rm[:, 2]}
        ok = True
        for q in CONS:
            R = Rproj[q] if q in Rproj else np.array([res["res_" + q].get(i, np.nan) for i in range(n)])
            # 規模が数値的に 0 の成分 (一様流域の y・z 運動量など) は自由流基準 × 節点の面積 (単位を合わせる)
            ref = {"ro": RI * UI, "roe": RI * UI * HI}.get(q, RI * UI * UI) * area
            # 運動量 3 成分は共通の規模 (面流束 3 成分の絶対和、V1d と同じ): 成分別では流れに垂直な成分の分母が痩せる
            s = (sc["roUx"] + sc["roUy"] + sc["roUz"]).copy() if q in ("roUx", "roUy", "roUz") else sc[q].copy()
            zero = s < 1e-6 * ref; s[zero] = ref[zero]
            r = np.abs(R) / s
            w = int(np.nanargmax(r))
            corner = (xyz[:, 0] >= 0.195) & (xyz[:, 0] <= 0.23) & (xyz[:, 1] <= 0.01)
            exitb = (xyz[:, 0] >= 1.195) & (xyz[:, 1] <= 0.09)
            nbad = int(np.sum(r > 1e-5))
            ok &= nbad == 0 and np.all(np.isfinite(R))
            top = np.argsort(-np.nan_to_num(r))[:3]
            w = int(np.nanargmax(r))
            print(f"  {q:5s}: max |R|/Σ|F| {r[w]:.3e} (x {xyz[w,0]:.3f} y {xyz[w,1]:.4f})、> 1e-5 の節点 {nbad}、角 {np.nanmax(r[corner]):.2e}、出口下端 {np.nanmax(r[exitb]):.2e}、非有限 {int(np.sum(~np.isfinite(R)))}、上位 {[(round(float(xyz[t,0]),3), round(float(xyz[t,1]),4), round(float(xyz[t,2]),3), f'{r[t]:.1e}') for t in top]}")
        print(f"  → {'全項目 ≤ 1e-5' if ok else '閾値超えあり'}")
        ok_all &= ok
    print("VERDICT: " + ("両状態とも全節点で ≤ 1e-5 → 残差を残した停止の仮説を棄却 (異なる離散平衡の候補として扱う)" if ok_all
                         else "少なくとも片方が閾値超え → 「二つとも離散定常解」を棄却し、非零残差が更新されない箇所を追う"))


if __name__ == "__main__":
    {"setup": setup, "compare": compare}[sys.argv[1]]()
```

## 参考: `case/58.farfield_verification/setup_v2d.py`

```
#!/usr/bin/env python3
"""V2d TP での境界の独立参照 (plan boundary-node-farfield-characteristic §6 V2d) の forge run を作る (AWS で実行)。
V2a と同じ薄板チャネル (Δx 5 mm、y・z 各 3 層 slip、左右端 farfield)、dual-time。自由流 = SERN 外気 (220 K、2851 Pa、Y_EXH 0)。
  --gas cpg | tp1 (単成分 AIR、NASA-9) | tp2 (EXH/AIR)
  --test contact : V2d-1。P・u 一様 (P∞、u = ±0.5 c∞)、中心 0.5 m (長領域も同じ x) にガウス (半値全幅 0.1 m) の T 600 K・Y_EXH 0.13 の塊。
                   塊は右端 (u > 0) から出ていく。評価点 = 右端から 5 セル内側 (x = 0.975 m)。場を 10 step ごとに書く (T・Y の時系列)。
  --test acoustic: V2d-2。内部 = 一様の T 600 K・Y_EXH 0.13 (P∞、u = 0.3 c_i)、中心 0.4 m に右向きの圧力パルス (振幅 1e-3 P∞、FWHM 40 mm)。
                   farfield の自由流は外気 (内外で γ・ρc が違う)。評価点 x = 0.8 m。--nopulse でパルスなし対照。
リミッタ基準値は自由流で固定 (V2a と同じ理由)。
  python3 setup_v2d.py RUN --gas G --test T --len 1|3 [--dir 1|-1] [--nopulse] [--cfl-ac 0.113] [--nsub 20]
"""
import argparse, math, os, shutil, subprocess, sys
import h5py, numpy as np, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from setup_v1 import SPECIES_SRC, RU, nasa9_h_mass  # noqa: E402

P0, T0 = 2851.0, 220.0


class Gas:
    def __init__(self, kind, db):
        self.kind, self.db = kind, db

    def R(self, Y):
        if self.kind == "cpg":
            return 1004.5 - 1004.5 / 1.4
        if self.kind == "tp1":
            return RU / self.db["AIR"]["MW"]
        return (Y / self.db["EXH"]["MW"] + (1 - Y) / self.db["AIR"]["MW"]) * RU

    def h(self, Y, T):   # 顕熱 (datum 298.15 K、setup_v1 と同じ規約)
        if self.kind == "cpg":
            return 1004.5 * T
        sp = {"AIR": 1.0} if self.kind == "tp1" else {"EXH": Y, "AIR": 1 - Y}
        return sum(w * (nasa9_h_mass(self.db[s], T) - nasa9_h_mass(self.db[s], 298.15)) for s, w in sp.items())

    def e(self, Y, T):
        return self.h(Y, T) - self.R(Y) * T if self.kind != "cpg" else 1004.5 / 1.4 * T

    def c(self, Y, T):
        if self.kind == "cpg":
            return math.sqrt(1.4 * self.R(Y) * T)
        cp = (self.h(Y, T + 0.01) - self.h(Y, T - 0.01)) / 0.02
        R = self.R(Y)
        return math.sqrt(cp / (cp - R) * R * T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--gas", choices=("cpg", "tp1", "tp2"), required=True)
    ap.add_argument("--test", choices=("contact", "acoustic"), required=True); ap.add_argument("--len", type=float, required=True)
    ap.add_argument("--dir", type=int, default=1); ap.add_argument("--nopulse", action="store_true")
    ap.add_argument("--left-hot", action="store_true", help="診断: 左端 farfield の自由流を内部と同じ高温状態にする (左からの接触面を作らない)")
    ap.add_argument("--cfl-ac", type=float, default=0.113); ap.add_argument("--nsub", type=int, default=20)
    ap.add_argument("--bin", default=os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff"))
    a = ap.parse_args()
    run = a.run
    os.makedirs(run)
    for f in ("species_db.yaml", "species_meta.yaml"):
        shutil.copy(os.path.join(SPECIES_SRC, f), run)
    db = yaml.safe_load(open(os.path.join(run, "species_db.yaml")))
    g = Gas(a.gas, db)
    dx, W, nyz = 0.005, 0.015, 3
    nx = int(round(a.len / dx))
    subprocess.run([sys.executable, os.path.join(HERE, "make_box_msh.py"), os.path.join(run, "chan.msh"), str(nx), str(nyz), str(nyz),
                    repr(a.len), repr(W), repr(W)], check=True)
    ro0 = P0 / (g.R(0.0) * T0); c0 = g.c(0.0, T0)
    Ti, Yi = 600.0, (0.13 if a.gas == "tp2" else 0.0)
    ci = g.c(Yi, Ti)
    if a.test == "contact":
        U = a.dir * 0.5 * c0; umax = abs(U) + max(c0, ci); tend = (0.475 + 0.15) / abs(U); probes = {"eval": 0.975 if a.dir > 0 else 0.025}
    else:
        U = 0.3 * ci; umax = U + ci; tend = 0.4 / (U + ci) + 0.2 / (ci - U) + 0.15 / (ci - U); probes = {"eval": 0.8, "mid": 0.4}
    dt = a.cfl_ac * dx / umax
    steps = int(math.ceil(tend / dt))
    if a.gas == "cpg":
        phys = "{thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}"
    else:
        spl = '["AIR"]' if a.gas == "tp1" else '["EXH", "AIR"]'
        phys = f'{{thermalMethod: 2, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4, pMin: 20.0, species: {spl}, speciesDBFile: "species_db.yaml", thermoHrefTemp: 298.15}}'
    outint = 10 if a.test == "contact" else steps
    open(os.path.join(run, "solverConfig.yaml"), "w").write(f"""mesh: {{discretization: "node", nodeWallDirichlet: 1, meshFileName: "chan.h5", valueFileName: "chan.h5"}}
gpu: 1
solver: "SLAU"
physProp: {phys}
time:
  unsteady: 1
  dualTime: 1
  nSubIterDualTime: {a.nsub}
  last: {{nStepOuter: {steps}}}
  deltaT: {{control: 0, dt: {dt!r}, cfl: 1.0, cfl_pseudo: 12.0, dt_min: 1e-12, dt_max: 1.0, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1}}
  outStepStart: 0
  outStepInterval: {outint}
  timeIntegration: 11
  nStepInner: 5
space: {{convMethod: 1, limiter: 2, pRef: {P0}, limiterRefLength: 1.0, limiterRoRef: {ro0!r}, limiterPRef: {P0!r}, limiterARef: {c0!r}}}
turbulence: {{model: "none"}}
initial: "uniform_p101325_u10"
""")
    fl = f"ro: {ro0!r}, Ux: {U!r}, Uy: 0.0, Uz: 0.0, Ps: {P0!r}" + (", Y0: 0.0, Y1: 1.0" if a.gas == "tp2" else "")
    kinds = {1: ("xmin", "farfield"), 2: ("xmax", "farfield"), 3: ("ymin", "slip"), 4: ("ymax", "slip"), 5: ("zmin", "slip"), 6: ("zmax", "slip")}
    flh = fl
    if a.left_hot:
        roh = P0 / (g.R(Yi) * Ti)
        flh = f"ro: {roh!r}, Ux: {U!r}, Uy: 0.0, Uz: 0.0, Ps: {P0!r}" + (f", Y0: {Yi!r}, Y1: {1 - Yi!r}" if a.gas == "tp2" else "")
    open(os.path.join(run, "bcondConfig.yaml"), "w").write("".join(
        f"{nm}: {{physID: {p}, kind: {k}, outputHDFflg: 0, ints: , floats: {{{(flh if p == 1 else fl) if k == 'farfield' else ''}}}}}\n" for p, (nm, k) in kinds.items()))
    yp = W / 3.0 + 0.001
    open(os.path.join(run, "probe.yaml"), "w").write("outStepInterval: 1\noutStepStart: 0\npoints:\n" + "".join(
        f"  {k}: {{x: {v + 0.001}, y: {yp}, z: {yp}}}\n" for k, v in probes.items()) + "surfaces:\n")
    subprocess.run([os.path.join(a.bin, "convertGmshToForge"), "chan.msh", "chan.h5"], cwd=run,
                   stdout=open(os.path.join(run, "convert.log"), "w"), stderr=subprocess.STDOUT)
    if not os.path.exists(os.path.join(run, "chan.h5")) or "writeInputH5: wrote /VIZMESH" not in open(os.path.join(run, "convert.log")).read():
        raise SystemExit(f"変換に失敗: {run}/convert.log")
    with h5py.File(os.path.join(run, "chan.h5"), "r+") as f:
        V = f["VALUE"]; n = V["ro"].shape[0]
        x = np.array(f["MESH/COORD"]).reshape(-1, 3)[:, 0].astype(float)
        if a.test == "contact":
            sig = 0.1 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            w = np.exp(-0.5 * ((x - 0.5) / sig) ** 2)
            T = T0 + (Ti - T0) * w; Y = Yi * w; P = np.full(n, P0); u = np.full(n, U)
        else:
            sig = 0.04 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            roi = P0 / (g.R(Yi) * Ti)
            pp = 0.0 if a.nopulse else 1e-3 * P0 * np.exp(-0.5 * ((x - 0.4) / sig) ** 2)
            P = P0 + pp * np.ones(n); u = U + pp / (roi * ci); Y = np.full(n, Yi)
            # 等エントロピー擾乱: T を P から (小振幅、線形)
            gi = ci * ci * roi / P0
            T = Ti * (P / P0) ** ((gi - 1.0) / gi)
        R = np.array([g.R(y) for y in Y]); ro = P / (R * T)
        e = np.array([g.e(y, t) for y, t in zip(Y, T)])
        vals = {"ro": ro, "roUx": ro * u, "roUy": 0.0 * ro, "roUz": 0.0 * ro, "roe": ro * (e + 0.5 * u * u)}
        if a.gas == "tp2":
            vals.update({"roY0": ro * Y, "roY1": ro * (1 - Y)})
        for kk, v in vals.items():
            if kk in V:
                V[kk][...] = v.astype(V[kk].dtype)
            else:
                V.create_dataset(kk, data=v.astype(V["ro"].dtype))
    open(os.path.join(run, "IC_FROM.txt"), "w").write(
        f"V2d {a.test} gas {a.gas}: left_hot {a.left_hot}、L {a.len} m、U {U:.5g} (dir {a.dir})、nopulse {a.nopulse}、dt {dt:.4e}、steps {steps}、"
        f"自由流 ρ {ro0:.6g} c {c0:.5g}、内部 T {Ti} Y {Yi} c {ci:.5g}\n")
    print(f"prepared {run}: steps {steps}, dt {dt:.3e}")


if __name__ == "__main__":
    main()
```

## 参考: `case/58.farfield_verification/eval_v2d.py`

```
#!/usr/bin/env python3
"""V2d の判定 (plan boundary-node-farfield-characteristic §6 V2d、判定 (b) = 短領域 vs 長領域 = 境界が加える誤差)。
  contact SHORT LONG : V2d-1。評価点 (右端から 5 セル内側) で
       max_t |P_short − P∞|/P∞ ≤ 1e-3、 max_t |T_short − T_long| ≤ 0.02 × (600 − 220) K、 max_t |Y_short − Y_long| ≤ 0.002
       (Y は 10 step ごとの res_*.h5 から評価点の節点値を読む。T・P はプローブ)
  nopulse SHORT      : V2d-2 パルスなし対照。評価点 (x 0.8 m) で max_t |P − P∞| ≤ 0.01 × 1e-3 P∞
  acoustic SHORT LONG: V2d-2。反射 = max_t |P_short − P_long| / 入射振幅 ≤ 0.05 (入射振幅 = max_t |P_long − P∞|)
"""
import glob, re, sys
import h5py, numpy as np

P0 = 2851.0


def probe(run, i=0):
    d = np.genfromtxt(f"{run}/point_probe_{i}.out", delimiter=",", names=True)
    return np.asarray(d["TotalTime"], float), np.asarray(d["P"], float), np.asarray(d["T"], float)


def y_series(run, xeval):
    fs = sorted((int(re.search(r"res_(\d+)\.h5$", p).group(1)), p) for p in glob.glob(run + "/res_*.h5") if re.search(r"/res_\d+\.h5$", p))
    with h5py.File(run + "/chan.h5") as m:
        xyz = np.array(m["MESH/COORD"]).reshape(-1, 3)
    node = int(np.argmin((xyz[:, 0] - xeval) ** 2 + (xyz[:, 1] - 0.005) ** 2 + (xyz[:, 2] - 0.005) ** 2))
    t, Y = [], []
    for st, p in fs:
        with h5py.File(p) as f:
            V = f["VALUE"]
            if "roY0" not in V:
                return None
            Y.append(float(V["roY0"][node] / V["ro"][node]))
            t.append(st)
    return np.array(t), np.array(Y)


mode = sys.argv[1]
if mode == "contact":
    s, l = sys.argv[2], sys.argv[3]
    ts, Ps, Ts = probe(s); tl, Pl, Tl = probe(l)
    tm = min(ts[-1], tl[-1]); k = ts <= tm
    dP = np.max(np.abs(Ps[k] - P0)) / P0
    dT = np.max(np.abs(Ts[k] - np.interp(ts[k], tl, Tl)))
    Tamp = np.max(Tl) - 220.0
    out = [f"max|P−P∞|/P∞ {dP:.3e} (≤ 1e-3)", f"max|ΔT| {dT:.3f} K (≤ {0.02 * 380:.1f} K、塊の到達振幅 {Tamp:.1f} K)"]
    ok = dP <= 1e-3 and dT <= 0.02 * 380.0
    ys, yl = y_series(s, 0.975), y_series(l, 0.975)
    if ys is not None:
        n = min(len(ys[0]), len(yl[0]))
        assert np.array_equal(ys[0][:n], yl[0][:n])
        dY = np.max(np.abs(ys[1][:n] - yl[1][:n]))
        out.append(f"max|ΔY| {dY:.2e} (≤ 0.002、塊の到達 Y {np.max(yl[1]):.4f})")
        ok = ok and dY <= 0.002
    print(f"{s} vs {l}: " + "、".join(out))
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
elif mode == "nopulse":
    s = sys.argv[2]
    ts, Ps, _ = probe(s)
    d = np.max(np.abs(Ps - P0))
    print(f"{s}: パルスなし max|P − P∞| {d:.4e} Pa (≤ {0.01 * 1e-3 * P0:.4e})")
    print(f"VERDICT: {'PASS' if d <= 0.01 * 1e-3 * P0 else 'FAIL'}")
elif mode == "acoustic":
    s, l = sys.argv[2], sys.argv[3]
    ts, Ps, _ = probe(s); tl, Pl, _ = probe(l)
    tm = min(ts[-1], tl[-1]); k = ts <= tm
    A = np.max(np.abs(Pl[tl <= tm] - P0))
    d = np.abs(Ps[k] - np.interp(ts[k], tl, Pl))
    print(f"{s} vs {l}: 入射振幅 {A:.4g} Pa、max|ΔP| {d.max():.4g} Pa (t {ts[k][d.argmax()]:.4e} s) → 反射率 {d.max() / A:.4%}")
    print(f"VERDICT: {'PASS' if d.max() / A <= 0.05 else 'FAIL'} (≤ 5 %)")
```

## 参考: `notes/reviews/2026-09-29-farfield-v2c-corner-and-v2d2-redesign-diagnose.md`

```
# codex 諮問 (diagnose): farfield-v2c-corner-and-v2d2-redesign

- **brief**: [`notes/reviews/briefs/2026-09-29-farfield-v2c-corner-and-v2d2-redesign.md`](../../notes/reviews/briefs/2026-09-29-farfield-v2c-corner-and-v2d2-redesign.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-09-29
- **commit**: `71442b34` (feature/sern-design)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **次は同じ低領域 B の格子・BC に、`run_0064_v2c_B` と `run_0094_v2c_B_fromC` の最終保存量をそれぞれ与え、初回の絶対残差を比較する A/B で、「残差を残した停止」と「異なる離散定常解の候補」を分ける。**
- **extra**: `case/58.farfield_verification/v2c_operator_ab.py`, `case/58.farfield_verification/eval_v2c.py`, `case/58.farfield_verification/setup_v2d.py`, `notes/reviews/2026-09-29-farfield-v2c-v2d-fail-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**Major 4 件、Minor 1 件**。run 数値はブリーフ記載値であり、原本による独立検証はできていない）。

| 重大度 | 諮問事項 | 採否・根拠・対案 |
|---|---|---|
| **Major** | V2c を C 共通初期場で判定し直す | **追加の隔離試験として採用。旧試験の訂正・置換は要再検証**。`fromC()` は保存量を座標対応で移しており、初期場を揃える操作として妥当（[v2c_operator_ab.py:145](/home/sano/work/forge-sern-design/case/58.farfield_verification/v2c_operator_ab.py:145)）。ただし「異なる離散定常解がある」は未証明。作用素 A/B が測ったのは**残差の差**であり、各状態の残差がゼロに近いかではない（同ファイル:121）。対案は旧 FAIL を保持し、「履歴依存を制御する追加試験」として §6 に事前登録すること。全線・0.02 の基準は維持する。 |
| **Major** | 全線 0.0210 を「反射約 2 %」と記録する／出口端を避ける／境界式を変更する | **現段階ではいずれも却下**。評価器は最終スナップショットの壁圧差を測るだけで、反射波を分離していない（[eval_v2c.py:23](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2c.py:23)、同:43）。記録できるのは「共通初期場試験でも、最終場の全線誤差が基準超過」。超過幅は `0.001 Δp ≈ 87.5 Pa` であり、この幅より十分小さい時間変動の確認が必要。対案は下記の残差 A/B を先行させること。出口延長は後続の原因切り分けにはなり得るが、反射到達前で評価線を切る案は検出対象そのものを外す。 |
| **Major** | `--left-hot` を正式な V2d-2 にする | **右端流出の隔離試験として採用。合格判定は要再検証**。コードは左端だけ高温にし、右端の冷たい自由流を保つため、内外物性差を持つ右端の検証目的は残る（[setup_v2d.py:102](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2d.py:102)）。旧配置の FAIL は別記する。対案は短・長の双方でパルスなし対照を通し、CPG／単成分 TP／多成分 TP について、入射窓で振幅を測り、反射到達窓で短長差を測ること。現評価器の全時間最大値による分母・分子は修正が必要（[eval_v2d.py:61](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2d.py:61)）。独立参照の精度確認も §6 のまま残す。 |
| **Major** | V2 未受入れのまま V3 に進む | **却下**。V3 は V1–V2 合格後という明示条件がある（[plan:187](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:187)）。対案は V2c の定常性・誤差帰属と、再設計した V2d-2 の受入れを先に終えること。幅系列 4 本ではこの未解決事項を代替できない。 |
| **Minor** | 作用素比較器の `VERDICT` を単独で根拠にする | **要修正**。幾何差・状態差・未対応面数は表示されるが、`ok` に入るのは残差差だけ。ゼロ規模も宣言された自由流基準ではなく `1.0` になる（[v2c_operator_ab.py:86](/home/sano/work/forge-sern-design/case/58.farfield_verification/v2c_operator_ab.py:86)、同:103、同:127）。対案は全項目・対応の完全性・有限性を判定に含めること。今回報告された「全項目差 0、未対応なし」が原本で確認できれば、その観測まで否定する理由はない。 |

結論: **次は同じ低領域 B の格子・BC に、`run_0064_v2c_B` と `run_0094_v2c_B_fromC` の最終保存量をそれぞれ与え、初回の絶対残差を比較する A/B で、「残差を残した停止」と「異なる離散定常解の候補」を分ける。**

第 1 仮説: **角の履歴差は、少なくとも一方が非零の離散残差を残す反復停滞であり、二つの定常解を示していない。** 確度: **低**
  
  根拠: ブリーフでは旧系列の残差低下が 2.4–2.5 桁でプラトー、新系列は未確認。`case/58.farfield_verification/run_0092_v2c_opab_low`／`run_0093_v2c_opab_high` の報告値 `|ΔR|/Σ|F| ≤ 7.7e-8` は、両者の `|R|` を制限しない。`eval_v2c.py:25` も最終場しか読まない。**停滞の直接証拠は未確認**であり、優先する理由は二解説を安く反証できるため。

  反証条件: 同一 B 作用素で、異なる両状態とも全保存量の局所残差が下記の診断閾値内なら、「閾値を超える残差を残した停止」という仮説を棄却する。ただし、それだけで厳密な二解・収束を認定しない。

第 2 仮説: **共通初期場でも残る下流壁圧差には、上面 `farfield` の斜入射への応答が含まれる。** 確度: **中**。現仕様も斜入射の限界を明記する（[methods/boundary.md:150](/home/sano/work/forge-sern-design/methods/boundary.md:150)）。ただし **0.0210 の出口処理への帰属、反射率への読み替えは未確認**。

第 3 仮説: **旧 V2d-2 の擾乱源は、左端から入る接触面の TP 保存形混合である。** 確度: **中**。`--left-hot` が接触面を消す操作であることはコードで確認したが、境界固有の寄与は未除外。周期内の段差から「同じ桁の 2.5 Pa」が出るだけでは原因帰属できない。周期配置に伴うもう一つの接触面との干渉、混合率・局所 EOS 復元圧力・波の発生位置を区別する必要がある。

判別 A/B: **変えるのは初期保存量だけ。各 1 step、判定は更新前の初回評価。**

- 格子・BC・数値設定・バイナリ・CUDA block size を `run_0064_v2c_B` と同一にする。
- A 入力＝`run_0064_v2c_B/res_6000.h5`、B 入力＝`run_0094_v2c_B_fromC/res_6000.h5`。いずれも `case/58.farfield_verification/` 配下。同一格子であることを確認して index コピーする。
- 全節点・全 5 保存量の**更新に渡る残差**を比較する。帳簿の `res_final` は壁射影等の後に取得される（[main.cpp:1555](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1555)）。角・出口交線の値も別集計する。
- 診断閾値は事前に `|R_q| / Σ面|F_q| ≤ 1e-5`。分母には境界面も含め、ゼロ規模には単位の合う非零の自由流基準を使う。**これは今回の診断閾値であり、収束基準の代用ではない。**

→ **少なくとも片方が超えるなら**「二つとも離散定常解」という説明を棄却し、非零残差が更新されない箇所を追う。  
→ **両方とも全項目で許容内なら**第 1 仮説をその精度で棄却し、異なる離散平衡の候補として履歴依存を扱う。既存系列の `check_convergence.py` と壁圧時系列の `check_quasisteady.py --series-csv` の判定は、別途必要。

やらない方がよいこと: **出口端の除外で合格化する、0.0210 に合わせて許容値を緩める、壁圧差を反射係数と呼ぶ、未確認の C を収束床に指定すること。** `--from-floor` は参照 run 自身の通常判定 PASS を要求する（[check_convergence.py:359](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:359)）。

呼び出し側の前提への異議: 観測は「初期場を変えると、6000 step 後の角圧力が別の値を保持した」。解釈は「少なくとも二つの定常解」「領域高さが解を選ぶ」「出口端差は反射と出口処理の重なり」であり、後者はまだ確定できない。局所作用素差の除外は、報告された入力状態・対象領域・精度に限定して受け入れ、同じ試験を繰り返す提案はしない。**plan 未反映（依頼どおり変更なし）。**

不足情報: 手元の対象 `run_*` は **0 件**（[case README:5](/home/sano/work/forge-sern-design/case/58.farfield_verification/README.md:5)）。AWS の実効 YAML、バイナリ識別、帳簿原本、全残差履歴、保存場時系列、メッシュ品質・収束・準定常性の VERDICT と判定区間が必要。したがって今回、**run 数値の再現確認も収束認定も行っていない**。
```

## 出力形式 (この形のまま)

```
結論: <次にやる一手を 1 文で>
第 1 仮説: <内容>  確度: <高/中/低>
  根拠: <ファイル:行 / run パスと数値>
  反証条件: <何が観測されたらこの仮説は誤りか>
第 2・第 3 仮説: <あれば 1 行ずつ>
判別 A/B: <変える設定 1 点、回す長さ、見る量>  → A なら … / B なら …
やらない方がよいこと: <呼び出し側が取りそうな誤った一手>
呼び出し側の前提への異議: <ブリーフの枠組み・除外判断・指標の定義で受け入れなかったものと理由。無ければ「無し」>
不足情報: <あれば>
```
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
