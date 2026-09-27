# 特性型の遠方境界 `farfield` (node)

## メタ

- **area**: `boundary`
- **status**: `draft`
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

- **やる**: node・密度ベース・SLAU/ROE 共通の境界半割面流束に渡す境界状態 (`bvar`) を作るカーネル、CPG と TP (単一・多成分)、
  SST の $k,\omega$、化学種、陽解法・block-DPLUR。runner の `side_far` / `top_out` を切替可能にする。
- **やらない**: cell 方式 (使わない方針)、凝縮・トレーサ・遷移モデルとの併用 (起動時に拒否)、動的格子、
  粘性の遠方境界流束 (現行の非壁境界と同じく 0)、陰解法の境界 Jacobian の厳密化 (現行の ghostless 対角のまま)。

## 3. 関連 docs と前提

- 理論は [`methods/boundary.md`](../../methods/boundary.md) 「特性型の遠方境界」。参照実装は SU2 `CEulerSolver::BC_Far_Field`
  (`.external/su2-src/SU2_CFD/src/solvers/CEulerSolver.cpp:4822`)。
- node の境界流束は種別によらず `bvar` (rob, Uxb/Uyb/Uzb, roeb, Psb) を R 状態として組む
  (`convectiveFlux_boundary_d.inc.cuh:113-118`、質量流束も R 状態から 193 行)。したがって本 BC は `bvar` を整合的に埋めればよい。
- 既存の `outlet_statPress` が同じ特性構成 (外向き Riemann + 内部エントロピー、TP は `gamma_cell[ic]`) を持つ (`boundaryCond_d.cu:576-629`)。
- 化学種・受動種・トレーサ・凝縮の境界は名前が `inlet_` で始まるかで分岐し、RANS は明示リスト (`ransBoundary_d.cu:254-290`、
  未知の種別は**何もしない**)。

## 4. 設計方針

**改訂 (2026-09-27、codex plan 段 NO-GO C1/M6/m1 を全件採用)**: 「境界状態 `bvar` を埋めれば既存の境界流束がそのまま使える」という初稿の前提は誤り
(下の 4.0)。`farfield` は **同一評価時点の境界状態から、保存量 5 成分・化学種・k/ω の面流束をすべて自前で組む** node 専用の境界にする。

### 4.0 既存経路の事実 (コードで確認、2026-09-27)

- 境界流束 `convectiveFlux_boundary_d` は、質量流束は境界状態 (R) から作るが、**流出時の運動量・エネルギーは内部速度・内部全エンタルピーで風上化**する
  (`convectiveFlux_boundary_d.inc.cuh:192` 以降)。構成状態の物理流束 $F(U_b)$ ではない (codex の数値例: エネルギー流束で −4.45 %)。
- node 境界半割面のスカラー移流は、流入・流出によらず**内部節点の値**を使う (`scalarTransport_d.cu:166` の `ext_is_self`、`passiveKernels_d.cuh` の `nodeBnd`)。
  既存の入口で組成・k/ω が入るのは節点ピンのおかげ。→ ゴーストを埋めても外の値は運ばれない。
- スカラー境界処理 (`applySpeciesBoundaries` 等) は対流流束より前に呼ばれる (`main.cpp:1458`) ので、そこで `massflux` を読むと前回評価の値になる。
- SST 全エネルギー (`sstEnergyIncludesK`) の境界 k は `inlet` 接頭辞のときだけ bvar、他は内部値 (`convectiveFlux_d.cu:413`)。

### 4.1 入力

`bcondConfig.yaml`: `{kind: farfield, floats: {ro, Ux, Uy, Uz, Ps, k, omega, Y0, ...}}` = 自由流 (`inlet_uniformVelocity` と同じキー)。多成分は `Y{s}`/`X{s}` 必須、RANS は `k`, `omega` 必須。
`valueTypesOfBC["farfield"]` を新設し、自由流の値 (type 1) と組み立てた境界状態 (別名の bvar) を分けて持つ (自由流を上書きしない)。

### 4.2 境界状態の構成 (frozen-γ 契約)

1 境界半割面 = 1 スレッド、`ic` = 境界節点 (内部状態)、$\hat{\mathbf n}$ = 外向き単位法線。

- **γ の契約**: 面ごとに**単一の $\gamma^*$ = 内部の $\gamma_i$** (CPG は `cfg.gamma`、TP は `gamma_cell[ic]` = frozen $\gamma_{mix}$) を使い、**音速・Riemann 不変量・エントロピー・$\rho,P$ の復元をすべて $\gamma^*$ で**行う:
  $c_i=\sqrt{\gamma^*P_i/\rho_i}$、$c_\infty=\sqrt{\gamma^*P_\infty/\rho_\infty}$ (自由流の音速も $\gamma^*$ で作る)、$s=P/\rho^{\gamma^*}$。
  熱力学的な内部エネルギーだけは TP の実物性で作る ($T_b=P_b/(\rho_b R_{mix}(Y_b))$ → `thermo_state_at_T`)。
  この近似の誤差は V1 では検出できない (一様流では $\gamma_i=\gamma_\infty$) ので、**V2b (内外の温度・組成が違う流入/流出)** で許容を決める。
  SERN 側方の実測範囲 (run_0986 の遠方面 51,143 節点、瞬時値): $Y_{EXH}$ 0–0.128、T 176–597 K、$\gamma$ 1.372–1.405。
- **分岐の決定表** (評価状態を明記、優先順に上から):

| 条件 (評価に使う量) | 境界状態 |
| --- | --- |
| 内部が非有限・$\rho_i\le0$・$P_i\le0$ | 起動中なら拒否、走行中は診断カウンタ + 自由流状態 (面数をログ) |
| $M_{n,i}=U_{n,i}/c_i \ge 1$ **かつ** $M_{n,\infty}=U_{n,\infty}/c_\infty \ge 1$ | 超音速流出: 全量内部 |
| $M_{n,i}\le -1$ **かつ** $M_{n,\infty}\le -1$ | 超音速流入: 全量自由流 |
| それ以外 (分類が食い違う場合を含む) | 亜音速の式: $R^+=U_{n,i}+2c_i/(\gamma^*-1)$、$R^-=U_{n,\infty}-2c_\infty/(\gamma^*-1)$、$U_{n,b}=(R^++R^-)/2$、$c_b=(\gamma^*-1)(R^+-R^-)/4$ |
| 亜音速で $R^+-R^-\le 0$ ($c_b\le0$、大膨張) | 診断カウンタ + $c_b=\max(c_b, 10^{-3}\min(c_i,c_\infty))$ ではなく**全量内部に落とす** (音速を捏造しない) |

  亜音速の続き: **流入/流出は $U_{n,b}$ の符号**で決める (SU2 は $U_{n,\infty}$。側方で $U_{n,\infty}=0$ だと常に流入扱いになり、境界に達したプルームのエントロピー・組成を自由流で置換する)。
  流出 ($U_{n,b}>0$): $s=s_i$、接線速度・組成・k/ω は内部。流入 ($U_{n,b}\le0$): $s=s_\infty$、接線速度・組成・k/ω は自由流。
  $\rho_b=(c_b^2/(\gamma^*s))^{1/(\gamma^*-1)}$、$P_b=\rho_bc_b^2/\gamma^*$。
- 単体試験 (ホスト側で同じ関数を呼ぶ): codex の反例 ($\gamma$ 1.4、$c_i=c_\infty=1$、$U_{n,i}=2$、$U_{n,\infty}=0.5$ → 亜音速式 $U_{n,b}=1.25$、$c_b=1.15$)、
  音速通過 ($M_{n,i}$ 0.99/1.01)、流向反転 ($U_{n,b}$ ±ε)、$c_b\le0$。

### 4.3 面流束 (farfield 専用分岐)

- 構成状態 $U_b=(\rho_b,\mathbf u_b,P_b,H_b,Y_b,k_b,\omega_b)$ の**物理流束**を面流束にする:
  $\dot m=\rho_bU_{n,b}|S|$、運動量 $\dot m\,\mathbf u_b + (P_b-p_{ref})\mathbf S$、エネルギー $\dot m H_b$、化学種 $\dot m Y_b$、$k,\omega$ も $\dot m k_b$ 等。
  (上の分岐で「流出なら内部の値、流入なら自由流の値」になっているので、これが風上化に当たる。)
- 実装: `convectiveFlux_boundary_d` に farfield 分岐を足すのではなく、**`farfield_flux_d` を別カーネル**にする (既存境界のコードを触らない = V0)。
  同じカーネルが `massflux[ip]` と、スカラー用の面値 (面ごとの $Y_b$, $k_b$, $\omega_b$) を書く。
- スカラー移流: node 境界半割面の `ext_is_self` / `nodeBnd` 経路に、**farfield 面だけ面値配列を読む**分岐を足す
  (面フラグ配列で判定。既存の面は従来どおり内部値 = ビット不変)。一次輸送 (`scalarTransport_d.cu`) と S3 経路 (`passiveKernels_d.cuh`、化学種) の両方。
  呼び出し順は対流流束 → スカラー移流なので、同じ評価時点の $\dot m$ と面値を使う (境界処理での `massflux` 読みはしない)。
- `sstEnergyIncludesK`: エネルギー流束の $H^*$ と圧力の $p^*$ に、**スカラー流束と同じ $k_b$** を使う。
- ピン (`scalarDirichletPin`) はしない。凝縮・トレーサ・受動種・遷移モデル・軸対称・周期との併用は初版では起動時に拒否。

### 4.4 その他

- ディスパッチ・読込 (`boundaryCond.{hpp,cpp}`): `farfield` を追加し、Y/X と k/omega の読込・検査を「入口または farfield」に広げる (`inletCornerWall` の再割当は入口専用のまま)。
- 陰解法: 境界半割面は現行の ghostless 対角 A⁺ のまま (初版)。安定性は V2 の陽解法/block-DPLUR 両方で確かめる。
- 既存境界の経路は変えない (V0)。
- runner: `evaluate.side_far_kind` / `top_out_kind` に `farfield` を追加し、`top_out_kind: outflow` が slip に落ちる不具合も直す。生成 YAML の実効 BC を照合する。

## 5. 実装ステップ

1. `boundaryCond.{hpp,cpp}`: `valueTypesOfBC["farfield"]`・読込・ディスパッチ・起動時の拒否。
2. `cuda_forge/boundaryCond_d.cu`: 境界状態の構成 (§4.2、ホスト/デバイス共通の関数にして単体試験から呼ぶ)。
3. `cuda_forge/convection/`: `farfield_flux_d` (§4.3)、`convectiveFlux_d.cu` の境界ループで farfield 面をこちらへ振る。
4. `cuda_forge/scalarTransport_d.cu`・`passiveKernels_d.cuh`・`speciesTransport_d.cu`・`ransBoundary_d.cu`: farfield 面の面値読み、k/ω の分岐。
5. `design/forge_design/evaluate/runner_sern3d.py`・`metrics/sern_momentum.py` (`OPEN_KINDS`)。
6. docs: `methods/boundary.md` (計画中 → 実装済み)、`procedures/recommended-settings.md` の境界の節。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~codex plan 段レビュー~~ **済 (2026-09-27 NO-GO C1/M6/m1、全件採用で §4/§6 を改訂)** | C1 → §4.3 スカラー面値経路、M2 → §4.3 `farfield_flux_d`、M3 → §4.2 決定表と単体試験、M4 → §4.2 γ 契約と V2b、M5 → §4.3 k 整合、M6 → §6 V1/V2 追加、M7 → §6 V3 再設計、m8 → §6 V0 | F |
| 1b | ~~codex plan 段 再レビュー (plan-2)~~ **済 NO-GO C0/M5/m2 (2026-09-27)** — 採否は §6.1。**検証範囲をユーザ判断待ち** | 改訂版 §4/§6 | F |
| 2 | 実装 (§5 の 1–4) | ビルド (AWS)、単体試験、V0 | O |
| 3 | 検証 V1–V2 (AWS、軽量) | §6 の合格条件 | O |
| 4 | SERN V3 (AWS) | §6 の合格条件。R4d の結論へ反映 | O |
| 5 | codex result 段 → accepted | | F |

## 6. 検証

すべて AWS で回す (ローカルは使わない)。合否は事前に固定する。判定ツールは `check_convergence.py` (判定区間を明記) と
`check_quasisteady.py`、原本 (`CONVERGENCE_VERDICT.txt`・判定出力) を run に残す。未収束 (NOT CONVERGED) の run は「感度診断」用途に限る。

- **V0 既存境界の不変**: farfield を含まない既存構成 (run_0971 設定) で、同一初期状態からの**初回面流束 `massflux` と状態ダンプ**
  (`FORGE_DUMP_MASSFLUX`) が変更前バイナリとビット一致。更新後の保存量は旧バイナリの反復実行 (3 回) の再現性幅以内。
- **V0u 単体試験** (§4.2 の決定表): 反例・音速通過・流向反転・$c_b\le0$ で期待どおりの分岐と有限値。
- **V1 自由流保持**: 一様流の 3D hex 直方体 (z 3 層以上)、全境界 farfield (または入口 1 面 + 他 farfield)、陰解法、2000 step。
  (a) CPG M 0.5、(b) TP 2 種 lump (SERN 外気) M 6、(c) (b) を面に対し 30° 傾ける。
  合格: 全節点で $|P/P_\infty-1|,|\rho/\rho_\infty-1|,|\mathbf u-\mathbf u_\infty|/|\mathbf u_\infty|\le10^{-5}$、$|Y-Y_\infty|\le10^{-6}$。
- **V2a 亜音速の微小音響パルス (反射率)**: 3D 薄板チャネル (z 3 層以上)、一様流 M 0.3、中央にガウス型の微小圧力パルス (振幅 1e-3 P∞)、
  一方の端を farfield、比較に 3 倍長の領域 (反射が評価時間内に戻らない)。評価点の圧力履歴から入射と反射を時間窓で分離し、
  **反射振幅/入射振幅 ≤ 0.05** (slip は ≈1、outflow も同時に測って記録)。SLAU と ROE、陽解法と block-DPLUR の 4 組合せで。
- **V2b 内外で状態が違う流入/流出**: 同じチャネルで、領域内部を $Y$・T・k・ω が自由流と異なる状態 (Y 0.13、T 600 K、k・ω を 10 倍) にして
  (i) 流出、(ii) 流入 (自由流が境界から入る向き) を回す。合格: 境界面流束と体積積分の時間変化の収支 (質量・全エネルギー・化学種・k) が
  相対 1e-4 以内で閉じる、(ii) で内部が自由流の値へ置き換わる (定常後 $|Y-Y_\infty|\le10^{-4}$)、NaN・床到達なし。
- **V2c 超音速の斜め衝撃波**: 3D 薄板、M 2.5、半角 10° ウェッジ、衝撃が上側境界に当たる配置。上側境界 (A) slip / (B) farfield / (C) 上方に 2 倍広げた
  slip (共通領域の格子は同一、`z_append` 型の追加)。評価線 = ウェッジ下流の壁面、定常後。合格: $\max|p_B-p_C|\le0.02\,\Delta p_{shock}$ かつ A は 0.1 Δp 以上
  (試験が反射を検出できること)。
- **V3 SERN 側方 (R4d)**:
  - **V3a 同一格子の対照**: g3・遠方面 2.50 H、同一新バイナリで `side_far` だけ slip ↔ farfield (top_out 等は同一、生成 YAML で照合)。run_0986 最終場から各 20000 step。
    R4d と同じ ε・D・窓条件。これは「BC の差」の測定で、合否でなく記録。
  - **V3b farfield の幅系列**: farfield のまま遠方面 2.50 / 3.42 / 4.35 H (`z_append`)。**合格 = farfield で隣り合う幅の D ≤ ε が 2.50→3.42 から成立する**
    (= farfield なら 2.50 H で領域独立)。slip の幅系列 (run_0986/0988/0989/0990) との比較は参考値 (両者とも無限遠との一致ではない)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan.md) | NO-GO, C1/M6/m1 | **全件採用** (C1・M2 はコードで再確認: `scalarTransport_d.cu:166` の `ext_is_self`、`convectiveFlux_boundary_d.inc.cuh` の流出時内部風上)。§4 を「面流束を自前で組む専用境界」に改訂、§6 に V0u/V2a–c/V3a–b。§5.1 #1 |
| plan | 2026-09-27 | [2026-09-27-boundary-node-farfield-characteristic-plan-2.md](../../notes/reviews/2026-09-27-boundary-node-farfield-characteristic-plan-2.md) | NO-GO, C0/M5/m2 | 構造 (専用流束 + 同時刻スカラー面値) は妥当。M1 (超音速切替で流束が不連続: 入力差 2e-6 で P_b 2.09 倍) → 採用予定: SU2 と同じく内部状態と構成状態を既存の近似 Riemann 流束 (SLAU) に渡して面流束を作る。M2–M5・m6–m7 (TP の独立参照解、SST ソース込みの収支、dual-time の音響試験、面流束ダンプでの収支) → 採用予定。**範囲をユーザ判断待ち (2026-09-27)** |

## 7. 影響範囲

- `solver_density_cuda/boundaryCond.{hpp,cpp}`、`cuda_forge/boundaryCond_d.cu`、`cuda_forge/speciesTransport_d.cu`、`cuda_forge/ransBoundary_d.cu`
- `design/forge_design/evaluate/runner_sern3d.py`、`design/forge_design/metrics/sern_momentum.py` (`OPEN_KINDS`)
- 既存ケース: 変更なし (新種別を書かなければビット不変、V0)
- docs: `methods/boundary.md`、`methods/index.md` (見出しのみ)、`procedures/recommended-settings.md`

## 8. 完了条件

- [ ] 関連 `methods/` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6 V0–V3)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/`
- [ ] [`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-27` — codex plan 段 NO-GO (C1/M6/m1) を全件採用し §4/§6 を改訂。
- `2026-09-27` — 初稿 (ユーザ「遠方境界入れたらすっきりかもね。やってみますか」)。`methods/boundary.md` に理論節を追加し、`outflow` の説明 (実装は全量コピー) を訂正。
