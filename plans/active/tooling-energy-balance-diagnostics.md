# エネルギー収支診断 (定常 node の同一状態における離散収支)

## メタ

- **area**: `tooling / architecture (output)`
- **status**: `draft`
- **related_docs**:
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md) — 残差組立てと時間積分の現在仕様
  - [`procedures/solver-settings.md`](../../procedures/solver-settings.md) — `output.level` / `extraFields`
- **related_plans**:
  - [`case-hypersonic-gap-heating-validation.md`](case-hypersonic-gap-heating-validation.md) — **発注元** (§4.8)。すきま壁の熱量を離散スキームの言葉で検算したい
  - [`../accepted/output-level-and-h0.md`](../accepted/output-level-and-h0.md)、[`../accepted/output-node-wall-surface-viz.md`](../accepted/output-node-wall-surface-viz.md) — **再利用する基盤** (重複ではない)
- **created**: `2026-09-19`
- **owner**: `sano`

## 1. 目的

**壁面の局所伝導寄与は既に取得できる** ([`viscousFlux_d.cu`:560](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) が
残差へ入れた寄与を `qwall = heatflux/sss` として保存し、[`output.cpp`:355](../../solver_density_cuda/output/output.cpp) が境界変数を出す)。
足りないのは、**(a) 壁 Dirichlet 拘束による熱授受 (拘束反力) と、(b) 領域境界を横切る全数値流束を、
同一状態で対応づけること**である。この 2 つが無いと「壁熱量が離散的に閉じているか」に答えられない
(場から積分した物理収支は SLAU の数値散逸・再構成・拘束操作を含まない)。

> 注意: フィールドが存在することと、その run で診断できることは別。旧 `case/24` の `run_isoT_condN_node` の
> 最終壁出力は上下とも `qwall` 全 17 点がゼロだった。

## 2. スコープ

- **やる (初版)**: **定常 node** の、**同一状態における離散収支**。
  収支式と拘束反力の定義、採取位相の固定、`res_roe` への全加算経路の採取、CV/面の HDF5 スキーマ、
  明示 opt-in、解除試験。
- **やらない (初版で受理しない = 診断要求時に明示的に拒否する)**: 非定常/dual-time、RK 各ステージ、
  `cell` 離散化、周期・軸対称 (§4.6 の単位・合算規約が閉じるまで)、任意幾何断面での CV 切り出し、
  収支に基づく自動判定ツール、運動量・化学種方程式そのものの診断
  (**ただしエネルギー式へ入る寄与は省略しない**、§4.3)。

## 3. 前提

- 出力登録は `variables.hpp` の `output_cellValNames` に無いと `extraFields` でも出ない (警告のみ)。
- `output.level>=2` は登録量を全出力し ([`output.cpp`:39](../../solver_density_cuda/output/output.cpp))、
  登録された CV 配列は一律確保される ([`variables.cpp`:282](../../solver_density_cuda/variables.cpp))
  → **静的登録は「要求していない run にもコストを課す」**。
- 面配列は `var.p` / `nPlanes` 管理で、CV 配列 (`var.c` / `nCells`) とは別系統
  ([`variables.cpp`:303](../../solver_density_cuda/variables.cpp))。

## 4. 設計方針

### 4.1 収支式と拘束反力 (codex M1)

面流束 $F$ を**外向き正**、体積ソース $S$ と拘束反力 $C$ を**流体への供給を正**とする。

$$R_i^{raw} = -\sum_f F_{if} + S_i,\qquad D_t(V_iE_i) = -\sum_f F_{if} + S_i + C_i$$

- [`nodeWallDirichlet_d.cu`:80](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu) は
  $roe=\rho(e(T_w,Y)+e_k)$ を上書きし、:90 で残差をゼロ化する。一方
  [`timeIntegration_d.cu`:1023](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) は
  **エネルギー補正をゼロにしても密度の更新は残す** → 壁温一定でも密度・組成変化で壁 CV の蓄積は動く。
- 定常の拘束行では $C_i = -R_i^{raw}$。**拘束前残差はすでに壁面の伝導流束を含む**ので、
  **実効壁熱量は「物理境界流束 + 拘束反力」から作る** (拘束前残差をそのまま壁熱量と呼ばない)。
- **壁別の実効熱量は集計式と帰属規則まで決める** (codex 2 巡目 M1)。角ノードは複数の `bcond.iCells` に
  重複して現れる ([`gmshReader.hpp`:2292](../../solver_density_cuda/mesh/gmshReader.hpp)) 一方、
  ゼロ化されるエネルギー残差は**その CV に 1 つだけ**
  ([`nodeWallDirichlet_d.cu`:170](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu)) なので、
  壁ごとに `iCells` 上で $C_i$ を積分すると**共有ノードの反力を二重計上する**。
  - 重複しない拘束 CV 集合 $W$ と、対応する物理壁面集合 $B_W$ に対して

    $$Q_{\mathrm{eff}}(W) = -\sum_{f\in B_W} F_f + \sum_{i\in W} C_i$$

    と定義する (**外向き $F$ と供給方向 $C$ をそのまま足さない**)。
  - **$C_i$ は CV ごとに一度だけ保存**し、複数壁が共有する反力は初版では**接合部の別勘定**にする
    (壁別に無理に配賦しない)。壁∩出口のような非壁面流束も保持する。
  - **第一内部列への供給は別定義**: 壁 CV と内部 CV の間の**全数値流束**。壁部分領域には接線輸送・ソースも
    あるので、$Q_{\mathrm{eff}}$ との無条件な同一視を禁止する。
  - 検証: **共有角を含む試験で、壁別勘定 + 接合部勘定の合計が領域全体に戻ること**。
- 状態上書き (温度ピン・no-slip) による $\Delta(VE)$ も記録する (過渡を扱うときに必要)。

### 4.2 採取位相 (codex M2)

- **残差組立ての最中に、同一状態のエネルギー・面流束・体積ソース・拘束前残差を 1 組で保存**する。
  [`main.cpp`:1583](../../solver_density_cuda/main.cpp) で残差を組み、:1610–1615 で DPLUR 補正、
  :1744 で更新後の場を出力するので、**素朴に退避配列を足すと「更新前の流束」と「更新後の `roe`」が
  同じファイルに入る**。
- 属性に `step` / 評価位相 / (将来) RK stage・dual-time subiteration / **空間残差か BDF 込みか**を書く。
- **定常 DPLUR の反復差分を物理的なエネルギー変化率に換算しない**。
- 出力時の再評価 (`assembleResidual` の再実行) は**しない** (状態ピンを含むため診断が計算を変える)。

### 4.3 `res_roe` への全加算経路 (codex M3)

**面流束と体積ソースを別々に**、かつ**各カーネルが実際に加算した値**を採取する (後処理で再計算しない)。

| 経路 | 実装 | 初版の扱い |
| --- | --- | --- |
| 対流 (SLAU) | [`convectiveFlux_slau_d.inc.cuh`:528](../../solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh) | **受理** (面流束) |
| 粘性・伝導・粘性仕事 | [`viscousFlux_d.cu`:320](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) | **受理** (面流束) |
| 化学種拡散のエンタルピー輸送 $\sum h_s J_s$ | [`speciesTransport_d.cu`:271](../../solver_density_cuda/cuda_forge/speciesTransport_d.cu) | **受理** (面流束)。非反応 TP が中核 |
| $k$ 拡散 (`sstEnergyIncludesK`) | [`viscousFlux_d.cu`:285](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) | 受理するなら**スキーマに対応を明記** (下記) |
| SST のエネルギーソース | [`ransSource_d.cu`:230](../../solver_density_cuda/cuda_forge/ransSource_d.cu) | **受理** (体積ソース、条件付き) |
| 軸対称の幾何ソース | [`axisymmetricSource_d.cu`:232](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu) | **拒否** (§4.6 の単位規約が未確定) |
| 体積力の仕事 | [`bodyForce_d.cu`:43](../../solver_density_cuda/cuda_forge/bodyForce_d.cu) | **受理** (体積ソース) |
| **化学反応熱** `res_roe += V·Qdot` | [`chemistry_d.cu`:108](../../solver_density_cuda/cuda_forge/chemistry_d.cu) (残差組立て中に [`main.cpp`:1402](../../solver_density_cuda/main.cpp) から) | **拒否** (初版は非反応に限る) |
| **陰的連成補正** | [`speciesTransport_d.cu`:502](../../solver_density_cuda/cuda_forge/speciesTransport_d.cu) ([`main.cpp`:1610](../../solver_density_cuda/main.cpp) から) | **採取対象外**。空間残差の物理ソースとして採ると誤る |
| 凝縮・壁モデルなど未検証の組合せ | — | **拒否** (設定条件で起動時に判定) |

- **採取区間を固定する**: 「**境界状態・物性・勾配の評価後**から、**壁残差射影の直前**まで」。
  これにより陰的連成補正 (`main.cpp`:1610 以降) は自然に外れる。
- **収支対象**: `sstEnergyIncludesK` を受理する場合、[`ransTransport_d.cu`:163](../../solver_density_cuda/cuda_forge/ransTransport_d.cu) の
  分割保持に従い「**残差は $roe+roK$ の式、保存配列 `roe` は平均流エネルギー**」という対応をスキーマに明記する。
- 初版の中核は**発注元に必要な「非反応 TP・層流 / 低 Re SST」**。拒否対象は**起動時に明示的に拒否**する
  (黙って一部だけ閉じない)。

### 4.4 データモデル (codex M4)

- **CV 診断と面診断で HDF5 スキーマを分ける**。面配列は `var.p`/`nPlanes` 管理なので
  `output_cellValNames` に足すだけでは出力できない。
- 面側の最小項目: `face_id` / owner・neighbor / 境界種別・`physID` / 向き / **使用した面積ベクトル** / 流束。
- **node の可視化用 primal 面と、残差を組む dual 面を混同しない**。
- §5 は「登録」より先に**データモデルと採取位置の確定**を置く。影響範囲に流束カーネル・`main.cpp`・
  変数確保/転送を含める。

### 4.5 opt-in の契約 (codex M5)

- **明示要求を初期化時に解決**し、有効時だけ登録・確保・採取・転送する。
- **`output.level: 2` だけでは有効化しない** (level 2 は登録量を全出力するため、静的登録だと巻き込む)。
- 確認: **ON/OFF で解が変わらないこと**、OFF で追加の確保・カーネル・同期が無いこと。

### 4.6 CV 選択・周期・単位 (codex M6)

- 対象領域は **「指定した solver CV の和集合」**と定義し、**所属が片側だけの dual 面**を領域境界にする。
  任意幾何断面で CV を切る機能は初版から外す。
- 周期: [`mesh.cpp`:699](../../solver_density_cuda/mesh/mesh.cpp) は合併前の部分体積を保存し、
  :710–717 で合併体積を全 member に複写する。さらに [`main.cpp`:1423](../../solver_density_cuda/main.cpp) は
  **壁拘束の後**に周期残差を合算・broadcast する → **出力 `volume` と各ノード値の単純積分は重複計上**。
  → root 単位の一回集計 / 部分 CV 単位の集計のどちらかに統一し、**対応表を出力**する。初版は受理しない。
- 単位: 平面 2D は単位スパン当たり `W/m`、3D は `W`、軸対称 method 0 は面積に半径を掛ける方式
  ([`variables.cpp`:505](../../solver_density_cuda/variables.cpp)) なので単位角度当たり量と全周換算を区別する。

## 5. 実装ステップ

1. **収支式と拘束反力の確定** (§4.1) — 何を「実効壁熱量」と呼ぶかを先に決める。
2. **採取位相の確定** (§4.2) — どの時点の状態を 1 組にするか。
3. **全加算経路の棚卸しと対応機能の確定** (§4.3) — 未対応なら拒否する条件も。
4. **データモデル (CV/面スキーマ) と opt-in 契約** (§4.4・§4.5) — ここまで決めてから登録・確保に触る。
5. **実装** (流束カーネルでの採取 → 転送 → 出力)。
6. **解除試験** (§6)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | 収支式と拘束反力 | §4.1。$R^{raw}$ / $C$ / 実効壁熱量、壁 CV 込み領域と第一内部列起点領域の別定義 |
| 2 | 採取位相 | §4.2。同一状態で 1 組、属性に位相。DPLUR 反復差分を変化率にしない |
| 3 | 全加算経路 | §4.3。面流束と体積ソースを別採取、`roe` か `roe+roK` か、未対応は起動時拒否 |
| 4 | CV/面スキーマと opt-in | §4.4・§4.5。面は `var.p` 系。level 2 では有効化しない |
| 5 | 解除試験 | §6。(1)–(6) + $\gamma_n$ 丸め限界の事前定義 + 受理経路ごとの非ゼロ試験 + 拒否の起動拒否試験 |
| 6 | 発注元との受け渡し | 発注元 ([case-hypersonic-gap-heating-validation](case-hypersonic-gap-heating-validation.md) §4.8) の離散収支評価は**本計画の解除試験が通ってから**開始 |
| 7 | **解除マイルストーン: dual-time** (発注: [boundary-conjugate-heat-transfer](../accepted/boundary-conjugate-heat-transfer.md) §4.3、CHT の V1 の前提) | 初版は dual-time を拒否している (§2) が、CHT 検証 V1 が dual-time を要する。追加で: **(a) $R^{raw}$ は壁残差射影より前の空間残差と固定** ([`main.cpp`:1814](../../solver_density_cuda/main.cpp) は空間残差の後に BDF 項を足すので、後から `res_roe` を読んでも復元できない)、**(b) 同一評価状態の BDF 係数・履歴・体積から $D_t(VE)$ を作り $C=D_t(VE)-R^{raw}$**、**(c) 蓄積が非零の試験で BDF1/BDF2 の符号と単位を検証** |
| 8 | **解除マイルストーン: 周期** (発注: 同上、CHT の V5 = 翼列の前提) | 周期は **root 単位の 1 回集計**とし、部分 CV と壁反力の対応を出力する。軸対称・壁モデルは引き続き対象外 (CHT 側もそれに合わせて合格対象外と宣言済み) |
| 9 | codex レビュー | plan 2 巡完了 (NO-GO → GO-with-changes)。実装後に `--stage result`。**#7・#8 は CHT plan の 2/3 巡目レビューで要求された追加項目**なので、実装前に本 plan の §2 スコープも更新する |

## 6. 検証 (codex M7 で全面改訂)

**符号**: ソース・拘束が無ければ $D_t(V E) = -\sum_f F$ (**負の流束和**)。実装は対流が owner に負、
粘性が owner に正を加えるので、採取時に符号規約へ揃える。

- **(1) 組立て恒等式**: 同一状態で、独立に保存した残差と「面流束 + ソース」を比較。
  **内部面の相殺・壁列・第一内部列・領域和**をそれぞれ検査する。
- **(2) 壁の拘束収支**: 既知の純伝導問題で、**物理境界流束・拘束反力・実効壁熱量を別々に**検査する。
  **全領域の正味ゼロだけでは上下壁を両方誤ってゼロにしても通る**。
- **(3) 実用途**: SLAU の対流・再構成・粘性仕事を含む node ケースと、発注元の**低 Re SST / TP 経路**。
  標準の `case/48` と整合させる。
- **(4) 定量判定 (codex 2 巡目 M3 で具体化)**:
  - **組立て試験の合格式**: CV $i$ に実際に加えた項を $a_{ij}$ として

    $$\left|R_i^{raw}-\sum_j a_{ij}\right| \le B_i,\qquad
      B_i \sim \gamma_{n_i}\sum_j |a_{ij}|,\quad
      \gamma_n=\frac{nu}{1-nu},\ u=2^{-24}$$

    **限界は $\sum|F|$ だけでなく体積ソース加算も含める** (体積力仕事などが漏れる)。
    保存時の丸めと集計誤差も含め、**領域和は FP64 で評価**する。
  - **「組立て誤差の範囲内」と「壁熱量を必要精度で検算できる」を分ける**。
    誤差限界が対象熱量の誤差予算を超えるときは、後者を**判定不能**とする (合格にしない)。
  - 解析解との差は**別の離散化誤差**として判定する。定常解の主張には `check_convergence.py`、
    熱量系列には `check_quasisteady.py` の VERDICT を貼る。
- **(5) 受理経路ごとの非ゼロ試験 (事前登録)**: 経路が**実際に発動したこと**も検査する。
  低 Re SST を回すだけでは [`ransSource_d.cu`:230](../../solver_density_cuda/cuda_forge/ransSource_d.cu) の
  条件付きソースを検証できず、一様組成の TP では
  [`speciesTransport_d.cu`:263](../../solver_density_cuda/cuda_forge/speciesTransport_d.cu) の種拡散が 0 になりうる。
  → **組成勾配 / $k$ 勾配 / 非平衡の $P_k-D_k$ / 体積力仕事**がそれぞれ非ゼロになる試験を割り当て、
  **拒否対象には起動拒否試験**を用意する。
- **(6) 純伝導の独立検証**: 上下壁それぞれの**符号・期待熱量・離散化誤差許容**を固定する。
  $C=-R^{raw}$ から作った恒等式が閉じるだけでは**熱量の正しさの独立検証にならない**。
- **既存 run を解除試験に流用しない**: `case/24.laminar_channel_bl/run_isoT_condN_node/` を再判定すると
  **`NOT CONVERGED (stalled/plateau)`** (`rms_roe` 4.42e-4, `rms_roUy` 8.27e-7, ともに rising)。
  新規 run はメッシュ品質・IC・段階起動・run 索引 (case README) の手順を踏む。
- **回帰**: 診断 OFF で既存ケースの場・残差・速度が不変。**ON/OFF で解が変わらない**ことも確認。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-19` | [`notes/reviews/2026-09-19-tooling-energy-balance-diagnostics-plan.md`](../../notes/reviews/2026-09-19-tooling-energy-balance-diagnostics-plan.md) | **NO-GO**, C0/M7/m1 | **全件採用**。M1→§4.1 (収支式・拘束反力・実効壁熱量の定義、密度更新は残る)、M2→§4.2 (採取位相を残差組立て中に固定、DPLUR 差分を変化率にしない)、M3→§4.3 (種拡散 $\sum h_sJ_s$・$k$ 拡散・SST/軸対称/体積力まで棚卸し、未対応は拒否)、M4→§4.4 (面は `var.p` 系で別スキーマ、primal と dual を混同しない)、M5→§4.5 (初期化時解決の opt-in、level 2 で有効化しない、ON/OFF で解不変)、M6→§4.6 (CV 和集合・周期の重複計上・単位規約、周期/軸対称は初版で受理しない)、M7→§6 (解除試験を 4 本に分割、符号を訂正、float 丸め限界、既存 case/24 run は未収束なので流用しない)、m1→§1 (壁の局所伝導寄与は `qwall` で取得可能。足りないのは拘束反力と領域境界の全数値流束の対応づけ) |
| plan (2 巡目) | `2026-09-19` | [`notes/reviews/2026-09-19-tooling-energy-balance-diagnostics-plan-2.md`](../../notes/reviews/2026-09-19-tooling-energy-balance-diagnostics-plan-2.md) | **GO-with-changes**, C0/M3/m0 (旧 M2/M5/M6/m1 は解消) | **全件採用**。M1→§4.1 ($Q_{\mathrm{eff}}(W)$ の集計式、**角ノードは複数 bcond に重複するが残差ゼロ化は CV に 1 つ**→二重計上の禁止、$C_i$ は CV ごとに一度、共有反力は接合部の別勘定、第一内部列への供給は別定義)、M2→§4.3 (**化学反応熱 `V·Qdot` の欠落を補い拒否**、陰的連成補正は採取対象外、採取区間を「境界状態・物性・勾配の評価後〜壁残差射影の直前」に固定、`sstEnergyIncludesK` は残差 $roe+roK$ / 配列 `roe` の対応をスキーマ化、経路ごとに有効条件・採取対象・拒否条件を表に)、M3→§6 (組立て合格式 $\gamma_n$、体積ソース加算も限界に含める、領域和は FP64、「組立て誤差内」と「熱量を必要精度で検算できる」の分離、受理経路ごとの非ゼロ試験と起動拒否試験、純伝導の壁別独立検証) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/` の流束・ソースカーネル (採取点)、`main.cpp` (採取位相)、
  `variables.{hpp,cpp}` / `output/output.cpp` (確保・転送・スキーマ)。
- 診断 OFF では追加の確保・カーネル・同期を持たない (§4.5 の契約)。

## 8. 完了条件

- [ ] §6 の (1)–(4) が通る (丸め誤差限界を事前定義したうえで)
- [ ] 診断 OFF の回帰が不変、ON/OFF で解が変わらない
- [ ] codex レビュー (`plan` 2 巡目 / `result`) を §6.1 に記録
- [ ] 発注元へ「離散収支評価を開始してよい」と伝える条件を満たす
- [ ] `status: done` にし `plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-09-19` — 初稿 (スタブ)。発注元のレビューで「離散収支の取得経路が無い」と指摘されたのを受けて起票。
- `2026-09-19` — codex plan レビュー 1 巡目 (**NO-GO**, C0/M7/m1) を全件採用して全面改訂。
  初版を**定常 node の同一状態における離散収支**に限定し、収支式と拘束反力・採取位相・全加算経路・
  CV/面スキーマと opt-in・CV 選択と周期/単位・解除試験 4 本立てを定義。周期/軸対称/非定常/cell は初版で受理しない。
- `2026-09-19` — codex plan レビュー 2 巡目 (**GO-with-changes**, C0/M3/m0) を全件採用。壁別実効熱量の集計式と
  共有角の別勘定、経路表への化学反応熱の追加と陰的連成補正の除外・採取区間の固定、解除試験の定量化
  (丸め限界 $\gamma_n$・経路別非ゼロ試験・起動拒否試験・純伝導の壁別独立検証)。**仕様確定 → 実装の順で進める**。
