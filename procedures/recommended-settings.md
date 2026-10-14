# 推奨解析設定 (現行レシピの正本)

**このファイルが「今どう設定するべきか」の正本** (2026-09-08 起票)。各キーの意味は
[`solver-settings.md`](solver-settings.md)、起動手順は [`divergence-and-startup.md`](divergence-and-startup.md)、
メッシュ/実行手順は [`calculation-workflow.md`](calculation-workflow.md) を参照し、ここには**選ぶべき値と理由・根拠 run**だけを書く。

- 各節の見出しに **現行 (YYYY-MM-DD)** を付ける。設定を変えたら日付を更新し、旧設定は末尾の
  「§9 旧設定 (superseded)」に移す。**日付の無い記述・§9 にある記述を新規 config に使わない**。
- 記憶や過去の run の config をコピーする前に、必ずこのファイルと突き合わせる (古い run の config には
  廃止キー・旧既定が残っている)。
- 新しい config を書くときは `.claude/skills/forge-config` (Claude) の手順どおり、解析種別 → 該当節 → 段階起動 → 検証の順で組む。
- **値の正本は機械可読のファイル** (2026-10-14、plan [tooling-generator-defaults-and-launch-gate](../plans/active/tooling-generator-defaults-and-launch-gate.md) §4.0・§4.1):
  機種と段ごとの推奨設定 (プロファイル) の値の正本は
  [`solver_density_cuda/tools/recipes/recommended_settings.yaml`](../solver_density_cuda/tools/recipes/recommended_settings.yaml)。
  **§1 は機種に依らない共通部分と、機種の決まっていない手書きの run の既定** (YAML の `common_node_ns_sst`・`common_node_euler`)。
  設計チェーンの機種ごとの値は [付録 A](#付録-a-機種と段ごとの推奨設定-生成区間) (YAML から生成した表) を見る。本文の値と YAML が
  食い違ったら YAML が正しい (本文を直す)。照合と宣言は §0.1。

## 0. 解析種別の選び方

| 解析対象 | 節 | 主な根拠 case |
|---|---|---|
| 平面/3D ノズル・ダクトの定常 NS (層流・SST) | §1 + §2 (+§3 TP) | case/16, case/26, case/40 |
| 非粘性 (Euler) 設計評価・MOC 比較 | §1 + §4 | case/23, 44, 45, 46 (design chain) |
| 軸対称 | §1 + §5 | case/23, 29, 41, 42 |
| 多成分 semi-perfect (TP)・凝縮・化学 | §3 | case/16 (TP+凝縮), 28, 44, 47 |
| 非定常 (LES/DES, 音響, 周期箱) | §6 | case/09, 18, 38, 39 |
| 変換メッシュ (Fluent/非直交) | §7 | case/24〜27 |
| 設計チェーンの機種 (ベル・風洞 axis-Mach・wall-driven・SERN など) | 付録 A (機種のプロファイル) | case/40, 41, 45, 46, 68 |

### 0.1 機種と段ごとの推奨設定 (プロファイル) — 現行 (2026-10-15)

- **プロファイル** = 機種と段 (`soft`・`mid`・`main`、SERN の層流暖機 `warm_lam` など) ごとの推奨設定。各キーは `fixed` (違えば FAIL)、
  `explicit` (明示すべきキー。照合は `fixed` と同じ。下の項)、`advisory` (推奨の範囲。外れたら WARN)、`forbid` (書けば FAIL) のどれか。
  どれにも無いキーは、問題 YAML の物理条件・入出力のキー
  (気体の物性・組成・本段の長さ・出力の間隔・`pRef` など、YAML の `problem_keys`) でなければ FAIL。
- **省略は FAIL** (2026-10-15、plan [tooling-generator-defaults-and-launch-gate](../plans/active/tooling-generator-defaults-and-launch-gate.md)
  §4.3 の 4 (a)・§5.1 #16): 宣言したプロファイル・段の `fixed`・`explicit` のキーを config に書いていなければ、出所 (生成器・手書き) によらず
  そのキーを名指しで FAIL にする。省略したキーの値は実行するバイナリの既定で決まる (`run_case.sh` は `FORGE_BIN` で別ビルドも回す) ので、
  **今のソースの既定がプロファイルの値と同じでも通さない** (2026-10-14 の版 1 は通していた)。省略そのものを試すときは
  `--deviate key=~ --reason "1 行の理由"` (「書かない」の逸脱) で宣言する。`RECIPE.json` の逸脱には実効値「未確定」が付く。
- **明示すべきキー** (`explicit`、plan §5.1 #16): 値に依存するのに設計チェーンの生成器が書いていなかったキー (`scalarGradient`・
  `implicitRelax`・`slauWallNormalChi`・`limiterScaled` など) を、機種・段ごとに旧い実効値 (書き写しに使ったソルバの既定、出所つき) で
  付録 A の表に「(明示すべき)」として載せてある。基準と出所は付録 A の「明示すべきキーの基準」(YAML の `explicit_keys`)。
  **今の生成器の config はこれらのキーの省略で FAIL する** (生成器に書かせるのは plan §5.1 #11。そのときこの一覧を
  「足してよいキーと旧い実効値」の例外表にする)。手書きの run でこの機種のプロファイルを宣言するなら、表の値を書く。
- **段取り 1 (2026-10-14)**: 設計チェーンの機種のプロファイルは**今の生成器が実際に書く設定の書き写し**
  (`recipes/regen_profiles.py` が代表の問題で全段の起動列を作って取る。手で写さない)。§1 と違うキーは付録 A の「§1 との差」に
  受入の状態「未」で載せてある。機種の生産レシピに揃えるのは段取り 2、§1 に寄せるのは受入の計算を事前登録してから (plan §6.4)。
- **手書きの run の照合と宣言**: `python3 solver_density_cuda/tools/recipe_settings.py check <run_dir> --profile P --stage S`
  (書かずに照合)、`... declare <run_dir> --profile P --stage S` (通れば `RECIPE.json` を書く)。プロファイルから外す値は
  `--deviate key=value --reason "1 行の理由"` (根拠 `--basis` は任意。キーを書かないなら `key=~`)、推奨に従わない試し計算は
  `--experiment --reason "..."` (外れたキーと理由だけを記録)。禁止キー・起動時に拒否されるキー・未知のキー・誤った節のキー・
  YAML の重複キーは、逸脱でも実験でも FAIL。一覧は `recipe_settings.py list`、中身は `recipe_settings.py show P`。
- 投入の関門 (`run_case.sh` の起動前の照合) はまだ無い (plan §5.1 #10)。宣言は今は記録で、強制はしない。

## 1. 共通の基本設定 — 現行 (2026-09-08)

```yaml
mesh: {discretization: "node", nodeWallDirichlet: 1, nodeInletCornerWall: 1, meshFileName: X.h5, valueFileName: X.h5}
solver: "SLAU"
space: {convMethod: 1, limiter: 2}                      # 本段。起動は convMethod 0 (§1.2)
                                                        # limiterScaled は既定 1 (無次元化 Venkatakrishnan, venkatK 0.05)
time:
  unsteady: 0
  last: {nStepOuter: N}
  deltaT: {control: 1, dt: 1.0e-5, cfl: C, cfl_pseudo: C, implicitRelax: 0.7, blockDPLUR: 1, dt_min: 1.0e-8, dt_max: 1.0, detectNaN: 1}
  timeIntegration: 11
  nStepInner: 4                                          # 2026-09-12: 5→4 (3D node SST の 12000 step 継続で 3 と 5 の残差経路が一致、余裕で 4)
output: {level: 1}                                      # 保存量 + 原始量 + h0 (2026-09-08〜)
```

- **`physProp` に `isCompressible` と `ro` は書かない**、`mesh.meshFormat` と `time.last.control` も書かない
  (§9.2 の段階移行キー。`time.deltaT.control` は別物で必要)。
- **離散化は node (median-dual) のみ** ([user-prefers-node-base]、2026-09-16 ユーザ決定)。cell の回帰対照も組まない。node のメッシュは
  **node 用に変換した h5** (`discretization: node` を書いた config で `convertGmshToForge`) を使い、2D は
  **平面メッシュ** (押し出し 2 ノード spanwise は 2 次 MUSCL の散逸が消えて発散)。
- 対流は SLAU。`convMethod: 1, limiter: 2` が本段の標準、`limiter: 0` は使わない、**`mesh.bndFirstOrder` は禁止**。
- **リミッタは `limiterScaled: 1` が既定** (2026-09-20 変更)。旧式 (`0`) は $\varepsilon^2=K^3|V|$ を次元のある $\Delta$ と
  比べており、**メッシュを拡大すると実質 OFF になり、変数ごとに効き方が桁違い**になる。
  修正版は $\Delta$ を変数ごとの基準で無次元化し、評価点を流束と揃える。`venkatK` の既定は **0.05**。
  根拠: 厳密 Riemann 解との L1 収束次数が 0.76/1.00 → **0.82/1.04**、近傍逸脱が数千 → **0**、
  座標 ×1024 でビット不変 ([`limiter-config-simplify.md`](../plans/active/limiter-config-simplify.md) §4.5)。
  **切り替えで答えは変わる** — ⑤ SERN の目的量で +0.54 %。旧値を再現したい run は `limiterScaled: 0` を明記すること。
  node 以外 / `convMethod` 対象外では**自動で 0 に落ちる** (警告を出す)。
- 定常は陰解法 `timeIntegration: 11` + `blockDPLUR: 1`。実効 CFL は `cfl_pseudo` (§solver-settings「CFL の定義」)。
  `cfl` は表示用なので同じ値を入れておく。**`nStepInner` は 4** (node NS の soft/mid 段は 10)。
  根拠 (2026-09-12, 3D node SST TP case/16 run_0410–0412): 本段の `nStepInner: 3` は 5 と残差経路が全列一致
  (到達 step も同じ) で 1 step が 11 % 速く、2 は発散 → 余裕を見て 4 をユーザ決定 (旧レシピ 5 は過剰反復)。
- **TP の温度反転は `physProp.thermoFloat: 1` が既定** (2026-09-12 ユーザ決定): float Newton + double 研磨 (収束まで最大 3 段, 通常 1 段) で
  基準バイナリと同じ (未収束プラトーの) 状態に留まり、12000 step 後の壁圧差 ≤1e-4・温度反転誤差 ≤4e-10·T、1 step −13 %。`thermoHrefTemp: 298.15` が前提で、datum 無し config では
  自動で 0 (double Newton) に落ちて警告が出る (plan performance-3d-node-sst-speedup)。
- **`cfl_pseudo` の目安**: node NS ノズル 4〜6 + `implicitRelax: 0.7` (上限は EOS 圧力床の洗浄で決まり、
  relax のみ有効: cfl 8 + relax 0.7 ≈ 3〜4 倍速 [implicit-cfl-ceiling-eos-floor])、Euler 設計評価 4、
  bump など易しい流れは 50 まで。TP 多成分の陰解法は 0.5〜2 から (§3)。
- `detectNaN: 1` は常時。EOS 床 `pMin/roMin/tMin` は既定 (1 Pa) のままで良いが、**無次元・低圧ケースは下げる**
  (Taylor-Green: 1e-6)。
- 出力は `output: {level: 1}` (既定)。勾配・リミッタ・診断が要る run だけ `level: 2` か `extraFields`。
  全温・全圧は `VALUE/h0` から `tools/total_quantities.py` で作る (AGENTS.md「出力と後処理の原則」)。

### 1.0a 壁隣接面の $\chi$ (`space.slauWallNormalChi`) — 現行 (2026-09-27)

**既定は auto = 実効 1** (node ∧ `nodeWallDirichlet: 1` ∧ SLAU/SLAU2 のとき。それ以外は 0。2026-09-26 ユーザ決定、
plan [convection-slau-wall-normal-chi-default](../plans/accepted/convection-slau-wall-normal-chi-default.md))。**書かない**のが推奨。
旧挙動 (2026-09-25 以前の結果) の再現だけ `slauWallNormalChi: 0` を明記する (演算はビット同一)。起動ログに実効値と解決理由が 1 行出る。

- **検証 (2026-09-27)**: 省略 (auto 1) と明示 0 を同じ起点から分岐し、SERN 2D 3 作動点の力係数 5 列 (差区間が R5n 帯内)、
  case/16 SST と凝縮の壁 p/p0 (L∞ ≤ 0.0012 %)・onset (差 0)、case/39 周期丘の下壁 $C_f$ (相対 L2 0.0002 %)・再付着点・継ぎ目比、
  case/40 軸対称ノズルの $\eta_{CF}$・$\dot m$ (±0.0001 %)・壁 p/p0 (0.38 %) がすべて許容内。3D 接続模型では flag 0 が壁 CV を排出して
  発散する構成を auto 1 が救う (起点から 8〜14 step で回復)。
- **周期・軸対称・凝縮も auto に含む** (上の追加域の検証で成立条件と許容を満たした)。回転周期は node では未対応 (起動エラー)。
- **既定変更をまたぐ run は途中から再開しない** (`scalarGradient` と同じ運用ルール、2026-09-27): 旧既定 0 で始めた段階起動 run を
  省略のまま新バイナリで再開すると、`stage_manifest.py` は同じ設定の起動を区別できず (起動記録との結び付けは plan
  [tooling-stage-manifest-launch-binding](../plans/active/tooling-stage-manifest-launch-binding.md) で後回し) 0 の段と 1 の段を 1 区間につなぐ。
  最初から回し直すか、`slauWallNormalChi: 0` を明記して続ける。**段の途中で切り替えない**。
- **0 に落とすかの判断材料** (既定化前は「1 にする条件」だったもの): 診断ツール `case/46.sern_design/cad/diag_wall_cv_budget.py` と
  `case/46.sern_design/cad/diag_applicability.py` (設定キーの検査 + 1 step `FORGE_DUMP_MASSFLUX` で対象 CV の全接続面を照合) で、
  `convMethod: 0` の起動区間に (i) 壁 CV の $\rho_w$ が 3 dump 以上単調減少し $\rho_w/\rho_i<0.1$、(ii) 全接続面の正味流出が 3 dump 以上持続、
  (iii) $\chi=0$ の壁隣接面を $\chi_n$ に置換すると補充が増える、の 3 つが成り立つ構成は **0 に落とさない** (0 では壁 CV が排出される)。
  検査に落ちる設定 (`slauContactFloor ≠ 0`、`convMethod: 1`、`lowMachPrecond ≠ 0` 等) は診断できないので、既定 (1) のまま使い、0 にする理由にはしない。
- **残る制限**:
  - **衝撃がランプに当たる点の壁圧への影響は未評価** (m6_on・m4_off とも登録条件で衝撃足を同定できず判定保留)。衝撃衝突のある構成で
    0/1 の差が問題になりうるときは、その構成で壁圧分布の比較を別途行う。
  - **軸対称ノズルの出口角 (壁∩出口) の壁圧は 0/1 で ±1 % 程度動く** (case/40 等温壁、他の壁点は ≤ 0.03 %)。出口角の壁圧を評価に使う場合は注意 (2026-09-27、plan §5.1 #12)。
  - 1 の run は SLAU の node カーネルのレジスタ上限に当たりうるので `FORGE_CUDA_BLOCKSIZE=128` (比較する run どうしで揃える)。
- **旧方針 (2026-09-25〜26、履歴)**: 「既定 0 のまま、1 にするのは 3D 側壁∩後端面接続構成だけ」「2D 生産で 1 にしたときの差は生産許容内だが、
  これは 2D を 0 のままにしてよい根拠であって 1 にしてよい根拠ではない」(plan [convection-slau-wall-normal-chi-usage-rule](../plans/accepted/convection-slau-wall-normal-chi-usage-rule.md))。
  2026-09-26 のユーザ決定で既定 1 に変わり、上の検証 (2026-09-27) で 2D・周期・軸対称・凝縮も含めた。

### 1.1 境界条件 — 現行 (2026-09-07)

- 入口: `inlet_Pressure` (Pt/Tt/組成/k/ω)。node NS では `mesh.nodeInletCornerWall: 1` を**変換時**に付ける
  (入口∩壁角の質量溜まりを根治)。入口境界層・全温分布・組成 (H2O) 分布は `ints: {inletProfile: 1}` +
  `inlet_profile_<physID>.csv` (`tools/gen_inlet_profile.py` で生成・照合, 手順は
  [`inlet-profile.md`](inlet-profile.md); 超音速入口 `inlet_uniformVelocity` は Tt を ρ/U/Ps に換算)。
- 出口: 亜音速は `outlet_statPress` に **逆流用 Pt/Tt を必ず添える**。超音速出口は `outflow`、または
  `outlet_statPress` の Ps を実出口圧に合わせる (node は壁列が常に亜音速なので Ps ≪ 実圧だと SST が出口列から
  unstart [node-supersonic-exit-outflow])。3D の出口角線で unstart が続く場合は出口バッファ (slip 延長, physID 別)
  + `mesh.wallDistExtraPhysIDs: [そのID]` + Ps 一致 (case/16 run_0226〜0228)。
  **⑤ SERN (case/46) は出口・`far_bottom` とも設計上つねに超音速なので、runner の `evaluate.outlet_kind` の既定を
  `outflow` にした (2026-09-23)**。この罠は**同じ case で 2 回踏んでいる**: run_0121 で出口 P 7.5 → 128 kPa に積み上がって発散
  (対策が `problem_r5_3d_sst_outflow.yaml` の個別 YAML に留まり既定へ反映されなかった)、接続模型 run_0430–0436 で
  出口の亜音速率 3.5 → 99.5 %・`far_bottom` の圧力 22 MPa・残差 +0.8 桁。**`outflow` に変えるだけで**同一起点・同一設定で
  亜音速率 4.4 %・逆流 0.1 %・全域 P>2e5 が 0 節点・残差 −1.3 桁になった (run_0435 対 run_0437)。
  **超音速出口に静圧を課さない**こと。
- 壁: NS は `wall` (断熱 no-slip)、等温は `wall_isothermal` (キーは `Ts`)。Euler は `slip`。
  **SST の壁は変換時に no-slip `wall` でないと wall_dist=0 になり ω が step 0 で発散**。
- 対称面/疑似 2D 側面は `slip`。

### 1.2 段階起動 — 現行 (2026-09-07, case/16 `run_user_profile.py` の手順)

1. **soft**: `convMethod: 0, limiter: 0`, `cfl_pseudo 0.5`, `nStepInner 10` (node NS), 2000 step。
2. **mid**: 同上で `cfl_pseudo 1.0`, 2000 step。
3. **本段**: `convMethod: 1, limiter: 2`, `cfl_pseudo 4〜6` + `implicitRelax 0.7`, 8000〜12000 step。

段間の引き継ぎは**同一メッシュなら index コピー** (`continue_run.py` / `--ic-index-from`)。**3D で `interp_field.py` を
同一メッシュ引き継ぎに使わない** ((x,y) 最近傍で z 列・座標一致ノードを混同 [interp-field-coincident-nodes-trap])。
メッシュを変えたときは `interp_field.py` (cross-mesh、原始変数補間)。一様 IC から超音速/SST を直接始めない。
初期 k/ω は非ゼロ (例 k=1, ω=1000; ω=0 は初手発散)。

### 1.3 検証と報告 (必須)

投入前 `check_mesh_quality.py` (AR ≤ 1000, skew ≤ 0.9)。まとめ前 `check_convergence.py` (全列) と
`check_quasisteady.py` (報告する量) の VERDICT を貼る。NaN チェック、run パスの明示、case README の run 一覧更新。

## 2. SST (RANS) — 現行 (2026-09-08)

```yaml
turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 2, katoLaunder: 1,
             wallTreatmentSST: 0, turbulentPrandtl: 0.9, kInit: 1.0, omegaInit: 1000.0}
```

> **訂正 (2026-09-17)**: 旧記述の `kInf` / `omegaInf` は**存在しないキー**で、書いても黙って無視され k=ω=0 の初期値になる
> (node SST が step 0 で壊れる既知の原因 [node-sst-init-and-ghostless-fixes])。実キーは `turbulence.kInit` / `omegaInit`。

- 既定で ON (キー省略時): `sstNodeWallKPin: 1`, `sstOmegaProdFromPk: 1`, `sstSigmaBlend: 1` (2026-09-08〜)。
  旧挙動を再現するときだけ 0 を明記。`sstIsotropicStress` / `sstEnergyKSource` は 0 のまま (分離型、非推奨)。
- **`sstEnergyIncludesK: 0` が既定で確定** (2026-09-08 ユーザ決定: roe に k を含めない E_m 形が安全。
  乱流の跳ねが T に伝播しない)。エネルギー整合の検証で必要なときだけ 1 (opt-in) にし、全温は h0 (k 込み属性) から。
- `dilatationCorrection: 2` はノズル BL を SU2 比 −16 % 薄くする (モデル形式差、ソルバ無罪)。SU2 と揃える A/B は 0。
- `katoLaunder: 1` はノズル喉のよどみ偽生産抑制用。平板など剪断層主体では 0 でも可 (A/B で無影響)。
- **壁処理: `wallTreatmentSST: 0` (低 Re 壁解像) を使う。`1` (SST 壁関数) は使わない** (2026-09-20 ユーザ方針)。
  **ソルバ既定も 0 に変更**し、`1` を指定すると起動時に既知欠損を列挙した警告が出る。
  設計チェーンの生成器で `evaluate.wall_treatment_sst` (既定 0) を読んで書くのは **SERN (`runner_sern`) だけ**。③ベル TP
  (`runner._solver_config_tp`) は `wallTreatmentSST: 0` を書き、`evaluate.wall_treatment_sst` が 0 以外なら止める。風洞の NS
  (`runner_wt._config_sst_node`、axis-Mach の NS も) は `wallTreatmentSST: 0` を直書き、③ベル CPG (`runner._solver_config`) は書かない
  (= ソルバ既定 0)。機種ごとの実際の値は付録 A (2026-10-14 訂正: 旧記述「生成器は `evaluate.wall_treatment_sst` を読む」は SERN だけの話だった)。

  **使わない理由** (積み上がった既知欠損):

  | 欠損 | 根拠 |
  | --- | --- |
  | Cf −6 % | 壁関数 P_k 規約の欠損 ([node-wallfunction-pk-convention-deficit]) |
  | 3D の角線ノードで u_τ=0 | 代表点が無い ([node-sst-wallfunction-utau-zero]) |
  | **壁モデル渦粘性に上限が無い** | `ν_t = ν(1/g−1)` が `g→0` で発散。低密度域で **4660 m²/s → k 2.06e9** となり 3D SST が死ぬ ([`tooling-nozzle-sern-3d.md`](../plans/active/tooling-nozzle-sern-3d.md) §4.25) |
  | case/40 の壁温 | y⁺1 低 Re と SU2 壁関数のみが根拠。**node 壁関数系列は撤回済み** |

  **帰結**: y⁺≈1 を満たすメッシュが要る。y⁺≈1 と AR ≤ 1000 が両立しないときは
  **壁関数に逃げず**、AR 緩和 (壁法線構造格子は ≤5000 可) か形状・解像度の見直しで解く。
  やむを得ず `1` を使う run は **README に理由を書く**こと。
  断熱壁の壁温出力は `sstThermalWallFunction: 3` (defect-flux) だが、これも `wallTreatmentSST: 1` 依存である。
- 3D node SST の角部加熱は k/ω 拡散の相対ゼロ割ガード (2026-09-08 修正) で解消済み。**3D SST の壁圧は実験より +9 %
  (側壁合流域の乱流 BL 過厚) が未解決** (case/16 run_0228)。定量比較には 3D 層流 (0.7 %) か 2D SST (+1.5 %) を使う。
- restart で `vis_turb` が再現されない (敏感な擬似衝撃波は位置が動く) [forge-sst-restart-nonfidelity]。

### 2.1 遷移モデル (γ–Re_θt) — 現行 (2026-09-22)

```yaml
turbulence: {model: "sst", scalarDiffusion: 1, dilatationCorrection: 0, katoLaunder: 0,
             wallTreatmentSST: 0, turbulentPrandtl: 0.9, transition: "lm2009"}
```

- **使う場面**: 層流域が長く、そこでの摩擦・熱伝達が目的量に効くとき (翼の負圧面前半、低 Re 平板)。既定は `none`。
  キーの意味と受付条件は [`solver-settings.md`](solver-settings.md)「遷移モデル」。
- **前提は局所 $y_1^+\le1$** (`check_wall_resolution.py` の VERDICT)。C3X は第一層 2 µm だと面積の 67 % が 1 を超えるので 1 µm メッシュを使う。
- **SST の収束場から継続する** (段階起動の最後の段)。$k$/$\omega$ を入口値に戻す必要は無い (T3A で同じ解)。
- **遷移位置が落ち着くまで長い**: T3A は cfl 5 で 50000 step、C3X は cfl 0.5 で 150000〜200000 step (粘性比 30 の正圧面は 60000 step で +18 %、200000 step で −22 % — **短い run は符号まで違う値を読む**)。`rms_roGamma` が下がっても遷移終了位置は動き続けるので、
  **報告量の時系列** (`cf_plate.py --series` / `h_series.py`) を `check_quasisteady.py --series-csv` に掛けて止める。
- **入口の $\omega$ (入口粘性比) の感度を必ず付ける**。入口から物体までの乱れの減衰は $\omega$ で決まり、実験報告は普通これを与えない。
  C3X では入口粘性比 1 / 10 / 30 / 100 で正圧面の $h$ 偏差が −39 / −33 / −22 / +25 % と動く。**実験に合う値を事後に選ばない**。
- `dilatationCorrection` との組み合わせは未検証。`katoLaunder: 1` は前縁よどみ点の $k$ 過大生成を抑え、C3X の正圧面の結果を大きく変える (粘性比 100 で +25 → −26 %)。どちらが正しいかは決まっていない。
- 検証の到達点: T3A 平板で同一メッシュの SU2 LM と遷移開始 −2.2 %・$C_f$ ±0.8 % ([`case/57`](../case/57.transition_flat_plate/README.md))。翼では遷移が実測より遅く急に起きる (C3X 負圧面 $s/S$ 0.32、実測は 0.2 から緩やかに)。

## 3. 多成分 semi-perfect (TP)・凝縮・化学 — 現行 (2026-09-07)

```yaml
physProp: {thermalMethod: 2, thermoHrefTemp: 298.15,
           species: [{name: "MIXDRY", lump: {"CO2": 0.1255, "O2": 0.0226, "N2": 0.8518}, basis: mole}, "H2O"], ...}
```

(設計チェーンの生成器が書く形: 2026-10-14 訂正。lump の係数はソルバが起動時に構成種から合成するので、合成済みの擬似種の
`species_db.yaml` も `speciesDBFile` も書かない。ソルバ内蔵で解決できない実種があるときだけ、その生エントリを
`species_db_external.yaml` に置いて `speciesDBFile` で渡す。lump 内の値は全桁のモル分率。上の数値は例で、生成器は
`evaluate.tp_species` (例 `split_h2o` = 乾き成分を MIXDRY にまとめて H2O を残す) から作る)

- **`thermoHrefTemp: 298.15` は必須** (絶対基準 h では χ_eos が桁違いになり TP×node 軸対称が発散 [isobutane-wt-semiperfect])。
- **LJ パラメータの出典 `physProp.ljSource` は書かない (既定 `[gri30, svehla1962]`)** — 現行 (2026-10-01, plan
  [thermophysics-solver-owned-species-db](../plans/active/thermophysics-solver-owned-species-db.md) §4.10 #14)。#14 前の run を同じ物性で
  再現・継続するときだけ `ljSource: [legacy_v1]` を明示する (旧既定 = §9)。既定変更の影響 (H2・OH・H・O・NO・CO の LJ が GRI 値に) の
  確認は plan #14-L2 で進行中。
- 輸送種の配置は `evaluate.tp_species` (`forge_design.gas.composition` の `SpeciesLayout`、例 `split_h2o` は乾き空気を lump MIXDRY に
  まとめ H2O を残す)。config には lump 記法 (`solver_species_config`・`physprop_species_flow`) で書く (上の例)。
- TP 陰解法の `cfl_pseudo` は 0.5〜2 から上げる (H2O 生成エンタルピーの増幅で上限が低い)。**`implicitRelax: 0.7` を付ければ 6〜8 まで可**
  (2026-09-16 case/44 va3 M4.19 node Euler 軸対称 TP 2 種 + 非平衡凝縮 `run_0181`–`0189`: cfl 6/8 + relax 0.7 は乾き一様場からの起動でも安定で場・残差床が cfl 2 と同じ;
  relax 無しの cfl 6 は軸列の EOS 床洗浄で発散 [implicit-cfl-ceiling-eos-floor と同じ機構]、relax 無しの cfl 4 は完走するが残差床が 2〜10 倍高い)。
  定常 precond × 多成分は `speciesPrecondDt: 1` (既定) — **ただしこのキーは化学ブランチ `feature/chemistry-finite-rate` 限定で、main / sern には未マージ** (2026-09-17 確認)。TP 亜音速 `outlet_statPress` の γ 混用は修正済。
- 凝縮: 平衡凝縮を選ぶなら `condensation: 1, condEquilibrium: 2` (EOS 拘束形、厳密 S=1) を推奨 (設定既定値は 0 = 非平衡)。蒸発は既定 ON。
- 凝縮 (2026-09-15, plan [condensation-source-limiter-steady](../plans/accepted/condensation-source-limiter-steady.md)): 非平衡の θ 律速は
  **`condLimiterMode: 1` (既定) で更新クランプ**になり、ソース残差の明示的な Δτ 依存 (θ×Δτ_loc) を除去した (旧 0 は残差に θ を掛け、大型ノズルで成長を 1/4 に絞っていた; case/44 で cfl 2 の解が旧 cfl 0.5 と一致)。cfl 間の固定点一致の検証状況は plan condensation-source-limiter-steady §9。
  凝縮 ON の定常 run は `output.level 2` の `condLim_<s>` が収束時に全域 ≈1、`condClampCorr_<s>` が 0 であることを確認する。
  RK 陽解法では自動で 0 に降格 (起動ログ `[condensation] condLimiterMode=`)。dual-time + S3 (`speciesFaceReconstruction 2`) では既定 `passiveFct 1` の保存的 FCT 補正が受動種の有界性を担う (`check_passive_budget.py` で収支 PASS を確認)。dual-time は既定 `passiveScalarScheme 1` ならモーメントに BDF 物理時間項が付き
  更新クランプ (mode 1) のまま有効 (2026-09-17, plan species-passive-scalar-unification); `passiveScalarScheme 0` のときだけ 0 に降格。上限は `condDgMaxStep` 5e-3 / `condDTmaxStep` 1 K。
  凝縮 run は h0 保存を確認する (面温度修正済み)。onset は実験より ~5 mm 下流 (case/16 2026-09-08 比較)。
- 受動スカラ (2026-09-17, plan [species-passive-scalar-unification](../plans/accepted/species-passive-scalar-unification.md)): トレーサ・凝縮モーメントは既定で化学種経路
  (`passiveScalarScheme 1`)。2 次面移流 (S3) を使うなら `speciesFaceReconstruction 2` + `speciesImplicitCoupling 1` の組で (coupling 0 + S3 は定常で発散)。
  **`implicitRelax 0.7` は定常 (擬似時間) の安定化として推奨**で、**dual-time の非定常計算には使わない** (§6 の dual-time 節): sub-iter の残差ノルムは下がるのに遅いモードが
  収束せず、同じ物理時刻の解が sub-iter 数に依存する (case/44 `run_0399`–`0404`: nSub 40→80 の ρ 差が緩和なしの 1.3e-5 に対し 6.3e-4; 2026-09-17,
  plan [species-passive-scalar-unification](../plans/accepted/species-passive-scalar-unification.md) §5.1 #26)。dual-time で安定化が要るときは `cfl_pseudo` を下げるか nSub を増やす。
  S3 は凝縮 onset を 0.2 r_t 程度下流に動かす (数値拡散減) ので、実験比較の基準を変えるときは明記する。
- 凝縮 (2026-09-10, plan [condensation-kantrowitz-gamma-twophase-sonic](../plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)):
  Kantrowitz 補正 (`condKantrowitz: 1`) の γ は凝縮種 (蒸気) の γ_v が既定 (`condKantrowitzGammaMode: 0`; 旧=1)。
  凝縮セルの音速は `condSonicModel` (未指定=自動): **TP carrier H2O・`condEquilibrium 0`・境界が
  `inlet_Pressure`/`outflow`/`wall`/`slip`/`periodic` のみ** のとき二相 frozen 音速 (1)、それ以外は旧の全蒸気音速 (0)。
  起動ログ `[condensation] condSonicModel=... (理由)` で解決値を確認する。`outlet_statPress`/`wall_isothermal` を含む構成と
  `condEquilibrium 1/2`・pure N2・CPG は未検証 (明示 1 は警告つきで可)。
- 凝縮 (2026-09-13, branch feature/condensation-air; plan [condensation-kantrowitz-carrier](../plans/accepted/condensation-kantrowitz-carrier.md) /
  [condensation-air](../plans/accepted/condensation-air.md)):
  - `condKantrowitz`: **省略時既定 0** (非等温補正なし; `solverConfig.cpp`) / 1=Kantrowitz (純蒸気形) / **2=Feder carrier 形** (キャリアの冷却込み, $\hat q=b-\tfrac12$) / 3=同+表面項 ($\hat q=b-\tfrac12-\ln S$)。
    Wysłouzil 参照設定 (case/16 run_0335 系) は 1 を明示する。2/3 は θ 167→3–4 で計算上の onset が mode 1 より 8.2 mm 上流に動く (**モデル間差**; 実験との位置差ではない) ため推奨は 1 のまま。診断出力は `condTheta_0`/`condLim_0` (`extraFields`, 種番号付き)。
  - `condFloat` (既定 1, 2026-09-13): 凝縮経路の float 実体 (凝縮 ON の速度; [methods/condensation.md](../methods/condensation.md) 実装 §9)。A/B や旧結果の再現には 0。
  - `condSigmaScale` (既定 1.0): 表面張力の定数倍感度。Wysłouzil mode 3 基準で σ×1.03 は onset +2.5 mm (収束 PASS+STEADY)、σ×0.97 は −2.2 mm だが**未収束の準定常参考値** (残差 plateau, h0 偏差 OSCILLATING; case/16 run_0354)。Tolman 補正の代替ではない。
  - 空気凝縮 (CPG carrier 形): `physProp: {cp: 1008.7, gamma: 1.4}` (二成分空気 R 288.19) + `condensation: {condModel: 0, nCondSpecies: 1, condVaporMassFraction: 0.7671}`。
    受付は SLAU × `thermalMethod 0` × 非平衡 × `condKantrowitz ≤1` × 境界 `inlet_uniformVelocity`/`outlet_statPress`/`slip` のみ (config 検査で拒否)。
  - N2 低温物性: `condN2LatentLowT: 1` / `condN2PsatLowT: 1` (既定, 70 K 未満の潜熱線形外挿と整合する飽和圧), `condN2LiquidCp: 2000` (J/kg/K; 1500/2000/2500 で onset 38.95/39.67/40.25 K = +500 で +0.58 K, −500 で −0.72 K の同符号)。
    0/0 は旧一式 (A/B 用)。ログの `[condensation] WARNING: CPG two-phase temperature inversion failed` が出た run は使わない (保存量を凍結して継続する退避)。
- 化学 (H₂): 定常陰解法は `speciesImplicitCoupling: 2` + `jacobianMode: 2` 必須、Cabra は dual-time 必須
  ([chemistry-finite-rate-direction], branch feature/chemistry-finite-rate)。

## 4. Euler 設計評価 (design chain) — 現行 (2026-09-04)

**値は機種ごとのプロファイル (付録 A) が正本** (2026-10-14 訂正: 旧記述「`runner.py` が生成する config が正本」は誤り。
③ベル CPG の `runner._solver_config` は `viscMethod 1` 固定で NS (層流・SST) しか作れない)。Euler の設計評価の config を作るのは
`runner_wt._config_euler_node` (axis-Mach Euler・wall-driven・帰還 `euler_loop`。① モード F の `runner_wt.prepare(euler=True)` は
cell の `_config_euler`) と、③ベル TP の `runner._solver_config_tp` (`evaluate.model: euler`)。今の生成器の値は
全域 2 次 (`convMethod 1, limiter 2`) + 陰解法 `timeIntegration 11, blockDPLUR 1`、本段 `cfl_pseudo 4` (axis-Mach は問題 YAML の
`evaluate.cfl_main`)、`nStepInner 5` (§1 の 4 と違う。③ベル TP は 4)、`lowMachPrecond 0`、壁 `slip`。§1 との差と受入の状態は付録 A。
δ* (排除厚) は CONTUR 初期壁 + 固定 Euler 基準・帯局所抽出 1 pass ([deltastar-production-recipe])、質量流量ゲート必須。

**軸対称ノズルの物理壁 (NS) — 現行 (2026-10-07)**: joint 壁 (`wall_repr: joint`) の物理壁の上流は **`geometry.pw_upstream: poly`**
(既定。新しい YAML・テンプレートには明記する): スロートより上流に δ_r を足さず、配管〜設計スロートを 5 次多項式 1 本にする。
`pw_ramp` は書かない (`poly` と併記すると例外)。寸法 (`spec.r_throat`) は出口径から (`solve_rt`) か物理スロート径から (`solve_rt_throat`)
決め、`spec.sizing` に残す。手順は [`nozzle-design-workflow.md`](nozzle-design-workflow.md) §1 ④・§2 ⑤、仕様は
`methods/design/overview.md`「上流の作り方」(plan [tooling-nozzle-upstream-poly-and-throat-sizing](../plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md))。

## 5. 軸対称 — 現行 (2026-08)

- `isAxisymmetric: 1`, `axisymMethod: 0` (r 重み方式) が既定。SU2 流 `axisymMethod: 1` は実装済みだが implicit 深収束未達。
- `space.pRef` に動作静圧 (hoop 源とゲージ整合、自由流誤差 1.6e7 倍改善)。`hoopAreaFromClosure` 利用可。
- node の軸ノードは通常 DOF + u_r=0 ピン (固定、キー不要。旧 `nodeAxisDirichlet` 等は起動時エラー)。
- 近軸の block-DPLUR 粘性対角は幾何修正済み (float で解ける)。

## 6. 非定常・LES/DES・周期箱 — 現行 (2026-09-17)

- 過渡 (音響) を含む閉じた系は **`unsteady: 1` 必須** (定常局所 dt は数 step で発散 [steady-localdt-acoustic-transient-instability])。
  陽解法は `timeIntegration: 4` (RK4) か 3 (RK3) を `unsteady: 1` で。dual-time は `timeIntegration 11 + dualTime 1`、
  物理 CFL ≲ 12 ([dualtime-subiter-divergence-fingerprint])。
- **dual-time の内部反復 — 現行 (2026-09-17)**: `cfl_pseudo: 12`〜`20`、`nSubIterDualTime: 10`〜`20`、**`implicitRelax` は書かない** (既定 1.0 = 緩和なし)。
  - 擬似時間の CFL に安定限界が見当たらない (物理時間項が対角に $V/(a\Delta t)$ を足すため)。case/44 で `cfl_pseudo` 2〜**160** を緩和なしで完走、NaN 0。
  - 収束解は `cfl_pseudo` に依存しない (nSub 80 同士で ρ 差 3.1e-4 = 同一 cfl の nSub 40 vs 80 差 2.3e-4 と同水準)。
    一方**必要な nSub は `cfl_pseudo` で決まる**: nSub 80 との ρ 差は `cfl_pseudo` 2 で nSub 10 → 1.2e-1 / 20 → 4.5e-3 / 40 → 2.3e-4 なのに対し、
    **`cfl_pseudo` 12 以上なら nSub 10 で既に 2.2e-4** (float のノイズ床)。壁時間は nSub 10 が 0.045 s/step、40 が 0.14、80 が 0.265 →
    **同じ解が 1.8〜3 倍速い**。根拠は case/44 `run_0459`–`0480` (生産 float, 湿り核生成場, convMethod 1 + S3 + FCT)。
  - **`implicitRelax` (定常の推奨 0.7) を dual-time に持ち込まない**: 安定性は物理時間項が担うので緩和の効果は無く、遅いモードの収束だけ遅らせる。
    しかも残差ノルムには出ない (2.8〜3.5 桁下がって見える) ので、同じ物理時刻の解が nSub 依存になる
    (case/44 `run_0399`–`0404`: nSub 40→80 の ρ 差が緩和なし 1.3e-5 に対し 0.7 で 6.3e-4)。
  - **sub-iteration が足りているかは残差の桁数でなく `nSub` を倍にして解が動かないかで見る** (float では残差床がノイズ床に当たるため)。
    安定化が要るときは緩和ではなく `cfl_pseudo` を下げるか nSub を増やす。
  - 確認は `python3 solver_density_cuda/tools/check_solver_config.py <run_dir>` (投入前) と、上の nSub 倍増比較 (精度が要る run)。
  - 限定: 検証は case/44 の 1 形状・物理 CFL 2・dt 8e-6 の範囲。物理刻みを大きく取ると対角が痩せるので擬似 CFL の余裕も減るはず。
- 解像 LES/DNS は `solver: "KEEP"` + `keepDissType/Coeff` (**トップレベル**キー; `space` 配下は無視される) で
  σ=0.05 (市松抑制) / 0.02 (解像 LES)、`WALE` は off + ES 散逸 ([keep-es-dissipation-status], [wale-inactive-fix])。
  一様流 U∞≠0 は `space.roRef/uRef` (KEEP, CPG) で機械精度保存。fdblend は音響モード成長のため opblend (増分フル c) を使う。
- 低マッハ (backstep 等) の市松は `lowMachPrecond: 2` (1 は不可)。
- 周期境界は node/cell とも実装・検証済 (DOF 同一視)。体積ソースは `volumePartial` (seam 2 倍バグ修正済)。

## 7. 変換メッシュ (Fluent / 非直交) — 現行 (2026-07)

`calculation-workflow.md`「変換メッシュを安定に回す実行レシピ」を厳守: 1 次で投入、初期場を入口流速に合わせる、
陰解法 `cfl_pseudo`≈1 `nStepInner`≈10、非直交は `pRef`、出口は `outlet_statPress` + 逆流 Pt/Tt。

## 8. メッシュ — 現行 (2026-09)

- AR ≤ 1000, skew ≤ 0.9 (`check_mesh_quality.py`)、y⁺≈1 が要るときは第一層と接線長のバランスで AR を守る。**壁法線の構造層に限り AR ≤ 5000 まで緩和可** (2026-09-12; `--ar-max 5000` / 問題 YAML `mesh.ar_max`, 台帳に明記)。冷却壁では y⁺ が ×5〜6 に上がるので `mesh.wall_first_frac_throat` (スロートだけ第一セルを詰める) と併用する。
- node 2D は平面メッシュ、3D は六面体押し出し可。3D の角線ノード (2 壁交線) は内部隣接ゼロだが受動で問題なし
  (2026-09-07 検証)。y₁ 0.6 µm × z₁ 2 µm の極端な異方角セルは旧バイナリで角部加熱を起こした (相対ガードで解消)。
- SST メッシュは壁 bcond を no-slip `wall` で変換 (wall_dist)。slip 延長壁は `wallDistExtraPhysIDs`。

## 9. 旧設定 (superseded) — 新規 config に使わない

### 9.1 削除したキー (2026-09-18, plan [config-key-pruning](../plans/accepted/config-key-pruning.md))

**書くと起動時エラーになる**。黙って無視されないので、古い config を再実行すると気づける。

| 削除したキー | 理由 | 移行先 |
| --- | --- | --- |
| `time.deltaT.lineDtWallRelief` | 壁半割面の λ 除外は発散した診断スイッチ | 既定経路 (何も書かない) |
| `time.deltaT.implicitRelaxSST` | 独立 3 件の A/B で効果なし | `implicitRelax` (SST も同じ値に従う) |
| `turbulence.C_DES_kw` / `C_DES_ke` / `wmlesNewtonTol` / `wmlesNewtonMaxIt` / `mesh.gradLSQDegenThresh` | 一度も使われず掃引の実測も無い内部定数 → コードに固定 | (固定値。変えたいときはコードを直す) |
| `turbulence.wmlesPrt` | Kader 原式への修正で壁法則が Pr_t を使わなくなり、**読む側が消えた**引数だった | (移行先なし。乱流 Prandtl 数は `turbulentPrandtl`、壁法則とは別物) |

### 9.2 第 2 陣で削除したキー (2026-09-18)

**書くと起動時エラーになる** (solver・変換器の両方)。2026-09-18 に S1 (任意化 + 警告) を経て S2 (拒否) へ移した。

| 削除したキー | 理由 | 移行先 |
| --- | --- | --- |
| `physProp.isCompressible` | solver はこの値をどこでも読まない (圧縮性は常に有効。非圧縮 SMAC 経路は現在のコードから到達できない) | 行ごと消す |
| `physProp.ro` | solver はこの値をどこでも読まない (密度は EOS で決まる) | 行ごと消す |
| `time.last.control` | solver はこの値をどこでも読まない (終了条件は `nStepOuter` のみ) | 行ごと消す。**`time.deltaT.control` は別物で、そちらは必須のまま** |
| `mesh.meshFormat` | 受け付ける形式が `hdf5` 1 つだけだった | 行ごと消す |
| `time.deltaT.detectNaNInterval` | トップレベル `detectNaNInterval` と同じものを読む 2 つ目の綴り | **トップレベルへ移す** (値は同じ意味) |

**互換の打ち切りと移行手順**: 2026-09-18 より前に作った run config は、そのままでは**再実行できない**。
移送は 1 コマンドで済む。

```bash
python3 solver_density_cuda/tools/migrate_solver_config.py [--dry-run] [PATH ...]
```

既定はこのツリーの `case/` 配下。**他の worktree は明示したパスでしか触らない** (各ブランチの作業物なので、
必要になったそのブランチで回すこと)。ツールは節を特定してから編集し、書き戻す前に
「元 − 消したキー」と厳密一致することを確かめる。**自動移送しない入力**は理由つきで報告する:
anchor / alias / merge key を含むもの (節どうしが同じ実体を共有しうる)、重複キーを含むもの
(solver の yaml-cpp は先勝ち、PyYAML は後勝ちで照合できない)、`meshFormat` が `hdf5` 以外のもの
(現行 solver では動かない旧入力で、キーを消すだけでは移送にならない)。

**削除していないもの (opt-in のまま残置)****削除していないもの (opt-in のまま残置)**: 過去に「不採用」と判定したスイッチでも、元 plan が「opt-in で残す」と
決めているものは残っている。`blockDPLURDiagCache` / `blockDPLURDqPack` / `mesh.primPack` (性能の再現用)、
`updateGuardAlpha`、`lowMachThornber`、`multispeciesRhoYCommonLimiter`、`passiveFctTolAbs`、
`turbulence.turbulentSchmidt` (`physProp.Sc_t` の別名)、`physProp.isAxisymmetric` / `axisymMethod` (deprecated 読み)。
**既定に採用しないことと、設定機能を廃止することは別の決定**であり、後者には残置決定を置き換える判断が要る。

存在しないキーを書いても solver は黙って無視するので、**投入前に `check_solver_config.py` を通す**。完全修飾パスで
照合し、**節の位置が違うキー** (トップレベルの `lowMachPrecond`、`physProp.speciesFaceReconstruction`、
`space.keepDissType` など。書いた場所が違うと黙って無視される) とどこにも無いキーを WARN、solver が起動時に拒否する
キーを FAIL で拾う。

| 旧設定 | 状態 | 代替 |
|---|---|---|
| `mesh.bndFirstOrder: 1` | **禁止** (粘性応力破壊・疑似 2D で全域に効く) | 段階起動 (§1.2) |
| 設計 YAML の `geometry.pw_ramp` (物理壁の上流に δ_r をランプで足す) | **2026-10-06 まで既定の作り方** (2026-10-07 からコードの既定は `pw_upstream: poly`。`pw_ramp` を書くなら `pw_upstream: ramp` の明示が必須、`poly` と併記すると例外) | `pw_upstream: poly` (キー無しと同じ)、`pw_ramp` を消す。旧壁の再現・case/45 の生産 (U4 で採否を決めるまで) だけ `pw_upstream: ramp` + `pw_ramp` |
| `solve_rt` (CFD 前) の δ = 未較正・未平滑化の `integral_bl` 直接 | **2026-10-06 まで** (生産の壁と出口半径が 2.4 mm ずれた) | 2026-10-07 から `prepare_ns` と同じ経路 (コード側、キー無し)。旧値の再現はできない |
| `nodeAxisDirichlet` / `nodeMidpointFx` / `nodeValueAtNode` / `nodeReconEdgeMidpoint` / `nodeAxisUrDirichlet` | 廃止 (2026-08-16, 書くと起動エラー) | node は固定スキーム (§1) |
| `turbulence: {LESorRANS: 2, RANSmodel: 1}` 旧キー体系 | 旧 config に残存 | `turbulence: {model: "sst", ...}` |
| `space.slauWallNormalChi` 省略 = 0 | **2026-09-25 まで既定** (2026-09-26 から auto = node+nodeWallDirichlet+SLAU で 1、plan [convection-slau-wall-normal-chi-default](../plans/accepted/convection-slau-wall-normal-chi-default.md)) | 旧結果の再現は `slauWallNormalChi: 0` を明記 |
| node の `mesh.scalarGradient` 省略 = gg | **2026-09-26 まで既定** (2026-09-27 から node の既定は lsq、plan [gradient-scalar-lsq-unification](../plans/accepted/gradient-scalar-lsq-unification.md) #6) | 旧結果の再現は `mesh.scalarGradient: gg` を明記。**切り替え日をまたぐ run は途中から再開せず最初から回し直す** (段の区間判定が gg 段と lsq 段をつなぐため、ユーザ決定 2026-09-27) |
| `sstOmegaProdFromPk: 0` / `sstSigmaBlend: 0` | 2026-09-08 まで既定 | 既定 1 (旧挙動が要るときだけ 0 明記) |
| `sstEnergyKSource: 1` / `sstIsotropicStress: 1` (分離型) | 非推奨 (境界未完備・離散保存せず) | 必要なら `sstEnergyIncludesK: 1` (opt-in) |
| 定常 + 陽解法 (`timeIntegration: 3`, `unsteady: 0`) | 非推奨 (局所 dt で不安定) | 陰解法 11 + blockDPLUR、または `unsteady: 1` |
| `thermoHrefTemp` 未指定 (絶対 h 基準) の多成分 TP | 発散要因 | `thermoHrefTemp: 298.15` |
| 同一 3D メッシュの段間引き継ぎに `interp_field.py` | 禁止 (z 列混同) | index コピー |
| 超音速 node 出口に `outlet_statPress` (Ps ≪ 実圧) | unstart の原因 | `outflow` か Ps 一致 |
| `output` 未指定 = 全量出力 | 2026-09-08 まで | `output.level` 既定 1 (全量は `level: 2`) |
| k/ω 拡散の絶対ゼロ割ガード 1e-12 [m³] | バグ (2026-09-08 修正) | 相対ガード (コード側、キー無し) |
| `wall_dist` を双対重心から測る変換 | バグ (2026-09-08 修正) | ノード座標 (コード側) |
| 内蔵種の LJ = #14 前の単一の値 (GRI-Mech 3.0 と Svehla 1962 の混在、`physProp.ljSource` 無し) | 2026-10-01 まで既定 | 既定 `[gri30, svehla1962]` (書かない)。旧 run の再現・継続だけ `ljSource: [legacy_v1]` |

## 付録 A. 機種と段ごとの推奨設定 (生成区間)

下の表は `solver_density_cuda/tools/recipes/recommended_settings.yaml` から `python3 solver_density_cuda/tools/recipe_settings.py --doc --write procedures/recommended-settings.md` で生成した (手で編集しない。試験 `solver_density_cuda/tests/unit/test_recipe_settings.py` が YAML との一致を確かめる)。設計チェーンの機種のプロファイルは 2026-10-14 の生成器の書き写し (§0.1)。共通プロファイル (`common_*`) は §1・§1.1・§1.2・§2 の値。

<!-- BEGIN GENERATED: recipe_settings.py --doc (手で編集しない。正本は solver_density_cuda/tools/recipes/recommended_settings.yaml) -->

推奨設定のファイルの版 2。プロファイルの一覧 (`python3 solver_density_cuda/tools/recipe_settings.py list`):

| プロファイル | 機種 | 適用条件 | 段 | 根拠 |
| --- | --- | --- | --- | --- |
| `common_node_ns_sst` | 機種の決まっていない手書きの node NS (SST・低 Re 壁) の定常 run | discretization node, time steady, physics ns_sst | soft → mid → main | §1・§1.1・§1.2・§2 2026-09-08 (§1・§2)・2026-09-07 (§1.1・§1.2) |
| `common_node_euler` | 機種の決まっていない手書きの node Euler (slip 壁) の定常 run | discretization node, time steady, physics euler | soft → mid → main | §1・§1.1・§1.2 2026-09-08 (§1)・2026-09-07 (§1.1・§1.2) |
| `bell_cpg_ns_sst` | ③ベル CPG (node SST・断熱壁、opt/driver の CPG 評価) | dimension 2d_axisym, discretization node, axisymmetric true, physics ns_sst, gas cpg, wall noslip_adiabatic, mesh structured, time steady | soft → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `bell_tp_ns_sst` | ③ベル TP (凍結組成の lump、node SST・低 Re・断熱壁) | dimension 2d_axisym, discretization node, axisymmetric true, physics ns_sst, gas tp, wall noslip_adiabatic, mesh structured, time steady | soft → mid → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `bell_tp_euler` | ③ベル TP Euler (slip 壁、(M) 特性曲線モードの評価) | dimension 2d_axisym, discretization node, axisymmetric true, physics euler, gas tp, wall slip, mesh structured, time steady | soft → mid → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `wt_modef_euler_cell` | ① 風洞 モード F (Euler、cell) | dimension 2d_axisym, discretization cell, axisymmetric true, physics euler, gas cpg, wall slip, mesh structured, time steady | soft → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `wt_feedback_euler_node` | ① 風洞 B10 の帰還 (feedback/euler_loop、Euler node、段階起動なし) | dimension 2d_axisym, discretization node, axisymmetric true, physics euler, gas cpg, wall slip, mesh structured, time steady | main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `axismach_euler_tp` | ① 風洞 axis-Mach Euler (TP、slip 壁) | dimension 2d_axisym, discretization node, axisymmetric true, physics euler, gas tp, wall slip, mesh structured, time steady | soft → mid → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `axismach_ns_sst_tp` | ① 風洞 axis-Mach NS (TP、SST 低 Re、等温壁、main_line なし) | dimension 2d_axisym, discretization node, axisymmetric true, physics ns_sst, gas tp, wall isothermal, mesh structured, time steady | soft → mid → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `axismach_ns_sst_tp_b0` | ① 風洞 axis-Mach NS (TP、SST 低 Re、等温壁、本段に main_line の B0) | dimension 2d_axisym, discretization node, axisymmetric true, physics ns_sst, gas tp, wall isothermal, mesh structured, time steady | soft → mid → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `walldriven_euler_cpg` | ① 風洞 wall-driven (Euler、node) | dimension 2d_axisym, discretization node, axisymmetric true, physics euler, gas cpg, wall slip, mesh structured, time steady | soft → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `sern2d_finite_base_frozen_tp_sst` | ⑤ SERN 2D 有限ベース (frozen TP、node SST、等温壁、既定の作動点) | dimension 2d_planar, discretization node, axisymmetric false, physics ns_sst, gas tp, wall isothermal, mesh structured, time steady | warm_lam → soft → mid → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |
| `sern3d_frozen_tp_sst` | ⑤ SERN 3D (frozen TP、node SST、等温壁、既定の作動点) | dimension 3d, discretization node, axisymmetric false, physics ns_sst, gas tp, wall isothermal, mesh structured, time steady | warm_lam → soft → mid → main | plans/active/tooling-generator-defaults-and-launch-gate.md §4.0 段取り 1 (今の生成器の書き写し)・§5.1 #6 2026-10-14 |

**読み方**: 各表の値は `fixed` (宣言したプロファイルと違えば FAIL、理由つきの逸脱で通す)、`〜` の範囲と「任意」は `advisory` (外れたら WARN)。「(明示すべき)」の値は `explicit` (plan §5.1 #16: 生成器が省略していた、値に依存するキーと旧い実効値。照合は `fixed` と同じ)。**`fixed`・`explicit` の省略は FAIL** (省略の実効値は実行するバイナリで決まるので、solver の既定との一致で代用しない。省略を試すなら `--deviate key=~ --reason "..."`、実効値は未確定として記録)。本段 (main) 以外の列の「=」は本段と同じ値、空欄はその段では照合しない (本段の長さ `time.last.nStepOuter` は問題 YAML のキー)。表に無いキーは、問題 YAML の物理条件・入出力のキー (`problem_keys`) か無視するキー以外は「プロファイルに無いキー」として FAIL (逸脱を宣言すれば通す)。禁止キーは書けば FAIL。

**明示すべきキーの基準** (YAML の `explicit_keys`、plan §5.1 #16): 各機種・段で、下の適用条件が成り立ち、生成器の config が書いていないキーを、その段の `explicit` に旧い実効値で載せる (`recipes/regen_profiles.py`)。旧い実効値の出所は書き写しの時のソースで確かめてある (変わると書き写しの台本が止まる)。

| キー | 適用条件 | 旧い実効値 | 出所 | 値に依存する理由 |
| --- | --- | --- | --- | --- |
| `time.deltaT.implicitRelax` | time.timeIntegration 11, time.dualTime 0 | 1.0 | solverConfig.cpp:370 | 定常の陰解法の緩和 (§1 は 0.7、省略は 1.0 = 緩和なし) |
| `time.deltaT.lowMachPrecond` | 常に | 0 | solverConfig.cpp:479 | 低マッハ前処理 (流束と擬似時間の両方を変える) |
| `time.deltaT.blockDPLUR` | time.timeIntegration 11 | 0 | solverConfig.cpp:391 | 陰解法の型 (block / scalar DPLUR) |
| `time.nStepInner` | time.timeIntegration 11 | 1 | solverConfig.cpp:531 | 陰解法の内部反復の回数 |
| `time.deltaT.detectNaN` | 常に | 0 | solverConfig.cpp:460 | NaN の検出で止めるか (§1 は常時 1) |
| `time.deltaT.lineImplicit` | time.timeIntegration 11, time.deltaT.blockDPLUR 1 | 0 | solverConfig.cpp:372 | 線陰解法 (axis-Mach の B0 の本段が書く) |
| `time.deltaT.lineDtDirectional` | time.timeIntegration 11, time.deltaT.blockDPLUR 1 | 0 | solverConfig.cpp:389 | 線方向の擬似時間刻み (axis-Mach の B0 の本段が書く) |
| `time.deltaT.lineDtDirectionalCap` | time.deltaT.lineDtDirectional {not: [0, 書かない]} | 0.0 | solverConfig.cpp:390 | 線方向の擬似時間刻みの上限 (lineDtDirectional が 0 でないときだけ効く) |
| `time.deltaT.implicitThermalJacobian` | time.timeIntegration 11, time.deltaT.blockDPLUR 1, physProp.viscMethod {not: [0]} | 0 | solverConfig.cpp:453 | エネルギー行の熱伝導のヤコビアン (粘性のときだけ。axis-Mach の B0 の本段が書く) |
| `space.limiterScaled` | mesh.discretization node, space.limiter {not: [0]} | 1 | solverConfig.cpp:692 | リミッタの評価点と無次元化 (既定 1 は 2026-09-20 から。切り替えで SERN の目的量が +0.54 %) |
| `space.venkatK` | mesh.discretization node, space.limiter 2, space.limiterScaled [1, 書かない] | 0.05 | solverConfig.cpp:706 | Venkatakrishnan の K (既定は limiterScaled で変わる) |
| `mesh.scalarGradient` | mesh.discretization node, applies.physics ns_sst | lsq | solverConfig.cpp:255-256 | k/ω・化学種の勾配 (node の省略は 2026-09-27 から lsq、それ以前は gg) |
| `mesh.scalarGradient` | mesh.discretization node, species_count {min: 2} | lsq | solverConfig.cpp:255-256 | k/ω・化学種の勾配 (node の省略は 2026-09-27 から lsq、それ以前は gg) |
| `space.slauWallNormalChi` | mesh.discretization node, solver [SLAU, SLAU2], applies.wall [isothermal, noslip_adiabatic], mesh.nodeWallDirichlet [1, 書かない] | 1 | solverConfig.cpp:731-740 | 壁隣接面の SLAU の χ (省略は auto。2026-09-26 から node ∧ nodeWallDirichlet 1 ∧ SLAU で 1、それ以前は 0。作用するのは wall・wall_isothermal の節点だけ) |
| `space.slauWallNormalChi` | mesh.discretization node, solver [SLAU, SLAU2], applies.wall [isothermal, noslip_adiabatic], mesh.nodeWallDirichlet {not: [1, 書かない]} | 0 | solverConfig.cpp:731-740 | 壁隣接面の SLAU の χ (省略は auto。2026-09-26 から node ∧ nodeWallDirichlet 1 ∧ SLAU で 1、それ以前は 0。作用するのは wall・wall_isothermal の節点だけ) |
| `mesh.nodeWallDirichlet` | mesh.discretization node, applies.wall [isothermal, noslip_adiabatic] | 1 | solverConfig.hpp:472 | node の壁の速度の Dirichlet (wall・wall_isothermal の節点だけに作用する) |
| `mesh.axisymMethod` | applies.axisymmetric true, physProp.axisymMethod [書かない] | 0 | solverConfig.hpp:432 | 軸対称の定式化 |
| `mesh.axisSegmentRWeight` | mesh.discretization node, applies.axisymmetric true, mesh.axisymMethod [0, 書かない], physProp.axisymMethod [書かない], mesh.axisRFloor [0, 書かない] | 1 | solverConfig.hpp:454 | 軸対称の r 重みを双対面の区間ごとに掛ける (既定 1 は 2026-10-10 から。node・axisymMethod 0・axisRFloor 0 のときだけ効く) |
| `physProp.thermoFloat` | physProp.thermalMethod 2, physProp.thermoHrefTemp {gt: 0} | 1 | solverConfig.hpp:627・solverConfig.cpp:1094-1101 | TP の温度反転 (既定 1 は 2026-09-12 から。thermoHrefTemp が無いと 0 に落ちる) |
| `turbulence.turbulentPrandtl` | applies.physics ns_sst | 0.85 | solverConfig.cpp:920 | 乱流プラントル数 (§2 は 0.9、省略は 0.85) |
| `turbulence.scalarDiffusion` | applies.physics ns_sst | 1 | solverConfig.cpp:842 | k/ω の拡散 |
| `turbulence.dilatationCorrection` | applies.physics ns_sst | 2 | solverConfig.cpp:843 | SST の膨張補正 (ノズルの境界層の厚さを変える) |
| `turbulence.katoLaunder` | applies.physics ns_sst | 0 | solverConfig.cpp:851 | Kato–Launder の生成項の修正 |
| `turbulence.wallTreatmentSST` | applies.physics ns_sst | 0 | solverConfig.cpp:875 | SST の壁の扱い (既定 0 は 2026-09-20 から) |
| `turbulence.sstNodeWallKPin` | applies.physics ns_sst, mesh.discretization node | 1 | solverConfig.cpp:852 | node の壁の k の固定 (既定 1 は 2026-09-08 から) |
| `turbulence.sstOmegaProdFromPk` | applies.physics ns_sst | 1 | solverConfig.cpp:853 | ω の生成を P_k から作る (既定 1 は 2026-09-08 から) |
| `turbulence.sstSigmaBlend` | applies.physics ns_sst | 1 | solverConfig.cpp:854 | SST の σ の混合 (既定 1 は 2026-09-08 から) |
| `turbulence.sstEnergyIncludesK` | applies.physics ns_sst | 0 | solverConfig.cpp:857 | 全エネルギーに k を含めるか (既定 0 を 2026-09-08 に確定) |

#### `common_node_ns_sst`

recommended-settings.md §1・§1.1・§1.2・§2。本段の cfl_pseudo は目安 4〜6 (advisory)。

| キー | soft | mid | main |
| --- | --- | --- | --- |
| `mesh.axisCentroidShift` | = | = | 任意 |
| `mesh.axisymMethod` | = | = | 0 (書くなら) |
| `mesh.discretization` | = | = | node |
| `mesh.isAxisymmetric` | = | = | 任意 |
| `mesh.nodeInletCornerWall` | = | = | 1 |
| `mesh.nodeWallDirichlet` | = | = | 1 |
| `output.level` | = | = | 1 / 2 のどれか |
| `physProp.thermCondMethod` | = | = | 任意 |
| `physProp.thermalMethod` | = | = | 0 / 2 のどれか |
| `physProp.thermoFloat` | = | = | 1 (書くなら) |
| `physProp.viscMethod` | = | = | 1 / 2 のどれか |
| `solver` | = | = | SLAU |
| `space.convMethod` | 0 | 0 | 1 |
| `space.limiter` | 0 | 0 | 2 |
| `space.limiterScaled` | = | = | 1 (書くなら) |
| `space.venkatK` | = | = | 0.05 (書くなら) |
| `time.deltaT.blockDPLUR` | = | = | 1 |
| `time.deltaT.cfl` | = | = | 任意 |
| `time.deltaT.cfl_pseudo` | 0.5 | 1.0 | 4〜6 (目安) |
| `time.deltaT.control` | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | 0.7 |
| `time.deltaT.lowMachPrecond` | = | = | 0 (書くなら) |
| `time.dualTime` | = | = | 0 |
| `time.last.nStepOuter` | 2000 | 2000 | 8000〜12000 (目安) |
| `time.nStepInner` | 10 | 10 | 4 |
| `time.timeIntegration` | = | = | 11 |
| `time.unsteady` | = | = | 0 |
| `turbulence.dilatationCorrection` | = | = | 2 |
| `turbulence.katoLaunder` | = | = | 1 |
| `turbulence.model` | = | = | sst |
| `turbulence.scalarDiffusion` | = | = | 1 |
| `turbulence.turbulentPrandtl` | = | = | 0.9 |
| `turbulence.wallTreatmentSST` | = | = | 0 |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

#### `common_node_euler`

recommended-settings.md §1 (cfl_pseudo の目安は Euler 設計評価 4)・§1.1 (Euler の壁は slip)。viscMethod 0 は定数粘性なので visc・thermCond も 0。§1 は soft・mid の nStepInner を node NS (10) にしか決めていないので任意。

| キー | soft | mid | main |
| --- | --- | --- | --- |
| `mesh.axisCentroidShift` | = | = | 任意 |
| `mesh.axisymMethod` | = | = | 0 (書くなら) |
| `mesh.discretization` | = | = | node |
| `mesh.isAxisymmetric` | = | = | 任意 |
| `mesh.nodeInletCornerWall` | = | = | 任意 |
| `mesh.nodeWallDirichlet` | = | = | 任意 |
| `output.level` | = | = | 1 / 2 のどれか |
| `physProp.thermCond` | = | = | 0.0 |
| `physProp.thermalMethod` | = | = | 0 / 2 のどれか |
| `physProp.thermoFloat` | = | = | 1 (書くなら) |
| `physProp.visc` | = | = | 0.0 |
| `physProp.viscMethod` | = | = | 0 |
| `solver` | = | = | SLAU |
| `space.convMethod` | 0 | 0 | 1 |
| `space.limiter` | 0 | 0 | 2 |
| `space.limiterScaled` | = | = | 1 (書くなら) |
| `space.venkatK` | = | = | 0.05 (書くなら) |
| `time.deltaT.blockDPLUR` | = | = | 1 |
| `time.deltaT.cfl` | = | = | 任意 |
| `time.deltaT.cfl_pseudo` | 0.5 | 1.0 | 4 (目安) |
| `time.deltaT.control` | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | 0.7 |
| `time.deltaT.lowMachPrecond` | = | = | 0 (書くなら) |
| `time.dualTime` | = | = | 0 |
| `time.last.nStepOuter` | 2000 | 2000 | 8000〜12000 (目安) |
| `time.nStepInner` | 任意 | 任意 | 4 |
| `time.timeIntegration` | = | = | 11 |
| `time.unsteady` | = | = | 0 |
| `turbulence.model` | = | = | none |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

#### `bell_cpg_ns_sst`

書き写し: generator `runner._solver_config + opt/driver の soft (別ディレクトリ)`、problem `case/40.nozzle_design_tool/problem_bell_top.yaml`、call `opt.driver.Campaign.evaluate (prepare・forge・collect を差し替え)`

| キー | soft | main |
| --- | --- | --- |
| `mesh.axisCentroidShift` | = | 1 |
| `mesh.axisSegmentRWeight` | = | 1 (明示すべき) |
| `mesh.axisymMethod` | = | 0 (明示すべき) |
| `mesh.discretization` | = | node |
| `mesh.isAxisymmetric` | = | 1 |
| `mesh.nodeWallDirichlet` | = | 1 |
| `mesh.scalarGradient` | = | lsq (明示すべき) |
| `physProp.thermalMethod` | = | 0 |
| `physProp.viscMethod` | = | 1 |
| `solver` | = | SLAU |
| `space.convMethod` | 0 | 1 |
| `space.limiter` | = | 2 |
| `space.limiterScaled` | = | 1 (明示すべき) |
| `space.slauWallNormalChi` | = | 1 (明示すべき) |
| `space.venkatK` | = | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | = | 1 |
| `time.deltaT.cfl` | 0.5 | 4.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | 4.0 |
| `time.deltaT.control` | = | 1 |
| `time.deltaT.detectNaN` | = | 1 |
| `time.deltaT.implicitRelax` | = | 1.0 (明示すべき) |
| `time.deltaT.implicitThermalJacobian` | = | 0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | 0 |
| `time.dualTime` | = | 0 |
| `time.last.nStepOuter` | 3000 |  |
| `time.nStepInner` | = | 5 |
| `time.timeIntegration` | = | 11 |
| `time.unsteady` | = | 0 |
| `turbulence.dilatationCorrection` | = | 2 |
| `turbulence.katoLaunder` | = | 1 |
| `turbulence.model` | = | sst |
| `turbulence.scalarDiffusion` | = | 1 |
| `turbulence.sstEnergyIncludesK` | = | 0 (明示すべき) |
| `turbulence.sstNodeWallKPin` | = | 1 (明示すべき) |
| `turbulence.sstOmegaProdFromPk` | = | 1 (明示すべき) |
| `turbulence.sstSigmaBlend` | = | 1 (明示すべき) |
| `turbulence.turbulentPrandtl` | = | 0.85 (明示すべき) |
| `turbulence.wallTreatmentSST` | = | 0 (明示すべき) |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_ns_sst`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `mesh.nodeInletCornerWall` | soft, main | 1 | 書かない (solver 既定 不明) | 未 |
| `space.limiter` | soft | 0 | 2 | 未 |
| `time.deltaT.implicitRelax` | soft, main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.last.nStepOuter` | soft | 2000 | 3000 | 未 |
| `time.nStepInner` | soft | 10 | 5 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |
| `turbulence.turbulentPrandtl` | soft, main | 0.9 | 書かない (solver 既定 0.85) | 未 |

#### `bell_tp_ns_sst`

書き写し: generator `runner._solver_config_tp (TP_STAGES)`、problem `case/68.bell_two_mode/problem_bell_tp_m_ns_l100_corr.yaml`、call `runner.run_staged_tp (prepare・forge を差し替え)`

| キー | soft | mid | main |
| --- | --- | --- | --- |
| `mesh.axisCentroidShift` | = | = | 1 |
| `mesh.axisSegmentRWeight` | = | = | 1 (明示すべき) |
| `mesh.axisymMethod` | = | = | 0 (明示すべき) |
| `mesh.discretization` | = | = | node |
| `mesh.isAxisymmetric` | = | = | 1 |
| `mesh.nodeInletCornerWall` | = | = | 1 |
| `mesh.nodeWallDirichlet` | = | = | 1 |
| `mesh.scalarGradient` | = | = | lsq (明示すべき) |
| `output.level` | = | = | 1 |
| `physProp.thermalMethod` | = | = | 2 |
| `physProp.thermoFloat` | = | = | 1 (明示すべき) |
| `physProp.viscMethod` | = | = | 2 |
| `solver` | = | = | SLAU |
| `space.convMethod` | 0 | 0 | 1 |
| `space.limiter` | 0 | 0 | 2 |
| `space.limiterScaled` |  |  | 1 (明示すべき) |
| `space.slauWallNormalChi` | = | = | 1 (明示すべき) |
| `space.venkatK` |  |  | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | = | = | 1 |
| `time.deltaT.cfl` | 0.5 | 1.0 | 4.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | 1.0 | 4.0 |
| `time.deltaT.control` | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | 0.7 |
| `time.deltaT.implicitThermalJacobian` | = | = | 0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | = | 0 (明示すべき) |
| `time.dualTime` | = | = | 0 |
| `time.last.nStepOuter` | 2000 | 2000 |  |
| `time.nStepInner` | 10 | 10 | 4 |
| `time.timeIntegration` | = | = | 11 |
| `time.unsteady` | = | = | 0 |
| `turbulence.dilatationCorrection` | = | = | 2 |
| `turbulence.katoLaunder` | = | = | 1 |
| `turbulence.model` | = | = | sst |
| `turbulence.scalarDiffusion` | = | = | 1 |
| `turbulence.sstEnergyIncludesK` | = | = | 0 (明示すべき) |
| `turbulence.sstNodeWallKPin` | = | = | 1 (明示すべき) |
| `turbulence.sstOmegaProdFromPk` | = | = | 1 (明示すべき) |
| `turbulence.sstSigmaBlend` | = | = | 1 (明示すべき) |
| `turbulence.turbulentPrandtl` | = | = | 0.9 |
| `turbulence.wallTreatmentSST` | = | = | 0 |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

#### `bell_tp_euler`

書き写し: generator `runner._solver_config_tp (TP_STAGES、evaluate.model: euler)`、problem `case/68.bell_two_mode/problem_bell_tp_m_euler_l100.yaml`、call `runner.run_staged_tp (prepare・forge を差し替え)`

| キー | soft | mid | main |
| --- | --- | --- | --- |
| `mesh.axisCentroidShift` | = | = | 1 |
| `mesh.axisSegmentRWeight` | = | = | 1 (明示すべき) |
| `mesh.axisymMethod` | = | = | 0 (明示すべき) |
| `mesh.discretization` | = | = | node |
| `mesh.isAxisymmetric` | = | = | 1 |
| `output.level` | = | = | 1 |
| `physProp.thermCond` | = | = | 0.0 |
| `physProp.thermalMethod` | = | = | 2 |
| `physProp.thermoFloat` | = | = | 1 (明示すべき) |
| `physProp.visc` | = | = | 0.0 |
| `physProp.viscMethod` | = | = | 0 |
| `solver` | = | = | SLAU |
| `space.convMethod` | 0 | 0 | 1 |
| `space.limiter` | 0 | 0 | 2 |
| `space.limiterScaled` |  |  | 1 (明示すべき) |
| `space.venkatK` |  |  | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | = | = | 1 |
| `time.deltaT.cfl` | 0.5 | 1.0 | 4.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | 1.0 | 4.0 |
| `time.deltaT.control` | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | 0.7 |
| `time.deltaT.lineDtDirectional` | = | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | = | 0 (明示すべき) |
| `time.dualTime` | = | = | 0 |
| `time.last.nStepOuter` | 2000 | 2000 |  |
| `time.nStepInner` | 10 | 10 | 4 |
| `time.timeIntegration` | = | = | 11 |
| `time.unsteady` | = | = | 0 |
| `turbulence.model` | = | = | none |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

#### `wt_modef_euler_cell`

書き写し: generator `runner_wt._config_euler (cell 固定) + run_staged の文字列置換`、problem `case/41.wind_tunnel_design/problem_m4.yaml`、call `runner_wt.main (prepare・forge を差し替え)`

| キー | soft | main |
| --- | --- | --- |
| `mesh.axisymMethod` | = | 0 (明示すべき) |
| `mesh.discretization` | = | cell |
| `mesh.isAxisymmetric` | = | 1 |
| `physProp.thermCond` | = | 0.0 |
| `physProp.thermalMethod` | = | 0 |
| `physProp.visc` | = | 0.0 |
| `physProp.viscMethod` | = | 0 |
| `solver` | = | SLAU |
| `space.convMethod` | 0 | 1 |
| `space.limiter` | = | 2 |
| `time.deltaT.blockDPLUR` | = | 1 |
| `time.deltaT.cfl` | 0.5 | 4.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | 4.0 |
| `time.deltaT.control` | = | 1 |
| `time.deltaT.detectNaN` | = | 1 |
| `time.deltaT.implicitRelax` | = | 1.0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | 0 |
| `time.dualTime` | = | 0 |
| `time.last.nStepOuter` | 3000 |  |
| `time.nStepInner` | 10 | 5 |
| `time.timeIntegration` | = | 11 |
| `time.unsteady` | = | 0 |
| `turbulence.model` | = | none |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_euler`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `mesh.discretization` | soft, main | node | cell | 未 |
| `space.limiter` | soft | 0 | 2 | 未 |
| `time.deltaT.implicitRelax` | soft, main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.last.nStepOuter` | soft | 2000 | 3000 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |

#### `wt_feedback_euler_node`

書き写し: generator `runner_wt._config_euler_node (euler_loop の cold のパス)`、problem `case/41.wind_tunnel_design/problem_m4_b10.yaml`、call `euler_loop の `_cfg(n, n // 4)` (n = 24000)`

| キー | main |
| --- | --- |
| `mesh.axisCentroidShift` | 1 |
| `mesh.axisSegmentRWeight` | 1 (明示すべき) |
| `mesh.axisymMethod` | 0 (明示すべき) |
| `mesh.discretization` | node |
| `mesh.isAxisymmetric` | 1 |
| `physProp.thermCond` | 0.0 |
| `physProp.thermalMethod` | 0 |
| `physProp.visc` | 0.0 |
| `physProp.viscMethod` | 0 |
| `solver` | SLAU |
| `space.convMethod` | 1 |
| `space.limiter` | 2 |
| `space.limiterScaled` | 1 (明示すべき) |
| `space.venkatK` | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | 1 |
| `time.deltaT.cfl` | 4.0 |
| `time.deltaT.cfl_pseudo` | 4.0 |
| `time.deltaT.control` | 1 |
| `time.deltaT.detectNaN` | 1 |
| `time.deltaT.implicitRelax` | 1.0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | 0 |
| `time.dualTime` | 0 |
| `time.nStepInner` | 5 |
| `time.timeIntegration` | 11 |
| `time.unsteady` | 0 |
| `turbulence.model` | none |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_euler`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `time.deltaT.implicitRelax` | main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |

#### `axismach_euler_tp`

書き写し: generator `runner_wt._config_euler_node → runner_axismach._apply_gas_to_config + run_staged (_cfg_set)`、problem `case/45.isobutane_m6_d155/problem_d155_euler_pin_G1_recal_mono.yaml`、call `runner_axismach.main (prepare・forge を差し替え)`

| キー | soft | mid | main |
| --- | --- | --- | --- |
| `mesh.axisCentroidShift` | = | = | 1 |
| `mesh.axisSegmentRWeight` | = | = | 1 (明示すべき) |
| `mesh.axisymMethod` | = | = | 0 (明示すべき) |
| `mesh.discretization` | = | = | node |
| `mesh.isAxisymmetric` | = | = | 1 |
| `mesh.scalarGradient` | = | = | lsq (明示すべき) |
| `physProp.thermCond` | = | = | 0.0 |
| `physProp.thermalMethod` | = | = | 2 |
| `physProp.thermoFloat` | = | = | 1 (明示すべき) |
| `physProp.visc` | = | = | 0.0 |
| `physProp.viscMethod` | = | = | 0 |
| `solver` | = | = | SLAU |
| `space.convMethod` | 0 | = | 1 |
| `space.limiter` | = | = | 2 |
| `space.limiterScaled` | = | = | 1 (明示すべき) |
| `space.venkatK` | = | = | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | = | = | 1 |
| `time.deltaT.cfl` | 0.5 | 1.0 | 2.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | 1.0 | 2.0 |
| `time.deltaT.control` | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | 1.0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | = | 0 |
| `time.dualTime` | = | = | 0 |
| `time.last.nStepOuter` | 3000 | 3000 |  |
| `time.nStepInner` | = | = | 5 |
| `time.timeIntegration` | = | = | 11 |
| `time.unsteady` | = | = | 0 |
| `turbulence.model` | = | = | none |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_euler`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `space.convMethod` | mid | 0 | 1 | 未 |
| `space.limiter` | soft, mid | 0 | 2 | 未 |
| `time.deltaT.cfl_pseudo` | main | 4 (目安) | 2.0 | 未 |
| `time.deltaT.implicitRelax` | soft, mid, main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.last.nStepOuter` | soft, mid | 2000 | 3000 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |

#### `axismach_ns_sst_tp`

書き写し: generator `runner_wt._config_sst_node → _apply_gas_to_config(viscous) + ns_stage_configs`、problem `case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw300.yaml`、call `runner_axismach._ns_configs → run_staged_ns(stages=full) (forge を差し替え)`

| キー | soft | mid | main |
| --- | --- | --- | --- |
| `mesh.axisCentroidShift` | = | = | 1 |
| `mesh.axisSegmentRWeight` | = | = | 1 (明示すべき) |
| `mesh.axisymMethod` | = | = | 0 (明示すべき) |
| `mesh.discretization` | = | = | node |
| `mesh.isAxisymmetric` | = | = | 1 |
| `mesh.nodeWallDirichlet` | = | = | 1 |
| `mesh.scalarGradient` | = | = | lsq (明示すべき) |
| `physProp.thermalMethod` | = | = | 2 |
| `physProp.thermoFloat` | = | = | 1 (明示すべき) |
| `physProp.viscMethod` | = | = | 2 |
| `solver` | = | = | SLAU |
| `space.convMethod` | 0 | 0 | 1 |
| `space.limiter` | = | = | 2 |
| `space.limiterScaled` | = | = | 1 (明示すべき) |
| `space.slauWallNormalChi` | = | = | 1 (明示すべき) |
| `space.venkatK` | = | = | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | = | = | 1 |
| `time.deltaT.cfl` | 0.5 | = | 1.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | = | 1.0 |
| `time.deltaT.control` | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | 1.0 (明示すべき) |
| `time.deltaT.implicitThermalJacobian` | = | = | 0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | = | 0 |
| `time.dualTime` | = | = | 0 |
| `time.last.nStepOuter` | 3000 | 3000 |  |
| `time.nStepInner` | 10 | 10 | 5 |
| `time.timeIntegration` | = | = | 11 |
| `time.unsteady` | = | = | 0 |
| `turbulence.dilatationCorrection` | = | = | 2 |
| `turbulence.katoLaunder` | = | = | 1 |
| `turbulence.model` | = | = | sst |
| `turbulence.scalarDiffusion` | = | = | 1 |
| `turbulence.sstEnergyIncludesK` | = | = | 0 (明示すべき) |
| `turbulence.sstNodeWallKPin` | = | = | 1 (明示すべき) |
| `turbulence.sstOmegaProdFromPk` | = | = | 1 (明示すべき) |
| `turbulence.sstSigmaBlend` | = | = | 1 (明示すべき) |
| `turbulence.turbulentPrandtl` | = | = | 0.9 |
| `turbulence.wallTreatmentSST` | = | = | 0 |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_ns_sst`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `mesh.nodeInletCornerWall` | soft, mid, main | 1 | 書かない (solver 既定 不明) | 未 |
| `space.limiter` | soft, mid | 0 | 2 | 未 |
| `time.deltaT.cfl_pseudo` | main | 4〜6 (目安) | 1.0 | 未 |
| `time.deltaT.implicitRelax` | soft, mid, main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.last.nStepOuter` | soft, mid | 2000 | 3000 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |

#### `axismach_ns_sst_tp_b0`

書き写し: generator `runner_wt._config_sst_node → _apply_gas_to_config(viscous) → _apply_main_line + ns_stage_configs (PRESTAGE_CONFIG)`、problem `case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml`、call `runner_axismach._ns_configs → run_staged_ns(stages=full) (forge を差し替え)`

| キー | soft | mid | main |
| --- | --- | --- | --- |
| `mesh.axisCentroidShift` | = | = | 1 |
| `mesh.axisSegmentRWeight` | = | = | 1 (明示すべき) |
| `mesh.axisymMethod` | = | = | 0 (明示すべき) |
| `mesh.discretization` | = | = | node |
| `mesh.isAxisymmetric` | = | = | 1 |
| `mesh.nodeWallDirichlet` | = | = | 1 |
| `mesh.scalarGradient` | = | = | lsq (明示すべき) |
| `output.level` |  |  | 1 |
| `physProp.thermalMethod` | = | = | 2 |
| `physProp.thermoFloat` | = | = | 1 (明示すべき) |
| `physProp.viscMethod` | = | = | 2 |
| `solver` | = | = | SLAU |
| `space.convMethod` | 0 | 0 | 1 |
| `space.limiter` | = | = | 2 |
| `space.limiterScaled` | = | = | 1 (明示すべき) |
| `space.slauWallNormalChi` | = | = | 1 (明示すべき) |
| `space.venkatK` | = | = | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | = | = | 1 |
| `time.deltaT.cfl` | 0.5 | 1.0 | 4.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | 1.0 | 4.0 |
| `time.deltaT.control` | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | 0.7 |
| `time.deltaT.implicitThermalJacobian` | 0 (明示すべき) | 0 (明示すべき) | 5 |
| `time.deltaT.lineDtDirectional` | 0 (明示すべき) | 0 (明示すべき) | 1 |
| `time.deltaT.lineDtDirectionalCap` |  |  | 50.0 |
| `time.deltaT.lineImplicit` | 0 (明示すべき) | 0 (明示すべき) | 1 |
| `time.deltaT.lowMachPrecond` | = | = | 0 |
| `time.dualTime` | = | = | 0 |
| `time.last.nStepOuter` | 3000 | 3000 |  |
| `time.nStepInner` | 10 | 10 | 5 |
| `time.timeIntegration` | = | = | 11 |
| `time.unsteady` | = | = | 0 |
| `turbulence.dilatationCorrection` | = | = | 2 |
| `turbulence.katoLaunder` | = | = | 1 |
| `turbulence.model` | = | = | sst |
| `turbulence.scalarDiffusion` | = | = | 1 |
| `turbulence.sstEnergyIncludesK` | = | = | 0 (明示すべき) |
| `turbulence.sstNodeWallKPin` | = | = | 1 (明示すべき) |
| `turbulence.sstOmegaProdFromPk` | = | = | 1 (明示すべき) |
| `turbulence.sstSigmaBlend` | = | = | 1 (明示すべき) |
| `turbulence.turbulentPrandtl` | = | = | 0.9 |
| `turbulence.wallTreatmentSST` | = | = | 0 |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_ns_sst`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `mesh.nodeInletCornerWall` | soft, mid, main | 1 | 書かない (solver 既定 不明) | 未 |
| `space.limiter` | soft, mid | 0 | 2 | 未 |
| `time.last.nStepOuter` | soft, mid | 2000 | 3000 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |

#### `walldriven_euler_cpg`

書き写し: generator `runner_wt._config_euler_node + run_staged の文字列置換`、problem `case/41.wind_tunnel_design/problem_m4_walldriven.yaml`、call `runner_walldriven.main (prepare・forge を差し替え)`

| キー | soft | main |
| --- | --- | --- |
| `mesh.axisCentroidShift` | = | 1 |
| `mesh.axisSegmentRWeight` | = | 1 (明示すべき) |
| `mesh.axisymMethod` | = | 0 (明示すべき) |
| `mesh.discretization` | = | node |
| `mesh.isAxisymmetric` | = | 1 |
| `physProp.thermCond` | = | 0.0 |
| `physProp.thermalMethod` | = | 0 |
| `physProp.visc` | = | 0.0 |
| `physProp.viscMethod` | = | 0 |
| `solver` | = | SLAU |
| `space.convMethod` | 0 | 1 |
| `space.limiter` | = | 2 |
| `space.limiterScaled` | = | 1 (明示すべき) |
| `space.venkatK` | = | 0.05 (明示すべき) |
| `time.deltaT.blockDPLUR` | = | 1 |
| `time.deltaT.cfl` | 0.5 | 4.0 |
| `time.deltaT.cfl_pseudo` | 0.5 | 4.0 |
| `time.deltaT.control` | = | 1 |
| `time.deltaT.detectNaN` | = | 1 |
| `time.deltaT.implicitRelax` | = | 1.0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | 0 |
| `time.dualTime` | = | 0 |
| `time.last.nStepOuter` | 3000 |  |
| `time.nStepInner` | = | 5 |
| `time.timeIntegration` | = | 11 |
| `time.unsteady` | = | 0 |
| `turbulence.model` | = | none |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_euler`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `space.limiter` | soft | 0 | 2 | 未 |
| `time.deltaT.implicitRelax` | soft, main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.last.nStepOuter` | soft | 2000 | 3000 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |

#### `sern2d_finite_base_frozen_tp_sst`

書き写し: generator `runner_sern._solver_config + run_staged の正規表現・文字列置換 (問題の opt: の梯子)`、problem `case/46.sern_design/problem_moo_frozen_tp_cycle3op.yaml`、call `runner_sern.main (prepare・forge・restart を差し替え)`

| キー | warm_lam | soft | mid | main |
| --- | --- | --- | --- | --- |
| `mesh.discretization` | = | = | = | node |
| `mesh.isAxisymmetric` | = | = | = | 0 |
| `mesh.nodeWallDirichlet` | = | = | = | 1 |
| `mesh.scalarGradient` | = | = | = | lsq (明示すべき) |
| `physProp.thermalMethod` | = | = | = | 2 |
| `physProp.thermoFloat` | = | = | = | 1 (明示すべき) |
| `physProp.viscMethod` | = | = | = | 2 |
| `solver` | = | = | = | SLAU |
| `space.convMethod` | 0 | 0 | = | 1 |
| `space.limiter` | = | = | = | 2 |
| `space.limiterScaled` | = | = | = | 1 |
| `space.slauWallNormalChi` | = | = | = | 1 (明示すべき) |
| `space.venkatK` | = | = | = | 0.05 |
| `time.deltaT.blockDPLUR` | = | = | = | 1 |
| `time.deltaT.cfl` | 0.1 | 0.1 | 0.1 | 0.5 |
| `time.deltaT.cfl_pseudo` | 0.1 | 0.1 | 0.1 | 0.5 |
| `time.deltaT.control` | = | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | = | 1.0 (明示すべき) |
| `time.deltaT.implicitThermalJacobian` | = | = | = | 0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | = | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | = | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | = | = | 0 |
| `time.dualTime` | = | = | = | 0 |
| `time.last.nStepOuter` | 4000 | 4000 | 4000 |  |
| `time.nStepInner` | = | = | = | 5 |
| `time.timeIntegration` | = | = | = | 11 |
| `time.unsteady` | = | = | = | 0 |
| `turbulence.dilatationCorrection` |  | = | = | 2 |
| `turbulence.katoLaunder` |  | = | = | 1 |
| `turbulence.model` | none | = | = | sst |
| `turbulence.scalarDiffusion` |  | = | = | 1 |
| `turbulence.sstEnergyIncludesK` |  | = | = | 0 (明示すべき) |
| `turbulence.sstNodeWallKPin` |  | = | = | 1 (明示すべき) |
| `turbulence.sstOmegaProdFromPk` |  | = | = | 1 (明示すべき) |
| `turbulence.sstSigmaBlend` |  | = | = | 1 (明示すべき) |
| `turbulence.turbulentPrandtl` |  | = | = | 0.85 (明示すべき) |
| `turbulence.wallTreatmentSST` |  | = | = | 0 |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_ns_sst`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `mesh.nodeInletCornerWall` | soft, mid, main | 1 | 書かない (solver 既定 不明) | 未 |
| `space.convMethod` | mid | 0 | 1 | 未 |
| `space.limiter` | soft, mid | 0 | 2 | 未 |
| `time.deltaT.cfl_pseudo` | soft | 0.5 | 0.1 | 未 |
| `time.deltaT.cfl_pseudo` | mid | 1.0 | 0.1 | 未 |
| `time.deltaT.cfl_pseudo` | main | 4〜6 (目安) | 0.5 | 未 |
| `time.deltaT.implicitRelax` | soft, mid, main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.last.nStepOuter` | soft, mid | 2000 | 4000 | 未 |
| `time.nStepInner` | soft, mid | 10 | 5 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |
| `turbulence.turbulentPrandtl` | soft, mid, main | 0.9 | 書かない (solver 既定 0.85) | 未 |

#### `sern3d_frozen_tp_sst`

書き写し: generator `runner_sern._solver_config + run_staged の正規表現・文字列置換 (問題の opt: の梯子)`、problem `case/46.sern_design/problem_3d_prod_m6on_g2.yaml`、call `runner_sern3d.main (prepare・forge・restart を差し替え)`

| キー | warm_lam | soft | mid | main |
| --- | --- | --- | --- | --- |
| `mesh.discretization` | = | = | = | node |
| `mesh.isAxisymmetric` | = | = | = | 0 |
| `mesh.nodeWallDirichlet` | = | = | = | 1 |
| `mesh.scalarGradient` | = | = | = | lsq (明示すべき) |
| `physProp.thermalMethod` | = | = | = | 2 |
| `physProp.thermoFloat` | = | = | = | 1 (明示すべき) |
| `physProp.viscMethod` | = | = | = | 2 |
| `solver` | = | = | = | SLAU |
| `space.convMethod` | 0 | 0 | = | 1 |
| `space.limiter` | = | = | = | 2 |
| `space.limiterScaled` | = | = | = | 1 |
| `space.slauWallNormalChi` | = | = | = | 1 (明示すべき) |
| `space.venkatK` | = | = | = | 0.05 |
| `time.deltaT.blockDPLUR` | = | = | = | 1 |
| `time.deltaT.cfl` | 0.2 | 0.2 | 0.2 | 0.25 |
| `time.deltaT.cfl_pseudo` | 0.2 | 0.2 | 0.2 | 0.25 |
| `time.deltaT.control` | = | = | = | 1 |
| `time.deltaT.detectNaN` | = | = | = | 1 |
| `time.deltaT.implicitRelax` | = | = | = | 1.0 (明示すべき) |
| `time.deltaT.implicitThermalJacobian` | = | = | = | 0 (明示すべき) |
| `time.deltaT.lineDtDirectional` | = | = | = | 0 (明示すべき) |
| `time.deltaT.lineImplicit` | = | = | = | 0 (明示すべき) |
| `time.deltaT.lowMachPrecond` | = | = | = | 0 |
| `time.dualTime` | = | = | = | 0 |
| `time.last.nStepOuter` | 2000 | 2000 | 2000 |  |
| `time.nStepInner` | = | = | = | 5 |
| `time.timeIntegration` | = | = | = | 11 |
| `time.unsteady` | = | = | = | 0 |
| `turbulence.dilatationCorrection` |  | = | = | 2 |
| `turbulence.katoLaunder` |  | = | = | 1 |
| `turbulence.model` | none | = | = | sst |
| `turbulence.scalarDiffusion` |  | = | = | 1 |
| `turbulence.sstEnergyIncludesK` |  | = | = | 0 (明示すべき) |
| `turbulence.sstNodeWallKPin` |  | = | = | 1 (明示すべき) |
| `turbulence.sstOmegaProdFromPk` |  | = | = | 1 (明示すべき) |
| `turbulence.sstSigmaBlend` |  | = | = | 1 (明示すべき) |
| `turbulence.turbulentPrandtl` |  | = | = | 0.85 (明示すべき) |
| `turbulence.wallTreatmentSST` |  | = | = | 0 |

禁止: `mesh.bndFirstOrder` (使用禁止 (粘性応力を壊し、疑似 2D では全域に効く; AGENTS.md))

照合しない: `time.deltaT.dt`、`time.deltaT.dt_min`、`time.deltaT.dt_max` (定常 (unsteady 0) では効かない (plan §3.2))

§1 (`common_node_ns_sst`) との差:

| キー | 段 | §1 | この機種 | 受入 |
| --- | --- | --- | --- | --- |
| `mesh.nodeInletCornerWall` | soft, mid, main | 1 | 書かない (solver 既定 不明) | 未 |
| `space.convMethod` | mid | 0 | 1 | 未 |
| `space.limiter` | soft, mid | 0 | 2 | 未 |
| `time.deltaT.cfl_pseudo` | soft | 0.5 | 0.2 | 未 |
| `time.deltaT.cfl_pseudo` | mid | 1.0 | 0.2 | 未 |
| `time.deltaT.cfl_pseudo` | main | 4〜6 (目安) | 0.25 | 未 |
| `time.deltaT.implicitRelax` | soft, mid, main | 0.7 | 書かない (solver 既定 1.0) | 未 |
| `time.nStepInner` | soft, mid | 10 | 5 | 未 |
| `time.nStepInner` | main | 4 | 5 | 未 |
| `turbulence.turbulentPrandtl` | soft, mid, main | 0.9 | 書かない (solver 既定 0.85) | 未 |

<!-- END GENERATED: recipe_settings.py --doc -->

## 変更ログ

- `2026-10-15` — §0.1 の照合の規則を直した: `fixed`・`explicit` の省略は出所によらず FAIL (今のソースの既定との一致で通さない)、省略を試すのは `key=~` + 理由の逸脱 (実効値は未確定として記録)。機種・段ごとの明示すべきキー (`explicit`、旧い実効値と出所つき) と基準の表を付録 A に追加、推奨設定のファイルを版 2 に (plan [tooling-generator-defaults-and-launch-gate](../plans/active/tooling-generator-defaults-and-launch-gate.md) §5.1 #16)。
- `2026-10-14` — 値の正本を機械可読の `solver_density_cuda/tools/recipes/recommended_settings.yaml` (機種と段ごとのプロファイル) に移し、§1 を「機種に依らない共通部分と機種の決まっていない手書きの run の既定」と位置づけた (冒頭・§0.1)。機種ごとの値の表を付録 A (`recipe_settings.py --doc` の生成区間) に置いた。古い記述 3 件を訂正: §4「`runner.py` が生成する config が正本」、§2「生成器は `evaluate.wall_treatment_sst` を読む」(SERN だけ)、§3 の化学種の例 (生成器の lump 記法に)。plan [tooling-generator-defaults-and-launch-gate](../plans/active/tooling-generator-defaults-and-launch-gate.md) §5.1 #5・#6。
- `2026-10-07` — §4 に軸対称ノズルの物理壁 (NS) の現行 (`geometry.pw_upstream: poly`、寸法は `solve_rt` / `solve_rt_throat`) を追加、§9 に `geometry.pw_ramp`・CFD 前の旧 `solve_rt` を旧設定として (plan [tooling-nozzle-upstream-poly-and-throat-sizing](../plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md))。

- `2026-10-01` — §3 に LJ の出典 `physProp.ljSource` (既定 `[gri30, svehla1962]`、旧 run の再現は `[legacy_v1]`) を追加、§9 に旧既定 (plan [thermophysics-solver-owned-species-db](../plans/active/thermophysics-solver-owned-species-db.md) §4.10 #14)。
- `2026-09-25` — §1.0a に `space.slauWallNormalChi` の適用規則 (3D 側壁接続のみ 1、新構成は診断可能性の検査 + 3 条件) を追加 (plan [convection-slau-wall-normal-chi-usage-rule](../plans/accepted/convection-slau-wall-normal-chi-usage-rule.md))。
- `2026-09-22` — §2.1 遷移モデル (γ–Re_θt, `turbulence.transition: lm2009`) のレシピを追加 (plan [turbulence-transition-lm2009](../plans/active/turbulence-transition-lm2009.md))。
- `2026-09-17` — §6 に dual-time の内部反復レシピを追加 (`cfl_pseudo` 12–20 + `nSubIterDualTime` 10–20 + 緩和なし; 擬似 CFL に安定限界が見つからず、必要な nSub は `cfl_pseudo` で決まる)。定常の `implicitRelax 0.7` は据え置き。投入前チェック `check_solver_config.py` を追加。

- `2026-09-08` — 初稿 (散在していた推奨値を集約。ユーザ要請「既定の解析設定を 1 か所に、最新/旧を明記」)。
