# node の既定を `slauWallNormalChi: 1` にする (三値 auto)

## メタ

- **area**: `convection`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/convection/theory.md`](../../methods/convection/theory.md) (SLAU の $\chi$、「既知の限界」「対策」「既定」節)
  - [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) §1.0a
- **related_plans**:
  - [`convection-slau-wall-normal-chi.md`](../accepted/convection-slau-wall-normal-chi.md) (実装と受入)
  - [`convection-slau-wall-normal-chi-usage-rule.md`](../accepted/convection-slau-wall-normal-chi-usage-rule.md) (適用規則。本 plan の決定で §4.1 を置き換える)
  - [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (R5o-chi: 3D 格子収束は委譲先)
  - [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) (設計チェーンの runner)
- **created**: `2026-09-26`
- **owner**: `CFD Dev`

## 1. 目的

**ユーザ決定 (2026-09-26)**: node では `nodeWallDirichlet: 1` と `slauWallNormalChi: 1` を既定にする。理由 (ユーザ): 傾向はほぼ変わらず
(2D 生産 m6_on で 5 列とも生産許容内、case/16・case/48 はノイズ同程度)、壁 CV の排出による発散を防げる。
**本 plan が検証する利益は「既知構成 (node 側壁接続) の壁 CV 排出を抑える」に限定する** (codex plan M7。ユーザ決定は運用方針を変えるが、検証範囲を広げる根拠ではない)。
`nodeWallDirichlet` の既定は既に 1 (`solverConfig.hpp:443`) なので、本 plan で変えるのは `slauWallNormalChi` の既定のみ。
2026-09-25 の `diagnostician` 判断「既定化しない」(usage-rule plan §4.1) をユーザ決定で置き換える。

## 2. スコープ

- **やる**: 三値 auto の導入と解決・起動エコー (§4.1–4.2)、`stage_manifest` の実効値化 (§4.3)、ツールと設計チェーンの追随 (§4.4)、docs (§4.5)、
  ビット不変の確認 (B0) と標準 node ケースの A/B (B1)。
- **やらない**: `SLAU_d` カーネルの変更 (面流束の式は前 plan のまま)。3D の格子収束・固定点 (sern-3d R5o-chi)。カウル衝撃足への影響の評価
  (usage-rule plan で判定保留、衝撃衝突構成で flag 1 を使うときの注意として docs に残す)。SLAU のレジスタ上限 (既存の別問題)。

## 3. 関連 docs と前提

- 前 plan の受入 (試験した条件で回帰許容内) と usage-rule plan §6.2 (2D 生産の帯判定、診断可能性、衝撃足は判定保留)。
- 実装の制約: `solverConfig.cpp:680-693` は `slauWallNormalChi: 1` を cell・`nodeWallDirichlet ≠ 1`・solver が SLAU/SLAU2 以外で**起動エラー**にする。
  → 既定は構成に応じて解決する必要がある。
- `stage_manifest.py` は前 plan #12 で「`"0"` ならキーを落とす (省略 ≡ 0)」にしている。既定を変えると省略の意味が変わる。
- 未検証域: 周期・軸対称の生産規模 (V6 は小規模のみ)、凝縮/二相、衝撃衝突点。

## 4. 設計方針 (2026-09-26 `diagnostician`)

### 4.1 三値と解決

`space.slauWallNormalChi` は **省略 (= auto、内部 −1) / 0 / 1**。全キー読込後に `resolveSlauWallNormalChi()` で
**auto → 1 iff `discretization == "node"` ∧ `nodeWallDirichlet == 1` ∧ `solver ∈ {SLAU, SLAU2}`、それ以外 → 0**。
明示 1 の検証エラー (cpp:680-693) は不変。明示 0 は旧挙動。cell・非 SLAU で auto → 0 は**静かな解決 + エコー** (既定なので警告は出さない)。

### 4.2 起動エコー (常に 1 行)

`'slauWallNormalChi' effective: 1 (auto: node+nodeWallDirichlet+SLAU)` / `0 (auto: cell)` / `0 (auto: solver=ROE)` /
`0 (auto: nodeWallDirichlet=0)` / `1 (explicit)` / `0 (explicit)`。manifest・`RUN_PROVENANCE`・`diag_applicability.py` はこの行を正本にする。

### 4.3 実効値の来歴 (codex plan M1 で改訂)

- **forge 自身が起動ごとに `<run_dir>/forge_launches.jsonl` へ 1 行追記**する (`main.cpp` の `appendLaunchRecord`):
  `{time, cfg_fnv, bcond_fnv, exe_size, exe_mtime, slauWallNormalChi, slauWallNormalChi_source}`。`cfg_fnv` は `solverConfig.yaml` の FNV-1a 64。
  run_case.sh は `forge_run.log`・`RUN_PROVENANCE.txt` を起動ごとに上書きし、段階起動の runner は残差しか退避しないので、ログ抽出には頼らない。
  `run_case.sh` は起動エコーの行を `RUN_PROVENANCE.txt` にも写す。
- **`stage_manifest`**: 「`"0"` ならキーを落とす」を廃止し、段の key に **`space.slauWallNormalChi.effective`** (推定値) と `cfg_fnv`・`chi_source` を持たせる
  (`manifest_version: 2`)。`segments()` は `forge_launches.jsonl` を `cfg_fnv` で段に結び付けて**確定値で上書き**する。
  **推定 (`inferred`) の段と確定の段は、値が同じでも連結しない** (由来不明の推定値を既知と自動連結しない)。
- **旧形式の移行**: `space.slauWallNormalChi: "1"` はその値、キー無しは当時の既定 0 (`legacy`)。
- `diag_applicability.py` は新エコー (`'slauWallNormalChi' effective: N (...)`) を読み、省略かつエコー無しは実効値未確定で「診断不能」。

### 4.4 ツールと設計チェーン

- `check_solver_config.py`: node + SLAU で省略なら INFO「既定が 1 (2026-09-26)。旧結果の再現は 0 を明記」。
- `migrate_solver_config.py`: **書き換えない** (黙って 0 を足すと旧/新の区別が消える)。
- **runner は `slauWallNormalChi` を書かない (auto)** (codex plan M2: 明示 1 を書くと品質検査用の cell 変換 `convertGmshToForge` が起動エラーになる)。
  問題 YAML `mesh.slau_wall_normal_chi: 0` のときだけ明示 0 を書く。runner 経路でも auto の除外 (§4.6) が効く。
- **設計 DB (codex plan M3)**: `metrics.json` と ledger 行に `slau_wall_normal_chi_effective` (起動記録から) と `flag_policy: "2026-09-26"` を書く。
  driver は `flag_policy` が campaign と一致する行だけ学習に使う。旧行 (フィールド無し) は commit hash で「旧既定 0」と分類できるものだけ旧扱い、
  分類できない行は**学習から除外**。旧台帳からの再開で除外件数をログに出す。

### 4.5 docs

- `recommended-settings.md` §1.0a: 「既定 1 (2026-09-26 ユーザ決定)。旧挙動の再現は 0 を明記。診断 3 条件は『0 に落とす/落とさない』の判断材料として残す」。
  §9 旧設定表に「省略 = 0 (〜2026-09-25)」を追加。
- `methods/convection/theory.md` の「対策」「既定」節 (本 plan 起票時に更新済み、実装完了で「実装完了までは既定 0」の注記を外す)。
- `solverConfig.hpp:515` のコメント。前 plan 2 本 (accepted) に「ユーザ決定で既定化 (本 plan)」の 1 行。

### 4.6 未検証域の扱い (codex plan M4・M5、ユーザ指示 2026-09-26 で改訂)

- ユーザ指示: 周期・軸対称・凝縮は**除外でなく、flag が発火するケースで検証して auto に含めたい**。候補は case/39 周期丘 (並進周期 + 滑りなし壁、SLAU RANS)、
  case/40 軸対称ノズル (SST + 等温壁)、case/16 wys 凝縮 (SST + 滑りなし + 凝縮)。§6 B1-p / B1-a / B1-c。
- **回転周期は対象外**: node の周期処理 (`periodicNode_d`) に速度ベクトルの回転が無い (cell の `boundaryCond_d.cu:1728` のみ)。回転周期の 90° セクタは別作業。
- **成立条件**: 各ケースで (a) `wall_flag` 面 > 0、(b) 初回 `FORGE_DUMP_MASSFLUX` で省略 vs 明示 0 の差がある面 ≥ 1。成立しない域は「検証済み」に数えず auto から除外する
  (発火しない試験は根拠にならない)。
- **除外の実装**: 不合格・判定不能・成立しない域は、auto の解決規則に条件を足して 0 にする (周期 bcond の存在 / `isAxisymmetric` / 凝縮モデル有効)。明示 1 は使える。
  **通常域 (接続模型・case/46・case/16) の不合格は「除外」でなく既定化の保留**。

## 5. 実装ステップ

1. codex plan 段。
2. `solverConfig.{hpp,cpp}` の三値化・解決・エコー。
3. `stage_manifest.py` の実効値化と試験の書き換え。
4. `check_solver_config.py` の INFO、runner の明示キーと問題 YAML キー。
5. B0 (ローカル、1 step) → B1 (case/05 ローカル、他は AWS)。
6. docs、codex result、accepted。

### 5.1 残作業 (優先順)

**計算資源**: B0 と case/05 はローカル (1 step・小規模)。B1 の case/36・case/44・case/46 は AWS (ユーザ指示: ローカルで大きな計算をかけない。
AWS は他セッションと共有なので起動前に `aws_instance.sh status`/`busy`、手動 stop はしない)。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| ~~1~~ (**済 2026-09-26**: GO-with-changes C0/M7/m1、全件採用) | codex plan 段 | | F |
| ~~2~~ (**済 2026-09-26**) | 三値化 + 解決 + エコー | `solverConfig.{hpp,cpp}`、kernel へは `cfg > 0` を渡す防御。確認: エコー 5 パターン (省略→1 auto / 明示 0 / 明示 1 / ROE→0 auto / nodeWallDirichlet 0→0 auto) と明示 1 のエラー 3 種 (ROE・nwd 0・値 2) | O |
| ~~3~~ (**済 2026-09-26**) | 来歴 (起動記録) + manifest + diag_applicability | `main.cpp` の `appendLaunchRecord` (FNV は Python と一致を確認)、`stage_manifest.py` v2 (`test_stage_manifest_wall_normal_chi.py` 15/15、`test_gate_bad_input.py` PASS)、`check_convergence.py` が起動記録を渡す、`run_case.sh` が RUN_PROVENANCE に写す、`diag_applicability.py` (試験 22/22) | O |
| ~~4~~ (**済 2026-09-26**: runner は元からキーを書かない (= auto)。`mesh.slau_wall_normal_chi: 0` で明示 0。2D で「品質検査 cell 変換 → node 変換 → config」が auto・明示 0 とも rc=0。`metrics.json` に実効値 (起動記録の最後) と `flag_policy`、driver は `flag_policy` 不一致の PASS 行を学習から除外し件数を表示。3D runner は 2D と同じ cell 置換経路で、実 run は B1 (i) で確認) | runner auto + 変換経路、ledger/driver | `design/forge_design/evaluate/runner_sern*.py` (キーを書かない、問題 YAML で明示 0)、`metrics.json`/ledger/driver (§4.4)。合格: 2D・3D とも「品質検査 cell 変換 → node 変換 → 計算」が通り起動記録に auto→1、旧台帳からの再開で除外件数が出る | O |
| ~~5~~ (**済 2026-09-26**: 判定に影響しない INFO として表示) | `check_solver_config` INFO | node + SLAU で省略なら INFO | O |
| ~~6~~ (**済 2026-09-27: PASS**、AWS・現行 `2fa3826c` / 実装前 `22976398`: (a) 省略 vs 明示 1 の massflux 180,007 面ビット同一、(b) 明示 0 (現行、`scalarGradient: gg` 明記 — 2026-09-27 に node 既定が lsq になったため条件を揃えた) vs 実装前ビット同一、(c) cell・ROE・nodeWallDirichlet 0 の省略は 0 (auto)。`case/46.sern_design/CHI_DEFAULT_B0.txt`、`case/48.flat_plate_cooled_m4/chi_default_b0.py`) | **B0** | §6 B0。ローカル 1 step | O |
| ~~7~~ (**完了 2026-09-27: (i)(ii)(iii) とも PASS** → 既定化は保留しない。(ii) SERN 2D (`problem_moo_frozen_tp_cycle3op.yaml`、基準場 `run_0963/0964_chidef_2d_*_base` = runner 段階起動 12000 step・GATES PASS、分岐 `run_0965`–`0968` 各 36000 step・500 step 毎): m4_off・m10_on とも 4 run GATES PASS (実効 chi 0/1・scalarGradient lsq を metrics で確認)、力係数 5 列すべて帯内 (m4_off C_T −0.0001・C_L +0.0011・C_M −0.021 / m10_on −0.0003・+0.0012・−0.026)、省略側残差 ≤ 2×、m6_on は usage-rule plan の帯内結果を再利用 (`case/46.sern_design/CHI_DEFAULT_B1ii_*.txt`)。(iii) case/16 SST (`run_user_profile.py --disc node --phys sst` + 出口 outflow の基準 `run_0962_chidef_sst_base`、分岐 `run_0965/0966_chidef_sst_*` 各 12000 step): 違う面 627、輪郭壁 p/p0 L∞ 0.0003 %、3 点 STEADY、残差 ≤ 2× (`case/16.nozzle_wys/CHI_DEFAULT_B1iii.txt`)。(i) PASS 2026-09-27: AWS の場ファイルが削除済みのため起点を作り直した (`run_0960` = `jm_c16.msh` から準備 + run_0430 と同設定の層流暖機 1800 step・明示 0、監視 153797 の ρ 1.711e-4 で当時の記録と一致)。`run_0961` (run_0437 の設定から chi キーだけ削除、6000 step): 監視 3 点が step 14/11/8 で 10ρMin 超え・以後維持、壁の床到達 0、NaN 0、初回面流束 741 万面と流束状態が明示 1 (`run_0962`、1 step) とビット同一。`CHI_DEFAULT_B1i.txt`・`cad/chi_default_b1i.py`。(ii) m6_on は usage-rule の帯内結果を再利用、m4_off・m10_on は走行中。(iii) 未) | **B1 通常域** (i) 接続模型の救済 (AWS) (ii) case/46 2D 3 作動点 (AWS) (iii) case/16 SST (ローカル) | §6 B1 | O (結論 F) |
| 8 (**完了 2026-09-27: 3 域とも PASS → auto から除外する域なし**。B1-p case/39 (`run_0952/0953_chidef_b1p_*`、12000 step): 成立条件 = 初回面流束の違う面 3791、周期 group の壁/非壁混在 0、量 = 下壁 Cf(x) 相対 L2 0.0002 %・x_r 差 0・継ぎ目バルク速度比の差 0.0001 %、全量 STEADY、省略側残差 ≤ 2× (`case/39.periodic_hills/CHI_DEFAULT_B1p.txt`)。B1-a case/40 (`run_0954/0955_chidef_b1a_*`): 違う面 371、壁∩軸ノード 0、η_CF・ṁ の差区間 ±0.0001 % 以内、輪郭壁 p/p0 L∞ 0.38 %、③ q_w は断熱壁のため対象外、STEADY、残差 ≤ 2× (`CHI_DEFAULT_B1a.txt`)。B1-c case/16 凝縮 (`run_0963/0964_chidef_b1c_*`、起点 run_0335/res_48000 `aa22ab18…`): 違う面 434、壁 p/p0 L∞ 0.0012 %、onset 22.508/22.508 mm (差 0)、g_exit 0.01085、STEADY (onset の step 0 = restart 入力は g 未計算で非有限のため除外、既存の凝縮系列と同じ流儀)、残差 ≤ 2× (`CHI_DEFAULT_B1c.txt`)。2026-09-27 着手。**B1-p の基準場を `run_0007_coarse_rans` から `run_0039_r1_gradfix_new_ext/res_800000` に変更** (投入前): run_0007 の config は現在黙って無視される旧キー (トップレベル `lowMachPrecond` 等) を含み、run_0039 は同じ粗さ (80×50×30) で 80 万 step・`--from-floor` PASS・全量 STEADY、`slauWallNormalChi: 0` 明示なので §6 の「PASS または plateau + STEADY の基準場」を既に満たす。抽出器 `case/39.periodic_hills/chidef_b1p.py` は結果を見る前に commit。B1-a の基準場は `run_0045/res_12000` (§6 の `res_24000` は通算 step の書き誤り、sha256 `15e47af0…`)) | **B1 追加域** p: case/39 周期丘 / a: case/40 軸対称 / c: case/16 凝縮 | §6 B1-p/a/c。小規模は手元、重ければ AWS | O (結論 F) |
| ~~9~~ (**済 2026-09-27**: `procedures/recommended-settings.md` §1.0a を「既定 auto = 実効 1、旧挙動は 0 明記、診断 3 条件は 0 に落とす判断材料」に書き換え + §9 に「省略 = 0 (〜2026-09-25)」、`methods/convection/theory.md` の「実装完了までは既定 0」注記を外し検証完了を記載、`solverConfig.hpp` コメント、前 plan 2 本 (accepted) の変更ログに 1 行) | docs (§4.5) | 4 ファイル。`check_plans.py` PASS | O |
| 10 | codex result → accepted | | F |

## 6. 検証

**判定基準はすべて測る前に固定する。許容は既存のものを再利用し、新設しない。**

### B0 ビット不変 (面流束、1 step。codex plan m8 で縮小)

| # | 比較 | 合格 |
| --- | --- | --- |
| (a) | node + SLAU で省略 (auto → 1) vs 明示 1 | `massflux` ビット同一 |
| (b) | node + SLAU で明示 0 vs 実装前 commit のバイナリ | ビット同一 |
| (c) | cell・非 SLAU・`nodeWallDirichlet: 0` | **設定解決の単体確認のみ** (エコーが 0 (auto) を出す。流束比較はしない: cell の回帰対照は組まない運用) |

### B1 通常域 (不合格 = 既定化保留)

収束場 (無ければ STEADY まで延長した場) から**省略 (auto 1) と明示 0 を分岐**、同 step 数。プラトーは「非有限・上昇を拒否したうえで許す」(前 plan :278)。
各 run の収束: `check_convergence` を併記し、**省略側の `rms_*` 末尾平均 ≤ 2× 明示 0 側、RISING 列なし、非有限なし**。

| # | ケース | 量と抽出器 | 合格 (出典) |
| --- | --- | --- | --- |
| (i) | **接続模型の救済** (codex M7): `run_0437` と同じ起点・設定で**省略** 6000 step (AWS) | 前 plan V1-a/e (監視 3 CV が step 200 までに $10\rho_{Min}$ を超えて維持)、V1-d (床到達 0)、NaN 0、**初回面流束が明示 1 とビット同一** | 前 plan V1 の基準そのまま |
| (ii) | case/46 2D 3 作動点 (設計チェーン、AWS) | runner の GATES、5 列 (`C_T`, `C_T_with_shear`, `C_L`, `C_L_with_shear`, `C_M`) の差区間 (usage-rule plan §4.2)、準定常は usage-rule plan §4.3 の閾値 | GATES PASS、差区間が R5n 帯内 (0.002 / 0.002 / 0.05)。m6_on は usage-rule plan §6.2 を再利用 |
| (iii) | case/16 SST (`run_user_profile.py --disc node --phys sst`、出口 `outflow`) | 前 plan V3 の式: $p_0$ 59070、`extract_wall_pp0.py` (`MESH/COORD`)、輪郭壁 $x\in[10,94]$ mm の $L^\infty$、3 点系列 | $L^\infty \le 0.5$ %、3 点系列 `--drift 0.001 --osc 0.0025` |

### B1 追加域 (不合格・判定不能・成立しない = その域を auto から除外)

成立条件: (a) `wall_flag` 面 > 0、(b) 初回 `FORGE_DUMP_MASSFLUX` で省略 vs 明示 0 の差がある面 ≥ 1。各 VERDICT ファイルに成立条件の記録を含める。

| | B1-p case/39 周期丘 | B1-a case/40 軸対称ノズル | B1-c case/16 凝縮 |
| --- | --- | --- | --- |
| 基準場 | `run_0007_coarse_rans` (12 万セル、SLAU) を 30000 step まで延長し PASS または plateau + 下記量 STEADY を確認してから分岐 (生産格子・DDES は使わない) | `run_0045_node_yp1_outletfix_cont/res_24000` (ALL STEADY) | `run_0335` (node 2D SST 凝縮) の最終場。無ければ段階起動で本段 24000 step |
| 分岐 | `restart_field.py` → 省略 / 明示 0、各 12000 step、`--out-interval 500` | 同 | 同 |
| 追加の成立条件 | 周期 seam の両メンバに `wall_flag` が揃う (前 plan V6 P6-per の検査、混在 group 0) | 壁∩軸ノード数を記録 (前 plan V6 P6-ax) | — |
| 量 | ① 下壁 $C_f(x)$ (`res_wall_*.h5` の `twall_x`、$\tfrac12\rho_bU_b^2$ は両 run 共通に flag 0 の値) ② 再付着点 $x_r$ (下壁 $C_f$ の負→正ゼロ交差、`MESH/COORD`) ③ seam 列/隣接内部列のバルク速度比 | ① 推力効率 $\eta_{CF}$ と $\dot m$ ② 壁 $p/p_0$ (輪郭壁、$p_0$ = 入口全圧) ③ 等温壁 $q_w(x)$ (`qwall`) | ① 壁 $p/p_0$ (`extract_wall_pp0.py`) ② onset = `center_g` が初めて 1e-3 を超える x ③ $g_{exit}$ (報告のみ) |
| 許容 (出典) | ① 相対 L2 ≤ 1 % (case/48 V3 の $C_f$ ±1 %) ② ≤ 下壁の局所ノード間隔 1 つ (SERN の位置規則) ③ ≤ 0.5 % (case/16 V3 の値を流用 — **出典は流用と明記**) | ① 差区間 ≤ ±0.1 % (SERN V3 $C_T$) ② $L^\infty$ ≤ 0.5 % (case/16 V3) ③ 相対 L2 ≤ 1 % (case/48 V3 $q_w$) | ① 3 点 + $L^\infty$ ≤ 0.5 % (前 plan V3 case/16) ② ≤ 中心線ノード間隔 1 つ |
| 準定常 (各 run) | $C_f$ 3 位置・$x_r$・seam 比: `--drift 0.002 --osc 0.005` | $\eta_{CF}$・$\dot m$: 0.0002/0.0005、壁 3 点 $p/p_0$: 0.001/0.0025、$q_w$ 3 点: 0.002/0.005 | 3 点 $p/p_0$・onset: 0.001/0.0025 |
| 不合格時 | 周期 bcond の存在を auto→0 の条件に | `isAxisymmetric` を auto→0 の条件に | 凝縮モデル有効を auto→0 の条件に |

- 新規に書く抽出器は case/39 の $C_f$ / $x_r$ / seam 比だけ。**結果を見る前に commit** する。

### 不合格・判定不能の対応表 (測る前に固定)

| 事象 | 処置 |
| --- | --- |
| B0 不一致 | 実装バグ。B1 に進まない |
| B1 (i) 救済の不合格 | **既定化保留** (既定 0 のまま `in_progress`) |
| B1 (ii) 帯外 / GATES FAIL、(iii) 帯外 | 既定化保留 |
| DRIFTING / TRANSIENT | 1 回だけ本段を倍に延長。それでも同じなら、通常域は保留・追加域は除外 |
| 省略側に RISING 列・非有限 | 通常域は保留・追加域は除外 |
| 差区間の部分重なり (判定不能) | 通常域は保留・追加域は除外 |
| 追加域の成立条件を満たさない | その域を auto から除外 |

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-26` | [2026-09-26-convection-slau-wall-normal-chi-default-plan.md](../../notes/reviews/2026-09-26-convection-slau-wall-normal-chi-default-plan.md) | **GO-with-changes**, C0/M7/m1 | **全件採用** (2026-09-26 `diagnostician` 判断)。M1 → forge 自身の起動記録 `forge_launches.jsonl` + manifest v2 (推定と確定を連結しない、旧形式移行) (§4.3)。M2 → runner はキーを書かない (auto) (§4.4)。M3 → metrics/ledger に実効値と `flag_policy`、driver は不一致行を学習から除外 (§4.4)。M4/M5 → ユーザ指示で周期・軸対称・凝縮は**発火するケースで検証** (case/39・case/40・case/16 凝縮)、成立条件と除外処置を事前固定 (§4.6、§6 B1 追加域・対応表)。M6 → case/36・case/44・case/05 を落とし、抽出器・許容・準定常を登録済みのものに (§6)。M7 → 接続模型の救済を受入項目に、利益の記述を限定 (§1、§6 B1 (i))。m8 → cell 等は設定解決の単体確認に (§6 B0) |

### 6.2 結果

(未実施)

## 7. 影響範囲

- `solver_density_cuda/input/solverConfig.{hpp,cpp}` (設定の解決のみ。カーネルは不変)。
- **既存の node + SLAU の config で `slauWallNormalChi` を省略しているもの**は、新バイナリで挙動が変わる (flag 1 になる)。旧結果の再現には `0` を明記する。
- `stage_manifest.py`・`check_solver_config.py`・設計チェーンの runner。
- docs: `methods/convection/theory.md`、`procedures/recommended-settings.md`、前 plan 2 本。

## 8. 完了条件

- [ ] `methods/convection/theory.md` の既定節を実装完了に合わせて更新 (「実装完了までは既定 0」の注記を外す)
- [ ] 実装・検証完了 (§6 の B0・B1)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-27` — B0 (ビット不変) PASS、B1 通常域 (i) 接続模型の救済・(ii) SERN 2D 3 作動点・(iii) case/16 SST、追加域 B1-p 周期丘・B1-a 軸対称・B1-c 凝縮すべて PASS (auto から除外する域なし)。起点は AWS の場ファイル削除のため作り直した (接続模型) / 変更した (case/39 は run_0039)。docs 更新。

- `2026-09-26` — codex plan 段 GO-with-changes (C0/M7/m1) を全件採用して §4/§5.1/§6 を改訂。ユーザ指示で周期・軸対称・凝縮は発火するケースで検証する。#2 (三値化)・#3 (起動記録・manifest v2・diag_applicability) を実装。
- `2026-09-26` — 初稿。ユーザ決定 (node の既定を `slauWallNormalChi: 1` に) を受けて起票。§4/§6 は `diagnostician` の設計。
