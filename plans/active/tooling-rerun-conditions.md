# rerun_conditions: 形状を固定したまま入口条件・背圧・組成だけ変えて回す run を作る

## メタ

- **area**: `tooling`
- **status**: `in_progress`
- **related_docs**:
  - [`procedures/nozzle-design-workflow.md`](../../procedures/nozzle-design-workflow.md) §3a (本ツールが置き換える手作業の手順)
  - [`.claude/skills/nozzle-design/SKILL.md`](../../.claude/skills/nozzle-design/SKILL.md)
  - [`methods/design/overview.md`](../../methods/design/overview.md) (対応入力・初期場変換・拒否条件の節を追加する)
- **related_plans**:
  - [`tooling-nozzle-cfd-pinned-initial-line.md`](../accepted/tooling-nozzle-cfd-pinned-initial-line.md) (case/45 の生産 run — 本ツールの最初の利用先)
  - [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (将来の campaign 層。本ツールは run 単位の小道具で、それとは独立)
- **created**: `2026-10-06`
- **owner**: `Claude (主セッション) / ユーザ`

## 1. 目的

作成済みのノズル (既存 run の `nozzle.h5` に物理壁の格子が入っている) で、**形状を変えずに**入口条件 (Pt・Tt・H2O 分率・乱流量)・
背圧・輸送種の定義を保った入口分率 (H2O) だけを変えた計算を回したい (乾き成分 lump の組成変更は v1 で拒否) (ユーザ 2026-10-06「単純に入口条件変えて計算回したい」)。
problem YAML を書き換えて runner/deltastar_loop に通すと `design_chain` が壁を作り直してしまうので、run ディレクトリを直接複製・
編集する手順 (§3a) が要る。手作業だと複製漏れ (`resolved_species_*.yaml` を忘れて forge が起動しない前例)・化学種の互換判断の誤り・
変更記録の欠落が起きやすいので、run の準備 (複製・条件の書き換え・初期場の投入・記録) をツールにする。

## 2. スコープ

- **やる**: run の準備まで (複製 → bcond/solverConfig の書き換え → 初期場 → 記録)。forge の起動・判定はしない (既存の `run_staged_ns`/`run_staged` と判定ツールに任せる)。
- **やらない**: 形状の変更・再設計 (§1・§2 の設計作業)、格子の変更 (cross-mesh)、凝縮 ON/OFF の切り替え (種構成が変わる — `convert_species_field` の既存手順)、複数条件のスイープ管理 (campaign 層)。

## 3. 関連 docs と前提 (観測事実)

- run ディレクトリの入力: `nozzle.h5` (格子 + 初期場、`meshFileName`/`valueFileName`)、`solverConfig.yaml`、`bcondConfig.yaml` (1 行 1 境界の flow 形式、
  例 `inlet: {physID: 1, kind: inlet_Pressure, …, floats: {Y0: 0.9142, Y1: 0.0858, Pt: 5500000.0, Tt: 1600.0, k: 1.0, omega: 18000.0}}`、
  `outlet: {… kind: outlet_statPress, floats: {Ps: 2237.0, Pt: 2237.0, Tt: 300.0}}`)、`species_meta.yaml`・`resolved_species_*.yaml` (化学種の解決記録)、
  `probe.yaml`、`prepare_info.json` (`nozzle_report` が `scale_m`・`mesh.ni`・`Md`・`x_E`・`ic_from`・`dstar_source` を読む)、`wall_*.csv`、`MESH_QUALITY.txt`。
- `solverConfig.yaml` の `physProp.species` は `[{name: MIXDRY, lump: {CO2: …, O2: …, N2: …}, basis: mole}, H2O]`。Pt・Tt に依存する基準値は config に無い (case/45 run_0094 で確認)。
- 同一格子の初期場: `solver_density_cuda/tools/restart_field.py SRC_res.h5 DST_input.h5` は保存量を index コピーし、`forge --resolve-species` で宛先の化学種を解決して
  互換ハッシュが一致しなければ**書き込まずに停止**する (lump の成分比を変えると不一致)。種構成が違うときは `convert_species_field.py --mode conserve|reinit`。
- 一様初期場から始めない (AGENTS.md)。細分格子 NS は本段 CFL 5 で発散する (plan cfd-pinned §5.1 #11a) → 段階起動。

## 4. 設計方針 (diagnostician 2026-10-06 → codex plan 段 M1〜m8 を diagnostician が採否して確定)

1. **置き場所と CLI**: `solver_density_cuda/tools/rerun_conditions.py`。
   `REF_RUN NEW_RUN [--res res_N.h5] [--Pt] [--Tt] [--Y NAME=v | --Y1 v] [--balance NAME] [--k] [--omega] [--Ps P | --keep-Ps] [--Tw T | --keep-Tw] [--steps] [--out-interval] [--cfl] [--scale-ic none|pt] [--forge BIN] [--dry-run]`。
   forge は起動しない。`--lump` は受理するが停止 (4 参照)。`--Y` の NAME→index は `solverConfig.yaml` の `physProp.species` から。`--steps`/`--out-interval` は `nStepOuter % outStepInterval == 0` を強制、
   `--steps`・`--out-interval`・`--cfl` は YAML 上の位置 (`time.last.nStepOuter`・`time.outStepInterval`・`time.deltaT.cfl`/`cfl_pseudo`) の値トークンだけを書き換え (`yaml_strict.replace_scalars`; コメント・引用符付きキー・コロン前の空白・flow/block の両形式)、読み直して要求値との完全一致と他の値の不変を検査し、合わなければ作成前に停止 (正規表現置換はコメントに当たり実効値を変えずに「変更済み」と記録していた — codex result 段 4 回目 #3; runner_axismach の `_cfg_set` も同じ関数)。指定値と参照の実効値は有限・物理範囲 (Pt・Tt・Ps・Tw > 0、k ≥ 0、omega > 0、cfl > 0、steps・out-interval > 0) を作成前に検査し、入口 Y は変更の有無によらず各成分 [0,1]・ΣY=1 (変更時 1e−12、参照のまま 1e−9) (codex result 段 #1)。`--res` 既定は `res_[0-9]*.h5` の最大番号 (参照 config の nStepOuter と違えば警告)。SRC が float64 なら `restart_field.py --keep-src-dtype`。
   **Pt を変えたら `--Ps` か `--keep-Ps` が必須** (無ければ停止し、参照の出口圧 P_exit_ref [参照 res の出口断面 (最終 x の節点列) の内部節点の静圧 P の中央値。壁・軸の BC 節点は除く。節点列が取れなければ出口 BC の節点と同じ x の節点で代替し、それも取れなければ停止せず `P_exit_ref: null` と警告を記録。`res_outlet_*` の Ps は課した値そのものなので使わない (2026-10-06 改訂)] に対する `Ps/(f·P_exit_ref)` を表示)。
   node の出口は壁列が常に亜音速で Ps を見るため、Pt だけ下げて Ps 据え置きにすると出口列の不安定要因になる (`boundaryCond_d.cu:610-625`、run_0094 の出口壁側 3 節点 M 0.06/0.03/0)。
2. **入力契約 (v1 の対応入力、作成前に検査)**: 単一の `inlet_Pressure` (floats に `Y{s}` 形式、`inletProfile` 無し)、`outlet_statPress`、`wall`/`wall_isothermal`、`axis`、`meshFileName == valueFileName == "nozzle.h5"` (run 内相対)。
   `X{s}` 形式・`inletProfile: 1`・複数 inlet・外部参照・未知の入力依存は NEW_RUN を作る前に拒否。**ファイル参照は値の拡張子・存在でなく設定キーの完全修飾パスで判定** (`rerun_conditions.FILE_REF_KEYS`、出典は solverConfig.cpp の文字列キーとファイルを開く箇所の grep): `mesh.meshFileName`/`valueFileName` は `nozzle.h5` だけ、`physProp.speciesDBFile` は `species_db_external.yaml` の run 内相対 (実在) か空文字列だけ、`physProp.chemistry.mechanismFile`・`conjugate.solid` は拒否、表に無いがキー名 (`*File`/`*Path`/`*Dir`) か値の拡張子でファイル参照と見えるものも拒否。暗黙のファイル依存 (`conjugate:` ブロック → `conjugate_state_<physID>.h5`、bcond の `inletProfile`/`wallProfile` → 境界分布 CSV) は値によらずキーの存在で拒否 (拡張子で探していたので `physProp.speciesDBFile: /tmp/x` を受理していた — codex result 段 5 回目 #2)。
   **重複キーを含む YAML も拒否** (solverConfig・bcondConfig・species_meta、全階層の flow/block 形式; PyYAML は後勝ち・solver の yaml-cpp は先勝ちで、`turbulence: {model: "sst", model: "none"}` を PyYAML は none と読み k・ω をスケール対象から外していた — codex result 段 3 回目 #1)。読み込みは共通ローダー `solver_density_cuda/tools/yaml_strict.py` (runner_axismach の段 config 検査も同じ)。species_meta は変更の有無によらず検査する (書き換えは組成変更時だけ; codex result 段 4 回目 #4)。
   **merge key `<<` も全階層で拒否** (solver の yaml-cpp は merge を展開せず未定義キーは既定値を使う — `solverConfig.cpp:26`・`:840`。PyYAML の展開値で検査すると `turbulence: {<<: {sstEnergyIncludesK: 1}}` を 1 と読み solver は 0、merge 内の重複キーも検査を迂回した — codex result 段 4 回目 #1)。
   **使用禁止・廃止キーは値によらず作成前に拒否** し、参照設定から削除してから再実行するよう案内する (`mesh.bndFirstOrder` [AGENTS.md]、recommended-settings §9.1/§9.2 の削除キー、§9 の node 廃止キー・旧乱流キー体系; 一覧は `rerun_conditions.BANNED_KEYS`)。§9 の旧既定・非推奨の値 (`sstOmegaProdFromPk: 0`・`sstSigmaBlend: 0`・`sstEnergyKSource: 1`・`sstIsotropicStress: 1`) は拒否せず警告 (codex result 段 4 回目 #2)。
   **Euler (全壁 `slip`) と扱うのは輸送が無効と確認できる設定だけ**: 乱流・遷移なし、`viscMethod: 0`・`visc: 0`・(`thermCondMethod: 0` なら) `thermCond: 0`・`physProp.transport` なし (viscMethod 0 の μ は定数 visc で粘性流束は毎反復評価、種拡散は `viscMethod ≠ 0` のときだけ — `gasProperties_d.cu:91`・`speciesTransport_d.cu:960`)。全壁 slip でも輸送が有効な入力は未対応として作成前に拒否 (codex result 段 3 回目 #2)。粘着壁があれば NS。
   **複製は許可リスト**: `nozzle.h5`・`nozzle.xmf`・`bcondConfig.yaml`・`solverConfig.yaml`・`species_meta.yaml`・`resolved_species_*.yaml`・`species_db_external.yaml`・`probe.yaml`・`prepare_info.json`・`wall_*.csv`・`target_axis_M.csv`・`delta_r_initial.*`・`MESH_QUALITY.txt`。
   系譜・ログ・VERDICT・series・report は持ち込まない。NEW_RUN が既にあれば失敗、途中失敗は NEW_RUN ごと削除。
3. **書き換え**: bcond は対象行の floats だけを正規表現で置換 (他の行はバイト一致)、書き換え後に YAML で読み直して検証。Y は**全種を書き** Σ=1 を 1e−12 で検査
   (forge の起動検査は 1e−3 で入口カーネルが黙って正規化する; `speciesDB.cpp:992-1004`・`boundaryCond_d.cu:1064-1069`)。吸収種は 2 種なら自動、3 種以上は `--balance` 必須。
   `--Ps` は出口の Ps と Pt を個別に照合し、どちらかが違えば両方をその値に揃えて変更として記録する (`changes` の `Ps`・`outlet_Pt`; Tt 300 は据え置き。参照 `Ps: 2237, Pt: 3000` に `--Ps 2237` で Pt の同期を省き変更なし扱いしていた — codex result 段 5 回目 #4)。等温壁 + `--Tt` は `--Tw`/`--keep-Tw` が無ければ停止。`species_meta.yaml` の `streams.inflow.Y_transport`・`Y` を同期 (`X` は実種 MW が取れるときだけ)。
4. **初期場**: (1) config から**必要保存量集合**を決める = `forge_species.required_conserved` (ro/roU/roe/roY*/roXi) + SST → roK・roOmega + 遷移 → roGamma・roReth + 凝縮 → rog_*・roQ*_* + トレーサ → roXi。
   SRC/DST の存在・shape・有限・ρ>0・0 ≤ roY/ro ≤ 1+1e−6・|ΣroY − ro| ≤ 1e−6·ro を restart_field の前に検査し、DST `/VALUE` にこの集合と `wall_dist` 以外があれば拒否 (codex M1)。
   (2) `restart_field.py REF_res NEW/nozzle.h5 --dst-run NEW [--forge]` (VERDICT OK 行必須)。
   (3) **`--scale-ic pt` のときだけ** (明示 opt-in) 必要保存量を全部 × f (= Pt_new/Pt_ref)。これは **T・U・Y・k・ω を保つ初期場変換**であって境界値問題の相似ではない (Ps・壁温は別の BC)。
   検査 `allclose(d_new, f·d_ref, rtol=1e−6, atol=0)` でゼロはゼロのまま (比で検査しない — ゼロ成分で NaN)。変換後の場にも (1) と同じ場の検査 (有限・ρ>0・Y の範囲) をかける (codex result 段 #1)。凝縮 block・Tt 変更・組成変更が同時にあればスケール不可 (指定されても停止)。
   (4) lump の組成変更は v1 で停止: restart_field (互換ハッシュ不一致)・convert_species_field conserve (実種ごとの ρY 保存検査)・reinit (ξ の出所なし)・forge (旧ハッシュ属性) の全経路が拒否する。
5. **記録**: `NEW_RUN/RERUN_CONDITIONS.json` (参照 run/res、変更前後の全値、必要保存量集合、scale_ic と f・スケール検査の結果、P_exit_ref と `Ps/(f·P_exit_ref)`、restart_field の VERDICT 行、forge の sha256、ツールの commit、`recommended_stages`)。
   `prepare_info.json` は幾何を据え置き、`ic_from` を参照 res で上書き、`rerun_of` を追加。
6. **入口 k・ω は変えない** (指定時のみ上書き; 設計チェーン自身も固定値、`runner_wt.py:280-285`)。
7. **`recommended_stages`**: 変更なし → `none` (参照の `cfl`・`cfl_pseudo` をそれぞれ保持し、推奨 cfl との照合はしない。照合は条件を変えて推奨 cfl があるときだけ、生成 config の実効 `cfl_pseudo` [定常で効く値; procedures/solver-settings.md「CFL の定義」] と `cfl` に対して行う — codex result 段 5 回目 #1); 条件を 1 つでも変えた → `full` (細分格子は本段 cfl 1); **Pt を変えた → `full` かつ本段 cfl 1、`--scale-ic pt` を推奨 (scale none なら警告)** — §6 (ii′)・A3 の結果 (2026-10-06)。**確定 (2026-10-06)**: 同条件の B3 (scale none + full + 本段 cfl 1) は延長後も STEADY 不達で、入口配管の壁際に逆流域 (≈ 800 節点) が残り A3 と別の状態 (出口 M 5.9940 vs 5.9921) に向かった — scale none は規定時間内に準定常に達せず、入口配管の壁際の逆流域が残った (長い過渡か別の定常解かは未区別)。「スケール IC + none」を Pt のみ変更の推奨に昇格するかは §6 (ii) の結果で別途判断 (§5.1 #8)。**推奨と生成 config の整合** (codex result 段 #3): `run_staged_ns(stages="full")` の本段は生成 config の cfl・nStepOuter で回るので、NS 参照で推奨 (Pt 変更: 本段 cfl 1・nStepOuter ≥ 60000) と食い違えば必要な引数 (`--cfl 1.0 --steps 60000`) を示して作成前に停止 (`--override-recommended` で明示的に通し、記録に残す)。§4.7 は NS の実測に基づくので Euler 参照 (§4.8 の対; `run_staged` で回す) には整合検査をせず警告のみ (Euler の実績は §6 (ii) 腕 E: none・cfl 2・6000)。**NS で Pt と Tt・組成・壁温を同時に変える複合条件** (凝縮 run の Pt 変更を含む) は起動・整定が未検証で推奨対象外 (full・本段 cfl 1・60000 と警告を記録)。`--scale-ic pt` の禁止条件は §4.4(3) の Tt・組成・凝縮だけで、壁温は禁止条件でない (Pt + Tw は受理; 記録に「禁止」と書くのは Tt・組成・凝縮のときだけ — codex result 段 3 回目 #4)。
8. **Euler 参照**: Pt のみの変更でも、δ_E 評価には同条件の Euler rerun を対で作る (旧 Pt の Euler 参照だと edge 帯の |傾き| 判定が変わり δ が −4 % 動く合成反例; `deltastar.py:243-250`)。旧条件参照の `mdot_ratio_vs_euler` は診断量として記録のみ。
9. **前提作業: `run_staged_ns` の改修** (codex M4、design 側コード、ソルバ数値は触らない): StageManifest 記録 (`stage_manifest.py` 既存 API)・段ごとの `residual_history_<tag>.csv` 保持・前段の 1 次化 (convMethod 1/2 → 0)・段終了ゲート (最終 res の必要保存量が有限・ρ>0、非有限なら次段へ進まず停止)。段の config 変更は YAML 上の位置で値を読み値トークンだけを書き換えて読み直し、起動前に各段の実効値 (convMethod・cfl・cfl_pseudo・nStepOuter・outStepInterval) を照合して違えば例外 (正規表現置換は `convMethod:  2`・指数表記・block 形式で黙って外れた; codex result 段 #2)。`run_staged` (Euler) も同じ方式 (段 config の書き換えと照合に加え、各段 [soft・mid] の restart 前の段終了ゲート・段ごとの `residual_history_<tag>.csv`・StageManifest も; codex result 段 3 回目 #3)。

## 5. 実装ステップ

1. 諮問 (diagnostician、ユーザ指定) で §4 の要判断 4 点と §6 を確定。
2. codex plan 段レビュー。
3. 実装 + 単体試験 (`solver_density_cuda/tools/test_rerun_conditions.py`)。
4. 検証 run (§6) と文書同期 (procedures §3a・skill)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4 の要判断と §6 の確定 — 判断: 2026-10-06 diagnostician「lump・Tt スケールは v1 から外す」→ 同日 codex plan 段 M1〜m8 を diagnostician が採否し全件採用 (Pt スケールは opt-in、条件変更は full 既定、検証は粗格子でツールのゲート・run_0117 系列で物理受入れ) | 完了 | F |
| 2 | `run_staged_ns` 改修 (§4.9) | `design/forge_design/evaluate/runner_axismach.py`: StageManifest・段ごとの残差履歴・`convMethod: [12]`→0・段終了ゲート。合格: 段ゲートの単体 (NaN 注入で次段へ進まない)、既存の design/tests が通る | O |
| 3 | `rerun_conditions.py` 実装と単体 | §4.1〜4.8、§6 単体 (a)〜(m)。合格: `test_rerun_conditions.py` FAIL 0 | O |
| 4 | `methods/design/overview.md` 追記 | 対応入力・初期場変換 (T・U・Y・k・ω を保つ)・拒否条件 | O |
| 5 | 検証 run (i)(ii)(iv) | §6、AWS、粗格子 run_0094 系列。case README 台帳に run を追記 (主セッション) | O |
| 6 | 文書同期 — 完了 2026-10-06 (procedures §3a・skill・methods の rerun 節; 推奨は (ii″) で確定) | `procedures/nozzle-design-workflow.md` §3a・skill `nozzle-design` | O |
| 7 | codex result 段レビュー — 1 回目 NO-GO (C0/M4/m2) → 修正 (a2d5c658・3c613995)、2 回目 NO-GO (C0/M3/m2) → 修正 (Euler CLI の __main__ 位置・NS/Euler の実効 config 分類・推奨の適用範囲・増分の定義・文書同期)、3 回目 NO-GO (C0/M3/m2) → 修正 (#1 重複 YAML キーの拒否 `yaml_strict`・#2 Euler は輸送が無効な設定に限定・#3 `run_staged` の段終了ゲートと段記録・#4 Pt + Tw の文言統一; #5 台帳は主セッション)、4 回目 NO-GO (C0/M3/m2) → 修正 (#1 merge key `<<` の全階層拒否・#2 使用禁止・廃止キーの拒否 `BANNED_KEYS`・#3 --steps/--out-interval/--cfl の構造ベース書き換えと読み直し照合 `replace_scalars`・#4 species_meta の常時検査・#5 §4.7 最終形の文言)、5 回目へ | `codex_review.py --stage result` | O |
| 8 | ~~「スケール IC + none」を Pt のみ変更の推奨に昇格するか~~ 決着 2026-10-06 (diagnostician): 昇格しない ((ii) で両腕とも none・cfl 5 で発散) | 完了 | F |
| 10 | ~~スケール IC の残置/推奨/削除~~ **確定 2026-10-06 ((ii″) B3 → (b))**: Pt 変更時は `--scale-ic pt` を推奨・stages full・本段 cfl 1、scale none には警告: **Pt 変更時は `--scale-ic pt` を推奨・stages full・本段 cfl 1**、scale none + Pt 変更には警告。scale none の本段 cfl 1 も B3 (run_0133/0139、§6 (ii″)) で検証済み (延長後も STEADY 不達) | 完了 (ツール・単体 (n)・§4.7) | F |
| 9 | (本 plan 外) lump 変更対応・物理受入れ (iii) | `convert_species_field` に lump 同名・輸送 Y 保持のモード; (iii) は最初の生産利用で | F |
| 11 | Tt・組成を変えた rerun の整定長 ((iv″) で流量・δ_E が上限 4 ブロックでも漸近値を出さず、ユーザ決定「1」で生産利用へ持ち越し) | (iii) の最初の生産利用で、細分格子・本段 cfl 1・60000 step のとき量が STEADY になるまでの長さを実測し §4.7 の Tt/Y 行に書く 。**追加 (2026-10-06 監査)**: 同条件の Euler 参照 (Tt・組成の変更) を先に STEADY まで整定させ (run_0134 は 6000 step で DRIFTING)、その参照で δ_E を再抽出してから NS の整定を判定する | F |

## 6. 検証 (事前登録; diagnostician 2026-10-06 の改訂版)

**参照の使い分け**: ツールのゲート (再現性・入力契約・full 経路) = 粗格子 `case/45.isobutane_m6_d155/run_0094_ns_c2pin_pass2_ext6k/res_6000.h5` 系列 (壁解像 FAIL 29.4 % は既知として記録)。
物理の受入れ = 壁解像 PASS の `run_0117_ns_recal_final_ext/res_60000.h5` 系列 ((iii)、最初の生産利用で実施し本 plan の完了条件には含めない)。

**準定常の判定規約 (全 run 共通)**: 系列は res_1000..res_N (res_0 は含めない)。`check_quasisteady.py --series-csv <run>/quantities_series.csv --series-cols <列> --tail 0.4 --min-snaps 10` (n=12 → 末尾 5 枚)。
- δ_E(x_F) [`delta_E`]: `--drift 5e-5 --osc 1e-4` (末尾幅 ≤ 0.01 %)、直前窓差 |mean(8000–12000) − mean(3000–7000)| ≤ 0.01 %。
- 出口コア M [`exitM_A`]: `--drift 3e-6 --osc 3e-6` (末尾幅 ≤ 1.8e−5)、直前窓差 ≤ 2e−5。
- 流量 [`mdot`、2π∫ρU_x r dr の中央値]: `--drift 1e-5 --osc 2e-5`。
- STEADY 以外は不合格・延長 1 回 (6000 step) まで、なお不達は保留。check_convergence は NOT CONVERGED (plateau) を記録のみ、DIVERGED は FAIL。`--from-floor` は使わない。

**単体** (`solver_density_cuda/tools/test_rerun_conditions.py`、fixture は run_0094 の入力):
(a) 無変更 → bcond・solverConfig バイト一致 (`--steps/--out-interval` 指定時はその 2 トークンだけ差)、restart_field の VERDICT OK、必要保存量が参照 res と `np.array_equal`、`species_hash` 一致;
(b) `--Pt 4.4e6 --Ps 1789.6` → inlet/outlet 行の差分は当該トークンだけ・YAML 再読込で一致・他の行バイト一致; (c) `--lump` → 非ゼロ終了・NEW_RUN なし; (d) inlet が `inlet_uniformVelocity` → 拒否;
(e) `--Y H2O=0.09` → Y0 = 0.91 (Σ=1 を 1e−12)、3 種 fixture で `--balance` なし → 停止、`species_meta.yaml` の Y_transport [0.91, 0.09]; (f) 等温壁 + `--Tt` のみ → 停止; (g) 凝縮 block + `--scale-ic pt` → 停止;
(h) restart_field を失敗させる → NEW_RUN が残らない; (i′) 必要保存量が SRC に欠ける → 拒否; (j) 3 種 + トレーサ fixture で `--scale-ic pt` → roY2・roXi も f 倍、ゼロ成分はゼロ;
(k) `inletProfile: 1`・X{s} 形式・inlet 2 本・`valueFileName` が別名 → 作成前に拒否; (l) `--Pt` のみ (Ps 指定なし) → 停止し `Ps/(f·P_exit_ref)` を表示; (m) 無変更 → `recommended_stages: none`、`--Tt` → `full`。
`run_staged_ns` 改修の単体: 段ゲートが非有限の res で次段へ進まない (soft 段の res に NaN を注入)。

**検証 run** (AWS):
- (i) 無変更の再実行 (run_0094 → `--steps 12000 --out-interval 1000`、`stages="none"`、cfl 5): NaN 0; δ_E STEADY (規約) かつ末尾 5 枚平均が run_0104 の 0.725280 に |Δ| ≤ 0.01 %;
  出口コア M 末尾平均が 5.999265 に |Δ| ≤ 2e−5 (〜1e−4 は保留、> 1e−4 は FAIL)。Euler 参照は run_0086 (同条件)。
- (ii) Pt 0.8 倍 (4.4e6 Pa、`--Ps 1789.6` = 0.8 × 2237、cfl 5・12000 step・`stages="none"` を実験として明示):
  - 腕 E (Euler 対参照): run_0086 → `--Pt 4.4e6 --Ps 1789.6 --scale-ic pt`、`run_staged(stages="none", cfl 2)` 6000 step。ゲート: 出口コア M 末尾平均が run_0086 に |Δ| ≤ 2e−5、流量が 0.8 × run_0086 に相対 1e−4 (状態変換が Euler 解を保つことの直接検証)。E が不合格なら (ii) 全体を保留し再諮問。
  - 腕 A (`--scale-ic pt`) / 腕 B (`--scale-ic none`)、δ_E は腕 E を参照に抽出。ゲート: 両腕 NaN 0・δ_E と M が STEADY; |A−B| が δ_E ≤ 0.02 %・M ≤ 5e−5。
    壁解像は記録のみ (既知 FAIL; y1+>1 割合・最大値が run_0094 の 29.4 %/13.4 より下がることを記録)。診断量 (合否外): `mdot_ratio_vs_euler` (run_0086 基準、予想 ≈ 0.7995)、δ_E(A)/0.725280 (予想 1.02〜1.10)、M の向き (予想: 低下)。
  - 解釈 (先に固定): B 不安定・A STEADY → 「スケール IC + none」を Pt のみ変更の opt-in 推奨として §4.7 に追記 (既定は full のまま); A 発散 → スケール機能を削除候補に; 両腕 STEADY で一致 → 記録のみ; 不一致 → 延長 1 回 → 保留。
- (iv) `full` 経路の結合試験: run_0094 → `--Tt 1500 --Y H2O=0.10 --keep-Ps`、`run_staged_ns(stages="full")` (soft 3000・mid 3000・本段 cfl 5・6000)。
  ゲート: `stage_manifest.json` に 3 段、各段の `residual_history_<tag>.csv` が有限、段ゲートで停止なし、本段 NaN 0、実効入力 (bcond の Tt/Y1/Y0、`species_meta` の Y_transport) が要求値、本段の δ_E・M が STEADY (値の合否なし)。
  (Tt 1500・H2O 0.10 は経路試験用の任意値で報告に使わない。)
- (iii) 物理受入れ (最初の生産利用で実施、本 plan の完了条件外; 担当 F で条件決定): 参照 run_0117 res_60000、同条件 Euler rerun (run_0114 から) + NS rerun (`stages="full"`、本段 cfl 1・60000 step、Pt 変更なら `--Ps` 必須)。
  ゲート: NaN 0; δ_E・出口コア M・流量が規約で STEADY; 壁解像 PASS; 値は報告 (参照値が無いので合否なし)。
- **方針が誤りと言える観測**: (i) の M 差 > 1e−4 または δ_E 差 > 0.01 %; (iv) で段ゲートが非有限を通した、または実効入力が要求値と違う; (ii) 腕 E の不合格。

- **(ii′) スケール IC の残置判定 (diagnostician 2026-10-06、(ii) 両腕発散を受けた事前登録)**: 起点 run_0094 res_6000、`--Pt 4.4e6 --Ps 1789.6`、
  `run_staged_ns(stages="full")` (soft 3000・mid 3000・本段 cfl 5・12000 step・1000 ごと = §4.7 既定)。腕 A2 `--scale-ic pt` (`run_0129_rerun_pt08_scale_full`) / 腕 B2 `--scale-ic none` (`run_0130_rerun_pt08_noscale_full`)。δ_E は腕 E (run_0120) を参照に抽出。
  ゲート: 全段で段ゲート停止なし・NaN 0; 本段の δ_E・出口コア M・流量が規約で STEADY; |A2−B2| が δ_E ≤ 0.02 %・M ≤ 5e−5。
  診断量 (合否外): `mdot_ratio_vs_euler` (run_0120 基準)、δ_E(A2)/0.725280、本段で |δ_E − 末尾平均| ≤ 0.01 % に初めて入って以後出ないスナップショット番号。
  解釈 (先に固定): (a) A2・B2 とも STEADY で一致 → scale-ic pt は **opt-in 残置** (「T・U・Y・k・ω を保つ IC 変換であって起動補助ではない」と §4.4 に明記)、§4.7 不変;
  (b) A2 STEADY・B2 が破綻または延長 1 回でも STEADY 不達 → §4.7 に「Pt 変更時は `--scale-ic pt` を推奨 (stages は full のまま)」、ツールは Pt 変更 + scale none に警告;
  (c) A2 が段を問わず破綻 → **スケール機能を削除** (`--scale-ic` CLI・`scale_fields`・単体 (g)(j)・§4.4(3)・methods 節・記録の scale 項目を撤去し、腕 E の結果は §9 に残す)。B2 も破綻なら Pt 0.8 の BC/粗格子側の問題として再諮問;
  (d) 両腕 STEADY だが不一致 → 延長 1 回 (6000) → なお不一致は保留・再諮問。
  補足 (今は回さない): (c) のうち soft/mid は通り本段 cfl 5 だけが出口角で破綻した場合は run_0098/0101 と同指紋なので、削除の前に本段 cfl 1・60000 step の A3 を 1 本だけ登録して回す (A3 STEADY → §4.7 を「Pt 変更の本段は cfl 1」とし (a)/(b) を A3 に適用; A3 も破綻 → 削除)。
- **(iv) の決着規約 (延長 run_0128 の結果を見る前に固定)**: 延長後 STEADY → 合格。不達のとき: `classify` の detail が単調・増分減衰 (漸近値あり) で NaN 0・残差 rising でなければ「経路合格・量は漸近中 (漸近値 X・最終比 ±Y %)」と記録して plan 完了の妨げにしない (経路試験で値は報告に使わない; 量の定常性は (iii) で確認)。OSCILLATING・残差 rising・漸近値なしの DRIFTING → 保留・再諮問。これ以上延長しない。
- **方針が誤りと言える観測 (追記)**: (ii′) で A2 が soft/mid 段で破綻 (スケール場が NS で不整合) または B2 のみ健全。

- **(ii″) B3 (diagnostician 2026-10-06、(ii′)(b) の同条件確認; 先に固定)**: run_0094 res_6000 → `--Pt 4.4e6 --Ps 1789.6 --scale-ic none`、`run_staged_ns(stages="full")`、本段 cfl 1・60000 step・5000 ごと (A3 と同一)、延長 1 回 (6000) まで (`run_0133_rerun_pt08_noscale_full_cfl1`)。δ_E は run_0120 参照。
  ゲート: 段ゲート停止なし・NaN 0; δ_E・M・ṁ が規約で STEADY; |A3−B3| が δ_E ≤ 0.02 %・M ≤ 5e−5。
  解釈: (a) B3 STEADY・一致 → scale-ic pt は opt-in 残置、§4.7「Pt 変更: full + 本段 cfl 1 (scale は任意、どちらも検証済み)」、scale none の警告は情報表示に格下げ; (b) B3 が破綻 (段を問わず) または延長後も STEADY 不達 → (b) 確定、警告維持; (d) 両方 STEADY だが不一致 → 保留・再諮問。B3 の結果が出るまで §4.7 の (b) は「暫定 (根拠 A3 のみ)」。
- **(iv″) 整定長の確定 (diagnostician 2026-10-06、(iv) 保留の決着手順; 先に固定)**: run_0128 から restart_field (ビット一致) で cfl 5・6000 step・500 ごとのブロックを最大 4 回 (本段累計 36000) (`run_0135`〜`run_0138_rerun_fullpath_blk{1..4}`)。Euler 参照 = run_0086 → `--Tt 1500 --Y H2O=0.10 --keep-Ps`、`run_staged(stages="none", cfl 2)` 6000 step (`run_0134_rerun_euler_tt1500`、δ_E 抽出用、値の合否なし)。
  各ブロック終了ごとに規約 (M drift 3e−6/osc 3e−6、ṁ drift 1e−5/osc 2e−5、δ_E drift 5e−5/osc 1e−4、`--tail 0.4 --min-snaps 10`) で判定し、3 量 STEADY で停止。
  診断量 (合否外): ṁ の到達予想 ≈ 17150 (Pt/√(R·Tt) スケーリング、BL 変化含まず; 1 % 以上外れて単調に進み続ければ実効入力を疑う)、本段 res の入口列の T・Y1 実効値 (1500 K・0.10)。
  解釈: STEADY 到達 → (iv) 合格、§4.7 に「Tt/Y 変更: full + 本段 cfl 参照、量が STEADY になるまで (粗格子で N step)」を実測で記載; 上限で単調・増分減衰 (classify が漸近値) → 「経路合格・量は漸近中」で記録し §5.1 に F 項目を残して plan は進める; 上限で線形 (漸近値なし) または OSCILLATING → 保留のまま、ツール外としてユーザ判断に上げる (plan は done にしない)。注 (codex result 段 #5): `check_quasisteady` は DRIFTING のとき漸近値を計算しない (漸近値の外挿は STEADY の経路でのみ行う) — 判定と外挿診断は別であり、「classify が漸近値を出さない」ことは (iv)・(iv″) の線形ドリフトの証拠にならない (単調・増分減衰かは増分の時系列で見る)。
- **§4.7 最終形 (2026-10-06、実測で確定)**: (凝縮 ON の run で Pt を変える場合は scale-ic pt が使えない [§4.4] ので「Pt と凝縮 run の複合条件」として full・本段 cfl 1・60000 step・未検証の警告とする。) 無変更 → none・参照 cfl (run_0119: δ_E +0.0003 %・出口 M −4e−6 で参照を再現); **Euler 参照 (条件を変えた対)**: Pt のみ → `run_staged(stages="none")`・cfl 2・6000 step (腕 E run_0120 STEADY); **Tt・組成の変更 → 推奨は未確立** (run_0134 は同じ設定で 6000 step では DRIFTING — 出口 M 6.022・直前窓差 0.021、流量の直前窓差 212; ツールは stages full を記録し警告、STEADY まで延長が必要) [2026-10-06 監査で訂正: 先に「run_0134 も STEADY」と書いたのは誤り]; Pt 変更 → full + 本段 cfl 1・60000 + `--scale-ic pt` (A3 run_0131/0132 STEADY; scale none は B3 run_0133/0139 で規定時間内に準定常に達せず入口壁際の逆流域が残った、cfl 5 の本段は A2/B2 とも発散); Tt/Y 変更 → full + 本段 参照 cfl (粗格子 cfl 5 で本段 36000 step まで延長しても流量・δ_E は DRIFTING、簡易予想値 ≈ 17150 との差は約 0.1 %、整定余量は未確定 (§5.1 #11)、出口 M は 30000 step 付近で STEADY — 整定長は (iii) の生産利用で確定、§5.1 #11)。
- **方針が誤りと言える観測 (追記)**: (iv″) で ṁ が到達予想から 1 % 超外れて単調に進む (実効入力の不整合); B3 と A3 が両方 STEADY で不一致 (IC 依存)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-06` | diagnostician (Fable、ユーザ指定) — ブリーフ [`notes/reviews/briefs/2026-10-06-rerun-conditions-plan.md`](../../notes/reviews/briefs/2026-10-06-rerun-conditions-plan.md) | C2/M6/m6 | 全件採用 → §4・§6 確定。C: lump の自動切替却下 (v1 停止)・Tt スケール却下。M: Pt スケール既定 ON (凝縮 run 禁止)・Y 全種書き Σ=1 1e−12・検証参照を run_0094 系列へ・複製許可リスト・記録の拡充・k/ω 自動調整しない |
| plan | `2026-10-06` | [`notes/reviews/2026-10-06-tooling-rerun-conditions-plan.md`](../../notes/reviews/2026-10-06-tooling-rerun-conditions-plan.md) | GO-with-changes, C0/M7/m1 | 採否は diagnostician (同日 2 回目) が判断し全件採用: M1 必要保存量集合とゼロ対応の検査、M2 Pt スケールは opt-in・条件変更は full 既定・Pt 変更時は Ps 指定必須 (前回の既定 ON を撤回)、M3 入力契約、M4 run_staged_ns 改修を前提作業に、M5 同条件 Euler 参照、M6 準定常の閾値を量ごとに事前登録、M7 物理受入れは壁解像 PASS の run_0117 系列に分離、m8 methods と台帳 → §4・§5.1・§6 |
| result | `2026-10-06` | [`notes/reviews/2026-10-06-tooling-rerun-conditions-result.md`](../../notes/reviews/2026-10-06-tooling-rerun-conditions-result.md) | NO-GO, C0/M4/m2 | 全件採用: M1 入力範囲・参照 ΣY・スケール後検査、M2 段 config を YAML 構造ベース、M3 推奨と生成 config の食い違いで停止 (NS)、M4 監査可能化 (原データ取得・QS_VERDICT・rerun_audit) と B3 の結論の限定、m5 DRIFTING で漸近値は計算されない旨、m6 文書同期; 監査で run_0134 (Euler Tt/Y) が DRIFTING と判明し §4.7 の Euler 行を訂正 |
| result | `2026-10-06` | [`notes/reviews/2026-10-06-tooling-rerun-conditions-result-2.md`](../../notes/reviews/2026-10-06-tooling-rerun-conditions-result-2.md) | NO-GO, C0/M3/m2 | 全件採用: M1 `runner_axismach.py` の `__main__` を全定義の後へ (Euler CLI が `_first_order` で NameError) + 回帰試験 (design/tests/run_mesh_params_tests.py (g))、M2 NS/Euler を実効 config (粘着壁・乱流) で分類し prepare_info は照合のみ、M3 Euler の none は検証条件 (scale-ic pt・Ps も同じ比) に限定・NS の Pt と Tt/組成の複合変更では scale を推奨しない、m4 増分の定義を訂正、m5 §5.1 #7・#11 と methods を同期 |

## 7. 影響範囲

新規ツール 1 本と試験。既存ツールは `design/forge_design/evaluate/runner_axismach.py` の `run_staged_ns`・`run_staged` を変えた (段の記録・段終了ゲート・前段の 1 次化 convMethod 1/2 → 0・段 config の構造ベース変更と起動前の実効値照合; 段の CFL・step 数・nStepInner は不変)。他の runner (`runner_wt`・`runner_sern`・`runner_walldriven`) の同種の正規表現置換は変えていない。ソルバは変えない。
- 2026-10-06 (result 段 4 回目の対応): 新規 `solver_density_cuda/tools/yaml_strict.py` (重複キー・merge key の拒否、値ノード位置の書き換え `replace_scalars`)。`runner_axismach._cfg_set` を `replace_scalars` に共通化。同じ型の文字列置換 (`cfl: 4.0, cfl_pseudo: 4.0` の literal、`convMethod: 1` だけの置換) は `runner_wt.py:406-407`・`runner_walldriven.py:120-121`・`runner_sern.py:914,930,959` に残る — 本 plan の範囲外 (rerun_conditions が生成する run は通らない) として残し、別途扱う。

## 8. 完了条件

- [ ] 諮問・codex レビュー 2 回を記録
- [ ] 単体試験と検証 run が §6 を満たす
- [ ] procedures §3a・skill を同期
- [ ] `status: done`、`plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-10-06` — 起票 (ユーザ依頼「rerun_conditions.py つくってほしい。diagnostician に諮ってね」)。§4 は草案、要判断 4 点を diagnostician に諮る。
- `2026-10-06` — diagnostician (ユーザ指定) の判断を全件採用し §4・§6 を確定 (§6.1)。lump 変更と Tt スケーリングは v1 対象外、Pt スケーリングは既定 ON (凝縮 run は禁止)、検証は粗格子 run_0094 系列。
- `2026-10-06` — codex plan 段レビュー (GO-with-changes, C0/M7/m1) の採否を diagnostician (同日 2 回目、ユーザ指定) に諮り全件採用。diagnostician は前回の「Pt スケール既定 ON + stages none」を撤回 (出口壁列が Ps を見るので Pt だけ下げると出口不安定要因、合成反例で旧 Pt の Euler 参照は δ を −4 % 動かす)。§4・§5.1・§6 を確定版に置換、status → in_progress。
- `2026-10-06` — §5.1 #2・#3・#4 を実装 (implementer)。#2 `runner_axismach.run_staged_ns`: 段ごとに `stage_manifest.json` (`S1_soft`/`S2_mid`/`R<i>_cfl<c>`/`main`) と `residual_history_<tag>.csv` を残す・前段の 1 次化を `convMethod: [12]\b`→0・段終了ゲート `stage_gate` (最終 res の必要保存量が有限・ρ>0 でなければ次段へ進まず RuntimeError、restart_field へ渡さない)。段の CFL・step 数は不変。#3 `solver_density_cuda/tools/rerun_conditions.py` (§4.1〜4.8) と `test_rerun_conditions.py`: 単体 (a)〜(m)・段ゲート・restart を模擬した作成経路で FAIL 0。ただし (a)(j) の**作成経路** (実 restart_field) は、ローカルに `--resolve-species` を持つ forge が無いため SKIP (理由つき; `--force-species` では通していない) — AWS など新しい forge (FORGE_BIN) で再実行が要る。#4 `methods/design/overview.md`「既存 run の条件変更 (rerun_conditions)」節。
- `2026-10-06` — §4.1 の P_exit_ref の定義を改訂 (implementer の指摘を主セッションが採用: `res_outlet_*` の Ps は課した値 2237 そのもので比較にならない)。参照 res の出口断面の内部節点 (壁・軸の BC 節点を除く) の P の中央値に変更し、`rerun_conditions.exit_pressure_ref` と単体試験を更新 (run_0094 res_6000 で P_exit_ref = 2241.998 Pa、95/97 節点; FAIL 0・SKIP 2)。
- `2026-10-06` — **単体 FAIL 0・SKIP 0 (AWS)**: ローカルでは --resolve-species を持つ forge が無く (a)(j) の作成経路を SKIP したが、AWS で `FORGE_BIN=~/forge-wallfit-bin/.../forge` を渡して全項目合格。**検証 run 投入** (`case/45.isobutane_m6_d155/run_rerun_validation.sh`、量の時系列 `rerun_series.py`): (i) `run_0119_rerun_ctrl`、(ii) `run_0120_rerun_euler_pt08` (E)・`run_0121_rerun_pt08_scale` (A)・`run_0122_rerun_pt08_noscale` (B)、(iv) `run_0123_rerun_fullpath`。§6 からの変更 (投入前): E と (iv) は 6000 step だと準定常判定の 10 枚に届かないので出力間隔を 500 に (判定規約 `--tail 0.4 --min-snaps 10` を守るため; step 数は不変)。
- `2026-10-06` — 検証バッチが Euler 対参照 (run_0086、壁 `kind: slip`) で「対応外の境界種別」により作成前に停止 (forge 未起動; run_0119 は作成済み)。§4.2 の入力契約に `slip` (Euler の滑り壁、壁温なし) を追加し、P_exit_ref の除外にも含めた。単体 (f′) 追加、FAIL 0。バッチは既存の run_0119 を作り直さずに再開。
- `2026-10-06` — 再開した検証バッチが Euler 対参照の作成で再停止: run_0086 の `nozzle.h5` に乱流モデルなしでも変換器が作る roK/roOmega の入れ物があり、§4.4 (1) の「必要保存量と wall_dist 以外は拒否」に当たった。乱流モデルなし (Euler・層流) の run に限り roK/roOmega を「使われない余りの量」として許し警告に記録 (スケールしない; 他の未知の量は従来どおり拒否)。単体 (f″) 追加、FAIL 0。
- `2026-10-06` — **検証 run の結果** (AWS、粗格子 run_0094 系列): **(i) 合格** — `case/45.isobutane_m6_d155/run_0119_rerun_ctrl` δ_E 末尾 5 枚 0.7252822 (run_0104 の 0.725280 に +0.0003 % ≤ 0.01 %)、出口コア M 5.999261 (5.999265 に −4e−6 ≤ 2e−5)、δ_E・M・流量 STEADY (§6 の閾値)。**(ii) 腕 E 合格** — `run_0120_rerun_euler_pt08` 出口コア M 5.999998 (run_0086 の 5.999998 に |Δ| < 1e−6)・流量 13344.97 = 0.8 × run_0086 の 16681.20 に相対 7.5e−7 (≤ 1e−4)、run_0120 自身は STEADY。**限定 (codex result-5 Major 3)**: 比較元 run_0086 は step 6000 の単一スナップショットしか残っていない (中間は AWS のディスク整理で削除) ので、run_0086 の準定常性は監査資料から判定できない — 主張は「基準スナップショットに対する 2 量の差が閾値内」に限る。**(ii) 腕 A・B は両方発散** — A `run_0121_rerun_pt08_scale` (scale あり) は step 334 で ro NaN、非有限 54 節点が x/r_t 94.7〜95.2・r/r_t 10.09〜10.10 (出口の壁際の角、T → 6000 K); B `run_0122_rerun_pt08_noscale` (scale なし) は step 120、x/r_t −12.5〜−12.0 (入口、T → 50 K の膨張)。いずれも §6 で「実験として明示」した stages none・cfl 5。登録解釈「A 発散 → スケール機能を削除候補に」に当たる → 諮問 (diagnostician)。**(iv)** `run_0123_rerun_fullpath` は経路のゲート (stage_manifest 3 段・段ごとの残差履歴・段ゲート停止なし・NaN 0・実効入力 Tt 1500/Y 0.9,0.1/Y_transport) 合格、本段 6000 step の出口コア M・流量は DRIFTING → 登録どおり 6000 step 延長 1 回 (`run_0128_rerun_fullpath_ext`)。集計スクリプトの欠陥: 発散 run (スナップショットなし) で `rerun_series.py` が例外 → `set -e` でバッチが残りの集計を飛ばした (手で再実行)。AWS は別セッションの 3D 変換の OOM で不通になり、ユーザ許可で停止・起動 (IP 100.54.10.165)。
- `2026-10-06` — diagnostician に諮った (エスカレーション条件 3): 「A 発散 → 削除候補」は発火したが、(ii) は両腕とも stages none・cfl 5 で発散し scale の良否を識別していない (B は入口 step 120、A は出口角 step 334; 出口 BC の ghost 構築はスケール不変 `boundaryCond_d.cu:614-624`、config に Pt 依存定数なし)。第 1 仮説 = 非固定点 IC からの本段 cfl 5 の安定限界超え (run_0098/0101 と同指紋、run_0031 は full なら Euler IC からでも cfl 5 完走) — 原因の確定は (ii′) の後。#8 は「昇格しない」で決着。残置/推奨/削除は (ii′) [scale-ic pt/none × full] の事前登録表で決める (§6)。(iv) の延長不達時の扱いを事前登録。
- `2026-10-06` — **(ii′) の結果**: 腕 A2 `case/45.isobutane_m6_d155/run_0129_rerun_pt08_scale_full` (scale-ic pt) は soft・mid 段を通過し**本段 (2 次 cfl 5) の step 468 で発散** (非有限 54 節点 x/r_t 94.7〜95.2・r/r_t 10.09〜10.10 = 出口の壁際の角、T 上限 6000 K; (ii) 腕 A と同じ位置); 腕 B2 `run_0130_rerun_pt08_noscale_full` (scale none) は soft・mid 段を通過し**本段の step 2 で発散** (x/r_t −12.5〜−12.0・r/r_t 5.7〜6.1 = 入口、T 50 K)。段記録 (stage_manifest 3 段・段ごとの残差履歴) は両腕とも正常。→ 登録 (c) の補足「soft/mid は通り本段 cfl 5 だけが出口角で破綻 → 削除の前に本段 cfl 1・60000 step の A3 を 1 本だけ」に該当 → `run_0131_rerun_pt08_scale_full_cfl1` (scale-ic pt、stages full、本段 cfl 1・60000 step・5000 ごと) を投入。**(iv) 延長** `run_0128_rerun_fullpath_ext`: 出口コア M 6.0274 → 6.0299 (末尾 5 枚、直前窓差 +0.0008)・流量も増加中、`classify` は DRIFTING (extremum at tail-end、漸近値なし) → 決着規約により**保留・再諮問** (延長はしない)。
- `2026-10-06` — **(ii′) の補足 A3 の結果と決着**: `case/45.isobutane_m6_d155/run_0131_rerun_pt08_scale_full_cfl1` (scale-ic pt、stages full、本段 cfl 1・60000) は完走・NaN 0、延長 1 回 `run_0132_rerun_pt08_scale_cfl1_ext` (6000 step・500 ごと、restart_field ビット一致) の末尾 5 枚で δ_E 0.7494484 (直前窓差 0.0005 %)・出口コア M 5.992124 (幅 3.5e−6)・流量 13341 (run_0120 比 0.99970) が §6 の閾値で全 STEADY (延長前は δ_E の直前窓差 0.0114 % が 0.01 % を僅かに超えたため登録どおり延長)。壁解像は粗格子の既知 FAIL (y1+>1 25.9 %、run_0094 の 29.4 % より下がった)。診断量: δ_E/0.725280 = 1.033 (予想 1.02〜1.10 内)、出口 M は低下 (予想どおり)。→ 登録 (c) の補足「A3 STEADY → §4.7 を『Pt 変更の本段は cfl 1』とし (a)/(b) を A3 に適用」: B2 は破綻 → **(b) Pt 変更時は `--scale-ic pt` を推奨 (stages full・本段 cfl 1)、scale none に警告**。ツール (`recommended_stages` に cfl 1、警告) と単体 (n) を追加、FAIL 0。注記: B 腕は cfl 5 でしか試しておらず「scale none + 本段 cfl 1」は未検証 — 推奨の根拠は A 腕が安定したことだけ。
- `2026-10-06` — diagnostician に諮った ((iv) 保留・再諮問): (iv) の DRIFTING は Tt/H2O 変更の正当で遅い整定 (NaN 0・plateau・ṁ は到達予想 ≈ 17150 まで残り ≈ 0.2 %) と判断、「cfl 1・60000」は累積擬似時間が同じで整定を進めないため却下 → cfl 5 のまま上限付き延長 (iv″) と同条件 Euler rerun を登録。(ii′)(b) は A3 (cfl 1) と B2 (cfl 5) の非同条件比較なので暫定とし、B3 (none + full + cfl 1) で確定 (ii″)。result 段レビューは (ii″)・(iv″)・#6 の後。
- `2026-10-06` — **(iv″) の結果とユーザ決定「1」**: Euler 参照 `case/45.isobutane_m6_d155/run_0134_rerun_euler_tt1500` (run_0086 → Tt 1500・H2O 0.10・keep-Ps、cfl 2・6000) 完走。ブロック延長 `run_0135`〜`run_0138_rerun_fullpath_blk1..4` (cfl 5・各 6000・restart_field ビット一致) の末尾 5 枚: 出口コア M 6.03169 → 6.03243 → 6.0326 → **6.03269 (blk4 で STEADY)**; 流量 17128 → 17132 → 17135 → 17137 (各 run 内の隣接 5 枚窓 [= 2500 step 離れた窓] の平均差 +2.9 → +1.4 → +0.95 → +0.82; 6000 step 離れたブロック末尾平均どうしの差は blk1→2 +4.26・blk2→3 +2.43・blk3→4 +2.05 — codex result-2 m4 の再計算); δ_E 0.73321 → 0.72946 → 0.72765 → 0.72651 (ブロック間の差 −0.00374 → −0.00182 → −0.00114)。「残り約 0.1 %」は漸近値でなく簡易予想値 ≈ 17150 (Pt/√(R·Tt) スケーリング) との差。流量・δ_E は上限 4 ブロックで DRIFTING (classify は漸近値を出さず) → 登録の「上限で線形 (漸近値なし) → 保留・ユーザ判断」に該当 → **ユーザ決定 2026-10-06「1」: ツールの検証としては「経路合格・量はゆっくり整定中 (流量は簡易予想値との差 約 0.1 %、ブロック間の増分は減衰しているが末尾 2 ブロックで鈍い; DRIFTING の記録から真の整定余量は確定できない)」と記録して先に進む。整定に要る長さは実際の生産利用 ((iii)、細分格子・cfl 1・60000 step) で確かめる** (ツールの欠陥ではなく条件変更後の流れの整定の遅さ = ソルバ側の性質; §5.1 に F 項目)。NaN 0。
- `2026-10-06` — **(ii″) B3 の結果 → (b) 確定**: `case/45.isobutane_m6_d155/run_0133_rerun_pt08_noscale_full_cfl1` (scale none、full、本段 cfl 1・60000) は完走・NaN 0 だが出口コア M・流量が DRIFTING (流量 13367 → 13721 → 13814 → 13433 → … → 13371 と大きく動いてから単調減少)、入口配管の壁際 (x/r_t −12.5〜−9.1、r/r_t 6.24〜6.48) に Ux < 0 の逆流域 854 節点 (A3 は 0)、δ_E 抽出は破綻 (1e9 — 逆流で縁帯が取れない)。延長 1 回 `run_0139_rerun_pt08_noscale_cfl1_ext` (6000 step) でも DRIFTING・逆流域 779 節点・出口 M 5.993977 (A3 5.992124 と 0.0019 違う別の状態)。→ 登録 (b)「B3 が延長後も STEADY 不達 → (b) 確定、警告維持」: **Pt 変更時は `--scale-ic pt` を推奨 (stages full・本段 cfl 1)、scale none に警告** を確定 (§4.7・§5.1 #10)。
- `2026-10-06` — codex result 段 (NO-GO) の指摘 1・2・3・5・6 を修正 (implementer): #1 指定値・参照値の有限性と物理範囲、入口 Y を常時検査 (参照 ΣY 1.0001 を拒否)、スケール後の場の検査; #2 `run_staged_ns`/`run_staged` の段 config を YAML 上の位置で書き換え起動前に実効値を照合 (空白 2・指数表記・block 形式の回帰); #3 NS 参照で推奨 (Pt 変更: 本段 cfl 1・≥ 60000) と生成 config が食い違えば停止 (`--override-recommended`)、Euler 参照は警告のみ・回し方を `run_staged` と表示; #5 §6 (iv″) に注; #6 methods・§5.1 #10・§7・plans/README。`test_rerun_conditions.py` FAIL 0・SKIP 2 (ローカル)、design/tests 4 本 FAIL 0。指摘 4 は未了。
- `2026-10-06` — codex result 段 (NO-GO, C0/M4/m2) の指摘 1・2・3・5・6 を実装担当が修正 (入力の範囲・参照 ΣY の常時検査・スケール後の field_problems; 段 config を YAML 構造ベースに・各段の実効値を起動前に照合 [run_staged_ns・run_staged]; 推奨と生成 config の食い違いで作成前に停止 [NS のみ、`--override-recommended`]; §6 (iv″) の注記; methods・§7・plans/README の同期)。主セッションの判断: 停止は NS 参照に限り、Euler 参照の推奨は実績 (腕 E・run_0134) から「stages none・cfl 2」として §4.7 に追加 (ツールの recommended_stages も同じ)。指摘 4 (監査可能な原データ) は主セッションが対応中。
- `2026-10-06` — **指摘 4 (監査可能化) と監査で見つかった訂正**: 検証 run 17 本の小さな成果物 (config・全段の残差履歴・stage_manifest・量の時系列 CSV・RERUN_CONDITIONS.json・判定ログ) と主要な最終場 (run_0119・0120・0132・0138・0139 の res + nozzle.h5) を AWS からローカル `/home/sano/work/forge/case/45.isobutane_m6_d155/` へ取得 (AWS は ユーザ許可で起動、IP 13.220.68.91)。`case/45.isobutane_m6_d155/rerun_audit.py` で各 run に `QS_VERDICT.txt` (§6 の判定規約のコマンドと出力) を書き、要約を `case/45.isobutane_m6_d155/_band_ab/rerun_audit.txt` に残した。**訂正**: Tt・H2O を変えた Euler 対照 `run_0134_rerun_euler_tt1500` は stages none・cfl 2・6000 step で **DRIFTING** (出口 M 6.0218・直前窓差 0.021、流量 16901・直前窓差 212) — 先に §4.7 に「run_0134 も STEADY」と書いたのは誤り (主セッションが確認せずに書いた) → §4.7 の Euler 行を「Pt のみ実績あり、Tt・組成は未確立」に直し、ツールも Tt/組成の Euler 参照には full と警告を記録するよう変更。**影響**: (iv″) の δ_E (run_0135〜0138) は定常でない Euler 場 (run_0134 の最終場) を参照に抽出していた — (iv″) の δ_E の整定判断は参照が動いていた分だけ割り引く (出口 M・流量は Euler に依存しない)。B3 の結論は「規定時間 (60000 + 6000 step) 内に準定常に達せず、入口配管壁際の逆流域 (779 節点) が残った」に限定し、「偽のはく離」とは断定しない (長い過渡と別の定常解は区別できていない; codex result Major 4)。
- `2026-10-06` — **codex result 段 2 回目 (NO-GO, C0/M3/m2) を反映**: Euler CLI (`runner_axismach.py`) の `__main__` を全関数定義の後へ (`run_staged` の `_first_order` が CLI で NameError だった; 回帰試験 run_mesh_params_tests (g) は runner_axismach と deltastar_loop の両方を検査)。rerun_conditions の NS/Euler 分類を実効 config に (粘着壁 wall/wall_isothermal → NS、全 slip → Euler、sst + 全 slip・粘着と slip の混在・prepare_info.viscous との食い違いは停止)。推奨の適用範囲: Euler の stages none は「Pt のみ・scale-ic pt・Ps = f·Ps_ref」に限定 (それ以外は full + 未検証の警告)、NS で Pt と Tt/組成/壁温/凝縮を同時に変えたときは scale-ic pt を推奨せず複合変更の未検証を警告。§6 (iv″) の増分の定義を訂正 (隣接 5 枚窓の差とブロック間の差を区別、「残り 0.1 %」は簡易予想値との差)。§5.1 #7・#11 更新。単体 (q) 5 項目追加、全試験 FAIL 0。
- `2026-10-06` — **codex result 段 3 回目 (NO-GO, C0/M3/m2) の指摘 1〜4 を修正** (implementer): #1 重複 YAML キーを全階層で拒否する共通ローダー `solver_density_cuda/tools/yaml_strict.py` を rerun_conditions (solverConfig・bcondConfig・species_meta) と runner_axismach の段 config 検査 (`_cfg_load`・`stage_gate`) に適用; #2 全壁 slip を Euler と扱うのは乱流・遷移なし・viscMethod 0・visc 0・thermCond 0 (thermCondMethod 0)・transport なしに限定し、輸送が有効なら作成前に拒否 (run_0086 は Euler 判定); #3 `run_staged` (Euler) の soft・mid 段に restart 前の段終了ゲートと `residual_history_<tag>.csv`・`stage_manifest.json` (段の CFL・step 数は不変); #4 scale-ic pt の禁止条件から Tw を外し、Pt との複合条件は「起動・整定は未検証で推奨対象外」に統一 (methods 同期)。単体 (r)(s)(t) と run_staged の ro=−1/Inf ゲート試験を追加: `test_rerun_conditions.py` FAIL 0・SKIP 2 (ローカル)、design/tests 4 本 FAIL 0。指摘 5 (case README) は主セッション。
- `2026-10-06` — **codex result 段 4 回目 (NO-GO, C0/M3/m2) の指摘 1〜5 を修正** (implementer): #1 `yaml_strict` が merge key `<<` を全階層で拒否 (`MergeKeyError`); #2 参照 config の使用禁止・廃止キー (`mesh.bndFirstOrder` は値によらず・§9 の削除キー等) を作成前に拒否、旧既定・非推奨値は警告; #3 `--steps`/`--out-interval`/`--cfl` を `yaml_strict.replace_scalars` (compose した木の値ノードの文字範囲だけ置換・読み直して要求値と他の値の不変を検査) に置換、runner_axismach `_cfg_set` も共通化; #4 species_meta を常時検査; #5 §4.7 最終形の Tt/Y 行を訂正。単体 (u#1〜#4) を追加: `test_rerun_conditions.py` FAIL 0・SKIP 2 (ローカル)、design/tests 4 本 FAIL 0。
- `2026-10-06` — **codex result 段 5 回目 (NO-GO, C0/M3/m3) の指摘 1・2・4 を修正** (implementer): #1 無変更 rerun は参照の cfl・cfl_pseudo をそれぞれ保持し推奨と照合しない (参照 `cfl 4.0, cfl_pseudo 5.0` を止めて `--cfl 4.0` を提案していた)、照合は条件変更で推奨 cfl があるときだけ生成 config の実効 cfl_pseudo と cfl に; #2 ファイル参照を設定キーの完全修飾パスで検査 (`FILE_REF_KEYS`; `speciesDBFile`・`chemistry.mechanismFile`・`conjugate.solid`、暗黙依存の `conjugate:`・bcond `wallProfile` も拒否、`speciesDBFile` は `species_db_external.yaml` だけ [以前は `probe.yaml` 等の許可リスト名も通した]); #4 `--Ps` で出口 Ps・Pt を個別照合し `outlet_Pt` を変更記録。単体 (q) を追加し `test_rerun_conditions.py` FAIL 0 (SKIP 2 はローカル forge 旧版)、design の 4 試験 FAIL 0。指摘 3・5・6 は主セッション。
- `2026-10-06` — **codex result 段 5 回目 (NO-GO, C0/M3/m3) を反映**: M1 (前回修正の退行: 無変更 rerun で `cfl ≠ cfl_pseudo` を拒否) — 無変更時は両値を保持し照合しない、推奨 cfl の照合は条件変更時の実効 cfl_pseudo/cfl に; M2 外部ファイル参照を設定キーの完全修飾パスで検査 (speciesDBFile は species_db_external.yaml のみ、mechanismFile・conjugate・inletProfile/wallProfile は拒否; 実装担当の自己点検で speciesDBFile が任意の実在ファイルを通していた件・wallProfile の素通り・conjugate の暗黙依存も修正; 全 run の solverConfig で誤検出 0 を確認); m4 `--Ps` は出口 Ps・Pt を個別に照合; M3 腕 E の主張を「基準スナップショットに対する差」に限定 (run_0086 は単一スナップショット); m5 methods・procedures の旧説明を修正; m6 発散の記録を訂正: **残差 CSV で最初に非有限になる step は run_0121 333・run_0122 119・run_0129 467・run_0130 1 (いずれも `rms_roe`; `rms_ro` は最後まで有限)**。§9 に書いた 334/120/468/2 は forge の detectNaN ログの「Non-finite value detected in 'ro' at step N」(場の ro の NaN、1 から数えた反復番号) の値で、残差とは別の量。場の NaN の位置・節点数 (出口壁際 54 節点など) は計算時に AWS 上の res_nan_*.h5 で調べた値で、その h5 と forge_run.log は監査資料に未取得 (取得には AWS 起動が要る)。
