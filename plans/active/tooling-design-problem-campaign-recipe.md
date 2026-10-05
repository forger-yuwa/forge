# 設計ツールの入力構造: problem (何を作るか) / campaign (どう探し・どう評価するか) / recipe (機種別の実行手順)

## メタ

- **area**: `tooling / optimization`
- **status**: `draft`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (§問題定義 YAML — 本 plan の移行中である旨を注記)
  - [`design/CAPABILITIES.md`](../../design/CAPABILITIES.md) (機種・メニューの対応表)
- **related_plans**:
  - [`tooling-nozzle-axismach-length-dv.md`](../accepted/tooling-nozzle-axismach-length-dv.md) (§2 の「ドライバは dv dict に対して汎用」はコードと矛盾 — 本 plan §5.1 #8 で訂正)
  - [`tooling-nozzle-deltastar-core-matched-euler.md`](../accepted/tooling-nozzle-deltastar-core-matched-euler.md) (δ\* 補正の生産レシピ = `axismach.deltastar_ns/v1` の中身)
  - [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) (SERN の多作動点 MOO — 後続で同じ契約に載せる)
- **created**: `2026-09-27`
- **owner**: `Claude (主セッション) / ユーザ`

## 1. 目的

forge_design を「ユーザと対話しながら、機種 (風洞 axismach / ベル+スラスタ / デュアルベル / SERN、将来は熱制約・FEM) ごとの
最適化フローを組める」基盤にする。現状は **YAML に書いた内容と実際に動くコードの間に検査できる約束事が無い**ため、
書いても無視される設定・黙って入る既定値・入口の名前で決まる物理が生じている (§3)。
完了時には、axismach で「problem (基準値) + campaign (探索・評価要求)」から、何が回るかを実行前に解決・表示・拒否でき、
「Euler で探索 → 上位 N 点を δ\* 補正 NS で評価 → 補正後で最終順位」が 1 つの campaign で回る状態にする。

## 2. スコープ

- **やる**:
  - スキーマ `forge.problem/v1` (何を作るか: 仕様・ガス・形状パラメータの基準値・作る物の物理条件) と
    `forge.campaign/v1` (探索: dv と範囲・目的・制約・探索法 / qualification: 上位の選び方・評価手順の列・最終順位の決め方) の骨格。
  - 機種 adapter の公開契約: 変数一覧 (名前・単位・独立/従属・対応モード)、評価手順 (recipe) 一覧と各 recipe の必須入力・出力・検査する VERDICT。
  - `resolve_evaluation(problem, candidate, recipe)`: 副作用なしで実行計画 (形状・メッシュ・実効 solverConfig/bcond・段・必須入力・検査予定 VERDICT・値の出所) を返す。
  - axismach の recipe 2 本: `axismach.euler/v1` (現 `prepare` + 段階起動) と `axismach.deltastar_ns/v1` (現 `deltastar_loop` の生産レシピ)。
  - 既存の不具合修正: `run_staged` の同一メッシュ段間引き継ぎを `restart_field.py` に、semiperfect の γ/cp の扱いの分離、dv の min/max の置き場所。
  - 利用例 1 本で端から端まで検証: case/44 va3 M4.19 で L_c を Euler 探索 → 上位 2〜3 点を δ\* 補正 → 補正後で最終順位。
- **やらない** (別 plan / 後続):
  - 機種非依存の汎用 stage DAG エンジン (codex 2026-09-27 諮問で却下: δ\* は壁を変えて再評価する反復で段の直列に乗らない)。
  - 既存 problem YAML の一括移行 (旧形式は読み続ける。新形式は新規 campaign から)。
  - ベル (`opt/driver.py`)・SERN (`opt/driver_sern.py`) の adapter 化 (axismach で契約を確かめてから、§5.1 #9 以降)。
  - 熱構造 FEM・CHT の recipe (CAPABILITIES §4 で 📋。契約上の置き場所だけ §4.5 で決める)。
  - MOO (サロゲート) を axismach に繋ぐこと (初回の探索法は grid。§5.1 #10)。

## 3. 関連 docs と前提 (観測事実, 2026-09-27)

根拠の詳細はブリーフ [`notes/reviews/briefs/2026-09-27-design-problem-schema.md`](../../notes/reviews/briefs/2026-09-27-design-problem-schema.md) の観測事実 1〜10。

- semiperfect でも `gas.gamma`/`gas.cp` を書かせ、省略すると黙って γ 1.4 / cp 1004.5 (`probdef.py:167-168`)。
  値は forge `physProp` への素通し (TP では無効)・出口一様性診断への引数 (`runner_axismach.py:581-582`; ただし `core_radius_traced` は出力の `sonic` から M を作り
  `gamma` を使っていない = 未使用引数, codex plan レビュー m7)・NS の ω 床初期化 (`runner_axismach.py:820-829`; **TP の `roe` から CPG の式で T を逆算しており不整合**。
  codex の `run_0509` 場での検算で 57 % のノードが 50 K 下限に張り付く, plan レビュー M1)・`cfd_gas: cpg` の CPG 定数・SERN `frozen_tp` の設計側 `GasCPG` で使われる。
- `dv.L_c.min/max` は axismach の探索に使われない (`probdef.py:221-226` の範囲検査と MOO ドライバだけが読む)。
  runner は自分で計算した許容範囲 (va3 では `[0.01, 10.94]`) で検査する。YAML の `max: 14.0` は設計上ありえない値のまま。
- MOO ドライバは機種ごとに別実装で runner を固定 import (`opt/driver.py:40,77`, `opt/driver_sern.py:41`)。風洞の L_c 最適化は
  problem YAML を読まない使い捨てスクリプト (`case/42.isobutane_wt/optimize_axislaw_A_shortest.py`, 条件・範囲がスクリプト内定数)。
- Euler/NS は YAML で表せず入口で決まる (`runner_axismach` CLI は常に Euler + slip、NS は `feedback.deltastar_loop` からだけ)。
- `--prepare-only` では YAML の `cfl_main` が効かず準備時 cfl 4.0 (`runner_axismach.py:477,620`)。同じ problem で入口により実効設定が変わる。
- `run_staged` は同一メッシュの段間引き継ぎに `interp_field.py` を使う (`runner_axismach.py:538-540`)。AGENTS.md 違反。
  case/44 run_0509–0511 では runner を prepare で止め、段階起動を case 内スクリプト `run_lumpX_staged.py` に置き換えた。
- 数値設定の正本 `procedures/recommended-settings.md` の §1/§3 (cfl 6+relax 0.7, nStepInner 4) と §4 Euler 設計評価 (cfl 4 / nStepInner 5) が
  食い違う (codex 諮問 Minor)。recipe は適用範囲と版を持たせ、どちらを採るかを recipe 側で明示する。
- 方向の判断: codex (diagnose) に諮った — [`notes/reviews/2026-09-27-design-problem-schema-diagnose.md`](../../notes/reviews/2026-09-27-design-problem-schema-diagnose.md)
  — 結論「汎用 DAG ではなく機種別 runner + 小さな共通契約。まず axismach に副作用なしの `resolve_evaluation` を作り Euler/δ\* NS の計画差分で境界を検証」。
  採否表は同記録の末尾。ユーザ承認 2026-09-27 (「現状方針でいきましょう」)。

## 4. 設計方針

### 4.1 3 つの置き場所

| 置き場所 | 持つもの | 持たないもの |
| --- | --- | --- |
| `problem.yaml` (`forge.problem/v1`) | 機種 `kind`、仕様 (Pt, Tt, M_design, r_throat …)、ガス (モデル・組成・入力基準)、形状パラメータの**基準値** (例 `L_c: 8.0`)、**作る物の物理条件** (壁の熱条件・乱流の有無など、NS 評価で必須になるもの) | 探索範囲 (min/max)、目的・制約、評価手順の選択、数値設定 |
| `campaign.yaml` (`forge.campaign/v1`) | `problem` への参照、seed、`search` (探索する dv と範囲・目的・制約・探索法・探索中の評価 recipe)、`qualification` (上位の選び方・評価手順の列と入力参照・最終順位の決め方) | 物理条件・数値設定の直書き |
| recipe (機種別コード, 名前 + 版 `axismach.euler/v1`) | 形状生成・メッシュ・段階起動・δ\* 反復・停止条件・数値設定 (推奨設定のどの節に従うかを明記)・検査する VERDICT | 探索範囲 |

- 1 点だけ回す (今回の再計算のような) ときは `search` の無い campaign (`evaluations` だけ) にする。
- 探索範囲には 2 種類ある: **ユーザの探索範囲** (campaign に書く) と**設計が成立する範囲** (adapter が候補ごとに計算して検査する、人は書かない)。
  設備の上限 (全長など) は campaign の制約に書く。
- 旧形式 (1 ファイルに dv の min/max を含む) は読み続ける。新形式では problem に min/max があれば**拒否**する (黙って無視しない)。

### 4.2 adapter の公開契約

機種 adapter (初回は axismach のみ) は次を返す関数を持つ:

- `variables()`: 変数名・単位・独立/従属・有効な `Lc_mode` などの対応モード。campaign はここに無い変数を探索できない。
  `Lc_mode: from_length` のとき `L_c` の探索は拒否 (従属量なので)。
- `recipes()`: recipe 名・版 → 必須入力 (他の評価結果の参照を含む)・problem に必須の物理条件・出力 (評価量と単位)・検査する VERDICT。
- `feasible_range(problem, name)`: 設計が成立する範囲 (L_c の許容範囲など)。

### 4.3 `resolve_evaluation(problem, candidate, recipe) -> ExecutionPlan`

- **副作用なし** (forge・メッシュ変換・ファイル書き込みを呼ばない)。
- 返すもの: 形状の識別子 (Euler 設計壁か δ\* 補正壁か + 生成元)、メッシュ条件、**実効の** solverConfig/bcond (壁 BC、CFL、段、γ/cp の出所)、
  必須入力とその充足状況、検査予定の VERDICT (種類・対象量・判定区間)、各値の出所 (problem / recipe 既定 / 推奨設定の節)。
- 矛盾は実行前に拒否する: NS recipe なのに壁の熱条件が無い、`from_length` なのに L_c を探索する、problem に min/max がある、必須入力の欠落 など。
- 実行時は解決済みの計画を run ディレクトリに `execution_plan.json` として保存する (来歴)。実行は計画だけを読む (CLI・case 内の後編集に依存しない)。
- **事前に決まるものと実行後に決まるものを分ける** (codex plan レビュー M2): campaign 開始時の検証では、入力を「充足」「将来生成される参照」
  (例: δ\* pass 1 の壁は pass 0 の NS 場から作られる) 「欠落」に分類し、欠落だけを拒否する。具体的な実行計画は recipe が**各 pass・各段の直前に確定**して保存する。
- 成果物の参照は run ディレクトリ名でなく **candidate・形状 ID・メッシュ・輸送種配置・エネルギー基準・特定 snapshot と内容ハッシュ**で持ち、不一致を拒否する
  (現状の `ic_from` は run ディレクトリの最大番号 `res_*.h5` を拾い、継続計算後は別の場を読む: `runner_axismach.py:805`)。
- 段ごとの実効設定と履歴は既存 `StageManifest` (`stage_manifest.json`) に保存し、**全段**を計画と照合する (最終 solverConfig だけを見ない)。

### 4.4 VERDICT の扱い

- 生の VERDICT (`NOT CONVERGED` など) は書き換えない。campaign の採否は「どの VERDICT を合格条件にしたか」を別に記録する。
- 要求された評価の未実行・判定ファイルの欠落は不合格。
- qualification の用途判定は生の VERDICT と別に `accepted / rejected / incomplete` の 3 値で持つ (codex plan レビュー M3)。
  δ\* のゲート (ṁ_NS/ṁ_E, 出口コア M) は **NS・Euler 両系列の定常性**を確認したうえで判定し、反復上限でも不合格なら `rejected` として最終順位から除く。
  判定ファイルの欠落・solver 失敗は `incomplete` (合格にしない)。全候補が不合格なら「勝者なし」を返す。
  最終順位は補正後の目的量・制約を**再評価**して決める (Euler の順位を流用しない)。既存 `run_pass` は solver 戻り値を表示して先へ進み収束検査も `check=False`
  (`deltastar_loop.py:156`) なので、recipe 側で判定を強制する。
- δ\* 補正後の壁は**別の形状**として扱い、後段 (凝縮・3D・FEM) がどちらの形状を評価するかは入力参照で固定する。

### 4.5 熱制約・FEM の置き場所 (今回は実装しない、契約だけ決める)

- 熱が**制約**なら探索中の評価 (`search.evaluator` か制約評価) に入れる。勝者だけに FEM をかけても熱制約付き最適化にはならない。
- 片方向 FEM は追加 recipe、双方向 CHT は界面反復を所有する別 recipe。

### 4.6 γ/cp の分離

- semiperfect では熱力学の本体は `gas` のモデル (NASA-9)。**TP の R は組成から直接取り、温度は組成・エネルギー基準 (`thermoHrefTemp`) と整合した EOS で復元する**
  (NS の ω 床初期化の CPG 逆算はこれで置き換える; codex plan レビュー M1)。
- 定数 γ/cp が本当に要る箇所 (`cfd_gas: cpg`・`frozen_tp` の設計側 `GasCPG`) は**近似用の参照定数**として別名で扱う。**γ* と cp(Tt) のように別温度の値を組み合わせない**
  (va3 で R が 292.59 vs モデル 285.27 = +2.6 %)。同一参照温度の γ/cp の組、または R と一方の定数で持つ。ユーザ記入値とモデル値が食い違えば警告。
- 暗黙の既定 (1.4 / 1004.5) への fallback は semiperfect では禁止 (CPG では従来どおり)。
- 出口一様性診断の `gamma` は未使用引数なので整理だけする (ガスモデル経由化は不要; plan レビュー m7)。
- 試験は 2 種類に分ける: **熱力学の整合性試験** (R・T の復元が EOS と一致) と、**不具合修正に伴う結果変化の評価** (ω 床修正後の NS が旧 run からどれだけ動くか; 「旧結果不変」を合格条件にしない)。

### 4.7 評価量の定義は 1 つの抽出関数で固定する (codex plan レビュー M4)

評価量ごとに**抽出位置・補間格子・積分重み・正規化・版**を固定し、探索・時系列判定・旧 run 比較を同じ関数で行う。現状は runner (`runner_axismach.py:572`,
軸 200 点補間・出口は x_E) と case 側 (`case/44.vitiated_air_wt/lumpX_series_csv.py:13,30`, 軸ノード上・出口は x_max−2 r_t) で定義が違い、
同じ `run_0509` の場で軸 M 目標差 max が 0.002129 vs 0.002120 (0.41 %) と食い違う。`run_0509` は**未収束 (残差 plateau) の準定常回帰参照**として扱う。

### 4.8 メッシュ条件は recipe ごとに固定し、壁解像をゲートにする (codex plan レビュー M5)

- Euler と NS でメッシュ条件 (特に `wall_first_frac`) は別物。現状 `prepare_ns` は problem の `mesh.wall_first_frac` を優先し (`runner_axismach.py:746-747`)、
  Euler 用 0.005 がそのまま NS に渡る (NS 既定 4.5e-5 の約 111 倍)。メッシュ条件は **recipe の持ち物**にし、problem の `mesh` は recipe ごとの節に分けるか recipe 既定を上書きする明示キーにする。
- δ\* 補正の各 pass で: メッシュ品質 (`check_mesh_quality.py`) 不合格なら投入停止、cross-mesh IC の適合性確認、低 Re SST の**局所 y₁⁺** (`check_wall_resolution.py`; ソルバ `ypls` は使わない, AGENTS.md 壁解像確認)
  が不足なら qualification を `rejected`。

### 4.9 recipe `axismach.contur_c2/v1` — M6 で確定した C2 方式 (codex 諮問 2026-10-05 で条件付き採用)

出典: `plans/active/verification-m6-axis-wave-mesh-su2.md` §5.1 #8f・#9、諮問 `notes/reviews/2026-10-05-m6-wall-fit-and-design-pipeline-diagnose.md`。
ユーザ要求 (2026-10-05): マッハ数を変えて同じ形状を作り直せること、軸長・出口径・入口径の決め打ちとパレート探索の両方。

- **手順** (recipe の内部): 逆 MOC → Euler 参照 → CONTUR 積分法の δ で物理壁 → NS → E 法 (`band_select="edge"`, 測定器) で出口 δ を測る → CONTUR の `cf_scale` (k_f) を出口 δ に合わせる →
  出口半径の仕様から r_t を解く → 最終 NS → (後段) 凝縮 ON の NS → `nozzle_report` (pptx)。
- **反復は recipe 内部に置く**: r_t 解き・k_f 較正・再測定・停止判定。campaign は候補の探索と評価の順序だけを持つ (§4.1 の分担と同じ)。
  r_t 比の −0.2 乗による δ の換算は**初期予測**であり、最終 NS で出口 δ を再測定して閉じる (残差が許容外なら反復、上限で未達なら `rejected`)。
- **寸法固定は r_t 解きだけでは閉じない**: `from_length` の `L_total` は r_t 単位で設計スロート起点 (`runner_axismach.py:386`)。実寸の全長・入口径を固定するなら、
  r_t を更新するたびに `L_total`・`r_inlet` も換算し直し、長さの起点 (設計スロートか物理スロートか) を problem のキーで定義する。
  recipe は寸法残差 (全長・出口半径・入口半径) と出口 δ 残差をともに検査する。現 `problem_d155_ns_c2final.yaml` の入口半径は 0.4980966 m でコメントの 0.5 m と合っていない (同期漏れの実例)。
- 凝縮評価は確定した物理壁を参照する後段 recipe として campaign に置く (壁を作り直さない)。
- 入力成果物の来歴 (Euler 参照・δ 表・NS snapshot) は §4.3 のハッシュ参照で持つ。
- **合格の分け方**: 「旧結果の再現」(同一入力・物性・バイナリで run_0051 系列を再現) と「修正後の設計資格」(用途上の絶対ゲート: 波 ≤ 0.01 %・オーバーシュート ≤ +0.035 %・出口コア M 6 ± 0.02 %、
  残差・準定常・壁解像の各 VERDICT) を別に判定する。main の物性変更や TP 初期化修正による差は再現ノイズに含めず「変更影響」として報告する。
  `run_0051` は残差 plateau・壁解像 FAIL・波/オーバーシュート DRIFTING なので**歴史的な回帰参照**であって、修正後の recipe にこれとの数値一致を課さない。
- 再現許容差は量ごとに max(同一条件 3 反復の最大差 × 3, 事前登録の絶対下限)。別 case のノイズ実測 (§5.1 #2 の va3) を M6 に流用しない。波・オーバーシュートの下限は別に決める。
  許容差が用途上の余裕を食い潰すなら判定不能とする。
- **実装順** (諮問の推奨を採用): M6 の旧入力・物性・壁・メッシュ・初期 snapshot・バイナリのハッシュ保存と再現幅の取得 → #4 評価量統一 → #3 熱力学修正 → #6・#7 の最小契約 →
  search 無しの C2 campaign → 壁表現の変更 (verification-m6 §5.1 #13) → 寸法固定・パレート探索。

## 5. 実装ステップ

1. スキーマと検証 (`design/forge_design/probdef.py` 周辺に `campaign.py` を新設)。旧形式の読み込みは維持。
2. axismach adapter の `variables()` / `recipes()` / `feasible_range()` (`evaluate/runner_axismach.py` から抽出)。
3. `resolve_evaluation` と `execution_plan.json`。
4. recipe `axismach.euler/v1` (段間は `restart_field.py`)、`axismach.deltastar_ns/v1` (`feedback/deltastar_loop.py` を呼ぶ)。
5. campaign 実行器 (grid 探索 → qualification → 最終順位、ledger)。
6. case/44 で端から端まで。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~plan 段 codex レビュー~~ | 完了 2026-09-27 (GO-with-changes, C0/M5/m2, 全件採用 → #2–#7 と §4.3/4.4/4.6/4.7/4.8/§6 に反映; §6.1) | O |
| 2 | ~~V0 ノイズ床の取得~~ | 完了 2026-09-27 (種 DB plan の基準として取得; 抽出は現行 case 側 `lumpX_series_csv.py`、#4 の統一抽出関数ができたら再抽出して併記): `run_0513`–`0515` の 3 反復で許容差 ṁ 1e-4 相対・M 1e-5・T 0.01 K・軸 M 目標差 1e-5 (下限が効く)。詳細は `thermophysics-solver-owned-species-db.md` §6 冒頭 | O |
| 3 | 熱力学の修正 | §4.6: TP の R を組成から・NS ω 床の T を EOS 整合で復元、γ/cp の参照定数化・fallback 禁止・未使用 `gamma` 引数整理。合格: 熱力学整合性の単体試験 (va3 で R 285.2704、IC 場で逆算 T と `res` の T の差 ≤0.01 K)、CPG / semiperfect / `frozen_tp` の読込試験。NS への影響は旧 run 比の変化量を記録 (不変を合格条件にしない) | O |
| 4 | 評価量の抽出関数 | §4.7: 軸 M 目標差・軸 M 出口・出口質量流束平均 M/T・ṁ を 1 モジュールに定義 (版付き)。runner の `collect` と case の系列スクリプトをこれに置換。合格: `run_0509` で両旧定義との差を記録し、新定義で `check_quasisteady --series-csv` を再判定 | O |
| 5 | ~~`run_staged` の段間引き継ぎ修正~~ | 完了 2026-09-27 (種 DB plan #3b で実施): `run_staged`/`run_staged_ns` の段間を restart_field (ビット一致、species ハッシュの継承つき) に。配管 run `case/44.vitiated_air_wt/run_0520_species_attrs_runner`。旧 runner 経路の run との本段の差は未計測 (丸め程度の見込み) | O |
| 6 | スキーマ骨格と検証 | `forge.problem/v1` / `forge.campaign/v1`。problem の min/max 拒否、未知キー拒否、recipe ごとの `mesh` 節。合格: 単体テスト (正例 2・負例 5 以上) | O |
| 7 | axismach adapter 契約 + `resolve_evaluation` | §4.2/§4.3 (入力の 3 分類、pass ごとの計画確定、snapshot ハッシュ参照、全段照合)。合格は §6 V1 | O |
| 8 | recipe 2 本と campaign 実行器 | §4.4 の 3 値判定・勝者なし・補正後の再評価と最終順位、§4.8 のメッシュ/y₁⁺ ゲート。合格は §6 V2 | O |
| 9 | length-dv plan の記述訂正 | `plans/accepted/tooling-nozzle-axismach-length-dv.md` §2 の「ドライバは汎用」に訂正注記 | O |
| 10 | ベル・SERN の adapter 化 | axismach で V1/V2 が通った後。方針は別 plan で | F |
| 11 | axismach へのサロゲート MOO 接続 | 探索法の追加。grid で足りない需要が出てから | F |
| 12 | result 段の解釈と codex result レビュー | V1/V2 の結果解釈を上位に諮ってから `--stage result` | F |
| 13 | M6 の旧基準を固定 (C2 recipe の前提) | `case/45.isobutane_m6_d155/run_0051_ns_final_c2` 系列の入力 (problem・δ 表・Euler 参照)・物性・壁・メッシュ・初期 snapshot・バイナリのハッシュを保存し、同一条件 3 反復で評価量 (§4.7 の関数) の再現幅を取る。合格: 3 反復が完走し量ごとの許容差を §6 に登録 | O |
| 14 | recipe `axismach.contur_c2/v1` (§4.9) | #3・#4・#6・#7 の後。search 無しの campaign で M6 を作り直す。合格: 寸法残差・出口 δ 残差の検査、3 値判定、用途上の絶対ゲート、旧基準との差を「再現差」「変更影響」に分けて報告 (§4.9)。合否条件の数値は #13 の後に §6 へ事前登録 | F |
| 15 | 寸法固定・パレート探索の campaign | #14 の後。寸法固定 (全長・入口径・出口径を実寸で固定、長さの起点を problem キーで定義) と、Euler で探索 → 上位だけ C2 で qualification | F |

## 6. 検証

事前に決める参照値・許容差 (結果を見てから変えない):

- **V0 ノイズ床** (§5.1 #2, 数値を変える作業の前に取る): 同一バイナリ・同一初期場・同一実効設定で 3 回反復し、§4.7 の抽出関数で量ごとに
  許容差 = max(反復差の最大 × 3, 絶対下限 ṁ 1e-4 相対・M 1e-5・T 0.01 K)。1 回の差をそのまま許容差にしない (codex plan レビュー m6)。
- **V1 (配線の検証, CFD 0 step)**: va3 M4.19 L_c8 の同じ候補に `axismach.euler/v1` と `axismach.deltastar_ns/v1` の実行計画を解決して項目比較する。
  合格: Euler 計画は slip 壁・Euler 設定・Euler メッシュ条件、NS 計画は no-slip + problem の壁熱条件 + NS メッシュ条件 (`wall_first_frac` を数値で照合) + δ\* 入出力に解決され、
  **説明できない設定差 0・無視された有効入力 0・共通層に入った物理の分岐 0**。pass 1 の壁が「将来生成される参照」として分類されること。
  負例: NS 計画で壁熱条件を消す → 拒否、problem に min/max → 拒否、`from_length` で L_c を探索 → 拒否、snapshot ハッシュ不一致 → 拒否。
  どちらかが CLI や case 内の後編集に依存したら不合格。
- **V1b (判定ロジック, CFD 0 step)**: 合成した評価結果で「Euler 順位が NS で逆転」「判定ファイル欠落 → incomplete」「全候補 rejected → 勝者なし」を試験する。
- **V2 (端から端まで)**: case/44 va3 M4.19 で campaign 1 本: `search` = grid L_c ∈ {7.0, 7.5, …, 10.0} を `axismach.euler/v1`、目的 = 軸 M 目標差 max 最小 (§4.7 定義)、
  `qualification` = 上位 2 点に `axismach.deltastar_ns/v1`、最終順位 = 補正後の同じ目的量。合格:
  - 各 Euler 評価: NaN 0、評価量が `check_quasisteady` STEADY (残差は plateau でも可、生の VERDICT はそのまま記録)。L_c8 の Euler 評価が `run_0509` と V0 許容差以内。
  - 各 δ\* 評価: メッシュ品質 PASS、局所 y₁⁺ ≤1 (面積割合と位置を報告)、NS・Euler 両系列 STEADY のうえでゲート |ṁ_NS/ṁ_E − 1| ≤ 0.3 %・出口コア M ±0.1 %
    (deltastar-core-matched-euler plan の値)。**ゲート不合格は `rejected` として記録され最終順位から除かれること** (記録しただけでは V2 合格にならない)。
  - 最終順位が補正後の目的量で決まり、少なくとも 1 候補が `accepted`。全候補 `rejected` なら V2 は不合格 (勝者なしの挙動は V1b で確認済みとする)。
  - すべての run ディレクトリに段ごとの `execution_plan.json` と `stage_manifest.json` があり、全段の実効設定が計画と差分 0。
- **単体**: スキーマ検証の正例/負例テスト、`resolve_evaluation` の決定性 (同入力で同出力)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-27 | [2026-09-27-tooling-design-problem-campaign-recipe-plan.md](../../notes/reviews/2026-09-27-tooling-design-problem-campaign-recipe-plan.md) | GO-with-changes, C0/M5/m2 | 全件採用。M1 (γ*/cp(Tt) の R 不整合・NS ω 床の CPG 逆算) → §4.6・#3 (R 292.59 vs 285.27 と `runner_axismach.py:822-829` を当方で再現)。M2 (事前計画と実行後成果物の区別・snapshot ハッシュ) → §4.3・#7。M3 (qualification の 3 値判定・勝者なし・補正後の再評価) → §4.4・#8・V1b/V2。M4 (評価量の定義差 0.41 %) → §4.7・#4。M5 (NS メッシュに Euler の `wall_first_frac` が渡る・y₁⁺ ゲート) → §4.8・V1/V2 (`runner_axismach.py:746-747` を再現)。m6 (V0 を 3 反復×安全係数に) → #2・V0。m7 (出口診断の `gamma` は未使用) → §3・§4.6 訂正。判断役 (codex) 自身の指摘で却下が無いため、採否の別途諮問は省略 |
| diagnose | 2026-10-05 | [2026-10-05-m6-wall-fit-and-design-pipeline-diagnose.md](../../notes/reviews/2026-10-05-m6-wall-fit-and-design-pipeline-diagnose.md) | C0 / Major 4 (本 plan 関係: B1 C2 recipe 条件付き、B2 順序変更、B3 再現と資格の分離、寸法固定の連成) | 全件採用 → §4.9・§5.1 #13–#15 |

## 7. 影響範囲

- `design/forge_design/probdef.py`, 新規 `design/forge_design/campaign.py` (仮), `evaluate/runner_axismach.py`, `feedback/deltastar_loop.py`
- `design/CAPABILITIES.md` (契約・recipe 一覧を追記)、`methods/design/overview.md` §問題定義 YAML、`.claude/skills/design-intake/SKILL.md` (インテークの出力を problem + campaign に)
- 既存 case の旧形式 YAML は変更しない

## 8. 完了条件

- [ ] 関連 `methods/design/overview.md` の現在仕様を更新済み
- [ ] 実装・検証完了 (本計画の §6 を満たす)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [ ] ファイルを `plans/active/` → `plans/accepted/` へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期

## 9. 変更ログ

- `2026-10-05` — ユーザ要求 (M6 設計の再現・別マッハ・寸法決め打ち/パレート) を受け、codex (diagnose) に諮って C2 方式を recipe `axismach.contur_c2/v1` として §4.9 に追加。実装順は旧基準の保存 → 評価量統一 → 熱力学 → 最小契約 → C2 → 壁表現 → 寸法固定/探索 (§5.1 #13–#15)。

- `2026-09-27` — codex plan 段レビュー (GO-with-changes, C0/M5/m2) を全件採用し §3・§4.3/4.4/4.6/4.7/4.8・§5.1・§6 を改訂。実装順は ①熱力学 ②pass ごとの計画と成果物参照 ③合否・最終順位 ④評価量と V0 ⑤メッシュ・壁解像ゲート (V0 取得は数値変更の前)。
- `2026-09-27` — 初稿。codex diagnose 諮問 (`notes/reviews/2026-09-27-design-problem-schema-diagnose.md`) の推奨とユーザとの対話
  (dv は campaign、Euler/δ\* NS の切り替えは campaign の recipe 選択、NS に要る物理条件は problem、探索後に上位だけ δ\* 補正 = qualification) を方針として確定。ユーザ承認。
