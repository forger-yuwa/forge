# メッシュの規模を抑える: 節点数の上限 (全経路・格子を作る前)・規模の事前の提示 (Euler の ni は NS と同じまま)

## メタ

- **area**: `tooling / design (格子の決め方)`
- **status**: `draft`  <!-- 2026-10-14 起票 (ユーザ「黙ってたら、50万点以上のメッシュのオイラーメッシュを作ってきた。aiがね。どうしたら対処できる？」→「あなたのおすすめでいいよ」) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「実寸の要求カードと成果物」の節 (格子の見積もり)
  - [`procedures/nozzle-design-workflow.md`](../../procedures/nozzle-design-workflow.md) §0・§1
- **related_plans**:
  - [`../accepted/verification-case45-euler-total-enthalpy.md`](../accepted/verification-case45-euler-total-enthalpy.md) (Euler 専用の格子 `mesh_euler`、全域 0.005)
  - [`../accepted/tooling-nozzle-wt-deliverables.md`](../accepted/tooling-nozzle-wt-deliverables.md) (要求カードの格子の見積もり `estimate_mesh`)
- **created**: `2026-10-14`
- **owner**: `sano`

## 1. 目的

AI が黙って 50 万点を超える Euler のメッシュを作った (ユーザ 2026-10-14)。原因の見立て (§3): 要求カードの生成器は NS の ni を壁際の AR の上限から見積もり (上限 6000)、
**Euler の ni に同じ値を写す**。完了時:

1. ~~Euler の ni が NS の ni から切り離され、Euler に要る軸方向の解像度で決まる~~ **取り下げ (2026-10-14、§4.1)**: V1b で Euler の解が ni に依存し (ni 600 で出口 M +3.8e-4 FAIL、
   900 は判定不能)、case/70 の 1 例から一般の上限を決められない (codex 諮問 mesh-size-guard-v1b・plan 段レビュー M1)。Euler の ni は明示が無ければ NS と同じまま。
2. 要求カードの生成器と axis-Mach の prepare (Euler・NS) で**実効の節点数を格子を作る前に**数え、上限を超えたら明示の許可と理由が無い限り止まる (AI が YAML を手で書いても prepare で止まる)。
   他の runner (`runner.py`・`runner_wt.py`・`runner_walldriven.py`・`euler_loop.py`) は対象外 (§5.1 #6)。
3. メッシュを作る前に節点数と所要の見込みをユーザに示し、上限を超えるならユーザの了承を取ってから許可を書くことが、AGENTS.md と skill の規則になる。

## 2. スコープ

- **やる**: 生成器の Euler の ni の決め方、prepare の節点数の上限、記録、試験、AGENTS.md・skill・手順書。
- **やらない**: NS の格子の決め方 (AR の制約と壁解像) の変更、格子の雛形の見直し (core grid、別セッション)。

## 3. 関連 docs と前提 (観測事実, 2026-10-14)

- 生成器 (`intake/wt_card.py`): NS の `mesh.ni` は `estimate_mesh` が「設計壁で壁際の AR ≤ 0.97 × ar_max (5000)」から見積もる (上限 `NI_CAP` 6000、`wt_card.py:67・727`)。
  第 1 層は y1+ 0.8 の見積もり。`mesh_euler.ni` は `None` にしてから生成器が `mesh.ni` を写す (`wt_card.py:562-563`、雛形 `axismach_final_v1.yaml` の「ni は生成器が mesh.ni を写す」)。
  nj は両方 97。したがって NS の ni が 6000 なら Euler も 6000 × 97 ≈ 58 万点。
- 節点数の上限はどこにも無い。品質の検査は AR ≤ 1000 (壁に沿った構造の層は 5000)・skew ≤ 0.9。
- 規模の例: 雛形の Euler 1600 × 97 ≈ 15.5 万、case/70 の NS 1800 × 97 ≈ 17.5 万 (第 1 層 1.15e-5・スロート 2.75e-6、AR の最大 約 4830)、M6 生産の NS 19.4 万。
- δ* の抽出 (`metrics/deltastar.py:465-473` `euler_q_at`) は Euler の ρU と壁の半径を NS の断面の x へ**列の間で線形**に補間する。Euler の列が NS より粗いと、
  補間の誤差 (曲率 × Δx²/8 程度) が δ* に入りうる — Euler の ni を NS に揃えていた理由である可能性 (§4.1 で確かめる)。
- **補間の成分は小さい (2026-10-14 実測、§6 V1a)**: Euler `run_0012b` (ni 1800) の列を 2・3・4 本おきに間引いて NS `run_0010b` を抽出
  (`case/70.m5_air_wt_demo/euler_ni_interp_ab.py`、Euler の解は同じで抽出が使う列だけを間引く): 使用域 [0, x_e] の max|Δδ_E| は
  1.7e-5・2.2e-5・2.7e-5 r_t (x ≈ 6〜11)、流量の比の変化 ≤ 1.5e-7。抽出の補間は ni を揃える理由にならない。残るのは Euler の解そのものの ni 依存 (§6 V1b)。
- 揃えた経緯: 雛形の注記「E4 は生産の G1 と同じ 2000 点」(`intake/recipes/axismach_explore_v1.yaml:48`)。case/45 の評価で NS の G1 (2000) に合わせただけで、
  verification-case45-euler-total-enthalpy は「2000 × 97 は今回の評価条件で、全ケースの汎用の既定にする根拠はない」と書いている。
- Euler の所要: 1790 × 97 で 3.56 ms/step、48,000 step で約 3 分 (`run_0013_euler_ns2`、手元の GPU)。

## 4. 設計方針

### 4.1 Euler の ni を切り離す

- **決定 (2026-10-14、codex 諮問 mesh-size-guard-v1b の C 案を採用)**: 切り離さない。明示が無ければ Euler の ni = NS の ni (今の生成器のまま)。
  明示の `overrides: {mesh_euler.ni: N}` は今どおり優先 (`wt_card.py:1497`; レビュー m4) し、その値にも §4.2 の上限を当てる。
  理由: V1b (§6) で Euler の解は ni に 2 次程度で依存し (出口 M の差 +1.53e-4 (900)・+3.79e-4 (600))、出口較正 (`Md_moc_offset`) は Euler から設計に入る
  (`runner_axismach.py:373-377`)。比の規則 (A 案) は事前の分岐の条件 (比 2 の合格) を満たさず、Δx の目標 (B 案) は最大間隔だけで局所の形・勾配を拘束できない。
  Richardson の見積もり (ni 1800 の M の格子誤差 約 5e-5) は仮説として残し、既定の根拠には使わない。
- 取り下げた案の記録: `mesh_euler.ni = min(mesh.ni, 2000)` (§6 V1b で棄却)。

### 4.2 節点数の上限 (安全網)

- `runner_axismach.prepare` (Euler) と `prepare_ns` (NS) で、**run dir を作る前・格子の配列を確保する前**に**実効の節点数**を数え、上限を超えたら止める (作った後の実際の節点数とも照合する)。
  実効の節点数は指定の ni × nj ではない: `axis_gap_frac` を書くと nj は導出値になる (`meshing/mesh2d.py:170-173`、例 ni 2000・nj 97・axis_gap_frac 0.005 で nj ≥ 201 = 40 万点以上)。
  そこで径方向の点数を求める処理を `mesh2d` の関数 (例 `node_count(prm)`) に切り出し、生成器 (`generate_axisym_mesh`)・事前の判定・要求カードの生成器の表示が同じ関数を使う (codex 諮問 mesh-size-guard、2026-10-14)。上限の既定: **2D 軸対称は Euler・NS とも 30 万点** (Euler の ni = NS の ni なので同じ値にする;
  雛形 1600〜2000 × 97 = 15.5〜19.4 万、M6 生産 19.4 万の 1.5 倍。2026-10-14 親の決定、ユーザには報告)。
  問題 YAML の `mesh_euler.max_nodes` / `mesh.max_nodes` に明示の値と `max_nodes_reason` (文字列) があれば、その値まで許す (理由が無ければ止める)。
  - **要求カードの生成器にも同じ判定を置く** (plan 段レビュー M3): `estimate_mesh` は AR の確認で NS の全格子を作る (`wt_card.py:871-875` → `_mesh_ar_full` `:716`)。
    その前に実効の節点数を数え、上限を超えたら全格子を作らずに止める (見積もりの節点数・ni・理由を示す)。カードの `overrides: {mesh.max_nodes: N, mesh.max_nodes_reason: "…"}`
    (Euler は `mesh_euler.*`) で許可する。
  - 点数の計算そのものも大きな配列を作らない (`_radial_fracs_capped` はリストを作る、`mesh2d.py:103`): 上限を超えた時点で数えるのをやめる実装にする
    (極小の `axis_gap_frac` でも巨大な確保をしない)。
  run の `prepare_info.json` と `MESH_QUALITY.txt` の近くに節点数と上限・理由を記録。3D は別の上限 (今回は 2D だけ; 3D は既定で上限なしのまま記録だけ)。

### 4.3 規模の事前の提示 (規則)

- AGENTS.md の「計算・実行ルール」と skill `nozzle-design`: メッシュを作る前 (生成器の出力・prepare の前) に、節点数と、手元の GPU での 1 万 step あたりの時間の見込みを示す。
  上限 (30 万) を超えるときは、ユーザの了承を取ってから `max_nodes` と理由を書く (了承なしに AI が書かない)。生成器の標準出力と YAML のコメントに節点数 (Euler・NS) を出す。

## 5. 実装ステップ

1. ~~§6 V1b~~ (済、§4.1 の決定の根拠)。
2. `mesh2d` の実効の節点数の関数 (上限で打ち切る)、生成器と prepare の上限・記録・試験。
3. AGENTS.md・skill・手順書・`methods/design/overview.md` (レビュー m4)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | plan 段の諮問とレビュー | 諮問は済み (判断: 2026-10-14・事前の判定は実効の節点数で、§6.1)。次は plan 段のレビュー | F |
| 2 | V1b: Euler の解の ni 依存 | **済 (2026-10-14)**: 比 3 は FAIL、比 2 は流量が判定不能 → 判断: 2026-10-14・切り離さない (C 案、§4.1)。延長の A/B は決定に要らないので回さない | O |
| 3 | 実装 | **済 (2026-10-14)**: `mesh2d.node_count`・`radial_node_count` (nj の決め方を `generate_axisym_mesh` と共有、上限で打ち切る)、`runner_axismach.MAX_NODES_2D`・`node_budget`・`check_node_budget` (run dir と設計チェーンの前に判定、作った後に照合、`prepare_info` に記録)、`wt_card` (全格子の前・ni の繰り上げごと・明示の Euler の ni)、試験 `run_mesh_size_guard_tests.py` 54 ok。V2・V3 合格 (§6) | O |
| 4 | 規則 | **済 (2026-10-14)**: AGENTS.md「メッシュの規模は作る前に示す」、skill `nozzle-design`「メッシュの規模」、`procedures/nozzle-design-workflow.md` §0 の表、`methods/design/overview.md` | O |
| 6 | 他の runner の格子生成 | `evaluate/runner.py:155`・`runner_wt.py:312`・`runner_walldriven.py:73`・`feedback/euler_loop.py:303` に判定が無い (axis-Mach の外の旧チェーン)。使う場面が出たら同じ `check_node_budget` を当てる | O |
| 7 | NS の nj の既定の食い違い (既存) | 生成器の `_mesh_ar_full` は nj の既定 81、`prepare_ns` は 97。雛形は両方 97 を書くので今は効かないが、雛形が nj を省くと生成器の数えが小さくなる (最後の網は prepare_ns)。既定を 1 か所に | O |
| 5 | result 段 | 諮問とレビュー | F |

## 6. 検証 (事前登録)

- **V1a 抽出の補間** (済、§3): Euler の列を間引いた抽出。max|Δδ_E| ≤ 2.7e-5 r_t (1/4)。
- **V1b Euler の解の ni 依存** (手元、Euler 3 本 + 抽出、各 約 3 分): case/70 の `problem_m5_air_c2b_euler.yaml` (r_t 0.024973561 m) の `mesh_euler.ni` だけを
  1800・900・600 に変えた問題で、`runner_axismach` の Euler を同じ段取り (soft + 本段 cfl 1・48,000 step、run_0013 と同じ 1 本通し) で回す
  (1800 は既存の run_0012 + 0012b の継ぎでなく 1 本通しで作り直し、3 本の段取りを揃える)。NS `run_0010b_ns_c2b_cont60k` の最終の場を 3 本それぞれに対して
  抽出 (`extract_and_merge(.., band_select="edge")`、今の既定)。
  見る量: 使用域 [0, x_e] の max|δ_E(ni) − δ_E(1800)|、出口コア M (`qs_series.py` の末尾 5 枚の平均) の差、流量の比の差。
  判定 (事前、ni 900 = 比 2 と ni 600 = 比 3 のそれぞれ): δ_E の差 ≤ 1e-4 r_t・出口コア M の差 ≤ 2.5e-4・流量の比の差 ≤ 1e-4。
  比 3 (= NS 6000 に対して Euler 2000) まで合格 → §4.1 の上限 2000 でよい。比 2 だけ合格 → 上限を NS の ni の 1/2 以上に (Euler の ni = max(min(mesh.ni, 2000), ⌈mesh.ni/2⌉)、
  このとき上限 25 万は ni 2577 まで)。比 2 も不合格 → 切り離さず、§4.2 の上限と §4.3 の提示だけを入れ、諮問に戻す。
  補足: 比べる 3 本はどれも残差が床で横ばい (run_0013 は NOT CONVERGED (stalled/plateau)) の Euler なので、出口コア M と流量は `check_quasisteady` で STEADY を確かめてから比べる。
  **結果 (2026-10-14)**: `run_0015_euler_c2b_ni1800`・`run_0016_euler_c2b_ni900`・`run_0017_euler_c2b_ni600` (NaN 0・品質 PASS・残差 NOT CONVERGED)。
  出口コア M (STEADY) 4.999365 / 4.999518 / 4.999744 → 差 +1.53e-4 (900)・+3.79e-4 (600)。流量 (`--drift 5e-5 --osc 5e-5` で 3 本とも OSCILLATING、1800 の振幅は相対 1.3e-4) の比の差
  −9.86e-5 (900)・−2.43e-4 (600)。max|Δδ_E| 6.0e-5 (900)・3.0e-5 (600) (単調でない)。**事前の判定: 比 3 は出口 M で FAIL、比 2 は δ_E・M 合格で流量は判定不能**
  (差が 1800 の振動の中)。差は 2 次の格子収束に乗る (M の差の比 2.48、流量 2.6; 2 次なら 2.67) → 解釈と既定の決め方を諮問 (`notes/reviews/briefs/2026-10-14-mesh-size-guard-v1b.md`)。
  **V1b の扱い (2026-10-14 諮問の結論)**: case/70 の条件に限った感度評価。ni 900 の流量は判定不能 (差 −9.86e-5 と許容 1e-4 の余裕 1.4e-6 は 1800 の振幅 相対 1.3e-4 より小さい)、
  δ_E も単一時刻なので合格とはしない。延長で判定可能にする A/B (諮問の提案) は §4.1 の決定 (切り離さない) に要らないので回さない。
- **V2 生成器**: 上限以下のカード (case/70 の c2b) は今と同じ生成物 (時刻の値を除く) で、YAML のコメントと標準出力に Euler・NS の節点数が出る。
  NS の ni が 6000 になるカード (第 1 層を薄くしたもの) は、**全格子を作る前に** (`_mesh_ar_full` を呼んだら失敗する差し替えで確かめる) 見積もりの節点数を示して止まる。
  カードに `mesh.max_nodes`・`mesh_euler.max_nodes` と理由を書けば通る (理由が無ければ止まる)。明示の `mesh_euler.ni` は保持され (上限以下・超過の両方)、超過なら止まる。
- **V3 上限**: 上限を超える問題で prepare が止まる (Euler・NS)、`max_nodes` と理由を書けば通る、理由が無ければ止まる、上限以下は今と同じ。
  **実効の節点数の A/B** (諮問 mesh-size-guard): Euler の ni 2000・nj 97・wall_first_frac 0.005・上限 25 万を固定し、`axis_gap_frac` だけを A: 未指定 / B: 0.005 に変える (forge 0 step)。
  合格: A は 194,000 点で通る、B は 40 万点以上と数えて**格子の生成 (`generate_axisym_mesh`) を呼ぶ前・run dir を作る前**に止まる (生成の関数を呼んだら失敗する差し替えで確かめる)。
  `node_count(prm)` が実際に生成した座標の数と一致すること (axis_gap_frac・axis_cap_frac・wall_first_frac_throat の各経路)。
- **V4 規則**: AGENTS.md と skill に §4.3 の規則があり、生成器の標準出力に節点数が出る。`methods/design/overview.md` に上限と許可の書き方。

**結果 (2026-10-14)**:
- **V3 合格**: `node_count` は 5 経路 (一様・wall_first_frac_throat・axis_cap_frac・axis_gap_frac 2 種) で実際の座標数と一致、HEAD 版の mesh2d と 9 通りで座標・接続・境界辺がビット一致。
  A/B (上限 25 万に差し替え): A (axis_gap_frac なし) 194,000 で通る、B (0.005) は 402,000 (ni 2000 × nj 201) と数えて生成の関数・設計チェーン・run dir の前に止まる。
  axis_gap_frac 1e-7 でも limit を渡せば 0.03 ms で返る。許可の値と理由の組 4 通り。
- **V2 合格**: case/70 の c2b カードは変更前の生成器と、節点数のコメント行と時間の値を除いて YAML・Euler の写し・card_log が一致。第 1 層 2e-6 / 5e-7 のカード (NS ni 6000 = 582,000 点) は
  `_mesh_ar_full` の前に止まる (見積もりの ni 9972、許可の書き方を示す)。mesh と mesh_euler の両方に許可を書けば通り、NS だけでは Euler で止まる。明示の `mesh_euler.ni` は保持 (1000 はそのまま、4000 は止まり許可で通る)。
- **V4 合格**: 上の文書。case/70 の `problem_m5_air_c2b_euler.yaml` を実際に `--prepare-only` → mesh 欄 n_nodes 174600・max_nodes 300000・source default、品質 PASS。
- 既存の試験: `run_wt_card_tests` (287 s)・`run_axismach_wall_tests`・`run_moc_exit_cut_tests`・`run_wt_card_profile_tests` ほか 12 本が通る。`run_wall_temperature_tests` V4 の prepare_info の比較から
  mesh の新しい 4 キーを外した (代わりに n_nodes = ni × nj を確かめる)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose (result の解釈) | 2026-10-14 | [2026-10-14-meshguard-fallback-result-diagnose.md](../../notes/reviews/2026-10-14-meshguard-fallback-result-diagnose.md) | 両 plan の accepted を保留: Major fallback V1 は無応答 2 回を補充で埋めた (事前の条件の読み替え)、却下 既存 YAML 7 本の停止は破壊でない、却下 V3 は道具の一覧で満たす、Minor mesh の完了条件の V1b の文言 | **採用** (Major と Minor): fallback は V1 を条件を固定して 1 回再試験、mesh は §8 の文言を直す |
| plan | 2026-10-14 | [2026-10-14-tooling-nozzle-mesh-size-guard-plan.md](../../notes/reviews/2026-10-14-tooling-nozzle-mesh-size-guard-plan.md) | GO-with-changes, C0/M3/m1: M1 比 3 の合格を全ケースの上限 2000 に一般化できない (配点は領域長・細分で決まる)、M2 準定常のゲートが許容差に対応していない (流量の振動が 1e-4 を覆う)、M3 要求カードの `estimate_mesh` が prepare の前に全格子を作る、m4 明示の `mesh_euler.ni` を黙って上書きしない | **全件採用** (諮問 mesh-size-guard-v1b): M1 → 既定を変えない (§4.1)、M2 → V1b は判定不能で閉じ、δ_E の合格を主張しない、M3 → 生成器にも格子を作る前の判定と打ち切る点数計算 (§4.2・V2・V3)、m4 → 明示の Euler の ni を保持し上限を当てる |
| diagnose (V1b の解釈) | 2026-10-14 | [2026-10-14-mesh-size-guard-v1b-diagnose.md](../../notes/reviews/2026-10-14-mesh-size-guard-v1b-diagnose.md) | 推奨 C (切り離さない; ガードと事前の提示を先に)。A (比の規則) は事前の分岐の条件を満たさない、B (Δx の目標) は最大間隔では局所を拘束できない。M2 は延長の A/B で判定可能か確かめる案。異議: δ_E の合格は保留、Richardson は既定の根拠にしない、case/45 の −0.0008 は ni 単独の証拠でない | **採用** (C 案、§4.1)。延長の A/B は決定に要らないので回さない (理由を V1b に) |
| diagnose | 2026-10-14 | [2026-10-14-mesh-size-guard-diagnose.md](../../notes/reviews/2026-10-14-mesh-size-guard-diagnose.md) | Major 1: 事前の判定を指定の ni × nj ですると `axis_gap_frac` の導出 nj (`mesh2d.py:170-173`) で過小評価し、上限を超えた格子を作ってから止まる | **採用** (親が `mesh2d.py:170-173` で確認): §4.2 を実効の節点数 (生成器と共通の関数) で格子の確保の前に判定、§6 V3 に axis_gap_frac の A/B |

## 7. 影響範囲

- 生成器の Euler の ni は変えない (§4.1)。prepare・prepare_ns・要求カードの生成器は上限 (2D 30 万点) を超えるときだけ止まる。設計チェーン本体・ソルバは変えない。
- **既存の問題 YAML で prepare_ns が止まるようになるもの** (全 310 件を走査、2026-10-14): case/44 `problem_va_R2_LU6_Lc8_ns_fine{,_cpg,_cpg_iso300,_iso300}.yaml` (4001 × 113 = 452,113)、
  case/45 `problem_d155_ns_finemesh3_pin.yaml` (435,000)・`problem_d155_ns_recal_final_cap0133.yaml`・`problem_d155_ns_recal_final_cap020.yaml` (2000 × 257 = 514,000)。
  回し直すときは `mesh.max_nodes` と理由 (ユーザの了承) を書く。main に入れるときに該当するセッションへ知らせる。
- `prepare_info.json` の mesh 欄に 4 キーが増える。要求カードの生成器の YAML と標準出力に節点数の行が 1 行増える (データ部分は不変)。

## 8. 完了条件

- §6 V1a (合格)・V1b (判定不能で閉じる: Euler の ni を切り離さない決定 (§4.1) の根拠であり、合格を主張しない; 追加の計算は不要)・V2〜V4 が合格、plan 段・result 段の codex レビューの行、AGENTS.md・skill・手順書の同期。

## 9. 変更ログ

- 2026-10-14: 起票。原因の見立て (Euler の ni が NS の ni を写す) と、抽出の補間が Euler の列を使うこと (§3) を記録。
- 2026-10-14: V1a (抽出の補間) を実測: 列を 1/4 にしても max|Δδ_E| 2.7e-5 r_t。V1 を V1a (済) と V1b (Euler の解の ni 依存、3 本) に分けた。
- 2026-10-14: 諮問 (mesh-size-guard) を採用: 事前の判定を実効の節点数 (axis_gap_frac の導出 nj を含む) で格子の確保の前に、§6 V3 に A/B。
- 2026-10-14: V1b を実施 (`run_0015`〜`0017`): 比 3 は出口 M で FAIL、比 2 は流量が判定不能。plan 段のレビュー (GO-with-changes, M3/m1) と合わせて諮問へ。
- 2026-10-14: 諮問 (mesh-size-guard-v1b) を採用: **Euler の ni は NS から切り離さない** (目的 1 を取り下げ)。上限は 2D で Euler・NS とも 30 万点、要求カードの生成器にも格子を作る前の判定、点数の計算は上限で打ち切る、明示の Euler の ni は保持。plan 段レビュー M1〜M3・m4 を全件採用。
- 2026-10-14: §5.1 #3・#4 済。V2〜V4 合格 (試験 54 ok、既存の試験も通る)。既存の YAML 7 本が上限で止まるようになる (§7)。残作業 #6 (他の runner)・#7 (nj の既定の食い違い)。
- 2026-10-14: result 段の諮問 (meshguard-fallback-result) を採用: 完了条件の V1b を「判定不能で閉じる (切り離さない決定の根拠)」に直した。既存の YAML 7 本が止まるのは目的どおり (却下側)、main に入れるときに該当セッションへ知らせる。
