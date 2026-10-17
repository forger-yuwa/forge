# 諮問: 格子の節点数の上限の承認 (§5.1 #8) の実装と V2a〜V2c の事前登録 (2026-10-18)

plan `plans/active/tooling-design-orchestrator.md` §4.6・§5.1 #8・§6 V2a〜V2c。関連: `plans/accepted/tooling-nozzle-mesh-size-guard.md` §4.2・§4.3。

聞きたいこと: 下の設計案 (1)〜(6) と事前登録の案 (V2a〜V2c) の穴。特に (Q1) 上限を問題 YAML に書く経路と承認の照合の二重化で
抜け道が残らないか、(Q2) V2a を「#8 の前の hand_mesh と後の hand_mesh + 承認の A/B」で決めてよいか (V6 の格子との比較にしない理由は下)、
(Q3) V2b の「雛形の格子の段」を recipe に入れずに、段の契約の関数として試してよいか。

## 1. 観測事実 (コードを読んで確かめた)

- `case/71.m55_combustion_wt/_mesh/hand_mesh.py:354〜357`: build は数えた節点数から上限を `ceil(n/1000)·1000` に自動で広げ
  (NS 442,520 点 → 443,000、Euler 466,570 → 467,000)、理由に固定の文「ユーザ決定「a」2026-10-15: …」を書く。承認の記録は読まない。
- `design/forge_design/evaluate/runner_axismach.py:791` `node_budget`: 上限は `max_nodes` が正の整数・`max_nodes_reason` が空でない文字列で
  あれば受け入れる (理由の中身・出所は見ない)。既定 `MAX_NODES_2D` = 300,000。`check_node_budget` が格子を作る前に照合。
- 要求カード (`design/forge_design/intake/wt_card.py:771`): v1 でも v2 でも、`overrides` の `mesh.max_nodes`・`mesh.max_nodes_reason`
  (と `mesh_euler.*`) で上限を広げられる (理由は自由な文)。v2 のカードには `mesh_limits {max_nodes_ns, max_nodes_euler}` (正の整数) があり、
  承認の記録の保護する条件に入る (`approval.py:309`、無ければ既定 `MAX_NODES_2D`)。今は `mesh_limits` を誰も照合しない
  (`approval.py:14` に「照合は §5.1 #8」)。
- recipe `axismach.wt_card_chain/v1` の格子の段 (`campaign/recipes/wt_card_chain.py:286`) は `hand_mesh.py build` と `check` を回し、
  その後 `approval.check` (カード・問題 YAML の承認との照合; 上限は見ない) を通す。recipe v1 は case/71 専用 (手組みの格子だけ)。
- V1a で、今のコードの準備は V6 (run_0014) と初期壁・格子の座標が違うと分かった (x_E 40.633257 vs 40.630797; 候補は #6e eabbab97、未確定)。
  したがって今の hand_mesh で作る格子は V6 の格子と座標で一致しない見込み。

## 2. 期待値と出典

- plan §4.6: 上限の承認はカードのユーザ欄に入れ、承認元に結び付ける。case/71 の承認は case/71 のカードだけに効く。Euler と NS の両方を照合。
  規則は mesh-size-guard §4.3 (上限を超えるときはユーザの了承を取ってから書く、AI が了承なしに書かない)。
- §6 V2a: 手組みの格子を包んだ段が case/71 の節点の座標と一致する。V2b: 雛形の格子で、別のカード (case/70 M5 の空気) が品質の関門を通る。
  V2c: 別のカード + 過去の承認の文で上限を広げようとすると格子を作る前に止まる。承認した上限を超える作り直しも止まる。case/71 のカードと承認では V2a と同じ格子。

## 3. 設計案

(1) **承認から上限を取る**: `approval.mesh_limits_for_card(approval, card)` — 承認の記録を検査して読み (sha256)、カードが承認の原本か、承認の記録を
指す試行の入力かを確かめ (`approval.check` と同じ規則)、`protected.mesh_limits` の `max_nodes_ns`・`max_nodes_euler` と承認元
(承認の記録の sha256・承認の user_ref・カードの sha256) を返す。結び付かなければ例外。

(2) **手組みの格子** (`hand_mesh.py build --approval <card_approval.json>`): 数えた節点数 (NS・Euler) を承認の上限と格子を作る前に比べ、
超えたら問題 YAML を書かずに止める (承認した上限を超える作り直しも同じ)。書く `max_nodes` は承認の上限そのもの、`max_nodes_reason` は
承認元 (記録のファイル名・sha256・user_ref)。自動の拡大と固定の文は消す。`--approval` 無しのときは入力の問題 YAML の上限 (`node_budget`;
明示が無ければ既定 30 万) を超えたら止める (自動で広げない)。recipe の格子の段は `--approval` を渡す。

(3) **問題 YAML の照合**: `approval.check` (問題 YAML を渡したとき) に、問題 YAML の `mesh`・`mesh_euler` の実効の上限 (`node_budget`; 明示が
無ければ既定) が承認の `max_nodes_ns`・`max_nodes_euler` 以下であることを足す (超えれば MISMATCH)。recipe の問題・格子の段と
判定の入口 (`approval.check_run`) で効く。

(4) **v2 のカード**: `overrides` の `mesh.max_nodes`・`mesh.max_nodes_reason`・`mesh_euler.max_nodes`・`mesh_euler.max_nodes_reason` を
誤りにする (上限は `mesh_limits` に書く)。生成器 (`wt_card`) は v2 のカードの `mesh_limits` を問題 YAML の `max_nodes` と理由
「カードの mesh_limits (承認の対象)」に写す (`mesh_limits` が無ければ既定のまま)。v1 のカードは今のまま。

(5) **雛形の格子の段の契約** (`forge_design/campaign/mesh_stage.py`): `check_template(problem, out_dir, limits)` — 生成器の問題 YAML の
`mesh`・`mesh_euler` から、`prepare_ns`・`prepare` と同じ構築 (設計の連鎖・積分法の初期壁・物理壁・`mesh_params`・`generate_axisym_mesh`・
msh・変換・`check_mesh_quality`) で NS と Euler の格子を作り (forge は起動しない)、記録 {問題 YAML の sha256、節点数と承認の上限、
格子の不変部分の sha256 (座標と接続)、品質の VERDICT}。上限は作る前に照合。recipe v1 には入れない (case/71 専用のまま; 雛形の recipe は後)。

(6) **試験** (`design/tests/run_mesh_limits_tests.py`、GPU なし): 下の事前登録。

## 4. 事前登録の案

- **V2a (手組みの格子を包んだ段の A/B)**: 同じ入力 (case/71 `card_m55_3coef.yaml` を版 2 にして `mesh_limits` {443000, 467000}
  (出所: case/71 のユーザ決定「a」2026-10-15)、問題 `problem_m55_p3_card.yaml`) で、A = #8 の前の `hand_mesh.py` (HEAD 0161e51c の版を一時に写す)
  の build、B = #8 の後の `hand_mesh.py build --approval`。期待: 問題 YAML の `mesh`・`mesh_euler` は `max_nodes_reason` 以外の全キーが一致
  (`max_nodes` も 443000・467000 で一致)、それぞれの問題 YAML から作った NS・Euler の格子の節点の座標 (`generate_axisym_mesh` の戻り) が
  ビットで一致。参考 (合否にしない): case/71 の `problem_m55_p3.yaml` (10-15 の版) の mesh 節との差を記録する (V1a の差の候補 #6e の影響を見る)。
  V6 の格子と比べない理由: V1a で今のコードの初期壁が V6 と違うと分かっており、包む前後の差と混ざるため。
- **V2b**: case/70 `card_m5_air_demo.yaml` を版 2 にして (mesh_limits なし = 既定 30 万) 生成 → `mesh_stage.check_template` → NS・Euler とも
  品質 PASS・節点数が上限以下。
- **V2c (負例)**:
  (a) 版 2 のカード (case/70) の `overrides` に `mesh.max_nodes: 443000`・`mesh.max_nodes_reason: "ユーザ決定「a」2026-10-15: …"` →
  カードの検査で誤り (生成しない)。
  (b) case/71 の版 2 のカードで `mesh_limits` を書かない (既定 30 万) 承認 → `hand_mesh.py build --approval` が問題 YAML を書かずに止まる。
  (c) `mesh_limits` {400000, 467000} の承認 (NS が要る 442,520 より小さい) → 止まる (承認した上限を超える作り直し)。
  (d) 承認の上限 30 万の承認に、`max_nodes: 443000` と過去の理由を書いた問題 YAML を `approval.check` → MISMATCH。
  (e) 別のカード (case/70) の承認を case/71 のカードの build に渡す → 承認とカードが結び付かず止まる。
  正例: case/71 のカードと {443000, 467000} の承認では V2a の B と同じ問題 YAML。
