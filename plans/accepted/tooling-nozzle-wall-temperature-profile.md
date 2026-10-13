# 風洞ノズルの設計チェーンに壁温分布 $T_w(x)$ を入れる (要求カード・問題 YAML・NS・積分法・成果物)

## メタ

- **area**: `tooling / boundary (design chain)`
- **status**: `done`  <!-- 2026-10-13 起票 (ユーザ「壁面温度分布って入れられる？入れられるようにしてほしいし、設計開始時に尋ねるようにしてほしい。実際にデモしてみてほしい」) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「壁の熱境界条件 (断熱 / 等温) と壁温影響の評価」節・「実寸の要求カードと成果物一式」節
  - [`methods/boundary.md`](../../methods/boundary.md) 「壁温分布の入力 (`wallProfile`)」(ソルバ側は実装・検証済み)
  - [`procedures/inlet-profile.md`](../../procedures/inlet-profile.md) (CSV の書式: 座標列 1 つなら 1D 線形補間・範囲外は端値)
- **related_plans**:
  - [`tooling-nozzle-isothermal-wall-chain.md`](../active/tooling-nozzle-isothermal-wall-chain.md) (`spec.wall_thermal` を単一ソースにした §4.1。§5.1 #6 の `wallProfile` はソルバ側として CHT plan へ移管済み)
  - [`boundary-conjugate-heat-transfer.md`](../active/boundary-conjugate-heat-transfer.md) (`ints: {wallProfile: 1}` の実装と検証 case/48 `run_0015_wallprofile`)
  - [`../accepted/tooling-nozzle-wt-deliverables.md`](tooling-nozzle-wt-deliverables.md) (要求カード `forge.wt_card/v1`・成果物一式。本計画はその入力と出力に壁温分布を足す)
- **created**: `2026-10-13`
- **owner**: `sano`

## 1. 目的

設計チェーンの壁温は一様の $T_w$ しか受けない (`probdef.Problem.wall_thermal`)。実機のノズルは冷却の配置でスロート付近が熱く、
入口と出口が冷たいことが多い。完了時の状態:

1. 要求カードと問題 YAML で $T_w(x)$ の表を与えられ、NS (ソルバの `wallProfile`) と積分法の初期壁 (`thermal_bc.Tw_table`) の両方が同じ表を読む (単一ソース)。
2. 設計開始時のインテークで「断熱 / 等温 (一様) / 等温 (分布)」を必ず尋ねる (skill `design-intake`・`nozzle-design`)。
3. 成果物 (Fluent の境界条件) に壁温分布が載る (Fluent の profile ファイル)。
4. case/70 の M5 デモで、壁温分布の NS が回り、**場の壁ノードの温度が与えた分布と一致する**ことを確かめる。

## 2. スコープ

- **やる**: `spec.wall_thermal` の分布の書式と検査、bcond の `ints: {wallProfile: 1}` と `wall_profile_<physID>.csv` の書き出し、積分法の `Tw_table`、
  要求カードのキー、格子の見積もり (局所の $T_w$)、成果物 (Fluent profile・境界条件表・Excel の条件)、skill 2 本、M5 デモ。
- **やらない**: 壁温を熱流束から決める連成 (CHT、`cht_loop`)、時間変化する壁温、`rerun_conditions.py` の分布対応 (今は一様の `--Tw` だけ; 後続)、
  SERN (`runner_sern.py`) の分布対応 (後続)、CONTUR の熱閉包の改良 (isothermal plan §5.1 #10〜#12 の領分)。

## 3. 関連 docs と前提 (観測事実, 2026-10-13)

- ソルバ: `wall_isothermal` に `ints: {wallProfile: 1}` を付けると `wall_profile_<physID>.csv` (座標列 + `Ts`) を壁ノードの座標で補間して
  per-face の `Ts` を埋める (`solver_density_cuda/boundaryCond.cpp:471`、`methods/boundary.md` 582 行〜)。case/48 `run_0015_wallprofile` で
  $T_w=300+100x$ が場の壁ノードの `VALUE/T` に最大差 0.068 K で入ることを確認済み。座標列 1 つなら 1D 線形補間・範囲外は端値。
- 温度を拘束する壁どうしが CV を共有して $T_w$ が食い違うと起動時に拒否 (`conjugateWall::checkWallTemperatureSharing`)。ノズルは壁が 1 本なので当たらない見込み。
- 積分法: `deltastar_integral.integral_bl(thermal_bc={"mode": "prescribed_temperature", "Tw_table": [[x_rt, Tw], ...]})` に対応済み (`deltastar_integral.py:132-140`)。
- 設計チェーン: `probdef.Problem.wall_thermal` は `{mode: adiabatic}` か `{mode: isothermal, Tw}` だけ (`probdef.py:60-74`)。bcond の壁の 1 行
  (`wall_bcond_line`)・積分法の `thermal_bc` (`wall_thermal_bc_integral`)・帳簿がここを読む。
- 要求カード (`intake/wt_card.py`): `wall: {thermal: isothermal, Tw}` だけ。格子の見積もり (Eckert の参照温度) は一様の $T_w$。

## 4. 設計方針

### 4.1 問題 YAML の書式

```yaml
spec:
  wall_thermal:
    mode: isothermal
    Tw_profile: {x_m: [-0.5, -0.1, 0.0, 0.1, 0.4, 1.2], Tw_K: [320, 380, 500, 420, 350, 330]}
```

- 座標は 2 通り (どちらか 1 つ):
  - **`s` (区間で正規化; 長さのパラスタ向けの既定の推奨)** (ユーザ 2026-10-13「長さをパラスタすることもあるから、温度分布は軸方向長さをスケールする必要ありそう」):
    入口直管 [入口端 → 縮流部の始まり] を s ∈ [−2, −1]、縮流部 [縮流部の始まり → 設計スロート] を [−1, 0]、拡大部 [設計スロート → 物理出口 (壁の終点)] を [0, 1]
    に区間ごとに線形で写す。スロートは常に s = 0、出口は常に s = 1。区切りの点はその設計で実際にできた壁から取る (出口は MOC の $x_F$ でなく壁の終点)。
    書式 `Tw_profile: {s: [...], Tw_K: [...]}`。
  - **`x_m` (実寸 [m]、原点は設計スロート)**: 要求カード・STEP・寸法ゲートと同じ座標。長さを変えても表は動かない。
  どちらも実寸の x [m] に写してから、積分法へは `x / scale` (実効の $r_t$)、ソルバの CSV へは x [m] で渡す。写しは**呼ぶたびに**その場の寸法で作る。
  「スロート」は設計スロート (設計壁の x = 0)。物理スロート (δ_r を足した壁の最小径) とは数十 µm 違う (case/70: −0.05 mm) が、壁温分布の位置としては無視できる。
- `Tw` (一様)・`Tw_profile` (表)・`Tw_profile_csv` (カードのみ) は厳密に 1 つ。座標列 (`s` か `x_m`) と `Tw_K` は**同じ長さの 1 次元の有限な数値列**、
  座標は狭義単調増加・2 点以上、`Tw_K` は正 (NaN・Inf は拒否)。CSV の相対パスはカードの所在が基準 (codex plan レビュー m5)。
- **壁の範囲を覆うこと** (codex plan レビュー M2): `s` の表は [−2, 1] を覆うこと。`x_m` の表は、**積分法を呼ぶたび** (寸法の逆算の試しの $r_t$ を含む)
  にその場の寸法で積分区間を覆うかを検査し、足りなければ必要な範囲を示して例外 (`np.interp` の端値クランプに頼らない)。最終の物理壁でも検査する。
- 補間は x の 1D 線形 (ソルバと同じ)。積分法も同じ線形補間で読む。**s の表は点を実寸の x に写してから x で線形補間する**ので、区切りの点 (s = −1・0) を
  またぐ区間では s で線形補間した値とは違う (区切りの点を表に含めれば一致する; カードとインテークで「区切りの点 −2・−1・0・1 を表に含める」と案内する)。

### 4.2 単一ソースの配管 (`probdef.Problem`)

- `wall_thermal` は `{"mode": "isothermal", "Tw": <代表値>, "Tw_profile": {...}}` を返す。代表値は表の x = 0 (設計スロート) の値 (bcond の `Ts` の既定値・帳簿用; 分布があるときソルバは CSV で上書きする)。
- **[実装で確定 2026-10-13]** 積分法の熱条件は `forge_design.wall_temperature.integral_thermal_bc(p, d, scale)` (s の表は設計の区切りの点 `d` が要るので Problem の
  メソッドにしない)。property `wall_thermal_bc_integral` は分布の問題では例外 (scale も区切りの点も知らないので; これを直接使う case スクリプト
  `case/70/c2_calibrate_kf.py`・`case/45/calibrate_contur.py`・`c2pin_solve.py` は分布の問題では止まる; §5.1 #7)。**`scale` は積分法を呼ぶその場の実効の $r_t$**
  (寸法の逆算 `solve_rt_throat` の反復中は試しの $r_t$; `deltastar_loop.py:134,185,197`)。`integral_delta_r` (`runner_axismach.py:1161-1172`) は解決した同じ `scale` を
  熱条件と `integral_bl` の両方に渡す。表を `spec.r_throat` で換算してキャッシュしない (codex 諮問 2026-10-13 の採用)。
  (履歴: 当初案は `wall_thermal_bc_integral_at(scale)` を Problem のメソッドにし、property は `spec.r_throat` で換算する互換としていた。実装では上のとおり
  モジュール関数にし、property は分布で例外にした — codex result レビュー m3 で本文の矛盾を整理。)
- 代表値 (bcond の `Ts`・帳簿) は表の座標 0 の値。**s の表は s = 0 を点として含むこと** (実寸の x に写してから補間するので、0 が点でないと代表値と実際のスロートの温度が
  食い違う; codex result レビュー m1 の再現例 s [−2,−1,1]・T [300,400,600] で代表値 500 K・実際 440 K)。**s の表は直管長 > 0** (直管なしは x_m の表; 同 m4)。
- `wall_bcond_line` は分布なら `ints: {wallProfile: 1}` を付ける。
- `prepare_ns` は分布なら run に**実寸に解決した表** `wall_temperature_profile.json` (元の座標の種類・元の表・x [m] の表・区切りの点) と
  `wall_profile_<wallID>.csv` (**空白区切り**の列 `x Ts`、x は m; ソルバは空白で分割して読む `boundaryCond.cpp:264` — カンマ区切りは起動で落ちる) を書き、`prepare_info.json` の `wall_thermal` に両方の sha256 を記録する。**run の中のこの 2 つが実行条件の正本**
  (codex plan レビュー M3): `run_staged_ns` は各段の起動前に CSV の sha256 を照合し、成果物 (`deliverables`) は run の JSON から壁温を読み (外の問題 YAML や代表 `Ts` で
  代替しない) CSV の sha256 を照合する。不一致は止める。
- 一様の $T_w$ と断熱の問題は、生成される config・CSV・`thermal_bc` が**変更前とバイト同一** (既存の run を変えない)。

### 4.3 要求カードとインテーク

- カード: `wall: {thermal: isothermal, Tw_profile: {s|x_m: [...], Tw_K: [...]}}` か `wall: {thermal: isothermal, Tw_profile_csv: <path>}` (列 `s,Tw_K` か `x_m,Tw_K`)。
  `Tw` との併用は拒否。生成器は表をそのまま `spec.wall_thermal.Tw_profile` に写す (s は無次元、x_m は実寸なので、どちらも r_t での割り直しは不要)。
- 格子の見積もり (`estimate_mesh`) は断面ごとに局所の $T_w(x)$ を使う。s の表の区切りの点は、見積もりでも**その設計の壁** (`geometry.L_pipe`・`L_U`・
  設計壁の終点 `chain["wall"].x_e` に r_t を掛けたもの) から取る (prepare_ns の `anchors_from_design` と同じ写し; 区切りの点を記録に残す)。
  **[2026-10-13 訂正]** 当初は「カードの長さから取る (差は長さの許容内)」としたが、デモのカード (物理スロート径から r_t を解く) で設計壁の終点が
  カードの拡大部長より +0.089 mm 長くなり、カードの長さで写した s の表が壁を覆わず止まった (§4.1 の「区切りの点はその設計で実際にできた壁から」と矛盾していた)。
- skill `design-intake`: 壁の熱条件を必ず聞く — 断熱 / 等温 (一様) / 等温 (分布: x [m] と T_w [K] の表、原点は設計スロート)。
  分布は「ノズル全長 (入口直管の始点〜出口) を覆う」ことを伝える。skill `nozzle-design` の対応表に行を足す。

### 4.4 成果物

- Fluent: `fluent/wall_temperature.prof` (Fluent の point profile 書式) を、**同梱する Fluent メッシュの壁面中心**で原表を評価して書く
  (Fluent の point profile の既定は境界面への最近傍の割り当てなので、壁ノードで書くと面の値が線形分布からずれる; codex plan レビュー M1 の実測で 362/1599 面が 0.1 K 超)。境界条件表の壁の行を
  「Thermal: Temperature = profile `wall-temperature` `temperature`」にする。**Fluent での読み込みは未確認** (Fluent が無い)。
- Excel: `conditions.csv` に壁温の表を載せ、CFD 壁面量表の `T_w_K` (壁ダンプの `Ts`) で分布が入ったことが見える。
- **[実装で確定 2026-10-13]** 成果物は run の JSON (実寸の表) を正本とし、`--problem` の YAML の分布が run と違っても run を使う (一様の YAML を渡したときだけ注記)。
  V3 の「原表」は run の JSON の実寸の表 (補助として元の s の表を独立に写した値とも照合)。prepare_ns 以外で作った run (正本の記録が無いのに bcond に `wallProfile` がある) は
  書き出しで止まる (安全側; §5.1 #7 の後続で `rerun_conditions.py` を対応させるときに正本を書く)。

## 5. 実装ステップ

1. `probdef.py` (書式・検査・`wall_thermal_bc_integral` の r_t・`wall_bcond_line`)、`runner_axismach.prepare_ns` (CSV と記録)、`deltastar_loop` の積分法の呼び出し (r_t を渡す)。
2. `intake/wt_card.py` (カードのキー・格子の見積もり)、`export/deliverables.py`・`export/fluent_msh.py` (profile ファイル・境界条件表)。
3. 試験、skill・methods の同期、M5 デモ。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~plan 段の諮問とレビュー~~ **完了 2026-10-13**: 諮問 (実効の r_t の契約) とレビュー (GO-with-changes, M4/m2) を全件採用 (§6.1)。判断: 2026-10-13 codex | F |
| 2 | ~~配管~~ **完了 2026-10-13** | `wall_temperature.py` (新規)・`probdef.py`・`runner_axismach.py`・`runner_wt.py`。試験 `run_wall_temperature_tests.py` 100 ok: V1b 最大差 5.7e-14 K (A・B、s・x_m; 旧半径の誤実装なら 22.4 K、B だけ範囲外は例外)、V1d の負例 6 種で止まる (`run_staged_ns` は forge の起動前)、V4 は case/45 生産 (断熱)・case/70 c2b (一様) の prepare-only の生成物が 138e83b2 とバイト同一。既存の試験 5 本も同数で合格 | O |
| 3 | ~~カード・見積もり・成果物~~ **完了 2026-10-13** | カード (`wt_card.py`: `Tw_profile`・`Tw_profile_csv`、見積もりの局所 T_w、s の区切りの点は設計壁から) と成果物 (`export/fluent_msh.py`・`deliverables.py`・`wall_tables.py`: Fluent の point profile を Fluent メッシュの壁面中心で、境界条件表・conditions.csv・xlsx シート `wall_temperature`・README、書き出しの入口で正本の照合)。試験: `run_wt_card_profile_tests.py` 87 ok (V1c は期待値と相対差 0、誤実装 (a) 代表温度・(b) 無次元の x・(c) カードの長さの区切り を検出)、`run_deliverables_walltemp_tests.py` 36 ok (V3: 1799 点 = 壁面数、最近傍の値と run の表の相対差 0 [元の s の表を独立に写した値とは 6.4e-11]、壁ノードで書く誤実装なら 585 面が 0.1 K 超; V1d 成果物側の負例 3 種で出力を作る前に止まる)。一様・断熱は 5a6dd76f と同一 (時刻の値を除く)。既存の試験 8 本も同数で合格 | O |
| 4 | ~~skill~~ **完了 2026-10-13** | `design-intake` 3b で壁の熱条件を必ず聞く (断熱 / 一様 / 分布 s か x_m、区切りの点を表に含める)、`nozzle-design` の対応表に壁温 (分布) の行 | O |
| 5 | ~~M5 デモ~~ **完了 2026-10-13** | `case/70.m5_air_wt_demo/run_0011_ns_walltemp_demo`: **V2 合格** (壁 1800 節点の `VALUE/T` と表の線形補間の最大差 1.5e-5 K、許容 0.1 K)。積分法が `Tw_table` を読んだことは `delta_r_initial.json` の `thermal_bc` で確認 (x は実効の r_t 0.0249699 m で割った値、終点 48.0615 = 1.20009 m / r_t)。warm start の元 run_0010 (一様 350 K)、cross-mesh、段 full (soft・mid が壁温の変化を吸収)、品質 PASS、NaN 0。収束は合否と別: `check_convergence --segment` NOT CONVERGED (main)、`check_quasisteady machmax,pmax` ALL STEADY、壁解像 PASS。成果物 PASS (Fluent profile 1799 点) | O |
| 6 | ~~result 段~~ **完了 2026-10-13** | 諮問 (物理壁への反映の確認が不足 [Major] → 採用、V1e 判定 A)。レビュー GO-with-changes C0/M0/m4 → 全件採用して対応 (§6.1)。判断: 2026-10-13 codex — V1e を足せば完了の主張を支える | F |
| 7 | (後続) `rerun_conditions.py`・SERN・較正スクリプトの分布対応 | 本計画の範囲外。`rerun_conditions.py` (一様の `--Tw` だけ)、SERN (`runner_sern.py` は `profile_ok=False`)、`wall_thermal_bc_integral` を直接使う case スクリプト (`case/70/c2_calibrate_kf.py`・`case/45/calibrate_contur.py`・`c2pin_solve.py`) は分布の問題では止まる (安全側)。分布で C2 の較正をするときに `integral_thermal_bc(p, d, scale)` へ移す | O |
| 8 | (ユーザ) Fluent での profile の読み込み | `fluent/wall_temperature.prof` を Fluent で読み、壁の Thermal で `wall-temperature` `temperature` を選べること・面の値が表と合うことの確認 (手元に Fluent が無い; deliverables plan §5.1 #10 と同じくユーザに依頼) | — |

## 6. 検証 (事前登録)

- **V1b 実効の r_t の A/B** (codex 諮問 2026-10-13; 負例は codex plan レビュー M2): 問題と表を固定し、`integral_delta_r` の `scale` だけを A = `spec.r_throat`・B = 1.2 倍にする (forge 不要)。
  積分法が実際に使った温度 `res_init["Tw"]` (`deltastar_integral.py:223`) と、原表を実寸座標 `scale × res_init["x"]` で評価した温度の最大差が **A・B とも ≤ 1e-8 K**
  (表は内部で 10 K 以上の差が出る勾配のものを使う)。`x_m` の表で、B だけ積分区間が表の外に出るように作った場合は**例外**になる (端値クランプで通らない)。
- **V1 単体**: 書式の正例 (表・CSV)・負例 (非単調・負の温度・`Tw` との併用・壁を覆わない表)、`thermal_bc` の `Tw_table` の x が `x_m / r_t` と相対 1e-15、
  bcond の行に `wallProfile: 1`、CSV が表と一致 (値の文字列の往復で相対 1e-12)。
- **V1e δ_r が物理壁まで伝わる** (codex result 段の諮問 2026-10-13 で追加; 事前登録は `case/70.m5_air_wt_demo/walltemp_wall_ab.py` の docstring):
  腕 A = デモの表、腕 B = スロートより下流 (s = 0.2・0.5・1) だけ +100 K。設計壁・scale・補正の設定は同じで、`--init-integral --prepare-only` を各 1 回。
  拡大部 x ∈ [1, x_e − 1] r_t で、積分法の δ_r の差 Δr_ref と保存された物理壁 (`wall_repr.json`) の半径の差 Δr_out を比べる。
  ε = max(1e-8 m, Σ (物理壁の許容 tol.r + 旧表現との差 d0) × scale)。**判定 A**: max|Δr_ref| > 10 ε かつ max|Δr_out − Δr_ref| ≤ ε。
  **結果 (2026-10-13): 判定 A** — max|Δr_ref| 0.155 mm (検査区間の終点 x_e − 1 = 1.175 m; 実際の出口 1.200 m では Δr_ref = Δr_out = 0.159 mm、参考)、
  max|Δr_out − Δr_ref| 7.1e-17 m、ε 1e-8 m (信号/ε 1.6e4)。**[2026-10-13 訂正、codex result レビュー m2]** 当初「出口 0.155 mm」と書いたのは検査区間の終点の値。出力 `case/70.m5_air_wt_demo/_walltemp_wall_ab/`。
- **V1c 格子の見積もり** (codex plan レビュー M4): スロートの温度が同じで上流・下流の温度が違う 2 つの表で、選んだ断面の局所の $T_w$ と丸め前の第 1 層厚が、
  独立に計算した期待値と相対 1e-8。代表温度で置き換える誤実装・無次元の x をそのまま使う誤実装が落ちることを確かめる。
- **V1d 正本の照合** (codex plan レビュー M3): CSV を書き換えた run・別の run の JSON を写した run・スロートの温度が同じで分布だけ違う表の JSON に差し替えた run の
  それぞれで、`run_staged_ns` の起動と成果物の書き出しが止まる (負例)。**外の問題 YAML (`--problem`) の壁温では run の正本を上書きしない** (分布が違っても run の表を使い、
  一様の YAML を渡したときだけ注記; §4.4) — **[2026-10-13 同期、codex result レビュー m3]** 当初の「別の問題 YAML で止まる」は §4.4 の決定と矛盾していた。
- **V2 M5 デモ (case/70) — 温度拘束の配管の検査**: 壁温分布 (例: 入口 320 K → スロート 500 K → 出口 330 K) の要求カードで、生成 → CFD 前の物理壁 (積分法が `Tw_table` を読んだことを
  `prepare_info.json`・`delta_r_initial.json` の記録で確認) → NS (warm start、短い run で可)。**合格: 最終 snapshot の場の壁ノードの `VALUE/T` が与えた分布の線形補間と
  最大差 ≤ 0.1 K** (node の壁は温度をピンするので配管の検査であって流れの妥当性ではない; codex plan レビュー m6)。warm start の元・引き継ぎの方法・壁温を変えた後の
  soft 段・格子品質・序盤と最終の有限性を記録し、収束の VERDICT は合否と分けて必ず記録する (case/48 `run_0015` の 0.068 K の実績に合わせた値; `bvar` の `Ts` の再表示ではなく場で見る)。NaN なし。収束・出口 M はこの検証の対象外 (記録だけ)。
- **V3 成果物**: profile ファイルを自前の読み込みで読み直し、**Fluent メッシュの壁面中心に最近傍で割り当てた値**が、独立に計算した原表の値 (面中心で線形補間) と
  相対 1e-9。境界条件表の壁の行が profile を指す。
- **V4 後方互換**: 一様 $T_w$ と断熱の問題 (case/45 の生産 YAML・case/70 の既存 YAML) の `prepare_ns` の生成物 (bcondConfig・`thermal_bc`) が変更前とバイト同一、既存の試験 (`run_wt_card_tests` ほか成果物の試験 6 本) が全部合格。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-13 | [2026-10-13-tooling-nozzle-wall-temperature-profile-plan.md](../../notes/reviews/2026-10-13-tooling-nozzle-wall-temperature-profile-plan.md) | GO-with-changes, C0/M4/m2 | **全件採用**: M1 Fluent の profile を壁面中心で (§4.4・V3)、M2 範囲を積分法の呼び出しごとに (§4.1・V1b の負例)、M3 run の中の表と CSV を正本に起動・再開・成果物で照合 (§4.2・V1d)、M4 格子の見積もりの数値試験 (V1c)、m5 入力の契約 (§4.1)、m6 V2 を配管の検査と明記し記録項目を追加 |
| result | 2026-10-13 | [2026-10-13-tooling-nozzle-wall-temperature-profile-result.md](../../notes/reviews/2026-10-13-tooling-nozzle-wall-temperature-profile-result.md) | GO-with-changes, C0/M0/m4 | **全件採用**: m1 s の表は s = 0 を点として含むこと (`parse_profile`、回帰試験に codex の再現例) / m2 V1e の値を「検査区間の終点」と「実際の出口 0.159 mm」に分けた (`walltemp_wall_ab.py`・§6) / m3 §4.2 の旧案を履歴として分け、V1d を §4.4 と同期、較正スクリプトと Fluent の実機確認を §5.1 #7・#8 に / m4 s の表は直管長 > 0 をカードで早期に拒否 (x_m を案内、回帰試験)。試験: `run_wall_temperature_tests` 102 ok・`run_wt_card_profile_tests` 89 ok・`run_wt_card_tests` 86 ok・`run_deliverables_walltemp_tests` 36 ok |
| diagnose (result 段の解釈) | 2026-10-13 | [2026-10-13-wall-temperature-profile-result-diagnose.md](../../notes/reviews/2026-10-13-wall-temperature-profile-result-diagnose.md) | Major 1: 検証が積分法の温度の読み取りで止まり、δ_r が初期の物理壁に反映されないことを検出できない | **採用**: §6 V1e (壁温だけを変えた prepare-only の A/B) を追加し、判定 A (差 7.1e-17 m ≤ ε 1e-8 m、信号 0.155 mm)。「V2 の NS を延ばして埋めない」も採用 |
| diagnose | 2026-10-13 | [2026-10-13-wall-temperature-profile-diagnose.md](../../notes/reviews/2026-10-13-wall-temperature-profile-diagnose.md) | 結論: 壁温表の無次元化には寸法逆算中の実効 `scale` を使う契約と A/B を足す。Major 1 | **採用**: §4.2 に契約 (`wall_thermal_bc_integral_at(scale)`)、§6 に V1b。「V2 の 0.1 K だけで配管全体を合格にしない」も採用 (V1b・V3 と組で判定) |

## 7. 影響範囲

- `probdef.py`・`runner_axismach.py`・`deltastar_loop.py` (設計チェーン本体、分布のときだけ動作が変わる)、`intake/wt_card.py`・`export/*`、skill 2 本、methods。
- ソルバは変えない。一様・断熱の問題は生成物がバイト同一。

## 8. 完了条件

- §6 V1〜V4 が合格、plan 段・result 段の codex レビューの行、skill と methods の同期。

## 9. 変更ログ

- 2026-10-13: 起票。ソルバ (`wallProfile`) と積分法 (`Tw_table`) は既に分布を受けられ、設計チェーンの `spec.wall_thermal` が一様しか受けないことを確認 (§3)。
- 2026-10-13: ユーザ指摘 (長さのパラスタで分布が伸縮すべき・スロートの座標は 0 前提か) を受け、区間で正規化した座標 `s` (直管 [−2,−1]・縮流部 [−1,0]・拡大部 [0,1]) を §4.1 に追加し、長さのパラスタ向けの推奨にした。実寸の `x_m` も残す。
- 2026-10-13: §5.1 #2 実装 (共通モジュールと配管)。CSV は空白区切り `x Ts` (§4.2 を訂正)、熱条件は `integral_thermal_bc(p, d, scale)` (§4.2 を訂正)、s の補間の意味を §4.1 に追記。
- 2026-10-13: §4.3 を訂正 (見積もりの s の区切りの点もその設計の壁から取る)。カードの長さで写すと、公差の内で設計壁の終点がカードの長さを超えたときに表が壁を覆わない (デモで +0.089 mm)。
- 2026-10-13: §5.1 #3 実装 (カード・見積もり・成果物)。成果物の正本の扱いを §4.4 に追記。デモの YAML `case/70.m5_air_wt_demo/problem_m5_air_walltemp.yaml` (r_t 0.0249699 m、設計壁の終点はカードの拡大部長より +0.089 mm、寸法 PASS)。
- 2026-10-13: §5.1 #5 (M5 デモ) — V2 合格 (`run_0011_ns_walltemp_demo`、最大差 1.5e-5 K)。成果物にも profile が載った。
- 2026-10-13: result 段の諮問を採用し V1e (δ_r → 物理壁) を追加・判定 A。
- 2026-10-13: result 段のレビュー (GO-with-changes, m4) を全件採用して対応し、status done → accepted へ移す。
