# 風洞ノズルの設計チェーンに壁温分布 $T_w(x)$ を入れる (要求カード・問題 YAML・NS・積分法・成果物)

## メタ

- **area**: `tooling / boundary (design chain)`
- **status**: `in_progress`  <!-- 2026-10-13 起票 (ユーザ「壁面温度分布って入れられる？入れられるようにしてほしいし、設計開始時に尋ねるようにしてほしい。実際にデモしてみてほしい」) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「壁の熱境界条件 (断熱 / 等温) と壁温影響の評価」節・「実寸の要求カードと成果物一式」節
  - [`methods/boundary.md`](../../methods/boundary.md) 「壁温分布の入力 (`wallProfile`)」(ソルバ側は実装・検証済み)
  - [`procedures/inlet-profile.md`](../../procedures/inlet-profile.md) (CSV の書式: 座標列 1 つなら 1D 線形補間・範囲外は端値)
- **related_plans**:
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) (`spec.wall_thermal` を単一ソースにした §4.1。§5.1 #6 の `wallProfile` はソルバ側として CHT plan へ移管済み)
  - [`boundary-conjugate-heat-transfer.md`](boundary-conjugate-heat-transfer.md) (`ints: {wallProfile: 1}` の実装と検証 case/48 `run_0015_wallprofile`)
  - [`../accepted/tooling-nozzle-wt-deliverables.md`](../accepted/tooling-nozzle-wt-deliverables.md) (要求カード `forge.wt_card/v1`・成果物一式。本計画はその入力と出力に壁温分布を足す)
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
  `case/70/c2_calibrate_kf.py`・`case/45/calibrate_contur.py`・`c2pin_solve.py` は分布の問題では止まる)。当初案: `wall_thermal_bc_integral_at(scale)` は分布なら `{"mode": "prescribed_temperature", "Tw_table": [[x_m / scale, Tw], ...]}`。**`scale` は積分法を呼ぶその場の実効の $r_t$**
  (寸法の逆算 `solve_rt_throat` の反復中は試しの $r_t$; `deltastar_loop.py:134,185,197`)。`integral_delta_r` (`runner_axismach.py:1161-1172`) は解決した同じ `scale` を
  熱条件と `integral_bl` の両方に渡す。表を `spec.r_throat` で換算してキャッシュしない (codex 諮問 2026-10-13 の採用)。既存の property `wall_thermal_bc_integral` は
  `spec.r_throat` での換算 (case スクリプト向けの互換)。
- `wall_bcond_line` は分布なら `ints: {wallProfile: 1}` を付ける。
- `prepare_ns` は分布なら run に**実寸に解決した表** `wall_temperature_profile.json` (元の座標の種類・元の表・x [m] の表・区切りの点) と
  `wall_profile_<wallID>.csv` (**空白区切り**の列 `x Ts`、x は m; ソルバは空白で分割して読む `boundaryCond.cpp:264` — カンマ区切りは起動で落ちる) を書き、`prepare_info.json` の `wall_thermal` に両方の sha256 を記録する。**run の中のこの 2 つが実行条件の正本**
  (codex plan レビュー M3): `run_staged_ns` は各段の起動前に CSV の sha256 を照合し、成果物 (`deliverables`) は run の JSON から壁温を読み (外の問題 YAML や代表 `Ts` で
  代替しない) CSV の sha256 を照合する。不一致は止める。
- 一様の $T_w$ と断熱の問題は、生成される config・CSV・`thermal_bc` が**変更前とバイト同一** (既存の run を変えない)。

### 4.3 要求カードとインテーク

- カード: `wall: {thermal: isothermal, Tw_profile: {x_m: [...], Tw_K: [...]}}` か `wall: {thermal: isothermal, Tw_profile_csv: <path>}` (列 `x_m,Tw_K`)。
  `Tw` との併用は拒否。生成器は表をそのまま `spec.wall_thermal.Tw_profile` に写す (実寸なので割り直し不要)。
- 格子の見積もり (`estimate_mesh`) は断面ごとに局所の $T_w(x)$ を使う。
- skill `design-intake`: 壁の熱条件を必ず聞く — 断熱 / 等温 (一様) / 等温 (分布: x [m] と T_w [K] の表、原点は設計スロート)。
  分布は「ノズル全長 (入口直管の始点〜出口) を覆う」ことを伝える。skill `nozzle-design` の対応表に行を足す。

### 4.4 成果物

- Fluent: `fluent/wall_temperature.prof` (Fluent の point profile 書式) を、**同梱する Fluent メッシュの壁面中心**で原表を評価して書く
  (Fluent の point profile の既定は境界面への最近傍の割り当てなので、壁ノードで書くと面の値が線形分布からずれる; codex plan レビュー M1 の実測で 362/1599 面が 0.1 K 超)。境界条件表の壁の行を
  「Thermal: Temperature = profile `wall-temperature` `temperature`」にする。**Fluent での読み込みは未確認** (Fluent が無い)。
- Excel: `conditions.csv` に壁温の表を載せ、CFD 壁面量表の `T_w_K` (壁ダンプの `Ts`) で分布が入ったことが見える。

## 5. 実装ステップ

1. `probdef.py` (書式・検査・`wall_thermal_bc_integral` の r_t・`wall_bcond_line`)、`runner_axismach.prepare_ns` (CSV と記録)、`deltastar_loop` の積分法の呼び出し (r_t を渡す)。
2. `intake/wt_card.py` (カードのキー・格子の見積もり)、`export/deliverables.py`・`export/fluent_msh.py` (profile ファイル・境界条件表)。
3. 試験、skill・methods の同期、M5 デモ。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~plan 段の諮問とレビュー~~ **完了 2026-10-13**: 諮問 (実効の r_t の契約) とレビュー (GO-with-changes, M4/m2) を全件採用 (§6.1)。判断: 2026-10-13 codex | F |
| 2 | ~~配管~~ **完了 2026-10-13** | `wall_temperature.py` (新規)・`probdef.py`・`runner_axismach.py`・`runner_wt.py`。試験 `run_wall_temperature_tests.py` 100 ok: V1b 最大差 5.7e-14 K (A・B、s・x_m; 旧半径の誤実装なら 22.4 K、B だけ範囲外は例外)、V1d の負例 6 種で止まる (`run_staged_ns` は forge の起動前)、V4 は case/45 生産 (断熱)・case/70 c2b (一様) の prepare-only の生成物が 138e83b2 とバイト同一。既存の試験 5 本も同数で合格 | O |
| 3 | カード・見積もり・成果物 | §4.3・§4.4。合格: §6 V1・V3 | O |
| 4 | skill | §4.3。合格: インテークの質問が skill に入る | O |
| 5 | M5 デモ | §6 V2 | O |
| 6 | result 段 | 諮問とレビュー | F |
| 7 | (後続) `rerun_conditions.py`・SERN の分布対応 | 本計画の範囲外 | O |

## 6. 検証 (事前登録)

- **V1b 実効の r_t の A/B** (codex 諮問 2026-10-13; 負例は codex plan レビュー M2): 問題と表を固定し、`integral_delta_r` の `scale` だけを A = `spec.r_throat`・B = 1.2 倍にする (forge 不要)。
  積分法が実際に使った温度 `res_init["Tw"]` (`deltastar_integral.py:223`) と、原表を実寸座標 `scale × res_init["x"]` で評価した温度の最大差が **A・B とも ≤ 1e-8 K**
  (表は内部で 10 K 以上の差が出る勾配のものを使う)。`x_m` の表で、B だけ積分区間が表の外に出るように作った場合は**例外**になる (端値クランプで通らない)。
- **V1 単体**: 書式の正例 (表・CSV)・負例 (非単調・負の温度・`Tw` との併用・壁を覆わない表)、`thermal_bc` の `Tw_table` の x が `x_m / r_t` と相対 1e-15、
  bcond の行に `wallProfile: 1`、CSV が表と一致 (値の文字列の往復で相対 1e-12)。
- **V1c 格子の見積もり** (codex plan レビュー M4): スロートの温度が同じで上流・下流の温度が違う 2 つの表で、選んだ断面の局所の $T_w$ と丸め前の第 1 層厚が、
  独立に計算した期待値と相対 1e-8。代表温度で置き換える誤実装・無次元の x をそのまま使う誤実装が落ちることを確かめる。
- **V1d 正本の照合** (codex plan レビュー M3): CSV を書き換えた run・別の問題 YAML・スロートの温度が同じで分布だけ違う表のそれぞれで、`run_staged_ns` の起動と
  成果物の書き出しが止まる (負例)。
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
