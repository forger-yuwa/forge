# case/45 の Euler で、スロート付近の全温が Tt を数百 K 超える — 発生する段・原因・設計の入力への影響

## メタ

- **area**: `verification`
- **status**: `draft`  <!-- 診断の段階。次の手の登録が揃ったら plan 段レビュー (§5.1 #4) を回して in_progress に -->
- **related_docs**:
  - `procedures/solver-settings.md` (Euler・node・SLAU・block-DPLUR の設定)
  - `solver_density_cuda/tools/total_quantities.py` (全温・全圧の求め方)
- **related_plans**:
  - [discretization-moc-axis-limit-and-corrector.md](discretization-moc-axis-limit-and-corrector.md) (V5・V5b・V5c でこの異常が見つかった。MOC の生産採用の判断は本 plan の切り分けを待つ)
  - [tooling-nozzle-throat-monotone-r2.md](../accepted/tooling-nozzle-throat-monotone-r2.md) (単調壁の採用の E′ も同じ種類の Euler の比較)
- **created**: `2026-10-07`
- **owner**: `Claude (Opus 5.5)`

## 1. 目的

case/45 の生産 Euler 格子 G1 の run (非粘性・すべり壁・断熱) で、スロート付近の壁際の数千〜1 万節点の全温が、入口の Tt = 1600 K を最大 +373 K 超える。V5c (MOC の plan) で、これが出力・後処理の換算ではなく保存された状態そのものにあることを確かめた。
本 plan は、この超過が**どの段で生じるか**、**原因は何か**、**設計チェーンの入力 (CFD でピン止めする初期線 run_0062・出口較正) にどの程度効くか**を切り分ける。MOC の変更の採否は扱わない (MOC の plan の側で、本 plan の結果を待つ)。

## 2. スコープ

- **やる**: case/45 の Euler の全エンタルピーの異常の、発生する段の特定・原因の切り分け・設計の入力への影響の見積もり。
- **やらない**: 全温のクリップ、出口較正の調整、再延長による救済 (原因が不明なまま入力や判定対象を変えない)。MOC の採否の判断。NS の場の評価 (影響の範囲の調査として NS の保存物を見ることは、必要になったら別に登録する)。

## 3. 関連 docs と前提

- **観測** (MOC の plan §9、2026-10-07):
  - 全温の最大 − Tt: run_0114 (旧壁の Euler の収束場、番号写像の IC の元) +241 K、腕 A run_0140 +334 K、腕 B run_0143 +334 K、腕 M run_0150 +373 K (いずれも本段 step 18000)。
  - ISEN (run_0153、等エントロピー IC から soft 3000 step → 本段 18000): 本段の初期場 (`res_0.h5` = soft 段の後) で最大 1915 K (軸上、x = 45.9 r_t)、本段の最後で +75 K、延長 (run_0160) で +343 K (壁から 7 層目、x = 0.39 r_t)。最大の位置は時間とともに移る。
  - 保存された P・T・h0 からの全温と、保存量から独立に復元した全温は、全 56 標本で 0.65 K 以内で一致 (V5c、`case/45.isobutane_m6_d155/_band_ab/moc_v5c_thermo_ab.json`)。
- **設定** (run_0150 の `solverConfig.yaml`): Euler、node、本段は 2 次 (convMethod 1)・limiter 2・SLAU・block-DPLUR (timeIntegration 11)・cfl 2・implicitRelax 0.7・nStepInner 5、soft 段は 1 次・cfl 0.5。TP (MIXDRY・H2O)、壁はすべり壁 (`slip`)。
- **諮問の注意** (`notes/reviews/2026-10-07-moc-v5c-next-step-diagnose.md`): 異なる run の最大値を並べても同じ節点の時間変化ではない。全温の超過から直ちにエネルギーの非保存や設計の誤差量を断定しない。すべり壁では SLAU の `wall_flag` 分岐 (`slauWallNormalChi`) は効かない (`mesh.cpp:979` は `wall`・`wall_isothermal` だけに立つ) ので、その切り替えは判別操作にならない。
- **IC の生成**: `design/forge_design/evaluate/ic.py:86`・`:97` が TP の内部エネルギーと運動エネルギーから `roe` を作る。段の引き継ぎ: `runner_axismach.py:885` が soft 段の最終場を入力に移し、`:886` で旧出力を消し、`:898` で本段を起動する。

## 4. 設計方針

切り分けは安い順に、1 回に 1 つの要因だけを変える。各段の A/B の登録 (§6) を結果を見る前に書き、結果を見たら次の手を上位に諮る。

- **E1 (最初の一手、CFD 0 step)**: ISEN の起動前の prep と soft 段の後の状態を、同じ保存量の復元 (V5c の経路 B、`moc_v5c_thermo_ab.path_b`) で比べ、超過が IC の生成で生じるのか、起動後に生じるのかを判別する。

## 5. 実装ステップ

1. `case/45.isobutane_m6_d155/euler_t0_stage_ab.py` (新規): E1 の比較。V5c の `path_b` (保存量 → float64 の原始量と全温) を流用し、ファイルには書かない。
2. E1 を AWS の保存物で回す (forge は起動しない)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~§4・§6 の諮問~~ 完了 | 判断: 2026-10-07 codex (diagnose) `notes/reviews/2026-10-07-moc-v5c-next-step-diagnose.md` — 別 plan に切り出し、E1 (起動前の prep 対 soft 段の後、0 step) を先に | F |
| 2 | ~~E1~~ 完了 (2026-10-07) | 判定は §9 | O |
| 3 | ~~E1 の解釈と次の手~~ 完了 | 判断: 2026-10-07 codex (diagnose) — 半径方向の配点だけを変える E2 を新規に回す。既存の比較は交絡している | F |
| 3b | E2 | §6 E2 の起動スクリプト・評価器 (soft 段の出力の退避を含む) を書き、AWS で回す。合格条件: §6 E2 の判定を出す | O |
| 3c | E2 の解釈と次の手 | 上位に諮る | F |
| 4 | codex plan 段レビュー | 次の手の登録が揃った時点で `codex_review.py <本 plan> --stage plan` | O |
| 5 | codex result 段レビュー | `--stage result` | O |

## 6. 検証

### 6.0 事前登録

- **E1 起動前と soft 段の後の比較** (2026-10-07、結果を見る前。諮問の判別 A/B):
  - 設定は変えず、比べる状態の取得時点だけを変える。追加の計算は 0 step。
  - A = 起動前: AWS の `case/45.isobutane_m6_d155/_prep_moc_v5_mocG1_isen/nozzle.h5` (`IC_MAP.json` の生成後の sha256 と照合する)。
  - B = soft 段の後: `case/45.isobutane_m6_d155/run_0153_euler_icdep_mocG1_isen/nozzle.h5`。保存量が本段の `res_0.h5` と一致するかも確かめ、起動の処理で差が入る場合は区別して記録する。
  - 両者の `ro`・`roU*`・`roe`・`roY*` を、同じ解決済みの熱物性・エンタルピーの基準で float64 に復元する (Euler なので k は足さない)。見る量: max\\|T₀ − 1600\\|、+1 K と +100 K を超える節点の数、最大の位置と壁からの層、h₀ − h_mix(1600 K, Y)。
  - 判定:
    - A が全節点で ±1 K 以内、かつ B が既報 (本段の res_0 の最大 1915 K = +315 K) を再現 → 「数百 K の超過を IC の生成だけで説明する」説を退け、起動後の処理 (soft 段の計算・段の引き継ぎ・本段の起動) を対象にする。
    - A にも +100 K 以上がある → 「起動後に初めて生じた」説を退け、IC の整合を先に扱う。
    - その中間 → 寄与の大きさを記録し、単独の原因とは認定しない。
  - この比較では、延長中の壁近傍の増大までは説明できたとしない。

- **E2 半径方向の配点への感度 (2026-10-07 登録、結果を見る前。諮問 `notes/reviews/2026-10-07-euler-t0-e1-next-diagnose.md`)**:
  - 変えるのは `mesh.wall_first_frac` だけ: A = 1.3e-5 (G1 と同じ)、B = 0.005 (1100 × 65 の Euler と同じ)。等比の配点は半径方向の全体を変える (nj = 97 で軸側の第 1 間隔は 0.088758 r_w → 0.018745 r_w) ので、**判別するのは「半径方向の配点への依存」であって、壁際の float32 の誤差とは即断しない**。
  - 固定: E1 の run_0153 と同じ壁 (`problem_d155_euler_pin_G1_recal_mono_moc.yaml`)・2000 × 97・軸方向の節点・BC・熱物性・MOC・出口較正・solver と変換器のバイナリ。両側とも同じ等エントロピー IC の生成手順で、起動前に全節点で \|T₀ − 1600\| ≤ 1 K を確かめる。既存の G1 の異常な場を片側だけに引き継がない。
  - 長さ: soft 3000 step (1 次、cfl_pseudo 0.5) + 本段 18000 step (2 次、cfl_pseudo 2)。既存の ISEN と対応させた固定の予算で、収束の予測ではない。soft 段の出力も残す (今の `runner_axismach.py:886` は段の終わりに `res_*` を消すので、その前に退避する)。
  - run: run_0161_euler_t0cluster_wff1p3em5 (A)、run_0162_euler_t0cluster_wff5em3 (B)。
  - 見る量: 同じ物理座標・η = r/r_w の領域 (軸域・コア・壁域。境界は結果を見る前に固定する: 軸域 η < 0.1、壁域 η > 0.9、その間をコア) ごとに、\|T₀ − 1600\| の最大・99 % 点、±1 K と ±100 K を超える節点の割合、h₀ − h_mix(1600 K, Y)、最大の節点の位置。壁からの同じ層番号だけでは比べない。soft 段の後と本段の 1000 step ごとに記録する。
  - 判定の前提: 本段だけの区間の `check_convergence` と、本段 6000〜18000 の全 13 枚と末尾 5 枚で、上の領域別の量の `check_quasisteady` を記録する。未収束の短い過渡の差を原因の判定に使わない。
  - 判定:
    - A で > 100 K の異常が再現し、B は全域で ±1 K 以内、かつ判定の前提を満たす → 「総点数だけが原因」(H3) を退け、半径方向の配点への依存を支持する (float32 の幾何の誤差の確定にはしない)。
    - 両側で異常が残り、領域別の異常量の差が 1 K 以内、かつ判定の前提を満たす → 「配点が支配的」を退ける。
    - A が再現しない・未収束・未整定・中間的な改善 → 判別不能。窓の変更や自動の延長で救済しない。
  - 実効のメッシュの軸側・壁側の間隔、全入力とバイナリのハッシュ、run のパスを記録する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-07` | [`notes/reviews/2026-10-07-moc-v5c-next-step-diagnose.md`](../../notes/reviews/2026-10-07-moc-v5c-next-step-diagnose.md) | M5 | 全件採用 (MOC の plan §6.1 にも記録)。本 plan の起票と E1 の登録 |
| diagnose (E1 の次の手) | `2026-10-07` | [`notes/reviews/2026-10-07-euler-t0-e1-next-diagnose.md`](../../notes/reviews/2026-10-07-euler-t0-e1-next-diagnose.md) (ブリーフ [`briefs/2026-10-07-euler-t0-e1-next.md`](../../notes/reviews/briefs/2026-10-07-euler-t0-e1-next.md)) | M4 | 採用: 既存の比較は交絡している (run_0120 は Pt 0.8 倍・壁 interp・別の IC で、solverConfig の一致だけでは対照にならない) → 条件をそろえた E2 を新規に (M) / wall_first_frac は半径方向の全体を変えるので「配点への感度」として登録 (M) / AR の PASS は双対の幾何の精度の証明ではない (M) / 出口較正と G1 の Euler の比較による判断は保証を保留、較正値は診断中は固定 (M)。倍精度のビルド・全温のクリップ・`Md_moc_offset` の調整はしない |

## 7. 影響範囲

- 診断スクリプトの追加だけ。ソルバ・設計のコードは変えない (原因が分かって修正が要る場合は、改めて plan を書く)。

## 8. 完了条件

- [ ] 発生する段と原因 (または原因の候補と未確定の範囲) を §9 に記録する
- [ ] 設計チェーンの入力 (初期線 run_0062・出口較正) への影響の見積もりを記録する
- [ ] codex レビュー 2 回 (plan / result) を §6.1 に記録
- [ ] `status: done` にして accepted へ移動し、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-07` — 起票 (MOC の plan の V5c と諮問 `notes/reviews/2026-10-07-moc-v5c-next-step-diagnose.md` を受けて)。E1 を登録。
- `2026-10-07` — **E1 の判定: 「数百 K の超過を IC の生成だけで説明する」説を退ける** → 起動後の処理 (soft 段の計算・段の引き継ぎ・本段の起動) を対象にする。スクリプト `case/45.isobutane_m6_d155/euler_t0_stage_ab.py` (commit f89f2835)、出力 `_band_ab/euler_t0_stage_ab.json` (AWS で実行)。
  - A (起動前の prep、等エントロピー IC。sha256 は IC_MAP.json の生成後の値と一致): 全温 1599.99〜1600.00 K、max\|T₀ − Tt\| 0.011 K。
  - B (soft 段の後、1 次・cfl 0.5・3000 step): 全温 1245.18〜1915.38 K、max\|T₀ − Tt\| 354.8 K (x = 52.4 r_t、軸上)。+1 K を超える節点 95,983、−1 K を下回る節点 91,475、+100 K を超える節点 707。超過だけでなく、Tt を下回る側にも同じ程度に広がる。
  - C (本段の res_0): 全温の統計は B と同じ。保存量は ro・roU* がビット一致、roe・roY0・roY1 は一致しない (本段の起動時に forge が何かを変える。差の大きさは記録していない)。
  - 解釈 (どの処理で、なぜ生じるか) は確定していない。次の手は上位に諮る (§5.1 #3)。
- `2026-10-07` — **ユーザ決定「1」(Euler の全温の異常の切り分けを続ける)**: 設計チェーンの入力 (初期線 run_0062・出口較正) に関わるので、影響の大きさを先に見積もる。MOC の生産採用は本 plan の結果を待つ。上流の多項式化の plan は並行して進める。- `2026-10-07` — **格子による違いの観測** (主セッション、別ツリー `/home/sano/work/forge/case/45.isobutane_m6_d155/` の保存物を `total_quantities.py --Tt 1600` で): 1100 × 65・wall_first_frac 0.005 の Euler (run_0047・0059〜0064・0086 は cfl 6、run_0120 は cfl 2) は全温の最大 − Tt が +0.07〜0.12 K で Tt + 1 K を超える節点 0。2000 × 97・wall_first_frac 1.3e-5 の G1 (run_0114、cfl 2) は +241 K・3,490 節点。run_0120 と run_0114 の `solverConfig.yaml` は cfl・step 数・出力間隔以外に差が無い。**CFD でピン止めする初期線の元 run_0062 は全温が正常 (+0.073 K)**。出口較正は G1 の run_0113+0114 由来。原因は確定していない (次の手を諮問中)。
- `2026-10-07` — **codex (diagnose) に諮った (E1 の次の手)**: `notes/reviews/2026-10-07-euler-t0-e1-next-diagnose.md` — 同じ壁・2000 × 97・同じ初期化で wall_first_frac だけを変える E2 を新規に回す。**訂正**: 上の「格子による違いの観測」で run_0120 を対照のように扱ったのは不適切だった。run_0120 は Pt を 0.8 倍にした試験で、壁は interp、IC の元も run_0114 と違う (一致していたのは solverConfig だけ)。観測 (1100 × 65 の Euler はすべて正常、G1 は異常) は事実だが、格子以外の条件と交絡している。また E1 の soft 段の後のずれ (全節点の約 96.6 % が ±1 K を外れ、最大は軸上の低温側) と、後期の壁際の高温側の超過を、同じ原因とする根拠はまだ無い。出口較正 (run_0113+0114 由来) は異常な場に依存していることまでは確認したが、出口 M の誤差の量・符号は未確定。較正値は診断中は固定する。
