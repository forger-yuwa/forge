# 収束を「残差の準定常・流量の収支・目的の量の準定常」で判定し、全場を出し続けずに監視する

## メタ

- **area**: `tooling` (判定の道具・停止規則・生成器の既定) と `architecture` (forge の出力の追加; 数値には触らない)
- **status**: `draft` (2026-10-15 起票)
- **related_docs**:
  - [`procedures/solver-settings.md`](../../procedures/solver-settings.md) の output の節 (新しい出力のキー)
  - [`AGENTS.md`](../../AGENTS.md) の「収束確認」(移行の後に書き換える対象)
  - `methods/design/overview.md` の停止規則・壁の熱流の節 (`iface_Qf_eff`)
- **related_plans**:
  - [`tooling-nozzle-card-local-mesh.md`](tooling-nozzle-card-local-mesh.md) — 起点 (§4.4 と §5.1 #3q・#3r・#3v・#3x・#3y・#4a・#4b をこの plan へ移した)。M5.5 の壁の設計はあちらに残る
  - [`tooling-convergence-and-wall-resolution-gates.md`](tooling-convergence-and-wall-resolution-gates.md) — 今の `check_convergence.py` と判定区間 (stage_manifest) の正本。この plan はその後継の判定を作り、区間の管理を引き継ぐ
  - `tooling-generator-defaults-and-launch-gate.md` (作業ツリー `forge-gate`) — 生成器の既定の全体。ここでは要求カードの雛形の NS の本段・出力・判定の既定だけを扱い、あちらの段取りと整合させる
- **created**: `2026-10-15`
- **owner**: `M6ノズルの設計` セッション (作業ツリー `forge-cardmesh` / `feature/card-local-mesh`)

## 1. 目的

収束の判定を「残差の低下の桁数」から「残差が準定常・流量の収支が取れる・目的の物理量が準定常」の 3 条件に替え (ユーザ 2026-10-15)、その判定に要る量を全場 (1 枚 約 100 MB) でなく小さい出力 (節点の状態・境界の質量流束) から取れるようにする。判定の基準はケースごとに書けて、評価の対象が無ければ判定不能を返す (勝手な台本で判定しない)。要求カードのノズルでは、本段をライン陰解法 B0 にし、全場は粗い間隔と止めた step だけにする。

## 2. スコープ

- **やる**: forge の出力の追加 (全場の別の間隔・協調停止・節点の状態の出力・境界の質量流束の積算; 出力だけ)、判定の共通の部品と道具 `check_settled.py`、基準のファイル `settle_criteria.yaml` とプリセット、停止規則 v2、要求カードの雛形の既定 (出力・本段 B0・リミッタの基準値の記録と継承)、呼び出し側 (run_case.sh・受け渡しの関門・設計チェーン) の移行、AGENTS.md と Stop フックの書き換え (試した後)。
- **やらない**: 数値の計算の変更。非定常の問題・周期の問題の判定 (未設定なら判定不能を返すだけ; 必要になったら量の種類を足す)。δ_E の時系列 (境界層の厚みの時系列は判定から外す、ユーザ了承)。ノズル以外のケースのプリセット (移行の範囲は §4.7)。

## 3. 前提 (観測事実)

- 今の `check_convergence.py` は系列の最大からの低下の桁数 (既定 3 桁) で判定する。収束した前の解から始めた run (case/71 `run_0010_ns_b0`) は始めの残差が小さく 2.3〜3.0 桁で NOT CONVERGED、`--from-floor` も過渡の山で落ちる (`case/71.m55_combustion_wt/_runs/run_0010_ns_b0/judge/`)。
- 停止規則 (`design/forge_design/evaluate/stop_rule.py`) は全場を毎出力読んで r (流量の収支)・θ_r・出口コア M を測り、止めた後に δ の 2 窓で全場 21 枚から δ* を抽出する (1 時間強; 抽出は高速化済み 15.9 s / 1 回 (12 コア))。止め方は SIGTERM (`stop_rule.py:784`)。
- 境界の h5 (`output.cpp:333`、`outputHDFflg`) の中身は `bc.bvar` = 境界条件の値で、流れの値ではない (出口の `Ps` は課した背圧; `output.cpp:425`、`boundaryCond_d.cu:570`)。音速・局所の組成も無い。壁の出力に `iface_Qf_eff` は無い (`output.interfaceDiag: 1` のときだけ作る、適用は定常・等温壁など; `conjugateWall.cpp:206`・`:265`)。
- point probe (`probe/point_probes.cu`) は重心の最近傍で T・P・U だけ。
- 本段 B0 (lineImplicit・方向別の擬似 dt・cap 50・itj 5・cfl 4・緩和 0.7) は M5.5 で point (cfl 1.0) より約 2.7 倍速く床に着き、同じ解に向かう (case/71 run_0007・run_0008、card plan §6 V6)。リミッタの基準値は forge が最初の段の開始の場の体積平均で自動に決め、基準長は格子の境界箱 (`main.cpp:1768`)。log は 10 桁。
- 試作 `solver_density_cuda/tools/check_settled.py` (残差だけ) は、codex の再現で、短い履歴・欠測・列が 1 つ・末尾で 10⁸ 倍の跳ねを SETTLED にする (下の §6.1 の plan 段レビュー M1・M2)。

## 4. 設計方針

ユーザの決定 (2026-10-15): 「残さが準定常になり、流量バランスが取れ、目的とする物理量が準定常担ったら収束と判定」「収束判定基準はケースによって変えられるようにしてほしい」「評価対象を設定していないと収束判定が不明と返す…あほなaiが勝手にスクリプトを作り出してスタックしないように」「出口境界と入り口境界と軸中心のh5を出し続けて…結果丸ごと出し続けて容量圧迫されるのがな」「境界層厚さの時系列チェックはできなくなる…最終だけあればね」「line implicitを規定化」。codex 諮問 2 回 (card plan の §6.1: m55-three-coef・m55-settle-and-monitor) と plan 段レビュー (card plan の plan-2) の採用をすべて入れた現行の仕様。

### 4.1 forge の出力 (出力だけ; キー無しは今とビット同一)

1. `time.fieldOutStepInterval` (省略時 = `outStepInterval`): 全場 `res_<step>.h5` だけに効く。境界の h5・probe・下の 3・4 は `outStepInterval`。最後の step は今どおり全場を書く。
2. **協調停止**: run dir の `STOP_REQUEST` を毎 step 確かめ、見つけたらその step の全場・境界の h5・節点の出力・流量の行・残差の行を書く。**全部の書き込みの成功を確かめた後だけ** `STOPPED_AT` (step・理由・書いたファイルの一覧と大きさ) を書いて終了コード 0。書き込み・通知の失敗は終了コード非 0 (監視側は REACHED にしない)。
3. **節点の状態の出力**: `output.nodeProbeFile` (節点の番号の表) で、`res_nodes_<step>.h5` に MESH/COORD と VALUE/ro・Ux・Uy・Uz・P・T・sonic (全場と同じ値・型)。
4. **境界ごとの質量流束の積算**: `output.boundaryMassflux: 1` で `boundary_massflux.csv` (step と境界ごとの Σ massflux、実寸 [kg/s] と生の和、流出が正)。値はソルバが流束の計算で使った境界の面の質量流束。
5. **壁の熱流** (等温壁のプリセットのとき): 生成器が `output.interfaceDiag: 1` を有効にし、壁の h5 の `iface_Qf_eff` を節点ごとに使う。適用できない壁 (断熱など) はプリセットの解決のときに明示して外す (黙って外さない)。

### 4.2 判定の基準のファイル `settle_criteria.yaml` (run の中に 1 つ)

- 判定の道具の入力は run の `settle_criteria.yaml` だけ。設計チェーンの run は問題 YAML の `evaluate.settle` (`preset` + 足し引き・閾値の上書き) を準備のときに解決して写し、出典と版を記録する。手で組む run は同じ形式を置く (ひな形を用意)。
- 中身: 判定区間 (本段の開始の step と stage_manifest の区間の識別; 区間の管理は gates plan の stage_manifest を引き継ぎ、本段の開始のときに記録する)、観測の間隔、窓 W と 2W、残差の**必要な列**と**非活動を許す列**、各量の種類・入力のファイル・対応の ID (physID・節点の表)・単位・閾値・適用条件。
- 未設定・誤用は合格にしない: ファイルが無い・目的の量が 1 つも無い → `INDETERMINATE (評価対象が未設定)`・終了コード非 0・置き場所と量の種類の一覧を出す。知らない種類・キー、閾値の欠けは例外。入力のデータが無い・点が足りない・step が単調でない・欠測が窓に入る → 判定不能。
- 量の種類 (最初の版): `residual`、`mass_balance`、`exit_core_mach`、`line_mach` (r/r_w と x の範囲で決まる線; 最大偏差・片振幅・位置ごとの変化)、`exit_flow_angle`、`wall_heat_nodes`。新しい量は道具に種類を足す (試験つき)。
- プリセット `wt_nozzle_v1` (要求カードの風洞ノズル): 上の 6 つ (等温壁でないなら `wall_heat_nodes` を外したと明示)。

### 4.3 判定 (共通の部品 `settle.py` と道具 `check_settled.py`)

- **残差**: 必要な列ごとに log10 を取り、末尾の窓 W (既定 2 万 step) と 2W の被覆・欠測・step の単調を先に検査。窓は端点込みで `min_points` 以上。傾き |s|·W ≤ log10(1.05)、振幅 p95 − p5 ≤ 0.1 桁、直前の窓からの持続的な上昇 ≤ log10(1.10)。**生の系列の末尾の跳ね・短い反復のバースト**を別の検査で落とす (閾値は基準のファイルに; 単発の急増・反復のバースト・持続的な上昇を負例にする)。非有限は DIVERGED。`OSCILLATING` 相当は不合格。
- **流量の収支**: `boundary_massflux.csv` から r = (ṁ_in − ṁ_out)/ṁ_in、窓の平均の絶対値 ≤ 5e−4・振れ ≤ 5e−3、ṁ_in・ṁ_out の変化 ≤ ṁ の 1e−4 (案)。
- **目的の量**: `res_nodes_*.h5` から量を作る関数は既存の抽出 (`nozzle_report.eta_line`・停止規則の出口コア M) と**演算の順まで同じ**にする (M を節点で作ってから補間する順に固定)。止まっている条件は設計の公差から決めた固定の尺度 (平均で割らない): 出口コア M の窓での変化 ≤ 公差の 1/10、線の最大偏差 ≤ 0.01 %・片振幅 ≤ 0.005 % の変化、線の位置ごとの差 ≤ 0.01 %、流れ角 ≤ 0.01°、壁の熱流は節点ごとの相対変化 (案)。値は試して決める (§6 V3)。
- 3 条件の AND で `SETTLED`。結果のファイル (`SETTLE_VERDICT.json`) に、基準の中身と sha256・起動の識別・格子・判定区間・終端の step・使ったデータの一覧を書く。

### 4.4 停止規則 v2

同じ部品で毎出力判定し、3 条件が保持 (当面 4 万 step; B0 での短縮は新しい条件で再生して検証した後に範囲を限って) 続いたら `STOP_REQUEST` を置く → forge が同じ step を書き終えて `STOPPED_AT` → その step を判定の対象に固定して再評価 (要求の後に条件が崩れたら REACHED にしない) → 止めた step の全場から δ_E を 1 回抽出。θ_r・δ の時系列はやめる。

### 4.5 要求カードの雛形の既定 (NS)

入口・出口・壁・軸の境界の出力と節点の出力と流量の積算を有効、`outStepInterval` 2,000・`fieldOutStepInterval` 20,000、本段 B0 (前の段 S1・S2 は point、Euler の cfl は別に保持)、`evaluate.settle: {preset: wt_nozzle_v1}`。B0 が NS の本段だけに効くことを試験する。

### 4.6 リミッタの基準値

最初の本段の開始時に forge の自動の値 (ρ・P・a と基準長 L_ref) を決め、forge が機械可読・往復できる精度で書いた値を run の記録と stage_manifest に残す。同じ設計の系列 (再開・小さな壁の変更) は継承して固定、新しい運転条件・物理問題は初回の本段で自動。`main_line.limiter_ref` は省略を「自動で決めて記録」の意味にし、検証 (`runner_axismach.py:1533`) を直す。

### 4.7 移行

共通の部品を先に作り、オンラインの監視 (停止規則)・終了後の判定 (`check_settled.py`)・`run_case.sh` の後処理 (`run_case.sh:66` は今の判定を呼ぶ)・受け渡しの関門 (`deltastar_loop.py:653`)・設計チェーンを同じ部品へ移す。旧のデータ・未設定のケース・入口と出口の無い周期の問題・非定常の問題は `INDETERMINATE (対象外)` を返し、今の `check_convergence.py` を残して使う (判定の正本は 2 つのまま、どちらを使ったかを結果に書く)。ノズルで検証が済んでも、全ケースの移行完了とはしない。AGENTS.md の「収束確認」と Stop フックは、ノズルでの検証の後に「設定のある run は `check_settled.py`、無い run は今の道具」と書き換える。

## 5. 実装ステップ

1. forge の出力 (#1): `solver_density_cuda/` の input/solverConfig、output、main.cpp。
2. 判定の部品・道具・基準のファイル (#2〜#4): `solver_density_cuda/tools/settle.py`・`check_settled.py`、`design/forge_design/evaluate/settle_presets.py` (案)。
3. 停止規則 v2 (#5): `design/forge_design/evaluate/stop_rule.py`。
4. 雛形の既定・リミッタの記録 (#6): `design/forge_design/intake/recipes/axismach_*`・`runner_axismach.py`。
5. 移行 (#7): `run_case.sh`・`deltastar_loop.py`・AGENTS.md・`hook_convergence_gate.py`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | forge の出力の追加 (§4.1 の 1〜4) | card plan #4a から移した (implementer が実装中、plan 段レビュー M4・M5 の要件を追加で伝えた)。合格: §6 V1 | O |
| 2 | 判定の部品 `settle.py` と `settle_criteria.yaml` の形式・プリセット `wt_nozzle_v1` | §4.2・§4.3。合格: §6 V2 (負例) | O |
| 3 | 量の種類 6 つと演算の順の固定 | §4.3。合格: §6 V2・V3 | O |
| 4 | 試作の残差の判定を直す | 試作 `check_settled.py` を §4.3 の契約に (被覆・欠測・必要な列・跳ねの検査)。合格: §6 V2 | O |
| 5 | 停止規則 v2 と協調停止の監視側 | §4.4。合格: §6 V4 | O |
| 6 | 雛形の既定・`interfaceDiag`・リミッタの基準値の記録と継承 | §4.5・§4.6。合格: §6 V5 | O |
| 7 | 呼び出し側の移行・AGENTS.md・Stop フック | §4.7。V1〜V5 の後。他のセッションに効くので調整板に書く | F |

## 6. 検証 (事前登録 2026-10-15)

- **V1 forge の出力**: (a) キー無しで変更前のバイナリと全場・境界の h5・残差がビット同一。(b) `fieldOutStepInterval` で全場が指定の間隔と最後だけ。(c) 節点の出力が同時の全場の同じ節点とビット同一。(d) 流量の積算が、同じ流束評価の時点の面ごとの massflux を独立に集計した値と一致 (符号・2π・軸端の有限な双対面・境界の角・周期面の除外を確かめる; 全場の折れ線の積分との差は診断として報告)。(e) 協調停止: 出力の step でない step の停止要求・書き込み失敗・通知失敗の各場合で、成功のときだけ `STOPPED_AT` と終了コード 0、失敗は非 0。(f) 出力を有効にしても全場の値がビット同一 (出力が数値を変えない)。
- **V2 判定の負例 (すべて SETTLED にならない)**: 11 点だけの短い履歴、窓に大きな欠測、必要な列が欠ける (`rms_ro` だけ)、末尾の 1 点だけ 10⁸ 倍、短い反復のバースト、持続的な上昇、活動中の乱流の列が全ゼロ、基準のファイルが無い、目的の量が無い、知らない種類・キー、入力のファイルが無い、総熱量は一定で局所の分布が変わる壁、位置が動くが最大偏差と片振幅が一定の線。正例: 2,000 ごと 0〜40,000 の一定系列は SETTLED (端点込みの点の数を合わせる)。
- **V3 再生**: 既存の run の系列 (case/71 run_0002 (未整定)・0005・0007・0008・0010 の残差と、停止規則の記録・全場から作った出口コア M・r・線の時系列) で、新しい判定が run_0002 を途中のどこでも合格させないこと、run_0008・0010 で最初に満たす step と旧の停止規則の到達 (5 万・4.6 万) の関係を記録 (閾値はこの結果で決める; 決めた後は動かさない)。
- **V4 停止規則 v2**: 小さいノズルの run で、協調停止で止めた step の全場・節点・流量・残差が揃う、オンラインの判定と終了後の再生の判定が同じ、要求の後に条件が崩れたら REACHED にしない。
- **V5 雛形の既定**: B0 が NS の本段だけに効く (前の段・Euler は今のまま)、ρ・P・a・L_ref の往復の精度と継承 (再開で同じ値)、出力の容量が 20 万 step で 1 GB 未満 (全場 10 枚 + 最後と小さい出力)。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-15 | [2026-10-15-tooling-nozzle-card-local-mesh-plan-2.md](../../notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan-2.md) (起点の card plan の §4.4 に対するもの) | NO-GO, C0/M7/m0 | 全件**採用**してこの plan に書き直した: M1 被覆・欠測・必要な列・step の単調の検査 (§4.2・§4.3、V2)、M2 跳ねとバーストの検査 (§4.3、V2)、M3 量ごとの出力の設定と `interfaceDiag` (§4.1 の 5・§4.2、V2)、M4 協調停止の書き込みの成功の確認 (§4.1 の 2、V1 (e); 実装中の #1 へ伝えた)、M5 流量の正解は面ごとの独立の集計・演算の順の固定 (§4.1 の 4・§4.3、V1 (d)・V3; #1 へ伝えた)、M6 共通の部品と移行の範囲・区間の管理の引き継ぎ (§4.2・§4.7)、M7 現行の仕様を 1 つに (この plan)。この plan で plan 段レビューをやり直してから #2〜#7 |

## 7. 影響範囲

forge の出力 (キー無しは今とビット同一)、判定の道具 (新しい道具を足し、今の道具は残す)、停止規則、要求カードの雛形の既定、AGENTS.md の収束確認と Stop フック (最後; 他のセッションに効く)。

## 8. 完了条件

V1〜V5 がすべて合格、§5.1 #7 の移行 (呼び出し側が共通の部品を使う、AGENTS.md と Stop フックの書き換え) が済み、result 段のレビューを通す。ノズル以外のケースのプリセットはこの plan の完了条件に入れない (対象外を返すことだけを確かめる)。

## 9. 変更ログ

- `2026-10-15` — 起票。card plan (`tooling-nozzle-card-local-mesh.md`) の §4.4 と §5.1 #3q・#3r・#3v・#3x・#3y・#4a・#4b、codex 諮問 2 回と plan 段レビュー (NO-GO、Major 7、全件採用) をまとめて現行の仕様に書き直した。
