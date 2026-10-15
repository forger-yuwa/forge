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
6. **随時の全場の書き出し** (2026-10-16 ユーザ「めったに出さないとしたら、フルのh5も流し切るまで手にはいらないのかな?それは良くないかも」): run dir の `SNAPSHOT_REQUEST` を協調停止と同じ毎 step の確認で見つけたら、その step の全場・境界の h5・節点の出力・流量の行を**止めずに**書き、書き込みの成功を確かめてから `SNAPSHOT_REQUEST` を消して `SNAPSHOT_DONE` (step・ファイルの一覧) を書き足す。失敗しても run は止めず、`SNAPSHOT_FAILED` を残す。`fieldOutStepInterval` (既定 20,000) の間でも、見たいときに全場を取れるようにする。

### 4.2 判定の基準のファイル `settle_criteria.yaml` (run の中に 1 つ)

- 判定の道具の入力は run の `settle_criteria.yaml` だけ。設計チェーンの run は問題 YAML の `evaluate.settle` (`preset` + 足し引き・閾値の上書き) を準備のときに解決して写し、出典と版を記録する。手で組む run は同じ形式を置く (ひな形を用意)。
- 中身: 判定区間 (本段の開始の step と stage_manifest の区間の識別; 区間の管理は gates plan の stage_manifest を引き継ぎ、本段の開始のときに記録する)、観測の間隔、窓 W と 2W、残差の**必要な列**と**非活動を許す列**、各量の種類・入力のファイル・対応の ID (physID・節点の表)・単位・閾値・適用条件。
- 未設定・誤用は合格にしない: ファイルが無い・目的の量が 1 つも無い → `INDETERMINATE (評価対象が未設定)`・終了コード非 0・置き場所と量の種類の一覧を出す。知らない種類・キー、閾値の欠けは例外。入力のデータが無い・点が足りない・step が単調でない・欠測が窓に入る → 判定不能。
- 量の種類 (最初の版): `residual`、`mass_balance`、`exit_core_mach`、`line_mach` (r/r_w と x の範囲で決まる線; 風洞ノズルの既定は r/r_w = 0.03 (2026-10-16 ユーザ決定、card plan §6 V5c); 最大偏差・片振幅・位置ごとの変化)、`exit_flow_angle`、`wall_heat_nodes`。新しい量は道具に種類を足す (試験つき)。
- プリセット `wt_nozzle_v1` (要求カードの風洞ノズル): 上の 6 つ (等温壁でないなら `wall_heat_nodes` を外したと明示)。

### 4.3 判定 (共通の部品 `settle.py` と道具 `check_settled.py`)

- **残差**: 必要な列ごとに log10 を取り、末尾の窓 W (既定 2 万 step) と 2W の被覆・欠測・step の単調を先に検査。窓は端点込みで `min_points` 以上。傾き |s|·W ≤ log10(1.05)、振幅 p95 − p5 ≤ 0.1 桁、直前の窓からの持続的な上昇 ≤ log10(1.10)。**生の系列の末尾の跳ね・短い反復のバースト**を別の検査で落とす (閾値は基準のファイルに; 単発の急増・反復のバースト・持続的な上昇を負例にする)。非有限は DIVERGED。`OSCILLATING` 相当は不合格。
- **流量の収支**: `boundary_massflux.csv` から r = (ṁ_in − ṁ_out)/ṁ_in、窓の平均の絶対値 ≤ 5e−4・振れ ≤ 5e−3、ṁ_in・ṁ_out の変化 ≤ ṁ の 1e−4 (案)。
- **目的の量**: `res_nodes_*.h5` から量を作る関数は既存の抽出 (`nozzle_report.eta_line`・停止規則の出口コア M) と**演算の順まで同じ**にする (M を節点で作ってから補間する順に固定)。止まっている条件は設計の公差から決めた固定の尺度 (平均で割らない): 出口コア M の窓での変化 ≤ 公差の 1/10、線の最大偏差 ≤ 0.01 %・片振幅 ≤ 0.005 % の変化、線の位置ごとの差 ≤ 0.01 %、流れ角 ≤ 0.01°、壁の熱流は節点ごとの相対変化 (案)。値は試して決める (§6 V3)。
- 3 条件の AND で `SETTLED`。結果のファイル (`SETTLE_VERDICT.json`) に、基準の中身と sha256・起動の識別・格子・判定区間・終端の step・使ったデータの一覧を書く。

### 4.4 停止規則 v2

同じ部品で毎出力判定し、3 条件が保持 (当面 4 万 step; B0 での短縮は新しい条件で再生して検証した後に範囲を限って) 続いたら `STOP_REQUEST` を置く → forge が同じ step を書き終えて `STOPPED_AT` → その step を判定の対象に固定して再評価 (要求の後に条件が崩れたら REACHED にしない) → 止めた step の全場から δ_E を 1 回抽出。θ_r・δ の時系列はやめる。

**実装で決めたこと (2026-10-16、#5)**: v2 は run に `settle_criteria.yaml` があり停止規則を要求するとき (起動前の検査 `preflight_v2` に通らなければ S1 の前に止める; 旧方式の fallback に落とさない)。完成の印は流量の行 s (forge は step s の h5 を全部閉じた後に書く) と残差の outer_end の行 s−1 の改行。最初の判定は開始 + 2W + 間隔。REACHED は「要求の種類が settled・forge の終了コード 0・STOPPED_AT の検査 (step の範囲・residual_step = step − 1・要求の照合の token・ファイルの大きさ) に通る・要求の後の格子点で崩れない・止めた step の再評価が SETTLED」のすべてのとき、それ以外は FAILED (理由の分類つき)、要求を置かずに上限まで回れば CENSORED。**オンラインの判定が DIVERGED、または直らない「未設定」の判定不能 (基準に無い残差の列が値を持つ等) なら、理由を書いた STOP_REQUEST で協調停止して FAILED にする** (上限まで GPU を使い切らない; 主セッションの判断で採用)。kill は要求の後 `stop_wait_s` (既定 900 s) を過ぎた・見張りの例外のときだけ (FAILED)。状態は `STOP_RULE_STATE.json` (method v1/v2)、δ_E の抽出は `stop_rule_v2/delta_E_<step>/` (帯の位置の感度 band_sensitivity.json つき)。

### 4.5 要求カードの雛形の既定 (NS)

入口・出口・壁・軸の境界の出力と節点の出力と流量の積算を有効、`outStepInterval` 2,000・`fieldOutStepInterval` 20,000、本段 B0 (前の段 S1・S2 は point、Euler の cfl は別に保持)、`evaluate.settle: {preset: wt_nozzle_v1}`。B0 が NS の本段だけに効くことを試験する。

### 4.6 リミッタの基準値

最初の本段の開始時に forge の自動の値 (ρ・P・a と基準長 L_ref) を決め、forge が機械可読・往復できる精度で書いた値を run の記録と stage_manifest に残す。同じ設計の系列 (再開・小さな壁の変更) は継承して固定、新しい運転条件・物理問題は初回の本段で自動。`main_line.limiter_ref` は省略を「自動で決めて記録」の意味にし、検証 (`runner_axismach.py:1533`) を直す。

### 4.7 移行

共通の部品を先に作り、オンラインの監視 (停止規則)・終了後の判定 (`check_settled.py`)・`run_case.sh` の後処理 (`run_case.sh:66` は今の判定を呼ぶ)・受け渡しの関門 (`deltastar_loop.py:653`)・設計チェーンを同じ部品へ移す。旧のデータ・未設定のケース・入口と出口の無い周期の問題・非定常の問題は `INDETERMINATE (対象外)` を返し、今の `check_convergence.py` を残して使う (判定の正本は 2 つのまま、どちらを使ったかを結果に書く)。ノズルで検証が済んでも、全ケースの移行完了とはしない。AGENTS.md の「収束確認」と Stop フックは、ノズルでの検証の後に「設定のある run は `check_settled.py`、無い run は今の道具」と書き換える。

### 4.8 plan 段レビュー (2026-10-15、GO-with-changes) で足した契約

1. **判定の格子と止めた step** (M1): W・2W・保持は**通常の出力の格子** (`outStepInterval` の倍数の step) の点だけで判定する (端点込みの点の数もこの格子で数える)。協調停止で書いた非出力の step の観測は、別の検査 (非有限なし・直前の通常の観測からの急変なし・最終値が各量の許容の内) だけに使う。残差の番号の対応を固定する (forge は場を `iStep+1` で書き、残差は `iStep` で記録し、残差は書き込みがバッファされる; `main.cpp:2571`) — 判定は「出力の step s の場」と「残差の行 s−1 まで」を同じ時点として扱い、必要な出力一式 (場・節点・流量・残差の行) がそろったことを各出力の完成の印 (例: 流量の行の後に書く小さな印、または協調停止の `STOPPED_AT` の一覧) で確かめてから判定する。
2. **受け渡しの関門の置き換え** (M3): 新方式 (`settle_criteria.yaml` のある run) では、今の `handoff_gate` (`deltastar_loop.py:1006`) の (a) 停止の状態 → 停止規則 v2 の REACHED、(b) `check_convergence --segment` → 残差の準定常、(c) 流量の保存の検査 → 流量の収支 (数値流束の和)、(e) δ_E・δ_target の 2 窓 → **最後の全場からの δ_E の抽出 1 回 + 帯の位置の感度 (y_b × 1.25 の差を記録; ユーザ了承で時系列は見ない)** に置き換える。(d) 局所 y₁⁺・(f) Euler の参照と ω・(g) 抽出の被覆・(h) 次の壁の健全性は**そのまま残す**。旧方式の run は今の関門のまま (分岐は `settle_criteria.yaml` の有無、結果にどちらを使ったかを書く)。上書きの規則は今と同じ (被覆・参照・次の壁は上書きで通さない)。
3. **リミッタの基準値を区間の鍵に** (M4): forge が確定した ρ・P・a・L_ref (自動で決めたなら決まった値) を stage_manifest の区間比較の hard key に入れる (今は入っていない: `stage_manifest.py:266`、codex の再現で limiterRoRef 1→10 でも区間 1 つ)。値が分からない (自動で記録が無い) 段は連結しない。自動の値を forge が書いた後に本段の記録を確定する順序にする。起動の識別は既存の plan `tooling-stage-manifest-launch-binding.md` の担当、ここは基準値の鍵だけ。
4. **壁の熱流の許容** (M5): 節点 i で |ΔQ_i| ≤ Q_abs,i + ε_rel·|Q_ref,i|、Q_abs,i = q_tol·A_i (熱流束の固定の公差 × 局所の面積)。q_tol・ε_rel はプリセットの値 (試して決める)。符号の反転は許容の内なら可。
5. **容量の条件** (m6): case/71 の規模 (442,520 節点、FP64 のビルド、`output.level` 1、全場 1 枚 約 100.8 MB) で 20 万 step、`fieldOutStepInterval` 20,000 なら全場は res_0・2 万ごと 10 枚・最後で最大 12 枚 (約 1.2 GB)。条件は「全場の枚数 ≤ ⌈上限/間隔⌉ + 2」と「全場以外の出力 (境界・節点・流量; 残差の CSV は除く) ≤ 300 MB (10⁶ bytes 単位)」にする。

## 5. 実装ステップ

1. forge の出力 (#1): `solver_density_cuda/` の input/solverConfig、output、main.cpp。
2. 判定の部品・道具・基準のファイル (#2〜#4): `solver_density_cuda/tools/settle.py`・`check_settled.py`、`design/forge_design/evaluate/settle_presets.py` (案)。
3. 停止規則 v2 (#5): `design/forge_design/evaluate/stop_rule.py`。
4. 雛形の既定・リミッタの記録 (#6): `design/forge_design/intake/recipes/axismach_*`・`runner_axismach.py`。
5. 移行 (#7): `run_case.sh`・`deltastar_loop.py`・AGENTS.md・`hook_convergence_gate.py`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~forge の出力の追加 (§4.1 の 1〜4・6)~~ (2026-10-16 実装済み・commit) | **結果** (implementer、試験は case/71 `run_0009_euler` の写し 466,570 節点・float・手元の GPU、一時 run は scratchpad): (b) 全場は fieldOutStepInterval と最後だけ・境界/節点/流量は outStepInterval ごと。(c) 節点の出力は全場と 15 時点でビット同一 (座標を含む)。(d) 流量の CSV は面ごとの massflux の書き出しを独立に集計した値とビット同一 (入口・出口・壁・軸、3 run)、面の集合は互いに素で和が境界の全 9,814 面、軸は厳密に 0、全場の折れ線の積分との差 ≤ 2.3e−7 (診断)。(e) 停止は出力の格子でない step (1507・311・265) で終了コード 0・全出力が揃う、全場/通知の書き込み失敗は終了コード 3 で STOPPED_AT なし、起動時の不正入力 12 種は終了コード 1。(g) 随時の書き出しは格子でない step (163・401) で書けて計算は続き、失敗は SNAPSHOT_FAILED 1 行で続き、次の要求は成功。**(a)・(f)・(g) のビット同一は判定不能**: 変更前のバイナリどうしが再実行で step 0 からビット一致しない (node の atomicAdd)。代わりに決定的な部分 (res_0・xmf・出力の名前の集合・最初の流束評価の面の massflux とその入力の状態) のビット一致と、600 step の差が変更前どうしの差と同じ桁 (RMS 相対 1e−6〜5e−5) であることを確かめた → **代替の検査で足りるかは上位に諮る (保留)**。未試験: 周期面の除外、RK・dual-time での停止、FP64 のビルド。未対応: 停止・随時の step で CHT の固体の状態を書かない (#1b)、`check_field_regress.py` が res_nodes を境界の出力と取り違える (#1b) | O |
| 1b | #1 の残り: 協調停止・随時の書き出しの step に CHT の固体の状態 (`conjugateWall::writeConjugateState`) を書く、`check_field_regress.py` が `res_nodes_*.h5` を境界の出力と取り違えないようにする、残差履歴を随時の書き出しの step で flush する | 合格: CHT の小さい run で停止の step の壁と固体の状態の step が一致、check_field_regress の既存の試験が通る | O |
| 1c | **V1 (f′)・(g′) の状態不変 A/B** (codex diagnose settle-v1-and-line-height の判別) | 同じ node・GPU・陰解法の状態を、通常の出力の step でも最後の step でもない位置で同期・退避し、A = 要求なし・B = `SNAPSHOT_REQUEST` あり (と、キー無し・キーあり) で**出力の処理だけ**を実行して、次の更新が読む device と host の状態を比べる (数値更新は 0 step)。既存の `pdeSnapshot`/`pdeRestoreMismatch` (`main.cpp:2247` 付近) を使い、今は省いている面の配列の被覆と全場の出力が更新する host の写しを足す。合格: 差 0 (試験した経路の直接の変更を除外; 将来の全軌道のビット一致までは主張しない)。差があれば非決定性だけで説明する説を退け、原因を調べる。#6 の main.cpp の変更が終わってから | O |
| 2 | ~~判定の部品 `settle.py` と `settle_criteria.yaml` の形式・プリセット `wt_nozzle_v1`~~ (2026-10-16 実装済み・commit) | `solver_density_cuda/tools/settle.py`・`settle_criteria.template.yaml`、`design/forge_design/evaluate/settle_presets.py` (`resolve`・`write_settle_criteria`・`probe_node_ids`)。決めたこと: 壁の面積は壁の h5 の COORD・CONNE から (軸対称は 1 rad あたり、`iface_Qf_eff` [W/rad] と単位を揃える); 節点の表は 1 行 1 整数; 線の節点は使う列とその両隣の列の軸・η を挟む対とその両側・壁、出口の列は全部; float のビルドは res_nodes の COORD が float32 なので `structured.mesh_file` (格子の h5) の座標を NODE_ID で引いて照合; required にも inactive_ok にも無い列が値を持てば判定不能; residual・mass_balance・目的の量 1 つ以上がそろわなければ「評価対象が未設定」; 出力の前提 (boundaryMassflux 等) は `preconditions` に記録。**格子点 0 の扱い (2026-10-16 決定)**: forge は本段の開始 (step 0) に res_0 しか書かず残差の行 −1・節点・流量が無いので、格子点 0 は観測なし → 最初に判定できる終端は 2W + 間隔 (2 万・2,000 なら 42,000)。開始の行で代用しない (開始直後の跳ねで上昇の検査が甘くなる)。V2 の正例はこの形に読み替える | O |
| 3 | ~~量の種類 6 つと演算の順の固定~~ (2026-10-16 実装済み・commit) | 試験 `run_settle_presets_tests.py` 37/37: run_0010 res_86000 から抜いた 5,335 節点 (全体の 1.21 %) で、出口コア M (= stop_rule の exit_core_M = nozzle_report の M_mean = 5.502840924572144)・線 (η 0.03 は maxabs 0.062304 %、η 0.1 は 0.094020 %)・位置ごとの d・流れ角 22 節点が全場からの抽出と**ビット同一**。η 0.03 は 747 列すべてで第 1 層の η と差 1e−9 未満、np.interp の対が列で (0,1)・(1,2) に分かれても両側の節点を含めるので一致。プリセットの既定の線の x は (x_E + 2, x_F − 1) r_t = 4.273〜7.900 m → case/71 は V5c と同じ 4.0〜7.95 m を問題 YAML の `x_range_m` で与える (#6) | O |
| 4 | ~~試作の残差の判定を直す~~ (2026-10-16 実装済み・commit) | `check_settled.py` を settle.py の CLI に書き直し (終了コード 0/1/2/3、引数・基準の誤りも 3)。`test_settle.py` 54/54 (V2 の負例・正例、design venv と system python3 の両方)、`test_gate_bad_input.py` PASS。V2 の負例は理由の分類 (データ不足・閾値の超過・未設定) まで検査。残差の部分の再生 (閾値は暫定): run_0005 は NOT_SETTLED (roUx の振幅 0.131 桁・roK の傾き 0.102 桁/窓と振幅 0.361 桁)、42k〜200k のどの格子点でも SETTLED にならない; run_0002 (未整定) は途中のどこでも SETTLED にならない (V3(a) の要件の向き)。run_0008・0010 の残差の履歴は AWS にあるので V3(a) の残りは後 | O |
| 5 | ~~停止規則 v2 と協調停止の監視側~~ (2026-10-16 実装済み; commit は #6 と同時 — `runner_axismach.py` を共有) | `stop_rule.watch_main_stage_v2`・`preflight_v2`、`runner_axismach._stop_rule_v2_for_main` (§4.4 の「実装で決めたこと」)。試験 `run_stop_rule_v2_tests.py` 56/56 (偽の forge): 非出力の step 80,001 で止めて REACHED (M1 の正例、再評価は 80,000 の窓 + 停止の観測の別検査 17 個)、オンラインの判定が終了後の `settle.judge`・`check_settled.py --end-step 80001` と verdict・量の sha256 まで一致、要求の後に崩れる・終了コード 3・STOPPED_AT のずれ・上限・要求を無視・非有限・未設定・残りの制御ファイル・旧方式の分岐。既存の停止規則・runner の試験は合格。δ_E の抽出は run_0010 res_86000 で保存済みの抽出と #3z の band_ab (0.16475590862106337 mm @ 2.2596 m) に**ビット同一** (17.9 s、8 並列)。本物の forge (#1 のバイナリ、run_0009_euler の写し、scratchpad、間隔 10・窓 50・保持 50・閾値は配管の確認用に緩めた): 160 で要求 → 非出力の step 198 で STOPPED_AT → 再評価 SETTLED → REACHED。残り (V4 を AWS の case/71 規模で): 判定 1 回の時間 (残差の CSV を毎回全部読む; 長ければ増分に)、要求から停止までの遅れ、δ_E の所要、stop_wait_s 900 s が FP64 の全場の書き出しに足りるか、容量 (停止の step の全場 1 枚を含む)。v2 の run は今の受け渡しの関門を通らない (#7 まで安全側) | O |
| 6 | ~~雛形の既定・`interfaceDiag`・リミッタの基準値の記録と継承~~ (2026-10-16 実装済み・commit) | forge が起動ごとに `limiter_ref.json` (ρ・P・a・L_ref を %.17g、cfg_fnv・自動の印; 出力だけ) を書く。`stage_manifest.py` (v3) の hard キーに `limiter.ref` (記録 > 4 つの明示 > limiterScaled 0 の明示 > unknown)、unknown の段は連結しない。`main_line.limiter_ref` は省略 = 自動・`{ro,p,a[,L]}`・`{from_run}` (cfg_fnv・stage_manifest と照合)。雛形 (`axismach_final_v1`・`axismach_explore_v1`) に B0 (`main_line`、基準値は自動)・`cfl_euler: 1.0`・`ns_output` (2,000・全場 20,000・4 境界の h5・流量)・`settle: {preset: wt_nozzle_v1}`、出力は本段の設定だけに書く。case/71 の例 `problem_m55_p4_settle.yaml` (p3_b0 + `ns_output` + `settle` の x 4.0〜7.95 m; リミッタは p3_b0 の固定値)。試験: `run_settle_defaults_tests.py` 69 項目・`test_stage_manifest_limiter_ref.py`・既存の設計チェーンの試験 (主セッションで `run_deliverables_walltemp` の xlsx のシートの順の期待を #3u に合わせて直した)。本物の `prepare_ns` で S1・S2 が run_0010 の記録と sha256 一致、PRESTAGE は run_0010 の設定とバイト同一。GPU の短い run: JSON と log の値が一致、`from_run` の継承で %.17g の文字列が一致、区間の鍵 (自動→同じ値の固定は 1 区間、ρ を変えると別、その後の自動も別)、NS 40 step で 4 境界の h5・res_nodes・流量・limiter_ref.json が書け check_settled が 6 量を計算できる。容量: 全場 12 枚 ≤ 12・全場以外 289.6 MB ≤ 300 MB (余裕 3.5 %; 等温壁の界面の診断が壁の h5 の約 45 %、随時の書き出しは数えない)。振る舞いの変化: 自動の基準値で起動し直す段は別の区間になる (Euler の mid と本段も別に; case/49・57・59 の gen_runs のように forge の前に段を記録する台本は自動の段が unknown で連結しない) | O |
| 6b | 系列のリミッタの基準値の自動の継承 (§4.6 の「同じ設計の系列は継承」の実装) | 今は `from_run` を YAML に明示したときだけ継承し、雛形の既定は自動なので、再開・次の壁の pass が起動ごとに違う基準値になる (case/71 の S2 の場からの自動の値は ρ 0.57614、run_0008/0010 の固定値は 0.58100)。決定 (2026-10-16 主セッション、§4.6 の方針の実装): NS の初期場を別の NS run から作るとき (`prepare_ns` の初期場の元・deltastar_loop の前の pass) は、その run の `limiter_ref.json` (computed) か固定値を自動で継承し、出所を記録する。`main_line.limiter_ref` の明示はそれに優先、新しい問題 (元の NS run が無い) は自動。合格: 試験で継承・明示の優先・元の記録が無いときの扱い (自動に戻さず止めて報告) **決定の修正 (2026-10-16、implementer が実装の前に止めた報告を受けて)**: 既存の固定値の run (case/45 の B0 の系譜・case/71 の p2/p3/p4、run_0010) は ρ・P・a の 3 つだけを明示し、L_ref (`limiterRefLength`) は forge が格子の境界箱から自動で決めていた (場によらない) → **(A) を採用**: 元の run の `limiter_ref.json` (cfg_fnv が本段の設定と一致・computed) があれば 4 つを継承、無ければ本段の設定の ρ・P・a の明示を継承し L_ref は新しい run でも格子から自動 (L も明示されていれば 4 つ)。新しい run の `limiter_ref.json` に L の実際の値が残るので区間の鍵は決まる (stage_manifest の「値が分かる」= 4 つとの違いはここに書いて残す)。あわせて: B0 の本段 (`main_line` が mapping) だけに効かせ point には効かせない; 元の run の NS/Euler は solverConfig (`turbulence.model`・`viscMethod`) と prepare_info の `chain` で決め、食い違い・決められないときは ValueError (`limiter_ref: auto` か明示で通る); 元の run の本段の設定は `MAIN_STAGE_STARTED` の `main_config` の sha256 (と #6 以後なら fnv) で確かめる; `limiter_ref.json` が壊れている・computed false は ValueError、cfg_fnv の不一致は固定値に回して理由を記録; stage_manifest は tag `main` の段と照らす; 入口は `prepare_ns(ic_from)` の 1 か所 (deltastar_loop・run_from_table・run_pass0_integral はここに集まる) なので deltastar_loop に配線は足さない | O |
| 7 | 呼び出し側の移行・AGENTS.md・Stop フック | §4.7。V1〜V5 の後。他のセッションに効くので調整板に書く 追加: §4.8 の 2 (受け渡しの関門の置き換え) | F |

## 6. 検証 (事前登録 2026-10-15)

- **V1 forge の出力**: (a) キー無しで変更前のバイナリと全場・境界の h5・残差がビット同一。(b) `fieldOutStepInterval` で全場が指定の間隔と最後だけ。(c) 節点の出力が同時の全場の同じ節点とビット同一。(d) 流量の積算が、同じ流束評価の時点の面ごとの massflux を独立に集計した値と一致 (符号・2π・軸端の有限な双対面・境界の角・周期面の除外を確かめる; 全場の折れ線の積分との差は診断として報告)。(e) 協調停止: 出力の step でない step の停止要求・書き込み失敗・通知失敗の各場合で、成功のときだけ `STOPPED_AT` と終了コード 0、失敗は非 0。(f) 出力を有効にしても全場の値がビット同一 (出力が数値を変えない)。(g) 随時の書き出し (§4.1 の 6): 出力の格子でない step の `SNAPSHOT_REQUEST` で全場ほかが書かれ run は続く・その後の場が要求なしの run とビット同一・書き込み失敗で run が止まらず `SNAPSHOT_FAILED`。
- **V1 の改訂 (2026-10-16、事後の改訂; codex diagnose settle-v1-and-line-height Major 1・2 採用)**: 上の (a)・(f)・(g) のビット同一は、変更前のバイナリ自身が再実行でビット一致しない (node の加算の順序) ことを確かめずに書いた条件で、**元の条件のままでは判定不能**として残す (PASS に書き換えない)。改訂した条件 (2026-10-16 以後に適用): (a′) キー無しでは数値計算の経路・既存の出力の仕様を変えない。再実行で決定的な項目 (res_0・xmf・出力の名前の集合・最初の流束評価の面の massflux とその入力の状態) はビット一致、非決定的な場は `check_field_regress.py` (step 600・factor 2、保存量・化学種・P・T・h0、変更前の 2 本の反復をノイズ床) で判定する。(f′) 出力の処理が、次の数値更新が読む状態 (device と host の両方) を変えないことを、数値更新 0 step の状態不変 A/B で直接確かめる (§5.1 #1c)。別 run の全場のビット一致は求めない。(g′) 随時の書き出しの成功・失敗の動作の条件は元のまま、数値への影響は (f′) と同じ状態不変の検査で確かめ、続けた場の比較は非退行の補助とする。既に得た 600 step の結果は参考資料として新しい条件の結果と区別する (codex が check_field_regress で事後に評価: キー無し 1.48・キーあり 1.66・随時の書き出し 1.14・書き出しの失敗の後 1.01 (ノイズ床比の最大)、いずれも PASS; 各 1 本で統計的な同等性の証明ではない)。
- **V2 判定の負例 (すべて SETTLED にならない)**: 11 点だけの短い履歴、窓に大きな欠測、必要な列が欠ける (`rms_ro` だけ)、末尾の 1 点だけ 10⁸ 倍、短い反復のバースト、持続的な上昇、活動中の乱流の列が全ゼロ、基準のファイルが無い、目的の量が無い、知らない種類・キー、入力のファイルが無い、総熱量は一定で局所の分布が変わる壁、位置が動くが最大偏差と片振幅が一定の線。正例: 2,000 ごと 0〜40,000 の一定系列は SETTLED (端点込みの点の数を合わせる; forge の実際の出力の形では格子点 0 に観測が無いので最初の判定は 42,000、§5.1 #2 の決定)。
- **V3 再生 (M2 で 2 つに分けた)**: (a) **残っているデータでの部分ごとの再生**: 残差 (case/71 run_0002 (未整定)・0005・0007・0008・0010 の履歴) と、停止規則の記録の出口コア M・r (旧の折れ線の積分 = 別定義の診断)・保存した V5c の 2 量の系列で、量ごとに判定する。run_0002 を途中のどこでも合格させないこと、run_0008・0010 で最初に満たす step と旧の停止規則の到達 (5 万・4.6 万) の関係を記録。流量の数値流束・線の位置ごとの変化・壁の熱流は既存の run に無いので、この段では**判定不能**として扱い、負例を落とせた数に数えない。(b) **新しい出力を持つ run での全条件の検証**: 新しい出力 (§4.1) を有効にした run (小さいノズルと case/71 規模を 1 本ずつ) で全条件を判定し、閾値を固定 (固定した後は動かさない)。その後、別の系列 (別の run か同じ run の別の区間) で正例・負例を確かめる。各試験で、不合格の理由が閾値の超過かデータ不足かを結果に書かせ、それも検査する。
- **V4 停止規則 v2**: 小さいノズルの run で、協調停止で止めた step の全場・節点・流量・残差が揃う、オンラインの判定と終了後の再生の判定が同じ、要求の後に条件が崩れたら REACHED にしない。**正例: 一定の系列を非出力の step (例 80,001) で止めて正しく合格する** (M1)。**移行の試験** (M3): V4・V5 の後、新方式の run を受け渡しの関門まで通し、(d)・(f)・(g)・(h) が今どおり効き、旧方式の run は今の関門のままであること。
- **V5 雛形の既定**: B0 が NS の本段だけに効く (前の段・Euler は今のまま)、ρ・P・a・L_ref の往復の精度と継承 (再開で同じ値)、**区間の鍵: 同じ値の継承は同じ区間・どれかの基準値の変更は別の区間・自動の値が分からなければ連結しない** (M4)、出力の容量は §4.8 の 5 の条件。壁の熱流 (M5): 一定のゼロは正例、許容の内のゼロの交差は正例、許容を超える局所の変化は負例。
- **V3(b)・V4 の実施 (事前登録 2026-10-16、codex diagnose `2026-10-16-settle-v3b-v4-preregistration-diagnose-settle-v3b-v4.md` を全件採用して確定)**:
  - 共通: AWS B、FP64 のビルド (`~/forge-cardmesh-fp64` = 91027fde + `flowFormat.hpp` の typedef を double、run_0010 のバイナリと同じ作り方)、`FORGE_CUDA_BLOCKSIZE=128`。問題 `case/71.m55_combustion_wt/problem_m55_p4_settle.yaml` (3 係数の壁 = run_0010 と同じ格子・B0・リミッタの固定値、`ns_output` 2,000・全場 20,000、`settle: {preset: wt_nozzle_v1, x_range_m: [4.0, 7.95]}`、線 η 0.03)。初期場は run_0010 と同じ (run_0008_ns_b0 res_90000 から補間)、Euler の参照 run_0009_euler res_12000。**1 本ずつ回す** (codex Minor 1)。
  - **閾値**: 今の暫定値 (`wt_nozzle_v1` 版 1) に**固定**して検証する。実測に合わせて較正しない (codex M2・M3: 流量の許容を実測に追随させると 1 % の不釣合いでも通る; 壁の熱流の 2 係数は最大 × 1.5 から一意に決まらない)。目的の量の許容は設計の公差から決めた値のまま (M1 の採用部分)。壁の熱流の確定値はユーザの熱設計上の絶対・相対の許容を別に要る (未定)。
  - **(1) `run_0012_ns_free` (V3(b)、先に)**: 停止規則なしで本段 200,000 step。終了後に全格子点で判定を再生し、**停止候補 t_s** = 停止規則 v2 が STOP_REQUEST を置く最初の格子点 (固定の閾値で 3 条件が保持 40,000 続いた点) を求める。**同じ軌道の t_s から t_s + 80,000 まで**で、(a) 出口コア M・線の各位置の d・出口の流れ角 (節点ごと)・壁の局所の熱流 (節点ごと) について **t_s の値からの差の最大**が固定の許容 (settle の各量の tol と同じ; 壁は q_tol·A_i + ε_rel·|Q(t_s)|) の内、(b) 各出力の格子点で全条件 (残差・流量の収支・目的の量) が SETTLED。**結果 A** (全期間で (a)・(b) とも内): この条件・80,000 step の範囲で早すぎる停止の懸念を退ける。**結果 B** (どれか超過): 今の停止の判定を反証 → 閾値を緩めず停止規則 (窓の幅の判定・保持) を見直す。t_s > 120,000 または停止候補が無い → 検証未成立。80,000 step は有限の検証範囲で、永久の整定の保証とは扱わない (codex M1: 窓ごとの幅の判定は 2 万 step あたり許容の 0.73 倍のドリフトを 80 窓すべて SETTLED にし、82k→200k の累積は許容の 4.29 倍になる人工の例)。あわせて、目的の量だけなら通るのに残差か流量で落ちる格子点の数を記録する (残差の条件が独立の情報を足すかの確認; AND に入っているだけでは足したと言わない)。
  - **(1) の判定の道具と細則 (2026-10-16、run_0012 の結果を見る前に固定)**: `solver_density_cuda/tools/settle_virtual_stop.py` (試験 `test_settle_virtual_stop.py` 41/41、両 Python)。t_s は停止規則 v2 の保持の規則 (`stop_rule._hold_step`・`watch_main_stage_v2` の要求の分岐) と同じ (偽の forge の 5 通りで要求の種類・step・保持の経過まで一致)。settled 以外の要求 (diverged・criteria_unset) で止まる場合は停止候補なし → 未成立。**A/B/未成立の順**: 基準なし・停止候補なし・t_s > 120,000 → 未成立 (先に決める) → span の中で正しく読めた格子点に (a) の超過か (b) の NOT_SETTLED・DIVERGED が 1 つでもあれば **B** (データの欠け・判定不能が同時にあっても B、欠けは記録; 超過はそれ自体が反例) → 超過が無く欠け・判定不能があれば未成立 → どれも無ければ A。判定不能: (b) の INDETERMINATE、(a) の値を作れない (データ不足・線の補間に軸の節点が入る・壁の熱流の非有限)。線は位置ごとの d (tol_pos) で判定し、maxabs・halfamp の差は記録だけ。壁の A_i は t_s の壁の h5 から。**追加の記録 (判定には使わない)**: 目的の量だけで保持した場合の停止候補 t_s′ と、t_s′ から span の同じ検証 (残差・流量の条件が早すぎる停止を防いだかの直接の比較; codex の「目的の量だけなら通る不適切な区間を残差が落とせるか」)。
  - **(2) `run_0011_ns_v2` (V4、(1) の結果 A の後)**: 同じ設定・固定の閾値で停止規則 v2 を本物で回す。測るもの: 判定 1 回の開始〜終了の秒、要求を置いた時刻・step から STOPPED_AT までの秒・step (今の記録 `judged_at` だけでは測れないので見張りの記録に足す)、`stop_wait_s` の足り具合、δ_E の抽出の所要、容量 (§4.8 の 5)、状態と理由、オンラインの判定と `check_settled.py --end-step` の一致。停止に到達した正例が実在することを要件にする (codex M4)。この run の停止の後の妥当性は (1) で見るもので、run_0011 では見ない。
  - **プリセット全体の `provisional: false` は保留**: 小さいノズル (別の条件) の確認をしていないので (codex M4)。小さいノズルの case の選定は後で (既存の候補の README・設定・判定の記録を読んで決める)。
  - 前提への異議の採用 (codex): run_0010 は新方式の正例として確定していない (旧の系列は ALL STEADY だが残差は NOT CONVERGED)。新しい数値流束・局所の壁の熱流を含む全条件の正解のラベルとして使わない。
- 判定不能は合格にしない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-15 | [2026-10-15-tooling-nozzle-card-local-mesh-plan-2.md](../../notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan-2.md) (起点の card plan の §4.4 に対するもの) | NO-GO, C0/M7/m0 | 全件**採用**してこの plan に書き直した: M1 被覆・欠測・必要な列・step の単調の検査 (§4.2・§4.3、V2)、M2 跳ねとバーストの検査 (§4.3、V2)、M3 量ごとの出力の設定と `interfaceDiag` (§4.1 の 5・§4.2、V2)、M4 協調停止の書き込みの成功の確認 (§4.1 の 2、V1 (e); 実装中の #1 へ伝えた)、M5 流量の正解は面ごとの独立の集計・演算の順の固定 (§4.1 の 4・§4.3、V1 (d)・V3; #1 へ伝えた)、M6 共通の部品と移行の範囲・区間の管理の引き継ぎ (§4.2・§4.7)、M7 現行の仕様を 1 つに (この plan)。この plan で plan 段レビューをやり直してから #2〜#7 |
| plan | 2026-10-15 | [2026-10-15-tooling-settle-convergence-and-monitor-plan.md](../../notes/reviews/2026-10-15-tooling-settle-convergence-and-monitor-plan.md) | GO-with-changes, C0/M5/m1 | 全件**採用** (§4.8 と §6 に反映): M1 判定の格子と止めた step・残差の番号の対応・出力一式の完成の印 (§4.8 の 1、V4 の正例 80,001)、M2 V3 を残っているデータの部分ごとの再生と新しい出力の run での全条件の検証に分ける (V3)、M3 受け渡しの関門の条件ごとの置き換えと残すもの (§4.8 の 2、V4 の移行の試験)、M4 リミッタの基準値を区間の鍵に・自動の値の確定の順序・launch-binding plan との分担 (§4.8 の 3、V5)、M5 壁の熱流の絶対 + 相対の許容 (§4.8 の 4、V5)、m6 容量の条件を測れる形に (§4.8 の 5) |
| diagnose | 2026-10-16 | [2026-10-16-settle-v1-and-line-height-diagnose-settle-v1-and-line-height.md](../../notes/reviews/2026-10-16-settle-v1-and-line-height-diagnose-settle-v1-and-line-height.md) (ブリーフ [briefs/2026-10-16-settle-v1-and-line-height.md](../../notes/reviews/briefs/2026-10-16-settle-v1-and-line-height.md)) | 結論: V1 の合格の確定を保留し、出力の処理の前後で状態が変わらないことを数値更新 0 step の A/B で確かめる。Major 4・Minor 2 | (M1 600 step の差が同じ桁だから無影響とは言えない) → **採用** (§5.1 #1c の状態不変 A/B)。(M2 V1 (a)・(f)・(g) は事後の改訂と明記し元の条件は判定不能で残す) → **採用** (§6 V1 の改訂)。(M3 η 0.03 でも格子によっては軸の節点が混ざる) → **採用**: `settle.py` の線の補間で軸の節点の重みを記録し 1e−6 を超えたら判定不能 (格子と線の定義の不整合)、試験 X20・X21 (test_settle 56/56); run_0010 の格子では重みの最大 3.3e−16。(M4 η 0.03 を試験部全体の一様性と扱わない) → **採用** (card plan §6 V5c の名前を「近軸の線の一様性」に、η 0.1 は比較用に残す)。(m 公差 0.1 %・0.05 % の継承は可、ただし同じ厳しさの保証ではない) → **採用**。(m η・x 範囲・補間・点数・公差を一緒に版管理、報告と監視の定義を揃える) → **採用** (card plan §5.1 #5c)。(前提への異議: 「軸の節点の癖」は原因の確定でない、η の固定は正規化した位置で実寸の半径でない) → **採用** (card plan の文言を直した) |
| diagnose | 2026-10-16 | [2026-10-16-settle-v3b-v4-preregistration-diagnose-settle-v3b-v4.md](../../notes/reviews/2026-10-16-settle-v3b-v4-preregistration-diagnose-settle-v3b-v4.md) (ブリーフ [briefs/2026-10-16-settle-v3b-v4-preregistration.md](../../notes/reviews/briefs/2026-10-16-settle-v3b-v4-preregistration.md)) | 案のままの閾値の確定は NO-GO、C0/M4/m1。結論: 閾値を今の暫定値に固定し、run_0012 の単独の実行で「仮想の停止の後の同じ軌道」を検証する事前登録に改める | 全件**採用** (§6「V3(b)・V4 の実施」を書き直した): M1 T_q を整定の正解にしない → 停止候補 t_s からの累積の変化で検証 (人工の例: 80 窓すべて SETTLED でも累積は許容の 4.29 倍)。M2 残差・流量を最大 × 1.5 で較正しない (流量の 1 % の不釣合いが通る) → 流量は固定、× 1.5 は使わない。M3 壁の熱流の 2 係数は同じ方法で決まらない → 今は固定、確定には熱設計の許容が要る。M4 同じ条件の再実行は独立の検証でない・停止の正例が要る → 自由継続の run を先に、V4 は閾値の固定の後、プリセットの provisional: false は小さいノズルまで保留。m1 同時実行は時間の測定を乱す → 1 本ずつ、判定の所要と停止の遅れを秒と step で記録。前提への異議 (run_0010 は正例として確定していない) も採用 |

## 7. 影響範囲

forge の出力 (キー無しは今とビット同一)、判定の道具 (新しい道具を足し、今の道具は残す)、停止規則、要求カードの雛形の既定、AGENTS.md の収束確認と Stop フック (最後; 他のセッションに効く)。

## 8. 完了条件

V1〜V5 がすべて合格、§5.1 #7 の移行 (呼び出し側が共通の部品を使う、AGENTS.md と Stop フックの書き換え) が済み、result 段のレビューを通す。ノズル以外のケースのプリセットはこの plan の完了条件に入れない (対象外を返すことだけを確かめる)。

## 9. 変更ログ

- `2026-10-16` — V3(b)・V4 の事前登録を codex diagnose (2026-10-16-settle-v3b-v4-preregistration-diagnose-settle-v3b-v4.md、NO-GO) を全件採用して書き直した: 閾値は暫定値に固定、run_0012 (停止なし) を先に単独で、t_s から 8 万 step の累積の変化で検証、V4 はその後。
- `2026-10-16` — #6 (雛形の既定・リミッタの基準値の記録と区間の鍵) を実装 (implementer)、#5 と一緒に commit。系列の継承の自動化を #6b に。
- `2026-10-16` — #5 (停止規則 v2) を実装 (implementer)。DIVERGED と直らない未設定の判定不能でも協調停止して FAILED にする形を採用 (§4.4)。
- `2026-10-16` — codex diagnose (2026-10-16-settle-v1-and-line-height-diagnose-settle-v1-and-line-height.md) を全件採用: V1 を事後の改訂として書き直し (元の条件は判定不能)、#1c の状態不変 A/B を追加、線の補間の軸の節点の重みの検査を settle.py に追加。
- `2026-10-16` — #2〜#4 (判定の部品・プリセット・CLI) を実装・commit (implementer)。test_settle 54/54・run_settle_presets_tests 37/37 (主セッションで再実行して確認)。格子点 0 は観測なしと決めた。
- `2026-10-16` — #1 (forge の出力) を実装・commit (implementer)。V1 の (b)・(c)・(d)・(e)・(g) の挙動は合格、(a)・(f)・(g) のビット同一は変更前のバイナリが再実行でビット一致しないため判定不能 → 決定的な部分の検査で代替、上位に諮る (保留)。
- `2026-10-16` — 随時の全場の書き出し (§4.1 の 6、ユーザ指摘) と線の既定 r/r_w = 0.03 (ユーザ決定) を足した。
- `2026-10-15` — plan 段レビュー (GO-with-changes、Major 5・Minor 1) を全件採用して §4.8・§6 に反映。
- `2026-10-15` — 起票。card plan (`tooling-nozzle-card-local-mesh.md`) の §4.4 と §5.1 #3q・#3r・#3v・#3x・#3y・#4a・#4b、codex 諮問 2 回と plan 段レビュー (NO-GO、Major 7、全件採用) をまとめて現行の仕様に書き直した。
