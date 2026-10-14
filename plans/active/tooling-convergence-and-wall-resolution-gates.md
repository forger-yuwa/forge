---
status: in_progress
owner: CFD Dev
created: 2026-09-19
related_docs:
  - methods/turbulence/implementation.md
  - procedures/verification/README.md
  - AGENTS.md
---

# 収束判定と壁解像判定の恒久対応 (リポジトリ全体)

## 1. 背景 — なぜケース内の修正では足りないか

`case/49` で、**2 種類の判定の甘さ**が実際に誤報告を生んだ (2026-09-19)。どちらもケース固有では
なく、**段階起動を使う全ケース・node 方式の全ケース**で同じことが起きる。

### 1.1 収束判定: 別の方程式・BC の履歴を混ぜていた

段階起動 (S0 全 slip → S1 層流 → S2 等温 → S3/S4 SST → S5 2 次ランプ → S6 本段) では、
**本段だけを `check_convergence.py` にかけると必ず「plateau」**になる (前段の収束場から
始まるので低下桁数が小さい)。これを避けようとして**全段を連結**して判定したが、
`check_convergence.py` は**系列全体の最大値を低下桁数の基準**に取るため、
**slip や層流の大きな過渡を本段の基準にしてしまう**。結果、実際には本段で 2 桁しか
落ちていない run が「起動全体で 3.3〜8.0 桁低下」と報告された。

### 1.2 壁解像: ソルバの `ypls` は $y_1^+$ ではない

`viscousFlux_d.cu` の低 Re 経路は
$$\texttt{ypls\_b} = \rho\,u_\tau\,d_{cc}/\mu$$
で、$d_{cc}$ は**ゴーストセル重心と内点セル重心の距離**である。
**node 方式では壁ノードが壁面上に乗るので $d_{cc}$ が退化**し、値が 1 桁以上小さく出る
(ソルバ自身、流束計算では `dcc` の退化を避けるために別経路を通している)。
`case/49` では「y⁺ 平均 0.04 だから低 Re 解像は十分」と報告したが、**第一内部ノード基準の
$y_1^+$ は平均 0.15〜0.30・最大 2.6〜4.4** で、リップ近傍では 1 を超えていた。

### 1.3 共通する構造

AGENTS.md で**守られているルール**(run パス明示・VERDICT 貼付・case README の run 一覧) は、
いずれも「**ツールがあり、応答に痕跡が残る**」形をしている。収束・準定常はこの形になっている
一方、**壁解像にはツールも貼付義務も無かった**。また収束ツールは**段階起動を知らない**ため、
使い手が誤った区間を渡せてしまう。

## 2. 目的

- 段階起動の run で、**同じ方程式・BC・空間離散化の区間**でしか収束合否を出せないようにする。
- 壁解像 ($y_1^+$) を**正しい定義で出すツール**を用意し、報告に VERDICT の貼付を義務づける。
- 準定常ツールが**欠損列・非有限値を合格にしない**ようにする。

## 3. 現状の確認 (実測)

| 事実 | 根拠 |
| --- | --- |
| 本段のみの判定は必ず plateau | `case/49` run_0101〜0103 の `CONVERGENCE_VERDICT.txt` |
| 全段連結は別 BC の過渡を基準にする | `check_convergence.py` の `analyze()` が系列全体の最大値を使う |
| node の `ypls` は 1 桁小さい | `case/49` run_0103: ソルバ 0.0221 に対し第一ノード基準 0.146 |
| `--series-csv` 経路が NaN を通す | `[1,1,1,1,1,NaN]` が `STEADY` を返す (codex 再現) |

## 4. 設計方針 (2026-09-19 codex plan レビュー NO-GO/M9 を受けて全面改訂)

**方針の柱**: ①**不正入力を合格にしない** (欠損・非有限・空検査・判定不能を明示的に落とす)、
②**履歴の同一性を名前でなくデータで保証する**、③**壁量は壁面ごとの局所量として出す**。

### 4.1 判定ツール共通の入力検証 (最優先。これ自体が実バグ)

codex が現行コードで再現した**誤合格**を先に塞ぐ。いずれも「ケース固有」ではない。

| ツール | 現行の穴 | 対応 |
| --- | --- | --- |
| `check_convergence.py` | **残差列が 1 つも無い CSV で `ok=True`**。`rms_ro` だけでも `ok=True` | 必須保存量列を確定し、**欠損・空検査は判定不能で非ゼロ終了** |
| 〃 | 19 点が 1 で**最終 1 点だけ 0** → `drop=inf` で合格 ([:153](../../solver_density_cuda/tools/check_convergence.py)) | 低下桁数を**末尾窓の代表値**で測り、末尾の 1 点ゼロで合格させない |
| `check_quasisteady.py` | CSV 経路は正しいが **HDF5 経路が今も `classify` を呼ぶ** ([:448](../../solver_density_cuda/tools/check_quasisteady.py)) | 全経路を `classify_series` に寄せ、**step の有限性・単調性**も検査 |
| `check_mesh_quality.py` | **体積ゼロ四面体が通常合格** (AR 1.732 / skew 0.500)。**NaN 座標 1 セルで `SOFT-PASS` exit 0** | 非有限座標・接続不正・ゼロ/負体積・不正 Jacobian は**外れ値許容の前に無条件拒否** |
| 呼び出し側 | `run_case.sh` が `|| true` で判定失敗を吸収、`sern_gates.py` は CLI を通さず `analyze()` 直呼び | **共通判定 API** を用意して両方をそこに寄せる |

参考にすべき既存実装: `check_field_regress.py` は欠損・非有限を事前拒否、
`check_passive_budget.py` は検証不能を `INDETERMINATE` にしている。**この扱いを共通化する**。

### 4.2 収束判定 — 段名でなく**設定の同一性**で区間を決める

**却下した当初案**: 段名プレフィックスで連結する / 共通列だけ残す。
- 段名は方程式・BC・メッシュ・組成・再構成の同一性を**保証しない**。
- `residual_history_all.csv` のような集約済みファイルが同じ命名に入り**二重連結**する。
- 共通列への縮退は、層流段と SST 段を選ぶと `rms_roK/rms_roOmega` が消え、
  **検出したかった誤合格をそのまま再導入する**。

**採用案**:
- 段ごとに **`stage_manifest.json`** (順序・履歴ファイル・restart 元・実効 config のハッシュ・
  メッシュ/BC/入口分布の識別子) を run 生成側が書く。
- `check_convergence.py` は**残差の定義が同じ区間だけ**連結する。CFL・反復数の変更は連結可、
  方程式・BC・離散化の変更は**不可**として拒否する。
- **`--segment main` による本段単独判定は許可する**。「本段だけなら必ず plateau」は一般命題として
  誤りで、本段で十分落ちていれば合格してよいし、落ちていなければ `NOT CONVERGED` が正しい動作。
- 既定では自動連結せず、複数段が曖昧なら**選択を要求**する。`--from-floor` の参照 run にも同一性検査。

### 4.3 壁解像 — 壁面ごとの**局所**第一内部点と**接線**応力

**却下した当初案**: 場の `wall_dist` の全域下位パーセンタイルを全壁に適用 / Sutherland 固定 /
「平均 ≤1 かつ最大 ≤5」。
- 全域の 1 つの距離では、**細かい壁群が 1 つあると粗い壁群の解像不足を隠す**。
- `wall_dist` は**最近接壁点への距離**であって対象壁面への法線距離ではなく、slip 壁も距離計算に入る。
- 粘性は定数 / Sutherland / 多成分 Chapman–Enskog+Wilke の**3 経路**がある (`gasProperties_d.cu`)。
- SST の壁解像は**局所**の壁間隔に感度があり、平均で合格させてはいけない。

**採用案**:
- 壁面要素 (境界半割面) ごとに**接続情報から第一内部点を対応付け**、局所法線方向の距離 $y_1$ を出す。
  node は節点座標、cell は実セル中心。`ransWallFunction_d.cu` の `Normal_Neighbor` 選択が再利用候補。
  角・斜交・複数壁の曖昧さは**検出して「評価不能」にする** (ゼロ距離や合格に変換しない)。
- $u_\tau$ は**接線 traction** から $u_{\tau}=\sqrt{|\boldsymbol\tau_{w,t}|/\rho_w}$ で組む。
  **ソルバの `utau` をそのまま使わない**: node の既定経路は `twall_*` だけを上書きし
  `utau`/`ypls` を更新しないため不整合が残る (`case/49` run_0103 の `cav_outer` で
  $|\tau_w|/(\rho u_\tau^2)$ は中央値 0.996 だが**最大 74.5**。最大値が出る高せん断域がまさに壊れている)。
- 粘性は run の輸送モデル・組成・壁面状態に対応するものを使い、**必要データが無ければ判定不能**。
- 判定: **低 Re は評価対象領域の局所 $y_1^+\le1$ を既定目標**とし、**超過面積と位置**を報告する。
  緩和は熱流束・摩擦の格子感度で裏付けた場合だけ認める。
  mode 1 は forge では **automatic treatment** (粘性低層〜対数層を接続) なので
  `30≤y_1^+≤300` を全点必須にしてはいけない。**低層/buffer/log の面積分布**を出し、
  対数層配置を要求する運用だけ別の明示条件にする。mode 2 の対象範囲も明記する。
- 壁群ごとの**面積加重統計・最大値の位置・評価できた面積割合**を出す。周期の重複点・軸境界・
  slip 境界は明示的に区別する。

### 4.4 ソルバ側

- `ypls_b` の**定義 (dcc 基準)** と mode ごとの違いをコメントと `methods/turbulence/implementation.md` に明記。
- **当面は値を変えない**。ただし理由は「同じ配列だから」ではない (572 行で既に mode 分岐している)。
  正しくは「**後処理で独立した `y1plus` を作る方が検証しやすい**」。将来ソルバ出力へ移すなら、
  壁近傍対応を初期化時に構築し、node 壁応力補正後の診断経路から**別配列**へ出す。

### 4.5 AGENTS.md への追記

- 「**壁解像確認 (必須)**」節: $y_1^+$ を報告する応答に新ツールの VERDICT 貼付を義務づけ、
  **ソルバ出力の `ypls` を壁解像の根拠にしない**ことを明記。
- 「収束確認」節: **段階起動では判定区間を明示する**ことを追記。

### 4.6 計算の後の判定を 1 本にして自動で回す (ユーザ 2026-10-15「計算後の収束判定は完全に機械化してほしい」「あほなAIでもね」; 2026-10-15 plan 段レビュー judge を全件採用して改訂)

今は `run_case.sh` が forge の後に `check_convergence.py` を**履歴全体**で回すだけで、判定区間・準定常・非有限の監査・物理性・壁解像・停止規則の状態は
エージェントが手で順に回している (所要 5〜10 分、漏れ・取り違え・区間の誤りの実例あり)。弱いモデルでも判定を誤らないよう、**選ぶ場面を作らずに**自動で回す。
`judge_run.py` は既存の判定ツールの**共通の集約層**で、設計の受理の条件 (deltastar_loop の関門・SERN のゲートなど) はそれぞれのゲートに残す。

1. **判定の仕様を実行の前に確定する** (`JUDGE_SPEC.json`、run の中): 量の定義・単位・抽出器・対象の領域・判定区間・窓・閾値・必須か・適用外の理由。
   設計チェーンの run は機種の既定から生成器が書く (ノズル: 出口コア M・質量の収支 (絶対の平均・振幅も)・θ_r など、停止規則の系列は定義・出典・区間が一致するときだけ再利用。SERN: C_T_with_shear・C_L・C_M)。
   手書きの run は `judge.yaml` (case か run に置く) から作る。どちらも無ければ準定常は「対象未設定」= INDETERMINATE (黙って省かない)。継続の run の `--from-floor` と非定常の扱いも仕様に書く。
   壁解像の適用は**実効の粘性と境界条件**で決める (`viscMethod 0` の定数粘性も粘性: `gasProperties_d.cu:91`)。適用になる run で壁の出力 (`outputHDFflg`) が無ければ**投入の関門で止める** (終わってからでは計算し直しになる)。
2. **確定の時点**: 単段の run は `run_case.sh` が forge の後に判定する。段階起動の run は、`run_case.sh` は段ごとの実行の結果と健全性だけを保存し、
   親 (runner) が段の履歴・manifest・停止の記録の確定・出力の整理を終えた後に総合の判定を**一度だけ**確定する (`runner_axismach` の段の記録 2414・2467 行の後、`stop_rule.postprocess()` の後)。
3. **入力の完全性**: 要求した対象と実際に評価した対象の完全一致を確かめる。段・列・量・壁・時刻の欠損、空集合、時刻の重複・逆行、要る窓の不足は INDETERMINATE。
   共通の判定ツール側も直す (`check_convergence.py` の欠けた段の `continue`・`check_quasisteady.py` の評価量ゼロで STEADY・step が進まない系列で STEADY・`check_wall_resolution.py` の step の代替; plan 段 judge の M3 の再現)。
   manifest が無いときの単段の扱いは、単段と確かめられるときだけ。
4. **状態と終了コード**: 原判定を保持した構造化の結果と、総合の真理値表を先に決める。総合は PASS / FAIL / INDETERMINATE / NOT_APPLICABLE (理由つき)。既知の FAIL は判定不能と併存しても隠さない。
   停止規則の状態は `REACHED / CENSORED / DIVERGED / EXEC_ERROR / DATA_ERROR / UNAVAILABLE / NOT_USED` をそのまま。solver の異常終了・正常な監視の停止・判定器の例外を区別する。
   `CONVERGENCE_VERDICT.txt` は収束の副判定として残し、`FORGE_STRICT_CONVERGENCE` の意味は変えない。起動の関門・`LAUNCH_ECHO`・残差のプロットは今のまま。
5. **鮮度**: 判定に実行の識別子・確定の状態・入力の識別 (config・manifest・系列・壁の出力・停止の記録・`JUDGE_SPEC.json` の sha) を記録し、原子的に書く。
   Stop フックは同じセッションの範囲で「判定済みか・今の入力に対応するか」を確かめる。**FAIL や判定不能を正しく報告して終えるのは許す** (PASS になるまで止める設計にしない)。
6. **監査の範囲**: 通常は AGENTS.md どおり残差の全履歴 + 最終の場 (+ 準定常の抽出で読む場はその場で監査、停止規則の体積場の監査の記録を再利用)。全スナップショットの完全な監査はオプション。
   所要の目標 (45 万点・出力 80 枚で 2 分以内) は出力量・媒体・キャッシュの条件を固定した性能試験として登録する。
7. 報告の規則 (AGENTS.md): 「収束した」「定常」と書く応答には `RUN_VERDICT.txt` の総合と各判定の行を貼る。

## 5. 実装

### 5.1 残作業 (優先順。codex の推奨順「①不正入力拒否と段情報 → ②局所 y1・接線応力・物性 → ③閾値と適用範囲 → ④呼び出し側移行と独立検証」に合わせた)

| # | 項目 | 状態 |
| --- | --- | --- |
| 1 | `check_convergence.py`: 残差列ゼロ/必須列欠損を**判定不能で非ゼロ終了**、低下桁数を**末尾窓**で測る | ✅ |
| 2 | `check_quasisteady.py`: HDF5 経路を `classify_series` に寄せ、step の有限性も検査 | ✅ |
| 3 | `check_mesh_quality.py`: 非有限座標・ゼロ/負体積・向き混在を**外れ値許容の前に無条件拒否** | ✅ |
| 4 | 不正入力テスト `test_gate_bad_input.py` (10 件) + 実データ回帰 (過去 PASS の 3 run・実メッシュ 3 種) | ✅ |
| 5 | 呼び出し側: 検証を `analyze()` 内に置いたので `sern_gates.py` は自動継承。`run_case.sh` は判定 rc を `CONVERGENCE_EXIT:` として記録し `FORGE_STRICT_CONVERGENCE=1` で伝播 | ✅ |
| 6 | 段情報 `stage_manifest.json` (`stage_manifest.py` の `StageManifest`) と `check_convergence.py --segment` ✅。hard キー (convMethod/limiter/turbulenceModel/wallTreatmentSST/viscMethod/thermalMethod/unsteady/dualTime/axisymmetric + bcond のハッシュ) が同じ連続区間だけ連結。**列が区間内で違えば連結しない** (共通列への縮退は誤合格を再導入するため)。case/49 の段構成で S5 ランプ以降+本段が 1 区間と正しく出る |
| 7 | `check_wall_resolution.py`: 壁面ごとの局所 $y_1$ (PLANES/STRUCT の接続から法線方向の第一内部点)・接線 traction からの $u_\tau$・run の `viscMethod` に応じた $\mu$・壁ダンプ節点と DOF の**厳密一致による対応づけ**。自己検査: $y_1$ が第一層厚と一致 | ✅ |
| 8 | 壁解像の小規模テスト `test_wall_resolution.py` (一様直交 / 非一様 = 壁ごとに別 y1 / 斜交 75° = 評価不能 / 角 = 2 面共有) ✅。⏳ 軸対称・周期の配置は未追加 |
| 8b | 判定は**超過面積割合** (`--over-frac`, 既定 2 %) で行う。最大値は幾何的特異点で収束しないので使わない ✅ |
| 9 | 実ケース展開 ✅: **`case/26`** (cell・壁関数) が設計 3 水準を再現 (y1 4.35/115/356 µm → y1+ 0.79/25.5/66.5)、**`case/48`** (冷却平板) が等温壁で断熱の 5 倍 (0.101→0.506) = 既知挙動と一致、**`case/49`** (node・低 Re) で `cyl_top` の全面 y1+>1 を検出。⏳ TP ケースは未 |
| 10 | AGENTS.md に「壁解像確認 (必須)」節 + 収束節に判定区間と判定不能の 2 項目 ✅ / `methods/turbulence/implementation.md` §3.7.x に `ypls`/`utau` の mode 別定義と node での退化を記載 ✅ / `viscousFlux_d.cu` に同趣旨のコメント ✅ (コード変更なし=解はビット不変) |
| 11 | `check_cavity_steady.py` を `classify_series` 経由に ✅ / `cavity_eval` の y1+ は参考値と明記し正式判定は共通ツールへ ✅ / ⏳ `stack_residuals.py` は #6 の段情報が入ってから統合 |
| 12 | 計算の後の判定を 1 本に (§4.6) | 順に: ① 確定の時点 (単段 / 段階起動、runner と `stop_rule.postprocess` の後) → ② `JUDGE_SPEC.json` (機種の既定・`judge.yaml`) と責務の分け方 → ③ 入力の完全性 (共通の判定ツールの欠損の扱いを直す) → ④ 壁解像の適用 (実効の粘性と BC、投入の関門で壁の出力を確かめる) → ⑤ 状態・終了コード・鮮度・Stop フック → ⑥ 統合試験 (欠損の反例、段階起動、監視の停止、solver の異常終了、判定器の例外、古い判定、定数粘性の壁、他のセッションの除外、関門で未起動、エコーとプロットの維持、終了コードの互換)。変更対象: tools (judge_run・check_*・run_case.sh・hook)、design の runner・stop_rule・設計ゲート、AGENTS.md | ⏳ 途中 (2026-10-15、implementer、前半: ③・本体・単段): ③ 70ef3d59 — `check_convergence --segment` は区間の段の履歴が無い・空・step の行が無いとき判定不能 (前段だけで PASS していた)、`check_quasisteady` は step の重複・逆行・評価量ゼロ (`--quantity asym` だけで適用不能)・`--series-cols` 空を INDETERMINATE (終了 1)、`check_wall_resolution` は M6 の 86dace29 と同じ規則 (指定の群・step が無ければ判定不能); 負例 7 件 (対照 3 件) を `test_gate_bad_input.py` に足し、旧コードで 7 件とも NG・新コードで OK (M6 の壁の受入 7 件も同じ形で入れた)。本体・単段 dd33fbc3 — `judge_run.py` (項目 solver・nonfinite・physical・convergence・stop_rule・quasisteady・wall_resolution、真理値表 FAIL > INDETERMINATE > PASS > NOT_APPLICABLE、終了 PASS 0・FAIL 1・INDETERMINATE 2・例外 3・段の記録 4、`RUN_VERDICT.json`/`.txt` を原子的に、`--check-fresh`)、`JUDGE_SPEC.json` 版 1 と `judge.yaml` の書き方 (docstring・procedures/calculation-workflow.md「計算の後の判定」)、`run_case.sh` は副判定の後に `judge_run.py --from-run-case --solver-rc` (段階起動 = manifest か `FORGE_JUDGE_DEFER=1` は `state: deferred`)。試験 `tests/unit/test_judge_run.py` ALL PASS (94 件)、既存の unit 6 本・test_gate_bad_input・design 3 本・check_plans は合格のまま。残り (後半): runner と `stop_rule.postprocess` の後の確定 (`FORGE_JUDGE_DEFER=1` を S1 と stages=none の本段にも渡す)、生成器の `JUDGE_SPEC.json`、投入の関門で壁の出力、Stop フック、AGENTS.md、性能試験 |
| 12a | ~~#12 前半の判断 7 点~~ | 判断: 2026-10-15 (主セッション): (1) 停止規則の対応 CENSORED (上限で未到達)→FAIL、UNAVAILABLE→NOT_APPLICABLE (他の判定は続ける) — 採用。(2) 段の記録 (deferred) は確定ではなく、既知の FAIL があれば FAIL と表示 — 採用。(3) 不明な solver の終了コードは INDETERMINATE — 採用。(4) 計算の開始の後に書いた `judge.yaml` は記録だけで結果に使わない (仕様は実行の前に確定) — 採用。(5) 壁は常に全 no-slip 壁 — 採用。(6) 物理性は ρ・P・T > 0 だけ → 静圧 ≤ 全圧 (`total_quantities.py` の h0 から) を後半で足す。(7) h5 の入力の識別は名前・大きさ・mtime_ns — 採用 (速さのため; 判定の記録に明記) | F |
| 13 | 判定と停止規則の量を全場の出力でなく注目の位置の出力から取る (後回し、ユーザ 2026-10-15「h5丸ごと出すのは容量的に大変だから、追って、着目する位置（出口分布？）を自動出力するようにしてほしい。probeなど？まだそこまで話はひろげなくていい」) | 停止規則 (3 断面の δ*・θ_r、質量の収支) と準定常の量に要る位置 (出口の断面の分布・3 断面の壁法線の線・入口と出口の流量) を `JUDGE_SPEC.json` から決めて、forge の probe・線・境界の出力で自動で書き、全場の `res_*.h5` の出力の間隔を広げる。今は 1,000 step ごとに全場 (M5.5 で 1 枚 約 50 MB) | ⏸ 後回し |

## 6. 検証

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 | 対応 |
| --- | --- | --- | --- | --- |
| plan (§4.6 judge) | 2026-10-15 | [2026-10-15-tooling-convergence-and-wall-resolution-gates-plan-judge.md](../../notes/reviews/2026-10-15-tooling-convergence-and-wall-resolution-gates-plan-judge.md) | GO-with-changes, C0/M6/m2 | **全件採用** → §4.6 を改訂: M1 確定の時点 (段階起動は親の確定の後)、M2 判定の仕様を実行の前に確定 (`JUDGE_SPEC.json`、責務は共通の集約層)、M3 入力の完全性 (共通ツールの欠損の扱いを直す)、M4 壁解像の適用は実効の粘性と BC (定数粘性も) + 投入の関門で壁の出力、M5 状態の真理値表と終了コード、M6 鮮度と Stop フック (FAIL・判定不能の報告は許す)、m7 監査の範囲は AGENTS.md どおり、m8 統合試験と影響範囲 (#12) |
| plan | `2026-09-19` | [`notes/reviews/2026-09-19-tooling-convergence-and-wall-resolution-gates-plan.md`](../../notes/reviews/2026-09-19-tooling-convergence-and-wall-resolution-gates-plan.md) | **NO-GO**, C0/M9/m0 | **全件採用**。M1 段名では同一性を保証できない・本段単独判定は許可すべき → §4.2 を書き直し / M2 共通列縮退は誤合格を再導入・`check_convergence` の空検査と最終 1 点ゼロの穴 → §4.1, §5.1 #1 / M3 `wall_dist` 全域パーセンタイルは局所距離でない → §4.3 / M4 `utau` と `twall` の不整合 (最大 74.5 倍) → §4.3 で接線 traction から組む / M5 Sutherland 固定は 3 経路ある物性と不整合 → §4.3 / M6 閾値が局所精度と automatic treatment を取り違え → §4.3 / M7 準定常の穴は HDF5 経路 → §5.1 #2 / M8 `check_mesh_quality` がゼロ体積・NaN 座標を通す → §5.1 #3 / M9 呼び出し側の迂回と独立検証の欠落 → §5.1 #5, #8, #9 |

### 6.2 合否条件

- **不正入力テスト**: 残差列ゼロ / 必須列欠損 / 最終 1 点だけゼロ / `[1,1,1,1,1,NaN]` /
  NaN 座標メッシュ / ゼロ体積四面体 が**すべて非ゼロ終了**する。
- **正常入力の回帰**: 既存ケースの判定が変わらない (単段 CLI の互換性)。
- **壁解像**: 解析的に $y_1$ が分かる小規模テストで一致し、`case/26` / `case/48` で妥当な値を返す。
  **`case/49` の 0.146/2.60 は `utau` 不整合の影響を受けるため合否条件にしない**。
- 判定区間・評価できた面積割合・判定不能の理由が**出力に残る**。

## 7. 影響範囲

- ソルバのコード変更は**無し** (コメントのみ)。解はビット不変。
- `solver_density_cuda/tools/` の 3 ファイルと `AGENTS.md`、`methods/turbulence/implementation.md`。

## 8. 完了条件

§5.1 が全て完了し、§6.2 の合否条件を満たし、codex の plan / result レビューを通ること。

## 変更ログ

- `2026-10-15` — plan 段 judge (GO-with-changes C0/M6/m2) を全件採用して §4.6・#12 を改訂。
- `2026-10-15` — ユーザ依頼「計算後の収束判定は完全に機械化してほしい」: §4.6 (判定を 1 本の道具 `judge_run.py` にまとめ、`run_case.sh`・設計チェーン・Stop フックから自動で回す) と §5.1 #12 を起案。

- 2026-09-19 (3): #8 (小規模テスト 4 配置)・#9 (case/26, case/48, case/49 展開)・#11 の一部を実装。
  判定を**超過面積割合**に変更。`case/49` 最細格子で **`cyl_top` が全面 y1+>1** (平均 1.22) と判明 —
  入熱 84 W とキャビティ 3 壁の 2 倍を受ける面が壁解像できていない (ソルバの `ypls` はここで 0 を返していた)。
- 2026-09-19 (2): §5.1 #1〜#5, #7, #10 の一部を実装。実測: `case/49` 全ヘキサ 3 格子で
  $y_1$ が第一層厚 (28.57/20.00/14.29 µm) と厳密一致、$y_1^+>1$ の面積が外筒壁で
  13.5 → 10.5 → **7.8 %**、平板は最細で **0 %**。**最大値はリップ特異点のため収束しない**
  (18.9 → 14.4 → 11.1) ので判定は超過面積割合で行う。過去に報告した「y⁺ 0.04」も
  「$y_1^+$ 平均 0.146」も誤りだった (前者は `dcc` 退化、後者は全域共通 $y_1$ と stale な `utau`)。
- 2026-09-19: 起票。`case/49` の result レビュー (codex NO-GO, C1/M10/m1) の M3・M4・M9 を
  ケース内で塞いだが、**同じ誤りが他ケースでも起きる構造**なのでユーザ指示により恒久対応に格上げ。
