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
8. **停止規則 v2 の run** (2026-10-18、`feature/card-local-mesh` への取り込みの後; codex 諮問
   [judge-v2](../../notes/reviews/2026-10-18-orchestrator-judge-v2-diagnose-judge-v2.md) を全件採用。作業は plan `tooling-design-orchestrator` §5.1 #5b、
   読み方の正本は同 plan §4.2「解析の成立の読み方」): `run_case.sh` と同じ道具の選択 (`settle.select_judge`) で、`check_settled.py` を選ぶ run の `convergence` は
   `check_settled` (照合つき・書き込みなし; PASS = 設定した量・区間で SETTLED、残差の低下は意味しない)。停止の時の再評価は `stop_rule.match_reeval` と
   停止の状態・`STOPPED_AT` で照合し、通らなければ判定不能 (判定し直しの FAIL は FAIL のまま併記)。`quasisteady` は同じ結果を criterion `settle` で共有し、
   仕様の量と基準の量の照合は省かない。v2 の状態 `FAILED` → FAIL。鮮度の入力は判定と入力の一覧を作る共通の処理から取り、今のファイルを計算し直す。

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
| 12 | 計算の後の判定を 1 本に (§4.6) | 順に: ① 確定の時点 (単段 / 段階起動、runner と `stop_rule.postprocess` の後) → ② `JUDGE_SPEC.json` (機種の既定・`judge.yaml`) と責務の分け方 → ③ 入力の完全性 (共通の判定ツールの欠損の扱いを直す) → ④ 壁解像の適用 (実効の粘性と BC、投入の関門で壁の出力を確かめる) → ⑤ 状態・終了コード・鮮度・Stop フック → ⑥ 統合試験 (欠損の反例、段階起動、監視の停止、solver の異常終了、判定器の例外、古い判定、定数粘性の壁、他のセッションの除外、関門で未起動、エコーとプロットの維持、終了コードの互換)。変更対象: tools (judge_run・check_*・run_case.sh・hook)、design の runner・stop_rule・設計ゲート、AGENTS.md | ⏳ 途中 (2026-10-15、implementer、前半: ③・本体・単段): ③ 70ef3d59 — `check_convergence --segment` は区間の段の履歴が無い・空・step の行が無いとき判定不能 (前段だけで PASS していた)、`check_quasisteady` は step の重複・逆行・評価量ゼロ (`--quantity asym` だけで適用不能)・`--series-cols` 空を INDETERMINATE (終了 1)、`check_wall_resolution` は M6 の 86dace29 と同じ規則 (指定の群・step が無ければ判定不能); 負例 7 件 (対照 3 件) を `test_gate_bad_input.py` に足し、旧コードで 7 件とも NG・新コードで OK (M6 の壁の受入 7 件も同じ形で入れた)。本体・単段 dd33fbc3 — `judge_run.py` (項目 solver・nonfinite・physical・convergence・stop_rule・quasisteady・wall_resolution、真理値表 FAIL > INDETERMINATE > PASS > NOT_APPLICABLE、終了 PASS 0・FAIL 1・INDETERMINATE 2・例外 3・段の記録 4、`RUN_VERDICT.json`/`.txt` を原子的に、`--check-fresh`)、`JUDGE_SPEC.json` 版 1 と `judge.yaml` の書き方 (docstring・procedures/calculation-workflow.md「計算の後の判定」)、`run_case.sh` は副判定の後に `judge_run.py --from-run-case --solver-rc` (段階起動 = manifest か `FORGE_JUDGE_DEFER=1` は `state: deferred`)。試験 `tests/unit/test_judge_run.py` ALL PASS (94 件)、既存の unit 6 本・test_gate_bad_input・design 3 本・check_plans は合格のまま。残り (後半): runner と `stop_rule.postprocess` の後の確定 (`FORGE_JUDGE_DEFER=1` を S1 と stages=none の本段にも渡す)、生成器の `JUDGE_SPEC.json`、投入の関門で壁の出力、Stop フック、AGENTS.md、性能試験。**実装 (2026-10-15、implementer、後半)**: ① 確定の時点 — 新しい `design/forge_design/evaluate/judge.py` の `Staged`・`finalize` を全 runner から呼ぶ。`runner.run_forge`・`runner_sern.run_forge`・停止規則の見張りの `run_case.sh` は `FORGE_JUDGE_DEFER=1` (S1・stages=none の本段も段の記録だけ)、親は段の記録・manifest・`STOP_RULE_STATE.json`・`stop_rule.postprocess` の後に系列の CSV を書いて `judge_run.py <run>` を 1 度だけ (段の途中の例外でも起動していれば; 関門で未起動なら確定しない)、実行中は `JUDGE_PENDING.json` (pid・host)。結線: runner_axismach (Euler `run_staged`・NS `run_staged_ns`、停止規則の forge の終了コードを渡す)・runner_wt・runner_walldriven・runner_sern (2D・3D・driver_sern)・runner (`run_staged_tp`・`continue_tp`・`runner.main`)・`continue_cpg` の後の起動 (opt/driver の soft と本段・opt/polish)・euler_loop (パス・継続チャンク)。② 仕様 — 生成器が最初の起動の前に `JUDGE_SPEC.json` を書く (`judge.spec_*`; 中身は `judge.py` の docstring の表): axis-Mach NS = 出口コア M・質量の収支 r (`criterion: absolute`、停止規則の |平均| 5e-4・振幅 5e-3)・θ_r × 3 (停止規則を使う run だけ)、系列 `judge_series.csv` は `stop_rule.json` の行を `STOP_RULE_STATE.json` の record_sha256 と一致するときだけ使い、それ以外は本段の res から同じ式で; Euler = 出口コア M・r (本段の出力 3 枚では「窓の不足」で INDETERMINATE とした — 軸 M などの報告量があるので理由つきの適用外にはしない); SERN = 目的量 (粘性 C_T_with_shear / 非粘性 C_T)・C_T・C_L・C_M (`force_history.csv`、sern_forces の 0.4/0.02/0.05; collect と同じ関数で確定の時点に書く); ベル TP = tp_gates (ii) の 5 列 (`thrust_history.csv`、1.5e-4・ṁ 1e-3); ベル CPG・WT・walldriven・euler_loop = snapshot の shock・machmax・pmax。同じ数値設定の継続 (polish のチャンク・euler_loop の継続チャンク; `stage_key` で判定) は `from_floor`。`judge_run.py` に `criterion: absolute` (STEADY / OSCILLATING / OFFSET) と `notes`。③ 投入の関門 — 壁解像が適用になる run (実効か仕様) で no-slip 壁の `outputHDFflg` ≠ 1 は FAIL (移行期・`FORGE_LAUNCH_GATE=warn` でも)。④ 物理性に静圧 ≤ 全圧 (`total_quantities.total_state`; 許容 `P ≤ P0 (1 + 1e-5)`、実測の丸めは case/57・case/56 の 25 run で最大 6.4e-7、case/45 run_0138 (TP) 1.0e-6; h0・物性が無ければ判定不能)。⑤ Stop フック — `RUN_VERDICT.json` の有無・`check_fresh`・`state: deferred` (親が生きていれば見逃す) を見る; FAIL・判定不能・判定器の例外は止めない; claim は保つ。⑥ AGENTS.md「収束確認」に `RUN_VERDICT.txt` が正本の 1 項、モデル分担の Stop フックの記述を同期、procedures/calculation-workflow.md「計算の後の判定」の後半の記述。⑦ `nozzle_report.QS_VERDICTS` に INDETERMINATE。⑧ 性能 (`tests/unit/perf_judge_run.py`、45 万点・出力 80 枚 4.60 GB・snapshot の量は 80 枚を読む・壁 1 枚・残差 16 万行): 1 回目 3.5 s・2 回目 3.0 s (≤ 120 s); WSL2 では `posix_fadvise` でページキャッシュを落としたと確かめられない (1 枚 56 MB の読み直し 0.023 s vs 0.016 s) ので**温かい状態**。⑨ 試験: test_judge_run ALL PASS (118 件; 静圧 ≤ 全圧・absolute・Stop フック 12 件)、test_launch_gate ALL PASS (壁の出力 10 件)、design/tests/run_judge_finalize_tests.py (新規 27 件、段階起動・例外・関門で未起動・系列・仕様・判定器の例外・古い判定) ALL PASS、run_stop_rule_tests に確定の判定 (監視の停止・CENSORED・異常終了) を足して ALL PASS、既存の unit 4 本 (recipe_settings・mesh_provenance・check_solver_config・restart_field_mesh_identity)・test_wall_resolution・test_gate_bad_input・design の run_launch_path・run_generator_recipe・run_stop_rule_settling・run_recipe_profiles (+ bell_tp・opt・h370_prelaunch・sern_design_path・walldriven・main_line・bell_cv_gate) は合格のまま (run_sern_gates_tests は HEAD a555e611 でも同じ IndexError で落ちる既存の失敗)。commit 411759bd・b2cbccac・51706a3a |
| 12a | ~~#12 前半の判断 7 点~~ | 判断: 2026-10-15 (主セッション): (1) 停止規則の対応 CENSORED (上限で未到達)→FAIL、UNAVAILABLE→NOT_APPLICABLE (他の判定は続ける) — 採用。(2) 段の記録 (deferred) は確定ではなく、既知の FAIL があれば FAIL と表示 — 採用。(3) 不明な solver の終了コードは INDETERMINATE — 採用。(4) 計算の開始の後に書いた `judge.yaml` は記録だけで結果に使わない (仕様は実行の前に確定) — 採用。(5) 壁は常に全 no-slip 壁 — 採用。(6) 物理性は ρ・P・T > 0 だけ → 静圧 ≤ 全圧 (`total_quantities.py` の h0 から) を後半で足す。(7) h5 の入力の識別は名前・大きさ・mtime_ns — 採用 (速さのため; 判定の記録に明記) | F |
| 12b | ~~#12 後半の判断 6 点~~ | 判断: 2026-10-15 (主セッション): (1) axis-Mach の Euler は本段の全場の出力が 3 枚で準定常の窓が足りず INDETERMINATE — 量は残す (not_applicable にしない)。全場でなく入口・出口の境界の出力 (1,000 step ごと、小さい) から出口コア M・質量の収支の系列を作るよう直す (#14)。(2) 閾値の方針: **その機種に受理のゲート (停止規則・SERN のゲート・ベル TP のゲート) があれば、判定の閾値はそのゲートと同じ** (PASS の意味を揃える): axis-Mach NS は θ_r・出口コア M を停止規則の θ_r と同じ (drift 5e−4・振幅 5e−3)、r は停止規則の絶対の閾値、SERN・ベル TP は今のまま。ゲートの無い機種 (CPG・WT・walldriven・euler_loop) は `check_quasisteady` の既定のまま (#14)。(3) CPG・風洞の snapshot の `shock` は既存のゲートに合わせて残す。(4) Stop フックの pid の再利用の見逃しは既知の限界として記録。(5) `from_floor` の参照が PASS でないと INDETERMINATE — 採用 (AGENTS.md の規則どおり)。(6) `launch_gate.py` の `note_initial_field` が `AREAS` に無い区分 '初期場' を使う既存の欠陥 → 直す (#14) | F |
| 14 | #12b の後始末 | (1) axis-Mach の Euler・NS の準定常の系列を境界の出力から作る、(2) axis-Mach NS の閾値を停止規則と同じに、(6) `note_initial_field` の区分を `AREAS` に。**実装 (2026-10-15、implementer)**: (6) 済 — 区分を既存の `来歴` に (`launch_gate.py` の `note_initial_field`)、負例 (初期場の来歴の記録が JSON でない run: 旧コードは関門の内部の例外 FAIL で起動せず 3 件 NG、新コードは `来歴` の NOTE だけで起動) を `test_launch_gate.py` に。(1)(2) は**未実装**、方針の穴として主セッションへ返した: (1) forge は境界の出力を全場と同じ `outStepInterval` で書く (`output.cpp` 321・339 行) ので Euler の本段の点数は増えない、axis-Mach の入口は `outputHDFflg: 0` (`runner_wt._bcond`)、入口の出力はゴーストの状態 (`inlet_Pressure_d`) で節点の値でない、境界の出力に `sonic` が無い (`mesh.hpp` の `bplaneValNames`)。(2) `classify` の窓は点数の割合・drift は最小二乗の傾き × 末尾の幅・fluct は全幅で、停止規則の D (20,000 step の窓の前後半の平均の差、2 万 step あたり)・amp (半幅) と定義が合わない (osc だけは 2 倍で対応) | O |
| 14a | ~~#14 の (1)(2) の判断~~ | 判断: 2026-10-15 (主セッション; 報告: 境界の出力は全場と同じ `outStepInterval` で書かれ、入口の出力は無く、境界の値は節点の値と定義が違う。`check_quasisteady` の drift は最小二乗の傾き × 末尾窓の点数の幅で、停止規則の D・amp と定義が違う): (1) Euler の本段の全場の出力の間隔を n/6 (12,000 step なら 2,000) に下げて点数を 6 枚にする (道 b; 境界の出力の別の間隔は #13 の範囲で後回し)。生成器の段の config が変わるので、プロファイルの書き写しと例外表を更新し、変わるのは出力の間隔の 1 行だけであることを試験で確かめる。(2) axis-Mach NS の θ_r・出口コア M・r の判定は、新しい criterion `stop_rule` で停止規則の `judge_rows` そのものを呼んで判定する (道 b; 定義を二重に持たない、閾値は `stop_rule.json` の settings から)。`check_quasisteady` の relative は機種のゲートが無いときだけ | F |
| 14b | #14a の実装 | (1)・(2)。**実装 (2026-10-15、implementer)**: (2) 済 — `judge_run.py` に criterion `stop_rule` (量の役 `stop_rule_role`: r・theta・relative)。停止規則の `stop_rule.judge_rows` そのもの (`design/forge_design/evaluate/stop_rule.py` を judge_run から import; この作業ツリーのファイルでなければ判定不能) で系列の最後の点の時点を判定し、窓・閾値・断面は run の `stop_rule.json` の settings (`stop_rule.thresholds` で検査; 無ければ判定不能)。出口コア M を θ_r と同じ式・閾値で判定するため `judge_rows` に `relative_keys` (既定 () = 停止規則の判定は不変) を足した。仕様に thresholds・window は書かせない (SpecError)、accept は [STEADY] だけ。停止規則の到達と並べた照合を `stop_rule_cross` (一致・不一致と理由) に記録し RUN_VERDICT.txt に 1 行 (総合は変えない)。道具の sha256 に stop_rule.py。`judge.spec_axismach`: 停止規則を使う NS は exit_core_M・r・θ_r × 3 を criterion stop_rule に、使わない NS・Euler は今のまま。試験: test_judge_run (+11: 到達・未到達で judge_rows と値・合否が同じ、照合の一致・不一致、窓の不足・settings 無し・stop_rule.json 無し・断面に無い列、仕様の誤り 6 種) ALL PASS、run_judge_finalize_tests (+11: (i) Euler 6 枚で STEADY・DRIFTING を判定でき 3 枚は窓の不足、(j) 記録 → judge_series.csv → criterion = judge_rows) ALL PASS、run_stop_rule_tests の本段の判定の確認を criterion stop_rule + 照合の一致に直して ALL PASS (偽の forge + 本物の見張りの REACHED と一致)、既存 (run_stop_rule_settling・run_launch_path・run_generator_recipe・run_recipe_profiles・test_launch_gate・test_recipe_settings・test_gate_bad_input・check_plans) は合格のまま。**(1) は未実装** — 方針の穴として主セッションへ返した: case/ の axis-Mach の問題 YAML 300 本はすべて `evaluate.outStepInterval` を明示している (12,000/4,000 が 213 本・24,000/8,000 が 15 本 = 3 枚、12,000/1,000 が 47 本 = 12 枚 ほか) ので、`runner_axismach.prepare` (987 行) の既定 `max(n // 3, 1)` を n/6 にしても既存の問題の出力の間隔は 1 つも変わらない。書き写しの代表 (`recipes/fixtures/axismach_euler_tp.yaml`) も 12,000/1,000 なので段の config・バイト同一の比較は変わらない (`regen_profiles.py --check` 一致) | O。**判断 (2026-10-15、主セッション)**: (1) Euler の出力の間隔は**変えない** — 3 枚になるのは case/41 の問題 YAML が `outStepInterval` 4,000 を明示しているからで、カードの雛形 (`axismach_final_v1`、書き写しの代表も) は 12,000 / 1,000 = 12 枚 (M5.5 の run_0001_euler も 12 枚)。生成器の既定 n/3 を変えても既存の問題には効かず、明示の値の上書きは採らない。古い case/41 の Euler は準定常が判定不能のまま (理由つき)。(2) `stop_rule.judge_rows` に `relative_keys=()` を足したこと (既定の動きは不変) は採用 |
| 15 | 停止規則 v2 の run の判定 (§4.6 の 8) | ⏳ 作業と合格の条件は plan `tooling-design-orchestrator` §5.1 #5b・§6 J (2026-10-18 の取り込みの後; 取り込みの前は v2 の run が無かった)。判断: 2026-10-18 codex (diagnose) judge-v2・全件採用 | O |
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
- (#12、2026-10-15) 計算の後の判定: `solver_density_cuda/tools/` の `judge_run.py`・`run_case.sh`・`check_*`・`launch_gate.py` (壁の出力の FAIL)・
  Stop フック `hook_convergence_gate.py` (判定の有無・鮮度・段階起動の確定)。
- (#12 後半) 設計チェーンの runner (`design/forge_design/evaluate/runner.py`・`runner_axismach.py`・`runner_wt.py`・`runner_walldriven.py`・
  `runner_sern.py`・`runner_sern3d.py`・`feedback/euler_loop.py`・`opt/driver.py`・`opt/polish.py`) と新しい `evaluate/judge.py`: 全段
  `FORGE_JUDGE_DEFER=1`・親の 1 度の確定・`JUDGE_SPEC.json`・系列の CSV。停止規則 `evaluate/stop_rule.py` (見張りの `run_case.sh` に
  `FORGE_JUDGE_DEFER=1`)。生成器の数値・段の config は変えない。
- 設計ゲート (停止規則の到達・受け渡しの関門・`sern_gates`・`tp_gates`) の受理の条件は変えない (`judge_run.py` は共通の検査だけ)。
  `report/nozzle_report.py` の `QS_VERDICTS` に INDETERMINATE。
- 振る舞いの変化: 粘性の run で no-slip 壁の `outputHDFflg` が 0 だと forge を起動しない。`runner.run_forge` を直接呼ぶ古い case の台本
  (例 case/16・case/41) の run は `state: deferred` のまま残るので、最後に `judge_run.py <run>` を回す (Stop フックが止める)。

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
