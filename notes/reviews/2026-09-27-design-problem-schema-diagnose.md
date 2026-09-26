# codex 諮問 (diagnose): design-problem-schema

- **brief**: [`notes/reviews/briefs/2026-09-27-design-problem-schema.md`](../../notes/reviews/briefs/2026-09-27-design-problem-schema.md)
- **date**: 2026-09-27
- **commit**: `02f871a8` (feature/gap-heating-precision)
- **codex**: effort `xhigh`, 5.8 min, rc=0
- **結論**: **まず axismach に、副作用なしで評価要求を実行計画へ解決する `resolve_evaluation(problem, candidate, recipe)` を一つ設け、既存 Euler/NS 経路を明示的に選択できる境界を検証する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **H2 を採用**：機種別 runner＋小さな共通契約 | ベルは固定の `DV_ORDER` と「長さは厳密値」という目的関数実装、SERN は別の変数・目的関数・多作動点処理を持つ（[driver.py:48](/home/sano/work/forge/design/forge_design/opt/driver.py:48)、[driver_sern.py:3](/home/sano/work/forge/design/forge_design/opt/driver_sern.py:3)）。共通化するのは**候補の受け渡し・評価要求・結果の形式**とし、物理と評価手順は runner が所有する。 |
| **Major** | **H1 の「機種非依存 stage DAG を先に作る」は却下** | δ* 経路は前回場から補正を抽出し、**物理壁を変更して**再評価する（[deltastar_loop.py:135](/home/sano/work/forge/design/forge_design/feedback/deltastar_loop.py:135)）。単純な精度向上の段列ではない。反復・緩和・停止条件は機種別 recipe に閉じ込め、外側は名前付き recipe の呼び出しに留める。 |
| **Major** | **「YAML の dv 名を変えれば MOO に乗る」は却下** | [length-dv plan:31](/home/sano/work/forge/plans/accepted/tooling-nozzle-axismach-length-dv.md:31) はコードと矛盾する。ベル runner の固定 import、固定変数順、固定目的関数がある（[driver.py:40](/home/sano/work/forge/design/forge_design/opt/driver.py:40)、[driver.py:77](/home/sano/work/forge/design/forge_design/opt/driver.py:77)）。対案は、adapter が公開する変数・評価量を campaign が選択する契約。import の差し替えだけでは足りない。 |
| **Major** | **明示的な評価 recipe を採用** | `prepare` は Euler/slip を生成する一方、NS は別関数・別入口である（[runner_axismach.py:461](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:461)、[runner_axismach.py:641](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:641)）。CLI 名に依存する物理選択を、検証可能な評価要求へ移す。 |
| **Major** | **semiperfect の γ/cp の単純削除は却下** | 既定値への暗黙フォールバックがあり、出口診断・NS 初期化で使用される（[probdef.py:167](/home/sano/work/forge/design/forge_design/probdef.py:167)、[runner_axismach.py:582](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:582)、[runner_axismach.py:822](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:822)）。対案は熱力学モデルと近似用参照定数の分離。SERN の `frozen_tp` は設計側で `GasCPG` を使うため、その用途も残す。 |
| **Minor** | **「CFL 4・内部反復 5 は一律に旧設定」は要再検証** | [recommended-settings.md:36](/home/sano/work/forge/procedures/recommended-settings.md:36) は 4 を示すが、Euler 専用節は CFL 4・内部反復 5 を指定する（[同:197](/home/sano/work/forge/procedures/recommended-settings.md:197)）。対案は適用範囲と版を持つ recipe。全機種への一括置換は支持しない。 |

結論: **まず axismach に、副作用なしで評価要求を実行計画へ解決する `resolve_evaluation(problem, candidate, recipe)` を一つ設け、既存 Euler/NS 経路を明示的に選択できる境界を検証する。**

第 1 仮説: **主因は YAML が一枚であることではなく、宣言された変数・物理・評価手順と、実際に呼ばれるコードの間に検査可能な契約がないこと。** 確度: 高  
  根拠: `KNOWN_TYPES` は入力の許可リストであり、runner の dispatch ではない（[probdef.py:161](/home/sano/work/forge/design/forge_design/probdef.py:161)）。さらに `--prepare-only` は YAML の `cfl_main` を段階実行へ渡す前に終了する。CLI の `--cfl` がなければ準備時は 4.0、通常実行では後から YAML 値を適用する（[runner_axismach.py:477](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:477)、[同:620](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:620)）。**同じ problem でも入口によって実効設定が変わる具体例**である。  
  反証条件: 明示的な評価要求と入力成果物をすべて与えても、既存 Euler/NS の計画を表現するために共通層へ機種固有の物理判断を追加する必要があるなら、この境界だけで十分という設計仮説は棄却する。

  推奨構造は **「problem に設計の基準値、campaign に探索と評価要求、機種別 recipe に実行手順」**。物理的なファイル分割は二つでよい。以下は新スキーマの骨格であり、現在の loader が受け付ける形式ではない。探索範囲も検証前の例である。

```yaml
# problem.yaml
schema: forge.problem/v1
name: va3
kind: wind_tunnel_axisym_axismach

definition:                 # kind ごとの schema で検証
  gas:
    model: semiperfect
    composition_basis: mole
    species:
      H2O: 0.0609135
      N2: 0.664860
      O2: 0.216072
      AR: 0.00797588
      CO2: 0.0490034
  spec:
    Pt: 1139000.0
    Tt: 1161.0
    r_throat: 0.205711
    M_design: 4.19
  parameters:
    L_c: 8.0                # 単位・独立性は adapter の変数定義が持つ
  geometry:
    R: 2.0
    L_U: 6.0
    Lc_mode: explicit
---
# campaign.yaml
schema: forge.campaign/v1
problem: problem.yaml
seed: 0

search:
  method: grid
  dv:
    L_c: {min: 7.5, max: 8.0, count: 3}
  evaluator: {recipe: axismach.design/v1}
  objectives:
    - {metric: length_m, direction: minimize}
  constraints:
    - {metric: geometry_feasible, equals: true}

qualification:
  select: {from: search.feasible, best: 3}
  evaluations:
    - id: euler
      recipe: axismach.euler/v1
      inputs: {design: candidate.design}
      options:
        mesh: {ni: 365, nj: 65}
        numerics_recipe: va3_euler_r1

    - id: corrected
      recipe: axismach.deltastar_ns/v1
      inputs:
        design: candidate.design
        euler_reference: euler.result
      options:
        initializer: contur
        numerics_recipe: va3_ns_r1
```

  この構造で固定すべき契約は次の四点。

- **変数の所有者を決める。** adapter が変数名・単位・独立／従属・対応モードを公開し、campaign が探索対象と範囲を指定する。problem には基準値だけを置く。`Lc_mode: from_length` では `L_c` を探索対象として受け付けない。ユーザの探索範囲と設計成立領域は別物で、候補ごとに両方を検査する。現在の bound 検査は値が矩形範囲内かしか見ず、設計側は別途許容窓を検査している（[probdef.py:221](/home/sano/work/forge/design/forge_design/probdef.py:221)、[runner_axismach.py:325](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:325)）。
- **評価 recipe は機種と能力を宣言する。** 未対応のキー・物理の組合せは実行前に拒否する。ケース固有スクリプトも同じ入出力契約で包める。soft/mid/main、δ* の反復、多作動点の実行は recipe 内部に置く。SERN の作動条件は problem、目的関数への重み付けは campaign が所有する。
- **成果物はファイル名以上の情報を持つ。** 設計点・物理形状・メッシュ・組成配置・エンタルピー基準を識別し、生成元 recipe とコード版を記録する。δ* 補正後は別の形状として扱う。後続の凝縮・3D・FEM がどちらの形状を評価するかを入力参照で固定する。同一メッシュ restart と cross-mesh 移植もここで区別する。
- **生の VERDICT と用途上の採否を分離する。** 結果契約には評価量・単位・作動点・判定区間・判定ファイルを残す。`NOT CONVERGED` を campaign の採用判断で `CONVERGED` に書き換えない。要求された評価の未実行・判定欠落も合格にしない。

  詳細解析を後段に置けるのは、それが確認用途の場合に限る。**熱制約が最適化の実行可能性を定義するなら、その判定は探索側の制約評価に接続する必要がある。** 勝者だけの FEM では「熱制約付き最適化」を実施したことにはならない。片方向 FEM は追加 evaluator、双方向 CHT は界面反復を所有する別 recipe とする。両者は現在も別の未実装機能として整理されている（[CAPABILITIES.md:57](/home/sano/work/forge/design/CAPABILITIES.md:57)）。

第 2 仮説: **現場の差し替え需要の一部は、汎用 workflow がないためではなく、既存 runner の不具合・規約追従不足による。** 確度: 高。`run_staged` の同一メッシュ移植は `interp_field.py`（[runner_axismach.py:538](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:538)）、ケース側は `restart_field.py`（[run_lumpX_staged.py:96](/home/sano/work/forge/case/44.vitiated_air_wt/run_lumpX_staged.py:96)）。これは **Major**。対案は recipe 内の移植処理を修正することであり、DAG 導入では直らない。

判別 A/B: **同じ候補・物理条件・参照成果物を固定し、`recipe` だけを Euler／δ* NS に切り替える、実行計画の比較を一つ行う。**

- A: `axismach.euler/v1`
- B: `axismach.deltastar_ns/v1`
- 長さ: 各一回の解決処理。**CFD は 0 step、forge・メッシュ変換・書き込みは呼ばない。**
- 見る量: 選択関数、壁 BC、実効数値設定、必須入力、成果物の形状識別、検査予定の VERDICT。既存経路との差分を項目単位で比較する。

**A は Euler/slip、B は NS・指定熱境界・δ* 入出力へ解決され、説明不能な設定差分 0、無視された有効入力 0、共通層への物理分岐追加 0 なら、第 1 仮説を支持する。どちらかが CLI やケース内の後編集に依存したら、この最小境界は不十分として棄却する。** この試験が確認するのは配線と責務分担であり、数値的妥当性ではない。

やらない方がよいこと: **先に汎用 DAG エンジンを作ること、axismach 専用へ縮退すること、全 YAML を一括移行すること。** 捨てるべきなのは、スクリプト内定数を仕様の正本にする運用、暗黙の設定無視、「全機種の dv と目的関数は既に汎用」という前提である。既存の `explicit`／`max` は保持決定があるため削除対象にしない（[length-dv plan:27](/home/sano/work/forge/plans/accepted/tooling-nozzle-axismach-length-dv.md:27)）。

呼び出し側の前提への異議:

- **Minor — 観測事実 7 の解釈。** モル分率入力を runner が質量分率へ変換して `Y{s}` を出すこと自体は誤りではない。[runner_axismach.py:151](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:151) は解決済み輸送種の質量分率を取得している。対案は入力基準・輸送種配置・出力基準を記録し、`X` 出力が必要な要求と、組成が物理的に正しいことを分けて検証すること。
- **Major — 「ALL STEADY」を全体の合格とする解釈は受け入れない。** 以下の保存済み判定は三件とも `NOT CONVERGED`、メッシュは `PASS`、選択された六量の系列は `ALL STEADY` だった。準定常判定に記載された系列は **7 行、step 0–24000**。対案は、この対象量と区間を含む別々の判定として保持すること。

| 根拠 run | 残差判定 | 派生量判定 |
|---|---|---|
| `case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX/` | [NOT CONVERGED](/home/sano/work/forge/case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX/CONVERGENCE_VERDICT.txt:2) | [ALL STEADY](/home/sano/work/forge/case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX/QUASISTEADY_SERIES_VERDICT.txt:10) |
| `case/44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX/` | [NOT CONVERGED](/home/sano/work/forge/case/44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX/CONVERGENCE_VERDICT.txt:2) | [ALL STEADY](/home/sano/work/forge/case/44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX/QUASISTEADY_SERIES_VERDICT.txt:10) |
| `case/44.vitiated_air_wt/run_0511_va3_M4.19_Lc8_eq_lumpX/` | [NOT CONVERGED](/home/sano/work/forge/case/44.vitiated_air_wt/run_0511_va3_M4.19_Lc8_eq_lumpX/CONVERGENCE_VERDICT.txt:2) | [ALL STEADY](/home/sano/work/forge/case/44.vitiated_air_wt/run_0511_va3_M4.19_Lc8_eq_lumpX/QUASISTEADY_SERIES_VERDICT.txt:10) |

不足情報: 原データは指定どおり読んでいないため、上記は保存済み判定の確認であり、抽出処理・判定区間の独立検証ではない。FEM の実際の入力・出力・許容値も未提示なので、その内部 schema はまだ固定できない。今回の推奨は診断上の仮説であり、**plan 未反映。ファイル変更・forge 実行は行っていない。**

## 当方の採否 (2026-09-27, 主セッション)

| 指摘 | 採否 | 理由 |
|---|---|---|
| H2 採用 (機種別 runner + 小さな共通契約) | 採用 (方向として) | ベル/SERN の DV_ORDER・目的関数が機種固有なのはコードで確認済み (driver.py:40,77 / driver_sern.py:41) |
| 汎用 stage DAG を先に作らない | 採用 | δ* は壁を変えて再評価する反復で、段の直列では表せない (deltastar_loop.py:135) |
| 「dv 名を変えれば MOO に乗る」の却下 | 採用 | length-dv plan:31 の記述は誤り。plan への訂正は起票時に行う |
| 評価 recipe を明示 (Euler/NS を入口名でなく要求で選ぶ) | 採用 | 今回 `--prepare-only` と自作段階起動に分けざるを得なかった実例がある |
| γ/cp は単純削除でなく「熱力学モデル」と「近似用参照定数」に分離 | 採用 | frozen_tp の設計側 GasCPG 用途を見落としていた |
| run_staged の同一メッシュ interp_field 使用は Major | 採用 | recipe 化とは独立に直すべき不具合 |
| recommended-settings の Euler 節 (cfl 4 / nStepInner 5) との整合 | 採用 | 今回 run_0509–0511 は §1/§3 の値 (cfl 6+relax 0.7, nStepInner 4) を使った。Euler 設計評価 §4 と食い違う点を README に注記する |
| 「ALL STEADY を全体合格としない」 | 同意 (既にそう報告済み) | case README に NOT CONVERGED と STEADY を別々に記載済み |
| 最初の一歩 = axismach に副作用なしの `resolve_evaluation` を作り Euler/δ* NS の計画差分を検証 | 保留 (ユーザ判断待ち) | 構造の方向はユーザが決める。plan 未反映 |
