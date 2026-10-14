**未成立 (incomplete)** — incomplete: 証拠不足: M5.5 の元カード・固定した設計壁・`m55_mesh_est` の再現台本と出力。ni 3,163、AR 4,850、補正後約45万点を追試するために必要。 / M6 補正係数の元となる run・step・壁量と、対応する収束／準定常 VERDICT、冷却倍率および比率分布の生成記録。

# codex レビュー: tooling-nozzle-card-local-mesh (plan)

- **plan**: [`plans/active/tooling-nozzle-card-local-mesh.md`](../../plans/active/tooling-nozzle-card-local-mesh.md)
- **stage**: `plan`
- **date**: 2026-10-15
- **commit**: `db47ee8d` (feature/card-local-mesh)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 4.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-nozzle-card-local-mesh-plan.codex.log` (codex: incomplete)
- **判定**: **NO-GO**, 指摘 C0/M7/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO

局所的な第一層・流れ方向密度を既存の `mesh2d` で生成する方向は妥当です。
ただし、第一層の補間保証に反例があり、既定変更で既存の成果物出力も欠落します。現案のまま実装へ進めません。
M5.5 の点数見積もりと M6 の補正係数は、提示された作業ツリーだけでは再現できませんでした。

1. **Major — §4.1 の下包絡の作り方では、補間後の第一層厚を保証できない。**

   根拠: [plan:53](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:53) は「実寸厚さの最小を取ってから局所半径で割る」方式ですが、[mesh2d.py:103](/home/sano/work/forge-cardmesh/design/forge_design/meshing/mesh2d.py:103) が log 線形補間するのは**半径で正規化した比率**です。

   既存関数で反例を実行しました。要求厚さを一定の 10⁻⁵、両端の半径を 1、2、中点を 1.5 とすると、中点の生成厚さ／要求厚さは **1.060660**、推定 y₁⁺ は **0.848528**。両端が要求以下でも中間で超過します。

   **対案:** 先に f_req(x) = c(x)·y₁,req(x)/r_w(x) を作り、その下包絡を表にする。表の節点だけでなく補間区間内を検査し、不適合区間を細分する。さらに、見積もりは設計壁、実際の NS 格子は物理壁なので、最終物理壁でも実寸厚さを再検査する必要があります。実生成は [runner_axismach.py:1814](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:1814) の物理壁を使います。

2. **Major — `wall_normal_layer` の既定化で、最終カードの既存成果物が欠落する。**

   根拠: [mesh_structure.py:233](/home/sano/work/forge-cardmesh/design/forge_design/export/mesh_structure.py:233) は、同じ列の x 座標のずれが全長の 10⁻⁹ を超える格子を拒否します。壁法線層はこの条件に抵触します。[deliverables.py:195](/home/sano/work/forge-cardmesh/design/forge_design/export/deliverables.py:195) は例外を捕捉して戻るため、`mesh_structure.json`・`rebuild_from_step.py` と再構築検証が作られません。これは既存仕様の [methods/design/overview.md:1342](/home/sano/work/forge-cardmesh/methods/design/overview.md:1342) にある成果物です。

   **対案:** 壁法線層を表現できる構造ファイルと STEP 再構築を、既定切り替えの前提作業に追加する。M5 の final カードで「生成→NS 準備→成果物→再構築座標・接続照合」まで検証する。現在の「`wt_card.py`・雛形・試験だけ」という影響範囲は不十分です。

3. **Major — M6 の表を「実測による保守的な補正係数」とする根拠が不足している。**

   根拠: [M6 の H の YAML:3](/home/sano/work/forge-cardmesh/case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml:3) は、上流の表を C/H の半径比で変換したものです。同ファイル [97 行](/home/sano/work/forge-cardmesh/case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw300_cgref_h.yaml:97) も、断熱 NS に冷却倍率を掛けた**見込み**と記しています。300 K の実測要求厚さそのものではありません。

   また、plan 自身の縮流部の比は **0.61〜2.2** なのに、採用係数は縮流部 **1.0**。比が 1 未満の場所を含む区間全体について、これを保守側とは断定できません。`x≈−2 r_t` の遷移位置・幅も未定義です。

   **対案:** 係数は「初期格子用の経験係数」と位置づけ、元の run・step・壁形状・物性・冷却倍率・比率の分布を固定する。縮流部を含め、補間後の c(x) が根拠とする比率を超えないことを検査する。M5 空気と M5.5 燃焼ガスへの既定適用は、実際の壁解像確認を必須条件にする。

4. **Major — V2 は現状では実行できず、AR 緩和の条件も定量化されていない。**

   根拠: 記載コマンドを実行すると、現行ツールは **`unrecognized arguments: --bcond`** で停止します。引数定義は [check_mesh_quality.py:198](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_mesh_quality.py:198) です。

   また、`--ar-max 5000` の PASS は全域の skew ≤ 0.9 を確認するだけです。「壁直交性 PASS」の判定器・閾値が未指定では、高 AR セルに必要な条件を保証できません。既に [runner.py:807](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner.py:807) に **AR > 1000 の全セルで skew ≤ 0.1** を確認する `mesh_ar_relax_check` があります。

   **対案:** 正しい品質コマンドと既存の AR 緩和検査を使い、変換後の全セルで違反数 0 を要求する。壁法線層の遷移域も対象に含める。薄い第一層については M6 と同じ `msh_digits: 17`、変換後の第一層厚の照合、使用した変換器・solver の精度と版も検証条件へ明記する。

5. **Major — V4 を「記録だけ」で完了にでき、未達・別時刻の壁量でも閉じられる。**

   根拠: [plan:81](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:81) と [96 行](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:96) は、壁解像未達でも計画を完了できる構成です。一方、停止規則の `REACHED` は収束・準定常判定とは別であると [stop_rule.py:57](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:57) が明記しています。

   さらに [check_wall_resolution.py:316](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_wall_resolution.py:316) は指定 step が無いと最新ダンプへ代替します。指定した壁群が一部欠落していても、存在する群だけがループ対象です。「全壁・停止 step」を引数に書くだけでは保証になりません。

   **対案:** 初回 V4 の失敗は再調整の契機として許し、**最終完了には再調整後の V4 PASS を要求**する。全 no-slip 壁の同一 step の出力存在を先に検査する。`check_convergence --segment` と、壁せん断・熱流束などの系列に対する `check_quasisteady` の VERDICT、未収束時の扱いを事前登録する。格子修正後は `interp_field.py` による restart、品質・節点数の再判定を手順に含める。

6. **Major — 見積もりと実生成の設定契約、既存 overrides の移行規則が抜けている。**

   根拠: [_ar_mesh_params:784](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:784) と実生成の [_mesh_params_from:738](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:738) に同じ設定を渡して確認したところ、前者では `wall_first_frac_table`・`x_density_table`・`wall_normal_layer` がすべて **None**、後者では指定値が保持されました。

   また、現行雛形には `wall_first_up_x0/x1` があり、新しい表と併存すると既存関数が実際に `ValueError` を返します。plan は併用禁止を認識していますが、ユーザーが明示した旧 overrides の扱いまでは決めていません。

   **対案:** 見積もり・節点数確認・品質確認・実生成で同じパラメータ構築を使う。雛形由来の旧キーは除去し、明示 overrides の競合は説明付きで拒否する。自動生成表と明示表の優先順位、`ni` 固定時の不適合処理も定義する。AR 制約は固定値 4850 ではなく、明示された `mesh.ar_max` に従わせる。

7. **Major — Euler を影響範囲外とする記述と、精度検証の省略が成立しない。**

   根拠: [wt_card.py:2267](/home/sano/work/forge-cardmesh/design/forge_design/intake/wt_card.py:2267) は NS の `ni` を Euler にコピーします。NS の `ni` を変える本変更は、Euler の格子も変えます。

   accepted の [mesh-size-guard:52](/home/sano/work/forge-cardmesh/plans/accepted/tooling-nozzle-mesh-size-guard.md:52) は、この結合を維持する決定を記録しています。同 [55 行](/home/sano/work/forge-cardmesh/plans/accepted/tooling-nozzle-mesh-size-guard.md:55) は、最大 Δx だけでは局所形状・勾配を拘束できない点も指摘しています。したがって「ni 2000 相当の h_cap 以下」は、解の精度の証明にはなりません。

   **対案:** Euler の半径方向設定を維持することと、Euler の `ni` が変わることを分けて記載する。M5 を使い、出口 M・流量・δ_r の基準格子との差と準定常性を検証する。M5.5 でも局所格子の細分比較と許容差を登録する。節点予算は Euler・NS 両方で確認する。

8. **Minor — 関連計画・再現資料の参照が、この作業ツリーでは解決しない。**

   根拠: [plan:11](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:11) の `tooling-generator-defaults-and-launch-gate.md`、[61 行](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:61) の `h370_wall_envelope.py` は見つからず、参照先の等温壁 plan にも `#36r` がありません。`case/71` と `m55_mesh_est/` もありません。

   **対案:** 対象 commit・正しい参照先と再現入力を固定し、`methods/design/overview.md` の更新を実装前の作業へ追加する。既存の表生成機能や節点数ガードとの重複実装は避け、この plan の責務をカードへの接続に限定する。

**推奨:** `local_v1` を採用する方向を維持し、まず **①正規化した要求値の下包絡、②成果物出力の対応、③設定の一元化と品質ゲート、④補正係数の根拠、⑤Euler／NS の受入条件** の順で計画を修正し、再レビューすることを推奨します。メッシャや solver の新方式を追加する必要はありません。今回の範囲では node の M5・M5.5 が適切で、cell・周期の solver 回帰を追加する必要もありません。

ファイル変更はありません。補間の反例、設定伝達の差、CLI の不成立を実行確認しました。CFD の収束・定常性を新たに認定した結果はありません。

証拠不足: M5.5 の元カード・固定した設計壁・`m55_mesh_est` の再現台本と出力。ni 3,163、AR 4,850、補正後約45万点を追試するために必要。
証拠不足: M6 補正係数の元となる run・step・壁量と、対応する収束／準定常 VERDICT、冷却倍率および比率分布の生成記録。

指摘数: Critical 0 / Major 7 / Minor 1
