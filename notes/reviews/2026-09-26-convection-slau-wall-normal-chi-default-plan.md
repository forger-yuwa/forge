# codex レビュー: convection-slau-wall-normal-chi-default (plan)

- **plan**: [`plans/active/convection-slau-wall-normal-chi-default.md`](../../plans/active/convection-slau-wall-normal-chi-default.md)
- **stage**: `plan`
- **date**: 2026-09-26
- **commit**: `eb2c9e30` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m1
- **focus**: (i) 三値 auto の解決規則と明示 1 の検証エラーの整合 (cell/非 SLAU で auto→0 を静かに解決してよいか)。(ii) manifest で『省略』の意味が 0→1 に変わることの旧 run 互換 (キー無し = 0) と、バイナリ更新をまたいだ継続 run の区間分離がエコー正本で担保されるか。(iii) 未検証域 (周期・軸対称・凝縮) を auto に含める判断と、§6 の A/B・不合格時の除外処置が測る前に固定されているか。(iv) 許容の再利用 (case/16 0.5 %、SERN R5n 帯、局所格子 1 つ) が各ケースの量に対して適切か。(v) 設計チェーンの provenance (runner が明示 1 を書く) で旧 DB エントリと区別できるか。(vi) 本 plan は 09-25 の diagnostician 判断『既定化しない』をユーザ決定で置き換えるもの — その根拠 (傾向不変・発散防止) が §6 のどの量で検証されるか。
- **extra**: `plans/accepted/convection-slau-wall-normal-chi-usage-rule.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

三値 auto と、対象外構成での静かな `auto→0` は妥当です。ただし、旧 run の実効値保存、runner のメッシュ変換、検証範囲に実装前の修正が必要です。  
m6_on の「5 列とも許容帯内」は再計算で確認できましたが、一般的な発散防止や未検証域への既定適用までは裏付けません。  
ファイルは変更していません。以下の指摘は **plan 未反映**です。

1. **Major — 起動エコーを追加するだけでは、旧 run 互換と区間分離を保証できない。**

   **根拠:** [stage_manifest.py:215](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:215) は追加時の YAML からキーを確定します。実際に [make_case.py:304](/home/sano/work/forge-sern-design/case/51.gap_turbulent_thd8233/tools/make_case.py:304) は全段の manifest を実行前に生成しています。一方、[run_case.sh:37](/home/sano/work/forge-sern-design/solver_density_cuda/tools/run_case.sh:37) は `RUN_PROVENANCE.txt` と `forge_run.log` を実行ごとに上書きし、SERN の [_archive_stage:742](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:742) は残差だけを保存します。

   この経路では、後から最新ログを読んでも過去段の実効値は復元できません。また、旧 manifest の明示 1 は `space.slauWallNormalChi: "1"` です。「新しい `.effective` キーが無いから 0」という移行は誤りになります。旧ログも無い省略 YAML に新規則を適用すると、旧実行の 0 を 1 と誤認します。

   **対案:** manifest を版管理し、実行前の予定値と実行後の確定値を分離する。各起動のログ・設定・バイナリ hash・残差を段 ID に結び付けて保存し、終了時に実効値を確定する。旧キーの明示値を移行し、既知の旧形式でキー無しの場合だけ 0 とする。由来不明の推定値は、既知の実効値と自動連結しない。

   同じ変更で `RUN_PROVENANCE` と `diag_applicability.py` の更新も §5.1 に入れるべきです。現行 [diag_applicability.py:65](/home/sano/work/forge-sern-design/case/46.sern_design/cad/diag_applicability.py:65) に提案エコーを入力した結果は **`{}`** で、省略時には旧既定 0 に落ちます。

2. **Major — runner の明示 1 が、品質検査用の cell 変換を壊す。**

   **根拠:** [runner_sern.py:523](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:523) と [runner_sern3d.py:141](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:141) は、生成 config の離散化だけを一時的に `cell` に変更します。[convertGmshToForge.cpp:26](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/convertGmshToForge.cpp:26) も `cfg.read()` を通るため、残った明示 1 が [solverConfig.cpp:684](/home/sano/work/forge-sern-design/solver_density_cuda/input/solverConfig.cpp:684) のエラーになります。

   **対案:** 品質検査用 config は明示 0、対象の node 計算用 config は解決済みの 0/1 として別生成する。§5.1 #4 の合格条件を「文字列に 1 が出る」から、**2D・3D とも品質検査→node 変換→計算用 config 復元まで通る**ことへ変更する。明示 1 の構成エラー自体は維持してよいです。

3. **Major — config への明記だけでは、設計 DB の旧結果との混在を防げない。**

   **根拠:** [driver_sern.py:73](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:73) は既存 `ledger.jsonl` を読み込み、[同:224](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:224) は `status == "PASS"` の全行を学習に使います。[同:156](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:156) の評価要約にはフラグの実効値がありません。

   したがって、既存 campaign を再開すると旧既定 0 と新既定 1 の評価が同じ応答関数として混在します。

   **対案:** 実効フラグと設定方針の版を `metrics.json`・台帳・campaign の識別情報に保存し、不一致の評価を同じ学習集合へ入れない。旧行は保存された設定・ログから分類し、不明な行は新 campaign で再利用しない。旧台帳からの再開試験を追加する。

4. **Major — B1 は周期・軸対称・凝縮で、変更対象の流束が発火することを保証していない。**

   **根拠:** B0(d) の KEEP は auto が 0 なので、周期境界で SLAU の変更を有効にした試験ではありません。前 plan の [V6:406](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:406) も **1 step** の試験です。

   case/44 の既存凝縮入力 [run_0024…/bcondConfig.yaml:3](/home/sano/work/forge-sern-design/case/44.vitiated_air_wt/run_0024_va_R2_LU6_Lc8_split_cond/bcondConfig.yaml:3) は `slip` です。[mesh.cpp:967](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/mesh.cpp:967) が `wall_flag` を立てるのは `wall` / `wall_isothermal` のみで、この入力の対象境界数は **0** でした。このケースをそのまま選ぶと、凝縮 A/B が通っても変更分岐は未検証です。

   **対案:** 基準 run・BC・起点 dump を先に指定し、各追加領域で「対象面数 > 0」「少なくとも一部の面で変更による質量流束差が検出される」を試験成立条件にする。周期は壁を持つ SLAU ケース、軸対称・凝縮は変更対象壁を持つケースで継続計算を行う。成立しない領域を「検証済み」に数えない。

5. **Major — 「不合格ならその域を除外」は、まだ事前固定された処置になっていない。**

   **根拠:** [対象 plan:132](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:132) は周期・軸対称・凝縮の除外しか定めていません。case/44 は複合条件なので、不合格時に軸対称と凝縮のどちらを除外するか未定です。case/36、通常の SERN 2D、残差悪化、非有限値、差区間の部分重なりにも対応がありません。

   さらに runner が常に明示 1 を書けば、**auto の除外を追加しても設計チェーンには効きません**。

   **対案:** 試験ごとに不合格・判定不能と除外条件の対応表を作る。複合ケースだけで判定するなら、その複合条件を保守的に除外し、原因を単独領域へ帰属させない。通常領域の不合格は既定化を保留する。runner の省略時方針も同じ適用条件から解決して明示値を出力し、ユーザーの明示指定と区別する。

6. **Major — §6 の許容値には、測定定義と合否規則の不足がある。**

   **根拠:** [対象 plan:124](/home/sano/work/forge-sern-design/plans/active/convection-slau-wall-normal-chi-default.md:124) 以下では、次が未確定です。

   - case/36・44 の壁圧 **0.5%** は、分母・対象壁・座標窓・空間ノルムが未定義。case/16 の既存条件は特定位置の壁圧系列であり、その数値だけでは別ケースの分布比較を定義できません。
   - 「局所格子 1 つ」は、位置抽出法と採用する間隔が未定義。現在の [check_quasisteady.py:73](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:73) は `CELLS/centCoords` を使い、[shock 抽出:89](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:89) は中心帯の圧力閾値です。node 値位置と衝撃位置の定義を確認せず転用できません。
   - `asym` は準定常判定だけで、flag 間の許容差がありません。凝縮 onset も抽出法・未検出時の扱いがありません。
   - 残差比 ≤2 は、出典の [旧条件:301](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:301) にあった **RISING 列なし**を落としています。
   - Sod の同 step 数は、適応時間刻みの場合に同じ物理時刻を保証しません。

   **対案:** ケース別に抽出器、値位置、正規化、評価窓、差区間、数値許容、判定不能時の扱いを登録する。node の位置量は `MESH/COORD` で抽出した系列を `check_quasisteady --series-csv` に渡す。Sod は同一物理時刻で厳密解と比較する。

   また「同一収束場から分岐」は実態に合わせて修正が必要です。保存済み SERN 残差を再実行した VERDICT は **両側 `NOT CONVERGED (stalled/plateau)`** でした。非有限値・上昇を拒否したうえで、プラトーを許すケースと目的量による受入条件を明記してください。

7. **Major — 既定化の利益である「発散防止」を、新既定の受入試験が測っていない。**

   **根拠:** B1 は主に既存場からの比較で、既知の壁 CV 排出構成を新既定で救済する試験がありません。既存の [APPLICABILITY_v2.txt:17](/home/sano/work/forge-sern-design/case/46.sern_design/_rule_usage/APPLICABILITY_v2.txt:17) は、対象面で流出 **`4.3455e-9` → 流入 `−1.1705e-7 kg/s`** となる局所機序を裏付けますが、一般的な発散防止の証明ではありません。

   一方、m6_on の再計算では `C_L` 差区間 **[0.000790, 0.001584]**、`C_M` **[−0.031365, −0.013462]** を含め、5 列すべて帯内でした。準定常 VERDICT は flag0 の `C_L`・摩擦込み `C_L`・`C_M` が **OSCILLATING**、他は **STEADY**。これは指定窓の回帰許容を支持します。

   **対案:** 既知の接続模型で **新省略値と明示 1 の同等性・壁密度回復・床到達数・NaN の有無**を受入項目に追加する。前 plan の V1 の期限・床・準定常基準を再利用でき、3D 格子収束の全面再実施は不要です。§1 の利益は「既知構成の壁 CV 排出を抑える」と限定する。ユーザ決定は運用方針を変更できますが、検証済み範囲を拡張する根拠にはなりません。

8. **Minor — B0(c) が現行の検証運用規則と不整合。**

   **根拠:** [procedures/verification/README.md:63](/home/sano/work/forge-sern-design/procedures/verification/README.md:63) は cell の回帰対照を組まない方針ですが、B0(c) は cell の GPU 流束比較です。

   **対案:** 本変更で必要な cell 保護は設定解決・明示 1 拒否・品質検査用変換の試験へ置き換える。cell の流束比較も必要と判断するなら、本 plan 限定の例外として理由を明記する。

**推奨は、三値 auto の方針を維持し、上記を反映してから実装することです。** 課題は既存実装と重複せず、既定変更として独立しています。面流束式を変えない方針も費用対効果は妥当です。実装前の優先順は、①実行単位の来歴と旧形式互換、②変換用 config と台帳分離、③発火する検証ケースと除外規則、④定量ゲートと既知の救済効果確認です。

指摘数: Critical 0 / Major 7 / Minor 1
