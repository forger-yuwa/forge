# codex レビュー: thermophysics-solver-owned-species-db (plan)

- **plan**: [`plans/active/thermophysics-solver-owned-species-db.md`](../../plans/active/thermophysics-solver-owned-species-db.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `192d0448` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m1
- **extra**: `notes/reviews/2026-09-27-solver-owned-species-db-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

目的の同定と「温度区切りの和集合による起動時合成」は妥当です。前身計画の未解決部分を扱っており、単なる重複ではありません。  
ただし、restart の記録契約、全種表の名前解決、化学種拡散、検証範囲を実装前に補う必要があります。

1. **Major — 解決済み記録を「どの保存場の物性か」に結び付ける契約が不足しています。**

   **根拠:** [plan:85](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:85) は run 内への記録と比較を規定しますが、個々の `res_*.h5`・restart 入力との対応、記録の上書き防止、起動前の宛先物性解決を規定していません。現状の [restart_field.py:36](/home/sano/work/forge/solver_density_cuda/tools/restart_field.py:36) は種署名を確認せず配列をコピーします。ソルバ本体も [main.cpp:1193](/home/sano/work/forge/solver_density_cuda/main.cpp:1193) で `valueFileName` を直接読み込みます。

   `compare_signatures` の単独再実行では、内蔵 N2 の係数を `+0.001` 変更しても不一致一覧は `[]` でした。V1 は必要ですが、比較関数だけ通っても実際の restart を保護したことにはなりません。

   **対案:** 使用した物性の内容ハッシュを各 HDF5 に格納し、不変の解決済み記録と対応させてください。宛先は起動前に同じ resolver で解決し、ソルバ直接読込・同一メッシュ restart・補間・設計 runner の全入口で照合します。V1 に「保存場だけコピー」「記録の取り違え」「外部 DB の変更」「起動前の宛先」を追加してください。輸送に使う構成実種の係数・LJ まで記録対象です。

2. **Major — 平均 LJ を廃止すると、既定の化学種拡散モデルが未定義になります。**

   **根拠:** [plan:91](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:91) は粘性・熱伝導しか扱っていません。しかし既定は `speciesDiffusionMethod = 1`（[solverConfig.hpp:579](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:579)）で、二元拡散係数も LJ を使います（[thermo_d.cuh:456](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:456)）。実際の面カーネルはその係数で種流束を作り、補正流束から `Σh_s J_s` をエネルギーへ加えています（[speciesTransport_d.cu:254](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:254)）。

   実種ごとの拡散速度が異なる場合、固定内部組成の lump はその変化を表現できません。μ・λ の実種展開だけでは、この閉じ方は決まりません。

   **対案:** 初回の lump は**共通 Schmidt 数による拡散に限定し、lump と kinetic 混合平均拡散の併用を入力で拒否**してください。kinetic 拡散が必要なケースは既存の `full` を使います。組成勾配を持つ試験で、種質量収支とエンタルピー拡散を検証する条件を追加してください。

3. **Major — 「CEA 全種表」を既存の名前解決へそのまま投入できません。**

   **根拠:** ローカル `thermo.inp` を現行パーサで読み取ると 2,012 エントリあり、`CO` は MW 28.0101、`Co` は MW 58.9332 という別種でした（[thermo.inp:2593](/home/sano/work/forge/.venv-cea/nasa_cea/thermo.inp:2593)、[同:4135](/home/sano/work/forge/.venv-cea/nasa_cea/thermo.inp:4135)）。一方、Python は種名を大文字化し（[composition.py:35](/home/sano/work/forge/design/forge_design/gas/composition.py:35)）、C++ の index 検索も大小文字を同一視します（[speciesDB.cpp:75](/home/sano/work/forge/solver_density_cuda/input/speciesDB.cpp:75)）。

   また、全種表には `H2O(L)` があり、LJ 表は 23 種分だけです。現行生成器は LJ 不明種へ N2 相当値を仮置きします（[cea_thermo_to_species_db.py:105](/home/sano/work/forge/solver_density_cuda/tools/cea_thermo_to_species_db.py:105)）。この動作を継承すると、「データがある」が「気相 EOS・輸送に使用可能」に化けます。

   **対案:** 大小文字を区別する canonical ID と、明示的な互換 alias を定義してください。相・熱力学の利用可否・輸送データの有無を分け、未対応相と輸送データ欠落は使用時に拒否します。現行 `AIR` は CEA の `Air` と別の互換擬似種として保持してください。N2 等の 6000–20000 K 区間を有効化することも外挿規約の変更なので、移設回帰から分離すべきです。

4. **Major — 新形式への移行対象から、実際に使われる後処理・restart reader が漏れています。**

   **根拠:** [total_quantities.py:119](/home/sano/work/forge/solver_density_cuda/tools/total_quantities.py:119) は `speciesDBFile` がなければ `species_db.yaml` を開き、[_TPGas:48](/home/sano/work/forge/solver_density_cuda/tools/total_quantities.py:48) は二温度域を前提とします。このファイルは plan §7 にありません。種変換器もこの `_TPGas` を import します。

   SERN の [_species_signature:571](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:571) にも独自の種名処理と `species_db.yaml` 必須条件があります。したがって、ソルバの V1–V4 が通っても、DB ファイル廃止後に全温・全圧算出や段階 restart が壊れ得ます。

   **対案:** Python の物性解決・区間評価・署名比較を共通 API に集約し、これらの reader を明示的に移行対象へ追加してください。**新 config で prepare → 段階 restart → `total_quantities` → 種変換**までを、生成 `species_db.yaml` なしで通すことを runner 切替の合格条件にします。

5. **Major — 区間可変化の検証が、実動する float 経路と datum 処理を覆っていません。**

   **根拠:** double の `SpeciesThermo` とは別に、float の `SpeciesThermoF` と専用評価関数があります（[thermo_d.cuh:64](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:64)）。datum は現在 `low[7]` と `high[7]` の二箇所へ焼き込み、その後 float 表へ変換します（[thermo_d.cu:69](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cu:69)）。

   和集合合成自体は支持します。区切りと上下限を変えた合成試験を既存 Python 評価式で検算し、100–20000 K の区間内サンプルで cp・h・s° の尺度正規化最大誤差はそれぞれ **5.98e−16／2.08e−15／6.67e−16** でした。ただし、これは float の面評価・温度反転を検証していません。

   **対案:** V2/V3 の `1e−12` は double の合成・評価に限定し、量別の絶対許容差も付けてください。全区間への datum 適用、境界温度そのものと両側、float mirror、e↔T 往復を別試験にします。温度反転には既存 [test_thermo_float.cpp:108](/home/sano/work/forge/solver_density_cuda/tools/test_thermo_float.cpp:108) の `errHyb/T < 3e−8` 等を継承できます。輸送種数・展開実種数・区間数の上限も分離し、超過を起動時に拒否してください。

6. **Major — 潜熱共通化の項目 #8 は、エンタルピー基準を指定しないまま実装可能な書き方です。**

   **根拠:** [plan:121](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:121) の「種 DB の評価」に対し、実動 DB は sensible 基準へ変更済みです。一方、[h2o_latent:250](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:250) は絶対基準の気相・液相 h を引きます。

   関数単体の数値検算では、H2O の h_abs(298.15 K) は **−13.423291 MJ/kg**、現行 L(300 K) は **2.438293 MJ/kg**。sensible の気相 h と絶対基準の液相 h を引く置換では **15.861584 MJ/kg** となり、現行クランプ後も **3.5 MJ/kg** です。これは将来実装の反例であり、現行ソルバの測定誤差ではありません。

   **対案:** #8 は独立した後続計画へ切り出してください。その計画で、気相の絶対 h 復元、液相との datum 整合、外部 H2O 上書き時の扱い、Python の二相 EOS 複製も含む整合を定義します。onset・g の変化記録だけでなく、datum を変えても L が不変である単体試験を必須にしてください。

7. **Major — V5 は準定常回帰として使えますが、汎用熱物性コア変更の合格条件としては不足しています。**

   **根拠:** `case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX/` の既存区間 CSV を再判定しました。判定区間は `S2_main`、最終 step 23999 です。

   - `check_convergence`: **NOT CONVERGED (stalled/plateau)**。末尾代表値は `rms_ro=2.78e−6`、`rms_roUy=6.72e−4`、`rms_roe=2.70`、種残差も停滞。
   - `check_quasisteady --series-csv`: 六つの報告量は、drift/fluct 閾値を各 0.2% として **ALL STEADY**。

   根拠ファイルは [CONVERGENCE_VERDICT.txt](/home/sano/work/forge/case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX/CONVERGENCE_VERDICT.txt)、[lumpX_series.csv](/home/sano/work/forge/case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX/lumpX_series.csv)、索引は [case README:783](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:783) です。

   また、参照する V0 の三反復と抽出関数は、別 plan の未完了作業です（[campaign plan:154](/home/sano/work/forge/plans/active/tooling-design-problem-campaign-recipe.md:154)）。「NS run 一本で変化量を記録」も合否を定義していません。

   **対案:** V0 の抽出関数・数値許容差・旧バイナリを**数値変更前に固定する依存関係**を追加してください。case/44 は準定常回帰として残し、前身計画の小型 node TP ケースを使う全残差 `PASS` 試験も追加します。NS は対象ケース・有効な輸送設定・保存収支・判定量を事前指定してください。メッシュ品質、IC と datum の整合、段階起動、判定区間も §6 に明記すべきです。

8. **Minor — 現在仕様ドキュメントの更新が実装後に置かれています。**

   **根拠:** [plan:108](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:108) と残作業 #9 は docs を最後に置いていますが、[AGENTS.md:164](/home/sano/work/forge/AGENTS.md:164) は仕様文書と計画を揃えてから実装する順序です。

   **対案:** resolver・記録スキーマ・区間評価・拡散制約の仕様追記を実装前へ移し、完了時の docs 同期とは別作業にしてください。

**推奨は、起動時合成方式を維持した段階移行です。** 上記の優先順で計画を補正し、最初の実装は「保存場に結び付いた解決済み記録と全 restart 経路の照合」に限定してください。そこを通過後に共通 DB、区間可変化、輸送、runner を順に切り替え、潜熱変更は後続計画で扱うのが妥当です。

ファイル変更なし。数値確認は既存判定ツールとメモリ内の物性検算のみで、新規 CFD は実行していません。指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
