# codex レビュー: thermophysics-cea-mole-fraction-species (plan)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `plan`
- **date**: 2026-09-15
- **commit**: `3a065ca8` (feature/sern-design)
- **codex**: effort `high`, 6.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M8/m2
- **focus**: §4.1 の YAML 入力仕様 (composition_basis / tp_species full|lumped|pseudo / species_db) と §4.3 の forge 側 X{s} 入力の妥当性、§6 の full≡lumped 等価性ゲート、見落としている forge 側の 2 種前提 (condGasSpecies, 後処理の Y1 決め打ち)
- **extra**: `design/forge_design/gas/semiperfect.py`, `design/forge_design/evaluate/runner_axismach.py`, `solver_density_cuda/tools/cea_thermo_to_species_db.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **GO-with-changes**

目的は妥当ですが、現計画のまま実装開始は勧めません。  
主な不足は、外部 DB の設計側への反映、種順序変更時の restart、凝縮種指定の整合性、検証ゲートです。  
主要な凝縮カーネルに２種固定は見つかりませんでした。ただし「forge 側は入力とログだけ」というスコープは不足しています。

## 指摘一覧

### 1. **Major — `species_db` を CFD に渡すだけでは、設計と CFD が異なる熱物性になる**

**根拠:** [semiperfect.py:176](/home/sano/work/forge/design/forge_design/gas/semiperfect.py:176) は種の存在・MWを、同ファイルの [cp/h 評価:247](/home/sano/work/forge/design/forge_design/gas/semiperfect.py:247) は係数を、すべて `SPECIES_NASA9` から取得します。[probdef.py:47](/home/sano/work/forge/design/forge_design/probdef.py:47) に外部 DB を渡す経路もありません。

したがって、§4.2 の出力処理だけを追加すると、外部 DB だけにある種は設計時に失敗し、既存種の係数を上書きすると MOC・IC と CFD が食い違います。また「外部 DB をコピー」と「`lumped` の擬似種を生成」の優先関係も未定義です。

**対案:** 問題読込時に、内蔵 DB に外部 DB を上書きした**解決済み DB**を１つ構築し、換算・`GasSemiPerfect`・擬似種生成・IC・出力の全経路へ渡してください。

- 外部パスは問題 YAML の所在基準で解決する。
- run には、使用種と生成擬似種を含む自己完結 DB を出力する。
- MW・係数・温度域・出典を検証する。
- 異なる温度区切りを持つ種の `lumped` は、２区間で厳密に表現できない場合に拒否する。

### 2. **Major — `condensing_species` と `condGasSpecies` の不整合を拒否する仕様がない**

**根拠:** [runner_axismach.py:139](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:139) は、既存 `condGasSpecies` があると自動決定を行いません。計画の `full` 順序は `[H2O,N2,O2,AR,CO2]` なので、旧設定の `condGasSpecies: 1` を残すと **N2 を凝縮対象として参照**します。

一方、[condensationSource_d.cu:30](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cu:30) の凝縮物性は `condModel` から作られ、種名からは決まりません。[dependentVariables_d.cu:128](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:128) には index の直接参照もあります。

**対案:** 種名を正本として、数値 index は生成値に限定してください。明示 index があれば一致を検査し、不一致はエラーにします。加えて、以下を入力段階で拒否してください。

- 凝縮種が `keep` に残っていない `lumped`。
- 混合擬似種 `MIX` を純凝縮種として扱う設定。
- 種名と `condModel` の不一致、範囲外 index。
- 擬似種名と実種名の衝突、重複種、未対応の凝縮種。

単成分時は [variables.cpp:56](/home/sano/work/forge/solver_density_cuda/variables.cpp:56) で `roY` が登録されないため、`full` が１種になった場合も別途定義が必要です。

### 3. **Major — ２種→５種の restart 手順が欠落している**

**根拠:** [runner_axismach.py:482](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:482) は IC 作成後に `interp_field.py` を実行します。同ツールは [68行](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:68) と [88行](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:88) で、`Y0`／`Y1` を種名の照合なしに同番号へコピーします。

旧 `[MIXDRY,H2O]` から計画の順序へ引き継ぐと、旧乾燥成分が H2O、旧 H2O が N2 に入ります。先に作った５種 IC も部分的に上書きされます。§6 の既存プロトコル再利用で踏む経路です。

**対案:** 種名と擬似種内組成に基づく restart 変換を、CFD 回帰より前に追加してください。

- `MIXDRY` の保存量を構成種へ分配する。
- H2O と凝縮モーメントを正しく移す。
- DB・エンタルピー基準が変わる場合は、元の温度から新しい `roe` を構成する。
- `ΣρY=ρ`、総水量、温度の保存を検査する。
- 通常の restart は、種順序・DB の不一致を黙って受け付けない。

### 4. **Major — §4.3 の入力位置・既定値・CSV の仕様が実装構造と合っていない**

**根拠:**

- [solverConfig.cpp:763](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:763): `initial` は組成 mapping ではなく**文字列**です。
- [main.cpp:954](/home/sano/work/forge/solver_density_cuda/main.cpp:954): DB 初期化は `cfg.read()` の後です。「species 読込直後」には解決済み MW がありません。
- [convertGmshToForge.cpp:36](/home/sano/work/forge/solver_density_cuda/mesh/convertGmshToForge.cpp:36): 変換器も境界を読みますが、先行する `thermo_init_db()` はありません。
- [boundaryCond.cpp:139](/home/sano/work/forge/solver_density_cuda/boundaryCond.cpp:139): 未指定 Y は `Y0=1`・他ゼロで補完されます。
- [gen_inlet_profile.py:223](/home/sano/work/forge/solver_density_cuda/tools/gen_inlet_profile.py:223): 分布生成では未指定種を境界の Y 比で補完します。

`X_H2O` が「生成器への入力列」なのか「forge が直接読む CSV 列」なのかも、§4.3 と §10 で確定していません。

**対案:** §4.3 を具体的な入力契約に書き直してください。

- GPU 初期化に依存しない host の DB 解決処理を用意する。
- `initial` の既存文字列形式を維持し、組成を追加する具体的 schema と HDF5 への反映箇所を定義する。
- X は double で検証・換算し、最後に `flow_float` 化する。
- 負値・非有限値・総和ゼロ・未知 index・X/Y 混在を拒否し、省略種の扱いを明記する。
- 今回の CSV は **生成器が `--X`／`X_NAME` を受け、forge 用には `Y{s}` を出力する仕様**に統一する。直接 X 補間は範囲外とする。

### 5. **Major — ParaView のヘルプ変更だけでは、N2 を水蒸気として後処理する**

**根拠:** [forge_filters.py:502](/home/sano/work/forge/solver_density_cuda/tools/paraview/forge_filters.py:502) と [522行](/home/sano/work/forge/solver_density_cuda/tools/paraview/forge_filters.py:522) は既定配列を `Y1` に固定し、[615行](/home/sano/work/forge/solver_density_cuda/tools/paraview/forge_filters.py:615) でそれを蒸気分圧に使用します。

計画の５種順序では `Y1=N2` です。配列は存在するので、警告なく誤った飽和度を出します。§4.4 の「ヘルプ文だけ更新」は目的を満たしません。

なお、[total_quantities.py:85](/home/sano/work/forge/solver_density_cuda/tools/total_quantities.py:85) は既に全種を読み、凝縮警告も `g_0` 基準です。こちらに計画が述べる `Y1` 決め打ちはありません。

**対案:** ParaView に run の設定を明示的に読み込む経路を追加し、種名から配列を解決してください。設定を得られない場合は配列選択を必須にし、`Y1` を自動採用しないこと。H2O を先頭・中間・末尾に置く順序入替試験を追加してください。

### 6. **Major — 「dry なら同一解」は適用条件が広すぎる**

**根拠:** `cp/h` の線形混合は、擬似種内部の組成比が固定なら成立します。しかし [speciesTransport_d.cu:246](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:246) は種ごとの拡散係数を評価し、[260行](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:260) は種エンタルピー拡散をエネルギーへ加えます。

したがって、反応なしでも組成分布・差動拡散・輸送物性の違いがあれば、`full` と `lumped` は一般には同一解になりません。`pseudo` １種は独立した組成輸送自体を持ちません。

**対案:** 等価性を次の条件に限定してください。

> 同じ解決済み DB・温度域処理・エンタルピー基準を使用し、擬似種内部比が空間的に一定の非粘性 frozen ケース。

その上で、熱物性の代数的等価性と、float32 CFD の許容差付き回帰を分けます。`condKantrowitz: 1` の比較方針はこの条件下では妥当です。`2/3` は等価性を要求せず、[carrier 集計:85](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:85) の種別和を独立計算と照合する試験にしてください。

### 7. **Major — CEA 直読みと内蔵転記の「係数一致 → run 同一」は現物で成立しない**

**根拠:** 同梱 `thermo.inp` を既存パーサで読み、`SPECIES_NASA9` と照合しました。

| 項目 | 内蔵転記 | CEA 直読み |
|---|---:|---:|
| H2O MW［kg/mol］ | 0.0180153 | 0.01801528 |
| AR 高温域 `a0` | 0 | 20.10538475 |

根拠は [semiperfect.py:52](/home/sano/work/forge/design/forge_design/gas/semiperfect.py:52)、[thermo.inp:5688](/home/sano/work/forge/.venv-cea/nasa_cea/thermo.inp:5688)、[AR 高温域:656](/home/sano/work/forge/.venv-cea/nasa_cea/thermo.inp:656) です。係数差の大きさが、そのまま熱物性誤差の大きさを意味するわけではありませんが、完全一致という前提は誤りです。

さらに [cea_thermo_to_species_db.py:124](/home/sano/work/forge/solver_density_cuda/tools/cea_thermo_to_species_db.py:124) の `--check` は **N2 の low だけ**を比較し、不一致でも失敗終了しません。

**対案:** §6.5 を全使用種の MW・両温度域・係数の照合へ変更してください。外部 DB の値をその経路の正本として設計・CFD 双方で使い、内蔵との差は物性値の許容差で評価します。「同一 run」を要求するのは、同じ DB を使う入力表現間に限定してください。

### 8. **Major — §6 は未収束の比較元を使い、追加する species 方程式も監視していない**

**根拠:** `check_convergence.py` を再実行し、保存済み判定と一致しました。

| 比較元 run | VERDICT | 代表的な最終残差 |
|---|---|---|
| `case/44.vitiated_air_wt/run_0126_va3_M4.19_Lc8_noneq_rerun_merged/` | `NOT CONVERGED (stalled/plateau)` | `rms_roUy=4.24e-4`、`rms_roe=5.27e-1` |
| `case/44.vitiated_air_wt/run_0131_va3_M4.19_Lc8_noneq_inletTt_cfl05/` | `NOT CONVERGED (stalled/plateau)` | `rms_roUy=3.68e-4`、`rms_roQ0_0=1.33e9` |

両 run の４スナップショットの `VALUE/*` に NaN/Inf はありませんでした。保存済みメッシュ品質は `VERDICT: PASS`。`check_quasisteady.py --quantity machmax,pmax` は両者 `STEADY` ですが、**これは onset・出口平均 g・質量流量の定常性を判定していません**。索引は [case README の run 一覧](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:566) です。

また [main.cpp:158](/home/sano/work/forge/solver_density_cuda/main.cpp:158) の残差列には species が追加されず、実際の CSV にも `rms_roY*` はありません。さらに node のみの計画は [検証手順:40](/home/sano/work/forge/procedures/verification/README.md:40) の共有境界コードに対する両離散化検証を満たしません。

**対案:** §6 に以下を合否条件として追加してください。

- 収束済みの小型 TP ケースを node/cell 両方に設け、X/Y・種順序・旧入力の回帰を行う。
- 全 `roY{s}` 残差を出力して `check_convergence.py` に通し、`ΣY`・負値も検査する。
- onset の閾値、出口 g の重み、質量流量の断面・相対誤差を定義し、その時系列を `check_quasisteady.py --series-csv` で判定する。
- case/44 の現行２本は未収束の回帰参考とし、収束済み解の一致の根拠にしない。
- [凝縮ソース律速の別 plan](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:16) とのコード世代を固定する。同じ CFL だけでなく、両経路の `condLim` を確認する。

### 9. **Minor — 既存機能との重複と SERN の対象範囲を整理すべき**

**根拠:** [frozen.py:24](/home/sano/work/forge/design/forge_design/gas/frozen.py:24) に既に `mole_to_mass()` があり、SERN はモル分率入力に対応しています。一方、[runner_sern.py:87](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:87) は `[EXH,AIR]` 固定で、[281行](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:281) の IC も２領域・２擬似種です。axismach の３ヘルパを一般化して済む共通経路ではありません。

**対案:** 今回の対象を `gas.model: semiperfect` の設計経路と明記し、既存 `frozen.mole_to_mass()` は共通換算関数へ委譲してください。SERN の `full` 化は、排気・外気の種集合を統合する別仕様として切り出すのが妥当です。CEA 変換器も既存機能として再利用します。

### 10. **Minor — 単体試験の絶対誤差とビット同一条件が過剰**

**根拠:** 計画例を既存関数で計算すると、`ΣX=0.99882478`、`Y_H2O=0.03769539643469918` で、換算例は正しいです。一方、200–3000 K を 1 K 刻みで評価した既存 full／split の最大絶対差は、

- `cp`: `1.59e-12 J/(kg·K)`
- `h`: `2.33e-9 J/kg`

でした。§6 の `1e-10` を絶対誤差と解釈すると、正しい既存混合も失敗します。出力は [runner_axismach.py:171](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:171) で８桁丸めされ、DB 生成には浮動小数点の加算順序も影響します。

**対案:** `rtol/atol` と単位を明記してください。例えば double の熱物性比較は `rtol=1e-12` に量別の絶対許容差を併用し、`h(Tref)=0` 近傍も検査します。YAML は解析後の値を比較し、ビット同一は同一の正規化済み入力を再出力する決定性試験に限定してください。

## 推奨

**解決済み DB・正規化済み組成・species 順序を一体で管理する共通処理を先に設計し、それを各入力・IC・後処理から使う方針に修正して進める**ことを推奨します。

実装前の修正優先順は次のとおりです。

1. §4.1–4.2：共通 DB、凝縮種の整合性、種変換 restart を確定する。
2. §4.3–4.4：実在する読込経路に合わせた X 入力仕様と、ParaView の配列解決を確定する。
3. §6：等価性の適用条件、CEA 差分、収束・species 残差・派生量定常性・node/cell のゲートを書き直す。
4. §5.1：上記の単体試験を CFD 回帰より前へ置き、SERN の範囲を整理する。

ファイル変更は行っていません。以上の提案は **plan 未反映**です。

指摘数: Critical 0 / Major 8 / Minor 2
