# codex レビュー: tooling-convergence-and-wall-resolution-gates (plan)

- **plan**: [`plans/active/tooling-convergence-and-wall-resolution-gates.md`](../../plans/active/tooling-convergence-and-wall-resolution-gates.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `be527007` (feature/sern-design)
- **codex**: effort `high`, 7.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M8/m1
- **focus**: case/49 の result レビュー (同日) の M3 (別 BC の履歴を連結した収束判定)・M4 (準定常ツールが欠損/NaN を通す)・M9 (y+ が第一内部ノードでなく ghost 間距離 dcc から出ている) を、リポジトリ全体の恒久対応に格上げする plan。**web 検索は不要、リポジトリのコードと既存 run だけで判断してほしい**。見てほしい点: (1) §4.1 の check_convergence.py 段階認識の設計が妥当か。既存の呼び出し側を壊さないか。既定挙動を変えるべきか明示オプション必須にすべきか。(2) §4.2 の check_wall_resolution.py の y1 の取り方 (場の wall_dist の正の最小値の下位パーセンタイル) が頑健か。node の双対重心ベース wall_dist の既知の癖を踏まえてより正しい y1 の定義があれば示してほしい。判定閾値の妥当性も。(3) §4.5 でソルバの ypls の定義を変えない判断 (mode 1/2 が同じ配列を使うため) が正しいか。node 用に正しい y1+ をソルバで出すべきならどこにどう実装するか (壁カーネルからは wall_dist[ic]=0 で第一内部ノードが直接引けない制約がある)。(4) 見落としている同種の甘い判定が repo に他にもないか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

課題の恒久対応は必要ですが、現案では異なる設定の履歴混合を防げず、壁解像も誤判定します。
特に「全域の `wall_dist` 下位値を第一層厚にする」設計は、既存runで誤りを確認しました。
以下を実装前に修正してください。

ファイルは変更していません。ローカルに `run_0101`〜`run_0103` はなく、計画の **0.146／2.60** は原データで検証できませんでした。以下の壁面数値は保存スナップショットの診断値であり、定常値という主張ではありません。

1. **Major — §4.1のプレフィックス指定と共通列化では、異なる問題の連結を防げない**

   **根拠:** [plan:65](/home/sano/work/forge/plans/active/tooling-convergence-and-wall-resolution-gates.md:65) は開始プレフィックスを指定するだけで、方程式・BC・空間離散化の同一性を検証しません。さらに68行の「共通列のみ」は、層流とSSTを混ぜた際に `rms_roK`／`rms_roOmega` を落とし、その問題を隠します。

   実際の生成側は段ごとに設定ファイルを上書きし、退避するのは残差・画像・VERDICT・ログだけです。[gen_runs.py:296](/home/sano/work/forge/case/49.plate_annular_cavity_m5/gen_runs.py:296)  
   現在の `solverConfig.yaml` から過去全段の設定を確定することはできません。

   **対案:** 各段の実効設定、BC・入口分布、メッシュ、方程式列、残差定義、実行順、restart元を保存する段階manifestを先に設計してください。連結は、それらの互換性を確認した連続区間に限定する。**必須残差列の欠損は拒否**し、共通列化は起動診断専用にするべきです。`--from-floor` の参照側にも同じ区間検査を適用してください。

   CFL変更は直ちに別問題にはなりません。「最後の数値設定変更」ではなく、離散方程式と残差尺度の互換性で境界を決めるべきです。

2. **Major — 本段単独の一律拒否は過剰で、CLIだけの変更では既存経路を守れない**

   **根拠:** 本段が同一設定の区間なら、それを検査すること自体は正当です。低下桁数不足という結果を受け入れればよく、「必ずplateauだから拒否」は原因と対処が逆です。

   正本ツールを再実行した結果、次の3runはいずれも全保存量判定で  
   **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`** でした。

   - `case/49.plate_annular_cavity_m5/run_0003_stageB_cpg/`
   - `case/49.plate_annular_cavity_m5/run_0004_stageC_cpg/`
   - `case/49.plate_annular_cavity_m5/run_0005_stageD_cpg/`

   呼び出し側にも差があります。[run_case.sh:23](/home/sano/work/forge/solver_density_cuda/tools/run_case.sh:23) は無指定でCLIを実行し、判定の終了コードを捨てます。一方、[sern_gates.py:83](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:83) は `analyze()` を直接呼びます。CLIだけに段階拒否を入れると、同じrunの扱いが分かれます。

   ファイル名も `case/16` は `soft/mid`、`case/49` は `S0…S5`。`_all`／`_turb` という連結済みファイルまで存在します。

   **対案:** 推奨する既定挙動は、**単一区間は従来どおり検査し、複数段では明示的な区間選択を必須にする。本段単独も選択可能にする**、です。ファイル名から自動連結しないこと。CLIとPython呼び出し側が共有する区間検証APIを設け、wrapper・設計ゲートも移行対象に入れてください。

3. **Major — §4.2の `y₁` は局所的な第一内部ノード距離ではない**

   **根拠:** 実データで次を確認しました。

   | 対象 | 距離 |
   |---|---:|
   | `case/49.plate_annular_cavity_m5/run_0006_hex_cpg/` の底面第一内部層、側壁から離れた363点 | 約 **80.0006 µm** |
   | 同runのキャビティ内 `wall_dist>1e-9` の下位0.1% | 約 **19.9984 µm** |

   計画の方法を底面へ適用すると、距離を約**4分の1**にします。構造化ヘキサでも、全壁の第一層厚が同じとは限りません。

   また、現行の [calcWallDistance_kdtree.cpp:133](/home/sano/work/forge/solver_density_cuda/input/calcWallDistance_kdtree.cpp:133) は利用可能ならnode座標を使います。「常に双対重心基準」という理解は現行コードには当たりません。ただし計算しているのは最近接**壁点**距離で、壁面への法線距離ではありません。古い値は [variables.cpp:698](/home/sano/work/forge/solver_density_cuda/variables.cpp:698) でそのまま読み込まれます。

   **対案:** 入力メッシュの接続・境界ID・node座標を読み、壁パッチごとに第一内部層を対応付けて、
   \[
   y_{1,b}=(\mathbf x_{I_b}-\mathbf x_{W_b})\cdot\mathbf n_{\mathrm{in},b}
   \]
   を求めてください。構造層では層接続を優先し、非構造では隣接候補と法線整合性を検査する。角部は壁パッチごとに扱い、対応不能・負距離・曖昧な対応は合格にしないこと。

   既存の [ransWallFunction_d.cu:136](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransWallFunction_d.cu:136) にCSR隣接から代表内部点を選ぶ実装があります。ただし曲面・角部で完全ではないことも [accepted plan:62](/home/sano/work/forge/plans/accepted/turbulence-node-wf-representative-point.md:62) に記録済みです。無条件に流用するのではなく、対応の妥当性を出力してください。

4. **Major — 距離以外の `uτ`・物性定義も不整合で、Sutherland固定ではリポジトリ共通にならない**

   **根拠:** [viscousFlux_d.cu:557](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:557) の `utau` は、接線成分へ射影する前のtraction全体から計算されます。その後、[1018行](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:1018) でnodeの `twall_*` が別カーネルにより上書きされますが、`utau` は更新されません。

   `case/49.plate_annular_cavity_m5/run_0003_stageB_cpg/res_cav_floor_9_24000.h5` では、水平底面について
   \[
   \frac{\sqrt{twall_x^2+twall_y^2}}{\rho\,utau^2}
   \]
   の中央値が **2.352**、最大が **14.33** でした。保存された `utau` と接線壁応力は同じ定義になっていません。

   また輸送モデルは一定粘性・Sutherland・混合気のkinetic theoryに分岐します。[gasProperties_d.cu:52](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:52)  
   `Ts` だけからSutherlandで再計算すると、`viscMethod: 0/2` を誤評価します。

   **対案:** 壁単位の定義を
   \[
   \tau_{w,t}=\|(\mathbf I-\mathbf n\mathbf n^\mathsf T)\mathbf t_w\|,\quad
   u_{\tau,w}=\sqrt{\tau_{w,t}/\rho_w},\quad
   y_1^+=y_1\sqrt{\rho_w\tau_{w,t}}/\mu_w
   \]
   と明記し、応力の取得経路を固定してください。物性は対応する壁状態の `vis_lam`／実際の輸送モデルに従う。必要データがない旧出力を推測で合格にしないこと。既存 `utau`・`ypls` は比較診断として残すべきです。

5. **Major — 壁解像の閾値とモード分類が、既存仕様を正しく表していない**

   **根拠:** §4.2の「平均≤1、最大≤5」は、局所的な解像不足を平均で埋め合わせます。例えば面積の90%が0.5、10%が5でも平均0.95でPASSです。この条件が熱流束や摩擦精度を保証する裏付けは計画にありません。

   また `wallTreatmentSST: 1` は純粋な対数則専用ではなく、粘性低層・バッファ層・対数層を接続するautomatic処理です。[theory.md:308](/home/sano/work/forge/methods/turbulence/theory.md:308)  
   `y₁⁺<30` を一律FAILにするのは仕様と矛盾します。現行推奨の粗い壁関数メッシュも30〜80であり、300までの無条件PASSを裏付けません。[recommended-settings.md:106](/home/sano/work/forge/procedures/recommended-settings.md:106)

   **対案:** 低Reの既定スクリーニングは局所 `y₁⁺≤1` とし、緩和は対象量の壁法線細分化感度で裏付けて明示してください。壁群別の面積加重平均・分位値・最大・超過面積率を出す。軸対称では回転面積を使うこと。

   automatic処理は適用レジームと検証済み条件を区別し、30未満を機械的に不合格にしない。内部mode 2はWMLESで、`wallTreatmentSST: 2` ではありません。[implementation.md:503](/home/sano/work/forge/methods/turbulence/implementation.md:503) 壁ごとの実効モデルを判定し、未対応モデルは未評価として返してください。`y₁⁺` 合格だけで熱的・空間的な十分解像を宣言してはいけません。

6. **Major — §4.3は既修正の入口を対象にし、実際に残るNaN通過経路を外している**

   **根拠:** [check_quasisteady.py:298](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:298) のCSV経路は、既に欠損列を拒否し、`classify_series()` で非有限値を拒否します。既存の [回帰試験:64](/home/sano/work/forge/design/tests/run_sern_gates_tests.py:64) もあります。

   書き込みを伴わない再現結果は次のとおりです。

   | 入力・経路 | 実際の判定 |
   |---|---|
   | `[1,1,1,1,1,NaN]` → `classify()` | **STEADY** |
   | 同系列 → `classify_series()`／CSV経路 | **NONFINITE** |
   | CSVの指定列欠損 | **ERROR** |
   | 値は一定、`step` にNaNまたは末尾重複 → `classify_series()` | **STEADY** |

   既存runでも、`case/49.plate_annular_cavity_m5/run_0005_stageD_cpg/cavity_series.csv` の `h_ref` は正本CLIで **`NONFINITE — 12/12 non-finite value(s)`** でした。

   一方、HDF5経路は抽出失敗をNaNにしてから、NaNを落とす `classify()` を呼びます。[check_quasisteady.py:364](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:364)  
   適用不能な要求量をskipし、評価対象ゼロでもSTEADYになる経路も残ります。

   **対案:** §4.3を**全入口で共有する厳格な判定関数への統一**に変更してください。値・時刻の有限性、時刻の順序と重複、空の要求量、抽出失敗、必須量欠損を検査する。絶対許容は単位付きで、相対許容との合成規則と振幅への適用も定義してください。CSV修正を新規成果として扱うのは不適切です。

7. **Major — 収束ツールにも、段階認識とは独立した「検査不足でPASS」が残る**

   **根拠:** [check_convergence.py:39](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:39) は存在する列だけを読みます。メモリ上のCSVで以下を再現しました。

   - `step,phase` だけの2行：`report={}`、**`ok=True`**。
   - `rms_ro` だけが4桁低下：他の保存量がなくても **`ok=True`**。
   - NS全5列を「99点は1、最後の1点だけ `1e-6`」にする：全列 **`drop=6.0dec flat`、`ok=True`**。

   最後の例は、低下桁数を末尾窓ではなく最終1点で決めるためです。[check_convergence.py:152](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:152)  
   また、`outer_end` がある場合、他phaseのNaNはこの読込経路から除外されます。

   **対案:** 設定から期待される残差列を検証し、空の検査を拒否する。履歴全体の健全性検査と、収束トレンドに使うphase選択を分けてください。低下量は末尾窓での持続を要求し、最終1点のゼロや小値だけでPASSにしないこと。

   必須列・全phase検査には既存の [passive_gate_common.py:311](/home/sano/work/forge/solver_density_cuda/tools/passive_gate_common.py:311) が参考になります。これらを修正する以上、§6.2の「非段階runの判定は変わらない」は、**正当な既存入力の互換性を保ち、誤PASSは変更する**に改める必要があります。

8. **Major — §5・§6に、正しさを独立に検証する試験と移行範囲が不足している**

   **根拠:** [plan:122](/home/sano/work/forge/plans/active/tooling-convergence-and-wall-resolution-gates.md:122) の主な正解値は、同じ下位パーセンタイル法を使う [cavity_eval.py:270](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:270) に依存しています。同じ値の再現は、その定義が正しい証明になりません。対象 `run_0103` も今回の作業ツリーにはありません。

   既存計画では、準定常ツールは `accepted/tooling-quasisteady-check.md`、代表内部点診断は `accepted/turbulence-node-wf-representative-point.md` にあり、SERN側でもCSV厳格化と残差判定変更が進んでいます。[SERN計画:820](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:820)  
   これらとの責務分担・呼び出し側移行が§5.1にありません。

   **対案:** 実装順を「データ契約・現在仕様更新 → 判定器 → 呼び出し側移行 → 実測回帰」に変更し、少なくとも次を合否表へ入れてください。

   - 同一設定のrestart連結と、同じ列を持つ別BCの連結拒否。
   - SST・化学種・凝縮列の欠損、派生CSVの二重取込み、参照床の不適合。
   - 解析的な距離が分かる平板、壁群ごとに厚さが異なる格子、曲面・角部・軸対称。
   - NaN・Inf・欠損・不正時刻を、CLI／Python／HDF5各経路で検査。
   - 冷却壁の標準ケース `case/48`、既存automatic処理の `case/26` を使う後処理回帰。

   各試験の期待VERDICT、数値許容、評価壁とstepを固定してください。今回は後処理中心なので、新規CFD投入を必須にする必要はありません。cellの定義分岐は設計に残しつつ、実計算の検証は [現行のnode限定規則](/home/sano/work/forge/procedures/verification/README.md:40) に合わせるのが妥当です。

9. **Minor — §4.5の「定義を変えない」結論は妥当だが、理由と説明が不正確**

   **根拠:** [viscousFlux_d.cu:570](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:570) は既にmode 0だけで `ypls_b` を更新します。mode 1の `ypls` は代表内部点の距離と物性を使います。[ransWallFunction_d.cu:187](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransWallFunction_d.cu:187)  
   配列を共有しているだけで、mode 0限定の変更が必ずmode 1/2へ影響するわけではありません。またdcc基準式の分母は正確には `mu_total` です。

   **対案:** 既存 `ypls` を維持し、モード別の距離・密度・粘性の定義を文書化する方針を推奨します。今回、正しい診断のためにソルバ変更を必須にする必要はありません。

   将来ソルバ出力へ入れるなら、既存値の上書きではなく別配列 `y1plus` を追加する。メッシュ初期化時に壁パッチから第一内部点への対応を作り、最終的な壁応力評価後に出力専用カーネルで計算すれば、`wall_dist[ic]=0` の制約を避けられます。

**推奨は、このplanを修正してから、既存ソルバの解を変えない共通後処理・判定基盤として実装することです。** 実装前の優先順は以下です。

1. §4.1を、設定の出自・必須列・共有APIで区間を検証する設計へ変更する。
2. §4.2を、壁パッチ別の実際の第一内部点距離・整合した接線応力・実物性へ変更する。
3. §4.3の修正対象を全入口へ広げ、指摘7の誤PASSも残作業に入れる。
4. モード別の判定条件と例外の根拠を定義する。
5. 独立した正解値を持つ試験、既存呼び出し側の移行、利用可能なrunを§5・§6に明記する。

ファイル変更禁止に従い、これらの提案は **plan未反映**です。

指摘数: Critical 0 / Major 8 / Minor 1
