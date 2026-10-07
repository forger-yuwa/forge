# codex レビュー: architecture-solver-host-memory (plan)

- **plan**: [`plans/active/architecture-solver-host-memory.md`](../../plans/active/architecture-solver-host-memory.md)
- **stage**: `plan`
- **date**: 2026-10-07
- **commit**: `7259d840` (feature/sern-design)
- **codex**: effort `high`, 3.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
削減対象の同定と R1・R2 は妥当です。R3 はホスト参照の保護と初期化順序を補強する必要があります。  
検証計画には未定義値の比較が含まれ、計測用の最小格子も実際には発散しています。現状の §6 では受入判定できません。

目的は未解決の課題を扱っています。既存の `architecture-residual-monitor-async` は同期削減、起点 plan の B4-5 は変換器の省メモリ化であり、本件との重複ではありません。`matrix::initMatrix` は自身のコンテナを構築する処理で、GPU の `block-DPLUR` 行列とは別です。デバイス側の確保・初期化・転送内容を維持する方針自体に、数値的な問題は見当たりません。

1. **Major — 転送時のガードでは、ホスト書込み前の未確保アクセスを防げない。**

   **根拠:** [`main.cpp:1148`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1148) の checkpoint 復元は、`var.c.at(nm)` を取得し、`v[i]` に書いてから H2D を呼びます。H から履歴変数が漏れると、計画の転送ガードへ到達する前に範囲外書込みになります。`map::at()` はキーの存在しか検査せず、空 vector を拒否しません。`readValueHDF5` に追加するガードも、この別経路には効きません。

   **対案:** GPU 経路のホスト参照を、変数名・期待長を検査する共通アクセサ、または各処理入口の検査に通してください。検査は添字アクセス・イテレータ演算・生ポインタ引渡しより前です。負例も「出力用変数を外す」だけでなく、**checkpoint 復元用変数を外して書込み前に停止する**試験を含めます。

2. **Major — H の確定時点と出力依存関係の解決方法が未設計。**

   **根拠:** 確保は [`main.cpp:1715`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1715) ですが、`FORGE_OUT_RESIDUALS`・`FORGE_RESID_SNAP` による出力登録は初期化終了後の [`main.cpp:3417`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:3417) です。また、[`output.cpp:117`](/home/sano/work/forge-sern-design/solver_density_cuda/output/output.cpp:117) の `h0` は出力名そのものではなく `Ht`・`k` を必要とし、同ファイル 217 行以降の `/CHECKPOINT` 履歴は通常の出力名リストとは別です。

   現在の出力リストを確保時に一度読むだけでは、これらを網羅できません。別実装で名前を再現すると、今後の出力変更で再び食い違います。

   **対案:** **変数・診断登録 → 出力と checkpoint の依存名確定 → H 構築 → 確保**という順序を §4・§5 に明記してください。出力側と確保側で同じ依存名生成関数を使用し、監査表には名前・有効条件・最初の利用箇所を記録します。`output/output.cpp` の変更は「必要なら」ではなく、この共通化の対象として扱うべきです。

3. **Major — 「同じ入力」で回帰する最小計測ケースが、既に発散している。**

   **根拠:** 保存された [`time_v_s050.txt:24`](/home/sano/work/forge-sern-design/notes/investigations/2026-10-07-forge-memlog/time_v_s050.txt:24) は `Exit status: 1` です。実体の [`mem_s050/forge_run.log:424`](/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad/forgemem/mem_s050/forge_run.log:424) には、step 9 の `roe` 非有限値検出があります。読取り検査でも `res_nan_9.h5` の `P`・`roe`・`sonic` に各 **4 個**の非有限値を確認しました。

   `check_convergence.py` を実行した結果は、3 ケースとも **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。`s050` の CSV は step 8 までなので、ツールはその後の発散を観測できていません。これは短いメモリ測定であり、プラトーの原因診断にも使えません。

   **対案:** 現在の記録は「発散前までのメモリ観測」と明記して保持し、数百 step の回帰基準から外してください。同じ本段設定に対し、妥当な初期場・段階起動を用意して正常終了する基準入力を固定します。入力、実行環境、バイナリ、全ログを恒久的な `run_*` と case README に残すことも必要です。

4. **Major — `res_0` の「全データセットがビット一致」は成立を保証できない。**

   **根拠:** [`variables.hpp:276`](/home/sano/work/forge-sern-design/solver_density_cuda/variables.hpp:276) では level 2 の出力に `res_roK`・`res_roOmega`・平均流残差が含まれます。一方、[`variables.cpp:365`](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:365) の確保時ゼロ初期化対象には平均流残差が含まれません。初期出力は [`main.cpp:3526`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:3526) で、通常の残差組立より先です。

   したがって、level 2 の初期出力には値がまだ定義されていない配列が含まれます。そのビット列を合格条件にすると、正しい変更も不合格になり得ます。

   **対案:** 初期出力では初期化済みの保存量・原始量・幾何量を比較し、残差や補正量は最初に定義される組立・更新後に比較してください。比較時点を変数ごとに固定し、未定義値を扱う既存問題は別途記録します。出力互換性にはデータセット集合・shape・dtype・属性も含め、再実行ばらつきの評価指標と許容式を事前に決めてください。

5. **Major — 構成名の列挙だけでは、変更対象の分岐と受入条件を覆えていない。**

   **根拠:** [`procedures/verification/README.md:23`](/home/sano/work/forge-sern-design/procedures/verification/README.md:23) は `case/08.bump` を node の二次精度回帰基準に使えないとしています。候補の `case/63`・`64` もこの checkout にはありません。さらに、計画 §6 は `restart`・`probe` を挙げるだけで、前述の **dual-time checkpoint 復元、環境変数診断、遷移モデル、実際の probe／CHT 出力**の検証条件を特定していません。収束・準定常の判定手順もありません。

   **対案:** 監査表の到達分岐に対応する、具体的な入力と確認成果物の表へ置き換えてください。2D node は標準の `case/36`、軸対称・多成分・凝縮は `case/44`、CHT は実在する `case/52.conjugate_slab` 等から選定します。特に dual-time は **連続実行と checkpoint 経由の分割実行**を比較し、履歴属性も確認します。

   短時間試験の合格は「回帰差が許容内」に限定し、全保存量・出力の NaN/Inf 検査を必須化してください。定常結果の「収束・一致」には両側の `check_convergence.py`、派生量の定常性には `check_quasisteady.py` の VERDICT を要求します。TGV のような非定常ケースには定常収束を要求せず、同じ物理時刻で比較します。

6. **Minor — 一時メモリを定数扱いしており、ピーク移動先の断定が早い。**

   **根拠:** [`memlog_table.txt:17`](/home/sano/work/forge-sern-design/notes/investigations/2026-10-07-forge-memlog/memlog_table.txt:17) の構造量設定後の `HWM − RSS` は **40／98／172 MiB**。再計算すると約 **148 B/節点 + 3.6 MiB**で、固定の「約150 MB」ではありません。実装も [`variables.cpp:567`](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:567) で面数・セル数に比例する作業配列を確保しています。

   同表の構造量設定時の傾き 2751 から計画の削減量 312 + 195 + 850 を引いても **1394 B/節点**であり、`readMesh` の 1365 より大きいままです。単純な差引きモデルでも、ピークが必ず `readMesh` に移るとは言えません。

   **対案:** 工程別の「傾き＋切片」の最大値として見積もり直してください。1000 万節点は暫定値とし、正常終了時の実測と共有環境の実効メモリ予算から更新します。≤1500 B/節点という改善目標自体は維持して構いません。

**推奨は、R1–R3 の方針を維持し、監査と受入条件を先に確定してから段階実装することです。** 実装前の修正優先順は次のとおりです。

1. 指摘 1・2：ホストアクセス保護、H の依存名生成と初期化順序を確定する。
2. 指摘 3・4・5：正常な基準入力、比較対象・時点、分岐別の試験と判定条件を固定する。
3. 指摘 6：工程別メモリモデルを修正し、R1 → R2 → R3 の各段階で効果を確認する。

ファイル変更・計算投入は行っていません。以上はレビュー上の修正提案であり、**plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
