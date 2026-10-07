# codex 諮問 (diagnose): hostmem-h2d-channel

- **brief**: [`notes/reviews/briefs/2026-10-07-hostmem-h2d-channel.md`](../../notes/reviews/briefs/2026-10-07-hostmem-h2d-channel.md)
- **plan**: [`plans/active/architecture-solver-host-memory.md`](../../plans/active/architecture-solver-host-memory.md)
- **date**: 2026-10-07
- **commit**: `415626bd` (feature/sern-design)
- **codex**: effort `xhigh`, 5.8 min, rc=0
- **結論**: **`done` は保留し、`c44dual_ckpt100` だけで、共通の GPU 読戻しデータを与えたホスト処理が同じ CUDA 操作列を生成するかを比較する A/B を先に行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 諮問事項 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 入力経路から R1–R3 を検証する | **採用** | `9c9f623c..93e55957` に CUDA ソースの差分はなく、`c_d`・`p_d` の確保長も維持されている（[variables.cpp:315](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:315)、[同:348](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:348)）。ホスト側の変更を切り分ける方向は妥当。ただし、以下の修正が必要。 |
| H2D・確保・カーネル名/grid/block で経路を網羅できる | **却下・Major** | 実際に、定数への転送（[convectiveFlux_d.cu:186](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:186)）、`cudaMemset`（[variables.cpp:326](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:326)）、D2D（[main.cpp:2719](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2719)）、`dt`・係数のカーネル引数渡し（[timeIntegration_d.cu:1420](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1420)）がある。**転送先・値・全引数・stream/context・同期依存関係まで含む操作列**にする。 |
| base 2 本で変わる H2D はデバイス由来として内容を免除する | **却下・Major** | 反復間の差からデータの由来は判定できない。また D2H 後にホストで処理して H2D する実経路がある（[calcGradient_d.cu:737](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:737)、[同:799](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:799)）。例えば base が読戻し値 `r`、new が `2r` を送っても、提案条件では両方とも免除される。**由来をコードで登録し、共通の読戻しデータに対するホスト処理を比較する。** |
| 同じ確保列と生バイトのハッシュで入力を同一視できる | **却下・Major** | H2D にはデバイスポインタの配列も含まれる（[boundaryCond_d.cu:911](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/boundaryCond_d.cu:911)）。さらに確保直後に初期化される配列は限定される（[variables.cpp:315](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:315)）。**ポインタを「確保 ID・世代・offset」で対応付け、読み取る領域の初期化を別途確認する。** 同じ確保要求は、同じアドレス・初期内容の保証ではない。 |
| 初期化＋1〜2 step の一致から、既存の差を非決定性と確定して `done` にする | **却下・Major** | 対象入力は100 stepで出力する（[solverConfig.yaml:18](/home/sano/work/forge-sern-design/case/66.hostmem_regression/inputs/c44dual_ckpt100/solverConfig.yaml:18)）。checkpoint のホスト処理も比較範囲に必要（[output.cpp:173](/home/sano/work/forge-sern-design/solver_density_cuda/output/output.cpp:173)）。また、同じ操作列でもホストの投入タイミングが変われば、浮動小数点 `atomicAdd` の実行順序・結果の出現頻度まで同じとは限らない。**観測した経路の同等性と、通常実行の結果分布を区別する。** |

結論: **`done` は保留し、`c44dual_ckpt100` だけで、共通の GPU 読戻しデータを与えたホスト処理が同じ CUDA 操作列を生成するかを比較する A/B を先に行う。**

第 1 仮説: **R1–R3 は、同じ GPU 読戻しデータを与えれば、GPU に渡す意味のある値と操作を変えていない。** 確度: **中**
  
根拠: 実差分では CUDA ソースと主要デバイス確保長が維持され、ホスト参照には長さ検査が追加されている。ただし、実際にロードされたデバイスコード・全引数・全転送内容の同等性は未確認である。
  
反証条件: 同じ外部入力・同じ読戻しデータに対して、最初の異なる転送値、転送先、カーネル引数、確保寿命、または操作順序が検出されること。

第 2 仮説: **ホスト処理の変更による差が、提案した記録項目または D2H 由来の免除に隠れる。** 確度: **低・未確認**。上表の未捕捉経路があるため除外できない。

第 3 仮説: **未初期化領域の読取り、または投入タイミングの変化が実行間変動に影響する。** 確度: **低・未確認**。`atomicAdd` の存在（[calcGradient_d.cu:251](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:251)）だけでは、今回の差の原因を確定できない。

判別 A/B: **比較因子は base/new のホスト実装だけ。対象は同じ `c44dual_ckpt100`、初期化から100 step・checkpoint 書出しまで。**

まず base の読戻しデータを1系列記録し、同じ検証用経路で base/new に再供給する。D2H は**完了後・ホスト利用前**に共通化する。GPU ポインタを含むデータは各実行の対応先へ変換する。これはホスト処理の試験なので、再供給中の計算場を物理解として評価しない。

- **A：全操作が一致** → 「この共通入力系列で、R1–R3 がホスト処理を通じて GPU 入力を変えた」という第2仮説を退ける。保証は通過した経路に限る。
- **B：操作が不一致** → 第1仮説を退け、最初の不一致を調べる。同じ読戻しを与えているため、その不一致を GPU の実行間変動では説明できない。
- 未対応 API・読戻しの対応失敗・記録欠落は**判定不能**。免除して A にしない。

事前登録する合格条件の文案は次とする。

> 固定した入力、環境、GPU、ビルド条件、および共通のデバイス読戻し系列に対し、初期化から指定した終了点まで、正規化した CUDA 操作列の不一致が0件である。比較対象は確保・解放と寿命、全方向の転送と範囲、symbol 更新、memset の値、カーネルの全引数・grid/block・shared memory、stream/context と同期依存関係を含む。D2H 由来の値を内容比較から一括除外しない。未対応操作・由来不明データ・記録欠落は合格にしない。合格の意味は「当該入力系列でのホスト側操作生成の同等性」に限定する。

D2H の扱いは、単純退避・復元なら**同一 run 内で読戻した値と再送値のビット一致**を要求する。実例は [main.cpp:2311](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2311) → [同:2325](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2325)。変換を伴う場合は共通入力で変換結果を比較する。読戻し値がカーネル引数や分岐へ流れる経路も対象にする。例えば適応 `dt` はこの形である（[setDT_d.cu:414](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/setDT_d.cu:414)）。今回の `c44dual` は `control: 0` なので、この例を今回の原因とは扱わない。

記録は **CUPTI の Runtime/Driver callback を基礎にする**。直接呼出しを個別にラップする方式より漏れを管理しやすい。ただし CUPTI が変数名・引数型・データ由来を自動で理解するわけではなく、その対応表は必要である。callback 引数は有効期間内に複製し、callback 内から追加の CUDA 転送を呼んでハッシュを取る設計は避ける。[NVIDIA Callback API](https://docs.nvidia.com/cupti/main/main.html#cupti-callback-api)、[CallbackData の寿命](https://docs.nvidia.com/cupti/api/structCUpti__CallbackData.html)。

`nsys` は操作数・stream・転送サイズ等の照合には使えるが、通常のトレースは転送内容や全カーネル引数の比較を代替しない。[NVIDIA CUDA GPU Trace](https://docs.nvidia.com/nsight-systems/AnalysisGuide/index.html#cuda-gpu-trace)。

やらない方がよいこと: **現案のまま30構成へ展開すること、base 2 本の不一致を根拠に内容比較を免除すること、全配列をゼロ初期化して一致させること、トレース下の合格を通常実行の頻度差解消と扱うこと。** 64bit ハッシュ一致も厳密なバイト一致の証明ではない。厳密一致を完了根拠にする箇所は実データで照合する。

呼び出し側の前提への異議: **「同じ入力・操作」という前提自体が、提案した記録では確認できない。** また、確認できても「§6.2・§6.3 の差は変更起因ではない」という因果の断定までは進めない。

一方、**用途に基づく許容差が、あらゆる受入れ方法で必須というわけではない。** デバイスコード、意味のある入力、ホスト処理、初期化・同期、および出力処理の同等性を十分に立証できれば、数値演算を変えない変更として受け入れる道はある。今回の省略付きトレースと1構成の A/B は、その立証を完成させない。現行 [plan §8](/home/sano/work/forge-sern-design/plans/active/architecture-solver-host-memory.md:251) の完了条件も満たさないため、**今回の提案だけでは `done` 不可**とする。

不足情報: 対象の `run_0194`〜`run_0205` の HDF5・残差CSV・個別 VERDICT はローカルに存在しない。確認できたのは、base/new とも4/6本が幅を超えたという[保存集計:16](/home/sano/work/forge-sern-design/case/66.hostmem_regression/fixedwidth_c44dual_ckpt100/RESULT.txt:16)までで、原本の再集計ではない。実バイナリのデバイスコード同等性、読取り前初期化の証拠、完全な経路台帳も不足している。収束・準定常性の主張はしない。

ファイル変更・forge 実行なし。**plan 未反映**。呼び出し側で §6.3・§8 に反映する診断として返す。
