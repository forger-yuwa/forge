# codex 諮問 (diagnose): farfield-v0-deterministic-design

- **brief**: [`notes/reviews/briefs/2026-10-03-farfield-v0-deterministic-design.md`](../../notes/reviews/briefs/2026-10-03-farfield-v0-deterministic-design.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-10-03
- **commit**: `f18bac06` (feature/sern-design)
- **codex**: effort `xhigh`, 9.9 min, rc=0
- **結論**: ?

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | A：固定小数点の影アキュムレータ | **却下（今回の第一選択として）**。各寄与の量子化と最終丸めを導入し、既存の逐次 `float` 加算とは別の演算になる。小さい新旧差を量子化で消す可能性もある。現行は各面寄与を直接加算している。[SLAU:661](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:661)。対案は、加算命令を保持して実行順だけ固定する。 |
| **Major** | B：double 積算後に float 化 | **却下**。「丸め境界だけが問題」は成立しない。float32 で表現できる寄与 `{2⁶⁰, −2⁶⁰, 1}` を double で順次加算して float32 に戻す小試験でも、順序によって **0 / 1** になった。今回の流束でこの条件が発生した証拠ではないが、決定性の保証には使えない。対案は固定順実行。 |
| **Major** | D：差分カーネルと流れ残差だけ決定化 | **却下**。変更されていないカーネルでも、更新に入る非決定性は残せない。例えば化学種拡散は `res_roY`・輸送対角・`res_roe` に加算する（`a6ceee0b:solver_density_cuda/cuda_forge/speciesTransport_d.cu:271–285`）。対案は、**対象設定で実行され、更新に影響する共有書込み**をすべて対象にする。 |
| — | C：静的確認／E：新ソースから farfield を除去 | **C の静的確認は補助として採用、C/E を主試験にする案は却下**。farfield 不在時の未確保経路は確認できたが、それだけで生成コード・更新全体の同一性は証明できない。[未確保経路:58](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:58)。比較対象は、実際の旧・新ソースにそろえる。 |
| **Major** | 既存 `FORGE_DUMP_LEDGER` を全場比較にそのまま使う | **却下**。面バッファは **200,000 面で打切り**、節点量は **1 要素ずつ `cudaMemcpy`** している。[容量:665](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:665)、[転送:706](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:706)。対案は元の面・節点 ID に対応するバイナリ配列の一括／チャンク出力。欠落・容量超過は診断失敗とする。 |
| **Major** | 旧コミットを再ビルドすれば配布済み旧と同一 | **要再検証**。回収原本にはバイナリ識別子はあるが、ビルド条件との対応を証明する情報がない。[原本:1](/home/sano/work/forge-sern-design/notes/investigations/2026-10-03-farfield-evidence/c46_V0B_REPRO_20261003.txt:1)。対案は、診断パッチ適用前の再ビルドを配布物に照合し、証明できなければ結論を「再ビルドした両ソースの比較」に限定する。 |

**結論:** A〜E ではなく、**同じ節点へ書き込む面を別の色に分け、色を固定順に実行する診断方式 F**を両版へ適用し、元の `float atomicAdd` を保持した完全な 1 outer step を比較する。

**第 1 仮説:** 非決定的な面寄与の積算と、旧 3 回の観測範囲を許容限界にした判定規則が、V0 後半 FAIL の主因である。**確度: 中**

- **根拠:** `case/46.sern_design/run_1032_v0b_{o,n}{1,2,3}/` の 200 step 比較では、`roUx` の `S_old / S_new / D` が `0.2778 / 0.3838 / 0.5125`、`roOmega` が `4.459e5 / 7.511e5 / 6.332e5`。版内変動と版間差が同じ桁にある。ただし、同分布の証明ではない。[判定原本:19](/home/sano/work/forge-sern-design/notes/investigations/2026-10-03-farfield-evidence/c46_V0B_REPRO_20261003.txt:19)
- 流れだけでなく、スカラー残差と輸送対角にも面単位の競合加算がある。[scalarTransport:181](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/scalarTransport_d.cu:181)
- **反証条件:** 入力と加算順を固定し、旧同士・新同士がビット一致する条件で、新旧の面寄与・残差・対角・更新量に再現する差が残ること。

**第 2 仮説:** 共有カーネルの変更、またはコンパイラによる生成コードの違いが、決定的な差を生んでいる。**確度: 低、未確認。**

**第 3 仮説:** 実バイナリのソース対応、実効設定、診断環境変数に交絡がある。**確度: 低、未確認。** `FORGE_DIAG_FACE_VEL_CELL` の無効化は、今回も実行時の証拠が必要。

**判別 A/B:** 変える要因は旧版／新版だけ。`run_0971` の同じ入力を固定し、両版で **1 outer step・元の `nStepInner: 5` を各 2 回、別プロセスで実行**する。実装と事前判定を次の形にする。

1. **面の色分けは CPU で一度だけ決定する。**  
   元の面 ID 順に走査し、書込み先節点が同じ面には異なる色を割り当てる。同じ色の中では、同じ配列要素へ複数スレッドが書かないことを検査する。書込み先に ghost が含まれるカーネルでは、それも含める。色分け表を両版で共用し、ハッシュを照合する。

   各カーネルの色を同じ CUDA stream 上で順番に起動する。既存の `atomicAdd(float*, float)` は残す。これなら節点ごとの加算順が色順に固定され、各加算の丸めも保持できる。global memory の float atomic 加算には subnormal の扱いもあるため、通常の `+=` への置換は避ける。[NVIDIA PTX 仕様](https://docs.nvidia.com/cuda/archive/11.4.4/parallel-thread-execution/index.html#atom)

2. **改修は全 `atomicAdd` ではなく、実行経路の面 scatter に絞る。**  
   最初に対象にするのは次の経路。

   - 内部 SLAU、既存境界対流。
   - 内部粘性、壁粘性。
   - SST／化学種のスカラー移流・拡散と輸送対角。
   - 化学種 Fick 拡散とエンタルピー拡散。
   - 実効設定で使う追加の面加算経路。

   SLAU は既に `FaceGeom.nLoopPlanes / loop_planes` を持つため、色ごとのリストを渡せる。[SLAU:61](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:61)。スカラーも既存の面リストを利用できる。直接 `ip`／`ib` を作る粘性・境界カーネルには、元 ID への間接参照を追加する。境界の `bvar[ib]` は元の `ib` で読む。

   **wrapper 全体を色ごとに呼び直さない。** 残差ゼロ化などは元の位置で一度だけ行い、個々の kernel launch を色ループに置き換える。対流 wrapper は冒頭で残差をゼロ化している。[初期化:318](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:318)

3. **既に固定順の部分は維持し、比較対象に含める。**  
   ブリーフの改修規模の見積りは過大である。対象版では、node の `gradLSQ=2` は節点ごとに隣接面を順次走査し、スカラー LSQ も同様。[NS 勾配:865](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:865)、[スカラー勾配:935](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:935)。block-DPLUR と化学種 DPLUR も節点ごとの gather である。

   ただし、**実効設定が GG なら GG の面加算も色分けする**。診断を簡単にするために LSQ へ設定変更してはいけない。起動処理・BC の共有書込み・更新に使う reduction も監査し、未対応の競合が残れば前提ゲートを通さない。

   また、`nStepInner: 5` は固定残差に対する **5 回の Jacobi sweep** であり、残差を 5 回作り直す意味ではない（`a6ceee0b:solver_density_cuda/main.cpp:1579–1629`）。SST 更新、化学種更新・再正規化、最終ピン、outer 更新まで実行する。

4. **比較は「最初に違う場所」が分かる配列で行う。**  
   起動後の幾何・接続順・実効設定・状態、勾配・リミッタ、面寄与、段別残差、輸送対角・ソース Jacobian、`dt_local`、block の対角・右辺、各 sweep の `ΔQ`、更新後 9 保存量を採取する。

   面データは元 ID で配置し、atomic カウンタで得た記録順を比較しない。出力はバッファを使い回して一括／チャンク転送する。比較するのは有効要素のビット列であり、HDF5 ファイル全体、未使用領域、時刻、ポインタ値ではない。

   **事前閾値:** 欠落・重複・範囲外・有効数値の NaN/Inf は 0。旧同士・新同士の全比較配列は **不一致要素数 0**。新旧間も **0 bit 差**を一致条件とする。1 ULP でもあれば不一致として位置を追う。ただし、それだけで物理的に重大な回帰とは判定しない。

5. **両版へ同じ診断パッチを適用する。**  
   一方だけの決定化は不可。旧は `2fa3826c`、新は `b0240cd7…` を生成した実際のソースと未コミット差分を確定させる。「`a6ceee0b` 系」だけでは不足する。

   配布済み旧との対応は、まず**パッチなし旧**を元の CUDA/compiler・オプション・依存ライブラリで再ビルドして照合する。成果物ハッシュ一致ならバイナリ同一性を示せる。不一致なら実行コード・定数・依存関係まで調べる必要があり、同じコミットという理由だけでは担保できない。対応が取れなければ、今回の診断は再ビルド版どうしに限定する。

→ **A＝版内も新旧間も一致:** この入力・1 更新・固定順実行について、決定的な差を棄却する。次は同じ診断を元の 200 step 区間へ延ばし、後続状態で初めて差が出ないか確認する。**これでも通常実行の V0 合格には自動昇格しない。**

→ **B＝版内は一致、新旧間だけ不一致:** 再実行変動だけの説明を棄却する。最初に異なる面寄与／段／sweep を特定し、その箇所のソースと生成コードへ調査を限定する。

→ **版内で不一致:** 診断不成立。未固定の加算、共有書込み、未初期化、出力欠落を調べる。新旧差の解釈には進まない。

**V0 の扱い:** 現行の `VERDICT: FAIL` は保持する。1 更新の一致を V0 合格にする経路はない。V0 の再判定には、実バイナリとの対応を確立した上で、生産実行に対する受入試験が別途必要。基準を再設計する場合も、今回の FAIL を消さず、新しい検証計画として事前に確定する。[plan §5.1 #5a](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:167)

**やらない方がよいこと:** 全 478 か所の一括置換、固定小数点／double 化の一致を元演算の一致と扱うこと、`nStepInner` を減らすこと、乱流・化学種を凍結すること、blocksize の固定だけで決定性を主張すること。warp 内でも同一アドレスへの atomic 実行順は保証されない。[NVIDIA CUDA 仕様](https://docs.nvidia.com/cuda/cuda-programming-guide/03-advanced/advanced-kernel-programming.html)

**呼び出し側の前提への異議:** 「全更新の決定化」と「ソルバ全体の全 atomic 改修」は別である。対象経路には既に固定順の勾配・陰解法がある。一方、ソース差分のない面加算を非決定のまま残す D では不十分。また、固定順実行で一致しても、通常実行の多 step 分布が同じとは証明できない。

**不足情報:** 対象 6 run の実体、入力・物性 DB の完全なハッシュ、実効設定、実行時環境、両配布バイナリのビルド記録がローカルにない。したがって、色分けすべき最終カーネル一覧と所要時間は未確定。今回確認した実測根拠は回収済み判定原本である。

ファイル変更・forge 実行はしていない。提案は **plan 未反映**。反映先は `plans/active/boundary-node-farfield-characteristic.md` §5.1 #5a・§6。
