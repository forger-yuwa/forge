# 出力の処理が書き換える host の写しと、その読み手の網羅 (settle plan §5.1 #1d)

- **日付**: 2026-10-16
- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md) §5.1 #1c・#1d、§6 V1 の改訂 (f′)・(g′)
- **発端**: codex diagnose [`2026-10-16-settle-run0012-result-diagnose-settle-run0012-result.md`](../reviews/2026-10-16-settle-run0012-result-diagnose-settle-run0012-result.md) の Major 4 (#1c の「全配列の差 0」は不合格のまま残し、契約 = 次の数値更新が読む状態の不変の合否は、変わった host の配列ごとに読み手を網羅して決める)
- **コード**: `feature/card-local-mesh` 199ceb18 + 診断 `namespace oid` の追記 (下の §2、既定 off)。float のビルド `.build-native/card4a/forge` (sha256 2e4687c8…)、FP64 のビルド (同じソースで `flowFormat.hpp` の typedef を double、scratchpad、sha256 c2536a2a…)

## 1. 結論

- #1c で腕 B・C (`SNAPSHOT_REQUEST`) だけが変えた host の配列 (Euler 43・B0 NS 58・CHT 57) を、`gpu: 1` の時間進行中の読み手で分けると、**(a) 数値の経路で読まれない 153、(b) 読む前に必ず上書きされる 5 (CHT の連成の壁だけ)、(c) 上書き前に読む 0** (計 158)。
- 「上書き前に読む」読み手は 2 つある (farfield の流束が読む自由流の `bvar`、CHT の `flux: q_compact` が読む壁温 `Ts`)。どちらの配列も device 側に書き手がなく、host と device が常に同じ値を持つので、出力の D2H は同じビットを写すだけである。4 ケースの実測で、出力が D2H しても値は変わらず (B.nDiff 0)、S0 の時点で host と device が一致した (s0twin 0)。
- farfield のケース (case/58 V2b oblique、TP + SST、全 6 面 farfield) で #1c と同じ A/B を回した。device は全腕でビット同一で、host で変わったのは CV の写し 21 個だけだった。farfield の host の `bvar` 96 個は、境界の出力の有無 (`outputHDFflg` 0/1) にも step (最初の境界の出力の前 63・後 163) にもよらず不変だった。
- FP64 のビルドでも 4 ケースの結果は同じだった (CHT だけ slip 境界の 2 配列が増えるが、分類は同じ (a))。
- したがって §5.1 #1d の合格条件 (「上書き前に読む」が 0、または値の不変性が全部示せる) を満たす。ただし V1 (f′)・(g′) の合否は上位の判断として残す (§9)。

## 2. 方法

1. **診断の追記** (`solver_density_cuda/main.cpp` の `namespace oid`、既定 off、数値の経路は変えない):
   - 報告の表に `s0twin` 列を足した。S0 (腕を回す前) の時点で host の写しが対応の device の配列と違う要素の数で、出力が D2H しない配列でも数える。
   - `FORGE_DIAG_OUTPUT_INVARIANCE_POISON_EXTRA=<名前,...>` を足した (`CONTINUE=poison` のときだけ使える)。名前で指定した host の配列も NaN で埋める陽性対照である。
2. #1c と同じ 3 ケースを現行のバイナリで回し直した (配列の数は #1c と同じ 43・58・57)。
3. host の写し (`var.c`・`var.p`・`bc.bvar`・`bc.bint`・`bc.diagVar`) の**読み出し箇所**を全ソースで grep した。対象は `.c[`・`.c.at`・`.c.find`・`.c.count`・`var.c` の参照と反復、`.bvar` のすべての形、`diagVar`、`->c`。`variables` を参照・ポインタとして保持する型も grep し、診断の `TpuState` だけだった。読み出し箇所ごとに、(i) 時間進行中に到達するか、(ii) 直前に D2H・再計算で上書きされるかを読んだ。間接の経路 (H2D・`cudaMemcpyToSymbol`・カーネルの値渡し) も、どれも host の読み出しから始まるので、この列挙に含まれる。H2D の全 114 箇所は送り元も確かめた (§4)。
4. farfield の A/B (§7) と FP64 のビルドでの再実行 (§8)。

## 3. 出力の処理が書き換える host の状態

`writeStepOutputs` (`main.cpp:1902-1944`) が呼ぶもののうち、host の状態を変えるもの:

| 書き換える状態 | 箇所 | 内容 |
| --- | --- | --- |
| `var.c[名前]` | `output.cpp:120-123` | 全場に書く量 (`effectiveOutputNames`) と `Ht`・`k` を D2H |
| `var.c[履歴]` | `output.cpp:227` | dual-time の `/CHECKPOINT` の履歴 (`roN` ほか) を D2H (dual-time だけ) |
| `bc.bvar[全名前]` | `output.cpp:350` → `mesh.cpp:122-134` | `outputHDFflg: 1` の境界の bvar をすべて D2H |
| `bc.diagVar` | `output.cpp:354` → `conjugateWall.cpp:327-340` | 界面の診断 14 量を作り直す (`interfaceDiag: 1` の壁だけ、`:208-212`) |
| `bc.planes_local`・`inodes_l2g`・`inodes_g2l` | `output.cpp:363` → `mesh.cpp:136-171` | 最初の 1 回だけ作る。読むのは `output.cpp:384-463`・`:568` だけ |
| `conjugateWall` の `firstInterior` の static キャッシュ | `conjugateWall.cpp:182-199` | 幾何と `cfg` だけから決まる。連成の壁は起動時 (`:634`) に作る。誰が先に作っても中身は同じ |

`monitorOutputs` の節点の出力・流量 (`monitorOutputs.cpp:225`・`:538`) と `passiveFctHistoryToHost` (`speciesTransport_d.cu:2164-2178`) は、局所の配列に D2H するだけで host の写しを変えない。point probe は probe の step だけで、診断は K がその格子に乗るのを拒否する。`writeConjugateState` はファイルを書くだけで、固体の状態は比べて不変だった。

## 4. 時間進行中の読み手 (`gpu: 1`)

| host の状態 | 時間進行中の読み手と、読む前の上書き | 分類 |
| --- | --- | --- |
| `var.c` | **数値の経路の読み手なし**。読むのは出力自身 (`output.cpp:196-211`・`:230`・`:263-264`、どれも自分の D2H の後)、診断 `gatherImplicitDiagSnapshot` (`main.cpp:855` で D2H → `:858-864`。`implicit_diag_logger` が有効なときだけ `:2152`)、`FORGE_PIN_DIAG` (`main.cpp:2240`・`:2260` で D2H → 読む)。host 分岐は `gpu ≠ 1` のときだけ通る: `conjugateWall.cpp:86-89` (gpu 1 は `:82-84` で D2H)、`monitorOutputs.cpp:227` (gpu 1 は `:225`)、`main.cpp:3375` (gpu 1 は device の集約 `:3351-3372`)、CPU の解法 (`update.cpp:9`・`:76`・`:118`・`:153`・`:189`、`gradient.cpp:9`、`dependentVariables.cpp:9` は gpu 1 で即 return、`solveNavierStokes.cpp` の関数は呼び出し元なし)。起動時だけのもの: `setInitial.hpp:11-456`、`variables.cpp:1003-1188`、`calcWallDistance_kdtree.cpp:193-203`、`main.cpp:1156-1161`・`:4072` | (a) |
| `bc.bvar` (連成・farfield 以外) | 出力自身 (`output.cpp:430-447`、`:350` の後) と、出力から呼ぶ `fillInterfaceDiagnostics` (`conjugateWall.cpp:226-270`、`output.cpp:350` の後)。数値の経路の読み手なし | (a) |
| `bc.bvar` (連成の壁 `conjugate: 1`) | `updateConjugateWalls` の q_eff 経路は `conjugateWall.cpp:930` で D2H → `:931` の `fillInterfaceDiagnostics` が `ifaceFw`・`ifaceRraw`・`ifaceRro`・`qwall`・`Ts` を読む (`:226-270`) → `:965` で `Ts` を読む | (b) |
| `bc.bvar["Ts"]` (連成の壁、`flux: q_compact`) | `conjugateWall.cpp:965` (非 fem) と `:745` (fem2d) が D2H なしで読む → `:992`・`:850` で書き、すぐ `:1006`・`:883` で H2D | (c) → §6 |
| `bc.bvar` (farfield の `ro`・`Ux`・`Uy`・`Uz`・`Ps`・`k`・`omega`・`Y{s}`) | `farfieldFlux_d_wrapper` (`convectiveFlux_d.cu:119-125`) が要素 [0] を読み、`inf` として値渡しでカーネルへ (`:144`)。D2H なし | (c) → §6 |
| `bc.diagVar` | 出力自身 (`output.cpp:450-466`、`:354` で作り直した後) と、`updateConjugateWalls` (`conjugateWall.cpp:932-941`、`:931` で作り直した後。fem2d は `:944` → `:720-731`)。`fillInterfaceDiagnostics` が早期 return する構成では誰も作らないので、`:935-940` で停止する | (b) |
| `bc.bint`・`var.p` | 出力は書かない (3 ケースとも変化なし) | — |
| 間接 (H2D) | 時間進行中に host の値を device へ送るのは `conjugateWall.cpp:883`・`:1006` (`Ts`) だけ。他の H2D は、起動時 (`mesh.cpp:79`・`:114`、`variables.cpp:1076-1188`、`main.cpp:1161`・`:1763`、`boundaryCond.cpp:412`、`setInitial.hpp:454`、`calcWallDistance_kdtree.cpp:203`)、ポインタ配列・定数・幾何 (`boundaryCond_d.cu:916`・`:1233`、`ransTransport_d.cu:235`、`residualMonitor_d.cu:201`、`speciesTransport_d.cu:1784`、`condensationTransport_d.cu:1376`、`calcGradient_d.cu:853`)、診断 (`main.cpp:2350`・`:2533`、`:1356-1521`) のどれか。`copyVariables_cell_plane_H2D_all`・`copyVariables_probes_H2D_all` は呼び出し元がない。`cudaMemcpyToSymbol` は `cfg` か定数を送る (例 `convectiveFlux_d.cu:183-198`) | — |

## 5. 変わった配列ごとの分類

列の意味: 名前は報告の名前、n は B・C の差の要素数と `s0twin` (すべて一致)、twin `=dev` は「写した値 = S0 の device の値」。

### Euler (case/71 `run_0009_euler` の写し、K 163): 43 個、すべて (a)

| 群 | 配列 | 分類 |
| --- | --- | --- |
| CV の写し (16) | Ht・P・T・Ux・Uy・Y0・Y1・k・omega・ro・roUx・roUy・roY0・roY1・roe・sonic | (a) |
| b1 `inlet_Pressure` (6) | Ps・Ts・Ux・ro・roUx・roe | (a) |
| b2 `outlet_statPress` (7) | Ts・Ux・Uy・ro・roUx・roUy・roe | (a) |
| b3 `slip` (8) | Ps・Ts・Ux・Uy・ro・roUx・roUy・roe | (a) |
| b4 `axis` (6) | Ps・Ts・Ux・ro・roUx・roe | (a) |

### B0 NS (run_0010 の nozzle.h5、K 163): 58 個、すべて (a)

| 群 | 配列 | 分類 |
| --- | --- | --- |
| CV の写し (20) | Ht・P・T・Ux・Uy・Y0・Y1・k・omega・ro・roK・roOmega・roUx・roUy・roY0・roY1・roe・sonic・vis_lam・vis_turb | (a) |
| b1 `inlet_Pressure` (6) | Ps・Ts・Ux・ro・roUx・roe | (a) |
| b2 `outlet_statPress` (9) | Ts・Ux・Uy・kb・omegab・ro・roUx・roUy・roe | (a) |
| b3 `wall_isothermal` (`interfaceDiag: 1`、連成なし) (15) | Ps・ifaceFw・ifaceRconv・ifaceRpre・ifaceRraw・ifaceRro・omegab・qwall・ro・roe・twall_x・twall_y・utau・ypls・diagVar | (a) (`updateConjugateWalls` は `conjugateWall.cpp:901` で return。読むのは出力だけ) |
| b4 `axis` (8) | Ps・Ts・Ux・kb・omegab・ro・roUx・roe | (a) |

### CHT (case/65 C1 の写し、fem2d・q_eff、K 263): 57 個、(a) 52・(b) 5

| 群 | 配列 | 分類 |
| --- | --- | --- |
| CV の写し (10) | Ht・P・T・Ux・Uy・ro・roUx・roUy・roe・sonic | (a) |
| b1 `inlet_uniformVelocity` (5) | Ps・Ts・ro・roUx・roe | (a) |
| b2 `outlet_statPress` (7) | Ts・Ux・Uy・ro・roUx・roUy・roe | (a) |
| b3・b4・b6 `slip` (8・8・6) | Ps・Ts・Ux・Uy・ro・roUx・roUy・roe (b6 は Uy・roUy なし) | (a) |
| b5 `wall_isothermal` 連成 (8) | Ps・ifaceRconv・ifaceRpre・ro・roe・twall_x・twall_y・utau | (a) (`conjugateWall.cpp:930` の D2H で上書きもされる) |
| b5 `wall_isothermal` 連成 (4) | ifaceFw・ifaceRraw・ifaceRro・qwall | (b) `:930` で D2H → `:226-270` で読む |
| b5 `wall_isothermal` 連成 (1) | diagVar | (b) `:931` で作り直す → `:932-941` で読む |

連成の更新は iStep 250・300 (`warmup` 100・`interval` 50)。K 263 の後に数値の経路が読むのは iStep 300 の `:930` の D2H の後である。

## 6. 上書き前に読む (c) 読み手と、値の不変性

どちらの読み手の配列も、4 ケースの「変わった配列」に入っていない。理由は、device 側に書き手がなく host と device が同じ値を持つことである。

**farfield の自由流** (`ro`・`Ux`・`Uy`・`Uz`・`Ps`・`k`・`omega`・`Y{s}` は type 1 = config の一様値):
- host の値は起動時に floats から入れ (`mesh.cpp:64-69`)、同じ値を H2D する (`mesh.cpp:76-80`)。
- device の書き手はない。`applyBconds` は farfield で何もしない (`boundaryCond.cpp:498`)。既存の境界流束は通らない (`convectiveFlux_d.cu:519-522`)。farfield の流束は自由流を値渡しで受け、残差・`massflux`・面の値配列だけを書く (`:144`)。RANS の境界は return する (`ransBoundary_d.cu:256-258`)。化学種の境界が `bvar_d` を読むのは入口だけ (`speciesTransport_d.cu:909-913`)。境界の勾配は `bvar` を参照しない (`calcGradient_d.cu:1265`)。
- type 0 のキー (`roUx`・`roUy`・`roUz`・`roe`・`Tt`・`Pt`・`Ts`) は device が未初期化のまま確保される。ただし `hv` は読まない。実測では host・device とも 0 で一致した。
- 実測 (§7): farfield の host の bvar 96 個 (6 境界 × 16) は、float・FP64、`outputHDFflg` 0/1、K 63/163 のどれでも B.nDiff 0・s0twin 0。

**連成の壁温 `Ts`** (`flux: q_compact` が D2H なしで読む):
- 連成の壁は `wall_isothermal` に限る (`conjugateWall.cpp:629-633`)。
- `wall_isothermal` の device のカーネルは `Tsb` を読むだけである (`boundaryCond_d.cu:308-447`、WMLES `wmlesWallModel_d.cu:175`)。`Tsb` に書くのは他の種別だけ: `boundaryCond_d.cu:102` (slip)・`:245` (wall)・`:684-689` (outlet)・`:898`・`:1215`・`:1462` (入口)・`:1618-1621` (outflow)・`:1756` (periodic)、`nodeWallDirichlet_d.cu:231` (`:249` で断熱の `wall` だけ)、`fluct_variables_d.cu:286`。
- host が書いたら直後に H2D する (`conjugateWall.cpp:850`→`:883`、`:992`→`:1006`、起動時の `boundaryCond.cpp:412`・`mesh.cpp:79`)。
- 実測: CHT の float・FP64 とも `host b5:Ts` は B.nDiff 0・s0twin 0。q_eff の run での確認で、q_compact の run は回していない (読み手は同じ配列で、書き手の議論は経路によらない)。

したがって、出力の D2H が写す値 = その時点の device の値 = (c) の読み手が読む値。

## 7. farfield の A/B (case/58 V2b oblique の写し)

`case/58.farfield_verification/setup_v2b.py --dir oblique --ek 0 --steps 600` を scratchpad に作った。条件は 24×12×12 の箱、全 6 面 farfield、TP (EXH/AIR) + SST、3 面から斜めに流入、内部は T 600 K・Y_EXH 0.13・k・ω 10 倍で始める。species は手元の SERN run の写しを使い、IC に species の印が無いので `FORGE_ALLOW_UNVERIFIED_SPECIES=1` で回した。`outStepInterval` は 100、`output: {boundaryMassflux: 1}` を足して腕 C も回した。

| run | ビルド | outputHDFflg | K | device (345 配列) | host で変わった配列 | farfield の host の bvar 96 個 |
| --- | --- | --- | --- | --- | --- | --- |
| `ff_h1_k063` | float | 1 | 63 (最初の境界の出力より前) | 全腕で同一 | 21 (CV の写しだけ、全部 `=dev`) | B.nDiff 0・s0twin 0 |
| `ff_h1_k163` | float | 1 | 163 | 全腕で同一 | 21 | 同上 |
| `ff_h0_k163` | float | 0 (case/58 の元の設定) | 163 | 全腕で同一 | 21 | 同上 |
| `ff64_k063` | FP64 | 1 | 63 | 全腕で同一 | 21 | 同上 |

続けた試験 (K 163 で S0 に戻してから 600 step):
- **毒入れ** `ff_poi1..5` (変わった 21 配列を NaN で埋める): 非有限なし。farfield の置換の計数は 0。`check_field_regress.py --repeat ff_ref1..5 --candidate ff_poi1..5 --step 600 --factor 2` は PASS (比の最大 1.10)。
- **陽性対照** `ff_posctl` (さらに `host b1:ro` を NaN で埋める): farfield のカーネルは非有限の外側状態を内部の状態で置き換える (`farfieldFlux_d.inc.cuh:111`)。そのため NaN は場に入らず、置換の計数が 1 回あたり 169 面 (b1 の全面) で増えた。計数の合計は 6253 = 37 × 169 で、K 163 から次の境界の出力 (200) までの 37 step にあたる。step 200 の出力の D2H が `host b1:ro` を device の値 (自由流) で上書きして毒が消えた。これは §6 の値の不変性と同じ仕組みである。
- 600 step の場は ref に対して PASS (ノイズ床の内) で、**場の比較ではこの読み手の毒を検出できなかった**。#1c の「毒を入れても NaN が出ず、差はノイズ床の内」は補助の証拠であり、網羅の代わりにはならない (§9)。

## 8. FP64 のビルドへの適用

- **コード**: §4 の分類は「どの関数が host の写しを読むか・その前に上書きがあるか」で決まり、`flow_float` の型によらない。出力・CHT・farfield・節点の出力の経路に typedef で分かれる分岐は無い。`sizeof(flow_float)` で分かれるのは診断 (`quantizeStateDiag.cpp:157`、`geomAbDiag.cpp:238`・`:243`) と報告だけである。D2H・H2D の wrapper も host の `std::vector<flow_float>` も同じ関数で、変わるのは要素の大きさだけ。§6 の「device に書き手が無い」も型によらない。
- **実測** (FP64、`FORGE_CUDA_BLOCKSIZE=128` が要る。既定 512 は `convectiveFlux_d.cu:494` で `too many resources requested for launch`):

| ケース | float で変わった配列 | FP64 で変わった配列 | 違い | device |
| --- | --- | --- | --- | --- |
| Euler `e64_k163` | 43 | 43 | 名前も同じ | 全腕で同一 (307) |
| B0 NS `n64_k163` | 58 | 58 | 名前も同じ | 全腕で同一 (329) |
| CHT `c64_k263` | 57 | 59 | `b6:Uy`・`b6:roUy` (slip。float では S0 で host = device) が増えた。分類は (a) | 全腕で同一 (311) |
| farfield `ff64_k063` | 21 | 21 | 名前も同じ。farfield の 96 個は不変 | 全腕で同一 (345) |

全ケースで、変わった配列は `s0twin` = B.nDiff = C.nDiff、twin `=dev` だった (`diagVar` は twin なし)。

## 9. 残った懸念・範囲外

1. **将来の読み手**: 網羅は現行のコードについてだけ言える。時間進行中に `var.c`・`bvar`・`diagVar` を host で読むコードを足すと前提が崩れる。機械的な見張り (lint・試験) は作っていない。
2. **dual-time と陽解法**: 出力の dual-time の CHECKPOINT (`output.cpp:217-260`) は host の `var.c` の履歴を書き換えるが、読み手は無い (§4)。A/B は定常の陰解法だけで、dual-time・RK の経路では回していない。`writeStepOutputs` と読み手は同じである。
3. **CHT の q_compact**: §6 はコードと q_eff の run の実測で示した。q_compact の run は回していない。
4. **診断が比べないもの**: 数値のモジュールのファイル内の static な device の作業配列 (出力の関数からは届かない)、`cfg` (出力の関数は const で受ける。point probe だけは probe の step 以外は何もしない)。
5. **毒入れの検出力**: §7 の陽性対照のとおり、ガードのある読み手 (farfield の置換) では、毒を入れても場の比較では見えない。#1c の毒入れの結果は網羅の代わりにならない。
6. farfield の run は species の印の無い IC を `FORGE_ALLOW_UNVERIFIED_SPECIES=1` で通した (診断専用)。

## 10. run の所在 (一時、scratchpad)

`/tmp/claude-1000/-home-sano-work-forge/e254baf2-8857-4243-99d2-53804f553ad9/scratchpad/hr/runs/` (セッション終了で消える。case/ の下には作っていない):
`e_k163`・`n_k163`・`c_k263` (float、#1c の再実行)、`ff_h1_k063`・`ff_h1_k163`・`ff_h0_k163` (float、farfield)、`ff_ref1..5`・`ff_poi1..5`・`ff_posctl` (float、続けた試験)、`e64_k163`・`n64_k163`・`c64_k263`・`ff64_k063` (FP64)。各 run の `output_invariance_<K>.txt` が報告の原本で、`changed.tsv` が変わった配列の抜き出しである。FP64 のビルドは `.../scratchpad/hr/fp64/build/forge`。
