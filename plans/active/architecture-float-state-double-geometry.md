# 状態を float にして、差を取る幾何の量を double の座標から作る

## メタ

- **area**: `architecture`
- **status**: `draft`
- **related_docs**:
  - `methods/architecture/overview.md` (§6.2a 幾何の量の精度)
  - `methods/time_integration/implementation.md` (commit の丸め、`qAccumulatorFP64`)
- **related_plans**: [`time_integration-line-implicit-speed`](time_integration-line-implicit-speed.md) §5.1 #23・§6.21・§6.22 (発端)、[`tooling-nozzle-isothermal-wall-chain`](tooling-nozzle-isothermal-wall-chain.md) (冷却壁の格子と FP64 のビルド)
- **created**: `2026-10-10`
- **owner**: Claude (Opus 5.5)

## 1. 目的

case/45 の冷却壁の M6 ノズル (第一層厚 / 半径 ≈ 1e-7〜1e-6) は、float の座標では第一層が数 ulp にしかならないので、FP64 のビルドで回している。
float のビルドは同じ構成で 1 step が FP64 の 0.64 倍 (20.5 対 32.0 ms、line-implicit-speed §6.21) だが、残差は約 10 倍で使えない。
座標の差を取る量を、読み込み時に double の座標から作って float で渡す。これで状態は float のまま、θ_r・Q_w の到達を FP64 と同等に保ち、1 step を 3 割ほど短くする。

## 2. スコープ

- **やる**:
  - メッシュの座標・面重心を double でも読み込み (ホスト)、面ごとの差の量 (節点間ベクトル、拡散の幾何係数、再構成のベクトル) を double で作って `flow_float` の配列に入れる。
  - 毎 step 絶対座標の差を取っているカーネルを、それらの配列を読むように書き換える (調査 [`notes/investigations/2026-10-10-geometry-precision-inventory.md`](../../notes/investigations/2026-10-10-geometry-precision-inventory.md) の §1)。
  - 読み込み時に作る量 (LSQ の係数、軸対称の closure、node の壁関数の代表距離、弱形式の等温壁の d1) を double の座標から作る。
  - case/45 で float の到達を FP64 と比べ、標準ケースの回帰を確かめる。
- **やらない**:
  - ~~変換器の float 経路の修正。冷却壁の格子は従来どおり FP64 の変換器 (`stod`・17 桁) で作る。~~ **2026-10-10 にスコープへ移した** (ユーザ決定、§4.6)。
  - `geom_float` と `flow_float` を別の型にする (調査 §0-1 のとおり、配列の表とカーネルの引数の作り直しになる)。
  - 状態の commit の精度 (`qAccumulatorFP64` の軸対称対応)。§6 の測定で丸めによる停滞が見えたときに別に立てる (§5.1 #11)。
  - cell 離散化の正しさの確認 (使っていない)。ただし cell の経路もビルドと既存の回帰は壊さない。

## 3. 関連 docs と前提

- 現状の精度の扱い: `methods/architecture/overview.md` §6.2a。
- 調査: `notes/investigations/2026-10-10-geometry-precision-inventory.md`。
- 冷却壁の格子と FP64 の経緯: メモ [cooled-wall-mesh-precision]、`plans/active/tooling-nozzle-isothermal-wall-chain.md`。
- 発端の測定: `plans/active/time_integration-line-implicit-speed.md` §6.21 (float 20.5 対 FP64 32.0 ms/step、float の残差は約 10 倍) と §6.22 (codex の諮問。「混成 typedef では試験にならない、差を double で作る案を先に」)。
- B0 の水準での 1 step の変化と float の刻み (2026-10-10、`case/45.isobutane_m6_d155/run_0387_dqulp_B0`、FP64 で B0 の最終状態から 10 step、9 → 10 step の変化):
  - |Δq| が float32 の刻みの半分に満たない節点は、ρ 17.7 %、ρu 8.9 %、ρv 3.5 %、ρE 11.9 %、ρk 3.4 %、ρω 5.9 %。
  - 比の中央値は 2〜31 (刻みより大きく揺れている)。

## 4. 設計方針 (2026-10-10、codex diagnose `notes/reviews/2026-10-10-float-geometry-design-diagnose.md` の採否を反映)

### 4.1 どこを直すか

調査 §5 の順。**毎 step カーネルが絶対座標の差を作っている箇所**が本体で、読み込み時の量・接続の生成はその次。
ただし「座標の差の丸めが壁際の誤差の過半を説明する」こと自体はまだ確かめていない (codex の第 1 仮説、確度 中)。全部を書き換える前に、§6 の **V0 (状態を固定した演算の A/B)** で判別する。

### 4.2 差の量の作り方

1. **double の入力を持つ** (`mesh.cpp`): 次を `geom_float` の配列とは別に double のホスト配列にも読む。HDF5 が double ならそのまま。
   - `/MESH/COORD`、`/PLANES/centCoords`、`/CELLS/centCoords`
   - **面ベクトル `surfVect` と面積 `surfArea`** (codex Major 2: 面ベクトルも丸めて読み、さらに float の半径を掛けているので、幾何の精度と離散の閉包を分けて検査する)
   - node 離散化の中心 (= 節点座標) も double で持つ。
2. **使う側ごとに、今の式を保った係数を作る** (codex Major 1: 単一の `ediff = |S|²/|e·S|` は却下)。今の式は次のように違うので、精度の変更と式の統一を混ぜない。
   - 粘性: 絶対値 + `1e-30` のガード (`viscousFlux_d.cu:151`)
   - スカラー拡散: 符号付き + `1e-6·|e|·ss` の相対ガード (`scalarTransport_d.cu:110`)
   - 受動種の FCT: 符号付き + 絶対ガード (`passiveFct_d.cuh:32`)
   - 面ごとに e = cc1 − cc0 (3 成分、double で差を取って float に入れる) を作り、ガードと絶対値の扱いは利用側の式のまま e から計算する (`ss²` を `S·S` に替えることも、同一性の確認の対象にする)。
   - 再構成: 辺中点の ±0.5e は `g_reconEdgeMid == 1 && ip < nNormalPlanes` の**内部辺に限る** (codex Minor)。境界面と辺中点でない経路は、pc − cc (double で作って float) と境界の `fx` を残す。
   - node の境界の半割面 (e ≈ 0) では、係数を一律に評価しない。block DPLUR はその粘性の対角を明示的に除いている (`timeIntegration_d.cu:946`)。退化面は「有限値・使わない」を検査する。
3. **読み込み時の量と接続**:
   - LSQ の係数 (`gradLSQ` 1/2)、軸対称の closure (double の面ベクトルと半径から作るが、流束に渡す float の面ベクトルの和との整合も検査する)、弱形式の等温壁の d1/d2、`delta_les` を double の入力から作る。
   - 壁関数の代表点と距離は、**境界面ごとの `(irep, y)`** として作る (codex Major 3: 角では同じ節点でも面によって代表点が違う)。WMLES も同じ。`conjugateWall` の合算した法線は流用しない。
   - **周期の相手の対応付け (`mesh.cpp:568`)、LSQ の継ぎ目の同値類、ラインの接続 (`mesh.cpp:1067`) も double の入力から作る** (codex Major 4: 接続が FP64 と違うと同じ比較にならない)。
4. **書き換えるカーネル**: 調査 §1c〜§1j。読むのは e・r と、利用側の式で作る係数。座標の差の式は消す。

### 4.2a 精度の取り決めと、使う側ごとの表 (codex plan 段 M2 の反映)

- **取り決め**: 置き換えるのは `cc1 − cc0` (と `pc − cc`) の引き算だけ。double の座標から引いて 1 回だけ `flow_float` に丸めた値 e32 = fl(cc1₆₄ − cc0₆₄) を面ごとに持ち、**その先の式・型・ガード・対象面・向きは今のまま**にする。
- FP64 のビルドでは、e は今と同じ double の引き算になるので、残差はビット単位で一致する (V2 で確かめる)。例外は、block DPLUR の `ST` に落としてから差を取る経路 (`implicitSolvePrecision 0` で ST = float)。新しい経路は fl(e64) を ST にするので値が変わる。左辺だけの変化で、不動点は変わらない。

| 使う側 | 場所 | 今の差 | 差の後の型 | ガード | 対象の面 |
| --- | --- | --- | --- | --- | --- |
| 粘性 | `viscousFlux_d.cu:113-124` | flow_float の cc の差 | flow_float | \|D\| + 1e-30 | 内部面 (node は壁節点と内部節点の辺を含む) |
| スカラー拡散 (k/ω など) | `scalarTransport_d.cu:105-111, 202-208` | flow_float | flow_float | 符号付き、1e-6·\|e\|·ss の相対 | node は境界の半割面を除く |
| 化学種の拡散 | `speciesTransport_d.cu:283-289` ほか | flow_float | flow_float | 実装時に表へ追記 | 同上 |
| 受動種の拡散・FCT | `passiveKernels_d.cuh:184-188`、`passiveFct_d.cuh:28-32` | flow_float | flow_float | 符号付き、1e-30 の絶対 | 同上 |
| block DPLUR の粘性の対角 | `timeIntegration_d.cu:952-963` | **ST にしてから差** | ST | \|e·S\| + 1e-30 | 内部面 (node の境界の半割面を除く、`:946`) |
| スカラー DPLUR・前処理版 | `timeIntegration_d.cu:587-602, 1327-1335` | flow_float | flow_float | 実装時に表へ追記 | 実装時に表へ追記 |
| 壁関数の代表点 | `ransWallFunction_d.cu:146-166`、`wmlesWallModel_d.cu:158-175` | flow_float の cc_I − cc_W | flow_float | 候補なし・同率の候補の扱いを今のまま | 境界面 ib ごと |
| 再構成 | `convection/*.inc.cuh` | flow_float の cc・pc の差 | flow_float | — | 辺中点は内部辺だけ、ほかは pc − cc |
| リミタ | `limiter_d.cu`・`limiterPeriodic_d.cuh`・`passiveLimiter_d.cuh` | flow_float | flow_float | — | 同上 |

- 段 ① では本番の配列を変えない。新旧の e と、それを使った係数・分岐の結果を**別の配列に記録して比べる**だけにする。
- **接続** (周期の相手・LSQ の同値類・ラインの接続) と壁の `(irep, y)` (候補なし・同率の候補を含む) の一致の確認は、段 ② の合格条件にする (段 ① では作らない)。

### 4.2b 軸対称の closure (codex plan 段 M3 の反映)

- 今の closure は、半径の重みを掛けた後の実際の `sx/sy` を足し上げている (`variables.cpp:637`)。これは自由流の保持のための設計である ([`axisymmetric-freestream-hoop-gauge`](axisymmetric-freestream-hoop-gauge.md))。
- **ソースに使う closure は、デバイスに渡す最終の面ベクトル (半径の重みを掛けて float に丸めた後の値) を double で足し上げ、最後に丸める**方式にする。生の double の幾何から作る closure は、診断用として別に持つ。
- 段 ② の必須の確認: `P ≠ pRef` の一様な静止場で、正規化した残差と偽の加速度が今より悪くならないこと (許容値は段 ② の前に事前登録する)。

### 4.2c V0 の評価の経路 (codex plan 段 M1 の反映、段 ① に含める)

通常の実行を 3 回回しても、面ごとの比較にはならない。粘性とスカラー拡散は面の流束を残差へ直接足し込み、`F1` はその直前に更新される (`main.cpp:1967`)。そこで、専用の評価の経路を 2 つの環境変数で作る。

1. **`FORGE_DIAG_GEOMAB_DUMP=<file>` (float のビルド)**
   - 1 step 目の組み立てで、粘性とスカラー拡散を呼ぶ直前まで通常どおり進める。
   - その時点の入力を保存する。対象は、原始量・勾配・物性 (`vis_lam`・`vis_turb`・熱伝導率・cp)・`F1`・`axisym_divU`・壁モデルの入力・`fx`・半径の重みを掛けた後の `sx..ss`・接続・float の座標・e32。
   - 同じ入力から、粘性とスカラー拡散 (k・ω) のカーネルを 2 回呼ぶ。旧腕は今の座標の差、新腕は e32 を使う。
   - **atomicAdd の前の面の流束** (運動量の応力 3 成分・エネルギー (熱と仕事)・k と ω の拡散) を、面の番号つきで別のバッファに書く。残差も別のバッファに貯め、本番の残差は変えない。
   - 書き出して終わる。境界条件・EOS・勾配・乱流モデル・commit はやり直さない。
2. **`FORGE_DIAG_GEOMAB_REF=<file>` (FP64 のビルド)**
   - 同じ状態から 1 step 目の組み立てを同じ位置まで進める。
   - 保存した入力 (float を double に広げたもの) で配列を上書きする。面ベクトルも `double(S32/ss32)` で上書きし、FP64 の面ベクトルとは入れ替えない。
   - そのうえで、double の座標の差 (FP64 のビルド自身のもの) で同じカーネルを呼び、面の流束を書き出して終わる。
3. **先に確かめること**: 診断の旧腕の面の流束を足し上げたものが、変更していないカーネルの残差への寄与とビット単位で一致すること。
4. **比較** (Python):
   - 面ごとの誤差 (旧腕 − 参照、新腕 − 参照) を量ごとにまとめる。量は、運動量 (3 成分の和)、エネルギー、k、ω。
   - ノルム: マスクの中の面の L2 を、参照の流束の L2 で割る。参照の L2 が 0 のマスクは比較から外す。
   - マスク (壁からの節点の番号で固定。case/45 は j = 120 が壁):
     - 第一層: 端点の片方が j ≥ 119 の辺
     - 壁から 3 層: j ≥ 117
     - 全域
   - 再実行の差は 3 回の実行の差の最大 (面の流束は原子加算の前なので、0 と見込む)。
   - 結果 A・B・判別不能の規則は §6 V0 のとおり。判別不能なら段 ② に進まず、上位に諮る。

### 4.3 既定の扱い

- **判断 (codex、採用)**: 案 A (切り替えなし、常に新しい経路) を最終の設計にする。新しい常設の opt-in は作らない。比較には固定した旧コミットのバイナリを使う。既存の opt-in は、この変更を理由に消さない。
- FP64 のビルドでは、ガードと分岐を変えない限り、丸めの順序の範囲で同一になるはず。分岐やガードを変えた場合は、丸めの差として許容しない (§6 V2)。

### 4.4 メモリと速さ

- 面ごとに float 3 (e) + 必要な経路だけ 6 (r0・r1) + 壁の境界面ごとに (irep, y)。2D の node で +25〜50 B/節点。
  - **訂正 (2026-10-10、段 ④ の実装の後)**: r0・r1 は全面に持ち、デバイスとホストの両方にある。e と合わせて 2D の node で約 72 B/節点、case/45 (57 万節点・115 万面) の float で、ホストとデバイスにそれぞれ約 41 MB。今は受け入れ、V5 でピークの使用量を測る (codex diagnose 2026-10-10)。
- 速さの目標: float の新しい経路で、FP64 の 1 step の 0.70 倍以下。採用の判断には、単価だけでなく到達までの総時間も含める (codex V5)。

### 4.5 状態の commit の精度

- 2026-10-10 の測定 (FP64 の 9 → 10 step の変化と float の刻みの比) は、float で要求された更新と実際に反映された更新の測定ではない (codex Major: H2 は要再検証)。float の run で `dq_requested` と `Q_after − Q_before` を別々に記録する (§6 V4)。
- `qAccumulatorFP64` の軸対称対応は別 plan とする (codex 採用)。ただし、V4 が不合格だったときに、すぐ commit のせいにはしない。軸対称の拒否だけを外すこともしない (拒否には実装上の理由がある、`main.cpp:3371`)。

### 4.6 格子の幾何の精度を変換器で保証する (2026-10-10、ユーザ決定「これはやるんでしょうね」)

- **問題**: float のビルドの変換器は、座標を `stof` で読み、双対の面・体積・重心・壁距離を `geom_float` (float32) で計算して HDF5 に float32 で書く (`gmshReader.hpp:506-511, 2343`)。その格子では double の写し (cc64 など) も float の値そのものになり、段 ①〜④ の利点が消える (§5.1 #13)。
- **方針** (両方やる):
  1. **変換器は、ビルドの型 (`flow_float`・`geom_float`) に依らず、幾何を double で計算して double で書く**。対象は、座標の読み取り、primal の幾何、双対の面ベクトル・面積・面重心、体積・CV の重心、境界の半割面、壁距離。
     - **やり方は案 (b)** (codex diagnose 2026-10-10 で採用、`notes/reviews/2026-10-10-converter-double-geometry-diagnose.md`)。変換器は `variables`・`cuda_forge` を link して共用の型で H2D をするので、案 (a) (変換器だけ double の型でコンパイル) は共用の型と関数の整合を保証できない。
     - **型の置き換えだけでは足りない**: `gmshReader` は `mesh` を継承し、座標を共用の `node` に入れ、双対を作った後も共用の `cell`・`plane` に戻し、出力はそれらを読む (`gmshReader.hpp:34, 516, 2224, 2240, 2343, 2423, 2465`)。そこで、座標・primal の幾何・双対の幾何・境界の半割面の**正本を変換器専用の double の配列・構造体に置き**、出力はそこから書く。共用の `mesh` へは、初期値の設定に要る分だけ一方向に丸めて渡し、丸めた値から幾何を計算し直して書くことはしない。
     - **壁距離**: 壁の点と評価点を `msh.nodes/cells/planes` から取り、`var.c["wall_dist"]` に入れて、それを書いている (`calcWallDistance_kdtree.cpp:137, 159, 175`、`gmshReader.hpp:2525`)。double の点群を受けて double の配列を返す計算の部分を切り出し、`/VALUE/wall_dist` はその配列から直接書く。ソルバ用の写しと HDF5 用の正本を分ける。
     - 共用の型と `cuda_forge` のビルドはそのまま。変換器から CUDA の依存を外すのは別の変更とする。
  2. **ソルバは起動時に、使う幾何のデータセット全体の型を確かめる** (座標だけでなく `/PLANES/surfVect`・`surfArea`・`centCoords`・`/CELLS/centCoords`・`volume` など、`mesh.cpp:238, 261, 265, 313`)。float32 が 1 つでもあれば、まとめて 1 行の警告「入力の時点で失われた幾何の精度は回復できない」を出す。止めはしない (既存の float32 の格子でも回せるように)。後から double に変換した格子は精度を回復しないので、型の確認を生成の精度の証明とは扱わない。
  3. **手順** (`procedures/calculation-workflow.md` の格子の変換の節): 格子は double の幾何で書く (1 の後はどのビルドの変換器でもよい) と明記し、run の準備で座標の型を確かめる。
- FP64 のビルドの変換器の出力は今と同じ値になること (型も値も変えない)。

## 5. 実装ステップ

codex の勧めどおり一度に変えず、4 段に分けて各段で V2 を通す。

1. **段 ①**: double の入力と係数の生成 (`mesh.cpp`・`mesh.hpp`・`variables.cpp`・`variables.hpp`)。この段では利用側を切り替えず、旧の係数と並べて比べられるようにする (V0・V1 の準備)。
2. **段 ②**: 前処理 (LSQ・周期の相手・壁の代表点・ラインの接続・closure・d1・`delta_les`)。
3. **段 ③**: 粘性・拡散・陰解法の対角 (`viscousFlux_d.cu`・`scalarTransport_d.cu`・`speciesTransport_d.cu`・`passiveKernels_d.cuh`・`passiveFct_d.cuh`・`timeIntegration_d.cu`・壁関数)。
4. **段 ④**: 再構成・リミタ (`convection/*.inc.cuh`・`limiter_d.cu`・`limiterPeriodic_d.cuh`・`passiveLimiter_d.cuh`)。
5. 検証 (§6) と文書の更新 (`methods/architecture/overview.md` §6.2a、`procedures/` の FP64 の手順)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 設計の諮問 | **判断: 2026-10-10 codex diagnose** — 案 A を目標に、まず係数の符号・ガード・境界での定義を直し、状態を固定した粘性・拡散の A/B を 1 回行って、座標の差の修復が誤差を減らすか判別する。指摘 (Major 8・Minor 1) は全件採用 (§4・§6 に反映) | F |
| 2 | codex plan 段 | 反映後の §4・§6 を `codex_review.py --stage plan` | F |
| 3 | 段 ① と V0・V1 | 段 ① を実装する (§4.2a の取り決め、本番の配列は変えない、§4.2c の評価の経路)。合格: float と FP64 のビルドが通る。FP64 で e が今の差とビット一致。診断の旧腕の和が変更していないカーネルの寄与とビット一致。V1 の係数の検査 (接続は段 ②)。そのうえで V0 を回し、結果 A なら段 ②〜④ へ、B か判別不能なら上位に諮る **完了 (2026-10-10、§6.2)**: 段 ① を実装し、既定の経路は再実行の差の範囲で同じ。V1 は e32 の誤差 ≤ 6e-8 (旧は最大 30 %)。V0 は結果 A (新旧の誤差の比 ≤ 1e-5) → 段 ②〜④ へ | O → F |
| 4 | 段 ② | §5 の 2。合格: §6.3 の 1〜5 **完了 (2026-10-10、§6.5、改めた基準ですべて PASS)**: 事前登録の判定は case/48 の cInt だけ FAIL (物差しの欠陥)。**判断: 2026-10-10 codex diagnose** (`notes/reviews/2026-10-10-float-geometry-stage2-result-diagnose.md`) — 無条件の合格は保留し、旧基準の FAIL を残して改めた基準 (§6.5 の 4a) で判定し直す。判定の台本の欠陥 (欠損を一致と数える等) を直して集計し直し、closure を使う経路の FP64 の新旧比較を足す (§6.5 の 4b・4c)。これらが通れば段 ③ へ | O → F |
| 5 | 段 ③ | 合格: §6.4 (V2 の実施形) と §6.5 の 4d の追加 **完了 (2026-10-10、§6.7、すべて PASS)**: 本番の経路が e32 を読む (粘性・スカラー・化学種・受動種・3 種の DPLUR の対角)。FP64 で値が変わるのは block DPLUR の isp 0 の左辺だけ (implementer の報告) | O |
| 6 | 段 ④ | 合格: §6.6 (事前登録 2026-10-10) **完了 (2026-10-10、§6.9、すべて PASS)** | O |
| 7 | float の同じ状態の残差 | V3。面の流束と集めた後の残差を分けて測る (§6.5 の H3、codex diagnose 2026-10-10) **PASS (2026-10-10、§6.10)**。解釈は諮問待ち | O → F |
| 8 | 到達 | V4。**§6.11 で実施形を事前登録し直した (2026-10-10)**。見張り (`m9_watch.py`) の出力の間隔・点数の検査と末尾 9 出力の保存を直してから回す。§6.12 の A/B を先に行う **結果 (2026-10-10、§6.16): 到達、フラグあり** — θ_r +1.3〜1.8 %、Q_w +1.4 %、残差の床 2.6〜3.9 倍。解釈は諮問待ち | O → F |
| 9 | 速さと総時間 | V5 **PASS (2026-10-10、§6.17): 比 0.636** (float 20.7・FP64 32.5 ms/step) | O |
| 10 | 標準ケースの回帰 | V6 (5 ケース) **2026-10-10 改訂 (§6 V6、codex diagnose)**: 定常のケースは float と FP64 の精度の比較を主軸にし、従来の物差しも残す。最初は case/48 (§6.21、事前登録して投入) | O |
| 11 | 状態の commit の精度 (条件付き) | V4 の記録で停滞が見えたら、`qAccumulatorFP64` の軸対称対応を別 plan に (拒否の理由に向き合う) | F |
| 12 | r0・r1 のメモリ (§6.9) | 面あたり +24 B (2D の node で約 48 B/節点) で §4.4 の見積もりを超える。使う経路だけ確保する・ホストの写しを捨てる、のどちらにするか (または受け入れるか) を決める **ユーザ決定 (2026-10-10)「いいんじゃね」: 受け入れる** (case/45 でホスト・GPU それぞれ約 41 MB、e と合わせて約 72 B/節点)。大きな格子で GPU のメモリが足りなくなったら、使う経路だけ確保する形に削る | F |
| 13 | 格子の HDF5 の精度 (2026-10-10、V6 の下調べで判明) | float のビルドの変換器は座標を `stof` で読み `geom_float` (float32) で書く (`gmshReader.hpp:506-511, 2343`)。float32 の格子では cc64 は float の値そのものなので、段 ①〜④ の利点が出ない (case/56 で新旧が再実行の差の範囲で同じだったのはこのため)。利点を得るには FP64 の変換器で書いた格子が要る (case/45 は `forge-wallfit-bin-fp64` の変換器、case/48 は AWS の `run_0046/0047` のみ)。変換器を常に double で幾何を書くようにするか、手順で FP64 の変換器を必須にするかを決める。V6 の既存の基準の run の格子は 5 ケースとも float32 **判断: 2026-10-10 ユーザ決定 — 両方やる** (§4.6: 変換器は常に double で幾何を書く、ソルバの警告、手順)。合格は §6.14 **完了 (2026-10-10、§6.15、すべて PASS)** | O |
| 14 | 精度の切り替えの A/B (§6.18) | 事前登録どおり 3 本を投入し判定する **結果 (2026-10-10、§6.19): 判別不能** (θ_r 3 断面は支持の条件を満たし、FP64 に切り替えると FP64 の値に戻る。Q_w は A の揺れで条件を外れた)。**ユーザ決定 (2026-10-10)**: float がなぜ違う所で止まるか (何が効いているか) を調べる。ただし hoop の修正 ([axisymmetric-freestream-hoop-gauge](axisymmetric-freestream-hoop-gauge.md) §4.5・§4.6) を入れてから。それまで case/45 の設計の計算は FP64 | F |
| 15 | float が違う所で止まる原因の切り分け (ユーザ決定 2026-10-10、hoop の修正の後) | (1) hoop の修正 (キー 1) を入れた格子で、float と FP64 の比較 (V4 と同じ停止規則) をやり直し、差が残るかを見る **(§6.20 に事前登録して投入した、2026-10-10。FP64 の腕は hoop plan §4.6 の 7 のキー 1 の腕を共用) 結果 (2026-10-11、§6.25): 到達・フラグあり。差は V4 と同じ向き・同じ桁で残った (θ_r +1.4〜1.9 %、Q_w +2.1 %)。(2) の最初の段は §6.26 (共通の Q32 の起点からの軌跡の A/B、codex diagnose 2026-10-11) に事前登録** (差の出る場所が収縮部で、hoop の偽の力も収縮部にあるため、先に混ざりを除く)。(2) 残るなら、float に依存する部分 (状態の commit・残差の足し込み・面の値・状態から原始量・SST のソースの大きな項の釣り合い) を、float のビルドの中で部分ごとに double にする A/B で切り分ける。設計 (§4・§6) は (1) の結果を見てから書き、上位に諮る **先行の証拠 (2026-10-10、`forge-faceh` の `notes/sessions/2026-10-10-return-viscous-jacobian-close.md` (ブランチ `feature/faceh-audit-viscjac-close`)、どちらも FP64 のビルド)**: (i) 面エンタルピーの評価の精度だけを変えた探索試験 (point、2000 step × 3 本、`run_0500`〜`0537`) は「double の評価で 10 % 以上の改善」を支持しない (NOT_SUPPORT、FP32 の状態では未判定)。(ii) LHS の座標の差を段 ② と段 ③ で変えた A/B (`~/forge-faceh-audit` の `run_0560`〜`0566`): 同じ凍結入力で残差は再実行の揺れの範囲で変わらず、K・D の変化は e による係数の変更から 5 ulp 以内で再現した (段 ③ は意図どおり LHS だけを変える独立な確認)。値 3・マスク 7 の破綻はどちらでも起きる。(2) の候補の「面の値」「LHS の幾何」を外す根拠にはしない (FP64 での結果で、float の状態の丸めとは別の問い) | F |
| 16 | 粘性 Jacobian の値 2・3 を float のビルドで使う場合 (条件付き) | 粘性 Jacobian plan (faceh ブランチで accepted) の値 2・3 は FP64 のビルドでだけ確かめた (FP32 は起動も未確認)。float の状態のビルドで使うなら、先にその plan の §6.17 の再開の条件の最初の試験 (FP32 での起動と U2・U3) を行う (`forge-faceh` の `notes/sessions/2026-10-10-return-viscous-jacobian-close.md` (ブランチ `feature/faceh-audit-viscjac-close`))。**更新 (2026-10-11、forge-8a の連絡)**: 粘性 Jacobian の上限 50 の比較 (plan time_integration-line-viscous-jacobian-dt §6.4、commit 4ead4918、faceh ブランチ) は分岐 2 で閉じた。~~値 0 を既定のまま、値 2・3 は診断用に残す~~ **訂正 (forge-8a、同日)**: dt plan は既定を決めない (観測だけ: 全部入りは上限 50 で生産より速くも良くもならなかった、止まる状態の残差 4〜5 倍・θ_r/Q_w 0.7 % ずれ、各 1 本)。既定の扱いは速度 plan §5.1 #9 の判断欄と関係するセッションに委ねる (ユーザ「何を既定にするかは、ほかのセッションも関わるから決められんかも」)。float のビルドで値 2・3 を使う予定は無くなったので、この項目は再開されるまで保留。再開の条件の 1 つは「残差の評価の精度を上げた経路 (面エンタルピーの double など) が既定になる」で、§6.40〜§6.46 の float の調査 (残差の側) の結論とつながる | O |
| 17 | 段ごとの更新の診断 `FORGE_DIAG_STAGE_DUMP` (§6.30) | 環境変数で有効にする診断。平均流の保存量 5 つについて、指定の step (既定 1) の組み立て後の残差 R・commit の基準 b・緩和とガードの後の要求の更新 d・commit の直後の状態 q と、境界条件・状態方程式の上書きの後の状態を、節点ごとに double で HDF5 に書く。数値を変えないこと (同じ入力で診断の有無の 2 step の結果がビット一致、FP64 と float の両方)。float と FP64 の両方でビルドし、§6.30 の 3 回ずつの 1 step を回す台本と判定 (`stage_an.py`) を用意する **実装 (2026-10-11、implementer)**: `cuda_forge/stageDumpDiag.{hpp,cpp}`、`main.cpp` (`advanceImplicitSteady` の頭と末尾・solve の直前のフック)・`update_d.cu` (commit のラッパのカーネルの直後のフック)。Q0・Q_asm (組み立ての後・solve の直前の ro..roe)・R (res_ro..res_roe)・b (roN..roeN)・d (block DPLUR は dq_block_old_k、implicitRelax は掛け済み)・q (commit の直後)・q_final (step の末尾) を節点ごとに double で書き、q == fl(b + d) を自己照合する。float・FP64 ともビルド可、`update_d.cu` の 8 カーネルの SASS は HEAD と一致 (環境変数なしでは数値の経路は同じ)。未対応の経路 (陽解法・dual-time・qAccumulatorFP64・updateGuardAlpha > 0) は書かずに標準出力に出す | O |
| 18 | qAccumulatorFP64 の軸対称 (診断用) と採用の mask の診断 (§6.34) | **判断 (2026-10-11、codex diagnose `notes/reviews/2026-10-11-qacc-axisym-state-precision-diagnose.md`)**: 呼び出し側の案 (作業用 Q から運動エネルギーを取り `roe = (float)Qacc` で上書き) は壁と軸の両方に属する節点で誤る (Major 1) ので採らない。**実装する契約**: (i) `enforceAxisSymmetry_d` に `double* qacc_roUy, double* qacc_roe` を足す (nullptr = OFF、OFF の演算列は 1 文字も変えない)。作業用 Q と Q_N への既存の射影はそのまま。Q_N を射影する**前**に、既存と同じ `flow_float` の式と密度のガード (`rN > 0`) で `ke_N = 0.5·roUyN²/rN` を作り、ON では `qacc_roe −= (double)ke_N`・`qacc_roUy = 0` を足す。作業用 roe を Qacc から上書きしない (ミラーへの書き戻しは commit が担う)。(ii) 軸対称の拒否は全面解除しない (Major 4)。環境変数 `FORGE_QACC_AXISYM_DIAG=1` があるときだけ通し、起動時に「診断用、rollout plan の S1b-① の合格条件は未検証」と出す。無ければ今の拒否のまま。(iii) 採用の mask の診断 `FORGE_DIAG_QACC_ADOPT=<接頭辞>`・`FORGE_DIAG_QACC_ADOPT_STEPS=<step,...>`: 指定の step の `updateVariablesOuter` (reconcile) の直前に Q と Qacc を host に写し、量ごとに Q ≠ (float)Qacc の mask・節点の座標・軸/壁/境界の flag を h5 に書く (host だけ、カーネルを変えない)。(iv) 単体試験 `tests/unit/test_qacc_axis.cu`: (a) rollout の反例 (ρ=1, ρu_r=0.125, ρE=2, dq=0 → 軸ピン → commit で Q・Qacc とも ρu_r=0, ρE=1.9921875、旧い作り (Qacc に当てない) では 0.125・2 に戻ることも示す) (b) 1 ULP 未満の残余の保持 (c) 残余ゼロでピン後の (float)Qacc が OFF の Q_N と**生のビット**で一致 (乱数 20 万セル) (d) 壁→軸→commit の順序で壁と軸に重なる節点 (ρ=1, ρu_x=1, ρu_r=0.125, ρE=2 → ON・OFF とも commit 後 ρE=1.9921875) (e) 軸でない節点は不変 (f) dq ≠ 0 (g) 残余を持つ複数 step。(v) `cuobjdump -sass` で `enforceAxisSymmetry_d` 以外のカーネルが HEAD と一致することを確かめる。合格条件: 単体試験の全項目 PASS、float・FP64 のビルド成功、SASS の照合。rollout plan の S1b-① (case/44 の回帰・改善ゲート) はこれで済ませない。**実装 (2026-10-11、implementer)**: 軸ピンを `cuda_forge/axisSymmetryPin_d.cuh` (`axisSymmetryPinCell`) に切り出し、ON では Q_N を射影する前の `ke_N` を同じ式・同じガードで作り直して `qacc_roe` から引き、`qacc_roUy = 0` とした (OFF の行は HEAD の文字のまま。共有すると OFF の SASS で負号の置き場所が変わるので作り直す形にした)。`main.cpp` の拒否は `FORGE_QACC_AXISYM_DIAG=1` のときだけ外し、毎回「診断用」と出す。`cuda_forge/qaccAdoptDiag.{hpp,cpp}` (h5 は量ごとの `/<q>/{mask,Q,Qacc}` と `/cells/{cc64,axis_flag,wall_flag,iso_wall_flag,bnode_flag}`、属性に内訳の数)。単体試験 `tests/unit/test_qacc_axis.cu` (a)〜(g) PASS (旧い作りと却下した案が反例で外れることも示した)、既存の `test_qacc_commit` PASS、float・FP64 のビルド成功。SASS は 292 関数のうち 290 が HEAD と一致、残りは `enforceAxisSymmetry_d` (シグネチャ) と `dependentVariables_d` (nvcc の揺れ、HEAD の再コンパイルの 2 通りの一方と一致)。**分かった性質**: 軸の節点では Qacc が `roeN − ke_N` を float に丸めずに持つので、残余ゼロから始めても dq ≠ 0 の 1 step で ON と OFF の ρE がビット一致しない (単体試験 (f) で 20 万中 22,683、ON は常に正本に最も近い float)。軸対称では「残余ゼロなら ON は OFF とビット同一」を補助ゲートに使わない | O |
| 19 | 状態の量子化の診断 `FORGE_DIAG_QUANTIZE_STATE` (§6.36) | 環境変数 `FORGE_DIAG_QUANTIZE_STATE=1` で、定常の陰解法の各 step の末尾 (`advanceImplicitSteady` の `applyNodeIsothermalWallPin` の後・`updateVariablesOuter` の前) に、保存量 9 つ (ρ・ρu・ρv・ρw・ρE・ρk・ρω・ρY_s すべて) を全 CV (`nCells_all`) で float32 に丸めて double に戻すカーネルを 1 本足す (読み書きは状態の配列だけ、ほかの補正を混ぜない)。**FP64 のビルド専用** (float のビルドでは何もしないので起動時に拒否)、定常の陰解法 (`timeIntegration 11`・`unsteady 0`) 以外と `qAccumulatorFP64 1` との併用も拒否。自己点検: 環境変数 `FORGE_DIAG_QUANTIZE_STATE_CHECK=<step,...>` (既定 1) の step で、`updateVariablesOuter` の後に Q と履歴 (Q_N・Q_M、ρk・ρω の N・M) を host に写し、float32 で表せない値の個数を量ごとに出す (0 であること)。ρY の履歴 (`roY_N`) は次の step の化学種の更新の頭で ρY から取られるので、その step の化学種の更新の直前にも ρY の個数を数える。合格条件: FP64 のビルド成功、float のビルド成功 (環境変数なしで数値の経路が同じ)、`cuobjdump -sass` で新しいカーネル以外が HEAD と一致 (FP64・float とも)、ローカルの短い確認は不可 (run は回さない) なので、自己点検の出力の書式を単体で確かめる。**実装 (2026-10-11、implementer)**: `cuda_forge/quantizeStateDiag.{hpp,cpp}`・`quantizeStateDiag_d.{cu,cuh}`、`main.cpp` のフック (丸め・(a)・(b)×2・init・finalize)。単体試験 `tests/unit/test_quantize_state.cu` PASS (20 万個でデバイス = (double)(float)x のビット一致・冪等・行の書式。`-ftz=true` では FAIL して試験が効くことを確認)。float・FP64 のビルド成功。SASS は新しいカーネル以外が HEAD と一致 (`dependentVariables_d` は HEAD のソースの再コンパイルの 2 通りの一方と一致、float・FP64 とも)。**未決**: (b) の roY_N は node の入口のピンの節点 (`speciesTransport_d.cu:136,145` が step の中で double で書く) で 0 にならない見込み。§6.36 を再開するときは、前提 3 をピンと ghost を除いた個数で見るかを run の前に決める | O |
| 20 | 状態の量子化の診断で丸める量を選ぶ `FORGE_DIAG_QUANTIZE_STATE_ONLY` (§6.40) | `FORGE_DIAG_QUANTIZE_STATE_ONLY=<量,...>` (例 `roe`) で丸める対象をその量だけに絞る。対象に無い名前は起動時に拒否、未設定なら #19 のまま。自己点検 (a) は選んだ量とその履歴だけを見る。自己点検の step で、丸めの直前・直後の選ばなかった保存量をビット比較し、変わった値の個数を出す (0 が期待、codex Minor 6)。合格条件: FP64・float のビルド、単体試験 PASS、デバイス側の SASS が #19 の commit と同じ (host だけの変更)。**実装 (2026-10-11、implementer)**: 名前の検査 (空の値・空の要素・対象に無い名前・重複・前後の空白を拒否)、履歴を選んだ量に絞る、自己点検の step で丸めの直前・直後の選ばなかった量の生のビット比較 (`[quantizeState] step N: 選ばなかった量の丸め前後の変化 … / 合計 0`)。`quantize` に `iStep` を足した (呼び出しは main.cpp の 1 か所)。単体試験 (5)(6) を足して PASS (変異を入れると FAIL することも確認)。SASS は float・FP64 ともデバイス側が HEAD と一致 (差は nvcc の揺れの `dependentVariables_d`・`cond_moment_update_limited_d` だけ、基準のソースの再コンパイルの 2 通りの一方と一致)。未設定では #19 と同じ | O |
| 21 | ノズル上流 (入口〜収縮部、M 0.016〜0.08) の低マッハ補正の検討 (ユーザの見立て 2026-10-11「低マッハ補正はノズル上流で何らか入れるべき」。方針は未決) | 使える選択肢と制約: `lowMachPrecond 1` (流束の散逸だけの Weiss–Smith 前処理、`precondEps` 既定 0.15、line implicit と併用可)、`lowMachPrecond 2/3` (完全な前処理、line implicit と併用不可)、`lowMachThornber 1` (再構成の補正、前処理と併用可)。SLAU は全速度向けの散逸を持つので重ねたときの挙動は未検証 (深いすきまで `lowMachPrecond 2` が非物理な貫通流を作った前例)。離散化を変えるので、(i) FP64 の解 (スロートの流量・境界層・θ_r) の変化と (ii) float と FP64 の差への効き方を別々に測る。§6.46 (pRef の長い試験) の後に、設計 (§4・§6) を上位に諮ってから試験を組む | F |
| 22 | 生産の陰解法 (B0) が細かい x の格子で途中の状態に止まる疑い (forge-25 の観察 2026-10-11、原因は未確認・上位にも諮っていない) と float の調査の関係 | forge-25 の記録 (feature/nozzle-core-grid、plan tooling-nozzle-core-grid §4.27〜§4.29、commit b1249258): 参照の格子 Gref (4719 × 160) で B0 は残差がプラトー (rms_ro 1.6e-5) で量は STEADY、同じ場から point (lineImplicit 0・lineDtDirectional 0・ITJ 0) に替えると残差は 12 万 step で 2.4 桁下がり続け θ_r −0.7 %・Q_w −0.44 % 動く。x を粗くした Grefx (3086 × 160) では B0 で 5 桁下がり point と 0.002 % 以内。G0 (4719 × 121、float の調査の格子) も B0 でプラトー。**float の調査への含意 (仮説、未検証)**: §6.26 以降の試験は同じ B0 での float と FP64 の比較なので精度の差の試験としては有効だが、FP64 の参照自体がプラトーの途中の状態かもしれない。float のずれが「B0 が減衰させない遅い方向に丸めの揺れが溜まる」現象である可能性がある。§6.46 の後に、point で float と FP64 を比べる試験を組むかを上位に諮る。**続報 (forge-25、同日、plan tooling-nozzle-core-grid §4.27〜§4.31、commit 12869fc7、判定 `_band_ab/core_grid/grefx_refpoint_judge.json`)**: Gref を point で通算 20 万 step 解くと残差は 2.4 桁下がって (rms_ro 7e-8) 横ばい、8 量 STEADY (2 万 step の窓の差 ≤ 0.01 %)。止まった値は B0 の窓より θ_r −0.57〜−0.63 %・Q_w −0.51 %・δ_loc −0.03〜−0.04 %。Grefx は B0 と point が 0.002 % 以内で同じで、止まった Gref との差は θ_r ≤ 0.05 %・Q_w 0.04 %。よって B0 は ni 4719 の格子で止まるべき値から θ_r・Q_w で 0.5〜0.6 % 離れた所に止まる (原因の設定は未確認)。G0 は未測定。§6.48 (G0 の point での float と FP64 の比較) はこの観察と直接つながる (G0 の point の値の参考にもなる)。**続報 2 (forge-25、同日、plan tooling-nozzle-core-grid §4.32・§4.33、commit d51af1d0、判定 `_band_ab/core_grid/line_cap_judge.json`)**: Gref の B0 の場 (run_0608 res_60000) から `lineDtDirectionalCap` だけを 50 → 1 にした A/B (各 2 万 step、run_0615・0616): 上限 1 は壁際の擬似 dt を 1/50 (point と同じ) に縮め、平均流の残差は上限 50 の 0.6〜4 % まで下がり、θ_r は point の止まった値へ動き始めた (差 0.63 → 0.29 %)。上限 50 はプラトーのまま。上限 1 は ω の残差が 3〜13 倍で check_convergence は末尾で RISING、Q_w は一時的に離れる向き。B0 の停滞は方向別の擬似 dt の伸ばしと関係している読み (原因は未確定、forge-25 は上位に諮る前)。**float の調査の次の候補**: §6.48 の段 2 の後に、ラインは残して上限 1 にした float と FP64 の比較 (方向別 dt の伸びの寄与を分ける) を上位に諮る。**続報 3 (forge-25、同日、plan tooling-nozzle-core-grid §4.34・§4.35、commit 1045449d)**: (1) 上限 (lineDtDirectionalCap) は setDT_d.cu でしか使われず、dt_local は残差の組立に入らない (更新・左辺の側だけ) → 上限・方向別 dt は右辺を変えず、B0 と point の止まる値の違いは反復の側。(2) Gref の B0 の場から上限 50/1 で 2032 step 回し最後の 33 step を毎 step 書き出すと、交互に揺れる節点は 75.5 万中 20〜38、残差の上位の節点は毎 step ほぼ同じ向き・大きさの更新 (lag-1 +0.9998)、dq_block_old と実際の ΔQ のずれ 1e-11、commit の n_clip・n_lost 0 → 揺れでも床でもなく「遅いモード」。(3) 残差の上位 0.1 % (755 節点) が ρ・ρE の残差二乗和の 67 % を担い、全部が縮流部 x < −2 r_t の主流側 (壁から 40〜160 節点)。ω の残差は壁の第一内部節点のスロート付近。(4) 仮説: 流れ方向に細かい縮流部で擬似時間のモードが抜けにくい (Grefx では B0 で 5 桁下がる)。**float の調査との対応 (仮説、未検証)**: float と FP64 の差は縮流部〜スロートで先に・最大に育つ (§6.38・§6.47) こと、point では float と FP64 が一致する (§6.49) ことと合わせ、「B0 の反復が縮流部の遅いモードを減衰させず、float の丸めの揺れがそのモードを押して別の所へ流す」読みが有力 (原因の確定は上位に諮ってから)。**続報 4 (forge-25、同日、plan tooling-nozzle-core-grid §4.36・§4.37、commit 13a5522a)**: 格子 Grefc = Gref の縮流部 x −8〜−2 r_t だけ流れ方向の密度を max(1, d/2) に並べ直した (列 589 → 337、ni 4467、区間の外はビット一致)。同じ B0 の場から B0 のまま各 6 万 step: 縮流部の単位体積あたりの残差の RMS は対照の ρ 0.027〜0.028・ρE 0.016・ρk 0.026・ρω 0.31〜0.33 倍、全域の rms_ro 4.3e-7。Grefc の量 (40k〜60k) は point で止まった Gref から θ_r −0.05 %・Q_w +0.04 %・M ≤ 0.001 (B0 の Gref は +0.6 %・+0.5 %)。機構は未確認。**float の調査の次の候補 (最有力)**: G0 の縮流部だけを同じように粗くした格子 (G0c) で、B0 のまま float と FP64 を比べる (陰解法を変えずに格子だけで FP64 の停滞と float の差の両方が消えるか)。**ユーザ決定 (2026-10-11)**: G0c (G0 の縮流部 x −8〜−2 r_t だけ流れ方向の密度を max(1, d/2)、ni 4467、forge-25 の `core_grid_mesh.py` の Grefc と同じ station) は AR の最大が 8,000〜9,000 になる見込みで、**float の診断試験の系列に限り AR ≤ 10,000 を許す** (forge-25 と同じ扱い、`check_mesh_quality.py --ar-max 10000`)。生産に使うかは別に決める。**訂正 (forge-25、2026-10-11、plan tooling-nozzle-core-grid §4.35 の訂正・§4.38、commit 2859b299)**: 続報 3 (3) の「主流側 (壁から 40〜160 節点)」は誤り。Gref は境界層の中に節点を多く入れているので、壁から 60〜100 節点目は物理的には乱流・温度境界層の縁 (冷却壁の温度境界層は縮流部で厚く、壁から 70〜80 節点目で 1600 K に届く。島はその温度の前線に重なり、μt/μ が 50〜150 から 8〜11 に落ちる所。M 0.02〜0.08)。モードの主役はエントロピー (s の RMS が圧力の 18 倍)。速度の変化は境界層の縁から主流へ斜めの縞になってスロートまで伸びる。G0c の区間の作り方には影響しない。**forge-8a の観測と仮説 (同日、セッション間のメッセージ。forge-8a の plan dt §6.4 にも記録、上位には諮っていない)**: G0 の生産の開始場 (case/45 `run_0571_lvc75_m5_b1` の res_0) で壁から 40 節点より内側の M は x/r_t −8 で 0.004〜0.041、−4 で 0.039〜0.078、−2 で 0.13〜0.23。case/45 は SLAU・`lowMachPrecond 0`・`lowMachThornber 0`。仮説: M ≈ 0.02 では擬似時間の dt が u + c で決まり、流れに乗るエントロピー・渦の誤差は 1 step に CFL·M ≈ 0.1 セルしか進まない (縮流部を流れ方向に粗くすると停滞が消えることと向きは合うが、大きさは説明しない)。ラインの陰解法・ITJ・`lineViscCoupling` 2/3 は `lowMachPrecond` ≥ 2 を起動時に拒否する (forge-8a の確認、本 plan では未確認)。#21 の材料。**§6.52 (G0c で B0 のまま float と FP64 を比べる) を事前登録 (2026-10-11)**。forge-25 の確認 (同日): Grefc は runner の経路で `msh_digits` 17 で書かれ、区間の外のビット一致は実際の Gref の入力と照合したもの (G0c の 10 桁の欠陥は本 plan の `make_g0c.py` だけ)。**続報 5 (forge-25、同日、plan tooling-nozzle-core-grid §4.39・§4.41)**: Gref・FP64・B0 で ITJ だけ 5 → 0 にした `run_0627_cg_gref_itj0` (AWS A、対照 `run_0619`): 縮流部 x∈[−8,−2] の残差の比 (ITJ 0 / 対照、窓 10k〜15k・15k〜20k) は ρ 0.97・0.94、ρE 0.84・0.83、集中の場所 (境界層の縁) も同じ → forge-25 の事前登録の結果 B (熱の Jacobian は停滞の原因でない)。全域の SST の残差は ITJ 0 で 1.6〜2.8 倍。SST の分離更新 (§4.39) と ITJ (§4.41) が外れ、停滞は平均流のブロックの反復に残る (forge-25 の読み)。実装の観察 (forge-25): `implicitRelax` は各 sweep の中で掛かる (`timeIntegration_d.cu` 1980 付近、dq_new = relax·dq)、無限回 sweep の解は (L/0.7 + N)⁻¹。**続報 6 (forge-25、同日、plan tooling-nozzle-core-grid §4.43〜§4.45、commit 7ba8ba1e)**: 内部反復 (nStepInner 5) の解き残しは原因でない。停滞場 (Gref・FP64・`run_0619` res_60000) から 1 外反復を nStepInner 41 で回し `FORGE_LINE_DUMP` で縮流部の 289 ラインを書き出した (`run_0628〜0631`、上限 50/1 を各 2 回): 5 回目の補正は 40 回目と島 (境界層の縁) で 0.04〜0.16 % しか違わず、後退誤差 1e-15・40 回目の線形残差 ≤ 7e-8。**訂正**: 以前の「上限 1 は停滞の場所の dt をほぼ変えない (0.977)」は全節点の中央値。島の列の高さ別では、B0 の擬似 dt は point の 50 倍 (壁から 1〜49 節点目)・13.6〜3.3 倍 (60〜79)・2〜1.3 倍 (80〜89)・1 倍 (90 以上)。残差の山は壁から 70〜84 節点目で、温度の前線 (1343 → 1592 K)・μt の落ち込みと、B0 の dt の拡大率が 5 倍から 1 倍へ落ちる帯が重なる (重なりは観察、因果は未確認)。forge-25 は次に左辺と右辺の Jacobian のずれ (十分解いた補正で実際の残差が予測どおり減るか) を測る案を諮っている。本 plan の §6.51 の前提 2 (dt の縮小) は壁側 j ≥ 100 (G0 で壁から 20 節点以内) と全域の中央値で見るので、この訂正で登録は変えない。**続報 7 (forge-25、同日、plan tooling-nozzle-core-grid §4.46〜§4.49、commit f11c7dcb)**: 停滞場 (Gref・FP64、SST 凍結) で B0 の補正 δ を 1 回当てたときの残差の応答: 停滞帯 (境界層の縁、壁から 60〜89 節点目) では δ が大きい (δρ/ρ_ref の RMS 3e-3〜1e-2) のに、実際の残差の変化の R₀ の向きへの射影 ⟨ΔR, R₀⟩/‖R₀‖² は 1 step あたり −0.1 % で、変化はほぼ R₀ と直交 (残差は減らずに向きを変える。**訂正 (forge-25、同日)**: ノルムの減り方は √(1 + 2p + r²) で、停滞帯の ρ では約 0.04 %)。左辺 (1 次 FVS とスカラーの粘性の対角) は 4〜29 % の減少を予測。右辺を 1 次の再構成に替えても、ずれ e のノルムは 90〜95 %、R₀ への射影は 83〜94 % が残る (**訂正 (forge-25、同日)**: 最初の「ずれの 2 割弱を説明」はベクトルの差のノルムを内訳として書いた誤り。境界層の外ではずれのほぼすべてが再構成)。`run_0619` の 2 万〜6 万 step の残差は平均が動かず ±20 % の不規則な振れ (時間尺度 2000〜10000 step) で、止まり方は有界な揺れ。残りの候補 (流束 SLAU と左辺の FVS の違い = 低マッハの散逸 / 粘性の左辺の近似 / μt) を forge-25 が切り分け中 (次は左辺の粘性スカラー対角の作用を差し引く診断、core-grid §4.50)。forge-25 は本 plan §6.51 の結果を core-grid §4.50 の補助情報に記録した。診断 `FORGE_LINE_DUMP_FINAL` (commit fb17dd5b、既定 off)。**続報 10 (forge-25、2026-10-12、plan tooling-nozzle-core-grid §4.60・§4.61、commit e4271699)**: 本 plan §6.53 の増幅を forge-25 の停滞と同じものと読む (1 step あたり約 0.37 % の増幅、G0c では粗くした区間のすぐ上流の細かいままの所に場所が移る = 流れ方向に細かい所で増幅する手がかり)。B0 にライン面の熱伝導の近傍結合だけを足す (`lineViscCoupling 3`・`FORGE_LVC_TERMS=2`、Gref、2 万 step) と停滞は緩和せず (縮流部の残差の比 0.89〜1.15)、止まる水準が上がり (rms_ro 3.8e-5 対 2.1e-5)、θ_r +0.86 %・Q_w +0.9 % と単調に動き続けた。左辺の 1 項の整合では増幅は止まらない読み。forge-25 は次に格子間隔ごとの増幅率 (同じ入力の 2 本の差の成長率) を測る。**続報 8 (forge-25、同日、plan tooling-nozzle-core-grid §4.51・§4.54、commit d71416a9)**: 停滞場で補正に対する残差の応答と左辺の予測のずれを面の種類で分けた: 流れ方向 (ライン外) の面の対流は左辺と右辺の応答がほぼ一致。運動量のずれの大部分は壁法線 (ライン面) の対流 (左辺の閉形式 FVS と右辺の SLAU の違い、‖E − Eₗ‖/‖E‖ は ρu 0.23・ρv 0.18)。連続の式のずれは壁法線の粘性スカラー対角。エネルギーは粘性対角で約 1/3、残りは未特定 (熱伝導の近傍結合が候補)。近似の左辺の誤差は壁法線の方向にあり流れ方向には無い。「流れ方向に細かいと止まる」とのつながりは未確定 (forge-25 の候補: 方向別 dt の拡大率が 5 → 1 に落ちる高さが流れ方向の格子間隔で動き、Gref では温度の前線と重なる)。forge-25 の判定は自身の閾値が厳しく形式上は判定不能 (諮問中)。B0 の左辺は `implicitSolvePrecision 0` で単精度で組まれる (要素の丸め 1e-7、本 plan §6.32 の感度試験と関係)。**続報 9 (forge-25、同日、plan tooling-nozzle-core-grid §4.57、commit 54fa80b6)**: エネルギー行の予測のずれの主要部分は、B0 の左辺に無い壁法線 (ライン面) の熱伝導の近傍結合。固定 δ への予測にその作用 (凍結係数) を足すと停滞帯の ρE のずれは 0.31 倍 (射影 0.14 倍、forge-25 の事前登録の結果 A)。停滞帯の左辺の読み違いは 3 成分とも壁法線の方向の近似から来る: 連続の式 = 粘性スカラー対角 (真の式に無い拡散)、運動量 = ライン面の対流 (FVS と SLAU の違い)、エネルギー = 粘性スカラー対角 (約 1/3) と熱伝導の近傍の欠落。「流れ方向に細かいと止まる」経路は未確定。固定 δ の予測の診断であり、forge-8a の「熱伝導の近傍を反復に入れると悪化」(`lineViscCoupling 3`・マスク 7) とは操作が違うので「入れれば治る」とは言えない | F |

## 6. 検証 (事前登録。数値は codex の提案する判定基準で、実測で裏付けた許容差ではない)

- **ビルド**: float と FP64 の両方が通る。
- **V0 状態を固定した演算の A/B** (段 ① の後、因果の切り分け。到達性能の合格ではない):
  - 同じ保存場から状態・勾配・物性・面ベクトルを一度だけ float にそろえ、3 つの腕で固定する。
    - 旧腕: 今の絶対座標の差。
    - 新腕: double の座標の差から作った float の係数 (符号・ガード・面の選択は変えない)。
    - 参照: 同じ float の状態・勾配・物性を double に広げ、double の座標から作った幾何で評価する。元の FP64 の状態を直接使って入力の丸めを混ぜることはしない。
  - 見る量: 面ごとの粘性応力・熱流束・k/ω の拡散と、その残差への寄与。第一層、壁から 3 層、全域を分けて集計する。ρ の残差は主判定に使わない。各腕 3 回。
  - **結果 A**: 熱・運動量・k/ω の各対象で、新旧の誤差の比 ≤ 0.5、かつ再実行の差を十分上回る改善 → 第 1 仮説を支持して段 ②〜④ へ進む。
  - **結果 B**: 誤差の比 > 0.5 → その量について「直接の座標の差が過半を説明する」を棄却し、上位に諮る。
  - 旧の誤差が再実行の差と同程度なら判別不能。
- **V1 係数と接続**:
  - 非退化面の |e| と、各用途の拡散係数が FP64 比 ≤ 1e-6。ベクトルは `‖e_new − e_ref‖/‖e_ref‖` で判定する。
  - 退化面は相対比較せず、有限値・使わない扱いを確かめる。
  - 周期の相手・LSQ の同値類・ラインの接続は参照と一致する。
  - closure: `Σ±S_device − A_closure` を、局所の `A_planar` と、打ち消し前の `Σ|S_device|` の両方で正規化して記録する。一様圧の試験は `P = pRef` と `P ≠ pRef` の両方 (前者だけではゲージが誤差を隠す、`axisymmetricSource_d.cu:51`)。
- **V2 FP64 の同一性** (各段の後):
  - 同じコミット系列・コンパイラ・フラグ・入力・実効設定で比べる。
  - 各保存量について、節点の残差の差を「参照の各面流束・源項の絶対値の和」で正規化し、L1 と最大値の両方で評価する。閾値は `max(1e-12, 旧版の再実行の差 × 10)`。
  - 面流束・残差・dq を別々に保存する。20 step 後の状態は、固定した物理の尺度で正規化して ≤ 1e-10 を補助条件にする。
  - 分岐・ガードを変えた場合は、丸めの差として許容しない。
- **V3 同じ状態の残差** (float 対 FP64):
  - すべての保存量・組成・乱流量を共通の Q32 にそろえ、FP64 側にも double(Q32) を渡す。比べるのは commit 前の残差の場。
  - 全域と壁際で、各保存量の正規化した誤差が `E_new ≤ 1.1 E_old + 再実行の差 × 10` (悪化しない)。
  - 旧の誤差が再実行の差の 10 倍を超える熱・運動量・乱流の対象では、少なくとも半分になること (改善)。保存量を 1 つのノルムに混ぜない。
- **V4 到達** (**2026-10-10 に §6.11 で置き換えた。以下は元の事前登録で、窓の条件は撤回** — 評価の目的を「定常の同等性」から「停止規則のもとでの比較」に変えたため。codex diagnose 2026-10-10) (float の新しい経路、case/45 の B0 の構成、`run_0183` の res_100000 から、2 出力連続の水準まで、2500 step ごと、最大 20 万 step):
  - θ_r (断面ごと) と Q_w (符号つき) を、両方の run の同じ設定の区間の末尾 2 万 step (9 点以上) で比べる。分母は参照の窓の平均。
  - 両方の窓のドリフトと振幅が閾値 (θ_r 0.05 %・Q_w 0.1 %) の 1/5 以下で、平均の差に両方の時間変動の幅を足しても閾値以内。単調なら漸近値も書く。
  - 抽出した系列 (末尾 2 万 step に切り出した窓) を `check_quasisteady.py --series-csv --tail 1.0` に渡す。θ_r は `--drift 0.0001 --osc 0.0001`、Q_w は `--drift 0.0002 --osc 0.0002` を明示する (閾値 0.05 %・0.1 % の 1/5。既定の drift 0.05 は 5 % で緩すぎる、codex plan 段 m1)。
  - `check_convergence` が NOT CONVERGED なら、結論は「指定した水準への到達の比較」に限る。PASS・対象量の STEADY・判定区間がそろわない限り、「FP64 と同等に収束」とは書かない。
  - FP64 側の参照の窓: B0 の到達の場は消えているので、FP64 の B0 の参照の run は同じ記録の形 (2500 step ごと、末尾 2 万 step の場を残す) で回し直す。
  - 記録: `dq_requested` と `Q_after − Q_before` の統計 (§4.5)、壁の熱流束の分布の差、step 数。
- **V5 速さ**: 専有・交互に 1000 step × 3 本。新しい float が FP64 の 0.70 倍以下。採用の判断は V4 までの総時間でも行う。
- **V6 標準ケースの回帰** (`case/08` は使わない。`procedures/verification/README.md:23` のとおり回帰の基準にならない):
  - **2026-10-10 の改訂 (ユーザの目標「double をやめたい」、codex diagnose で採否)**: 定常のケースは、新しい変換器 (double の幾何) の同じ格子・同じ初期保存量で float と FP64 を並べる精度の比較を主軸にする (最初は case/48、§6.21)。新しい float と新しい FP64 の一致だけでは両者に共通の退行を見逃すので、下の表の物差し (Sod の厳密解、TG の保存則と時刻歴、凝縮の成分の収支、段 ①〜④ で通らない経路の回帰) は残す。非定常のケースを「同じ所に止まるか」で評価しない。

| 対象 | 主な物差し |
| --- | --- |
| `case/05.sod_shock_tube` | 同じ物理時刻の厳密解に対する ρ・u・p の L1、逸脱・正値性 |
| `case/09.Taylor-Green` | 周期の質量・運動量・エネルギーの収支、KEEP のエネルギーの履歴、受動種の収支 |
| `case/48.flat_plate_cooled_m4` | SST、壁の熱流束・摩擦・θ の分布、全残差 |
| `case/44.vitiated_air_wt` | 軸対称・TP・凝縮の成分の収支、出口の量、凝縮量 |
| `case/16.nozzle_wys` の既存の 3D node 構成 | 角・壁の代表点・3 方向の変位、熱流束と全残差 |

  - 時系列か準定常の窓で、新しい float の参照からの誤差が `旧 float × 1.05 + 再実行の差 × 10` 以下。Sod は FP64 より厳密解を優先する。ほかの FP64 は比較の参照であって真値ではない。
  - 未実行の壁関数・WMLES などは、係数の検査だけで実行の回帰を済ませたことにしない。
  - 各ケースの基準の run・固定する物理時刻 / 窓・再実行の差は、V6 を回す前に表に書き足す。
  - **2026-10-10 の下調べ** (Explore): 既存の基準の run の格子は 5 ケースとも float32 の HDF5 (float のビルドの変換器)。そのままでは V6 は「新しい経路が壊していないか」の回帰の確認にしかならず、利点は測れない (§5.1 #13)。候補の run: case/05 `run_newdef_run_0020_limA_prod` (200 step、Toro の厳密解との L1)、case/09 `run_0174_r2_keep_eul_new` (KEEP 非粘性、AWS のみ) と `run_0156_passiveG_fct_step_seam` (受動種 FCT)、case/48 `run_0025_B_tw300_y3_fx05` (`cooled_plate_eval.py`、FP64 の格子は AWS の `run_0046/0047`)、case/44 `run_0490_limA_va3` (4000 step の再開、凝縮)、case/16 `run_0300_corner_node3d_sst_coarse` (3D node、変換し直しが要る)。V6 の規模と、FP64 の格子で回すかどうかは V4・V5 の結果を見てから決める (plan 未確定)。

### 6.2 段 ① の結果 (2026-10-10)

- **実装** (implementer): `mesh` に double の幾何の写し (`coord64`・`surfVect64`・`surfArea64`・`planeCent64`・`cellCent64`・`cc64`)、`variables` に面の配列 `ge_x/ge_y/ge_z` (e32 = fl(cc64[ic1] − cc64[ic0])、全面)。
  - V0 の評価の経路 `cuda_forge/geomAbDiag.{hpp,cpp}` を追加した。`FORGE_DIAG_GEOMAB_DUMP` / `_REF` で起動する。`viscousFlux_d` と SST k/ω の `scalar_diffusion_multi_d` に、省略できる引数 `ge_*`・`faceFlux` を足した。本番は nullptr。
  - AWS でビルドした (float sha a874b748…、FP64 sha 6059b99c…、`~/forge-fgeom-{f32,fp64}`)。台本は `case/45.isobutane_m6_d155/fg1.sh`、集計は `fg1_an.py`・`t1mag.py`。
- **既定の経路の同一性**:
  - 新旧のバイナリで 20 step: FP64 case/45 は `run_0390〜0392_fg_id64_*`、float case/56 は `run_0079〜0081_fg_id32_*`。
  - ビット単位では判定できない。node でも、原子加算の順序で、旧どうしの再実行でも場が 1,684 万バイト違う。
  - 新旧の差は、再実行の差と同じ桁だった。
    - 場の最大の相対差: FP64 で 5e-8〜7e-7 (再実行 9e-8〜5e-7)、float で 3e-7〜9e-6 (再実行 3e-7〜3e-5)。
    - 残差履歴の相対差の最大: FP64 で 8.9e-6 (再実行 3.4e-6)、float で 7.3e-3 (再実行 4.3e-3)。
  - FP64 のビルドで e32 と今の差 fl(cc1) − fl(cc0) が違う面は 0 (取り決めのとおり)。
- **V1** (B0 の最終状態、`run_0393_fg_v0_dump`):
  - ‖e − e64‖/‖e64‖ は、新 (e32) で最大 5.95e-8 (99 % 点 5.25e-8)。
  - 旧 (float の座標の差) では最大 0.299 (99 % 点 0.062) で、1e-3 を超える面が 185,855。
  - dump の e64 と REF の e は完全に一致した。
- **V0** (B0 の最終状態、float の旧腕・新腕 `run_0393`・`run_0394` (再実行) と FP64 の参照 `run_0395`、面の流束の相対 L2 誤差):

| 量 | 第一層 旧 / 新 | 壁から 3 層 旧 / 新 | 全域 旧 / 新 |
| --- | --- | --- | --- |
| 粘性の運動量 | 2.19e-2 / 9.8e-8 | 2.00e-2 / 9.8e-8 | 6.34e-3 / 9.7e-8 |
| 粘性のエネルギー | 3.01e-2 / 8.6e-8 | 2.72e-2 / 8.6e-8 | 8.82e-3 / 1.06e-7 |
| 熱流束 | 3.04e-2 / 8.2e-8 | 2.78e-2 / 8.3e-8 | 1.12e-2 / 9.6e-8 |
| k の拡散 | 1.29e-2 / 8.6e-8 | 1.15e-2 / 8.0e-8 | 4.43e-3 / 8.1e-8 |
| ω の拡散 | 1.00e-1 / 8.1e-8 | 1.00e-1 / 8.1e-8 | 1.00e-1 / 8.1e-8 |

  - 新旧の誤差の比は、すべて ≤ 1e-5。再実行の差 (面の流束) は 0。**結果 A**: 第 1 仮説 (座標の差の丸めが壁際の粘性・拡散の誤差の過半を説明する) を支持する。
  - 規則どおり、段 ②〜④ に進む。
  - 診断の照合で、旧腕と本番の残差は 61 万件で一致しない。旧腕の再実行どうしでも同程度 (61.3 万) 一致しないので、原子加算の順序の差で、面の流束の不一致は 0 である。

### 6.3 事前登録: 段 ② の合格条件 (2026-10-10、段 ② の実装の前)

- **対象**: LSQ の係数 (`gradLSQ` 1/2。case/45 は 2 = 前処理の係数 `cInt`)、周期の相手の対応付けと LSQ の継ぎ目の同値類、ラインの接続、壁の代表点 `(irep, y)` (境界面ごと)、軸対称の closure、弱形式の等温壁の d1/d2、`delta_les`。
  - どれも読み込み時に double の入力から作る。カーネルが毎 step 読む値の型は `flow_float` のまま。
- **合格条件**:
  1. **接続の一致** (新しい float のビルド 対 FP64 のビルド):
     - ラインの接続 (`line_offsets`・`line_cells`) が case/45・case/56 で一致する。
     - 周期の相手と LSQ の継ぎ目の同値類が case/39 で一致する。
     - 壁の代表点 irep が、壁関数を使うケース (case/48 の格子) の全境界面で一致する。
     - 古い float のビルドとの不一致の件数も記録する。
  2. **係数の精度**:
     - LSQ の係数 `cInt` の、FP64 に対する相対差 (係数ごと、`|c_new − c64| / max(|c64|, 1e-30·max|c64|)`) の 99.9 % 点が ≤ 1e-5 (float の保存の範囲)。旧 (丸めた座標の差) の値も記録する。
     - 壁の代表距離 y の相対差 ≤ 1e-6。d1/d2 も同じ。
  3. **closure の非劣化** (軸対称、case/45 の格子、一様な静止場 P = 1.0e5 Pa・T = 300 K・u = 0、`pRef` = 0 なので P ≠ pRef):
     - 1 step 目の `res_roUy` を `P·A_planar` で割った値の最大と RMS が、新しい float で旧 float の 1.1 倍以下。FP64 では新旧の差が再実行の差の 10 倍以下。
     - 10 step 後の最大の |u| (偽の加速度) も、新しい float で旧 float の 1.1 倍以下。
  4. **FP64 の同一性 (V2)**: 新旧の FP64 で case/45 の B0 を 20 step。場の最大の相対差が、旧どうしの再実行の差の 10 倍以下 (段 ① の確認と同じ形)。
  5. 段 ① の既定の経路の確認 (float の case/56 で 20 step) を、段 ② の後も同じ形で回す。新旧の差は、再実行の差の 10 倍以下。
  - **3' (2026-10-10 追加、3 の結果を見た後・3' の実行の前に登録)**: 3 の run は case/45 の既定 (`hoopAreaFromClosure 0`・`axisRFloor 0`) で回していた。この設定では hoop ソースは `A_planar` を使い、closure は使われていない (段 ② のダンプの `closure_active = 0`)。したがって 3 は、段 ② の closure の変更を通っていない。
    - そこで、`mesh.hoopAreaFromClosure: 1` を足した同じ一様な静止場で、3 と同じ判定をもう一度行う。float は旧・新・旧の再実行の 3 本、FP64 は新 1 本 (記録)。
    - 合格条件は 3 と同じ: float の新が旧の 1.1 倍以下 (1 step 目の `res_roUy/(P·A_planar)` の最大と RMS、10 step 後の最大の |u|)。
    - 3 の判定 (PASS) は既定の経路の非劣化として残し、closure の変更の確認は 3' で行う。

### 6.4 事前登録: 段 ③ の合格条件 (2026-10-10、段 ③ の実装の前・段 ② の結果を見る前)

- **対象**: 毎 step の座標の差 `cc1 − cc0` を e32 (`ge_x..ge_z`) に置き換える。
  - 粘性の内部面 (`viscousFlux_d.cu`)
  - スカラー拡散のすべての量 (k/ω・トレーサ・遷移・凝縮モーメント)
  - 化学種の拡散
  - 受動種の拡散と FCT
  - 陰解法の粘性の対角: スカラー DPLUR (`timeIntegration_d.cu:587`)、block DPLUR (`:952`、e32 を `ST` にする)、前処理版 (`:1327`)
  - 式・型・ガード・対象面・向きは §4.2a のとおり変えない。
- **FP64 で値が変わる経路**: block DPLUR は `implicitSolvePrecision 0` (既定) だと FP64 のビルドでも `ST = float` になる。旧は float に落とした座標の差、新は fl(e64) なので値が変わる。
  - 左辺だけの変化で、不動点は変わらない (§4.2a の例外)。
  - そのため、FP64 の同一性は `implicitSolvePrecision 1` (ST = double) で判定し、既定の 0 は記録だけにする。
- **合格条件**:
  1. float と FP64 のビルドが通る。
  2. **FP64 の面流束の同一性**: FP64 のビルドで `FORGE_DIAG_GEOMAB_DUMP` を case/45 の B0 の状態 (`run_0252` の res_60000、§6.2 の V0 と同じ状態) で回す。旧腕 (座標の差) と新腕 (e) の面の流束 (粘性 6 成分・k・ω) がビット単位で一致する。
  3. **FP64 の同一性 (V2 の実施形)**: case/45 の B0 を `implicitSolvePrecision 1` で 20 step 回す。新 (段 ③) 対 旧 (段 ②) の場の最大の相対差が、各量で旧どうしの再実行の差の 10 倍以下。
     - V2 の補助条件「20 step 後の状態 ≤ 1e-10」は判定に使わない。旧どうしの再実行の差がすでに 9e-8〜5e-7 ある (§6.2、原子加算の順序) ためで、段 ① と段 ② の実施形 (再実行の差との比) に合わせる。この変更は段 ③ の結果を見る前に決めた。
  4. **float の既定の経路**: float の case/56 で 20 step。新旧の差が、再実行の差の 10 倍以下。
     - case/56 の格子の HDF5 は float32 なので、cc64 は float の値そのものになる。e32 = fl(cc64[ic1] − cc64[ic0]) は float の引き算 (正しく丸められる) と同じ値なので、一致するはずである。
  5. **記録** (合否に使わない):
     - FP64 の既定 (`implicitSolvePrecision 0`) の 20 step の新旧の差。
     - float の case/45 B0 の 20 step の新旧の差 (改善の向きの変化)。
     - どちらも NaN/Inf が無いことは必須。
- 化学種・受動種・遷移・凝縮モーメントの経路には、case/45・56 で通らないものがある。実行の回帰は V6 (case/09・44・16・48) で見る。§6 V6 の注記のとおり、係数の検査だけで済ませたことにしない。

### 6.5 段 ② の結果 (2026-10-10)

- **実装** (implementer): 読み込み時の量と接続を double の写しから作る。
  - LSQ の係数: gradLSQ 1/2 とも、変位は e32 に向きを付けたもの (`calcGradient_d.cu` の `lsq_disp`)
  - 周期の相手と、LSQ の継ぎ目の同値類 (`mesh.cpp` の `periodicNearest<double>`)
  - ラインの接続 (`implicitLineTopology`、節点座標は `coord64`)
  - node の壁関数・WMLES の境界面ごとの (irep, y): `cuda_forge/wallRepPoint_d.cu` で起動時に 1 回作る。毎 step の選び直しをやめた
  - d1/d2 (`conjugateWall.cpp` の `buildFirstInterior`)、`delta_les`
  - 軸対称の closure: デバイスの最終の面ベクトルを double で足して、最後に 1 回丸める
  - 診断 `FORGE_DIAG_GEOM_STAGE2_DUMP` (`cuda_forge/geomStage2Dump.cpp`)
- **ビルド**: float sha `25aed33a…`、FP64 sha `1b8590e8…` (AWS `~/forge-fgeom2-{f32,fp64}`)。
- **run**: case/45 `run_0396〜0408`、case/56 `run_0082〜0084`、case/39 `run_0956・0957`、case/48 `run_0046・0047` (各 case の README)。
  - 台本は `case/45.isobutane_m6_d155/fg2.sh` → `fg2b.sh` → `fg2c.sh` (準備の不備で 2 回回し直した)、`fg2d.sh` (3')。集計は `fg2_an.py`・`fg2d_an.py`。

| 条件 | 結果 | 判定 |
| --- | --- | --- |
| 1 ラインの接続 case/45・56 | float と FP64 で一致 (旧 float も一致) | PASS |
| 1 周期の相手・継ぎ目の同値類 case/39 | 一致 (1581 + 4131 面、5559 group) | PASS |
| 1 壁の代表点 irep case/48 | 一致 (1001 面) | PASS |
| 2 cInt 99.9 % 点 ≤ 1e-5 | case/45 7.3e-6、56 5.7e-8、39 5.8e-8、**48 0.0722** | **case/48 FAIL** |
| 2 壁の y ≤ 1e-6 | 最大 2.2e-7 (case/45。旧 float 0.279) | PASS |
| 2 d1/d2 ≤ 1e-6 | 4 ケースとも FP64 と完全に一致 (旧 float は case/45 で 0.279) | PASS |
| 3 closure の非劣化 (既定の経路) | float 新 / 旧 0.95〜1.00、FP64 の新旧の差は再実行の差の 0.89 倍 | PASS (closure は通らない) |
| 3' closure の非劣化 (`hoopAreaFromClosure 1`) | float 新 / 旧 最大 0.91・RMS 0.87、\|u\| 1.01 | PASS |
| 4 FP64 の同一性 | 新旧の差 / 再実行の差の最大 1.30 | PASS |
| 5 float の既定の経路 case/56 | 同 1.25 | PASS |

- **case/48 の cInt の FAIL の中身**: 直交格子で真値が 0 の成分 (x 方向の隣への係数の y 成分など) が、x・y とも約 50 % ある。FP64 にもその成分に丸めの雑音があり (0 ちょうどは 1 % 未満)、係数ごとの相対差の分母に入っている。
  - その成分を除くと、99.9 % 点・最大とも 8.4e-8。
  - incidence のベクトルの相対差 ‖c − c64‖/‖c64‖ も最大 8.6e-8 (旧 float 3.5e-5)。
  - 事前登録の判定は FAIL のまま記録する。改めた基準 (4a) の判定は別に書く (下の解釈)。
- **段 ② の外の観察** (静止場の run から):
  - case/45 の既定 (`hoopAreaFromClosure 0`・`pRef 0`) では、一様な圧力で、軸の近く (j 2〜4) に `3.3e-3·P·A_planar` の偽の半径方向の力が出る。FP64 でも float でも同じ値で、10 step で偽の速度 0.257 m/s になる。
  - `hoopAreaFromClosure 1` にすると、FP64 の残差は 8.7e-10 まで、偽の速度は 6.4e-4 m/s まで下がる。
  - B0 の残差の床の最大も j 3・x/r_t −1.8 にある (line-implicit-speed §6.18)。結び付けるかどうかは諮問中。
  - float の壁際 (j 118・119) の静止場の残差は 0.36〜0.69 で、段 ② の前後でも closure の有無でもほとんど変わらない (FP64 は 7e-5)。
- **解釈**: codex (diagnose) に諮った: `notes/reviews/2026-10-10-float-geometry-stage2-result-diagnose.md` — 段 ② の無条件の合格は保留。旧基準の FAIL を残し、改めた基準で判定し直してから段 ③ へ。次の計算は B0 の保存場から `hoopAreaFromClosure` 0/1 の 2000 step の短い比較 (line-implicit-speed §6.23)。

| 重大度 | 指摘 | 採否 | 対応 |
| --- | --- | --- | --- |
| Major | case/48 の FAIL は物差しの問題 (成分の相対差の分母が真値 0 の成分を拾う) | 採用 | 元の 0.0722・FAIL を残し、改めた基準 (4a) の判定を別に書く。実装全体に欠陥が無いことの証明とはしない |
| Major | ベクトルの相対差に変えるなら条件を足す | 採用 | 4a |
| Major | 3・3'・4・5 は限定的に有効。closure を使う経路の FP64 の新旧比較が無い。closure を使わないと配列を作らないので、ダンプの「float と FP64 の差 0」は精度の証拠にならない | 採用 | 4b。§1.2 の「A_closure の差 0」の記録は証拠から外す |
| Major | 判定の台本の欠陥: `neq` が欠損を `None` にし、呼び出し側が 0 (一致) と数える。本番の配列との一致の属性が表示だけ。3' の非有限が合否に入らない | 採用 | 4c。保存済みの成果物で集計し直す |
| Major | H3 (float の壁際の静止場の残差は原子加算が原因で段 ③・④ では減らない) は未分離 | 採用 (要再検証) | V3 で面の流束と集めた後の残差を分けて測る (§5.1 #7 に追記) |
| Major | H4 (軸の近くの偽の力が B0 の床を作る) は未検証。静止場は半径の運動量、床は密度の残差で、最大の位置も違う (x/r_t −3.02 と −1.79)。closure の切り替えは周方向の粘性応力に掛ける面積と軸方向の補正も変える | 採用 (要再検証) | line-implicit-speed §6.23 に、同じ保存場からの 0/1 の比較を事前登録して回す |
| — | 段 ③ の V2 を 20 step 後の状態の比較だけで代えない。同じ入力での本番の面流束・commit 前の残差・dq の照合を残し、isp 0 も係数・ガード・対象面の変更を検査する | 採用 | 4d |

- **改めた基準と追加の確認** (2026-10-10、段 ② の結果を見た後の再評価。以後の段の係数の比較の基準として固定する):
  - **4a LSQ の係数**: 4 ケースそれぞれで次をすべて満たす。
    - 非ゼロの incidence のベクトルの相対差 ‖c32 − c64‖/‖c64‖ の最大 ≤ 1e-6。
    - 全係数が有限。
    - c64 がゼロベクトルの incidence (境界の面 `ip ≥ nNormalPlanes` を含む) では、c32 もゼロベクトル。
    - 線形の再現の行列 M_i = Σ_j c_ij ⊗ d_ij (d は同じ倍精度の変位を両方に使う) の差 ‖M32 − M64‖_F / ‖M64‖_F の最大 ≤ 1e-6 (各節点。z を切り捨てた疑似 2D では z の行・列は両方 0)。
  - **4b closure を使う経路**: case/45 に `hoopAreaFromClosure: 1` を足してダンプ (float・FP64) し、デバイスの `A_closure` を、デバイスの面ベクトルの double の和 (`closure_sum64`) を 1 回丸めた値と照合する (float・FP64 とも不一致 0)。さらに、静止場 (3') の FP64 の旧・旧の再実行を足し、FP64 の新旧の差 (1 step の `res_roUy`) が再実行の差の 10 倍以下であること。
  - **4c 判定の台本**: 欠損・形の不一致・非有限・本番の配列との不一致 (`lines_device_equal_new ≠ 1`、`lsq_production_vs_rebuild_mismatch ≠ 0`) は判定不能 (= FAIL) にする。3' の非有限も合否に入れる。保存済みの成果物で集計し直す。
  - **4d 段 ③ の V2 の追加** (§6.4 に足す): FP64 (`implicitSolvePrecision 1`) で、同じ入力 (B0 の保存場) から 1 step、旧・新・旧の再実行。commit 前の残差 (`res_ro`・`res_roUx`・`res_roUy`・`res_roe`・`res_roK`・`res_roOmega`) と dq = Q₁ − Q₀ の各量の最大の差が、再実行の差の 10 倍以下。isp 0 は、implementer の報告と差分の読みで、変えたのが `ST` に入れる値だけ (係数の式・ガード・対象面は不変) であることを確かめて記録する。

- **やり直しの判定の結果** (2026-10-10、`case/45.isobutane_m6_d155/fg2_rejudge.py` → `fg2_rejudge.json`、4b の run は `fg3.sh` の `run_0415〜0418`): **すべて PASS**。
  - 4c: 4 ケース × 2 ビルドで、ライン (lineImplicit のとき)・LSQ の本番の配列 = ダンプの作り方、double の写しあり。欠損・非有限を判定不能にしたうえで 1・2・3・3'・4・5 も PASS (値は上の表と同じ)。
  - 4a: ベクトルの相対差の最大 1.2e-7 (case/45。旧 float 0.34)・1.0e-7 (56)・5.9e-8 (39)・8.6e-8 (48。旧 float 3.5e-5)。ゼロベクトルの incidence (境界の面、9680・1494・16446・2285 本) は float も 0。線形の再現の行列の差の最大 7.9e-8・9.5e-8・4.6e-8・7.4e-8 (旧 float 0.24・9.5e-8・4.6e-8・1.2e-5)。元の基準の case/48 の FAIL (0.0722) は記録として残す。
  - 4b: closure を使う経路で、デバイスの `A_closure` は Σ±S の double の和を 1 回丸めた値と全節点で一致 (float・FP64)。FP64 の `A_closure` は、生の double の面ベクトルと面の重心の半径から作った closure とも差 0。静止場 (closure あり) の FP64 の新旧の差は 5.7e-10 で、再実行の差 5.1e-10 と同じ。
  - → **段 ② は合格** (改めた基準。§5.1 #4 を完了にし、段 ③ へ)。

### 6.7 段 ③ の結果 (2026-10-10)

- **実装** (implementer): 本番の経路で、毎 step の座標の差 `cc1 − cc0` を e32 (`ge_x..ge_z`) にした。
  - 粘性の内部面 (`viscousFlux_d.cu`)
  - スカラー拡散 (融合版と、トレーサの単独版。k/ω・遷移・凝縮モーメント・化学種のスカラーの全経路)
  - 化学種の拡散 (`species_diffusion_d`)、二相拡散と監査、tpfd 診断
  - 受動種の拡散と FCT (`passiveKernels_d.cuh`・`passiveFct_d.cuh`)
  - 3 種の DPLUR の粘性の対角 (`timeIntegration_d.cu`): スカラー版・block 版 (`ST` に e を入れる)・前処理版。`(ic0 == ic) ? e : −e` で向きを合わせる。
  - V0 の診断の旧腕 (ge が nullptr のときの座標の差) は残した。
- **ビルド**: float sha `d4991d4b…`、FP64 sha `129de3f4…` (AWS `~/forge-fgeom3-{f32,fp64}`)。
- **run**: case/45 `run_0419〜0428`、case/56 `run_0085・0086`。台本 `case/45.isobutane_m6_d155/fg4.sh`、集計 `fg4_an.py` (`fg4_verdict.json`)。旧は段 ② のバイナリ。

| 条件 (§6.4・§6.5 の 4d) | 結果 | 判定 |
| --- | --- | --- |
| 1 ビルド | float・FP64 とも通った | PASS |
| 2 FP64 の面流束の同一性 (GEOMAB_DUMP、`run_0419`) | 旧腕と新腕の面の流束のビット不一致 0 (粘性 6,881,028・スカラー 2,293,676 値)。非有限は境界面 9680 面の未評価の印だけで、両腕で同じ | PASS |
| 3 FP64 isp 1 の 20 step | 新旧の差 / 再実行の差の最大 1.77 | PASS |
| 4 float case/56 の 20 step | 同 2.81 | PASS |
| 4d FP64 isp 1 の 1 step (commit 前の残差・dq) | 同 1.50 | PASS |

- **記録** (合否に使わない):
  - FP64 の isp 0 (既定) の 20 step の新−旧: 場の最大の相対差 3e-6〜2.5e-5 (再実行の差 2e-9〜5e-7)。block DPLUR の左辺に入る値が fl32(cc1) − fl32(cc0) から fl32(e64) に変わったため (§4.2a の例外)。
  - float の case/45 B0 の 20 step の新−旧: 場の最大の相対差 1e-2 (ro・roUx・roe・roOmega)、P 1.5e-4。20 step 後の `rms_roK` 7.6e-4 → 2.7e-4、`rms_roOmega` 6.6e3 → 1.5e3 (旧 → 新)。ほかの rms は同じ桁。
  - どの run も非有限なし。
- **isp 0 の差分の読み** (§6.5 の 4d): 変えたのは `ST` に入れる値 (`static_cast<ST>(ccx[other]) − static_cast<ST>(ccx[ic])` → `static_cast<ST>(±ge_x[ip])`) だけで、係数の式・ガード・対象面 (`!(isNode && !has_nbr)` の除外) は変えていない (`timeIntegration_d.cu` の段 ③ の差分、implementer の報告と一致)。
- **段 ③ は合格**。解釈 (float の k・ω の残差の低下が到達に効くか) は V3・V4 で見る。

### 6.6 事前登録: 段 ④ の合格条件 (2026-10-10、段 ④ の実装の前・段 ③ の結果を見る前)

- **対象**: 再構成とリミタが毎 step 作る座標の差 (段 ③ の implementer の目録、2026-10-10)。
  - 再構成: `convection/convectiveFlux_{roe,hlle,slau,keep}_d.inc.cuh`、境界の `convectiveFlux_boundary_d.inc.cuh` (pc − cc)、`AUSM_d.cu`、`convection/legacy/*.inc.cuh`
  - リミタ: `limiter_d.cu`、`limiterPeriodic_d.cuh`、`passiveLimiter_d.cuh`
  - 辺中点の ±0.5e は `g_reconEdgeMid == 1 && ip < nNormalPlanes` の内部辺に限る (§4.2a)。それ以外 (境界面・cell の双対面重心) は、面の両側の pc − cc を double の写しで引いて 1 回だけ丸めた配列 (r0・r1) を新しく作って読む。
  - 式・型・ガード・対象面・向きは変えない。
- **合格条件**:
  1. float と FP64 のビルドが通る。
  2. **FP64 での値の同一性**: FP64 では 0.5·e と r0・r1 は今の差とビット単位で同じになるはず。段 ① と同じく、FP64 のビルドで新しい配列 (r0・r1) と今の差 `pc − cc` の不一致が全面で 0 (起動時の確認を診断の環境変数で出す)。
  3. **FP64 の同一性 (V2 の実施形)**: §6.4 の 3 と同じ (case/45 の B0、`implicitSolvePrecision 1`、20 step、新 = 段 ④・旧 = 段 ③・旧の再実行、各量で ≤ 10 倍)。§6.5 の 4d と同じ 1 step の commit 前の残差と dq の照合も行う。
  4. **float の既定の経路**: case/56 で 20 step、新旧の差が再実行の差の 10 倍以下。case/56 の格子は float32 なので、0.5·e と r0・r1 は float の引き算と同じ値になるはず。
  5. 記録: float の case/45 の B0 の 20 step の新旧の差。どちらも NaN/Inf が無いこと。
- 段 ④ の後は V3 (同じ状態の残差、float 対 FP64)、V4・V5・V6 へ進む。V3 では、面の流束と集めた後の残差を分けて測る (§6.5 の H3 の指摘)。

### 6.8 事前登録: V3 の実施形 (2026-10-10、段 ④ の結果を見る前)

- **状態**: B0 の状態 (`run_0252` の res_60000、FP64) の保存量・組成・乱流量を float に丸めた Q32 を、double のまま `nozzle.h5` に書く。float のビルドも FP64 のビルドも同じ値を読む (FP64 側は double(Q32))。
- **腕** (case/45、B0 の構成、1 step、`output.level 2`):
  - 旧 float: 段 ① の float のバイナリ (`a874b748…`、本番の経路は幾何の変更の前と同じ) を 2 回
  - 新 float: 段 ④ の float のバイナリを 2 回
  - 参照: 段 ④ の FP64 のバイナリを 2 回 (`implicitSolvePrecision 1`。残差は isp に依らないが、そろえる)
- **比べる量**: 1 step 目の出力の commit 前の残差 `res_ro`・`res_roUx`・`res_roUy`・`res_roe`・`res_roK`・`res_roOmega` (Q32 から計算したもの)。
- **ノルム**: 量ごとに E = ‖R32 − R64‖₂ / ‖R64‖₂。領域は全域と壁際 (j ≥ 117、case/45 は j = 120 が壁) の 2 つ。量を 1 つのノルムに混ぜない。
- **再実行の差** N: float の 2 本の差と FP64 の 2 本の差の大きい方を ‖R64‖₂ で割ったもの (腕・量・領域ごと)。
- **合格条件**:
  - 悪化しない: すべての量・領域で E_new ≤ 1.1·E_old + 10·N。
  - 改善: 熱 (`res_roe`)・運動量 (`res_roUx`・`res_roUy`)・乱流 (`res_roK`・`res_roOmega`) のうち、E_old > 10·N の量・領域では E_new ≤ 0.5·E_old。
- **記録**: 量・領域ごとの E_old・E_new・N の表。最大の誤差の位置 (j と x/r_t)。float の静止場の壁際の残差 (§6.5 の観察 B・D) と同じ場所かどうか。面の流束の誤差は粘性・k/ω について V0 の経路 (GEOMAB の DUMP/REF) で測れるので、集めた後の残差の誤差のうち面の流束で説明できない分を、`res_roUx`・`res_roUy` の差として分けて書く (§6.5 の H3)。

### 6.9 段 ④ の結果 (2026-10-10)

- **実装** (implementer): 再構成 (Roe・HLLE・SLAU・KEEP・境界) とリミタ (r1・fused5・周期の 2 段・G1 の検査・受動種) が、e と面の配列 r0・r1 を読む。
  - r0・r1 は `gr0_*`・`gr1_*` = fl(planeCent64[ip] − cc64[ic0/ic1])。全面が対象で、`variables.cpp` の `fillGeomDiffR` が作る。
  - リミタが ic1 側から見る辺中点は `0.5·(0 − e)` (符号付きゼロも従来と同じにするため)。
  - 起動時の照合 `FORGE_DIAG_GEOM_STAGE4_CHECK=1` を足した。
  - 本番で起動されない legacy (`AUSMp_d` など) と、ビルド対象外の `AUSM_d.cu` は触っていない。
- **メモリ**: r0・r1 で面あたり +24 B (float)。2D の node で約 48 B/節点、e と合わせて約 72 B/節点。デバイスとホストの両方に持つ。**§4.4 の見積もり (+25〜50 B/節点) を超える**。減らすなら、使う経路だけ確保するか、ホストの写しを捨てる (未決定、§5.1 #12)。
- **ビルド**: float sha `eac9dc38…`、FP64 sha `5b4b0169…` (AWS `~/forge-fgeom4-{f32,fp64}`)。単体試験 4 本 (`test_periodic_limiter`・`test_passive_scalar`・`test_passive_fct`・`test_farfield_flux`) は PASS。
- **run**: case/45 `run_0429〜0435`、case/39 `run_0958`、case/56 `run_0087〜0089`。台本 `case/45.isobutane_m6_d155/fg5.sh`、集計 `fg5_an.py`。旧は段 ③ のバイナリ。

| 条件 (§6.6) | 結果 | 判定 |
| --- | --- | --- |
| 1 ビルド | float・FP64 とも通った | PASS |
| 2 FP64 での値の同一性 (起動時の照合) | case/45・39 (周期)・56 の全面で、e・ic1 側の e・r0・r1・境界の r0 のビット不一致 0 | PASS |
| 3 FP64 isp 1 の 20 step | 新旧の差 / 再実行の差の最大 1.38 | PASS |
| 3 (§6.5 の 4d と同じ形) 1 step の commit 前の残差と dq | 同 1.24 | PASS |
| 4 float case/56 の 20 step | 同 1.37。照合では内部面の e・r0・r1 の不一致 0 (境界面の 83〜85 面はゴーストの中心を double で鏡映するための差で、段 ④ の経路では使わない) | PASS |
| 5 記録: float case/45 B0 の 20 step | 新−旧 の場の最大の相対差 4e-6〜1.3e-4、rms はほぼ同じ。照合では内部面の 99 % 以上で e・r0・r1 が旧の float の差と違う (最大の相対差 r1 で 5.0) | — |

- **段 ④ は合格**。

### 6.10 V3 の結果 (2026-10-10、§6.8 の事前登録)

- run: `run_0440・0441` (旧 float = 段 ①)、`run_0442・0443` (新 float = 段 ④)、`run_0444・0445` (FP64 の参照 = 段 ④)、`run_0446_v3_dump32`・`run_0447_v3_ref64dump` (粘性・k/ω の面の流束)。台本 `fg6.sh` (`q32.py` で Q32 を作る、8 本とも同じ sha `66ecd4ba…`)、集計 `v3_an.py` (`v3_verdict.json`)。
- **判定: PASS** (全量・全領域で悪化なし、改善を要求した 10 件すべてで E_new ≤ 0.5·E_old)。

| 量 | 全域 E_old → E_new | 壁際 (j ≥ 117) E_old → E_new | N (全域 / 壁際) |
| --- | --- | --- | --- |
| `res_ro` | 0.39 → 0.028 | 12.0 → 0.33 | 8e-4 / 4e-7 |
| `res_roUx` | 0.60 → 0.074 | 210 → 1.20 | 1.1e-2 / 0.28 |
| `res_roUy` | 0.48 → 0.14 | 19.9 → 1.42 | 2.6e-2 / 0.26 |
| `res_roe` | 0.45 → 0.044 | 680 → 0.42 | 1.3e-3 / 1.9e-4 |
| `res_roK` | 19.9 → 0.028 | 1420 → 0.018 | 4.6e-3 / 5.8e-4 |
| `res_roOmega` | 9.8e4 → 0.37 | 9.8e4 → 0.37 | 4.5e-3 / 4.5e-3 |

- E は ‖R32 − R64‖₂/‖R64‖₂。B0 の状態は床の近くで ‖R64‖ が小さいので、E が 1 前後になる量・領域がある。
- **残った誤差の内訳 (記録)**: 粘性の面の流束の誤差を節点に集めた大きさは ‖R64‖ の 5e-6〜2e-3 で、残差全体の誤差 E_new (0.04〜1.4) のごく一部。~~残りは粘性の面の流束以外にある~~ **撤回 (2026-10-10、codex diagnose)**: この比較の FP64 の REF は、原始量・勾配・物性・面ベクトルを float 側の入力で上書きするので、測ったのは「同じ入力での演算の差」だけで、入力の丸めの差は消えている。k・ω の拡散も集計していなかった。入力の丸めを候補に戻し、§6.12 の A/B で確かめる。
- 最大の誤差の位置 (新): `res_roOmega` は j 119・x/r_t −0.91 (壁の 1 つ内側、スロートの少し上流)。運動量・エネルギーは入口の近く (x/r_t −11〜−8) の軸寄り〜壁際。
- **解釈**: codex (diagnose) に諮った: `notes/reviews/2026-10-10-float-geometry-v3-v4-diagnose.md` — V4 の前に、同じ Q32 から FP64 の流束の直前の入力を「自分で作った値 / float 側から移した値」で切り替える 1 step の A/B (§6.12) を行い、粘性・ω の拡散を除外した解釈を直す。
  - 段 ①〜④ で、座標の差に由来する誤差は桁で減った (相対の改善の判定は採用)。残りの発生箇所は未確定。
  - `res_roOmega` の E_new 0.37 は、固定した Q32 での作用素の差であって、反復の後の残差の床でもその下限でもない。床の値は V4 の末尾の系列で見る。

| 重大度 | 指摘 | 採否 | 対応 |
| --- | --- | --- | --- |
| Major | 「残りは粘性の面の流束以外」は REF が入力を上書きするので言えない。k・ω の拡散も集計していない | 採用 | 上の記録を撤回。§6.12 の A/B |
| Major | V3 の集計器が、参照のノルム 0 を除外・面の流束の非有限を 0 に置換・内訳の集計の失敗でも総合 PASS (評価件数 0 で PASS を合成データで再現) | 採用 | `v3_an.py` を直して集計し直す (既知の未評価の境界面だけを除き、必要な量・領域の欠損は判定不能。残差の改善の判定と内訳を分ける) |
| Major | E_new 0.37 から ω の床を予測しない | 採用 | V4 の末尾の系列で見る |
| Major | V4 を「水準への到達 + 到達時のフラグ」に直すのは目的を限って採用。発散・不完全なデータは許さない | 採用 | §6.11 に事前登録し直す |
| Minor | 「FP64 の B0 自身が旧 V4 を満たせない」は根拠不足 (B0 の最初の到達のドリフトは +0.007/−0.004/+0.009 %) | 採用 (呼び出し側の誤り) | 変更の理由は「不可能だから」でなく「評価の目的を定常の同等性から停止規則のもとでの比較に変える」と書く (§6.11) |
| Minor | r0・r1 のメモリは今は受け入れ、見積もりを直してピークの使用量と V5 を測る | 採用 | §4.4 を訂正、§5.1 #12 |


### 6.11 事前登録: V4 の実施形 (2026-10-10、§6 の V4 を置き換える。V4 を回す前)

- **目的**: 登録した停止規則 (水準) のもとで、float の新しい経路が FP64 と同じように止まるかを比べる。定常解の同等性や、物理の精度 0.05 % を保証する判定ではない。
- **腕**: 段 ④ と同じコードの系列 (commit `56b1a6fe`、段 ④ + 既定 off の commit の診断) の float と FP64。同じ格子・同じ出発点 (`run_0183` の res_100000)・B0 の実効設定 (cfl 4・方向別 Δτ・上限 50・`implicitThermalJacobian 5`・緩和 0.7・ライン陰解法)。2 本を同時に投入する (GPU は共有)。両方に `FORGE_DIAG_COMMIT_LOSS=2500` を付ける。
- **停止**: 2500 step ごとに出力し、水準 (|欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step) を 2 出力連続で満たした **2 つ目**の出力で止める。上限は両腕とも 20 万 step。参照との差や準定常の判定を理由に延長しない。
- **比較が成り立つ条件** (満たさなければ判定不能):
  - 到達を保存した系列から計算し直して確かめる。出力の間隔が 2500 step でそろい、末尾 2 万 step に 9 点以上ある。
  - 必要な残差の列 (`rms_ro`・`rms_roUx`・`rms_roUy`・`rms_roe`・`rms_roK`・`rms_roOmega`・`rms_roY*`) がそろい、系列と残差に非有限がなく、発散していない。
  - 見張りは末尾 9 出力の場を消さずに残す (今の `m9_watch.py` は最新 5 つと到達だけを残し、出力の間隔と点数を確かめていないので、それを直してから使う)。
- **判定**:
  - float が上限までに水準に入らなければ FAIL (到達しない)。
  - 各腕の到達時の値で比べる。分母は FP64 の到達時の値の絶対値。
    - θ_r(40/70/94) の各断面の差 ≤ 0.05 %
    - 符号つきの Q_w の差 ≤ 0.1 %
    - 末尾 2 万 step の各 `rms_*` の中央値 ≤ FP64 の 2 倍 (FP64 の中央値が 0 の列は、float も 0 でなければフラグ)
    - 到達時の 2π·Σ|res_ro| ≤ FP64 の 2 倍
  - 超えたものはフラグ (改善・悪化は断定しない)。
  - **合格の名前**: 「登録した水準へ到達し、停止時の比較でフラグなし」。
- **記録**: 到達の step 数、`check_convergence` と `check_quasisteady` (θ_r は `--drift 0.0001 --osc 0.0001`、Q_w は `--drift 0.0002 --osc 0.0002`、`--series-csv --tail 1.0`) の VERDICT はそのまま併記する (NOT CONVERGED・非定常は許すが、発散と不完全なデータは許さない)。commit の診断の統計、壁の熱流束の分布の差、共有 GPU での実際の経過時間 (参考。採用に使う総時間は V5 の専有の単価・出力・起動の費用から出す)。

### 6.12 事前登録: V3 の残差の差のうち、流束の直前の入力の丸めで説明できる分 (2026-10-10、codex diagnose、実行前)

- **問い**: V3 の ω の残差の差 ΔR = R32 − R64 (壁際 j ≥ 117) のうち、拡散の面の流束に入る入力 (原始量・勾配・物性・F1・幾何) の float の丸めで説明できる分はどれだけか。
- **組み方** (同じ Q32、段 ④ の FP64 のバイナリ、各 1 step・各 2 回):
  - A: FP64 で `FORGE_DIAG_GEOMAB_DUMP` (流束の直前の入力を自分で作る)
  - B: FP64 で `FORGE_DIAG_GEOMAB_REF` (float の DUMP `run_0446` の入力を移す。既存の `run_0447` に再実行を 1 本足す)
  - 入力の差の寄与 C = Σ±(F_B − F_A) を、本番と同じ符号 (ic0 に +F、ic1 に −F) で節点に集める。ω の拡散を必ず含める (粘性も記録する)。
- **判定** (ω、壁際 j ≥ 117): r = ‖ΔR − C‖₂ / ‖ΔR‖₂。ΔR は V3 の `run_0442` − `run_0444`。必要な内部面の欠損・非有限は判定不能。
  - **支持** (入力の丸めが過半を説明): r ≤ 0.5、かつ再実行の差 (float・FP64 の各 2 本、A・B の各 2 本) が ‖ΔR‖ の 10 % 以下 → 「粘性・拡散以外」の除外を撤回する。
  - **棄却**: r ≥ 0.9 (同じ再現性の条件) → 入力の丸めが過半を説明するという仮説を棄却する。演算・足し込みの丸めが残るが、足し込みが原因とまでは決めない。
  - それ以外は判別不能。どれも長期の残差の床の判定にはしない。
- **記録**: k についても同じ r。float の演算の差まで含めた C_all = Σ±(F_32 − F_A) (F_32 は `run_0446` の新腕) での r。

### 6.13 §6.12 の結果 (2026-10-10)

- run: A = `run_0448_v3_self64a`・`run_0449_v3_self64b` (FP64 の GEOMAB_DUMP、入力を自分で作る)、B = `run_0447_v3_ref64dump`・`run_0450_v3_ref64dump2` (FP64 の GEOMAB_REF、float の入力を移す)。段 ④ の FP64 (`5b4b0169…`)、Q32 は V3 と同じ (sha `66ecd4ba…`)。台本 `fg7.sh`、集計 `fg7_an.py` (`fg7_verdict.json`)。
- **判定 (主判定 ω、壁際 j ≥ 117): 棄却**。r = ‖ΔR − C‖/‖ΔR‖ = 1.010 (‖C‖/‖ΔR‖ = 0.44 だが ΔR と向きがそろわない)。再実行の差 / ‖ΔR‖ = 1.2e-2 (≤ 0.1 で判定できる)。
  - → 拡散の面の流束に入る入力の丸めは、ω の残差の差の過半を説明しない。
- **記録**:
  - float の演算の差まで含めた C_all でも r_all = 1.004 (ω)。拡散の面の流束全体の差 (入力 + 演算) でも説明できない。
  - k (壁際): r = 0.83 (判別不能)、r_all = 0.69、‖C‖/‖ΔR‖ = 0.57。全域では r = 0.94。
- **解釈 (未確定、諮っていない)**: ω の残差の差は拡散以外 (ソースの大きな生成と散逸の釣り合い、対流、足し込み) にある。j 119 の ω は 6〜8e10 1/s で、大きな項どうしの釣り合いで残差が決まる (line-implicit-speed §6.22)。どれかは決めていない。V4 の末尾の系列で、float の ω の床が実際にどうなるかを見る。
- 合わせて、V3 の集計器を codex の指摘どおり直して集計し直した (`v3_an.py`、参照のノルム 0 は判定不能、内部面の非有限は判定不能、内訳と総合を分ける)。結果は変わらず PASS。

### 6.14 事前登録: §4.6 (変換器の幾何の精度) の合格条件 (2026-10-10、実装の前)

- 対象の格子: case/45 の `nozzle.msh` (node、軸対称)、case/48 の `fp_y1_3um.msh` (node、平面)、case/39 の丘 (3D、周期) の 3 つ。どれも、それぞれの run の設定 (`solverConfig.yaml`・`bcondConfig.yaml`) で変換する。
- **参照の固定**: 旧 FP64 の変換器は段 ④ の FP64 のビルド (`~/forge-fgeom5-fp64`、commit `56b1a6fe`、sha を記録) とする。新旧とも同じコンパイラ・浮動小数点のフラグ・依存ライブラリ (`HAVE_KDTREE` の有無を含む)・壁距離の探索の経路 (環境変数)・入力 (`.msh` と YAML の sha256) で変換する (codex diagnose 2026-10-10)。
- **比較器**: 必要なデータセットの欠損・形の不一致・非有限は FAIL。符号付きゼロを含むビット単位の比較。両側に同じ NaN があっても合格にしない。対象は `/MESH/*` (`RENUMBER_PERM` があればそれも)・`/PLANES/*`・`/CELLS/*`・`/DUAL/*`・`/BCONDS/*`・`/VIZMESH/*`・`/VALUE/wall_dist`。
- **1 FP64 の変換器は変わらない**: 新しい FP64 の変換器と旧 FP64 の変換器の出力で、既存のデータセットの型・形・値・属性が一致する。不一致が出たら先に変更箇所を特定する。根拠なく許容差に緩めない。
- **2 float の変換器が double の幾何を書く**: 新しい float の変換器の出力の幾何のデータセット (座標・primal と双対の幾何・境界の半割面・壁距離) が binary64 で、新しい FP64 の変換器の出力とビット単位の不一致 0 件。接続・番号・境界の所属・可視化の接続 (整数の表) も完全に一致する。float の初期値・境界条件の実数値までビルド間で同じにすることは求めない。
- **3 ソルバの警告**: float32 の格子 (例: case/56 の `run_0073` の `mesh.h5`) で起動すると警告が 1 行出る。新しい double の格子では出ない。座標だけ double で他が float32 の混在の格子 (試験用に作る) でも警告が出る。
- **4 ソルバが新しい格子を読む**: 新しい float の変換器で作った case/48 の格子をソルバ (float) が読み、double の写し (`coord64`・`cc64` など) に保持できる (段 ② のダンプ `FORGE_DIAG_GEOM_STAGE2_DUMP` で `has_geom64` と値を確かめる)。
- **5 互換性**: 旧の float32 の格子 (case/56) で新しい float のバイナリを 20 step (§6.4 の 4 と同じ形)。新旧の差 ≤ 再実行の差の 10 倍。これは読み込みの互換性の試験で、収束の証拠ではない。
- 判別の最初の段 (codex): case/48 の `fp_y1_3um.msh` と YAML を固定し、ビルドの型 (float / FP64) だけを変えて各 1 回変換する (ソルバは回さない)。1 件でも不一致があればビルドの型からの独立性を棄却し、最初に違う生成の段を追う (許容差で通さない)。
- 変換器の HDF5 はバイト単位では再現しない (同じバイナリでも md5 が違う。メモ converter-h5-not-byte-reproducible) ので、比べるのはデータセットの値。

### 6.15 §6.14 の結果 (2026-10-10)

- **実装** (implementer、案 (b)): 変換器は幾何の正本を `gmshReader::geo64` (`mesh/convGeom64.hpp`、double) と double の双対の作業配列に持ち、共用の `mesh` へは 3 か所 (節点の読み込み・makeMesh の最後・replacePrimalWithDual) で一方向に丸めて渡す。`writeInputH5` は幾何を正本から binary64 で書く。壁距離は `calcWallDistance64` (double) で計算し、`/VALUE/wall_dist` を double で書く (ソルバ用の写しは従来どおり丸める)。ソルバは幾何のデータセットのどれかが float32 なら起動時に警告を 1 行出す (`mesh.cpp` の `warnFloat32Geometry`)。比較器 `solver_density_cuda/tools/compare_mesh_h5.py`。手順は `procedures/calculation-workflow.md` に追記。
- **ビルド**: AWS `~/forge-fgeom6-{f32,fp64}`。旧 FP64 の変換器は `~/forge-fgeom5-fp64` (56b1a6fe)。台本 `case/45.isobutane_m6_d155/fg8.sh`、集計 `fg8_an.py` (`fg8_verdict.json`)。入力 (sha256 の先頭): case/45 `nozzle.msh` 9a891968…、case/48 `fp_y1_3um.msh` c60a7dae…、case/39 `hill_des_80x50x30.msh` d5b54520…。
- **結果: すべて PASS**。

| 条件 | 結果 |
| --- | --- |
| 1 新 FP64 対 旧 FP64 (全件) | case/45 101・case/48 113・case/39 148 データセットでビット不一致 0、属性も一致 |
| 2 新 float 対 新 FP64 (`--geometry-only --require-float64`) | 幾何 (binary64) と整数の表がビット不一致 0 (case/45 32・case/48 37・case/39 42 データセット) |
| 3 ソルバの警告 | float32 の格子 (case/56) と混在の格子 (`surfVect` だけ float32) で出る。新しい double の格子 (float・FP64) で出ない |
| 4 新しい格子を float のソルバが読む | 段 ② のダンプで has_geom64 = 1、FP64 のダンプとの比較で cInt のベクトルの相対差の最大 8.6e-8、irep の不一致 0、d1 の差 0 |
| 5 互換性 (case/56、float32 の格子、20 step) | 新旧の差 / 再実行の差の最大 1.67 |

- 記録: 変換器の最大常駐メモリは case/39 で +6.6 % (float)、+1.3 % (FP64)。float の変換器の出力ファイルは大きくなる (case/48 で 26.6 → 約 40 MB)。`fluent_h5_to_forge.py` は既定で幾何を float32 で書く (範囲外、ソルバが警告する)。case/48 の `gen_mesh.py` の `CONV_CFG` には現行の変換器が止まる廃止キーが 4 つある (直していない)。
- **§5.1 #13 は完了**。

### 6.16 V4 の結果 (2026-10-10、§6.11 の事前登録)

- run: `run_0451_v4_f32` (float)・`run_0452_v4_fp64` (FP64)。どちらも段 ④ + commit の診断 (`56b1a6fe`)、`run_0183` の res_100000 から B0 の構成、同時に投入 (GPU を共有)。台本 `case/45.isobutane_m6_d155/v4.sh`、判定 `v4_judge.py` (`_band_ab/cold_pair/v4_judge.json`)。
- 比較が成り立つ条件: 満たした (出力の間隔 2500・末尾 2 万 step に 9 点・非有限なし・発散なし。到達は系列から計算し直して見張りと一致)。
- **判定 (事前登録のとおり): 到達、フラグあり**。float は 140,000 step、FP64 は 125,000 step で水準に入った。

| 量 (各腕の到達時) | float | FP64 | 差 (分母 FP64) | 上限 |
| --- | --- | --- | --- | --- |
| θ_r(40) | 0.0793136 | 0.0778965 | +1.819 % | 0.05 % |
| θ_r(70) | 0.116610 | 0.114963 | +1.432 % | 0.05 % |
| θ_r(94) | 0.129414 | 0.127768 | +1.288 % | 0.05 % |
| Q_w | 9.42662e6 W | 9.29339e6 W | +1.434 % | 0.1 % |
| 2π·Σ\|res_ro\| | 9.158 | 2.340 | 3.91 倍 | 2 倍 |
| 末尾 2 万 step の rms の中央値の比 | ρ 3.45・ρu 1.15・ρv 1.91・ρE 2.58・ρk 3.58・ρω 1.38・ρY 3.45 | — | ρ・ρE・ρk・ρY がフラグ | 2 倍 |

- **記録**:
  - `check_convergence` は両方 CHECK FAILURES (NOT CONVERGED)、`check_quasisteady` (θ_r・Q_w の系列、閾値 1e-4・2e-4) は両方 NOT ALL STEADY。止まった値は定常解ではない。
  - commit の診断 (float、到達時): S_lostfrac = Σ|dq_req − 実際の更新|/Σ|dq_req| (丸めで切り上がった分も数えるので「消えた割合」ではない、codex 2026-10-10) は ρ 0.22 %・ρu 0.44 %・ρv 0.15 %・ρE 0.18 %・ρk 0.28 %・**ρω 4.9 %**・ρY 0.17 %。更新が丸ごと消えた節点は ρ で 56k/561k。FP64 は 1e-11 の水準。
  - **軌跡**: θ_r(70) の差は 2 万 step まで 0.00 %、2.5 万 −0.002 %、3 万 +0.23 %、3.5 万 +0.74 %、4 万 +1.13 %、5 万 +1.59 %、6 万 +1.68 % で、その後 1.45〜1.52 % で止まった。Q_w の差は 2.5 万 step から先に開いた (+0.12 %)。
  - **場の違いの場所** (到達時の場、float res_140000 対 FP64 res_125000): 収縮部 (x/r_t −15〜−1) の中間の半径で |U| の差が最大 45 % (|U| ≥ 1 m/s で正規化)、T 15 %、μt 最大 30 倍 (x/r_t −3.4・j 37)。喉より下流は |U|・T とも 1 % 未満。壁の 1 つ内側 (j 119) の T は収縮部で平均 +0.4〜0.8 K、下流で −0.1〜−0.2 K。
- **解釈 (未確定、諮る)**: float の新しい経路は、case/45 の B0 で FP64 と 1.3〜1.8 % 違う所で止まる。設計チェーンの許容差 (0.05〜0.1 %) を大きく超えるので、**今のままでは case/45 の設計の計算に float は使えない**。違いが 2.5 万〜4 万 step で急に開き、収縮部の流れの状態が違うことから、(i) float の丸め (残差・commit) で止まる点がずれた、(ii) 収縮部の遅いモードが別の状態に移った (どちらの精度でも止まりうる状態が 2 つある)、のどちらかと見ている。判別していない。

### 6.17 V5 の結果 (2026-10-10)

- run: `run_0453〜0459_v5_*` (D64 F32 D64 F32 D64 F32 D64、各 1000 step、専有の GPU、V4 と同じバイナリ `56b1a6fe`・診断 off)。台本 `case/45.isobutane_m6_d155/v5.sh`、記録 `v5.log`。
- **判定: PASS**。float 20.69・20.70・20.68 ms/step、FP64 32.42・32.31・33.04・32.32 ms/step。比 0.636 (基準 ≤ 0.70)。
- GPU のメモリの最大 (nvidia-smi の使用量): float 2254 MiB、FP64 2938 MiB (§5.1 #12 の記録)。
- 総時間の目安 (codex): V4 の到達の step で換算すると 0.636 × 140,000/125,000 = 0.71 (約 29 % 短い)。ただし V4 で float は FP64 と違う所で止まったので、**同じ品質に届く費用の比較ではない**。

### 6.18 V4 の解釈の諮問と、事前登録: 精度の切り替えの A/B (2026-10-10、実行前)

- codex (diagnose) に諮った: `notes/reviews/2026-10-10-float-v4-and-hoop-review-diagnose.md` — 同じ float の到達の場から「float で続ける / FP64 に切り替える」の対照つきの A/B を 1 つ行う。case/45 の設計の計算は当面 FP64 のままにする。

| 重大度 | 指摘 | 採否 | 対応 |
| --- | --- | --- | --- |
| Major | 相互の再開 (float の場を FP64 で・FP64 の場を float で) だけで「状態が 2 つ」を判定しない (起点も精度も違う 2 本では、再開の後の過渡と精度の効果を分けられない) | 採用 | 下の A/B (同じ起点、float の継続を対照) |
| Major | S_lostfrac から commit を主因としない (切り上げも数える) | 採用 | §6.16 の記述を訂正。領域別・符号つき・時系列で見るのは A/B の後 |
| Major | V4 の判定器: 残差の有限性を末尾の窓だけで見ている、判定ツールの結果が総合に入らない、Q_w が NaN でも REACH を返しうる | 採用 (要再検証) | 保存した成果物で、全期間の非有限・ツールの失敗・系列の欠損を確かめ直す (§5.1 #8) |
| — | `qAccumulatorFP64` の拒否だけを外さない (軸の基準の射影と `sstEnergyIncludesK` の増分にも対応が要る。ω は別の float の足し込み)。残差の double の足し込みを真因の確認なしに実装しない。`pRef` を純粋な丸めの A/B として扱わない (閉性の欠損があるので ΔpRef·(ΣW_y − A) が変わる) | 採用 | 手当の候補はこの A/B の後に決める |

- **第 1 仮説** (codex、確度 中): float に依存する作用素・更新の差が、収縮部の遅い過渡を偏らせ続けている (commit だけが主因とは絞れない)。対立仮説: 精度の変更への応答より、共通の長い過渡・停止規則の影響が大きい。状態が複数あって丸めがその選択を変えた可能性も残る (排他ではない)。
- **組み方**: 起点は `run_0451_v4_f32` の res_140000 (float の到達の場)。`cold_cfl.py prep … --field-from run_0451_v4_f32` (`restart_field.py`、同じ格子・ビット一致、FP64 側は float の値をそのまま広げる)。B0 の構成、V4 と同じバイナリ (`56b1a6fe`)。格子・hoop・`pRef`・CFL・リミタの基準は変えない。
  - A: float で続ける。A2: A の再実行 (再実行の差)。B: FP64 に切り替える。
  - 各 **30,000 step 固定** (水準で止めない)、2,500 step ごとに出力 (`res_ro` 入り)。3 本を同時に投入する。
- **判定** (θ_r(40/70/94) と符号つきの Q_w の 4 量。D = V4 の到達時の差 float − FP64 = +1.819・+1.432・+1.288・+1.434 %。動き Δ = 30,000 step 後 − 起点):
  - **支持** (精度に依存して差が保たれている): 4 量すべてで、B が FP64 の値の向きに 0.5·|D| 以上動き、|Δ_A| ≤ 0.1·|D|、かつ |Δ_B − Δ_A| が再実行の差 |Δ_A2 − Δ_A| の 10 倍を超える。
  - **棄却** (この期間で精度の変更は支配的でない): 4 量すべてで |Δ_B − Δ_A| ≤ 0.1·|D|、かつ A と B がともに 0.5·|D| 以上同じ向きに動く。
  - それ以外 (両方とも留まる場合を含む) は判別不能。短い期間で動かないことから「状態が 2 つ」とは言わない。
- **記録**: 固定した収縮部の領域 (x/r_t ∈ [−5, −1)、j 20〜60) の |U|・T・μt の A と B の差の時系列、各腕の `check_convergence`・`check_quasisteady` (この 30,000 step の区間)、commit の診断 (`FORGE_DIAG_COMMIT_LOSS=2500`、A のみ)。

### 6.19 §6.18 の結果 (2026-10-10)

- run: `run_0460_ab_f32cont` (A)・`run_0461_ab_f32cont2` (A2)・`run_0462_ab_fp64switch` (B)、起点は `run_0451_v4_f32` の res_140000、各 30,000 step、同時に投入。台本 `case/45.isobutane_m6_d155/v6ab.sh` (1 回目は `restart_field.py` の座標の照合で準備に失敗し、直して回し直した)、判定 `v6ab_judge.py` (`_band_ab/cold_pair/v6ab_judge.json`)。
- **判定 (事前登録のとおり): 判別不能**。θ_r の 3 断面は「支持」の条件を満たしたが、Q_w で A の動き (+0.227 %) が 0.1·|D| (0.143 %) を超え、4 量すべてはそろわなかった。
  - **訂正 (2026-10-10、codex diagnose `notes/reviews/2026-10-10-float-v6-redesign-diagnose.md` の Minor)**: Q_w は A の静止の条件だけでなく、効果と再実行の差の条件も満たしていない (|Δ_B − Δ_A| = 1.340 % に対し 10·|Δ_A2 − Δ_A| = 2.250 %)。「判別不能」は変えない。A の静止の条件だけを緩めれば支持になる、とは読まない。

| 量 | D (V4 の差) | Δ_A | Δ_A2 | Δ_B | 支持の条件 |
| --- | --- | --- | --- | --- | --- |
| θ_r(40) | +1.819 % | −0.052 | −0.039 | −2.183 | 満たす |
| θ_r(70) | +1.432 % | +0.059 | +0.010 | −1.663 | 満たす |
| θ_r(94) | +1.288 % | +0.089 | −0.026 | −1.425 | 満たす |
| Q_w | +1.434 % | +0.227 | +0.002 | −1.113 | A の動きが 0.1·|D| を超える (B は 0.78·|D| 戻った) |

- **記録**:
  - θ_r は、B が 1 万 step から FP64 の値の向きへ動き始め、3 万 step で V4 の FP64 の値の位置を少し越えた。A・A2 は ±0.1 % 以内。
  - Q_w の A と A2 の差は 0.225 % で、Q_w の再実行の揺れがこの程度ある (0.1·|D| と同じ桁)。
  - 収縮部 (x/r_t ∈ [−5, −1)、j 20〜60) の A と B の差の平均は、|U| 1.5e-2・T 1.2e-2・μt 1.26 (相対)。
  - 3 本とも `check_convergence` は CHECK FAILURES (NOT CONVERGED)。
- **解釈 (未確定、諮っていない)**: float の B0 が FP64 と 1.3〜1.8 % 違う所で止まるのは、float の精度 (状態・残差・commit の丸め) がその状態を保っているためと見ている。同じ場を FP64 で続けると 3 万 step で FP64 の値に戻り、float で続けると留まる。事前登録の判定は「判別不能」のまま残す (Q_w の条件が再実行の揺れと同じ桁で、判定の作りに余地があった)。

### 6.20 事前登録: hoop の修正 (キー 1) の格子での V4 のやり直し (2026-10-10、§5.1 #15 (1)、実行前)

- **目的**: V4 (§6.16) の float と FP64 の差 (θ_r +1.3〜1.8 %、差の大きい場所は収縮部) に、hoop の偽の半径方向の力 (同じ収縮部の軸の近く、[axisymmetric-freestream-hoop-gauge](axisymmetric-freestream-hoop-gauge.md) §4.10 で修正) が混ざっていたかを分ける。
- **腕**: float = `run_0484_fl1_f32_k1` (`~/forge-fgeom7-f32`)、FP64 = `run_0483_hp7_k1` (`~/forge-fgeom7-fp64`、hoop plan §4.6 の 7 のキー 1 の腕をそのまま使う)。同じソース (commit `0b4dff4e`)・同じ新しい格子 (`_conv9/c45/new64/nozzle.h5`)・同じ出発点 (`run_0183` の res_100000 を `restart_field.py` で移す、9 量)・B0 の実効設定・キー 1。V4 との違いは格子 (区間ごとの r 重み・半割面の端点)・キー 1・commit の診断なし。GPU はキー 0 の腕 (`run_0482_hp7_k0`) と合わせて 3 本で共有する。台本 `case/45.isobutane_m6_d155/fl1.sh`。
- **停止・比較が成り立つ条件・判定**: §6.11 のとおり (式と閾値は同じ)。判定 `fl1_judge.py` (`v4_judge.py` の腕だけを替えたもの)。
- **解釈の区分** (到達して比較が成り立ったとき。θ_r の 3 断面で見る):
  - (a) フラグなし → V4 の差は hoop の偽の力と float の組み合わせで生じていた。#15 (2) の部分ごとの A/B は行わない。
  - (b) 3 断面とも V4 と同じ符号で、大きさが V4 の差の 1/2 以上残る → hoop は主因でない。#15 (2) の設計に進む (上位に諮る)。
  - (c) それ以外 → 部分的。内訳を見てから上位に諮る。
  - Q_w は再実行の揺れが 0.225 % ある (§6.19 の A と A2) ので、Q_w だけのフラグでは区分を決めない。
- **記録**: §6.11 の記録に加えて、V4 の差との並び (`fl1_judge.py` の末尾)。
- **解釈の区分の訂正 (2026-10-11、結果の前、codex diagnose `notes/reviews/2026-10-11-g0-fluctuation-controls-diagnose.md` の Major 4)**: 上の (a)・(b) の因果の読み方 (「差が消えれば hoop と float の組み合わせが原因」「半分以上残れば hoop は主因でない」) は外す。この比較は V4 と格子・ビルド・停止の時点も違うので、hoop の寄与だけを取り出せない。§6.11 の数値の判定 (フラグ) だけを記録し、原因は書かない。あわせて、新格子・fgeom7 では θ_r・Q_w の出力ごとの揺れが V4 より大きい (hoop plan §4.13)。停止規則の「到達」は揺れの小ささを保証しない (平均の ±0.1 % を交互に取る合成の系列でも `REACH` になる、codex が `m9_watch.py` の判定関数で確認)。FP64 の腕 `run_0483` が上限 20 万 step までに到達しなければ、§6.11 のとおり判定不能とする (閾値を緩めない・延長しない)。

- **判定器の訂正 (2026-10-10、結果を見る前、codex diagnose の Major 6)**: `fl1_judge.py` (と元の `v4_judge.py`) は残差を末尾 2 万 step しか見ず、系列の Q_w が NaN でも到達と判定し差の比較も素通りした (合成入力で確かめた)。また総合の判定を正式ツールより前に確定し、`check_quasisteady` に run_dir を渡して組み込みの量まで判定させていた。`fl1_judge.py` を直した: 全出力の θ_r・Q_w・欠損の有限性、残差の全期間・全列の有限性、正式ツールが実行できないときの判定不能、`check_quasisteady` は系列だけを渡す。閾値と区分は変えていない。
  - 直した判定器で V4 (`run_0451`・`0452`) を再評価した (`_band_ab/cold_pair/v4_judge_v2.json`): 判定と数値は §6.16 と同じ (到達、フラグあり)。NaN を仕込んだ合成入力では判定不能を返す。

### 6.21 事前登録: case/48 の float と FP64 の比較 (2026-10-10、codex diagnose、実行前)

- **諮問**: codex (diagnose) に諮った: `notes/reviews/2026-10-10-float-v6-redesign-diagnose.md` — V6 を精度の比較に寄せるのは採用するが、回帰の責務の全面的な置き換えは却下。case/48 の精度の A/B を、判定器と実効設定を直してから §6.20 を待たずに回す。採否は §6.1 の表。
- **目的**: ユーザの目標「double をやめたい」(2026-10-10) に対し、case/48 (平面の冷却平板 M4、SST、等温壁 300 K) で、精度だけを変えたときに float と FP64 の判定量が許容内に入るかを見る。
  - 解釈の範囲 (Major 2): 許容外の差が出れば「軸対称・TP・ライン陰解法・収縮部がなくても精度差が現れる」まで言える。許容内なら「この case/48 の構成では許容外の差は残らない」まで。どちらでも「SST・低 Re 壁が原因」「軸対称が原因」とは確定しない (case/48 は平面化のほかに CPG・点陰解法・圧力勾配なし・格子・壁解像も case/45 と違う)。
- **腕**: float ×2 (`run_0057_pv_f32a`・`run_0058_pv_f32b`、`~/forge-fgeom7-f32`)、FP64 ×2 (`run_0059_pv_f64a`・`run_0060_pv_f64b`、`~/forge-fgeom7-fp64`)。同じソース (commit `0b4dff4e`)・同じ double の幾何の格子 (`case/45/_conv614/c48/new64/m.h5`)・同じ入力 HDF5 (`run_0025_B_tw300_y3_fx05/res_48000.h5` の 7 保存量を `restart_field.py` で 1 回だけ移し、4 本で共有。FP64 はその値を正確に広げて読む)。各 48,000 step 固定、2,000 step ごとに出力。台本 `case/48.flat_plate_cooled_m4/c48_prec.sh`。
- **実効設定** (Major 3): `run_0048` (B 300 K、生産 SST、`blockDPLUR 1`・`cfl_pseudo 2`) から、`mesh.scalarGradient` を `gg` (S2 の双子の名残り) から `lsq` (node の既定) に明示した。ほかは複製元のまま。起動時の実効値 (`scalarGradient`・`slauWallNormalChi`・`limiterScaled`・`venkatK`・`implicitSolvePrecision`) を記録する。起点の `run_0025` から設定が変わるので、新しい判定区間として扱う。`run_0025` 自体は `NOT CONVERGED (stalled/plateau)` であり、「float で収束した場」ではない。
- **判定量** (Major 4): `cooled_plate_eval.py` の `station()`・`plate_integrals()` をそのまま呼ぶ `tools/c48_prec_series.py` で、各出力の x = 0.3/0.6/0.9 の θ・δ*・Cf・q_w と CD・HF (x ∈ [0.002, 0.998] の 200 点の台形積分) を CSV にする (同じ double の座標・同じ節点列)。縁の点 je と縁の高さも記録する。`--series` (公差 1〜1.5 %) と閉合は判定に使わない。
- **判定** (Major 5、`tools/c48_prec_judge.py`):
  - 前提 (満たさなければ判定不能): 4 本とも完走、残差の全期間・全列が有限、系列が 24 点そろって有限、正式ツールが実行できる。
  - 窓は 24,000〜48,000 step の 13 点。各 run の窓を `check_quasisteady.py --series-csv --tail 1.0` で判定 (θ・δ* は `--drift 0.0001 --osc 0.0001`、ほかは 0.0002)。4 本すべてで全量が STEADY でなければ判別不能。48,000 step で静定しなくても窓を選び直さない。
  - 各量で μ32・μ64 = 2 本・全窓点の平均、D = |μ32 − μ64|/|μ64|、b = 各精度内で平均からの最大偏差、E = (b32 + b64)/|μ64|。τ = 0.0005 (θ・δ*)、0.001 (Cf・q_w・CD・HF)。精度間の比較の暫定の許容差であって、物理の正しさの保証ではない。
  - 全量で D + E ≤ τ → 「登録窓の対象量は許容内」。いずれかで D − E > τ → 「許容外の精度差を検出」(D − E ≥ 0.005 なら「0.5 % 以上の明瞭な差」と追加表示)。それ以外は判別不能。4 本の変動幅による判定であって、統計的な信頼区間とは呼ばない。
  - `check_convergence` は記録する。NOT CONVERGED なら「同じ定常解に収束」「FP64 を全面的にやめられる」とは書かない。
- **記録**: 実効設定の起動時の出力、je の範囲、各 run の準定常の詳細、`check_convergence`。再実行の差の物差しとして `run_0028/0029` (CHT)・`run_0030/0031` (層流)・S2 (gg/lsq の設定差) の数値は使わない (codex の異議を採用)。

### 6.22 §6.21 の結果 (2026-10-10)

- run (AWS `~/forge-wallfit/case/48.flat_plate_cooled_m4/`): float `run_0057_pv_f32a`・`run_0058_pv_f32b`、FP64 `run_0059_pv_f64a`・`run_0060_pv_f64b`。4 本とも 48,000 step を完走し、残差の全期間が有限、系列 24 点がそろった。起動時の出力で `scalarGradient: lsq (explicit)`・`slauWallNormalChi 1`。`check_convergence` は 4 本とも CHECK FAILURES (case/48 の既知の残差の床)。判定 `c48_prec_judge.json`。
- **判定 (事前登録のとおり): 判別不能** (窓の対象量が 4 本すべてで STEADY でない)。FP64 の 2 本は全量 ALL STEADY。float の 2 本は閾値 1e-4 の組で NOT ALL STEADY: θ(0.3)・δ*(0.3) が OSCILLATING (両方)、δ*(0.6) が DRIFTING (`run_0057` だけ。縁の点 je が窓の中で 74〜76 を行き来した)。閾値 2e-4 の組 (Cf・q_w・CD・HF) は 4 本とも ALL STEADY。
- 量ごとの D・E: δ*(0.6) は D 3.5e-5・E 8.0e-4 で判別不能、ほかの 13 量は D + E ≤ τ。D の最大は δ*(0.3) の 4.0e-5。
- 記録: 窓の中の揺れの幅は、float が θ・δ*・q_w で相対 0.7〜2.0e-4、FP64 が 1e-7 程度 (3 桁小さい)。FP64 の 2 本は窓の平均で 1e-14 まで一致した。float の窓の平均と FP64 の差は、全量で 7e-5 以下。
- **codex (diagnose) に諮った**: `notes/reviews/2026-10-10-c48-precision-result-diagnose.md` — 追加の run は行わず、共通の double の座標に揃えた保存済みの場で、縁の点 je を動的と固定で比べる後処理の A/B を 1 回だけ行う。採否は §6.1 の表。
  - **訂正 (Major 1、コードと実測で確かめた)**: 抽出器は各結果ファイルの `/MESH/COORD` で判定量を作っていた。float のビルドの結果ファイルは座標を float32 で書く (共通の double の座標と最大 2.97e-8 m の差、節点の並びは同じ)。§6.21 の「共通の double の座標・同じ節点列」を満たしていないので、§6.23 で再集計する。上の判定はその前の値として残す。
  - **判別不能の理由の訂正 (Major 2)**: 準定常の条件だけではない。δ*(0.6) は D + E = 8.35e-4 = 1.67τ で、変動幅を含む比較でも判別不能になる。「準定常の条件と D + E の条件の 2 つを満たさない」と記録する。
  - **書いてよい範囲 (Major 3)**: 「この 4 本・登録の窓・指定の 14 量の抽出値では、窓の平均に % 級の差は観測されなかった」まで。定常的な精度差が小さいとも、SST・commit・残差の集積・面の演算を #15 の候補から外すとも言わない。#15 の候補の順も、これだけでは変えない。
  - 揺れの要約の訂正 (Minor): float の揺れは θ・δ*(0.3)・q_w で相対 0.7〜2.0e-4 (FP64 の約 1000 倍)、`run_0057` の δ*(0.6) は例外で 8.7e-4 (FP64 の約 5400 倍)。ここでの「時間の揺れ」は定常計算の反復の履歴の変動。

### 6.23 事前登録: case/48 の再集計と縁の点 je の後処理の A/B (2026-10-10、codex diagnose、実行前、ソルバの計算 0 step)

- **再集計** (Major 1・4): `tools/c48_prec_series.py --coords-from _prec/mesh_init.h5` で、4 本の窓の場から共通の double の座標で判定量を作り直す (各結果ファイルの座標が、共通の座標をその型に丸めたものと全節点で一致することを確かめてから置き換える)。判定器 `tools/c48_prec_judge.py` は前提の検査を強める (系列 24 点の全体の有限性、`check_convergence` が実行できること。どちらかを欠けば判定不能)。§6.21 の規則・閾値・窓は変えずに判定し直す。§6.22 の値は消さずに並べる。
- **je の A/B** (第 1 仮説「`run_0057` の δ*(0.6) の大きな揺れは、縁の点 je の整数の切り替えで主に増幅されている」、確度 中):
  - A: 登録どおり、各場から je を選ぶ (再集計の系列)。
  - B: x = 0.6 の je を全 4 本・全時刻で 76 に固定する (縁の値 ρe・ue も同じ節点から取り直す)。ほかの断面と量は A と同じ。
  - 主判定: `run_0057` の δ*(0.6) の相対変動幅 (最大 − 最小)/平均 の比 B/A ≤ 0.5 → 「切り替えが変動幅の過半を説明する」を支持。B/A > 0.5 → 第 1 仮説を棄却。
  - 記録: B での `check_quasisteady` (同じ閾値) と D・E。**B の結果で §6.21 の判定を合格に変えない** (抽出の感度の診断)。x = 0.3 の OSCILLATING は je と無関係 (je は一定) なので、別に残る (第 2 仮説「float の反復の状態に変動が残る」、確度 中)。
- これで case/48 は止める (窓の延長・閾値の緩和はしない)。#15 は `case/45` の `run_0483_hp7_k1`・`run_0484_fl1_f32_k1` (§6.20) の結果で進める。

### 6.24 §6.23 の結果 (2026-10-10、計算 0 step)

- `c48_prec_series.py --self-test`: je を固定しない `station_je` は `station()` と 4 本とも全量で一致。共通の double の座標への置き換えは、4 本とも全節点で対応を確かめた (float の結果ファイルの座標は、共通の座標を float32 に丸めたものと一致)。
- **再集計** (`c48_prec_judge_dbl.json`、窓 13 点、共通の double の座標): 判定は §6.22 と同じ **判別不能** (準定常: float の 2 本で θ(0.3)・δ*(0.3) が OSCILLATING、`run_0057` の δ*(0.6) が DRIFTING。D + E: δ*(0.6) が D 3.5e-5・E 8.0e-4 で τ を超える)。D の値は §6.22 から 3 桁目が動く程度 (例 Cf(0.3) 7.75e-6 → 7.78e-6)。座標の精度は判定を変えなかった。
- **je の A/B** (`c48_je_ab.json`): **支持**。`run_0057` の δ*(0.6) の相対変動幅は A (動的 je、74〜76) 8.67e-4 → B (je 76 に固定) 6.39e-5 で、B/A = 0.074 (≤ 0.5)。θ(0.6) は 0.761。ほかの 3 本は je が窓の中で 76 のままなので A = B。
  - B での判定 (`c48_prec_judge_dbl_je76.json`、記録): δ*(0.6) は D 2.6e-6・E 3.8e-5 で D + E ≤ τ。総合は準定常 (x = 0.3) のために判別不能のまま。**§6.21 の判定は変えない**。
- **まとめ (書いてよい範囲)**: case/48 の §6.21 は判別不能。理由は (a) float の反復の履歴に、θ・δ*・q_w で相対 1e-4 の揺れが残ること (FP64 は 1e-7。x = 0.3 では je は一定なので抽出の切り替えでは説明されない。第 2 仮説は未検証のまま)、(b) δ*(0.6) の変動幅の大半が縁の点 je の切り替えによること (B/A 0.074)。この 4 本・登録の窓・指定の 14 量の抽出値では、float と FP64 の窓の平均に % 級の差は観測されなかった (最大 7e-5)。case/45 の候補 (SST・commit・残差の集積・面の演算・軸対称・ライン陰解法・TP・収縮部) はどれも外さない。
- case/48 はここで止める。#15 は `run_0483_hp7_k1`・`run_0484_fl1_f32_k1` (§6.20) の結果で進める。

### 6.25 §6.20 の結果 (2026-10-11)

- run: float `run_0484_fl1_f32_k1` (到達 120,000 step)、FP64 `run_0483_hp7_k1` (到達 145,000 step、res_145000 の sha256 `217c7fc6…`)。判定 `fl1_judge.py` (結果の前に直した版、§6.20 の末尾) → `_band_ab/cold_pair/fl1_judge.json`。前提の検査 (全期間の系列・残差の有限性、正式ツールの実行) はすべて通った。
- **判定 (§6.11 の規則のとおり): 到達、フラグあり**。

| 量 | float | FP64 | 差 | V4 (§6.16) の差 |
| --- | --- | --- | --- | --- |
| θ_r(40) | 0.079589 | 0.078111 | +1.892 % | +1.819 % |
| θ_r(70) | 0.117096 | 0.115199 | +1.647 % | +1.432 % |
| θ_r(94) | 0.129852 | 0.128021 | +1.431 % | +1.288 % |
| Q_w | 9.4725 MW | 9.2764 MW | +2.114 % | +1.434 % |

- 残差の末尾 2 万 step の中央値の比 (float/FP64): `rms_ro` 1.09、`rms_roUx` 1.29、`rms_roUy` 1.13、`rms_roe` 1.16、`rms_roK` **2.83** (フラグ)、`rms_roOmega` 1.27、`rms_roY*` 1.09。2π·Σ|res_ro| の比 1.02。V4 (2.6〜3.9 倍) より小さい。
- 記録: 両腕とも `check_convergence` は CHECK FAILURES、θ_r・Q_w の系列は NOT ALL STEADY。末尾の出力ごとの揺れは float で θ_r(40) 0.092 %・Q_w 0.194 %、FP64 で 0.050 %・0.077 % (hoop plan §4.13)。到達時の差 (1.4〜2.1 %) は、この揺れより 1 桁大きい。
- **書いてよい範囲** (§6.20 の訂正のとおり): 「hoop の修正の格子・ビルドでも、停止時の float と FP64 の θ_r・Q_w の差は V4 と同じ向き・同じ桁で残った」まで。hoop の寄与を取り出す比較ではないので、原因は書かない。#15 の次の段 (float のビルドの中で部分ごとに double にする A/B) の設計は、上位に諮ってから書く。

### 6.26 事前登録: 共通の Q32 の起点からの float と FP64 の軌跡の A/B (2026-10-11、codex diagnose、実行前、§5.1 #15 (2) の最初の段)

- **codex (diagnose) に諮った**: `notes/reviews/2026-10-11-float-cause-step2-design-diagnose.md` — 共通の Q32 の初期状態から、全領域・固定 30,000 step の float と FP64 の軌跡の A/B を、float の再実行 1 本付きで行う。一部だけ double にする A/B と解析領域の短縮は先にしない。採否は §6.1。
  - SST の commit の丸めを原因に絞らない (S_lostfrac は切り上がりも数え、場所も符号も分からない)。SST は重点的に観測するが、発生源を限らない。
  - V3 の PASS は旧実装からの改善であって、残差の一致ではない (ω の残差の誤差は参照のノルム比 0.37 が残る、§6.10)。commit の前の作用素の差も候補に残す。
  - 解析領域の短縮は、この結果で差の発生域を見てから、切る面の条件 (法線のマッハ数・逆流・壁際の亜音速域・出口条件) と、保持する領域の指標が全領域の場合からほとんど変わらないことを確かめてから考える。
- **起点**: `run_0483_hp7_k1` の res_145000 (到達の場、sha256 `217c7fc6…`)。新格子の `nozzle.h5` に `restart_field.py` で 9 保存量を移し (ビット一致)、**保存量 9 つだけを float32 に丸めて double に戻した共通の入力** (Q32) を 1 つ作る。幾何と `wall_dist` は変えない。3 本ともこの入力を読む (FP64 には double(Q32) を、float には Q32 をそのまま)。
- **腕**: A = FP64 (`~/forge-fgeom7-fp64`)、B = float (`~/forge-fgeom7-f32`)、B′ = B の再実行 (再実行の差の対照)。同じソース (commit `0b4dff4e`)・格子・境界条件・B0 の実効設定・キー 1・環境変数 (`FORGE_CUDA_BLOCKSIZE=128`・`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500` を 3 本とも)。インスタンス B で 3 本同時。台本 `case/45.isobutane_m6_d155/tr.sh`、判定 `tr_an.py`。
- **期間・出力**: 各 30,000 step 固定、500 step ごとに出力し全出力を残す。到達で止めない。結果を見て延ばさない。`output: {level: 1, extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, omg_prod, omg_dest, omg_cross, omg_trans, omg_axisym]}`。
- **主判定** (θ_r(40/70/94)。Q_w は独立した副判定として記録し、主判定には混ぜない):
  - 尺度 D = §6.25 の停止時の差 (+1.892, +1.647, +1.431) % (真値ではなく効果量の尺度)。窓は W1 = 20,000〜25,000、W2 = 25,000〜30,000 step (各 11 点)。各腕・各窓の平均の、起点 (step 0 の出力、なければ最初の出力) からの移動を ΔA・ΔB・ΔB′ (D と同じ分母 = A の起点の値で割った %) とする。
  - **精度依存を支持**: 両窓・3 断面とも、ΔB と ΔB′ がともに ΔA から float の向き (D の符号) へ 0.5|D| 以上離れ、|mean(ΔB, ΔB′) − ΔA| が |ΔB − ΔB′| の 10 倍を超え、かつ |ΔA| ≤ 0.1|D|。→ 共通の過渡だけという説明を退け、「float に依存する作用素・更新の差が偏りを作る」(第 1 仮説) を支持する。SST の commit の証明にはしない。
  - **共通の過渡を支持**: 両窓・3 断面とも、|mean(ΔB, ΔB′) − ΔA| ≤ 0.1|D| で、ΔA と mean(ΔB, ΔB′) が同じ向きに 0.5|D| 以上動く。→ 第 1 仮説の「この期間で支配的」を棄却し、共通の過渡・restart・揺れ (第 3 仮説) を優先する。
  - それ以外 (両方とも動かない、再実行の差が大きい、窓の間で反転する、非有限・欠損) は判別不能。係数は暫定の診断の閾値で、信頼区間ではない。
- **場の記録** (判定には使わない): ρk・ρω・k・ω・μt・ρ・P・T・|U| を、領域 (収縮部全体 x/r_t ∈ [−15, −1)、その内部 x/r_t ∈ [−5, −1)・j 20〜60、軸の近く j 0〜8、壁際 j 117〜120、下流 x/r_t ≥ 1) ごとに、起点の場の領域の RMS で正規化した差の L2 で追う。局所の差の立ち上がりは、「‖B − A‖ が ‖B − B′‖ の 10 倍を超える出力が 3 回続いた最初の区間」とする (500 step より細かい先後は分からないとする)。先に離れた量は原因の確定ではない。あわせて ω の収支 (`omg_*`) の差、commit の診断の集計、全期間の非有限の検査、区間を明示した `check_convergence` と系列の `check_quasisteady` を記録する。
- 所要の見込み: 専有の単価で float 2 本と FP64 1 本の 30,000 step は合計約 37 分の GPU の仕事量 (同時投入で専有の単価は保たれないので、実際はこれより長い)。出力は 3 本で 20 GB 前後。

### 6.27 §6.26 の結果 (2026-10-11)

- run (AWS B `~/forge-wallfit/case/45.isobutane_m6_d155/`): A `run_0487_tr_fp64`・B `run_0488_tr_f32`・B′ `run_0489_tr_f32b`。共通の入力 `_tr/q32_init.h5` (sha256 `ba859e02…`、保存量 9 つを float32 に丸めて 455 万個の値が変わった)。3 本とも 30,000 step を完走 (RUN_RC 0)、起動時の出力でキー 1 と追加の出力 11 項目を確認。1 回目の投入は起点の場の隣に化学種の記録が無く `restart_field.py` が拒否した (台本の誤り、記録を置いて回し直した)。判定 `tr_an.py` → `_band_ab/cold_pair/tr_judge.json`。
- **判定 (事前登録のとおり): 精度依存を支持** (θ_r の 3 断面・2 窓の 6 条件すべて)。

| 窓 | 断面 | ΔA | ΔB | ΔB′ | \|mean(ΔB,ΔB′) − ΔA\| | 10·\|ΔB − ΔB′\| | 0.5\|D\| |
| --- | --- | --- | --- | --- | --- | --- | --- |
| W1 | θ_r(40) | −0.056 % | +1.577 % | +1.439 % | 1.564 | 1.378 | 0.946 |
| W1 | θ_r(70) | +0.012 % | +1.229 % | +1.120 % | 1.163 | 1.091 | 0.824 |
| W1 | θ_r(94) | +0.000 % | +0.995 % | +0.908 % | 0.951 | 0.866 | 0.716 |
| W2 | θ_r(40) | −0.121 % | +1.883 % | +1.776 % | 1.950 | 1.064 | 0.946 |
| W2 | θ_r(70) | −0.056 % | +1.512 % | +1.426 % | 1.525 | 0.853 | 0.824 |
| W2 | θ_r(94) | −0.054 % | +1.267 % | +1.194 % | 1.285 | 0.728 | 0.716 |

- 副判定 Q_w (主判定に混ぜない): ΔA +0.085/+0.169 %、ΔB +2.816/+2.749 %、ΔB′ +2.524/+2.400 % (W1/W2)。float は FP64 から D の向きへ 2.4〜2.7 % 離れたが、10·|ΔB − ΔB′| (2.92・3.49) を超えないので、副判定の支持の条件は満たさない。
- **書いてよい範囲**: 同じ状態から出発すると、FP64 はその場に留まり (θ_r の移動 ≤ 0.12 %)、float の 2 本は 2〜3 万 step で以前の float の停止時の差の 6〜13 割まで戻る。**float に依存する作用素・更新の差が偏りを作る** (第 1 仮説、支持)。SST の commit の証明ではなく、どの演算かは未特定。
- **場の記録** (判定には使わない): 「‖B − A‖ > 10‖B − B′‖ が 3 出力続いた最初の step」の早い順は、P・|U| の軸の近く (2,500〜3,000)、ω の壁際 (3,500)、ρk・k・ω の収縮部の内部 (7,000)、ρk・k の下流 (8,500〜9,000)、ρω の収縮部の内部 (9,000)、k の収縮部全体 (11,000)、ρ の下流 (12,500)、ρω の壁際・下流・収縮部全体 (14,000〜15,500)、ρk の壁際・|U| の下流 (22,000)。
  - **注意**: 出発の直後は B と B′ が同じ入力からほとんど離れていないので、‖B − B′‖ が小さく、この判定は早い時刻で立ちやすい。P・|U| の軸の近く・ω の壁際は step 30,000 では ‖B − A‖ ≲ ‖B − B′‖ で、精度の差とは言えない。
  - step 30,000 でも ‖B − A‖ が ‖B − B′‖ の 3 倍以上なのは、k・ρk・ω・ρω の収縮部 (内部・全体)・下流・壁際と、ρ・|U| の下流。収縮部の P・T・|U| は float どうしの差より離れない。
- 記録: 3 本とも `check_convergence` は CHECK FAILURES、20,000〜30,000 step の θ_r・Q_w の系列は `check_quasisteady` で ALL STEADY (既定の閾値)。commit の診断 (float、step 30,000、内部) は ρ 0.24 %・ρu 0.41 %・ρk 0.26〜0.46 % (FP64 は 1e-11 の水準)。
- **次の段は諮ってから** (§5.1 #15、hoop plan §4.13・§4.14 の揺れの結果とまとめて)。

### 6.28 事前登録: SST の輸送の更新を止めた軌跡の A/B (2026-10-11、codex diagnose、実行前)

- **codex (diagnose) に諮った**: `notes/reviews/2026-10-11-next-step-fluct-and-float-diagnose.md` — §6.27 の前提 (全期間の有限性・ツールの実行) を確かめたうえで、`FORGE_FREEZE_TURB=1` の FP64・float・float 再実行を比べ、「SST の更新を止めても精度差が残るか」を見る。採否は §6.1。
  - §6.27 の確かめ直し (Major 3、`tr_recheck.py` → `_band_ab/cold_pair/tr_recheck.json`): 3 本とも全期間の残差 14 列と全 61 出力の全データセットが有限、正式ツールは実行できていた → **§6.27 の判定は確定**。ただし W1 の効果/再実行の差は 10.7〜11.3 倍で、10 倍の基準に対する余裕は小さい (codex の再計算)。`check_quasisteady` の ALL STEADY は既定の閾値 (ドリフト 5 %・揺れ 10 %) なので、揺れが小さいことの根拠にしない。
  - **定義の訂正** (Major 1): `FORGE_FREEZE_TURB=1` は SST の更新のブロックを飛ばすだけで (`main.cpp:2380`)、k = ρk/ρ・ω = ρω/ρ は毎回計算し直され、壁の ω・ρω の上書きも残り、μt は密度・速度勾配・物性で変わる。「μt を固定する試験」ではなく「**SST の輸送の更新を止める試験**」とする。μt が一定であることを合格の条件にしない。
- **組み方**: §6.26 と同じ共通の入力 `_tr/q32_init.h5` (sha256 `ba859e02…`)・格子・キー 1・実効設定。3 本とも `FORGE_FREEZE_TURB=1` (`FORGE_DIAG_COMMIT_LOSS=500` も付ける。§6.26 の `FORGE_OMEGA_BUDGET` は付けない)。A = FP64 `run_0491_fz_fp64`、B = float `run_0492_fz_f32`、B′ = float の再実行 `run_0493_fz_f32b`。各 30,000 step 固定・500 ごと。窓 W1 = 20,000〜25,000・W2 = 25,000〜30,000 step。台本 `fz.sh`、判定 `fz_an.py` (インスタンス B)。
- **判定** (θ_r の 3 断面で主判定、Q_w は副判定として同じ式で記録): Δ は §6.26 と同じ (窓の平均 − 自分の起点、A の起点で割った %)。凍結した FP64 自身の動きを差し引くため、δB = ΔB − ΔA、δB′ = ΔB′ − ΔA を比べる (|ΔA| ≤ 0.1|D| は要求しないが ΔA は記録する)。D = (+1.892, +1.647, +1.431) %。
  - **差が残る**: 両窓・3 断面とも δB・δB′ が D の向きへ 0.5|D| 以上、かつ |mean(δB, δB′)| > 10|δB − δB′| → 「観測した規模の偏りには SST の更新が要る」(第 1 仮説) の必要性を棄却し、「SST の更新を止めても残る経路 (平均流・物性・μt の再評価など) だけで偏りが育つ」(第 2 仮説) を支持。
  - **差が抑えられる**: 両窓・3 断面とも max(|δB|, |δB′|) + 10|δB − δB′| ≤ 0.1|D| → 「SST の更新を止めても従来の規模の差が育つ」説を退ける。結論は「この精度差の発達には SST の更新・フィードバックが要った」まで。**SST の中の演算が発生源だとは言わない**。
  - それ以外 (中間・窓の間の不一致・欠損・非有限) は判別不能。期間を延ばさない。閾値は診断用の暫定値。
  - 判定の前に、全期間の残差の全列と全出力の全データセットの有限性、正式ツールの実行を確かめる (欠ければ判定不能)。
- **凍結が効いた証拠** (記録): 境界の節点 (軸 j 0・壁 j 120・入口 i 0・出口 i 末) を除く内部の ρk・ρω が、起点と step 30,000 で変わらないこと (変われば凍結が効いていないとして判定不能)。k・ω・μt の変化は別に記録する。
- 次の段 (部分的な double 化) は、同じ入力での「残差・要求の更新・実際の更新」のどこに説明力のある差があるかを見てから決める (codex: S_lostfrac は符号を持たないので累積の偏りを示さない)。

### 6.29 §6.28 の結果 (2026-10-11)

- run (AWS B): A `run_0491_fz_fp64`・B `run_0492_fz_f32`・B′ `run_0493_fz_f32b`。3 本とも 30,000 step を完走 (RUN_RC 0)、全期間の残差と全出力が有限、`check_convergence` は 3 本とも実行でき CHECK FAILURES。判定 `fz_an.py` → `_band_ab/cold_pair/fz_judge.json`。
- **凍結の確認**: 3 本とも、内部の ρk・ρω は起点と step 30,000 で全節点が同じ (変化 0)。k・ω の相対 L2 の変化は A 9.8e-4・2.9e-4、B 3.5e-3・3.3e-4、B′ 3.1e-3・3.2e-4 (ρ の変化による)。μt の変化は起点の出力 (res_0) の μt が計算前の値なので測れない (記録の不備)。
- **判定 (事前登録のとおり): 判別不能**。「差が残る」は 6 条件のうち W2 の θ_r(40) の 1 つだけ成立、「差が抑えられる」は 0。

| 窓 | 断面 | ΔA | ΔB | ΔB′ | δB | δB′ | \|mean δ\| | 10\|δB − δB′\| | 0.5\|D\| |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| W1 | θ_r(40) | −0.116 % | +1.055 % | +0.812 % | +1.171 | +0.927 | 1.049 | 2.436 | 0.946 |
| W1 | θ_r(70) | −0.096 % | +0.870 % | +0.641 % | +0.966 | +0.737 | 0.852 | 2.288 | 0.824 |
| W1 | θ_r(94) | −0.126 % | +0.675 % | +0.492 % | +0.801 | +0.619 | 0.710 | 1.823 | 0.716 |
| W2 | θ_r(40) | −0.028 % | +1.000 % | +1.061 % | +1.028 | +1.089 | 1.059 | 0.611 | 0.946 |
| W2 | θ_r(70) | −0.026 % | +0.848 % | +0.949 % | +0.873 | +0.974 | 0.924 | 1.011 | 0.824 |
| W2 | θ_r(94) | −0.055 % | +0.696 % | +0.802 % | +0.752 | +0.857 | 0.804 | 1.058 | 0.716 |

- 副判定 Q_w: δB −0.06/+0.08 %、δB′ +0.56/+0.47 % (W1/W2)。float どうしの差が大きく、向きもそろわない。
- **記録 (判定ではない)**: SST の輸送の更新を止めても、float の 2 本は D の向きに θ_r で 0.6〜1.2 % 動いた (凍結なしの §6.27 の δ の約 5〜8 割)。凍結した FP64 は ±0.13 % 以内。一方で float どうしの差 (|δB − δB′|) は W1 で 0.18〜0.24 % と、凍結なし (§6.27、0.07〜0.14 %) より大きい。「差が抑えられる」(SST の更新が要った) は明らかに当たらない。「差が残る」は再実行の差の 10 倍の条件で外れた。
- **codex (diagnose) に諮った**: `notes/reviews/2026-10-11-float-after-freeze-diagnose.md` — 判別不能は採用。採否は §6.1。
  - **訂正 (Major 1・Minor)**: 「差が残る」が外れたのは 10 倍の条件だけではない。W1 では B′ の δ が 3 断面とも 0.5|D| に届かない (0.927 < 0.946、0.737 < 0.824、0.619 < 0.716)。上の「float の 2 本は D の向きに 0.6〜1.2 % 動いた」は Δ と δ の混同で、正しくは「凍結した FP64 の動きを差し引いた float の θ_r の差 δ は +0.62〜+1.17 ポイント」(自分の起点からの ΔB・ΔB′ は +0.49〜+1.06 %)。「凍結なしの約 5〜8 割」は記述の比であって寄与率ではない。**SST の更新の必要性 (第 1 仮説) も、残る経路だけで育つこと (第 2 仮説) も、採否は未確定**。
  - **確かめ直し (Major 2、`fz_recheck.py` → `_band_ab/cold_pair/fz_recheck.json`)**: 抽出した θ_r・Q_w の系列は 3 本・61 出力とも有限で点数がそろい、3 本の格子の座標は同じ。20,000〜30,000 step の系列の準定常は 3 本とも NOT ALL STEADY (θ_r は閾値 1e-4、Q_w は 2e-4、凍結の後の移り変わりの途中)。判定は判別不能のまま。
  - float どうしの差 (B − B′) は W1 の +0.18〜+0.24 ポイントから W2 の −0.06〜−0.11 ポイントへ全断面で符号が反転した (Major 4)。判定や閾値を変える材料にしない。

### 6.30 事前登録: 共通の Q32 の状態からの最初の 1 回の更新を段ごとに測る診断 (2026-10-11、codex diagnose、実装の前)

- **目的**: float と FP64 の最初の 1 回の平均流の更新の差が、commit の前の演算 (残差・要求の更新) と commit の丸めのどちらで主に生じるかを分ける。分かるのは差の発生の段までで、30,000 step 後の θ_r の偏りの主因は確定しない。
- **組み方**: 共通の入力 `_tr/q32_init.h5`・格子・キー 1・実効設定・環境変数を固定し、両方で `FORGE_FREEZE_TURB=1`。A = FP64、B = float。各 1 step を 3 回 (毎回同じ Q32 から始める。連続した 3 step にしない)。
- **記録** (新しい診断 `FORGE_DIAG_STAGE_DUMP`、§5.1 #17。数値は変えない): 平均流の保存量 5 つ (ρ・ρu・ρv・ρw・ρE) について、組み立て後の残差 R、実際の commit の基準 b、緩和・ガードを当てた後の要求の更新 d、commit の直後の状態 q。double で書き出す。
  - **時点の訂正 (実装の報告、2026-10-11)**: 境界条件・EOS の床・TP の roe の組み直し・壁と軸のピンは、commit の後ではなく、step の頭の組み立ての前処理 (`main.cpp:1891-1925`) で作業配列 ro..roe に当たり、commit の基準 roN には入らない (軸のピンが roUyN・roeN を射影する分は例外)。その効果は同じ step の Q_asm − b に出る (Q_asm を足して記録する)。q_final − q は commit の後の等温壁のピンだけ (SST を止めているので)。
- **比較**: 共通の起点 Q₀ に対し a = q − Q₀、e = q − b − d。精度の間の差は Δa = Δb + Δd + Δe (double で照合する)。領域を固定する: 収縮部の内部 (x/r_t ∈ [−5, −1)・j 20〜60)、軸の近く (j 0〜8)、壁の近く (j 117〜120)。保存量ごとに符号つきの和と L2 を出し、異なる保存量を 1 つのノルムに混ぜない。
- **判定** (各領域・各保存量。差は両精度の中の再実行の差の最大の 10 倍を超えることを要求):
  - 結果 A: ‖Δe‖ ≤ 0.1‖Δa‖ → 最初の更新の差は主に commit の前で生じる (第 1 仮説を支持、commit の丸めを主因とする説を退ける)。
  - 結果 B: ‖Δb + Δd‖ ≤ 0.1‖Δa‖ → commit の丸めを優先する。
  - 中間・領域の間の不一致・差が小さい → 判別不能。10 %・10 倍は診断用の暫定の基準。
- **run と台本** (2026-10-11 追記、投入の前): インスタンス B、バイナリ `~/forge-fgeom8-{fp64,f32}` (commit `ba9e889a`)。FP64 = `run_0494_st_f64a`〜`run_0496_st_f64c`、float = `run_0497_st_f32a`〜`run_0499_st_f32c`。診断の有無の確認として、各精度で診断なしの 1 step を 2 本 (`_stage/off_*`) 回し、res_1 の保存量の「診断あり − なし」が「なしどうし」以下であることを確かめる。前提: 6 本の Q0 が double で一致、経路が lineImplicit、commit の自己照合の不一致 0。台本 `stage.sh`、判定 `stage_an.py`。
- やらない: この結果の前に平均流だけの double の commit を実装すること、`qAccumulatorFP64` の軸対称の拒否だけを外すこと (commit の基準の射影を FP64 の正本に反映する必要がある、`main.cpp:3396`)、30,000 step を延ばして合格の窓を探すこと。

### 6.31 §6.30 の結果 (2026-10-11)

- run (AWS B): FP64 `run_0494`〜`0496`・float `run_0497`〜`0499` (各 1 step、`~/forge-fgeom8-*`、commit `ba9e889a`、sha256 FP64 `049cfbfd…`・float `9e70a27d…`)、診断なしの確認 `_stage/off_*`。10 本とも完走。記録は 6 本とも経路 lineImplicit・commit の自己照合 (q == fl(b + d)) の不一致 0・SST 凍結。
- **診断の有無の確認: 合格**。res_1 の保存量の「診断あり − なし」の最大は FP64 7.3e-5 (なしどうし 3.4e-4)、float 2.0 (なしどうし 2.0、ρE の絶対値)。
- **判定 (事前登録のとおり): 判定不能**。前提の「6 本の Q0 が double で一致」を満たさない。ρ・ρu・ρv・ρw は 6 本とも入力とビット一致だが、**ρE だけが step 1 の頭で入力から書き換わっていた** (起動時の初期化で ρE を状態方程式から組み直す往復を、それぞれの精度で計算するため)。FP64 は 243,765 節点で中央値 1.8e-16・最大 9.8e-8 (相対)、float は 105,643 節点で中央値 7.4e-8・最大 5.5e-4。相対の大きいところは ρE がほぼ 0 の節点 (壁際 j 113〜115 の下流、内部エネルギーの負の値と運動エネルギーが打ち消し合う) で、絶対値の変化は成分の丸め 1 つぶん。毎 step の組み立ての前処理での ρE の変化 (Q_asm − b) は FP64 最大 4.3e-10・float 1.5e-7。commit の後の上書き (q_final − q) は等温壁の 4,719 節点だけ。
- **記録 (登録外、`stage_an_ownq0.py` → `_band_ab/cold_pair/stage_ownq0_record.json`)**: 各 run 自身の Q0 を起点にした同じ分解では、3 領域 × 4 保存量の 12 組すべてが「差が小さい」(‖Δa‖ が精度の中の再実行の差の 1.5〜3 倍で、10 倍に届かない)。内訳の比 ‖Δbb + Δd‖/‖Δa‖ は 0.83〜1.00 (差はほぼ要求の更新 d の側、commit の丸め e は小さい) だが、揺れの範囲内なので判定に使わない。残差 R の精度の間の差は、軸の近くで ‖ΔR‖/‖R‖ が 1.1〜1.2 (ρ・ρE)・0.6 (ρu)、壁の近くで 0.06〜0.63、収縮部の内部で 1e-3〜2e-2 (R の再実行の差は測っていない)。
- **書いてよい範囲**: 共通の入力からの最初の 1 回の更新の精度の差は、atomicAdd の順番による再実行の揺れと同じ桁で、この 1 step の診断では差の発生の段を分けられない。ずれ (§6.27) は多くの step にわたる積み重ねで育つ。初期化で ρE を組み直す往復が精度ごとに違うこと自体は、丸めの水準 (相対の大きいところは ρE ≈ 0 の節点)。
- 次の段は諮ってから。

### 6.32 事前登録: block の系の組み立ての精度の感度試験 (2026-10-11、codex diagnose、実行前)

- **codex (diagnose) に諮った**: `notes/reviews/2026-10-11-double-solve-trajectory-diagnose.md` — 追加の比較は採用。ただし「線形の解き全体の double 化」ではなく「block の系の組み立ての精度の感度試験」と定義し直す。採否は §6.1。
  - **範囲の訂正 (Major 1)**: `time.deltaT.implicitSolvePrecision: 1` は `implicit_defect_correction_block_d<ST>` の演算 (対角・ライン外の近傍の寄与・右辺の組み立て) を double にするが、`diag`・`rhs`・`K` の保存は `flow_float` のまま (`timeIntegration_d.cu:952, 1184`)。ライン Thomas の分解と代入は別のカーネルで、既定で double (`:2339, 2417`)。
  - **結論の範囲 (Major 2)**: 「登録の窓の偏りを、この設定の変更で抑えられた/抑えられなかった」まで。発生源・必要性・固定点の偏りは言わない。
- **腕**: §6.26 の float の腕 B = `run_0488_tr_f32`・B′ = `run_0489_tr_f32b` を対照に、同じバイナリ (`~/forge-fgeom7-f32`)・入力 (`_tr/q32_init.h5`)・環境変数 (`FORGE_CUDA_BLOCKSIZE=128`・`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500`)・設定で、`implicitSolvePrecision` だけを 0 → 1 にした C = `run_0590_ds_f32`・C′ = `run_0591_ds_f32b`。各 30,000 step・500 ごと、インスタンス B。A = `run_0487_tr_fp64` は参照の軌跡。台本 `ds.sh`、判定 `ds_an.py`。
- **前提** (Major 4、総合の判定の前に確かめる。欠ければ判定不能): C・C′ の全期間の残差・全出力の有限性と 61 出力、`check_convergence` の実行。C と B の設定の差が `implicitSolvePrecision` の 1 行だけ。バイナリの sha256 が §6.26 の float と同じ。初期化後の保存量 (res_0) が B・B′・C・C′ でビット一致 (§6.31 の ρE の組み直しは同じ float のビルドなら同じはず)。
- **判定** (Major 3、θ_r の 3 断面・W1・W2、Δ と D は §6.26 と同じ): mB = mean(ΔB, ΔB′)、mC = mean(ΔC, ΔC′)、N_B = |ΔB − ΔB′|、N_C = |ΔC − ΔC′|。
  - 結果 A「抑制」: max(|ΔC − ΔA|, |ΔC′ − ΔA|) + 10N_C ≤ 0.1|D|、かつ sign(D)·(mB − mC) ≥ 0.5|D|、かつ |mB − mC| > 10·max(N_B, N_C) → 「block の系の組み立ての丸めが偏りを増幅しており、この部分の double 化だけで登録の窓の偏りを 9 割以上抑えられる」(第 1 仮説) を支持。
  - 結果 B「抑制不十分」: C・C′ がそれぞれ A から D の向きへ 0.5|D| 以上、かつ |mC − ΔA| > 10N_C → 第 1 仮説の「この変更だけで十分」を棄却。
  - 両窓・3 断面で揃わなければ判別不能。Q_w は副判定として記録。暫定の診断の基準で、信頼区間ではない。
- やらない: ライン Thomas の double 化を新しく実装すること、差が残っただけで左辺を候補から外すこと、差が縮んだだけで根治とすること、窓の延長・選び直し。

### 6.33 §6.32 の結果 (2026-10-11)

- run (AWS B): C `run_0590_ds_f32`・C′ `run_0591_ds_f32b` (各 30,000 step、RUN_RC 0)。判定 `ds_an.py` → `_band_ab/cold_pair/ds_judge.json`。前提はすべて満たした: 61 出力・全期間の有限性・`check_convergence` の実行 (CHECK FAILURES)、`run_0488` との設定の差は `implicitSolvePrecision` の 1 行、バイナリの sha256 は §6.26 の float と同じ、初期化後の保存量 (res_0) は B・B′・C・C′ でビット一致。
- **判定 (事前登録のとおり): 判別不能**。結果 A (抑制) は 6 条件中 0。結果 B (抑制不十分) は 6 条件中 5 で成立し、W2 の θ_r(40) だけ |mC − ΔA| 1.860 が 10·N_C 1.87 をわずかに下回って外れた。

| 窓 | 断面 | ΔA | ΔB | ΔB′ | ΔC | ΔC′ | mB − mC | N_B | N_C |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| W1 | θ_r(40) | −0.056 % | +1.577 % | +1.439 % | +1.487 % | +1.398 % | +0.066 | 0.138 | 0.089 |
| W1 | θ_r(70) | +0.012 % | +1.229 % | +1.120 % | +1.136 % | +1.097 % | +0.058 | 0.109 | 0.040 |
| W1 | θ_r(94) | +0.000 % | +0.995 % | +0.908 % | +0.911 % | +0.892 % | +0.050 | 0.087 | 0.018 |
| W2 | θ_r(40) | −0.121 % | +1.883 % | +1.776 % | +1.832 % | +1.645 % | +0.091 | 0.106 | 0.187 |
| W2 | θ_r(70) | −0.056 % | +1.512 % | +1.426 % | +1.467 % | +1.333 % | +0.069 | 0.085 | 0.134 |
| W2 | θ_r(94) | −0.054 % | +1.267 % | +1.194 % | +1.217 % | +1.117 % | +0.064 | 0.073 | 0.101 |

- 副判定 Q_w: mB − mC は +0.33 (W1)・+0.17 (W2) で、N_B・N_C (0.25〜0.44) と同じ桁。
- **記録 (判定ではない)**: block の系の組み立てを double にしても、float の 2 本は B とほぼ同じように FP64 の参照から離れた (mC は mB の 94〜97 %)。この変更で登録の窓の偏りが抑えられる兆しは無い。発生源・必要性は言わない (§6.32 の範囲)。
- 次の段は諮ってから (残る候補: 残差の評価 (流束・勾配・リミタ)、平均流の状態の commit・保存の丸め、状態方程式の往復)。

### 6.34 事前登録: 平均流の状態を FP64 で積む float の軌跡 (qAccumulatorFP64 の感度試験) (2026-10-11、codex diagnose、実装の前)

- **問い**: §6.27 の float の偏りは、平均流の保存量の積み方 (commit のたびに float32 に丸める) を変えると抑えられるか。**感度試験であって、発生源の特定ではない** (codex diagnose 2026-10-11 Major 2: Qacc は 1 ULP 未満の増分の吸収だけでなく、切り上げ・切り捨てを含む累積の丸め全体を変える。要求の更新・EOS の再構成・残差は float のまま)。
- **記録 (登録外、2026-10-11 に読んだ)**: §6.26 の腕の commit の診断 (`FORGE_DIAG_COMMIT_LOSS=500`) で、float の `run_0488_tr_f32` は 500 step ごとの各時点で、内点 561,323 節点のうち要求した更新が丸ごと消えた節点が ρ 53,853〜81,879・ρu 23,896〜38,682・ρv 7,601〜15,843・ρE 30,538〜47,338 (k・ω・ρY も数千〜数万)。FP64 の `run_0487_tr_fp64` は全量で 0。更新量の合計に対する実際の比は float でも 0.9999 (消えた量は合計の 0.01 % 程度)。消えた節点の位置は見ていない。これが偏りの原因だとは言わない。
- **実装**: §5.1 #18 の契約 (codex Major 1 の修正を含む)。前提の単体試験が通るまで投入しない。
- **腕** (~~インスタンス B~~ **インスタンス A** (2026-10-11 run の前に変更: B は他の 2 セッションの run 3 本でメモリの空きが 7.4 GB だった。A と B は同じ型 g5.xlarge・同じ AMI。A・B・B′ は B で回した記録で、機械の違いも対照 F が拾う)、3 本同時、各 30,000 step・500 ごと): §6.26 と同じ共通の入力 `_tr/q32_init.h5` (sha256 `ba859e02…`)・格子・キー 1・実効設定・環境変数 (`FORGE_CUDA_BLOCKSIZE=128`・`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500`)。バイナリは #18 の commit を A でビルドした float (`~/forge-fgeom9-f32`、`qa_build.sh`)。基準の float (`~/forge-fgeom7-f32`) も A にある同じ commit のビルド (SASS の照合はこれと)。
  - E = `run_0592_qa_f32`、E′ = `run_0593_qa_f32b`: `time.deltaT.qAccumulatorFP64: 1` (B との設定の差はこの 1 語)。環境変数に `FORGE_QACC_AXISYM_DIAG=1` と採用の mask (`FORGE_DIAG_QACC_ADOPT`、step 1 と 30,000) を足す。
  - F = `run_0594_qa_off`: 同じ新しいバイナリで `qAccumulatorFP64: 0` (バイナリの変更の対照)。
  - 台本 `case/45.isobutane_m6_d155/qa.sh`、判定 `qa_an.py` → `_band_ab/cold_pair/qa_judge.json`。
- **前提** (総合の判定の前に確かめる。一つでも欠ければ**判定不能**。結果を見て対照を差し替えない、codex Major 3):
  1. 単体試験 `test_qacc_axis.cu` の全項目 PASS。新しいバイナリと §6.26 の float のバイナリで `enforceAxisSymmetry_d` 以外のカーネルの SASS が一致し、ホスト側の差 (初期化・引数・実効設定の出力) が #18 の範囲だけ。
     - **追記 (2026-10-11、実装の後・run の前)**: nvcc 12.0 は同じソース・同じフラグでも一部のカーネルで 2 通りの SASS を出す (implementer が `dependentVariables_d`・`cond_moment_update_limited_d` で再現。HEAD のソースを 3 回再コンパイルして 2 通り)。そこで、一致しない関数は、基準 (`~/forge-fgeom7-f32`) のソースのその `.cu` を同じフラグで 3 回再コンパイルした SASS のどれかと一致すれば「コード生成の揺れ」として一致とみなす。それでも外れれば判定不能。台本 `qa_build.sh` (ビルド・単体試験・照合)・`qa_sass_cmp.py` → `_band_ab/cold_pair/qa_pre.json`。揺れとみなしたカーネルの数値の差は対照 F (前提 4) が拾う。
  2. E・E′・F の全期間の残差と全 61 出力が有限、`check_convergence` が実行できる。初期化後の保存量 9 つ (res_0 の ρ..ρE・ρk・ρω・ρY0・ρY1) が B・E・E′・F でビット一致。
  3. E・E′ の両方で、採用の mask (step 1 と 30,000) が軸の節点・内点で 0、採用は等温壁の節点の ρE (と壁ピンの書く量) だけ。monitorInterval ごとの採用数が (等温壁の節点数 × monitorInterval) 以下。**2 時点の観測であり、全期間で軸・内点の採用が 0 だったとは言わない** (codex Minor 2)。格子の軸と壁の flag が重なる節点の数も記録する。
  4. F が B・B′ の範囲に入る: 両窓・3 断面の 6 組すべてで min(ΔB, ΔB′) − N_B ≤ ΔF ≤ max(ΔB, ΔB′) + N_B。
- **判定** (θ_r の 3 断面、W1 = 20,000〜25,000・W2 = 25,000〜30,000、Δ は §6.26 と同じ、D = (+1.892, +1.647, +1.431) %、k = 10。k は信頼係数でなく暫定の診断基準。§6.33 の僅差の未達を理由に下げない、codex Minor 1):
  - mB = mean(ΔB, ΔB′)、mE = mean(ΔE, ΔE′)、N_E = |ΔE − ΔE′|、N_BF = max(ΔB, ΔB′, ΔF) − min(ΔB, ΔB′, ΔF)。
  - **結果 S (抑制)**: max(|ΔE − ΔA|, |ΔE′ − ΔA|) + 10 N_E ≤ 0.1|D|、かつ sign(D)(mB − mE) ≥ 0.5|D|、かつ |mB − mE| > 10 max(N_BF, N_E)。
  - **結果 R (抑制不十分)**: E・E′ がそれぞれ A から D の向きへ 0.5|D| 以上、かつ |mE − ΔA| > 10 N_E。
  - 両窓・3 断面で揃わなければ判別不能。**延長・窓の選び直し・倍率の変更はしない**。Q_w は副判定として記録。
- **解釈の範囲**:
  - S: 「平均流の保存量の累積の方法を変えると、登録の窓の参照 (A) からの偏りを抑えられた」。吸収が原因とは言わない。A は保存済みの参照の軌跡で、再実行の差は測っていないので、A に近いことを「FP64 一般と同等」とは呼ばない。
  - R: 「この変更だけでは抑えられない」。残差の評価・要求の更新の生成・float のままの k/ω・ρY の更新のどれかは決まらない。R を根拠に左辺・増分の生成を候補から外さない。
  - 呼び出し側の「状態の保存か残差の評価か」の二択は採らない (codex)。
- **やらないこと**: commit の軸の分岐 (`update_d.cu:345` の nullptr) を有効にすること、軸ピンの直前に Qacc を Q へ合わせて残余を捨てること。
- **所要**: 実装・単体試験・ビルド 1〜2 時間、run 約 1 時間 (B で 3 本同時)、ディスク 約 10 GB (判定後に途中の場を消す)。

### 6.35 §6.34 の結果 (2026-10-11)

- run (AWS **A**): E `run_0592_qa_f32`・E′ `run_0593_qa_f32b` (`qAccumulatorFP64: 1`、`FORGE_QACC_AXISYM_DIAG=1`)・F `run_0594_qa_off` (同じバイナリで 0)。各 30,000 step、RUN_RC 0。バイナリ `~/forge-fgeom9-f32` (commit `4d45f5f6`、sha256 `a9ce00b1…`)。判定 `qa_an.py` → `_band_ab/cold_pair/qa_judge.json` (前提 1 の記録は `qa_pre.json`)。
- **前提はすべて満たした**:
  1. 単体試験 `test_qacc_axis`・`test_qacc_commit` は A でも PASS。SASS は 283 関数のうち 282 が基準 (`~/forge-fgeom7-f32`、sha256 `79321591…` = §6.26 の float と同じ) と一致し、違うのは `enforceAxisSymmetry_d` (シグネチャ) だけ。nvcc の揺れは出なかった。最初の照合は A に `cuobjdump` が無く空の SASS どうしを「不一致なし」と通していたので、空の入力を不合格にする修正を入れてやり直した (`qa_sass_cmp.py`、CUDA 13.2 の `cuobjdump`・`nvdisasm` を venv に入れた)。
  2. 3 本とも 61 出力・全期間で有限、`check_convergence` は実行でき OVERALL: CHECK FAILURES (未収束、この区間の判定に使わない)。res_0 の 9 保存量は B・E・E′・F でビット一致。設定の差は E・E′ が `qAccumulatorFP64` の 1 語、F は 0 行。
  3. 採用の mask (step 1・30,000、E・E′ とも) は等温壁の節点の ρE だけ (軸・内点・壁の外・ほかの量は 0)。等温壁 4,719 節点、軸と壁の重なり 0。1 step ごとの採用数は E 3,978〜4,550・E′ 3,890〜4,516 (上限 4,719 × 1 以内)。**2 時点の観測であり、全期間で軸・内点の採用が 0 だったとは言わない**。
  4. F は 6 組すべてで B・B′ の範囲に入った (W1 θ_r(94) は下限 0.821 に対し 0.823 で、ぎりぎり)。
- **判定 (事前登録のとおり): 判別不能**。結果 S は 6 条件中 0。結果 R は 6 条件中 0 で、外れた理由はすべて同じ: 「E・E′ がそれぞれ A から D の向きへ 0.5|D| 以上」は 6 組とも成立したが、「|mE − ΔA| > 10 N_E」が成立しない (E と E′ の差 N_E が 0.20〜0.47 ポイントと大きく、10 N_E が 2.0〜4.7 になる)。

| 窓 | 断面 | ΔA | ΔB | ΔB′ | ΔF | ΔE | ΔE′ | mB − mE | N_B | N_BF | N_E |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| W1 | θ_r(40) | −0.056 % | +1.577 % | +1.439 % | +1.325 % | +1.686 % | +1.302 % | +0.014 | 0.138 | 0.252 | 0.384 |
| W1 | θ_r(70) | +0.012 % | +1.229 % | +1.120 % | +1.025 % | +1.281 % | +1.026 % | +0.021 | 0.109 | 0.204 | 0.254 |
| W1 | θ_r(94) | +0.000 % | +0.995 % | +0.908 % | +0.823 % | +1.028 % | +0.832 % | +0.022 | 0.087 | 0.172 | 0.196 |
| W2 | θ_r(40) | −0.121 % | +1.883 % | +1.776 % | +1.713 % | +2.134 % | +1.663 % | −0.069 | 0.106 | 0.169 | 0.471 |
| W2 | θ_r(70) | −0.056 % | +1.512 % | +1.426 % | +1.363 % | +1.685 % | +1.327 % | −0.037 | 0.085 | 0.148 | 0.358 |
| W2 | θ_r(94) | −0.054 % | +1.267 % | +1.194 % | +1.130 % | +1.407 % | +1.107 % | −0.026 | 0.073 | 0.137 | 0.300 |

- 副判定 Q_w: 両窓で結果 R の条件が成立 (|mE − ΔA| 2.63・2.31 > 10 N_E 1.98・1.16)。
- **記録 (判定ではない、codex diagnose 2026-10-11 の範囲に直した)**: 登録の W1・W2 では、E・E′ の A からの差の平均は基準の float (B・B′) の 97.68〜103.56 % だった。大幅な抑制は観測されなかったが、事前登録の判定は S 0/6・R 0/6 で判別不能。**平均流の丸めの因果的な関与は除外できない** (qAccumulatorFP64 は累積の正本を double にするが、演算に渡す Q は float のままで、状態の表現・残差の入力の丸めまで除いた試験ではない)。有限期間の応答の記述であり、収束解の一致や固定点の偏りの確定ではない (E・E′・F とも `check_convergence` は CHECK FAILURES)。E と E′ の再実行の差 (N_E 0.20〜0.47) は B と B′ (0.07〜0.14) より大きい (N_E/N_B 2.27〜4.43)。F の合格範囲は B・B′ の最小〜最大の左右に N_B を足した許容範囲。§6.29・§6.33・§6.35 はいずれも総合判定が判別不能で、「どの介入も効かなかった」とは書かない (§6.29 は効果量の条件も一部未達)。
- 後始末: 3 本とも途中の場を消した (残したのは res_0・20,000・25,000・30,000、採用の mask、ログ・残差履歴は gzip、VERDICT)。系列は `qa_judge.json` にある。
- 次の段は諮ってから (エスカレーション条件 3・7)。

### 6.36 事前登録: FP64 のビルドで step の末尾に状態を float32 に丸める試験 (十分性の試験) (2026-10-11、codex diagnose、実装の前)

- **問い** (codex diagnose `notes/reviews/2026-10-11-qacc-result-next-step-diagnose.md`): FP64 のビルドで、step の末尾に保存量 9 つを float32 に丸めて戻すだけで、§6.27 の float と同じ向き・同じ規模の偏りが出るか。**十分性の試験**であり、「状態の保存が float の真因」「float の演算は無関係」のどちらも言わない。
  - 第 1 仮説 (確度 低): step の末尾の量子化だけでは観測の規模の偏りを再現せず、step の中の低精度の演算・途中の丸めとの組み合わせが要る。
  - 第 2 仮説 (未確認): 9 保存量の末尾の量子化だけでも同じ向き・同じ規模の偏りを生む。
- **実装**: §5.1 #19 (`FORGE_DIAG_QUANTIZE_STATE`)。
- **腕**: 同じ新しい FP64 のバイナリ (`~/forge-fgeom10-fp64`、#19 の commit を A でビルドし、B へは**同じファイルを写す** (nvcc のコード生成の揺れを持ち込まないため))。
  - 対照 A\* = 量子化 OFF × 5 (`run_0595_sq_off1`〜`run_0599_sq_off5`)、介入 G = 量子化 ON × 5 (`run_0620_sq_on1`〜`run_0624_sq_on5`)。各 30,000 step 固定・500 ごと。各反復は同じ起点から独立に開始する。
  - **反復を 5 にした理由** (codex の案は 3。結果を見る前に変更): B・B′・F の 3 本の float の再実行の標準偏差を G のばらつきと仮定して区間の半幅を見積もると、3 反復では G が float を完全に再現しても W1 の 3 断面で L < 0.5|D| (半幅 0.43〜0.64、効果 0.91〜1.50) になり、判定はほぼ確実に判別不能になる。5 反復なら半幅 0.14〜0.26 で 6 組とも L ≥ 0.5|D| になりうる。ただし G のばらつきが E・E′ 並みに大きければ (N_E 0.20〜0.47)、5 反復でも W1 で外れうる (見積もりの半幅 0.29〜0.58)。
  - 共通の入力 `_tr/q32_init.h5` (sha256 `ba859e02…`)・格子・キー 1・実効設定は §6.26 の FP64 の腕 (`run_0487`) と同じ (`implicitSolvePrecision` は 0 のまま。FP64 のビルドにも `ST=float` の経路があるが、既存の FP64 の基準の演算経路を維持する、codex Minor)。環境変数は §6.26 と同じ (`FORGE_CUDA_BLOCKSIZE=128`・`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500`) に、G だけ `FORGE_DIAG_QUANTIZE_STATE=1` を足す。
  - 実行環境: A に off1・on1・off2・on2・off3、B に on3・off4・on4・off5・on5 (A と B は同じ型 g5.xlarge・同じ AMI・同じ GPU で、各機械に両方の腕を置く)。ディスクのため、run の途中で step 1〜19,500 の res を消す (判定は step 0 と 20,000〜30,000 だけを使う)。旧 A (`run_0487`) は参考で、新しい A\* の反復は省略しない。
  - 台本 `sq.sh`、判定 `sq_an.py` → `_band_ab/cold_pair/sq_judge.json`。
- **前提** (一つでも欠ければ判定不能):
  1. 10 本とも全期間の残差と使う出力 (step 0・20,000〜30,000 の 22 出力) が有限、`check_convergence` が実行できる (VERDICT を記録)。`check_quasisteady` は θ_r の系列で `--drift 0.0001 --osc 0.0001` を記録する (判定に使わない)。
  2. res_0 の 9 保存量が 10 本でビット一致 (旧 A の res_0 との一致も記録)。設定は 10 本とも `run_0487` と同じ (差 0 行)。
  3. 介入の実効性: G の 5 本は、使う全出力で 9 保存量の全 CV が float32 で表せる値。自己点検 (step 1 と 30,000) で履歴の配列の表せない値が 0。A\* の 5 本は step 30,000 で表せない値が多数ある (割合を記録)。
  4. 新しい A\* が旧 A と同じくほぼ動かない: 5 本の Δ がすべての組で |Δ| ≤ 0.1|D|。
- **判定** (θ_r の 3 断面、W1 = 20,000〜25,000・W2 = 25,000〜30,000、Δ = (窓の平均 − 自分の step 0) / |A\*1 の step 0| × 100、D = (+1.892, +1.647, +1.431) %):
  - 独立な標本は **run ごとの窓の平均** (窓の中の 11 出力や W1・W2 を反復に数えない)。e = sign(D)·(mean(ΔG) − mean(ΔA\*))。腕ごとの標本分散から Welch の差の区間を作り、6 組に Bonferroni の補正をした全体 95 % の両側区間 [L, U] (各組 t の分位 1 − 0.05/12、自由度は Welch–Satterthwaite。両腕の分散が 0 なら区間は点)。独立性と run の平均の近似正規性を仮定する (5 反復で分布を精密に推定できるとは扱わない)。
  - **結果 A**: 6 組すべてで L ≥ 0.5|D| → 第 1 仮説を退け、「この量子化の介入だけで、観測の差の少なくとも半分の規模を生める」。必要性や float の真因は確定しない。
  - **結果 B**: 6 組すべてで [L, U] ⊂ [−0.1|D|, +0.1|D|] → 第 2 仮説を退け、次の対象は途中の丸め・演算との組み合わせ。
  - それ以外は判別不能。**有意差が出ないだけでは結果 B にしない**。倍率・窓を変えて救済しない。Q_w は同じ式の副判定として記録する。
- **やらないこと** (codex): §6.35 を結果 R に読み替えること、異なる介入の run をプールすること、最初から平均流だけ・k/ω/ρY だけの腕に分けること、qAccumulatorFP64 の削除・既定の変更。
- **所要**: 実装・ビルド 1〜2 時間、run は A と B で 5 本ずつ同時に約 2 時間 (FP64 の case/45 は RSS 約 2.1 GB/本、A の空き 14 GB・B の空き 12 GB)、ディスクは 1 本 約 3〜4 GB (途中の res を消す)。

### 6.37 事前登録: case/45 の陰解法の設定一式を入れた case/48 の float と FP64 の有限期間の再現判定 (2026-10-11、codex diagnose、実行前)

- **背景** (ユーザ決定 2026-10-11): case/45 は 1 本 1〜2 時間で切り分けが遅すぎる。偏りが再現する小さいケースを先に作る。**§6.36 (case/45 で 10 本) は保留** (#19 の実装は小さいケースでも使う)。
- **位置づけ** (codex diagnose `notes/reviews/2026-10-11-c48-line-implicit-reproducer-diagnose.md`): 「line implicit 単独の試験」ではなく、**case/45 の陰解法の設定一式 (`lineImplicit 1`・`lineDtDirectional 1`・上限 50・`implicitThermalJacobian 5`・`cfl 4`・`cfl_pseudo 4`) を 4 本に共通で固定した精度の A/B**。両腕の違いはビルドの精度だけ。case/48 の §6.21〜§6.24 は正式には判別不能で、「case/48 は偏りなし」とは言わない (言えるのは指定の窓の平均に % 級の差を観測しなかったことだけ)。case/45 と case/48 の違いは軸対称・TP 多成分・陰解法の設定のほか、収縮部・圧力勾配・格子・壁解像も含む。
  - 第 1 仮説 (確度 低): この設定一式が case/48 でも float の正方向の θ の偏りを増幅する十分条件になる。
  - 第 2 仮説 (確度 低、未確認): case/45 の収縮部・圧力勾配・軸対称・TP・格子との組み合わせが要る。
- **腕** (AWS A、4 本同時、各 48,000 step・2,000 ごと): §6.21 と同じバイナリ (`~/forge-fgeom7-{f32,fp64}`)・格子・入力 (`_prec/mesh_init.h5`。保存量 7 つはもともと float32 で表せる値 (丸めで変わる値 0 個を確認))・設定 (`run_0057`) に上の設定一式を足す。float ×2 = `case/48.flat_plate_cooled_m4/run_0061_li_f32a`・`run_0062_li_f32b`、FP64 ×2 = `run_0063_li_f64a`・`run_0064_li_f64b`。台本 `c48_li.sh`、判定 `tools/c48_li_judge.py` → `c48_li_judge.json`。
- **抽出**: `c48_prec_series.py --steps 0:48000:2000 --coords-from _prec/mesh_init.h5` (共通の double の座標、節点の対応の検査つき)、加えて縁の点 je を各 run の step 0 の値に固定した抽出 (`--fix-je`、`prec_series_fixje.csv`)。
- **前提** (欠ければ判定不能):
  1. 4 本とも RUN_RC 0・48,000 step、残差の全行・全列が有限、系列 25 点 (2 種) が有限、`check_convergence` が実行できる (VERDICT は記録)。
  2. 線: 4 本のログの `[lineImplicit] lines=… covered CVs=…` が同じで、被覆率 ≥ 95 % (評価断面の境界層が線に載ることの代用。直接の確認は無い)。方向別 dt のログは記録する。
  3. 初期化後: step 0 の θ(0.3/0.6/0.9) の精度間の相対差がすべて ≤ 0.0005 (超えたら更新による増幅の判定を保留)。res_0 の保存量 7 つの精度間の最大相対差は記録する (§6.31 で初期化の ρE の組み直しが精度で違った実績がある)。
- **判定** (主判定は θ(0.3)・θ(0.6)・θ(0.9)、窓 W1 = 24,000〜36,000・W2 = 36,000〜48,000。結果を見て窓・閾値を変えない):
  - Δ_r(s) = (θ_r(s) − θ_r(0)) / θ_FP64(0)。e = (float の 2 本・窓の全点の Δ の平均) − (FP64 の同じ平均) (符号つき。case/45 と同じ向きが正)。E = 各精度の窓の全点の、その平均からの最大偏差の和 (観測の範囲で、信頼区間ではない)。
  - **結果 A (再現)**: 同じ θ の断面で両窓とも e − E ≥ 0.005。その断面の窓で je が step 0 から変わった run があれば、je を固定した抽出でも A が保たれること。→ 「case/45 と同じ向きの 0.5 % 以上の精度依存の偏りを、小さいケースの有限期間で再現」。軸対称・TP が必須という説を退ける。真因が case/45 と同じとは確定しない。以後の切り分けの試験台にする。
  - **結果 B (再現せず)**: 3 断面・両窓とも |e| + E ≤ 0.0005 → この case・入力・48,000 step では設定一式だけで十分という仮説を退ける。「line 系は無関係」とは読まない。次の候補は **case/45 の物理・BC・陰解法の設定と近壁の解像を保った粗格子版** (既知の症状へ構成を近づける探索の順であって、原因の順位ではない)。
  - それ以外は判別不能。逆符号の大差や q_w・Cf・δ* だけの大差は記録し、A にも B にも読み替えない。
  - 「定常差」「収束解の偏り」とは呼ばない (有限期間の応答)。記録: `check_quasisteady` (窓 24,000〜48,000、θ・δ* は 1e-4、Cf・q_w・CD・HF は 2e-4)、§6.21 の旧判定器の結果 (`c48_li_judge_old.json`、判定には使わない)。
- **所要**: 4 本同時で 10〜30 分の見込み。

### 6.38 記録: float と FP64 の差が x のどこから育つか (既存の出力の解析、2026-10-11、判定ではない)

- **きっかけ**: ユーザの見立て「スロート上流の亜音速域で問題が起きているのでは。低マッハ向けの補正が足りないのでは」。新しい run は回さず、§6.27 の軌跡の run (AWS B の `case/45.isobutane_m6_d155/run_0487_tr_fp64`・`run_0488_tr_f32`・`run_0489_tr_f32b`、共通の Q32 の起点、500 step ごとの全出力) で、各節点の「自分の step 0 からの変化」の float − FP64 を A の step 0 の値で割った量を、step 1,000・2,000・5,000・10,000・20,000・30,000 で見た。道具 `c45_diffmap.py` (中間 `_band_ab/cold_pair/c45_diffmap.npz`、B)・`c45_diffmap_sum.py` (x ごとの要約 `c45_diffmap_sum.npz`)・`c45_diffmap_loc.py` (上位の節点と割合)。
- **観測** (B = `run_0488`。B′ も同じ傾向):
  - step 1,000 では、下の指定の時点・正規化・閾値で、入口〜収縮部 (x/r_t < −2) だけが閾値を超えた (下流も 0 ではなく閾値未満。例: 超音速部の中心の |ΔU|/U の最大 3e-6)。|ΔU|/U > 1e-3 の節点の割合は 入口 (x/r_t < −8) 23 %・収縮部 (−8〜−2) 7 %・スロート付近 (−2〜2) 0 %・超音速部 0 %。|ΔT|/T > 1e-4 は 28 %・10 %・0 %・0 %。中心の列の |ΔU|/U の最大は 入口・収縮部 2.1e-2、スロート付近 1.2e-4、超音速部 3e-6。
  - step 10,000 では差が下流へ広がる (|ΔT|/T > 1e-4 の割合: 収縮部 71 %・スロート付近 75 %・超音速部 47 %)。step 30,000 では超音速部も 80 %。θ_r の差も step 1,000 では x/r_t ≤ −4 にだけあり (x/r_t −8 で −1.2 %、−2 より下流は 0.001 % 未満)、step 30,000 では x/r_t 40〜94 で +1.5〜2.1 %。
  - 差の最も大きい所は x/r_t ≈ −7・壁に近い列 (j 100〜106、壁は j 120)、A の step 0 のマッハ数 0.005〜0.007。ΔT/T は step 1,000 で 3.3e-3、10,000 で (x/r_t −6.9・j 56、M 0.015) 0.20、30,000 で 0.37。|ΔU|/U の上位は壁の 1 つ内側 (j 119、M 0.0004〜0.0009) で 0.1 (step 1,000) → 1.5 (step 30,000)。
  - A の step 0 のマッハ数: 入口 x/r_t −12 で軸 0.016、収縮部 −6 で 0.036、−4 で 0.076、スロート x/r_t 0 で 0.91。
- **訂正 (codex diagnose 2026-10-11 `notes/reviews/2026-10-11-lowmach-origin-next-test-diagnose.md`、Major 1・2・5)**:
  - 「入口で生まれて下流へ運ばれる」は観測を超えている。書いてよいのは「指定の保存時点・正規化・閾値では、上流域で先に閾値を超え、後に下流域でも超えた」まで。局所の初期値で割るので低速域の相対差は大きく見えやすく、定常の局所擬似時間の step の先後だけでは物理的な移流の経路も決まらない。原因・伝播の経路は仮説として分ける。
  - 「ρE は生成エンタルピーを含む」は誤り。case/45 は `thermoHrefTemp: 298.15` で、種ごとの基準エンタルピーを除く (`thermo_d.cu:91`)。ρE は x/r_t < −2 で −5.4e6〜1.37e7 J/m³ (冷却壁の近くで負)。
  - **場の値での見積もり** (`run_0487` の step 0、組成固定の局所線形化で δp ≈ (γ−1)·ulp(ρE)/2、γ−1 = 0.274 の代表値): 丸めの圧力換算の中央値は 入口・中心 0.14 Pa・入口の壁側 0.03 Pa・差の最大の所 (x/r_t ≈ −7、j 100〜106) 0.07 Pa に対し、動圧 ½ρ|U|² の中央値は 559 Pa・63 Pa・180 Pa。比 δp/q は中央値 2e-5〜7e-4、最大 0.15 (壁のすぐ横の最も遅い節点)。**「ρE の刻みが動圧と同じ桁」は壁のすぐ横の節点を除いて当たらない** (前の版の代表値の見積もりは誤り)。隣り合う節点の間の圧力差は動圧より小さいので、これで仮説が退くわけでもない。
  - 「低マッハ前処理では直らない」とは言えない (前処理は失われたビットを戻さないが、圧力差に対する流束の応答や増幅を変えうる。`lowMachPrecond` 1/2 は質量流束の圧力差の項を変える。`lineImplicit` と `lowMachPrecond ≥ 2` は併用できない)。効果は未検証。
- **仮説 (未検証、原因とは書かない)**: 入口〜収縮部の低マッハ域は小さな乱れを増幅しやすく (float 同士でも差が育つ)、float では状態の丸め (ρE を含む) がその乱れの源になる。ρE の丸めだけで十分かは §6.40 で確かめる。直し方の候補に順位は付けない (codex)。

- 次は codex に諮り、この仮説を安く確かめる試験を決める (候補: 収縮部〜スロート直後だけの小さい格子での精度の A/B、FP64 で ρE だけを float32 に丸める試験)。

### 6.39 §6.37 の結果 (2026-10-11)

- run (AWS A): float `case/48.flat_plate_cooled_m4/run_0061_li_f32a`・`run_0062_li_f32b`、FP64 `run_0063_li_f64a`・`run_0064_li_f64b` (各 48,000 step、RUN_RC 0、系列の抽出 2 種とも rc 0)。判定 `tools/c48_li_judge.py` → `case/48.flat_plate_cooled_m4/c48_li_judge.json`。
- **前提はすべて満たした**: 4 本とも残差・系列 25 点が有限、`check_convergence` は実行でき OVERALL: CHECK FAILURES (この判定には使わない)。線は 4 本とも 1,001 本・被覆 102,102/106,080 (96.2 %)・最大長 102。step 0 の θ の精度間の相対差は 1.6e-7〜7.6e-7 (≤ 0.0005)。res_0 の保存量の精度間の最大相対差は ρ・ρu 0、ρE 1.5e-7、ρk 5.1e-8、ρω 2.1e-7 (初期化の組み直しの差)。
- **判定 (事前登録のとおり): 結果 B (再現せず)**。3 断面・両窓とも |e| + E ≤ 0.0005。|e| は 0〜1.3e-5、|e| + E は 1.1e-5〜9.2e-5 (θ(0.3) が最大)。

| 量 | 窓 | e (float − FP64) | E | \|e\| + E |
| --- | --- | --- | --- | --- |
| θ(0.3) | W1 / W2 | −1.3e-5 / +2e-6 | 8.0e-5 / 8.7e-5 | 9.2e-5 / 9.0e-5 |
| θ(0.6) | W1 / W2 | +5e-6 / +4e-6 | 2.6e-5 / 2.9e-5 | 3.1e-5 / 3.3e-5 |
| θ(0.9) | W1 / W2 | 0 / +1e-6 | 1.1e-5 / 1.4e-5 | 1.1e-5 / 1.4e-5 |

- 記録: θ 以外 (δ*・Cf・q_w・CD・HF) にも |e| − E ≥ 0.005 の大差は無い (最大は δ*(0.6) の |e| + E 8.3e-4、E が大きい)。逆符号の大差も無い。θ(0.6) の窓で je が step 0 から変わったのは float の 2 本 (結果 B なので je 固定の確認は要らない)。`check_quasisteady` (窓 24,000〜48,000、θ、1e-4): float の θ(0.3) が OSCILLATING (±3.2e-8 = 0.008 %)、ほかは STEADY (記録 `/home/sano/.claude/jobs/e254baf2/tmp/c48_li_qs.txt` は一時ファイル、判定には使わない)。方向別 dt の実際の伸びはログに出ておらず確かめていない (設定の読み込みの行だけ)。§6.21 の旧判定器は系列の点数 (25 点、step 0 を含む) が想定の 24 点と違って判定不能 (記録用で、判定には使わない)。
- **書いてよい範囲** (codex の §6.37 の範囲): この case・入力・48,000 step では、case/45 の陰解法の設定一式だけでは float の θ の偏りは再現しない。「line 系は無関係」とは読まない。§6.38 の観測 (差は入口〜収縮部の低マッハ域で生まれる) と合わせると、case/48 に無い要素 (低マッハ域・収縮・軸対称・TP) の方へ探索を進める。

### 6.40 事前登録: FP64 のビルドで ρE だけを step 末尾に float32 に丸める短い試験 (入口〜収縮部の初期の応答、2026-10-11、codex diagnose、実行前)

- **問い** (codex diagnose `notes/reviews/2026-10-11-lowmach-origin-next-test-diagnose.md`): 同じ FP64 のバイナリで step 末尾に ρE だけを `double(float(ρE))` にすると、§6.38 で観測した入口〜収縮部の初期の温度差の主な成分が出るか。**ρE 単独の末尾の量子化の十分性**に問いを限る (全 CV の末尾の量子化で、上流域だけが発生源だという試験ではない。陽性でも float の実際の原因・必要性とは呼ばない)。
  - 第 1 仮説 (確度 低): ρE の末尾の量子化だけで、入口〜収縮部の初期の温度差の主な成分を生じさせられる。
  - 第 2 仮説 (未確認): ほかの保存量・原始量・物性・残差などの途中の丸めとの組み合わせが要る。
  - 第 3 仮説 (未確認): 上流域は微小な摂動を増幅し、観測の一部は初期化の差や再実行の揺れを表す。
- **腕** (~~AWS B~~ **AWS A、1 本ずつ順に** (2026-10-11 run の前に変更: B は他の 2 セッションの run 4 本で空きが 5 GB だった。A でも forge-25 の 2 本と共存するので、各 run の前に available ≥ 4.5 GB を確かめて 1 本ずつ C1・E1・F1・C2・… の順。FP64 のバイナリは B でビルドしたものを A へそのまま写した (sha256 `74954e7b…`)。参照の A・B・B′ の step 0・1,000・1,500・2,000 の場も B から A へ写した)、各 2,000 step 固定・step 0 と 500 ごとに出力、共通の入力 `_tr/q32_init.h5`・§6.26 と同じ設定・環境変数 (`FORGE_CUDA_BLOCKSIZE=128`・`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500`)、全格子):
  - 対照 C = FP64・量子化 OFF × 3 (`run_0595_rq_c1`〜`run_0597_rq_c3`)
  - 介入 E = 同じ FP64 のバイナリで `FORGE_DIAG_QUANTIZE_STATE=1`・`FORGE_DIAG_QUANTIZE_STATE_ONLY=roe`・`FORGE_DIAG_QUANTIZE_STATE_CHECK=1,2000` × 3 (`run_0598_rq_e1`・`run_0599_rq_e2`・`run_0620_rq_e3`)
  - 陽性対照 F = float (`~/forge-fgeom7-f32`、§6.26 の B・B′ と同じバイナリ) × 3 (`run_0621_rq_f1`〜`run_0623_rq_f3`)。短い期間・この指標で既知の症状が出るかの確認だけに使い、因果を判別する比較は C 対 E の一つ。
  - C・E のバイナリは #20 の commit を B でビルドした FP64 (`~/forge-fgeom11-fp64`)。1 回に C_i・E_i・F_i の 3 本を同時に回し、3 回に分ける。§6.36 の run 番号の割り当て (`run_0595〜0599`・`run_0620〜0624`) は §6.36 の保留中はこちらに使い、§6.36 を再開するときに取り直す。
  - 台本 `rq_build.sh`・`rq.sh`、判定 `rq_an.py` → `_band_ab/cold_pair/rq_judge.json`。
- **参照の型** (結果を見る前に固定、新しい run で選び直さない): 評価域は x/r_t < −2 の内点 (軸 j = 0・壁 j = 120・入口 i = 0 を除く。境界の節点は別に記録)。共通の double の座標 (`nozzle.h5` の `/MESH/COORD`) で各 run の節点の対応を確かめる。既存の A (`run_0487`)・B (`run_0488`)・B′ (`run_0489`) から、z_{r,s} = (T_{r,s} − T_{r,0}) / T_{A,0}、v_s = mean(z_{B,s}, z_{B′,s}) − z_{A,s} (s = 1,000・1,500・2,000)。
- **指標**: 各 run の a_r = mean_s ⟨v_s, z_{r,s}⟩ / ⟨v_s, v_s⟩ (s = 1,000・1,500・2,000)。この定義で既存の float − FP64 の平均の応答が 1。独立な標本は **run ごとの a_r** (節点数・時点を反復に数えない)。e_E = mean(a_E) − mean(a_C)、e_F = mean(a_F) − mean(a_C)。腕ごとの標本分散から Welch の区間、2 つの比較で全体 95 % に補正 (各比較 t の分位 1 − 0.05/4、自由度 Welch–Satterthwaite)。
- **判定** (0.5・0.1 は診断用の暫定の閾値):
  - **結果 A**: e_F の区間の下限 ≥ 0.5、かつ e_E の区間の下限 ≥ 0.5 → 「ρE の末尾の量子化だけで、登録した初期の温度差の成分を少なくとも半分の規模で生じさせられる」。第 2 仮説の「その成分には他の丸めも必須」を退ける。
  - **結果 B**: e_F の区間の下限 ≥ 0.5、かつ e_E の区間全体が [−0.1, +0.1] の中 → 第 1 仮説を退ける。
  - 中間・逆符号・陽性対照の不成立は判別不能。3 反復で検出力が足りる保証は無く、**有意差が無いことを B と扱わない**。この判定から下流の θ_r・定常解・全場の再現は主張しない。
- **前提** (欠ければ判定不能): C・E の res_0 の保存量 9 つがビット一致。E の自己点検 (step 1・2,000) で ρE とその履歴 (roeN・roeM) の表せない値が 0、選ばなかった量の丸め前後の変化が 0。9 本とも全残差・保存した場が有限。共通の座標と節点の対応。F の初期化の差 (res_0 の ρE の C との差) は記録する。`check_convergence`・`check_quasisteady` の結果は記録する (有限期間の応答の試験なので、収束・準定常の合否とは分ける)。絶対値の分布・速度・圧力は副記録。
- **やらないこと** (codex): 4 腕を一度に投入して「離れなければ棄却」と読むこと、2,000 step の上流の応答を下流の θ_r の原因の確定に使うこと、直し方の候補に順位を付けること。
- **所要**: 実装・ビルド 30 分〜1 時間、run は 3 回 × 数分〜10 分。

### 6.41 §6.40 の結果 (2026-10-11)

- run (AWS A、1 本ずつ): C `run_0595_rq_c1`〜`run_0597_rq_c3` (FP64 OFF)、E `run_0598_rq_e1`・`run_0599_rq_e2`・`run_0620_rq_e3` (FP64 + ρE だけ量子化)、F `run_0621_rq_f1`〜`run_0623_rq_f3` (float)。各 2,000 step、9 本とも RUN_RC 0。判定 `rq_an.py` → `_band_ab/cold_pair/rq_judge.json`。
- **前提** (`check_quasisteady` の記録は判定器が取っていなかったので後から補った。下の §6.42 の記録を見よ): C・E の res_0 の 9 保存量がビット一致、E の 3 本とも自己点検 (step 1・2,000) で ρE と roeN・roeM の表せない値 0・選ばなかった 8 量の丸め前後の変化 0、9 本の全残差・保存した場が有限、座標が一致。B での単体試験は PASS (`rq_build.sh` が PASS の行の書式を読み違えて False と書いたので、`rq_pre.json` を手で訂正した記録あり)。
- **判定 (事前登録のとおり): 結果 B** — ρE 単独の末尾の量子化では、登録した初期の温度差の成分は出ない (第 1 仮説を退ける)。e_E = +0.006、区間 [−0.007, +0.018] (a_E 0.218〜0.225、a_C 0.2164〜0.2166)。陽性対照は成立: e_F = +0.91、区間 [+0.57, +1.25] (a_F 1.04〜1.23)。
- **記録 (物理量、`rq_dT.py`、判定ではない)**: step 2,000 の温度を FP64 の C1・C2 の平均と比べる (x/r_t < −2 の内点 83,310 節点):

| 腕 | RMS ΔT | 最大 \|ΔT\| | \|ΔT\| > 1 K の節点 |
| --- | --- | --- | --- |
| C (FP64 の再実行) | 0.002〜0.003 K | 0.05〜0.08 K | 0 % |
| E (FP64 + ρE の量子化) | 0.16〜0.19 K | 4.4〜5.5 K | 0.56〜0.58 % |
| F (float) | 3.7〜4.5 K | 119〜134 K | 5.4〜5.7 % |
| 参考: 既存の A / B・B′ | 0.006 K / 4.3 K | 0.16 K / 126 K | 0 % / 5.5〜5.7 % |

- **書いてよい範囲**: ρE の step 末尾の量子化は、FP64 の再実行の約 60 倍 (RMS) の温度差を生むが、float の差の RMS の約 4 % にとどまる。float の入口〜収縮部の初期の温度差の主な成分は、ρE 単独の量子化では生じない。他の保存量の丸め、原始量・物性・残差の float の演算のどれが要るかは決まらない (codex の §6.40 の範囲)。
- 入口〜収縮部の静圧 (`run_0487` の step 0): 入口 x/r_t < −10 で 5.4989e6〜5.4995e6 Pa (幅 約 550 Pa)、収縮部で 5.317e6〜5.4995e6 Pa。float の 5.5e6 の刻みは 0.5 Pa。`space.pRef` (対流流束の圧力項を (p̃ − pRef)·S で組む基準静圧、既定 0 = オフ) は case/45 では設定していない。

### 6.42 事前登録: FP64 のビルドで保存量 9 つを全部 step 末尾に float32 に丸める短い試験 (G9 対 C、2026-10-11、codex diagnose、実行前)

- **codex diagnose** (`notes/reviews/2026-10-11-pref-and-full-quantize-diagnose.md`) の採否: 全件採用。§6.41 から「主な成分は残差の演算側」とは書かない (未試験の 8 保存量と相互作用が残る、F と C は初期化後の保存量も違う、FP64 のビルドにも `ST=float` の経路がある)。**pRef の試験はこの登録から外す** (内部面・境界面・hoop を同時に変えるので通り道を分けない、Cp と C の同等性に事前の許容幅が要る、上限 ≤ 0.5 の式は大きな逆符号まで「抑制」に数える穴があった)。pRef の誤差の見積もり (「数百 Pa」「2〜数倍」) は入口付近だけ・未測定として取り下げる (p − pRef は収縮部全体で約 −182,000〜500 Pa)。
- **問い**: 既存の FP64 の演算経路に、保存量 9 つ (ρ・ρu・ρv・ρw・ρE・ρk・ρω・ρY0・ρY1) の step 末尾の量子化を加えるだけで、§6.40 で登録した入口〜収縮部の初期の温度差の成分を生じさせられるか (十分性)。
  - 第 1 仮説 (確度 低): 再現できない。第 2 仮説 (未確認): 他の保存量の丸めや相互作用を加えると再現できる。第 3 仮説 (未確認): float 固有の初期化の差の増幅が応答の一部を担う (C 対 F の比較だけでは除外されない)。
- **腕** (AWS A、1 本ずつ、各 2,000 step・500 ごと): G9 = `case/45.isobutane_m6_d155/run_0706_rp_g9a`・`run_0707_rp_g9b`・`run_0708_rp_g9c`。§6.40 の C と同じ FP64 のバイナリ (`~/forge-fgeom11-fp64`、sha256 `74954e7b…`)・入力 (`_tr/q32_init.h5`)・設定・環境変数に `FORGE_DIAG_QUANTIZE_STATE=1`・`FORGE_DIAG_QUANTIZE_STATE_CHECK=1,2000` を足す (`FORGE_DIAG_QUANTIZE_STATE_ONLY` は付けない = 全 9 本)。`pRef` は変えない。C は既存の `run_0595〜0597` を再利用、陽性対照 F は既存の `run_0621〜0623`。台本 `rp.sh` (G9 だけを回す形にする)、判定 `rg_an.py` → `_band_ab/cold_pair/rg_judge.json`。run 番号 `run_0700〜0705` は pRef の試験のために取っておく。
- **C の再利用の条件** (欠ければ判定不能): G9 と C でバイナリの sha256・入力の sha256・実効設定 (run_0595 との差 0 行)・res_0 の 9 保存量がビット一致。
- **G9 の自己点検** (量の名前まで検査し、「合計 0」の文字列だけで通さない。許す場所・時点は結果を見る前にこう固定する):
  - (a) step 1・2,000 の `updateVariablesOuter` の直後 (量子化の直後): Q・Q_N・Q_M の 9 量すべて (ro, roUx, roUy, roUz, roe, roK, roOmega, roY0, roY1 とその N・M) で表せない値 0。
  - (b) step 1 の次の step の化学種の更新の頭 (`speciesUpdateOuter` の直後) の roY{0,1}N・roY{0,1}M: 全 CV で 0、または 0 でないなら「入口のピン (scalarDirichletPin 1) を除く内点」で 0 (入口のピンの節点と ghost は step の中で double で書かれるので許す。§5.1 #19 の未決をこの形で決める)。step 2,000 の (b) は次の step が無いので無い。
- **指標・判定** (§6.40 の定義を固定): 評価域・参照の型 v_s・時点 1,000/1,500/2,000・run ごとの a。e_G9 = mean(a_G9) − mean(a_C)、e_F = mean(a_F) − mean(a_C) の 2 つの比較、Welch の区間、各比較の t の分位 1 − 0.05/4 (この 2 比較への補正で、探索全体の 95 % の保証ではない)。
  - **結果 A**: L_F ≥ 0.5 かつ L_G9 ≥ 0.5 → 第 1 仮説を退ける。「この介入だけで登録した成分の少なくとも半分を生じさせられる」まで。
  - **結果 B**: L_F ≥ 0.5 かつ [L_G9, U_G9] ⊂ [−0.1, 0.1] → 第 2 仮説を退ける。「9 保存量の末尾の量子化だけでは不十分」まで。残差の演算への局在化はしない。
  - それ以外は判別不能。延長・閾値の変更で救済しない。
- **記録**: 温度の K の表 (`rq_dT.py` に G9 を足す)、`check_quasisteady` の結果 (C・E・F・G9 の全 12 本、§6.41 の記録の不足を補う。有限期間の試験なので合格条件にしない)、`check_convergence`。収束判定は CHECK FAILURES で、定常解・下流の θ_r・全場の一致・原因の確定へは広げない。
- **所要**: 3 本 × 約 3.5 分。

### 6.43 §6.42 の結果 (2026-10-11)

- run (AWS A、1 本ずつ): G9 `case/45.isobutane_m6_d155/run_0706_rp_g9a`・`run_0707_rp_g9b`・`run_0708_rp_g9c` (各 2,000 step、RUN_RC 0)。判定 `rg_an.py` → `_band_ab/cold_pair/rg_judge.json`。
- **前提はすべて満たした**: G9 と C の FP64 のバイナリ (sha256 `74954e7b…`)・入力・設定 (run_0595 との差 0 行)・res_0 の 9 保存量が一致。自己点検 (a) は step 1・2,000 とも Q・Q_N・Q_M の全量で表せない値 0 (名前まで検査)。(b) は全 CV で roY{0,1}{N,M} が 244 節点ずつ表せない値だったが、入口のピン (121 節点) と ghost で、入口のピンを除く内点では 0 (事前に許した形)。
- **判定 (事前登録のとおり): 結果 B** — 9 保存量の step 末尾の量子化だけでは不十分 (第 2 仮説を退ける)。e_G9 = −0.013、区間 [−0.024, −0.001] (a_G9 0.201〜0.206、a_C 0.2164〜0.2166)。陽性対照 e_F = +0.91 [+0.57, +1.25] (§6.41 と同じ値の再掲)。残差の演算への局在化はしない (codex)。
- **記録 (物理量、`rq_dT.py`)**: step 2,000 の温度の C1・C2 の平均との差 (x/r_t < −2 の内点): G9 は RMS 0.28〜0.31 K・最大 5.4〜6.5 K・1 K を超える節点 1.2〜1.3 % (E は 0.16〜0.19 K・4.4〜5.5 K・0.6 %、float は 3.7〜4.5 K・119〜134 K・5.4〜5.7 %)。9 保存量を全部丸めても温度差の RMS は float の約 7 %。
- **記録 (`check_quasisteady`、§6.41 の不足を補う、`qs_record.py` → `_band_ab/cold_pair/qs_record.json`)**: C・E・F・G9 の 12 本とも OVERALL: NOT ALL STEADY (2,000 step の有限期間の試験なので合格条件にしない)。
- 書いてよい範囲: 保存量 9 つの step 末尾の float32 への量子化 (FP64 の演算経路のまま) では、入口〜収縮部の初期の温度差の主な成分は生じない。step の中の何 (原始量・物性の float での格納、流束・勾配・リミタ・粘性・ソースの演算、初期化の差の増幅) が要るかは決まらない。
- 次は諮ってから (候補: pRef の試験 (別に登録)、FP64 で原始量 (P・T・速度) を EOS の直後に float32 に丸める試験)。

### 6.44 事前登録: pRef (対流流束の圧力の基準値) の短い試験 (2026-10-11、codex diagnose、実行前)

- **codex diagnose** (`notes/reviews/2026-10-11-after-g9-next-diagnose.md`) の採否: 全件採用。案 A (pRef) を先に、案 B (原始量の量子化) は今回作らない (挿入位置の後で境界条件・`turbulent_viscosity` が上書きするので一括の試験にならない、G9 対 Gp は 2 因子が同時に変わる)。Cp と C の a が近いことを作用素の同等の根拠にしない (a は一方向の射影)。
- **問い**: `space.pRef: 0 → 5,499,000 Pa` の 1 つの変更だけで、§6.40 で登録した入口〜収縮部の初期の温度差の成分を十分に抑えられるか。第 1 仮説 (確度 低): 抑えられる。第 2 仮説 (未確認): 圧力などの格納・EOS・再構成・面の値の計算の丸めが重要で、後段の pRef の差し引きでは抑えられない。第 3 仮説 (未確認): ほかの演算・幾何の丸めや初期化の差との相互作用が要る。
- **腕** (AWS A、1 本ずつ、各 2,000 step・500 ごと、§6.40 と同じバイナリ・入力・設定・環境変数に `space.pRef: 5499000.0` だけを足す): Cp = FP64 × 3 (`case/45.isobutane_m6_d155/run_0700_rp_cp1`〜`run_0702_rp_cp3`、1 本目で `FORGE_DIAG_HOOP_CLOSURE` の閉性の診断も取る (読むだけ))、Fp = float × 3 (`run_0703_rp_fp1`〜`run_0705_rp_fp3`)。既存の C (`run_0595〜0597`)・F (`run_0621〜0623`) を、各精度でバイナリ・入力・環境・初期化後の 9 保存量が一致し、設定の差が `space.pRef` だけのときに限り再利用する。**Cp を先に回し、前段の関門を通ったときだけ Fp を回す**。台本 `rpref.sh`、関門 `cp_gate.py` → `_band_ab/cold_pair/cp_gate.json`、判定 `rpref_an.py` → `rpref_judge.json`。
- **前段の関門 (FP64 の同等性、外れたら前提不成立。幅を広げて救済しない。値は診断用の提案の許容幅で、実証済みの精度ではない)**:
  1. **作用素 (幾何)**: pRef による残差の変化は ΔR = pRef·(Σ±S_x, Σ±S_y − A_planar) (case/45 は `hoopAreaFromClosure: 0`、hoop は A_planar)。Cp1 の閉性の診断 (FP64、デバイスの最終の面ベクトル、境界の半割面を含む) から各実 CV で g = √((Σ±S_x)² + (Σ±S_y − A_planar)²) / (Σ|S_x| + Σ|S_y| + A_planar) を作り、全 CV で g ≤ 1e-12。軸の節点では hoop を当てず半径方向の行は 0 にする (拘束行) ので、そこは幾何の閉性として見る。
  2. **有限期間**: mean(a_Cp) − mean(a_C) の 95 % Welch 区間全体 ⊂ [−0.002, +0.002] (既存の C の a の最大差 0.000174 の 10 倍程度)。step 500・1,000・1,500・2,000 の Cp_i − C_j の全 9 組で、温度差 RMS ≤ 0.05 K かつ最大 |ΔT| ≤ 1 K を、上流域 (x/r_t < −2 の内点) と全域で別々に満たす。
- **主判定** (評価域・参照の型・時点・a は §6.40 のまま): e_Fp = mean(a_Fp) − mean(a_Cp)、d = [mean(a_F) − mean(a_C)] − e_Fp (4 腕の対比、SE² = Σ v_g、ν = (Σ v_g)² / Σ v_g²/(n_g − 1))。新しい 2 比較 (e_Fp と d) に各 t の分位 1 − 0.05/4。既存の F の陽性は過去の確認として扱い、探索全体の 95 % の保証とは混同しない。
  - **結果 A**: CI(e_Fp) ⊂ [−0.5, +0.5] かつ L_d ≥ 0.5 → この pRef の変更で、参照の尺度で 0.5 以上減らし、残りも ±0.5 以内。第 1 仮説を支持。「残差の加算が真因」とは言わない。
  - **結果 B**: L_Fp > 0.5、U_Fp < −0.5、または U_d < 0.5 → この変更だけでは抑制の目標を満たさない (逆符号の大差は抑制と呼ばない)。「圧力の丸めは無関係」とは言わない。
  - それ以外は判別不能。延長・窓の変更はしない。
- **記録**: 温度の K の表 (`rq_dT.py` に Cp・Fp を足す)、各 run の `check_convergence`、実効設定・入力・座標の対応。既定値の変更・本番への採用は、この試験だけでは決めない。
- **所要**: 6 本 × 2〜3.5 分。

### 6.45 §6.44 の結果 (2026-10-11)

- run (AWS A、1 本ずつ): Cp `case/45.isobutane_m6_d155/run_0700_rp_cp1`〜`run_0702_rp_cp3` (FP64 + `space.pRef: 5499000.0`)、Fp `run_0703_rp_fp1`〜`run_0705_rp_fp3` (float + 同じ pRef)。各 2,000 step、6 本とも RUN_RC 0、`check_convergence` は CHECK FAILURES (有限期間の試験なので判定に使わない)。関門 `cp_gate.py` → `_band_ab/cold_pair/cp_gate.json`、判定 `rpref_an.py` → `rpref_judge.json`。
- **前段の関門は通った**: 閉性 g の最大 2.6e-16 (99 % 点 1.2e-16、1e-12 を超える CV 0 / 570,999)。mean(a_Cp) − mean(a_C) = −0.000002、95 % 区間 [−0.00037, +0.00036]。Cp − C の温度差は全 9 組・全時点で RMS の最大 0.0054 K (上流域)・0.0020 K (全域)、|ΔT| の最大 0.137 K。前提 (設定の差は `space.pRef` だけ、バイナリ・入力が §6.40 と同じ、res_0 の 9 保存量の一致、有限) もすべて満たした。
- **判定 (事前登録のとおり): 判別不能**。e_Fp = +0.092、区間 [+0.046, +0.138] (⊂ [−0.5, +0.5] は成立)。d = e_F − e_Fp = +0.821、区間 [+0.494, +1.147] で、**L_d ≥ 0.5 に 0.006 届かない** (pRef なしの float の 3 本の a が 1.04〜1.23 とばらつき、SE が大きい)。結果 B の条件 (L_Fp > 0.5・U_Fp < −0.5・U_d < 0.5) にも当たらない。閾値・窓の変更で救済しない。a: C 0.2164〜0.2166、Cp 0.2163〜0.2166、F 1.04〜1.23、Fp 0.294〜0.317。
- **記録 (物理量、`rq_dT.py`、判定ではない)**: step 2,000 の温度の FP64 (C1・C2 の平均) との差 (x/r_t < −2 の内点):

| 腕 | RMS ΔT | 最大 \|ΔT\| | \|ΔT\| > 1 K の節点 | 差の最大の所 (x/r_t≈−7, j 100〜106) の平均 ΔT |
| --- | --- | --- | --- | --- |
| Cp (FP64 + pRef) | 0.004 K | 0.09 K | 0 % | 0.000 K |
| F (float) | 3.7〜4.5 K | 119〜134 K | 5.4〜5.7 % | −0.91〜−1.04 K |
| Fp (float + pRef) | 1.5〜1.6 K | 37〜41 K | 3.3〜3.6 % | −0.005〜−0.007 K |

- **書いてよい範囲**: 正式には判別不能。記録として、この pRef の変更で、登録した型への射影 (e) は 0.91 から 0.09 へ下がり、入口〜収縮部の温度差は RMS で約 6 割・最大で約 7 割減った。~~型と直交する差 (RMS 約 1.5 K) は残る~~ **訂正 (codex diagnose `notes/reviews/2026-10-11-pref-result-long-run-diagnose.md` Major): 温度差の RMS は 1.5〜1.6 K 残る。直交成分の大きさは未測定** (a は正規化した温度変化の射影係数の 3 時点の平均で、step 2,000 の絶対温度差の RMS と同一視できない)。「残差の加算 (圧力流束の足し込みの桁落ち) が真因」とは言わない (pRef は内部面・境界面・hoop を同時に変える)。下流の θ_r・Q_w への効果は未測定。
- 次は諮ってから (解釈の確定と、下流の θ_r を見る 30,000 step の試験の組み方) → §6.46。
- **判定器の穴 (codex Major)**: `rpref_an.py` の res_0 の照合は、存在する保存量だけを比べていて 9 量の欠損を確実には拒否しない。入力のハッシュ・座標の対応の証跡も判定記録に無い。§6.46 の判定器では 9 キーの完全一致を先に検査し、照合結果を記録に残す。

### 6.46 事前登録: pRef の長い試験 (下流の θ_r・Q_w、30,000 step、2026-10-11、codex diagnose、実行前)

- **codex diagnose** (`notes/reviews/2026-10-11-pref-result-long-run-diagnose.md`) の採否: 全件採用。§6.45 は判別不能のまま (L_d 0.494 を丸めて合格にしない)。短い試験の再試行 (反復を増やす) は挟まない。旧 A・B・B′ と新しい run を比べるには同じ実バイナリと陽性対照が要る (別の試験の無介入の F は W1・θ_r(40) で旧 B・B′ の平均より約 0.18 ポイント低かった)。「2 本とも満たす」は再実行の誤差への保証ではなく、この試験の診断基準として登録する (統計的な同等性・一般の再現性の証明とは呼ばない)。過去の判定は変えない。
- **問い**: `space.pRef: 0 → 5,499,000 Pa` の 1 点の変更だけで、§6.27 で登録した下流の θ_r の有限期間の差を参照の尺度 D の 0.1 倍以内に抑えられるか。第 1 仮説 (確度 低): 抑えられる。第 2 仮説 (未確認): 上流の射影の減少だけでは下流の目標を満たさず、他の演算・初期化の差との相互作用が残る。第 3 仮説 (未確認): 旧対照との実行条件の差が見かけの改善を作る (今回の陽性対照で点検する)。
- **腕** (AWS A、4 本同時、30,000 step 固定・500 ごとに出力 (step 1〜19,500 の場は判定に使わないので途中で消す)、共通の入力 `_tr/q32_init.h5`、§6.26 と同じ設定・環境変数 (`FORGE_CUDA_BLOCKSIZE=128`・`FORGE_OMEGA_BUDGET=1`・`FORGE_DIAG_COMMIT_LOSS=500`)、バイナリは旧 A・B・B′ の実ファイル `~/forge-fgeom7-f32` (sha256 `79321591…`)・`~/forge-fgeom7-fp64` (`7e95e5bb…`)):
  - F0 = float・pRef なし ×1 (`case/45.isobutane_m6_d155/run_0709_lp_f0`、陽性対照、設定は `run_0488` と同じ)
  - Bp = float + `space.pRef: 5499000.0` ×2 (`run_0710_lp_bp1`・`run_0711_lp_bp2`)
  - Ap = FP64 + 同じ pRef ×1 (`run_0712_lp_ap`)
  - 起動の順と同時実行の条件を固定する: F0 → Bp1 → Bp2 → Ap の順に、各起動の前に A の available ≥ 4.5 GB を確かめて 3 分ずつずらして起動し、4 本を同時に走らせる (30 分待っても空かなければ中止し、判定不能として記録する)。台本 `lp.sh`、判定 `lp_an.py` → `_band_ab/cold_pair/lp_judge.json`。
- **照合** (欠ければ判定不能、結果は判定記録に残す): バイナリの sha256、入力の sha256、設定の差 (F0 は `run_0488` と差 0 行、Bp は `run_0488` と `space.pRef` の 1 行だけ、Ap は `run_0487` と `space.pRef` の 1 行だけ)、初期化後の保存量 9 つの**キーの完全一致と値のビット一致** (F0・Bp は旧 B の res_0、Ap は旧 A の res_0)、座標の一致、全期間の残差と使う出力 (step 0・20,000〜30,000) の有限性。
- **指標**: Δ = 100 × (各 run の窓の平均 − その run の step 0) / |旧 A の step 0|。W1 = 20,000〜25,000、W2 = 25,000〜30,000、D = (+1.892, +1.647, +1.431) % (θ_r の 3 断面)。独立な反復は run で、窓の中の 11 点や 6 組を反復に数えない。旧 A・B・B′ の系列は `tr_judge.json`。
- **前提**: 全 6 組で |ΔAp − ΔA| ≤ 0.1|D|。F0 が旧 B・B′ の範囲を N_B だけ両側へ広げた区間に入り、かつ sign(D)·(ΔF0 − ΔAp) ≥ 0.5|D| (6 組すべて)。欠ければ判別不能。
- **判定**:
  - **結果 A**: 全 6 組で、Bp の 2 本がそれぞれ |ΔBp − ΔAp| ≤ 0.1|D| → 「今回の 2 本・登録の窓では目標を満たした」。第 2 仮説の「この条件でも不十分」を退ける。**0.1|D| 以内は「観測した差を 90 % 減らした」と同じ意味ではない** (D は別の停止時点の尺度)。
  - **結果 B**: 全 6 組で、Bp の 2 本がそれぞれ sign(D)·(ΔBp − ΔAp) ≥ 0.5|D| → 第 1 仮説の「この変更だけで十分」を退ける。
  - それ以外は判別不能。各組の残りの差と、無介入 (F0) との差も数値で残す。
  - Q_w は独立した副判定 (θ_r が A でも Q_w の改善は主張しない)。
- **記録**: `check_convergence` と θ_r・Q_w の系列 (20,000〜30,000) の `check_quasisteady` (1e-4) の VERDICT。未収束なら有限期間の応答として報告する。
- **やらないこと**: 長い試験の陽性を「残差の加算が真因」「double を廃止できる」へ広げること。既定値の変更・本番への採用をこの試験だけで決めること。
- **所要**: 4 本同時で 1〜1.5 時間 (A の forge-25 の run と共有)。メモリ 計 約 9 GB、ディスク 1 本 約 2.5 GB (途中の場を消す)。

### 6.47 §6.46 の結果 (2026-10-11)

- run (AWS A、4 本同時、各 30,000 step、RUN_RC 0): F0 `case/45.isobutane_m6_d155/run_0709_lp_f0` (float・pRef なし)、Bp `run_0710_lp_bp1`・`run_0711_lp_bp2` (float + `space.pRef: 5499000.0`)、Ap `run_0712_lp_ap` (FP64 + pRef)。判定 `lp_an.py` → `_band_ab/cold_pair/lp_judge.json`。
- **照合はすべて通った**: バイナリ・入力の sha256 が §6.26 と同じ、設定の差は F0 0 行・Bp と Ap は `space.pRef` の 1 行だけ、res_0 の 9 保存量はキーの完全一致と旧 B・旧 A とのビット一致、座標の一致、有限性。`check_convergence` は 4 本とも CHECK FAILURES。`check_quasisteady` (θ_r・Q_w、20,000〜30,000、1e-4) は 4 本とも DRIFTING (Ap の θ_r(40) は drift 0.0 %/窓・fluct 0.1 %、float の 3 本は drift 0.4〜1.0 %/窓)。
- **判定 (事前登録のとおり): 判別不能**。理由は前提の F0 (陽性対照) の範囲: W2 の 3 断面で ΔF0 (+2.167・+1.715・+1.426 %) が旧 B・B′ の範囲を N_B だけ広げた上限 (+1.989・+1.597・+1.340 %) を超えた。前提の Ap (全 6 組で |ΔAp − ΔA| ≤ 0.1|D|) と、F0 の「Ap から 0.5|D| 以上」は満たした。
- **記録 (判定ではない)**: Bp の 2 本は全 6 組で結果 B の条件 (sign(D)(ΔBp − ΔAp) ≥ 0.5|D|) を満たした。

| 窓 | 断面 | ΔA | ΔB / ΔB′ | ΔF0 | ΔAp | ΔBp1 / ΔBp2 |
| --- | --- | --- | --- | --- | --- | --- |
| W1 | θ_r(40) | −0.056 % | +1.577 / +1.439 % | +1.620 % | −0.039 % | +1.417 / +1.700 % |
| W1 | θ_r(70) | +0.012 % | +1.229 / +1.120 % | +1.229 % | −0.001 % | +1.105 / +1.303 % |
| W1 | θ_r(94) | +0.000 % | +0.995 / +0.908 % | +0.981 % | −0.013 % | +0.892 / +1.045 % |
| W2 | θ_r(40) | −0.121 % | +1.883 / +1.776 % | +2.167 % | −0.048 % | +1.799 / +1.891 % |
| W2 | θ_r(70) | −0.056 % | +1.512 / +1.426 % | +1.715 % | −0.005 % | +1.440 / +1.551 % |
| W2 | θ_r(94) | −0.054 % | +1.267 / +1.194 % | +1.426 % | −0.024 % | +1.199 / +1.308 % |
| W1 / W2 | Q_w (副) | +0.085 / +0.169 % | +2.816・+2.524 / +2.749・+2.400 % | +2.872 / +2.827 % | +0.155 / +0.177 % | +2.460・+1.940 / +2.144・+1.825 % |

- **書いてよい範囲** (codex diagnose `notes/reviews/2026-10-11-after-pref-long-diagnose.md` で訂正): 正式には判別不能。**登録した抑制の目標の達成は確認できず、今回の 2 本には大きな差が残った**。W2 の θ_r の差は F0 − Ap が +1.450〜2.215 %、Bp − Ap が +1.223〜1.939 % で、対応する断面で約 8〜17 % 小さい (F0 は 1 本・前提不成立なので比較は記録)。~~pRef なしの float とほぼ同じ、型の減少は下流の差の減少につながらなかった~~。Q_w の Ap からの差は全窓で +1.648〜2.304 ポイント (上の表の値は各 run の起点からの変化で、Ap の変化を引いていない)。副判定で主張はしない。~~下流の差の通り道は短い試験の型とは別~~ **訂正**: 初期の射影の縮小と後期の差の残存だけでは、別の成分の成長と、残った同じ成分の遅れた増幅を分けられない。`lp_diffmap.py` の領域別の絶対差 (`_band_ab/cold_pair/lp_diffmap.txt`、step 20,000〜30,000 で収縮部〜スロート付近の |ΔT|/T の中央値が 1e-2 程度で最大、pRef ありで F0 より約 2 割小さい) は、場所・時点・正規化を付けた観測であって、モードや因果の経路を同定していない。
- 関連: forge-25 の観察 (§5.1 #22)。採用するのは「B0 の参照が未収束の途中の状態かもしれない」まで (Gref の観察で G0 とは壁法線の解像度が違い、point 化は `lineImplicit`・`lineDtDirectional`・ITJ の組を同時に変える)。H-slow (遅い方向への溜まり) と H-res (残差の偏り) は排他的でない。「FP64 の揺れは 1e9 倍小さい」は演算精度の比で、軌跡の揺れの実測ではない。
- 次は諮ってから (候補: Bp と Ap の 20,000〜30,000 の場で残る差の場所を調べる (run 不要)、point で float と FP64 を比べる)。

### 6.48 事前登録: point の陰解法に固定した float と FP64 の比較 (案 P、2026-10-11、codex diagnose、実行前)

- **ユーザの選択 (2026-10-11)**: 案 S (収縮部〜スロートの小さい格子) より案 P を先に (案 S は x/r_t ≤ +2 で切っても節点が 57 % 残り、速さは 1.8 倍程度と分かったため)。案 S・E (面エンタルピーの double)・低マッハ補正 (§5.1 #21) は同時に始めない (codex)。
- **問い** (codex diagnose `notes/reviews/2026-10-11-after-pref-long-diagnose.md`): 陰解法を point に固定しても、float と FP64 の大きな差 (§6.27) が残るか。「B0 の構成がこの大きな精度差に要るか」を調べる。第 1 仮説 (確度 低): point でも十分に進めれば大きな差が残る (B0 の停滞だけでは説明できず、初期化や step の中の演算の精度差が関わる)。第 2 仮説 (確度 低): B0 の反復の特性が精度差を大きく増幅しており、point では大差を保たない (forge-25 の観察が動機、G0 で未確認)。第 3 仮説 (未確認): 精度ごとの初期化の差の増幅が一部を担う。
- **腕** (AWS B、4 本同時): 変える因子はビルドの精度だけ。両精度とも §6.26 の設定から `time.deltaT` の `lineImplicit`・`lineDtDirectional`・`lineDtDirectionalCap`・`implicitThermalJacobian` の 4 キーを除く (= point、既定 0。forge-25 の point と同じ書き方)。ほか (blockDPLUR 1・cfl 4・cfl_pseudo 4・implicitRelax 0.7・nStepInner 5・空間の離散化・BC・pRef 0) と環境変数は §6.26 のまま。共通の入力 `_tr/q32_init.h5`、G0 の全格子、旧バイナリ `~/forge-fgeom7-{fp64,f32}`。FP64 ×2 (`case/45.isobutane_m6_d155/run_0713_pt_d1`・`run_0714_pt_d2`、FP64 の再実行の幅も測る)、float ×2 (`run_0715_pt_f1`・`run_0716_pt_f2`)。500 step ごとに出力し、判定に使わない場 (step 1〜19,500、30,500〜99,500) は途中で消す。台本 `rpt.sh`、判定 `rpt_an.py` → `_band_ab/cold_pair/rpt_judge_s1.json`・`rpt_judge_s2.json`。
- **期間** (最大 120,000 step を事前登録): まず 30,000 step で段 1 の判定 (窓 W1 = 20,000〜25,000・W2 = 25,000〜30,000)。段 1 で結果 A なら run を止めて終える。結果 A でなければ 120,000 step まで続けて段 2 の判定 (窓 W3 = 100,000〜110,000・W4 = 110,000〜120,000)。12 万 step そのものを収束の保証にはしない。差が小さいだけで抑制と読まない。
- **指標**: θ_r(40/70/94) の Δ = 100 × (窓の平均 − 自分の step 0) / |旧 A の step 0| (§6.26 と同じ)。D = (+1.892, +1.647, +1.431) % は尺度としてだけ使う。N_d = |ΔD1 − ΔD2|、N_f = |ΔF1 − ΔF2|、N = max(N_d, N_f)。独立な反復は run で、窓の中の点を反復に数えない。
- **判定**:
  - **結果 A (大差が残る)**: その段の両窓・3 断面で、全組 (F_i − D_j、i, j ∈ {1, 2}) の min sign(D)·(ΔF_i − ΔD_j) − N ≥ 0.5|D| → 第 1 仮説を支持し、「B0 の構成が大差に必須」を退ける。面エンタルピー (や他の特定の演算) が原因とは言わない。
  - **結果 B (大差が残らない)**: 段 2 の両窓 (W3・W4)・3 断面で、全組の max |ΔF_i − ΔD_j| + N ≤ 0.1|D|、かつ 4 本とも `check_convergence` が PASS、かつ θ_r の系列 (100,000〜120,000) の `check_quasisteady` (1e-4) が STEADY → 第 1 仮説の予測を退け、第 2 仮説を支持する。B0 のどの設定が効いたかは未確定。
  - 進行の不足・未収束・ドリフト・反復の間のばらつき・中間の差は判別不能。
- **前提と記録**: 各精度の中で res_0 の 9 保存量がビット一致 (キーの完全一致)、精度の間の初期化の差は記録。設定の差は `run_0487`・`run_0488` からの 4 キーの削除だけ (**追記、run の開始直後・判定の前**: `nStepOuter` の 30,000 → 120,000 も必然の差として許す。登録の時点で書き漏らした)。バイナリ・入力の sha256。全期間の残差・使う出力の有限性。Q_w、上流・スロートの温度差 (`lp_diffmap.py` 相当)、全残差、局所擬似時間の刻みを記録する。
- **所要**: 段 1 は 4 本同時で 約 1 時間、段 2 まで進めば 約 4〜5 時間。ディスクは B の空き 147 GB に対し最大 約 35 GB。
- **やらないこと**: この試験を理由に opt-in を削除したり既定値を変えたりしない。

### 6.49 §6.48 の段 1 の記録 (30,000 step、2026-10-11、判定は段 2 へ持ち越し)

- run (AWS B、4 本同時、走行中): FP64 `case/45.isobutane_m6_d155/run_0713_pt_d1`・`run_0714_pt_d2`、float `run_0715_pt_f1`・`run_0716_pt_f2` (point に固定、最大 120,000 step)。段 1 の判定 `rpt_an.py --stage 1` → `_band_ab/cold_pair/rpt_judge_s1.json`。前提 (設定の差は 4 キーの削除と nStepOuter、res_0 の精度の中のビット一致、有限性) は満たした。初期化の差 (D1 と F1 の res_0、最大相対差): ρE 5.5e-4・ρω 3.0e-7・ρY 5〜6e-8、ほかは 0。`check_convergence` は 4 本とも CHECK FAILURES (段 1 の時点、判定には使わない)。
- **段 1 の判定 (事前登録のとおり): 判別不能 → 段 2 へ持ち越し** (結果 A の条件 min sign(D)(ΔF − ΔD) − N ≥ 0.5|D| は成立しない)。
- **記録 (判定ではない)**: point では float と FP64 の θ_r の Δ が両窓・3 断面で表示の桁 (0.001 ポイント) まで一致した。

| 窓 | 断面 | ΔD1 / ΔD2 | ΔF1 / ΔF2 |
| --- | --- | --- | --- |
| W1 (20,000〜25,000) | θ_r(40) / (70) / (94) | −1.135 / −0.819 / −0.680 % (2 本とも) | −1.135 / −0.819 / −0.680 % (2 本とも) |
| W2 (25,000〜30,000) | θ_r(40) / (70) / (94) | −1.485 / −1.124 / −0.952 % | −1.485 / −1.124 / −0.952 % |
| W1 / W2 | Q_w (副) | +0.334 / +0.339 % | +0.326〜0.327 / +0.328 % |

- 同じ起点から B0 では float だけが +0.9〜1.9 % 動いた (§6.27) のに対し、point では float と FP64 が同じ向き・同じ量だけ動く (θ_r −0.7〜−1.5 %、まだ動いている)。動きの向きは forge-25 の観察 (§5.1 #22: Gref で B0 の止まる値は point の収束値より θ_r +0.6 % 程度大きい) と同じ。
- **書いてよい範囲 (段 2 の前、上位にまだ諮っていない)**: この起点・この 30,000 step では、point に固定すると float と FP64 の θ_r の差は検出できない大きさ。§6.27 の大差は B0 の構成のもとで生じた。B0 のどの設定 (ライン・方向別 dt・ITJ・上限) が要るかは未確定。段 2 (120,000 step、収束・準定常の付帯条件) の正式の判定を待つ。

### 6.50 §6.48 の段 2 の結果 (120,000 step、2026-10-11)

- run (AWS B、4 本同時、各 120,000 step、RUN_RC 0): FP64 `case/45.isobutane_m6_d155/run_0713_pt_d1`・`run_0714_pt_d2`、float `run_0715_pt_f1`・`run_0716_pt_f2`。判定 `rpt_an.py --stage 2` → `_band_ab/cold_pair/rpt_judge_s2.json`。前提は満たした (段 1 と同じ)。残差は起点から 1.9 桁下がった (rms_ro 3.7〜4.3e-7、worst は rms_roOmega)。
- **判定 (事前登録のとおり): 判別不能** — θ_r の差の条件 (両窓・3 断面で max|ΔF − ΔD| + N ≤ 0.1|D|) は満たした (max|F − D| + N は 0.026〜0.034 ポイント、0.1|D| は 0.143〜0.189) が、付帯条件の `check_convergence` (4 本とも CHECK FAILURES) と θ_r の `check_quasisteady` (100,000〜120,000、1e-4、4 本とも DRIFTING、drift 0.1 %/窓) を満たさない。

| 窓 | 断面 | ΔD1 / ΔD2 | ΔF1 / ΔF2 | max\|F − D\| + N |
| --- | --- | --- | --- | --- |
| W3 (100,000〜110,000) | θ_r(40) / (70) / (94) | −1.466 / −1.234 / −1.140 % | −1.492〜−1.493 / −1.258〜−1.260 / −1.163〜−1.164 % | 0.029 / 0.027 / 0.026 |
| W4 (110,000〜120,000) | θ_r(40) / (70) / (94) | −1.396 / −1.170 / −1.080 % | −1.427〜−1.428 / −1.200〜−1.202 / −1.110〜−1.111 % | 0.034 / 0.033 / 0.032 |
| W3 / W4 | Q_w (副) | −0.448 / −0.510〜−0.511 % | −0.464 / −0.524 % | 0.017 / 0.014 |

- **書いてよい範囲 (codex diagnose `notes/reviews/2026-10-11-point-result-diagnose.md` で訂正)**: point に固定すると、**登録した窓 (20,000〜30,000・100,000〜120,000) で** float と FP64 の θ_r の差は小さかった (W4 の max|F − D| + N は 0.032〜0.034 ポイント、符号つきの平均の差は −0.030〜−0.031 ポイントで B0 の正方向の差と符号も違う)。30,500〜99,500 の場は消したので全期間の上限ではない。「point では同じ解へ向かう」とは言わない (4 本とも CHECK FAILURES・θ_r は DRIFTING)。ただし両精度とも θ_r はまだ動いており (W3 → W4 で +0.07 ポイント、向きが反転して戻りつつある)、収束・準定常の条件を満たさないので、「point では精度差が無い」とは確定しない。差は時間とともにわずかに増えている (W1・W2 の 0.000 → W3 0.026〜0.029 → W4 0.032〜0.034)。
- §5.1 #22 (forge-25 の観察: B0 の停滞は縮流部の主流側の遅いモード、上限・方向別 dt は右辺を変えない。**場所は後に「境界層の縁」と訂正された**、§5.1 #22) と合わせた読みは仮説のまま。次は上位に諮る (result の解釈と、上限 1 での float と FP64 の比較などの次の手)。

### 6.51 事前登録: 方向別 dt の上限だけを 50 → 1 にした float と FP64 の比較 (案 C1、2026-10-11、codex diagnose、実行前)

- **codex diagnose** (`notes/reviews/2026-10-11-point-result-diagnose.md`) の採否: 全件採用。§6.50 の言い方を訂正 (上)。「B0 の遅いモードを丸めの揺れが押す」は仮説のまま (point 化はライン・方向別 dt・ITJ を同時に変えた、Gref の局所の自己相関と場所だけでは G0 の精度差が同じモードとも駆動が揺れとも確定しない)。C1 は旧 B0 の結果と比べず、上限 50 の陽性対照を同じ実行条件で取り直す。共通の Q32 の入力と初期化後の状態の一致は別 (精度の間の res_0 の ρE の最大相対差 5.5e-4)。「上限・方向別 dt は右辺を変えない」は、観測した停止状態が右辺ゼロの同じ解だという証明ではない。効かなかった介入を一括して候補の除外の根拠にしない。
- **問い**: B0 で観測した 30,000 step の大きな精度差に方向別 dt の拡大が関わり、上限 1 でその期間の差を小さくできるか。第 1 仮説 (確度 中): できる。第 2 仮説 (未確認): ライン結合・ITJ・精度依存の演算・初期化と反復の相互作用が残り、上限だけでは不十分。第 3 仮説 (未確認): dt の縮小は差の発達を遅らせるだけで長期の差を除かない。
- **腕** (AWS B、**設定の差は `lineDtDirectionalCap` 50 対 1 だけ**、それぞれ FP64 ×2・float ×2、計 8 本): 共通の入力 `_tr/q32_init.h5`、§6.26 の B0 の設定 (`lineImplicit 1`・`lineDtDirectional 1`・ITJ 5・cfl 4・緩和 0.7) とバイナリ (`~/forge-fgeom7-{fp64,f32}`)・環境変数、出力は §6.26 と同じ量に `dt_local` を足す (両方の上限で同じ)。各 30,000 step・500 ごと (step 1〜19,500 の場は消す)。同時に 4 本ずつ 2 回、各回に上限 50・1 と両精度を 1 本ずつ置き、起動の順も両方の回で同じ (上限 50 FP64 → 上限 1 FP64 → 上限 50 float → 上限 1 float、2 分ずつずらす)。
  - 1 回目: `case/45.isobutane_m6_d155/run_0717_c50_d1`・`run_0718_c1_d1`・`run_0719_c50_f1`・`run_0720_c1_f1`
  - 2 回目: `run_0721_c50_d2`・`run_0722_c1_d2`・`run_0723_c50_f2`・`run_0724_c1_f2`
  - 台本 `rcap.sh`、判定 `rcap_an.py` → `_band_ab/cold_pair/rcap_judge.json`。
  - **実行の記録 (2026-10-11)**: `rcap.sh` は 1 回目 (10:54〜12:13 UTC) の後、設定の差を記録する `diff` の終了コード 1 を失敗と受け取って止まった (中止の記録も残らない形)。不具合を直し (`; true`)、2 回目と判定だけを回す `rcap2.sh` で続ける。2 回目は同じ B で、同じ起動順・同じ 4 本同時 (投入の前に G0c の変換の終了を待つ)。1 回目との間が約 1.5 時間空いた。
  - **判定器の修正 (2026-10-11、結果の前・2 回目の走行中に B に置いた。codex diagnose `2026-10-11-g0c-precision-test-diagnose.md` の不足情報)**: `rcap_an.py` の設定の照合は各精度の 1 回目どうしだけだった。8 本すべてを `run_0717_c50_d1` と比べ、上限 50 は差なし・上限 1 は上限の行だけ、に広げた。判定の式は変えていない。
- **指標**: W1 = 20,000〜25,000、W2 = 25,000〜30,000。Δ・D・N は §6.48 と同じ定義で、N は上限ごとに求める。
- **前提** (欠ければ判別不能):
  1. 上限 50 (陽性対照) で両窓・3 断面とも min sign(D)(ΔF − ΔD) − N ≥ 0.5|D|。
  2. 実効の dt が縮んだ: step 20,000 の `dt_local` の、縮流部 (x/r_t ∈ [−8, −2)) の壁側 (j ≥ 100) と全域で、上限 1 の FP64 の上限 50 の FP64 に対する比の中央値 < 0.5 (縮んだ程度は記録)。
  3. 入力・実バイナリの sha256、各精度の中の res_0 の 9 保存量のビット一致、全期間の残差と使う出力の有限性。
- **判定**:
  - **結果 A**: 上限 1 で両窓・3 断面とも max|ΔF − ΔD| + N ≤ 0.1|D| → 第 1 仮説を支持し、「この期間でも上限の変更だけでは不十分」を退ける。長期の解の一致や遅いモードの同定とはしない。
  - **結果 B**: 上限 1 でも両窓・3 断面とも min|ΔF − ΔD| − N ≥ 0.5|D| → 第 1 仮説を退ける。符号も記録し、逆符号の大差を改善と扱わない。
  - 中間・混在・陽性対照の不成立は判別不能。`check_convergence` と θ_r・Q_w の `check_quasisteady` を添え、未収束なら有限期間の応答として扱う。
- **やらないこと**: point の延長や C2 を同時に始めること、30,000 step の小差を恒久的な抑制と呼ぶこと、この結果だけで double の廃止・既定値の変更・opt-in の削除を決めること。
- **所要**: 2 回 × 約 1 時間 (B)。ディスク 計 約 25 GB。

### 6.51 の結果 (2026-10-11、判定 `_band_ab/cold_pair/rcap_judge.json`)

- **前提**: すべて満たした。バイナリ・入力の sha は §6.26 と同じ、8 本の設定の差は上限の行だけ (2 回目の組も照合)、各精度の中の res_0 の 9 保存量はビット一致、全期間の残差と使う出力は有限。step 20,000 の `dt_local` の比 (上限 1 / 上限 50、FP64 の 1 回目どうし) の中央値は縮流部の壁側 0.020 (10 % 点〜90 % 点 0.0200〜0.0200)、全域 0.161。
- **判定**: **結果 A** (上限 1 で両窓・3 断面とも max|ΔF − ΔD| + N ≤ 0.1|D|)。上限 50 の陽性対照は成立 (min sign(D)(ΔF − ΔD) − N = +0.77〜+1.52 ポイント ≥ 0.5|D| = 0.72〜0.95)。

| 上限 | 窓 | 断面 | ΔD1 / ΔD2 | ΔF1 / ΔF2 | N | max\|F − D\| + N |
| --- | --- | --- | --- | --- | --- | --- |
| 50 | W1 | θ_r(40) / (70) / (94) | −0.063 / +0.057、−0.010 / +0.073、−0.016 / +0.045 % | +1.387 / +1.486、+1.080 / +1.169、+0.878 / +0.945 % | 0.120 / 0.089 / 0.067 | 1.669 / 1.268 / 1.027 |
| 50 | W2 | 同上 | −0.070 / +0.016、−0.021 / +0.046、−0.029 / +0.033 % | +1.622 / +1.657、+1.307 / +1.355、+1.095 / +1.151 % | 0.086 / 0.067 / 0.061 | 1.812 / 1.444 / 1.241 |
| 1 | W1 | 同上 | −1.082、−0.779、−0.649 % (2 本とも) | −1.083、−0.779、−0.649 % (2 本とも) | 0.000 | 0.000 / 0.000 / 0.000 |
| 1 | W2 | 同上 | −1.400、−1.055、−0.897 % (2 本とも) | −1.400、−1.056、−0.897 % (2 本とも) | 0.000 | 0.000 / 0.001 / 0.000 |
| 50 / 1 | W1・W2 | Q_w (副) | 50: +0.139〜+0.258 %、1: +0.290〜+0.296 % | 50: +2.113〜+2.582 %、1: +0.273〜+0.274 % | — | 50: 1.990 / 2.617、1: 0.018 / 0.023 |

  Δ は旧 A (`run_0487`) の step 0 で割った各 run の step 0 からの変化 [%]。
- `check_convergence`: 8 本とも OVERALL: CHECK FAILURES ABOVE。`check_quasisteady` (θ_r・Q_w、20,000〜30,000): 8 本とも DRIFTING。step 30,000 の rms_ro は上限 1 で 2.2e-7・3.2e-7 (2 回目、2.0〜2.2 桁の低下)、上限 50 で 2.7e-5・4.1e-5 (2 回目、横ばい)。
- **書いてよい範囲**: G0 の B0 で、方向別 dt の上限だけを 50 → 1 にすると、登録した窓 (20,000〜30,000 step) で float と FP64 の θ_r の差は 0.001 ポイント以下になった (上限 50 では同じ条件で +0.9〜+1.7 ポイント)。第 1 仮説 (30,000 step の大きな精度差に方向別 dt の拡大が関わり、上限 1 でその期間の差を小さくできる) を支持し、「この期間でも上限の変更だけでは不十分」を退ける。長期の解の一致・遅いモードの同定・float の丸めの通り道の特定とはしない。両精度とも未収束 (CHECK FAILURES・DRIFTING) なので有限期間の応答として扱う。
- **記録 (判定外)**: 上限 1 では両精度とも θ_r が B0 の起点から −0.65〜−1.40 % 動き、同じ起点からの point (§6.49 の W1・W2: −1.135 / −0.819 / −0.680 %、−1.485 / −1.124 / −0.952 %) の 94〜95 % の位置にある。上限 50 の FP64 は起点に留まる (|Δ| ≤ 0.073 %)。forge-25 の Gref の観察 (上限 1 で θ_r が point の止まった値へ動き出す、§5.1 #22 続報 2) と向きが同じ。
- 判定台本の版: 設定の照合を 8 本に広げた版 (sha 807a68d0) が回った。前提 2 が欠けたときにも Δ を記録する修正 (結果の前に作った) は判定の開始に間に合わなかったが、前提はすべて満たしたので結果は変わらない。
- 解釈 (上の範囲を超える読み) は、§6.52 の結果とまとめて上位に諮る。

### 6.52 事前登録: 縮流部を粗くした格子 G0c で、B0 のまま float と FP64 を比べる (2026-10-11、codex diagnose、実行前)

- **codex diagnose** (`notes/reviews/2026-10-11-g0c-precision-test-diagnose.md`、ブリーフ `notes/reviews/briefs/2026-10-11-g0c-precision-test.md`) の採否: 全件採用。
  - Major 1: 旧 G0 の結果だけでは陽性対照が足りない → §6.51 の上限 50 の 4 本 (`run_0717_c50_d1`・`run_0721_c50_d2`・`run_0719_c50_f1`・`run_0723_c50_f2`) で大差が再現したことを前提 1 にする。G0c の比較相手は G0c の FP64。
  - Major 2: `g0c_prep.sh` は失敗しても続いて `.done` を作った → 各段で止める形に直し、作り直した (下)。
  - Major 3: `make_g0c.py` が比べたのは作り直した G0 で、実際の入力 `_tr/q32_init.h5` ではない → `g0c_check.py` で実際の入力と照合した。**照合で欠陥が見つかった** (下)。
  - Major 4: 起点からの増分の差と、窓の値そのものの差は別 → 結論は「θ_r の増分の差」に限り、値の差は別の条件 (下) を満たすときだけ書く。精度の間の res_0 の差は記録する。
  - Minor 1: D の出典 (§6.25 から持ち込んだ尺度) と正規化を固定 → D_c = D_old × S_old / S_c。N は 2 反復で観測した幅で、信頼区間とは呼ばない。
  - Minor 2: 大差は絶対値で判定し、符号は別に記録する。
  - あわせて: 同じ B で §6.51 の 4 本に G0c の 4 本を足さない (同時実行の条件を変えない)。投入は §6.51 の終了と陽性対照の成立の後。
- **G0c の準備の記録 (2026-10-11)**:
  - **欠陥と修正**: `make_g0c.py` は `mesh2d.write_msh41_2d` に座標の桁数を渡さず、既定の有効数字 10 桁で .msh を書いていた。生産の書き出し (`runner_axismach.py:1377`) は問題 YAML の `mesh.msh_digits` (冷却壁は 17) を使う。実際の G0 の入力との照合で、区間の外の全節点の座標が相対 ±5e-10 (10 桁の丸めの幅) ずれ、壁際 (j 120) では局所の最小の格子間隔の 0.13 % に当たった。`make_g0c.py` を `msh_digits` (17、assert) で書くように直し、作り直した。10 桁の版 (`_g0c_d10/`、入力の sha 5bd456f340f77ede) は使わない (B に残し、§6.52 の判定後に消す)。
  - 作り直した G0c: `.msh` sha 285f0fa4ecc0ed1c (`make_g0c.py 96197821`、forge-25 の `core_grid_mesh.py`)。区間 x/r_t −7.99〜−2.00 の列 589 → 337 (ni 4467)、つなぎ目の dx の比 0.996・0.999。
  - `g0c_prep.sh` (失敗で止まる版、記録 `g0c_prep.log`): 変換 (`~/forge-fgeom7-fp64` の変換器) rc 0。閉性 (`hoop_verify_conv.py check --rw`) の r 重みの max E_x 9.2e-10・E_y 2.3e-9、基準超え 0。品質 AR 最大 8,704.6・スキュー最大 0.457 で `--ar-max 10000` は PASS (既定 1000 では FAIL、AR > 1000 のセル 9.47 %。ユーザ決定により診断に限り許す)。`interp_field.py` は 9 保存量を移した (最近傍距離 0 の節点 499,851 / 540,507)。Q32 の丸めで変わった値 0。
  - **実際の G0 の入力との照合** (`g0c_check.py` → `_band_ab/cold_pair/g0c_check.json`): 区間の外 (前 100 列・後 4,030 列) は座標・9 保存量ともビット一致 (VERDICT: 照合 OK)。区間の中の最近傍の距離は最大 0.0170・99 % 点 0.0100・中央値 0.0029 r_t。壁法線の行 j が同じ行に写った節点は 99.89 %。違う行に写った 44 節点 (区間の 0.11 %) はすべて壁側の j 106〜120 で、壁から 1〜2 行離れた節点の値を受け取った。両方の精度が同じ起点から始まるので比較は成り立つが、起点の作り方が G0 と違う点の一つとして記録する。
  - G0c の共通の入力: `_g0c/q32_init.h5` sha 6f8fcda6d97482e3。
- **問い**: 縮流部を流れ方向に粗くした格子では、B0 のままでも 30,000 step の float と FP64 の θ_r の差が小さくなるか。第 1 仮説 (確度 中): 縮流部の格子と B0 の反復の性質の組合せが精度の差の増幅に関わり、G0c では登録した窓で増分の差が小さくなる。第 2 仮説 (未確認): 格子は差の育つ速さを変えるだけで、30,000 step より後に大差が育つ。第 3 仮説 (未確認): 格子の移植・初期化と精度の相互作用が応答を変える。
- **腕** (AWS B、4 本同時、**設定は §6.51 の上限 50 の `run_0717_c50_d1` の `solverConfig.yaml` とバイト一致** (B0・出力に `dt_local`)、バイナリ・環境変数も §6.51 と同じ、入力だけ G0c): 各 30,000 step・500 ごと (step 1〜19,500 の場は消す)。起動は FP64 → float → FP64 → float の順に 2 分ずつずらす。
  - `case/45.isobutane_m6_d155/run_0725_g0c_d1`・`run_0726_g0c_d2` (FP64)、`run_0727_g0c_f1`・`run_0728_g0c_f2` (float)
  - 台本 `rg0c.sh` (§6.51 の `rcap.done` を待ち、`rg0c_an.py --pc-only` で陽性対照が成立し、available ≥ 11 GB のときだけ投入)、判定 `rg0c_an.py` → `_band_ab/cold_pair/rg0c_judge.json`。
- **指標** (θ_r(40/70/94) が主、Q_w は副で記録だけ): W1 = 20,000〜25,000、W2 = 25,000〜30,000。
  - S_c = `run_0725_g0c_d1` の step 0 の |θ_r| (断面ごと、全腕に共通)。Δ = 100 × (窓の平均 − 各腕の step 0) / S_c。
  - D_c = D_old × S_old / S_c。D_old = (+1.892, +1.647, +1.431) % (§6.25 から持ち込んだ尺度)、S_old = 旧 A (`run_0487`) の step 0 の |θ_r|。
  - e_ij = ΔF_i − ΔD_j、N = max(|ΔD1 − ΔD2|, |ΔF1 − ΔF2|)。
  - 値の差: V_ij = 100 × (F_i の窓の平均 − D_j の窓の平均) / S_c、N_v = 精度の中の 2 本の窓の平均の差の大きい方 (同じ正規化)。
  - 断面の抽出は `cold_series.snapshot` と同じ量を ni = 節点数 / 121 で取る (`XC.mesh_info` は G0 の ni に固定のため使わない)。G0 の run では `cold_series.snapshot` とビット一致を確認済み (`run_0717_c50_d1` の step 0・20,000)。
- **前提** (欠ければ判別不能):
  1. 陽性対照: §6.51 の上限 50 の 4 本で、両窓・3 断面とも min sign(D)(ΔF − ΔD) − N ≥ 0.5|D| (§6.51 の前提 1 と同じ式、`rg0c_an.py` が 4 本から直接計算する)。
  2. G0c の準備: `g0c_prep.log` が中止なく終わり、入力の sha が 6f8fcda6d97482e3、`g0c_check.json` が照合 OK。
  3. 抽出の検査: `run_0725_g0c_d1` と `run_0717_c50_d1` (G0、同じ FP64 のバイナリ・同じ設定) の step 0 の θ_r(40/70/94) の相対差 ≤ 1e-9 (断面は区間の外で、入力はビット一致)。
  4. 入力・実バイナリの sha256、4 本の設定が `run_0717_c50_d1` とバイト一致、各精度の中の res_0 の 9 保存量のビット一致、全期間の残差と使う出力の有限性。精度の間の res_0 の差は記録だけ。
- **判定** (θ_r の 3 断面・両窓):
  - **結果 A (大差が残る)**: min|e| − N ≥ 0.5|D_c| → 第 1 仮説の有限期間の予測を退ける。符号も記録する。
  - **結果 B (差が小さい)**: max|e| + N ≤ 0.1|D_c| → 第 1 仮説の有限期間の予測を支持する。書いてよいのは「G0c では登録した窓で θ_r の増分の差が小さい」まで。**値の差**も max|V| + N_v ≤ 0.1|D_c| を満たすときだけ「θ_r の値の差が小さい」と書く。
  - 中間・断面の間の混在・前提の不成立は判別不能。`check_convergence` と θ_r・Q_w の `check_quasisteady` を添え、未収束なら有限期間の応答として扱う。
- **解釈の範囲**: 結果 B でも、格子が原因・遅いモードを同定した・double は不要、とは言わない (格子を変えると離散の作用素・局所の dt・初期場の移植も同時に変わり、寄与を分けられない)。
- **やらないこと**: §6.51 の走行中に G0c の run を足すこと。G0c を生産に使うこと (AR の緩和は診断に限る)。結果を見てから閾値・窓・正規化を変えること。
- **所要**: 約 1.2 時間 (B、4 本同時)。メモリ 計 約 9 GB、ディスク 約 12 GB (途中の場を消す)。

### 6.52 の結果 (2026-10-11、判定 `_band_ab/cold_pair/rg0c_judge.json`)

- **前提**: すべて満たした。バイナリ・G0c の入力の sha は登録と同じ、4 本の設定は `run_0717_c50_d1` とバイト一致、各精度の中の res_0 の 9 保存量はビット一致 (精度の間の res_0 の最大相対差は ρE 5.5e-4・ρω 3.0e-7、ほかは 0 か 1e-7 以下、記録だけ)、全期間の残差と使う出力は有限。陽性対照 (§6.51 の上限 50 の 4 本から直接計算) は min sign(D)(ΔF − ΔD) − N = +0.77〜+1.52 で成立。抽出の検査: step 0 の θ_r (40/70/94) は G0 と相対差 0.0。S_c = S_old なので D_c = D_old。
- **判定**: **判別不能** (結果 A の min|e| − N ≥ 0.5|D_c| に 3 断面・両窓とも届かず、結果 B の max|e| + N ≤ 0.1|D_c| も満たさない)。

| 窓 | 断面 | ΔD1 / ΔD2 (FP64) | ΔF1 / ΔF2 (float) | N | min\|e\| − N | max\|e\| + N | 符号つきの平均の差 F − D | 0.5\|D_c\| |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| W1 | θ_r(40) | −1.362 / −1.360 % | −0.185 / −0.489 % | 0.304 | +0.567 | 1.481 | +1.024 | 0.946 |
| W1 | θ_r(70) | −1.060 / −1.059 % | −0.176 / −0.392 % | 0.215 | +0.452 | 1.099 | +0.776 | 0.824 |
| W1 | θ_r(94) | −0.921 / −0.920 % | −0.199 / −0.369 % | 0.170 | +0.381 | 0.891 | +0.636 | 0.716 |
| W2 | θ_r(40) | −1.416 / −1.415 % | −0.106 / −0.410 % | 0.304 | +0.702 | 1.613 | +1.158 | 0.946 |
| W2 | θ_r(70) | −1.143 / −1.143 % | −0.090 / −0.333 % | 0.243 | +0.566 | 1.297 | +0.931 | 0.824 |
| W2 | θ_r(94) | −1.019 / −1.019 % | −0.117 / −0.330 % | 0.213 | +0.475 | 1.116 | +0.796 | 0.716 |
| W1 / W2 | Q_w (副) | −1.557〜−1.602 % | −1.217〜−1.408 % | 0.126 / 0.012 | +0.088 / +0.181 | 0.467 / 0.218 | +0.277 / +0.200 | 1.045 |

  Δ は G0c の FP64 の step 0 (= G0 と同じ値) で割った各 run の step 0 からの変化。値の差 max|V| + N_v も 0.89〜1.61 で小差の条件を満たさない。
- `check_convergence`: 4 本とも OVERALL: CHECK FAILURES ABOVE。`check_quasisteady` (θ_r・Q_w、20,000〜30,000): 4 本とも DRIFTING。
- **残差の推移 (記録、判定外)**: rms_ro は step 2,000 まで両精度がそろって下がる (FP64 7.5e-6、float 8.2〜8.6e-6)。その後 FP64 は 1e-6 前後 (step 30,000 で 7.5e-7・1.2e-6、step 0 から 3.6〜3.8 桁) まで下がるが、float は step 5,000 までに 2.0〜2.3e-5 へ上がり、30,000 まで 1.8〜2.1e-5 で止まる (G0 の上限 50 の停滞の水準 2.6〜4.2e-5 と同じ桁)。G0 の上限 1 は両精度とも 2〜4e-7 までそろって下がる。step 1〜19,500 の場は消したので、float が離れた場所は分からない。
- **書いてよい範囲 (事前登録の範囲)**: G0c の B0 では、登録した窓で float と FP64 の θ_r の増分の差は小さくならなかった (符号つきの平均で +0.64〜+1.16 ポイント、float の 2 本の間の幅 0.17〜0.30)。ただし事前登録の結果 A の閾値には届かないので、形式上は判別不能。FP64 は G0c で B0 の停滞を抜けて、θ_r が同じ起点からの point (§6.49) と同じ向きに動いた (同じ窓の point に対して W1 で 120〜135 %、W2 で 95〜107 %)。float の動きは FP64 の 8〜36 % にとどまった。解釈は codex に諮る (`notes/reviews/briefs/2026-10-12-cap-g0c-result.md`)。
- **codex diagnose** (`notes/reviews/2026-10-12-cap-g0c-result-diagnose.md`) の採否 (§6.51 と §6.52 をまとめて諮った): 全件採用。
  - 正式判定は変えない (§6.51 の θ_r の max|F − D| + N は 0.000547 ポイント以下、§6.52 は 0.891〜1.613 ポイントで小差でも大差でもない、を codex が JSON から再計算して再現)。
  - 書いてよいのは「上限 1 は登録窓の θ_r の精度差を小さくした。G0c の粗格子化だけでは登録した小差条件を満たさなかった」まで。上限 1 の Q_w には max|F − D| + N が 0.023 ポイント残るので、θ_r の数値を全量に広げない。12 本とも CHECK FAILURES・DRIFTING なので、収束解の一致・FP64 の収縮性・float の定常的な誤差の水準は主張しない。
  - (Major) 「精度の差は B0 の中でだけ育つ」「上限 1 では反復が収縮的」は測っていない。上限 1 は局所 dt を大きく縮めており、差の発達を遅らせただけの可能性も残る (第 3 仮説)。採用する言い方は「この起点・登録窓では、方向別 dt の拡大が精度差の発達に関与した」まで。
  - (Major) 「float の停滞の水準は格子より B0 と float の組で決まる」は、G0 と G0c の rms_ro が同じ桁というだけで同じ局所モードとは言えない → 初期場への依存を試す仮説に限る。
  - (Major) 「移植直後の過渡が原因で B0 の性質ではない」という対立の立て方は却下 (両方が同時に関わりうる、共通の入力でも初期化後の精度の間の ρE の最大相対差 5.5e-4)。問いを「移植直後の大きな過渡が再上昇に必要か」に変える。
  - (Major) 出力の頻度だけを変えた再実行 (ブリーフの N3) では判別できない → 初期場だけを変える A/B に直す (§6.53)。
  - (Major) ブリーフの「これまでに外したもの」の書き方の訂正: `qAccumulatorFP64` (§6.35) と pRef (§6.45・§6.47) は正式には判別不能。step 末尾の量子化 (§6.41・§6.43) で退いたのは「G0・2,000 step の初期の温度差を状態の丸めだけで再現できる」という十分性。plan 本文の記述は元からこの範囲。
  - やらない: 12 万 step の N1 を今すぐ出す、N4 を原因の近道と扱う、既定値の変更・opt-in の削除・G0c の生産採用。N2 (上限の中間値) は生産設定の候補探しとして有用だが今回の問いを分けないので後。

### 6.53 事前登録: G0c・B0 で初期場だけを変える A/B (2026-10-12、codex diagnose、実行前)

- **問い**: G0c・B0 (上限 50) の float の残差の再上昇 (§6.52 の結果: step 2,000 で 8e-6 → 5,000 で 2e-5) に、最近傍の移植の直後の大きな過渡は必要か。第 1 仮説 (確度 中): 必要でない (残差の低い継続場からでも 6,000 step 内に再発する)。第 2 仮説 (確度 低): 移植・初期化の過渡と B0 の相互作用が早い再上昇を誘発した。第 3 仮説 (確度 低、§6.51 について): 上限 1 は差の発達を遅らせただけ (この試験では扱わない)。
- **腕** (AWS B、変えるのは初期場だけ。G0c・B0・上限 50・バイナリ (`~/forge-fgeom7-{fp64,f32}`)・境界条件・環境変数は §6.52 と同じ。設定は `run_0717_c50_d1` から `nStepOuter` 6,000・`outStepInterval` 250 の 2 行だけを変えたものとバイト一致):
  - 旧起点: `_g0c/q32_init.h5` (sha 6f8fcda6d97482e3)。
  - 継続場: `run_0725_g0c_d1/res_30000.h5` (G0c の FP64、30,000 step、rms_ro 7.5e-7。「収束場」とは呼ばない) の 9 保存量を `restart_field.py --keep-src-dtype` (同一メッシュ) で `_g0c/q32_init.h5` の複製に写し、float32 に丸めて double に戻した `_ic/cont_q32.h5` (両精度に共通、sha は `ric.log` に記録)。
  - 各起点で FP64 ×2・float ×2、計 8 本、各 6,000 step・250 ごと。4 本ずつ 2 回、各回に両起点・両精度を 1 本ずつ (起動の順は旧 FP64 → 継続 FP64 → 旧 float → 継続 float、30 秒ずつ)。
  - 1 回目: `case/45.isobutane_m6_d155/run_0729_ic_old_d1`・`run_0730_ic_cont_d1`・`run_0731_ic_old_f1`・`run_0732_ic_cont_f1`。2 回目: `run_0733_ic_old_d2`・`run_0734_ic_cont_d2`・`run_0735_ic_old_f2`・`run_0736_ic_cont_f2`。
  - 台本 `ric.sh`、判定 `ric_an.py` → `_band_ab/cold_pair/ric_judge.json`。
- **指標**: `residual_history.csv` の rms_ro の窓平均。V1 = 4,000〜5,000、V2 = 5,000〜6,000。閾値は §6.52 の記録をもとにした**この試験用の値** (codex の提案): 高 = 1.5e-5、低 = 3e-6。rms_ro は既知の再上昇の検出に使うだけで、収束の判定ではない。
- **前提** (欠ければ判別不能):
  1. (a) 旧起点で再上昇が再現: 旧起点の float 2 本が両窓 ≥ 1.5e-5、旧起点の FP64 2 本が両窓 ≤ 3e-6。
  2. (b) 継続場は起動直後の残差が小さい: 継続場の各 run の step 1〜250 の rms_ro の平均 ≤ 0.1 × 旧起点の 4 本の同じ平均の最小。
  3. (c) 継続場で精度を切り替えた直後に大差が無い: 継続場の各 float の step 1〜250 の平均 ≤ 3 × 継続場の FP64 の同じ平均の最大 (超えたら初期化と反復中の増幅を分けられない)。
  4. (d) バイナリ・入力の sha、8 本の設定、各起点・各精度の中の res_0 の 9 保存量のビット一致、完走、残差の有限性。精度の間の res_0 の差は記録。
- **判定**:
  - **結果 A**: 継続場の float 2 本が両窓 ≥ 1.5e-5 かつ継続場の FP64 2 本が両窓 ≤ 3e-6 → 移植直後の大きな過渡を必要条件とする読み (第 2 仮説の限定版) を退ける。
  - **結果 B**: 継続場の 4 本がすべて両窓 ≤ 3e-6 → 「初期場によらず 6,000 step で戻る」という第 1 仮説の予測を退け、初期場への依存を優先する (長期の再発までは否定しない)。
  - それ以外 (中間・反復の間の不一致) は判別不能。
- **記録 (判定外)**: 全残差列の窓平均、各起点の float と FP64 (1 本目どうし) の温度差の最大値・位置 (x/r_t、壁からの節点数)・RMS の 250 step ごとの推移、`check_convergence`。
- **やらないこと**: θ_r の小差や double の要否をこの 6,000 step で判定すること。結果を見てから閾値・窓を変えること。
- **所要**: 4 本同時で 1 回 10〜15 分 × 2 回、判定 10 分 (B)。ディスク 最大 約 28 GB (250 step ごとの場、判定後に途中の場を消す)。

### 6.53 の結果 (2026-10-12、判定 `_band_ab/cold_pair/ric_judge.json`)

- **判定**: **判定不能** (前提 (a) が欠けた)。旧起点の float の再上昇は再現した (V1・V2 で 1.8〜2.3e-5 ≥ 1.5e-5) が、旧起点の FP64 の V1 (4,000〜5,000) が 3.7e-6・3.9e-6 で登録の上限 3e-6 を超えた (V2 は 1.1e-6・1.4e-6)。前提 (b)・(c)・(d) は満たした。閾値は変えない。
- rms_ro の窓平均:

| 起点 | run | step 1〜250 | V1 (4k〜5k) | V2 (5k〜6k) |
| --- | --- | --- | --- | --- |
| 旧起点 FP64 | `run_0729_ic_old_d1` / `run_0733_ic_old_d2` | 1.19e-4 / 1.19e-4 | 3.71e-6 / 3.94e-6 | 1.08e-6 / 1.45e-6 |
| 旧起点 float | `run_0731_ic_old_f1` / `run_0735_ic_old_f2` | 1.19e-4 / 1.19e-4 | 2.25e-5 / 2.06e-5 | 1.83e-5 / 2.29e-5 |
| 継続場 FP64 | `run_0730_ic_cont_d1` / `run_0734_ic_cont_d2` | 7.60e-7 / 7.60e-7 | 1.05e-6 / 1.24e-6 | 1.54e-6 / 8.29e-7 |
| 継続場 float | `run_0732_ic_cont_f1` / `run_0736_ic_cont_f2` | 8.07e-7 / 8.08e-7 | 2.90e-5 / 2.76e-5 | 3.06e-5 / 2.93e-5 |

  継続場の側だけを見ると結果 A の型 (float ≥ 1.5e-5、FP64 ≤ 3e-6) だが、前提が欠けたので結果 A とは書かない。全残差列も同じ向き (継続場の float は ρ・ρu・ρv・ρE・ρY が FP64 の 15〜30 倍、ρk・ρω は同程度)。8 本とも `check_convergence` OVERALL: CHECK FAILURES ABOVE。精度の間の res_0 の最大相対差は ρE で旧起点 5.5e-4・継続場 4.9e-5、ほかは 3e-7 以下。
- **記録 (判定外): float と FP64 (1 本目どうし) の温度差の推移**:

| step | 旧起点: 最大 \|ΔT\| (x/r_t, 壁から) / RMS | 継続場: 最大 \|ΔT\| (x/r_t, 壁から) / RMS |
| --- | --- | --- |
| 250 | 1.0 K (−7.45, 52) / 0.016 K | 1.9 K (−7.28, 58) / 0.028 K |
| 1,000 | 1.6 K (−7.33, 51) / 0.036 K | 3.0 K (−7.45, 51) / 0.069 K |
| 1,250 | 3.3 K (−7.87, 71) / 0.054 K | 4.3 K (−8.41, 69) / 0.095 K |
| 1,500 | 8.8 K (−7.70, 71) / 0.098 K | 10.7 K (−8.08, 70) / 0.163 K |
| 1,750 | 18.5 K (−7.18, 72) / 0.19 K | 29.2 K (−7.79, 70) / 0.35 K |
| 2,000 | 33.4 K (−6.68, 72) / 0.37 K | 62.6 K (−7.35, 70) / 0.78 K |
| 2,500 | 77.9 K (−5.60, 72) / 1.04 K | 136.3 K (−7.78, 67) / 2.27 K |
| 3,000 | 107.9 K (−5.38, 73) / 1.68 K | 194.0 K (−7.74, 64) / 3.05 K |
| 4,000 | 113.6 K (−5.00, 77) / 2.50 K | 200.3 K (−7.58, 64) / 3.98 K |
| 6,000 | 130.8 K (−6.05, 70) / 3.49 K | 277.9 K (−7.58, 59) / 8.12 K |

  両方の起点で、差は step 1,000 まで 1〜3 K で止まり、step 1,250 前後から 250 step ごとに約 2〜2.5 倍で育ち、step 3,000 以降は 100〜300 K で頭打ちになる。場所は縮流部の上流側 (x/r_t −8.4〜−4.8) の壁から 51〜77 節点目 (G0c の 121 節点の数え方)、forge-25 の停滞帯 (境界層の縁・温度の前線) と同じ高さ。継続場では FP64 の残差が起動から 6,000 step まで 1e-6 前後のままなので、差の育ちは FP64 の過渡に乗ったものではない。
- **事後の追加の記録 (登録外、同じ場から後処理だけ)**: 同じ精度の再実行どうしの温度差も、同じ場所・同じ速さで育つ。
  - 継続場の FP64 どうし (`run_0730` 対 `run_0734`): step 250 で 0.016 K (入口)、2,000 で 0.15 K (x/r_t −7.8、壁から 70)、2,500 で 0.74 K、3,000 で 2.6 K、4,000〜6,000 で 7〜9 K (x/r_t −4.8〜−3.5、壁から 75〜77) で頭打ち (RMS 0.13〜0.17 K)。旧起点の FP64 どうしもほぼ同じ (最大 4.6〜9.1 K)。
  - 継続場の float どうし (`run_0732` 対 `run_0736`): 2,000 で 15 K、3,000 で 80 K、4,000〜6,000 で 110〜300 K (RMS 1.2〜2.5 K)。旧起点の float どうしも 90〜170 K。
  - 継続場で、step 1,000〜2,000 の FP64 どうしの差の育ち (250 step あたり 2.5〜2.7 倍) は、float と FP64 の差の育ち (2.1〜2.7 倍) と同じ。場所も同じ (x/r_t −8.6〜−7.8、壁から 67〜71)。同じ時点の大きさは FP64 どうしの約 1/500。
  - 差の量の順 (step 1,250〜2,000、継続場の float 対 FP64): T・ρ・dt_local (壁から 69〜70) と U・k・ω (壁から 58〜77) がそろって育ち、P の相対差は 1e-6〜5e-5 と小さい。
  - **読み (仮説、未検証、諮る)**: 縮流部の上流の境界層の縁に、B0 (上限 50) のもとで小さな乱れを増幅するモードがあり、FP64 の再実行の差 (起源は未確認。GPU の加算順か、単精度で組む左辺か) も float と FP64 の差も同じモードに乗って育つ。FP64 では最大 7〜9 K で飽和して残差は 1e-6 前後、float では 100〜300 K の振れになり残差は 3e-5。§6.51 の上限 1 では FP64 の再実行の幅 N も 0.000 だった。
- **codex diagnose** (`notes/reviews/2026-10-12-ic-result-next-diagnose.md`、ブリーフ `notes/reviews/briefs/2026-10-12-ic-result-next.md`。後処理の追加の記録を足して諮り直した) の採否: 全件採用。
  - 正式判定は判定不能のまま (codex が JSON から再計算: 前提 (a) は float 側は成立・FP64 側だけ不成立、前提 (b) の比 0.00679、(c) 1.063)。「継続場の float 2 本は低い初期残差から再上昇した」は観測として書いてよい。
  - (Major、却下) 「線形の増幅なら FP64 も約 6,000 step 遅れて離れるはずなので、有限の大きさの閾値が必須」という推論: 実効の擾乱の比を機械精度の比 1e9 と置く根拠が無い (FP64 にも単精度の経路がある、精度の間で反復の写像が同じとは限らない)。有限期間の増幅を作業仮説として残し、閾値を必須としない。
  - (Major、要再検証) 「同じモード」「飽和」「2 つの状態」: 継続場の float 対 FP64 の RMS は 3,000 → 6,000 step で 3.05 → 8.12 K (2.66 倍) で頭打ちの根拠にならない。最大差の位置も動く。「近い領域で差の増幅を観測した」までに限る。上の「事後の追加の記録」の読みもこの範囲に限る。
  - (Major、要再検証) 「float の右辺が系統的にずれる」という局在化: 末尾の量子化の陰性の結果は右辺・左辺・原始量の格納・初期化を区別しない。
  - (Minor、採用) 判定器: 必要なファイルが欠けたら判定不能にする、前提の欠けは条件ごとに出す (`ric_an.py` は前提 (a) の FP64 側の欠けを「再上昇が再現しない」と書いていた)。次の判定器 `rqs_an.py` から直した。
  - 次の手は M1 (§6.54)。M2 (上限 1) は G0 で観測済みで増幅の抑制と進行の遅れを分けにくいので後。M3 (最初に離れる量) は因果の順を決めないので原因と読まない。M4 (1 回だけの摂動) は K1・K2 の判別にならない。

### 6.54 事前登録: 継続場からの step 末尾の量子化の十分性 (2026-10-12、codex diagnose、実行前)

- **問い**: G0c・B0 の継続場から float が FP64 から離れる差 (§6.53) を、FP64 のビルドで step 末尾に保存量 9 つを float32 に丸めるだけで再現できるか。第 1 仮説 (確度 中): できない (step 内の演算・格納の精度の差が差を供給し、境界層の縁で増幅される)。第 2 仮説 (確度 低): 継続場では末尾の量子化だけで十分な駆動になる (§6.43 の G0・2,000 step の陰性は当てはまらない)。第 3 仮説 (確度 低): 精度に依存する初期化の差が主な種 (継続場でも res_0 の ρE の最大相対差 4.9e-5)。
- **腕** (AWS B、変えるのは量子化の OFF/ON だけ。陽性対照に同じソースの float を添える。G0c・B0・上限 50、共通の入力 `_ic/cont_q32.h5` (sha 3baafa3d23d31e5e)、設定は `run_0717_c50_d1` から `nStepOuter` 3,000・`outStepInterval` 250 の 2 行だけを変えたものとバイト一致、環境変数は §6.53 と同じ):
  - C (FP64・量子化 OFF): `~/forge-fgeom11-fp64` (commit b0303547、sha 74954e7be5ef38a8)。`case/45.isobutane_m6_d155/run_0737_qs_c1`・`run_0740_qs_c2`。
  - Q (FP64・保存量 9 つを量子化): 同じバイナリ + `FORGE_DIAG_QUANTIZE_STATE=1`・`FORGE_DIAG_QUANTIZE_STATE_CHECK=1,2750,3000` (`_ONLY` は付けない)。`run_0738_qs_q1`・`run_0741_qs_q2`。
  - F (陽性対照、float): 同じ commit b0303547 の既定 (float) のビルド `~/forge-fgeom11-f32` (`rqc_build.sh` で作る。FP64 のビルドとの差は `flowFormat.hpp` の typedef 4 行だけ)。`run_0739_qs_f1`・`run_0742_qs_f2`。
  - 3 本ずつ 2 回 (各回に C・Q・F を 1 本ずつ、C → Q → F の順に 30 秒ずつ)。各 3,000 step・250 ごと。台本 `rqs.sh`、判定 `rqs_an.py` → `_band_ab/cold_pair/rqs_judge.json`。
- **指標**: 評価域 = x/r_t ∈ [−9, −6] かつ壁から 40〜90 節点 (G0c の 121 節点)。窓 W1 = step 2,000・2,250・2,500、W2 = 2,500・2,750・3,000。E_run(W) = 評価域・窓の中の (T_run − T̄_C)² の平均の平方根 [K] (T̄_C は C 2 本の同じ時点の平均)。N(W) = 同じ域・窓の (T_C1 − T_C2)² の平均の平方根 [K] (C の再実行の差。全差を使って保守側に取る)。差の最大 |ΔT| も K で残す。
- **前提** (欠ければ判別不能): 陽性対照 = 両窓で min E_F ≥ 10 N。バイナリ・入力の sha、6 本の設定、完走・有限、C と Q の res_0 の 9 保存量と座標のビット一致、Q の自己点検 (step 1・2,750・3,000 の Q・Q_N・Q_M で表せない値 0、step 1 の化学種の履歴は全 CV か入口のピンを除く内点で 0。step 3,000 の化学種の履歴は次の step が無いので数えない)、C と F に量子化の起動の記録が無いこと、必要なファイルがそろうこと (欠けたら判定不能)。
- **判定**:
  - **結果 A**: 両窓で min E_Q − N ≥ 0.5 × max E_F → 第 1 仮説の「末尾の量子化だけでは不十分」を退ける。この介入の十分性までで、float の実際の主因・同じモードの証明とはしない。
  - **結果 B**: 両窓で max E_Q + N ≤ 0.1 × min E_F → 第 2 仮説をこの起点・3,000 step に限って退ける。第 1・第 3 仮説は残る。
  - 中間・反復の間の不一致・陽性対照の不成立は判別不能。閾値・窓・評価域は事後に変えない。
- **やらないこと**: 最初に離れた量を原因と読むこと、この 3,000 step で θ_r の差や double の要否を判定すること。
- **所要**: 3 本同時で 1 回 5〜7 分 × 2 回 + ビルド 10〜15 分 (B)。ディスク 約 11 GB。

### 6.54 の結果 (2026-10-12、判定 `_band_ab/cold_pair/rqs_judge.json`)

- **前提**: すべて満たした。float のビルド `~/forge-fgeom11-f32` (commit b0303547、sha 1f9ff32232da8823、FP64 のビルドとの差は typedef 4 行)、FP64 74954e7be5ef38a8、入力 3baafa3d23d31e5e、6 本の設定は登録どおり、C と Q の res_0 の 9 保存量・座標はビット一致。Q の自己点検: step 1・2,750・3,000 の Q・Q_N・Q_M で表せない値 0、step 1 の化学種の履歴は全 CV で 244 × 4 個 (入口のピンを除く内点で 0)。陽性対照は両窓で成立 (min E_F ≥ 10 N)。評価域の節点 3,378。
- **判定**: **結果 B** (両窓で max E_Q + N ≤ 0.1 × min E_F)。

| 腕 | E (W1 / W2) [K] | 最大 \|ΔT\| (W1 / W2) [K] | rms_ro の窓平均 (2,000〜2,500 / 2,500〜3,000) |
| --- | --- | --- | --- |
| C (FP64) `run_0737`・`run_0740` | 0.011 / 0.039 (N = 0.022 / 0.077) | — | 8.0e-7 / 8.7〜9.0e-7 |
| Q (FP64 + 末尾の量子化) `run_0738`・`run_0741` | 0.42・0.34 / 1.23・1.47 | 3.8・4.0 / 14.6・16.7 | 6.9〜7.0e-7 / 1.3〜1.5e-6 |
| F (float) `run_0739`・`run_0742` | 18.5・18.3 / 27.3・26.8 | 131・132 / 187・187 | 1.8e-5 / 2.9〜3.0e-5 |

- 6 本とも `check_convergence` OVERALL: CHECK FAILURES ABOVE。C1 と F1 の res_0 の最大相対差は ρE 4.9e-5、ρω 3.0e-7、ほかは 6e-8 以下 (記録)。
- **書いてよい範囲**: G0c・B0 の継続場から 3,000 step では、step 末尾に保存量 9 つを float32 に丸めるだけでは、float の温度差の 2〜5 % (評価域の RMS) しか再現しなかった。第 2 仮説 (継続場では末尾の量子化だけで十分) をこの起点・3,000 step に限って退ける。第 1 仮説 (step 内の演算・格納の精度の差が差を供給する) と第 3 仮説 (初期化の差が種) は残る。
- **記録 (判定外)**: Q の差は FP64 の再実行の差 N の 15〜20 倍あり、W1 → W2 で約 3.5 倍に育つ (最大 |ΔT| 4 → 15〜17 K)。Q の rms_ro も W2 で C の 1.5〜1.7 倍に上がり始めている。3,000 step より後に Q も float と同じように離れるかは測っていない (差の発達が遅いだけの可能性は残る)。
- **有限性の確認範囲 (codex の指摘で補った、2026-10-12)**: `rqs_an.py` は選んだ時点の T だけを見ていた。6 本の残差の必須列 (rms_ro・ρu・ρv・ρE・ρk・ρω・ρY0・ρY1) がそろい全 9,000 行が有限、保存した場 13 時点 × 35 量がすべて有限であることを別に確かめた。
- **codex diagnose** (`notes/reviews/2026-10-12-quant-result-next-diagnose.md`、ブリーフ `notes/reviews/briefs/2026-10-12-quant-result-next.md`) の採否: 全件採用。
  - 結果 B と「この起点・評価域・3,000 step では末尾の量子化だけで差を再現できない」は採用 (codex が JSON から (max E_Q + N) / min E_F = 0.024・0.058 を再計算)。
  - (Major、却下) 「残りの大部分は残差の演算の誤差」という寄与率の読み: Q と F は供給される摂動も反復の経路も違うので、RMS の比を原因の内訳にできない。右辺・途中の格納・陰解法・初期化を候補に残す。
  - (Major、却下) 「左辺は両方のビルドで単精度だから候補から外れる」: `ST=float` でも取り込む γ・dt_local・ρ・速度・音速・Ht・乱流粘性は各ビルドの計算結果で、組み立てた左辺は同じとは限らない (`timeIntegration_d.cu:829`、commit b0303547 でも同じ)。組み立ての精度が同じことと組み立てた左辺が同じことを区別する。
  - (Major、採用) 次の P1 (pRef) に FP64 + pRef の対照 Cp と G0c での閉性の検査を足す (§6.55)。
  - (Minor、採用) 有限性の検査範囲 (上で補った。次の判定器 `rpc_an.py` は全量・全時点・必須列を見る)。
  - やらない: P2・P4 を先に回して 3,000 step 離れなければ生産に使えると判断すること (上限は平均流と SST の擬似 dt を縮めるので、短期の小差では抑制と遅延を分けられない)。P3 の部分 double 化の実装 (対象を選ぶ根拠が弱い)。pRef と `lowMachPrecond`・Thornber・左辺の修正を同時に試すこと。

### 6.55 事前登録: 継続場からの pRef の十分性 (2026-10-12、codex diagnose、実行前)

- **問い**: G0c・B0 の継続場から float が離れる差 (§6.53・§6.54) を、圧力の基準値 `space.pRef` (SLAU の面の圧力・運動量流束・hoop から差し引く。`convectiveFlux_slau_d.inc.cuh:618`、`axisymmetricSource_d.cu:52`) だけで少なくとも半減できるか。試せるのはこの圧力のゲージの変更の十分性。第 1 仮説 (確度 中): 半減できる (G0 の §6.45 で float + pRef の温度差 RMS は 1.5〜1.6 K、無介入 3.7〜4.5 K、正式は判別不能)。第 2 仮説 (確度 中): pRef が変えない原始量・物性・残差・陰解法の内部の精度差と B0 の相互作用が大きく、pRef だけでは半減できない。第 3 仮説 (確度 低): 精度に依存する初期化の差が主な種。
- **腕** (AWS B、変えるのは `space.pRef` の 0 / 5,499,000 だけ。G0c・B0・上限 50・`_ic/cont_q32.h5`・§6.54 と同じバイナリ fgeom11 (FP64 74954e7be5ef38a8・float 1f9ff32232da8823)・環境変数を固定。設定は `run_0717_c50_d1` から nStepOuter 3,000・outStepInterval 250 の 2 行、p の腕はさらに `space: {pRef: 5499000.0, ` だけ):
  - C (FP64) `case/45.isobutane_m6_d155/run_0743_pc_c1`・`run_0747_pc_c2`、Cp (FP64 + pRef) `run_0744_pc_cp1`・`run_0748_pc_cp2`、F (float) `run_0745_pc_f1`・`run_0749_pc_f2`、Fp (float + pRef) `run_0746_pc_fp1`・`run_0750_pc_fp2`。各 3,000 step・250 ごと。4 本ずつ 2 回 (各回に C → Cp → F → Fp、30 秒ずつ)。
  - 関門の診断 run `run_0751_pc_hoop` (FP64 + pRef、1 step、`FORGE_DIAG_HOOP_CLOSURE`)。本番の 8 本には診断の環境変数を付けない。
  - 台本 `rpc.sh`、関門 `rpc_gate.py` → `_band_ab/cold_pair/rpc_gate.json`、判定 `rpc_an.py` → `_band_ab/cold_pair/rpc_judge.json`。
- **指標** (評価域・窓は §6.54 と同じ): E0_i = (T_Fi − T̄_C) の RMS、Ep_i = (T_Fpi − T̄_Cp) の RMS、N = max(RMS(C1 − C2), RMS(Cp1 − Cp2))、G = Cp_i と C_j の全組の差の RMS の最大 [K]。
- **前提** (欠ければ判別不能): 両窓で min E0 ≥ 10 N かつ G ≤ 0.05 × min E0。G0c の実際の面ベクトルと hoop の面積で閉性 g ≤ 1e-12 (§6.44 と同じ式、旧 G0 の結果は流用しない。通らなければ本番を回さない)。バイナリ・入力・設定の照合、完走、残差の必須列がそろい全行が有限、保存した全場の全量が有限、同じ精度の中で pRef の有無によらず res_0 の 9 保存量と座標がビット一致。
- **判定** (この試験用の提案の閾値、過去の判定は変えない):
  - **結果 A**: 両窓で max Ep + N + G ≤ 0.5 × min E0 → 第 1 仮説の有限期間の予測を支持し、「pRef だけでは半減できない」を退ける。
  - **結果 B**: 両窓で min Ep − N − G > 0.5 × max E0 → 第 1 仮説を退ける。圧力に関わる誤差全般を除外する根拠にはしない。
  - 中間・反復の間の不一致・前提の不成立は判別不能。`check_convergence` と評価域の平均温度の `check_quasisteady` を添え、有限期間の効果と収束を分ける。
- **やらないこと**: pRef の成否から「低マッハの補正の不足が真因」と判定すること、pRef と他の介入を同時に入れること、生産の設定を変えること。
- **所要**: 関門 1 分、4 本同時で 1 回 7〜9 分 × 2 回、判定 10 分 (B)。ディスク 約 15 GB。

### 6.55 の結果 (2026-10-12、判定 `_band_ab/cold_pair/rpc_judge.json`、関門 `rpc_gate.json`)

- **関門**: G0c の閉性 g の最大 2.58e-16 (99 % 点 1.21e-16)、1e-12 を超える CV 0 / 540,507 (`run_0751_pc_hoop`、FP64 + pRef、1 step)。
- **前提**: すべて満たした (両窓で min E0 ≥ 10 N、G ≤ 0.05 × min E0、バイナリ・入力・設定・res_0・座標の照合、残差の必須列と全行・全場の全量の有限性、完走)。pRef は 4 本の起動の記録で `'pRef' in 'space': 5.499e+06`。
- **判定**: **判別不能** (W1 は結果 A の条件を満たしたが、W2 は max Ep + N + G = 13.80 K に対し 0.5 × min E0 = 13.63 K で 0.17 K 届かない。結果 B の条件はどちらの窓も満たさない)。

| 窓 | E0 (F 対 C) | Ep (Fp 対 Cp) | N | G | 最大 \|ΔT\| F / Fp |
| --- | --- | --- | --- | --- | --- |
| W1 (2,000〜2,500) | 19.17・18.28 K | 6.31・4.33 K | 0.055 K | 0.043 K | 134・132 / 60・42 K |
| W2 (2,500〜3,000) | 27.34・27.26 K | 13.54・9.51 K | 0.153 K | 0.109 K | 189・195 / 114・94 K |

- rms_ro の窓平均 (2,000〜2,500 / 2,500〜3,000): C 8.0e-7 / 8.7〜8.9e-7、Cp 7.9e-7 / 8.4〜8.7e-7、F 1.7〜1.9e-5 / 2.9e-5、Fp 4.6〜6.0e-6 / 1.3〜1.7e-5。8 本とも `check_convergence` OVERALL: CHECK FAILURES ABOVE。評価域の平均温度の `check_quasisteady` (後半) は 8 本とも STEADY (平均は動かず、差は局所)。
- **書いてよい範囲**: G0c・B0 の継続場から 3,000 step で、pRef は float の評価域の温度差を W1 で 67〜77 %、W2 で 50〜65 % 減らした (FP64 自体は pRef でほとんど動かない、G ≤ 0.11 K)。ただし事前登録の「両窓で半減」には W2 で 0.17 K 届かず、正式は判別不能。pRef を入れた float の差と残差も W1 → W2 で約 2〜3 倍に育っており、pRef は差を小さくするが float の離れを止めてはいない。
- 解釈・次の手は codex に諮る。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-10 | [2026-10-10-architecture-float-state-double-geometry-plan.md](../../notes/reviews/2026-10-10-architecture-float-state-double-geometry-plan.md) | GO-with-changes, C0/M3/m1 | 全件採用: (M1) V0 の評価の経路を段 ① に含める (§4.2c、atomicAdd の前の面の流束、FP64 側は保存した入力で上書きして同じカーネル)、(M2) 精度の取り決めと使う側ごとの表 (§4.2a)、接続の一致は段 ② へ、(M3) closure はデバイスの最終の面ベクトルを double で足し上げ、一様な静止場の非劣化を段 ② の必須に (§4.2b)、(m1) V0 のノルム・マスクと V4 の `--tail 1.0`・量ごとの drift/osc を固定 |
| diagnose | 2026-10-10 | [2026-10-10-float-v6-redesign-diagnose.md](../../notes/reviews/2026-10-10-float-v6-redesign-diagnose.md) | C0/M6/m1 | 全件採用。M1 V6 は精度の比較を主軸にし、Sod の厳密解・TG の保存則・凝縮の収支・未検証の経路の回帰は残す (§6 V6)。M2 case/48 の解釈の範囲を限る (§6.21)。M3 `scalarGradient: lsq` を明示し実効設定を記録 (§6.21)。M4 判定量の系列を正式ツールに渡す抽出器 `c48_prec_series.py` (§6.21)。M5 D・E・τ の判定 (§6.21)。M6 `fl1_judge.py` を結果の前に直し、V4 を再評価 (§6.20 の末尾、判定は同じ)。m1 §6.19 の Q_w の不合格の理由を訂正 |
| diagnose | 2026-10-10 | [2026-10-10-c48-precision-result-diagnose.md](../../notes/reviews/2026-10-10-c48-precision-result-diagnose.md) | C0/M4/m1 | 全件採用。M1 抽出の座標が float32 だった (実測で確認) → §6.23 で共通の double の座標で再集計。M2 判別不能の理由は準定常と D + E の 2 つ (§6.22)。M3 書いてよい範囲を窓の平均の記述に限る (§6.22)。M4 判定器の前提の検査を強める (§6.23)。m1 揺れの要約に例外を明記 (§6.22)。追加の run はしない、je の後処理の A/B を 1 回 (§6.23) |
| diagnose | 2026-10-11 | [2026-10-11-g0-fluctuation-controls-diagnose.md](../../notes/reviews/2026-10-11-g0-fluctuation-controls-diagnose.md) | C0/M4/m1 (うち本 plan に関わるのは M4) | 採用: §6.20 の因果の読み方を外し、数値のフラグだけを記録する (§6.20 の訂正)。ほかは hoop plan §4.13 |
| diagnose | 2026-10-11 | [2026-10-11-float-cause-step2-design-diagnose.md](../../notes/reviews/2026-10-11-float-cause-step2-design-diagnose.md) | C0/M4/m1 | 全件採用。M1 SST の commit に原因を絞らない。M2 V3 の PASS を残差の一致と読まない (作用素の差を候補に残す)。M3 軌跡の診断は共通の Q32 の起点で (§6.26)。M4 一部だけ double にする A/B は後回し。m1 所要の見積もりを直した。解析領域の短縮は後回し (§6.26) |
| diagnose | 2026-10-11 | [2026-10-11-next-step-fluct-and-float-diagnose.md](../../notes/reviews/2026-10-11-next-step-fluct-and-float-diagnose.md) | C0/M4/m1 | 全件採用。M1 凍結の定義を「SST の輸送の更新を止める試験」に直す。M2 差が消えたときの結論を「SST の更新・フィードバックが要った」までに限る。M3 §6.27 を全期間の有限性とツールの実行で確かめ直した (確定)。M4 G1・G1x の揺れを格子の性質と言わない (hoop plan 側)。m1 キーを戻す試験は骨格を採用し後回し (hoop plan 側) |
| diagnose | 2026-10-11 | [2026-10-11-float-after-freeze-diagnose.md](../../notes/reviews/2026-10-11-float-after-freeze-diagnose.md) | C0/M4/m1 | 全件採用。M1 §6.29 の失敗の理由と Δ/δ の書き方を訂正、両仮説とも未確定。M2 系列の有限性・座標・準定常を確かめ直した (`fz_recheck.py`)。M3 次は最初の 1 回の更新の段ごとの診断 (§6.30、§5.1 #17)、double の commit は保留。M4 B − B′ の符号反転を判定に使わない。m1 記録の文を直した |
| diagnose | 2026-10-11 | [2026-10-11-double-solve-trajectory-diagnose.md](../../notes/reviews/2026-10-11-double-solve-trajectory-diagnose.md) | C0/M4 | 全件採用。M1 `implicitSolvePrecision` の範囲は block の系の組み立て (Thomas は既に double、保存は float) → 感度試験と呼ぶ。M2 結論の範囲を限る。M3 再実行の差を含む判定の式に直した。M4 同じ float のビルドの B 対 C を主の対照とし、ハッシュ・設定・初期化後の保存量を前提で確かめる (§6.32) |
| diagnose | 2026-10-11 | [2026-10-11-qacc-axisym-state-precision-diagnose.md](../../notes/reviews/2026-10-11-qacc-axisym-state-precision-diagnose.md) | C0/M4/m2 | 全件採用。M1 軸ピンの案を直す (Q_N を射影する前の運動エネルギーを Qacc にだけ当て、作業用 Q は上書きしない。壁と軸の重なる節点の反例を単体試験に入れる、§5.1 #18)。M2 S/R を感度試験の言い方に限る (§6.34)。M3 旧い A・B・B′ の流用は F の 6 条件と初期化後の照合が通るときだけ、外れたら判定不能・S の効果対ばらつきは B・B′・F の幅で。M4 軸対称の拒否は全面解除せず診断用の環境変数で通す。m1 k = 10 のまま。m2 採用の mask は 2 時点の観測と書く。§6.31 の判定不能の正式な理由は Q0 の不一致 (ブリーフの要約の誤り、plan の記述は元から正しい) |
| diagnose | 2026-10-11 | [2026-10-11-qacc-result-next-step-diagnose.md](../../notes/reviews/2026-10-11-qacc-result-next-step-diagnose.md) | C0/M5/m1 | 全件採用。M1 §6.35 から平均流の丸めを候補から外さない (記録の文言を直した)。M2「3 つの介入がどれも効かなかった」と書かない。M3 step 末尾の丸めは「末尾の 9 保存量の量子化の十分性の試験」に限る。M4 G が離れても必要性、離れなくても演算側とは言わない。M5 異なる介入の run をプールしない、腕ごとの反復の分散で区間を作る (§6.36)。m1 FP64 のビルドにも ST=float の経路がある (設定は変えない)。反復は codex の 3 から 5 に増やした (検出力の見積もり、§6.36) |
| diagnose | 2026-10-11 | [2026-10-11-c48-line-implicit-reproducer-diagnose.md](../../notes/reviews/2026-10-11-c48-line-implicit-reproducer-diagnose.md) | C0/M5/m1 | 全件採用。M1 line implicit 単独の試験でなく設定一式を固定した精度の A/B と読む。M2 §6.21 の判定器は使わず有限期間の再現判定を別に登録 (§6.37)。M3 主判定は符号つきの θ (q_w だけ・逆符号は記録)。M4 線の被覆・接続と方向別 dt を確認・記録。M5 初期化後の精度間の差を別に検査。m1 共通の double の座標を明示し対応を検査。「case/48 は偏りなし」とは言わない |
| diagnose | 2026-10-11 | [2026-10-11-lowmach-origin-next-test-diagnose.md](../../notes/reviews/2026-10-11-lowmach-origin-next-test-diagnose.md) | C0/M5/m1 | 全件採用。M1 §6.38 を「指定の閾値で上流が先に超えた」に直す。M2 生成エンタルピーの前提は誤り (`thermoHrefTemp` で除く)、場の値で ρE の丸めと動圧を比べ直した (§6.38 の訂正)。M3 問いを ρE 単独の十分性に限る (§6.40)。M4 符号と場所を保つ指標・初期化後の一致を前提に (§6.40)。M5「前処理では直らない」と言わない。m1 選んだ量に自己点検を追従させる (§5.1 #20) |
| diagnose | 2026-10-11 | [2026-10-11-pref-and-full-quantize-diagnose.md](../../notes/reviews/2026-10-11-pref-and-full-quantize-diagnose.md) | C0/M3/m2 | 全件採用。M1「主な成分は残差の演算側」と書かない。M2 pRef の陽性を「積・足し込みが主な通り道」と読まない、pRef はこの登録から外す。M3 pRef の判定式の穴 (逆符号の大きな偏りを抑制に数える)。m1 pRef の誤差の見積もりを取り下げ。m2 §6.41 の check_quasisteady の記録の不足を §6.42 で補う。G9 対 C だけを登録 (§6.42) |
| diagnose | 2026-10-11 | [2026-10-11-after-g9-next-diagnose.md](../../notes/reviews/2026-10-11-after-g9-next-diagnose.md) | C0/M3 | 全件採用。結論: 案 A (pRef) を先に、FP64 の同等性を作用素 (閉性 ≤ 1e-12) と有限期間 (a の差 ⊂ [−0.002, 0.002]、温度差 RMS ≤ 0.05 K・最大 ≤ 1 K) の 2 段で確かめてから Fp (§6.44)。M1 Cp と C の a の近さを作用素の同等の根拠にしない。M2 案 B の挿入位置では境界条件・`turbulent_viscosity` が上書きする。M3 G9 対 Gp は 2 因子が同時に変わる。案 B は今回作らない |
| diagnose | 2026-10-11 | [2026-10-11-pref-result-long-run-diagnose.md](../../notes/reviews/2026-10-11-pref-result-long-run-diagnose.md) | C0/M4 | 全件採用。§6.45 は判別不能のまま。M1「直交する差 1.5 K」を訂正。M2 長い試験は旧バイナリの実ファイルと陽性対照 F0 で (§6.46)。M3「2 本とも満たす」は再現性の保証と呼ばない診断基準。M4 判定器の 9 キーの完全一致と照合の記録。短い試験の再試行は挟まない |
| diagnose | 2026-10-11 | [2026-10-11-after-pref-long-diagnose.md](../../notes/reviews/2026-10-11-after-pref-long-diagnose.md) | C0/M4/m1 | 採用 (§6.47 の訂正)。M1「pRef は下流の差を減らさなかった」→「抑制の目標の達成は確認できず大きな差が残った (8〜17 % 小さい)」。M2「別の通り道」は要再検証。M3 forge-25 の観察は「B0 の参照が途中の状態かもしれない」まで。M4 案 P (point に固定した精度の A/B、最大 12 万 step を事前登録) を推す。m1 Q_w の表記。案 S・E・低マッハ補正を同時に始めない。ユーザは最初に案 S を希望したが、節点が 57 % 残る (速さ 1.8 倍) と分かったうえで**案 P を先に選んだ** (§6.48) |
| diagnose | 2026-10-11 | [2026-10-11-point-result-diagnose.md](../../notes/reviews/2026-10-11-point-result-diagnose.md) | C0/M3/m1 | 全件採用。M1「point では同じ解へ向かう」と言わない (§6.50 を訂正)。M2「遅いモードを揺れが押す」は仮説のまま、まず上限だけを変える。M3 C1 は上限 50 の陽性対照を同じ条件で取り直す (§6.51)。m1「1/50」は保存した窓に限る・符号も違う |
| diagnose | 2026-10-11 | [2026-10-11-g0c-precision-test-diagnose.md](../../notes/reviews/2026-10-11-g0c-precision-test-diagnose.md) | C0/M4/m2 | 全件採用 (§6.52)。M1 陽性対照は §6.51 の上限 50 の 4 本。M2 `g0c_prep.sh` を失敗で止める形に。M3 実際の G0 の入力と照合 → 10 桁の .msh の欠陥を見つけて作り直した。M4 増分の差と値の差を分ける。m1 D_c = D_old × S_old / S_c。m2 大差は絶対値で。不足情報: `rcap_an.py` の設定の照合を 8 本に広げた (§6.51) |
| diagnose | 2026-10-12 | [2026-10-12-cap-g0c-result-diagnose.md](../../notes/reviews/2026-10-12-cap-g0c-result-diagnose.md) | C0/M5 | 全件採用 (§6.52 の結果の末尾)。正式判定は変えない。M1 H1 は「方向別 dt の拡大が精度差の発達に関与」まで。M2 H2 は初期場への依存の仮説に限る。M3 H3 の対立の立て方を却下、問いを移植の過渡の必要性に。M4 N3 を初期場の A/B に直した (§6.53)。M5 外した介入の書き方を限定 |
| diagnose | 2026-10-12 | [2026-10-12-ic-result-next-diagnose.md](../../notes/reviews/2026-10-12-ic-result-next-diagnose.md) | C0/M3/m1 | 全件採用 (§6.53 の結果の末尾)。正式判定は判定不能のまま。M1 閾値を必須とする推論を却下。M2 同じモード・飽和・2 状態は言わない。M3 右辺への局在化は言わない。m1 判定器の欠落・条件ごとの出力を直した。次は §6.54 (末尾の量子化の十分性、陽性対照つき) |
| diagnose | 2026-10-12 | [2026-10-12-quant-result-next-diagnose.md](../../notes/reviews/2026-10-12-quant-result-next-diagnose.md) | C0/M4/m1 | 全件採用 (§6.54 の結果の末尾)。結果 B は採用。M1 寄与率の読みを却下。M2 左辺を候補から外す読みを却下。M3 P1 に Cp と G0c の閉性を足す (§6.55)。m1 有限性の検査範囲を補った |

## 7. 影響範囲

- すべてのケースのカーネル (粘性・拡散・再構成・リミタ・陰解法の対角・壁関数)。float のビルドの回帰の値が変わる (案 A の場合)。
- メモリ: 面ごとに float 5〜11 個。

## 変更ログ


- 2026-10-10: 起票 (draft)。ユーザ「float化よろ」。調査と §4・§6 の案。
- 2026-10-10: codex diagnose (`notes/reviews/2026-10-10-float-geometry-design-diagnose.md`) を全件採用して §4・§5・§6 を書き直し。段 ①〜④ に分け、まず V0 (状態を固定した演算の A/B) で因果を判別する。
- 2026-10-10: codex plan 段 (GO-with-changes、C0/M3/m1) を全件採用。§4.2a (取り決めと表)・§4.2b (closure)・§4.2c (V0 の評価の経路) を追加し、V4 の判定ツールの引数を固定。
- 2026-10-10: 段 ① を実装・確認 (§6.2)。V0 は結果 A (座標の差の丸めが壁際の誤差のほぼすべて: ω の拡散で旧 10 % → 新 8e-8)。段 ②〜④ へ進む。
- 2026-10-10: §6.3 に段 ② の合格条件を事前登録した (接続の一致、係数の精度、closure の非劣化、FP64 の同一性)。
- 2026-10-10: 段 ② を実装・測定 (§6.5)。事前登録の判定は case/48 の cInt だけ FAIL (真値 0 の成分を分母に拾う物差しの欠陥)。closure は case/45 の既定では使われていなかったので 3' を足した (PASS)。codex diagnose を全件採用: 旧基準の FAIL を残し、改めた基準 (4a)・closure を使う経路の FP64 の新旧比較 (4b)・判定の台本の修正 (4c)・段 ③ の V2 の追加 (4d) を登録。段 ③ を実装 (未測定)。
- 2026-10-10: 段 ② の判定をやり直した (`fg2_rejudge.py`)。改めた基準 4a・closure を使う経路 4b・台本の修正 4c ですべて PASS。段 ② を完了とした。
- 2026-10-10: §6.6 に段 ④ の合格条件を事前登録した。
- 2026-10-10: 段 ③ を実装・測定 (§6.7)。FP64 の面流束はビット一致、FP64 isp 1 と float case/56 の同一性・1 step の残差と dq は再実行の差の範囲で PASS。float の case/45 では 20 step で `rms_roOmega` が 6.6e3 → 1.5e3。
- 2026-10-10: §6.8 に V3 の実施形を事前登録した (段 ④ の結果を見る前)。
- 2026-10-10: 段 ④ を実装・測定 (§6.9、すべて PASS)。V3 を回した (§6.10、PASS: 座標の差に由来する誤差は桁で減り、残りは粘性の面の流束以外)。
- 2026-10-10: codex diagnose (`notes/reviews/2026-10-10-float-geometry-v3-v4-diagnose.md`) を全件採用。V3 の内訳の主張を撤回し (§6.10)、§6.12 に入力の丸めの A/B を事前登録。V4 を §6.11 で事前登録し直した (停止規則のもとでの比較、旧の窓の条件は撤回)。§4.4 のメモリの見積もりを訂正。
- 2026-10-10: §6.12 を回した (§6.13): 棄却 (拡散の入力の丸めは ω の残差の差を説明しない、r = 1.01)。V3 の集計器を直して集計し直した (PASS のまま)。見張り `m9_watch.py` に `--interval`・`--min-points`・`--keep-tail` を足し (合成の系列で、3 点は DATA_ERROR・9 点 + 2 出力連続で到達・8 点は未到達を確認)、V4 (§6.11) を投入した。
- 2026-10-10: V6 の下調べ (既存の格子は 5 ケースとも float32) を記録し、§5.1 #13 (格子の HDF5 の精度) を足した。
- 2026-10-10: ユーザ決定で §4.6 (変換器は常に double で幾何を書く・ソルバの警告・手順) を足し、§6.14 に合格条件を事前登録した。
- 2026-10-10: codex diagnose (`notes/reviews/2026-10-10-converter-double-geometry-diagnose.md`) を全件採用して §4.6 と §6.14 を補強した (案 (b)、double の正本を変換器専用に持つ、壁距離の切り出し、警告は幾何のデータセット全体を見る、参照の固定と比較器の条件)。§2 のスコープを同期した。
- 2026-10-10: §4.6 (変換器は常に double で幾何を書く・ソルバの警告・手順) を実装し、§6.14 をすべて PASS (§6.15)。
- 2026-10-10: V4 を回した (§6.16)。float は到達したが、θ_r +1.3〜1.8 %・Q_w +1.4 % でフラグ。今のままでは case/45 の設計の計算に使えない。
- 2026-10-10: V5 PASS (§6.17、比 0.636)。V4 の解釈を codex に諮り全件採用 (§6.18): case/45 の設計は当面 FP64、精度の切り替えの A/B を事前登録。S_lostfrac の記述を訂正。
- 2026-10-10: §6.18 の A/B の 1 回目は準備で止まった (`restart_field.py` が float の res の float32 の座標と double の格子を完全一致で照合していた)。`restart_field.py` を、SRC の精度が低いときは DST を SRC の精度に丸めて照合するよう直し (同じ格子なら一致する)、回し直した。
- 2026-10-10: §6.18 の A/B を回した (§6.19): 判別不能 (θ_r は FP64 に切り替えると戻り、float では留まる。Q_w は揺れで条件を外れた)。
- 2026-10-10: ユーザ決定: float が違う所で止まる原因を調べる。hoop の修正を入れてから (§5.1 #15)。
- `2026-10-10` — §6.20 に #15 (1) (hoop の修正の格子での V4 のやり直し) を事前登録し、float の腕を投入した (FP64 の腕は hoop plan の `run_0483_hp7_k1`)。
- `2026-10-10` — ユーザ決定: r0・r1 のメモリの増加を受け入れる (§5.1 #12)。
- `2026-10-10` — codex diagnose (V6 の改訂、C0/M6/m1) を全件採用。§6.21 に case/48 の float と FP64 の比較を事前登録し、`fl1_judge.py` を結果の前に直した (V4 の再評価は同じ判定)。§6.19 の Q_w の理由を訂正。
- `2026-10-10` — 粘性 Jacobian plan の返却を反映した (§5.1 #15 の先行の証拠、#16 の条件)。
- `2026-10-10` — §6.21 (case/48 の float と FP64 の比較) を回した (§6.22): 判別不能 (float の時間の揺れ 1e-4 が準定常の閾値を超えた)。解釈は諮問中。
- `2026-10-10` — §6.21 の結果の解釈を codex に諮り全件採用 (§6.22 の訂正)。§6.23 に再集計と je の後処理の A/B を事前登録した (計算 0 step)。
- `2026-10-10` — §6.23 を回した (§6.24): 再集計も判別不能 (座標の精度は判定を変えない)、je の A/B は支持 (B/A 0.074)。case/48 は止める。
- `2026-10-11` — §6.20 の解釈の区分を結果の前に訂正した (因果の読み方を外す)。
- `2026-10-11` — §6.20 を判定した (§6.25): 到達・フラグあり、差は V4 と同じ向き・同じ桁。
- `2026-10-11` — #15 (2) の設計を codex に諮り全件採用。§6.26 に共通の Q32 の起点からの軌跡の A/B を事前登録した。
- `2026-10-11` — §6.26 を回した (§6.27): 精度依存を支持。同じ状態から float だけが以前の float の止まり所へ戻る。場では k・ω の差が先に・大きく出る (記録)。
- `2026-10-11` — §6.27 を確かめ直して確定。§6.28 に SST の輸送の更新を止めた軌跡の A/B を事前登録した (codex diagnose)。
- `2026-10-11` — §6.28 を回した (§6.29): 判別不能。SST を止めても float は D の向きに 0.6〜1.2 % 動く (記録)。解釈は諮る。
- `2026-10-11` — §6.29 の解釈を codex に諮り全件採用 (訂正と確かめ直し)。§6.30 に最初の 1 回の更新の段ごとの診断を事前登録し、実装を §5.1 #17 に入れた。
- `2026-10-11` — §5.1 #17 (段ごとの更新の診断 `FORGE_DIAG_STAGE_DUMP`) を実装した (implementer)。float・FP64 ともビルド可、`update_d.cu` の SASS は HEAD と一致。§6.30 の記録の時点を実装に合わせて訂正。
- `2026-10-11` — §6.30 を回した (§6.31): 判定不能 (初期化で ρE が精度ごとに組み直され Q0 がそろわない)。登録外の記録では 1 step の差は再実行の揺れと同じ桁で、段を分けられない。
- `2026-10-11` — §6.32 に block の系の組み立ての精度の感度試験 (`implicitSolvePrecision` 0 → 1) を事前登録した (codex diagnose)。
- `2026-10-11` — §6.32 を回した (§6.33): 判別不能 (抑制 0/6、抑制不十分 5/6)。組み立てを double にしても float はほぼ同じように離れる (記録)。
- `2026-10-11` — §6.34 に平均流の状態を FP64 で積む float の軌跡 (qAccumulatorFP64 の感度試験) を事前登録した (codex diagnose、C0/M4/m2 全件採用)。
- `2026-10-11` — §5.1 #18 を実装した (implementer。qAccumulatorFP64 の軸対称の診断用の経路・採用の mask の診断・単体試験 `test_qacc_axis` (a)〜(g) PASS)。§6.34 の前提 1 に nvcc のコード生成の揺れの扱いを追記 (run の前)。
- `2026-10-11` — §6.34 を回した (§6.35、AWS A): 判別不能 (抑制 0/6、抑制不十分 0/6。後者は 6 組とも再実行の差の 10 倍の条件だけで外れた)。平均流の状態を FP64 で積んでも float は B とほぼ同じに離れた (98〜104 %)。
- `2026-10-11` — §6.35 の記録を codex の範囲に直し、§6.36 に FP64 のビルドで step 末尾に状態を float32 に丸める十分性の試験を事前登録した (各 5 反復、Welch の区間判定)。§5.1 #19 を追加。
- `2026-10-11` — §6.37 に case/48 で case/45 の陰解法の設定一式を入れた float と FP64 の有限期間の再現判定を事前登録した (ユーザ決定: 小さいケースで再現させる。§6.36 は保留)。
- `2026-10-11` — §5.1 #19 (状態の量子化の診断) を実装した (implementer。単体試験 PASS、SASS は新しいカーネル以外 HEAD と一致)。§6.36 は保留のまま。
- `2026-10-11` — §6.38 に既存の軌跡の出力から float と FP64 の差の発生場所を記録した (入口〜収縮部の低マッハ域で生まれて下流へ広がる。仮説で原因とは書かない)。
- `2026-10-11` — §6.37 を回した (§6.39、AWS A、case/48): 結果 B (再現せず)。θ の float と FP64 の差は 3 断面・両窓で |e| + E ≤ 9.2e-5。
- `2026-10-11` — §6.38 を codex の指摘で訂正 (場の値で ρE の丸めと動圧を比べ直し、前の見積もりの誤りを記録)。§6.40 に ρE だけを丸める短い試験を事前登録、§5.1 #20 を追加。
- `2026-10-11` — §5.1 #20 (`FORGE_DIAG_QUANTIZE_STATE_ONLY`) を実装した (implementer。単体試験 PASS、デバイス側の SASS は HEAD と一致、未設定では #19 と同じ)。
- `2026-10-11` — §6.40 の実行の場所を B から A (1 本ずつ) に変えた (run の前、B が満杯のため)。B での単体試験は PASS (台本が PASS の行の書式を読み違えていたので記録を訂正)。
- `2026-10-11` — §6.40 を回した (§6.41、AWS A): 結果 B (ρE 単独の量子化では入口〜収縮部の初期の温度差の主な成分は出ない、float の差の RMS の約 4 %)。
- `2026-10-11` — §6.42 に保存量 9 つ全部の末尾の量子化 (G9 対 C) の短い試験を事前登録した。pRef の試験は外して別に登録する (codex)。
- `2026-10-11` — §6.42 を回した (§6.43、AWS A): 結果 B (9 保存量の末尾の量子化だけでは不十分、温度差の RMS は float の約 7 %)。
- `2026-10-11` — §6.44 に pRef の短い試験を事前登録した (FP64 の同等性の 2 段の関門の後に float)。
- `2026-10-11` — §6.44 を回した (§6.45、AWS A): 関門は通過、判定は判別不能 (L_d 0.494 が 0.5 にわずかに届かない)。記録として pRef で射影 e は 0.91 → 0.09、温度差の RMS は約 6 割減。
- `2026-10-11` — §6.45 の文言を訂正、§6.46 に pRef の長い試験 (30,000 step、F0・Bp ×2・Ap) を事前登録した。
- `2026-10-11` — §5.1 #21 にノズル上流の低マッハ補正の検討を追加した (ユーザの見立て、方針は未決)。
- `2026-10-11` — §5.1 #22 に forge-25 の観察 (B0 が細かい x の格子で途中の状態に止まる疑い) と float の調査への含意を記録した。
- `2026-10-11` — §5.1 #16 を更新 (forge-8a: 粘性 Jacobian の dt plan は分岐 2 で閉じた、既定は決めない (同日の訂正)。再開の条件は残差の評価の精度を上げた経路が既定になること)。
- `2026-10-11` — §6.46 を回した (§6.47、AWS A): 判別不能 (陽性対照 F0 が W2 で旧 B・B′ の幅を超えた)。記録として pRef の下でも float の下流の θ_r は FP64 から +0.9〜1.9 % 離れた (Bp は全 6 組で結果 B の条件)。
- `2026-10-11` — §6.47 の記述を codex の指摘で訂正 (抑制の目標は未達で大きな差が残った、8〜17 % 小さい。「別の通り道」は要再検証)。
- `2026-10-11` — §6.48 に point の陰解法に固定した float と FP64 の比較 (案 P) を事前登録した (ユーザの選択)。
- `2026-10-11` — §6.48 の run を AWS B で始めた (08:20 UTC、B は容量不足で 3 回目の起動の試行で起動)。判定器の設定の照合に nStepOuter の必然の差を足した (判定の前)。
- `2026-10-11` — §5.1 #22 に forge-25 の続報 (Gref で B0 の止まる値は point の収束値から θ_r −0.6 %・Q_w −0.5 %) を記録した。
- `2026-10-11` — §6.48 の段 1 (30,000 step) を記録 (§6.49): 判別不能で段 2 へ持ち越し。記録として point では float と FP64 の θ_r が 0.001 ポイントまで一致。
- `2026-10-11` — §5.1 #22 に forge-25 の続報 2 (方向別 dt の上限 1 で B0 の停滞が解け始める) と、float の次の候補 (上限 1 での float と FP64) を記録した。
- `2026-10-11` — §5.1 #22 に forge-25 の続報 3 (B0 の停滞は縮流部の主流側の遅いモード、上限は右辺を変えない) と float の調査との対応 (仮説) を記録した。
- `2026-10-11` — §6.48 の段 2 (120,000 step) を記録 (§6.50): 判別不能 (差の条件は満たしたが収束・準定常の付帯条件を満たさない)。point では float と FP64 の θ_r の差は 0.03 ポイント程度。
- `2026-10-11` — §6.50 の言い方を訂正、§6.51 に上限 50 対 1 の float と FP64 の比較 (案 C1、8 本) を事前登録した。
- `2026-10-11` — §5.1 #22 に forge-25 の続報 4 (縮流部の流れ方向を粗くするだけで B0 の停滞が消える) と float の次の候補 (G0c で B0 の float と FP64) を記録した。
- `2026-10-11` — ユーザ決定: G0c は float の診断試験に限り AR ≤ 10,000 を許す (§5.1 #22)。
- `2026-10-11` — §6.51 の 1 回目の後に台本が diff の終了コードで止まった件を記録、2 回目は `rcap2.sh` で続ける。G0c の .msh を `make_g0c.py` で作った (区間 589 → 337 列、区間の外は G0 とビット一致)。
- `2026-10-11` — G0c の試験の設計を codex に諮った (`notes/reviews/2026-10-11-g0c-precision-test-diagnose.md`、全件採用)。実際の G0 の入力との照合で `make_g0c.py` の .msh が 10 桁だった欠陥を見つけ、17 桁で作り直した (照合 OK)。§6.52 に G0c で B0 のまま float と FP64 を比べる試験 (`run_0725〜0728`) を事前登録した。§6.51 の判定器の設定の照合を 8 本に広げた。§5.1 #22 に forge-25 の訂正 (遅いモードの場所は境界層の縁) と forge-8a の低マッハの観測を記録した。
- `2026-10-12` — §6.51 の結果 (結果 A) と §6.52 の結果 (判別不能、G0c では FP64 は停滞を抜けるが float は 2e-5 へ戻る) を記録し、まとめて codex に諮った (`notes/reviews/2026-10-12-cap-g0c-result-diagnose.md`、全件採用)。§6.53 に G0c・B0 で初期場だけを変える A/B (`run_0729〜0736`) を事前登録した。
- `2026-10-12` — §6.53 の結果 (判定不能、継続場からも float だけが再上昇、差は縮流部の上流の境界層の縁で育つ。同じ精度の再実行どうしの差も同じ場所で育つ) を記録し、codex に諮った (`notes/reviews/2026-10-12-ic-result-next-diagnose.md`、全件採用)。§6.54 に継続場からの末尾の量子化の十分性の試験 (`run_0737〜0742`) を事前登録した。
- `2026-10-12` — §6.54 の結果 (結果 B、末尾の量子化は float の差の 2〜5 %) を記録し、codex に諮った (`notes/reviews/2026-10-12-quant-result-next-diagnose.md`、全件採用)。§6.55 に継続場からの pRef の十分性の試験 (`run_0743〜0751`) を事前登録した。§5.1 #22 に forge-25 の続報 10 を記録した。
