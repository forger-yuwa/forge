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
  - 変換器の float 経路の修正。冷却壁の格子は従来どおり FP64 の変換器 (`stod`・17 桁) で作る。
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
- 速さの目標: float の新しい経路で、FP64 の 1 step の 0.70 倍以下。採用の判断には、単価だけでなく到達までの総時間も含める (codex V5)。

### 4.5 状態の commit の精度

- 2026-10-10 の測定 (FP64 の 9 → 10 step の変化と float の刻みの比) は、float で要求された更新と実際に反映された更新の測定ではない (codex Major: H2 は要再検証)。float の run で `dq_requested` と `Q_after − Q_before` を別々に記録する (§6 V4)。
- `qAccumulatorFP64` の軸対称対応は別 plan とする (codex 採用)。ただし、V4 が不合格だったときに、すぐ commit のせいにはしない。軸対称の拒否だけを外すこともしない (拒否には実装上の理由がある、`main.cpp:3371`)。

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
| 3 | 段 ① と V0・V1 | 段 ① を実装する (§4.2a の取り決め、本番の配列は変えない、§4.2c の評価の経路)。合格: float と FP64 のビルドが通る。FP64 で e が今の差とビット一致。診断の旧腕の和が変更していないカーネルの寄与とビット一致。V1 の係数の検査 (接続は段 ②)。そのうえで V0 を回し、結果 A なら段 ②〜④ へ、B か判別不能なら上位に諮る | O → F |
| 4 | 段 ② | 合格: V1 (接続の一致)・V2 | O |
| 5 | 段 ③ | 合格: V2 | O |
| 6 | 段 ④ | 合格: V2 | O |
| 7 | float の同じ状態の残差 | V3 | O |
| 8 | 到達 | V4 | O |
| 9 | 速さと総時間 | V5 | O |
| 10 | 標準ケースの回帰 | V6 (5 ケース) | O |
| 11 | 状態の commit の精度 (条件付き) | V4 の記録で停滞が見えたら、`qAccumulatorFP64` の軸対称対応を別 plan に (拒否の理由に向き合う) | F |

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
- **V4 到達** (float の新しい経路、case/45 の B0 の構成、`run_0183` の res_100000 から、2 出力連続の水準まで、2500 step ごと、最大 20 万 step):
  - θ_r (断面ごと) と Q_w (符号つき) を、両方の run の同じ設定の区間の末尾 2 万 step (9 点以上) で比べる。分母は参照の窓の平均。
  - 両方の窓のドリフトと振幅が閾値 (θ_r 0.05 %・Q_w 0.1 %) の 1/5 以下で、平均の差に両方の時間変動の幅を足しても閾値以内。単調なら漸近値も書く。
  - 抽出した系列 (末尾 2 万 step に切り出した窓) を `check_quasisteady.py --series-csv --tail 1.0` に渡す。θ_r は `--drift 0.0001 --osc 0.0001`、Q_w は `--drift 0.0002 --osc 0.0002` を明示する (閾値 0.05 %・0.1 % の 1/5。既定の drift 0.05 は 5 % で緩すぎる、codex plan 段 m1)。
  - `check_convergence` が NOT CONVERGED なら、結論は「指定した水準への到達の比較」に限る。PASS・対象量の STEADY・判定区間がそろわない限り、「FP64 と同等に収束」とは書かない。
  - FP64 側の参照の窓: B0 の到達の場は消えているので、FP64 の B0 の参照の run は同じ記録の形 (2500 step ごと、末尾 2 万 step の場を残す) で回し直す。
  - 記録: `dq_requested` と `Q_after − Q_before` の統計 (§4.5)、壁の熱流束の分布の差、step 数。
- **V5 速さ**: 専有・交互に 1000 step × 3 本。新しい float が FP64 の 0.70 倍以下。採用の判断は V4 までの総時間でも行う。
- **V6 標準ケースの回帰** (`case/08` は使わない。`procedures/verification/README.md:23` のとおり回帰の基準にならない):

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

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-10 | [2026-10-10-architecture-float-state-double-geometry-plan.md](../../notes/reviews/2026-10-10-architecture-float-state-double-geometry-plan.md) | GO-with-changes, C0/M3/m1 | 全件採用: (M1) V0 の評価の経路を段 ① に含める (§4.2c、atomicAdd の前の面の流束、FP64 側は保存した入力で上書きして同じカーネル)、(M2) 精度の取り決めと使う側ごとの表 (§4.2a)、接続の一致は段 ② へ、(M3) closure はデバイスの最終の面ベクトルを double で足し上げ、一様な静止場の非劣化を段 ② の必須に (§4.2b)、(m1) V0 のノルム・マスクと V4 の `--tail 1.0`・量ごとの drift/osc を固定 |

## 7. 影響範囲

- すべてのケースのカーネル (粘性・拡散・再構成・リミタ・陰解法の対角・壁関数)。float のビルドの回帰の値が変わる (案 A の場合)。
- メモリ: 面ごとに float 5〜11 個。

## 変更ログ

- 2026-10-10: 起票 (draft)。ユーザ「float化よろ」。調査と §4・§6 の案。
- 2026-10-10: codex diagnose (`notes/reviews/2026-10-10-float-geometry-design-diagnose.md`) を全件採用して §4・§5・§6 を書き直し。段 ①〜④ に分け、まず V0 (状態を固定した演算の A/B) で因果を判別する。
- 2026-10-10: codex plan 段 (GO-with-changes、C0/M3/m1) を全件採用。§4.2a (取り決めと表)・§4.2b (closure)・§4.2c (V0 の評価の経路) を追加し、V4 の判定ツールの引数を固定。
