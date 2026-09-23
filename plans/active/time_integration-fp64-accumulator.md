# 保存量アキュムレータの倍精度化 (commit の丸めで定常解に到達できない問題)

## メタ

- **area**: `time_integration`
- **status**: `draft`
- **related_docs**:
  - [`methods/time_integration/implementation.md`](../../methods/time_integration/implementation.md) の「commit の丸め — 定常解の到達限界を決める」
- **related_plans**:
  - [`plans/active/case-hypersonic-gap-heating-validation.md`](case-hypersonic-gap-heating-validation.md) §4.7-6d-2 (発見の経緯と実測)
- **created**: `2026-09-23`
- **owner**: gap-heating セッション

## 1. 目的

定常陰解法の commit `Q += dq` が、$|dq| < \tfrac12\,\mathrm{ULP}(Q)$ になった時点で丸めに落ち、
**残差が残っているのに場が動かなくなる**。保存量のアキュムレータだけを倍精度にして、
この到達限界を実用上の問題でなくする。

## 2. スコープ

**含む**: 保存量 5 変数 (`ro`, `roUx`, `roUy`, `roUz`, `roe`) の**保持と commit の精度**。

**含まない**:

- 流束・勾配・リミッタ・ヤコビアンの計算精度 (float32 のまま)
- 幾何 (`geom_float`) の精度 — **無関係と実証済み** (§3 の (c))
- 化学種 `roY_s` / 乱流 `roK`, `roOmega` / 凝縮モーメント — 同じ問題を持ちうるが、本計画では扱わない
  (同じ機構なので横展開は容易。まず流れ 5 変数で効果と副作用を確かめる)
- `unsteady`/dual-time 経路 — 定常経路で確立してから

## 3. 関連 docs と前提

前提は `case/56.gap_tp1187` で実測した (詳細は
[case plan](case-hypersonic-gap-heating-validation.md) §4.7-6d-2、run パスもそこ)。

- **(a) 機構**: `update_d.cu:249` の `ro[ic] = roN[ic] + d0` で `d0` = `dq_block_old_0` は
  1 step の密度変化そのもの。深部 ($z/W>10$, 16607 CV) では $\langle\rho\rangle$=0.017755 に対し
  1 ULP = 1.86e-9、$\langle|dq|\rangle$ = 2.97e-10 = **0.159 ULP**。
  **全 CV が $\tfrac12$ ULP 未満**で、float32 加算を模擬すると動く CV は 0.03 %、
  実効 $\langle d\rho\rangle$ は意図の **0.3 %**、実測の傾きと比 0.6。
- **(b) 症状**: 定常解ならゼロであるべき $\lvert\dot m\rvert = |\int\rho U_y dx|$ が
  **べき乗則** $\propto\mathrm{step}^{-0.23}$ でしか減らない (step 倍加で 15 %)。
  同じ場を倍精度ビルドで継続すると**幾何級数** (25k step ごとに 1/2.81) になり、
  300k step で 4.54e-7 → **4.56e-12** (SU2 の床 8.59e-12 を下回る)。
- **(c) 交絡のうち 2 つは否定した。ただし「commit 丸めが唯一の律速」は未証明** (codex plan M2)。
  `run_0020_double` は `flow_float` と `geom_float` を**ともに**変えているので、状態保持だけでなく
  流束・EOS・残差・線形解法・SST 更新も変わっている。**閉性誤差 0 は、非一様流の流束評価・総和・勾配の
  丸めが 0 であることを意味しない**。言えるのは「commit 丸めは実在する」までで、「それを直せば十分」は未検証:
  - **幾何精度**: 深部 CV の **99.8 % が双対面の閉性誤差ちょうど 0** (中央値・90 %・99 % 点すべて 0)。
    閉性由来の偽の質量ソースは 0 で、実測の $\langle res\_ro\rangle$ = 3.42e-11 を説明しない。
  - **ブロックサイズ**: 倍精度ビルドは `FORGE_CUDA_BLOCKSIZE` 512→128 を要したが、
    float32 のまま 128 にしても 100k step で 1.40 % 減 (対照 512 は 1.88 %)。倍精度の 50k で 87.5 % とは別世界。
- **(d) 線形解法側ではない**: `nStepInner` 4→16 で減衰倍率 1.00 / 1.18 / 0.83 (**±20 % 以内**)。
  `lineImplicit: 1` はライン構築に成功 (643 本・CV の 92.5 %・最大長 107) するが、
  **最大長 107 は幅方向の列数 105** で、止まっている軸モード (191 セル) に届かず、
  壁時計あたり 1 桁悪化した。どちらも $dq$ を精緻にするだけで、その $dq$ が表現できない。
- **(e) スケール変更は無効**: ULP は値に比例するので、保存量を定数倍しても $|dq|/\mathrm{ULP}$ は不変。

## 4. 設計方針

### 4.1 採る案: アキュムレータのみ倍精度 (案③)

**保存量 5 変数を `double` で保持し、commit を倍精度で行う。** 流束・勾配・ヤコビアンは float32 のまま。

| | 実効 ULP ($\rho$=0.0178) | 解像できる $\lvert\dot m\rvert$ の下限 | SU2 の床 8.6e-12 に対し |
| --- | --- | --- | --- |
| 現状 float32 | 1.86e-9 | 2.85e-6 | **不足** (実測と一致) |
| 案② Kahan (float32) | 2.78e-17 | 4.24e-14 | 200 倍下 |
| **案③ FP64 アキュムレータ** | **3.47e-18** | (5.30e-15) | — |

> **注意 (codex plan M2)**: 右 2 列は**表現分解能から引いた上界**であって、到達できる収束床の保証ではない。
> 状態を保持できても **FP32 の残差評価に別の限界が残る**。実際 `run_0020_double` (全体 FP64) でも
> 最終的に 1e-11 台で頭打ちになった (§6 の G2 参照)。**この数値を合格条件に使わない**。

**コスト**: FP64 加算は 1 step あたり $N_{cell}\times5$ 回だけ。65,194 CV で 3.26e5 回、
RTX 3060 の FP64 (FP32 の 1/32、≈0.2 TFLOPS) で **1.63 µs/step**。現状 1.63 ms/step に対し **+0.1 %**。
律速は流束計算 (数千 flop/CV) で、そこは float32 のまま動く。
**「倍精度が速い GPU」は不要** — 全部を FP64 にする場合の話であって、本計画はそうしない。

### 4.2 案② (Kahan 補償加算) を採らない理由

性能要件としては案② でも足りる (上表)。採らないのは:

- `c = (t - ro) - y` は代数的には 0 に見えるので、**結合則を許す最適化 (`-ffast-math` / `-use_fast_math`) で
  丸ごと削除される**。黙って無効化され、しかも症状は「なんとなく収束が遅い」なので気づきにくい。
- 防ぐには `volatile` か関数単位の最適化抑制が要り、**後から読む人が理由を追えない形**になる。
- コスト差 (+0.1 % vs ~0 %) はこのリスクに見合わない。

ただし **FP64 が無い/極端に遅い環境へ移す可能性が出たら案② を再検討する**。その判断材料として
上表と本節を残す。

### 4.3 状態の正本をどこに置くか — **影アキュムレータ** (2026-09-23 確定、`diagnostician`)

> **撤回 1**: 初稿の「`Q` の保持型だけ `double` にし、読み出しは float32 へ落とす」は**効かない**。
> commit は `Q += dq` でなく **`Q = Q_N + dq`** で、`Q_N` は step 末尾に `roN[ic] = ro[ic]` で
> **float32 コピー**される (`update_d.cu:39`)。実測: $dq$=2.97e-10 を 100 回足して累積 **0.0000e+00**。
>
> **撤回 2** (#13 の誤り): 「`dependentVariables` の無条件書き戻しが毎 step 下位ビットを落とす」も
> **定常陰解法では誤り**。順序は `assembleResidual` (`main.cpp:1641` → 内部 `:1391` で
> `dependentVariables` が `ro`/`roe` を書く) → `blockDPLURSolve` (`:1668`) →
> commit `ro = roN + dq` (`:1671`) で、**`roN` は書き戻しを受けていない前 step の値**。
> したがって書き戻しは commit に上書きされ**状態に残らない**。役割は残差評価に使う原始量
> (U/P/T/Ht/sonic) を床済み $\rho$ と整合させることだけ (`main.cpp:1797-1800` のコメントが明記)。
> → **`dependentVariables` は触らない**。触ると TP 二相の `twophaseFail` 分岐 (`:209,211`) や
> `condEquilibrium==2` の `rog` 射影 (`:154,255`) を壊す。

**採る設計**: FP32 の `Q` / `Q_N` は**演算用にそのまま残す**。**内点 CV (`nCells`) だけの FP64 配列
`Qacc` 5 本を正本**として加える。

| 段 | 処理 |
| --- | --- |
| commit | `Qacc += dq` (**FP64・in-place**) → `Q = (float)Qacc` |
| step 末尾 (`updateVariablesOuter` に**融合**) | reconcile: `Q != (float)Qacc` のセルだけ `Qacc = (double)Q` |

reconcile は「**FP32 の writer が書き換えたセルはその値を採用し、それ以外は残余を保持する**」という規則である。
壁ピン・周期ミラー・化学種の `roe += droe` など、commit 以外の writer と自然に両立する。

**なぜこれか**:

- **OFF 経路はカーネル 1 本も変わらない** → G1 (後方互換) の切り分けが効く。
  案 A (`Q` の型変更) は `c_d["ro"]` **112 箇所**と `std::map<std::string, flow_float*>` に波及してこれを失う。
- **残余ゼロなら ON は OFF とビット同一**。$f_{32}(a+b) = f_{32}(f_{64}(a)+f_{64}(b))$ が最近接丸めで成り立つ
  (20 万サンプルで **100.0000 %** 一致を実測)。→ **ON 1 step と OFF 1 step の比較が G1 の安価な補助ゲート**になる。
- 100 step の模擬: 現行構造 **0.0**、影アキュムレータ **2.9700e-08** (= $N\cdot dq$、ミラーは 16 ULP 動く)。

**`Q_N` は廃さない**。定常では step 冒頭 `roN == ro` なので FP64 側は in-place 1 組で等価だが、
float32 の `roN` は `updateGuardScale` (`update_d.cu:246`)・RK (`timeIntegration_d.cu:396`)・
dual-time (`update_d.cu:407`, `implicitCorrection_d.cu:53`)・軸対称 (`axisymmetricSource_d.cu:321`) が読む。
**FP32 `roN` は `(float)Qacc` のミラーとして残す**。dual-time 対応時は `QaccN/QaccNN` を shift する形で横展開する。

### 4.4 切替の既定と後方互換

- `time.qAccumulatorFP64` (0/1、**初期版は既定 0**)。

> **訂正 (codex plan M5)**: 初稿は「OFF が静かに収束しない状態だから既定 1」と書いたが、
> **対応経路を検証する前に既定を変えてはいけない**。`Q` 配列は陽解法・dual-time でも共有され、
> commit も別実装 (`update_d.cu:464`)。§2 で対象外とした SST も平均流の `roe` を書く。
> **既定 ON は、未対応経路に黙って適用されることを意味する。**

- **対応する GPU・時間積分・物理モデルの組合せを明示**し、**非対応で明示 ON なら起動時に拒否**する
  (黙って劣化させない)。要求値と実効値の両方をログと `RUN_PROVENANCE.txt` に残す。
- 既定を 1 に変えるかは、対応経路の検証が済んでから別途判断する。

**v1 の対応範囲** (`diagnostician` 2026-09-23):

| | 経路 |
| --- | --- |
| **対応** | `timeIntegration 11 && unsteady 0`、block/scalar DPLUR、CPG/TP 単相、node/cell、平面 |
| **起動時拒否** | unsteady/dual-time、陽解法、**node 軸対称** (`axisymmetricSource_d.cu:321-322` が `roN` を直接書く)、**`sstEnergyIncludesK=1`** (`ransTransport_d.cu:175` が毎 step 全 SST セルの `roe` を書くので累積が全域で消える) |
| **カウンタつきで許容** | 等温壁ピン (`main.cpp:1802`、壁ノードのみ)、化学種の `roe += droe` (`speciesTransport_d.cu:312`)、node 周期ミラー (`main.cpp:1753`) |

**採用セル数のカウンタ**を `monitorInterval` ごとにログへ出す (reconcile が何セルで発火したか)。
これが想定外に多ければ、その writer が累積を消している。

### 4.4a runtime 型切替と restart (**codex plan M4**)

初稿は「5 配列の確保を倍精度にする」としか書いていなかったが、**それでは runtime 切替にならない**:

- 状態は `std::map<std::string, flow_float*>`、host 側も `vector<flow_float>` (`variables.hpp:22`)。
  確保・転送もこの型を前提 (`variables.cpp:312,347`)。
- **出力は `vector<flow_float>`** (`output.cpp:156`)、**保存量の読込みは `vector<geom_float>`**
  (`variables.cpp:683`)。→ **device だけ倍精度にすると restart で蓄積した下位ビットを失う。**

**設計に含めること**:

1. 型付きの専用格納領域と、その初期化・転送・解放
2. **FP64 checkpoint** (下位ビットを落とさない出力)
3. 旧 FP32 入力からの初期化経路 (既存 run からの継続)
4. ゲート: **連続実行と中断再開が一致すること** (§6 の G4)

### 4.5 メモリ

影アキュムレータは **40 B/CV** (FP64 5 本、内点のみ)。65k CV で +2.6 MB、**1000 万 CV の 3D で +400 MB**。
(初稿の「+200 MB」は `Q` の型を変える案 A の数字で、本設計では倍になる。)
GPU メモリが逼迫する規模では `qAccumulatorFP64: 0` で従来に戻せる (§4.4)。

## 5. 実装ステップ

1. `methods/time_integration/implementation.md` の現在仕様を更新する (**済**: 「commit の丸め」節)。
2. 本計画を書き、**codex の plan 段レビュー**を受ける (**次**)。
3. `time.qAccumulatorFP64` を `solverConfig` に追加し、起動ログに出す。
4. `Q` 5 本の device 確保を倍精度にし、commit カーネル 2 本 (`applyBlockImplicitCorrection_d` /
   `applyScalarImplicitCorrection_d`) を倍精度で書く。読み出し側は float32 へ落とす。
5. §6 の検証を回す。
6. 化学種・乱流・凝縮モーメントへの横展開は**別計画**に切る (本計画では扱わない)。

### 5.1 残作業 (優先順)

| # | 担当 | 項目 | 内容 |
| --- | --- | --- | --- |
| 1 | F | ~~codex plan 段レビュー~~ **済 (NO-GO 全件採用)** | §6.1 |
| 2 | F | ~~変換点の妥当性~~ **済** | §4.3 で影アキュムレータに確定。`Q` の型は変えないので 112 箇所の変換は不要になった |
| 4 | F | ~~**FP64 正本の定義**~~ **済 (2026-09-23 `diagnostician`)** | §4.3。`Qacc` 5 本 (内点のみ) を正本、`Q`/`Q_N` は FP32 のまま。**案 A (型変更) も案 B (EOS 書き換え) も却下** |
| 13 | O | ~~書き込み側の棚卸し~~ **済・結論は訂正** | §4.3 の撤回 2。定常経路では `dependentVariables` の書き戻しは commit に上書きされ**状態に残らない**。`periodicNode` も名前依存の grep で漏れていた writer がある (`periodicBroadcastFromRoot_d`) |
| 11 | O | ~~`mdot_decay.py` に区間分離~~ **済** | 膝検出で減衰区間と床を分離 |
| **S0** | O | ~~**G0 の単体テスト**~~ **PASS (2026-09-23)** | `solver_density_cuda/cuda_forge/qAccumulator_d.cuh` (commit/reconcile の device 関数) と `tests/unit/test_qacc_commit.cu`。`nvcc --expt-relaxed-constexpr -I. -o test_qacc_commit tests/unit/test_qacc_commit.cu` で単体ビルドできる。<br>**(a)** `Qacc-Q0` = 2.9700e-08 = $N\,dq$ (相対差 **0.00e+00**)、ミラーは **16 ULP** 動いた / **(c)** OFF 経路は **0.0000e+00** (現行の症状を再現) / **(b)** 壁ピンは**ピンから 0.478 ULP しか離れない**・採用 25/100 step / **(d)** 残余ゼロの 1 step は 20 万サンプルで **ON/OFF がビット一致** (差 0 件)。**VERDICT: PASS**<br>⚠ 初回は 3 件 FAIL したが**いずれもテストの期待値の誤り**だった: (c) の −7.93e-11 は `(float32)RO0` と `RO0` の**表現差**を累積と取り違えたもの、(b) の「毎 step 採用」は誤りで、**残余が ½ ULP を超えたときだけ reconcile が発火する** (理論値 ~32 回) のが正しい挙動 = 「ピンから離れられない」という性質そのもの |
| **S1** | O | 配線 | `time.qAccumulatorFP64` キー、`Qacc` 確保 (`nCells` のみ)、block/scalar 両 commit に FP64 分岐、**reconcile を `updateVariablesOuter_d` に融合** (新規 launch を増やさない)、起動時拒否リスト、採用セル数カウンタ。→ **G1** |
| **S2** | O | FP64 checkpoint | 既存 `/CHECKPOINT` の契約 (`output/output.cpp:167-194` 書き、`main.cpp:983-1059` 読み・layout 不一致は拒否) に `qacc_ro…` 5 本を double で追加。**読込は `vector<double>`** (現行 `readValueHDF5` は `vector<geom_float>`=float なので使い回せない)。→ **G4** |
| **S3** | F | **#12 の判別 A/B** | 下記 |
| **S4** | O | G3 (`case/36`・`48`・`44`・`09`) と G5 | 量ごとの許容値を事前登録。scalar/block 両 commit |
| 9 | F | codex result 段レビュー | `done` にする前 |
| 10 | O | 横展開の別計画起票 | 化学種 `roY_s` / 乱流 `roK`,`roOmega` / 凝縮モーメント。dual-time は `QaccN/QaccNN` の shift |
| 14 | O | 調査 1 件 | 起動時拒否に `speciesImplicitCoupling==2` を入れるか (`speciesEOSFinalCommit` が `ro` も書くか)。v1 は CPG 単成分なので急がない |

### 5.2 判別 A/B (S3、**合格条件を事前登録**)

`run_0014` の最終場から、**影アキュムレータ ON・float ビルド・`FORGE_CUDA_BLOCKSIZE=128`・他設定完全同一**で
**100k step** (`run_0022_qacc_f32`)。見る量は `mdot_decay.py` の**減衰区間の e 折り**と 100k での低下倍率。

| | 条件 | 結論 |
| --- | --- | --- |
| **A** | e 折りが **3.23e4 の ×2 以内**、100k で **10 倍以上**低下 | **commit 丸めが律速**。`geom_float` と流束精度が**同時に除外**される |
| **B** | 100k で **2 倍未満** (`run_0021` の 1.4 % と同級) | FP32 の**流束・残差評価**に別の律速がある → plan を「残差累積の FP64 化」へ広げる |

**注**: G2 の「3 桁以上」は e 折り 3.23e4 なら $\ge$2.3e5 step 要る。**100k の A/B とは別物**なので、
A/B の合格条件は上記 (1 桁) で登録し、G2 は 300k 以上で判定する。

## 6. 検証

> **改訂 (2026-09-23, codex plan M3/M6)**: 初稿の G1–G5 は (a) 参照 run 自身が合格しない値を
> 合格条件にしていた、(b) 「ビット不変」の成立条件を誤っていた、(c) `case/08.bump` を回帰基準に
> 使おうとしていた (現行手順が「まだ基準に使えない」と明記)。全面的に書き直す。

### 6.0 まず通すゲート (これが通らない設計は先へ進めない)

| ゲート | 内容 | 合格条件 (**事前登録**) |
| --- | --- | --- |
| **G0** | **微小増分の保持** | $\lvert dq\rvert$ = 0.1 ULP 相当の一様な増分を人工的に与え、$N$ step 後の状態変化が $N\cdot dq$ に**比例して増える**こと (現行構造では 0 のまま。§4.3 の表) |

### 6.1 本体

| ゲート | 内容 | 合格条件 (**事前登録**) |
| --- | --- | --- |
| **G1** | 後方互換 | `qAccumulatorFP64: 0` で既存 run の `res_*.h5` が**数値配列として一致**。<br>**「ビット不変」とは書かない** — 成立するのは*同じ FP32 入力に対する同じ演算*までで、蓄積で状態が変われば流束もリミッタも変わる (M6) |
| **G2** | 主検証 (`case/56.gap_tp1187`) | **減衰区間と末尾区間を分けて**判定する。**300k step 以上で判定** (e 折り 3.23e4 なら 3 桁に $\ge$2.3e5 step 要る。100k の A/B とは別物)。<br>① **減衰区間**: $\lvert\dot m\rvert$ が幾何級数で 3 桁以上落ちる<br>② **末尾区間**: 床の値と変動幅を報告し、`check_quasisteady` で `STEADY` または `OSCILLATING` (平均±振幅で報告)<br>③ **全保存量の `check_convergence` VERDICT** と**対象量の `check_quasisteady` VERDICT** を併記<br>④ SU2 比較も**両者の判定つき**で行う |
| **G3** | 標準検証ケース | **`case/08.bump` は使わない** (`procedures/verification/README.md:23` が「まだ回帰の基準には使えない」と明記)。代わりに **node の SST = `case/36`・`case/48`**、**軸対称 TP・凝縮 = `case/44`**、**共有の周期経路 = `case/09`**。量ごとの許容値を**事前に**登録する。**scalar / block の両 commit** を試験対象にする |
| **G4** | restart 一貫性 | **連続実行と中断再開が一致**すること (FP64 checkpoint が無いと下位ビットを失う → §4.3 #2) |
| **G5** | 速度 | 1 step のコスト増が **+3 % 以内を実測**。<br>⚠ 初稿の「+0.1 %」は FLOP モデルで**誤り** (codex M7)。65k CV で 1.63 ms/step = **25 ns/CV は launch+sync 律速**で、別カーネルを足すと +1〜3 %。→ **commit と `updateVariablesOuter` に融合し新規 launch を増やさない**ことを設計要件にする |
| **G6** | メモリ | 増分が設計どおり (正本の置き方が決まってから数える) |

### 6.2 `run_0020_double` は「参照値」でなく「上界の参考」

**初稿の G2 参照値は撤回する。** 再判定した実測:

| 確認 | 結果 |
| --- | --- |
| `\|mdot\|` @ 1.9M (初稿が引用) | 4.5634e-12 |
| **`\|mdot\|` @ 2.0M (最終)** | **1.4798e-11** (最後の 100k で **3 倍に戻った**) |
| `mdot_decay.py` 全区間 | **べき乗則** (残差 RMS: べき乗 0.837 / 指数 1.043) |
| `check_convergence.py` 全 400k | **`NOT CONVERGED (stalled/plateau)`** |
| `check_quasisteady.py` (`gap_series.csv`) | `zW422`・`U_rms_deep` **STEADY** / `zW844`・`zW1411` **DRIFTING** |
| `mdot` 時系列の `classify_series` | **DRIFTING** (drift 137.5 %/tail) |

→ **「SU2 の床 8.59e-12 を下回った」は撤回**。1.9M の谷を床として報告したもので、
AGENTS が戒める「過渡を定常値として報告」に該当する。

**区間分離した正しい記述** (`mdot_decay.py` の膝検出、#11 で実装):

| 区間 | 内容 |
| --- | --- |
| 減衰区間 (先頭 16 点) | $\lvert\dot m\rvert\propto\exp(-3.10\times10^{-5}\,\mathrm{step})$、e 折り **3.23e4 step**、**4.59e4 倍**落ちた |
| 床 (末尾 2 点) | **1.234e-11 ± 2.5e-12** (振れ幅 1.50 倍)。**最小値 4.56e-12 を床と呼ばない** |
| float32 (`run_0014`) | **膝が検出されない** = 床に届いていない。べき乗則のまま、e 折り 4.98e6 step |


なお、この `DRIFTING` 判定は「初期の指数的減衰が存在しない」という意味ではなく、
**床付近まで一括 fit する現在のゲートが不適切**ということである (→ `mdot_decay.py` に
減衰区間と末尾区間の分離を入れる、§5.1 #11)。

### 6.1 レビュー記録 (codex)

| stage | 日付 | 記録 | 判定 | 指摘 | 対応 |
| --- | --- | --- | --- | --- | --- |
| 診断 | `2026-09-23` | (セッション内 `diagnostician`、結論は §4.3・§5.2) | 案 A/B とも却下 → **第 3 案 (影アキュムレータ)** | **全件採用・実測で確認**。**#13 の結論を訂正**: 定常経路では `dependentVariables` の書き戻しは commit に上書きされ状態に残らない (`main.cpp` の順序 1641→1391→1671 を自分で確認)。$f_{32}(a+b)=f_{32}(f_{64}(a)+f_{64}(b))$ を 20 万サンプルで検査し **100.0000 % 一致** → 残余ゼロなら ON は OFF とビット同一。**G5 の「+0.1 %」は FLOP モデルで誤り** (25 ns/CV は launch 律速) → 「+3 % 以内を実測」+ 融合要件に訂正。メモリも 40 B/CV で **+400 MB** に訂正。対応範囲の拒否リスト (`sstEnergyIncludesK=1` は全域で累積が消えるので拒否) を §4.4 に |
| plan | `2026-09-23` | [2026-09-23-time_integration-fp64-accumulator-plan.md](../../notes/reviews/2026-09-23-time_integration-fp64-accumulator-plan.md) | **NO-GO**, C0/M7/m1 | **全件採用・実測で確認**。**M1→§4.3 全面改訂**: commit は `Q = Q_N + dq` で `Q_N` は step 末尾に float32 コピー (`update_d.cu:39`)。`Q` だけ倍精度にしても **100 step の累積が 0** (自分で再現)。書き込み側 5 箇所の棚卸しを #13 に。**M2→§3・§4.1 の断定を撤回** (`run_0020` は `geom_float` も変えており「commit 丸めが唯一の律速」は未証明。3 本比較を #12 に)。**M3→§6 全面改訂・G2 参照値を撤回**: 1.9M の 4.56e-12 は谷で、**2.0M では 1.48e-11 に戻る** (自分で再現)。「SU2 の床を下回った」は過渡を定常値と報告したもので撤回。減衰区間と末尾区間の分離を #11 に。**M4→§4.4a 新設** (型切替・FP64 checkpoint・restart 一貫性ゲート G4)。**M5→既定を 1 から 0 へ**、非対応経路は起動時に拒否。**M6→G3 から `case/08.bump` を外し** `case/36`・`48`・`44`・`09` へ、「ビット不変」の主張も成立条件つきに訂正 |

## 7. 影響範囲

- `solver_density_cuda/input/solverConfig.{cpp,hpp}` — キー追加
- `solver_density_cuda/variables.{cpp,hpp}` — `Q` 5 本の型
- `solver_density_cuda/cuda_forge/update_d.cu` — commit カーネル 2 本と `updateGuardScale`
- `Q` を読む全カーネル — float32 へ落とす変換を挟む (§5.1 #2 で洗い出す)
- **既存 run の再現性**: 既定 ON なので、既定のまま回すと過去 run と**ビット一致しない**。
  再現には `qAccumulatorFP64: 0` を明記する。`RUN_PROVENANCE.txt` と合わせて追えるようにする。

## 8. 完了条件

- [ ] 関連 `methods/time_integration/` の現在仕様を更新済み (**済**)
- [ ] 実装・検証完了 (本計画の §6 を満たす)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [ ] ファイルを `plans/active/` → `plans/accepted/` へ移動

## 9. 変更ログ

- `2026-09-23` — 初稿。`case/56.gap_tp1187` の実測 (§3) から案③ を選定。
  案② (Kahan) は性能要件を満たすが `-ffast-math` で消えるリスクを理由に見送り (§4.2)。
