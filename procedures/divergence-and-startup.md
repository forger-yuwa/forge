# 発散の主因と安定起動の手順 (新規計算 / 発散時に必ず参照)

新しい計算を投入するとき、および計算が発散 (NaN / 残差爆発) したときは、まずこの文書の
チェックリストと起動手順に従うこと。経験上、発散の大半は**物理やメッシュではなく投入設定**
が原因であり、ここを直すと回ることが多い (Fluent 変換メッシュ fan / StaticMixer も、当初の
発散はすべて投入設定が原因だった。実証は [`procedures/calculation-workflow.md`](calculation-workflow.md) の Fluent 節)。

## 0. まず確認するチェックリスト

- [ ] **初期値は妥当か**。静止場 (U=0) に速度/圧力入口で流れを当てていないか。初期場の流速の
      向き・大きさを入口流れとおおむね一致させたか (例: 入口 `Ux=50` → 初期場も `U≈(50,0,0)`)。
- [ ] **出口 BC が流れに合っているか**。亜音速出口に `outflow` (超音速向け) を使っていないか →
      **`outlet_statPress`** を使う。`outlet_statPress` には**逆流用 `Pt`/`Tt` を必ず設定**したか
      (`Ps` だけだと逆流発生時に `Pt=0`→密度0→NaN)。
- [ ] **最初から難しい条件にしていないか** (2 次移流 / 乱流モデル / no-slip 壁を初手から)。
      まず易しい条件で収束させ、引き継ぎ計算で段階的に上げる (下記「段階起動」)。
- [ ] **乱流モデル使用時、ω が全域 0 になっていないか**。初期 ω/k に妥当な非ゼロ値を入れる
      (ω=0 は初手で発散する)。
- [ ] **非直交メッシュで `space.pRef` を入れたか** (動作静圧。free-stream の桁落ち対策)。
- [ ] **メッシュは妥当か**。多面体セルは現状不可 (tet/hex/prism/pyramid のみ)。
      非常に歪んだセルがあると 2 次移流が overshoot しやすい。

## 1. 発散の主因と対策

### (A) 初期値が良くない
- **静止場に速度/圧力入口を当てると危険** (起動衝撃で発散)。初期場の流速ベクトルを入口流れに
  合わせる。難しければ「易しい条件で一度回した結果」を初期場に使う (引き継ぎ計算)。
- 圧力入口 (`inlet_Pressure`) は初期場の総圧と指定 `Pt` がかけ離れていると過渡が暴れる。

### (B) 超音速向けの境界条件を亜音速に使う
- **`outflow` は超音速出口向け** (全量外挿)。亜音速出口には **`outlet_statPress`** (静圧指定)。
- `outlet_statPress` の**逆流分岐は `Ptb`/`Ttb` から逆流ガスの静的状態を作る**。`Ps` だけだと
  `Pt=0`→ρ=0→NaN。混合・剥離・再循環など出口逆流が起こりうる流れでは `Pt`/`Tt` 必須。
- 入口も同様: 亜音速で `inlet_uniformVelocity` (全量固定=超音速向け) を使うと音響反射しやすい。
- **逆も同じくらい危ない**: **超音速入口に `inlet_Pressure` (亜音速の全条件入口) を使ってはいけない**。
  この閉包は**内点の Mach から等エントロピーで Ps・ρ を決める** (`boundaryCond_d.cu:1058`) ので、
  超音速入口 (全特性が流入) では内点が入口を支配してしまい、**`ρ↓ → U↑ → Ps↓ → ρ↓` の正帰還**で
  **入口の 1 節点だけが十数 step で枯れて発散する**。1 次では減衰するので**2 次にしたときだけ落ちる**
  という紛らわしい出方をする。実例: `case/08.bump` (入口 M=1.65 なのに `inlet_Pressure`) は
  `convMethod: 0` なら 5.1〜5.3 桁収束するが 1/2 次では step ~200 で NaN。
  入口を `inlet_uniformVelocity` に替えると 5000 step 完走した
  (plan [convection-node-wall-reconstruction](../plans/active/convection-node-wall-reconstruction.md) §4.34)。
  **切り分け方**: NaN 節点が 1 個で、その座標が入口境界上なら真っ先にこれを疑う。
  **⚠ `inlet_Pressure` 側にも実装の不整合がある** (2026-09-20): 速度を新しい音速に整合させないので
  **指定した Tt/Pt を再現しない** (反例: γ1.4/cp1005/T_i 100 K/|U_n| 100 m/s、指定 Tt 293.15・Pt 100000 →
  返る境界状態は Tt 284.232・Pt 89751)。亜音速で使う場合もこの誤差を承知しておくこと
  (修正は plan [convection-node-wall-reconstruction](../plans/active/convection-node-wall-reconstruction.md) §5.1 W1l)。
  亜音速は `inlet_Pressure` / `inlet_Pressure_dir` を推奨 (forge の亜音速 `inlet_uniformVelocity`
  は非反射化済みだが、圧力入口の方が安定なことが多い)。

### (C) 最初から難しい条件で計算する
- **2 次移流 (`convMethod` 1/2)**・**乱流モデル (`turbulence.model` ≠ none)**・**no-slip 壁 (`wall`)** を
  初手から使うと、未発達な場で overshoot して発散しやすい。
- → **段階起動** (下記 §2)。

### (D) 時間積分の選択 (発散主因ではないが収束効率に効く)
- 定常は陽解法 (`timeIntegration: 3`)・陰解法 (`timeIntegration: 11`, `blockDPLUR: 1`) のどちらも
  使える。陰解法 (block-DPLUR) は大きな実効 CFL を取れて定常収束が速い/頑健なことが多いので、
  難しいケースの収束には陰解法が有利。陽解法でも回る。

### (E) 非直交メッシュの free-stream 桁落ち
- 非直交セルでは一様流の圧力流束 Σ(P·s) が float32 で厳密相殺せず偽の運動量源になる。
  `space.pRef` に動作静圧を入れて差分化する (詳細: [`.github/plans/convection-freestream-preserving-flux.md`](../plans/active/convection-freestream-preserving-flux.md))。

## 2. 推奨する段階起動 (易→難の引き継ぎ計算)

新規ケース、特に複雑形状/非構造メッシュでは、いきなり最終条件で回さず段階的に上げる:

1. **第1段 (起動・易)**: `convMethod: 0` (1 次風上) / 層流 (`model: "none"`) / 壁は `slip` か
   gentle / `outlet_statPress`(+逆流 Pt/Tt) / 入口に合う一様初期場 / 非直交なら `pRef`。
   時間積分は陰解法 (`timeIntegration: 11`, `cfl_pseudo`≈1, `nStepInner`≈10) が定常収束に有利
   (陽解法 `timeIntegration: 3` でも可)。これで残差を数桁落として安定な場を作る。
2. **第2段 (精度↑)**: 第1段の収束場を**引き継ぎ** (その `res_*.h5` を初期場にする)、
   `convMethod` を 2 次 (1=2ndUp / 2=MUSCL) に上げる。
3. **第3段 (物理↑)**: さらに引き継ぎ、`model: "sst"` (RANS) や no-slip 壁を入れる。
   **乱流の初期 k/ω は妥当な非ゼロ値**を全域に入れる (例 k=4, ω=500; ω=0 は初手発散)。
4. 各段で `detectNaN: 1` を付け、残差全列が下がることを確認してから次段へ。

引き継ぎ計算の作り方: 前段の収束 `res_N.h5` の `/VALUE/{ro,roUx,roUy,roUz,roe}` を、メッシュ
入りの h5 (`/MESH`,`/PLANES`,`/CELLS`,`/BCONDS` を持つもの) の `/VALUE` に上書きし、RANS なら
`/VALUE/roK`,`/VALUE/roOmega` を妥当値で作る。`mesh.meshFileName`/`valueFileName` をその h5 に向ける。

## 3. 発散したときの診断手順

1. `time.deltaT.detectNaN: 1` で実行 → 最初に非有限が出た step で `res_nan_<step>.h5` が出る。
2. `residual_history.csv` の**全 `rms_*` 列**を見て、どの保存量が先に発散したか確認する
   (`rms_ro` だけで判断しない)。発散が何 step 目から始まったか、初手 (step 0-2) か後期かを見る。
3. `res_nan_*.h5` で**非有限セルの位置**を特定する (入口/出口/壁/内部のどこか)。出口域なら逆流、
   入口域なら起動衝撃/入口 BC、壁際なら no-slip/乱流、を疑う。
4. 上の §0 チェックリスト・§1 主因に照らして設定を直し、§2 段階起動からやり直す。
5. メッシュ起因が疑われるなら閉性・体積・向きを確認する (Fluent 変換は
   [`procedures/calculation-workflow.md`](calculation-workflow.md) の Fluent 節参照)。

## 4. 関連文書

- [`procedures/calculation-workflow.md`](calculation-workflow.md) — 計算準備・実行・Fluent 取り込みと実行レシピ。
- [`procedures/solver-settings.md`](solver-settings.md) — `convMethod`/`limiter` 等の意味。
- [`.github/plans/convection-freestream-preserving-flux.md`](../plans/active/convection-freestream-preserving-flux.md) — `pRef` (free-stream 保存)。
