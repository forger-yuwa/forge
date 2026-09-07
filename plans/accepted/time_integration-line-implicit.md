# line-implicit (壁法線ライン block-Thomas の DPLUR 組み込み)

## メタ

- **area**: `time_integration`
- **status**: `done` (実装済・本ケースでは上限不変 = 負寄りの結果。opt-in 残置)
- **related_docs**:
  - [`methods/time_integration/implementation.md`](../../methods/time_integration/implementation.md) (「line-implicit」節)
- **related_plans**:
  - [`time_integration-update-positivity-guard.md`](../accepted/time_integration-update-positivity-guard.md) (負の結果 — 上限の真因が defect-correction 反復不安定と特定した前段)
- **created**: `2026-09-02`
- **owner**: ユーザ指示 (2026-09-02)

## 1. 目的

NS 陰解法の cfl_pseudo 上限 (素で 2–3、implicitRelax 0.7 併用で 8) は、高 AR 壁 BL の
壁法線方向結合を point-DPLUR が反復で解き切れないことによる defect-correction 反復不安定が
律速 (case/45 run_0015_cflsweep で特定)。**壁法線ライン上の結合を block 三重対角の直接解
(Thomas) に昇格**し、反復に残る lag をライン外 (流れ方向・弱結合) だけにして CFL 上限を
桁で引き上げる (FUN3D/DLR-Tau 系の標準構成)。

## 2. 設計 (実装済み構成)

- **ライン構築** (`mesh::buildImplicitLines`, host・初期化時 1 回): 壁 CV を種に
  「初手 = 最近接隣 (壁法線 = 最短エッジ)、以降 = 前進方向と cos>0.7 で最も揃う未訪問隣」
  の greedy 鎖。長さ <3 は不採用。構造化 TFI では j 列を再構成する。載らない CV は
  point fallback (`line_prev/next = -1`)。
- **sweep カーネル** (`implicit_defect_correction_block_d`): ライン CV は
  (i) ライン面の dq_old lag をスキップ (sdq=0、対角 A⁺ は従来どおり)、
  (ii) loop==0 でライン面の近傍行列 K を単位ベクトル 5 本の列抽出で保存
  (`accumulate_split_jacobian_cf` を dummy diag で呼ぶ — 点経路と厳密に同一の Jacobian)、
  decouple 行 (壁運動量/軸 roUy/等温壁 roe) は K 行を 0 化、
  (iii) 点解せず diag (loop==0)・rhs (毎 sweep) を保存。非ライン CV は完全に従来経路
  (lineImplicit=0 なら全 CV 従来経路 = 挙動不変)。
- **Thomas カーネル** (`lineThomas_d`, 1 ライン = 1 スレッド, 内部 double):
  $D_k\Delta Q_k - K^{prev}_k\Delta Q_{k-1} - K^{next}_k\Delta Q_{k+1} = rhs_k$ を
  前進消去 ($\tilde D_k = D_k - K^{prev}_k W_{k-1}$, $W_k=\tilde D_k^{-1}K^{next}_k$,
  $y_k=\tilde D_k^{-1}\tilde b_k$) + 後退代入で解き、implicitRelax を掛けて dq_new を上書き。
  LU 失敗ラインは前回反復値を保持 (保険)。
- **ドライバ**: 各 sweep の直後・swap の前に Thomas を起動。スイープ回数/周期ミラー等は不変。
- config: `time.deltaT.lineImplicit` (0/1, 既定 0)。blockDPLUR==1・timeIntegration==11 専用、
  lowMachPrecond>=2 とは併用不可 (起動時 fatal)。

v1 の割り切り: 1 スレッド/ライン (ワープ内要素並列は後続最適化)・因子分解は sweep 毎に
再計算 (loop 間キャッシュは後続)・非構造 tet のライン同定はスコープ外 (point fallback)。

## 3. 検証計画

1. 回帰: lineImplicit=0 で挙動不変 (同一バイナリ 2 回走ノイズ水準)。
2. ライン構築診断: M6 ノズルメッシュ (1250×97 構造化) で本数・被覆率・最大長。
3. 一致性: lineImplicit=1, cfl1–2 warm restart → point と同水準の残差床・NaN 0。
4. 上限: cfl {8,16,32,64,128} (CPG)・{8,16,32} (TP)、relax=1。
5. 収束レース: coarse IC 4000 step で最速設定 vs cfl8+relax0.7 (現最速)。

## 4. 結果 (2026-09-02) — **CFL 上限は上がらず (+収束 ~15%/step のみ)**

case/45 M6 NS (run_0015_cflsweep probes7–10):

- ライン構築は構造化 TFI で完璧 (1250 本 × 97 セル・被覆 100 %)。
- **cfl 上限**: line 化しても **積 cfl×relax ≈ 6–8 で point と同一** (line 素: cfl2 OK/8 発散、
  line+r0.7: 8 OK/16 発散、line+r0.5: 16 OK/64 発散)。
- **理由**: cfl8-line の発散は**出口域 (x>95) の全断面**で起きる — 律速モードは壁法線でなく
  **streamwise (ライン外) の defect-correction lag**。壁法線剛性は点法でも巨大な粘性対角が
  実質吸収しており、古典的な「高 AR = 壁法線律速」の期待はこのケースには当てはまらなかった。
- **収束速度**: 同一設定 (cfl8+r0.7, coarse IC 4000 step) で line は point 比 rms_roUx/roOmega
  とも **−14 %** (1 step あたり)。Thomas コストはほぼ中立 (82 s/4000step)。
- **今後の本丸**: streamwise 第 2 ライン族の交互解 (ADI 流) か LU-SGS の流下順序付け。
  v1 の 1 族機構はその土台になる。

## 5. 実装時に踏んだ罠 (再発防止)

- **lu5_solve のピボット適用バグ**: getrf 形 (行交換が L 部分にも及ぶ格納) に対し、交換と
  前進代入をインタリーブすると**誤解を返す** (数値照合で確定)。LAPACK getrs 流に
  「①全行交換 (LASWP) → ②L 前進代入 → ③U 後退代入」の順で書くこと。既存 `solve_5x5` は
  消去中に rhs を同時処理する形なので無事 (この罠は分解と解を分離した時だけ現れる)。
- **device printf の引数上限**: 引数 ~32 個超は**末尾引数がゴミ表示** (float 配列から
  1e227 等の "不可能な値" が出たら疑う)。診断 printf は 12 引数程度に分割する。
- **host の `var.c["ccx"]` は未充填 (全ゼロ) のことがある** — ライン構築の座標は
  `msh.nodes[ic].coords` を正とした。
- 検証ハーネス側: 「NaN 行が residual_history に残らない即死」があるため、完走判定は
  first_nan 無し **かつ** 最終 step 到達で行う (これを怠り初回スイープを誤判定した)。

## 変更ログ

- 2026-09-02: 起票・実装 (mesh ライン構築 / sweep カーネル K 抽出・保存 / lineThomas / config)。
- 2026-09-02: lu5 ピボットバグ修正・検証完了。上限不変 (+15 %/step) で done。診断 env
  (`FORGE_LINE_MAXLEN`/`FORGE_LINE_NOOP`/`FORGE_LINE_DEBUG_POINT`) は残置。
- 2026-09-02: **壁解像 DDES での A/B (case/39 periodic hills ny160, run_diag_lineimp_ctrl/on)** —
  本番 dual-time 構成 (dt 2e-6, nSub20, cfl_pseudo1, ir0.5) の発達乱流場 300 step 同一窓。
  持続サブ反復収束は roe +0.12 桁/20subiter・roUx +0.04 桁のみ (ctrl@20 相当到達 subiter 17-19)、
  step 単価は 2.44 倍 → **DDES では不採用**。ωバースト統計不変。リスタート過渡の回復のみ顕著 (roe 終端 1/18)。
  解釈: v1 は対流 K のみで粘性隣接結合が line 行列に無く (`viscous_diag` の対角集中のみ)、
  pseudo-dt の粘性スペクトル半径 (`setDT_d.cu` の 2ν/Δn 項) を割り引く根拠が作れない。
  DDES で効かせるには ①粘性非対角ブロック (まずスカラー αI 近似) の追加 + line 方向粘性 CFL の段階割引
  (lineViscousDtRelief) と ② K 抽出コストの削減 (サブ反復間の K 凍結等) の両方が前提。
