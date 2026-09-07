# node (median-dual) 2D の傾斜壁・出口コーナー不安定 (ノズル系で初手〜百 step 発散)

## メタ

- **area**: `boundary / discretization`
- **status**: `in-progress`  <!-- 課題1・2 解決済み (nodeAxisDirichlet 実装/検証済)。課題3 (細壁) と出口コーナー再評価が残 -->
- **related_docs**:
  - `methods/discretization.md` / `methods/boundary.md`
- **related_plans**:
  - [`tooling-nozzle-phase0-foundation.md`](tooling-nozzle-phase0-foundation.md) (発見元: node での再計算が全滅)
  - [`../accepted/turbulence-node-inlet-dirichlet-conserved.md`](../accepted/turbulence-node-inlet-dirichlet-conserved.md) (入口側の類似修正 — 出口版が欠けている疑い)
  - [`../accepted/diffusion-node-boundary-real-distance.md`](../accepted/diffusion-node-boundary-real-distance.md) / [`../accepted/turbulence-node-sst-wallfunction.md`](../accepted/turbulence-node-sst-wallfunction.md)
- **created**: `2026-08-03`
- **owner**: `sano` (調査: Claude 自走)

## 1. 目的 (現象)

forge_design の③ベルノズルメッシュ (構造化 quad 221×65、壁クラスタ AR≦19、品質 PASS、
**cell モードでは同一 msh が完走**) を node で回すと、あらゆる設定で序盤発散する。
2D node の検証実績 (平板 case/26・backstep case/18・case/36) は全て**軸平行壁**であり、
傾斜壁 + 超音速沿い流れ + 超音速流出コーナーを持つノズル系は node 初適用で未踏だった。

## 2. 調査結果 (2026-08-03, 全て 400 step 以下の再現 run at scratchpad)

### 切り分けマトリクス (単独では全て無効)

| 変更 | 結果 |
| --- | --- |
| SST → laminar / convMethod 1→0 / cfl 4→0.5 / axisym→平面 | いずれも step 2–4 で NaN |
| outlet_statPress→outflow / wall→slip / inlet_Pressure→速度入口 | 同上 (残差爆発値がビット同一 = 種は境界種別に非依存) |
| bndFirstOrder 0→1 / implicitRelax 1.0→0.3 / 壁ノード IC 速度ゼロ化 | 同上 |
| 双対幾何検査 (体積・面積・法線外向き) | 異常なし |

### 段階起動 (stage A: convMethod 0 + cfl_pseudo 3 + implicitRelax 0.5) で層が分離

| 構成 | 発散 step | 種の位置 |
| --- | --- | --- |
| 既定設定 (2次, relax 1, cfl 4) | 2–4 | 傾斜壁 (収縮部/ベル) の壁・壁隣接ノード。step 1 で壁 slip ノードの ΔroUx ~ O(10³) |
| stage A + SST (wf 0/1 とも) | 11–15 | **出口∩壁コーナーの outlet ノードで roOmega → 1e19–1e20** (流れ残差は健全 1.7e-4) |
| stage A + laminar | 128 | ベル中腹の壁隣接帯 (x≈47mm)。それまで残差は順調に低下 (1e-4) |
| 参考: 準直線壁 (勾配≦5°) 対照 | 5–6 (SST) | 壁隣接。傾斜が緩いと遅延 → 勾配依存性あり |

### リグレッション対照 (2026-08-03 追記)

「以前は node が動いていた」との整合確認: 検証済み node ケース case/26 平板 SST
(run_node_sst_final の入力一式、turbulence ブロックのみ新 API に翻訳) を**現行バイナリ**で
500 step 再実行 → **NaN なし・残差順調低下 (1.6e-3→2.5e-4)**。すなわち**ブランチの
リグレッションではなく**、軸平行壁で検証された node が、ノズル構成 (傾斜壁・超音速流出
コーナー) という未検証領域で破綻している。なお既存 node run の solverConfig は旧
turbulence API のため現行バイナリではそのまま再実行できない (要翻訳) 点にも留意。

### 結論 (二層の独立した node 問題)

1. **壁境界の運動量不安定 (傾斜壁×沿い超音速流)**: 設定を柔らかくすると遅延するが根治しない
   (2→15→128 step)。壁隣接ノードで成長するモードで、slip でも no-slip でも発生。既知の
   未修正バグ「node slip + 接線密度勾配の市松スプリアス流」(node-slip-spurious-flow) と
   同族の可能性。壁半割面流束の傾斜壁での取り扱い (弱形式置換・実距離 over-relax の適用域)
   を疑う。
2. **出口境界の SST スカラー爆発 (コーナー)**: roK/roOmega が出口∩壁コーナーから成長。
   入口側は [turbulence-node-inlet-dirichlet-conserved] で「Dirichlet 保存量整合+残差除外」
   が入っているが、**出口側に同等の処置が無い**。

## 2.5 追補 (2026-08-03 後半): 実証レシピ発見と原因の確定 — §2 の初期仮説を大幅更新

ユーザ指摘「以前ノズルでも node をやっていた」に基づき **case/29.bell_vs_conical の node
キャンペーン (run_0019〜0041: 層流壁レシピ・CFL 掃引・implicit 掃引・SST 実走)** を発見。
§2 の「ノズル系は node 初適用」は**誤り**だった。以後の対照実験で真因を確定:

- **リグレッション対照 2**: case/29 `run_0038_node_sst_wc` (node SST ノズル) を現行バイナリで
  verbatim 再実行 (turbulence ブロックのみ新 API 翻訳) → 3000 step 健全 (rms_ro 7.9e-5→1.45e-6)。
  **リグレッションではない**。
- **真因 = 壁第一セルの細かさ**: 同一レシピで `wall_first_frac` 1e-3 (第一セル 10–30μm,
  y+~15–50) は全設定で発散、**5e-3 (50–150μm) で完全収束** —
  `case/40/run_0006_bell_node_warm`: **VERDICT PASS / ALL PASS** (全列 3.0–4.6 桁低下,
  η_CF=0.9809 vs cell 0.9790)。case/29 の壁寄せ (Bump0.4) メッシュが粗めだったことと整合。
- **実証レシピ (node ノズル SST)**: explicit RK3 + cfl 0.1 + convMethod 0 (1 次) +
  `bndFirstOrder/nodeWallDirichlet/axisCentroidShift` + katoLaunder + **収束場からの
  warm start 必須** (冷間 IC は粗壁でも発散)。forge_design runner に既定として実装済み。
- なお §2 の「出口コーナー/傾斜壁」は症状の**発現場所**であり、独立バグとしての扱いは
  下記課題の解決後に再評価する。

### 残る解決すべきソルバ課題 (本 plan の実体・優先順)

1. **node × SST × 陰解法 (blockDPLUR)**: 収束済み node 場からの引き継ぎでも step 4 で発散
   (`run_0007_bell_node_2nd_imp` 相当)。case/29 でも implicit 成功は層流のみ (run_0035 cfl2–5
   で rms_ro 2e-9)。→ 最適化ループの評価速度に必須 (explicit cfl0.1 は ~10 倍遅い)。
2. **node × 2 次精度 (convMethod 1)**: 収束場からでも step 7–10 で発散 (cfl 0.3/0.15 とも)。
   case/29 の node SST も全て 1 次だった。→ 生産精度に必須。
3. **壁第一セル細分の頑健性**: y+~1 (low-Re) 相当の壁間隔が node で回らない。→ q_peak
   壁解像 (親計画 §4.6 ①) に将来必要。

## 2.6 追補 (2026-08-03 深夜セッション): 残課題 2 点の真因確定 — §2.5 の課題定義を再構成

case/40 で再現 run (run_0009〜0014) + 感度実験 (scratchpad E1〜E3) を実施し、両課題の実態を確定した。

### 課題 2 (2 次精度) — 解決: `bndFirstOrder: 1` の欠落が原因だった

- 発散を記録した run_0007 世代の設定は **`bndFirstOrder` なし** (当時の runner は無効キー
  `nodeWallViscGradFlux` を出力しており、`bndFirstOrder` は実証レシピ commit af0801b0 で導入)。
- 反実仮想で確定: 2 次+explicit+`bndFirstOrder` なし → step 8–10 で発散 (記録と一致)。
  **`bndFirstOrder: 1` あり → 12000 step 完走・全列 2.5–2.8 桁低下** (`run_0014_node_2nd_expl_long`)。
- すなわち「node × 2 次」は境界隣接 CV の再構成勾配だけが不安定源で、現 runner 既定
  (`bndFirstOrder: 1` = 境界 1 次化) で解決済み。

### 課題 1 (SST 陰解法) — 真因は node 軸対称の軸行崩壊 × SST 生産シート

- 「step 4 で発散」は再現せず (旧設定の 2 次との交絡)。実態は **遅発性**: 1 次+陰解法
  cfl_pseudo 2 で深く収束 (rms_ro 3e-9 @step 2000) した後、**roK が e-folding ~2000 step で
  指数成長し step ~7900 で発散** (`run_0013_node_imp_long`)。2 次+陰解法も同様 step ~10600
  (`run_0012_node_2nd_imp_full`)。cfl_pseudo 1 で成長率 ~1/4 に鈍るが残存 (E2)。
  `implicitRelaxSST` 0.5 は発散時期を変えない (E1) = スカラー更新の大きさは無関係。
- **モードの実体** (E2 スナップショット差分 + 場の検分):
  1. **軸ノード行 (y=0) がベル部 x≈0.031–0.048 で真空まで過膨張**: ro が隣接行の 1/10〜1/50、
     有効圧力が負 → EOS 床 P=1 Pa にピン。**これは「PASS」の explicit run (run_0006/0008/0014)
     でも慢性的に存在** (~20 ノード)。cell は同域で完全に健全。laminar 陰解法でも治癒せず悪化
     (E3) = **SST 以前の基底スキーム欠陥** (node 軸対称の軸半 CV が radial 圧力平衡から
     デカップルし準 1D 的に過膨張)。case/29 でも SST run で同傾向 (minP 2e3、laminar は健全 3.4e4)
     = 膨張の強い case/40 (Pt 4 MPa, ε9) で顕在化が激しい。
  2. 崩壊した軸行と第一内点行の間に**偽せん断** (dUxdy ~ -6e3) が立ち、SST が第一行に
     **k ~ 1.3e4 (周囲の ~1000 倍) の生産シート**を形成 (μt/μ ~320, Pk>Dk)。explicit では有界の
     まま定在するが、陰解法の大 pseudo-CFL でこの平衡が弱不安定化し指数成長 → 発散。
- **修正方針 (本セッションで実装)**: **node 軸対称の軸ノードを解かず、radial 隣接ノードからの
  対称 Dirichlet に置換する** (`nodeAxisDirichlet`)。軸対称条件 ∂q/∂r=0, u_r=0 の 1 次離散化。
  代表点選択は forge SST 壁関数の代表点 (SU2 `FindNormal_Neighbor` 移植) と同じ機構、実装様式は
  forge 既存の `nodeWallDirichlet`/入口スカラーピンと同パターン。軸 CV の体積は O(Δr²) で保存への
  影響は無視できる。実装は §3.1。
- **SU2 との関係 (2026-08-04 訂正 — 当初「SU2 と同型」と書いたのは不正確)**: SU2 も同じ
  頂点中心 median-dual で**軸ノードを通常 DOF として解く** (`BC_Sym_Plane` = 鏡映 flux +
  法線運動量残差の射影のみ。Dirichlet 置換はしない)。SU2 で軸が壊れないのは、SU2 の軸対称が
  **planar 幾何 + 1/y ソース項方式** ([flow_sources.cpp `CSourceAxisymmetric_Flow`]:
  `Coord_i[1]>EPS` で `yinv·Volume·(…)`、軸上はソース 0) であり、軸 CV の体積・面積が普通の
  planar 値を保つため。forge は **r 重み幾何 (B 流儀)** で軸半 CV の体積・面積が r̄ に比例して
  消えるため軸で離散平衡が悪条件化する — `nodeAxisDirichlet` はこの forge 固有の脆弱性への
  **SU2 より強い対症**である。根治 (軸 CV を解けるようにする) は r 重み幾何の軸極限の精査
  (将来課題) に持ち越す。

## 2.7 追補 (2026-08-04): `bndFirstOrder` をレシピから撤去 — 課題 2 の種も軸だった

ユーザ指摘「境界隣接 CV の一律 1 次化は精度上危険 (壁境界層第一 CV を 1 次にすると η_CF の
壁摩擦積分を損なう)。やめてほしい」を受け再調査した結果:

- **課題 2 の発散種の再特定**: §2.6 で「壁隣接第一内点」とされていた 2 次発散
  (`bndFirstOrder` なし, iso_2nd_nobfo) の `res_nan` を精査すると、実際の NaN 種は
  **軸帯 (x≈0.048–0.051, 軸ノード+第一内点行) の roOmega** であり壁ではなかった。
  つまり**課題 2 も課題 1 と同根 = 軸行真空化**であり、`bndFirstOrder` は軸行の再構成を
  1 次化することで病変をマスキングしていただけ。
- **切り分け**: `gradLSQ: 2` (境界 bvar 閉包 LSQ) でも軸修正なしでは step 9 発散 (勾配法は
  無関係)。`nodeAxisDirichlet: 1` があれば **`bndFirstOrder` なしの全域 2 次で安定**
  (explicit 300 step / implicit cfl4 12000 step とも NaN なし)。
- **全域 2 次の挙動**: L7 は入口〜収縮部の微小リミットサイクル (ro 相対 ~1e-5 の場の揺れ) で
  rms_ro ~1e-7 プラトー (cell 全域 2 次 run_0002 のプラトーと同格)、計量は quasisteady
  **ALL STEADY**。L6 は **PASS (converged)**、L9 は rms_ro 9.2e-9 (実質床)。
- **処置**: runner の node 既定から `bndFirstOrder` を撤去 (全域 2 次)。ソルバ側フラグは
  診断用に残置 (既定 0)。η_CF (全域 2 次): L6/L7/L9 = **0.9822 / 0.9896 / 0.9942** (単調、
  cell 比 ~+1.1%)。検証 run: `run_0021_node_2nd_imp_nobfo` (主検証) /
  `run_0022_node_runner_nobfo` (E2E 基準) / `run_0023_bell_L6_nobfo` / `run_0024_bell_L9_nobfo`。
- 残課題 (追加): 入口/収縮部の 2 次リミットサイクルの解消 (limiter 挙動の精査) — 優先度低
  (計量 STEADY で実用上支障なし)。

## 2.8 追補 (2026-08-04): 壁 T の流れ方向エントロピー市松 — **解消済み (2026-08-12 実測確認)**

> **状態 (2026-08-12)**: 本節の見出しは長らく「未解決」のままだったが、**現行の生産構成では市松は
> 実質消滅している**。case/40 の `res_wall_3_12000.h5` からベル部 (x>10 mm, 93 点) の交番二階差分
> $|\Delta\Delta T|$ を実測した結果:
>
> | run | 世代 | 壁 T 平均 | \|ΔΔT\| mean / max | \|ΔΔP\|/P |
> | --- | --- | --- | --- | --- |
> | `run_0029_node_adiabfix` | 断熱リーク修正直後 (2026-08-05) | 1195.9 K | **43.5 / 173.6 K** | 0.639% |
> | `run_0055_node_yp30_dofonly_tawoff` | 現行 mode 0 | 1122.9 K | **0.32 / 9.9 K** | 0.066% |
> | `run_0066_node_yp30_taw_defectflux_prod` | **現行生産 mode 3** | 1414.0 K | **0.25 / 10.4 K** | 0.066% |
>
> 43.5 K → 0.25 K (~170 分の 1)、圧力側も 1 桁改善。以下「残る容疑」「次ステップ (項別計装)」は
> **不要になった**ため実施しない。
>
> **帰属は未確定 (交絡あり)**: run_0029 と run_0055/0066 はノード数同一 (14365) だが壁第一間隔が
> 3.8e-5 → 7.5e-5 m (bell y+ mean 44 → 101) と異なる。効いたのが §2.9 の双対 CV 幾何 float32 桁落ち
> 修正 (`3cdddb9e`) なのか、壁第一間隔が太くなってスリバー双対 CV が解消したのかは同一メッシュ A/B
> を取らない限り分けられない。**運用上は解決** (壁温は mode 3 で SU2 実状態と 10 K 一致し生産使用可)
> だが、y+≈1 級の細壁メッシュに戻す際は再測定すること。

ユーザ報告「壁面温度が流れ方向に振動」の調査結果:

- **現象**: 壁ノード列に沿う node-to-node の **等圧エントロピー市松** (P は滑らか 0.6%、
  T/ρ が逆相で ~12%/9%、ベル部 |ΔΔT| mean ~206 K)。**完全に静的** (8000→12000 step で
  0.02 K しか動かない)。
- **切り分け**:
  - warm start (cell→node 最近傍 interp) が壁列に ~520 K の市松を刻印する (res_0 で確認)。
    しかし **IC を等圧接線平滑化しても (res_0 で 5 K)、12000 step 後に run_0022 と
    ビット同一の 206 K 市松が再生成** → 種でなく**定常解そのものが市松を持つ**
    (決定論的アトラクタ)。
  - laminar でも発生し**むしろ増悪** (mean 426 K / max 6381 K) → SST 無関係の基底スキーム欠陥。
  - cell は同条件で ~32 K (1 桁小)。
- **機構仮説**: 壁ノードは u=0 で接線移流がなく、等圧エントロピーモードは接線熱伝導のみが
  平滑化手段。定常解に市松が残る = 壁半 CV のエネルギー収支に市松を汲む/伝導が見えない
  離散不整合がある (W–W′ 接線双対面の伝導・W–I 面の交互フラックス経路の監査が必要)。
  既知の node 壁市松ファミリ ([node-slip-spurious-flow]・近壁 dP/dy 振動) の一員。
- **対処状況**: warm start 平滑化 (`runner._smooth_wall_entropy`) は IC 清浄化として実装済み
  (定常市松は不変)。η_CF への影響は小さい (市松は等圧で、壁摩擦・推力積分は P/τ 経由) が、
  壁温・熱流束を使う評価では現状の node 壁 T は信頼できない。
- **根本原因の大半を特定・修正 (2026-08-05, Codex レビュー指摘の検証)**: node の断熱壁は
  **厳密断熱になっていなかった**。`viscousFlux_wall_d` の node 分岐が壁面伝導熱流束を
  「壁ノードのセル勾配 ∇T·S」で評価しており、法線成分がゼロでないため断熱壁から毎ステップ
  熱が漏れていた (cell は ghost 差分 Ts_g=Ts_c で厳密 0)。**断熱壁 (kind: wall) の伝導熱流束を
  明示的に厳密 0 化** (cell はビット不変、等温壁/WMLES 経路は不変)。
  反実仮想 (`run_0029_node_adiabfix`, この一変更のみ): **壁 T 市松 mean 206→43.7 K (−79%)、
  max 931→174 K**。壁温平均 1631.5→1195.5 K — 旧値は断熱境界の熱漏れで汚染されており、
  修正後は SU2 の 1428 K に対して低め側。η_CF 0.9896→0.9884 (−0.12%)。**残る ~44 K 市松は
  cell 並み (32 K) の水準**で、残余監査 (W–I 対流/τ·u/W–W′) は修正後ベースラインで続行する。
- **SU2 反実仮想で残る壁温差を特定 (2026-08-05, run_0030/0031)**:
  - SU2 low-Re 最終場へ `STANDARD_WALL_FUNCTION` だけを追加すると、第一内点は forge に接近
    (`U=925 vs 843 m/s`, `k=9.93e3 vs 8.54e3`, `μ_t/μ=8.3 vs 5.6`)。それでもベル壁温の
    面積重み平均は SU2 1422 K、forge 1193 K のまま。
  - コード差は明確で、SU2 `CNSSolver::SetTau_Wall_WF` は断熱壁を
    Crocco–Busemann の回復温度へ直接更新する。forge `wallTreatmentSST=1` は壁応力・SST だけを
    壁関数化し、**熱側は壁解像のまま**。forge 第一内点へ同じ
    $T_{aw}=T_1+\mathrm{Pr}^{1/3}U_{t,1}^2/(2c_p)$ (`Pr=0.72`) を適用すると面積重み 1418 K、
    SU2 壁関数 1422 K と 4 K 差 (点 RMS 10 K)。実値1193 Kとの差 ~225 K はほぼ全量が
    **壁面温度に対する圧縮性熱壁関数の欠落**で説明される。
  - 同時実施の forge low-Re D2 は SU2 low-Re と壁温1 K差・η 0.25%差へ揃うが、双方とも
    `y+≈19–40` でlow-Reとして未解像なので、`wf_pk` 欠陥や物理真値の確証にはならない。
    実際、SU2 standard wall function は x=30 mm 第一内点 `μ_t/μ=6.88`, T=1015 K と
    forge wf=1 cell (`3.80`, 997 K) に近いまま、Crocco–Busemann 閉包で壁温1427 Kを得る。
    よって「BL全体のμ_t不足だけが約240 K差を作る」という帰属は棄却する。
  - `wf_pk=ρu_τ^4/ν·g(1-g)` は Reichardt 則の
    `ν_t/ν=1/g−1`, `∂u/∂y=(u_τ²/ν)g` を代入した平衡生産
    `ρν_t(∂u/∂y)^2` と代数的に一致する。y+≈29で low-Re 生産より小さいこと自体は、
    未解像low-Re解を正本にしない限り欠陥の証拠ではない。
  - よって `y+≈23–138` の現メッシュで forge の壁温を定量使用してはならない。修正候補は
    (A) automatic SST に断熱回復/熱壁関数を追加、(B) `y+<1` の壁解像格子で三者比較。
  - 別交絡: forge CPG `viscMethod:1` は Sutherland μ に対し λ=0.0257 固定なので、ベル第一内点の
    分子 Pr は1.50–1.87。SU2 を同じ固定 λ にした run_0031 は壁温を +15 K 上げただけで主差ではないが、
    高温 air の物性としては不適切。constant-Pr 伝導率オプションが必要。
  - 旧 cell 基準 (run_0002/0017) もベル壁温面積重み ~1767 K、`y+` 平均~5であり、熱閉包/物性問題は
    node 専用ではない。旧値「Tt 超え=直ちに熱漏れ」とする説明は、実効 Pr>1 の指定物性も交絡するため撤回し、
    熱漏れの確証はコード上の非零境界熱流束と run_0029 反実仮想に置く。
- **準定常ツールの注意 (2026-08-05)**: `check_quasisteady.py` は `res_*.h5` に
  `res_outlet_2_*.h5` を混入し、拡張子 `.h5` の `5` まで step に連結するため、表示上の
  `7 snapshots, steps 5..2120005 / ALL STEADY` は不正。さらに既定量は η・ṁ・壁温を検査しない。
  run_0029 の η/ṁ/壁温は正しい本体3時点 (4000/8000/12000) を直接再抽出し定常を確認した。
- **エネルギー収支監査 (2026-08-04, 実験ベース)**: run_0025 収束場 (市松定在) からの摂動
  実験で項別に切り分けた:
  - **W–W′ 接線伝導は生きている**: κ×10 → 市松振幅 206→87 K (×0.42)。壁列を見ていない
    わけではない。ただし 1/κ スケーリングでない = ポンプ項も κ とともに変わる/基底場の
    変化と交絡。
  - **粘性フラックスの勾配補正項 (非直交補正) は無罪**: `VISC_DIAG_NOGRAD` 診断ビルド
    (面平均 GG 勾配項を全部 0 化、over-relaxed 直接項のみ) で市松 204 K = 不変。
  - **μ(T) Sutherland フィードバックは一部寄与**: viscMethod 0 (定数 μ) で 206→134 K (×0.65)。
  - 1 次精度でも発生 (run_0006 世代 191 K) = MUSCL 再構成は不要条件。laminar で増悪 = SST 無関係。
  - **残る容疑**: ① W–I 面の対流フラックス (SLAU の面平均/風上が交互の ρ_W・H_W を運ぶ)、
    ② 壁境界半割面 (弱形式) のエネルギー閉包、③ W–I 面の τ·u 仕事 (Uxf=面平均速度)。
  - **次ステップ (未実施)**: 項別計装 — res_roe への寄与を「伝導直接項/対流/τ·u 仕事/
    境界半割面」別に壁ノードで分解出力する診断配列を足し、市松モードへの射影 (±交互和) で
    ポンプ項を確定する。
- **W–W′ (壁ノード間) フラックスの現状 (コード確認)**: W–W′ 双対面は対流 (SLAU)・粘性とも
  **主ループの通常内部面**として処理される (特別扱いなし。`Tau_Wall`/`Qw_Wall` の再スケールは
  片端のみ壁の xor 条件で W–W′ を除外、壁応力エッジカーネルも壁-壁エッジを明示 skip)。
  u=0 Dirichlet 下では W–W′ の対流は 1 次で質量/エネルギー≈0 (圧力のみ)、2 次では面重心が
  壁からオフセット (双対面はエッジ中点→primal セル重心) するため再構成速度が非零になり
  実質的な交換が生じる。粘性は完全応力+伝導 (over-relaxed 直接項が市松を直接見る)。

## 2.9 追補 (2026-08-11): 課題 3 (細壁 y+≈1) 解決 — 真因は双対 CV 幾何の float32 桁落ち

case/40 の y+≈1 low-Re 三者基準解キャンペーン (`run_0032`–`run_0035`) で、node ×
`wallTreatmentSST: 0` × 壁第一間隔 ~1 μm (`wall_first_frac 1e-4`) が warm start・段階起動に
かかわらず step 0–4 で roK/roOmega 爆発する事象を切り分け、**独立バグ 2 件**を特定・修正した。

1. **双対 CV 重心/体積の float32 shoelace 桁落ち (`gmshReader.hpp buildMedianDual` 2D)** —
   高 AR スリバー CV (~176 μm × ~0.7 μm) では絶対座標の交差積 ~1e-4 に対し面積 ~1e-10 の
   6 桁相殺となり、float32 では体積が %級・**重心が ~100 μm 級に誤り CV 外に出る** (壁 CV の
   重心が壁の外 134 μm 等)。この `centCoords` 汚染が
   (a) node の `wall_dist` (壁点=壁 CV 重心, 照会点=CV 重心の kdtree) を桁で破壊
   (真値 1.4 μm の第一内点に 172–1227 μm を付与) → SST F1/F2・壁 ω ピン (`wall_y_eff`) が
   毒され BL 中層〜出口帯で ω→1e25 発散、
   (b) `interp_field.py` の dst 照会点を汚染し cross-mesh IC を市松化していた。
   **修正**: サブ四角形 shoelace を「ノード A 相対座標 + double」で計算 (厳密演算では従来と
   同値の精度修正)。5e-3 世代への影響は重心 ≤75 μm (median rel 3e-6)・体積 ≤0.33%・
   wall_dist 最大 33% (近壁, 絶対 ≤88 μm) の改善方向の修正で、既存 run の結論は不変。
2. **`interp_field.py` の dst 照会点** — node 入力 h5 で `CELLS/centCoords` を優先していた
   ため 1. の汚染を受け、さらに (汚染がなくても) cell→node の warm start では x 半間隔ずれ×
   壁半径成長により「上流列の壁第一セル値と下流列の BL 内部値を行ごとに交互に拾う」ω 市松
   IC (×300) を作っていた。**修正**: 値数=節点数 (node) なら `MESH/COORD` を照会点にする。

y+≈1 の node 段階起動レシピ (確定): warm start は **node 場から** (cell 場からの interp は
半間隔ずれの staircase が残るため非推奨) → stage A = **陰解法 cfl_pseudo 0.5・convMethod 0・
nStepInner 10・3000 step** (explicit は壁 ω~1e10 の剛性で不可) → 本計算 = 全域 2 次+陰解法
cfl_pseudo 4。結果: `run_0034_node_yp1_lowre` 全列 2.8–3.5 桁低下 (rms_ro 3.4e-9)・
ALL STEADY、cell (`run_0033`) と壁温 2.3 K 差・τ_w ~1% 差で整合。5e-3 世代で見えた
node wf=0 の ṁ −4.5% 異常 (§2.8 系) は y+1+幾何修正後は 0.26% 差に消滅。

## 2.10 追補 (2026-08-11): node 出口列の系統的不整合 — **node η_CF は出口積分アーティファクトで +1.3% 過大** (未解決・最優先)

ユーザ指摘「node の出口 Neumann がうまくいっていなさそう (わずかに不連続)」の定量化
(y+1 メッシュ `run_0037_node_yp1_prt09`、cell `run_0036` / SU2 `run_0035` を対照):

- **現象**: 出口列 (x=70 mm の境界ノード行) の場が内部 2 列からの線形外挿に対し
  **P で平均 19% (最大 24%)・T で平均 15%** 低下する不連続を持つ。**超音速コア (M≈4.1)
  でも P が 1 列で −12%** 下がる (32.9→28.9 kPa) — 超音速流出で背圧情報は上流に届かない
  ため純粋に数値欠陥。cell は同条件で平均 0.9% (最大 2.3%)、**SU2 (同じ頂点中心
  median-dual) は −0.2%** と清浄 → forge node 固有。
- **η_CF への直撃**: `thrust_metrics` は outlet の `iCells` = この汚染列で積分するため、
  **node η=0.9905 は出口列アーティファクト**。内部列 (x=67.8–69.3 mm) で trapz 再積分すると
  node η ≈ **0.970–0.977** = cell 内部列 0.973–0.976・SU2 出口 0.9796 と整合。
  **従来「node 2 次化の離散差」と帰属していた node−cell の η +1.1〜1.3% オフセット
  (README の dv 感度表・run_0016 以降の全 node η) はこの出口欠陥のアーティファクトであり、
  帰属を訂正する**。真の η(y+1 low-Re) ≈ 0.975–0.980。
- **切り分け済み**:
  - 全域 1 次 (convMethod 0, 診断専用 scratchpad D6) でも近壁行の sag −20% は残存、
    コア行は −12%→−4% に減少 → 2 次再構成 (出口ノードの片側勾配×limiter) は増幅要因だが
    種ではない。半割面の対流流束は `convectiveFlux_boundary_d` で F(bvar)、outlet_statPress
    の bvar は超音速分岐で node 自身の状態 = F(W) 純風上のため対流としては自己整合。
  - 残る容疑: **出口半割面の粘性流束の扱い** (近壁で τ_xy·S_x が大きい; sag が壁側ほど
    大きいことと整合) と、W–I 主ループ面の 2 次再構成の境界側閉包。次セッションで
    半割 CV の離散収支を項別に監査する。
- **当面の運用**: node run の η_CF は出口積分値をそのまま使わず、**内部列積分または
  cell/SU2 の値を正とする**。壁温・τ_w・ṁ は出口列に依存しないため node 値を使ってよい
  (ṁ も内部列で 1.2903–1.2988 と出口 1.2957 の差 ~0.3% に留まる)。
  → **§2.11 で根治済み (2026-08-11)。本注記の運用は不要になった。**

## 2.11 追補 (2026-08-11): 出口列欠陥の根治 — 真因は `outlet_statPress` の bvar `Psb` 未更新 (超音速分岐でも背圧を流束・勾配閉包が読む)

§2.10 の残容疑 (出口半割面の粘性流束 / W–I 面の 2 次境界閉包) はどちらも**種ではなかった**。
収束場 `run_0037/res_12000.h5` から出口ノード半割 CV の離散収支を SLAU 2 次込みで項別再構成
(質量収支 1e-9 で厳密再現) して真因を直接特定した:

- **真因**: `outlet_statPress_d` は ghost には超音速分岐で `P[ig]=P[ic]` を正しく書くが、
  **bvar `Psb[ib]` への書き戻しがコメントアウトされたまま** (`Psb` を規定背圧の持ち場として
  流用していたため書けなかった)。他の全 BC カーネルは `Psb` を毎ステップ動的に書いており、
  outlet だけが例外だった。その結果 node 弱形式で `Psb` を読む 2 経路が超音速流出でも
  **規定背圧 20 kPa** を使っていた:
  1. `convectiveFlux_boundary_d` の `p_tilde = Psb` — 出口半割面の圧力流束が内部 33–44 kPa
     に対し 20 kPa となり、CV が (P_int−P_back)·S_x の偽の運動量欠損を受ける (流束スケール比
     0.8–2.5%)。収束場はこの欠損とちょうど釣り合う位置 = P が沈んだ状態で静定していた
     (バグ式では res_x が 0.00% で閉じ、正しい式では欠損項がそのまま残差に現れることを確認)。
  2. GG/LSQ 勾配の境界閉包 (`calcGradient_b_d` 等) の `Pf = Psb`・`Tf = Tsb`
     (`Tsb` も `Psb` から計算され ro=内部×P=背圧の不整合温度) — 出口ノードの dP/dT 勾配を
     汚染し、2 次再構成経由で sag を増幅 (1 次で残存・2 次で増幅の観測と整合)。
  - cell モードは主ループが ghost (`P[ig]` 正) を使うため無傷、SU2 も ghost 方式で清浄 —
    三者の切り分け結果と完全に整合する。粘性の境界半割面欠落は cell/node 同構造であり
    判別要因ではなかった (コア行 M≈4.1 の −12% は粘性では説明不能)。
- **修正 (最終形, 2026-08-11 ユーザレビューで改訂)**: 当初実装 (commit e9b5e4f5) は
  「規定背圧をスカラー引数化し `Psb` を動的書き込みに正常化」だったが、**`Psb` は
  `inletProfile` CSV 等で面ごとの分布を与え得る入力であり、スカラー化・上書きのどちらも
  分布指定を壊す**とのユーザ指摘で以下に改訂 (後続 commit):
  - `Psb` は**入力専用** (outlet カーネルは読むだけ。分布指定が完全に保持される)。
  - node × `outlet_statPress` の消費側を **interior (ic) 参照**に切替: 弱形式流束の
    `p_tilde = Ps[ic]`、勾配境界閉包 (GG/LSQ) の P/T ジャンプ = 0 (zero-normal-gradient 閉包)。
    per-bcond フラグで outlet のみ。ghost (ig) は使わない (node の ghost 撤廃方針
    [discretization-node-boundary-ghostless] と整合)。入口/壁/slip は従来どおり bvar
    (入口は bvar 圧が BC の効き筋、壁/slip は Psb=P[ic] で同値)。
  - **亜音速の背圧アンカーは p_tilde でなく特性構成 bvar (ronew/Vn_exit) の質量流束が担う**
    ことを実測で確認 (下記 run_0029)。`Tsb` は出口では入力でないため、適用圧 Pnew と
    整合した派生出力として書く (勾配閉包はもう読まない)。
  - cell モードは境界カーネル skip + ghost 主ループのため全変更の影響外 (ビット不変)。
- **検証** (いずれも 12000 step 陰解法, 出力先は case/40 README の run 一覧参照):
  - `run_0043_node_yp1_outletfix` (run_0037 と同一入力の A/B): 出口列 sag **mean 19.4% → 0.39%**
    (max 24.3% → 1.6%, 超音速行 0.40%) で cell 0.9%・SU2 −0.2% と同水準。
    **η_CF 出口積分 0.9754 = 内部列 0.9747–0.9755** (基準達成)、ṁ 列間ばらつき消滅 (1.2958–1.2964)。
    quasisteady ALL STEADY。η は生産値帯 0.978±0.003 の下端内 (SU2 y+1 出口 0.9796, cell 内部列
    0.973–0.976 と整合)。
  - `run_0045_node_yp1_outletfix_cont` (+12000 step 継続): 場は不変 (出口 P 最大 6e-4) = 準定常。
    残差床は run_0037 比 ~10 倍高い低レベルリミットサイクル (rms_roUx 3.5e-5, rms_ro 1.1e-8) だが
    絶対レベルは十分低く派生量 ALL STEADY。
  - `run_0044_node_yp30_outletfix` (run_0040 と同一入力の A/B): sag 13.6% → 0.17%、
    **η 出口積分 0.9835 (アーティファクト) → 0.9687 = 内部列 0.9679–0.9687**、SU2 壁関数
    (run_0042: 0.9673) と +0.14%。全列 2.8–5.1 桁低下。
  - 非退行: 壁温は出口隣接 3 列を除き y+1 で ≤0.85 K・y+30 で ≤4.6 K 差 (出口隣接列の変化は
    汚染の解消方向)。ṁ 不変 (1.2864/1.2959)。
- **検証 (最終形 ic 参照実装, run_0050/0051/run_0029)**:
  - `run_0050_node_yp1_outletic` (run_0037 と同一入力): sag mean 0.40% / max 1.59%、
    η 出口積分 0.9754 = 内部列 0.9747–0.9756、ṁ 1.2959 — **前実装 run_0043 と同一の根治**。
  - `run_0051_node_yp30_outletic` (run_0040 と同一入力): sag 0.18%、η 0.9687 = 内部列、
    SU2 wf +0.14% — run_0044 と同一。
  - `case/26 run_0029_node_outletic_subsonic` (M0.2 平板 node, 収束場から 12000 step):
    **出口 P 平均 97250.6 Pa (規定 97250, ±3 Pa) を保持** = 亜音速の背圧アンカーが
    p_tilde=P[ic] でも質量流束経由で機能。残差は warm start から悪化なし・NaN なし。
  - 分布保持: `Psb` をどのカーネルも書かなくなったため構造的に保証。
- **運用変更**: node の η_CF は今後**出口積分値をそのまま使ってよい**。§2.10 の「内部列を正と
  する」運用と README 注記は解除。y+1 low-Re の node η_CF = **0.975**、y+30 wf の node
  η_CF = **0.969** が修正後の正値。

## 2.12 追補 (2026-09-05 → 本ブランチ 2026-09-07 移植): TP 亜音速出口の特性構成が γ を混用し偽逆流/出口列加圧を課していた — 修正済

`outlet_statPress_d` の亜音速分岐 (SU2 BC_Outlet 型特性構成) が `c_int = sonic[ic]` (TP の真の音速, γ_mix) と
`c_exit = √(ga_cfg·P_exit/ρ_new)` (config 定数 γ) を混用していた。P_int = P_exit でも c_int ≠ c_exit となり、
Riemann 不変量 $V_{n,exit} = U_n + 2(c_{int}-c_{exit})/(\gamma-1)$ が数十 m/s ずれる。超音速出口 (全量外挿) では
通らず CPG では γ が同一なので無害 — **TP × 亜音速出口だけが壊れていた**。

- **修正** (feature/chemistry-finite-rate `876ec8b4` の `boundaryCond_d.cu` 部分のみを本ブランチへ移植): TP では
  特性構成の γ をセル局所 γ_mix(T) (`var.c_d["gamma"]`) に統一し、`c_int` も同じ γ で √(γP/ρ) と再構成。
  超音速判定は従来どおり `sonic[ic]`。CPG・超音速はビット不変。
- **本ブランチでの指紋 (case/16 run_0198, node × TP split × SST, 2026-09-07)**: 出口壁ノード (u=0 ピン) は常に亜音速
  なのでこの分岐を通り、出口列の P が 2 kPa 指定を無視して 45〜57 kPa へ上昇 → ノズルが unstart (全域 M≈0.3,
  P≈Pt)。cell (run_0191/0196) は出口が全面超音速で無害だった。`nodeInletCornerWall` は入口角 P>Pt 暴走を消したが
  (run_0195 → run_0198)、この出口欠陥が残っていた。修正後の再計算は run_0199。

## 2.13 追補 (2026-09-07): node の出口壁列は常に亜音速 — `outlet_statPress` の Ps が実出口圧より桁違いに低いと SST で unstart (case/16)

case/16 (Wyslouzil, 平面, Pt 59 kPa, 出口 M1.8 / P≈10.4 kPa, Ps 指定 2 kPa) の node SST が出口側から unstart した (run_0198–0201)。
100 step 再現で、最終列に壁向き Uy (上流列の 5–10 倍) と +5 % の P バンプが**層流でも**あり (bounded)、SST では最終列の k が
上流の 5–8 倍に膨れて BL が閉塞→剥離前線が x≈94 mm から上流へ伝播→全域 unstart、と特定。TP/CPG・dilatation・壁関数・
`nodeWallDirichlet`・出口 k/ω・壁第一層 (0.6/2.4/8 µm) は無関係 (全て同一の失敗)。

- **機構**: cell では出口面が全面超音速で背圧ゴーストは通らないが、node では壁ノード (u=0) と壁列が常に亜音速分岐を通る。
  Ps=2 kPa に対し実出口圧 10.4 kPa なので Riemann 構成のゴーストが Vn_exit≈+365 m/s の強い膨張を毎 step 課し、最終列の
  半 CV に壁向き速度/P バンプを作る。case/45 (Ps 2237 ≈ 実出口 2245 Pa) では同じ最終列バンプが +1.2 % に収まり無害。
- **対策 (2 通り、どちらも有効を確認)**: (a) 出口を `outflow` (全量外挿) にする (run_0210/0211/0212/0213) — 超音速出口の標準扱い。
  (b) `outlet_statPress` のまま Ps を実出口圧に合わせる (run_0214, case/45 流)。**node の超音速ノズルでは (a) を既定**にする
  (`run_user_profile.py` / 3D も同じ)。
- **残課題**: 最終列の壁向き Uy・P バンプ自体 (§2.10 の残り) は `outflow` でも僅かに残る可能性があり、出口積分量を使うときは
  最終列を避ける (§2.11 の方針どおり)。

## 3. 修正方針 (次セッションの作業項目)

(§2.6 で課題定義を再構成済み。旧 2./3. の仮説 — 壁 k/ω 行 Jacobian・壁隣接 MUSCL — は
反実仮想で棄却: 発散種は壁でなく軸であり、2 次は bndFirstOrder で解決済み。)

### 3.1 `nodeAxisDirichlet` 実装設計 (2026-08-03)

- **config**: `mesh.nodeAxisDirichlet` (int, 既定 0 = 従来ビット不変)。有効条件は
  `discretization==node && isAxisymmetric==1`。runner の node 既定で 1 を出す。
- **代表点** (`mesh.cpp`): 軸ノード A ごとに内部双対面の相手ノード I のうち radial 方向
  cos = (y_I−y_A)/|x_I−x_A| 最大の非軸ノードを選び `axis_rep_d` に格納 (SU2 Normal_Neighbor・
  SST 壁関数代表点と同型)。
- **状態ピン** (`axisymmetricSource_d.cu::enforceAxisMirror_d`): `assembleResidual` 冒頭
  (enforceWallNoSlip の直後・dependentVariables の前) で毎ステージ
  ro/roUx/roUz/roK/roOmega[A] ← [I]、roUy[A]=0、roe[A] ← roe[I] − 0.5·roUy[I]²/ro[I]
  (radial KE 除去)。以後の派生量・境界・勾配・フラックスは全てピン後状態を見る。
- **残差除外** (`zeroAxisAllResiduals_d`): assembleResidual 末尾で軸ノードの res_ro〜res_roe
  (+RANS 時 res_roK/res_roOmega) を 0 化 (rms 汚染防止 + explicit 経路の状態固定)。
- **陰解法整合** (`timeIntegration_d.cu` block-DPLUR): axis_flag が渡された場合の decouple を
  従来の roUy 1 行から**全 5 行単位行化**へ拡張 (現状 axis_flag は常に nullptr 渡しなので
  ビット不変)。起動側で `nodeAxisDirichlet==1` のとき `msh.axis_flag_d` を渡す。
  scalar-DPLUR 経路は状態ピンで実質担保 (dq が入っても次ステージ冒頭で上書き)。
- **保存性**: 軸 CV は解かなくなる。V ∝ r̄Δr Δx (r̄~1e-4 m) で全体質量への影響 <1e-6。
  EOS 床で既に非保存だった現状より悪化しない。
- **やらない**: species/凝縮スカラーの軸ピン (対象ケースに無し、必要時に拡張)、
  cell モード・平面 2D への影響 (ガードで不変)。

### 3.2 検証項目

1. 軸行の健全化: 12000 step 後に軸 P ≈ 第一内点 P (床ピン 0 ノード)、第一行 k シート消滅。
2. 課題 1 再現の解消: 1 次+陰解法 cfl_pseudo 2 が 12000 step で roK 成長なし。
3. 最終: node SST **2 次+陰解法** 12000 step `check_convergence.py` PASS + η_CF が cell
   (run_0002: 0.9790) と ~1%。
4. 非退行: 実証レシピ (explicit 1 次) run が同等以上に収束、case/29 `run_0038` verbatim 健全、
   cell / `nodeAxisDirichlet=0` はビット不変。
5. 課題 3 (細壁間隔) は本 plan では扱わず後継へ。

## 4. 影響・当面の運用

- **node は実証レシピ (1 次・explicit cfl0.1・warm start・粗め壁) で運用可能** — forge_design
  runner の node 既定に実装済み。2 次+陰解法が要る生産評価 (Phase 2 最適化) は課題 1–2 の
  解決までは cell を併用。
- 関連 run: `case/40.nozzle_design_tool/` README の run 一覧 (run_0005 発散記録 /
  run_0006 node 収束 PASS / run_0007 2次・陰解法の発散記録)。
- `interp_field.py` を node 対応に改良 (入力 h5 は CELLS/centCoords 優先、node res は
  ノード座標直用)。

## 5. 変更ログ

- `2026-08-03` — 起票。調査マトリクス完了 (§2)。
- `2026-08-03` — ユーザ指摘で case/29 node ノズルキャンペーンを発見し §2.5 で原因確定:
  リグレッションでも未踏領域でもなく、**壁第一セル過細 + レシピ不一致 (warm start /
  explicit / 1 次)**。実証レシピで node SST ノズル初の VERDICT PASS (run_0006)。残課題を
  「SST 陰解法 / 2 次精度 / 細壁間隔」の 3 点に再定義。
- `2026-08-03 (深夜)` — **課題 1・2 解決**。§2.6 で真因確定 (課題2 = bndFirstOrder 欠落設定の
  固有問題 / 課題1 = node 軸対称の軸行真空化 → 偽せん断 → SST k シート → 陰解法で遅発性発散)。
  §3.1 の `nodeAxisDirichlet` を実装 (solverConfig / mesh axis_rep / enforceAxisMirror +
  zeroAxisAllResiduals / block-DPLUR 軸 5 行 decouple / runner 既定)。**検証**: 軸床ピン
  21→0 ノード・k シート 1.3e4→3.4 に消滅 (case/40 run_0015)。**最終目標達成 =
  run_0016/0020 (2 次+陰解法 cfl_pseudo 2/4) 12000 step `check_convergence` ALL PASS +
  quasisteady ALL STEADY、12000 step ≈ 20 秒 (explicit 比 ~20 倍)**。η_CF=0.9907 は cell
  0.9790 と +1.2% (katoLaunder 起因でないことを run_0017 で確認、帰属は node 2 次化の離散差 —
  Phase 2 の Rao 照合で判定)。L スイープ node 取直し (run_0018/0019: 0.9835/0.9957) で単調応答
  維持。非退行: OFF 経路は atomicAdd ノイズ床内で不変、case/26 平板 node 500 step は表示桁一致。
  **残課題 = 課題 3 (細壁間隔 y+~1) と、§2 の出口コーナー・傾斜壁症状の再評価** (軸修正後に
  再現するかの確認) — 本 plan は active 継続。
- `2026-08-04` — ユーザ指摘 2 点を反映 (§2.7)。① `bndFirstOrder` の一律境界 1 次化は精度上
  不適切としてレシピから撤去: 課題 2 の発散種を再特定したところ**壁でなく軸帯**であり、
  `nodeAxisDirichlet` があれば全域 2 次が安定 (run_0021/0022 = 12000 step NaN なし・ALL STEADY,
  η=0.9896; L6 は PASS, L スイープ 0.9822/0.9896/0.9942 単調)。`gradLSQ: 2` 単独では防げない
  ことも確認 (勾配法は無関係)。② §2.6 の「nodeAxisDirichlet は SU2 と同型」を訂正: SU2 は
  同じ頂点中心だが planar+ソース方式で軸ノードを通常 DOF として解く。forge の r 重み幾何が
  軸半 CV を悪条件化しており、本フラグは SU2 より強い forge 固有の対症 (根治は将来課題)。
- `2026-08-11` — **課題 3 (細壁 y+≈1) 解決** (§2.9)。真因 = 双対 CV 重心/体積の float32
  shoelace 桁落ち (`gmshReader.hpp`, A 相対座標+double で修正) と `interp_field.py` の
  dst 照会点 (node は `MESH/COORD` を使うよう修正)。y+≈1 node low-Re が段階起動
  (陰 cfl0.5 1次 3000 step → 陰 cfl4 2次) で安定完走 (`case/40 run_0034`)。
- `2026-08-11 (2)` — **node 出口列欠陥 (§2.10) 根治** (§2.11)。真因 = `outlet_statPress_d` が
  bvar `Psb` を超音速分岐で更新せず、node 弱形式の境界圧力流束 `p_tilde` と GG/LSQ 勾配境界
  閉包が超音速流出でも規定背圧を読んでいた。半割 CV 離散収支の項別再構成で定量確定
  (バグ式で res 0.00% 閉合 / 修正式で欠損項が残差に現れる)。修正 = `Psb` の動的化
  (規定値は `P_exit_cfg` 引数へ分離)。検証 = run_0043/0044/0045: sag 19.4%→0.39% (y+1)・
  13.6%→0.17% (y+30)、η 出口積分=内部列一致、SU2 と +0.14% (y+30)。cell ビット不変。
  §2.10 の「内部列を正とする」運用を解除。
