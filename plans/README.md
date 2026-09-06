# forge 計画・設計判断インデックス (`plans/`)

`plans/` は forge の**変更単位 (変更計画・設計判断・検証ログ)** の文書を置く。`methods/` が「現在の仕様・解説」なのに対し、
`plans/` は「**なぜその設計にしたか / これから何を変えるか**」を担う。

計画の作り方は [`_template.md`](_template.md) を起点に、命名規約 `<area>-<short-slug>.md` で配置する。
運用フロー全体は [`../AGENTS.md`](../AGENTS.md) の「## 開発フロー」節を参照。

> **凡例 (旧パス)**: これらの計画は以前 GitHub 配下の `.github/plans/` に置かれていた (現在の root `plans/` とは別物)。
> 各計画の完了ログや起票履歴に出てくる `.github/plans/<name>.md` という記述は、**現在の root
> `plans/{active,accepted,archived}/<name>.md`** を指す (移行済み)。履歴の文言はそのまま残してあるので、
> 当時の記録として読むこと。

ライフサイクルは**フォルダ移動**で表す (ファイル名は固定し、移動でリンクを壊さない)。

- `active/` — 検討中・実装/検証中 (`draft` / `in_progress`)。**いま手を動かしている計画はここだけ**。
- `accepted/` — 実装・検証済みで、**いまも現在挙動を支配している設計判断** (`done`)。仕様は `methods/` 側、ここには「なぜ」を残す。
- `archived/` — superseded / 棄却 / 役目を終えた計画 (`superseded`)。後継を本文に明記する。

各計画の詳細 (背景・検証 run・結論) は当該 `.md` 本文を参照。本インデックスは航行用の一覧に徹する。

## active (検討中・進行中)

| Plan | area | 概要 |
| --- | --- | --- |
| [chemistry-finite-rate-h2.md](active/chemistry-finite-rate-h2.md) | `thermophysics / chemistry` | **有限速度化学 (H₂ 燃焼・ノズル化学非平衡)** (2026-09-04 起票): 種ブロック point-implicit + sensible datum 反応熱陽注入。Phase 0 (CEA スクリーニング・熱力学 DB ツール・Jachimowski YAML) 完了、Phase 1 ソース項から実装 |
| [architecture-bndfirstorder-removal.md](active/architecture-bndfirstorder-removal.md) | `architecture / boundary` | **`mesh.bndFirstOrder` の廃止**。使用禁止をルール化済 (粘性応力を破壊 + 疑似 2D で全域に効く)、コード削除が残タスク |
| [axisymmetric-freestream-hoop-gauge.md](active/axisymmetric-freestream-hoop-gauge.md) | `axisymmetric` | 軸対称 hoop ソースの自由流保持 (pRef ゲージ整合 + 離散閉性面積、倍精度不要) |
| [architecture-node-option-consolidation.md](active/architecture-node-option-consolidation.md) | `architecture` | node の `node*` オプション削減 — 整合セットを既定化し旧規約系フラグを撤去する手順 |
| [architecture-node-centroid-value-position.md](active/architecture-node-centroid-value-position.md) | `architecture` | node-centered の centCoords を「値の位置 (ノード座標)」に統一し、双対重心/軸半径を分離する |
| [convection-freestream-preserving-flux.md](active/convection-freestream-preserving-flux.md) | `convection` | 対流流束の free-stream 保存 (基準静圧差分) |
| [convection-keep-cb-pressure-correction.md](active/convection-keep-cb-pressure-correction.md) | `convection` | 高周波圧力欠陥 δp^HF 駆動の mass-flux 補正 (Rhie–Chow 型市松キラー, σ_min 独立で丘頂 spanwise 鋸歯を減衰) |
| [convection-keep-diss-lowmach-precond.md](active/convection-keep-diss-lowmach-precond.md) | `convection` | KEEP ES 散逸レイヤの低マッハ前処理化 (圧力ジャンプ増強で市松と解像物理を分離) |
| [convection-keep-revive-node.md](active/convection-keep-revive-node.md) | `convection / turbulence` | KEEP スキーム復活 (modern API・cell/node) + node WALE で LES/ILES |
| [discretization-lsq-gradient.md](active/discretization-lsq-gradient.md) | `gradient / architecture` | node-centered 最小二乗 (LSQ) 勾配 |
| [discretization-median-dual-3d.md](active/discretization-median-dual-3d.md) | `architecture` | 3D median-dual (M4): 3D 双対生成 + periodic 双対面 (3D node DDES/LES の起点) |
| [axisymmetric-su2-source-formulation.md](active/axisymmetric-su2-source-formulation.md) | `axisymmetric / discretization` | 軸対称を SU2 流 (planar 幾何+1/y ソース) に切替える `axisymMethod: 1`。node 軸行真空化の根治 (nodeAxisDirichlet 撤去) |
| [boundary-node-nozzle-wall-outlet-stability.md](active/boundary-node-nozzle-wall-outlet-stability.md) | `boundary / discretization` | node ノズル安定性。課題1/2/3 解決済 (§2.6-2.9)、**node 出口列不整合も根治 (§2.11, outlet Psb 動的化 — η 出口積分が正値に)**。壁 T エントロピー市松も解消を実測確認 (§2.8, 43.5→0.25 K)。残: 帰属の同一メッシュ A/B のみ |
| [discretization-node-boundary-ghostless.md](active/discretization-node-boundary-ghostless.md) | `boundary` | node-centered 境界のゴースト撤廃 (段階導入: まず壁, 次に流入出/slip) |
| [tooling-cloud-gpu-env.md](active/tooling-cloud-gpu-env.md) | `tooling / infra` | クラウド GPU 環境 (AWS EC2 spot): Docker 構築・計算投入・速度計測基盤。MPI 化の実行基盤を先行整備 |
| [tooling-nozzle-design-tool.md](active/tooling-nozzle-design-tool.md) | `tooling / optimization` | 超音速ノズル設計ツール親計画 (5 機種・forge 評価器・サロゲート MOO・帰還エンジン・確認 CFD メニュー・AI 対話問題定義)。フェーズごとに子 plan を起票 |
| [tooling-nozzle-phase3-windtunnel.md](active/tooling-nozzle-phase3-windtunnel.md) | `tooling / optimization` | ↑の Phase 3 子 plan (2026-08-13 起票): ①軸対称風洞 — モード F (中心線マッハ全権。**B8 は 2026-08-16 に生産の座を axis-Mach チェーンへ譲渡**、記録として維持)・**帰還エンジン v1 (逆MOC+δ\*経験式) / v2 (Euler帰還・凍結マップ)**・一様性 εM/εθ・自己一貫/Sivells 照合。後半 = co-kriging MOO + 凝縮確認メニュー |
| [tooling-nozzle-walldriven-chain.md](active/tooling-nozzle-walldriven-chain.md) | `tooling / optimization` | ①風洞の **wall-driven チェーン** (2026-08-15 起票、外部推奨方針の評価・採用): 円弧廃止 U→T→D 多項式壁 + CFD 特性抽出 + D 下流 wave-cancellation MOC + δ\*/B-spline。$M_c(x)$ は診断量。W1–W3 済、W4 は平面近似止まり (軸対称厳密化 2 経路を記録)。**W5 以降は保留** — 主線は [tooling-nozzle-axismach-chain.md](accepted/tooling-nozzle-axismach-chain.md) へ (2026-08-15 ユーザ判断) |
| [tooling-nozzle-sern-chain.md](active/tooling-nozzle-sern-chain.md) | `tooling / optimization` | ⑤ **SERN チェーン** (2026-09-04 起票、branch `feature/sern-design`): 燃焼器出口 starting line + 2D 平面最大推力理論の key point $(M_c,\theta_c,\dot m_c/\dot m)$ を dv にした逆 MOC + カウル角/長 + forge 2D RANS 多作動点評価 + MOO。旧 ⑤ (壁圧 Bézier dv・局所帰還・3D FFD in-loop) は撤回。出典: [sern-design-method-survey.md](../notes/investigations/sern-design-method-survey.md) |
| [turbulence-iddes-sst.md](active/turbulence-iddes-sst.md) | `turbulence` | SST-DDES / SST-IDDES 実装計画 |
| [turbulence-sst-thermal-flux-model.md](active/turbulence-sst-thermal-flux-model.md) | `turbulence / boundary` | SST 壁関数のエネルギー流束モデル置換 (Kader q_w)。等温壁×粗メッシュの熱負荷予測と T_aw 強閉包の前提 (in_progress: 平板合格 ±7%・Kader T⁺ 原式修正済。残 = T⁺ 圧縮性補正 [ベル +87% 実測]) |
| [turbulence-wmles-wall-stress.md](active/turbulence-wmles-wall-stress.md) | `turbulence / boundary` | WMLES 用代数壁応力モデル (Reichardt + Kader)。既存 SST 壁関数資産 (Normal_Neighbor / AddTauWall) を流用し τ_w/q_w で壁粘性流束を置換 |
| [turbulence-reichardt-gap-residual.md](active/turbulence-reichardt-gap-residual.md) | `turbulence / boundary` | 壁解像プロファイルが $y^+\approx30$ で Reichardt 相関を約 5% 下回る件。**SU2 (全残差収束) も 0.948 で forge と 0.1–0.2% 一致**するため forge 固有性・自由流乱流・補間手法は除外済み。**壁法則スイープで Spalding は逆に 1.020** (法則間の広がり 7%) = **壁法則選択のモデル形式不確かさが観測差と同程度**。ただし forge/SU2 とも SST なので**共通 SST 誤差は未除外**、Spalding の正しさも未証明。**forge の実装バグではないが壁関数の低摩擦バイアスには寄与する** (壁法則段 1.0425)。残る候補は有限 $Re$ ($Re_\theta$ 3600–5800) / 入口・前縁の扱い |

## accepted (現役の設計判断)

| Plan | area | 概要 |
| --- | --- | --- |
| [tooling-docker-image-split.md](accepted/tooling-docker-image-split.md) | `architecture / tooling` | **Docker イメージの多段分割とビルド時間短縮 (2026-09-06)**: `Dockerfile.cuda.dev` / `.cloud` を多段 `Dockerfile.cuda` (base / cloud / dev) へ統合し、apt・gmsh・ParaView を工程ごとのレイヤに分割 (gmsh 失敗で apt ~2 GB からやり直す状態を解消)。`.dockerignore` でコンテキスト 1.058 GB→6.6 kB、`tools/build.sh` に CUDA arch 1 種絞り込み + ccache (フルビルド 2:25→1:07、ccache 温 3.9 s)。**`libboost-dev` はどちらの Dockerfile にも未宣言で、ParaView / matplotlib の依存に相乗りしていただけ**と判明 (base に明示)。残: buildx 導入・`boost::split` 置換 |
| [tooling-nozzle-deltastar-core-matched-euler.md](accepted/tooling-nozzle-deltastar-core-matched-euler.md) | `tooling / design / boundary layer` | **排除厚さ更新 (2026-09-04 完了)**: CONTUR 運動量積分の初期壁 + 固定 Euler 基準・帯局所抽出 (BL 直外の帯で q_NS/q_E を線形フィット) の半径方向固定点反復。旧相関はスロート δ\* を 3〜12 倍過大に与え NS 質量流量 +0.8〜3.7 %・試験部 M −0.2〜−0.7 % の主因だった。case/45 M6: ṁ 比 0.9999・出口コア M +0.01 % を Md トリムなしで達成 (run_0020–0023)、case/42 M5・case/44 も設定不変で合格。Md トリム/law 帰還/x_lo は廃止 |
| [discretization-median-dual-2d-facevect-precision.md](accepted/discretization-median-dual-2d-facevect-precision.md) | `discretization` | **2D 双対幾何+CW 判定の桁落ち除去** (ローカル原点+double, 2026-08-31): 真犯人は makeMesh の float shoelace CW 判定 (スリバーで境界 surfVect 誤反転)。閉性 1e-7・第一セル 2.4 μm (y+1 真値) 変換可能に |
| [design-isobutane-wt-m5-sweep.md](accepted/design-isobutane-wt-m5-sweep.md) | `design campaign` | イソブタン風洞 **M5** R×L_U×L_c 27 点スイープ完了 (2026-08-30): quintic 不成立→knot MK2.5・L_c≥14、パレート 9 点、**推奨 R3/L_U9/L_c14** (4.37 m, dM 0.043 %M_d)。凝縮後段確認 +29 K で dry 確定 |
| [design-isobutane-m6-d155.md](accepted/design-isobutane-m6-d155.md) | `design campaign` | **M6 出口径 1.55 m 最短全長スタディ完了** (2026-08-31, case/45): 最短ロバスト R2/L_c39.7/M_K2.7 → **スロート→出口 7.337 m** (r_t 76.4 mm)。粘性トリム Md6.0144 で **NS dry 出口軸 M 6.0008**・δ\* 2巡固定点。**Tt 1600 K は凝縮不可避** (S≈14, 凝縮 ON 出口 M 5.927・軸勾配)、凝縮フリーは Tt≳1820 K |
| [condensation-equilibrium-eos.md](accepted/condensation-equilibrium-eos.md) | `condensation` | **平衡凝縮の EOS 拘束形** (`condEquilibrium: 2`, 2026-08-19 実装・検証済): 湿り度 g を dependentVariables で $(T,g)$ 同時反転 (括弧付き Newton) して `rog` に射影、輸送凍結。飽和セルで厳密 S=1・過冷却 0 (緩和形は onset 帯 ~2 r_t で S≤1.18)。下流は緩和形と固定点一致 (case/44 run_0098 vs 0097: |ΔM|≤1e-3) |
| [condensation-equilibrium.md](accepted/condensation-equilibrium.md) | `condensation` | **平衡凝縮** (2026-08-18 実装・検証済): 各セルで $p_v=p_{sat}(T)$ の $g_{eq}$ へ緩和 (`condEquilibrium`, 平衡専用 ΔT/Δg 律速)、モーメントはソース 0。診断 `condTsat_<s>` ($T_{sat}(p_v)$)。case/44 新条件 6 点で S→1 (onset 後 ~1 r_t)、出口凝縮率 9–45 %、M4.19 は CEA 平衡 8 % と整合 |
| [condensation-evaporation.md](accepted/condensation-evaporation.md) | `condensation` | 非平衡凝縮に**蒸発** (2026-08-18 実装・検証済, 既定 ON): S≤1 で統一駆動力 (Kelvin 項は正帰還回避で既定 off, 質量収支不変) の負成長を $r_{30}$・λ スケール (半減/step 律速) で適用し $r_{30}<2$ nm で一括消滅。0-D HK 解析解一致、Euler 圧縮帯で S 0.67→0.87/−4.4 K、Wyslouzil 回帰同一、NS 高温壁 (920 K) で T>400 K の液相 0。残: 二温度・亜臨界・液滴乱流拡散 |
| [tooling-nozzle-semiperfect-gas.md](accepted/tooling-nozzle-semiperfect-gas.md) | `tooling / optimization` | **semi-perfect gas (NASA-9/CEA, frozen) を MOC に実装** (2026-08-17): forge 内蔵 DB と同一係数、ν/A/A*/ρV をテーブル化、MOC カーネルは γ の位置にガスモデルを受ける。燃焼ガスは γ 1.31→1.38 で一定 γ の CPG は出口径 8% 誤る。**イソブタン風洞 (5 MPa/1000 K/M4.2/出口 1 m) の R/L スタディ**: 推奨 R=3, L_U=6, L_c=max (全長 4.35 m, 誤差 0.26%)、R=5 設計不能、L_c=8 は R≥2 でゲート外、L_U は出口一様性に効く。**申し送り: TP×node 軸対称は forge 側で発散** (axisymMethod:1 は CPG でも出口軸コーナーで発散 [case/41 run_0076])。run case/42 0001–0019 |
| [tooling-nozzle-axismach-knot-law.md](accepted/tooling-nozzle-axismach-knot-law.md) | `tooling / optimization` | axis-Mach **A6 = knot 付き区分 C² 軸 M 則** (2026-08-16): 単一 quintic の $L_c$ 上限 (M6/R3 で 17.9) では終端特性線の戻り 55 $r_t$ に対し膨張がスロート 2 $r_t$ 内に押し込められ θ_w>μ_w で逆 MOC が fold → M6 全点 REJECT。knot (DOF $L_c$, $M_K$) で $L_c$ 30–60 が成立、軸点列スロート側細分 `axis_dx0` も追加。**M6 イソブタン (Tt 1550 K, 出口 1 m) スタディ 11 点全て 0.03–0.08% $M_d$**、推奨 R 2–3 / L_c 45–50 / L_U 12 (全長 ≈ 6 m)。副産物: TP `inlet_Pressure` BC の静圧参照ブレンドを撤去 (残差床 3e-2→9e-7)。run case/42 0039–0061 |
| [tooling-nozzle-tp-split-h2o-condensation.md](accepted/tooling-nozzle-tp-split-h2o-condensation.md) | `tooling / condensation` | 燃焼ガス TP を **H₂O 以外=擬似種 `MIXDRY` + H₂O 独立種** に分割 (`evaluate.tp_species: split_h2o`) し既存 H₂O 凝縮 (Kw+HK) を有効化 (2026-08-17)。dry で pseudo と同一 (M4.2 軸 M 差 1.7e-5)、Wyslouzil で既存 CPG 凝縮を再現 (cell; node は非粘性で cell と一致)、イソブタン M4.2 H₂O 5 % で凝縮 (g 0.0126, 出口 M 4.20→3.94)。**申し送り: node × 平面壁クラスタ no-slip は壁 T>Tt で unstart** |
| [tooling-nozzle-axislaw-onepoint.md](accepted/tooling-nozzle-axislaw-onepoint.md) | `tooling / optimization` | axis-Mach **D 案 = 端点アンカー固定 + 内部補間点 1 点 (ξ_P, η_P) の C⁴ 区分 5 次軸 M 則** (2026-08-16): 単調可能領域は細い尾根 (L_c=60 でほぼ消滅)。μ−θ 余裕で A を上回る候補は全て内部プラトー (M'' 反転 3) 由来で、A と同じ単峰条件では全 L_c で A に劣る (μ−θ<0, θ_max 21–27°)。トポロジ余裕 (min sin 角) も A の半分。**D 不採用、CFD なし、生産 A 維持** |
| [tooling-nozzle-topology-crossing-speedup.md](accepted/tooling-nozzle-topology-crossing-speedup.md) | `tooling / optimization` | 逆 MOC 非隣接特性線交差検査を AABB sweep 化。旧判定と一致のまま600/1200/2400点で75/149/454倍、候補process並列込み $M_K$ 感度を32.4 sで完走 |
| [tooling-nozzle-shortest-robust-axis.md](accepted/tooling-nozzle-shortest-robust-axis.md) | `tooling / optimization` | M6/R3の軸則Aを品質制約下で$x_F$最小化。600点探索+境界のみ1200点で、最短合格$L_c=39.05,M_K=2.76$を確定 |
| [tooling-nozzle-axismach-length-dv.md](accepted/tooling-nozzle-axismach-length-dv.md) | `tooling / optimization` | **全長 dv 化 `Lc_mode: from_length`** (2026-08-31): dv を $L_c$ からスロート→物理出口全長 `L_total` ($=x_F$) に替え、軸 M 極大点 $x_E$ は終端特性線込み設計パスから逆算 (解析初期推定+固定点反復、実測 1–2 パス)。離散 MOC の $x_F$ ノイズ床 (n500 で ~0.05 $r_t$) の発見と tol 設計。explicit/max は不変・往復 bit 同一。テスト `run_axismach_length_tests.py` |
| [tooling-nozzle-axislaw-smoothness.md](accepted/tooling-nozzle-axislaw-smoothness.md) | `tooling / optimization` | axis-Mach 軸 M 則の滑らかさ比較 (2026-08-16): A (現行 knot, $M'''$ が knot でジャンプ) vs B (C⁴ monotone B-spline) vs C (非負 $d\nu/dx$ B-spline)。**A・B が hard gate 合格、C は特性線網の恒常的 fold で不合格**。B は $M'''$ を5倍改善するが θ_max/margin が悪化し CFD 誤差が A の 4.4 倍 — 軸則の滑らかさは壁 fairness に直結せず、**生産は A を維持**。run case/42 0051(A)/0062(B) |
| [tooling-nozzle-axismach-physical-throat.md](accepted/tooling-nozzle-axismach-physical-throat.md) | `tooling / optimization` | axis-Mach **A13/A14** (2026-08-17, ユーザ+Codex レビュー採択): **A13 = 上流履歴込み δ\* + 真の幾何スロート探索 + 上流 Hermite 再生成** (`PhysicalNozzleWall`) — 採用 (NS overshoot +0.30→+0.03%、残差は設計側の谷のみ)。**A14 = 制約付き LSQ B-spline を本流表現として A/B** — 曲率振動 50 倍改善だが壁角乖離 0.34° で軸 M 0.24→0.95% と悪化、**CFD で棄却** (補間壁維持、コード残置)。run_0072–0075 |
| [tooling-nozzle-axismach-viscous-deltastar.md](accepted/tooling-nozzle-axismach-viscous-deltastar.md) | `tooling / optimization` | axis-Mach **A12 = 粘性 δ\* 補正** (2026-08-16): 物理壁 = inviscid (cplus) + δ\* 法線オフセット。RANS チェーン (coarse 中継 43 s → y+1.5 SST 162 s) を `prepare_ns`/`run_staged_ns` に実装。**相関 δ\* で十分** (CFD 抽出との差 ~10% は軸 M に無影響、δ\* は 1 反復で固定点 [比 0.995–1.010])。無帰還 RANS 軸 M **0.533% Md** (B8 系 1.2% の半分以下)。残差は δ\* 非表現 (law 側帰還が future work)。run_0069–0071 |
| [discretization-moc-axisymmetric-source-term.md](accepted/discretization-moc-axisymmetric-source-term.md) | `tooling / optimization` | **MOC の軸対称源項 $\sin\theta/r$ を点自身で評価する (A11, 2026-08-16)**。旧実装は「2 点のうち軸から遠い方」で評価しており台形則の片端が置換され、**全体が 1 次精度**だった。修正で放射源流の厳密解誤差が 10–40 倍減・収束次数 1.02→1.7。**A5 の「質量流束リーク」の正体**で、`n_axis_inv: 2000` を要求した根拠が消えた (流線壁と C⁺ 閉包の $r_F$ が 5 桁一致)。CFD: 旧生産比 ‖ΔM‖∞ 0.353→**0.224%** かつ設計 31.0→2.6 s。近軸は `AXIS_LIMIT_FRAC` でガード(生産経路では発火 0)。外部レビュー指摘を本セッションで独立検証。run_0056/0057 |
| [tooling-nozzle-moc-wall-unit-process.md](accepted/tooling-nozzle-moc-wall-unit-process.md) | `tooling / optimization` | axis-Mach チェーン **A10** (2026-08-16): 逆問題に欠けていた**壁点の閉包**を導入 — 壁 = 各 C⁺ 線上で累積質量流束が $\dot m^*$ に達する点 (古典 CONTUR/Zucrow–Hoffman 流)。**三角充填の反対角線がそのまま C⁺ 線**なので追加計算ゼロで読め、Delaunay も流線 ODE も不要になる。厳密解 (放射源流) で 1 次収束・バイアス床なし。**CFD: n=4000 で ‖ΔM‖∞ 0.353→0.299% $M_d$ かつ設計 31.0→10.9 s** (同一 n=2000 なら 2.6 s で実質同等 = 12 倍速)。生産既定を cplus + n4000 へ変更。run_0053–0055。未検証 WIP の `_CPlusMarch` は削除。残件: 出口 ε_M が ~2 倍 |
| [tooling-nozzle-axismach-throat-characteristic.md](accepted/tooling-nozzle-axismach-throat-characteristic.md) | `tooling / optimization` | axis-Mach チェーン **A8** (2026-08-15 起票・同日完了): 逆 MOC の初期値線を $M=1.05$ 縦線から**スロート特性線** (CONTUR 流) へ差し替え、壁を T から全域 MOC 出力にした (骨接放物線区間が不要に)。**到達不能帯 1.533 → 0.0019 $r_t$ で消滅**し、**アンカー更新なしの 1 パスで 0.5% $M_d$ ゲート通過** (‖ΔM‖∞ 0.353%・うねり 0.00043・出口 ε_M 0.020%、いずれも縦線対照と親計画 3 パスより良い)。run_0049/0050。あわせて充填ベクトル化 + Delaunay 共有で設計 1 回 400 s → 33 s (壁はビット一致) |
| [tooling-nozzle-axismach-chain.md](accepted/tooling-nozzle-axismach-chain.md) | `tooling / optimization` | ①風洞の **axis-Mach チェーン = 現行生産チェーン** (2026-08-15 起票、2026-08-16 に B8 を正式に置き換え): Hall–KL 遷音速 (`HallThroat`, CONTUR App.A 系列・検算済) + **5次 Hermite 軸 Mach law (自由度 $L_c$ のみ — B10-c を構造的に解消)** + 逆 MOC + CFD アンカー更新反復。**pass2/3 で ‖ΔM‖∞ 0.451/0.463% Md ≤ 0.5% ゲート** (B8 の 1/3 パス数)・出口 ε_M 0.05%/ε_θ 0.012°・**軸うねり 0.0019 = 円弧 R=2 の 1/11 (接合こぶ消滅)**。知見: 単調ゲートは $L_c$ **上限**を縛る (窓 9.57 → x_E≈9.6)・逆 MOC 質量流束リークは離散化 (生産 n_axis=2000)。run_0044–0047。A7 (粘性) は後続 plan |
| [tooling-nozzle-moo-loop.md](accepted/tooling-nozzle-moo-loop.md) | `tooling / optimization` | ノズル設計 Phase 2 (2026-08-13 完了): ③ベル **TOP 幾何 dv** の MOO ループ — 自作 EHVI (MC 照合 ≤0.1%・ZDT1 99.96%)・バッチドライバ (2 段起動・物理ゲート)・50 評価 28 分で前線 19 点 (ポリッシュ 19/19 settled)・**Rao 照合合格** (チャート点が前線上 +0.05%) |
| [tooling-nozzle-phase0-foundation.md](accepted/tooling-nozzle-phase0-foundation.md) | `tooling / optimization` | ノズル設計 Phase 0 (2026-08-13 完了): 問題定義 YAML・区分構成ジオメトリ・TFI→msh4.1→forge メッシュ・バッチ評価 CLI・推力メトリクス。E2E 基準 = run_0071/0072。Kliegel–Levine 高次は③2巡目へ後送 |
| [turbulence-sst-su2-taw-coupling.md](accepted/turbulence-sst-su2-taw-coupling.md) | `turbulence / boundary` | SST 断熱壁の熱的閉包 **mode 3 defect-flux を正式採用** (2026-08-11): Couette 恒等式 q+τu=0 に基づく保存的 W–I 全エネルギー流束 $H_T(T_{aw}-T_W)$ で T[W]→Taw を残差の解として実現 (壁温 1418K プラトー・SU2 実状態 10K 一致・η 不変)。mode 2 (SU2 overlay) は未採用の失敗記録。follow-up: ω irep ピン正式化・壁 μt 監査 |
| [turbulence-sst-adiabatic-taw-fluxmodel.md](accepted/turbulence-sst-adiabatic-taw-fluxmodel.md) | `turbulence / boundary` | SST 断熱壁 T_aw の W-I 流束注入初回試行は**実測発散で棄却** (壁ノードが EOS 床まで異常冷却; root cause は未確定で壁半 CV 全体が調査対象)。現行 = output-only fallback (Tsb/Taw_diag 出力専用, 生産 baseline)。後継: [turbulence-sst-su2-taw-coupling.md](accepted/turbulence-sst-su2-taw-coupling.md) (mode 3 defect-flux 採用で決着) |
| [architecture-node-boundary-gradient-dof-only.md](accepted/architecture-node-boundary-gradient-dof-only.md) | `architecture / discretization` | node 境界勾配 (GG/LSQ) から bvar を排除し owner-state のみに統一 (node outlet P/T interior 化の一般化)。outlet 非退行・cell 不変を検証済。turbulence-sst-adiabatic-taw-fluxmodel (未完了) の前提 |
| [architecture-median-dual-3d-double-geometry.md](accepted/architecture-median-dual-3d-double-geometry.md) | `architecture / discretization` | 3D median-dual 幾何の Newell ローカル原点化+境界蓄積 double 化 (堅牢化)。**監査結論: 3D は元から double 演算で 2D のような実害なし** (露出見積もりを訂正)。wall_dist 定義は 2D と一貫 (双対重心間距離) |
| [architecture-limiter-negative-skip-fix.md](accepted/architecture-limiter-negative-skip-fix.md) | `architecture` | `limiter: -1` (off) が早期 return を素通りし Venkatakrishnan フル計算 (KEEP では未使用) に落ちるバグ修正。KEEP 系 run 全体で ~20% 高速化・挙動不変 |
| [turbulence-sst-omega-crossdiff-jacobian.md](accepted/turbulence-sst-omega-crossdiff-jacobian.md) | `turbulence / time_integration` | SST ω 交差拡散の point-implicit Jacobian (sstCrossDiffJac): dual-time サブ反復の ω 収束改善 (1.5x→2.7x)、収束解は不変 |
| [timeint-bodyforce-massflow-control.md](accepted/timeint-bodyforce-massflow-control.md) | `time_integration` | 体積力の質量流量一定制御 (Benocci & Pinelli 圧縮性版, γ=0.02): 周期丘 DDES 駆動。case/39 で検証済 |
| [boundary-cell-periodic-conservation.md](accepted/boundary-cell-periodic-conservation.md) | `boundary` | cell 全周期境界の非保存バグ修正: device `bint_d[partnerCellID]` 未転送で ghost が誤値→運動量注入。setPeriodicPartner 直後に H2D コピーで根治。TGV で cell/node 一致を確認 |
| [thermophysics-eos-positivity-floor-config.md](accepted/thermophysics-eos-positivity-floor-config.md) | `thermophysics` | EOS 正値化フロア (pMin/roMin/tMin) の config 化。無次元・低圧ケース (Taylor-Green) で既定 1.0 Pa フロアが場を破壊する問題を解消、既定値据え置きでビット不変 |
| [boundary-inlet-profile.md](accepted/boundary-inlet-profile.md) | `boundary` | 入口分布プロファイル: CSV テーブルで inlet bvar を非一様化 (x/y/z 1D 線形 or xyz 最近傍)、壁法則 helper |
| [turbulence-node-wall-function-coverage.md](accepted/turbulence-node-wall-function-coverage.md) | `boundary` | node SST 壁関数の生産置換を第一内層ノードにも適用 (近壁 k 暴走修正、cell 不変・x_R が SU2 整合) |
| [architecture-axisym-axis-singularity.md](accepted/architecture-axisym-axis-singularity.md) | `architecture` | 軸対称 近軸の数値問題 (軸中心 k スパイク) の根本原因特定 |
| [architecture-axisym-nozzle-geometry.md](accepted/architecture-axisym-nozzle-geometry.md) | `architecture` | architecture-axisym-nozzle-geometry |
| [architecture-axisym-sst.md](accepted/architecture-axisym-sst.md) | `architecture` | architecture-axisym-sst — 軸対称 SST 幾何項 子計画 |
| [architecture-axisymmetric.md](accepted/architecture-axisymmetric.md) | `architecture` | architecture-axisymmetric — 軸対称ソルバ Phase 1 |
| [architecture-detect-nan.md](accepted/architecture-detect-nan.md) | `architecture` | detectNaN — NaN/Inf 検知・自動ダンプ診断オプション |
| [architecture-rans-sst.md](accepted/architecture-rans-sst.md) | `architecture` | architecture-rans-sst — RANS SST 親計画 |
| [architecture-residual-monitor-async.md](accepted/architecture-residual-monitor-async.md) | `architecture` | 残差モニタの device 常駐化 (per-step host 同期の除去) |
| [architecture-perphase-profiling-hotspot.md](accepted/architecture-perphase-profiling-hotspot.md) | `architecture` | per-phase 計測駆動の host 律速最適化 (13.4→5.9s/2000step, −56%)。残レバー (CUDA Graph/renumbering) は不採用理由ごと記録 |
| [boundary-outlet-characteristic-backflow.md](accepted/boundary-outlet-characteristic-backflow.md) | `boundary` | 静圧固定流出 BC の特性ベース化と逆流統一 |
| [convection-keep-es-dissipation.md](accepted/convection-keep-es-dissipation.md) | `convection` | KEEP 用 entropy-stable 散逸レイヤ (scalar/matrix, CPG/TP/多成分, σ較正済) |
| [convection-keep-diss-recon-jump.md](accepted/convection-keep-diss-recon-jump.md) | `convection` | KEEP ES 散逸の再構成ジャンプ化 (keepDissJump: 市松無傷で解像スケール保護, node LES 推奨構成) |
| [turbulence-wale-fix.md](accepted/turbulence-wale-fix.md) | `turbulence` | WALE 不活性バグ (壁なしメッシュ wall_dist≡0) + Sd テンソル誤式の修正; 64³ TGV では WALE off の ILES+ES が最良と判明 |
| [turbulence-sigma-model.md](accepted/turbulence-sigma-model.md) | `turbulence` | σ-model SGS (Nicoud 2011, LESmodel:2) 追加; 64³ TGV では WALE 同様過散逸で ILES+ES が最良のまま |
| [convection-multispecies-contact-pressure.md](accepted/convection-multispecies-contact-pressure.md) | `convection` | 多成分 TP contact 圧力振動: S2 (`speciesFaceReconstruction:1`, ρ-limiter face 組成) を production 採用。S3/中心補間/PEP は棄却 |
| [convection-slau2-lowmach.md](accepted/convection-slau2-lowmach.md) | `convection` | SLAU2 圧力流束 (低マッハ圧力–速度カップリング是正) |
| [diffusion-node-boundary-real-distance.md](accepted/diffusion-node-boundary-real-distance.md) | `diffusion` | node 境界半割面 拡散の実距離 over-relax 化 (∇·S 弱形式の置換) |
| [diffusion-node-wall-viscous-distance.md](accepted/diffusion-node-wall-viscous-distance.md) | `diffusion` | node-centered 粘性壁 flux の法線距離修正 (mirror-ghost 退化対策) |
| [diffusion-turbulent-thermal-conductivity.md](accepted/diffusion-turbulent-thermal-conductivity.md) | `diffusion` | 乱流熱伝導 (turbulent thermal conductivity) の追加 |
| [diffusion-viscous-shear-flux.md](accepted/diffusion-viscous-shear-flux.md) | `diffusion` | 粘性せん断フラックスの修正 (内部面の法線項 + 壁面 no-slip 抗力) |
| [discretization-median-dual.md](accepted/discretization-median-dual.md) | `architecture` | cell/node (median-dual) 両対応化の親計画 (M1–M3 完了)。現役の入口は子 plan (median-dual-3d / ghostless / centroid-value-position / lsq-gradient) |
| [discretization-node-wall-implicit-dirichlet.md](accepted/discretization-node-wall-implicit-dirichlet.md) | `boundary / time_integration` | node-centered 壁 no-slip の implicit Jacobian Dirichlet 化 (SU2 DeleteValsRowi 相当) |
| [output-node-wall-surface-viz.md](accepted/output-node-wall-surface-viz.md) | `architecture` | node 壁サーフェス出力の修正 (primal 面トポロジ退避 + h5 スキーマ + 壁関数 twall の ρu_τ² 化) |
| [condensation-nonequilibrium.md](accepted/condensation-nonequilibrium.md) | `condensation` | 非平衡凝縮 (4 モーメント)。N2=Arthur Fig.2 1–2% / H2O=Wyslouzil Fig.3 ~5% 一致。TP<200K は棚上げ (CPG が正) |
| [thermophysics-multicomponent-tpgas.md](accepted/thermophysics-multicomponent-tpgas.md) | `thermophysics` | 多成分 thermally-perfect gas ソルバ機能 |
| [tooling-quasisteady-check.md](accepted/tooling-quasisteady-check.md) | `architecture` | 準定常確認ツール check_quasisteady.py と運用ルール (残差収束 ≠ 量の定常化) |
| [thermophysics-species-implicit-coupling.md](accepted/thermophysics-species-implicit-coupling.md) | `time_integration` | 化学種の陰解法緩和整合 (matched-relaxation scalar-DPLUR) |
| [time_integration-explicit-pointimplicit-sst.md](accepted/time_integration-explicit-pointimplicit-sst.md) | `time_integration` | 陽解法 RK のスカラー (k/ω) 源項 point-implicit 化 |
| [time_integration-general-eos-jacobian.md](accepted/time_integration-general-eos-jacobian.md) | `time_integration` | block-DPLUR の一般EOSヤコビアン化 (TP gas の陰解法律速の根治) |
| [gpu-implicit-plan.md](accepted/gpu-implicit-plan.md) | `time_integration` | GPU 陰解法基盤 (block-DPLUR / 軸対称 / SST point-implicit / scalar 版 / dual-time BDF2)。親計画 |
| [time_integration-implicit-stable-cfl.md](accepted/time_integration-implicit-stable-cfl.md) | `time_integration` | 陰解法 (block DPLUR) の安定 cfl_pseudo 引き上げ |
| [time_integration-lowmach-preconditioning.md](accepted/time_integration-lowmach-preconditioning.md) | `time_integration` | 低マッハ前処理 (Weiss–Smith) — 密度ベース経路の散逸スケール是正と陰解法固有値前処理 |
| [time_integration-line-implicit.md](accepted/time_integration-line-implicit.md) | `time_integration` | **line-implicit (壁法線 block-Thomas)** 実装 (2026-09-02, `lineImplicit` opt-in): 構造化で 1250 本/被覆100%。**M6 ノズルでは cfl 上限不変** (律速は streamwise lag と判明、壁法線説は否定) — 収束 −14 %/step のみ。lu5 ピボット罠 (LASWP 先行必須)・device printf 引数上限罠を記録。次の本丸は streamwise 第 2 ライン族/LU-SGS 順序 |
| [time_integration-line-implicit-viscous-v2.md](accepted/time_integration-line-implicit-viscous-v2.md) | `time_integration` | **line-implicit v2 (factor/solve 分離 + `lineKFreeze` + 粘性結合 + dt 割引)** (2026-09-02): v1 の 2.44× コストの主犯は Thomas 毎 sweep 再分解 → 分離+凍結で 1.32×。壁法線律速は音響で粘性結合は僅差。**DDES dual-time では pseudo-CFL 上限を引き上げ (point cp2 発散 / line cp8 安定)、cp4+nSub13 で ctrl の 0.88 倍時間・同品質**。M6 定常 (streamwise 律速) とは律速モードが違う |
| [time_integration-update-positivity-guard.md](accepted/time_integration-update-positivity-guard.md) | `time_integration` | **陰的更新の正値性ガード (負の結果, 2026-09-02)**: commit 時の局所 under-relax (`updateGuardAlpha`, opt-in) は CFL 上限を上げない — 上限の真因は defect-correction 反復不安定で床 NaN は終端症状。lowMachPrecond=2 併用も不成立 (SST チャネルで爆発)。生産は cfl8+implicitRelax0.7 |
| [time_integration-scalar-dplur-axisym-source.md](accepted/time_integration-scalar-dplur-axisym-source.md) | `time_integration` | scalar DPLUR の軸対称ソース Jacobian 整合 (TP 対応) |
| [turbulence-enhanced-wall-treatment.md](accepted/turbulence-enhanced-wall-treatment.md) | `turbulence` | SST Enhanced (Automatic / y⁺ 非依存) Wall Treatment |
| [turbulence-kato-launder.md](accepted/turbulence-kato-launder.md) | `その他 (turbulence)` | SST 生産項 Kato–Launder 補正 (`katoLaunder`) |
| [turbulence-node-sst-wallfunction.md](accepted/turbulence-node-sst-wallfunction.md) | `turbulence` | node-centered SST 壁関数の代表点修正 + τ_w 付与 (ic=壁ノード u=0 トラップ) |
| [turbulence-node-inlet-dirichlet-conserved.md](accepted/turbulence-node-inlet-dirichlet-conserved.md) | `turbulence` | node 入口 Dirichlet の保存量整合+残差除外 (rms_roOmega プラトー根治、cell 不変、x_R が実験に接近) |
| [turbulence-sst-freestream-init.md](accepted/turbulence-sst-freestream-init.md) | `turbulence` | SST 初期乱流の freestream 初期化 (kInit / omegaInit) |
| [verification-passive-pseudoshock-control.md](accepted/verification-passive-pseudoshock-control.md) | `verification` | 多孔壁パッシブコントロールの逆解析 (case/36)。2D で中心軸波状 ~30% 低減=ショックレス化を定性再現 (Phase 1 完了) |
| [turbulence-node-wf-representative-point.md](accepted/turbulence-node-wf-representative-point.md) | `turbulence / boundary` | node 壁関数の代表点 (Normal_Neighbor) 診断。$\omega$ 側を出し切っても両ケースの符号が逆 (case/26 不足 / case/40 過剰) なので $u_\tau$ の入口側へ。**幾何 / 代表点速度 / 壁解像基準の同一位置サンプリング**の 3 分離が先、アルゴリズム変更は後 **完了 (2026-08-13)**: 候補 A (法線化) は凍結場で $\tau_w$ 比 0.610 と逆効果のため不採用、候補 B (自由流乱流) は棄却、SU2 クロスチェックで forge 固有性を除外。残件は Reichardt 差 plan へ |
| [turbulence-node-wf-omega-source.md](accepted/turbulence-node-wf-omega-source.md) | `turbulence / boundary` | node 壁関数経路の低摩擦解。$\omega$ 介入 3 種 (pin / limiter bypass / E3) を分離検証したが**両ケースのゲートを同時に満たせず** (case/26 は 0.84–0.86 止まり、E3 は case/40 で 1.1414)。E3 は不採用・実験経路として残置。後継 = 代表点診断 plan |

## archived (superseded / 終了)

| Plan | area | 概要 |
| --- | --- | --- |
| [tooling-nozzle-moc-flux-closure-wall.md](archived/tooling-nozzle-moc-flux-closure-wall.md) | `tooling / optimization` | 逆 MOC の壁を「流線積分」から「断面の質量流束閉包」へ替える案 (A9)。**CFD で棄却** (2026-08-16): 同一場で `wall_mode` だけ振ると ‖ΔM‖∞ 0.353→0.507% $M_d$・出口 ε_M 0.020→0.057% と悪化 (run_0050 vs run_0052)。実装は正しい (一様流 1.7e-7) が、**「解像度非依存」はバイアス床**で、放射源流の厳密解では流線が細かい場で下回る (6.1e-4 vs 1.4e-3)。MOC 場は離散レベルで質量保存を満たさないため、閉包を変えても不整合は壁半径から壁形状へ移るだけ。`mdot_ratio_moc` が循環指標になる罠も記録。副産物: `dx_wall` は 0.02 でよい |
| [diffusion-node-scalar-nonortho-limit.md](archived/diffusion-node-scalar-nonortho-limit.md) | `diffusion` | node-centered flux の dcc に node 座標を使う (双対 CV 重心由来の見かけ非直交を排し SST omega 爆発を根治) |
| [precision-mixed-axisym.md](archived/precision-mixed-axisym.md) | `precision` | 混合精度 (iterative refinement) で軸対称 近軸の陰解法を root-fix する |
| [turbulence-sst-thermal-wall-function.md](archived/turbulence-sst-thermal-wall-function.md) | `turbulence / boundary` | SST 壁関数の熱的閉包 (Crocco 型 T_aw「弱閉包」旧版)。**superseded**: 境界勾配が当時 bvar を読んでいたため場に触れる経路が実在し、かつ実データでは壁ノード実状態がほぼ動いていなかった (「4K 一致」は出力 Taw 診断値同士)。後継: `accepted/turbulence-sst-adiabatic-taw-fluxmodel.md` |
