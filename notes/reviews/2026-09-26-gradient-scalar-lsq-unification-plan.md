# codex レビュー: gradient-scalar-lsq-unification (plan)

- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **stage**: `plan`
- **date**: 2026-09-26
- **commit**: `56067325` (feature/sern-design)
- **codex**: effort `high`, 6.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m2
- **focus**: (1) cInt 共有と差分形の汎用 LSQ gather で k/ω・Y・ξ・モーメントを置き換える設計の穴 (境界・ghost を読む側・周期 gather の前提・軸対称・壁ピン節点) (2) opt-in → 本 plan 内で既定化する切り替えと provenance/manifest (3) §6 S0–S3 の合否が測る前に固定されていて検出力があるか (特に S2 の『既存ゲート + 格子対』規則と、ケースの起点・量の選び方) (4) 既定設定で生きているのは k/ω だけという前提
- **extra**: `plans/accepted/boundary-node-periodic-gradient-fix.md`, `methods/gradient.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

`cInt` を共有する差分形 LSQ は支持します。既存の周期修正とは対象が分かれており、cell を GG に残す方針も現行ルールと整合します。  
ただし、SST の利用先・壁 `F1` の説明、S0/S2 の判定条件、既定化時の履歴管理に修正が必要です。現状の検証仕様のまま実装を始めるべきではありません。

1. **Major — SST の影響評価に、コードと矛盾する前提がある。**

   **根拠:** [ransBoundary_d.cu:77](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransBoundary_d.cu:77) は低 Re 壁で `k=0` を代入します。[ransSource_d.cu:340](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:340) の F1 式では、通常 SST・有限入力なら `arg1_c=0`、したがって **`F1=0`** です。CPU float32 で式を再評価しても `arg1_c=0, F1=0` でした。§4.2 と S2 の「壁ノード F1 が両経路とも 1」は誤った期待値です。

   また、[ransSource_d.cu:259](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:259) は `axisymMethod=1` で `dKdy/dOmegady` を半径方向の移流拡散ソースに直接使用します。「`CDkw/F1` だけ」ではありません。拡散面係数も単純な F1 平均ではなく、各側の `σ(F1)μt` を含む係数の補間です（[scalarTransport_d.cu:116](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/scalarTransport_d.cu:116)）。

   **対案:** F1 自体の仕様変更は混ぜず、現行式に基づく期待値へ訂正する。壁・第一内層で `k/ω`、勾配、`F1`、実際の拡散係数を同一状態で照合する。軸対称は `axisymMethod` を明記し、method 1 の半径方向ソースを作用素試験に追加する。

2. **Major — 周期 group の状態同値性が、勾配を読む時点で保証されていない。**

   **根拠:** 初期化の [main.cpp:1257](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1257) は `periodicMirrorNSState` を呼びますが、その追加配列には `roK/roOmega` が含まれません。SST の専用ミラーは更新後に行われます。また、勾配前の壁 BC は各部分 stencil から求めた `wall_y_eff` で壁 `omega` を再設定します（[ransBoundary_d.cu:153](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransBoundary_d.cu:153)）。非対称 stencil で局所最短距離が `h` と `2h` なら、同じ粘性でも壁 `omega` は4倍違います。保存量の更新後ミラーだけでは、この再設定後の同値性を保証できません。

   さらに [periodicNode_d.cu:173](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:173) の NS・`Y`・受動種 gather は回転周期を除外しません。「登録は現行どおり」と「回転周期では gather しない」は、そのままでは両立しません。

   **対案:** 初回・restart・BC 適用後の**勾配入力 primitive**について group 同値性を検査する。壁∩継ぎ目では合併 stencil に基づく壁距離など、BC の共通値を決めてから同期する。スカラー gather の適用条件は係数合併と同じ関数に統一する。非対称壁 stencil と、意図的に非同期にした restart を S0 に入れる。

3. **Major — S0 は参照作用素と浮動小数点の判定仕様が一致していない。**

   **根拠:** [gharness.py:185](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/gharness.py:185) の `lsq_merged_ref` は `root` ごとに行列を合併して全 member に返します。§4.3 の軸対称×周期は**片側 LSQ**なので、この参照をそのまま使えません。

   また、局所 gather は非 atomic でも、周期合算は [periodicNode_d.cu:57](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:57) の `atomicAdd` です。NS とスカラーで同じ加算順になる保証はなく、相殺する成分の差を「2 ulp以内」とは保証できません。float32 の部分和 `1e8, −1e8, 1` は加算順だけで `1` と `0` になります。

   壁で非零定数 `k/ω` を焼いて通常 BC を通す試験も、壁ピンが定数場を壊します。既存ハーネスはその stencil を除外しています。

   **対案:** S0 を「BC 前の純粋作用素」と「BC 後の実入力」に分ける。前者で全節点の定数保存を検査し、後者は実際に読んだ場を参照へ渡す。軸対称×周期には非合併参照を用い、非線形場で識別する。ビット一致は局所 gather に限定し、周期合算は部分和の絶対値和に基づく誤差限界で判定する。モーメントを含む **N=5以上**でチャンク境界も検査する。

   NS と同じ結果であることは、ジッタ系列の次数未達を独立に解決した証拠にはなりません。§6末尾は「NS と同じ既知制約を継承」とするのが妥当です。

4. **Major — S2 の `--from-floor` は、作用素切り替えの収束判定として不適切。起点の成立確認も不足している。**

   **根拠:** [check_convergence.py:321](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:321) は、末尾平均だけでなく**全期間ピーク**が参照床の既定1.5倍以下かを判定します。GG の固定点から LSQ へ切り替えて残差が跳ねれば、後で収束しても失敗が残ります。参照自体の通常判定 PASS も必須です（同ファイル:353）。

   case/48 の [README.md:43](/home/sano/work/forge-sern-design/case/48.flat_plate_cooled_m4/README.md:43) は、後続の標準 run も `NOT CONVERGED (stalled/plateau)` と記録しています。これを「収束場」として自動的に床参照にできません。指定された `run_0005_B_tw300_y3`、case/40 の `run_0045_node_yp1_outletfix_cont` は、この checkout にはありません。AWS 上の所在・入力・履歴は未確認です。

   **対案:** 実装前に起点の所在、snapshot、メッシュ、BC、実効設定、バイナリを固定し、基準の VERDICT を確認する。各方式への切り替え後は別の離散化区間として収束を判定し、`--from-floor` はその方式自身の PASS 場からの継続に使う。指定 step 数は初回予算とし、未達なら同一設定で延長する。case/48 のプラトーを許容するなら、既存ゲートとの関係を明示して事前に別の合否を定義する。

5. **Major — 「格子対で差が縮めば合格」は、欠陥と離散化差を識別できない。**

   **根拠:** [plan:93](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:93) では、格子対、細分方向、縮小率、誤差床が未指定です。例えば実装誤差が `Ch` なら、誤った実装でも差は縮みます。逆に壁方向だけを細分して流れ方向由来の差が残れば、正しい実装を停止します。case/48 には壁方向の系列と流れ方向の系列が別々にあります（[README.md:36](/home/sano/work/forge-sern-design/case/48.flat_plate_cooled_m4/README.md:36)、:39）。

   「既存ゲート」も部分的です。case/48 の手順書は `Cf` に加えて `2St/Cf`、エネルギー閉合などを扱いますが、S2 の合格列は `Cf` だけです。case/16 の「0.5%」「onset 1間隔」も正規化と間隔が未定義です。準定常ツールの既定値は drift 5%、oscillation 10%であり、表の0.1〜1%比較には粗すぎます（[check_quasisteady.py:543](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:543)）。

   **対案:** 既存の物理ゲートは独立した必須条件として維持する。格子対での縮小は診断に限定し、許容差超過を自動合格にしない。事前に格子ファイル、細分方向、比較位置・ノルム・正規化、必要な縮小率、反復床を固定する。準定常の drift・振幅も各比較許容より十分小さく設定する。

   今回再実行した `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/r1_series_1p6M.csv` の判定は **`DRIFTING`**、`dF1_inf` は末尾平均 **0.03341、drift 0.3%/tail**でした。これを後に LSQ が STEADY になったというだけで「GG 経路の性質」と一般化せず、そのケースでの観測に限定してください。

6. **Major — 既定値の調査は正しいが、`Y/ξ/モーメント` の検証を縮小する根拠にはならない。**

   **根拠:** `speciesFaceReconstruction=0` と [main.cpp:1453](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1453) の条件分岐は確認できました。したがって「既定では、SST 有効時の `k/ω` が対象」は正しいです。

   しかし、過去の opt-in 使用は [case/16 README.md:381](/home/sano/work/forge-sern-design/case/16.nozzle_wys/README.md:381) の非一様組成・トレーサ run や :384 の凝縮 S3 run に記録されています。ローカルの config 検索結果から使用実績なしとは判断できません。

   指定された [compare_s3_wall_pp0.py:18](/home/sano/work/forge-sern-design/case/16.nozzle_wys/compare_s3_wall_pp0.py:18) は保存済み壁圧・凝縮率 CSV の描画です。`ξ` の検証ではなく、勾配が実際の面再構成へ届いたかも判定しません。定常 run では、今回影響する `Pface` を読む dual-time FCT も作動しません。

   **対案:** 凝縮 case は維持し、非一様 `Y/ξ` の既存 `run_0476`／後続修正版の構成を追加する。`thermalMethod=2`、複数種、SLAU、SFR=2、対応する coupling を固定し、面値・更新後の有界性、組成和、補正収支を測る。FCT は既存の小規模 dual-time 試験を LSQ で再実行する。作用素を意図的にゼロ化・誤登録した版を落とせることも確認する。

7. **Major — 既定化の provenance は記録だけでは足りず、段履歴と設計 DB の混在防止が必要。**

   **根拠:** [main.cpp:90](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:90) 自体が、`RUN_PROVENANCE` は起動ごとに上書きされるため追記型の起動記録が必要だと説明しています。

   その先例にも注意が必要です。[stage_manifest.py:148](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:148) は同じ config hash に最後の起動値だけを対応させます。同じ YAML の旧既定0・新既定1を持つ二段をメモリ上で与えると、現行 `segments()` は **`[['old', 'new']]`** と誤連結しました。この仕組みをそのまま流用してはいけません。

   設計側の学習データ選別は現在 `flag_policy` を見ています（[driver_sern.py:234](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:234)）。スカラー勾配の記録だけ追加しても、GG/LSQ の評価を分離できません。

   **対案:** 各段を一意の起動記録へ結び付け、バイナリ識別子と実効 `scalarGradient` を保存する。旧 manifest の欠落は `gg` とし、新しい記録の欠落とは区別する。設計評価・台帳・学習選別にも実効方式を含める。旧省略／新省略／明示 gg／明示 lsq／cell 実効 gg の対応と、同一 YAML を異なるバイナリで再起動した場合を、既定化前の必須試験にする。

8. **Minor — 共有アクセサの初期化・寿命と ghost 配列の契約が未記載。**

   **根拠:** 現行 `cInt` は [calcGradient_d.cu:959](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:959) で初回 NS 勾配計算時に作られ、再構築条件は `nCells` だけです。汎用化後に独立ハーネスから呼ぶ場合や、同じ節点数の別メッシュを扱う場合の条件が曖昧です。また受動種 wrapper は現在 ghost を含めてゼロ初期化しています（[speciesTransport_d.cu:1268](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1268)）。

   **対案:** アクセサは読み取り専用とし、初期化済みメッシュとの対応を検査する。境界 incidence は係数0を掛ける前に明示 skip し、ghost 出力の初期化契約を維持する。単一メッシュ・単一プロセス専用なら、その制約を明記する。

9. **Minor — S3 の5%は、現在の記述では既定化の条件になっていない。**

   **根拠:** [plan:103](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:103) は5%超過でも「項目を立てる」だけなので、未最適化のまま既定化できます。測定環境、反復回数、初期化・出力時間の扱いも未指定です。

   **対案:** AWS native・同一GPU・同一ブロックサイズ・同一状態でウォームアップ後の測定窓を固定し、反復してばらつきを記録する。5%を既定化の保留条件にする。`REG` に加え spill と実測時間を残す。

**推奨は、共有差分形 LSQ を維持し、上記を plan に反映してから opt-in 実装へ進むことです。** 実装前の優先順は **SST の事実訂正 → 周期入力の同値性 → S0 の修正 → 起点・収束判定 → S2 の定量化とスカラー検証 → provenance/DB 分離**です。既定化は修正版の検証と result レビューを通した後に限定してください。

ファイルは変更していません。以上の提案は **plan 未反映**です。新規 CFD は実行せず、コード、保存済み判定記録、case/39 の準定常再判定と小規模な数値反例で確認しました。

指摘数: Critical 0 / Major 7 / Minor 2
