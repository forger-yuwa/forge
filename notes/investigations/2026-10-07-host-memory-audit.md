> 2026-10-07、plan `architecture-solver-host-memory.md` §5.1 #2 の監査 (implementer、読むだけ、HEAD 9c9f623c)。行番号は 9c9f623c 時点。

PLAN: plans/active/architecture-solver-host-memory.md §5.1 #2 (監査、読むだけ)
変更: なし (どのファイルも編集していない。commit もしていない)
ビルド: 未実施 (読むだけの項目なので)
確認: grep/目視による監査のみ。作業ツリーは /home/sano/work/forge-sern-design (HEAD 9c9f623c)。対象は forge_cppfiles (CMakeLists.txt:70-85) と cuda_forge/・probe/ のホスト側コード

## 0. 結論 (要点)

- GPU 経路でホストの `var.c` を**読み書きする**のは 5 か所だけ: `readValueHDF5`・dual-time の checkpoint 復元 (main.cpp:1148)・出力 (output.cpp)・環境変数で有効になる 2 つの診断 (FORGE_IMPLICIT_DIAG_CSV、FORGE_PIN_DIAG)・lineImplicit (main.cpp:3503)。それ以外 (probe、壁出力、CHT、残差・NaN の集計、初期化で使う入口分布・壁距離・種/凝縮/受動種) はデバイスか自前の局所バッファで済んでいて、`var.c` を経由しない。
- GPU 経路でホストの `var.p` に触るのは `allocVariables` の resize と 2 つの診断だけ。そのうち **R2 で挙動が変わる箇所が 1 つある**: main.cpp:2303 の `pdeSize(s.var.p, k, 0)` (FORGE_DIAG_PSI_DUALEVAL)。ホストの p が空になると、面の配列が黙って退避・復元の対象から外れる。R2 では fallback を `msh.nPlanes` に直す必要がある。
- GPU 経路にはホストの `c[...]` の operator[] (暗黙生成) が無い。operator[] の使用は CPU 専用の経路と変換器だけ。`c_d[...]` の operator[] は多いが、触るのはデバイス側のマップだけ。

## 1. ホスト `var.c` 参照の表 (GPU 経路で到達するもの)

凡例: W→H2D = ホストに書いてから H2D する経路 (書込みが転送より前)。保護「at」= `map::at` でキーの有無だけを見ている (長さは見ない)。

| # | 名前 (式・リスト) | 有効になる条件 | 最初の利用 (ファイル:行) | 読み/書き | 現状の保護 | H に入れる理由 |
|---|---|---|---|---|---|---|
| A | ro, roUx, roUy, roUz, roe, wall_dist, roK, roOmega | 常に | variables.cpp:787-810 で書く、818-820 で wall_dist を読み書き、827 で H2D | W→H2D | at | 初期場の読込 |
| B | roY{s}, Y{s} (+ ro を読む) | nSpecies≥2 | variables.cpp:839-865 で書く (853・859 で ro を読む)、870 で H2D | W→H2D | at | 初期場の読込 |
| C | roGamma, roReth | transition: lm2009 かつ入力に両方ある | variables.cpp:880-883 | W→H2D | at | 初期場の読込 (入力に無いときは書かない。ただし出力の base には常に入る) |
| D | roXi, Xi (+ ro) | physProp.tracer | variables.cpp:891-913 | W→H2D | at | 初期場の読込 |
| E | ro{g,Q2,Q1,Q0}_s と g_s, Q2_s, Q1_s, Q0_s (+ ro) | condensation (nCondSpecies≥1) | variables.cpp:923-939 | W→H2D | at | 初期場の読込 |
| F | roN, roUxN, roUyN, roUzN, roeN, roKN, roOmegaN + roY{s}P + (passiveScalarScheme 1 のとき) roXiP, <cons>P | unsteady 1 かつ dualTime 1 | 復元: main.cpp:1148-1149 で書き、1153 で H2D。書出し: output.cpp:227 で D2H、230 で読む | W→H2D と R | at だけ (**M1 の経路**) | checkpoint の復元と書出し |
| G | level 0 の出力: ro..roe, roK, roOmega, roY{s}, <cons>, roXi, roGamma, roReth | 常に | output.cpp:122 で D2H、207 で読む (最初の到達は main.cpp:3526 の writeInitialOutputs) | R | **無し** (207 は長さを見ずに `begin()+nCells` で読む) | 出力 |
| H | level 1 の追加分: P, T, Ux, Uy, Uz, k, omega, sonic, vis_lam, vis_turb, wall_dist, Y{s}, <prim>, Xi, gammaTr, reTheta, gammaEff。h0 用の **Ht, k** | output.level≥1 (既定は 1) | output.cpp:83-88、117 (`c.count` なので常に真)、121-122 で D2H、263 | R | 無し / at | 出力と h0 の依存名 |
| I | level 2: env による除去の後の output_cellValNames 全部。登録で増える分も含む (roYraw{s} [FORGE_SPECIES_RAW_DIAG, variables.cpp:86]、chemQdot/chemTau [chemistry]、lm* と src_jac_gamma/reth [transition かつ level 2]、d{Xi,prim}d*・limiter_*・passive*Corr_* [受動種]、cond*_s [凝縮]) と env で残る分 (wi_* [FORGE_WI_FORCE_DIAG]、wf_sprod [FORGE_WF_OMEGA_SOURCE]、wf_g [FORGE_WF_CLOSURE_DIAG]、rep_* [FORGE_WF_REP_DIAG]、omg_* [FORGE_OMEGA_BUDGET]。除去は variables.cpp:301-356) | level 2 | output.cpp:75-76 | R | 無し | 出力 |
| J | extraFields (登録されている任意の名前: wall_y_eff, dY{s}d*, res_*, dq_*, roN ほか) | output.extraFields | output.cpp:91、101-107 (受け付ける条件は `c.count \|\| c_d.count`) | R | count | 出力 |
| K | res_ro, res_roUx, res_roUy, res_roUz, res_roe, dq_block_old_0..4 | FORGE_OUT_RESIDUALS≠0 (**level 2 のときだけ出力される**。危険 8) | 登録: main.cpp:3417-3421 (**確保の後**) | R | 無し | 出力 |
| L | res_{ro,roUx,roUy,roUz,roe}_m, dq_{ro,roUx,roUy,roUz,roe}_new | FORGE_RESID_SNAP が設定されている (値は問わない) かつ level 2 | 登録: main.cpp:3426-3429 (確保の後)。デバイス側の退避は 2699-2707 | R | 無し | 出力 |
| M | ro, Ux, Uy, Uz, sonic, vis_turb, dt_local | FORGE_IMPLICIT_DIAG_CSV が空でない (main.cpp:721) かつ timeIntegration 11 (2123) | main.cpp:847 で D2H、850-856 で読む | R | at | 診断の CSV |
| N | scalarDirichletPin, res_roY{k}, res_roXi, res_<cons> / Y{k}, Xi, <prim> | FORGE_PIN_DIAG≠0 かつ node かつ dual-time 経路 (2677, 2756) | main.cpp:2210-2216 / 2230-2238 | R | at | 診断の標準出力 |
| O | ccx, ccy, ccz | lineImplicit 1 | main.cpp:3503 (`.data()` を渡す)。node ではノード座標を使うので読まない (mesh.cpp:1076-1087) | R (cell で nodes<nCells のときだけ) | at | 現状維持 (GPU 経路では**未充填の 0**。危険 5) |
| P | G〜L と同じ集合 (res_nan_ ダンプ) | detectNaN 1 かつ非有限値 | main.cpp:2855 → output.cpp:330 | R | 出力と同じ | 出力と同じ関数を通る |
| — | var.c 全体の capacity | FORGE_MEMLOG | main.cpp:143 | capacity だけ | 0 長でも無害 | 不要 (縮めた後の値が出るだけ) |
| — | c は長さだけ参照 (fallback nCells_all)、**p も長さだけ参照 (fallback 0)** | FORGE_DIAG_PSI_DUALEVAL | main.cpp:2293-2303 | 長さだけ | — | c: 不要 (fallback で同じ値)。**p: R2 で fallback を nPlanes に** |
| — | T, thermCond, cp, vis_turb, ro, roe | interfaceDiag 1 / conjugate | conjugateWall.cpp:59-72。gpu では c_d から局所 vector へ写す (ホストは c_d が無いときだけ、長さ検査つき) | — | 長さ検査あり | 不要 |

**H に要らない (ホストの c を経由しない) ことを確かめた経路**:
- probe: point_probes.cu:247-260 がデバイス上で集めて自前の var_d へ写す (T, P, Ux, Uy, Uz)。
- 壁出力: `bc.bvar` だけ (output.cpp:426-461)。
- 残差・NaN の集計: `makeDeviceResidualReducer`、residualMonitor_d.cu:42/65、main.cpp:2819-2831。
- 入口分布・壁温分布: bvar を書いて H2D (boundaryCond.cpp:404-414)。
- 共役伝熱の初期化: msh と bvar だけ。
- 種/凝縮/受動種の init_d と収支ログ: デバイス。
- FCT の G/H (`passiveFctHistoryTo/FromHost`): 自前の配列。
- 物性の probe (FORGE_TRANSPORT_PROBE、main.cpp:1304-1409) と TP 診断 3 種 (main.cpp:2910-3377): c_d・p_d だけ。
- 速度揺らぎ (fluct): デバイス。

**env の全一覧** (getenv を全部拾った)。名前の集合を変えるのは次の 10 個だけ:
- 登録・除去を変える: FORGE_WI_FORCE_DIAG, FORGE_WF_OMEGA_SOURCE, FORGE_WF_CLOSURE_DIAG, FORGE_WF_REP_DIAG, FORGE_OMEGA_BUDGET, FORGE_SPECIES_RAW_DIAG
- 出力名を追加する: FORGE_OUT_RESIDUALS, FORGE_RESID_SNAP
- ホストで読む: FORGE_IMPLICIT_DIAG_CSV, FORGE_PIN_DIAG

残りはデバイスか局所バッファしか使わない: FORGE_DUMP_{PREGATHER, SCALARGRAD, MASSFLUX, LEDGER*, FARFIELD*}, FORGE_LINE_*, FORGE_CONTACT_*, FORGE_DIAG_{TP_*, PSI_*, SU2*, FACE_VEL_CELL}, FORGE_TPD3_*, FORGE_FREEZE_*, FORGE_TAW_WALL_MUT_SCALE, FORGE_NUTFLOOR_*, FORGE_DT_OUTLET_*, FORGE_AXIS_DIAG_ALPHA, FORGE_FACE_THERMOY, FORGE_VISC_WALL_DIAG, FORGE_WF_{INV_LAW, CLOSURE_LAW, LIMITER_BYPASS}, FORGE_TRANSPORT_*, FORGE_KERNEL_SYNC, FORGE_CUDA_BLOCKSIZE*, FORGE_PROFILE*, FORGE_MEMLOG, FORGE_ALLOW_UNVERIFIED_SPECIES。
- FORGE_WALLDIST_BRUTE は変換器専用。

### CPU 専用・到達しない参照 (除外とその根拠)

| 参照 | 除外の根拠 |
|---|---|
| update.cpp:14-204 (`v.c[...]`) | 各関数の先頭 9/76/118/153/189 行で gpu==1 なら *_d_wrapper を呼んで return |
| dependentVariables.cpp:15-32 | :9 で gpu return |
| gradient.cpp:35-61 | :9 で gpu return、しかも gradientGauss の呼び出し元が無い |
| variables.cpp:505-510 (setStructuralVariables の CPU 部) | :484 で gpu==1 なら _d に分岐して return |
| main.cpp:790-800, 827-835, 2837-2844 | gpu 分岐の else 側 (785/822/2814)。827 は [[maybe_unused]] |
| setStructualVariables.cpp:26-30 | 呼び出し元が無い |
| convectiveFlux.cpp | 全体がコメント |
| input/setInitial.hpp, input/calcWallDistance_kdtree.cpp:176-187 | 変換器 (convertGmshToForge.cpp:75) からしか呼ばれない。main.cpp:27 は include だけ |
| mesh/gmshReader.hpp:2442・2659-2748 | 変換器専用。main.cpp:29 は include だけ |
| copyVariables_cell_plane_{H2D,D2H}_all, copyVariables_plane_{H2D,D2H} | 呼び出し元が無い |

補足: gpu: 0 は applyBconds が起動を拒否する (boundaryCond.cpp:478) ので、CPU 経路は実質動かない。

## 2. 初期化の順序 (現状)

1. main.cpp:1594 cfg.read (outputLevel・extraFields・unsteady などが確定する。以後これらを変更する箇所は無い。grep で確認)
2. 1600-1624 種 DB・熱力学・化学
3. 1630 readMesh
4. 1638 initMatrix
5. 1642 readBcondConfig (bvar)
6. 1647-1688 cfg の解決 (condLimiterMode を書き換える、二相の判定)
7. 1692-1698 入口分布・壁温分布・CHT (bvar と msh だけ)
8. 1702 registerSpecies (中の FORGE_SPECIES_RAW_DIAG は variables.cpp:86)
9. 1705-1706 registerCondensation / registerTwoPhaseVaporResidual
10. 1709 registerTracer
11. 1712-1713 registerTransition (診断場は outputLevel≥2 のとき)
12. **1715 allocVariables**
    - 内部 301-356 で env の除去 (output_cellValNames からも外す)
    - 357-385 でホストを resize (0 で埋める)、cudaMalloc、memset、sstF1=1
    - 387-396 で面配列
13. 1719-1724 種・凝縮・受動種の init_d (デバイスのポインタ表)
14. **1732 readValueHDF5**: ホストに書く → H2D (表の A〜E)
15. 1737-1768 メッシュ写像の H2D と周期の対応
16. 1771 setStructuralVariables_d (局所 malloc から H2D。ホストの c には書かない)
17. 1776-1817 デバイスの初期化
18. 1822-1825 FP64 アキュムレータ
19. **1829 initDualTimeHistory**: ホストに書く → H2D (表の F)
20. 1833 pprobes.init (probe.yaml を読む)
21. 1843-1894 リミッタの基準値 (局所 D2H)

main() 側:

22. 3402 ImplicitDiagLogger の ctor (FORGE_IMPLICIT_DIAG_CSV を読む。initializeSimulation より前)
23. 3404 initializeSimulation
24. 3406/3409 物性 probe で早期終了
25. **3417-3431 FORGE_OUT_RESIDUALS / FORGE_RESID_SNAP の出力名を登録 (確保の後)**
26. 3434-3494 qAcc の検査
27. **3503 lineImplicit** (ホストの ccx)
28. 3509 残差ロガー
29. 3513-3522 TP 診断で早期終了
30. **3526 writeInitialOutputs** (D2H してホストで読む)
31. ループ: 出力・NaN ダンプ・implicit-diag・pin-diag・PSI_DUALEVAL

### §4 の順序へ組み替えるときに動かす処理

- (a) main.cpp:3417-3431 の登録を 1713 の直後 (確保の前) へ移す。
- (b) variables.cpp:301-356 の env 除去を allocVariables から切り出し、H を作る前に実行する。
- (c) H = 共通関数 (§3) を 1714 で作る。必要な情報 (cfg・登録の状態・getenv) はこの時点で全部そろっている。
- readValueHDF5・initDualTimeHistory・lineImplicit は動かさなくてよい (H に名前が入っていればよい)。

### 動かすと挙動が変わりうるもの (順序依存)

- (i) (a) は output_cellValNames に push_back する。除去対象 (wi_* など) と名前が重ならないので、**並び順は今と同じになる**。変わるのは標準出力の行の位置だけ。
- (ii) FORGE_RESID_SNAP の登録条件 (`e!=nullptr`, main.cpp:3426) と使用条件 (`atoi(e)>=0`, 2699) は食い違っている。登録条件は今のまま写すこと。
- (iii) effectiveOutputNames の `static bool warned` (output.cpp:108)。初期化時に呼ぶと、警告の出る位置が前に動く。
- (iv) extraFields の受付条件 `c.count||c_d.count` (output.cpp:104)。c_d は実行中に `c_d["..."]` で新しいキーが増えうるので、H を作る時点と出力の時点で結果が食い違う可能性がある。共通関数では「登録済み (cellValNames)」で判定することを勧める (c と c_d のキーは登録の時点で一致しているので、通常の入力では等価)。
- (v) **長さ依存**: H2D は nCells_all 全体を写すので、ゴースト部にはホストの 0 が書かれる (A〜E と F の 1153)。H の配列は **nCells_all 長で 0 初期化を保つこと**。nCells 長に縮めるとビット一致しなくなる。
  - 推測: F のゴーストを 0 で上書きしていることが、連続実行との差 (case/09 run_0164 の roUy 7e-11) に効いているかは未確認。
- (vi) variables.hpp に長さのメンバを足すと構造体レイアウトが変わる。**クリーンビルドが要る**。

## 3. 出力・checkpoint の依存名を作る共通関数 (案)

```cpp
// output/outputFieldNames.hpp (新規)。出力側・確保側・復元側で共用する
struct OutputFieldPlan {
  std::list<std::string> solution;   // /VALUE と XDMF の順 (現 effectiveOutputNames。重複も今のまま残す)
  bool h0 = false; std::list<std::string> h0Deps;   // level>=1 なら {"Ht","k"}
  std::list<std::string> checkpoint;                // dual-time の履歴 (/CHECKPOINT)
};
OutputFieldPlan outputFieldPlan(const solverConfig&, const variables&);         // 入力: cfg・登録済みの変数 (env 登録は済ませておく)
std::list<std::string> dualTimeHistoryNames(const solverConfig&, const variables&);
std::list<std::string> initialValueNames(const variables&);                     // readValueHDF5 が書く A〜E
std::set<std::string>  hostCellSet(const solverConfig&, const variables&);      // H = 上の和 ∪ M (env と tI11) ∪ N (env と node と dual-time) ∪ O (lineImplicit)
void registerOutputDiagnostics(variables&);   // 現 main.cpp:3417-3431 をそのまま移す (getenv の条件も同一)
void applyEnvGatedRemovals(variables&);       // 現 variables.cpp:301-356 をそのまま移す
```

置き換える箇所:
- output.cpp:68-112 (→ solution)
- output.cpp:117-122 (→ h0 / h0Deps)
- output.cpp:217-226 (→ checkpoint)
- main.cpp:1105-1110 (→ dualTimeHistoryNames)
- variables.cpp:826 と種/遷移/トレーサ/凝縮の名前 (→ initialValueNames)
- main.cpp:842-844, 2206-2209, 2225-2228, 3503 の名前リスト (→ hostCellSet の分岐)

アクセサ: `variables::hostCell(name)` を設ける。キーが無い、または長さが hostCellLen_ (allocVariables・allocVariablesConverter で nCells_all を設定) と違えば、名前と長さを出して exit する。これを次で通す:
- copyVariables_cell_{H2D,D2H}
- output.cpp:207/230/263
- variables.cpp:787-931
- main.cpp:850-856/1148/2211-2238/3503

output.cpp:196 の `for (auto& v : var.c)` は h5 の書込み順を変えないため残し、名前で絞った後に hostCell を通す。

## 4. 試験の割り当て (§6 構成表の元)

入力の所在: 下の run はすべて **/home/sano/work/forge/case/...** (別セッションのワークツリー) にある。そこには書き込まず、入力を自分の新しい run_* へ複製して回す前提。sern-design 側の git status で変更が出ている case/09・16・46・48・53・54・58 も、入力の複製だけで使う。

| 表の行 | §6 の構成 | 入力の候補 (case / run / config のキー) | 確認する成果物と判定 |
|---|---|---|---|
| A, G, H, probe, 壁出力 | 2D node (標準) | case/36 `run_sym_H_2up_node` (node, SST, tI 11, 既定の level 1、detectNaN 1)。res_80000 から restart_field.py で再開。probe.yaml に 2〜3 点を足す | res_*.h5 のデータセット集合・shape・dtype・属性 (h0_includes_k を含む)、残差 CSV、壁 h5、point_probe_*.out |
| M | 同上 + env | case/36 に FORGE_IMPLICIT_DIAG_CSV=diag.csv。cell 版は case/20 `001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke` | diag.csv がビット一致 |
| R2 の pdeSize | 同上 + env | case/36 に FORGE_DIAG_PSI_DUALEVAL="2,4" | psi_dualeval.csv と log `snapshot: X of Y device arrays` が一致 (X が減ったら fallback 漏れ) |
| F, I, J, K, L, B, D, N (経路だけ) | dual-time の checkpoint | case/09 `run_0160_passiveG_fct_ckpt100` / `run_0164_…restart100_fixed` / `run_0162_…cont200` (node 周期、species [N2,O2]、tracer exhaust、passiveScalarScheme 1、FCT、BDF2、**level 2**、extraFields [volume, limiter_Xi])。100 step、そこから 100 step の再開、連続 200 step を変更前後で回す。env の変種: FORGE_OUT_RESIDUALS=1 と FORGE_RESID_SNAP=0、FORGE_SPECIES_RAW_DIAG=1、FORGE_PIN_DIAG=1 | /CHECKPOINT の 14 データセットと属性 (layout・nHistoryValid・totalTime)、log `history restored … 10 levels` と `(G/H/mEff) restored`、res_*_m・dq_*_new・roYraw0/1 があること、分割と連続の差 |
| B, D, E, F (<cons>P), N (入口あり) | 軸対称・多成分・凝縮 | case/44: 定常生産の最新 active run と、dual-time 版 `run_0376_passiveG_s3venkat_relax07_double_bdf2_dt8e-6_nsub40` 型 (axisym node、species [MIXDRY,H2O]、condensation 1、scheme 1、level 2、inletProfile) を 100+100 に分割。**要確認**: 現行バイナリが ic.h5 を受け付けるか (種の記録の照合、main.cpp:1202-1248) | 種とモーメントのデータセット、/CHECKPOINT/<cons>P、pin-diag の行 |
| CHT (pullField, interfaceDiag) | 共役伝熱 | case/52 `run_0007_fxhalf` (node, interfaceDiag 1, conjugate) | res_<wall>_*.h5 (iface* を含む) と CHT 状態の CSV がビット一致 |
| cell | cell モード | case/36 `run_sym_H_2up_cell` (cell, SST、上の node と同じ物理)、case/20 `001.test/run_slau` (cell RK3、追跡対象)、case/13 `run_slau` (cell)。cell の dual-time: case/20 `001.test/run_case04_les_unsteady_dualtime` (cell LES dual-time、res_wall あり) | 出力・残差 (cell は CONNE を cells から組む) |
| C, I (lm*) | 遷移モデル | case/57 `run_0014_t3a_lm_unitcheck` (node, LM2009, level 2, run_0011 の場から。roGamma を読む経路)。入力に roGamma が無い経路は SST だけの場 (例 `run_0012_t3b_sst` の最終場) から LM を始める | lm* のデータセット、log `[variables] transition … read from input / will be initialised` |
| O, J (level 1 の extra) | line-implicit と extraFields | case/56 `run_0019_lineimplicit` (node, lineImplicit 1, extraFields [thermCond, vis_lam, vis_turb])。`run_0027_s6_f32_b` (extraFields に res_roUz など output_cellValNames に無い名前)。case/48 `run_0903_absorb2` ([roN.., dq_block_old_*]) | 出力の集合、ライン構築の log |
| I (env で登録が変わる) | その他の環境変数診断 | case/26 `run_0086_optin_base` (extraFields [wf_pk, wi_ftan, wi_fnrm, wf_g, vis_turb])。env 無し (名前が消えて警告で無視される) と、FORGE_WI_FORCE_DIAG・WF_CLOSURE_DIAG・OMEGA_BUDGET・WF_REP_DIAG=1 の両方 | データセット集合、警告行 |
| メモリ・本番出力 | SERN 3D | AWS g3 (run_1068 の最終場)・g4。output.level と extraFields は AWS の solverConfig.yaml で確認する (ローカルには無い) | §6 のメモリ・残差・出力 |
| P | NaN ダンプ | 2026-10-07 の s050 縮小格子は step 9 で res_nan_9.h5 を出す。ただし入力は別セッションの scratch にあり、消えうる | res_nan_*.h5 のデータセット集合 |
| ガード (負例) | 負例 2 つ | (1) case/36 で H から "P" を外す → main.cpp:3526 の D2H より前で名前つきで停止すること。(2) case/09 の再開で H から "roN" を外す → main.cpp:1148 の書込みより前で停止すること | 停止メッセージ |
| (共通) | ビルド | convertGmshToForge は variables.cpp を共有し、copyVariables_cell_H2D を使う (setInitial.hpp:454, calcWallDistance_kdtree.cpp:187) | 変換結果をデータセット単位で比較 (h5 はバイト単位では再現しない) |

## 5. 危険

1. **output.cpp:196-207**: `var.c` を回し、出力名に一致した配列から `begin()+nCells` を長さ検査なしで写している。H に入れ忘れると、その前の D2H (0 長なので何もしない) を通った後で範囲外読みになる。黙ってゴミを出すか落ちる。
2. **variables.cpp:437-472**: copyVariables_cell_H2D/D2H の要素数はホスト配列の長さ。0 長だと黙って何もしない (呼び出し元の無い _all 版も同じ)。
3. **map::at しか無い書込み経路**: variables.cpp:787-939 (readValueHDF5)、**main.cpp:1148-1149 (checkpoint の復元)**。0 長の配列に `v[i]` で書くとヒープを壊す。読みの経路: main.cpp:850-856, 2211-2238、output.cpp:230/263。
4. **main.cpp:2293-2303**: pdeSize はホストの長さをデバイス配列の長さとして使っている。c は fallback が nCells_all なので R3 では変わらない。p は fallback が 0 なので、**R2 で面配列が黙って退避から外れる**。R2 の修正に含める必要がある。
5. **main.cpp:3503**: ホストの ccx は GPU 経路では一度も書かれない (setStructuralVariables_d は局所 malloc から H2D するだけ、variables.cpp:579-713)。node はノード座標を使う。cell で nodes.size()<nCells のときは 0 座標を読み、それ以外はセル番号でノード座標を引いている。既存の潜在不具合 (範囲外。触らない)。R3 では lineImplicit のとき H に入れて現状を保つ。
6. **全変数を回すコード**:
   - output.cpp:196 (上の 1)
   - main.cpp:143-144 (capacity だけなので無害)
   - main.cpp:2302-2304/2326/2406 (c_d を回す。2302・2303 は上の 4)
   - variables.cpp:357/387 (確保)
   - convertGmshToForge.cpp:79 (変換器)
7. ゴースト部の 0 の H2D (§2 の (v))。
8. **FORGE_OUT_RESIDUALS / FORGE_RESID_SNAP は level<2 では効かない**。output.cpp:75-96 の level 0/1 の base に入らないため。printf は「出力に追加した」と言う。さらに level 2 では、res_ro/res_roUx/res_roUy/res_roe が output_cellValNames (variables.hpp:278) と重複するので、XDMF の Attribute が二重になり (output.cpp:296-298)、D2H も二重になる。既存の動作なので、共通関数でも重複を含めて今のまま再現する。修正は範囲外。
9. writeH0 (output.cpp:117) は `c.count` で判定している。R3 でもキーは残るので常に真で、ガードにならない。level≥1 では Ht と k を必ず H に入れる。
10. variables.cpp:381-384 の sstF1=1 の H2D がループの中にあり、sstF1 より後ろの変数の数だけ繰り返される。デバイスの最終状態は同じ。allocVariables を組み替えるときに、ついでに直さないこと。
11. readValueHDF5 はファイル側のデータセット長を検査していない (variables.cpp:796-810)。既存。範囲外。

## 未了・保留

- §6 の構成表への組み込みと、plan への記入は親に任せる (plan は編集していない)。
- case/44 の dual-time 入力が現行バイナリで起動するかは未確認。SERN の output 設定は AWS でしか確認できない。
- 方針との不整合は無い。補足が 2 点ある:
  - §4.2 の「ホストの p を読むのは CPU 経路だけ」は**ほぼ正しい**が、例外として main.cpp:2303 (診断。長さだけを参照) がある。
  - §4.3 の「H には probe が読む名前を入れる」は、**実装上は不要** (probe はデバイスで集める)。

## 振る舞いの変化

無し (読むだけ)。
