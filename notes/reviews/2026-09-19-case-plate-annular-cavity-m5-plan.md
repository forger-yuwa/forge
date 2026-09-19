# codex レビュー: case-plate-annular-cavity-m5 (plan)

- **plan**: [`plans/active/case-plate-annular-cavity-m5.md`](../../plans/active/case-plate-annular-cavity-m5.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `836d4d80` (feature/sern-design)
- **codex**: effort `high`, 6.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C2/M8/m1
- **focus**: §3.1 自由流条件の妥当性, §4.2 境界条件と領域寸法, §4.3 メッシュ (VL 13層10µm・すきま解像・節点数予算), §4.5 段階起動 S0-S6, §4.6-4.7 深いキャビティの収束性と定常RANSの限界, §6 検証ゲート (参照解が無い中で何をもって妥当と言うか)
- **extra**: `case/37.pintle_nozzle/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

現状のまま実装には進めません。3D restart が深さ座標を失う実装、成立しない熱収支ゲート、定常 RANS を時間平均場とみなす根拠の誤りがあります。  
ケース自体の意義と既存機能の利用方針は妥当ですが、計画を修正して再レビューすべきです。ファイル変更・新規計算投入は行っていません。

`plans/README.md` と関連する `accepted/` を確認した範囲では、本ケースと同じ課題を解決済みの計画はありません。`case/48` の利用は検証ケース選定として適切です。ただし、その実績の適用範囲を超えています。node、SLAU、低 Re SST、block-DPLUR の組合せは既存構造と整合し、周期・軸対称は本ケースには該当しません。

1. **Critical — Stage A→B の restart が深さ方向の場を破壊する。**

   **根拠:** [interp_field.py:34](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:34) は node 座標を `[:, :2]` に切り、cell 経路も同様です。そこで作った KD-tree を [同:159](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:159) で使用します。[plan:244](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:244) の cross-mesh restart は、キャビティ深さを表す z を無視します。

   メモリ内の最小再現では、同じ `(x,y)` の深さ −50 mm／0 mm に置いた供給点に対し、−49 mm／−1 mm の移植先が**両方とも同じ供給点**を選びました。同一メッシュだけを index コピーにしても解消しません。

   **対案:** 3 座標による移植を先行必須作業にする。壁・内部点の取り違え、移植距離、温度・圧力・保存量の再構成、新メッシュの `wall_dist` 保持を検証してから Stage B に進む。「原始変数補間」という記述も、実際の保存量転送と合わせる。

2. **Critical — エネルギー収支の検査対象が誤っている。**

   **根拠:** [plan:278](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:278) は「等温壁 4 面の総入熱」と「環状開口のエンタルピー流束」を比較します。しかし `cyl_top` はキャビティ内部の検査体積を囲まず、外部流に直接接する面です。この入熱を足した収支は閉じません。

   また低 Re 壁では、[viscousFlux_d.cu:521](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:521) が温度勾配から熱流束を計算します。実ファイル `case/48.flat_plate_cooled_m4/run_0005_B_tw300_y3/res_wall_4_48000.h5` の `VALUE/qwall` は **1001 点すべて 0** でした。`outputHDFflg: 1` だけでは必要な熱流束を取得できません。

   **対案:** キャビティ収支は `cav_outer + cav_floor + cyl_side` と環状開口で閉じ、`cyl_top` は別に報告する。開口では移流全エンタルピーに加え、分子・乱流熱伝導と粘性仕事を含める。熱流束抽出は低 Re の勾配再構成を明示し、node 温度ピンによる離散収支との関係、法線・符号・面積重み・半領域からの換算を先に検証する。

3. **Major — §4.7 の cavity 分類と、定常 RANS の正当化が誤っている。**

   **根拠:** [plan:213](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:213) の「20」は深さ／すきま幅であり、流れ方向開口長／深さではありません。対称面の各スリットでは後者は `2.5/50=0.05` です。さらに、せん断層が開口を跨ぐものは通常 **open cavity**、床へ再付着するものが **closed cavity** です。計画の説明は逆です。[NASA の cavity heating 研究](https://ntrs.nasa.gov/api/citations/20090007685/downloads/20090007685.pdf?attachment=true)

   定常擬似時間計算の収束は物理的振動の不在を証明せず、半割は反対称モードを禁止します。「振動すれば内部温度が高くなる」という偏りの方向も、本形状については未検証です。

   **対案:** 定常・半割 RANS を**条件付きの初期評価**に位置付ける。実際の時間平均場を完了条件に残すなら、物理時間による非定常確認と対称性制約の感度を、本番移行前の必須ゲートにする。定常計算が収束した場合も省略しない。

4. **Major — 面重心による境界分類では、指定の壁を分類できない。**

   **根拠:** [plan:125](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:125) は円筒面を重心の `r≈25`／`22.5` で分類します。しかし半円筒側面の面重心半径は `2R/π` で、**15.915 mm／14.324 mm** です。重心は円筒面上にはありません。流用元も [mesh_salome.py:124](/home/sano/work/forge/case/37.pintle_nozzle/cad/mesh_salome.py:124) で実際に面重心を使用しています。

   また `x=-200 mm` で CAD 面を分割しなければ、連続した平面を `runup` と `plate` に別タグ付けできません。総表面積一致では誤タグを検出できません。

   **対案:** CAD 生成時に境界の切替位置で面を分割する。円筒は曲面種別・軸・半径、平面は支持平面と範囲で分類し、グループ別面積、未分類ゼロ、重複ゼロ、法線方向まで検査する。

5. **Major — S0–S6 に起動の前提と段間の移行判定が不足している。**

   **根拠:** [plan:186](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:186) の全壁 `slip` からの起動には、変換時の壁指定がありません。[recommended-settings.md:71](/home/sano/work/forge/procedures/recommended-settings.md:71) は、SST 用の壁距離を作るため変換時の no-slip 壁指定を必須としています。非構造 float32 用の `space.pRef` も設定例から欠落しています。

   自由流一様＋壁から 300 µm のランプでは、すきま中央の初期速度はほぼ M5 のままです。これは「深部は準静止」という想定と整合せず、内外の静圧が同じでも大きな速度過渡が発生します。固定 step 数・NaN なしだけでは次段への移行を保証できません。

   **対案:** 最終壁タグで変換して壁距離を保持し、実行時に S0 の BC を切り替える。`pRef≈5475 Pa`、外部流と深部低速場を滑らかにつなぐ IC、各段の粘性・乱流・`nStepInner` を明示する。段間移行は全残差の下降、圧力・温度・速度範囲、キャビティ質量・エネルギーの変化を条件にし、固定 step 数は初回評価時点とする。

6. **Major — メッシュ品質 PASS と格子感度が、必要な解像度を保証していない。**

   **根拠:** 13 層・stretch 1.22 の総厚 **0.55746 mm** は正しいです。ただし、計画自身が予想する外部 BL 厚 4–7 mm の一部しか prism で覆いません。残りの BL、開口せん断層、prism→tet 遷移の法線解像度が未指定です。[plan:132](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:132)

   また [check_mesh_quality.py:70](/home/sano/work/forge/solver_density_cuda/tools/check_mesh_quality.py:70) の 3D 指標は辺長比と面内角です。単位辺長に対し厚さ `1e-6` の人工 tet を評価すると、**AR=1.414、skew=0.250、体積=1.667e-7** でした。体積が潰れた sliver を、このゲートだけでは排除できません。

   **対案:** 符号付き体積・scaled Jacobian／sliver 指標・dual CV の正体積と閉性を追加する。`y₁⁺` は壁ノードではなく第一内部点で評価し、10→20 µm の粗化だけでなく細化側を含める。開口・外部 BL・すきまコアも細化対象とする。節点数上限は暫定予算として、SST 本段の GPU・host メモリ双方で確定する。

7. **Major — 自由流・流入境界層・領域寸法を、根拠以上に確定している。**

   **根拠:** [plan:46](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:46) の高度 20 km は、M5 と壁温からは決まりません。CPG としての表の主要数値は概ね整合し、再計算した `Taw=1187.55 K`、指定 k/ω による `μt/μ≈10.16` も妥当です。しかしこれは、選んだ作動条件の算術が合うという意味です。

   200 mm 助走で得る BL と全乱流の仮定は、キャビティ流入状態を直接決めます。また [plan:112](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:112) のマッハ角は微小擾乱の角度であり、有限強度の圧縮波や亜音速 BL を介した影響まで排除しません。`cyl_top` の熱条件が侵入深さを変えないという断定にも根拠がありません。

   **対案:** 20 km・全乱流・助走長・`cyl_top=500 K` を解析の仮定として明記する。開口直前の `U,T,k,ω,δ*,θ` を評価し、流入 BL と上面・側面・出口位置の感度を必須化する。回復温度ゲートは、キャビティ・前縁・出口の影響を避けた指定区間に限定する。

8. **Major — TP 感度の物性定義と restart 手順が未成立。**

   **根拠:** [plan:226](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:226) の「1300 K で cp が約 +8 %」は小さすぎます。repo の [NASA-9 係数](/home/sano/work/forge/design/forge_design/gas/semiperfect.py:30) を用いた乾燥空気組成の計算では、`cp(298.15)=1004.73`、`cp(1300)=1188.31 J/(kg K)`、**+18.27 %、γ=1.3185** でした。同じ `T∞,U∞` に対する TP 全温は約 **1222.9 K** です。これはキャビティ温度誤差の予測ではなく、入口熱力学の計算結果です。

   `MIXDRY` は内蔵種ではなく、DB の指定が必要です。[speciesDB.cpp:209](/home/sano/work/forge/solver_density_cuda/input/speciesDB.cpp:209) は未定義種を拒否します。さらに CPG の `roe` を、298.15 K datum の TP に直接コピーできません。

   **対案:** 乾燥空気組成と専用 `speciesDBFile` を固定する。CPG→TP は `P,T,U,Y` を保持して新 EOS で密度・内部エネルギー・保存量を再構成する。輸送物性を揃えた TP 比較を Stage A に前倒しし、生産格子の感度評価に入る前に正本 EOS を確定する。

9. **Major — 準定常判定が、目的量の変化を見逃す。**

   **根拠:** [plan:274](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:274) の CSV 列には侵入深さがなく、`mdot_mouth` の定義もありません。符号付き正味流量は定常の底閉じキャビティではゼロに近づくため、交換量ではありません。`T_floor` が壁温なら、[nodeWallDirichlet_d.cu:79](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:79) により常に 500 K です。

   ツール既定値は drift 5 %・変動幅 10 %です。[check_quasisteady.py:487](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:487) の実関数に、500→510 K と単調上昇する 10 点系列を渡すと **`STEADY`** でした。500 K を分母にすると、侵入による小さな温度上昇の変化が隠れます。

   **対案:** 流入量・流出量・正味不釣合いを別列にし、交換率を定義する。侵入深さの温度閾値・周方向集約方法、壁から離れた温度測点を指定する。温度上昇 `T−Tw` と絶対 K 許容を併用し、反復誤差を格子感度 5 % より十分小さくする。`T≥500 K` は入口が 216.65 K の開いた圧縮性流れに対する一般的下限ではないため、無条件の棄却基準から外す。

10. **Major — 先行検証の到達点を過大に引用し、未収束時の停止条件がない。**

    **根拠:** 今回 `check_convergence.py` を実行した結果は以下でした。

    | `case/48.flat_plate_cooled_m4/` 配下 | VERDICT |
    |---|---|
    | `run_0004_A_ad_y3_cont/` | `NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)` |
    | `run_0005_B_tw300_y3/` | 同上 |
    | `run_0007_B_tw300_y6/` | 同上 |
    | `run_0008_B_tw300_y12/` | 同上 |
    | `run_0011_Bplain_tw300_y3/` | 同上 |

    例えば `run_0004` の最終 `rms_roK=5.30e-5`、`rms_roOmega=1.69e2` はまだ下降中です。[先行 README:33](/home/sano/work/forge/case/48.flat_plate_cooled_m4/README.md:33) も断熱温度の発達途中を記録し、SU2 比較には δ* 約 +4 %、断熱 Cf の VD-II 比 0.87–0.92 などが残っています。[plan:283](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:283) の一括した「±3 %/±1 %」は成立しません。

    **対案:** 先行実績を量・モデル・VERDICT ごとに限定して引用する。`PASS` と準定常性を別ゲートとして維持し、反復数を増やしても改善しない場合の中止・診断条件を決める。収束済み restart の `--from-floor` は、通常判定で PASS した参照がある場合だけ使用する。内部整合と格子感度から「実流れに対して検証済み」とは結論しない。

11. **Minor — 現在仕様への参照と物性定数が一致していない。**

    **根拠:** [plan:8](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:8) の `methods/boundary/index.md` は存在せず、現行文書は `methods/boundary.md` です。関連づけた [wall-function coverage plan:30](/home/sano/work/forge/plans/accepted/turbulence-node-wall-function-coverage.md:30) は 3D を対象外としており、3D 角線問題の直接の根拠ではありません。Sutherland 定数も計画の `273.15/110.4` に対し、[実装:59](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:59) は `273.0/111.0` です。

    **対案:** 実在する仕様文書・該当コード・run へ参照を直し、自由流表を実装の物性定数から生成する。

**推奨は、実装着手前に計画を改稿して再レビューすることです。** 優先順は、①3D restart と正しい検査体積・熱流束評価、②CAD 分類と起動条件、③物理モデルの適用限界と定量ゲート、④TP・格子・領域感度です。評価ツールを Stage A より前に準備し、Stage A でこれらの成立性を確認してから生産計算へ進む構成にしてください。

本回答はレビュー提案であり、指定どおり **plan 未反映**です。

指摘数: Critical 2 / Major 8 / Minor 1
