# codex 諮問 (diagnose): c48-line-implicit-reproducer

- **brief**: [`notes/reviews/briefs/2026-10-11-c48-line-implicit-reproducer.md`](../../notes/reviews/briefs/2026-10-11-c48-line-implicit-reproducer.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `8f957a85` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **case/48 に提案の陰解法設定一式を固定した float×2・FP64×2 の比較を行い、符号・時間的持続・再実行の変動を含む「有限期間の再現判定」を新規登録する。**
- **extra**: `case/48.flat_plate_cooled_m4/c48_prec.sh`, `case/48.flat_plate_cooled_m4/tools/c48_prec_judge.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| Major | 小さい case/48 に陰解法設定一式を移す案：**採用**。ただし「line implicit 単独の試験」という解釈は却下 | `lineDtDirectional` は擬似時間刻み、`implicitThermalJacobian: 5` はエネルギー行も変更する（[solver-settings.md:75](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:75)、[同:98](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:98)）。**設定一式を固定した精度 A/B**と位置付ける。 |
| Major | §6.21 の判定器をそのまま再現判定に使う：**却下** | [c48_prec_judge.py:95](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/tools/c48_prec_judge.py:95) は14量のどれかが非定常なら総合を判別不能にする。今回の目的である「有限期間に精度依存の偏りが育つか」と、定常解の比較を分ける。旧判定は残し、新しい再現判定を別に事前登録する。 |
| Major | 「θ か q_w のどれかに絶対差0.5%以上」で case/45 を再現：**要再検証** | case/45 の支持された主判定は**正方向の θ_r**であり、Q_w は副判定の支持条件を満たしていない（[plan §6.27](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:725)）。主判定は符号付きの θ とし、q_w だけの差や逆方向の差は別の精度感度として記録する。 |
| Major | 壁から線を作れるので介入は有効：**要再検証** | [mesh.cpp:1287](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.cpp:1287) は長さ3未満の鎖を不採用にする。線数・被覆率・長さだけでなく、評価断面の境界層内部が線に載ること、両精度の接続が同じことを確認する。方向別 dt の実際の伸びも必要。 |
| Major | 共通 Q32 入力なら起動後も同一状態：**要再検証** | [plan §6.31:801](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:801) に、共通入力でも初期化で ρE が精度別に変わった実績がある。入力7保存量の float32 表現可能性と、初期化後の保存量・評価量の差を別々に検査する。 |
| Minor | 既存抽出台本の流用：**条件付き採用** | [c48_prec.sh:31](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/c48_prec.sh:31) は `--coords-from` を付けない。共通 double 座標を明示し、節点対応検査を必須にする。現行判定器には `--arms` があるため、新しい4本を明示できる。 |

結論: **case/48 に提案の陰解法設定一式を固定した float×2・FP64×2 の比較を行い、符号・時間的持続・再実行の変動を含む「有限期間の再現判定」を新規登録する。**

第 1 仮説: line 系の結合と方向別擬似時間刻みを含む設定一式が、case/48 でも float の正方向の θ の偏りを増幅する十分条件になる。確度: **低（未検証）**。  
　根拠: `lineDtDirectional` は内部 line 面を CFL 評価から除き、上限50は平均流・SST双方の刻みに効く（`procedures/solver-settings.md:75–76`）。偏りの増幅経路としては候補になるが、case/48 での実効刻み・偏りの実測はない。  
　反証条件: 接続・設定の実効性を確認したうえで、下記の両窓・全 θ 断面が小差の条件を満たす。この場合に退けられるのは、**この case・入力・48,000 step で設定一式だけが十分という仮説**まで。

第 2 仮説: case/45 の収縮部・圧力勾配・軸対称・TP・格子条件との組合せが必要である。確度: **低、未確認**。現資料から軸対称と TP の優先順位は決められない。

判別 A/B: **両腕間で変えるのはビルド精度だけ**。提案の `lineImplicit: 1`、`lineDtDirectional: 1`、上限50、`implicitThermalJacobian: 5`、`cfl/cfl_pseudo: 4` は全4本に共通とする。

- 各48,000 step、2,000 stepごとに出力。判定窓は事前に W1＝24,000〜36,000、W2＝36,000〜48,000 と固定する。
- θ(0.3/0.6/0.9) を主判定とする。各 run の自分の起点からの変化を、共通の FP64 起点値で正規化する。各窓で、その変化の精度間平均差を符号付き e とし、各精度内の全窓点の平均からの最大偏差の和を E とする。E は観測範囲であり信頼区間ではない。
- **結果A：同じ θ 断面で両窓とも e − E ≥ 0.005。** → 「case/45 と同方向の0.5%以上の精度依存の偏りを、小さいケースの有限期間で再現」とする。軸対称・TPがこの現象に必須という説を退ける。真因が case/45 と同じとは確定しない。
- **結果B：全3断面・両窓とも |e|＋E ≤ 0.0005。** → この条件・期間では設定一式だけで θ の偏りを再現できず、第1仮説の限定された十分性を退ける。
- それ以外は判別不能。逆符号の大差や q_w だけの大差は別に記録し、AにもBにも読み替えない。

共通 double 座標を使い、`je` の履歴も保存する。主判定を満たす θ が `je` の切替えに依存する場合は、事前に定めた共通起点の `je` に固定した抽出でも判定が保たれることを要求する。入力が共通でも初期化後の差が大きい場合は、更新による増幅の判定を保留する。

全期間の残差・保存出力の有限性、`check_convergence`、上記区間の `check_quasisteady` の VERDICT は必須で記録する。非定常でも有限期間の差は診断できるが、**「定常差」「収束解の偏り」とは呼ばない**。旧§6.21の判定・閾値は変更しない。

やらない方がよいこと: 「再現条件に届かなかった」を「line 系は無関係」に読み替えること、結果を見て窓や閾値を変えること、軸対称・TPの順を原因の確度として決めること。結果Bなら、次候補は**case/45 の物理・BC・陰解法設定と近壁解像を保持した粗格子版**を優先する。これは既知の症状へ構成を近付ける探索順であり、原因の順位ではない。

呼び出し側の前提への異議: **「case/48 は偏りなし」は受け入れない。** 正式記録は「判別不能」で、限定して言えるのは指定窓の平均に%級の差を観測しなかったことだけ（[plan §6.24](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:683)）。また、両 case の違いは列挙された3群だけではなく、収縮部・圧力勾配・格子・壁解像も含む（同§6.21）。§6.29・§6.33・§6.35から SST・左辺・平均流の丸めを除外することもできない。

不足情報: 対象の case/48 の run、入力 HDF5、判定 JSON、指定バイナリはローカルに無く、実測の再検査はできなかった。起動後の初期状態、line の被覆・接続、実効 dt、バイナリの同一性は投入側で確認が必要。**ファイル変更・forge 実行は行っていない。本提案は plan 未反映であり、呼び出し側が登録してから投入する。**
