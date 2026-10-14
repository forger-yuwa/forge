# codex 諮問 (diagnose): float-v4-and-hoop-review

- **brief**: [`notes/reviews/briefs/2026-10-10-float-v4-and-hoop-review.md`](../../notes/reviews/briefs/2026-10-10-float-v4-and-hoop-review.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-10
- **commit**: `2ee223fd` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.7 min, rc=0
- **結論**: **同じ float 到達場から「float 継続／FP64 切替」の対照付き A/B を一つ行い、case/45 の設計計算は当面 FP64 を維持してください。**
- **extra**: `plans/active/axisymmetric-freestream-hoop-gauge.md`, `notes/reviews/2026-10-10-axisymmetric-freestream-hoop-gauge-plan.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。対象 HEAD は `2ee223fd`。`run_0451_v4_f32`・`run_0452_v4_fp64` と判定成果物はこの作業ツリーに存在しないため、以下の実測値と VERDICT は **plan の報告値であり、独立検証済みではありません**。

| 指摘 | 重大度 | 採否 | 根拠・必要な対応 |
|---|---|---|---|
| hoop M1：全域 FAIL の原因確認 | Major | **採用。ただし誤差の分け方を修正** | [plan:86](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:86) の FAIL を維持する。入力の binary64 座標を厳密に固定し、端点生成・区間積・集約・`A_planar` の誤差を分離する。「座標の桁落ち」と先に命名して閾値を緩めない。 |
| hoop M2：拘束と独立な閉性検査 | Major | **採用** | [nodeWallDirichlet_d.cu:56](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:56) で壁の運動量残差をゼロにするため、出力残差だけでは壁 CV を検査できない。最終デバイス面ベクトルの閉性を全 CV で測り、流れの残差とは分ける。 |
| hoop M3：float32 の非劣化 | Major | **採用。定量条件の追記が必要** | 面を float32 に丸めた長方形の合成例で、閉性欠損/A = **0.4551915** を再現した。軸面の float32 二乗和によるノルムの **0 化**も再現した。double 生成だけでは保証できない。 |
| hoop M4：非静止場の採用ゲート | Major | **採用。「非物理値」も明記** | [plan:115](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:115) の「合否に使わない」は既定化の条件として不足する。同一バイナリ・格子・保存量で比較し、非有限・非物理・発散・上限までの未到達を不合格とする。量の差による既定化保留条件も事前登録する。 |
| hoop m1：適用除外・境界形状の検査 | Minor | **採用** | [dual_segment_closure.py:36](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/dual_segment_closure.py:36) は quad 専用。tri・曲がった境界・複数マーカ角、およびデータセットが存在しても無効になる分岐を追加する。 |
| hoop m2：現在仕様の説明更新 | Minor | **採用** | [methods:399](/home/sano/work/forge-integ-1005/methods/axisymmetric/implementation.md:399) の float32 起因の説明と、今回の FP64 にも残る集約誤差を分ける。 |
| float：相互 restart だけで二状態を判定する | Major | **却下。対照付きに変更** | [float plan §6.16](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:553) は両腕とも未収束。起点も精度も違う二本では、restart 後の自然な過渡と精度変更の効果を分離できない。同一起点の float 継続を対照にする。 |
| float：`S_lostfrac` から commit を主因とする | Major | **要再検証** | [commitLossDiag.cu:67](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/commitLossDiag.cu:67) は Σ\|dq−実更新\|/Σ\|dq\|。切り上げも数えるため「消えた割合」ではない。領域別・符号付き・時系列で調べる。 |
| float：V4 の比較成立条件の検証 | Major | **要再検証** | [v4_judge.py:47](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v4_judge.py:47) は残差の有限性を末尾窓だけ検査し、[同:86](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v4_judge.py:86) の収束ツール結果は表示だけで総合判定に反映しない。全期間の異常・ツール失敗・系列欠損を判定不能へ接続する。 |

M1 の高精度再計算では、**同じ入力座標から一貫した双対多角形を作った参照**と、**現行コードが丸めた端点を使う参照**を分けてください。内部は `G−M`、境界は `0.5·(B−A)` なので、端点の共有が浮動小数点上で一致しない可能性があります（[gmshReader.hpp:1492](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:1492)、[同:1656](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:1656)）。`A_planar` も独立に照合し、**符号付き欠損の差を各段階に帰属**させます。高精度でも閉じなければ、閾値改訂ではなく接続・向き・所属の修正です。

M3 は、軸・j 2〜8・収縮部内部・壁際・壁・入口角を固定して、E_x/E_y の最大と RMS、静止場残差、10 step 後の速度それぞれに **新 ≤ 1.1×旧＋10×再実行差**を事前登録する案を推奨します。これは非劣化条件であり、設計精度の保証ではありません。幾何検査は最終デバイス値を double 以上で集計し、流束の float 加算誤差とは分けます。軸床は double の幾何で分類し、床適用後のベクトルから double でノルムを計算して最後に丸め、`ss>0`・有限性・法線長を確認してください。

結論: **同じ float 到達場から「float 継続／FP64 切替」の対照付き A/B を一つ行い、case/45 の設計計算は当面 FP64 を維持してください。**

第 1 仮説: **float に依存する作用素・更新の差が、収縮部の遅い過渡を継続的に偏らせている。commit 単独主因とは絞れない。**　確度: **中**

- 根拠: V3 は同じ Q32・commit 前でも残差差を報告しており、ω の相対差は 0.37。V4 の報告では θ_r の差が 2.5万〜6万 step に成長しています。一方、両腕の VERDICT は **`NOT CONVERGED`／`NOT ALL STEADY`** です。したがって「二つの定常解」より、精度に依存する未収束軌道として扱うのが妥当です。
- 反証条件: 同一起点で精度を切り替えても、下記の期間・尺度で float 継続との差が十分小さく、両者が共通に動くなら、**現在の差を維持する支配的要因が精度である**という仮説を棄却します。過去の分岐を丸めが誘発した可能性までは否定できません。

第 2 仮説: **精度変更への応答より、共通する長い過渡・停止規則の影響が大きい。** 未収束のため未除外。

第 3 仮説: **複数の持続する流れ状態があり、丸めがその選択を変えた。** 未確認。(i) と (ii) は排他的ではありません。

判別 A/B: **起点を `case/45.isobutane_m6_d155/run_0451_v4_f32/res_140000.h5` に統一し、A＝float 継続、B＝FP64 切替。変更はビルド精度だけ。**

- 同一メッシュの保存量・乱流量・組成を `restart_field.py` で移し、FP64 側は float 保存値を正確に拡張する。メッシュ再生成、hoop、`pRef`、CFL、リミタ基準の変更を混ぜない。
- 各 **30,000 step 固定**、2,500 step 間隔。既存の到達判定で途中停止しない。θ_r 三断面・符号付き Q_w に加え、固定した収縮部領域の速度・温度・乱流量の差を追う。
- 判断尺度は各量の V4 到達差 D とする。**B だけが旧 FP64 値の方向へ 0.5D 以上移動し、A の移動が 0.1D 以下、かつ再実行差を十分上回る**なら、精度依存の維持機構を支持する。逆に、**A/B 差が 0.1D 以下で両者が共通に 0.5D 以上動く**なら、当該期間で精度変更が支配的という第 1 仮説を棄却し、共通の過渡を優先する。
- **両者とも留まる場合は判別不能**。短期の不変から二状態を認定してはいけません。元の二択を、どんな結果でも短期 A/B だけで決着できるという前提は受け入れません。収束・準定常 VERDICT は継続区間について別途記録します。

やらない方がよいこと:

- **`qAccumulatorFP64` の拒否だけを外すこと。** 軸の基準状態射影に加え、`sstEnergyIncludesK` の増分にも対応が必要です（[main.cpp:3395](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:3395)）。現行 accumulator は平均流五保存量だけで、4.9% が報告された ω は別の float 加算です（[update_d.cu:292](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:292)、[同:418](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/update_d.cu:418)）。
- **残差の double 加算を真因確認なしで実装すること。** 面値自体の丸め、生成・散逸項の評価、状態変換の誤差は加算だけでは直りません。
- **`pRef` を純粋な丸め対策の A/B と扱うこと。** 現在は閉性欠損があるため、ゲージ変更は厳密演算でも離散圧力残差を変えます。[axisymmetricSource_d.cu:51](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:51) と対流項から、自由な半径 DOF では変更分が ΔpRef·(ΣW_y−A) になります。
- **V5 の速さで精度フラグを相殺すること。** 仮に単価比 0.70 でも、報告された到達 step を使えば総時間比は 0.70×140/125＝**0.784**。約22%短縮の見込みですが、現状は同等品質に到達する費用の比較ではありません。

呼び出し側の前提への異議: **「止まった点＝定常状態」「`S_lostfrac`＝消失した更新率」「ω の拡散入力の除外＝V4 全体で拡散を除外」の三点を受け入れません。** 合成計算では更新が切り上がって `n_lost=0` でも `S_lostfrac=1/3` になりました。また §6.13 の除外範囲は、特定 Q32・壁際・ω 残差差の過半を説明するかに限定されます。

不足情報: **両 run の実効設定・入力とバイナリの照合記録、全期間の異常検査、判定区間付き VERDICT、収縮部の時系列と再実行差、V5 の結果。** V4 判定器は合成系列の Q_w=NaN でも `REACH` を返し得ることを確認したため、「比較成立」は保存成果物で再確認が必要です。ただし、これが実 run で発生したとは主張しません。

ファイル変更・forge 起動は行っていません。提案は **plan 未反映**です。
