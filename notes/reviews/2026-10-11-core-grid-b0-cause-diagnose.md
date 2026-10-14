# codex 諮問 (diagnose): core-grid-b0-cause

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-b0-cause.md`](../../notes/reviews/briefs/2026-10-11-core-grid-b0-cause.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `d51af1d0` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.9 min, rc=0
- **結論**: **S1 を修正し、同一の保存量から `lineDtDirectionalCap` だけを 50／1 に変えて、更新前の残差ベクトルを比較する A/B を先に行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **S1 を修正し、同一の保存量から `lineDtDirectionalCap` だけを 50／1 に変えて、更新前の残差ベクトルを比較する A/B を先に行う。**

第 1 仮説: 擬似 dt の拡大が平均流・SST の連成反復の減衰を悪化させ、非ゼロ残差を持つ状態に停滞させている。**確度: 中**
  
  根拠:
  
  - `case/45.isobutane_m6_d155/run_0615_cg_gref_cap50/` 対 `run_0616_cg_gref_cap1/` は、plan の記録上、上限だけの変更。壁際の dt 比は約 1/50、後半の平均流残差比は 0.006〜0.041。ただし `check_convergence` は **A: NOT CONVERGED、B: NOT CONVERGED（RISING）**。これは dt への感度の証拠であり、安定化完了の証拠ではない。[plan §4.33](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:704)
  - コードではライン面を CFL の集計から除外し、全面の値による上限を課す。[setDT_d.cu:166](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:166)、[同:194](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:194)
  - SST は平均流のライン解法とは別の `applySSTPointImplicit` で更新され、同じ `dt_local` の影響を受ける。したがって、まだ「流れ方向の結合だけが原因」とは特定できない。[main.cpp:2378](/home/sano/work/forge-coregrid/solver_density_cuda/main.cpp:2378)

  反証条件: **同一の実効状態・境界処理・係数に対し、上限だけで更新前の残差が再現誤差を超えて変わること。** この場合、「右辺は不変で反復だけの問題」という第 1 仮説を退ける。ただし、反復の問題も併存する可能性は残る。

第 2 仮説: 上限変更が dt 依存のソース処理などを通じ、実効残差を変えている。**確度: 低・未確認**。`dt_local` の検索には凝縮ソースへの受け渡しも現れるが、この run で有効かは確認できていない。検索の先頭 30 件だけで右辺への依存を除外してはいけない。

第 3 仮説: 残差は同じでも、有限回の線形反復・更新制限・丸めなどで補正が十分に適用されず、動かなくなっている。**確度: 低・未確認**。「揺れの平均がずれる」以外にも、右辺を変えずに停滞する経路はある。

判別 A/B: **上限 50／1、時間更新前の残差評価を各 1 回。**

- 起点は `case/45.isobutane_m6_d155/run_0614_cg_gref_point_ext2/res_80000.h5`。同一メッシュへの保存量の index コピーを一度行い、その入力を複製する。両腕とも B0 の構成とし、変更は `time.deltaT.lineDtDirectionalCap` の 1 点だけ。バイナリ・環境変数・BC・リミッタの実効基準値も固定する。
- **更新回数は 0 回**。既存の実行経路を使うなら各 1 step に限定し、最初の `outer_begin` と同じ時点の残差を採取する。コードの順序は `assembleResidual` → `logResidualSnapshot` → `setDT` → 更新であり、`outer_begin` は更新前である。[main.cpp:2323](/home/sano/work/forge-coregrid/solver_density_cuda/main.cpp:2323)、[同:2350](/home/sano/work/forge-coregrid/solver_density_cuda/main.cpp:2350)
- ただし残差組立の前処理では壁・軸の状態射影などを行う。入力ファイルの一致だけで済ませず、**実際に残差評価へ渡る状態の一致と、設定に応じた dt の有効化**を確認する。[main.cpp:1880](/home/sano/work/forge-coregrid/solver_density_cuda/main.cpp:1880)
- 比較対象は平均流・SST・化学種の全 `res_*`。各成分について `max|R₅₀−R₁|`、差の L2 ノルム、最大差の節点位置を出す。全 RMS は補助値とする。同設定での再評価差を ε と測り、差が 10ε を超えて再現するかを判別基準とする。決定的に再現する経路ならビット一致を基準にできる。

  → **A：残差ベクトルが再現誤差内で一致**なら、「この参照状態で、上限が直接右辺を変えている」を退ける。反復・補正適用側の診断を優先する。ただし、B0 と point の全設定にわたる右辺同一性や、他の状態での同一性までは証明しない。  
  → **B：状態は一致するが残差が再現誤差を超えて異なる**なら、純粋な反復側だけの説明を退け、最初に差が生じる残差組立処理を追う。差の存在だけで既存の 0.6 % の偏り全体を説明したとはしない。  
  → 状態が一致しない場合は試験不成立。初期化・境界処理の差を先に切り分ける。

やらない方がよいこと: **S2・S3 を今すぐ追加すること。** まず右辺の比較を成立させる。point の場から B0 で離れたことだけでは、右辺変更と反復不安定性を区別できない。この点は既に §4.32 で採用済みであり、同じ判別の誤りを戻さない。

呼び出し側の前提への異議:

- **Major — S1 の「全 RMS が同じなら右辺変更を退ける」は却下。** ノルムが等しくても残差の符号・空間分布は異なり得る。根拠はブリーフの S1。対案は上記の成分別残差ベクトル比較。
- **Major — 「流れ方向の結合は陽だから不安定」は要再検証。** ライン外結合は内側反復で遅れた補正を使う構造で、単純な陽時間積分と同一ではない。また dt は SST にも効く。対案は「連成反復の減衰不良」を第一仮説とし、方向と変数への帰属を保留する。[methods/time_integration/implementation.md:67](/home/sano/work/forge-coregrid/methods/time_integration/implementation.md:67)
- **Major — point の参照を厳密な定常解として扱うことは却下。** §4.31 の記録は、point 通算 0〜200k の残差が **NOT CONVERGED**、160k〜200k の派生量が **8 量 STEADY**。対案は「量の資格を満たした参照状態」と呼ぶ。そこから動くこと自体は不安定性の証明にならない。[plan §4.31](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:660)

不足情報: AWS の実行バイナリと確認したソースの対応、実効設定、更新前の残差ベクトル、既存出力でその時点を採取できるか。指示に従い生の run データは読んでおらず、数値と VERDICT は plan の記録を引用した。ファイル変更・forge 実行はしていない。**plan 未反映**。呼び出し側が §4.34 に登録する診断案である。
