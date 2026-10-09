# codex 諮問 (diagnose): line-viscous-thermal-weak-mode

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-thermal-weak-mode.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-thermal-weak-mode.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `bc39d64a` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.9 min, rc=0
- **結論**: **次は (a) の全自由行へのスカラー対角復元を診断専用で実装し、値2と同じ壁拘束を保った σ=0／1 の上限なし2000 step A/B を行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 「ρE 94%だから温度主体の弱いモード」 | **要再検証／Major** | [linedump_analyze.py:61](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/linedump_analyze.py:61) は保存量を `ρ・ρc・ρc²` で尺度化している。**ρE 成分には密度変化とエネルギー基準の寄与も含まれる**。弱い右特異ベクトル自体を `δρ/ρ・δu/c_ref・δT/T_ref` に変換して判定する。 |
| §6.2 の SVD を完了したという扱い | **要再検証／Major** | [plan:164](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:164) は拘束自由度の消去を要求するが、[解析コード:44](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/linedump_analyze.py:44) の `dec` はスカラー追加の除外にしか使われず、SVD は拘束行を含む全行列に適用される。壁で `δ(ρE)=e_wδρ` を代入し、固定運動量も消去した系で再評価する。 |
| 「host と CUDA が合ったので実装経路に誤りなし」 | **限定して採用／Major** | 保存 JSON は最大相対差約 `4e−12` を示す。ただし [解析コード:92](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/linedump_analyze.py:92) は全保存量をまとめた最大値比。確認できたのは**採取された系の解法との整合**であり、D/K の組立・物性入力・実残差との整合ではない。成分別誤差と線形後退誤差を併記する。 |
| 次の候補 (a)〜(d) | **(a) を診断用途に採用／Major** | 入口側3ラインでは σ=1 により保存済み補正ノルムが σ=0 の `0.0395〜0.0465` 倍になる。一方、温度だけを原因とする証拠はない。**全自由行への対角復元だけ**を変えた試験にする。(b)(c) の先行と、(d) の系統終了は却下。 |
| 等温壁拘束と、壁側の熱伝導 K を消す処理 | **式として採用、実データでの除外は保留／Major** | [timeIntegration_d.cu:1144](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1144) と [block_dplur_jacobian_d.cuh:138](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:138)。静止壁・固定組成では拘束により δT_w=0 となるため、熱伝導 K の作用もゼロになる。対案は拘束を変更することではなく、採取状態で拘束残差と δT_w を確認すること。 |
| 速度改善の次候補 | **double の逆行列保存を次候補に採用、C は後回し／Major** | [timeIntegration_d.cu:2072](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2072) は各 sweep・各節点で LU 代入を行う。まず部分ピボット付き double LU から逆行列を作り、代入だけを行列ベクトル積に替える独立した opt-in 実験とする。**高速化は未確認**。既定 double 維持、逆行列化と float 化を分けるという前回 M1・M2 は引き続き採用。 |
| 補正値の表記 | **修正を採用／Minor** | [解析コード:77](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/linedump_analyze.py:77) の δT は緩和前。JSON の `implicitRelax=0.7` を掛けると、入口3ラインの sweep 0 は **5.99／6.92／5.47 K**。また列2183でほぼ不変なのは最大密度・温度補正であり、保存量補正ノルムは **0.448倍**、最大速度補正は **0.137倍**。量を区別して記載する。 |

結論: **次は (a) の全自由行へのスカラー対角復元を診断専用で実装し、値2と同じ壁拘束を保った σ=0／1 の上限なし2000 step A/B を行う。**

第 1 仮説: **ライン面の全行スカラー対角を除いたことで、薄層近似と実残差の不一致に対する減衰が不足し、初回から大きな補正を返している。全行への対角復元で短期増幅を抑えられる可能性が高い。** 確度: **中**

根拠:

- [timeIntegration_d.cu:961](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:961) では、値2がライン面のスカラー対角を薄層 D/K に置き換える。
- `case/45.isobutane_m6_d155/run_0288_dump_lvc2/` に対応する[保存解析 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/linedump_run_0288_analysis.json)を再集計すると、入口3ラインで σ=1／0 の最大補正比は、密度 **0.037〜0.043**、速度 **0.033〜0.036**、温度 **0.041〜0.045**。これは特定の保存量だけに限らない強い感度である。
- 同3ラインでは sweep 0→4 の補正ノルム増加は **約1.7〜2.8%**。初回の大補正を「後続 sweep で初めて生じる増幅」で説明する必要はない。
- 提供された `run_0261_vn1_lvc2/` の29 step非有限、`run_0265_diag_linedir_cap50_lvc2/` の成長、point Δτ では急激に壊れないという記録とも整合する。ただし、これらの run 本体は手元では再確認できない。

反証条件: **σ=1 が初回補正を抑えても登録した短期安定条件を満たさなければ、「全行対角復元だけで十分」という仮説を棄却する。** 成功しても「熱伝導が真因」「Jacobian が正しい」とは判定しない。

第 2 仮説: **実残差との不一致、特に対流・ライン外結合・粘性仕事・勾配項が外側の更新を不安定化する。** 確度: 中。薄層モデルであることは [block_dplur_jacobian_d.cuh:94](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94) に明記されている。初回の5 sweepだけでは、以後の状態更新に伴う増幅を除外できない。

第 3 仮説: **壁フラグ・物性・状態の実入力が想定と異なる。** 確度: 低、未確認。κ の組立式は仕様と整合するが、[timeIntegration_d.cu:974](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:974) に渡った面の値は保存 JSON から検証できない。現時点で主因には据えない。

判別 A/B:

**変更する一点は、値2のライン面スカラー対角への係数 σ=0／1。**

- **A:** 現行の値2、σ=0。
- **B:** 同じ薄層 D/K に、既存の `scalar_visc_line` 相当を全自由行の対角へ加える、σ=1。壁・軸の拘束は両側で同じにする。
- 同じ新バイナリ、同じ `run_0183_ns_coldmesh_tw300_ext/res_100000.h5` 由来の保存量、キー5、方向別、上限なし、`cfl_pseudo=4`、緩和0.7、同じ sweep 数。新しい別々の run に最大2000 step。**上限50の組は同時に追加しない。**
- 初回に D/K・各 sweep の RHS・緩和前後の補正・分解失敗件数を保存する。σ=1 の初回解が、同一入力の host 再解と整合することを先に確認する。
- 全残差を毎 step確認し、序盤は局所の `δρ/ρ・δu・δT` と非物理値も追う。SVD を併記する場合は拘束消去を実施し、弱いモードも原始量へ変換する。

**事前の分岐:**

- **Aで既知の早期増幅が再現し、Bが非有限・非物理値・分解失敗なし、全残差の最大／開始≤10、末尾500 stepの傾き×500<0.1桁を満たす:** 「全行対角復元がこの条件の短期安定化に十分」を支持する。
- **Bがいずれかを満たさない:** 十分性を棄却する。初回補正だけ小さくなったことを成功とはしない。
- **Aが再現しない:** 比較は判別不能。バイナリ・実効設定・restart の対応を確認する。

これは既存 V-n1 の短期判定を使う試験であり、収束判定ではない。両側の `check_convergence` の VERDICT を別に記録する。V-n2へ進む際は、目的量の `check_quasisteady` と共通の終了条件までの総壁時計を要求する。

やらない方がよいこと:

- **［Major］エネルギー行だけを正則化すること。** 保存量のρE割合を根拠に行を選べない。全行の同一スカラーはエネルギー基準の変更に対して整合するが、エネルギー行だけの追加は密度との連動を変える。
- **［Major］熱伝導 K の除去を「値0と同じ温度処理」と呼ぶこと。** 値2には粘性仕事の D/K と異なる壁拘束も残る。今回の候補(c)は原因を一意に分離しない。
- **［Major］速度改善を上記A/Bへ混ぜること。** 逆行列保存は D/K と組立・保存精度を固定し、LU・逆行列・W・前進代入を double に保つ別実験にする。判定は、緩和前のスケーリング済み後退誤差、独立解との差、拘束残差、分解失敗件数で行う。例えば **η≤1e−11、物理尺度で正規化した補正差≤1e−8、新規失敗0**を実装前の候補基準として登録する。性能には逆行列生成費用を含め、native・専有GPU・profilerなしの反復計測と、同じ品質までの総時間を使う。
- **［Major］Bの3実装が遅かったことから、Thomas高速化全体を終了すること。** [速度 plan:122](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:122) の「全ラインが同時に走り、直列連鎖で決まる」は実行時間だけでは確定しない。失敗した実装の不採用と、律速原因の確定を分ける。

呼び出し側の前提への異議:

**「減衰を失って大補正が出る」という見立ては整合するが、「熱伝導の最小固有モードへのほぼNewton補正」という限定は支持できない。**

ここで求めたのは、質量項・対流・境界を含む非対称な全行列の特異ベクトルであり、熱伝導作用素の固有ベクトルではない。また、ρE主体と温度主体は異なる。コードの式から作れる反例として、静止・δT=0の密度変化でも `δ(ρE)=eδρ` となり、`e=400000 J/kg、c_ref=300 m/s` なら保存量ノルムの **95.2%がρE**になる。これは run の測定値ではなく、分類指標の限界を示す代数的反例である。

したがって、「密度90%」という旧指標で不合格だった記録は残してよいが、**密度を含む連成機構の物理的除外には使えない**。壁から6〜12節点付近という位置も、δTモードの形状、局所の質量項・輸送係数・格子幅との対応がなく、理由を特定できない。

TP の e の基準も、まずはこの**診断指標の依存性**として扱うべきで、solverの不整合と断定しない。壁拘束と温度微分が同じ基準を使えば、基準のずれだけでδTは変わらない。

不足情報:

`run_0288_dump_lvc2_linedump/` の生の D/K・状態・RHS、対象の case/45・case/52 の run 本体、実効 YAML、残差CSV、VERDICT、バイナリ識別情報はこの作業ツリーにない。今回直接確認できたのはコードと保存解析 JSON である。速度案Aの20 step比較は JSON の8量すべてが再実行差の3倍以内と再確認できたが、速度の計測原本は未確認。

ファイル変更・forge起動なし。**plan 未反映**。呼び出し側で、粘性 plan §6.2・§6.3・§6.1、速度 plan §5.1 #7・§6へ採否と事前判定を反映すること。
