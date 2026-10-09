# 諮問: `lineViscCoupling: 2` (薄層の粘性・熱伝導の Jacobian) が方向別の擬似 dt で即座に発散する理由と、次の一手

日付 2026-10-09。諮問先 codex (diagnose)。作業ツリー `/home/sano/work/forge-integ-1005` (commit 42dd9858 以降)。
エスカレーション条件: 2 (発散)、3 (事前登録の比較の不合格・判定不能)、4 (plan に無い修正の前)、6 (cuda_forge の数値の変更の前)。
plan: `plans/active/time_integration-line-viscous-jacobian.md` (§4.1・§4.3・§6・§6.1 を全文読むこと)。前回の諮問: `notes/reviews/2026-10-09-line-viscous-jacobian-and-v0-diagnose.md`。
実装: `solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh` の `accumulate_thinlayer_visc_jacobian`、`timeIntegration_d.cu` の `implicit_defect_correction_block_d` の
`isLineFace && lineViscCoupling == 2` の分岐と等温壁の行 (`lineViscCoupling == 2` で [−e_w,0,0,0,1])、wrapper の `FORGE_BDPLUR_ARGS`。
バイナリ: 5ab83056 + typedef double (sha256 6631a87e…、`~/forge-linevisc-fp64`)。

## 観測事実

- **U-J (host の単体照合、`solver_density_cuda/tools/test_line_visc_jacobian.cpp`)**: PASS。中心差分との差 2.3e-8、零空間 δρ(1,u,v,w,E) 2.5e-16、等温壁の拘束 1.5e-16、
  8 節点ラインの Thomas と密行列の差 1.6e-15、float と double の差 1.3e-5。
- **U2 (case/52 の流体の板、純伝導、静止、一様 20×16、ライン 21 本 × 17 節点、ライン + 方向別)**: L∞ < 0.05 K に入る step と 5000 step の L∞
  — 値 0 cfl 5: 未到達・0.198 K (`run_0008`)、値 2 cfl 5: 未到達・0.094 K (`run_0009`)、値 0 cfl 50: 非有限 (`run_0010`)、値 2 cfl 50: 950 step・1.4e-6 K (`run_0011`)。
  → 熱伝導の項のラインの配線は効いている (値 0 が壊れる cfl で値 2 は収束)。事前登録の「250 step 以内」は不合格。静止流なので運動量・仕事の項は試していない。
- **U3 (せん断)**: 上壁を z 方向に動かす版 (`run_0012〜0015`) も、側面を周期にして x 方向に動かす版 (`run_0016〜0019`) も、壁の節点の速度が 0 のまま (u が立たない)。
  node の壁の速度は `nodeWallDirichlet` で no-slip に固定する作りで、bcond の Ux/Uz は node では効かない (`boundaryCond_d.cu` の `wall_isothermal_d` は cell の ghost 用)。→ 判定不能。
- **V-n1 (case/45、run_0183 の res_100000 から、同じ新バイナリ、キー 5・方向別・上限なし・cfl 4、2000 step)**:
  A = 値 0 (`case/45.isobutane_m6_d155/run_0260_vn1_lvc0`): 残差の最大/開始 ρ 572・ρv 942・ω 280 まで振れて回復、末尾 500 step の傾き × 500 は −0.13〜−0.34 桁 (k だけ +0.04)。
  B = 値 2 (`run_0261_vn1_lvc2`): **1 step 目から爆発**、rms_ro 8.0e-6 → 1.8e-4 (1) → 3.9e-3 (2) → 2.6e-2 (3)、rms_roUy 1.0e-3 → 6.0e-2 → 1.3 → 6.2、29 step で非有限。
  NaN の場所: ρ が列 60〜70 の全層 (0〜120) で非有限、T は列 65 の層 52 で 50 K (床)・層 53 で 4327 K (ライン方向の隣で交互)。
- **切り分け (計測)**: `run_0264_diag_lineonly_lvc2` (ラインだけ、Δτ は point、キー 5・値 2、100 step) は `run_0201` (ラインだけ・値 0・キー 0) と残差の推移がほぼ同じ
  (step 99 の残差/開始: ρ 0.585 vs 0.586、ρv 0.412 vs 0.414)。`run_0265_diag_linedir_cap50_lvc2` (方向別 + 上限 50 + キー 5 + 値 2、300 step) は成長
  (step 99: ρ 98・ρv 194・k 7000 倍、step 299: ρ 674 倍)、同条件の値 0 (`run_0223`) は step 299 で ρ 0.55 倍。
- **1〜3 step の場の変化 (`run_0266_diag_s3_lvc0` / `run_0267_diag_s3_lvc2`、毎 step 出力)**: 1 step 目の res_1 では T・Ux・Uy の変化が両方とも 0 (出力の原始変数が更新前の値の可能性、未確認)。
  ρ の変化の最大: 値 0 は 9.3e-3 kg/m³ (列 37・壁から 1 層)、値 2 は 0.29 kg/m³ (列 12・5 層)。P の変化の最大: 値 0 は 195 Pa、値 2 は 3619 Pa (列 33・壁の節点)。
  2 step 目: P の変化 1.5e3 vs 9.1e4 Pa (値 2 は列 33 の壁)、T 0.58 vs 6.9 K、3 step 目: P 3.0e3 vs 5.3e5 Pa。値 2 の変化は入口寄りの列 12〜70 の壁の近くに集まる。

- **U3 の組み直し (Poiseuille、事後に登録、`case/52.conjugate_slab/run_0020〜0023_u3p_*`)**: 側面を周期、両壁 300 K で止め、体積力 434.4 N/m³ (中央 300 m/s、+22.5 K)。
  5000 step の速度の L∞ (中央の速度に対する比): 値 0 cfl 5 は 0.276、値 2 cfl 5 は 0.259、値 0 cfl 50 は非有限、**値 2 cfl 50 は 0.121** (どれも許容差に未到達 = 不合格)。
  値 2 は値 0 が壊れる cfl 50 で安定 (運動量の結合の符号が逆なら壊れるはず、という間接の証拠)。ただし cfl を 10 倍にしても速さは 1.6 倍
  (1 step あたりの減衰 2.7e-4 → 4.3e-4)。U2 の熱 (値 2 cfl 50) は 5.6e-3/step で、運動量は 13 倍遅い。
  呼び出し側の見立て (未確認): この板は縦横比 1.6 なので、ライン外 (x 方向) の面に従来どおり入るスカラー 2α_x (全行) が、ライン方向の滑らかなモードの固有値 (約 0.04 α_y) の 40 倍で、ライン解の応答を抑えている。
  ノズルの近壁 (縦横比 約 4000) では α_x/α_y = 1/AR² で無視できる。
- **切り戻し (別 plan だが同じ問い)**: `case/45.isobutane_m6_d155/run_0262_ns_coldmesh_tw300_cutback_point` (run_0252 = directional の最終場を point cfl 4 で 40000 step) は
  θ_r(40/70/94) が −0.096 / −0.082 / −0.069 % 動き (事前登録の 0.05 % を超えて不合格)、θ_r(70・94) は末尾で TRANSIENT-UNSETTLED。Σ|res_ro| は 2.4 → 0.49 kg/s。
  point の延長 `run_0263_ns_coldmesh_tw300_cfl4_ext4` は θ_r が +0.03〜0.04 % 動き ALL STEADY (単調 +0.014 %/2 万 step)。両者の終値の差は θ_r +0.03 %・Q_w −0.05 %。
  詳細は `plans/active/time_integration-implicit-thermal-jacobian.md` §6.0 の末尾。

## 呼び出し側の仮説 (確かめていない)

ライン方向に一様な「速度と温度を保ったまま密度だけが変わる」モード δQ = δρ(1, u, v, w, E) (圧力の水準の変化) は:
(i) 対流の項 (流束の差) では内部でほぼ 0 の応答 (一様な摂動の流束の差)、(ii) 新しい粘性の項ではちょうど零空間、(iii) 方向別の Δτ では V/Δτ も小さい。
値 0 では全行のスカラー 2α (各面) がこのモードを抑えていた。値 2 は物理として正しく外したので、ラインの線形系がこのモードでほぼ特異になり、ライン解が大きな Δρ を返す
(流れ方向の結合は lag なので、ライン間の圧力の水準をそろえる力は 1 step の中では働かない)。Δτ への依存 (point なら無害、上限 50 で成長、上限なしで即座) と整合する。
ただし壁の節点の境界の半割面 (`has_nbr=false`) の対流の A⁺ が壁の端を固定しているはずで、本当に特異に近いかは計算していない。運動量・仕事の項の誤りも U3 が成立しないので除外できていない。

## 問い

1. この仮説は観測と整合するか。別の機構の候補は (特に運動量・仕事の項の誤りや、等温壁の行の拘束 [−e_w,0,0,0,1] の影響)。仮説を確かめる最も安い計測は (例: 1 本のラインの線形系を host で組んで最小特異値を出す、
   列 12〜70 の 1 step 目の ΔQ を δρ(1,u,v,w,E) に射影する)。
2. 対策の候補の妥当性と順序: (a) 値 2 にスカラー 2α I (全行) を残し、その上に薄層の D と K を足す (キー 5 と同じ考え、ただし人工の減衰が戻る)、(b) 連続の行だけスカラーを残す、
   (c) 圧力の水準を決める別の正則化 (例: 連続の行に σ·V/Δτ_point)、(d) 上限 R との併用だけを使う、(e) 値 2 は不採用で閉じる。
3. せん断・仕事の項の試験 (U3) を node でどう組むべきか。壁を動かす経路がないので、周期の側面 + 体積力の Poiseuille (u と T の解析解) で代わりにしてよいか (体積力が単位体積か単位質量かは未確認)。
4. 事前登録の V-n1 の判定 (B の不合格) の書き方。plan の §6 の U2 の不合格 (950 step)・U3 (Couette 判定不能、Poiseuille 不合格) の扱い。U3 を縦横比の大きい板で組み直す意味はあるか。
5. 切り戻しの解釈: directional (上限 50) の最終場と point の解は「同じ不動点」と言えるか、言えないなら何を足せば判定できるか。directional を本線の加速に使ってよい条件。

## 読んでよいもの

- 上記 plan と前回の諮問、`plans/active/time_integration-implicit-thermal-jacobian.md` §6.0 (run_0221〜0252)、`methods/time_integration/implementation.md` (line-implicit・v2・熱伝導 Jacobian)
- `solver_density_cuda/cuda_forge/timeIntegration_d.cu`、`block_dplur_jacobian_d.cuh`、`viscousFlux_d.cu`、`nodeWallDirichlet_d.cu`、`boundaryCond_d.cu`
- `case/45.isobutane_m6_d155/README.md`・`case/52.conjugate_slab/README.md` (run の一覧)、`case/52.conjugate_slab/u1_thermjac.py`
