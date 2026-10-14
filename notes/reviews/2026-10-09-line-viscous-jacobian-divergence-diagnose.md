# codex 諮問 (diagnose): line-viscous-jacobian-divergence

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-jacobian-divergence.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-jacobian-divergence.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `d61892c0` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.5 min, rc=0
- **結論**: **次は `lineViscCoupling: 2` の初回の実ライン行列を取り出し、従来のスカラー対角を「戻さない／戻す」だけの host 再解 A/B で、密度主体の弱いモードが大補正を生むか確認する。**
- **extra**: `case/52.conjugate_slab/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 「粘性の零空間だから、ライン全体もほぼ特異」 | **要再検証／Major** | 粘性部分の零空間は確認できるが、全行列には質量項、ライン外の A⁺、境界半割面、軸対称ソースが残る。[timeIntegration_d.cu:841](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:841)、[872](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:872)、[1046](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1046)。**実際の行列と RHS で、弱いモードが励起されているか測る。** |
| 等温壁の行 `[−e_w,0,0,0,1]` | **採用／Major：原因からの除外は不可** | 静止壁・固定組成では ΔT_w=0 と整合する。ただし値 0→2 は、この拘束も同時に変更するため、粘性結合だけの比較ではない。[timeIntegration_d.cu:1128](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1128)。対案は、次の A/B では壁拘束を両側で同一に保つこと。 |
| 対策 (a)〜(e) | **(a) を診断用途に限り採用／Major** | 失った全行スカラー対角を戻す操作は、壁拘束・粘性 D/K を保ったまま減衰除去の効果を調べられる。[timeIntegration_d.cu:951](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:951)。**恒久採用は保留**。(b)(c) は連続行だけ変更する根拠不足、(d) の R=50 は提供記録の `run_0265` で成長、(e) の機能削除は早計。 |
| V-n1・U2・U3 の結果判定 | **不合格／判定不能の記録を採用／Major** | [plan §6](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:105) の基準では、V-n1 B は不合格、U2 は **950>250 step で速度要件不合格**、Couette は試験不成立、Poiseuille は提供値 **速度誤差12.1%>0.1%で不合格**。閾値を変更して救済せず、解析解照合・短期安定・残差収束を別々に記録する。 |
| Poiseuille による U3 の代替 | **採用／Major：仕事の正しさの証明は未了** | `bodyForce` は **N/m³**。運動量に fV、エネルギーに f·u V を加える実装で、提示された解析解と整合する。[bodyForce_d.cu:39](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/bodyForce_d.cu:39)。ただし有限で走ったことだけでは仕事の Jacobian を検証できない。速度・温度の両方の許容差到達が必要。 |
| directional と point は同じ不動点に到達 | **現時点では却下／Major** | 保存系列の再判定で切り戻し側の θ_r(70・94) が `TRANSIENT-UNSETTLED`。また切り戻し変化は登録許容を超える。対案は、両側を同じ point 設定で評価し、全残差・局所残差・目的量の共通終了条件を満たしてから比較すること。 |

結論: **次は `lineViscCoupling: 2` の初回の実ライン行列を取り出し、従来のスカラー対角を「戻さない／戻す」だけの host 再解 A/B で、密度主体の弱いモードが大補正を生むか確認する。**

第 1 仮説: **全行スカラー対角の除去により、方向別 Δτ で密度・圧力主体の補正への抑制が不足し、初回から過大な補正を返す。** 確度: **中**。「ライン方向に一様な圧力モードでほぼ特異」という限定した機構は未確認。

  根拠:

- 値 2 は従来のスカラー対角を薄層 D/K に置換する。その薄層行列は δQ=δρ(1,u,v,w,E) を消す。[timeIntegration_d.cu:952](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:952)、[block_dplur_jacobian_d.cuh:117](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:117)。
- **提供記録では** `case/45.isobutane_m6_d155/run_0261_vn1_lvc2/` が29 stepで非有限、`run_0265_diag_linedir_cap50_lvc2/` は上限50でも成長、`run_0264_diag_lineonly_lvc2/` は100 stepでは同様の爆発を示さない。この Δτ 依存とは整合する。
- ただし、壁法線に T・u が変わる場では δρ一定でも δp=RTδρ は一定でない。各節点で粘性零空間に属することと、対流・境界を含むライン全体の零空間に属することは別である。

  反証条件: **無次元化した実行列の弱いモードが密度主体でない、実 RHS がそのモードをほとんど励起しない、または大補正が最初のライン解ではなく後続 sweep の RHS 増幅で初めて現れる場合、「密度モードの準特異性が初発原因」という説明を棄却する。**

第 2 仮説: **実残差と前処理行列の差、特にライン外の lag・対流結合・壁拘束が反復を増幅する。** 確度: 中。ライン外は sweep ごとの lag であり、行列が非特異でも反復安定性は保証されない。[timeIntegration_d.cu:880](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:880)。対案は各 sweep の RHS と補正を分けて保存すること。**［Major］**

第 3 仮説: **粘性仕事などの薄層近似と実残差の不一致、または CUDA への組み込み誤り。** 確度: 低〜中、未除外。U-J は選んだ薄層モデルとの照合であり、実応力は勾配項も含む。[viscousFlux_d.cu:205](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:205)。対案は実行列の host 解と CUDA 補正の照合を今回の計測に含めること。**［Major］**

判別 A/B:

**変更は、ライン面から失ったスカラー対角に掛ける係数 σ=0／1 の一点。非線形計算は初回1 step分の採取だけとし、比較自体は同じ凍結状態・同じ RHS の host 再解で行う。**

- **A:** 現在の値 2 の行列 M。
- **B:** M に、各ライン面の従来係数 `2ν_eff·δ/dcc` を全5行の対角へ戻す。壁・軸の拘束行は従来どおり最後に上書きする。粘性 D/K、等温壁拘束、Δτ、物性、RHS は A と同一。
- 対象は提供記録で初発位置となった列12・33・65を含む**実際のライン全長**。壁半割面、ライン外 A⁺、軸対称項を省略しない。初回の各 sweep の RHS・補正・`line_fail_d` を採取する。
- 固定した物理尺度で行・列を無次元化し、拘束自由度を消去して、弱い特異モードと実 RHS の投影を調べる。生の保存量行列の条件数だけでは判定しない。
- 補正は保存量から  
  `δu=(δ(ρu)−uδρ)/ρ`、  
  `δT={δ(ρE)−u·δ(ρu)+(½|u|²−e)δρ}/(ρc_v)`  
  を作って評価する。出力時刻の疑いがある `VALUE/T`・`Ux` の差は使わない。

**事前の分岐:**

- **A の大補正が密度主体の弱いモードに集中し、B でその応答が抑えられるなら:** 第1仮説を支持する。診断の目安は、無次元補正の密度部分が二乗ノルムの90%以上、B の密度補正ノルムが A の1/10以下。この場合も「粘性仕事まで正しい」「長期安定」とは判定しない。
- **その構造が見えないなら:** 狭い意味の準特異仮説を棄却する。host 解が CUDA と合わなければ実装経路、合っていて後続 sweep のみ増幅するなら第2仮説へ進む。
- **B が単に全補正を小さくしただけなら:** 正則化への感度を示しただけで、圧力モード説の支持とはしない。

やらない方がよいこと:

- **［Major］(b)(c) の連続行だけの正則化を先に実装すること。** 現在の零空間は5保存量が連動した方向であり、連続行だけが問題だという証拠はない。根拠は上記の共通 Jacobian。まず全行対角の診断 A/B に限定する。
- **［Major］U3 の遅さをライン外スカラーだけで説明し、高 AR にすれば検証が済むとすること。** ライン外には音響対角も残る。[block_dplur_jacobian_d.cuh:66](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:66)。[U3生成条件](/home/sano/work/forge-integ-1005/case/52.conjugate_slab/u1_thermjac.py:32)の静止一様状態での概算では、接線運動量に作用する音響対角 `cΔy≈0.217` は、両側スカラー合計 `4νΔy/Δx≈0.00385` の**約56倍**。これは run 実測ではないが、スカラーだけの律速説明は不十分である。高 AR 試験は後日の別試験として意味があるが、元の U3 不合格を置き換えない。
- **［Major］値2を本線の加速へ採用すること。** V-n1が不合格であり、V-n2へ進む条件を満たさない。値0の directional も、現状は point 仕上げを要する加速候補に留める。採用には、仕上げ込みで同じ終了条件までの壁時計短縮を示す必要がある。

呼び出し側の前提への異議:

**観測として受け入れるのは「提供記録上の発散・Δτ感度」であり、「準特異」「運動量の符号誤りを除外」「同じ不動点」は解釈である。**

**［Major］切り戻しは同じ不動点の証明になっていない。** 保存済み `series_run_*.json` に `check_quasisteady.py` の `classify_series` を書き込みなしで適用し、末尾2万 step・5点・`drift=0.0005` で次を再確認した。

| 根拠 run（`case/45.isobutane_m6_d155/` 配下） | 判定区間 | VERDICT |
|---|---|---|
| `run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2` | 40000〜60000 | θ_r 3点・後処理 Q_w：`STEADY` |
| `run_0262_ns_coldmesh_tw300_cutback_point` | 20000〜40000 | θ_r(40)・Q_w：`STEADY`、θ_r(70・94)：`TRANSIENT-UNSETTLED` |
| `run_0263_ns_coldmesh_tw300_cfl4_ext4` | 20000〜40000 | 全量 `STEADY`、ただし単調変化あり |

`run_0263` のツール外挿は θ_r(40/70/94)=**0.0777773／0.114823／0.127632**。これは漸近値の推定であり、収束証明ではない。`run_0252` の残差判定は plan 記録で **`NOT CONVERGED (stalled/plateau)`、0〜59999**。`run_0262/0263` の全残差 VERDICT は直接確認できていない。

さらに保存された局所残差統計では、`run_0252` 最終→`run_0262` 最終で Σ|res_ro| は **2.417→0.493** に下がる一方、`RMS(res_ro/V)` は **5.84×10⁴→4.57×10⁵** に増える。**総和の改善を局所残差全般の改善とは呼べない。** 対案は共通の局所ノルムと分布を終了条件に含めること。

**［Minor］数値の転記も修正が必要。** 保存系列では切り戻しの後処理 Q_w の開始比変化は **+0.00932%**で、ブリーフの +0.001% ではない。θ_r(70) は **−0.08094%**。合否は変わらないが、系列から再生成する。

不足情報:

対象の case/45 `run_0260〜0267`、case/52 の run 本体はこの作業ツリーに存在しない。したがって残差CSV、元HDF5、実効YAML、各 `CONVERGENCE_VERDICT.txt`、実行バイナリの対応は未検証。U-J の PASS 数値も提供記録として扱った。今回必要なのは、**初回の実 D・Kprev・Knext・各 sweep の RHS／補正・境界フラグ・Δτ・物性とバイナリ識別情報**である。

ファイル変更・forge 起動なし。**plan 未反映**。呼び出し側で本診断の採否と A/B を、ライン粘性 plan §4・§6・§6.1 に反映すること。
