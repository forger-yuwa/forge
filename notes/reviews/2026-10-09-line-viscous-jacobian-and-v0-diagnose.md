# codex 諮問 (diagnose): line-viscous-jacobian-and-v0

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-jacobian-and-v0.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-jacobian-and-v0.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `dc3bf1cc` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 7.7 min, rc=0
- **結論**: **薄層近似と壁拘束を明文化し、行列作用の単体照合を通したうえで、次のCFDは同一バイナリ・同一起点の `lineViscCoupling: 0 / 2`、方向別dt・上限なし・2000 stepのA/Bに絞る。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 判定・重大度 | 根拠と対案 |
|---|---|---|
| §4.1 の式・符号 | **採用。ただし近似の定義を修正［Major］** | 仮定した流束 `τ=βPΔu`、`q=κΔT` に対して、提示された D・K の符号と仕事の微分は正しい。非零速度・異なる密度・エネルギー基準のオフセットを含む100組の独立モデルで、数値微分との相対誤差は最大 `3.64e−15`。ただし、実装の勾配を凍結した二点微分は **βI** であり、**βP ではない**。[viscousFlux_d.cu:205](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:205)。P は転置・発散項の法線方向微分も近似的に取り込むものとして採用し、「勾配項の微分を全部除いた厳密 Jacobian」とは記述しない。 |
| 壁の K をゼロにする扱い | **条件付き採用［Major］** | **粘性の追加分だけ**について、実際に速度を固定する壁では運動量・仕事、温度を固定する壁では熱伝導の隣接寄与をゼロにする判断は妥当。ただし現在の等温壁行は、キー5では `Δroe=0` であり、温度固定の接線拘束 `Δroe=e_wΔρ` とは異なる。[timeIntegration_d.cu:1091](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1091)。後段の温度ピンを組み込んだ近似として明記し、対流 K を含む壁拘束全体が整合したとは扱わない。 |
| §4.3 の実装範囲 | **修正して採用［Major］** | 現在のカーネルには局所 `vis_lam` が渡らず、`cfg.visc` を使用する。また `thermCond/cp/fx` は熱伝導キーのビット1がないと `nullptr`。[timeIntegration_d.cu:1536](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1536)、[1562](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1562)。残差と同じ μ_f を使うには入力追加が必要。値2単独でも必要な配列を渡し、未対応の `heatCorrSU2`・弱形式壁などは明示的に制限する。 |
| U0・U2・V-n1・V-n2 | **現状の合格条件は却下［Major］** | U2の静止純伝導では、P・せん断・粘性仕事・高速TPの交差項を検査できない。U0も値1の不変性を検査していない。[line-viscous plan:86](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:86)。行列作用の単体照合、非物理値・線形solve失敗・末尾の成長率を追加する。V-n2の速度判定は「同じstepの量＋step単価」から、**同じ合格状態までの壁時計**へ変更する。 |
| V0 の「判定不能」 | **回帰の帰属不能として採用。合格への読み替えは却下［Major］** | JSONでは旧同士・新同士も不一致。ただし再実行対照はST0のみ、比較は8変数のみで、当初のST0/1・全 `VALUE/*` を覆わない。[v0_compare.py:16](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v0_compare.py:16)、[64](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v0_compare.py:64)。旧基準の不合格と原因判別不能を残し、代替基準は別版として登録して独立データで判定する。 |
| FREEZE_TURBで「両方が寄与」と確定 | **要再検証［Major］** | plan自身が `roK` の最大変化2.9%、壁・入口ピン除外未確認と記載する。[thermal plan:199](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:199)。まず非拘束節点の `roK/roOmega` 凍結を確認する。確認できても、分かるのは更新停止への感度であり、μ_tを含むSSTの作用全体を除いた試験ではない。 |
| `run_0224` の同じ不動点への到達 | **現時点では却下［Major］** | 保存済み時系列に判定関数を適用すると、末尾2万stepの θ_r・Q_w はすべて `DRIFTING`。後述のとおり、point参照にも判定窓の標本不足がある。残差上昇の場所を確認し、両側の準定常性と共通基準での残差をそろえてから判断する。 |
| Q_w を実際の壁入熱として使用 | **要再検証［Major］** | `cold_series.py` はCPG比較用の `cold_xcheck.reduce_fields` を呼び、熱伝導率を固定 `cp=1360`、`Pr=0.72` とSutherland粘性から再構成している。[cold_series.py:28](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_series.py:28)、[cold_xcheck.py:88](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:88)。現状は共通の後処理指標として扱い、TP残差と整合する物性・壁熱流束で確認するまで、MW値や熱収支一致の根拠にしない。 |

結論: **薄層近似と壁拘束を明文化し、行列作用の単体照合を通したうえで、次のCFDは同一バイナリ・同一起点の `lineViscCoupling: 0 / 2`、方向別dt・上限なし・2000 stepのA/Bに絞る。**

第 1 仮説: **方向別dtで露出する平均流の不安定には、温度・せん断・粘性仕事のライン内結合の不足が寄与している。** 確度: **中**。

  根拠:

- 現行の値0は、粘性を全5行のスカラー対角に置き、ラインの隣接粘性結合を持たない。[timeIntegration_d.cu:950](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:950)。
- 対流Kのせん断・エントロピー固有値は法線速度Vに依存し、V=0ではその直接結合が消える。[block_dplur_jacobian_d.cuh:66](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:66)。
- 提供記録では、キー5で初期の成長は弱まるが、上限なしの `run_0221` は残差増大の登録条件を外している。ただし、この事実だけで粘性の非対角不足を主因とは確定できない。

**連続行の減衰を除くことは、キー1の不整合と同一ではない。** 提案する粘性行列は、各節点の速度・温度を保つ密度補正 `δQ=δρ(1,u,v,w,E)` を全行で消す。独立モデルでもこの零空間を確認した。キー1の「行0〜3には人工減衰を残し、行4だけ別の応答にする」構造とは異なる。

ただし、人工減衰を除いた結果、既存の対流・境界・ライン外結合の誤差が成長する可能性は残る。物理的な零空間を持つことは、反復法の安定性保証ではない。

**粘性仕事の注意点:** 実際の応力を `τ_res=βΔu+τ_grad` と書くと、勾配を凍結しても仕事の微分には `τ_grad·δū` が残る。提示式はこれも薄層流束で置き換えたモデルである。より残差に沿わせるなら、`δτ≈βPδΔu` としつつ、積の微分を

```text
δ(τ·ū) ≈ (βPδΔu)·ū + τ_res·δū
```

とする。どちらを実装するかを先に固定し、対応する流束モデルに対して単体照合する。

  反証条件: **値2が単体照合を通っても、値0で再現する同じ局所モードの成長を止められなければ、「今回の薄層ライン結合の追加だけで上限なしを安定化できる」という十分性の仮説を棄却する。** 粘性・熱伝導全般の関与まで棄却しない。

第 2 仮説: **壁の密度更新とエネルギー拘束、またはライン外の近似結合が残る律速である。** 確度: 中〜低。壁拘束の不整合はコードで確認できるが、今回の成長への寄与は未確認。

第 3 仮説: **SST分離更新が平均流の成長を増幅している。** 確度: 低〜中。FREEZE_TURBへの感度は報告されているが、凍結の実効性確認が未完了である。

判別 A/B:

**変更は `lineViscCoupling: 0 → 2` の1点。** `run_0183` の同じ保存量から、新しい同一バイナリ、`implicitThermalJacobian: 5`、`lineDtDirectional: 1`、`lineDtDirectionalCap: 0`、`cfl_pseudo: 4`、同じ精度・BC・緩和・sweep数で最大2000 step。旧 `run_0221` は参考とし、対照Aも新バイナリで作る。

- 全残差を毎step、局所帳簿を毎step、場を100〜200 stepごとに保存する。対象は既知の縮流部近壁と、全域の残差上位領域。
- 短期合格は、非有限・非物理値なし、線形solve失敗0、全残差最大が固定した開始尺度の10倍以内、末尾500 stepで残差と局所振幅の成長が検出されないこと。ゼロに近い列には事前に絶対尺度を定める。
- **Aで成長を再現し、Bで止まるなら:** 今回の薄層結合で短期安定化できないという説を退け、第1仮説を支持する。ただし、P・熱伝導・仕事・人工減衰除去のどれが効いたかまでは分離できない。
- **Bでも同じモードが育つなら:** この変更単独で十分という説を退ける。SSTや境界拘束が残っているため、値2そのものの一般的な無効性とはしない。
- Aが再現しない場合は判別不能。入力・実効設定・バイナリの対応を確認する。

SSTが分離したままでも、この**十分性の試験**として上限なしを試す意味はある。2000 stepの合格は短期安定の判定に限り、高速化・収束・定常解一致とは呼ばない。

やらない方がよいこと:

- **V0の事後データを、そのまま新しい合格基準にも検証データにも使うこと。** `V0_repeat.json` の「一致0/3」は不一致0/3ではない。また、3本から作った3組、3×3の9組は独立標本ではない。
- **「旧同士の最小〜最大に旧新差が入る」を採用すること。** 代替基準は、全対象変数の固定尺度によるRMS・最大差、必要な局所指標、許容差を先に登録する。再標本化するなら組ではなく独立runを単位にし、差の上側信頼限界が許容差内かで判定する。ST0/1を両方含める。許容差を正当化できなければV0は保留のままにする。
- **V0未解決のまま `run_0191` を正式な後退判定の対照にすること。** 同一新バイナリ内の診断は進めてよい。旧runとの正式比較は改訂V0を独立データで通した後とし、それまでは参考比較に限定する。
- **`run_0211` の発散を、密度への人工拡散だけに帰属させること。** 値1はK追加と対角半減を同時に行う。一定補正に対するライン内の減衰も大きく失うため、観測だけでは分離できない。
- **場所を確認せず `run_0224` 系列を無条件に延長すること。** 次の延長判断より先に、保存済み場で `res_ro` と運動量残差の局所ノルム・寄与領域・移動／固定・振幅成長を調べる。総和の欠損減少は、局所残差の相殺増大とも両立する。

呼び出し側の前提への異議:

**`run_0224` は、pointの値に近づいた未収束過渡であり、同じ不動点に到達した証拠ではない。** 保存JSONに対して、`check_quasisteady.py` の `classify_series` を書き込みなしで実行した。θ_rの登録閾値に合わせてドリフト閾値を `0.0005`、窓を末尾2万stepとした結果は次のとおり。

| 根拠run | 判定区間 | 結果 |
|---|---|---|
| `case/45.isobutane_m6_d155/run_0224_ns_coldmesh_tw300_linedir_tj5_cap50_ext/` | 40000〜60000、5枚 | θ_r(40/70/94): **DRIFTING**。両端差は +1.373 / +1.413 / +1.434%。後処理Q_wも同じ閾値では **DRIFTING** |
| `case/45.isobutane_m6_d155/run_0217_ns_coldmesh_tw300_cfl4_ext3/` | 180000〜200000、3枚 | **TRANSIENT-UNSETTLED: only 3 snapshot(s) (<4)**。両端ドリフトの登録条件は満たすが、この窓では標本不足 |

これは保存済み抽出系列の判定であり、元HDF5の再検査ではない。`run_0224` の残差については、提供記録の **`NOT CONVERGED`、区間0〜59999、`rms_ro/rms_roUy` が `RISING`** を採用するが、原CSVは確認できていない。

「同じ不動点」と判断するには、両側で十分な標本を持つ準定常窓、全残差の判定、局所残差・場の発達、整合した壁熱量が必要である。切り戻し2〜4万stepも固定長だけで合格にせず、point側の遅い緩和を検出できる窓と許容差を登録する。元の点解法が数十万stepを要した以上、短い切り戻しで動かなかったことだけでは十分でない。

不足情報: 対象の `run_0183/0191/0203/0211/0217/0222〜0248/0252` の実runディレクトリは、この作業ツリーにない。残差CSV、元HDF5、実効YAML、各VERDICT、バイナリ照合は直接再検証できていない。確認できた実測資料は [V0_repeat.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/V0_repeat.json) と各抽出系列JSONである。run索引は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md)。

**ファイル変更・forge起動なし。plan未反映。** 呼び出し側がライン粘性planの§4・§6と、熱伝導planの§6.0・§6.1へ採否を反映する診断結果である。
