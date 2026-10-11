# codex 諮問 (diagnose): line-viscous-lvc3-result

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-lvc3-result.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-lvc3-result.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `d88d717d` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.0 min, rc=0
- **結論**: **次は同じバイナリ・同じ初期場で `lineViscCoupling: 0` を固定し、`implicitThermalJacobian: 5／7` だけを変える最大2000 stepのA/Bを一組行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 「全行のスカラー対角復元だけで十分」を棄却 | **採用** | [plan §6.5:202](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:202) の記録では、`run_0301_vn1b_lvc3` は122 stepで非有限。事前登録の不合格に該当する。ただし今回は残差原本・VERDICTを独立確認できていない。 |
| 「薄層 D/K を足すこと自体が不安定化の原因」と断定 | **却下／Major** | [timeIntegration_d.cu:1146](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1146) で、値2・3は壁の拘束行も強制変更する。**薄層結合と壁拘束の効果が未分離**。同じバイナリの値0でキー5／7を比較する。 |
| キー7で悪化しなければ「壁拘束は除外」 | **限定して採用／Major** | 同じ箇所から、値0のキー5→7は拘束行だけの変更と確認できる。ただし除外できるのは**壁拘束変更単独で破綻する説**であり、薄層 D/K との相互作用は残る。 |
| ρvの先行増幅から P の4/3成分・フープ項を主因とする | **却下／Major** | [フープ項の組立:1058](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1058) は密度・エネルギーとも連成する。ρvの増幅は温度・圧力補正の結果でも起こり得る。まず壁拘束を分離する。 |
| 値2系統を閉じ、加速＋point仕上げを採用 | **現時点では却下／Major** | [thermal plan:251](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:251) では、切り戻しの登録基準は不合格、判定区間20000〜40000のθ_r(70・94)は `TRANSIENT-UNSETTLED`。**壁拘束のA/Bを一組だけ済ませ、値2・3は本線不使用のまま保留**とする。その後、仕上げ込みの共通終了条件までの壁時計を評価する。 |

結論: **次は同じバイナリ・同じ初期場で `lineViscCoupling: 0` を固定し、`implicitThermalJacobian: 5／7` だけを変える最大2000 stepのA/Bを一組行う。**

第 1 仮説: **今回の早期破綻には薄層 D/K の変更、またはそれと壁拘束との相互作用が必要で、壁拘束変更だけでは再現しない。** 確度: **中。ただし帰属は未確認。**

根拠:
- [組立コード:961](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:961) では、値3はスカラー対角を保持したまま薄層結合を入れる。それでも、記録上は `case/45.isobutane_m6_d155/run_0301_vn1b_lvc3/` が122 stepで破綻している。**対角を失ったことだけでは説明できない**。
- [薄層関数:94](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94) は実残差の厳密な微分ではなく、近似前処理行列である。対角の追加だけから外側反復の安定性は保証できない。
- 一方、壁の行は静止壁で ΔT_w=0 を課す整合した式であり、式そのものを誤りとする根拠はない。ただし、その変更が対流結合を通じて不安定化する可能性は残る。

反証条件: **値0＋キー5が2000 step非有限なしで進む一方、値0＋キー7が同区間で破綻すれば、「薄層 D/K が破綻に必要」という限定を棄却する。** 壁拘束変更だけで破綻を起こせることになる。ただし、値3と同じ局所機構かは別途の問題である。

第 2 仮説: **壁拘束変更単独が、壁密度とエネルギー補正の連動を変え、早期破綻を起こす。** 確度: **低、未確認**。[拘束行:1140](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1140) が根拠となる候補で、今回のA/Bで直接調べられる。

判別 A/B:

- **A:** `lineViscCoupling: 0`、`time.deltaT.implicitThermalJacobian: 5`。
- **B:** 同じ設定で `time.deltaT.implicitThermalJacobian: 7`。
- 両側とも `run_0183_ns_coldmesh_tw300_ext/res_100000.h5` から同一メッシュの保存量をビット一致で継ぐ。同じFP64バイナリ、方向別、上限0、`cfl_pseudo: 4`、緩和0.7、5 sweep、最大2000 step。既存runを上書きせず、新しいrunを二つ作る。実効YAMLの差が上記キーだけであることと、バイナリのSHA256を保存する。
- 全 `rms_*` を毎step記録する。序盤200 stepは壁と隣接層の δρ/ρ・δu・δT、拘束残差、最初の非有限・非物理値のstepと場所を保存する。200 stepごとの場だけでは、122 step以前の発生順序を追えない。

**事前の分岐:**

1. **Aは2000 step有限、Bは途中で非有限または非物理値:** 壁拘束変更単独での不安定化を支持し、第1仮説の「薄層変更が必要」を棄却する。200 step以内なら既存の早期破綻と時間尺度も重なるが、同じ原因とはまだ断定しない。
2. **両側とも2000 step有限・非物理値なし:** この条件・期間での「壁拘束変更単独による破綻」を棄却する。薄層変更と、その壁拘束との相互作用は残す。**過渡場の一致や収束を意味しない。**
3. **Aも破綻:** 基準が再現していないため帰属不能。旧 `run_0260` を対照として補わず、バイナリ・実効設定・restartを確認する。

全残差の最大／開始、最初の増幅時刻、末尾500 stepの傾きを併記し、両側で `check_convergence.py` のVERDICTと判定区間を保存する。**今回の因果判別に「最大／開始≤10」を流用しない**。対照の `run_0260` 自体が記録上ρ572倍まで振れており、その閾値は既存の短期安定合格基準として別に扱う。

やらない方がよいこと:

- **［Major］P→Iを先に試すこと。** [Pの使用箇所:114](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:114) と[仕事項:146](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:146)から、この変更は運動量だけでなく粘性仕事の線形化も変える。意味があるのは「薄層モデル中の法線補正への感度」の検査であり、改善してもフープ項の誤りを証明しない。実施する場合もD・K・仕事項で一貫して変更する必要がある。
- **［Major］壁側Kの消去を、拘束行とは独立した犯人として扱うこと。** [薄層関数:123](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:123)・[同:138](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:138)の項は、静止壁のΔu_w=0・ΔT_w=0が満たされれば作用がゼロになる。まず実入力の壁フラグ・拘束残差を確認する。
- **［Major］今回の不合格を受けて値2・3を調整し続けること。** 上記一組で帰属の範囲を記録した後は系統を保留し、加速＋point仕上げの評価へ戻る。保留はコード削除や既定変更を意味しない。

呼び出し側の前提への異議:

- **［Minor］値0＋キー5は「スカラーだけ」ではない。** [timeIntegration_d.cu:1009](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1009)で、既に熱伝導の自己側Jacobianを持つ。値0→3の比較は、そこへ主に運動量・仕事のDと近傍K、壁拘束変更を導入する比較として記述すべきである。熱伝導Dを丸ごと新規追加した比較ではない。
- **［Minor］改訂SVDも「原始量の尺度で最弱モードを求めた」わけではない。** [解析コード:83](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/linedump_analyze.py:83)は保存量尺度の縮約行列にSVDを適用し、得たベクトルを原始量へ変換している。δT 81〜83%はそのベクトルの成分表示として扱い、外側反復の不安定固有モードや熱伝導原因の証明には使わない。

不足情報:

対象の `run_0260/0261/0288/0300/0301/0302` の原本、実効YAML、残差CSV、VERDICT、生のライン行列はこの作業ツリーにない。保存済み `linedump_run_0288_analysis.json` は旧版で、改訂した原始量割合も含まない。直接再確認できたhost／CUDA差は旧JSONの最大約 **3.9e−12**までで、値3の **≤5e−13** と改訂SVDの数値はplanの記録に依存する。

ファイル変更・forge起動なし。**plan未反映**。呼び出し側で粘性plan §5.1・§6.5・レビュー記録へ採否と上記の判別条件を反映すること。
