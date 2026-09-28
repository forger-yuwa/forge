# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-28
- **commit**: `30f71498` (feature/gap-heating-precision)
- **codex**: effort `high`, 6.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m1
- **extra**: `../forge-sern-design/methods/boundary.md`, `../forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
目的の同定と専用境界流束の構成は妥当で、既存の accepted plan に同等機能はありません。  
CPG 試作の主要結果を再現しました。ただし、SST の状態定義と化学種輸送の検証経路を実装前に修正してください。

添付と対象 worktree の plan・試作コードが一致することを確認しました。ホストで再実行した音響反射率は、Δx 5→1.25 mm で M 0 が **0.501→0.0791 %**、M 0.3 が **0.1517→0.0403 %**、高温内部が **0.1547→0.0416 %**。独立掃引 14,787 点も float32/float64 とも非有限 0、真空置換 1,017、HLL 退避 0 を再現しました。これはホスト試作の検証であり、forge 本体・TP・SST の合格を意味しません。

1. **Major — `sstEnergyIncludesK` の右状態が、混合後のスカラー状態と矛盾する。**

   **根拠:** [plan:84](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:84) は `k` も原始変数として混合し、流入スカラーに混合後の値を使う規定です。一方、[plan:107](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:107) は HLLC に渡す右側のエネルギー・圧力を **`k∞`** で作るとしています。

   計画の混合帯反例に `k_i=0.14`、`k∞=0.01` を与えると、混合後は `k_R=0.075` です。§4.3 の HLLC 式を熱力学音速による Davis 波速で直接評価すると、エネルギー流束は `k_R` 使用で **−17.5268**、`k∞` 使用で **−16.9307**。約 **3.4 %** 異なり、質量流束も変わります。後段で同じ質量流束を使っても、この状態の不一致は修復できません。

   **対案:** §4.3 を **L=`k_i`、R=混合・置換後の `k_R`** に統一してください。`E_R*`、`p_R*`、移流する `k_R` を同じ状態から生成し、V0u の混合帯試験では種の面値だけでなく、運動量・全エネルギー・`ρk` の流束も照合すること。

2. **Major — 一次化学種輸送と S3 の両経路を受け入れる検証指定がない。**

   **根拠:** [plan:105](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:105) は両経路の変更を約束しています。しかし §6 には `speciesFaceReconstruction` を切り替える試験指定がありません。実装では [speciesTransport_d.cu:870](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:870) で別カーネルへ分岐し、S3 の [passiveKernels_d.cuh:27](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:27) は現在、境界面で内部組成を使います。一次経路だけ直した実装でも、設定次第で全試験を通過できます。

   **対案:** V2b の流入・流出を `speciesFaceReconstruction: 0 / 2` の双方で実施すると明記してください。S3 は `convMethod≥1`、`speciesImplicitCoupling:1` を指定し、実効設定を保存すること。各経路で流入組成、Σ種流束＝質量流束、残差との離散恒等式を確認してください。HLLC 単体のホスト／CUDA 一致だけでは、この接続部分を検証できません。

3. **Minor — 法線方向の試験結果を、一般の音波に対する非反射性へ拡張している。**

   **根拠:** [methods/boundary.md:150](/home/sano/work/forge-sern-design/methods/boundary.md:150) は「外へ出ていく音波は…反射しない」と断定しています。しかし V2a は法線方向のパルスで、V1(c) は傾いた**一様流の保持**です。斜入射音波の試験ではありません。

   本案の線形極限 `p′−ρc u′ₙ=0` に静止媒質中の斜入射平面波を代入すると、圧力反射係数の絶対値は `(1−cosθ)/(1+cosθ)`。45° で約17 %になります。これは線形化からの評価であり、今回の forge 実測値ではありません。多次元非反射条件では接線方向の波数も関与します。[Giles の原論文](https://people.maths.ox.ac.uk/~gilesm/files/aiaa90.pdf)

   **対案:** 説明を「検証した法線方向の小振幅音波で低反射」に限定し、斜入射・境界に沿う平均流での限界を記載してください。V3 の幅系列による採否は維持し、法線音響試験の合格だけで領域独立性を保証しないこと。

**推奨は、現方式を維持して上記を修正し、V0h を通してから CUDA 実装へ進むことです。** 優先順は ①SST の右状態統一、②一次／S3 の検証条件固定、③非反射性の説明範囲の修正です。既存の ghostless A⁺ を近似前処理として残す方針、TP ホスト検証を先行する順序、V1–V2 と V3 の判定を分離する構成は妥当です。

ファイルは変更していません。指摘・推奨は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 1
