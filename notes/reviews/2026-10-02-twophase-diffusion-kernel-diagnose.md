# codex 諮問 (diagnose): twophase-diffusion-kernel

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-diffusion-kernel.md`](../../notes/reviews/briefs/2026-10-02-twophase-diffusion-kernel.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `ecbb4915` (feature/species-transport)
- **codex**: effort `high`, 12.0 min, rc=0
- **結論**: **カーネル編集を保留し、輸送と相変化が釣り合う状態を使って、P1 と非分割更新の固定点保存をホストで判別する。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**カーネル着手は NO-GO**）

| 論点 | 重大度・採否 | 根拠と推奨 |
|---|---|---|
| ① 補正項の面 z | **Major／風上を採用** | [ハーネス:650](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:650)を再検算し、算術平均では蒸気 −6.1653e−3、風上では +3.2442e−4。補正流束 −Σj⁰ の流出側の z を使う。ただし、風上でも BE 反復は上限 5000 に到達しており、**非線形反復の収束まで保証したとは扱わない**。 |
| ② ソース・θ・commit | **Major／P1 を却下、非分割の全残差方式を採用** | [設計メモ:100](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:100)の P1 は、輸送とソースが釣り合う固定点も動かす（下記）。輸送・ソース・BDF を含む **R_v = R_w − R_g** から前処理し、蒸気・液の制限と commit を一体で定義する。現行の液用 θ をそのまま移植しない。[現行コード:78](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:78)は更新済み総水分を参照し、`avail == 0` では枯渇制限も掛からない。P2 のクランプ依存、P3 の物理時間分割への置換も採らない。 |
| ③ 総水分の FCT | **Major／単純な `qDiff` 拡張を却下** | [passiveFct_d.cuh:339](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/passiveFct_d.cuh:339)は **成分ごとに異なる α** を選ぶ。総水分と液を別々に制限しても、蒸気非負・ΣρY=ρ は保証できない。推奨は、全気相種と液を連成して制限し、総水分を蒸気＋液から戻す設計。Q の実現可能性、補正した拡散流束のエンタルピー、BDF 履歴まで整合させる。完成するまでは二相拡散との dual-time 併用を明示的に拒否する。 |
| ④ float32 許容・停止条件 | **Major／許容の分離は採用、4ε 床の採用は要再検証** | S5 を再実行すると、16ε 床は保存誤差 **7.168e−6**、4ε 床は **4.999e−7 だが反復上限 4000 に到達**。[be_solve:419](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:419)は上限到達でも状態を返し、[判定側:752](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:752)は収束成功を必須にしていない。上限返却を FAIL にし、累積保存 ≤1e−6 は維持する。一方、S1 の入力丸めを含む流束許容は別に定義する。 |
| ⑤ g=0 のビット一致 | **Major／面ごとの旧式分岐を却下、0 step probe を採用** | ①の反例そのものが **両セル g=0**。[設計メモ:39](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:39)の分岐は、そこで負値を生む旧式へ戻る。風上化すると g=0 でも旧式とは代数的に異なる。推奨は、TP 凝縮 ON の新作用素には乾燥面でも風上を適用し、この経路の旧式ビット一致要求を改訂する。凝縮 OFF・Euler、および別変更である物性 G0 は維持し、実際の面流束・対角・エネルギー寄与を atomicAdd 前の probe で比較する。 |

結論: **カーネル編集を保留し、輸送と相変化が釣り合う状態を使って、P1 と非分割更新の固定点保存をホストで判別する。**

第 1 仮説: 最優先の設計リスクは、ソース由来の擬似反復中の総水分変化ではなく、**輸送とソースを異なる前処理で分けて固定点を変えること**。確度: **高**  
　根拠: P1 を θ=1・床なしに縮約すると、R_w=0、R_gᵗʳ=−S、R_g=0 の正しい固定点でも、Δρv=S/P_v−S/P_g となる。前処理対角 P_v=3、P_g=2、S=0.02 では、ホスト検算で **ΔρY_w=−0.0033333333**。非分割の `(R_w−R_g)/P_v` は **0**。[設計メモ:97](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:97)  
　反証条件: 提案する実際の更新写像が、この状態で再正規化・クランプ前から全増分ゼロを保ち、輸送・ソース双方が非ゼロの固定点でも成立すること。

第 2 仮説: float32 の「4ε 床なら保存合格」は、収束未確認の上限返却により合格して見えている可能性がある。確度: **高**。S5 の保存値だけでは「解き切った」を証明できない。

第 3 仮説: dual-time では、独立した FCT 制限が蒸気非負とエンタルピー整合を壊す主要因になり得る。確度: **中**。コード上の独立 α は確認済み、新経路の実測はない。

判別 A/B: **変更するのは蒸気増分の組み立て方だけ**。ρ=V=1、ρY_w=0.3、ρg=0.1、R_wᵗʳ=0、R_gᵗʳ=−0.02、S=+0.02、P_v=3、P_g=2、Q の全残差=0 を与え、float32／float64 で各 **1 擬似更新**。A=P1、B=ソース込みの全残差方式とし、同じ limiter・commit を通す。見る量は補正前の全残差、Δρv、Δρg、ΔρY_w、補正量。  
　→ **A が −0.003333… 動き、B が不動なら P1 を棄却**。  
　→ **B も動くなら、残差変換から limiter・commit までの実装が固定点を保存していない**。B の合格条件は、このゼロ残差入力で全増分・補正量がゼロになること。上記の代数検算は実施済みだが、実際の更新写像による確認は未実施。

やらない方がよいこと: float64 の `ALL PASS` を根拠に CUDA 実装へ進むこと、4ε に替えるだけで停止問題を解決済みにすること、g=0 分岐で旧式一致と新式の非負保証を両立したと扱うこと。

呼び出し側の前提への異議:

- **「θ が固定点でも効けばソース相殺が消える」は不正確**。全残差ゼロに正の θ を掛けてもゼロである。問題は残差の分割、床・射影による見かけの停止、丸めによる停止を区別していないこと。
- S1 の 1 mm 条件は、格納済み float32 入力を float64 で評価しても **最大流束 3.082e−10 > 許容 2.250e−12**。入力量子化と演算誤差を分ける必要がある。ただし、これは累積保存許容を緩める根拠にはならない。
- **Minor**: [plan:129](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:129)の D 指定を、[ハーネス:464](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:464)では二元 D に読み替えている。D_k 固定の作用素試験を残し、混合平均係数の試験と区別するべき。
- S1〜S6 は再実行したが、その `ALL PASS` は **判定対象の単体条件の合格**に限る。CFD の収束・準定常性は評価していない。

不足情報: ソース込み更新写像の参照実装、連成 FCT の制限・履歴・エネルギー仕様、反復上限返却を除外した float32 の受け入れ結果。**plan 未反映**。ファイルは変更しておらず、呼び出し側で §4.2・§6・§5.1 へ反映が必要。
