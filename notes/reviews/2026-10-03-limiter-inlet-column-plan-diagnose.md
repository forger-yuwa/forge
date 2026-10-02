# codex 諮問 (diagnose): limiter-inlet-column-plan

- **brief**: [`notes/reviews/briefs/2026-10-03-limiter-inlet-column-plan.md`](../../notes/reviews/briefs/2026-10-03-limiter-inlet-column-plan.md)
- **plan**: [`plans/active/limiter-inlet-column-oscillation.md`](../../plans/active/limiter-inlet-column-oscillation.md)
- **date**: 2026-10-03
- **commit**: `1826a143` (feature/species-transport)
- **codex**: effort `high`, 5.8 min, rc=0
- **結論**: **参照値を固定した同一 restart 対照を用意し、`time.deltaT.cfl_pseudo: 2 → 1` だけを変える A/B を、毎 step の ψ 測定を組み込んで先に行うことを推奨します。**
- **extra**: `plans/active/limiter-config-simplify.md`, `methods/limiter.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**現状の §4・§6 は要修正です。物理的な真因はまだ確定できません。**

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | H4 の説明：**却下** | block-DPLUR は状態・残差を読み、補正量 `dq` を反復します。状態は sweep 中固定です（[timeIntegration_d.cu:697](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/timeIntegration_d.cu:697)、同ファイル:782、1594）。`nStepInner` の変更は線形解法の反復精度・更新写像の感度試験であり、「内反復中の ψ 更新」の試験ではありません。**外反復更新と動的 ψ の結合**に書き換えてください。 |
| **Major** | ε̂・丸め仮説の除外：**要再検証** | [eps_probe.py:13](/home/sano/work/forge-species/notes/investigations/2026-10-03-twophase-onoff/eps_probe.py:13) は `√(volume/dz)` と `|∇q|h` を使用。実装は平面なら `√volume`、実際の増分は **∇q·(xⱼ−xᵢ)/2** です（[limiter_d.cu:362](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:362)、:379、:492）。大きい勾配ノルムは、小さい面方向増分や極値差の丸めを除外しません。**ψ を決める面の δm・δ±・ε̂²**で検査し直す必要があります。 |
| **Major** | §4.3 の主因／否定基準：**却下** | [対象 plan:89](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:89)。反例として、全節点の残差振幅が 1/10 になっても入口占有率は 71.4 % のままです。現基準はこれを「主因でない」とします。**入口領域の絶対残差ノルムと全域ノルムを主指標**にし、占有率は位置の説明に限定してください。改善しても「介入が効いた」と「真因を特定した」は区別すべきです。 |
| **Major** | I1：測定優先は**採用**、周期・因果判定は**要再検証** | [対象 plan:67](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:67) の 5 step 保存では高周波が折り返され、2 step 周期を一意に識別できません。また状態と勾配は保存時点が異なります（[limiter_d.cu:258](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:258)）。**連続する毎 step 保存**に変更し、ψ と残差の相関だけで因果方向を決めないでください。 |
| **Major** | restart の条件：**要再検証** | 対象 plan:57 は参照値が `auto`、:62 は異なる初期場からの restart。親 plan は A/B・継続で参照値の固定継承を要求しています（[limiter-config-simplify.md:57](/home/sano/work/forge-species/plans/active/limiter-config-simplify.md:57)）。**参照 4 値を元 run の精度で固定し、同じ restart から無変更対照も走らせる**必要があります。丸めたブリーフの数値を転記しないでください。 |
| **Major** | 24000 step と §6：**要再検証** | 対象 plan:85、121–122。提示された L0 の 48000 step の緩和実績に対し、半分の期間で「効かない」は判定できません。**48000 step を最初の評価点とし、下降中なら否定せず延長**してください。全保存量の `check_convergence` に加え、比較する壁圧・壁温について `check_quasisteady` と許容差の事前指定が必要です。 |

結論: **参照値を固定した同一 restart 対照を用意し、`time.deltaT.cfl_pseudo: 2 → 1` だけを変える A/B を、毎 step の ψ 測定を組み込んで先に行うことを推奨します。**

第 1 仮説: 入口近傍の強い制限と外反復の更新が結合し、現在の残差床を維持する振動を作っている。 **確度: 中。ただし run 未確認の暫定仮説。**

  根拠: コードで確認した ψ の決まり方は次のとおりです。

- 近傍集合は「自分＋内部双対面を共有する実ノード」。境界半割面と ghost は除外されます。壁ノードを除外する条件はなく、隣接していれば固定された壁速度も極値に入ります（[limiter_d.cu:319](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:319)、:326、[nodeWallDirichlet_d.cu:29](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:29)）。
- `scaled=1` の評価点は各内部エッジの中点。成分ごとに全対象面の ψ の最小値を採ります。**半 CV だから再構成距離をさらに半分にする処理ではありません**（`limiter_d.cu:362–396`）。
- `inlet_Pressure` は owner の法線 Mach を参照し、指定 Pt/Tt から境界状態を作ります。流入方向は法線方向であり、速度ベクトル全体の単純外挿ではありません（[boundaryCond_d.cu:1059](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/boundaryCond_d.cu:1059)、:1122–1136）。
- **入口境界半割面の流束自体は ψ を使いません。** 境界状態から `mdot` と圧力流束を作ります（[convectiveFlux_boundary_d.inc.cuh:182](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:182)、:192–207）。疑うべき結合は、この境界流束と、入口ノードにつながる内部面の再構成です。

  反証条件: 時間変動が十分落ち着いた比較で CFL を半減しても入口の ψ 変動・絶対残差ノルムが変わらなければ、**「CFL 2 の更新幅が床を大きく維持している」という限定した仮説**を棄却します。時間積分との結合全般を除外するわけではありません。

第 2 仮説: 壁ノードを含む片側近傍の極値差・面方向増分が小さく、float32 の変動で制限を決める面が切り替わる。**確度: 中〜低、未確認。** `|∇q|h` の中央値では除外できません。

第 3 仮説: 入口 BC と内部面離散化の組合せ自体に、更新幅を下げても残る局所的不整合がある。**確度: 低、未確認。** 現時点で BC の欠陥と断定する証拠はありません。

判別 A/B:

- **変更点は `time.deltaT.cfl_pseudo` の 2／1 だけ。** `nStepInner=5`、K、BC、メッシュ、参照値、バイナリを共通にします。IC は双方とも `case/16.nozzle_wys/run_0524_floor_dry_L1/res_48000.h5` の同一メッシュ保存量コピーです。
- **まず双方 48000 step。** 序盤と末尾に連続 200 step の毎 step 出力を残します。24000 step は中間確認に限定します。CFL 1 側が下降中なら延長し、単純な擬似時間の目安として 96000 step まで見ますが、それ自体を十分条件にはしません。
- 入口 2 列の実ノード集合 S を固定し、`R入口 = √Σᵢ∈S(res_ro,i²)`、全域ノルム、各成分の ψ 変動、全 `rms_*` を比較します。床の判定は複数の末尾区間で再現することを条件にします。
- **A：CFL 1 で `R入口` と全域密度残差が双方 1/5 以下、ψ 変動も明瞭に減衰** → 更新幅への強い依存を支持し、「空間離散化だけで決まる不変の床」を棄却。
- **B：十分な期間後も残差・ψ 変動が双方 2 倍以内** → 上記の「CFL 半減で大幅に消える更新幅起因」を棄却。
- 中間・下降継続は判定保留です。**安価な短 run だけで、必ず真因を二分できるという条件は、今回の証拠では満たせません。**

やらない方がよいこと: **入口境界面をさらに 1 次化する、ψ 凍結を追加する、K=0.2 が効かなかっただけで H3 を除外すること。** 許可範囲のコードには入口列だけの ψ 固定機能は見つかりません。`badReconFallback` は非物理再構成、`FORGE_CONTACT_1ST` は組成センサによる介入で、入口列指定の代用にはなりません（[convectiveFlux_slau_d.inc.cuh:108](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:108)、:263）。H1 を直接分離する config だけの試験は確認できません。

呼び出し側の前提への異議:

- **観測と解釈を分ける必要があります。** 「入口 2 列に残差が集中」「2400 step 離れた ψ が違う」は観測の報告ですが、「ψ が揺れ続け、それが床を作る」は未証明です。
- L0/L1 は ε の式だけでなく評価点・増分の計算法も違います。20 倍差を K の問題へ直結させられません（`limiter_d.cu:344–384`）。
- **Minor・要再検証:** 「ψ=0」は正確なゼロか表示丸めかを確認してください。正の ε̂²・有限で正常な入力では新式の分子は正です（[limiterFunctions_d.cuh:43](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiterFunctions_d.cuh:43)）。近ゼロと厳密なゼロでは疑う原因が変わります。

不足情報: `run_0522`〜`run_0525` はローカルに存在せず、残差履歴・HDF5・`CONVERGENCE_VERDICT.txt` を独立検証できませんでした。必要なのは実効 config、参照 4 値、メッシュの実ノード対応、バイナリ識別情報、各 VERDICT と判定区間です。**「dry_L0 ALL PASS」「凝縮と無関係」は今回は独立確認済みとして扱っていません。ファイル変更・forge 起動は行っておらず、提案は plan 未反映です。**
