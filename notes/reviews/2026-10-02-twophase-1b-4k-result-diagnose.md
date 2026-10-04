# codex 諮問 (diagnose): twophase-1b-4k-result

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-1b-4k-result.md`](../../notes/reviews/briefs/2026-10-02-twophase-1b-4k-result.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `69dade5d` (feature/species-transport)
- **codex**: effort `high`, 5.0 min, rc=0
- **結論**: **二相拡散ONのまま `condensation.condTwoPhaseRelax` だけを1→0.5に変える診断A/Bを、共通初期場から実施する。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **二相拡散ONのまま `condensation.condTwoPhaseRelax` だけを1→0.5に変える診断A/Bを、共通初期場から実施する。**

第 1 仮説: **[Major] 蒸気・液・QをDPLURから点対角更新へ置き換えた系で、緩和不足による振動・過大更新が生じ、θ制限と再正規化を繰り返している。** 確度: **中**

  根拠:
  
  - 水のDPLUR更新結果は更新前の値へ戻され、蒸気・液・Qは別の点対角増分で更新される。設定上の `speciesImplicitCoupling: 1`／`passiveImplicitCoupling: 1` は、この部分のDPLUR更新を保証しない。根拠: [main.cpp:2148](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2148)、[condensationTransport_d.cu:836](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:836)、[twoPhaseDiffusion_d.cuh:156](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:156)。
  - S3の既存安定性実績はDPLURとの組合せに依存する。[solver-settings.md:154](/home/sano/work/forge-species/procedures/solver-settings.md:154)。その実績を新経路へ移せない。
  - ブリーフの `case/16.nozzle_wys/run_0506_twophase_ab_B/` は、初回θ制限なし→約200 stepで残差上昇→長い停滞という経過。実ソースの1D試験にも、緩和1では周期運動、0.5では残差基準達成という先例がある。ただしCFDへの外挿は仮説に留める。[設計メモ:471](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:471)。
  
  反証条件: 残差を支配するセルが継続して「ρv=0かつδρv<0」または「ρg=0かつδρg<0」にあり、更新がθ=0で遮断され、振動ではなく同じ向きの非ゼロ残差を保持していること。その場合、主因は次の境界停止仮説へ移る。

第 2・第 3 仮説:

- **第2 [Major・中] 非負境界でのθ=0停止。** 上記条件では共通θが0となり、水・液だけでなくQの更新も止まる。[twoPhaseDiffusion_d.cuh:171](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:171)。ただし該当セルと残差寄与は未確認。対案は、θの制限理由と同一セルの継続性を診断A/Bで記録すること。
- **第3 [Major・低〜中] 異なる前処理による密度・化学種増分の不整合を、再正規化が持続的に補っている。** 液・Qにも実際に係数を掛けるため、θ=0でも更新全体が不変とは限らない。[speciesTransport_d.cu:2155](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2155)。対案は同じA/Bで補正前後の変化を分離すること。`max|f−1|=8.53e−3`だけでは原因と結果を区別できない。

判別 A/B:

- **両側とも二相拡散ON**。診断A=`condTwoPhaseRelax: 1`、診断B=`0.5`。共通初期場は `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/res_48000.h5`。BC・CFL・S3・物性・バイナリ・制限閾値は同一。各2000 step、新規runへ出力する。
- 全残差履歴、終了時独立監査の**絶対残差と許容比の両方**、更新数で割ったθ制限頻度、`Qcut`、再正規化ゲートを見る。両側共通の診断計測として、末尾200更新のθ=0セルについて、セルID・ρv/ρg・制限前増分・制限理由・Q残差寄与・再正規化前後差を残す。
- **改善の事前基準**: 診断BのQ2/Q1/Q0各独立残差が、各自の初期値と診断A終了値の両方の1/10以下。末尾のθ=0頻度・Qcut・再正規化補正も診断Aの1/10以下で、非有限・負値なし。
  - **A停滞／Bが上記達成** → 第1仮説を支持。「変更不能な恒久停止だけが原因」は却下。ただし更新写像全体の健全性や正式収束までは証明しない。
  - **両側で同じセルの境界停止が残差を支配** → 第2仮説を優先。正の緩和係数では、ゼロの保有量と負方向増分から生じるθ=0を直接解除できない。
  - **両側未達だけでは判別不能**。0.5の失敗を、点対角前処理一般の問題の棄却に使わない。この限界を省いた二者択一はできない。

やらない方がよいこと:

- θに正の下限を入れる、再正規化を外す、受入閾値を緩めること。非負性や更新契約を変えて症状を隠す。
- `passiveImplicitRelax`だけを変えて新経路の緩和試験とすること。直接効くのは `condTwoPhaseRelax`。[condensationTransport_d.cu:853](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:853)。
- **投入済み8000 stepは、事前登録どおり再発・計算長依存を確認する意味がある。** ただし改善がなければ、16000／32000へ機械的に延ばすことを次の診断にしない。現時点で「延長しても絶対に改善しない」とも断定しない。

呼び出し側の前提への異議:

- **[Major] 「θ=0が2.8〜2.9万セル」は読み違い。** 区間内のセル・更新数であり、固有セル数でも恒久停止数でもない。集計とリセットは [condensationTransport_d.cu:802](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:802)、[同:879](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:879)。対案は更新数による正規化とセルIDの追跡。
- **[Major] 監査自己検査の一致は、離散式の一致の証拠であり、収束や反復安定性の証拠ではない。** ブリーフ上、`run_0505_twophase_ab_A/`・`run_0506_twophase_ab_B/` はともに `[twophase-audit] VERDICT: NOT CONVERGED`、再正規化ゲートもFAIL。二相拡散の物理的影響は判別不能。
- `rms_roOmega`の上昇だけでは乱流起因と判定しない。二相更新による場の変化への応答でも説明できる。

不足情報: 実効 `cfl_pseudo`・`scalarCflMax`、θ=0の制限理由と継続セル、残差最大セルとの重なり、8000 stepの判定結果。禁止されたrunファイルは読んでおらず、実測値はブリーフに依拠した。主要4ソースは報告ビルドの `6a7a865d` と現HEADで差分なし。**ファイル変更・forge起動なし。plan未反映であり、本提案の反映先は `plans/active/condensation-two-phase-transport.md` §5.1 #1b。**
