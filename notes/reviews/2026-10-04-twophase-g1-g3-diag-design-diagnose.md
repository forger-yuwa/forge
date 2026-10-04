# codex 諮問 (diagnose): twophase-g1-g3-diag-design

- **brief**: [`notes/reviews/briefs/2026-10-04-twophase-g1-g3-diag-design.md`](../../notes/reviews/briefs/2026-10-04-twophase-g1-g3-diag-design.md)
- **plan**: [`plans/active/condensation-two-phase-default.md`](../../plans/active/condensation-two-phase-default.md)
- **date**: 2026-10-04
- **commit**: `4d4f0ec4` (feature/species-transport)
- **codex**: effort `high`, 6.4 min, rc=0
- **結論**: **次の一手は、更新へ進まない D1 専用診断で `run_0520` の同一状態に OFF／ON の面作用素を適用し、G1 の符号と絶対誤差を判別すること。**
- **extra**: `plans/active/condensation-two-phase-transport.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表。**D1 は条件付き採用、D2 は現案のままでは却下**。Critical はありません。以下はコードに基づく設計レビューであり、G1・G3 の合格判定ではありません。

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Minor | D1 の `tp_face_flux<double>` | **採用**。既存監査も、格納状態と本番の float 係数を double に上げて再評価している。[speciesTransport_d.cu:2332](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2332)。G1 の演算誤差監査には足り、Python 別実装を追加する必要はない。ただし「式・係数・幾何まで独立に検証した」とは呼ばない。 |
| Major | OFF 面計算の関数切り出しと G0 | **要再検証**。演算順、FMA、レジスタ配置、atomic の実行順が不変とは保証できない。G0 は「ノイズ床×2」の回帰試験であり、ビット不変の証明ではない。[default plan:93](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:93)。atomic 前の面流束・エネルギー・対角を旧演算とビット比較し、通常計算の回帰は別に行う。 |
| Major | D2 の wrapper 後の残差差分 | **却下、記録位置を修正**。OFF は `speciesTransport` 内で移流と拡散、ON は `condensationTransport` 内で移流と二相拡散を実施する。[speciesTransport_d.cu:840](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:840)、[condensationTransport_d.cu:397](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:397)。wrapper 内部に記録点を置き、本番の面流束・ソース加算項も直接保存する。 |
| Major | 残差と更新補正を同じ収支へ加算 | **却下、二つの収支に分離**。面流束と `S·V` は時間当たり、補正 `Δq·V` は更新当たり。さらに更新は DPLUR・緩和・制限を通る。[condensationTransport_d.cu:1287](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:1287)。作用素の収支と更新写像の収支を別々に閉じる。 |
| Major | 0 step で G3 の補正まで取得 | **却下**。更新しなければ、その更新の再正規化・commit 補正は発生しない。一方、組立前の `condensationPrimitive` 自体が保存量をクランプ・射影する。[main.cpp:1837](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1837)、[condensationTransport_d.cu:302](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:302)。0 step は作用素収支まで。更新補正は別の一更新診断として扱う。 |
| Major | 現在の G3 の合否尺度 | **要再検証**。液質量流束差の 10 % を、次元の異なる Q の収支や更新補正へ直接適用できない。また残差を含めて恒等式が閉じても、残差自身が十分小さいとは限らない。[default plan:98](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:98)。成分・単位・比較する液流束差の定義を実装前に固定する。 |

問いへの具体的な回答は次のとおりです。

**1. D1 の double 参照と出力契約**

同じ `tp_face_flux` の double 評価でよい。ただし、入力は float で作った `z`・蒸気分率・流束を昇格するのではなく、**格納した `ρ、ρY、ρg、ρQ` から double で差・除算・正規化をやり直す**。既存の `TpFaceInT<double>` 経路はこの形です。

出力には提案済みの項目に加え、`f、h[k]、L、vis_lam、vis_turb、Sm、up0`、面の向き、処理対象／skip 理由を含める。`TpFaceIn` 全体を保存すれば参照入力を再現できます。

分子蒸気流束は補正後の `j_v = j_v⁰ − z_up,v·Σj⁰` を出し、補正前と補正項も別に残す。`Jv` は乱流分を含むので、そのまま分子流束と呼べません。[twoPhaseDiffusion_d.cuh:86](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:86)

絶対誤差尺度には両端の差し引き前の大きさと、気相内補正の全種和を含める。蒸気では `ρY_w − ρg` の相殺も対象です。液については、例えば

`A_l = |ct·geo| (|ρg₀/ρ₀| + |ρg₁/ρ₁|)`

を用い、`|Jl_float − Jl_double| ≤ 8ε₃₂ A_l` とする。分子蒸気側も対応する被演算項から尺度を明文化し、結果を見て許容を広げない。

なお、A/B の係数が同じという前提には確認が要ります。OFF と ON は面組成の乗算順が異なり、定数 Schmidt 経路では D の分母がそれぞれ `ρ_f` と `ρ_g,f` です。[speciesTransport_d.cu:284](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:284)、[speciesTransport_d.cu:2087](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2087)。**各経路の実係数を保存・照合すること**。ON の係数を OFF へ強制注入して「本番 OFF」とするのは不可です。

**2. OFF の数値不変性**

切り出しだけで不変と判断しない。次を分けて確認します。

- 同じ格納入力、コンパイラ、最適化条件で、旧演算と新演算の **atomic 前**の `Jc、q、diag` を面ごとにビット比較する。
- 診断評価の前後で、本番の保存量・残差・係数配列をバイト比較する。
- 通常計算は G0 型の反復比較で回帰を確認する。この結果は「ノイズ床内」と報告する。

対象には湿った TP の OFF 経路を含める必要があります。不活性クラスだけの G0 では、今回触る湿潤面計算を保証できません。

**3. D2 の記録位置と閉じる式**

最低限、次の位置を区別します。

1. `speciesTransport` 内のゼロ初期化後、移流後、OFF 拡散後。
2. 化学ソース後、最初の `speciesPinResidual` 前後。
3. `condensationTransport` 内の受動種移流後、ON 二相拡散後。
4. `condensationSource` 前後。
5. `passivePinResidual`、二度目の `speciesPinResidual` の各前後。
6. 最終残差確定後。

入口ピンは実際に残差をゼロ化するため、その差を境界拘束の寄与として残す必要があります。[speciesTransport_d.cu:153](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:153)

ただし、**残差差分だけを面流束・ソースの正本にしない**。小さい加算項は float の加算で消えます。また、同じ残差スナップショットから全項を差分して足すだけでは、望遠鏡和が閉じることしか確かめられません。面流束は atomic 前、ソースは実際に加算する `SQn·V、Sg·V` を保存し、残差との差を組立誤差として測ります。ソースの蒸発・凝縮・早期退出・double フォールバックをすべて覆う必要があります。[condensationSourceKernels_d.cuh:374](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:374)、[condensationSourceKernels_d.cuh:500](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:500)

制御体積 Ω を節点 CV の集合、境界流束 F を外向き正とすると、作用素の収支は

`ΣΩ R_final = −F_adv,∂Ω − F_diff,∂Ω + S⁺Ω + S⁻Ω + B_pin,Ω + R_other,Ω + E_assembly`

です。`S⁻ ≤ 0`。既存 `res_*` は面流束と `S·V` を積算済みなので、**残差に再び体積を掛けない**。

蒸気は各段で、double に上げてから

`R_v = R_w − R_l、F_v = F_w − F_l、S_v = S_w − S_l`

と作る。凝縮だけなら `S_w = 0`。float の `res_roYv` を double に上げるだけで済ませない。

更新の収支は別に、

`ΣΩ V(q_after − q_before) = ΣΩ V·δq_update + Σa ΣΩ V·Cq,a + E_update`

とする。`δq_update` は実際の DPLUR・緩和・増分制限を通った更新、`Cq,a` は各補正操作の符号付き格納差です。`δq_update = ΔτR/V` と置いてはいけません。

補正は各操作について `double(q_after) − double(q_before)` を保存し、符号付き積分と絶対量積分を両方出す。再正規化中の負値床、commit、passive floor、上限クランプ、射影、液滴消滅、境界上書きを含める。液滴消滅は κ ゲート外でも収支からは落とせません。ON の一時的な水更新と `twoPhaseHoldWater` による巻き戻しを二重計上しないことも必要です。[main.cpp:2325](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2325)

既存の `[twophase-corr-gate]` は全域集計で、OFF では出力されません。局所・符号付き収支の代用にはなりません。[condensationTransport_d.cu:1005](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:1005)

**4. 0 step と ON 最終場**

G1 と G3 の**作用素部分**は同じ 0 step 診断で取得できます。`nStepOuter: 1` という値だけを「更新なし」の保証にせず、診断専用の分岐で組立・出力後に明示終了させてください。通常経路では組立後に DPLUR と状態更新へ進みます。[main.cpp:2230](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2230)

組立前処理のクランプ・境界上書きも前後差を記録し、**A/B は前処理を一度だけ通した同一状態から分岐**させます。

G1 の主判定入力は事前登録どおり `run_0520` の最終場。G3 は OFF 最終場を OFF、ON 最終場を ON で評価する必要があります。ON 最終場で D1 も出してよいですが、OFF 共通入力の主判定と混ぜない。

既存の `twoPhaseAudit_d_wrapper` を無副作用の診断としてそのまま挿入するのは不可です。モーメント残差をゼロ化してソースだけを再組立する処理があります。[speciesTransport_d.cu:2401](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2401)。専用バッファを使うか、完全な退避・復元が必要です。

**5. 面・境界の扱い**

- `J` は「セル 0 に入る向き」が正で、面積を含む面積分値です。外向きへ符号変換し、面積を二重に掛けない。[twoPhaseDiffusion_d.cuh:16](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:16)
- node 境界半割面の **拡散は skip、移流は残る**。境界面を一括除外すると入口・出口の収支が消えます。[speciesTransport_d.cu:2258](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2258)
- Ω は同じ節点 ID の CV 集合で固定し、片端だけが Ω に属する面を境界として集計する。ghost の体積は含めない。
- 周期は今回の包絡外。未対応なら明示拒否し、集約後の残差を全 member で二重積算しない。
- 壁向き判定には壁距離差と面の向きを使う。上下壁を `y` の符号だけでまとめない。G1 の乱流マスク・壁距離帯・面積重みは実行前に固定する。

結論: **次の一手は、更新へ進まない D1 専用診断で `run_0520` の同一状態に OFF／ON の面作用素を適用し、G1 の符号と絶対誤差を判別すること。**

第 1 仮説: 共通入力上で、ON の気相組成駆動の分子蒸気輸送と液の乱流輸送が、OFF にはない壁法線方向の相別輸送を作る。確度: **中**。  
根拠: OFF は総水分の差を駆動力にする一方、ON は `z_v` と液分率の差を別々に使う。[speciesTransport_d.cu:312](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:312)、[twoPhaseDiffusion_d.cuh:86](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:86)。実際の輸送方向・大きさは未確認。  
反証条件: 固定した壁距離帯で、ON の分子蒸気が外向き・液乱流が壁向きにならない、または差が事前の絶対誤差尺度以下。

第 2・第 3 仮説: 第 2 仮説は、ON/OFF で異なる再正規化が局所差へ寄与していること。確度は中、寄与量は未確認。OFF は液・Q を変更せず、ON は共通係数で変更する。[speciesTransport_d.cu:194](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:194)、[speciesTransport_d.cu:2159](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2159)。第 3 仮説は追加しない。

判別 A/B: `case/16.nozzle_wys/run_0520_twophase_off_16k/res_16000.h5` を共通入力とし、変更は拡散作用素だけ、A=OFF／B=ON、時間更新 0 回。A の総水分流束が G1 の 5 % 条件を満たし、B に予告した輸送方向が誤差尺度を超えて現れれば第 1 仮説を支持。現れなければ第 1 仮説を棄却。double 照合に失敗した場合は診断未成立とする。

やらない方がよいこと: 残差スナップショットの差分だけで G3 合格とすること、更新補正を任意の擬似時間刻みで割って物理ソースへ足すこと、G0 のノイズ床内一致をビット不変と呼ぶこと、G1 成功だけで温度差や正味凝縮量の機構まで確定すること。

呼び出し側の前提への異議: 「0 step で更新補正も得られる」「既存の全域補正ゲートから局所収支を復元できる」「収支恒等式が閉じれば残差の物理的影響も小さい」は受け入れない。恒等式の閉鎖誤差と、残っている残差・拘束・補正の寄与は別に評価する必要がある。

不足情報: 対象の `run_0520`・`run_0521` はこの環境に存在せず、保存場・実効 config・収束／準定常 VERDICT は再確認できていない。G1 の固定面マスクと、G3 の液流束差・成分別許容の定義も必要。**ファイル変更・forge 実行はしていない。plan 未反映**。反映先は `plans/active/condensation-two-phase-default.md` §5.1 #4・§6 と、親 plan §5.1 #4k。
