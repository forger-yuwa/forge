# codex 諮問 (diagnose): uj-mp-result

- **brief**: [`notes/reviews/briefs/2026-10-10-uj-mp-result.md`](../../notes/reviews/briefs/2026-10-10-uj-mp-result.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `26d51643` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 5.0 min, rc=0
- **結論**: **U-J の限定付き PASS を採用し、次は同一の凍結状態で、原入力からの host 再組立と CUDA の実際の係数・D/K・壁行を照合する一組だけを行う。**
- **extra**: `notes/reviews/briefs/2026-10-10-uj-mp.txt`, `notes/reviews/2026-10-10-uj-colwise-result-diagnose.md`, `plans/accepted/time_integration-line-viscous-jacobian.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| §6.6 の総合 PASS と「丸めの仮説を支持」 | **採用** | CSV の3560列を再集計し、分類を再計算した。B は全列合格、非有限0、最大誤差 **2.4082e−9**。旧不一致13列も全列合格、最大誤差 **1.0001e−12**。保存された [VERDICT:33](/home/sano/work/forge-faceh/notes/reviews/briefs/2026-10-10-uj-mp.txt:33) と整合する。元の §6.3 の FAIL は保存する。 |
| 共通関数の不整合を、この試験の範囲で退ける | **採用・限定付き** | 対象は指定標本・凍結物性・薄層モデル・host の double。製品の係数生成、CUDA の演算、実際の壁行置換、実残差との整合は未検証。特に壁試験は共通関数の K と拘束方向の照合であり、製品の行置換を実行していない。[試験:266](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:266)、[製品:1144](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1144)。 |
| 既存 dump だけで製品経路を照合できる | **却下・Major** | 状態・`cp`・`gamma` は既にあるが、面係数・実際に読んだ速度・物性・面接続・寄与の内訳が不足する。[dump:1575](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1575)。**原入力、共通関数へ渡した引数、組立前後を追加採取する。** |
| η を追加すれば実残差の方向微分も解像できる | **却下・Major** | η は最終出力の量子化しか測らない。内部には float の γ・音速や輸送物性評価が残る。[熱力学:186](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/dependentVariables_d.cu:186)、[輸送物性:111](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/gasProperties_d.cu:111)。**入力摂動の実現、複数幅での比例性、再評価ノイズも独立に検査する。** |
| η の非有限は必ず判別不能になる | **却下・Minor** | `apply_eta` は `C_UNRES` を `nf` より先に返すため、既に未解像の rand 列では η の非有限が1%枠に紛れ込める。[試験:203](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:203)。`C_NF_J` の優先を維持し、η の非有限を `C_UNRES` より先に扱う。今回の B は全列 raw PASS なので、今回の結論は変わらない。 |
| 旧不一致13列の η の最小値は5.3e−6 | **訂正・Minor** | 標本51・D・列3が **4.991996e−6**。[出力:19](/home/sano/work/forge-faceh/notes/reviews/briefs/2026-10-10-uj-mp.txt:19)。範囲を **5.0e−6〜2.9e−4** に直す。全列が閾値を超えるという結論は維持。 |

結論: **U-J の限定付き PASS を採用し、次は同一の凍結状態で、原入力からの host 再組立と CUDA の実際の係数・D/K・壁行を照合する一組だけを行う。**

第 1 仮説: **共通関数の前後、すなわち製品の係数生成・精度・配置・拘束処理に、意図した薄層作用素との不整合が残っている。** 確度: **中。ただし実際の不一致は未確認。**

根拠: 製品は `implicitSolvePrecision=0` で `ST=float` を使い、**座標を float に変換してから差を取る**。[精度選択:1746](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1746)、[座標差:952](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:952)。一方、粘性残差は `flow_float` の座標差を使う。[残差:122](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/viscousFlux_d.cu:122)。FP64 ビルドでも両者の幾何係数が同じ精度とは限らない。今回の U-J は β・κ・法線を直接与えるので、この経路を通らない。**これを発散原因と断定する証拠はまだない。**

反証条件: 採取した対象で、原入力から独立に求めた面係数、製品精度を再現した各寄与、K の格納先、拘束前後の行列が、事前登録した丸め許容内で全件整合すること。その対象について第1仮説を退ける。

第 2 仮説: **製品組立は意図どおりだが、薄層近似・凍結物性・省略した勾配項と実残差との相違が、大きな Δτ で問題になる。** 確度: 中。P が実残差の厳密微分でないこと自体は [共通関数:94](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94) に明記されているが、増幅との因果は未確認。

第 3 仮説: **有限振幅、ライン外結合、SST の分離更新が増幅を作る。** 確度: 低・未確認。今回の U-J はこれらを試していない。

判別 A/B: **変える一点は行列を作る経路。A＝CUDA の実組立、B＝同じ凍結状態・同じ設定から独立に組む host 参照。ソルバ設定は変えない。**

- **長さ・対象:** `case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/res_100000.h5` を起点とする新規診断 run で、最初の factor と5 sweepを含む **1 step**。値3・マスク7・ISP0を維持する。既存の5ラインを全節点・両向きで採取し、壁・軸・内部節点を含める。これは収束試験ではない。
- **追加出力:** 面番号、両端番号、`line_prev/next`、生の座標・面積ベクトル・`fx`、両端の `Ux/Uy/Uz`・`vis_lam`・`vis_turb`・`thermCond`、`Prt`、実効フラグと精度。さらに実際に渡した β・κ・法線・f_i・状態、薄層の各寄与、拘束直前と直後の D/K を残す。既存の状態 dump は再利用する。
- **照合を二段に分ける:** まず原入力から係数を独立に計算し、**誤った β を host と CUDA が共有して合格することを防ぐ**。次に製品の `ST=float`、加算順序、格納時の変換を再現して、面寄与から最終行列まで照合する。時間項・値3のスカラー・対流・薄層・軸対称項を区別し、大きな総 D の差だけで判定しない。
- **判定:** 接続番号・格納先・拘束フラグは完全一致。数値は列別に尺度化し、float の丸めと相殺を含む許容を採取前に登録する。double 前提の `1e−12` は使わない。許容を説明できない差は判別不能とする。

**A が B と再現性のある許容外差を持つなら第1仮説を支持し、最初に食い違う段階を直す。全段が整合するなら、その対象の第1仮説を退け、第2仮説の実残差応答へ進む。** 組立の不一致を見つけても、それだけで発散原因まで確定しない。

実残差応答は、この照合の通過後に行う。最小構成は、一つの共通状態・一つの補正方向について **Q の再評価2回、±h・±h/2・±h/4 の6回、実補正後の1回、計9評価**。同じ残差評価経路を使い、次を要求する。

- BC・EOS 処理後にも意図した摂動が入っていること。現在の `jprobe` は起点の一致をρだけで検査しているため、全保存量と拘束を照合する。[`jprobe.py`:63](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/jprobe.py:63)。
- 差分幅間の再現 **≤1%**、再評価ノイズと出力量子化の各比 **≤1e−3**を暫定の事前登録案とし、片側応答も確認する。η が小さくても、内部の量子化による局所的な平坦部を「真の微分」と認定しない。
- 一回の更新の増幅を論じるなら、選んだラインだけでなく**全域の実際の平均流補正**を使う。SST・組成を何で固定したかを明示する。現行 `jprobe` の有限振幅試験は、選択ラインだけへの摂動である。[`jprobe.py`:85](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/jprobe.py:85)。

この9評価で有効な幅が得られる保証はない。得られなければ測定不成立で止め、近似作用素の良否を判定しない。

やらない方がよいこと: **熱伝導 K の削除、共通関数の修正、ISP・CFL・緩和の同時変更、長い発散再現 run の追加、η だけを根拠にした旧方向微分結果の合格化。** 「全部入りを成立させる」という方針には、係数生成から監査する順序が合う。

呼び出し側の前提への異議: **「単体試験が通ったので、残るのは近似手法の問題」とまでは言えない。** 製品経路は未監査である。また「丸めの仮説を支持」は採用するが、13列すべての誤差機構を粘性仕事による相殺に一本化しない。

不足情報: 実 run の HDF5・生の行列 dump はこの作業環境にない。保存された `lvcfh_judge.json` の45ゲート通過と4本の `DIVERGED` は確認したが、生データの再監査はしていない。U-J はコード確認と CSV の独立再集計までで、再ビルド・再実行は未実施。ファイル変更・forge 起動なし。**plan 未反映**。呼び出し側で §6.6 の採否と §5.1 の次項を更新する必要がある。
