# codex 諮問 (diagnose): uj-colwise-result

- **brief**: [`notes/reviews/briefs/2026-10-10-uj-colwise-result.md`](../../notes/reviews/briefs/2026-10-10-uj-colwise-result.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `b196885b` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.0 min, rc=0
- **結論**: **元の FAIL を保存し、参照差分の演算精度だけを変える host A/B を新規登録して、微分の不一致か参照計算の誤差かを判別する。製品経路の診断はその後に進める。**
- **extra**: `solver_density_cuda/tools/test_line_visc_jacobian.cpp`, `notes/reviews/briefs/2026-10-10-uj-colwise-output.txt`, `solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **元の FAIL を保存し、参照差分の演算精度だけを変える host A/B を新規登録して、微分の不一致か参照計算の誤差かを判別する。製品経路の診断はその後に進める。**

採否表:

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 「FAIL（桁落ちによる）」と原因まで確定して記録する | **却下・Major** | 再現指標は二つの差分が近いことしか検査していない。[`classify`:98](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:98)。**登録結果は原文の FAIL、解釈は「参照差分の丸め誤差が第一候補。実装への反証としては未確定」**と別記する。 |
| `long double` 調査を正式な合格の代用にする | **却下・Major** | [調査コード:225](/home/sano/work/forge-faceh/solver_density_cuda/tools/probe_line_visc_jacobian_ld.cpp:225) は差分幅係数も `1e-6` → `1e-7` に変更し、一幅しか検査しない。零列36本は `−1` として評価を省略し、最大値更新も NaN を落としうる。**高精度参照による新規試験は採用**するが、今回は `long double` 単独より十分余裕のある多倍長演算を推奨する。 |
| 第2仮説を退け、直ちに製品経路へ進む | **要再検証・Major** | 不一致13列を取り除いても、未解像は rand **56/2000＝2.8%**、edge **2列**、拘束 **70列**。[保存出力:3](/home/sano/work/forge-faceh/notes/reviews/briefs/2026-10-10-uj-colwise-output.txt:3)。[登録:139](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:139) に照らして PASS にはならない。**新規試験で解像を確認してから進む。** |
| 不合格列の証拠は全件保存できている | **却下・Minor** | [出力処理:339](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:339) の `shown < 40` により、合格以外141列のうち101列を省略している。対案は全件の機械可読出力。元の結果は保存し、補足記録として追加する。 |

第 1 仮説: **今回の微分照合の FAIL は、中心差分の丸め誤差を判定器が実装の不一致と分類したもの。** 確度: **高**
  
根拠: rand 標本73・K・列4は、double で誤差 **1.796e−5**、再現指標 **0**。[保存出力:22](/home/sano/work/forge-faceh/notes/reviews/briefs/2026-10-10-uj-colwise-output.txt:22)。一方、`long double` の保存値は **2.032e−10**。[調査出力:740](/home/sano/work/forge-faceh/notes/reviews/briefs/2026-10-10-uj-colwise-ld-probe.txt:740)。また、この列の非零成分は κγ/(c_pρ) であり、[共通関数:144](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:144) はその式を実装している。再集計でも旧不一致13列の高精度側誤差はすべて `1e-6` 未満だった。

反証条件: **参照差分の数値誤差が判定許容より十分小さいことを確認したうえで、同じ列に両差分幅で `1e-6` 超の不一致が残ること。** これは「丸め誤差だけで説明できる」を反証する。

第 2 仮説: 共通関数または拘束写像に小さい系統誤差が残る。確度: **低・未排除**。交換対称性・零空間・等温壁については既存試験が支持するが、列別の正式判定は未完了。

判別 A/B: **共通関数の J・標本・物性・拘束方向・差分幅を固定し、独立な流束評価と差分計算の演算精度だけを変える。A＝double、B＝100桁程度の多倍長演算。** 全3560列を host で一巡し、時間発展は行わない。

- 差分幅は元登録の `h = 1e-6·max(|q_c|, 1e-3)` と `h/2` を維持する。温度固定方向も元の二幅を使う。
- 非零列は両誤差 `≤1e-6`、再現 `≤1e-7`。零列は元の尺度で両誤差 `≤1e-10`、再現 `≤1e-11`。非有限検査と未解像の扱いも維持する。
- 量子化の指標として、各幅で `η_c = max_r{[ulp(R⁺_r)+ulp(R⁻_r)]/(2h·s_c)}` を記録する。**微分そのものと R の ulp は比較しない。比較対象は差分信号 `2h·s_c`。** 新規試験では `η_c > 1e-8` を未解像とする。ただし、これは最終出力の量子化指標であり、内部演算の誤差上限ではない。内部の相殺も確認できるよう、熱伝導項と仕事項の大きさを併記する。
- **Aで旧不一致が再現し、Bが解像条件込みで合格** → 第1仮説を支持し、対象列について許容を超える実装不一致を退ける。製品経路の照合へ進む。
- **Bでも解像した不一致が残る** → 丸め誤差だけという説明を棄却し、第2仮説を優先する。
- Bも未解像なら判別不能。合格へ読み替えない。

やらない方がよいこと: 元登録の閾値や VERDICT の変更、この FAIL だけを理由にした共通関数・熱伝導 K の修正、今回の host 試験から case/45 の破綻原因を断定すること。

呼び出し側の前提への異議: **「登録上の不一致＝微分実装への反証」という解釈規則は強すぎる。** 参照差分の解像を保証できていない。また「高精度で全列が一致」は、正確には「評価した非零3524列の一幅で許容内」である。製品の物性入力・面重み・K配置・壁行置換は別経路であり、今回の結果では除外できない。[製品側呼び出し:970](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:970)。

不足情報: 高精度参照の二幅での誤差・再現性・零列評価・内部非有限検査が不足している。既存 `tlvj` バイナリの再実行は終了コード1、保存出力と全文一致を確認したが、再ビルドはしていない。forge・収束／準定常判定は今回実行していない。ファイル変更なし、**plan 未反映**。
