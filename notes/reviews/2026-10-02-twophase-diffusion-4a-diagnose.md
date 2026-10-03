# codex 諮問 (diagnose): twophase-diffusion-4a

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-diffusion-4a.md`](../../notes/reviews/briefs/2026-10-02-twophase-diffusion-4a.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `499050e6` (feature/species-transport)
- **codex**: effort `high`, 5.8 min, rc=0
- **結論**: **採用予定の全残差変換・limiter・commit をホストの保存試験へ接続し、更新経路だけを替える A/B を先に行う。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**#4 のカーネル着手は、現時点では NO-GO**）

以下の行番号は現 HEAD `499050e6`。ファイル変更・forge 起動・禁止された run 出力の閲覧は行っていない。

| 論点 | 重大度・採否 | 根拠と対案 |
|---|---|---|
| Q1：非分割 B の limiter・commit | **Major／全残差方式は採用、実装移行は要再検証** | S7 を再実行し、P1 の ΔρY_w=−3.3333302e−3、B の全増分・補正ゼロを再現した。ただし、**S3/S5 は B の limiter を通らない**。[更新分岐:503](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:503)は `jacobi_update` / `line_update` を呼び、`vl_limit_commit` は別経路。また [line_update:383](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:383)は直接の蒸気流束残差を使い、提案の「総水分全残差−液全残差」と float32 の演算順も異なる。**採用予定の更新写像そのものを保存試験に接続する**。 |
| Q2：停止則・持ち越し・許容 | **Major／元の BE と許容 ≤1e−6 を維持し、6ε を次の判別試験に固定する** | S5 を1000更新で再測定すると、6ε・持ち越しなしは `mean` **5.6867e−7**、`upwind` **4.3847e−7**、両者とも上限到達0・最大反復5。停止則だけでは不可能、とは言えない。一方、6ε の一般的保証もない。**残差持ち越しと許容緩和は今回採らない**。[carry:459](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:459)は残差全体を次ステップへ加え、格納時の丸めだけを補償していない。 |
| Q3：緩和か種間結合の前処理か | **Major／初版候補には共通増分の緩和0.5を採用** | ブリーフの反例では、緩和なし5000回で未達、0.5で232–247回。種間結合を前処理へ追加する必要性までは示していない。混合平均 D・組成正規化・空間非対角も反復へ影響する。[ハーネス:507](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:507)の緩和は更新後の共通緩和であり、[CUDA:409](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:409)の sweep 内緩和と同一視しない。**新経路の増分に作用する位置を明記し、dual-time 併用拒否は維持する**。 |

結論: **採用予定の全残差変換・limiter・commit をホストの保存試験へ接続し、更新経路だけを替える A/B を先に行う。**

第 1 仮説: 現在の最大の穴は、B の固定点試験と保存試験が異なる更新写像を検証しており、両者の合格を合成して実装可と判断していること。確度: **高**  
　根拠: [ハーネス:383](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:383)、[503](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:503)、[989](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:989)。追加のメモリ上の試験では、S5 に `vl_limit_commit` と `dg_max=1e−4` だけを接続すると、1000更新で **5.8731e−7、上限0、制限作動993回、最大反復83**。これは肯定材料だが、全残差変換・温度制限・相変化ソースを含む B 全体の検証ではない。  
　反証条件: 採用予定の更新写像を通して、制限が実際に作動した状態でも、原方程式の残差基準と保存・非負条件を同時に満たすこと。

第 2 仮説: θ=0 による「更新停止」を、原方程式の固定点と誤認する余地がある。確度: **中**。 [ハーネス:950](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:950)では ρv=0・δρv<0 で共通 θ が0になる。S7 の S=0.01 の例も残差非零のまま停止する。ただし、実際の相変化ソースでこの条件へ到達するかは未確認。

第 3 仮説: 持ち越しによる改善は残差評価精度の変更だけでは説明できないが、「格納丸めの偏りが主因」との特定は未完了。確度: **中**。[run_be_series:702](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py:702)では持ち越しと float64 評価が同時に有効になる。float64 評価だけの S5 再測定は **6.7364e−6、上限0**で不合格だった。

判別 A/B: **S5 の更新経路だけを変更**する。A=既存の `line_update`、B=採用予定の全残差変換＋`vl_limit_commit`。両者とも float32、風上補正、6ε、持ち越しなし、同じ前処理・初期値で **1000物理更新、各更新上限3000回**。B の制限閾値は先に固定する。見る量は再正規化前の全保存量誤差、原方程式残差、上限回数、最小蒸気量、θ、補正量。  
　→ **A合格・B不合格なら更新写像の問題**。**両者合格ならこの保存上の懸念を当該試験で除外**する。両者不合格なら停止則の問題へ戻す。更新量ゼロを残差合格の代用にしない。

やらない方がよいこと: **累積保存を合わせるために持ち越しを実装すること、未確認の種間結合を真因として前処理を大型化すること。**

呼び出し側の前提への異議: **「残差ゼロなら動かない」は「動かない状態では残差ゼロ」の証明ではない。** また「上限到達なし」だけでは原方程式を解いた証拠にならない。現行 S3/S5 の既定補正は `mean` であり、採用予定の `upwind` との区別も必要。

不足情報: 実際の相変化ソースを含む B の反復試験、ρ の更新に伴う項とソース Jacobian の接続、制限作動後の原方程式残差。**plan 未反映**。呼び出し側で §4.2・§5.1 #4a・§6 に反映すること。
