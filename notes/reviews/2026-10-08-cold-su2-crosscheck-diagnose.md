# codex 諮問 (diagnose): cold-su2-crosscheck

- **brief**: [`notes/reviews/briefs/2026-10-08-cold-su2-crosscheck.md`](../../notes/reviews/briefs/2026-10-08-cold-su2-crosscheck.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-08
- **commit**: `5eb46ef0` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.3 min, rc=0
- **結論**: **最初に、TP→CPG 変換後の実効設定で熱伝導則だけを比較する CFD 0 step の A/B を行い、定 Prandtl 数への復元を確認してください。**
- **extra**: `case/45.isobutane_m6_d155/cold_pair.py`, `case/45.isobutane_m6_d155/run_0012_su2_sst/sst.cfg`, `procedures/su2-cross-check.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**CPG での SU2 照合は採用しますが、現状の事前登録は修正が必要です。特に、記載どおりの TP→CPG 書換えでは熱伝導モデルが一致しません。**

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **CPG への設定変換：要修正** | TP 設定の生成時に `thermCondMethod` を削除しています（[runner_axismach.py:152](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:152)）。省略時は `0`（[solverConfig.cpp:957](/home/sano/work/forge-integ-1005/solver_density_cuda/input/solverConfig.cpp:957)）なので、`viscMethod: 1` に戻すだけでは熱伝導率は定数のままです。**`physProp.thermCondMethod: 1`・`prandtlLam: 0.72` を明示し、TP 専用の `transport` も外す**必要があります。Sutherland の粘性を合わせるだけでは足りません。 |
| **Major** | **「素 SST＝SU2 と同じ形」：要再検証** | [plan §4.3:98](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:98) の `sstOmegaProdFromPk: 0`・`sstSigmaBlend: 0` という対応表は採用できません。SU2 v8.5 の V2003 は制限後の P_k から ω 生産を作り、拡散係数もブレンドします（[SU2 公式ソース](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/include/numerics/turbulent/turb_sources.hpp#L880)）。対応する forge の設定はこの部分では **両方 `1`** です（[ransSource_d.cu:227](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransSource_d.cu:227)、[ransTransport_d.cu:110](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransTransport_d.cu:110)）。入口 k/ω、壁条件、エネルギー中の k の扱いまで実効条件表を作ってください。旧 SU2 cfg の出口圧も `1389 Pa`、TP 準備側は `2237 Pa` なので、そのままの流用は不可です。 |
| **Major** | **Q1：壁・格子を維持した CPG 比較：範囲を限定して採用** | 同じ形状の CPG 問題を両コードで解く比較として妥当です。ただし **TP の冷却応答全体の正しさは証明できません**。定比熱化により h(T)、温度反転、組成・種輸送、混合気の μ・λ、およびそれらと冷却・乱流との相互作用が検証対象から外れます（[solver-settings.md:233](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:233)）。合格時の結論は「この CPG 条件・格子・比較量でコード間差が許容内」。不合格時も forge の誤りと即断せず、条件・モデル・離散化の不一致として扱ってください。 |
| **Major** | **Q2：共通 y_b の δ_loc・θ_r：条件付き採用** | 共通区間の比較量としては妥当ですが、**境界層全体の厚さを代表する保証はありません**。保存された `theta_*.npz` では、冷却 `run_0183` が選んだ y_b は断熱 `run_0181` の **1.235／1.211／1.223 倍**（x=40／70／94）。`theta_diag` は y_b の点値を ρ_e・u_e に使うため、帯の変更は積分範囲と正規化を同時に変えます（[cold_pair.py:575](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair.py:575)）。これは切断誤差の証明ではありませんが、断熱 TP の帯を無条件に外縁とは呼べません。固定した観測量として定義し、CPG 各腕でも外縁の位置・帯への感度を確認してください。δ_E の検証を完了した扱いにはできません。 |
| **Major** | **Q3：R 2%・厚さ/C_f 3%：目標として採用、実績による裏付けは却下** | 旧比較の θ は `(ρu/ρ_eu_e)(1−ρu/ρ_eu_e)` を積分した代理量です（[compare_bl_su2.py:82](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/compare_bl_su2.py:82)）。今回の θ_r は `ρu(u_e−u)` であり、**「前回 θ≤0.4%」をその精度根拠にできません**。3% は暫定の照合目標とし、時間・抽出・写像の不確かさを含めて判定してください。各腕3%から比2%は保証されないので、両条件を残します。冷却の解き方を検証するなら q_w を記録だけに落とさず、既存の [§4.4](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:112) の **5%目標**と熱収支確認も維持すべきです。 |
| **Major** | **Q4：TP 場からの初期化：採用／「収束済み」と現ゲート：却下** | 初期推定に未収束場を使うこと自体は可能ですが、`run_0181/0183` は収束済みではありません。下記の判定が根拠です。同一節点対応で保存量を CPG の EOS に整合させ、TP の `roe` はコピーしないこと。ρ・P を維持すると T が変わるため、壁温ピン後の状態も確認が必要です。初期化経路が違っても比較はできますが、**NaN/RISING なし＋5枚 STEADY だけでは高残差停滞を排除できません**。これは [§5.1 #22](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:204) で既に却下されています。残差上限・保存収支・壁解像を登録し、θ_r・C_f・断熱壁温も含む全比較量の定常性を確認してください。延長は各コードのゲートで決め、同じ反復数には揃えません。 |
| **Major** | **Q5：素 SST と生産 SST の両方：採用／差を dilatation 単独に帰属：却下** | 素 SST の対はコード間比較、生産 SST の対は生産条件への接続に必要です。ただし `(dilat, KL)=(0,0)→(2,1)` は2因子変更です。[ransSource_d.cu:127](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransSource_d.cu:127) と [同:134](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/ransSource_d.cu:134) で別々に作用します。現行4本から言えるのは **「生産 SST 設定一式の感度」**までです。補正を物理的に必要と考えることと、比較用に同じモデルへ揃えることは両立します。 |

観測事実の確認として、保存 NPZ に `check_quasisteady.classify` を実行しました。判定窓は各 run の **80000〜100000、5枚全部、drift・osc 各0.1%**です。

- `case/45.isobutane_m6_d155/run_0181_ns_coldmesh_ad/`：δ_E の x=40／70／94 は **`STEADY`**。保存済み収束判定は **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。
- `case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext/`：δ_E は x=40 が **`TRANSIENT-UNSETTLED`**、70／94 が `STEADY`。Q_w は **`DRIFTING`**、窓内変化約2.5%。保存済み収束判定には **`rms_roOmega … RISING (divergent)`** があります。
- 保存された総合判定も **`判定不能 (ゲート不成立)`**です。「冷却効果約25%」は現在の抽出値であり、確定した定常応答ではありません。

根拠成果物は [_band_ab/cold_pair/](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/)、run 索引は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:25) です。

結論: **最初に、TP→CPG 変換後の実効設定で熱伝導則だけを比較する CFD 0 step の A/B を行い、定 Prandtl 数への復元を確認してください。**

第 1 仮説: **記載された変換案をそのまま適用すると、forge は定数熱伝導率、SU2 は定 Prandtl 数となり、冷却比較に物性の不一致が混入する。** 確度: **高**。  
根拠: 上表の設定生成・既定値に加え、[gasProperties_d.cu:97](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/gasProperties_d.cu:97)。同じ式を Python で評価すると、cp=1360・Pr=0.72 の λ は **300 K で0.034886、1500 K で0.099507 W/(m·K)**。テンプレートの定数 `0.0257` を残した場合、それぞれ **26.3%／74.2%低い**値になります。これは物性差であり、流れ場の誤差量ではありません。  
反証条件: 実際の生成済み設定では既に `thermCondMethod: 1` が有効で、対象温度範囲の μ・λ が SU2 指定則と一致していること。

第 2 仮説: **断熱 TP 由来の固定 y_b が、冷却・CPG の外縁を代表せず、δ_loc・θ_r の比較に基準位置依存を残す。** 確度: **中**。帯位置の差は確認済み、積分量への寄与は未確認。

第 3 仮説: **未静定の熱場・乱流場が、現在の TP 冷却比に残っている。** 確度: **中**。上記 VERDICT が根拠で、最終的な比への寄与量は未確認です。

判別 A/B:

- **A**：案どおり変換した実効設定。**B**：そこから **`physProp.thermCondMethod` だけを `1`** に変更。
- **追加 CFD は0 step**。同じ温度点・cp・Sutherland 定数・Pr で、設定が選ぶ μ(T)・λ(T) を評価し、SU2 指定則と照合する。
- **A が不一致、B の最大相対差が1e−6以下なら第1仮説を支持**。両方が一致すれば「今回の変換に熱伝導復元漏れがある」という仮説は棄却。B も不一致なら他の物性設定が未整合であり、CFD 比較は判定不能とする。

やらない方がよいこと: **`viscMethod: 1` だけで物性が揃ったとすること、旧 θ の実績を新 θ_r の保証にすること、不一致を直ちに forge のバグとすること、素／生産 SST の差を dilatation 単独の効果と呼ぶこと。**

呼び出し側の前提への異議: **「収束した TP 解」「2キーだけで同じ SST」「5枚 STEADY なら比較可能」は受け入れません。** 一方、既実施の座標 A/B は保存記録で最大 **0.059756%**、登録した1%未満を確認できました。この座標不整合を主因候補へ戻す理由はありません。

不足情報: 新しい6 run の実効設定・物性変換と初期化の実装、元 run の全場 HDF5・全残差系列、使用バイナリと対応ソースが手元にありません。したがって物性復元漏れは**案を記載どおり実施した場合の指摘**であり、既に新 run で発生したとの認定ではありません。

ファイル変更・forge 起動はしていません。**plan 未反映**。反映先は `plans/active/tooling-nozzle-isothermal-wall-chain.md` の **§4.3、§5.1 #26、§6 の新しい CPG 照合条件**です。
