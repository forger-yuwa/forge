# codex 諮問 (diagnose): transport-implementation-design

- **brief**: [`notes/reviews/briefs/2026-09-27-transport-implementation-design.md`](../../notes/reviews/briefs/2026-09-27-transport-implementation-design.md)
- **plan**: [`plans/active/thermophysics-solver-owned-species-db.md`](../../plans/active/thermophysics-solver-owned-species-db.md)
- **date**: 2026-09-27
- **commit**: `f89c73bb` (feature/gap-heating-precision)
- **codex**: effort `xhigh`, 5.1 min, rc=0
- **結論**: GPU 実装へ進む前に、CPU 上で「全 CEA」と「種別選択モデル」の参照値・入力定数・合格条件を分けて固定する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 対象 | 採否 | 根拠・対案 |
|---|---|---|
| **Major：参照値とモデル選択の対応** | **要再検証** | `mixing_ab.py:23–40` は単成分値も全て CEA。採用済みの `iapws_cea` に同じ参照値・許容差は使えない。**全 CEA の再現試験と、選択したモデルの実装一致試験を分ける。** |
| (a) `viscMethod: 3` を追加して旧 `2` を残す | **却下** | [plan:137](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:137) の「現行混合則をやめる」「種別指定必須」に、旧経路を残す例外はない。**`2` を置換し、指定不足は起動時エラー**。旧 run の再現は保存した旧バイナリで行う。 |
| (b) 選択を config、係数を species データに置く | **採用** | `physProp.transport` は実種ごとの必須指定、DB は評価に必要なデータを保持する。lump は展開後に同じ実種の分率を合算する。異種モデル間の ηᵢⱼ は、**両種が `kinetic` なら二元 CE、それ以外は CEA 相互作用データ、欠損時は CEA 剛体球近似**という対称な規約を提案する。これは追加の設計提案であり、混在モデルの精度は未検証。 |
| (c) GPU で double 評価、float 保存 | **初回実装として採用** | 現行は内部 float（[thermo_d.cuh:379](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:379)）。精度変更として扱う。μ・λ を同時評価し、ηᵢⱼ は各非対角組を一度だけ計算、φᵢᵢ＝ψᵢᵢ＝1。セルと壁が同じ評価関数を使うこと。壁にも独立した呼出箇所がある（[wmlesWallModel_d.cu:53](/home/sano/work/forge/solver_density_cuda/cuda_forge/wmlesWallModel_d.cu:53)）。 |
| (d) 新経路だけ記録・ハッシュを拡張 | **採用** | 条件は `viscMethod: 3` ではなく、新しい輸送指定の使用とする。既存 TP の正規化本文・旧記録は維持し、新スキーマで輸送ブロックを追加する。種別係数だけでなく、**解決後の ηᵢⱼ の出所・係数、混合則の版、接続温度、外挿規約、展開行列**を含める。根拠は [plan:109](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:109) の内容照合方針。 |

実装順と合格条件は、次の3段階を推奨する。誤差は特記しない限り相対誤差。

1. **CPU resolver・単成分・独立参照・記録。** 全 CEA の既存16状態だけを FCEA2 と比較し、μ・frozen λ とも **≤0.1%**。合成モデルは同じ MW・係数・接続規約を使う独立評価と **≤1e−12**。500/700 K の左右極限は、値の相対差 **≤1e−12**、無次元勾配 T·d ln f/dT の差 **≤1e−10**。未指定種の拒否漏れ **0件**、既存 TP 本文の変更 **0 byte**。
2. **混合・lump・GPU・壁への接続。** 同じ入力を使い、double 出力 **≤1e−12**、float 保存値 **≤1e−5**。重複実種を含む lump と full、種の列挙順変更、セルと壁の同一 T・組成を検査する。CPU 参照は実装側の混合関数を共有しない。
3. **NS 統合確認。** 別 plan の気相組成処理を固定してから、新しい run で輸送モデル変更だけを評価する。NaN/Inf **0件**、同一数値設定区間の `check_convergence` **PASS**、事前指定した報告量の `check_quasisteady` **STEADY**を要求する。準定常判定は事前に `--tail 0.4 --drift 0.01 --osc 0.01` と固定する。旧結果との不変性は合格条件にしない。

結論: GPU 実装へ進む前に、CPU 上で「全 CEA」と「種別選択モデル」の参照値・入力定数・合格条件を分けて固定する。

第 1 仮説: **Major：素案の参照値をそのまま使うと、正しい H2O 合成モデルを不合格にする。** 確度: 高  
  根拠: [h2o_blend_and_lambda.py:24](/home/sano/work/forge/notes/investigations/2026-09-27-cea-vs-forge-properties/h2o_blend_and_lambda.py:24) の式と、`trans.inp:291,294` の係数を用いた独立スカラー計算では、400 K・純 H2O は次の値になる。CFD の測定値ではない。

| 物性 | CEA | `iapws_cea`（400 K では IAPWS） | CEA 比 |
|---|---:|---:|---:|
| μ [Pa·s] | 1.32788714104e−5 | 1.33545407126e−5 | **+0.569847%** |
| λ [W/(m·K)] | 0.0270414202095 | 0.0264314431570 | **−2.255714%** |

  さらに、[mixing_ab.py:22](/home/sano/work/forge/notes/investigations/2026-09-27-cea-vs-forge-properties/mixing_ab.py:22) の H2O MW は 18.01528 g/mol、現行 DB は 18.0153 g/mol（`forge_species_v1.yaml:98`）。600 K・X_H2O＝0.5 では、この差だけで μ が −8.79e−8、λ が +3.40e−8 相対変化する。**既存 `mixB` を無変更で 1e−12 基準にすることもできない。**  
  反証条件: FCEA2 の0.1%基準が全 CEA 指定専用であり、合成モデルには別参照が用意され、double 比較の MW も一致していると確認できれば、この指摘は解消する。

第 2・第 3 仮説: 追加しない。

判別 A/B: **CFD 0 step**、400 K・純 H2O に固定し、出所だけ **A＝`cea`／B＝`iapws_cea`** に変える。μ・λ と解決済みモデル名を出す。→ A が CEA 値、B が上表の IAPWS 値にそれぞれ **1e−12以内**で合えば、原因は合格基準の取り違え。A/B が同値なら、モデル指定・分岐が効いておらず、合格基準以前の実装不良。

やらない方がよいこと: 合成モデルを FCEA2 に近づけるために係数・接続温度を調整すること。また、case/16 の変化だけで単成分評価・混合則・気相組成処理をまとめて検証したことにしない。

呼び出し側の前提への異議: 「double 評価は現行と同じ」は誤り。全 CEA の16状態による成功も、混在モデルや低温 H2O の検証には拡張できない。今回確認したのはコード・既存比較記録・独立スカラー計算であり、CFD の収束認定はしていない。

不足情報: case/16 の基準 `run_*`、判定区間・VERDICT、採用する報告量、H2O の冪外挿開始温度が未提示。Warnatz・極性補正の正式採用に必要な独立参照も不足している。**ファイル変更なし、plan 未反映**。呼び出し側で §4.3c・§5.1 #5t・§6 に反映すること。
