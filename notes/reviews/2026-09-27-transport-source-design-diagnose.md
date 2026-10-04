# codex 諮問 (diagnose): transport-source-design

- **brief**: [`notes/reviews/briefs/2026-09-27-transport-source-design.md`](../../notes/reviews/briefs/2026-09-27-transport-source-design.md)
- **plan**: [`plans/active/thermophysics-solver-owned-species-db.md`](../../plans/active/thermophysics-solver-owned-species-db.md)
- **date**: 2026-09-27
- **commit**: `78ba4361` (feature/gap-heating-precision)
- **codex**: effort `xhigh`, 5.6 min, rc=0
- **結論**: **承認済みの物性解決順を維持し、次は N2–H2O の固定組成・CFD 0 step の A/B で、CEA の frozen 混合則を採用する仕様を確定する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**修正して進める。最優先の指摘は Major：CEA の混合則が未確定なこと。**

| 対象 | 採否 | 根拠と対案 |
|---|---|---|
| ユーザ指定フィット → CEA → LJ | **採用** | 承認済みの優先順を維持する。ただし μ・λ ごとに係数、単位、適用区間、範囲外規約を解決し、記録・互換性ハッシュに含める。現 plan のハッシュ列挙は LJ までである（[plan:109](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:109)）。 |
| (a) H1：混合則を現状のままにする | **却下〔Major〕** | forge は λ にも粘性用の φ を使う（[thermo_d.cuh:448](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:448)）。CEA は相互作用粘性 ηᵢⱼ から φ を作り、λ には別の ψ を使う（[cea2.f:5613](/home/sano/work/forge/.venv-cea/nasa_cea/cea2.f:5613)）。**CEA の frozen 混合則全体を採用する。** |
| (b) H2O 低温域・H2 | **要再検証** | フィット下限は 373.2 K（[trans.inp:291](/home/sano/work/forge/.venv-cea/nasa_cea/trans.inp:291)）。下記の検算から、単純な LJ 切替も無条件外挿も承認できない。低温域は μ・λ 両方の根拠と接続規約が揃うまで対応済みとしない。IAPWS も名称だけで 200 K の妥当性を保証できない。[公式の適用域説明](https://iapws.org/technical-guidance/release/viscosity) |
| (c) D は LJ を継続 | **採用** | 現行 D は Ω⁽¹˒¹⁾、μ は Ω⁽²˒²⁾ を使う（[thermo_d.cuh:379](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:379)、[456](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:456)）。Sc・Le が変わること自体は不整合ではない。ηᵢⱼ から D を推定する道はあるが、衝突積分比の追加仮定が必要。CEA 内部にも `1.1` の近似がある（[cea2.f:5660](/home/sano/work/forge/.venv-cea/nasa_cea/cea2.f:5660)）。今回は D の変更を混ぜない。 |
| (d) lump を実種へ展開 | **採用** | [plan:129](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:129) の全実種展開を維持し、展開後に CEA 混合を一度だけ行う。同じ実種が複数 lump に現れる場合は分率を合算する。 |
| (e) CEA 本体を合否基準にする | **採用、比較対象を限定** | **同一 T・同一気相組成の μ と frozen λ** を比較する。平衡反応寄与込みの λ は別物（[cea2.f:5635](/home/sano/work/forge/.venv-cea/nasa_cea/cea2.f:5635)、[5738](/home/sano/work/forge/.venv-cea/nasa_cea/cea2.f:5738)）。単成分試験だけでは混合則の誤りを検出できない。 |

結論: **承認済みの物性解決順を維持し、次は N2–H2O の固定組成・CFD 0 step の A/B で、CEA の frozen 混合則を採用する仕様を確定する。**

第 1 仮説: **種別 μ・λ の置換だけでは CEA 基準を満たさず、混合則による数 %〜10 % 級の差が残る。**　確度: **高**

  根拠: 許可された `trans.inp` の係数と両コードの式を Python の double で評価した。**600 K、X_H2O = X_N2 = 0.5、単成分値は両方とも CEA** とした結果：

| 量 | 現行の混合則 | CEA の混合則 | 現行 / CEA − 1 |
|---|---:|---:|---:|
| μ [Pa·s] | 2.57333454e−5 | 2.90293601e−5 | **−11.3541 %** |
| λ [W/(m·K)] | 0.0453553669 | 0.0504103657 | **−10.0277 %** |

  これは **CFD 実測でも FCEA2 実行結果でもなく、コードと係数からの独立検算**。根拠は [trans.inp:290](/home/sano/work/forge/.venv-cea/nasa_cea/trans.inp:290)、[394](/home/sano/work/forge/.venv-cea/nasa_cea/trans.inp:394)、[thermo_d.cuh:409](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:409)、[cea2.f:5613](/home/sano/work/forge/.venv-cea/nasa_cea/cea2.f:5613)。

  対案は、相互作用データがある組では  
  `φᵢⱼ = 2Mⱼμᵢ / [(Mᵢ+Mⱼ)ηᵢⱼ]`、  
  `ψᵢⱼ = φᵢⱼ{1 + 2.41(Mᵢ−Mⱼ)(Mᵢ−0.142Mⱼ)/(Mᵢ+Mⱼ)²}`  
  を使うこと。無い組は CEA と同じ Wilke 相当の推定を使う。**`V3C0` でも ηᵢⱼ は ψ を介して λ に効く。** [NASA RP-1311 §5.2.1](https://ntrs.nasa.gov/api/citations/19950013764/downloads/19950013764.pdf?attachment=true)

  反証条件: 下記の全状態で、FCEA2 の同一組成の値に現行混合則が相対 0.1 % 以内で一致し、上の CEA 式の評価が一致しない場合。本検算の対応付けを撤回して調べ直す。

第 2・第 3 仮説: **追加しない。低温 H2O の適切な代替モデルは未確認。**

判別 A/B: **変えるのは混合則だけ。** 単成分フィット、MW、T、組成を固定し、A＝現行 φ 共用、B＝CEA の ηᵢⱼ・φ・ψ とする。`T = 400, 600, 1000, 2000 K`、`X_H2O = 0, 0.1, 0.5, 1` の **16 状態を各 1 回、CFD 0 step** で評価し、FCEA2 の μ・frozen λ と比較する。CEA 側で組成が変わっていないことも確認する。

  → **A が全点 0.1 % 以内なら**「現行混合則では不足」を棄却する。  
  → **A が失敗し B が全点合格なら**「種別フィットの置換だけで十分」を棄却し、B を採用する。両方失敗なら単位・組成・参照出力の対応が未解決であり、GPU 実装へ進まない。

  実装・検証順は **①この A/B → ②共通 resolver と単成分評価 → ③CEA 混合・lump 展開・CUDA 評価 → ④低温モデル確定後の NS 検証** とする。②では上書き・単位換算・区間端・記録を検査し、③では独立 double 基準への相対誤差を double ≤1e−12、float ≤1e−5 と事前設定する。④の CFD 結論には指定どおり収束・準定常 VERDICT を要求する。

やらない方がよいこと: **H2O の 373.2 K 未満を黙って LJ に切り替えること。** 同温度の検算は CEA μ＝1.22841825e−5、現行 LJ μ＝1.61600880e−5 Pa·s で、**LJ 側が 31.552 % 大きい**。一方、CEA 最低区間の外挿は約 218.64 K 以下で `dμ/dT < 0` になる。連続性だけを理由に外挿を正当化することも避ける。

呼び出し側の前提への異議: **「CEA の種別係数を使う」と「CEA の混合物輸送を再現する」は同義ではない。** また、CEA 自体は最低区間を下へ外挿する実装になっている（[cea2.f:5466](/home/sano/work/forge/.venv-cea/nasa_cea/cea2.f:5466)）。したがって低温で FCEA2 と一致しても、その温度での物理的妥当性の証明にはならない。H2O の差を双極子補正不足に帰する説明も、現時点では未検証のまま扱う。

不足情報: **200–373.2 K の希薄水蒸気について、μ・λ 両方の参照データと許容誤差、接続規約が不足している。** 混合則の判断は先に進められる。ファイル変更・forge 起動・禁止された run データの読み取りは行っていない。**plan 未反映**。呼び出し側で §4.3b・§4.4・§6 V4 に反映すること。
