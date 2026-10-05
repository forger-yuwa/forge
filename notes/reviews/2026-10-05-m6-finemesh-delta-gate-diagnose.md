# codex 諮問 (diagnose): m6-finemesh-delta-gate

- **brief**: [`notes/reviews/briefs/2026-10-05-m6-finemesh-delta-gate.md`](../../notes/reviews/briefs/2026-10-05-m6-finemesh-delta-gate.md)
- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **date**: 2026-10-05
- **commit**: `ca705054` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.5 min, rc=0
- **結論**: **②を保留し、`run_0094` の最終場を共通起点として、粗格子を固定した実効 CFL 5 対 1 の A/B を一度行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| P：FAIL のまま②へ進む | **却下・Major** | 同一形状での比較と①後のゲートは [plan:93](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:93) に明記されている。1.077 % は登録値を超えており、壁解像 PASS や②の再較正ではこの不合格を解消できない。②は保留する。 |
| 「粗格子が δ_E を約1 %過小評価した」 | **要再検証・Major** | 観測されたのは、格子・CFL・起動履歴が異なる二つの計算の差。細格子を真値とする根拠はない。まず粗格子上で CFL の寄与を分離する。 |
| Q：第三水準で確認する | **条件付き採用・Major** | CFL を揃えても FAIL なら、生産化への推奨経路は Q。ただし `nj=97` のまま第一セルを半分にする案は修正が必要。[mesh2d.py:67](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:67) の読み取り実行では、下流第一セル比 `1.3e-5→6.5e-6` に対して軸側間隔が **0.08876→0.09614 r_w** に増える。`nj` も増やし、比較対象域で実際に細分される第三水準を事前登録する。 |
| R：ni と第一セル厚の寄与を分離する | **次手として却下** | 寄与分離はできても、生産候補格子の誤差が小さいことは示せない。先に CFL の交絡を除く方が今回のゲート判断に直結する。 |
| δ_E の準定常判定は完了した | **要再検証・Major** | [plan:159](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:159) の `STEADY` は `machmax`・`pmax` に対する判定。δ_E 自身の VERDICT はない。提示された細格子末尾5枚の幅は約 **0.00436 %** と小さいが、粗格子側の同じ抽出定義による時系列と、登録した3窓の判定も必要。 |

結論: **②を保留し、`run_0094` の最終場を共通起点として、粗格子を固定した実効 CFL 5 対 1 の A/B を一度行う。**

第 1 仮説: **δ_E の差の主成分は格子分布変更による空間離散化・抽出値の変化である。ただし粗格子側だけの誤差とは断定できない。** 確度: **中**

根拠: 提示値では、以下の差は 0.007814 r_t、**+1.077374 %** である。

- `case/45.isobutane_m6_d155/run_0094_ns_c2pin_pass2_ext6k/res_6000.h5`：δ_E = 0.725282。
- `case/45.isobutane_m6_d155/run_0103_ns_finemesh_pass_cfl1/res_60000.h5`：δ_E = 0.733096。

細格子の提示された末尾変動はこの差より十分小さい。一方、格子変更は壁第一セルだけでなく、流れ方向分布と断面全体に及ぶ（[mesh2d.py:106](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:106)、[同:129](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:129)）。

反証条件: 粗格子の CFL 1 側が準定常化した後、細格子との差が時間変動・抽出の不確かさの範囲まで消えるなら、「差の主成分は格子変更」という仮説を棄却する。差が1 %以下になるだけなら、主因仮説の棄却ではなく**ゲート不合格への CFL の寄与**が確認されたことになる。

第 2 仮説: **残差プラトー下の到達状態が CFL または残存過渡に依存している。** 確度: **中、未確認**。報告された収束 VERDICT は `NOT CONVERGED` であり、CFL 非依存を仮定してよい状態とは確認できていない。

第 3 仮説: **平滑化・有効断面の選択が格子間の δ_E 差に寄与している。** 確度: **低、未確認**。`extract_and_merge` は単調性ガードで λ を最大1000倍まで変更する（[deltastar_loop.py:85](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:85)）。同じ `band_select="edge"` だけでは実効的に同じ抽出条件とは限らない。現在の読み手が未緩和の `delta_E` 列を使うことは確認できた（[c2pin_solve.py:22](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/c2pin_solve.py:22)）。

判別 A/B:

- **共通起点**：上記 `run_0094/res_6000.h5`。新規の二つの run に `solver_density_cuda/tools/restart_field.py` で保存量を引き継ぐ。格子・物理壁・r_t・k_f・BC・物性・バイナリを固定する。
- **変更は実効 CFL だけ**：A = `time.deltaT.cfl_pseudo: 5`、B = `1`。表示用 `cfl` も揃える。2次・`limiter: 2`・`implicitRelax: 0.7`・`nStepInner` は元 run と同じ値を保持し、段階起動を挟まない。実効キーの根拠は [solver-settings.md:35](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:35)。
- **長さ**：A は12000 step・1000 stepごと保存、B は60000 step・5000 stepごと保存。累積 CFL を比較上の目安として揃えるが、到達状態が等しい証拠にはしない。
- **見る量**：全 `rms_*` と NaN/Inf、同じ Euler 参照・x_F・未緩和列による δ_E 時系列。δ_E は `check_quasisteady.py --series-csv` で判定する。λ・`held`・出口近傍の有効断面・評価点の被覆も照合する。`NOT CONVERGED` をそのまま記録し、未収束の参照床で PASS に読み替えない。
- **事前の解釈**：
  - **CFL 5/1 の δ_E 差が十分小さく、CFL 1 同士でも1 %超が時間変動込みで残る** → CFL が今回の不合格を説明する説を退け、Qへ進む。
  - **CFL 5/1 で再現可能な差が出る** → CFL 非依存の前提を棄却する。粗・細とも CFL 1 の結果で格子比較を判定し直す。
  - **対照Aも動き続ける、δ_E が `DRIFTING`、または比較幅が1 %をまたぐ** → 判定保留。格子誤差へ帰属せず、②へ進まない。

**判別精度には注意が必要**。提示された細格子値を固定すると、粗格子 δ_E が **0.7258376**、現在値から **+0.07661 %** 動くだけで差は1 %になる。したがって「CFL差 ≤0.1 %だから無視」は不適切である。比較する末尾窓・代表値・時間変動幅の計算法を投入前に固定し、差の評価幅全体が1 %以下なら合格、全体が1 %超なら不合格、またぐなら保留とする。

やらない方がよいこと:

- 1.077 %を丸めて1 %扱いする、または不合格を記録するだけで②へ進むこと。
- 壁解像 PASS を格子収束の証拠にすること。
- ②③で変更した物理壁の結果を、旧形状の `run_0094` と比較して格子ゲートを判定すること。
- 現在の `run_finemesh_pin.sh` をそのまま再開すること。まだ `run_0098` を入力とし、①から②へ判定なしで進む（[run_finemesh_pin.sh:13](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_finemesh_pin.sh:13)、[同:28](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_finemesh_pin.sh:28)）。

呼び出し側の前提への異議:

**Q1：P は採らない。CFL の交絡を除いても不合格なら Q を採る。** 第三水準の格子仕様・比較対・抽出条件・判定窓を計算前に固定し、細→第三水準で既存の **≤1 %** を満たすことを②への条件とする。「粗→細の半分以下」は今回の非一様な格子系列では収束次数から導けないため、追加の合否条件にしない。この合格も、選んだ格子対での感度確認であり、真値に対する誤差保証ではない。

**Q2：CFL の確認は必要。** 前回の `run_0101/0102` が示したのは即発散の CFL 依存であり、δ_E の準定常値の CFL 非依存ではない。前回の発散診断をやり直す必要はない。

**Q3：③の1 %条件は削除・緩和せず、①の同一形状比較で満たす前提条件として明文化する。** 元の `run_0094→0103` の **FAIL は履歴に残す**。新しい比較の合否は別に記録する。③では変更後の壁に対する δ_E/δ_C、出口 M、波、壁解像などを独立に判定し、再較正による閉包達成を格子感度の合格に代用しない。

不足情報: 対象 run はローカルにないため、数値と VERDICT は依頼文・plan の報告値として扱った。実効 config とバイナリ hash、粗格子の未緩和 δ_E 時系列、両格子の δ_E 自身の準定常 VERDICT、抽出時の λ・有効断面情報が不足している。

ファイル変更・forge 実行は行っていない。**plan 未反映**。呼び出し側で本提案を `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` §5.1 #11・§6 に事前登録する。
