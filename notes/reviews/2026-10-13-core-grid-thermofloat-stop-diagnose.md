# codex 諮問 (diagnose): core-grid-thermofloat-stop

- **brief**: [`notes/reviews/briefs/2026-10-13-core-grid-thermofloat-stop.md`](../../notes/reviews/briefs/2026-10-13-core-grid-thermofloat-stop.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-13
- **commit**: `2e50dcfe` (feature/nozzle-core-grid)
- **codex**: effort `high`, 3.0 min, rc=0
- **結論**: Gref の B0 停滞場から、旧バイナリを共通にして `physProp.thermoFloat` 1／0 を各2本・60,000 step 回す一つの A/B を行い、実行間差と時間変化を超えて準定常の窓平均が変わるかを判定する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: Gref の B0 停滞場から、旧バイナリを共通にして `physProp.thermoFloat` 1／0 を各2本・60,000 step 回す一つの A/B を行い、実行間差と時間変化を超えて準定常の窓平均が変わるかを判定する。

第 1 仮説: 熱物性の評価経路が、B0 の準定常状態の θ_r[65,75] を再現性をもって相対 0.1 % 以上変える。確度: 低（局所的な経路差は確認済みだが、長時間の量への影響は未確認）。
  
根拠:

- [plan §4.95](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1332) では、経路切替で微小補正への非比例成分が GI で 0.058 倍、GB で 0.0048 倍になった。同一 Q0 の縮流部残差にも 0.80 % の差がある。ただし、いずれも停止先の差を直接測った結果ではない。
- [dependentVariables_d.cu:182](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/dependentVariables_d.cu:182) 以降では cp・h の評価、γ の float 除算、音速の `sqrtf` が経路で変わる。さらに [同:209](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/dependentVariables_d.cu:209) で `roe` を再構成している。切替の作用は「Newton の精度」や「初期の微小な種」だけに限定できない。

反証条件: 下記の資格を満たし、両判定窓で主指標の経路差に不確かさを足しても 0.1 % 未満なら、この起点・B0・観測期間における上記の限定仮説を退ける。一般的な無影響の証明にはしない。

第 2 仮説: B0 自体の増幅により、経路間差も同一経路の実行間差と区別できなくなる。確度: 中。[plan §4.77](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1116) の Grefc では、同一入力間の差 N が摂動ありとの差 S の 1.4 倍以上になった。ただし、今回の Gref の θ_r・Q_w のばらつきへの換算は未確認。

第 3 仮説: 経路切替は過渡の速さを変えるが、十分長時間後の準定常値は変えない。確度: 低。[plan §4.27](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:606) では、point の θ_r が行き過ぎて戻る過渡が記録されており、短時間の移動方向では区別できない。

判別 A/B:

- **構成**: 起点は `case/45.isobutane_m6_d155/run_0608_cg_gref_cont_k1/res_60000.h5`。T1a・T1b は `thermoFloat: 1`、T0a・T0b は `0`。各経路の複製を必須にする。point の新規腕は加えない。
- **バイナリ**: 4本とも `~/forge-fgeom7-fp64`、提示された SHA 接頭辞 `7e95e5bb` を選ぶ。今回必要な切替は既存機能なので、過去の B0 と同じバイナリを使う。実行前に完全な SHA と実効設定を記録する。
- **共通条件**: B0 のライン陰解法・方向別 dt 上限50・ITJ 5・CFL 4・緩和0.7、SST・BC・メッシュ・物性 DB を固定する。同一メッシュ用 `restart_field.py` で入力を一度作り、4本へ複製する。保存量9量のビット一致に加え、warm start 用 T を含む入力ファイル全体の同一性も確認する。キーは必ず `physProp` 配下に置き、実効の 0／1 を確認する。
- **長さ**: 各60,000 step を連続実行する。20,000 step は安全確認・途中経過に限り、効果判定に使わない。原則1本ずつ、T1a → T0a → T0b → T1b の順で回す。予算内で資格が得られなければ判定不能とし、自動延長しない。
- **観測**: 主指標は θ_r[65,75]、補助は θ_r[84,94]・Q_w・縮流部残差の保存量別系列。既存と同じ空間抽出を使い、1,000 step 間隔で量を得る。時間窓は **20k〜40k と40k〜60k** に固定する。point の通算160k〜200kの参照値は、百分率の分母と移動方向の記述にだけ使う。

**事前の資格・閾値**:

各量について、D_w を「窓 w の T0 2本平均 − T1 2本平均」の百分率とする。

- F は、両経路・両窓で測った同一経路の複製間差の最大値。
- 各 run・窓の時間変化 u は、前半／後半10kの平均差と、単調系列なら漸近値と窓平均との差の大きい方とする。漸近値を安定して推定できなければ資格不足。
- U は各経路で最大の u を取り、その2経路分を足す。H = F + U を判定の余裕幅とする。**これは2複製による運用上の幅であり、統計的な信頼区間ではない。**

4本とも非有限値・非物理値がなく、対象量が各判定窓で `check_quasisteady.py` の判定ロジックによる **STEADY**、かつ主指標の H ≤ 0.025 % を満たすことを資格にする。既存の θ_r・Q_w 抽出系列に適用し、平板専用の `--quantity theta` で代用しない。全保存量の `check_convergence.py` の VERDICT も別に併記する。

→ **結果 A**: 両窓で D_w の符号が同じ、かつ |D_w| − H ≥ 0.1 %。  
「この起点・B0 では、熱物性経路が主指標の準定常平均を実行間差以上に変えた」と判断する。point に近づくことは条件にしない。

→ **結果 B**: 両窓で |D_w| + H < 0.1 %。  
「この条件・期間では、主指標に0.1 %以上の経路依存を検出しなかった」と判断する。補助量にも同じ判定を行い、差があれば「熱物性経路は効かない」と総括しない。

→ **それ以外**: 判定不能。特に、実行間差が大きい、DRIFTING、窓ごとに符号が変わる場合を結果 B に含めない。

やらない方がよいこと:

- 2万 step で point に0.1 %近づいたことを「停止先が変わった」と読むこと。
- T0 だけ複製し、T1 の実行間差を未測定のまま経路効果を判定すること。
- 古い B0 の継続結果を新しい対照の代わりにすること。
- この試験だけで `thermoFloat` の既定を変更したり、Newton 単独を原因と断定したりすること。

呼び出し側の前提への異議:

1. **Major —「実行間差以内なら効いていない」は却下。**  
   [plan:1116](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1116) の実行間差は、Grefc・短時間・温度差などの結果であり、今回の長時間の積分量の床ではない。**対案**: 両経路を複製し、差を検出できない場合と、0.1 %未満に限定できる場合を分ける。

2. **Major —「point に近づくこと」を経路効果の定義にする案は却下。**  
   [plan:663](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:663) の参照 point は **`NOT CONVERGED (プラトー)`**。[同:664](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:664) の **STEADY** は派生量の資格である。**対案**: 経路効果は方向によらず判定し、point への接近は補足にする。結果 A でも「異なる収束解が存在した」とは言わない。

3. **Major —「固定のかたよりとしては小さい見込み」は要再検証。**  
   [plan:1336](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1336) の 7.6／23.3 は、残差差と一つの補正方向への応答のノルム比である。これだけでは残差差の方向、弱く減衰する方向への投影、θ_r の感度を拘束できない。**対案**: 相対状態変化1e−6を θ_r の差0.6 %の上限評価として使わず、今回の A/B で判定する。

4. **Minor —「GPU の加算順の丸めだけが種」は要再検証。**  
   [plan:1117](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1117) は発生源を未確認としている。**対案**: 「同一入力から生じる実行間差」と記述し、発生源を固定しない。

不足情報: Gref・両経路の長時間の実行間差、`thermoFloat: 0` の実測速度、対象バイナリの完全な SHA と実効設定は未確認。禁止された run データは読まず、既存数値と VERDICT は指定 plan の記録を根拠にした。ファイル変更・forge 実行は行っていない。**plan 未反映**であり、呼び出し側が `plans/active/tooling-nozzle-core-grid.md` §4.95 に続く事前登録と §5.1 に採否を反映する必要がある。
