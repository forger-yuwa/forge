# codex 諮問 (diagnose): sern-r8-stage0-thermo

- **brief**: [`notes/reviews/briefs/2026-09-30-sern-r8-stage0-thermo.md`](../../notes/reviews/briefs/2026-09-30-sern-r8-stage0-thermo.md)
- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **date**: 2026-09-30
- **commit**: `4b4b4f0c` (feature/sern-design)
- **codex**: effort `high`, 4.2 min, rc=0
- **結論**: **実 run の保存済み物性を使う演算経路 A/B を1件追加し、(0) の FAIL を保持したまま丸め差として受理できることを記録する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 判断 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| (0) を「丸めによる FAIL」と記録する | **採用** | 下記の独立再計算で、提示された係数差を再現できた。ただし元の `1e-15` 判定は FAIL のまま残す。 |
| (1) の PASS だけで段 (i) を合格にする | **却下・Major** | CFD のノイズ床内という結果だけでは、物性差の原因を特定できない。事前条件は [plan §5.1 R8](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:432)。実 run の物性記録に対する下記の照合を追加し、丸め差として受理する理由を別記する。 |
| 差は最大5 ulp | **訂正・Minor** | 提示された `1560.6275407330086 − 1560.627540733007` は **7 ulp**。ulp 上限を観測値に合わせて後付けする判断も避ける。 |
| 段 (ii) へ直ちに進む | **要再検証・Major** | 6 run の一次記録がこの環境に無く、`GATES PASS` を独立確認できない。物性照合と既存判定記録の確認後に進む。追加 CFD は現時点では不要。 |

結論: **実 run の保存済み物性を使う演算経路 A/B を1件追加し、(0) の FAIL を保持したまま丸め差として受理できることを記録する。**

第 1 仮説: **質量分率・モル分率の往復、再正規化、逐次加算の違いに、係数合成の相殺が重なった丸め差である。** 確度: **高**

根拠:
- 旧処理は `MW = 1 / Σ(Y/MW)`、重み `Y·MW_mix/MW` で係数を合成する。[composition.py:258](/home/sano/work/forge-sern-design/design/forge_design/gas/composition.py:258)
- 新 config のモル分率自体も、質量分率から逆変換した値。C++ はそれを再正規化し、`Σx·MW` と `Σx·a` を逐次加算する。[composition.py:493](/home/sano/work/forge-sern-design/design/forge_design/gas/composition.py:493)、[speciesDB.cpp:442](/home/sano/work/forge-sern-design/solver_density_cuda/input/speciesDB.cpp:442)、[speciesDB.cpp:271](/home/sano/work/forge-sern-design/solver_density_cuda/input/speciesDB.cpp:271)
- ローカルの `m6_on` 入力で旧式と C++ の逐次演算を Python で再計算すると、**旧 `1560.627540733007`／新 `1560.6275407330086`、18係数中10係数の差、最大相対差 `1.0198562351110776e-15`** を再現した。forge は起動していない。
- 同じ入力を80桁で合成した参照値に対し、係数誤差を `Σ|x·a|` で割った最大値は旧・新とも約 `3.50e-16`。問題の係数の相殺倍率 `Σ|x·a| / |Σx·a|` は **17.77**。一律の係数相対誤差 `1e-15` は、この相殺を考慮していない。
- 再計算した物性を50–7000 Kの1007点で評価した EXH の差は、cp 最大相対差 `8.46e-16`、顕内部エネルギー最大絶対差 `7.45e-9 J/kg`。**これは演算再現の結果であり、AWS の実バイナリ・実 run の検証ではない。**

反証条件: 実 run の入力・解決済み記録を使うと係数差を再現できない、または構成種の MW・生係数・温度区切り・datum に別の差が見つかること。

第 2 仮説: `/tmp/r8chk` と実際の B 群で、入力・DB・実行バイナリが異なる。**未確認、確度: 低**。一時 prepare の結果だけでは実 run との結び付きが不足する。

第 3 仮説: 無し。

判別 A/B: **変えるのは合成の演算経路だけ。CFD は0 step。** 実 run の同一構成種データから、A＝旧 Python の換算・合成経路、B＝新 C++ の再正規化・逐次合成経路を再計算し、旧 `species_db.yaml` と実 B run の `resolved_species_*.yaml` に照合する。MW・18係数の差を再現し、温度区切り・datum も確認する。補助確認として cp 相対差、および顕内部エネルギー差を `cp·max(T,298.15)` で割った値が、上記温度点で `1e-12` 以下かを見る。これは元の判定の置換ではない。

→ **保存値の差を再現できれば** H1 を支持し、入力取り違え仮説を退ける。既存6 run の判定記録も確認できれば、例外受理を記録して段 (ii) へ進む。  
→ **再現できなければ**「丸めと説明済み」という判断を撤回し、DB・config・バイナリの対応を調べる。

やらない方がよいこと: `1e-15` を観測値に合わせて緩めること、最大7 ulpをそのまま新しい合格上限にすること、この診断のために6本の CFD を再実行すること。

段 (ii) では、**新バイナリ・同一形状・同一メッシュ・同一 BC の段 (i) B 群**を比較基準にする。4力係数と摩擦寄与、壁別入熱、μ・λの変化、各量の定常性を記録する。`gas.transport` は4種だけでは足りず、m6_on の**全11構成種**に必要である。[composition.py:697](/home/sano/work/forge-sern-design/design/forge_design/gas/composition.py:697) 「合否なし」は変化量の大小についてであり、NaN・収束・定常性の確認を省略する意味にはしない。

呼び出し側の前提への異議: **`GATES PASS` と残差収束は同義ではない。** ゲートは `require_residual_pass` により残差 PASS の必須性を切り替えるため、実際の残差 VERDICT と4量の定常性 VERDICT、判定区間を併記する必要がある。[sern_gates.py:290](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:290)

不足情報: `case/46.sern_design/run_1002_r8A_1`～`run_1007_r8B_3` の実 config、解決済み物性、6本それぞれの残差・定常性判定記録、末尾6000 stepの個別平均値。この環境には対象 run が無いため、(1) PASS は依頼元の報告として扱った。

**plan 未反映**。依頼どおりファイルは変更していない。反映先は `plans/active/tooling-nozzle-sern-chain.md` §5.1 R8・§6.1。
