# codex 諮問 (diagnose): m55-three-coef

- **brief**: [`notes/reviews/briefs/2026-10-15-m55-three-coef.md`](../../notes/reviews/briefs/2026-10-15-m55-three-coef.md)
- **plan**: [`plans/active/tooling-nozzle-card-local-mesh.md`](../../plans/active/tooling-nozzle-card-local-mesh.md)
- **date**: 2026-10-15
- **commit**: `98ff121f` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 6.7 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-m55-three-coef-diagnose.codex.log` (codex: ok)
- **結論**: Δx重み・端点2拘束の3係数壁を採用し、既存のB0設定とリミッタ固定値を継承した新規NSを1本、2係数B0を対照として検証する。
- **extra**: `case/71.m55_combustion_wt/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0・Major 4・Minor 2）。ファイル変更・forge起動は行っていない。

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 3係数・端点2拘束 | **採用** | 保存済みの`run_0008_ns_b0/judge/extract/delta_r_next.csv.gz`から、提示されたc₀・c₁・δ₀と最大誤差0.5985 mmを再現した。自由な3係数最小二乗より、スロート・出口を固定する案を推奨する。 |
| 最小二乗の重み | **修正して採用／Major 1** | 提示係数は列ごとの等重み。x≥0の3,156点中2,024点がx≤2 mにあり、格子密度が較正を左右する。**x≥0全域を対象にΔx重み**で近似する。x≥2 mへの追加の重み付けは、影響の伝播を確認する前には採らない。 |
| 追加較正の条件 | **現案を却下／Major 2** | [plan:75](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:75)はV5bだけを再較正の条件にしている。出口が通ってV5cだけ落ちる場合を扱えない。**V5b・V5c両方を判定し、再較正可能な形状誤差が残る場合だけ追加1回**とする。 |
| リミッタ基準値 | **(b)を採用／Major 3** | [main.cpp:1792](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:1792)のとおり基準値は作用素を変える。最初の本段の開始時に一度だけ決め、再開・比較・形状更新へ継承する。(a)・(c)への置換を今回の形状変更と同時に行う案は却下する。 |
| B0既定化の範囲 | **NS本段への限定を採用／Major 4** | 雛形2本・NS本段だけでよく、S1・S2はpoint、SERNは対象外。ただし`cfl_main`を削除するだけではEulerが既定CFL 4へ変わる。[runner_axismach.py:1001](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:1001)、[同:1274](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:1274)。Euler用の実効設定を別途保持する。 |
| 数値・解釈の記述 | **訂正／Minor 1・2** | 提示係数のx≥3 mでの最大誤差は、保存CSVでは**0.2422 mm**であり「0.22 mm以内」ではない。またδ_E/δ_Cの非一定性だけではc₁の必要性を証明できない。2係数でもδ_E/δ_C = c + δ₀/δ_Cは変化する。根拠は再較正後も残る分布誤差とする。 |

較正は次の形で固定すると再現しやすい。実寸でC(x)=δ_C(x)、E(x)=δ_E(x)、出口をL、端点値をCₜ・Cₑ・Eₜ・Eₑとする。

- a = (Eₑ − Eₜ)/(Cₑ − Cₜ)
- b(x) = xC(x) − LCₑ(C(x) − Cₜ)/(Cₑ − Cₜ)
- c₁ = Σwᵢbᵢ[Eᵢ − Eₜ − a(Cᵢ − Cₜ)] / Σwᵢbᵢ²
- c₀ = a − c₁LCₑ/(Cₑ − Cₜ)、δ₀ = Eₜ − c₀Cₜ

wᵢはx≥0区間の台形積分重みとする。保存CSVでの再計算は以下。

| 重み | c₀ | c₁ [1/m] | δ₀ [mm] | 全域最大誤差 [mm] | x≥3 m最大誤差 [mm] |
|---|---:|---:|---:|---:|---:|
| 提示値：列ごと等重み | 0.674323 | 0.0214886 | −0.305388 | 0.5985 | 0.2422 |
| 推奨：Δx重み | 0.677948 | 0.0210363 | −0.305721 | 0.6118 | 0.1979 |

これは**旧NSへの近似誤差**であり、新しい壁のNSの誤差ではない。今の壁のNSを較正元にするのは局所的な更新として妥当だが、固定点反復の収縮性は未確認である。

x≈1 mの残差0.6 mmは、下流の試験部に影響し得る。局所の厚みだけでなく、その傾きと特性線の到達位置が効くため、単純に「2 mm→0.2 mmだからMach偏差も1/10」とは推定できない。事前評価には既存の[CFD特性線追跡](/home/sano/work/forge-cardmesh/design/forge_design/geometry/cminus_cfd.py:80)と、ガスモデルのPrandtl–Meyer感度を使い、残差の壁角が試験部へ届くかを見る。ただし、これは見積もりでありV5c合格の代用ではない。

リミッタは**(b)で初回だけ決定し、以後固定**を推奨する。今回のM5.5は既に固定値があるため、ρ=0.5810042584、P=186732.713、a=341.8295597を継承する。新規案件ではS2から本段へ渡す状態を使うが、単なる保存場の平均を「forgeと同一」とはしない。forgeは境界条件・壁温拘束の適用後に平均しているため、その状態と`CELLS/volume`を揃える必要がある。[main.cpp:1731](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:1731)、[同:1782](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:1782)

この方式は「一度決めた解析系列の再現性」を保証するもので、任意の初期場から同じ値が出る方式ではない。値・由来・入力場のhashを保存し、`stages="none"`では既存値を継承する。

基準値の影響はゼロとは言えない。[実装](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/limiter_d.cu:394)では、次元付きの平滑化項はq_ref²(Kh/L_ref)³に相当する。ブリーフのrun_0005→0007の値なら、この項は密度で約+1.78%、圧力で+0.475%、速度で+7.89%変わる。**この割合を定常Mach数の変化へ直接換算することはできない**。K=0.05だけを理由に影響を無視しない。

結論: **Δx重み・端点2拘束の3係数壁を採用し、既存のB0設定とリミッタ固定値を継承した新規NSを1本、2係数B0を対照として検証する。**

第1仮説: 2係数壁の途中の過大補正が試験部のMach過大の主因であり、3係数化でV5cを改善できる。確度: **中**。  
　根拠: `case/71.m55_combustion_wt/run_0008_ns_b0`の保存抽出では、途中の過大補正が約2 mm。Δx重みの3係数近似では最大0.612 mm、x≥3 mでは0.198 mmまで減る。一方、実流のV5cは最大偏差0.2247%、半振幅0.07879%で許容外。  
　反証条件: 新NSで実際のδ_in−δ_Eを十分減らしても、定常化後のV5cが改善しない場合、この誤差を主因とする説明を棄却する。

第2仮説: 壁形状変更に伴うδ_Eの変化が大きく、旧NSへの射影では新NSの補正誤差を減らせない。確度: **低・未確認**。

第3仮説: 残る偏差には、壁補正以外の設計場・離散化・参照場の誤差が支配的に寄与する。確度: **低・未分離**。

判別A/B: **対照は既存の2係数B0、介入は3係数への壁補正規則の変更だけ**とする。B0・精度・リミッタ基準値・格子生成規則を固定し、新しい理想形状のEuler参照を作る。新NSはrun_0008の最終場からcross-mesh restartする。

- 長さは既存の窓20,000・保持40,000・上限200,000 stepを維持する。90,000到達は予測に使わず、短い過渡で判定しない。対照の再実行は不要。
- 新NSのδ誤差最大値が旧値から半減以上し、V5cも両量が半減以上なら、**結果A：過大補正が主因という仮説を支持**する。製品合格は別途V5a〜V5dすべてで判定する。
- δ誤差が半減以上してもV5cが同程度改善しなければ、**結果B：単純な過大補正主因説を退け、第3仮説を戻す**。
- δ誤差自体が減らなければ、第2仮説側の結果である。未収束・未定常なら判別不能とする。
- 追加較正は最大1回。ただし、前回更新で誤差が縮小し、残差を同じ3係数空間でさらに減らせる場合だけ許す。V5cだけが落ち、既に形状残差が小さい場合は機械的に再投入しない。

やらない方がよいこと: 3係数化とリミッタ基準値方式の変更を同時に行うこと、近似誤差だけからV5c合格を宣言すること、出口Mだけを見て追加較正を決めること。

呼び出し側の前提への異議: **B0採用とV6合格は別**である。θ_r差は登録許容を超え、δ_E比較は未実施なので「同じ解を確認した」とは書けない。ユーザ決定によるB0採用は維持してよい。また、本段の判定には保存済みの`check_convergence_main.txt`を使う。`--segment`の保存記録はS1・S2を指しており、本段の根拠に使えない。

確認したVERDICTは以下。判定原文は[run_0008のjudge](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0008_ns_b0/judge/check_convergence_main.txt:2)、run索引は[case README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md)。

```text
check_convergence_main.txt:
  本段履歴、last step 92488
  PASS (converged)
  OVERALL: ALL PASS

check_quasisteady（保存CSVから今回再実行）:
  本段出力 2000..90000、45点
  停止規則の7量: OVERALL: ALL STEADY
  V5cの4量:     OVERALL: ALL STEADY
```

不足情報: 新しい壁でのδ_E応答とV5結果は未取得。これは上記A/Bで取得する対象であり、現時点で合格や反復収束を保証しない。既存B0の残差PASSは保存判定原文による確認で、元の残差履歴からの再実行ではない。

証拠: 足りる
