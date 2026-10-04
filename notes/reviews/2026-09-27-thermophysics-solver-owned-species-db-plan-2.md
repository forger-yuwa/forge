# codex レビュー: thermophysics-solver-owned-species-db (plan)

- **plan**: [`plans/active/thermophysics-solver-owned-species-db.md`](../../plans/active/thermophysics-solver-owned-species-db.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `ef83bdd5` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: 前回レビュー (notes/reviews/2026-09-27-thermophysics-solver-owned-species-db-plan.md) 後の改訂 3 点に集中: §4.4 lump の二元拡散係数を Blanc の法則で作り kinetic 混合平均拡散と併用可にした (前回 M2 の拒否を撤回)、§4.8/§5.1 #10/V7 潜熱は液相を気相と同じ datum でシフトし差で L を作る (前回 M6 の後続 plan 切り出しを撤回)、§4.8 末尾/§10 N2 統一は未決。撤回の妥当性と検証条件の穴を見てほしい
- **extra**: `notes/reviews/2026-09-27-thermophysics-solver-owned-species-db-plan.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

目的と起動時合成の方向は妥当で、前身 plan の未解決部分を扱っています。  
Blanc による併用禁止の撤回と、同一 datum による H2O 潜熱共通化は支持します。ただし、拡散の近似範囲と二相系の検証契約は補正が必要です。  
N2 統一は保留を維持し、今回の実装完了条件から切り離すことを推奨します。

1. **Major — Blanc の係数集約と、補正後の拡散流束の等価性を区別する必要があります。**

   **根拠:** [plan:111](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:111) の式は、lump 外の種 i について、混合平均係数の分母を集約する式として成立します。しかし実装は [thermo_d.cuh:471](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:471) の係数を質量分率勾配に掛け、さらに [speciesTransport_d.cu:271](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:271) で `Jᵢ* = Jᵢ − YᵢΣJ` と補正します。補正項は lump 内各種の流束にも依存します。

   現行の係数・式を double でメモリ内検算しました。va3 の乾き組成、300 K、101325 Pa、`Y_H2O=0.04`、乾き成分の内部比を固定した組成勾配では、次になります。**CFD 実測ではありません。**

   | 比較量 | 値 |
   |---|---:|
   | `full` の H2O 混合平均係数 | 2.289509×10⁻⁵ m²/s |
   | Blanc の H2O–MIXDRY 係数 | 2.289509×10⁻⁵ m²/s |
   | `full` の補正後 H2O 流束を ρ∇Y_H2O で除した値 | 2.278015×10⁻⁵ m²/s |
   | lump の補正後流束との差 | **+0.5045%** |

   したがって、係数の集約は厳密でも、最終的な外部種流束まで厳密に保存するわけではありません。混合平均拡散では勾配・速度基準も区別が必要です。[Cantera の公式仕様](https://www.cantera.org/3.1/python/transport.html)も、それぞれの係数を分けています。

   また、既存 SERN は `EXH` と `AIR` の**複数 lump**を許します（[前身 plan:135](/home/sano/work/forge/plans/accepted/thermophysics-cea-mole-fraction-species.md:135)）。改訂式は「実種 i 対 lump」だけで、lump 同士や構成実種が重なる場合を定義していません。

   **対案:** 併用可は維持し、**固定内部組成を仮定した縮約拡散モデル**と明記してください。lump 同士の係数、共有実種の自己拡散、補正流束、エネルギー流束まで式を確定させます。V4 は係数一致だけでなく、補正後の `Jᵢ*` と `ΣhᵢJᵢ*` を独立参照計算で検査し、`full` との差は近似誤差として記録すべきです。

2. **Major — 同じ datum を引くだけでは、外部 H2O DB と液相の基準整合は保証できません。**

   **根拠:** [plan:133](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:133) の共通シフトは正しいです。120–400 K の10,001点で、評価後の気液 h に同じ定数を引く検算では、参照温度 0／298.15／350 K の L は同一でした。前回の「15.9 MJ/kg」反例は、この方針には当たりません。

   ただし外部 DB は残り、現行 reader は係数を読み込むだけで絶対エンタルピーの基準を検証しません（[speciesDB.cpp:173](/home/sano/work/forge/solver_density_cuda/input/speciesDB.cpp:173)）。

   例えば外部気相 H2O の h に定数 **+100 kJ/kg** が入っていると、気相の sensible h は基準引き後に同じでも、CEA 液相との L は **+100 kJ/kg** 変わります。それでも V7 の「参照温度を変えて L 不変」は通ります。**datum 不変性は、気液データの正しい組合せを証明しません。**

   **対案:** 外部 DB の絶対基準と気液ペアの契約を定義し、整合を確認できない H2O 上書きは凝縮有効時に拒否してください。解決済み記録には液相係数・MW・延長規約・気相と共有するシフトを明記し、液相だけ変更しても restart が拒否される試験を V1 に追加します。V7 には datum 不変性とは別に、既知の気液差を検査する試験が必要です。

3. **Major — V7 は、実際に使われる潜熱経路と二相 EOS を十分に検証していません。**

   **根拠:** 潜熱の利用先は `h2o_latent` だけではありません。

   - [condensationTables_d.cuh:106](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTables_d.cuh:106)：既定の `condFloat: 1` が使う潜熱表の生成。
   - [condensationSourceF_d.cuh:283](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:283)：表と double 退避の切替。
   - [condensationEOS_d.cuh:184](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:184)：`dL/dT` に依存する二相熱容量・音速。
   - [convert_species_field.py:72](/home/sano/work/forge/solver_density_cuda/tools/convert_species_field.py:72)：気液係数を再ハードコードした Python 潜熱。`:119` で保存エネルギーの再構築に使用。

   V6 の経路試験が乾き場なら、Python 側の旧潜熱を残しても通ります。また、現行外挿との差を再計算すると、新 L は120 Kで **−2387.19 J/kg**、150 Kで **−465.92 J/kg**。湿り場の変換には実際に影響します。

   **対案:** #10 に表生成・範囲外退避・面流束・二相反転・Python `CondEOS` の移行を明示してください。V7 を次まで拡張します。

   - double の datum 不変性。
   - `condFloat=0/1` の L、`dL/dT`、区切り両側・表範囲外の一致。
   - `g>0` での保存エネルギーの datum 変換後に、T・P・二相音速が保たれること。
   - 湿り場の種変換・restart 前後の T とエネルギー整合。

   許容差は既存 [test_cond_float.cpp:64](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_float.cpp:64) の L 相対 `2e−6`、微分基準、[同:195](/home/sano/work/forge/solver_density_cuda/tests/unit/test_cond_float.cpp:195) の湿潤反転試験を継承できます。V7 の `1e−12` は double 試験に限定すべきです。

4. **Major — V4 の保存収支と V7 の onset・g に、実行可能な合否条件がありません。**

   **根拠:** [plan:185](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:185) の「許容差は判定ツールの既定」は、ツール名も収支の定義もありません。`check_convergence` は収支監査ではなく、[check_passive_budget.py:2](/home/sano/work/forge/solver_density_cuda/tools/check_passive_budget.py:2) は受動種の補正収支用です。化学種と `ΣhᵢJᵢ` のエネルギー収支を代用できません。

   V7 も変化量の記録だけです。`check_quasisteady` の標準量に onset・g はありません（[check_quasisteady.py:377](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:377)）。

   既存の `case/44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX/` を再判定すると、次でした。

   - 判定区間 `S2_main`、最終 step 23999：**`NOT CONVERGED (stalled/plateau)`**
   - 保存済み系列に対する準定常判定、drift/fluct 各0.2%：**`ALL STEADY`**
   - `g_max`：全7時点で **0**

   この run は V2 の係数参照には使えますが、潜熱の湿潤回帰には使えません。V7 が入口 Tt 分布を指定した方向は正しく、具体的な run と抽出条件の固定が残っています。索引は [case README:783](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:783) です。

   **対案:** 実装前に、V4 の収支式・正規化・数値許容差・実行コマンドを確定してください。node の開境界ではピンによる交換も含めます。V7 は `g>0` の基準 run、onset の定義、g の報告量を固定し、その系列を `--series-csv` で判定します。V5(iii) の Schmidt 試験だけでは、新しい kinetic 拡散の検証を代替できません。

5. **Major — V3 の連続性条件は、現行係数を保持する条件と両立しません。**

   **根拠:** [plan:182](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:182) は h の連続性にも V2 と同じ絶対 `1e−6 J/kg` を要求します。しかし [semiperfect.py:32](/home/sano/work/forge/design/forge_design/gas/semiperfect.py:32) 以降の現行係数を1000 Kで両区間から評価した段差は、次のとおりです。

   - N2：**1.38339×10⁻⁴ J/kg**
   - H2O：**1.89641×10⁻² J/kg**
   - CO2：**2.72470×10⁻³ J/kg**

   同じ datum を両区間に適用しても、この段差は消えません。現行係数を正しく移す実装が失格になります。

   **対案:** 「構成種の区間評価との一致」と「新たな段差を導入しないこと」を別条件にしてください。既存の段差を基準として保持し、datum 適用による段差の増分を検査します。連続化のために係数を修正する作業は今回に混ぜないでください。

6. **Minor — N2 の保留は妥当ですが、将来検証先の `case/28` は誤っています。**

   **根拠:** [plan:166](/home/sano/work/forge/plans/active/thermophysics-solver-owned-species-db.md:166) は空気凝縮の onset を `case/28` 系で確認するとしています。しかし [case/28 README:1](/home/sano/work/forge/case/28.cutler_coaxial_jet/README.md:1) は He／空気同軸ジェットです。空気凝縮の承認済み検証は [condensation-air.md:30](/home/sano/work/forge/plans/accepted/condensation-air.md:30) の **case/34 Arthur** です。

   **対案:** 検証先を訂正し、#12 は「判断待ち」として H2O 完了の阻害条件にしないでください。今回は N2 の現行 L・飽和圧・CPG 経路を保持し、共通ディスパッチ変更による無影響だけ確認します。

**推奨は、縮約モデルの契約を明記したうえでの段階移行です。** 実装前の優先順は、①拡散の最終流束と複数 lump の定義、②外部 DB の気液基準契約、③湿潤・float・Python を含む V7、④V4/V7 の定量ゲート、⑤V3 の段差条件修正です。その後、保存場との内容照合 → 共通 DB → 区間可変化・輸送 → runner → H2O 潜熱の順に進めてください。H2O を別 plan に戻す必要はありません。

確認対象は commit `e1e1a26a`。ファイル変更・新規 CFD 実行なし。数値確認は既存系列の再判定とメモリ内検算です。指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
