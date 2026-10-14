# 諮問ブリーフ: ③ベル (M) の Rao 生成器の実装で見つかった方針の穴

- 日付: 2026-10-11
- 依頼元: 主セッション (Opus)。エスカレーション条件 4 (plan に無い修正の前) と 3 (事前登録の TIC 切り替えの目安に形式上該当)。このセッションの諮問は diagnostician (Fable) — ユーザ指示
- 作業ツリー: `/home/sano/work/forge-thruster` (HEAD `56dc8013`、生成器は WIP で commit 済み)
- 対象: plan `plans/active/tooling-nozzle-bell-two-mode.md` §4.4 (S1 の導出)・§5.1 #5・§6 V1。前回の諮問 `notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md` (構成手順 4 の「|c| > 1 は解なし」を書いたのは前回のあなた)

## 1. 観測事実 (主セッションで試験の再実行を確認: `run_rao_axisym_tests.py` ALL PASS、未実施 6 件)

- 実装: `design/forge_design/gas/isentrope.py` (`Isentrope`)、`design/forge_design/geometry/rao_axisym.py` (11 関数 + `gtilde`・`cplus_closure`・`edge_function`)、`design/tests/run_rao_axisym_tests.py`。
- 通った試験: T1b (縁条件と `theta_lip` の差 3.9e-16 rad、DE の C⁺ 適合式の残差 1.3e-15・符号反転の対照 4.9)、T2 (包絡線の恒等式、CPG 4〜6e-11・EXH 0.8〜4.5e-10)、T1d (DE 1.3e-14、網 4.3e-12)、accessor (定比熱の TP と CPG の差 1.8e-14)。
- **穴 1 (|c| 規則)**: c = Φ₂ cos μ / V は DE 上で cos(θ−μ) に等しい (主セッションも代入で確認) ので c ≤ 1 は恒等的。c = 1 は θ = μ の点で、そこで dM = 0 (M の極値) になり、M でパラメータ化した線が折り返すだけ。r でパラメータ化すると、方程式系 [V cos(θ−μ)/cos μ = Φ₂, ln(rρV² sin²θ tan μ) = ln Φ₃] の Jacobian は θ = μ で det = −(∂F₁/∂M)·2 cot θ ≠ 0 で正則に続く (実装者の計算)。両方のパラメータ化が成り立つ範囲で r と M の DE は x_E・r_E・T が 10 桁一致、θ = μ を越えるときの Φ₂・Φ₃ のずれ ≤ 5e-12 (実装者の scratch、生成器には未実装)。
  - 規則のままだと、円弧 kernel・EXH・R_t 1・ε 50・Hall n 161 で L_k 60〜90 % が INFEASIBLE (中の根なし)、100 % だけ成立 (θ_B 29.76°・θ_max 30.71°・θ_E 5.26°、推力の 2 会計の差 −1.39e-4 @ n 161 / −4.4e-5 @ n 321、閉合は n 倍で 4 倍縮む)。§8-4 の TIC 切り替えの目安「5 本のうち 2 本以上で根が括れない」に形式上当たる。
  - r パラメータの参考値 (scratch): θ_B 27.2° → L 29.54・θ_E 2.2°、30.08° → L 22.06・5.6° (ここまでは θ = μ を越えない)、31.5° → L 19.86・7.2°、32.95° → 18.21・8.6°、34.4° → 16.96・9.7°、35.8° → 16.02・10.7° (越える)。70〜90 % は θ_B 31〜36° で成立。60 % は θ_B ≈ 38° 以上が要るが、その付近に別の枝 (θ_B ≳ 37° で θ_E ≈ 40°・L 6.3〜7.2) が混ざるため未確認。C2 族 (c·R_t² = −0.5) も同じ型。
- **穴 2 (T1a)**: Rao チャートの照合は R_t 0.382 が要る (case/29 `mesh/make_nozzle.py:41-48`: 80 % ベルの digitize、ε 10 は θn 34.0°/θe 11.8°、ε 50 は 38.5°/7.5°。L = K(√ε−1)R_t/tan15°、R_d 0.382 R_t。チャートの γ はどの資産にも記録が無い)。Hall の初期線では R_t 0.382 の前進 MOC が n 161 で前線 25 の壁点 θ+μ ≥ 90° になり通らない、n 81 は網が粗く B の探索域が空 (`HallThroat` は R < 0.5 で警告を出す範囲)。参考: case/40 run_0076 では ε 9・固定長の CFD 最適 θn が 24.5〜28.5° で、digitize の 33.8° より低い。
- **穴 3 (T6)**: n 641 は Hall の R_t ≤ 1 のガード (`moc_forward.py:136-139`) が拒否。n 81 は kernel 場が組めない。
- **実装上の選択 (方針の変更ではないが確認を求められている)**: 外側の根を c = dκ/ds で取る (κ_B と 1 対 1)、kernel 場の B の先は同じ曲線で延長して C2 で直線につなぐ、BD は網の論理座標で ξ = ξ_B の C⁻、制御面の会計の「軸 → D」は D を通る C⁺ で閉じる (θ_B 36° の BD は x ≈ 72・M ≈ 10 で軸に着き網が粗い)、Goursat の DE は E の 8 % 先まで、`GasSemiPerfect` の線形補間の ν 表 (差 4.3e-6 rad) は使わない、円弧 100 % では B の後に壁角がさらに約 1° 増える (κ_W(B+) +0.075)。
- 実行時間: 円弧の 1 設計 約 190 s (n 161・321)。C2 は c ごとに場を作り直すので L_k 1 本あたり約 10 倍の見込み。

## 2. 問い (推奨を 1 つずつ)

1. **|c| 規則**: DE を r (または弧長) でパラメータ化し θ = μ で σ を反転して通す、「dr ≤ 0」と「θ ≤ 0」は解なしのまま、という修正で正しいか。T2 (包絡線の恒等式) を折り返しの先でも試すべきか。この修正の後の TIC 切り替えの判定 (§8-4) は、規則の修正前の結果 (4/5 本が根なし) を数えないでよいか。
2. **60 % の別の枝**: θ_B ≳ 37° の別の枝 (θ_E ≈ 40°) の扱い (除外の規則、例: θ_E の上限や L の範囲)。
3. **T1a**: R_t 0.382 で Hall が通らない。(a) CFD ピンの初期線で T1a を回す、(b) 照合の R_t を変えて別の参照を使う、(c) T1a を CPG の古典 Rao の別の定量参照 (文献の輪郭の表など) に置き換える、から 1 つ。チャートの γ 不明の扱いも。
4. **T6**: n 641 が使えない中で分解能の試験をどう組むか (n 161/241/321 など)。
5. 実装上の選択で問題のあるものがあれば。

## 3. 読んでよいもの

- `plans/active/tooling-nozzle-bell-two-mode.md` (§4.4・§5.1・§6 V1・§8)
- `notes/reviews/2026-10-11-bell-s1-rao-tp-derivation.md` (前回)
- `design/forge_design/geometry/rao_axisym.py` (`rao_line` 540 行付近・623 行の |c| 判定、`solve_thetaB`・`solve_kappaB`)、`design/forge_design/gas/isentrope.py`、`design/tests/run_rao_axisym_tests.py`
- `design/forge_design/geometry/moc_forward.py` (`hall_initial_line` 136 行付近)、`design/forge_design/geometry/transonic.py` (`HallThroat`)
- `case/29.bell_vs_conical/mesh/make_nozzle.py:30-60`、`case/29.bell_vs_conical/README.md`
- 必要なら python で生成器を呼んでよい (1 設計 約 3 分。ファイルは書かない)。禁止: `*.log`, `res_*.h5`, `*.pdf`, `plans/README.md`
