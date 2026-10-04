# codex 諮問 (diagnose): m6-axis-wave-numerics

- **brief**: [`notes/reviews/briefs/2026-10-04-m6-axis-wave-numerics.md`](../../notes/reviews/briefs/2026-10-04-m6-axis-wave-numerics.md)
- **plan**: [`plans/active/verification-m6-axis-wave-mesh-su2.md`](../../plans/active/verification-m6-axis-wave-mesh-su2.md)
- **date**: 2026-10-04
- **commit**: `bc1c84d1` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.2 min, rc=0
- **結論**: **Mach 抽出を修正したうえで、同一ビルドの粗格子継続と `axis_gap_frac=0.0125` の細分格子を比較する A/B を最優先とし、判定対象を「山への半径方向離散化の寄与」に限定する。**
- **extra**: `case/45.isobutane_m6_d155/README.md`, `plans/accepted/tooling-nozzle-deltastar-core-matched-euler.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | Mach 抽出：**修正を採用** | [axis_wave_compare.py:23](/home/sano/work/forge/case/45.isobutane_m6_d155/axis_wave_compare.py:23) は `gamma` が無いと TP にも γ=1.27354 を使う。実際の `run_0038/res_12000.h5` は `Mach`・`gamma` が無く、`sonic` がある。x≈70 の軸 M はスクリプトで **6.2616**、`sonic` 使用で **6.0161**。**M=√(Ux²+Uy²)/sonic に統一**する。山は修正後も残る。 |
| **Major** | 腕 A の実験：**採用**／「物理」の判定：**却下** | [plan:107](/home/sano/work/forge/plans/active/verification-m6-axis-wave-mesh-su2.md:107) の 0.10 %pt は、修正後の山 0.2325 %pt の **43%**。この範囲の不変性では格子収束を示せない。細分で増大する場合も物理の証明にはならない。対案は、結論をまず**半径方向離散化への感度の有無**に限定する。 |
| **Major** | 収束ゲート：**却下** | [plan:95](/home/sano/work/forge/plans/active/verification-m6-axis-wave-mesh-su2.md:95) の「非 DIVERGED」は未収束も通す。今回再実行した `--from-floor run_0038` は **REFUSED**。参照自身が通常判定を PASS する必要がある（[check_convergence.py:353](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:353)）。**PASS と目的量の STEADY を別々に要求**し、不成立なら診断途中とする。 |
| **Major** | 旧 A0 をそのまま基準にする：**要再検証** | [旧 config:17](/home/sano/work/forge/case/45.isobutane_m6_d155/run_0038_ns_final_rt77p02/solverConfig.yaml:17) は `limiterScaled` 未指定。現行は既定 1、2026-09-20 に変更された（[solverConfig.cpp:619](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:619)）。旧実行バイナリは未確認。**同一バイナリ・実効設定で粗格子も再実行し、それを比較基準にする**。変化を単なる「ノイズ床」として閾値に吸収してはいけない。 |
| **Major** | CPG の腕 B を本件の決着条件にする：**却下** | [plan:65](/home/sano/work/forge/plans/active/verification-m6-axis-wave-mesh-su2.md:65) ではガスに加えて SST の生産項も変わる。波の発生・伝播条件が変わり、固定した x∈[60,80] で山が無くても TP の仮説を棄却できない。**今回の第一手から外す**。壁の局所差し戻しも、消えた場合に証明するのは「壁への感度」であり、物理波と数値応答の区別ではない。 |
| **Major** | 「SU2 の既知の軸バグ」への帰属：**却下** | SU2 は軸上で軸対称ソースをゼロにする実装だが、それだけで今回の 2% 差の原因とは断定できない（[SU2 v8.5 ソース](https://raw.githubusercontent.com/su2code/SU2/v8.5.0/SU2_CFD/src/numerics/flow/flow_sources.cpp)）。対案は**原因未同定のソルバ間差として保持**すること。forge 側の正しさの証拠には使わない。 |
| **Minor** | SU2 の圧力抽出：**修正を採用** | [axis_wave_compare.py:33](/home/sano/work/forge/case/45.isobutane_m6_d155/axis_wave_compare.py:33) は `Energy` から ρk を引いていない。SU2 の原始量復元は k を差し引く（[SU2 v8.5 実装](https://raw.githubusercontent.com/su2code/SU2/v8.5.0/SU2_CFD/src/variables/CNSVariable.cpp)）。修正して再評価する。ただし既存最終出力の軸差は修正後も **2.1386〜2.7178%**で、これが 2% 差の主因ではない。 |

結論: **Mach 抽出を修正したうえで、同一ビルドの粗格子継続と `axis_gap_frac=0.0125` の細分格子を比較する A/B を最優先とし、判定対象を「山への半径方向離散化の寄与」に限定する。**

第 1 仮説: **現在の粗い半径方向離散化が山の振幅を有意に増幅している。** 確度: **中**  
  根拠: `case/45.isobutane_m6_d155/run_0038_ns_final_rt77p02/nozzle.h5` の実測では、x=69.970 r_t で第一軸間隔は **0.07518 r_w = 0.7372 r_t**、Δx=**0.09042 r_t**。ブリーフの 0.68 r_t とは異なる。`sonic` で修正した山の高さは軸上 **0.23247 %pt**。ただし、粗い格子と軸への局在だけでは原因の確定にならない。  
  反証条件: 同一実効設定で細分しても、山の高さ・位置・半径方向形状が後述の時間変動幅を超えて変わらなければ、「この粗さが振幅を大きく増幅する」という仮説を棄却する。

第 2 仮説: 上流壁差に起因する波が軸へ集束している。確度: **中〜低**。親 plan の比較は壁の影響を示唆するが、発生区間と物理的な伝播の帰属は未確認。

判別 A/B: **変える設定は半径方向間隔上限だけ。**

- **A：粗格子**を `restart_field.py` で継続。
- **B：上限 0.0125 r_w**へ細分し、`interp_field.py` で開始。
- 壁点列、x 分布、TP、SST、境界条件、実行バイナリ、実効数値設定を固定する。旧 `run_0038` ではなく、**今回の A を B の基準**にする。
- 各 12000 step を初回予算、1000 step 間隔で保存。12000 step を自動的な終了条件にはしない。全残差の PASS と、`b(0)`・山位置・固定 η のプロファイルの準定常判定を要求する。
- 時間変動の許容幅を事前に **b の末尾変動幅 ≤0.005 %pt**とする。A/B 差が **0.10 %pt 以上の減少**なら「粗格子による増幅」を支持し、**0.02 %pt 以下かつ形状も維持**なら第 1 仮説を棄却する。増大・中間値・未収束なら、物理と命名せず判定保留。

→ **大きく減るなら半径方向離散化の寄与が重要／ほぼ変わらないなら第 1 仮説は後退。ただし後者でも、x 方向誤差や共通の軸処理誤差は除外できない。**

やらない方がよいこと: **CPG/SU2 比較や壁差し戻しを先行して、その結果だけで「物理／数値」を決着させること。** また、ε が大きいときに許容閾値を広げて「物理」に通すこと。

呼び出し側の前提への異議:

- **観測と解釈を分ける。** 山の存在は抽出修正後も確認できる。一方、「上流壁差への物理的集束応答」は未検証の解釈。
- **`ALL STEADY` と収束は別。** 今回の再判定は以下。判定区間は各 run の本段 0〜11999 step。  
  `case/45.isobutane_m6_d155/run_0038_ns_final_rt77p02`：**NOT CONVERGED**。  
  `case/45.isobutane_m6_d155/run_0037_euler_rt77p02`：**PASS**。  
  `run_0038` の修正後 `b(0)` は 4000/8000/12000 step で **0.23103 / 0.23165 / 0.23247 %pt**。`check_quasisteady.py::classify_series` を読み取り専用で適用すると **STEADY**（drift 1%、osc 3%、全4点・末尾3点で判定）。これは山の準定常性を支持するが、残差未収束を解消しない。
- **Euler 誤差は完全には相殺しない。** `max−median` は非線形なので、共通分母でも `b(A)−b(B)` から厳密には消えない（[抽出処理:55](/home/sano/work/forge/case/45.isobutane_m6_d155/axis_wave_compare.py:55)）。山指標と併せて **100·(M_B/M_A−1) の直接差**を見る。
- **「軸側だけ」の範囲が広い。** 現行 `_radial_fracs_capped` の実測では、上限 0.025 は η≲**0.691**、0.0125 は η≲**0.835** の格子を変える（[mesh2d.py:87](/home/sano/work/forge/design/forge_design/meshing/mesh2d.py:87)）。結果は軸境界条件単独ではなく、コア内の波の伝播解像度も含む。

不足情報: 旧 `run_0038` と次の A/B のバイナリ識別情報・実効設定、SU2 の軸量の時系列判定が不足している。**plan 未反映**。本回答は呼び出し側が §4・§6 に反映するための諮問結果であり、ファイルは変更していない。
