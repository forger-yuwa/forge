# 非平衡凝縮モデルの後続課題 (carrier 補正・空気凝縮の引き継ぎ)

## メタ

- **area**: `condensation`
- **status**: `draft`
- **related_docs**:
  - `methods/condensation.md` (現在仕様: 核生成・成長・Kantrowitz/Feder 補正・空気 CPG carrier 形・N2 低温物性)
- **related_plans**: 親 [condensation-kantrowitz-carrier.md](../accepted/condensation-kantrowitz-carrier.md), [condensation-air.md](../accepted/condensation-air.md), [condensation-kantrowitz-gamma-twophase-sonic.md](condensation-kantrowitz-gamma-twophase-sonic.md)
- **created**: `2026-09-13`
- **owner**: `sano`

## 1. 目的

2026-09-12/13 の 3 つの凝縮 plan (γ_v・二相音速 / Feder carrier 形 / 空気凝縮) で「後続」に送った課題を 1 か所に集め、
親 plan を accepted へ移した後も残作業の正本が消えないようにする (codex result レビュー 2026-09-13 の要求: 後続を active plan にリンク)。
親 plan の §5.1 で「後続」と書かれた行は本 plan の表を正本とする。

## 2. スコープ

- **やる**: 下表の課題の設計・検証 (各項目は着手時に本 plan を `in_progress` にし、§4 に設計方針を書いてから codex plan レビュー)。
- **やらない**: 親 plan で検証済みの実装の再検証。

## 3. 関連 docs と前提

- Feder carrier 形と σ 感度: [condensation-kantrowitz-carrier.md](../accepted/condensation-kantrowitz-carrier.md) §9 (Wysłouzil node run_0350–0356)。
- 空気 CPG carrier 形と N2 低温物性: [condensation-air.md](../accepted/condensation-air.md) §9 (case/34 run_0014–0035)。
- 実験差の現状: Wysłouzil は mode 1 が壁圧 −4.5 % (21 mm)、mode 2/3 は +15〜+18 % (計算上の onset が 8 mm 上流 = モデル間差)。
  Arthur は壁 cond/dry が 3–4 in で実験より ~9 % 過大 (物性修正で解消せず、原因未切り分け)。

## 4. 設計方針

各課題の着手時に記述する (現時点は未着手)。

## 5. 実装ステップ

着手時に記述する。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 (出典) |
| --- | --- | --- |
| 1 | 成長率 dr/dt の α 感度 | Hertz–Knudsen の質量適応係数 α=1 (上限) を `condAccommodation` キーにし、Wysłouzil mode 3 で α=0.5/0.1 の応答を見る。**onset だけでは J と dr/dt を区別できない**ので観測量は onset・壁圧上昇幅 (Δp/p₀ の最大と勾配)・出口液滴径 (r30, Q1/Q0)・実験の SAXS 液滴径 (Wysłouzil の SAXS 文献を取得) を併記して J/dr/dt の寄与を分離する (carrier plan #5b・§5 (6), ユーザ指摘 2026-09-12) |
| 2 | 分圧スイープ | Wysłouzil の p_v0 0.5 / 0.26 kPa の dry node 場を作り mode 1/2/3 の分圧応答 (実験 Fig.3 との比較) (carrier #5) |
| 3 | `condKantrowitz` 既定値の決定 | **グローバル省略時既定は現行 0 (補正なし)**。これを維持するか 1 または 2/3 に変えるかを #1–#2 と実験一致で判断する。Wysłouzil 参照設定 (case/16 run_0335 系) が 1 を明示していることとグローバル既定は別 (carrier #6) |
| 4 | σ 0.97 系の h0 変動 | case/16 run_0354 の中心線 h0 偏差 (平均 1.200 ± 0.015 kJ/kg, 4 点で上昇傾向) が限界サイクルか drift か: 保存間隔を密にして判定し `check_quasisteady.py` に h0err を統合 (carrier #7) |
| 5 | Arthur 3–4 in の壁圧 ~9 % 過大の原因切り分け | 候補: 核生成 J (CNT×Iland), 成長 dr/dt (Goodheart, α), 壁圧の抽出位置 (壁セル列 vs 静圧孔), Arthur 記号の読み取り。J/dr/dt 各 ×0.5/×2 の感度と Fig.2 再デジタイズ (air #10) |
| 6 | CPG 二相音速 | γ_2φ + 実 Ht + 一般 EOS 固有系を CPG 二相 (N2/空気) へ; 流束 FD 照合 (double/float32) (air #7, sonic plan) |
| 7 | 混合液 (露点線) モデル | O2/N2 理想溶液の露点線、2 成分凝縮 (air #8) |
| 8 | 条件の拡張 | Longshot 級 (M 10–14, Ṗ 小) と Daum & Gyarmathy の複数条件で理論線・実験点との比較 (air #9) |
| 9 | 境界試験の実装経由化・流束収支 | slip ghost/bvar を kernel 経由で通す試験、node/cell の質量・エネルギー流束収支 (air #6b) |
| 10 | `check_quasisteady.py` 統合 | onset/壁圧比 (dry 参照場を要する case 固有量) を `--quantity` に (air #6c) |
| F-cf1 | node 周期の凝縮モーメント輸送が非保存 (2026-09-13, condensation-float-speedup §5.1 #12b): case/09 `run_0053` (run_0052 + N2 CPG `condEquilibrium 2`, Q0–Q2 はソース 0) で Σ ρQ_n V が 1 step +5.6e-6、20 step +1.1e-3 (融合前後の両バイナリで同値)。処理別収支 (移流 / dual-time 部分反復 / clamp / periodic seam) で原因を切り分ける | 未着手 |
| F-cf2 | cell 凝縮 run の軸 h0 非保存が 0.1 % 規約を超える (2026-09-13, condensation-float-speedup §5.1 #10): case/34 cell 空気 434 J/kg (1.4e-3 of 3e5), N2 cell 482 J/kg (1.6e-3); node は 129 / 272 J/kg で規約内。基準バイナリ側の既存挙動 (float 化での変化は ≤6 J/kg)。cell の二相面エンタルピー/clamp の収支を確認する | 未着手 |
| F-cf3 | 凝縮 ON 3D 発達場の A10G 52.3 ms/step (目標 ≤50, dry 33.8) の残り: `condensation_source_f_d` 9.6 ms と `dependentVariables_d` の湿潤セル分 +4.8 ms は warp 分岐 (湿潤 22.5 %) と二相反転の double 研磨が主で、ブロックサイズは効かない。候補: 湿潤セルだけを別起動 (インデックス圧縮) で反転・ソース評価、src_jac の摂動評価 2 回→1 回、`cond_primitive_multi_d` の 2 起動統合 (1.3 ms) | 未着手 |
| F-cf4 | N2/空気の σ・ρ_l が 45 K で凍結しており、生産域に入っている (2026-09-14, methods/condensation.md §8b 既知の限界 1): case/34 の onset は T=39.7 K。相関を延長すると σ +9 % (39.2 K)、onset 点で J が 1.2e-6 倍 → onset が数 mm 下流へ動く見込み。`condSigmaScale` と同形の opt-in にして onset 感度と同時に測る | 未着手 |
| F-cf5 | H2O の潜熱が 373.15 K 超で増加する (§8b 既知の限界 2): h_l のクランプで実質 c_p,l=0 となり dL/dT の符号が反転。L(400)=2.320 MJ/kg (真値 2.183)、L(600)=2.711 (本来は T_c=647 K で 0)。凝縮しない温度域だが g>0 のセルが高温へ飛ぶと潜熱が過大。h_l を臨界点まで伸びる相関へ差し替える | 未着手 |
| F-cf6 | H2O の h_v が二重ソース (§8b 既知の限界 3): EOS は種 DB (定 cp 外挿)、`h2o_latent` は生の多項式。200 K 未満で暗黙の h_l が 150 K で 0.47、120 K で 2.39 kJ/kg ずれる (L の 9e-4)。運転域 200–240 K では実質ゼロ。種 DB 経由に統一するか、差を許容と明記するか決める | 未着手 |
| **F-cf7** | **非平衡ソースの θ 律速 (`dT_max` 1 K/step, `dg_max` 5e-3/step) が `dt_local` に比例するため、定常局所時間刻みの収束解が擬似 CFL に依存する** (2026-09-15, case/44 va3 M4.19 入口 Tt 分布 Euler): 大型ノズル (dt_local 2.6e-5 s) では潜熱 ΔT/step ≈3.4 K > 1 K で θ≈0.25 (内側) / 0.55 (壁 2 ノードは dt_local 半分) となり、成長が内側で 4 倍絞られ「壁第一層だけ液相が速く増える」偽の壁異常と凝縮完了の 3–4 r_t 遅れを生む。A/B: cfl 2 → 1 → 0.5 で θ 0.25 → 0.5 → 1、出口 g 平均 0.437 → 0.574 → 0.584 %、出口 M 4.089 → 4.052 → 4.050 (`run_0127`/`0130`/`0131`, 図 `figs/va3_inletTt_cfl_ab.png`)。Wysłouzil (case/16 run_0335, dt_local 小) は ΔT/step 最大 0.86 K で無影響。**修正 (実装済, plan condensation-source-limiter-steady)**: θ を残差に掛けず更新量の同率クランプに限定 (`condLimiterMode 1`)、蒸発は一様 ṙ 形。平衡緩和形 (`condEquilibrium 1`) はスコープ外で据え置き (輸送との釣り合いが Δτ 依存なのは既知の制約) | **実装済** (plan [condensation-source-limiter-steady](condensation-source-limiter-steady.md), 2026-09-15/16; codex result レビュー 2 回目の残件対応中) |
| F-cf8 | **凝縮モーメント (rog, roQ0..2) に dual-time の物理時間項が無い** (2026-09-15, codex plan レビュー M4 of condensation-source-limiter-steady): BDF 残差・対角・時間レベルシフトは平均流 5 変数と roK/roOmega のみで、モーメントはサブ反復ごとに N/M を現在値へコピーし定常と同じ point-implicit 更新。物理 Δt を小さくしてもモーメントがその物理時間で積分される保証がない。物理履歴・BDF 残差・対角を整備し、物理 Δt 半減とサブ反復数変更で検証する | 未着手 |
| F-cf9 | 旧経路 `condLimiterMode 0` (残差に θ を掛ける・λ スケール蒸発) の削除時期 (2026-09-16, plan condensation-source-limiter-steady §4.2-6): 回帰 (Arthur/Wysłouzil/case/44) が新既定で揃った後、A/B 用途が無ければ削除。RK 陽解法・dual-time は自動降格で旧経路に依存しているので、F-cf8 (物理時間項) と RK の更新クランプ試験を先に済ませる | 未着手 |
| F-cf10 | **低 cfl_pseudo での凝縮 onset 域の過渡減衰が擬似時間比より大きく遅い** (2026-09-16, plan condensation-source-limiter-steady §9 交差 restart `case/44 run_0177`–`0180`): cfl 0.5 系の成長遅れ (g L1 ~1.5e-3) は cfl 2 の restart では 6000–8000 step で消えるが、cfl 0.5 では 1 世代 24000 step あたり ~1e-4 しか減衰しない (擬似時間比 4 倍を大きく超える差)。固定点は Δτ 非依存 (交差 restart で確認済) なので実用上は「凝縮 run は cfl 2 で回す」で足りるが、点陰解法 (block-DPLUR + point-implicit モーメント) の低 cfl 側の減衰率がなぜ落ちるか (モーメントの src_jac/transport_diag と流れの結合、壁ノードの半 dt) は未同定 | 未着手 |

## 6. 検証

各課題の着手時に記述する。判定の共通ルール: 同一バイナリ反復 ≥3 run の全ペア最大をノイズ床にし `diff_res.py --tolfile --factor 2` の exit code で判定、
凝縮 run は `check_convergence.py` (凝縮列込み) + `compare_condfix.py --series` / `onset_analysis.py --series`、実カーネル照合は `verify_theta.py`。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan 免除 | `2026-09-13` | — | — | 起票のみ (draft; 設計方針未記述)。各課題の着手時に §4 を書いて `--stage plan` を回す |

## 7. 影響範囲

着手時に記述する。

## 8. 完了条件

表の各行が個別 plan に分離されるか、実装・検証済みで decision が書かれること。

## 9. 変更ログ

- `2026-09-13` — 起票 (codex result レビュー ② の要求: carrier / air の後続を active plan にリンク)。
- `2026-09-16` — F-cf9 追加 (旧 condLimiterMode 0 の削除時期)。
- `2026-09-16` — F-cf10 追加 (低 cfl での onset 域過渡の遅い減衰; 固定点は非依存)。
- `2026-09-15` — F-cf8 追加 (凝縮モーメントの dual-time 物理時間項, codex M4)。
- `2026-09-15` — F-cf7 追加: θ 律速の dt_local 依存 (定常解が擬似 CFL に依存) を case/44 入口 Tt 分布 run の CFL A/B で確定。修正方針を記載、実装は未着手。

## 10. 未確定事項

- なし (課題ごとに着手時に書く)。
