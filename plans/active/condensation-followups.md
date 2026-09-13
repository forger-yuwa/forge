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
| F-cf5 | H2O の潜熱が 373.15 K 超で増加する (§8b 既知の限界 2): h_l のクランプで実質 c_p,l=0 となり dL/dT の符号が反転。ずれは 473 K で +27 % (2.460 対 真値 1.940 MJ/kg)、573 K で +89 % (2.657 対 1.404)、T_c=647 K では 0 であるべきところ 2.807。1000 K 以上は [1.5,3.5] MJ/kg クランプの上限に張り付く。**現行ケースは全て 373 K 未満で場には出ていない** (Wyslouzil 全温 ~300 K, Arthur は N2) が、高全温の風洞 (case/44 は T_t=1161 K) で H2O 凝縮を入れると踏む。影響は蒸発の潜熱吸収過大 (過冷却しすぎ) と Newton の dL/dT 符号反転。対処: Watson 型 L(T)=L(373.15)[(T_c−T)/(T_c−373.15)]^0.38 (蒸気表と 473 K 1.5 %, 573 K 1.6 % 一致) で臨界点まで伸ばす。物性表 (373–1200 K) が変わるので case/16 の回帰が要る | 未着手 (2026-09-14 定量化) |
| ~~F-cf6~~ | ~~H2O の h_v が二重ソース~~ **決着 (2026-09-14)**: 種 DB ポインタを持ち回らずに、`h2o_latent` 内の h_v を種 DB と同じ外挿規約 (有効域端の cp 一定で線形外挿, `h2o_gas_h_mass`) に揃えることで解消。係数は元から同一なので、これで同じ温度に対する h_v が 1 通りになる。あわせて 273.15 K 未満の c_p,l を ±0.5 K 差分 (液相フィットの有効域外を踏んでいた) から解析形 4228.27 J/(kg·K) に変更。L の変化は運転域 200–240 K で ≤1.4e-6、表の下端 120 K で 8.4e-4。単体検査 `tests/unit/test_cond_air.cpp` (h) を追加 | 完了 |
| F-cf8 | **上記 2 件の場レベル回帰が未実施** (2026-09-14): 単体 7 本は全 PASS、物性差は解析的に運転域で ≤1.4e-6 (L) と押さえたが、case/16 の凝縮 run での場の比較はしていない (AWS インスタンスが停止中で ssh 不可)。次に AWS が上がったときに `perf_regress.py` / `cmp` でノイズ床比較を回す | 未着手 |
| F-cf9 | **蒸気 c_p,v の出どころ** (2026-09-14 決着・実装済): Kirchhoff の傾き L'=c_p,v−c_l の c_p,v は凝縮する蒸気の値であって混合気の値ではない。pure-condensible CPG (気相=凝縮種) のときだけ `physProp.cp` を `CondPropOpts::gasCp` 経由で渡し、CPG carrier (空気, physProp.cp=1008.7 は混合気) と TP では N2 の 1038.8 を使う。既存 run は case/34 の N2 が `cp: 1038.8`、空気が carrier なので**どちらもビット不変**。残件: 二相 EOS 側は液相の顕熱に混合気の c_v+R_w (空気 1014.7) を当てており、EOS 上の実効 c_l が意図した 2000 でなく 1976 になる (1.2 %)。これは CPG carrier 形の定式化の近似で、直すなら液相の顕熱を凝縮種の c_p で持つ必要がある | 実効 c_l の 1.2 % は未着手 |
| F-cf7 | **CPG × `condGasSpecies` が未ガードで、ソース項と EOS が食い違う** (2026-09-14 発見): ソース kernel の carrier 判定は `condGasSpecies>=0 || condVaporMassFraction>0` なので p_v=ρ(Y_w−g)R_wT を使うが、CPG の二相 EOS (`dependentVariables_d.cu` の else 枝) は `condVaporMassFraction>0` しか見ないため pure 扱い (p=(1−g)ρR_mixT, e=e_v+gR_mixT−gL) になる。`solverConfig.cpp` の検証も `condGasSpecies>=0` に thermalMethod 2 を要求していない。case/16 の旧 fig3 run (run_0049–0052, 幾何が不一致で使用停止) だけが該当し、生産 run は全て TP。Y_H2O=0.011 で p が最大 0.6 %、e が 0.24 % (T にして 0.5 K) ずれる。対処案: `condGasSpecies >= 0` なら `thermalMethod == 2` を要求する検証を足す (1 行)。または CPG 枝でも `condGasSpecies` を見て carrier 形にする | 未着手 |

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

## 10. 未確定事項

- なし (課題ごとに着手時に書く)。
