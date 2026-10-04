# 諮問: farfield 境界と SERN R8 の到達点・生産切替の妥当性・残作業の優先順位 (2026-10-01)

関連 plan:
- `plans/active/boundary-node-farfield-characteristic.md` (§4 設計、§5.1 #1–#4d、§6 V0–V3、§6.1 レビュー記録。status in_progress)
- `plans/active/tooling-nozzle-sern-chain.md` の R8 行 (種 DB 追従、完了扱い)
- `plans/active/tooling-nozzle-sern-3d.md` の R4d 行 (側方領域)
コード: `feature/sern-design` 2a17b4bc (species-transport を 237357fd で直接マージ済み、main には未統合)。
run: AWS `case/58.farfield_verification` (V0–V2)、`case/46.sern_design` (V3・R8、一部は `~/forge-r8` 側)。一覧は各 case README。
過去の諮問: `notes/reviews/2026-09-29-farfield-*-diagnose.md` (3 件)、`2026-09-29-farfield-v3-result-diagnose.md`、`2026-09-30-sern-r8-stage0-thermo-diagnose.md`。すべて全件採用。

## 1. これまでの結果 (観測事実)

### 1.1 farfield 境界そのもの (V0–V2)
- 合格: V0 (既存 BC のビット不変)、V0u/V0k (単体、GPU 一致 2.7e-8)、V0u(v) 起動拒否、V1 (自由流保持 ≤ 8.7e-6)、V2a (音響反射 0.10 %、SLAU2 0.12 %、dt/nSub の時間精度 OK)、V2b (保存収支 恒等式 ≤ 1e-9・収支 ≤ 6e-9)、V2e (陽解法・SST dual-time)、V2f (局所逆流で ṁ<0 面が外気組成)、V2d-1 (接触波 3 物性)。
- 未達・判定不能 (FAIL は保持):
  - V2c 斜め衝撃波: 旧形状 (ランプ→水平) は B vs C 0.087 Δp (圧縮角、A=B でありソルバの角の非一意、作用素は同一、絶対残差 A/B で凸角近傍に非零残差が残る)。凸角を無くした形状 v2 は slip の対照 A が発散、B・C は 2.5 桁プラトーで事前の必要条件未達 → 判定不能。同一初期場での境界の効果は壁 0.008 Δp (出口端 1 節点 0.021)。
  - V2d-2 TP 音響: 元配置は左端流入の接触面 (TP 保存形混合) で擾乱 2.5 Pa → FAIL 保持。右端隔離配置 (左端高温) で対照・反射 (0.06–0.9 %) は合格。**時間精度は dt 8.8e-7→4.4e-7 で差 1.13 %、4.4e-7→2.2e-7 で 2.18 % と dt 細分で悪化** (float32 の時間微分の丸め床の疑い、振幅 2.2 Pa / 2851 Pa)。
  - 独立参照 1D (`ref1d_euler_tp.py`): Δx 5 mm・1.25 mm とも参照自身が §6 の収束条件未達 (手元 CPU 制約で中断、記録用で合否外)。
- 途中で見つけて直した: GPU 上の構造体代入で組成・k/ω が 0 になる欠陥 (V0k 追加)、リミッタ基準値の自動決定が領域長に依存して比較を汚すこと (比較では固定)、帳簿の段の前後差の読み方、プローブ出力桁。

### 1.2 SERN での結果 (V3、g3、M6、生産設定)
- 側方 2.50 H で slip → farfield: C_L 0.06259 → 0.06182、C_M −1.6442 → −1.6212 (D > ε、記録)。
- farfield の幅系列 2.50/3.42/4.35 H: 2.50 と両広幅で全 4 量 D ≤ ε (C_L 8e-5、C_M 2.2e-3)、広幅同士 3.4e-5・7.4e-4 → 事前規則で必要幅 2.50 H (試験系列内、限定付き)。C_M の 2.50 H の平均差は広幅と 1.6e-3 (ε の 33 %)。
- 初期場履歴 A/B: 全量 D ≤ 0.2ε → 履歴依存を棄却。
- 格子差 G (g3→g4、同じ新バイナリ・固定基準値): farfield で C_L −0.00052・C_M +0.0159、G + D は §8 総許容内 (slip では C_L 0.0021・C_M 0.062 で超過)。slip の G は新旧で同じ (+0.00139/−0.0407 vs +0.00133/−0.0389) → G の改善は farfield 化による。slip→farfield の C_M の動きは g3 +0.023、g4 +0.080 (機序未同定)。

### 1.3 R8 (種 DB 追従)
- species-transport を直接マージ (衝突 5 件解消)、外気 lump AIR→AMB、lump 記法、ランナーの `gas.transport` 対応 (作動点ごとに実種へ絞る)。
- 段 (i) (改名 + lump、輸送は旧): 2D m6_on 3 反復ずつで 4 量とも差がノイズ床内 → PASS。照合 (0) は EXH 係数最大 7 ulp で事前基準 1e-15 超え → 演算経路 A/B で旧・新とも保存値をビット再現 → 丸めとして記録 (FAIL は保持)。
- 段 (ii) (種ごとの輸送、記録のみ): 2D で C_T −0.003 %、C_T_with_shear +0.019 %、摩擦 +0.64 %、C_L −0.029 %、C_M +0.023 %。
- 依頼文と実装の食い違い: EXH 構成種 H2/OH/NO/H/O/CO はソルバ内蔵でなく (`legacy_builtin: design` のみ)、生係数を `species_db_external.yaml` で渡している。

### 1.4 生産への切替 (ユーザ決定 2026-10-01)
- SERN 生産 YAML 8 件に `gas.transport` (種ごとの輸送)、3D 生産 YAML 7 件に `side_far_kind: farfield`。
- 組み合わせの初回確認 (`~/forge-r8/case/46.sern_design/run_1011`→延長 `run_1012`、run_0997 の最終場から `--force-species` で移行): GATES PASS・窓条件 OK・置換 0。run_0997 (旧輸送) からの変化: C_T +0.004 %、C_T_with_shear −0.003 %、**C_L +0.23 % (D 1.7e-4)、C_M −0.22 % (D 4.1e-3 = ε の 82 %)**。2D の約 0.02 % より一桁大きい (理由未同定)。

## 2. 残っていること
1. farfield plan の未達: V2c 判定不能、V2d-2 時間精度、独立参照の収束。plan は in_progress。result 段の codex レビュー未実施。
2. 3D で輸送切替による C_L・C_M の変化が 2D より一桁大きい理由。
3. slip→farfield の C_M 変化が格子で違う (g3 +0.023、g4 +0.080) 理由。
4. species-transport・sern-design の main への統合 (両ブランチとも main と大きく乖離)。
5. 依頼元セッションへの「内蔵種」食い違いの連絡。
6. 生産 run の旧場 (AIR 名・属性なし) は `--force-species` で 1 回移行が要る。移行後の出力も未検証属性 (species_input_unverified=1) のまま。

## 3. 問い
1. 1 章の解釈に誤り・言い過ぎはあるか (特に V3 の必要幅の限定、G の改善の帰属、R8 段 (i) の合格判断)。
2. 生産への 2 つの切替 (側方 farfield、種ごとの輸送) は、この証拠で妥当か。足りない確認があれば何か (例: C_M が ε の 82 % 動いたことの扱い、3D での輸送切替の格子依存)。
3. 2 章の残作業の優先順位と、farfield plan を accepted にするための最小の道筋 (V2c・V2d-2 を別 plan に切り出すべきか等)。
4. 6 (移行後の出力が未検証属性のまま) は問題か。生産の継続 run で毎回許可が要る運用になっていないか。
