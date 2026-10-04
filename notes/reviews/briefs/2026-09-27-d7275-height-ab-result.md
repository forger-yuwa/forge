# 諮問ブリーフ: D-7275 試験 26 の上端高さ A/B (T4-0b-H) の解釈と G15 収束ゲートの次の一手 (AGENTS 条件 7・4)

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` §4.12、§6 G15、§5.1 #61 の末尾
- 事前登録と結果: `case/60.flatplate_d7275_m7/acceptance.json` の **T4-0b-H** (`outcome` に数値を全部記載)、同 T4-0b-WALL
- 台帳: `case/60.flatplate_d7275_m7/README.md` (run_0012〜0014 の行)
- 前回の諮問: `notes/reviews/2026-09-27-d7275-top-corner-shock-diagnose.md` (この A/B の設計はこの提案どおり)
- 集計スクリプト: `case/60.flatplate_d7275_m7/tools/height_ab.py`、初期化 `tools/stack_init.py`、メッシュ生成 `case/56.gap_tp1187/gen_mesh.py --top-layers`
- 判定ツール: `solver_density_cuda/tools/check_convergence.py` (`--from-floor`・`--segment` の仕様)

## 1. 観測事実 (run は AWS。数値は acceptance.json の T4-0b-H outcome と同じ)

- A = `case/60.flatplate_d7275_m7/run_0013_t26_h05` (H 0.5、run_0012 最終場から restart_field、12 量ビット一致)、
  B = `run_0014_t26_h079` (A の全節点を座標完全一致で含み上に 8 層、H 0.794506 m。共通節点はビット一致コピー、追加節点は入口自由流)。各 20k、同一バイナリ (全域 FP64)・同一設定。
- 窓間の中央値変化 (12–16k → 16–20k): A ≤ 1.7 %、B ≤ 2.1 % → 事前登録どおり 20k で判定。
- 旧上端域 (x 2.6–2.8、y 0.45–0.50) の √Σres² B/A: ro 6.6e-4、roUx 7.6e-4、roUy 9.2e-4、roe 5.3e-3、roK 8.8e-4、roOmega 7.0e-4。
- 点プローブ (10 step 毎、末尾 4k) の P 振幅/平均: 旧上端直下の 3 点で A 3.8e-5〜1.5e-4、B は 0 (出力が 6 桁なので分解能 ~1e-6 以下)。
- B の残りの流れ残差は局在しない。自由流域 y 0.03–0.8 m に x 全域で薄く分布 (res_ro の二乗和の 43 % が y 0.5–0.7、23 % が 0.7–0.75、31 % が 0.03–0.5)。最大は (2.156, 0.7498) の 1.3e-7 (A の最大 3.9e-5 の 1/300)。
- 全域 √Σres² B/A (最終場): ro 0.032、roUx 0.042、roUy 0.022、roe 0.10、roK 0.20、roOmega 0.87。履歴 16–20k の √Σres² B/A: ro 0.025、roUx 0.032、roUy 0.017、roe 0.080、roK 0.17、roOmega 0.995。
- 比較域 (1.07–2.55 m) の流れ残差は B が A の 1.4–2.0 倍 (絶対値は A の旧上端の 1/100 以下)。
- **rms_roOmega は両腕 0.6 前後で横ばい** (B/A 0.995)。res_roOmega の二乗和の 99.99999 % が壁の**第一内部節点列 (y = 3 µm)** にあり、平板全長の個別の節点 (x 0.72、0.90、1.26、2.75、2.78 m など、上位 5 点で 40 %) に散在。壁上 (y = 0) は 0 %。
- 参考: Cary の PASS run (`case/59.flatplate_cary_m6/run_0005_tw02_re027_fine/CONVERGENCE_VERDICT.txt`) も rms_roOmega 最終値 3.0e-1 だが、段階起動の開始値 2e8 から 8.7 dec 低下なので PASS。
- B 自身の check_convergence: `NOT CONVERGED (stalled/plateau)`。この run の開始値から流れ・化学種 1.7–2.3 dec、roK 5.1 dec、roOmega 1.8 dec。
- 比較量: 比較域の壁 q_w の B/A−1 は最大 1.5e-6。R_A 平均 1.313 は同じ。St 10 点は両腕 ALL STEADY。

## 2. 期待値と出典

- 事前登録の読み (T4-0b-H decision): 旧上端域の各未達流れ・化学種残差 (√Σres²) が A の 1/10 以下で、局所振幅も減衰し、別の場所へ同水準の残差が移っていない → 上端近接説を支持。1/10 は診断上の支持条件に限り、G15 の代替合格条件にしない。
- G15: B 自身の同一設定区間で check_convergence PASS・St の check_quasisteady・壁解像。未収束の run_0012 を合格済みの参照床に使わない。

## 3. 実施済みの操作

前回までの A/B 5 本 (壁距離・出口・緩和・slip 延長・壁延長) と本 A/B。記録からの逸脱: 点プローブは T・P・U のみ (ρ・k・ω の時系列なし)。

## 5. 仮説・呼び出し側の見立て (未確認)

- 機械的な読みは「支持」。
- G15 の障害: (i) B は収束に近い場から始めたので、この run の開始値からの低下桁数が構造的に小さい。`--from-floor` は収束済みの参照 run を要するが、該当するものがない。(ii) rms_roOmega の床 ~0.6 は壁第一内部列の個別節点に由来し、上端とは無関係。Cary でも同水準の床があり、Cary は大きな開始値から測っているので PASS した可能性がある (未確認)。
- 次の一手の候補: (a) B の領域 (壁を出口まで、H 0.7945 m) で**一様初期値から段階起動をやり直し** (run_0001 と同じ段階起動 + 本段)、判定区間を `--segment` で最後の同一設定区間にする。所要は run_0001 と同程度。(b) B を延長して流れの床が下がるか見る (いまは横ばい)。(c) ω の第一内部列の床を別途調べる。

## 6. 聞きたいこと

1. 機械的な読み「支持」を確定してよいか。比較域の流れ残差が 1.4–2.0 倍になったこと、roOmega が不変なことは読みにどう入れるべきか。
2. G15 の収束ゲートの次の一手を 1 つ ((a)/(b)/(c)/その他)。(a) なら事前に何を固定するか。
3. rms_roOmega の床 (壁第一内部列の個別節点) の扱い: Cary の PASS とどう整合させるか。判定区間の取り方の問題か、別の問題か。
