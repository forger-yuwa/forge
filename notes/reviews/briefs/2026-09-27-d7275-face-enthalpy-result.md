# 諮問ブリーフ: D-7275 試験 26 の SLAU 面エンタルピー精度 A/B (T4-0b-HF) の解釈と、G15 へ向けた次の一手 (AGENTS 条件 7・1・6)

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` §4.12、§6 G15、§5.1 #61 の末尾
- 事前登録と結果: `case/60.flatplate_d7275_m7/acceptance.json` の **T4-0b-HF** (`outcome` に数値を全部記載、`stop_rule` あり)
- 台帳: `case/60.flatplate_d7275_m7/README.md` (run_0019〜0022)
- パッチ: `case/60.flatplate_d7275_m7/tools/hface_double.patch`、集計 `tools/hf_ab.py`
- 前回の諮問: `notes/reviews/2026-09-27-d7275-thermofloat-result-diagnose.md`
- 面エンタルピー float 化の出所: `plans/accepted/performance-3d-node-sst-speedup.md` §4.2-2

## 1. 観測事実 (AWS。acceptance.json の T4-0b-HF outcome と同じ数値)

- A = 面エンタルピー float (`~/forge56-double`、4687c3c 全域 FP64)、B = 同じ面状態 Y・T を double 係数・double 評価 (`~/forge56-hface` = 同じソース + パッチ 2 行)。両腕 thermoFloat 0、run_0018 の最終場からビット一致、各 5k + 延長 5k。
- 自由流域 (y 0.03–0.79) の √Σres² B/A (領域別時系列、窓中央値): 5k で ro 1.16e-3、roUx 8.7e-4、roUy 1.25e-3、roe 3.4e-4、roY 1.16e-3、roOmega 4.9e-3、roK 0.12。10k で同じ水準 (ro 1.17e-3、roe 3.5e-4、roOmega 5.7e-3、roK 0.16)。
- 窓間変化 (10k): A ≤ 0.8 % (横ばい)、B −13.5〜−30 % (まだ下降)。**延長の再開直後に両腕の残差が跳ねた** (B の rms_ro 6.8e-12 → 最大 5.5e-9 → 末尾 6.1e-12、A 6.15e-9 → 8.2e-9)。restart_field は保存量ビット一致。
- 全域 RMS 履歴 B/A (延長の 4–5k): ro 1.3e-3、roUx 1.0e-3、roUy 1.2e-3、roe 4.1e-4、roY 1.3e-3、roK 0.23、roOmega 0.53。
- 第一内部列 (y = 3 µm) の B/A (10k): roUy 0.072、ro 0.18、roY 0.19–0.33、roK 0.24、roe 0.34、roUx 0.41、roOmega 0.58。B の ω は横ばい (~80)。
- check_convergence: A NOT CONVERGED (全列横ばい)、B NOT CONVERGED (流れ・化学種 falling、roOmega stalled/plateau)。いずれも延長区間だけで判定。
- 点プローブ (1 step 毎、出力 6 桁): 両腕とも振幅は分解能 ~1e-6 以下で差は読めない。
- 比較量: 比較域 q_w の B/A−1 最大 4.4e-7、R_A 平均 1.313。

## 2. 期待値と出典

- 事前登録 (T4-0b-HF decision): B の自由流域の流れ・化学種残差が A の 1/10 以下かつ局所振幅減 → 支持 / 0.8–1.2 倍 → 棄却 / 3–4k と 4–5k の領域別指標の中央値が 5 % 超動けば一度延長、それでも動けば判別保留。stop_rule: 棄却・保留なら手当てを連鎖させず『G15 未達の探索結果として他項目へ進む』をユーザに提案。

## 5. 呼び出し側の見立て (未確認)

- 自由流域の床の主因は SLAU 面エンタルピーの float 評価 (量子化) — B/A ~1e-3 が 5k・10k で一貫。字面の『判別保留』は、B が A からさらに離れる方向に動いていること (と再開の跳ね) による。
- 第一内部列の ω 床 (~80、壁第一内部列の個別節点) は別の原因 (面エンタルピーでも thermoFloat でも消えない)。trans ≈ dest ≈ 6e8 の相殺の残り。
- G15 へ: 面エンタルピー double の B バイナリで**一様初期値から段階起動 + 本段**を回せば流れ・化学種は収束判定に届く可能性があるが、ω は届かない見込み。面エンタルピーの精度を選べるようにするには cuda_forge の数値カーネルの変更 (既定 float、opt-in double など) と plan が要る。

## 6. 聞きたいこと

1. 読みを「支持 (自由流域の流れ・化学種の床の主因は SLAU 面エンタルピーの float 評価)」と確定してよいか。事前登録の『それでも動けば保留』はこの動き方 (B がさらに下がる・再開の跳ね) でも保留と読むべきか。
2. stop_rule との関係: 支持なら次の一手は何か。(a) 面エンタルピー精度を opt-in にする plan を起こして (条件 1・6)、G15 用に B 設定で一様初期値から回す、(b) ω の第一内部列の床を別途調べる、(c) G15 を未達のまま他項目へ進む提案をユーザに出す。
3. ω の床 (相殺の残り ~80、全域 RMS 0.13–0.24) を G15 の判定でどう扱うべきか。
