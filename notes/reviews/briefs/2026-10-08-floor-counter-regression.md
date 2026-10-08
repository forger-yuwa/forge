# 諮問: 床事象のカウンタ (出力専用) の回帰で TP の時間発展が判定不能 — 受入れの物差し (2026-10-08)

関連 plan: `plans/active/tooling-sern-te-wake-grid.md` §5.1 #2 (回帰の事前登録: 決定的な比較はビット一致、組立残差・時間発展は `procedures/verification/README.md` の「振る舞いを変えないはずの変更の回帰」)。エスカレーション条件 3 (事前に書いた比較が判定不能)。実装 commit 0ffc1a21。

## 観測事実 (実装担当の測定、スクリプト scratch `regress_floor_events.py`、記録 `regress_*/REGRESSION*.txt`・`frozen_width.json(.sha256)`)
- 変更: EOS カーネル (`dependentVariables_d.cu`) に出力専用の計数 (`fe.acc == nullptr` なら計数しない分岐)、終了時の監査 (コピー上)、記録の書き出し。数値は変えない設計。
- **決定的な比較: 3 構成すべてビット一致** (TP 初期場から・TP 発達場 [HEAD で 2000 step の場を restart_field で継続]・CPG c52cht)。対象: res_0 の全データセット (初期化の EOS の後)、res_1 で HEAD 5 反復の全要素が一致したデータセット (step 1 冒頭の EOS の出力。TP では T・Ux・Uy・Uz・h0・k・omega・roUz・sonic・vis_lam・vis_turb・wall_dist) と壁 689 節点を除いた P (等温壁の P は commit 後の ρ でピンし直されるので HEAD 同士でも 1〜50 節点割れる)。カウンタ無効と有効 2 反復のすべてで一致。
- 小型 A/B (計測の閉じ方、事前登録): TP 59 件 OK (親が再実行して VERDICT PASS)、CPG 50 件 OK / 1 SKIP (CPG の pMin 1 Pa では圧力床だけを起こせない)。監査の前後で 223 本の cell 配列がビット不変。
- **時間発展** (HEAD 3 反復で幅 = 2 × HEAD 間の最大差を凍結 → 有効 3 反復で評価): CPG 200 step PASS (HEAD の追加 3 反復も幅内 = 判別力あり)。TP 初期場から 200 step: 変更後が Ux で幅超過 (36.5 > 34.7)、HEAD の追加反復も T・Ux・roUy で超過 → 判定不能。TP 発達場から 200 step: 変更後が omega で超過、HEAD の追加 r6 も P・ro・roe で超過 → 判定不能。TP 発達場から 20 step: 変更後は全量が幅内、HEAD の追加 r6 が rms_roUx で超過 → 厳密には判定不能。幅・本数は増やしていない。
- 既知: 「新旧差 ≤ 旧 3 回の範囲」型の規則は同分布でも高い確率で FAIL する (メモリ repeat-range-as-limit-rule-trap)。forge は atomicAdd の加算順で同一入力でも再実行で揺れる。
- 用途: この後の g3/g4 の判定区間 (§5.1 #4) は**両腕とも同じカウンタ入りのバイナリ**で回す。旧 run (run_1078/1079 = 変更前のビルド) との数値の直接比較は予定していない (run_1079 の最終場を初期場に使うだけ)。

## 問い
1. 出力専用の変更として、決定的な比較のビット一致 (3 構成) + CPG の時間発展 PASS + 小型 A/B で受け入れてよいか。TP の時間発展が判定不能のまま進めてよい条件。
2. 足りないなら、TP で判別力のある物差し (頻度比較・多反復の分布比較・より決定的な比較の拡張 [例: 組立残差の atomicAdd を避けた比較]) を何で事前登録するか。
3. カウンタ無効時のカーネルの同一性を機械語で示す (SASS 比較) ことは有効か (分岐が増えるので一致しない可能性が高い)。
