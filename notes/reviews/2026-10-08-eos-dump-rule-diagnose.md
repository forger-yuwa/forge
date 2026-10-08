# codex 諮問 (diagnose): eos-dump-rule

- **brief**: [`notes/reviews/briefs/2026-10-08-eos-dump-rule.md`](../../notes/reviews/briefs/2026-10-08-eos-dump-rule.md)
- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-08
- **commit**: `24a843f4` (feature/sern-design)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **有限性規則を EOS の読み書きに合わせて改訂し、(b) の診断専用再生で「本当に床下の凍結入力」に対する旧版／新版 EOS 1 回の直接比較を行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 有限性を「読む入力の pre」「書き込み先の post」に限定する | **採用。ただし EOS 単体の非干渉判定に限定する。** TP の `gamma`・`cp` は既存値を読まず上書きする（[dependentVariables_d.cu:251](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:251)）。現比較器は引数配列の pre/post 全部に有限性を要求している（[compare_eos_dump.py:162](/home/sano/work/forge-sern-design/solver_density_cuda/tools/compare_eos_dump.py:162)）。読み書き集合をコードから固定し、**壁・ghost を含む全書き込み先の比較、未変更入力の不変検査、差分 0 バイトは維持**する。旧 `INVALID` は保存し、改訂規則の再判定を別記する。 |
| Major | ghost の NaN を無害として別件化する | **別件化は採用、無害の認定は却下。** EOS が読まないことと、ソルバ全体が読まないことは別。初期化では EOS 直後に `gasProperties` が呼ばれ（[main.cpp:1851](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1851)）、別の設定経路には ghost を含め `cp_array` を読む処理がある（[gasProperties_d.cu:84](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/gasProperties_d.cu:84)、同:94・104）。今回の `viscMethod: 2` でその経路による汚染を確認したわけではない。**EOS 比較規則の変更は認めるが、初期化全体の健全性は別途参照順序を監査する。** |
| Major | B 枝の受入れ方法 | **(b) を採用。(c) は却下、(a) は今回の受入れ試験には採らない。** 初期化 EOS が床を適用してから基準状態へ保存する（[main.cpp:1846](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1846)、同:1876）。初期化位置へのダンプ移動は、境界ピン後という登録条件を変え、ghost の組成未設定も持ち込む。**診断専用の再生経路で境界ピン後の凍結状態を読み込み、EOS 1 回で終了する。** |
| Major | 同じ規則で AWS の g3 へ進む | **改訂規則での EOS ダンプ比較は採用。g3/g4 の時間積分による検証への移行は要再検証。** 登録入力には g3 継続場が明記されている（[plan:84](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:84)）。直接の床下比較と g3 の比較を通してから移行する。TP 時間発展の回帰は「判定不能」のまま残す。 |

結論: **有限性規則を EOS の読み書きに合わせて改訂し、(b) の診断専用再生で「本当に床下の凍結入力」に対する旧版／新版 EOS 1 回の直接比較を行う。**

第 1 仮説: **現在の TP の `INVALID` は数値差ではなく、書き込み専用配列の初期内容にまで有限性を要求した比較規則による。** 確度: **高**  
  根拠: 保存されている `tp_abB_old_r1`／`tp_abB_new_r1` の HDF5 を独立に再比較した。EOS 引数の post は **0 / 5,388,516 バイト差**、その他を含む全配列も pre/post とも差分 0。非有限は pre の ghost `gamma`・`cp` 各 **2082 点**だけで、実節点と post 全配列は有限だった。保存比較記録も同じ結果を示す（[cmp_tp_abB_old_r1_vs_new_r1.txt:22](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-eos-dump/cmp_tp_abB_old_r1_vs_new_r1.txt:22)）。  
  反証条件: 対象経路の EOS が既存の `gamma`・`cp` を上書き前に読むこと、または読む入力／書き込み後の出力に非有限や再現する旧新差が見つかること。

第 2・第 3 仮説: **ghost の組成未設定で初期化 EOS が非有限を作り、境界条件が主要状態量だけを埋め直している。** 確度: **中・発生直後は未確認**。入力保存量は実節点までしか読み込まれず（[variables.cpp:801](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:801)）、境界条件は初期化 EOS より後。実測の ghost `Rmix` は pre で全点 0。ただし、発生箇所の確定には初期化直後の記録が必要。

判別 A/B: **変えるのは EOS の実装だけ：旧版／新版カウンタ有効。各腕 2 回、時間更新 0 回、EOS 1 回。**

- 境界ピン後の共通状態から、指定内部節点の `roe` だけを事前登録した床下値にした入力を固定する。両腕へ同一の配列・入力 `T`・組成・物性 DB・設定を読み込む。
- 読み込み後、EOS 前に再ピンや別の EOS を挟まない。実際の pre で床下条件と登録済みの余裕を検証する。
- 合格条件は、読む入力の有限性・同一性、全書き込み先の post の有限性・差分 **0 バイト**、未変更入力の不変性、腕内再現性。新版の対象 EOS 行で指定節点の温度床事象 **1 件**を確認する。

→ **結果 A：全条件を満たす**なら、試験した床下状態での数値干渉仮説を退ける。  
→ **結果 B：同一入力から再現する出力差がある**なら、「今回の問題は比較規則だけ」という説明を退け、最初の差を調べる。入力不一致・非有限出力・床下条件不成立は試験不成立とする。

やらない方がよいこと: **B 枝という入力名だけで床通過を検証済みにすること。** 対象節点 18259 は実測で pre の `T = 50 K`、post の `T = 50.000011444 K`。記録も初期化の温度床 **1 件・補正量 449.740526**に対し、step 1 の温度床は **0 件**である（[SUMMARY.json:886](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-eos-dump/SUMMARY.json:886)）。ghost を比較対象から外したり、NaN をゼロ置換して合格させたりもしない。

呼び出し側の前提への異議: **「全 7 入力が改訂規則でビット一致」と「登録した床上下の試験を完了」は同義ではない。** 現結果が示すのは、初期化を通過した状態に対する一致であり、床下入力に対する EOS の全出力一致ではない。

不足情報: g3 の比較結果と、床下入力を直接与えた比較結果。保存比較記録 42 組は確認したが、元ダンプが残る旧新 4 組だけを独立に再計算できた。**plan 未反映**。呼び出し側で `plans/active/tooling-sern-te-wake-grid.md` §5.1 #2・§6.1 に規則改訂、直接比較、適用範囲を記録する。
