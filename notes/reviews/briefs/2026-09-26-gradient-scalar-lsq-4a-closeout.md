# 諮問ブリーフ: スカラー LSQ 統一 Phase 1 — #4a hook の commit 可否・ダンプ非干渉の判定規則・S0/S1 を閉じて S2 へ進むか (2026-09-26)

plan: `plans/active/gradient-scalar-lsq-unification.md` (§5.1 #4a–#4c・#5、§6 S0/S1/S1 (2) 訂正/S2)。
引き継ぎ: `notes/sessions/2026-09-26-gradient-scalar-lsq-handoff.md`。
ワークツリー `/home/sano/work/forge-sern-design`、ブランチ `feature/sern-design`、HEAD `eb8fb6a4` (実装は `7ebf6d7d`)。

決めてほしいこと (3 点、それぞれ結論と理由を 1 段落で):

1. **#4a 出力 hook の diff は commit してよいか** (欠陥があれば `ファイル:行` で)。
2. **「ダンプ有効で res 不変」の判定規則** — 事前に書かれていなかった。親の案 (下の仮説 H2) で良いか、別の規則か。事後規則であることは plan に明記する。
3. **§6 S1 (2) の訂正 (#4b = A を条件に事前に書いたもの) を適用して S0/S1 を PASS で閉じ、S2 に進んでよいか**。S2 着手前に足すべき確認があれば挙げる。

## 読んでよいファイル

- `plans/active/gradient-scalar-lsq-unification.md` (全文)
- hook の diff: `notes/sessions/gradient-scalar-lsq-unification-4a-hook.patch` (作業ツリーの未 commit 差分と同一。`git diff solver_density_cuda` でも同じ)
- 判定結果: `case/09.Taylor-Green/_g0_lsq_seam/{S1_hook4a.txt,S1.txt,DUMPCHECK_case48.txt,PREGATHER_4b.txt,S0e_case39.txt}`、`S0*_*.txt` は VERDICT 行だけ
- 判定スクリプト: `case/09.Taylor-Green/_g0_lsq_seam/{s1_nointerference.py,pregather_check.py,r3_compare.py}`
- 周辺コード (該当関数だけ): `solver_density_cuda/cuda_forge/periodicNode_d.cu` (`periodicGradientGather_d_wrapper`・`periodicGather1ToRoot_d`)、`ransTransport_d.cu` の `ransGradient_d_wrapper`、`speciesTransport_d.cu` の `speciesGradient_d_wrapper`・`passiveGradient_d_wrapper`、`main.cpp:1290-1310,1460-1505`、`output/output.cpp:30-80`
- 巨大ファイル (`main.cpp` 全体、`calcGradient_d.cu` 全体、`*.h5`) は読まない

## 観測事実

### hook (#4a) の内容 (diff 要約)

- `preGatherDump(tag, nCells, names, grads)` (`periodicNode_d.cu` 末尾): env `FORGE_DUMP_PREGATHER` が無ければ即 return (同期も D2H も無し)。有れば tag ごとに**最初の 1 回だけ** (static `std::set`) 勾配を `cudaMemcpy` (同期) で読み、raw float `[nVar][nCells][3]` と `.names` を書く。
- 呼び出し位置: `main.cpp` の 2 箇所の `periodicGradientGather_d_wrapper` 直前 (tag `init`・`loop1`、NS 6 量 × 3 成分と、gg 経路で同 gather に登録される dY・dξ/モーメント)。`ransGradient` の lsq 分岐 (gather 直前) と gg 分岐 (div_vol 後・gather 直前) に tag `rans.loop1`。`speciesGradient`・`passiveGradient` の lsq 分岐 (wrapper 内合併の直前) に `species_lsq.loop1`・`passive_lsq.loop1`。rans/species/passive の wrapper は `assembleResidual` からしか呼ばれない (`main.cpp:1473,1474,1501`) ので tag の `loop1` は正しい。
- `preGatherDumpMain` の登録条件 (`speciesFaceReconstruction >= 1 && !lsq`、受動種は `passiveScalarScheme == 1`) は `periodicGradientGather_d_wrapper` の登録条件と同じ。divU はダンプしない (plan #4a の「NS 18 本」どおり)。
- 出力: `variables.hpp` に `extraOnly_cellValNames = {"wall_y_eff"}`、`registerSpecies` が `dY{s}d{x,y,z}` を足す。`output.cpp` は `output.extraFields` に書かれたときだけ末尾に足す (level 2 の既定出力にも入らないので既存 run の `res_*.h5` の中身は変わらない)。
- env 無しでも毎 step 残るホスト側コスト: rans/species/passive の呼び出し位置で `std::vector<std::string>` と initializer list を組み (`var.c_d[...]` の map 参照を含む) から `getenv` で return。GPU 側の追加処理は無し。

### 測定 (AWS g5、共有インスタンス。バイナリ: 旧 `712141ec…` = 7ebf6d7d 以前、新 `a398afdc…` = 1a334275 + hook patch)

- **S1 (1) gg 旧 vs 新、hook 入りバイナリで再実行** (`S1_hook4a.txt`、4 ケース × 旧 3・新 3 本): case48・case16・axi・tgv すべて PASS (訂正規則: 両側 3 本のノイズ対)。
- **S1 (2) NS 勾配・リミタ gg vs lsq**: case48・case16・axi は 48 配列すべて 15 対でビット一致。**tgv は dUxdz (res_1) が FAIL**: gg 同士 0 不一致、lsq 同士 n 4 / 7.45e-9、gg–lsq n 32 / 2.98e-8。不一致は 3 member 以上 (4 member) の group のみ。同じ規則の片側条項「旧同士がビット一致の配列は旧新もビット一致」による FAIL。
- **#4b (事前に書いた A/B、`PREGATHER_4b.txt`)**: tgv gg 3・lsq 3 本で gather 前の NS 18 配列 (init・loop1) は **6 本すべてビット同一**。gather 後の不一致節点 (res_0 dUx 124、res_1 dUx 232、全て 4 member) は **全点が member 部分和の float32 順列和の集合に含まれる**。→ **A**。
- **ダンプ非干渉 (ダンプ無し gg 3 本 vs ダンプ有り gg 3 本、S1 の規則を流用)**:
  - tgv (`PREGATHER_4b.txt` 末尾): 168 配列中 141 は 15 対ビット一致。**FAIL 2 配列**: res_0 dUxdz (無し同士 0、**有り同士 32**、無し–有り 32 / 2.98e-8)、res_1 dUxdy (無し同士 0、有り同士 4、無し–有り 4 / 5.96e-8)。その他 25 配列 PASS。
  - case48 (`DUMPCHECK_case48.txt`): 152 配列中 135 ビット一致。**FAIL 1 配列**: res_1 roe (無し同士 0、**有り同士 3**、無し–有り 2 / 0.0156)。同じ roe は S1 の旧同士でも n 2 / 0.0156 (`S1_hook4a.txt` case48)。case48 は周期無しでも k/ω 勾配・残差が run 間で揺れる (GG の面ループ atomicAdd)。
  - いずれの FAIL も「無し同士がビット一致 → 無し–有りもビット一致」の片側条項で、有り同士のノイズは無し–有りの不一致数以上。
- **S0** (`S0*_*.txt`): tgv 系 5 変種・box_slip・channel・jitter32・axi・axi_m1 の a/b/c(/d) すべて PASS。Y 5 種 (#4c) は channel・box_slip・tgv・jitter32 で S0-a Y・S0-c dY/dξ・ΣdY すべて PASS。S0-e は channel (k・ω・wall_y_eff) と case/39 起点で group 内差すべて 0 (記録)。

## 期待値と出典

- hook は「既定 off・出力専用・数値不変」(plan §5.1 #4a、`diagnostician` 判断)。
- S1 (1) の規則と「旧同士ビット一致 → 旧新ビット一致」の片側条項は §6 S1 (測る前に 2026-09-26 訂正)。
- §6 S1 (2) の訂正 (#4b = A を条件、測定後だが決定的な段は厳しくする方向で事前に書いた): NS 勾配・リミタは **gather 前の局所配列でビット一致 (必須)**、gather 後は 2 member 以下ビット一致・3 member 以上は順列和の集合内。「旧同士ビット一致なら旧新もビット一致」は **gather を含む配列には適用しない**。
- 「ダンプ有効で res 不変」の合否規則は**事前に無かった**。

## 再現条件

- ハーネス `case/09.Taylor-Green/_g0_lsq_seam/` (`s1_nointerference.py`、`pregather_check.py`)、AWS `~/sglsq/s1h/`・`~/sglsq/pg/`。block `FORGE_CUDA_BLOCKSIZE=128 FORGE_CUDA_BLOCKSIZE_SMALL=128`、1 step (res_0・res_1)。

## 仮説 (親の案)

- **H1 (commit 可)**: env 無しでは GPU の命令列も同期も変わらない (早期 return) ので数値不変。S1 (1) 再実行 4 ケース PASS がその裏付け。ホストの毎 step の小コストは S3 の測定に入るが無視できる見込み。
- **H2 (判定規則)**: ダンプは `cudaMemcpy` の同期で kernel の実行タイミングを変えるので、atomicAdd を通る配列の順序ノイズが変わるのは想定内。規則は 2 段: (a) **決定的な段**: gather 前の局所配列 (NS 18 配列) がダンプ有り 3 本の間、およびダンプ無しの lsq/gg とビット同一 (#4b (i) で既に 6 本同一を示した)。必要なら `FORGE_DUMP_MASSFLUX` の step 1 面流束もダンプ有無でビット一致を見る。(b) **場**: 両側 3 本のノイズ対 (無し同士 3 + 有り同士 3) による S1 の訂正規則。ただし片側条項は atomicAdd を通る配列に適用しない (§6 S1 (2) 訂正と同じ理由)。この規則で tgv・case48 の 3 FAIL はいずれも PASS (無し–有りの不一致数・最大差 ≤ 有り同士の値 × 2)。**事後規則と明記**する。
  - 懸念: 片側条項を外すと、決定的な配列 (atomicAdd を通らない) に本当の差が出ても両側ノイズで隠れうる。対策案: 片側条項は「無し同士・有り同士の**両方**がビット一致の配列」にだけ適用する (これなら tgv・case48 の 3 件は有り同士が不一致なので対象外)。
- **H3 (S0/S1 を閉じる)**: #4b = A なので訂正を適用すれば tgv S1 (2) は PASS (gather 前 6 本ビット同一、gather 後は全点が順列和の集合内)。S0 は全変種 PASS、#4c PASS。→ S0/S1 を PASS で閉じ、S2 (§6 表 4 ケース + case/16 双子 2 + FCT smoke) に進む。S2 の前提作業 #2e (起点の sha256 照合) と #2f (収束規則スクリプト) は S2 の最初にやる。

## 潰した候補の証拠

- 「lsq 経路が NS の入力を壊している (#4b の B)」: gather 前 NS 18 配列が gg 3・lsq 3 本でビット同一 (`PREGATHER_4b.txt` (i))。
- 「gather 以外の非決定源 (B′)」: gather 後の不一致全 356 節点×成分が順列和の集合内、全 group でも集合外 0 (`PREGATHER_4b.txt` (ii))。
