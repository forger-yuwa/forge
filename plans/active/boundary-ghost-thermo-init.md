# ghost の熱力学状態・物性の初期化と参照順序

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md)
- **related_plans**:
  - 起点: [`tooling-sern-te-wake-grid.md`](../accepted/tooling-sern-te-wake-grid.md) §5.1 #2・#8 (EOS の凍結入力の比較で発見、2026-10-08〜10)
- **created**: `2026-10-11`
- **owner**: `CFD Dev (Claude セッション: SERN 3D)`

## 1. 目的

step 1 の EOS 直前に ghost の状態の一部が非有限になっている。(a) TP の全 ghost で gamma・cp が NaN (EOS は書くだけで読まない)、(b) farfield (BCONDS/10、`side_far`) の ghost 51,143 点 (g3) で roe・T・P・Ht・sonic が NaN (EOS は roe・T を読む)。
実節点の結果は 20000 step を通じて有限だが、**無害とは認定していない** (EOS が読まないことはソルバ全体が読まないことではない)。発生箇所を確定し、ghost を読む処理の参照順序を監査する。

## 2. スコープ

- **やる**: 初期化直後の記録で発生箇所の確定、ghost の gamma・cp・cp_array・roe・T を読む処理 (gasProperties ほか) の洗い出しと参照順序の監査、必要なら初期化の修正 (数値の振る舞いを変える修正は上位に諮る)。
- **やらない**: farfield の流束の方式の変更。

## 3. 関連 docs と前提

推定 (未検証、codex diagnose 2026-10-10): 初期化の EOS は境界条件より先に走り (`main.cpp:1846/1857`)、入力の保存量は実節点まで (`variables.cpp:801`)、farfield の状態の境界処理は何もしない (`boundaryCond.cpp:498`)、組成が全成分 0 だと温度反転の研磨で cp = R = 0 の分母保護が 0 になる (`thermo_d.cuh:742`)。roY は後で Neumann コピー、主要状態と RANS の ghost は充填されない。
記録: `notes/investigations/2026-10-08-eos-dump/` (g3 の比較・ghost の index 対応)。初期化の `gasProperties` が ghost を含めて `cp_array` を読む経路 (`gasProperties_d.cu:84/94/104`)。

## 4. 設計方針

未定 (§5.1 #1 の測定の後に書き、上位に諮る)。

## 5. 実装ステップ

1. 初期化直後・境界条件の後・step 1 の EOS 直前の ghost の状態を記録して発生箇所を確定。
2. ghost の各配列を読む処理の一覧 (カーネル名・読む時点) を作る。
3. 修正の要否の判断 (上位)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 発生箇所の確定 | §5 の 1。小さい 3D (farfield あり) と 2D TP で | O |
| 2 | 参照順序の監査 | §5 の 2 | O |
| 3 | 修正の要否 | 上位に諮る | F |

## 6. 検証

§4 を決めてから書く。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

TP・farfield を使う全ケース (ghost の状態)。

## 8. 完了条件

発生箇所と参照順序が確定し、修正するか無害の根拠を示すかが決まり、codex レビューを経ること。

## 未確定事項

- 定常陰解法で EOS の床が保存量に commit されない (作業状態と保存量の乖離) ことの是非 — `tooling-sern-te-wake-grid.md` から追跡先として写した (2026-10-11)。同じ初期化・EOS の順序の論点として扱う。
- 実節点への影響の有無 (現状 20000 step で実節点は有限だが、fmaxf 等で NaN が隠れる経路の有無は未確認)。

## 9. 変更ログ

- `2026-10-11` — 起票。`tooling-sern-te-wake-grid.md` §5.1 #8 から切り出し。
