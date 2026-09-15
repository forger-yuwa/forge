# Forge Verification Cases

この文書は、コード変更時にどのケースで確認するかを整理したものです。

## 標準検証ケース

標準検証ケースは次の 3 つとする。

- `case/08.bump`
- `case/20.naca_ml`
- `case/13.nozzle_H`
- `case/48.flat_plate_cooled_m4` (等温壁 × 低 Re SST の冷却超音速平板: [48-flat-plate-cooled.md](48-flat-plate-cooled.md))

## 既定の検証先

ユーザーから明示がなければ、まず `case/20.naca_ml/001.test/run_slau` を既定の検証先として使う。

このケースを既定とする理由:

- `run_case.sh` が用意されており、メッシュ生成から `forge` 実行までを一通り流せる。
- 出力ステップ数が比較的短く、変更確認の起点にしやすい。
- まず 1 本の標準ケースで異常の有無を見てから、必要に応じて他ケースへ広げやすい。

## ケースの使い分け

- `20.naca_ml`: 既定の確認先。明示がなければ最初に使う。
- `08.bump`: line profile 比較を含む既存の回帰確認フローが必要なときに使う。低マッハ・高マッハを陽解法/陰解法で流し、y≈0.3 の数値解の一致と収束度合いの非劣化を自動判定する検証一式 (`case/08.bump/verify/run_verification.sh`) を備える。
- `13.nozzle_H`: ノズル系の定常計算や、別系統の圧縮性流れで追加確認したいときに使う。
- `09.Taylor-Green`: **低散逸スキーム (純粋 KEEP) の運動エネルギー・エントロピー保存性**と**周期境界の保存性 (全運動量保存)** を見る非定常ケース。対流スキーム (特に KEEP)・周期境界・node/cell 共有コードの変更時に使う。

## 個別解説

- `20.naca_ml`: `procedures/verification/20-naca-ml.md`
- `08.bump`: `procedures/verification/08-bump.md`
- `13.nozzle_H`: `procedures/verification/13-nozzle-h.md`
- `09.Taylor-Green`: `procedures/verification/09-taylor-green.md` (KEEP の KE・エントロピー保存性)

## 運用ルール

- **node / cell の両方で検証する (必須)**: `mesh.discretization` には `cell` (既定) と `node` (median-dual) の 2 系統があり、対流流束・境界処理・主ループ対象 plane などが系統で分岐する。**両系統が通る変更 (対流/拡散/勾配/境界/時間積分などの共有コード) を検証するときは、cell ベースと node ベースの双方で実行し、双方の収束 (`check_convergence.py`) と場の妥当性を確認する。** 片方だけで「問題なし」と判断しない。
  - 同一ケースを `discretization: "cell"` と `discretization: "node"` の 2 run で回して比較する。
  - node 系統が未対応の構成 (3D median-dual 未実装 / periodic 未検証 / 一部ケースの近壁チェッカーボード発散など) でそのケースを node で回せない場合は、**その旨と理由を明記**したうえで cell のみで判断してよい (黙って片方を省略しない)。
  - **三者比較が必要なときは SU2 を基準に加える**: node / cell / SU2 をサンプリングライン (既定は下端から ~25% 高さ) 上の `P`/`T`/速度で突き合わせる。手順は [`procedures/su2-cross-check.md`](../su2-cross-check.md) の「ライン比較プロトコル」に従う。
- `run_*` の複製実行と `residual_history.png` の生成は、[AGENTS.md](../../AGENTS.md) の共通ルールに従う (検証実行でも同じ)。
- 生成される `residual_history.csv`、`res_*.h5`、壁面出力、実行ログは、その新しい `run_*` に保存して比較可能な状態を維持する。
- 既定ケースで十分に異常を切り分けられない場合は、変更内容に応じて他の標準ケースも追加で使う。
- 数値結果を変えないはずの変更では、既存ベースラインとの差分が不必要に増えていないことを確認する。
- アルゴリズム変更などで差分が意図的に出る場合は、不一致を失敗と決めつけず、どの量がどれだけ変わったかを報告する。
