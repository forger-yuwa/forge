---
name: forge-config
description: forge の solverConfig.yaml / bcondConfig.yaml を新規に組む・見直すときの手順。解析種別から現行レシピ (procedures/recommended-settings.md) を選び、廃止キー・旧既定を混入させずに段階起動と検証まで揃える。「設定どうする」「config 作って」「この run の設定は古い?」と言われたときに使う。
---

# /forge-config — 解析設定を現行レシピで組む

正本は [`procedures/recommended-settings.md`](../../../procedures/recommended-settings.md)。**記憶や過去 run の config を
そのままコピーしない** (廃止キー・旧既定が残っている)。手順:

1. **解析種別を決める** (§0 の表): 定常 NS/SST、Euler 設計評価、軸対称、TP/凝縮/化学、非定常 LES/DES、変換メッシュ。
2. **該当節の「現行 (日付)」ブロックから config を組む**。日付の無い記述と §9 (superseded) の設定は使わない。
   キーの意味に迷ったら `procedures/solver-settings.md`。
3. **境界条件** (§1.1): 出口の種類 (亜音速 statPress + 逆流 Pt/Tt / 超音速 outflow か Ps 一致)、壁 (no-slip 変換で wall_dist)、
   node の `nodeInletCornerWall`。
4. **段階起動** (§1.2): soft → mid → 本段。同一メッシュは index コピー、cross-mesh は `interp_field.py`。
   一様 IC から超音速/SST を直接始めない。初期 k/ω は非ゼロ。
5. **投入前チェック**: `check_mesh_quality.py` の VERDICT、バイナリ鮮度 (`find ... -newer build/forge`)、
   `output.level` (既定 1; 診断が要るときだけ 2)、新しい `run_NNNN_<slug>` ディレクトリ。
6. **報告**: run パス、`check_convergence` / `check_quasisteady` の VERDICT、case README の run 一覧。

既存 config の点検を頼まれたら、§9 の表と突き合わせて「廃止キー / 旧既定 / 非推奨の組み合わせ」を列挙し、
現行値への置換案を出す (勝手に書き換えない)。

レシピを変えたら `recommended-settings.md` の該当節の日付を更新し、旧値を §9 に移す (plan 反映と同じく「決めた時点で書く」)。
