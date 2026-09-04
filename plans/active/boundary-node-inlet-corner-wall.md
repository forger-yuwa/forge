# node 入口∩壁コーナーの半割面所有 (nodeInletCornerWall)

## メタ

- **area**: `boundary / discretization`
- **status**: `in_progress`
- **related_docs**:
  - `methods/discretization.md` §7.2 (D)
  - `procedures/solver-settings.md` (`mesh.nodeInletCornerWall`)
- **related_plans**: `chemistry-finite-rate-h2.md` (Phase 3 Burrows–Kurkov で発覚), `discretization-node-boundary-ghostless.md`
- **created**: `2026-09-04`
- **owner**: `CFD Dev`

## 1. 目的

node モードで入口と no-slip 壁を共有する角ノードの質量蓄積 → 圧力暴走 (case/47 run_0002–0005, 0009–0012) を根治し、
上流壁を slip に逃がさずに壁 BL 付き入口 (`inletProfile`) を使えるようにする。

## 2. スコープ

- **やる**: `buildMedianDual` (2D) で `inlet_*` 境界エッジの壁ノード側半割面を壁 bcond に帰属 (`mesh.nodeInletCornerWall: 1`, 変換時)。
- **やらない**: 3D (`buildMedianDual3D`) は未対応 (同じ問題があれば同型で追加)。出口∩壁 (流入なし)、solver 実行時の切替。

## 3. 関連 docs と前提

`methods/discretization.md` §7.2 (A) 壁優先所有 → マルチマーカ (ow=ib) への変更経緯 (§7.0 の軸∩入口角のため)。本 plan はマルチマーカを維持しつつ
**壁ノード × 入口** の組合せだけ壁帰属に戻す。

## 4. 設計方針

`gmshReader::inletCornerWall` (converter が `cfg.nodeInletCornerWall` から設定)。壁ノード集合は wall/wall_isothermal の iPlanes から。
半割面ベクトル・面積加重重心とも壁 bcond に合算するので閉性チェックは不変。

## 5. 実装ステップ

1. `mesh/gmshReader.hpp` (所有ロジック), `input/solverConfig.*` (キー), `mesh/convertGmshToForge.cpp` (受け渡し) — 済 (2026-09-04)
2. case/47 で全壁 no-slip + 入口 BL プロファイル (`run_0014` 以降) が 2 次・cfl 2 で安定することを確認 — 進行中

## 6. 検証

- case/47 run_0011 系 (全壁 no-slip + profile) の再現: 変換時 `nodeInletCornerWall: 1` で 2 次・cfl 2 が NaN なしで 20000 step 完走、角ノード P が 1 atm 近傍。
- 回帰: `nodeInletCornerWall: 0` (既定) で変換したメッシュは従来とビット一致 (case/35 run_0050 相当)。

## 7. 影響範囲

`gmshReader.hpp`, `solverConfig.*`, `convertGmshToForge.cpp`。既定 0 で無影響。

## 8. 完了条件

- [x] methods §7.2 (D) 記述
- [ ] case/47 検証
- [ ] status done → accepted へ移動

## 9. 変更ログ

- `2026-09-04` — 初稿・実装。
