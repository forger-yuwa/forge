# 出力量の絞り込み (`output.level`) と全エンタルピー `h0` の出力

## メタ

- **area**: `architecture / output`
- **status**: `done` (2026-09-08 実装・検証済。後続 2 件は残作業表)
- **related_docs**: `procedures/solver-settings.md` (output), `AGENTS.md` (出力と後処理の原則)
- **related_plans**: `accepted/turbulence-sst-energy-includes-k.md` (k 込み全温の後処理ルールの出所)
- **created**: `2026-09-08`
- **owner**: `CFD Dev`

## 1. 目的

`res_*.h5` が ~66 配列 (2.08M 節点で 680 MB/ファイル) を毎回書いており、うち勾配 27 本・リミッタ 5 本・SST/DES 診断 20 本超は
リポジトリ内のどのスクリプトも読んでいない (2026-09-08 調査: 参照 0 件)。後処理で導出できる量はソルバから出さない原則にし、
既定を「保存量 + 原始量 + h0」に絞る。あわせて全温・全圧の後処理が k を含めるかどうかで迷わないよう、
ソルバが全エンタルピー `h0` (sstEnergyIncludesK のときだけ +k、属性 `h0_includes_k`) を書く。

## 2. スコープ

- **やる**: config `output: {level, extraFields}`、`output.cpp` の絞り込みと `h0` 合成出力、docs/AGENTS.md のルール。
- **やらない**: 既存スクリプトの `T + u²/2c_p` 置換 (h0 逆算の共通関数 `tools/total_quantities.py` は後続)。境界面出力 (bvar) の絞り込み。

## 3. 関連 docs と前提

`variables.hpp` の `output_cellValNames` が出力候補、`output.cpp` の `writeSolutionH5_XDMF` が書く。スクリプトの読取実態:
ro/P/Ux/Uy/Uz/T/sonic が 30〜70 件、roU*/roe/roK/roOmega がリスタート・IC ツール、vis_lam/wall_dist が 15 件 (cf_node, check_quasisteady 等)、
k/omega/vis_turb が 3〜6 件、volume 5 件 (case/09 保存則)、cfl/ducros 各 1 件、他は 0 件。

## 4. 設計方針

| level | 出力 |
|---|---|
| 0 | `ro roUx roUy roUz roe roK roOmega roY*` + 凝縮モーメント保存量 (リスタート最小) |
| 1 (既定) | 0 + `P T Ux Uy Uz k omega sonic Y*` + `vis_lam vis_turb wall_dist` + `h0` |
| 2 | 従来 (全部) |

`extraFields` で個別追加。順序は `output_cellValNames` のまま。`h0 = Ht (+ k)` は出力時に host で合成 (Ht/k を D2H)。
XMF にも同じ名前を書く。D2H コピーも絞るので出力時間も減る。

## 5. 実装ステップ

1. `solverConfig.{hpp,cpp}`: `outputLevel`, `outputExtraFields` (済)。
2. `output/output.cpp`: `effectiveOutputNames`, h0 合成, XMF (済)。
3. docs: `procedures/solver-settings.md` output 節, `AGENTS.md` (済)。
4. 検証 (§6) → commit。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~検証 (§6)~~ | 済: case/16 `run_0313`–`0316` (level 1: 23 データセット 5.2 MB / level 2: 71, 13.2 MB / level 0+extra: 11 / energyK で `h0_includes_k`=1)。`check_quasisteady` 動作。3D 2.08M 換算 level 1 ≈290 MB (level 2 681 MB; メッシュ部 98 MB) |
| 2 | ~~`tools/total_quantities.py` (h0 → T0/P0)~~ | 済 (2026-09-08): CPG は閉形式 (平板 run_0028 で T+u²/2c_p と 7e-5 K 一致)、TP は species_db の NASA-9 (thermoHrefTemp datum) を Newton で逆算・P0 は s° 差 (case/16 run_0313/0316 で h0_includes_k 0/1 とも動作)。凝縮は凍結組成の気相逆算のみ (潜熱項未対応・警告) |
| 3 | ~~case/09 の保存則スクリプトの volume 読み替え~~ | 済: 4 本に `_cell_volume()` (res に無ければ solverConfig の meshFileName の `CELLS/volume`)。case/36 `analyze_ducros.py` は既存 run 専用で、新 run は `extraFields: [ducros]` が要る旨を注記。case/45 `sweep_cfl_implicit.py` の `cfl` は h5 フィールドでなく設定値 (対応不要) |

## 6. 検証

- level 2 のファイルが従来と同じデータセット集合であること。level 1 のデータセット数・サイズ、`h0` と属性の存在。
- 既存ツール (`check_convergence.py`, `check_quasisteady.py`, `cf_node.py`, `restart_field.py` / index コピー) が level 1 の res で動くこと。
- `sstEnergyIncludesK: 0/1` で `h0_includes_k` 属性が 0/1、値が Ht (+k) と一致。

## 7. 影響範囲

- 既定 level 1 では `dt_local cfl ducros volume limiter_* d*d* Pk_diag Taw_diag res_roK res_roOmega src_jac_* transport_diag_* wf_* roK_wf delta_les l_des fd_shield rd_des fe_iddes thermCond axisym_divU` が出なくなる。
  読む解析は `output: {level: 2}` か `extraFields` を config に足す (case/09 保存則, `sweep_cfl_implicit.py`, `analyze_ducros.py`)。

## 8. 完了条件

- [x] 関連 docs 更新済み (solver-settings, AGENTS)
- [x] 実装・検証完了 (§6)
- [x] `status: done`、§9 に変更ログ
- [x] `plans/active/` → `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-08` — 初稿・実装 (ユーザ要請「後処理で出せるものはソルバから出さない」「全温・全圧の k 込みルールを浸透させる」)。
- `2026-09-08` — 検証 (case/16 run_0313–0316) 完了、accepted へ。残: `tools/total_quantities.py`、case/09 保存則スクリプトの volume 読み替え。
- `2026-09-08` — `tools/total_quantities.py` 追加、case/09 保存則スクリプトの volume 読み替え、docs 反映。残作業なし。
