# 引き継ぎ: ⑤ SERN 設計チェーン (branch `feature/sern-design`, case/46)

- **日付**: 2026-09-05 / **branch**: `feature/sern-design` (origin に push 済み、最新は `git log` で確認)
- **正本**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md) (方針・変更ログ)、
  [`case/46.sern_design/README.md`](../../case/46.sern_design/README.md) (run 一覧と各節のまとめ)、
  [`notes/investigations/sern-design-method-survey.md`](../investigations/sern-design-method-survey.md) (方針の根拠)、
  結果ページ (図はクリック拡大): https://claude.ai/code/artifact/cb4ff402-3ea8-4c43-a131-f9d5b8fe6254
- **メモリ**: `sern-chain-status`, `interp-field-coincident-nodes-trap`, `matplotlib-japanese-font-noto-cjk`

## 1. 何ができているか (S0–S7 の状態)

| 段 | 状態 | 実体 |
| --- | --- | --- |
| 方針 | 壁圧 Bézier dv を撤回。**平面最大推力理論の key point (M_c, f) を dv にした逆 MOC** | survey §3–4, plan §4 |
| S1 逆設計 | 済 (解析検算 19/19)。kernel は key point で打ち切り | `design/forge_design/geometry/moc_sern.py`, `rao_planar.py`, `sern_geometry.py`, `design/tests/run_sern_moc_tests.py` |
| S2 メッシュ 2D | 済 (2 バンド構造 quad、カウルはスリット = 2 重ノード、TE 共有、板厚オプション) | `meshing/mesh_sern.py`, `tests/run_sern_mesh_tests.py` |
| S3 評価 2D | 済 (node/cell、Euler/SST、作動点 `--op`、δ\* オフセット、stage 間 restart は index コピー) | `evaluate/runner_sern.py`, `metrics/sern_forces.py` |
| S4 照合 | 済: Euler vs MOC (C_T 0.05 %)、NASA TM X-71972 の傾向 (カウル長・角) 再現 | run_0002–0007 |
| S5 δ\* | 済: **Euler 基準**の質量流束欠損で抽出 → 一発オフセット (揚力・モーメントが設計値に戻る) | `metrics/sern_deltastar.py`, run_0011/0013 |
| S6 MOO | 済: `opt/driver_sern.py` (LHS→Kriging→NSGA-II+EHVI、多作動点束ね、C_M 制約、剥離指標) | run_0010 (Euler, バグ入り), run_0017 (SST, バグ入り), **run_0019 (SST 3 作動点, 設計点固定, 正)** |
| S7 3D | 一部: hex 生成・runner・quad 力積分は済。**外側空間なし = 2D と一致 (run_0023)**。**外側空間あり: Euler 成立 (run_0027)、SST は後縁 3 重点で ω 発散 (run_0028) 未解決** | `meshing/mesh_sern3d.py`, `evaluate/runner_sern3d.py`, `case/46/plot_sern3d.py` |

主要な数値 (スモーク設計: M_in 2.5, ランプ 15°, カウル 5°/1H, M_c 3.9, f 0.45, p_ext/p_in 0.05):
MOC C_T 0.9666 / Euler cell 0.9660 / node SST C_T(p) 0.9691, 摩擦 −0.0065 / 3D 押し出し SST 0.9699 /
3D 外側あり Euler (加速点 M∞3.5) C_T 0.932 vs 2D 0.976 (横方向膨張でランプ圧 1/3、揚力反転)。

## 2. 使い方 (design/ で)

```bash
PY=.venv-opt/bin/python
# 単体テスト
PYTHONPATH=. $PY tests/run_sern_moc_tests.py ; PYTHONPATH=. $PY tests/run_sern_mesh_tests.py
# 2D 1 作動点 (段階起動 soft→本段、force は metrics.json)
PYTHONPATH=. $PY -m forge_design.evaluate.runner_sern ../case/46.sern_design/problem_smoke_sst_node_t0.yaml ../case/46.sern_design/run_NNNN_<slug> [--op accel] [--wall-offset off.json]
# MOO (多作動点; spec.operating_points と opt: を YAML で)
PYTHONPATH=. $PY -m forge_design.opt.driver_sern ../case/46.sern_design/problem_moo_sst_node_3op.yaml ../case/46.sern_design/run_NNNN_moo --n-doe 10 --n-iter 2 --batch 2
# 3D (mesh3d: 節が要る; op 指定可)
PYTHONPATH=. $PY -m forge_design.evaluate.runner_sern3d ../case/46.sern_design/problem_smoke_euler_node_3d_accel_slip.yaml ../case/46.sern_design/run_NNNN_3d --op accel --soft-steps 1500
# 図
$PY ../case/46.sern_design/plot_sern_run.py <run_dir> <problem.yaml>      # 2D: Mach 場 + 壁圧 (CFD vs MOC)
$PY ../case/46.sern_design/plot_pareto.py <campaign_dir> "<title>"          # MOO パレート
$PY ../case/46.sern_design/plot_sern3d.py <run_dir> <problem.yaml>         # 3D: ランプ圧力マップ + 対称面 Mach
```

- YAML の形: `case/46.sern_design/problem_*.yaml` を複製する。`spec.external` = **設計点** (逆設計はこれで固定)、
  `spec.operating_points[]` = CFD 条件 (external / inflow の上書き、weight)。`evaluate.model: euler|sst`、
  `mesh.discretization: node|cell`、`evaluate.top_out_kind: outlet|slip`、`geometry.mode: keypoint|straight`。
- 新 run は必ず `run_NNNN_<slug>` を新規作成し、README の run 一覧表に行を足す。

## 3. 踏んだ罠 (再発防止)

1. **座標一致ノード (スリット) と最近傍補間**: `interp_field.py` は双子ノードを同じ元に写す → stage 間は index コピー
   (`runner_sern.restart_by_index`, `runner_sern3d.warm_from_same_mesh` は状態量のみ。wall_dist を上書きしない)。
2. **twall の符号**: cell = 流体に働く力、node = 壁に働く力 (`sern_forces` で離散化ごとに切替)。
3. **作動点ごとの再設計バグ** (修正済): 逆設計は `spec.external` で固定。run_0010/0017 は作動点間で形状不一致。
4. **δ\* の単独抽出**は膨張扇を欠損と誤認 (20 倍過大) → Euler 基準 (`deltastar_sern_vs_euler`) を使う。
5. **3D**: 入口面のノズル幅外は `inlet_ext` / ランプ∩側壁ノードは内外 2 重 (2 種入口の共有ノード禁止) / index IC で
   カウル外面ノードは外部流 / ランプ延長 `top_out` は slip 推奨 (静圧出口は流れ平行で不安定) / M∞6 では幅外ランプ角部の
   真空膨張で 2 次が落ちる (加速点で評価)。壁面出力 CONNE は 3D で [5, n0..n3]、2D で [2, 2, a, b]。
6. **node の rms_roOmega** は本段で 3e18 一定になる (壁ノード ω ピン留めの診断値)。収束判定は力係数の STEADY と NaN なしで。
7. **図**: matplotlib は `~/.fonts/NotoSansCJKjp-Regular.otf` を addfont。Droid Sans Fallback は豆腐。上付き C⁺/C⁻ は ASCII に。
8. **pkill -f** に自分のコマンド文字列が引っかかる (自殺)。パターンは分割して書く。

## 4. 未解決・次にやること (優先順)

> **正本は plan の §5.1 残作業表** ([`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md))。
> 以下は 2026-09-05 時点の写しで、項目 1 は §4.10 で方針決着済み (実装未着手)。優先順が変わったら plan 側を直す。

1. **作動点をサイクル値で定義し直す** — **方針決着済 (plan §4.10)**: アンカー = NASA TM X-71972 TABLE 1 + CEA2。
   M∞6 巡航は M_in 1.67 / p_in 101 kPa / NPR 35.4 / γ 1.18 (初版 2.5 / 20 kPa / 20 / 1.4)。低 NPR の正体は燃料遮断。
   残実装 = `operating_points[].gas` 上書き + YAML 差し替え + MOO 再取得。CEA 入出力は `case/46.sern_design/cea/`。
2. **剥離評価にはランプ後縁側の外部流 (機体上面からの回り込み・ベース) が要る**: 今は後縁の先を静圧出口/slip で閉じている
   ので RSS/FSS は信用できない。2D にブロックを 1 つ足す (mesh_sern に 3 バンド目) のが次の作業。
3. **3D SST の後縁 3 重点** (カウル後縁∩側壁後縁で ω → inf): 候補 = 後縁に板厚 (2 本の剥離線を分離)、
   `FORGE_FREEZE_TURB=1` で soft 段を凍結乱流で通す、局所の ω 上限。Euler は同格子で成立 (run_0027)。
4. **側壁長 L_sw を dv/仕様に**: 3D 効果 (推力 −4.5 %, 揚力反転) の主因。
5. **亜音速外部流 (静止・地上試験) のチャンバー構成**: case/23 方式。超音速外部流なら不要。
6. **カウルの輪郭設計** (Lv & Xu 2021 はカウルにも最適化理論): 今は直線固定。
7. **音速給気 (ラムジェット/TBCC) のスロート接続** (`sonic_throat` モード、Hall + kernel MOC で M 1.3–1.5 まで対称内部ノズル)。
8. 非一様入口 (回転流 MOC)、凍結 γ → 有限速度化学との接続。
9. CFD 領域の下側深さ (`mesh.bot_depth` 3H) は超音速外部流なら 1H 程度に縮める。

## 5. ディスク

`case/46.sern_design/` は run 30 本で数 GB。破棄予定 (README の表で印) は削除してよい: run_0001, 0008, 0009, 0012, 0014, 0020,
0021, 0022, 0025, 0026, 0028。MOO の各点 run (`run_0010/0017/0019/*_<op>/`) は最終 res だけ残してある。
