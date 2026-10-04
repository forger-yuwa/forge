# 引き継ぎ: 壁隣接面の SLAU $\chi$ (`slauWallNormalChi`) と SERN 出口 BC (2026-09-23)

**正本は plan 2 本**。本文書は写しとポインタ。残作業を増やすときは plan §5.1 に書く。

- [`plans/accepted/convection-slau-wall-normal-chi.md`](../../plans/accepted/convection-slau-wall-normal-chi.md) — 対策本体 (`in_progress`)
- [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md) §4.13.1 / §5.1 B1d–B1f — 接続模型の診断

## 1. 何が起きていたか (1 分で)

⑤ SERN の接続模型 (全ヘキサ・有限厚の板、245 万節点) が起動しなかった。**独立した 2 つの欠陥**があった。

**① 壁節点の痩せ**: node の壁節点は $u=0$ に固定され移流の流入経路が無く、SLAU の圧力差項だけが補充源。
その項を切る $\chi$ の判定が**速度の大きさ**基準なので、隣の内点の**接線**速度が大きいと $\chi=0$ になり補充が消える。
実測: 壁 58 Pa 対 内点 3545 Pa の 61 倍差でも質量が戻らず、床を割って NaN。

**② 超音速出口への静圧指定**: `outlet_statPress` が node の壁列・後流の**亜音速ノード**に背圧を課し、そこから圧力が育つ。
**同じ case で 2 回目** (run_0121 で一度踏み、対策が個別 YAML のコメントに留まって既定へ反映されなかった)。

## 2. 対策と現状

| | 内容 | 状態 |
| --- | --- | --- |
| ① | `space.slauWallNormalChi` (opt-in, **既定 0**)。壁隣接面の**質量流束の $\chi$ だけ**を面法線 Mach で組む | 実装済 (`5a4886ea`)。**検証は残あり** |
| ② | `outlet_kind` の既定を `outflow` に (SERN の 3 runner)。`procedures/recommended-settings.md` に追記 | **完了** (`48aefa15`) |

## 3. 2×2 の結果 (効果の分離)

| | `outlet_statPress` | `outflow` |
| --- | --- | --- |
| **flag 0** | `run_0432`: step 815 で NaN | `run_0440`: **step 815 で NaN** (排出の軌跡が 5 桁一致) |
| **flag 1** | `run_0435`: 完走するが出口が破綻 (残差 +0.8 桁) | `run_0437–0439`: **通算 66000 完走・`ALL STEADY`** |

→ **2 つは独立。両方の手当てが要る。**

## 4. 検証の到達点

| 試験 | 結果 |
| --- | --- |
| V0 単体 | **ALL PASS** (`cad/test_diag_wall_cv_budget.py`)。カーネル照合・面反転・等状態・$\Delta p$ 比例 |
| V0 config | 誤用 4 経路すべて起動時に停止 |
| V2 (flag 0 の無害性) | **PASS** — ノイズ床比 0.86〜1.36 |
| V1-a/e | step 8–13 で $10\rho_{Min}$ 超 (期限 200) |
| V1-b/c | **PASS** — 観測/予測 **1.000** (3 点)、$\Sigma\dot m$ が 5 桁縮小 |
| V1-d | **`ALL STEADY`** (通算 66000 step) |
| V3 冷却平板 | $C_f$ 0.000 % / $q_w$ 0.011 % |
| V3 SERN 2D | $C_T$ **0.031 %** (許容 ±0.1 %) |
| V3 衝撃足 | 壁圧 L2 差 **0.842 %** (許容 1 %)、位置ずれ 0 |
| CFL 固定点 | 累積 CFL を揃えて **0.01 % 一致** |

**未実施 (§5.1 #9/#10)**: 格子感度 (3 水準)、case/16 の V2・V3。どちらも**メッシュ生成から**で数時間規模。

## 5. 既定化の条件 (§5.1 #11)

**opt-in のまま**にしてある。既定にするには次が要る。

- #7 周期・軸対称の小規模試験、#9 格子感度、#10 case/16
- **$C_L$ 0.433 % / $C_M$ 0.316 % / 衝撃足の壁圧 0.842 % の正否を独立に判定** (格子収束・別スキーム・実験/文献のいずれか)。
  これらは**ノイズではない** (各側 3 本の反復で床 0.005–0.07 %、差は床の 5–10 倍)。

## 6. 今日踏んだ罠 (同じ型が 3 回出た)

**「比較するものが揃っているか」の確認不足**:

1. **診断ツールが変更前の $\chi$ を計算していた** (codex result M1)。「予測と一致して PASS」が無効になった。
   → ツールは**カーネルと項ごとに一致する単体試験**で縛る。`--wall-normal-chi` の指定忘れも同じ事故になる。
2. **効果の分離を 1 点比較で済ませた** (M2)。2 変数なら **2×2 を揃える**。
3. **CFL 固定点の比較で収束途中の場を並べた**。8〜16 % の差が出たが、累積 CFL を揃えたら 0.01 % に収まった。
   → **反復数でなく累積 CFL を揃える**。

その他:
- `outlet_kind` の既定は SERN 3 runner のみ。**波及先は case/46 だけ** (`PHYS_SERN` 専用、風洞系は `runner_wt.py` の自前 bcond)。
- XMF は**同名 `.h5` を相対参照**する。改名するとパラビューで開けない。ディレクトリで分ける。
- AWS → 手元の転送は約 0.6 MB/s。h5py の gzip で 340 MB → 145 MB (2.3 倍) に圧縮してから送る。
- ローカル GPU (RTX 3060) は A10G 向けのブロックサイズでカーネルが起動できない。`FORGE_CUDA_BLOCKSIZE=64` を両バイナリ共通で指定する。

## 7. 環境

- AWS `i-0b1a5e0b8dc152f00`。**今日 3 回自動停止した** (4 分間隔の ssh keepalive を張っても落ちる。原因未特定)。
- 作業ツリーは **`/home/sano/work/forge-sern-design`** (`feature/sern-design`)。
  `/home/sano/work/forge` は別セッションが別ブランチで使用中なので触らない。
- ParaView 用の場: `/home/sano/work/forge/case/46.sern_design/_paraview/` (3D outflow / 2D flag0,1 / 出口破綻時の比較)。
