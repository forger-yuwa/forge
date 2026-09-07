# 入口分布の与え方 (全温 / 全圧 / 組成 / k・ω / 超音速入口の ρ,U,Ps)

入口境界の値を**一様でなく座標の分布として**与える手順。対象は「超音速ノズルの亜音速入口に全温分布・
H2O 分布を与える」「超音速入口 (全量固定) に全温分布を与える」「入口境界層 (壁法則) を与える」など。
機構の仕様は [`methods/boundary.md`](../methods/boundary.md) の「入口分布プロファイル」、設計判断は
[`plans/accepted/boundary-inlet-profile.md`](../plans/accepted/boundary-inlet-profile.md)。
Claude は skill `forge-inlet-profile` でこの手順を辿る。

## 0. 何ができるか (2026-09-08 時点)

forge の入口カーネルは **面ごとの境界値 `bvar`** を読む。`inletProfile: 1` を付けた入口では起動時に
`inlet_profile_<physID>.csv` を face 重心座標で補間して `bvar` に書くので、**カーネル無改修で分布入口**になる。

| 入口種別 | 主な用途 | CSV に書ける列 (分布にできる量) | 備考 |
| --- | --- | --- | --- |
| `inlet_Pressure` | **亜音速入口** (ノズル上流のプレナム・整流部)。超音速ノズルでも入口は亜音速なのでこれ | `Tt` `Pt` `Y0..Y{n-1}` `k` `omega` | 超音速流入は非対応 (`Un_c<0` 分岐のみ)。Tt/Pt から内面マッハで静的状態を組む |
| `inlet_uniformVelocity` | **超音速入口** (全量固定: ρ, U, Ps)。亜音速では ρ/Ps を内面で上書き | `ro` `Ux` `Uy` `Uz` `Ps` `Y0..` `k` `omega` | **`Tt` 列は無効** (kernel が読まない)。全温分布は ρ/U/Ps に換算して与える (§2) |
| `inlet_fluctVelocity` | 合成乱流入口 | `Ux0` `Uy0` `Uz0` `ro` `Ps` `Y0..` `k` `omega` | 平均場のみ分布化 |
| `inlet_Pressure_dir` | 方向指定の全圧入口 | `Tt` `Pt` `k` `omega` | **多成分非対応**: カーネルが組成 `Yb` を受け取らず `sp[0]` 単成分の熱物性で状態を組む。`Y{s}` 列は化学種 Dirichlet だけに効き境界状態と不整合になるので、多成分では `inlet_Pressure` を使う |

- **組成 `Y{s}`** は `physProp.species` の並び順 (index) で列名にする (`Y1` = 2 番目の種)。指定した種以外は
  残り (1−ΣY) を受け持つ。**多成分 TP (`thermalMethod: 2`) のときだけ有効**。
- **k / ω** は per-face Dirichlet で入る (node は境界ノードをピン)。
- **凝縮モーメント** (液相) は入口で常に 0 (dry Dirichlet)。分布化は不可。
- **補間**: 座標列が 1 つ (`y` など) なら **その軸で 1D 線形補間** (範囲外は端値クランプ)。座標列が 2〜3 つ
  (`y z`, `x y z`) なら **最近傍** (2D 分布は表を十分細かくする)。
- **restart でも毎回読む** (`bvar` は res に保存されない)。**CSV は run ディレクトリごとに置く**。
  `inletProfile: 1` のまま CSV が無いと起動時にエラー終了する (黙って一様には戻らない)。段階起動で run を
  複製するときは CSV も一緒にコピーする。
- 既知の制約: 時間変化する分布は不可 (合成乱流は `inlet_fluctVelocity`)。1D 線形は座標 1 軸のみ。

## 1. 亜音速入口 (`inlet_Pressure`) に Tt / H2O 分布を与える — 標準手順

1. **run を複製**して新しい `run_NNNN_<slug>` を作る (IC は同一メッシュなら index コピー)。
2. `bcondConfig.yaml` の入口に `ints: {inletProfile: 1}` を付ける。`floats` の一様値はそのまま残す
   (CSV に書かない量はこの一様値のまま; 換算のデフォルトにも使う)。
3. **CSV を生成する** ([`tools/gen_inlet_profile.py`](../solver_density_cuda/tools/gen_inlet_profile.py)):

   ```bash
   # 式で与える (座標の式。使える変数: x y z r, 一様値 cfg_Tt cfg_Pt cfg_Y_H2O ..., numpy 関数)
   python3 solver_density_cuda/tools/gen_inlet_profile.py gen --run case/XX/run_NNNN --physID 1 \
       --axis y --range -0.0127 0.0127 --n 101 \
       --Tt "286.65 + 30*exp(-(y/0.004)**2)" \
       --Y  "H2O=0.005 + 0.012*exp(-(y/0.006)**2)" --plot
   # 測定表で与える (空白区切り, ヘッダ先頭に座標列; 列名は Tt Pt Y_H2O k omega ...)
   python3 solver_density_cuda/tools/gen_inlet_profile.py gen --run case/XX/run_NNNN --physID 1 --table measured.csv
   ```

   出力は `run/inlet_profile_<physID>.csv` (`--plot` で `inlet_profile_<physID>.png`)。gen は各列の min/max と
   bcond 一様値を表示するので、**ここで桁と符号を確認**する (Tt が K か、Y が 0〜1 か)。
   `--Ps` と `--M` を与えると等エントロピーで `Pt` に換算する (TP は NASA-9 の h/s° で解く)。
4. **起動ログで確認**: `forge_run.log` に
   `[applyInletProfiles] physID=1 kind=inlet_Pressure: set 3 quantities from inlet_profile_1.csv (1D interp, 101 rows, 118 faces). applied: Tt Y0 Y1`
   が出ること。`set N` は CSV ヘッダの量列数、**`applied:` が実際に境界値に反映された列**で、その種別に無い列名は
   `IGNORED (not an input quantity of this kind):` に出る (1 列も反映できなければエラー終了)。多成分では各 face の
   0≤Y_s≤1・ΣY_s=1 も検査され、外れるとエラー終了する。この行が無ければ `inletProfile: 1` が付いていない。
5. **結果と照合** (`verify`): 境界ノード (node) / 第 1 セル (cell) の場を CSV の目標と比べる。

   ```bash
   python3 solver_density_cuda/tools/gen_inlet_profile.py verify --run case/XX/run_NNNN --physID 1 --plot
   ```

   `Tt` は `VALUE/h0` 由来の T0 (`total_quantities.py`) と比較する。**cell は第 1 セルの値なので境界値そのものでは
   なく**、境界層内 (低速部) は初期組成が洗い流されるまで目標と差が残る (過渡)。
6. 通常の収束・準定常確認 (`check_convergence.py` / `check_quasisteady.py`) と case README の run 一覧更新。

### 検証済みの例 (2026-09-08)

| run | 内容 | 結果 |
| --- | --- | --- |
| `case/16.nozzle_wys/run_0317_inletprof_tt_h2o/` | node TP (MIXDRY+H2O) SST, `inlet_Pressure`, Tt +30 K Gaussian, Y_H2O 0.005+0.012 Gaussian, 1000 step | 入口 T0 は目標と max 0.3 K 差 (幅 30 K の 1 %)。Y_H2O は **node ピン修正後**に 1e-7 で目標一致 (修正前は一様値のまま = 本手順書の動機) |
| `case/16.nozzle_wys/run_0318_inletprof_h2o_cell/` | cell 版、Y_H2O 分布のみ, 600 step | 第 1 セルの Y_H2O が目標に追従 (境界層内は過渡で初期値との中間) |
| `case/46.sern_design/run_0089_inletprof_supersonic_tt/` | node CPG Euler, `inlet_uniformVelocity`, Tt 2025+200 K Gaussian・M 2.5・Ps 20 kPa を ρ/U/Ps に換算, 600 step | 入口ノードの ρ・U が目標と 1e-5 で一致 (換算の往復誤差 1e-12 K) |

## 2. 超音速入口 (`inlet_uniformVelocity`) に全温分布を与える

超音速入口は (ρ, U, Ps) を全量固定するので、**Tt(y) は静的状態に換算**して与える。`gen` に `--Tt --M` と
`--Ps` (または `--Pt`) を渡すと、CPG は閉形式、TP は NASA-9 で $h(T_t)-h(T_s)=\tfrac{1}{2}M^2 a(T_s)^2$ を解いて
`ro Ux Uy Uz Ps` 列を書く (方向は `--dir` か bcond の速度方向)。組成分布 (`--Y`) も同時に換算に入る。

```bash
python3 solver_density_cuda/tools/gen_inlet_profile.py gen --run case/46/run_NNNN --physID 1 \
    --axis y --range 0.0 0.1 --n 51 --Tt "2025 + 200*exp(-((y-0.05)/0.02)**2)" --M 2.5 --Ps 20000 --plot
```

gen は換算の自己検証 (出力から Tt, M を復元した誤差) を表示する。`--set NAME=EXPR` で `ro`/`Ux`/`Ps` を
直接書くこともできる。

## 3. 入口境界層 (壁法則) を与える

[`tools/gen_inlet_walllaw.py`](../solver_density_cuda/tools/gen_inlet_walllaw.py) が Reichardt 合成則で
チャネル/片側 BL の $u(y)$ を書く (`--physID --ylo --yhi --Uc --nu --mode channel|bl`)。速度入口
(`inlet_uniformVelocity`) 用。全圧入口で境界層を表したいときは `--Pt` の分布 (全圧欠損) を `gen` で与える。

## 4. つまずきどころ

- **ログの `IGNORED` に列が出る**: その入口種別に無い列名 (例 `inlet_uniformVelocity` に `Tt`) は反映されない
  (`bc.bvar` に存在する量だけ反映)。§0 の表の列名を使う。`set N` はヘッダの列数であって反映数ではない。
- **node で Y が一様のまま** (2026-09-08 以前のバイナリ): 化学種の入口 Dirichlet が ghost にしか書かれず
  境界ノードに効かなかった。修正済み (境界ノードをピン + 残差除外 `speciesPinResidual_d_wrapper` + 連成陰解法
  `speciesImplicitCoupling: 1/2` の Jacobi sweep でもピン行を δ=0 に拘束)。古いバイナリで回した run は
  Y 分布が入っていないので再実行する。
- **座標の単位**: CSV は変換後メッシュの座標 (m) で書く。Gmsh の ScalingFactor 後の値。
- **2D 分布 (y z)**: 最近傍なので、表の刻みは面の刻みより細かくする。
- **段階起動**: soft → mid → 本段で run を複製するたびに CSV を一緒にコピーする。`CONTINUED_FROM.txt` に
  「inletProfile CSV 継承」と書いておく。
- **物理的な不整合**: 入口で Tt を上げると同じ Pt でも密度が下がり質量流量が変わる。設計比較では ṁ を併記する。
