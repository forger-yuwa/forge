# 引き継ぎ: ⑤ SERN 全ヘキサ接続模型の CFD 起動 (2026-09-22, session_014aoVV4K167qqfTw7TojyJc → 次セッション)

**正本は plan** [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md) (§4.10–§4.13、§5.1 の B1b/B1d、§6.1)。
本文書は写しとポインタ。残作業を増やすときは plan §5.1 に書く (ここには書かない)。

## 1. いま何がどこまでか (1 分で) — **2026-09-22 夜 更新**

- **格子は完成** (plan §4.11–§4.12): 全ヘキサ接続模型 `case/46.sern_design/cad/hex_junction_model.py`。
  採用設定 → 2,446,572 節点、`check_mesh_quality.py` PASS (AR max 676)。`cad/demo_out/jm_c16.msh` (330 MB, git 追跡外。AWS にもある)。
- **B1d ①②③ は完了** (plan **§4.13.1** が正本)。run_0430 (IC 修正) / run_0431 (断熱壁 A/B) / run_0432 (判別 A/B, 継続 + 毎 step プローブ)。
  - **IC 修正は発散時刻を動かすだけ**で排出は残る (run_0427 step 1720 NaN → run_0430 は 1800 完走だが打ち切りが早かっただけ、run_0432 で step 2615 NaN)。
  - **熱的拘束は無関係** (断熱壁でも同率で抜ける)。**閉性 (B1b) も無関係** (当該 CV で 4.2e-6)。
  - **面別 SLAU 収支**: 流出面は $\chi$=0 で圧力差項が消える W↔I 面。補充源は壁↔壁の $\chi$=1 面で、**壁列を伝って枯れる**。
  - **判別 A/B = 分岐 B (構造的)**: 内点の $V_n$ が 2600 step 不変 (276→287 m/s)。**起動側の工夫では消えない**。
  - **副産物**: **EOS 床が定常陰解法の commit で捨てられる** (`update_d.cu:248` が床前の `roN` 起点)。
    床後は流束が $\rho_f=\rho_{Min}$ 固定で評価され、**指数減衰が一定シンクに変わり有限 step で負密度に達する**
    (パラメータ・$\Delta t$ 非依存の検証で予測/実測 0.93–1.07)。plan §5.1 **B1e** に独立の欠陥として起票。
- **まだ書いていないこと**: 「機序確定」。テーパは候補外 (codex M5 と `diagnostician` の両方が同じ結論)。

## 2. 最初にやること

1. **AWS を起動してもらう** (`i-0b1a5e0b8dc152f00`。stop/start で IP が変わるのでユーザに聞く。手元に AWS 資格情報は無い)。
   接続したら `sudo swapon /swapfile; sudo swapon /swapfile2` と `cd ~/forge-sern && git pull --ff-only`。
   **自動停止に注意** (2026-09-22 は `ssh … 'sleep 7200' &` を張っていても落ちた。長い解析の前にまず短い確認を済ませる)。
   ディスクは 2026-09-22 時点で約 10 GB 空き (run_0423–0429 の `sern.h5`/`sern.msh` は削除済。config・log・残差履歴は残してある)。
2. **積み残しの確認 (数分)**: `res_nan_815.h5` で最初に負・NaN になった**節点 ID**。
   一定シンクなら先に 0 を切るのは $\rho$ が最小の **153880** で、**153797 ではない**。plan §4.13.1 の「未了」。
3. **次の run は Roe の判別** (下の §3)。

### 投入コマンド (AWS 上) — SLAU → ROE の判別 A/B

```bash
cd ~/forge-sern/case/46.sern_design
mkdir run_0433_3d_junction_diag_roe && cd run_0433_3d_junction_diag_roe
for f in sern.h5 bcondConfig.yaml probe.yaml species_db.yaml species_meta.yaml solverConfig_main.yaml prepare_info.json MESH_QUALITY.txt; do
  cp ../run_0432_3d_junction_diag_ext4k/$f .; done
sed -E 's/solver: "SLAU"/solver: "ROE"/' ../run_0432_3d_junction_diag_ext4k/solverConfig.yaml > solverConfig.yaml
diff <(sed -E 's/solver: "[A-Z]+"//' solverConfig.yaml) <(sed -E 's/solver: "[A-Z]+"//' ../run_0432_3d_junction_diag_ext4k/solverConfig.yaml) \
  && echo "solver 以外は run_0432 と一致"
cd ~/forge-sern
python3 -c 'import sys; sys.path.insert(0,"design")
from forge_design.evaluate import runner_sern as R2
R2.restart_by_index("case/46.sern_design/run_0430_3d_junction_diag_icfix/res_1800.h5",
                    "case/46.sern_design/run_0433_3d_junction_diag_roe/sern.h5")'
solver_density_cuda/tools/run_case.sh case/46.sern_design/run_0433_3d_junction_diag_roe
```

`run_0432` の `sern.h5` は既に res_1800 の場で上書きされている (restart 済) ので、**コピー後に必ず restart_by_index をやり直す**
(でないと run_0432 の初期場と同じに見えて実は同じ、という混同が起きる。上のコマンドは run_0430 の res_1800 から入れ直している)。

**判定** (plan §4.13.1 の「次手」): **床前区間** ($\rho_w > \rho_{Min}$ ⇔ probe の $P$ > 33.8 Pa) の $\rho_w$ 減衰率で見る。**NaN の step は使わない**。
- A: 減衰が止まる/大幅に鈍る → $\chi$=0 が本質 → **(ii) 壁隣接面の $\chi$ を面法線 Mach で評価**を plan §4 に起票
  (`cuda_forge/` の流束カーネル変更 = AGENTS.md 条件 6。**plan §4 起票 + codex plan 段が先**。実装から入らない)。
- B: 同率で抜ける → **(i) 隅の $\Delta x_1$ を粘性層尺度に振る** ($t_{base}$ 固定) か Dirichlet 側の設計を起票。

## 3. 解析の道具

**`case/46.sern_design/cad/diag_wall_cv_budget.py`** (commit `41508544`)。`PLANES/STRUCT`・`surfVect`・`CELLS/volume` と dump の原始量から、
1 次 SLAU の**面別質量流束**を項ごとに (移流項 / 圧力差項 $-\chi(P_R-P_L)/\hat c$ / $\chi$ / $\hat c$ / 両側 $V_n$) 再計算する。
式は `convectiveFlux_slau_d.inc.cuh:533-564` と項ごとに照合済み。

```bash
python3 ../cad/diag_wall_cv_budget.py sern.h5 --summary                 # dump ごとの ro/P 最小値・位置・床到達数
python3 ../cad/diag_wall_cv_budget.py sern.h5 --auto 3 --faces          # 最終 dump で ro 最小の 3 節点と面別内訳
python3 ../cad/diag_wall_cv_budget.py sern.h5 --ids 153797 --faces-all  # 節点 ID 固定で全 dump の内訳
```

**注意 (§4.13.1 の「まだ書かないこと」)**: この収支が**ソルバの収支と一致するのは $\rho\ge\rho_{Min}$ かつ $\rho R T_w \ge p_{Min}$ の間だけ**。
それ以降は dump (commit 後・pin 後) とソルバが使う状態 (pin → 床・クランプ後) が別物になる。
**ソルバ前処理を足した再収支は未実装** (plan §5.1 B1d の「未了」)。

**毎 step の時系列は probe が速い** (全場 dump は 1 本 325 MB)。`probe.yaml` に `outStepInterval: 1` と
`points: {名前: {x:, y:, z:}}` を書くと `point_probe_<i>.out` に `Step,TotalTime,T,P,Ux,Uy,Uz` が出る。
座標は**双対重心** (`CELLS/centCoords`) を与える。壁節点は $T$ がピンされるので $\rho = P/(R T_w)$ で読める ($R$ は dump から較正)。

## 4. 判断が要ること (ユーザ or Fable)

- **テーパは候補外になった** (codex M5 と `diagnostician` が同じ結論: 3 重隅が残り、$t_{base}$ と $\Delta x_1$ が交絡する)。
  残る二択は plan §5.1 B1d の (i) 隅の $\Delta x_1$ / (ii) 壁隣接面の $\chi$ の定義で、**Roe の判別 A/B の結果で選ぶ**。
- **(ii) に進む場合は `cuda_forge/` の流束カーネル変更** = AGENTS.md 条件 6。plan §4 起票 + codex plan 段 + `diagnostician` が先。
- **B1e (EOS 床が定常 commit で捨てられる)** をいつ直すか。**全ケースに効く**ので既定変更は慎重に。
  まず既存 opt-in `time.deltaT.updateGuardAlpha` (既定 0) で足りるかを確認する。
- 変換器の単精度 (B1b) をいつ直すか (別セッションの変換器作業と衝突しないよう着手前に確認)。
  **今回の排出の原因ではない**ことは確認済み (当該 CV の閉性 4.2e-6)。

## 5. 既知の罠 (今日踏んだもの)

- AWS の自動停止 (§2)。AWS の `~/forge-sern` は shallow clone、`solver_density_cuda/build` は `.build-native/relwithdebinfo` へのリンク。gmsh は未導入 (格子は手元で切って scp)。
- `check_mesh_quality.py` の六面体の体積判定は今日 頂点 Jacobian に変えた (`3c73ace8`)。旧判定は曲面壁の薄いセルを誤って FATAL にしていた。
- 共有ワークツリー: 別セッション (case/51・53・54 の CHT/遷移) が同じブランチで作業中。**commit は必ず `git commit -- <paths>`** (index 経由だと相手の staged を巻き込む。今日 1 回やった)。
- `centCoords` は双対重心。壁節点の位置判定は `MESH/COORD` を使う。
- run_junction_model.py の IC は接続模型の寸法 (H, ZW, TSW, TC, LSW, LCOWL) を定数で持つ。形状を変えたら合わせる (codex M3 の注意)。
- **`pgrep -f "build/forge"` の待ちループは自分自身にマッチする** (bash -c のコマンド行にその文字列が入るため無限ループ)。
  待つなら `pgrep -f "[b]uild/forge"` にするか、親側で `run_in_background` を使う。
- **暖機段の config は `run_0427` の `solverConfig.yaml` をコピーすれば足りる** (`run_staged` が書くものそのもの)。
  sed で組み直すと数値設定が微妙に変わって IC 以外の差分が入る (codex m6 の「比較対象を固定する」)。
- **forge バイナリは意図的に `c21fdd19` ビルドのまま使った** (sha256 `eec7f767…`)。run_0427 と揃えるため。
  別セッションの `transition_d.cu` 変更を取り込まずに済む (`turbulence: none` なので不活性)。再ビルドするなら比較の基準を引き直すこと。

## 6. 今日の commit (feature/sern-design)

`645f25a4` B0 設計 → `e580114c` codex 2 回目採用 → `1c0296af` 接続模型 + `check_dual_closure.py` (別セッションの case/51 を巻き込み) → `ded72e17` 旧メッシャ自由端修正 (R5r)
→ `3c73ace8` 16 µm 化 + 品質ツール修正 → `995a519a`/`c21fdd1x` ドライバと IC → `a10a5755`…`2153b6d3` 診断記録 → `13ea76ca` レビュー採用 + IC 連続化
→ `5964312c` 引き継ぎ → `41508544` 面別収支ツール + run_0430/0431 の行。

(付録の使い捨てスクリプトは `cad/diag_wall_cv_budget.py --summary` に置き換わったので削除した。)
