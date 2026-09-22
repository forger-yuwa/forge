# 引き継ぎ: ⑤ SERN 全ヘキサ接続模型の CFD 起動 (2026-09-22, session_014aoVV4K167qqfTw7TojyJc → 次セッション)

**正本は plan** [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md) (§4.10–§4.13、§5.1 の B1b/B1d、§6.1)。
本文書は写しとポインタ。残作業を増やすときは plan §5.1 に書く (ここには書かない)。

## 1. いま何がどこまでか (1 分で)

- **格子は完成** (plan §4.11–§4.12): 全ヘキサ接続模型 `case/46.sern_design/cad/hex_junction_model.py` (19 ブロック/断面 × 3 区間、有限厚の側壁・カウル板、
  鈍頭の後端面)。採用設定 `--set H1=1.6e-4 H1_END=2.5e-3 HX=0.05 HX_FAR=0.08 NZ=31` → 2,446,572 節点、`check_mesh_quality.py` PASS (AR max 676, 緩和なし)。
  出力 `case/46.sern_design/cad/demo_out/jm_c16.msh` (330 MB, git 追跡外。AWS の同じパスにもある)。
- **CFD は起動できていない** (plan §4.13): run_0423–0429 の 7 run。IC の罠 2 つは潰した。残る発散は**側壁後端面の壁節点で ρ が 1〜1.5 %/step で抜ける**現象
  (200 step 刻みの診断で順序を確認: step 800 に後端面 × カウル上面 × 内側の縁の 3 重の隅から始まり、床を割って NaN)。
- **codex 3 回目レビュー (NO-GO M5/m1, 全件採用)**: 「Dirichlet 壁 CV の排出が真因」は**仮説に格下げ**、「テーパ推奨」は**撤回**。
  次は**鈍頭のまま** ① IC 修正 (済) → ② 問題 CV の面別 SLAU 流束を場から計算 → ③ 断熱壁 A/B。記録 `notes/reviews/2026-09-22-tooling-sern-mesh-blocking-plan-2.md`。
- **ユーザ決定 (済)**: 第一層は 16 µm・200〜300 万節点で検証できる規模を優先 (§4.12)。テーパにするかは未決 (codex の指摘で判断材料が足りない)。

## 2. 最初にやること

1. **AWS を起動してもらう** (ユーザに IP を聞く。stop/start で IP が変わる。`i-0b1a5e0b8dc152f00`, `ssh -i ~/.ssh/test.pem ubuntu@<IP>`)。
   接続したら **swap を張り直す** (`sudo swapon /swapfile; sudo swapon /swapfile2`) と `cd ~/forge-sern && git pull --ff-only`。
   **自動停止対策**: 作業中は `ssh -i ~/.ssh/test.pem ubuntu@<IP> 'sleep 7200' &` を手元で張っておく (30 分ログインなしで stop する。今日 2 回落ちた)。
   ディスクは 14 GB 空き。足りなければ run_0423–0429 の `sern.h5`/`sern.msh` (各 1.1 GB) を消してよい (診断は plan に記録済み、破棄予定)。
2. **run_0430 / run_0431 を投入** (下のコマンド。手元 `13ea76ca` の修正 IC が AWS に pull されていることを確認)。
3. 終了後 (各 10〜15 分) に §3 の解析。

### 投入コマンド (AWS 上)

```bash
cd ~/forge-sern/case/46.sern_design
python3 cad/run_junction_model.py problem_3d_prod_m6on_wallres.yaml cad/demo_out/jm_c16.msh run_0430_3d_junction_diag_icfix --prepare-only
cd run_0430_3d_junction_diag_icfix
# 暖機段の config を手で作る (run_staged が書くものと同じ: 層流・1 次・CFL 0.2) + 200 step 刻み出力
sed -E 's/nStepOuter: [0-9]+/nStepOuter: 1800/; s/outStepInterval: [0-9]+/outStepInterval: 200/; s/convMethod: 1/convMethod: 0/; s/cfl: [0-9.]+, cfl_pseudo: [0-9.]+/cfl: 0.2, cfl_pseudo: 0.2/' solverConfig_main.yaml \
 | python3 -c 'import sys,re; s=sys.stdin.read(); print(re.sub(r"turbulence: \{model: \"sst\"[^}]*\}", "turbulence: {model: \"none\"}", s), end="")' > solverConfig.yaml
grep -o 'turbulence: {[^}]*}\|convMethod: [0-9]\|cfl: [0-9.]*\|nStepOuter: [0-9]*\|outStepInterval: [0-9]*\|nodeWallDirichlet: [0-9]' solverConfig.yaml
# 断熱壁 A/B (運動量固定だけを残す。codex M2)
cd .. && mkdir run_0431_3d_junction_diag_adiabatic && cd run_0431_3d_junction_diag_adiabatic
for f in sern.h5 probe.yaml species_db.yaml species_meta.yaml solverConfig.yaml prepare_info.json MESH_QUALITY.txt; do cp ../run_0430_3d_junction_diag_icfix/$f .; done
sed -E 's/wall_isothermal, /wall, /; s/, Ts: [0-9.]+//' ../run_0430_3d_junction_diag_icfix/bcondConfig.yaml > bcondConfig.yaml
cd ~/forge-sern
nohup solver_density_cuda/tools/run_case.sh case/46.sern_design/run_0430_3d_junction_diag_icfix > /tmp/rc0430.log 2>&1 &
nohup solver_density_cuda/tools/run_case.sh case/46.sern_design/run_0431_3d_junction_diag_adiabatic > /tmp/rc0431.log 2>&1 &
```

`build/forge` を直接呼ぶと PreToolUse フックが止める (run_case.sh 経由のみ)。case README の run 一覧に 2 行足す (run_0429 の行の下)。

## 3. 解析 (B1d ②③)

- **ρ min の時系列**: 各 `res_*.h5` の `VALUE/ro,P` の最小値と位置 (`CELLS/centCoords` は双対重心 = 壁節点は面から Δx/4 だけ内側に写る。
  x = 1.2006 H の節点は後端面の壁節点そのもの)。今日のスクリプトは plan §4.13 の表を作ったもの (この文書末尾に再掲)。
- **問題 CV の面別 SLAU 流束** (codex M1): 前回の最悪節点は (x/H 1.2006, y/H 0.009, z/H 1.0025–1.0038) と (1.2006, 1.293–1.316, 1.0046)。
  `PLANES/STRUCT` (`check_dual_closure.py` の `parse_struct` で owner/neighbor)、`PLANES/surfVect`、`CELLS/volume` と各 dump の原始量から、
  1 次 (L/R = 節点値) の SLAU 質量流束を面ごとに再計算: 移流項 $\rho_w\rho_i/(\rho_w+\rho_i)\,u_n$ 相当と圧力差項 $-\chi(P_R-P_L)/\hat c$、$\chi$ の値、両側音速。
  SLAU の式は `solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:533` 付近。**χ が 0 になっているか**が問い。
  併せて壁節点の ρ の dump 間の減少量と、面流束の和 × 双対体積の比較 (床補正が入っていれば差が出る。codex M4)。
- **断熱 A/B (run_0431)**: 同じ節点群の ρ 履歴を等温 (run_0430) と並べる。差が無ければ熱的固定は無関係、運動量固定側の仮説が残る。
- **合否は「完走」で言わない** (codex M5): 全壁タグの床到達数 (ρ ≤ roMin 1e-4 / P ≤ pMin 20)、固定 CV 群の密度履歴、逆流域。
- 結果は plan §4.13 に「訂正」の後ろへ追記し、B1d の内容欄を更新する。**真因を書く前に**: 主セッションが Fable でなければ `diagnostician` に諮る (AGENTS.md の条件 4)。

## 4. 判断が要ること (ユーザ or Fable)

- 鈍頭で起動できる手立てが見つからなければ、**テーパ** (t → t/5 を 5t で絞る) か **後端面の扱いのソルバ側変更** (後端面だけ弱い壁、など) の選択。
  codex M5: テーパを試すなら鈍頭対照を残し t_base と Δx₁ を独立に振る。
- 変換器の単精度 (CV 閉性 max 1.7e-4、B1b) をいつ直すか (別セッションの変換器作業と衝突しないよう着手前に確認)。

## 5. 既知の罠 (今日踏んだもの)

- AWS の自動停止 (§2)。AWS の `~/forge-sern` は shallow clone、`solver_density_cuda/build` は `.build-native/relwithdebinfo` へのリンク。gmsh は未導入 (格子は手元で切って scp)。
- `check_mesh_quality.py` の六面体の体積判定は今日 頂点 Jacobian に変えた (`3c73ace8`)。旧判定は曲面壁の薄いセルを誤って FATAL にしていた。
- 共有ワークツリー: 別セッション (case/51・53・54 の CHT/遷移) が同じブランチで作業中。**commit は必ず `git commit -- <paths>`** (index 経由だと相手の staged を巻き込む。今日 1 回やった)。
- `centCoords` は双対重心。壁節点の位置判定は `MESH/COORD` を使う。
- run_junction_model.py の IC は接続模型の寸法 (H, ZW, TSW, TC, LSW, LCOWL) を定数で持つ。形状を変えたら合わせる (codex M3 の注意)。

## 6. 今日の commit (feature/sern-design)

`645f25a4` B0 設計 → `e580114c` codex 2 回目採用 → `1c0296af` 接続模型 + `check_dual_closure.py` (別セッションの case/51 を巻き込み) → `ded72e17` 旧メッシャ自由端修正 (R5r)
→ `3c73ace8` 16 µm 化 + 品質ツール修正 → `995a519a`/`c21fdd1x` ドライバと IC → `a10a5755`…`2153b6d3` 診断記録 → `13ea76ca` レビュー採用 + IC 連続化。

## 付録: ρ min 時系列のスクリプト (AWS の run ディレクトリで)

```python
import h5py, numpy as np, glob, re
m=h5py.File("sern.h5","r"); cc=m["CELLS/centCoords"][...].reshape(-1,3)/0.1
fs=sorted(glob.glob("res_[0-9]*.h5"),key=lambda f:int(re.findall(r"\d+",f)[0]))+sorted(glob.glob("res_nan_*.h5"))
for fn in fs:
    V=h5py.File(fn,"r")["VALUE"]; ro=V["ro"][...]; P=V["P"][...]; ok=np.isfinite(ro)
    r2=np.where(ok,ro,np.nan); i=np.nanargmin(r2); p2=np.where(np.isfinite(P),P,np.nan); j=np.nanargmin(p2)
    print("%-16s ro_min %.3e at %s  P_min %8.3g at %s  ro<5e-4: %d  nonfinite %d"%(fn,r2[i],np.round(cc[i],4),p2[j],np.round(cc[j],4),(r2<5e-4).sum(),(~ok).sum()))
```
