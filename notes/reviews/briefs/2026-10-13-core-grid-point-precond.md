# 諮問ブリーフ: 低マッハの完全前処理 (`lowMachPrecond` 2/3) が Gref で効くかの安い試験 — point 陰解法での組み方

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.80・§4.81、背景は §4.26〜§4.31・§4.79)
- **date**: 2026-10-13
- **諮問の条件**: AGENTS.md「エスカレーション条件」1 (§4 の検証の新設)・4 (登録に無い run の前)
- **ユーザの依頼と決定**: 元の問い「なぜ流れ方向に細かいとライン陰解法の反復が遅くなるのか、つきとめてくれ」。ユーザの見立て「低マッハ補正が欠けたピースではないか」。§4.81 で `lowMachPrecond: 1` (右辺だけ) は B0 で発散した。そこで私が「(a) 安い確認: point 陰解法で `lowMachPrecond: 2` を回し、同じ起点の point の既存 run と比べ、完全な前処理がこのケースで効くかを先に確かめる / (b) ライン経路に Γ を実装する」を示し、ユーザは **(a)「安いけいさんでいいよ」** を選んだ。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.f64`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
grep -n '^### 4\.2[6-9]\|^### 4\.3[01]\|^### 4\.8[01]' plans/active/tooling-nozzle-core-grid.md   (§4.26〜§4.31・§4.80・§4.81 を sed で 30 行ずつ)
sed -n '225,245p' methods/time_integration/theory.md                 (前処理の 2 面)
sed -n '250,290p' solver_density_cuda/cuda_forge/setDT_d.cu          (前処理の Δτ′)
sed -n '1276,1300p;1345,1360p;1380,1395p;1436,1470p' solver_density_cuda/cuda_forge/timeIntegration_d.cu   (前処理の block カーネル、node の行の扱い)
sed -n '1383,1394p' solver_density_cuda/input/solverConfig.cpp        (lowMachPrecond ≥ 2 の起動条件)
sed -n '15,30p;100,125p;180,195p' notes/investigations/backstep-lowmach-checkerboard.md
```

## 1. 前提 (観測事実)

- case/45 の Gref (4719 × 160 = 75.5 万節点、node・軸対称・SLAU・SST・冷却壁 300 K、組成は一様 Y0 = 0.9142 (ばらつき 3e−13))。停滞帯は縮流部 x/r_t −6〜−3.25・壁から 60〜100 節点目、M 0.02〜0.08。
- B0 (ライン陰解法・方向別 dt 上限 50・ITJ 5) は Gref で停滞する (平均流の残差がプラトー、§4.22・§4.35)。
- point 陰解法 (lineImplicit 0・lineDtDirectional 0・ITJ 0、他は B0 と同じ cfl_pseudo 4・implicitRelax 0.7・nStepInner 5・`lowMachPrecond: 0`) を同じ起点 (`run_0608` res_60000) から回した `run_0611_cg_gref_point` (6 万 step、AWS A、FP64 バイナリ sha 7e95e5bb): rms_ro は 1.8 桁下がったが rms_roOmega は 29 倍に上昇、θ_r は −1.1 % まで動き続けた (§4.27)。続き `run_0613` (+6 万)・`run_0614` (+8 万) で通算 20 万 step に自己変化の資格を満たした (§4.29・§4.31)。つまり point は停滞しないが遅い。
- §4.81: B0 に `lowMachPrecond: 1` (右辺の c′ 散逸だけ) を足すと、停滞帯から step 10 で発散 (step 4 以降 1 step あたり約 3.3 倍)。左辺は FVS の c のまま。backstep (case/18) の `=1`・cfl 1 の発散と同じ型。backstep を直したのは `=2` (右辺の c′ + 左辺の Γ の時間項 + Δτ′、cfl_pseudo 1)。`=3` (左辺だけ、解は不変) も backstep の残差のプラトーを `=2` と同じ程度に破った。
- `lowMachPrecond ≥ 2` はライン陰解法と併用できない (`main.cpp:3446`)、ITJ ≠ 0 とも併用できない。point (blockDPLUR 1・timeIntegration 11) では使える。
- node の前処理カーネルは境界節点を凍結させる欠陥があった (境界半割面の粘性対角・軸と壁の拘束の行)。修正 (別ブランチの e224ae11) をこのブランチに取り込んだ (commit a60d8eb5、`lowMachPrecond` < 2 の経路は不変)。それを入れた FP64 バイナリを AWS A でビルド中 (`~/forge-fgeom7-fp64-lmp2`、元は d0474498 の diag5 と同じソース + この差分)。`run_0611` のバイナリ (7e95e5bb) とはソースが違う (診断のコードは既定 off)。
- 同じ入力から回しても GPU の加算順で実行ごとに差が出る (§4.77)。ビット比較はできない。
- 予算: point の 2 万 step は 1 本で約 10 分、4 本同時で 40 分程度 (AWS A、メモリ 1 本 約 2.5 GB、空き 15 GB)。

## 2. 私の案 (縛られなくてよい)

- 腕 (すべて `run_0608` res_60000 から restart_field、新しいバイナリ、point、各 2 万 step・2500 ごと): C = `lowMachPrecond: 0` (対照)、P3 = `=3` (左辺だけ、解は不変)、P2 = `=2` (左右とも)。cfl_pseudo は全腕 4 (前処理の腕は Δτ′ が setDT で拡大される)。P2 を cfl_pseudo 1 にした腕を足すか迷っている。
- 主指標: 縮流部 x/r_t ∈ [−8, −2] の体積重み付き残差 q (ρ・ρE) の比 P/C、窓 10k〜15k・15k〜20k (§4.40 以来と同じ定義)。全保存量の残差 (とくに ρω) と縮流部の外。
- 読み方の案: P3 または P2 が両窓で C の 1/3 以下 → 低マッハの時間項の前処理がこのケースの反復を速める → ライン経路への実装を検討する価値あり。両方とも 0.8 以上 → 低マッハの前処理を入れても point の遅さは変わらない。

## 3. 判断してほしいこと

推奨を 1 つに絞り、事前の読み方つきで:
1. 腕の組み (どの `lowMachPrecond`、cfl、対照を新しく回すか `run_0611` を使うか、長さ)。
2. 主指標と結果 A/B の閾値、資格 (前処理が効いていることの確認の仕方を含む)。
3. この安い試験の結果を、B0 (ライン) の停滞と「ライン経路に Γ を実装するか」の判断にどこまで持ち込めるか (持ち込めない点も)。
