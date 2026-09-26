# 諮問ブリーフ: V6′ の result 段の解釈 — 帯内 G-if を測った後 (2026-09-26)

AGENTS.md エスカレーション条件 **7** (result 段の解釈を確定する前)。前回の諮問
[`notes/reviews/2026-09-26-v6p-verdict-brief-diagnose.md`](../2026-09-26-v6p-verdict-brief-diagnose.md)
(結論「判定不能」・次の一手 ①〜④) を **plan §5.1 #100/#101 で全件採用**し、その一手を実施した結果を諮る。

## 0. 読んでよいファイル (巨大ファイル禁止)

- plan: `plans/active/boundary-conjugate-heat-transfer.md` は**丸ごと読まない** (1000 行)。§6 V6′ の行 (`grep -n "^| \*\*V6′" `)、
  直後の「V6′ の登録の修正」ブロック、§5.1 #100 / #101 だけ。
- 評価器: `case/58.conjugate_slot/eval_v6p.py` (今回 (e) を差し替え)、`case/58.conjugate_slot/sens_v6p.py` (新規)
- 判定ツール: `solver_density_cuda/tools/check_cht_interface.py` (今回 `--band-y` を追加)、`solver_density_cuda/tools/check_quasisteady.py`
- ソルバ: `solver_density_cuda/conjugateWall.cpp` (`appendInterfaceNodeLog` と `updateFem2dWall`)。差分は `git show c1265aeb`
- run (いずれも `case/58.conjugate_slot/` 配下): `run_0012_v6p_{df5_i50,df20_i50,df5_i200,df20_i200}_nlog/` の
  `solverConfig.yaml` / `conjugate_history.csv` / `conjugate_gate.json` / `CHT_INTERFACE_VERDICT.txt` /
  `CHT_INTERFACE_BAND_VERDICT.txt` / `QUASISTEADY_band.txt` / `v6p_band.json` / `CONVERGENCE_VERDICT.txt` /
  `conjugate_iface_nodes_5.csv`。**`conjugate_iface_log_5.csv` (52 MB) と `v6p_band_series.csv` は numpy で必要な列だけ読むこと**。
  `res_100000.h5` は 32641 節点 (必要な配列だけ)。
- case README: `case/58.conjugate_slot/README.md`

## 1. 事前登録 (前回から変更なし)

問題・漸近解・帯・合否 (a)〜(f)・「合否は全域 FP64 ビルドだけで取る」・G-if 閾値 (`eps_abs_Wm2` 23.5 / `eps_rel` 1e-3 /
`dT_K` 1e-2 / `tol_solid` 1e-2 / `n_consec` 80、範囲は帯)・準定常 (帯平均 $q_{w1}$・$T_{w1}$ が STEADY、流体残差は
NOT CONVERGED になりうるので**界面ゲート + 準定常で判定**、`--drift`/`--osc` は比較許容の 1/5 以下)・
感度 (2×2 `Df_scale` 5/20 × `interval` 50/200、帯内 $T_{w1}$ 差 ≤0.5 % of 固体上昇・$q$ 差 ≤0.5 %)。
**合否を取る run は `Df_scale` 5 / `interval` 50** (旧 `run_0010`)。

**前回の諮問後に変えたもの (採用した指摘の実装。閾値・帯の決め方は変えていない)**:

1. (e) を $Q_{\rm sol}=(K_su-b_s)_{\rm iface}$ (run の `solid.h5` から `Fem2DOperator.assemble_full`) と壁ダンプ
   `iface_Qf_eff` の符号つき差 / 集中辺長に差し替え (`eval_v6p.py`)。旧 run で codex の数値を完全再現
   (0.0025 / 0.0391 / 0.0369 / 0.4607 %)。
2. ソルバに `conjugate.node_log: 1` (界面全節点の更新前物理残差 $r_i$・$Q_{f,i}$・$\Delta T_i$・$T_{w,i}$ を毎更新 CSV に。
   host 側の出力のみ・更新式は不変)。`check_cht_interface.py --band-y YTOP YBOT` が帯内節点で ①②③、④は全域で判定。
   帯は `eval_v6p.py` が決めて `v6p_band.json` に書いた値を**そのまま**渡した。
3. 準定常は `eval_v6p.py` が節点ログから作る `v6p_band_series.csv` (毎更新、帯平均と**節点ごと**の $T_w-300$・$Q_f/A$) を
   `check_quasisteady.py --series-csv ... --tail 0.5 --drift 0.001 --osc 0.001` (前回 codex が使った値) で判定。
4. 感度は `sens_v6p.py` (共通帯 = 各 run の `v6p_band.json` の交わり、節点ごと、最終スナップショット)。

## 2. 再現条件

4 run とも AWS g5 (A10G)・**全域 FP64 ビルド** (`c1265aeb`、`flowFormat.hpp` を double に置換したツリー、
`FORGE_CUDA_BLOCKSIZE=128` — 既定 512 では `convectiveFlux_d.cu`:321 が `too many resources requested for launch`)。
入力 (`mesh.h5` = `run_0009` の最終場・`solid.h5`・config) は `run_0010` / `run_0011_*` と**バイト同一**で、
**差は `node_log: 1` の 1 行だけ**。100k step、更新は 1960 回 (interval 50) / 490 回 (interval 200)。

## 3. 観測事実

**NaN/Inf**: 4 run とも残差履歴・全スナップショットに無し。最終場 $\rho$ 0.0697 以上、$P$ 9830–13486 Pa、$T$ 294–500 K。
**流体残差** (`check_convergence.py`): 4 run とも `NOT CONVERGED (stalled/plateau)` (低下 0.0–0.3 桁)。

**再現性**: `run_0012_v6p_df5_i50_nlog` vs 旧 `run_0010_v6p_fp64` の共通帯で温度差 max **0.0001 %**、$q$ 差 **0.0021 %**。

### 3.1 評価器 (a)〜(f) (各 run 自身の帯、帯内 max)

| run | 帯 | (a) | (b) | (c) | (d) | **(e)** | (f) | VERDICT |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **df5_i50** (合否 run) | 211 行 −4.166…−19.779 mm 15.6 W | 0.1529 | 0.0080 | 0.2072 | 0.1992 | **0.0011** | 0.0000 K | **PASS** |
| df20_i50 | 212 行 15.7 W | 0.1408 | 0.0089 | 0.1966 | 0.2043 | 0.0386 | 0.0000 | PASS |
| df5_i200 | 212 行 15.7 W | 0.1428 | 0.0081 | 0.2008 | 0.2086 | 0.0370 | 0.0000 | PASS |
| df20_i200 | 230 行 −3.297…−19.779 mm 16.5 W | 0.1045 | 0.0991 | 0.1165 | 0.0233 | **0.5050** | 0.0000 | **FAIL (e)** |

### 3.2 帯内 G-if (`CHT_INTERFACE_BAND_VERDICT.txt`、末尾 80 更新) と全域 G-if (報告項目)

| run | 帯内 ① [W/m²] | ② | ③ [K] | ④ (全域) | 帯内 VERDICT | 全域 ① / ② | 全域 VERDICT |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **df5_i50** | **0.240** | **2.26e-5** | **2.66e-5** | 1.4e-9 | **PASS** | 32.8 / 2.32e-3 | NOT CONVERGED |
| df20_i50 | 2.14 | 1.95e-4 | 6.0e-5 | 1.4e-9 | PASS | 32.2 / 2.28e-3 | NOT CONVERGED |
| df5_i200 | 3.06 | 2.80e-4 | 3.4e-4 | 1.4e-9 | PASS | 33.2 / 2.35e-3 | NOT CONVERGED |
| df20_i200 | **23.88** | **2.03e-3** | 6.7e-4 | 1.4e-9 | **NOT CONVERGED** | 779 / 5.46e-2 | NOT CONVERGED |

### 3.3 準定常 (`QUASISTEADY_band.txt`)

| run | 帯平均 $T_w-300$ | 帯平均 $q$ | 局所 非 STEADY ($T$ / $q$) | OVERALL |
| --- | --- | --- | --- | --- |
| **df5_i50** | STEADY | STEADY | **0 / 0** (211 節点) | **ALL STEADY** |
| df20_i50 | STEADY (単調増加、漸近 94.187、最終比 +0.001 %) | STEADY | 15 / 0 | NOT ALL STEADY |
| df5_i200 | STEADY (単調増加、+0.003 %) | STEADY | 15 / 0 | NOT ALL STEADY |
| df20_i200 | STEADY (単調増加、漸近 94.1986、+0.074 %) | STEADY | 8 / 7 | NOT ALL STEADY |

**非 STEADY の節点はすべて帯の上端 (lip 側)**: df20_i50 と df5_i200 は $y$ −4.119 … −4.800 mm の連続 15 節点 (末尾半分で
+0.049 … +0.120 K、単調増加)、df20_i200 は −3.297 … −3.610 mm (温度 +0.05 … +0.09 K、$q$ +5 … +11 W/m²)。
深部は全 run で STEADY。

### 3.4 上端節点の時系列 (私が手で出した。ツール経由ではない — 節点ログを直接読んだ)

基準 = df5_i50 の最終壁温。差 $T_w-T_{\rm base}$ [K]:

| run | 節点 ($y$) | step 25k | 50k | 75k | 100k |
| --- | --- | --- | --- | --- | --- |
| df20_i50 | −4.119 mm | −0.2146 | −0.1423 | −0.0630 | −0.0225 |
| df20_i50 | −3.297 mm | −0.6935 | −0.3537 | −0.1193 | −0.0352 |
| df5_i200 | −4.119 mm | −0.2103 | −0.1364 | −0.0588 | −0.0206 |
| df5_i200 | −3.297 mm | −0.6796 | −0.3361 | −0.1101 | −0.0319 |
| df20_i200 | −4.119 mm | −0.2861 | −0.2569 | −0.2372 | −0.2184 |
| df20_i200 | −3.297 mm | −0.7332 | −0.7470 | −0.7247 | −0.6590 |

### 3.5 感度 (`sens_v6p.py`、共通帯 211 節点 −4.166…−19.779 mm、基準 df5_i50、最終スナップショット)

| 比較先 | 温度差 max / 平均 [% of 94.18 K] | $q$ 差 max / 平均 [% of $q_*$] | 登録許容 |
| --- | --- | --- | --- |
| df20_i50 | 0.0231 / 0.0037 | 0.0163 / 0.0024 | 0.5 % |
| df5_i200 | 0.0212 / 0.0034 | 0.0161 / 0.0018 | 0.5 % |
| df20_i200 | 0.2186 / 0.0607 | 0.0530 / 0.0432 | 0.5 % |

## 4. 仮説 (私の読み。確定していない)

**H1**: 合否 run (5/50) は登録した全条件 ((a)〜(f)・帯内 G-if・帯平均と局所の準定常) を満たすので **V6′ は PASS**。
流体残差 NOT CONVERGED と全域 G-if NOT CONVERGED は、登録で「流体残差は NOT CONVERGED になりうる」「全域 G-if は報告項目」と
した範囲に入る。

**H2**: 他の 3 run の局所非 STEADY・df20_i200 の (e) FAIL と帯内 G-if 不合格は、**連成の緩和が遅いことによる未収束**であって
**固定点の $D_f$/interval 依存ではない**。根拠: (i) 更新式は残差補正形で固定点は $D_f$ に依らない (`conjugateWall.cpp` の
コメント「$r$ が $D_f$ に依らない」)、(ii) 実効緩和の遅さは $D_f\times$interval に比例するはずで、4 倍遅い 20/50 と 5/200 が
**ほぼ同じ挙動** (同じ 15 節点・差 0.0231 vs 0.0212 %・§3.4 で同じ減衰列)、16 倍遅い 20/200 はずっと遅い、(iii) §3.4 で
20/50・5/200 の基準との差は 25k step ごとに約 0.4 倍に**単調に縮んでいて基準へ向かう**。
ただし (ii)(iii) は**私が手で読んだもの**で、20/200 を伸ばして基準に達するかは**測っていない**。

## 5. 諮りたいこと (それぞれ 1 つに絞ってほしい)

1. **V6′ の総合を PASS と書いてよいか** (PASS / 条件つき / FAIL / 判定不能)。特に、**感度 2×2 の比較 run が局所で静定していない**
   ことと **df20_i200 が (e) と帯内 G-if で不合格**であることが、感度の登録条件 (差 ≤0.5 %) を満たしたと書く妨げになるか。
   前回あなたは「局所静定と帯内 G-if が揃わないので収束解の感度が合格したとは判定しない」とした。今回、**合否 run は揃った**が
   比較 run は揃っていない。
2. **H2 の判別 A/B を 1 つ**。私の案: df20_i200 を `conjugate_state_5.h5` から再開して更新回数を 5/50 と揃える (追加 1470 更新 =
   294k step) か、または df20_i50 / df5_i200 を 100k step 延長して局所 ALL STEADY に達したときの基準との差を測る。
   事前に書く合格条件の案: 「延長後に局所 ALL STEADY・帯内 G-if PASS・(e) ≤0.1 %、かつ基準との差が温度 ≤0.5 % / $q$ ≤0.5 %」。
   **これが必要か**、必要なら**どちらを・どの長さで**。
3. **(e) の df20_i200 0.5050 % の扱い**。私は「未収束の連成の不釣合いがそのまま出ている」と読んでいる (帯内 G-if ① 23.9 W/m² は
   $q_*$ の 0.51 % で (e) と同じ大きさ)。**(e) は合否 run だけに課す条件か、感度の各 run にも課す条件か**が登録で曖昧。
4. **見落とし**。特に (i) 準定常の `--tail 0.5` は前回あなたが使った値で、登録には無い (登録は drift/osc の 1/5 だけ)、
   (ii) 合否 run の帯内 G-if は ①が許容の 1 %・③が 0.3 % と余裕が大きいが、全域 ① 32.8 > 23.5 は lip 行が支配 —
   帯の上端が lip の影響を受けていないと言える根拠は帯の決め方 (Pe・線形残差・勾配比) だけ、
   (iii) 局所準定常を $T_w-300$ (≈94 K) で割った drift 0.1 % = 0.094 K は、(a) の許容 0.471 K の 1/5 に一致する。

**両論併記で逃げず、推奨を 1 つに絞ること。根拠は `ファイル:行` か上の数値で示すこと。**
