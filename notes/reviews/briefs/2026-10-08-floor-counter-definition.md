# 諮問: 毎更新の床補正のカウンタの定義のやり直し (定常陰解法では床が保存量に commit されない) (2026-10-08)

関連 plan: `plans/active/tooling-sern-te-wake-grid.md` §4「床の判定」・§5.1 #2・§6 #1、起点 `plans/active/convection-zero-thickness-edge-reconstruction.md` §5 の 3。codex plan レビュー `notes/reviews/2026-10-08-tooling-sern-te-wake-grid-plan.md` M1・M2。
エスカレーション条件 1 (§4・§6 の定義の変更)。

## 観測事実 (実装担当のコード読み、編集なし。SERN 生産設定 = node・TP・unsteady 0・timeIntegration 11・blockDPLUR 1・nStepInner 5・pMin 20、例 `case/46.sern_design/run_0446_3d_gs050_flag1/solverConfig.yaml`)
`advanceImplicitSteady` (`main.cpp:2636`) の 1 外側 step:
1. EOS は step の冒頭に 1 回 (`main.cpp:2032`、`assembleResidualPre` の中) で保存量の配列をその場で書き換える: `dependentVariables_d.cu:76` `ro_temp = max(ro, roMin)`、`:207` `ro[ic] = ro_temp`。TP は `:172` の温度反転の各反復で `DEPVAR_TMIN` (50 K) にクランプし、その T から `:209` で `roe = ρ·(e_mix(T) + ek)` を作り直す (T が 50 K に当たると作業状態の roe が e(ρ,Y,50 K) 相当まで持ち上がる)。**TP の圧力床 (`:203`) は P だけを書き換え roe には触れない**。CPG 単相は tMin (`:290`) が T だけ、pMin (`:291`/`:297`) は roe も作り直す。
2. 床後の作業状態から流束・`setDT`・block-DPLUR の係数を組む。
3. commit は `ro = roN + d0`, `roe = roeN + d4` (`update_d.cu:301-308`、`main.cpp:2491`)。`roN` は前の step 末の `updateVariablesOuter` (`main.cpp:2647`) で取った commit 直後 = **床前**の値 → **床による補正は捨てられる** (FP64 アキュムレータ `update_d.cu:285-299` も同様)。
4. 等温壁ピン (`main.cpp:2644`) → `updateVariablesOuter` で `roN = ro` (床前) → 出力 (`main.cpp:2661`)。`res_*.h5` の ro/roe は commit 直後の床前、T/P はその step 冒頭の EOS の値 (床後) で、1 更新ずれる。既存の `floor_gate()` (`sern_gates.py:200`) の `T <= 50` は最後から 1 つ前の更新結果を EOS にかけた値。
- 陽解法 RK も基準 `roN`/`roM` を EOS の前に取る (`timeIntegration_d.cu:465`、`update_d.cu:126`) ので捨てられる。**dual-time だけは `ro += dq` をその場で加算 (`update_d.cu:517-521`) し、床後の作業状態が基準になる (commit される)**。メモリの記録 (eos-floor-not-committed-steady) と一致。
- 帰結: plan の定義「更新後の保存量の内部エネルギーが e(ρ, Y, T_min) を下回り、最終状態を床へ補正した事象」を定義どおり数えると、定常陰解法・陽解法では構造上 0 件 → ゲートが中身のない合格になる (A′ の 517160 のように床に張り付いた run も通る)。

## 代わりの候補 (実装担当の案)
- **A. EOS の入力時点での床該当** (`dependentVariables_d` に出力専用の atomic 計数): 実節点 (`ic < nCells`) と ghost を分ける。温度床 = 反転の最終値が `DEPVAR_TMIN` かつ `roe 再構成 − roe 入力` が丸め上限を超える (補正量 Δ(ρE)、反復中のクランプは数えない)。圧力床 = P の生値 < pMin (補正量 ΔP、TP では Δ(ρE) は定義できない)。密度床 = ρ < roMin (Δρ)。
  意味は「直前の更新で commit された状態 (壁ピン後) が床の領域にあり、その step の流束と係数を床後の作業状態で組んだ」。床が commit されないので張り付いた節点は毎 step 数えられ続ける (既存の floor_gate を毎 step に広げたもの)。穴: 最後の step の commit 結果は EOS にかからない (終了時の読み取り判定が要る)、開始・restart 直後の最初の EOS は初期場を見る (印で区別)。費用ほぼ 0。
- **B. commit 直後の読み取り専用カーネル** (`updateVariablesOuter` の後・出力の前): commit 済みの保存量を e(ρ,Y,T_min)・pMin・roMin と比べる。定義の前半をそのまま評価でき最後の更新も含む。費用は 1 節点・1 step あたり NASA-9 の評価 1 回分。化学種の commit は流れの commit より後なので Y の読み位置に注意。壁の滑りなしピンは step の冒頭にしか当たらないので壁節点では A と見る状態が少し違う。
- **C. 定義どおり (commit された床補正)**: 定常陰解法・陽解法では構造上 0 件、意味を持つのは dual-time だけ。ゲートに使えない。
- 実装担当の見立て: A の計数 + 終了時の B 相当 1 回が最も安価で趣旨 (判定区間のどの更新でも床に触れていない) に合う。

## 問い
1. 新しい定義をどれにするか (A / B / A + 終了時の 1 回 / 別案)。§6 #1 の合格条件 (判定区間の全実節点・全更新で 0 件) の読み方。
2. 決めるべき細部: (i) TP の圧力床は保存量を補正しないので補正量を Δ(ρE) にするか ΔP にするか、(ii) 判定区間の端 (最初の EOS は初期場・restart の状態、最後の commit は EOS にかからない)、(iii) dual-time では A と C が一致する — 生産以外の経路で数え方を分けるか、(iv) 「丸め上限」の具体値、(v) 床に**近い**状態 (例 T ≤ T_min + 1 K) を別に数えるか。
3. 定常陰解法で床が commit されないこと自体 (作業状態と保存量の乖離) は、この plan の受入れで問題として扱うべきか (別 plan か)。
