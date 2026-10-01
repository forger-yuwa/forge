# 諮問 (diagnostician): 共役平板 C の方針決定 — リミッタ A/B の後どう進めるか

ユーザ指示 (2026-10-01): 「diagnostician に諮って方針決めて、計算進めて」(ユーザ就寝中。判断は本諮問で確定させ、主セッションが実行する)。
AGENTS.md 条件 1 (plan §4・§6 の方針変更)・3 (事前登録の判定が判定不能)・4 (plan に無い run の前)。

## 読んでよいファイル
- plan: `plans/active/boundary-cht-conjugate-flat-plate.md` (§3・§4.1〜§4.5・§5.1・§6)。発注元 `plans/accepted/boundary-cht-conjugate-benchmarks.md` (§4.3・§4.5〜§4.7・§6、C の登録条件と許容)。
- 判定出力: `case/65.conjugate_flat_plate/AB_JUDGE_run0011_0012.txt`、`AB_JUDGE_CFL_run0015_0016.txt`、`AB_JUDGE_LIM_run0017_0018.txt`、`AB_JUDGE_LIMEXT_run0017_0019.txt`、各 run の `CONVERGENCE_CHECK.txt`・`solverConfig.yaml`。
- 判定器: `case/65.conjugate_flat_plate/ab_judge*.py`、`ab_extract.py`。過去の codex 諮問 `notes/reviews/2026-10-01-conjugate-plate-*-diagnose.md`。
- `methods/limiter.md`、`procedures/solver-settings.md` (limiter / convMethod)、`procedures/recommended-settings.md`。
- ソルバ: `solver_density_cuda/cuda_forge/limiter_d.cu` 等。

## 観測事実 (AWS FP64 `86115cb1`、すべて run_0007_c1_n64 step 600000 起点、C1 n64)
1. B1 (壁温の動的更新 on/off、run_0011/0012): 停滞は同じ (残差比 0.996〜0.999)。
2. B2′ (cfl_pseudo 2 vs 0.5、壁温固定、run_0015/0016): 前縁・上流 slip の 7 step 周期の振動は 0.5 で消えた。残差は 0.43〜0.51 倍で横ばい、後縁直後の下流 slip の P 変動は残った (0.84)。登録判定は中間。
3. B3 (limiter 2 vs 0 at cfl 0.5、壁温固定、run_0017/0018、12000 step): limiter 2 では残差二乗和の 95〜99 % が後縁帯 (x/L 0.98〜1.1、y≤0.05L、節点の 8 %)、上位は x/L 1.004〜1.017 の y=0 と第 1 層。同じ場所で limiter_Uy の時間標準偏差 0.38〜0.41、limiter_P 0.23〜0.29。limiter 0 で後縁帯の割合は 1 % 未満になり下降開始。登録判定は中間 (rms_roUx が切替過渡で 13 倍→1.4 倍)。
4. B3 延長 (run_0019、limiter 0 累計 48000 step): B/A rms_ro 0.003、rms_roUy 0.0035、rms_roe 0.003、下流 slip P 0.004 (いずれも 1/10 を大きく下回り減少中)。**rms_roUx だけ 0.286→0.244→0.213 (4000 step ごとに約 13 % 減)**。連結 check_convergence: rms_ro/roUy/roe 3.7 桁低下で下降中、rms_roUx 2.7 桁 still converging → NOT CONVERGED。check_quasisteady (下流 P・界面 q) ALL STEADY。残差は 97〜99 % が前縁帯・後縁帯以外に分散。登録判定は判定不能。最終場 NaN/Inf なし。
5. 本番 run_0005〜0010 (C1/C2 × n16/32/64、limiter 2・cfl 2、600000 step): 主判定 (独立参照解との界面温度・熱流束) は 6 本 PASS、流体収束 NOT CONVERGED (stalled)、G-if NOT CONVERGED (n64 は前縁直後 678 W/m²、n16/n32 は後縁 ~6 W/m²)、④ res_solid は丸め床 (2^-32 の整数倍、#1c 未着手)。
6. C の完了条件 (§6、維持と決定済み): 6 本すべてが流体収束・準定常・G-if・G-cons・メッシュ品質・参照解自己検査・主判定に合格、格子差減少、C2 の固体効果。事後の窓除外・時間平均化・閾値緩和はしない。limiter 0 を C 専用の事後改訂設定とする余地はあるが未採用 (採用時は理由・適用範囲を書き、旧結果は判定不能のまま残す)。
7. 計算資源: AWS g5 (1 GPU、他セッションと共有)。n64 で ~20 ms/step 単独、2 本並走で ~45〜55 ms/step。本番 600000 step × 6 本は n64 で単独 ~3.3 h/本。

## 選択肢 (主セッション案)
1. limiter 0・cfl 0.5 を C 専用の事後改訂設定として採用し、本番 6 本 (連成あり、wall 連成 Df_scale 5、warmup 5000) を一様流 IC (または旧 run の最終場から) 600000 step で回し直して §6 で判定。G-if ④ の丸め床は #1c が未着手なので、その扱いも要決定。
2. 遅い rms_roUx 成分の原因を先に調べる (低マッハの遅い全体モード? 出口 BC?)。
3. C を「前提ゲート不合格、判定不能・未完了」で閉じる。

## 質問
- どれを採るか (組み合わせ可)。事前登録すべき判定条件・run 構成・IC/再開元・長さ・完了条件との関係を、中位モデルが迷わず実行できる粒度で示してほしい (触るファイル・合格条件の VERDICT・回す run)。
- ④ res_solid (丸め床) を C 再判定でどう扱うか。
- 一晩 (〜8 時間、GPU 1 枚共有) で回せる範囲で優先順を。
