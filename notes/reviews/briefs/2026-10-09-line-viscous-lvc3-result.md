# 諮問: 全行のスカラー対角を戻した値 3 も発散した。次の一手 (切り分けか保留か)

日付 2026-10-09。諮問先 codex (diagnose)。条件 4・7。plan: `plans/active/time_integration-line-viscous-jacobian.md` (§6.0〜§6.5 を全文)。前回の諮問: `notes/reviews/2026-10-09-line-viscous-thermal-weak-mode-diagnose.md`。

## 観測事実 (plan §6.5)

- 同じバイナリ (6d49738e + typedef double) で、run_0183 の res_100000 から キー 5・方向別・上限なし・cfl 4・緩和 0.7・sweep 5:
  値 2 (`case/45.isobutane_m6_d155/run_0300_vn1b_lvc2`) は 29 step で非有限 (step 1: ρ 23 倍・ρv 58 倍、step 2: ρ 487 倍・ρv 1270 倍)。
  値 3 = 値 2 + ライン面のスカラー 2ν_eff δ/dcc を全行の対角 (`run_0301_vn1b_lvc3`) は 122 step で非有限 (ρv が先: step 1 で 2.6 倍、5 で 11 倍、20 で 34 倍、50 で 1080 倍、ρ は step 50 で 489 倍)。
  値 3 の 1 step 目の書き出し (`run_0302_dump_lvc3`) は host の再解と CUDA の差 ≤ 5e-13。
- 参考 (別のバイナリ 5ab83056、ただし値 0 の経路はこの間変えていない): 値 0 (`run_0260_vn1_lvc0`) は同じ条件で 2000 step 生き残り、序盤の残差/開始は step 1〜100 で 0.5〜2 倍、最大 ρ 572 倍で回復。
- 値 2・3 が値 0 と違う点: (i) ライン面の薄層の D/K (運動量は βP∂u/∂Q、エネルギーは κ∂T/∂Q + 仕事、連続 0)、(ii) 強制の等温壁の行 [−e_w,0,0,0,1] (値 0 は単位行 = Δ(ρE)_w = 0)、
  (iii) 隣が等温壁・速度の Dirichlet のときの K の消去。値 3 は値 0 の対角をすべて含む。
- 改訂した解析 (拘束の消去・原始量) で値 2 の弱いモードは入口寄りで δT 81〜83 %・δρ/ρ 14〜16 %。

## 問い

1. 値 3 の不合格で、「薄層の D/K を足すこと自体が不安定を強める」と言ってよいか。等温壁の行 (ii) の変更が主因の可能性は。単因子で切り分けるなら、
   値 0 + キー 7 (= キー 5 + 等温壁の拘束の行) の 2000 step を同じバイナリで回すのが適切か (値 0 と同じ結果なら (ii) は外れ、(i)(iii) が残る)。
2. ρv (半径方向の運動量) が先に育つことは、薄層の P = I + ⅓ n̂n̂ᵀ (法線成分 4/3) や軸対称のフープの項との不整合を示唆するか。P を I にする単因子の試験の意味は。
3. ここで値 2 の系統を保留 (記録して閉じる) し、本線は「方向別 + キー 5 + 上限 50 で加速し point で仕上げる」の評価に戻るべきか。続けるなら、どの順で何を回すか (それぞれの事前の判定を含めて)。

## 読んでよいもの

上記の plan と前回の諮問、`solver_density_cuda/cuda_forge/timeIntegration_d.cu` (値 2/3 の分岐・等温壁の行)、`block_dplur_jacobian_d.cuh`、`case/45.isobutane_m6_d155/README.md`、`plans/active/time_integration-implicit-thermal-jacobian.md` §6.0 (run_0221〜0263)。
