# 諮問 (結果の解釈): R5h 判別 A/B の結果 — 分岐 A (再構成依存) と読んでよいか、次に何を測るか

前回の諮問: notes/reviews/2026-09-27-sern3d-r5h-te-coldspot-diagnose.md (同一保存場から convMethod だけ変える A/B を指定)。

## 実施 (2026-09-27、AWS、バイナリ 2fa3826c、g3 格子)
起点: `case/46.sern_design/run_0971_3d_g3_chidef_cont16k/res_16000.h5` を `restart_field.py` で 4 run に写した (9 量ビット一致)。
設定は run_0971 本段 (SLAU・limiter 2 Venkat limiterScaled 1 venkatK 0.05・SST・frozen_tp・等温壁 1000 K・blockDPLUR・CFL 0.25・chi 1・lsq) から
convMethod / nStepOuter / nStepInner / outStepInterval / extraFields (res_ro, res_roUx, res_roUy, res_roe, res_roK, limiter_*, dt_local, volume) のみ変更。
| run | convMethod | step | nStepInner |
| --- | --- | --- | --- |
| run_0974_r5h_ab_A1 | 1 (現行) | 1 | 1 |
| run_0975_r5h_ab_B1 | 0 | 1 | 1 |
| run_0976_r5h_ab_A100 | 1 | 100 | 5 |
| run_0977_r5h_ab_B100 | 0 | 100 | 5 |

## 観測事実
- 共通状態: A1 と B1 の res_0 は互いにビット一致。SRC との差は roe・roOmega が相対 1.9e-7 (restart 読込時の丸め。ro, roU, roY は一致)。
- 低温点 = 節点 517159 (0.12000, −0.010534, 0.097192)、出力 T 105.37 K。
- EOS 検算 (NASA-9 lump、h は 298.15 K 基準): 全域 3000 点サンプルで出力 T との差 max 0.10 K / median 1.3e-4 K。低温点は出力 105.373 K に対し保存量から 105.692 K (差 0.32 K、近傍 204 節点の最大もこの点)。T・P とも床 (tMin 50 K 相当, pMin 20 Pa) より十分上 (P 9367 Pa)。
- 注意: res_N の出力 T は最後の更新前の値 (1 step run で ΔT_出力 = 0 だが保存量は動いている)。以下の T は保存量から EOS で出し直した値。
- 1 step (nStepInner 1) の低温点:
  | | A (conv 1) | B (conv 0) |
  | --- | --- | --- |
  | res_ro | +4.96e-10 | +2.08e-06 |
  | res_roUx / res_roUy | +9.8e-09 / +2.7e-08 | +9.67e-03 / −1.12e-03 |
  | res_roe | −2.38e-04 | +6.36 |
  | res_roK | +1.28e-05 | +1.59e-02 |
  | Δro / Δroe (1 更新) | +2.1e-07 / −0.094 | −1.61e-03 / +2766 |
  | T (EOS) 更新後 | 105.692 (不変) | **116.281 (+10.6 K)** |
  | limiter_ro / P / T / Ux (この節点) | **0.000 / 0.269 / 0.000 / 0.281** | 1 / 1 / 0 / 1 (不使用) |
  | dt_local / volume | 1.08e-8 / 2.73e-11 | 同じ |
- 100 step (EOS 再計算の低温点 T / 近傍 204 節点の出力 T 最小 / 全域の出力 T 最小):
  A: 0〜100 step すべて 105.7 / 105.4 / 105.4 (定常のまま)
  B: 0: 105.7 → 10: 192.1 → 50: 361.6 → 100: **534.5** (単調増加・未飽和) / 全域最小は 10 step で 177.9 K (プルーム遠方) に移る
- 同じ状態・同じ勾配 (lsq) なので、A と B の残差差は対流流束の差だけ (粘性・ソースは共通) — と考えているが、ソルバで確認はしていない。
- 面別の流束分解は取っていない (出力にない。取るにはダンプ追加が要る)。

## 仮説 (呼び出し側)
- 分岐 A を支持: 共通状態で再構成差が低温点を加熱側へ動かす更新差 (res_roe +6.36 vs −2.4e-4) を与え、EOS 書き換え・床は関与しない。
- 低温点自身の ro・T の再構成は既に 1 次 (limiter 0) なので、差は**隣接節点側の面外挿**から来ている可能性。

## 問い
1. これで分岐 A (再構成依存、EOS 単独原因の除外) と記録してよいか。足りない証拠は何か (面別分解が必須か)。
2. 次の 1 手: (i) 面別の対流流束ダンプ (低温点の双対面ごとの質量・エネルギー流束を A/B で) を足す、(ii) 近傍節点のリミッタ・勾配を見る、(iii) 有限厚後縁で消えるか確かめる (blocking plan)、のどれを先にするか。あるいは別案。
3. 1 次の B で 534 K まで単調に上がり続けることは、2 次の定常解の低温点が「再構成が作った偽の定常」であることを意味するか、それとも 1 次の過渡に過ぎず結論できないか。

## 読んでよいファイル
- plans/active/tooling-nozzle-sern-3d.md §4.35・§4.46・§5.1 R5h
- notes/reviews/2026-09-27-sern3d-r5h-te-coldspot-diagnose.md
- solver_density_cuda/cuda_forge/ の対流流束・リミッタ・再構成 (node SLAU 経路) の該当関数のみ
巨大な h5・ログは読まない。
