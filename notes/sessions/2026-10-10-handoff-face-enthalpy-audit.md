# 引き継ぎ ③: point の仕上げにおける残差評価の精度の監査 (面エンタルピーの double) (2026-10-10)

共通ルールは [`2026-10-10-handoff-common.md`](2026-10-10-handoff-common.md)。正本は `plans/active/time_integration-implicit-thermal-jacobian.md` の **§5.1 #5** (この plan はこのセッションが担当してよい)。

## なぜ今やるか

元のセッションは、case/45 を float のビルドで回せるようにしている (`plans/active/architecture-float-state-double-geometry.md`)。上位の見立ては次の 2 つで、どちらも未確認。

- float では、残差の評価の丸めが、収束の床や θ_r の収束先に効くかもしれない。
- TP の面エンタルピー `thermo_h_mix` の float の評価は、同じ状態のエネルギー行の微分に強く効いた (粘性ヤコビアン plan §6.15: ‖J_t p‖ が 97 % 変わった)。

**FP64 のビルドでも、面エンタルピーは float の係数で評価している** (`thermo_d.cuh:133, 1079` 付近)。この A/B は、float 化で「どこを double に残すべきか」の手がかりになる。

## A/B の設計 (上位の諮問で決まっている)

- 設計の出典: `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md` の「判別 A/B」。plan §5.1 #5 にも写しがある。
- 同じ point の仕上げの保存状態・同じバイナリ・同じ設定から、腕を 2 つ回す。
  - A: 今の float の面エンタルピー
  - B: `FORGE_DIAG_FACE_H_DOUBLE=1` だけを変える
- 各 2000 step。全残差を毎 step 記録し、開始・1000・1500・2000 step の状態を保存する。
- 末尾 500 step の全残差の水準と傾きを見る。保存した状態を、共通の残差評価の精度でも評価し直す (表示の残差の定義が違うだけの改善を分ける)。
- 判定:
  - B だけで、共通の評価のエネルギーの残差が 10 % 以上下がり、再評価のノイズを十分上回る → この期間の停滞への寄与を支持する。
  - 自分の方式の表示の残差だけが変わり、共通の評価では 10 % 未満 → 主要因説を支持しない。
  - まだ減衰中なら、残差の床の判定は不能とする。自動で延長はしない。

## 段取り

- 出発の状態: point の仕上げの run の最終の場。例 `case/45.isobutane_m6_d155/run_0354_m9_L5cut` (res_40000)、または `run_0263_ns_coldmesh_tw300_cfl4_ext4` (res_40000)。AWS 上にある。
- 準備: `cold_cfl.py prep <src> <run> --field-from <上の run>` (point の構成。`--line` を付けない)。
- バイナリは `lineM_fp64` (共通ルールのとおり)。
- run 番号は case/45 の 05xx。
- 「共通の残差評価の精度で評価し直す」方法: 保存した状態を、A の設定の 1 step だけの run (更新前の残差を記録) に通すのが簡単。
  - 元のセッションの `omg.sh` (状態から 1 step で残差の場を出す) が参考になる。
- 事前登録 (§6 に判定と run の表を書く) と、codex の plan 段を、run の前に通す。
