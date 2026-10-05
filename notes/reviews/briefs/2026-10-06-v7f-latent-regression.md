# 諮問ブリーフ: V7(f) (潜熱変更 #10 の湿潤回帰) を影響評価で閉じてよいか

エスカレーション条件 4 (plan §6 に事前登録した検証の方法を変える)。plan: `plans/active/thermophysics-solver-owned-species-db.md` §5.1 #10・§6 V7(f)。

## 観測事実
- #10 (2026-09-28): H2O の潜熱を「気相 H2O と液相 H2O(L) の同 datum のエンタルピー差」に変更。V7 (a)–(e′) PASS 済み。
  ΔL (新 − 旧) = 150 K −465.92 J/kg、120 K −2387.19 J/kg、**200 K 以上は 3e-15 (相対、丸め)** (plan #10 行の記録)。
  旧実装は H2O 気相係数の再ハードコードで 200 K 未満を多項式外挿していた; 新は種 DB の c_p(200 K) 一定外挿 → 差は 200 K 未満だけ。
- 現行の H2O 凝縮の代表場で、液のある節点 (g > 1e-8) の温度 (2026-10-06 AWS で計測):
  case/16 Wysłouzil (run_0567 OFF・run_0561 ON・run_0482 の res_48000): 最低 207.5 K、200 K 未満 0 %;
  case/44 va3 dual-time 凝縮 (G0 入力 ic.h5・run_0524 res_200): 最低 228.0 K、200 K 未満 0 %。
- N2 凝縮 (case/34) は #10 でビット一致 (V7 記録)。
- V7(f) の事前登録: 「湿潤回帰の基準 run (正確な run パス・バイナリ・判定区間)・onset の抽出定義・g の報告量・記録する変化量を実装前に固定し (不変を合格にしない) (候補: case/44 va3 入口 Tt 分布 noneq の run_0127 系)、系列を check_quasisteady --series-csv で判定」。
  これは実装前に固定するはずだったが固定されないまま #10 が実装された (2026-09-28)。旧バイナリの場は互換性ハッシュで拒否され、`convert_species_field.py --src-latent legacy-v0` で移行が要る。

## 提案
V7(f) を「影響範囲の評価」で閉じる: (1) ΔL(T) 表 (120–400 K) と、(2) 現行の全 H2O 凝縮ケースで液のある節点の温度分布 (上の計測) を記録し、
「200 K 以上の液だけを持つ現行ケースの結果は #10 で変わらない (L の差 3e-15 相対; float 表の丸めの範囲)」と結論する。200 K 未満の液を持つケースが将来出たら、そのとき回帰を回す (残作業に条件付きで残す)。
回帰 run (旧バイナリ + 場の変換) は、差が丸めの範囲と分かっているので回さない。

## 問い
1. この閉じ方は事前登録の趣旨 (「不変を合格にしない」= 変化量を測る) に反するか。反するなら最小の回帰は何か (例: case/16 の float 表の L・dL/dT の旧新差を表のビットで比べる、等)。
2. float 経路 (condFloat 1 の潜熱表) で 200 K 以上の表値が旧新で変わる可能性 (表の構築の丸め) を、run なしで確かめる方法。
3. 他に #10 の影響を受ける経路 (dL/dT の二相熱容量・音速、Python 変換) で、200 K 以上でも差が出うるものはあるか。

読んでよいファイル: plan 上記、`solver_density_cuda/cuda_forge/condensationProperties_d.cuh` (`h2o_latent`)、`condensationTables_d.cuh`、`condensationSourceF_d.cuh`、`condensationEOS_d.cuh`、`solver_density_cuda/tools/convert_species_field.py`、`tests/unit/test_cond_latent_pair.cu`。
