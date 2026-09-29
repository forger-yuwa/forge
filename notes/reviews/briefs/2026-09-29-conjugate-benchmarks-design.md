# 諮問ブリーフ: 固体が効く共役熱伝達の検証 (A 厚肉管・C 共役平板) の設計 (2026-09-29)

AGENTS.md エスカレーション条件 **1** (plan §4・§6 を新規に書く)。plan: `plans/active/boundary-cht-conjugate-benchmarks.md` (初稿、§6 は骨子)。

## 背景 (事実)
- 前身 `plans/accepted/boundary-cht-axisymmetric-graetz.md` は固体 k_s = 100 W/mK (ガスの約 1750 倍) で壁をほぼ等温 (ずれ 33.5 mK) にしたので、固体の熱抵抗・軸方向伝導が効く共役は試していない (ユーザ指摘で主張を限定済み)。
- 固体 `fem2d` の外側境界は Robin (h, T_c) のみ (`solver_density_cuda/conjugate/solidFem2d.hpp`)。一様熱流束は h 1e-4・T_c ≈ 1e7 で代用する案。
- 文献 (Faghri & Sparrow 1980、Barozzi & Pagliarini 1985、Luikov 1974) は手元に無い。判定は独立参照解で行い、文献値は入手後の任意照合。
- 前身の道具: `case/63.graetz_cht/{graetz_ref.py, eval_graetz.py, compare_nu.py, temp_reproduce.py}`、評価の作法 (期待節点集合・u_it・全節点準定常・対照差し引き)。起動は `cfl_pseudo` 2。

## 諮りたいこと (推奨を 1 つに絞る)
1. **A の条件**: r_o = 2 mm (厚さ = R)、k_s/k_f = 10 と 100、Pe 720 のまま。壁の軸方向伝導が効き、かつ上流予熱が測れる組み合わせとして妥当か。文献と比べやすい無次元パラメータ (壁厚比・伝導比・Pe) にすべきか。
2. **C の流れ**: 参照解の流速を Blasius にすると forge の流れ (前縁・有限 Re・圧縮性) との差が混ざる。(a) Blasius + 窓で除外、(b) forge の流速場を参照解に入れる (流れの独立性を捨てて伝熱だけを検証)、のどちらを主判定にすべきか。
3. **許容差の作り方** (§6 の骨子の数値) と、「固体が効いていることの確認」の定量化。
4. **一様熱流束の Robin 代用**の数値的な危うさ (行列の条件・界面の D_f 収束・FP64)。
5. 有限 M の散逸と対照差し引き (前身の限界) を今回どう扱うか。M を下げる方がよいか。
6. 見落とし。
