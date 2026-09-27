# 諮問ブリーフ: 軸対称 CHT の流れあり検証 (Graetz 問題) の設計と事前登録 (2026-09-27)

AGENTS.md エスカレーション条件 **1** (plan §4 設計方針・§6 検証計画を新規に書く)。
対象 plan: `plans/active/boundary-cht-axisymmetric-graetz.md` (初稿。§4・§6 を読むこと)。

## 0. 読んでよいファイル

- plan 本体 (上)。親 plan `plans/accepted/boundary-cht-axisymmetric-fem2d.md` は §4 と §6 の V-ax2 だけ (`grep -n`)
- 基準解ツール `case/63.graetz_cht/graetz_ref.py`
- 入口 BC: `solver_density_cuda/cuda_forge/boundaryCond_d.cu` の `inlet_uniformVelocity_d` (743 行付近)
- 界面診断: `solver_density_cuda/conjugateWall.cpp` の `fillInterfaceDiagnostics` (158 行〜、`iface_q_eff` の定義 250 行付近)
- 流用元: `case/61.conjugate_annulus/{axcht.py,run_0005_annulus_r32s16_w20k/solverConfig.yaml}`

## 1. 観測事実 (実施済み、2026-09-27、commit 前の作業ツリー)

- `python3 case/63.graetz_cht/graetz_ref.py --selftest` → VERDICT PASS:
  march (Crank–Nicolson、nr 100/200/400) の Nu(x⁺=0.3) 3.65611/3.65662/3.65675、|Nu_∞−3.656793|/3.656793 = 1.2e-5。
  独立に解いた固有値級数 (FV 一般化固有値問題 nr=4000、60 項) と march nr400 の差: x⁺ 1e-3…0.1 で −0.0007…−0.0012 %。
  代表値 x⁺ 1e-3 / 3e-3 / 1e-2 / 3e-2 / 0.1 → 10.130 / 7.043 / 4.916 / 3.894 / 3.658。
- ellip (軸方向伝導あり、上下流断熱区間 x⁺ 0.01 ずつ、nr 60):
  Pe 720 で nx 300→600 の差は x⁺ 1e-3 で 1.4 %、3e-3 で 0.45 %、1e-2 で 0.1 %、0.1 で 0.0003 % (**x 方向は未収束**)。
  nx 600 で Pe 360 / 720 / 1440 の Nu: x⁺ 3e-3 → 7.0184 / 7.0287 / 7.0322、1e-2 → 4.9097 / 4.9116 / 4.9122、0.1 → 3.65617/3.65614/3.65613。
- まだ forge の run・格子は無い。

## 2. 設計 (plan §4 の要約) と私の判断

- R=1 mm、M 0.05、Re_D 1000、Pr 0.72、Pe 720、定数物性 (μ 4.086e-5、k 5.70e-2)、上流断熱 10 mm/加熱 172.8 mm (x⁺ 0.12)/下流断熱 10 mm。
- 入口 `inlet_uniformVelocity` + `inletProfile` の放物速度 (亜音速はエントロピー固定、ρ,p は内側リーマン不変量)。IC は Poiseuille + 線形 p + 300 K。
- 基準 = 古典 Graetz (march)。窓 x⁺ ∈ [3e-3, 0.1]。
- 有限 M の散逸・膨張冷却は ΔT=0 の対照 run を差し引いて消す (流れ固定なら T に線形)。ΔT 5/10 K の 2 水準で非線形 (密度変化) を見る。
- 共役: 固体殻 t=0.5 mm、k_s=100 (k_s/k_f≈1750、UWT に近づける)、外面 Robin h=1e8。対照に非連成の等温壁 run。
- 許容: 主判定 max_W |Nu/Nu_G−1| ≤ 2 %、x⁺∈[0.03,0.1] で ≤ 1 %。格子 16/32/64 で単調減少。連成 vs 等温壁 ≤ 0.3 %。

## 3. 仮説 (未検証の見立て)

- H1: 流れがあっても `iface_q_eff` (壁 CV の残差から壁流束を引き、蓄積 e_w R_ρ を引く) は壁熱流束の保存形の値になり、Graetz と格子誤差の範囲で合う。
- H2: 散逸の寄与は対照差し引きで 0.1 % 未満まで消える。
- H3: N_r=32 (Δr=31 μm、x⁺ 3e-3 の熱境界層に ~9 セル) で窓全体 2 % 以内。

## 4. 諮りたいこと (それぞれ推奨を 1 つに絞る)

1. **条件の選び方**: M 0.05・Re 1000・ΔT 10 K は妥当か。散逸を対照で消す方式と、M をもっと下げる方式のどちらを主にすべきか
   (M を下げると低マッハで SLAU・陰解法の収束が悪くなる懸念。`lowMachPrecond` は深いすきまで非物理な貫通流を出した前歴があり使わない予定)。
2. **許容差 2 % / 1 % / 0.3 %** は事前登録として妥当か。根拠の薄い数字があれば、どう決めるべきか (例: 格子 3 水準の Richardson 外挿値に対する許容にする等)。
3. **k_s/k_f ≈ 1750** の連成で、親 plan の D_f・interval の設定が安定する見込み。危ないなら k_s を下げる (壁温非一様を受け入れて評価で実 T_w を使う) べきか。
4. **入口 BC**: `inlet_uniformVelocity` の亜音速分岐は T を直接固定しない (エントロピー固定)。Graetz の「入口温度一様」を満たすか、
   `inlet_Pressure` (Pt,Tt) にして上流の発達区間を長く取るべきか。
5. **軸を含むこと**: 親 plan は軸を含まない格子だけで検証した。軸ノードの扱い (u_r=0 ピン、r→0 の r 重み) が混合平均・Nu に効く懸念と、切り分けに要る追加 run はあるか。
6. **見落とし**: 評価量 (T_w に `iface_Tw_bc` を使う、x⁺ を実測 ṁ で再計算) の穴、node 壁列の T/ρ 市松・fx うねりの影響、判定区間の取り方。

**両論併記で逃げず、推奨を 1 つに絞ること。根拠は `ファイル:行` か上の数値で示すこと。**
