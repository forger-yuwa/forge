# 諮問ブリーフ: 軸対称 `fem2d` の設計 (§4) と検証計画 (§6) (2026-09-27)

AGENTS.md エスカレーション条件 **1** (plan §4 設計方針・§6 検証計画を新規に書く)。

## 読んでよいファイル (巨大ファイル禁止)
- 本 plan: `plans/active/boundary-cht-axisymmetric-fem2d.md` (全文、短い)
- 親 plan: `plans/accepted/boundary-conjugate-heat-transfer.md` は**丸ごと読まない**。§4.3 (`grep -n "### 4.3"` から 50 行)、§4.4 冒頭、§5.1 #95 の行だけ
- `methods/boundary.md` の「共役熱伝達 (CHT)」節 (`grep -n "軸対称" methods/boundary.md` の前後)
- ソルバ: `solver_density_cuda/conjugateWall.cpp` (`initConjugateWalls`, `initSolidFem2d`, `fillInterfaceDiagnostics`, `updateFem2dWall`)、
  `solver_density_cuda/conjugate/solidFem2d.cpp`、`solver_density_cuda/variables.cpp`:495-640 (r 重み幾何)、
  `solver_density_cuda/tools/solid_fem2d.py`
- 軸対称の例: `case/23.axi_nozzle/run_cpg_air_axisym/solverConfig.yaml`

## 観測事実 (2026-09-27)
- 軸対称メッシュは $z\equiv0$ で `initSolidFem2d` の平面ガードを素通りする。`variables.cpp`:530 が `axisymMethod: 0` で面ベクトル・面積・体積を $r$ 倍
  (デバイス配列のみ。ホストの `msh.planes[].surfArea` は平面のまま、`variables.cpp`:506)。
- `iface_Qf_eff` = `ifaceRraw − Fw − e_w R_ρ` は r 重み残差由来 → [W/rad]。`iface_q_eff` はそれをホストの平面面積で割る (`conjugateWall.cpp`:229-238)。
  `fem2d` の固体は平面の集中辺長 [m] で組む。→ **軸対称では荷重/熱流束が節点ごとに r 倍ずれたまま黙って連成していた** (コード読みで確定。数値 run は無し)。
- 対処 (commit `1a389fef`): `initConjugateWalls` で軸対称を全モード拒否。ローカル 1 step で rc=1・固有メッセージを確認、平面はガード通過。

## 設計 (plan §4 の要約) と検証計画 (§6)
plan を読むこと。要点: 単位は流体に合わせてラジアンあたり、固体弱形式に r 重み (剛性は要素重心 r で厳密、Robin 辺は線形 r の consistent 形、界面集中量 $A_i^r=\sum L(2r_i+r_j)/6$)、
平面ではビット同一の経路、`iface_q_eff` の分母を流体面の r 重み面積に、受理は `axisymMethod: 0` + `fem2d` のみ。
検証: V-ax0 拒否の負例 / V-ax1 固体単体 (C++↔Python 1e-12、厚肉円筒殻の解析解で収束率 ≥1.8・16 層 ≤0.1 %、平面のまま組むと FAIL することの確認) /
V-ax2 同心円環 (静止ガス + 厚肉殻 + 外面 Robin、軸を含まない、$T_w$ ≤0.5 % of 降下・G-cons・連成保存・G-if・準定常・層数感度) / V-ax3 平面回帰 (ビット一致)。

## 仮説・懸念 (私の)
1. 流体側の $r$ 重みは**面重心の $r$** (`r_face`) を掛ける一次近似で、固体の $A_i^r$ (形状関数重み) と定義が違う。荷重は積分済みを渡すので保存性は保つが、
   **熱流束の換算 ($D_f=g_fA_i^r$ や報告値) で差が出る**。$D_f$ は固定点に効かないので許容と読んでいるが、報告値は両方明記とした。
2. `iface_q_eff` の分母変更は CHT 以外の軸対称 `interfaceDiag` 利用 (ノズル壁熱流束の報告) の値を変える — 過去の報告が $r$ 倍誤っていた可能性 (§8-1)。
3. 同心円環は静止ガスで流体残差が機械ゼロ近傍になり、擬似時間の緩和が遅い (case/52 は低圧で対処)。軸対称で同じ手が効くかは未確認。

## 諮りたいこと (推奨を 1 つに)
1. §4 の設計に**誤りや見落とし**はあるか (特に単位 [W/rad] の一貫性、$D_f$・G-if の $A_i$ の選択、軸上節点の扱い、平面ビット同一の保証方法)。
2. §6 の合格条件は**事前登録として十分か** (閾値・検出力・比較相手)。足す/削る項目を挙げる。
3. §8-1 (既存の軸対称 `interfaceDiag` 報告値の誤り可能性) を本 plan で扱うべきか、別件にすべきか。
4. 拒否ガード (済) の範囲は妥当か (全モード拒否。`interfaceDiag` 単独の軸対称利用は止めていない)。
