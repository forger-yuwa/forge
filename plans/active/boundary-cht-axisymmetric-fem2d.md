# 軸対称の固体伝導 (`fem2d`) — ノズル壁の共役熱伝達

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) 「共役熱伝達 (CHT)」— 軸対称は現状拒否 (2026-09-27) と、本計画の設計 (「軸対称の `fem2d` (計画中)」)
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md) — 軸対称の $r$ 重み幾何 (`axisymMethod: 0`)
- **related_plans**:
  - [`../accepted/boundary-conjugate-heat-transfer.md`](../accepted/boundary-conjugate-heat-transfer.md) — **親** (平面 2D の `fem2d`、界面契約・G-if・G-cons)。§5.1 #95 を本計画に切り出した
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) — ノズル等温壁チェーン (壁温を仮定している側)
- **created**: `2026-09-27`
- **owner**: `sano`

## 1. 目的

`mode: fem2d` を**軸対称 (`mesh.isAxisymmetric: 1`, `axisymMethod: 0`)** で使えるようにし、
**ノズル壁の共役熱伝達** (壁厚 + 外表面に熱伝達率 $h$ と $T_c$ を貼る = 固体の `hole` Robin) を解けるようにする。
いま `case/40` 等は壁温を**仮定**して生産値を出しているので (壁温 1400 ± 15 K を仮定)、それを**解**にするのが目的。
完了時: 円筒殻の解析解 (対数抵抗) に一致する軸対称 CHT がソルバ内で回り、平面 2D の結果は変わらない。

## 2. スコープ

- **やる**:
  - 固体の弱形式に $r$ 重みを入れる (要素剛性・Robin 辺・界面の集中量)。C++ (`conjugate/solidFem2d.cpp`) と参照オラクル (`tools/solid_fem2d.py`) の両方
  - 軸対称の起動時検査 (`axisymMethod: 0` のみ受理、固体節点 $r\ge0$、界面の集中量 $A_i^r>0$)
  - `iface_q_eff` (熱流束の診断・`local1d` の入力) の面積を軸対称で $r$ 重みにする — **現状は平面面積で割っており $r$ 倍ずれる** (§3)
  - 検証: 固体単体の円筒殻解析解、ソルバ内連成の同心円環 (case/52 の軸対称版)、平面 2D の回帰
- **やらない**:
  - `axisymMethod: 1` (SU2 流 planar + $1/y$ ソース) — 流体の界面荷重の単位が別物になるので**拒否のまま**。必要になったら別 plan
  - `local1d` の軸対称 (薄肉の対数抵抗化) — **拒否のまま**。ノズル壁は `fem2d` で扱う
  - ノズル生産ケース (`case/40` 等) への適用 — 検証後に別 plan (壁厚・外表面 $h$ の出典決めが要る)
  - 3D (`fem2d` の 3D、親 plan §5.1 #96)、放射、dual-time

## 3. 関連 docs と前提 (2026-09-27 に実測・コードで確認)

- **軸対称は黙って間違っていた → 拒否を入れた (本計画 §5.1 #1、済)**。軸対称メッシュも $z\equiv0$ なので
  `initSolidFem2d` の平面ガード (`conjugateWall.cpp`) を素通りする。`axisymMethod: 0` では
  `variables.cpp`:530 が面ベクトル・面積・体積を $r$ 倍するので、流体の残差 `ifaceRraw` と積分済み荷重 `iface_Qf_eff` は
  **[W/rad]** ($2\pi$ を掛けない)。一方
  - `fem2d` の固体は平面 [W/m] の集中辺長で組む → **荷重が節点ごとに $r$ 倍ずれる**
  - `iface_q_eff` = `iface_Qf_eff` / `msh.planes[ip].surfArea` の分母は**ホストの平面面積** (`variables.cpp`:506 で $r$ を掛けるのはデバイス配列だけ) → **$q$ が $r$ 倍ずれる**。`local1d` はこれを使う
  → `initConjugateWalls` で**全モードの軸対称を拒否**した。ローカル 1 step で拒否を確認 (case/58 入力に `isAxisymmetric: 1`、rc=1)、平面はガードを通過。
- 親 plan の記述「軸対称の $r$ 重みは host・double で行う」(`methods/boundary.md`) は**未実装だった** (訂正済み)。
- 流体の $r$ は `y` 座標 (軸は $x$ 軸、$r=y\ge0$)。`axisRFloor` 以下の面は $r$ を床に上げる。

## 4. 設計方針

### 4.1 単位 — ラジアンあたりに揃える

流体の $r$ 重み幾何は $2\pi$ を掛けない**ラジアンあたり**なので、固体も [W/rad]・[m²/rad] で組む。
界面荷重 `iface_Qf_eff` は**そのまま**渡す (親 plan の「積分済み荷重を面積で割って掛け直さない」を維持)。

### 4.2 固体の弱形式 ($r=y$)

$$\int_\Omega k\,\nabla T\cdot\nabla v\,r\,dA \;+\; \int_{\Gamma_R} h\,(T-T_c)\,v\,r\,ds \;=\; \sum_i Q_{f,i}\,v_i$$

線形三角形・要素内 $k_e$ = 節点 $k_s(T)$ の平均 (現行規約) として

- **剛性**: 勾配が要素内一定なので $K^e_{ij}=k_e\,\bar r_e\,A_e\,\nabla N_i\cdot\nabla N_j$ ($\bar r_e$ = 要素重心の $r$)。**厳密**。
- **Robin 辺** ($r$ は辺上で線形): $M^e=\frac{hL}{12}\begin{pmatrix}3r_a+r_b & r_a+r_b\\ r_a+r_b & r_a+3r_b\end{pmatrix}$、
  荷重 $b^e=\frac{hT_cL}{6}\,(2r_a+r_b,\;r_a+2r_b)$。**厳密**。
- **界面の集中量**: $A_i^r=\sum_{e\ni i}\frac{L_e}{6}(2r_i+r_j)$ [m²/rad] (= $\int N_i\,r\,ds$)。
  使い道は親 plan と同じく**熱流束への換算と $D_f=g_f A_i^r$ だけ** (荷重には使わない)。
- **最終積分に $r$ を掛けるだけでは足りない** (親 plan §4.3 の codex M6)。要素剛性と Robin 辺の両方に入れる。
- 平面 ($r\equiv1$ と置いたもの) と**同じコード経路**にし、`isAxisymmetric == 0` ではビット同一にする
  (重みを掛ける分岐を入れず、平面では重み配列を作らない)。

### 4.3 界面熱流束の診断 (`iface_q_eff`)

軸対称では分母を**流体面の $r$ 重み面積** $\sum_f |S_f|\,\max(r_f, r_{\rm floor})$ (`variables.cpp` と同じ床) にする。
$Q$ 側 (`iface_Qf_eff`) は変えない。**平面では変えない**。
流体面の $r$ 重み面積と固体の $A_i^r$ は定義が違う (半割面の重心 $r$ vs 一次の形状関数重み) ので、
**熱流束の比較は両方を明記して行う** (G-if の ① は固体側 $A_i^r$、壁熱流束の報告は流体側)。

### 4.4 起動時検査

受理: `isAxisymmetric: 1` かつ `axisymMethod: 0` かつ `mode: fem2d`。
拒否: `axisymMethod: 1`、`mode: local1d`、固体節点に $r<0$、界面節点に $A_i^r\le0$ (軸に接する界面)。
固体内部の軸上節点 ($r=0$) は許す (要素重心は $r>0$ なので剛性は正定値のまま。Robin 辺が軸上だけにある構成は零空間になるので既存の「Robin 辺 1 本以上」検査に加え $\sum hLr>0$ を見る)。

### 4.5 参照オラクル

`tools/solid_fem2d.py` の `Fem2DOperator` に `axisym=False` 引数を足し、同じ式で組む。
`test_solid_fem2d_cpp.py` に軸対称の組立て・求解の突き合わせを足す (§6 V-ax1)。

## 5. 実装ステップ

1. `conjugateWall.cpp` — 軸対称の拒否 (**済**、2026-09-27)
2. `conjugate/solidFem2d.{hpp,cpp}` — 重み ($r$) を持たせ、剛性・Robin・集中量に適用。`solid_fem2d_tool` に軸対称フラグ
3. `tools/solid_fem2d.py` / `test_solid_fem2d_cpp.py` — オラクルと突き合わせ
4. `conjugateWall.cpp` — 拒否を §4.4 の受理条件に置き換え、`fillInterfaceDiagnostics` の分母を §4.3 に
5. `case/59.conjugate_annulus/` — 同心円環の検証ケース (メッシュ生成・固体生成・評価器)
6. `methods/boundary.md` を「計画中」から現行仕様へ

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | **軸対称の拒否** (済 2026-09-27) | `initConjugateWalls` で `isAxisymmetric == 1` を全モード拒否。ローカル 1 step で rc=1 を確認、平面はガード通過 | O |
| 2 | **§4 / §6 を上位に諮る** | AGENTS.md 条件 1。codex (diagnose) → codex plan レビュー | F |
| 3 | 固体の $r$ 重み (C++ + オラクル) | §5 ステップ 2–3。合格: §6 V-ax1 | O |
| 4 | 受理条件と `iface_q_eff` の分母 | §5 ステップ 4。合格: §6 V-ax0 | O |
| 5 | 同心円環ケースと判定 | §5 ステップ 5。合格: §6 V-ax2 | O |
| 6 | 平面の回帰 | 合格: §6 V-ax3 | O |

## 6. 検証 (事前登録。結果を見る前に書いた — 2026-09-27)

合否はすべて **全域 FP64 ビルド**で取る (親 plan #94 と同じ理由: float32 の深部精度床)。run は AWS。

- **V-ax0 起動時検査**: `axisymMethod: 1` / `mode: local1d` / 固体節点 $r<0$ / 軸に接する界面 がそれぞれ拒否され (rc≠0 と固有のメッセージ)、
  受理構成は起動する。負例は `case/59` の入力を 1 か所ずつ変えて 1 step。
- **V-ax1 固体単体** (流体なし、`solid_fem2d_tool` と Python):
  - (a) 組立て: C++ と Python の $K$・$b$ の相対差 ≤ 1e-12、同じ荷重での解の差 ≤ 1e-9 K
  - (b) **円筒殻の解析解**: 内面 $r_1$ に一様熱流束 $q_1$ (界面荷重として与える)、外面 $r_2$ に Robin ($h$, $T_c$)。
    $T(r)=T_c+q_1r_1\left[\dfrac{1}{h\,r_2}+\dfrac{\ln(r_2/r)}{k}\right]$ (ラジアンあたりの熱量 $q_1r_1$ が一定)。**厚肉** $r_2/r_1=2$ (平面近似との差 ≈ 30 %) で、
    半径方向 4/8/16 層の誤差が**収束率 ≥1.8**、16 層で内面温度誤差 ≤ 0.1 % of 殻の温度降下。
    **平面のまま組んだ場合に FAIL することも確かめる** (検出力の確認)。
  - (c) 収支: 界面入熱 = Robin 持ち去り (相対 ≤ 1e-12)
- **V-ax2 ソルバ内連成 — 同心円環** (`case/59.conjugate_annulus`、case/52 の軸対称版):
  静止ガス ($k_f$ 一定、低圧、重力なし) を $r_a\le r\le r_b$ に置き、内壁 $r_a$ を等温 $T_a$ (非連成)、外壁 $r_b$ を**共役壁**、
  固体殻 $r_b\le r\le r_c$、外面 Robin ($h$, $T_c$)。軸方向の両端は slip (1 次元性)。**軸を含まない** (軸の特異性を混ぜない)。
  解析解 (ラジアンあたり) $Q'=\dfrac{T_a-T_c}{\ln(r_b/r_a)/k_f+\ln(r_c/r_b)/k_s+1/(h r_c)}$、$T_w=T_c+Q'\left[\ln(r_c/r_b)/k_s+1/(hr_c)\right]$。
  **厚肉** ($r_c/r_b$ = 2 程度) にして平面近似との差を 10 % 以上にする。合格:
  - $|T_w-T_{w,*}|\le$ **0.5 % of 固体の温度降下** (全界面節点の max)、事前予測 ≤0.2 % (親 plan V1 のソルバ内 0.152 %)
  - G-cons (固体入熱 = Robin 持ち去り = 流体の内壁放熱) ≤ 0.5 %、連成保存 $|Q_{\rm sol}-Q_f|/A_i^r\le$ 0.1 % of $q_*$
  - G-if (`check_cht_interface.py`、全域、登録値 `eps_abs` / `eps_rel` 1e-3 / `dT_K` 1e-2 / `n_consec` 80) PASS
  - 壁温 $T_w$ と $q$ が `check_quasisteady.py` で STEADY (`--drift/--osc` 0.001)
  - 流体残差は静止場で初期から機械ゼロ近傍になりうる (case/52 と同じ) ので、`check_convergence` の NOT CONVERGED は報告項目
  - **感度**: 固体の半径方向 8/16 層で $T_w$ 差 ≤ 0.1 % of 降下
- **V-ax3 平面の回帰**: 変更後のバイナリで `case/52` の V1 ソルバ内 (`run_0003` 設定) と `case/58` の基準 (`run_0012_v6p_df5_i50_nlog` 設定、10k step) を回し、
  **固体・界面の出力が変更前バイナリとビット一致** (平面では重み経路に入らない設計なので、atomicAdd の揺らぎを除けば一致するはず。
  一致しなければ同一バイナリ 2 回のノイズ床と比べ、床の内側なら合格とする)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/conjugate/solidFem2d.{hpp,cpp}`、`conjugate/solid_fem2d_tool.cpp`、`conjugateWall.cpp`
- `solver_density_cuda/tools/solid_fem2d.py`、`test_solid_fem2d_cpp.py`、`check_cht_interface.py` (必要なら)
- 新規 `case/59.conjugate_annulus/`
- 平面 2D の既存ケース (case/48・52・53・54・58) は §4.2 の設計で不変 (V-ax3 で確認)
- `methods/boundary.md`、`methods/index.md`

## 8. 未確定事項

1. `iface_q_eff` の分母を変えると、**軸対称で `interfaceDiag: 1` を使っている既存の診断** (CHT 以外、例えばノズル壁熱流束の報告) の値が $r$ 倍変わる。
   既存 run の報告値が誤っていた可能性がある — 影響範囲の洗い出しが要る (該当 run を grep で列挙し、誤っていたなら撤回の要否を決める)。

## 9. 変更ログ

- `2026-09-27` — 初稿 (親 plan §5.1 #95 から切り出し、ユーザ指示「起票していいよ」)。軸対称が**黙って間違っていた**ことをコードで確定し、全モードで拒否 (§5.1 #1)。§4/§6 は上位に諮る前の草稿。
