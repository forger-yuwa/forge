# ③ベル: 特性曲線で辿るモード (TP の Rao) と形状を直接変えるモード (C2 B-spline) の比較

## メタ

- **area**: `tooling / optimization`
- **status**: `draft`  <!-- 2026-10-11 起票。同日 codex diagnose で NO-GO → 指摘を反映した改訂版。実装前に --stage plan のレビューを受ける -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (実装着手前に「ベルの 2 モード」節を足す — §5 ステップ 0)
  - [`notes/investigations/nozzle-top-internal-shock-diagnosis.md`](../../notes/investigations/nozzle-top-internal-shock-diagnosis.md) (TOP の内部衝撃波。Rd を変えても強さはほぼ不変)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) §3 確定事項 5〜11、§5.1 #3・#4
  - 兄弟: [`tooling-nozzle-ambient-separation.md`](tooling-nozzle-ambient-separation.md) (**別セッションの担当**、2026-10-11 ユーザ決定。本計画は外部領域と剥離評価が使える前提で進め、海面作動の比較 §6 V4 はその受理後。推力の会計は同じ定義)
  - 資産元: [`../accepted/tooling-nozzle-moo-loop.md`](../accepted/tooling-nozzle-moo-loop.md) (Phase 2)、[`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) §4.3 (平面の縁条件 `rao_planar.py`)・§4.5 (凍結 TP `gas/frozen.py`)、[`../accepted/tooling-nozzle-cfd-pinned-initial-line.md`](../accepted/tooling-nozzle-cfd-pinned-initial-line.md) (CFD ピンの初期線)
- **created**: `2026-10-11`
- **owner**: `sano`

## 1. 目的

③ベルの壁を決める**設計手順**として、(M) 特性曲線法で Rao の最大推力輪郭を作る手順と、(D) 壁の形を少数の幾何変数で直接振る手順を、
同じ評価器 (TP のロケット排気・同じ物理寸法・同じ格子規約・同じゲート・同じ推力定義) で比べる。比べるのは「完成した設計手順」で、
設計点 (真空・付着流) の $\eta$–$L$、形状が作れない点の率、較正・探索を含む総費用を比べ、外部領域の受理後は海面作動の剥離と内部衝撃波も比べる。
結果で、ベル (と④デュアルベルの基部) の設計モードを決める。

## 2. スコープ

- **やる**
  - ③の評価器を TP のロケット排気にする (基準ケース: LOX/LH2・燃焼室圧 10 MPa・O/F 6・面積比 $\varepsilon$ 50、CEA の燃焼室平衡組成を凍結、断熱壁)
  - 両モード共通の C2 の壁 (スロートの曲率は 1 つの値 $1/R_t$) と、共通の寸法決定 (§4.2)
  - (M) 凍結 TP の熱力学 ($h(T)$, $s(T,p)$) で導いた軸対称 Rao 条件と生成器、曲率連続の kernel、端点を任意に指定できる当てはめ、δ\* の一発補正
  - (D) 拡大部の 5 次 B-spline (内部の形を 2 変数で振る)
  - 比較の run、TOP (Phase 2) の参照点の再評価
- **やらない**
  - 海面作動の比較の実装 (兄弟 plan の外部領域を待つ。§6 V4 は登録だけ)
  - 冷却壁 (当面断熱、親 §3-10)、平衡流 (凍結のみ)
  - 壁静圧分布を指定する逆 MOC、TIC、軸マッハ分布 (M 側の追加の族として後で検討)

## 3. 関連 docs と前提

- **Phase 2 の実績**: TOP ($\theta_n,\theta_e,L$) の MOO で 77 評価、Rao 80 % のチャート点が同じ長さの最適 $\eta$ の 0.133 % 下 (判定帯 ±0.15 %)。固定長で $\theta_n$ 24.5〜28.5°・$\theta_e$ 8〜11° が平坦な台地。評価器は CPG・`viscMethod: 1` (空気の Sutherland)・SST・断熱壁・出口面で領域を打ち切り (`evaluate/runner.py:90,103,124`)。
- **既存の部品**: 軸対称・semi-perfect の MOC 単位過程 `geometry/moc_kernel.py`、Hall の遷音速級数 `geometry/transonic.py`、CFD ピンの初期線 `feedback/cfd_initial_line.py` (凍結源の形の照合は U→T 5 次 Hermite が前提)、位置と壁角の同時当てはめ `geometry/wall_axismach.py:999` `joint_fit_wall` (**始点の拘束は $r'=x_0/R,\ r''=1/R$ に固定** — `wall_axismach.py:1061`)、δ\* の積分法 `feedback/deltastar_integral.py` (粘性係数は空気の Sutherland、温度は回復係数の近似)、RANS からの δ\* 抽出 `metrics/deltastar.py`、凍結 TP `gas/frozen.py`、平面の縁条件 `geometry/rao_planar.py` (自前導出扱い — 同ファイル 12 行)。
- **輸送物性**: `viscMethod: 1` は空気の Sutherland。燃焼ガスには `viscMethod: 2` (種ごとの輸送物性、`physProp.transport` 必須、`thermalMethod: 2` のみ — [`procedures/solver-settings.md`](../../procedures/solver-settings.md) の `physProp.viscMethod`) を使う。
- **Rao の軸対称の条件**: Cain 2010 式 4.1–4.3 (軸対称形)、Shyne & Keith 1990 式 11–15 (平面版)。`papers/nozzle_design/` (git 管理外)。CPG の導出は TP の最適性の保証にならない (codex 2026-10-11 B2)。
- **内部衝撃波**: case/40 で TOP の Rd を 0.4→1.5 にしても内部の ΔM は −0.43→−0.40 でほぼ不変。C2 にしても内部衝撃波が消える保証は無い。

## 4. 設計方針

### 4.1 評価器 (両モード共通)

- **ガス**: CEA で燃焼室平衡組成を求めて凍結し、NASA-9 の擬似種 1 つ (`EXH`) の TP (`thermalMethod: 2`、`thermoHrefTemp: 298.15`)。輸送物性は `viscMethod: 2` で CEA 形 (lump は構成実種ごと)。全温は CEA の燃焼室温度。
- **壁**: 断熱、SST 低 Re ($y_1^+\le1$、`wallTreatmentSST: 0`)。SST の項は推奨レシピ。
- **作動点 (第 1 段)**: 真空の設計点。領域は出口面までで、出口は超音速の流出。
- **推力**: ノズル内壁・入口面・出口面で閉じた検査体積 CV から 2 つの会計を作る。壁会計 = 入口面の (運動量流束 + $(p-p_a)$ の面積積分 + 粘性応力) + 内壁の $(p-p_a)$ と摩擦の軸方向積分。出口面会計 = 出口面の (運動量流束 + $(p-p_a)$) の符号付き積分。2 つの差を CV 残差とする (兄弟 plan §4.4 と同じ定義。本計画はこの定義を自前で持つ)。$\eta_{vac}=F_{vac}/F_{1D,vac}(\varepsilon)$ ($F_{1D}$ は同じ凍結 TP で $\varepsilon$ まで等エントロピー膨張した 1D 推力)。
- **ゲート**: Phase 2 の物理ゲート ($\dot m/\dot m_{1D}$、$\eta$ の範囲) + `check_convergence` / `check_quasisteady` (`STEADY` 以外は判定不能) + $\eta$ の頭打ちの polish ($|\Delta\eta|<3\times10^{-4}$) + CV 残差 ≤ 0.1 %。
- **内部衝撃波の指標**: 診断ノートの方法 (各 $x$ 断面で Mach の最大降下の位置を追う) で最大の $|\Delta M|$ を記録する。
- **Euler の設定** (V2 で使う): `viscMethod: 0`、`visc: 0.0`、`thermCond: 0.0`、`turbulence.model: none`、物理壁はすべて `slip`。

### 4.2 共通の壁の構成と寸法 (C2)

- 入口直管 → 収縮部 (U→T 5 次 Hermite、スロートで $r''=1/R_t$。①の axis-Mach チェーンと同じ部品) → スロート → 拡大部。拡大部の作り方だけがモードで違う。
- **寸法の spec** (両モード共通、キャンペーン内で固定): スロート半径 $r_t$ (物理単位)、$R_t$、収縮比、**物理壁の** $\varepsilon$ と $L$。値は §8-1。
- 曲率は $\kappa=r''/(1+r'^2)^{3/2}$ と $r''$ を区別して扱う (スロートでは $r'=0$ なので一致)。

### 4.3 (D) 直接モード

拡大部 $r(x)$, $x\in[0,L]$ を 5 次 B-spline にし、係数は「拘束を満たし $\int(r''')^2dx$ を最小にする」線形問題で一意に決める。拘束は

- スロート: $r=1,\ r'=0,\ r''=1/R_t$
- 最大壁角の点 $x_m$: $r'=\tan\theta_{max},\ r''=0$
- 出口: $r=\sqrt\varepsilon,\ r'=\tan\theta_e$

**設計変数 = $(L,\ \theta_e,\ \theta_{max},\ x_m/L)$ の 4 個**。成立条件: $r'\ge0$、壁角の極大が 1 つ、曲率の上限。破れたら INFEASIBLE。この壁がそのまま物理壁 (δ\* 補正なし)。

### 4.4 (M) 特性曲線モード (TP の Rao)

- **初期線**: スロートの C⁻ 特性線。場の出所は **CFD ピン** (同じ収縮部・$R_t$・ガスの Euler 場を 1 本回して凍結) を既定、Hall は検算。凍結源の形の照合を本計画の収縮部と $R_t$ に合わせて広げる。
- **kernel**: スロートから終点 B までの膨張部。曲率がスロートの $1/R_t$ から連続に変わる曲線族で、B での曲率 $\kappa_B$ をパラメータに持つ。古典の円弧 kernel も同じ生成器で選べるようにする (回帰試験用)。
- **Rao 条件 (TP)**: 凍結 TP の $h(T)$・$s(T,p)$ を使い、質量流量と長さの拘束付きで推力の変分をとって制御面 DE の条件を導く (S1、§5.1 #4)。比熱比を陽に含まない形に書けるかは導出の結果で判断し、仮定しない。与えた $(L,\varepsilon)$ に対して B と D を求め、領域 BDE を MOC で埋め、質量流量で壁流線 BE を取る。
- **C2 接続**: 符号付きの関数 $g(\kappa_B)=\kappa_{Rao}(B^+;\kappa_B)-\kappa_B$ の根を求める ($\kappa_{Rao}(B^+)$ は kernel を $\kappa_B$ で作ったときの Rao 輪郭の出だしの曲率)。**各評価で B・D を解き直し**、$g$ の符号を挟めない候補は「C2 の成立未確認」として台帳に残し、評価に回さない。BE の MOC 壁点は、**始点の $r,r',r''$ と終点の $r,r'$ を任意に指定できる**よう `joint_fit_wall` を拡張した当てはめで B-spline にする。B を跨いで当てはめ直す方式は採らない。
- **δ\* と物理寸法**: 積分法の粘性係数を評価器と同じ $\mu(T)$ に替え、摩擦係数の倍率 $k_f$ を RANS 抽出で 2〜3 形状に合わせて固定する (温度の回復の近似と較正の誤差を記録)。非粘性壁を法線にずらしたあと、**物理壁の $\varepsilon$・$L$・スロート曲率を拘束した当てはめ**で最終壁を確定する (相似倍率では $\varepsilon$ が変わらないので使わない)。帰還ループは回さない。
- **設計変数 = $L$ の 1 個**。Rao 族は長さごとに非粘性の最大推力輪郭が 1 本なので MOO はせず、$L$ を掃引する。

### 4.5 比較の手順

- 両モードを §8-1 の同じ寸法 spec で回す。比較用の**共通の長さ $L_k$ を 5〜7 点**決め、(M) はその $L_k$ で、(D) は MOO (LHS 40 点 + EHVI infill 20〜30 点) の後に各 $L_k$ で $L$ を固定した $(\theta_e,\theta_{max},x_m)$ の再最適化 (Phase 2 の固定長最適化と同じ) を行って、同じ $L_k$ の値どうしを比べる (前線の補間で比べない)。
- TOP (Phase 2 の $\theta_n,\theta_e,L$) を 3 点だけ新しい評価器で回し、過去の結果との橋にする。
- 比べる量: $\Delta\eta(L_k)=\eta_D(L_k)-\eta_M(L_k)$、形状が作れない点の率、総費用 (評価回数・時間と、(M) の CFD ピン・$k_f$ 較正の run を含める)、内部衝撃波の指標、(M) の $g$ の根の残差と当てはめ残差。

## 5. 実装ステップ

0. `methods/design/overview.md` に「ベルの 2 モード」節を書く (§4 が plan 段レビューを通った後、実装前)
1. 共通の推力定義 (兄弟 plan §4.4)・寸法 spec・輸送物性を確定
2. 評価器の TP 化 (CEA 組成・`EXH`・`viscMethod: 2`・断熱・低 Re 格子)
3. (D) の壁生成と単体試験
4. (M) の Rao 条件の導出 (S1) → 生成器 → 単体試験 (§6 V1)
5. CFD ピンの凍結源と Hall 検算
6. 当てはめの拡張、δ\* の積分法の $\mu(T)$ 化と $k_f$ の較正、拘束付きの物理壁の確定
7. 比較の run (case/68、新設) と判定

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~§4・§6 を上位に諮る~~ | **判断: 2026-10-11 codex (diagnose) [記録](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) — NO-GO。C1・B1〜B4 を採用して改訂** (採否は §6.1)。着手順 TP 化 → D → M は採用 | F |
| 1b | 寸法 spec (仮値) で plan 段レビュー | §8-1 の仮値で codex `--stage plan` (2026-10-11 実施) | F |
| 2 | 評価器の TP 化 | §4.1。合格: 1D の凍結 TP の出口温度が CEA の凍結計算と相対 0.1 % 以内、Phase 2 の TOP 1 点を旧設定 (CPG・`viscMethod: 1`) のまま回して旧 $\eta$ と 3e-4 以内 (配管の回帰)、CV 残差 ≤ 0.1 % | O |
| 3 | (D) の壁生成 | §4.3。合格: 拘束の残差 ≤ 1e-10、C2、同一入力でビット同一、成立条件の負例試験 | O |
| 4 | (M) TP の Rao 条件の導出 | §4.4。凍結 TP の $h(T)$・$s(T,p)$ で拘束付き変分。平面縮約が `rao_planar.py` と一致するかも確かめる | F |
| 5 | (M) 生成器と単体試験 | §6 V1 | O |
| 6 | CFD ピンと Hall 検算 | 合格: 凍結源の照合 PASS、Hall との初期線の差を記録 (合否なし) | O |
| 7 | 当てはめの拡張・δ\*・物理壁 | §4.4。合格: 拡張した当てはめが旧 `joint_fit_wall` の拘束を与えたときビット一致、物理壁の $\varepsilon$・$L$ が spec と 1e-6 以内 | O |
| 8 | 比較の run と判定 | §6 V2・V3。case/68 の README に run 一覧 | O (実行) / F (解釈) |
| 9 | 海面作動の比較 | 兄弟 plan の外部領域が受理されてから。§6 V4 に参照・判定量を書いてから回す | F |

## 6. 検証

**判定基準は run の前にここで固定する。** codex diagnose (2026-10-11) の指摘を反映した案で、plan 段レビューを受けて確定する。

- **V1 (M の生成器)**:
  - 回帰: CPG・円弧 kernel で古典 Rao を再現 (Rao チャートの $\theta_n,\theta_e$ と ±1°、チャートの比熱比に合わせる)。TP 生成器を定比熱にすると CPG 生成器と 1e-8 で一致。
  - 保存: 制御面 DE を通る流量がスロートの流量と相対 1e-4 以内。制御面の推力と「入口 + 壁」の推力が相対 1e-4 以内。
  - **TP の停留性**: $c_p$ が実際に変わる条件 (基準ケースの EXH) で、拘束 ($\dot m$、$L$、$\varepsilon$) を保つ壁の摂動族を作り、MOC で求めた推力の一次変分が 0 であること (摂動振幅 $\epsilon$ に対して推力変化が $O(\epsilon^2)$ になる)。任意の形状に勝つことを最適性の証明として使わない。
  - C2 接続: 採用した候補で $|g(\kappa_B)|\le10^{-6}$ (曲率は $1/r_t$ 単位)。
- **V2 (評価器の整合)**: (M) の非粘性壁を §4.1 の Euler 設定で回し、推力が MOC の値と相対 0.2 % 以内 (粗い整合確認)。
- **V3 (比較、事前登録の分類)**:
  - 各 $L_k$ で $\Delta\eta(L_k)$ の不確かさ $U_k$ = 格子 (2 形状×2 水準の差の大きい方) + 反復 (polish の閾値 3e-4 を両モードぶん) + 探索 ((D) の固定長再最適化で、最良と 2 番目の差)。
  - **$|\Delta\eta|-U_k>1.5\times10^{-3}$ → 「差あり」** (符号で優位なモードを記録)。**$|\Delta\eta|+U_k\le1.5\times10^{-3}$ → 「1.5e-3 の工学的同等幅の中で同等」**。他は判定不能。1.5e-3 は Phase 2 の Rao 照合の判定帯と同じ値を工学的な同等幅として使う。
  - 併記 (合否なし): 形状が作れない点の率、総費用、内部衝撃波の指標、$g$ と当てはめの残差。
- **V4 (海面作動)**: 兄弟 plan の外部領域の受理後に、参照・判定量 (剥離位置、海面の推力、内部衝撃波) をここに書いてから回す (未登録)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-11` | [`2026-10-11-thruster-nozzle-child-plans-diagnose.md`](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) (兄弟 plan と一括) | NO-GO, C1/M4 (本 plan 分: 共通 C1・B1〜B4) | **C1 採用** (兄弟 plan §4.4 の閉じた CV を共通定義に)。**B1 採用** ($g(\kappa_B)$ を明示、各評価で B・D を解き直す、当てはめの端点拘束を一般化 — `wall_axismach.py:1061` の始点固定を確認、$\kappa$ と $r''$ を区別)。**B2 採用** (凍結 TP の $h,s$ で導出、TP の一次変分試験、CPG は回帰)。**B3 採用** (相似倍率では $\varepsilon$ が不変なのは自明 — 拘束付き当てはめで物理壁を確定、「完成した設計手順の比較」、総費用を比較)。**B4 採用** (Euler を明示、共通 $L_k$ で比較、$U_k$ の構成、「差なし」→「同等幅内で同等」)。前提への異議 (`viscMethod: 1` は空気の Sutherland で定数ではない) を採用し §3 を訂正。0.1〜0.2 % の予測は根拠不足として撤回 |

## 7. 影響範囲

- 新規: `design/forge_design/geometry/rao_axisym.py` (仮)、`geometry/wall_bell_c2.py` (仮)、`case/68.bell_two_mode/` (新設)
- 変更: `evaluate/runner.py` (TP・CEA 組成・`viscMethod: 2`・低 Re 格子・Euler 設定)、`geometry/wall_axismach.py` (当てはめの端点拘束の一般化。既定の挙動は不変)、`feedback/deltastar_integral.py` ($\mu(T)$)、`feedback/cfd_initial_line.py` (照合の形)、`probdef.py`、`design/CAPABILITIES.md`、`methods/design/overview.md`
- ソルバは変更しない予定

## 8. 未確定事項

1. **寸法 spec** — **仮値で plan 段レビューに出す (2026-10-11)**: スロート半径 $r_t=0.14$ m (推力 1000 kN 級、$A_t=F/(p_c C_F)$)、$R_t=1.0\,r_t$ (0.4 $r_t$ [TOP の Rd] は感度として 1 組)、収縮比 3、比較の長さ $L_k$ = 等面積比 15° 円錐の長さの 60・70・80・90・100 %。ユーザが別の値を決めたら差し替える
2. (M) の kernel の曲線族と、B・D の反復の組み方 (S1 で導出)
3. 輸送物性の出所 (CEA の trans.inp の種がそろうか、lump の扱い)

## 9. 変更ログ

- `2026-10-11` — codex diagnose (NO-GO) を反映して改訂: 推力を兄弟 plan の閉じた CV に、Euler 設定の明示、`viscMethod: 2` (燃焼ガスの輸送物性)、TP の Rao を凍結 TP の $h,s$ で導出し一次変分で検査、C2 接続を符号付き $g(\kappa_B)$ の根探しに、当てはめの端点拘束の一般化、δ\* 後の物理壁を拘束付き当てはめで確定 (相似倍率は撤回)、比較を共通 $L_k$ で行い不確かさ $U_k$ を構成、「完成した設計手順の比較」と総費用。
- `2026-10-11` — 起票 (親計画 §5.1 #3・#4)。§4・§6 は上位に諮る前の下書き。
