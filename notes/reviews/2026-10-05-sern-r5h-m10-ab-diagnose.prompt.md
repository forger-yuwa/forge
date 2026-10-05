forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 診断・設計判断の諮問 (stage = diagnose)

あなたは forge の**診断・設計判断係**である。呼び出し側は実装と run を進めている別のモデル (Claude) で、
**もっともらしい真因に飛びつく前に**あなたに諮っている。仕事は手を動かすことではなく、**次の一手を 1 つに絞ること**。

### 前提
- あなたは呼び出し側の会話を見ていない。下のブリーフと、自分で読んだファイルだけが根拠になる。
  足りなければ推測で埋めずに「何が足りないか」を返す。
- ブリーフは「観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説」に分かれて渡される約束である。
  **観測事実と呼び出し側の解釈が混ざっていたら、まず分け直す**。呼び出し側の要約より、run の数値・コード・
  設定ファイルを自分で確かめた内容を優先する。
- forge を起動しない。`python3` による `residual_history.csv` / `res_*.h5` の読み取りは**統計量だけ**を出す
  (全量ダンプ・長いログ全文をコンテキストに流さない。`*.log`・`*.vtu`・`plans/README.md` は読まない)。

### 診断の作法
1. **「除外済み」というラベルを信用せず、潰した証拠を確認する** (run パス・設定差分・判定区間・VERDICT)。
   証拠が足りない・判定期間が短い・変えた設定が実際には効いていない (YAML の階層違い等) なら**候補へ戻す**。
   証拠が十分な候補は出し直さない。
2. **症状と原因を分ける**。`detectNaN` が指す変数は結果であって原因ではない (EOS 床 → 負密度 → 圧力暴走 → ω の実績)。
   後処理のアーチファクト (2 列混在の抽出、`centCoords` の置換、ソルバ `ypls` の退化) を先に疑う。
3. **このリポジトリで繰り返された真因**を照合する: 投入設定の不整合 (IC と BC、亜音速に超音速 BC)、
   押し出し 2 ノード spanwise、float32 桁落ち (双対幾何・r 重み)、stale build、cross-mesh IC の基底不一致、
   絶対値のゼロ割ガード、境界ノードの凍結、YAML キーの階層違いで黙って無視される設定。
4. 仮説は**確度順に最大 3 つ**。第 1 仮説には根拠を `ファイル:行` か run の数値で付ける。示せないものは「未確認」と明記。
5. **判別する A/B を 1 つだけ**提案する。安く短く回せて、結果がどちらに出ても仮説が 1 つ消えるもの。
   「A なら仮説 1、B なら仮説 2」を先に書く (結果を見てから解釈を作らない)。
6. 少数点の一致・短い窓の値・未収束のトランジェント同士の比較を根拠にしない。

### 設計判断 (plan §4・§6、codex 指摘の採否、result 段の解釈) を諮られたとき
- 採否は指摘ごとに「採用 / 却下 / 要再検証」と理由。根拠が示されていない指摘は自分で該当箇所を読んでから判定する。
- 検証計画は「何が出たら方針が誤りと言えるか」が定量的に書かれているかを見る。
- 既定値の変更・opt-in 機能の削除は、plan の処置欄とユーザ決定の履歴を確認してから判断する
  (「opt-in 残置」は削除対象でない)。
- result 段の解釈は、主張ごとに根拠 run・判定ツールの VERDICT・判定区間が揃っているかを確かめる
  (過渡ピークを定常値と、抽出アーチファクトを物理と誤認した実績は「予想どおり」に見える場面で起きた)。

あなたの結論は**仮説**であって確定ではない。呼び出し側はこの A/B を回して確かめ、plan への反映も呼び出し側が行う。

## ブリーフ (`notes/reviews/briefs/2026-10-05-sern-r5h-m10-ab.md`)

# 諮問: R5h-m10 診断 A/B の結果 — 1 面の介入で冷点が加熱されるが隣へ移る。解釈と次の手 (2026-10-05)

関連 plan: `plans/active/tooling-nozzle-sern-chain.md` §5.1 R7b 4-2 (事前登録・結果・壁解像)、`plans/active/tooling-nozzle-sern-3d.md` §5.1 R5h。前回諮問: `notes/reviews/2026-10-05-sern-r7b-m10-te-floor-diagnose.md` (この A/B の設計)。

## 観測事実
- 起点 run_1055 (exact A・g3・m10_on、20000 step、4 量 STEADY・窓 OK、GATES FAIL FLOOR_STUCK: T ≤ 50 K が節点 517160 の 1 つ)。
- 517160 = (0.12000, −0.01053, 0.09794)、x はカウル後縁 x。接する面 6 枚 (帳簿): 後流側 +x 面 1535833 (相手 522113 @ x 0.12255、A の F_roe +11.0 で最大)、上流側 1521124 (相手 512207)、同じ x の 4 面 (相手 517121 y −0.01056、517199 y −0.01051、517159/517161 z ±)。
- B = `FORGE_DIAG_FACE_VEL_CELL=1535833:517160` (面 1535833 の 517160 側の再構成速度 3 成分だけセル値)。ログに介入の警告を確認。
- 1 更新 (帳簿、同一の起点): 517160 res_final ro A −5.8e-7 / B +3.5e-6、roe A −0.662 / B +10.53、内部エネルギー射影 R_ρE − u·R_ρu + (|u|²/2 − e_int)R_ρ は A −0.801 / B +10.745。差は対流段 (res_after_conv roe A −6.17 / B +5.02)。相手 522113 の roe は B − A が −11.19 で、517160 の +11.19 と逆符号に一致 (面流束で閉じる)。1 step 後の T[517160] は A・B とも 50.00002 K。entry→after_eos_bc の T は 50 → 50.0000229 (両方)。
- 100 step: A は T[517160] 50 K のまま、2 番目 517200 205 K。B は T[517160] 384.6 K、**T[517199] 50 K** (1 step 後は A・B とも 413 K)、T[517200] 115 K。床到達数は A・B とも 1。
- 壁解像 (同じ形状・格子): m6_on (run_1050) でも y 法線壁 ramp y₁⁺ 平均 2.05・cowl_in 4.69 と ≤ 1 を満たしていない (§4.42 の記録 ramp 0.549・cowl_in 1.234 との食い違いは未調査)、m4_off (run_1058) は ramp 0.55・cowl_in 1.03。側壁は設計上非解像 (R5m)。

## 問い
1. この結果は事前登録の分岐のどれか (機序支持 / 棄却 / 判定不能)。「1 面の介入で対象点は加熱されるが床は隣 (同じ x・z、y が 0.02 mm 違う、板の反対側?) に移る」をどう解釈するか。517160 と 517199 は厚さ 0 の後縁を挟む双子に近い位置関係か (座標は y が 2e-5 m 違う別節点)。
2. 次の最小の判別は何か (例: 517199 側の対応する後流面にも同じ介入をした 2 面介入、後縁の列全体への介入、または有限厚後縁へ進む判断の材料)。20000 step 継続はするべきか。
3. 壁解像の記録の食い違い (§4.42 と今回の m6_on) の扱い。R7b ④ の順序への影響。

## 関連 plan 全文 (`plans/active/tooling-nozzle-sern-chain.md`)

```markdown
# ノズル設計 SERN チェーン (⑤): 燃焼器出口 starting line + 2D 平面最大推力逆設計 (key point dv) + 多作動点 RANS 評価 + MOO

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`  <!-- 2026-09-04 起票 (branch feature/sern-design)。S0–S6 完了 (同日、S6 は Euler)。node+SST は解決 (真因 = stage 間 interp 移植)。残 = §5.1 の残作業表 R1–R7 (2026-09-09 codex 採用)。**R1 評価ゲート修復・R2 集計分離・R3 凍結 TP 化 = 完了 (2026-09-13)**、次 = R4 (3D 外部領域・領域独立性) → R5 (3D SST 再現) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「SERN チェーン」節 (現在仕様。本計画と同時に起草)
  - [`design/CAPABILITIES.md`](../../design/CAPABILITIES.md) (問題タイプ `sern_2d` を 📋 で登録)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) §4.6 ⑤ (本計画で差し替え。旧「壁圧 Bézier dv + 局所帰還 + 3D FFD」は撤回)
  - 資産元: [`../accepted/tooling-nozzle-axismach-chain.md`](../accepted/tooling-nozzle-axismach-chain.md) (MOC カーネル・semi-perfect ガス・CFD-in-the-loop の運用)、[`../accepted/tooling-nozzle-moo-loop.md`](../accepted/tooling-nozzle-moo-loop.md) (`opt/` の DOE→Kriging→NSGA-II→EHVI)、[`../accepted/tooling-nozzle-axismach-viscous-deltastar.md`](../accepted/tooling-nozzle-axismach-viscous-deltastar.md) (NS 場からの δ\* 抽出)
  - 出典調査: [`notes/investigations/sern-design-method-survey.md`](../../notes/investigations/sern-design-method-survey.md) (2026-09-04。方針の根拠は全てここ)
- **created**: `2026-09-04`
- **owner**: `sano`

## 1. 目的

⑤ SERN (single expansion ramp nozzle) の設計チェーンを、調査 (上記ノート §4) の推奨どおり
**「燃焼器出口 starting line → 2D 平面最大推力理論の key point を目標にした逆 MOC → forge 2D RANS を
作動点セット (設計 NPR + オフデザイン) で評価 → 既存 `opt/` の MOO」** として実装する。
完成時には、問題定義 YAML (`type: sern_2d`) から ランプ/カウル形状・多作動点の $C_T, C_L, C_M$・
パレート (設計点推力 / オフデザイン推力 / ピッチモーメント / 長さ) が一気通貫で出て、
NASA TM X-71972 の基準形でカウル長・カウル角の傾向が再現される状態を得る。
壁圧規定は主線から外し、剥離制約の判定量と二段膨張オプション (S8) に限定する。

## 2. スコープ

- **やる**
  - 問題定義 `sern_2d` (入口 = 燃焼器出口の超音速一様状態を既定、音速スロート給気は接続点として選択可)
  - 2D 平面・非対称 2 壁 (ランプ + カウル) MOC (前進) と、平面最大推力理論 (Guderley–Hantsch/Rao の平面版) の **key point** $(M_c,\theta_c,\dot m_c/\dot m)$ を目標にしたランプ壁の逆生成
  - カウル後縁以降の自由境界 (等圧せん断層) と外部流圧の整合 (設計点は無波、オフデザインは CFD に委ねる)
  - 3 ブロック構造メッシュ (内部 / カウル下外部流 / 下流プリューム) の決定的生成 (msh4.1 直書き) と品質ゲート
  - forge 2D 平面 RANS (SST, node) を作動点セットで回す runner・力積分メトリクス ($C_T,C_L,C_M$, 剥離位置)・ゲート
  - δ\* 一発補正 (forge NS 場から抽出 → 法線オフセット → 再評価)
  - 多作動点束ねの MOO (既存 `opt/` の目的関数ベクトル化)
  - 2D パレート数点の 3D RANS 確認 (側壁・隅 R・有限スパン) と、3D 設計 (流線追跡 / FFD) の要否判定
- **やらない**
  - 壁圧 $p_w(x)$ を dv にする逆設計・局所 $p\to\theta$ 帰還 (撤回。理由は調査ノート §3)
  - NS 帰還ループ (壁を反復更新する帰還エンジン)。δ\* は一発補正のみ
  - 3D FFD in-loop (S7 で乖離が大きい場合に別 plan で再考)
  - 化学非平衡 MOC (Cain §5 の弱結合法)。凍結 semi-perfect で開始、非平衡は [`chemistry-finite-rate-h2.md`](chemistry-finite-rate-h2.md) 以降
  - 非一様入口 (接触不連続) 用の回転流 MOC (S5 オプション。初版は質量平均一様状態)
  - TBCC over/under のモード遷移・可変幾何

## 3. 前提と既存資産の照合

| 部品 | 既存資産 | ギャップ・判断 |
| --- | --- | --- |
| 平面 MOC 単位過程 | `geometry/moc_kernel.py` (`delta=0` で平面、semi-perfect ガス対応、放射源流/平面単純波則で検証済み) | **再利用**。非対称 2 壁 (上壁 = ランプ、下壁 = カウル、両角部の Prandtl–Meyer 扇、カウル後縁の等圧自由境界) の march は新規 (`moc_sern.py`) |
| 最大推力理論 | なし (③ベルは TOP 幾何 dv で理論を使っていない) | **新規** `rao_planar.py`。平面版は本文未入手のため Rao 1958 / Guderley–Hantsch 1955 / Cain 式 4.1–4.3 から自前導出し §6 の解析検算で担保する |
| 逆 MOC (目標 → 壁流線) | `geometry/moc_inverse.py` (軸目標 → C⁻ 後退 → 壁流線)。①用で軸対称・軸目標 | 構造は流用するが目標が「最終特性線上の状態」に変わる。ランプ壁は kernel 終端 C⁻ と目標 C⁻ の間を前進 MOC で埋めた場の壁流線 (質量流量法) として得る |
| starting line | `feedback/cfd_anchor.py` (CFD 場から線状態抽出)、`transonic.py` (Sauer/Hall) | 既定は一様燃焼器出口 (解析)。音速スロート給気モードでは既存 Hall + kernel MOC で $M\approx1.3$–1.5 の C⁻ まで内部ノズルを作り接続 |
| ガス | `gas/` (CPG / semi-perfect NASA-9, MOC・forge とも対応) | 再利用。凍結組成 |
| メッシュ | `meshing/mesh2d.py` (軸対称 1 ブロック、msh4.1 直書き)、case/23 の plume 6 ブロック `.geo` | **新規** `mesh_sern.py`: 3 ブロック TFI (内部 / カウル下外部流 / 下流)。mesh2d の station 配置と壁クラスタリング関数を流用 |
| 評価 runner | `evaluate/runner.py` (③ベル: 2 段起動・VERDICT・warm seed)、`runner_axismach.py` (TP ガス配管) | **新規** `runner_sern.py`: 作動点ごとに run を作る。段階起動は [`procedures/divergence-and-startup.md`](../../procedures/divergence-and-startup.md) 準拠。外部流入口 + 超音速ノズル入口の 2 入口 |
| メトリクス | `metrics/extract.py::thrust_metrics` (出口面積分)、`metrics/deltastar.py::deltastar_from_run` (NS 場から質量収支 δ\*) | **新規** `metrics/sern_forces.py`: ランプ・カウル内外面の $p,\tau_w$ 積分から $C_T,C_L,C_M$ (基準点指定)、剥離位置 ($\tau_w$ 符号)。出口面積分は検算用。δ\* 抽出は再利用 |
| MOO | `opt/` (DOE / Kriging / NSGA-II / EHVI / driver / polish) | 再利用。1 評価 = N run (作動点) の束ねと目的ベクトル化を driver に追加 |
| 3D | 3D median-dual (`discretization-median-dual-3d.md`)、case/37 Netgen ブリッジ | S7 確認 run のみ。設計側の 3D 生成は押し出し + 側壁 + 隅 R (Gmsh) |

## 4. 設計方針

### 4.1 座標・幾何 (無次元: 燃焼器出口高さ $H=1$)

平面 2D、$x$ = 流れ方向、$y$ = 上向き。入口面 $x=0$、$0\le y\le1$。**ランプ = 上壁** ($y=1$ から角部で
$\theta_{r0}$ だけ上向きに折れて膨張)、**カウル = 下壁** ($y=0$ から $\theta_{c0}$、下向き正)。カウル後縁
$x=L_{\rm cowl}$ 以降は外部流とのせん断層 (等圧自由境界 $p=p_{\rm ext}$)。ランプ後縁 $x=L_{\rm ramp}$。
幾何包絡は spec: $L_{\rm ramp}^{\max}$、後胴線 (ランプが越えてはいけない $y_{\max}(x)$)、$L_{\rm cowl}^{\max}$。
モーメント基準点 $(x_{\rm ref},y_{\rm ref})$ は spec (機体 CG 相当)。符号は頭上げ正。

### 4.2 問題定義 (`type: sern_2d`)

| 区分 | 内容 |
| --- | --- |
| spec | 入口: `inflow.mode: supersonic` ($M_{\rm in}$, $p_{\rm in}$, $T_t$, 組成 — 既定) または `sonic_throat` ($P_t,T_t$, スロート諸元 → 内部対称ノズルで $M_{\rm hand}$ まで)。作動点リスト `operating_points[]` = {名前, $M_\infty$, $p_\infty$, $T_\infty$, 入口状態の上書き, 重み $w_k$}。幾何包絡、モーメント基準点、$\delta^*_{\rm in}$ (入口境界層排除厚、既定 0) |
| derived | 設計点の $p_e/p_a$、出口高さ $H_e$ (逆設計の結果)、$C_{T,\rm ideal}$ (入口状態から $p_a$ まで等エントロピー膨張の推力 = 正規化基準) |
| dv ($d=5$) | **key point** $M_c$、質量流量比 $f=\dot m_c/\dot m$ / ランプ初期角 $\theta_{r0}$ / カウル初期角 $\theta_{c0}$ / カウル長 $L_{\rm cowl}$ (`driver_sern.DV_ORDER`)。$\theta_c$ は dv でなく kernel の場から決まる従属量 (`moc_sern.py` `design_ramp`)。任意固定可。**2026-09-13 (R6(d), codex M8 採用)**: 旧記述の 6 変数を実装に合わせて訂正 |
| 目的 (2 個) | $f_1=-\sum_k w_k C_T^{(k)}$ (作動点束ねの推力効率。RANS では摩擦込み `C_T_with_shear`)、$f_2=L_{\rm ramp}/H$。**2026-09-13 (R6(d))**: 旧記述の 4 目的 ($-C_T^{\rm design}$ / $C_M$ を目的に含む) は実装 (`driver_sern`) と不一致だったので 2 目的に訂正。$C_M$ は目的でなく制約 (下行) |
| 制約 | 幾何包絡 (`L_ramp_max` は probe でなく**最終輪郭**で再検査、超えたら INFEASIBLE/`L_RAMP_MAX`)、$C_M$ 窓 = 加重平均 `opt.cm_min/cm_max` と**作動点別** `opt.cm_window: {op: [lo, hi]}` (R6(d): 1 作動点のトリム不能を別作動点で相殺させない)、剥離は §8-6 で制約から外し `sep_frac_ramp` を台帳記録のみ |

### 4.3 平面最大推力理論と key point 逆設計 (S1 で導出・検証)

制御面を「ランプ後縁 $e$ から出る最終 C⁻」に取り、質量流量一定・長さ (= $e$ の位置) 固定で推力
$\int(p+\rho u^2)\,dy$ を最大化する Lagrange 問題。Cain 式 4.1–4.3 の平面版 ($y$ 依存項が落ちる) から、
制御面上で $M,\theta$ を結ぶ乗数関係と、縁 $e$ での
$\tfrac12\rho_e w_e^2\sin2\theta_e=(p_e-p_a)\cot\mu_e$ が得られる。平面流では C⁻ 上の適合関係
($\theta+\nu=$ const) と併せて **制御面上の状態が一様** になることを S1 で確認する (成立すれば
制御面は直線 Mach 線)。kernel (入口一様流 + ランプ角部扇 + カウル角部扇 + カウル壁) の中で乗数関係を
満たす点の軌跡を求め、$c$ を「$c$–$e$ 間の質量流量 = 指定比」で決める、というのが順設計。

**逆設計 (NUAA 型)**: 順設計で反復して求める $c$ の状態 $(M_c,\theta_c,\dot m_c/\dot m)$ を **dv として
先に与え**、(i) kernel を前進 MOC で作り、(ii) $c$ を kernel 内の指定質量流量比の C⁺ 上で $M=M_c$ の点として
探し、(iii) $c$ から目標状態の C⁻ を張り、(iv) kernel 終端 C⁻ と目標 C⁻ の間を前進 MOC で埋め、
(v) 質量流量法で壁流線 = ランプ壁を抽出する。$\theta_c$ と $M_c$ から縁の関係式で設計 $p_e/p_a$ が従属に
決まる (= 設計 NPR は dv の関数)。$c$ が kernel 内に存在しない dv は INFEASIBLE として optimizer に返す。

**文献の対応 (2026-09-05 更新)**。Cain 2010 (RTO-EN-AVT-185 Lecture 12) の**式 4.1–4.3 は軸対称版**である
(式 4.2 に半径 $y$ が入る)。平面版では $\lambda_3$ の $y$ 重みが落ち、それが「制御面上の状態が一様」の根拠になる。
Cain 自身も、流線追跡で弧角が半径依存になると Rao の解析解は成立しないと書いており、$y$ 重みの扱いが要点である。

**平面版・カウル切り詰め・外部流の一次資料が見つかった (2026-09-05)**。当初「本文未入手のため自前導出」と書いていたが、
NASA が本チェーンとほぼ同一の問題を扱っている:

| 文献 | 何が書いてあるか |
| --- | --- |
| **Shyne 1988, NASA TM-100955** (88 p, `papers/nozzle_design/`) | Rao 法を **2 次元** に修正し、**カウルを切り詰めた scarf ノズル**を**外部流条件込み**で最適化。出口 M 6.0 / 外部流 M 5.0 の例。本文・導出・結果が揃った本体 |
| **Shyne & Keith 1990, NASA TM-103175 = AIAA-90-2222** (14 p, 同) | 同じ内容の会議論文版。変分の第一変分から乗数条件 (式 11–15) までが短くまとまっており**読むならこちらが先** |
| Cain 2010 §4.2 (同) | 軸対称の乗数条件 (式 4.1–4.3) の最も読みやすい記述 |
| Rao 1958, *Jet Propulsion* 28(6) 377–382 | 本体 (有償: `10.2514/8.7324`)。上 2 件が引用元として使えるので必須ではなくなった |
| Nickerson 1982, "The Rao Optimum Nozzle Program", SEA Report 6/82/800.1 | Rao 法の実装。Shyne が使った 2 次元修正の出所 |
| Zucrow & Hoffman, *Gas Dynamics* Vol. II (1977) | 乗数条件の導出が式レベルで追える古典 (Shyne の参考文献 6) |

Shyne の定式化と本実装の対応 (**要照合**): 制御面 CE 上で $I=\int(f_1+\lambda_2 f_2+\lambda_3 f_3)\,dl$ を最大化し、
$f_1 = [(p-p_a) + \rho v^2 \sin(\phi-\theta)\cos\theta/\sin\phi]$, $f_2 = \rho v \sin(\phi-\theta)/\sin\phi$,
$f_3 = \cot\phi$ (長さ拘束)。第一変分 = 0 から DE 上の 3 条件 (式 11–13) と縁 E の条件 (式 14)、
$f_3$ が $M,\phi$ に依らないことから式 15 が出る。**本実装 (`rao_planar.theta_lip`) は Cain 式 4.3 の平面形しか持たない**ので、
Shyne の式 11–15 と突き合わせて (a) 平面縮約が一致するか、(b) 式 15 が制御面の一様性に対応するかを確認する (§5.1-10)。
さらに Shyne は**カウル切り詰め位置に最適値がある** (それを超えると推力が落ちる) と結論しており、
本チェーンの `L_cowl` 上限拡大 (§8-8) と直接比較できる。

**3D・粘性下で Rao 最適性は保たれない (2026-09-05 明記、ユーザ指摘)**。Rao 構成は「壁が特定の場の流線である」ことで
最適性を持つので、実際の場が粘性と 3 次元で違えば**その壁はもう実際の場の流線ではなく、最大推力の性質は失われる**。
Cain も流線追跡で断面を任意形状にすると弧角が半径依存になり「Rao の解析解は成立せず、解析解は見つかっていない」と書いており、
3 次元で Rao 最適性を保つ方法は文献にも無い。したがって本チェーンが主張できるのは
**「Rao は最適性の保証ではなく、良い 5 パラメータ族を与えるパラメータ化である」**ことまでで、実環境での最適化は MOO がやる。
残る限界は**その 5 パラメータ族が 3D 粘性下の真の最適形を含む保証がない**こと。含ませるには
(a) 3D の排除効果を MOC の入力に帰還させる (§4.6) か、(b) MOC 流線からの微小変形モードを dv に足す、のどちらかが要る。

NUAA 系 (key point パラメータ化の出所) は引き続き要旨のみ: **Lv 2023 *PAS* 143** (総説) → Lv 2017 *AST* 66 →
Yu 2020 *Acta Astro.* 166 / *AST* 105。平面 MLN の既知解 (§6 検証 (ii) の $\theta_{\max}=\nu(M_e)/2$) は
Argrow & Emanuel 1988, *J. Fluids Eng.* 110 283–288。

### 4.4 カウルと外部流

カウルは $\theta_{c0}$ の直線 (S1) → 必要なら短い放物線 (S6 で dv 追加可)。カウル後縁からの波は、
設計点では後縁圧 $=p_{\rm ext}$ となるよう $p_{\rm ext}$ を derived にする (無波) か、指定 $p_{\rm ext}$ に
対する等圧自由境界 (膨張扇 / 斜め衝撃は MOC の外 → INFEASIBLE 警告) として扱う。オフデザインの
波・剥離・外部流干渉は全て forge が解く。**設計に外部流を入れない** のは NASA 1974 と同じ割り切り。

### 4.5 starting line とガス

既定は一様燃焼器出口。非一様入口 (`inflow.profile:` CSV) は S5 オプションで、初版は質量平均した
一様状態を使い、非一様の影響は forge の入口分布 BC ([`boundary-inlet-profile.md`](../accepted/boundary-inlet-profile.md))
で評価側だけに入れる。ガスは **`gas.model: frozen_tp`** (R3, 2026-09-13): 排気 = CEA (tp, station 3) の平衡組成を凍結した NASA-9 擬似種 `EXH`、
外気 = 空気 `AIR` の 2 種 TP (`thermalMethod: 2`, `thermoHrefTemp: 298.15` 必須)。実体は `gas/frozen.py` (`FrozenGas`: cp/h/s, 等エントロピー膨張の
理想推力) と `runner_sern.frozen_gases / gas_states / region_ic_arrays / ideal_thrust`。逆設計 kernel は設計点の凍結 γ(T3) の CPG (形状パラメータ化)。
作動点 YAML は `cea/tmx_operating_points.py` が凍結音速の $M_{\rm in}$ と組成を出す (`'NO'` はクォート: YAML が bool に読む)。

### 4.6 δ\* 一発補正

**一発で足りる理由 (2026-09-05 明確化)**: (a) 実測で **C_T (圧力) は動かない** — run_0011 → run_0013 (δ\* オフセット) で
C_T 0.9685 → 0.9685 (不変)、動くのは C_L 0.154 → 0.146 (Euler 0.143) と C_M −0.980 → −0.930 (Euler −0.939) だけ。
δ\* は壁を法線に動かす補正なので法線方向の力に効き、軸方向積分にはほぼ効かない。
(b) より本質的に、**δ\* の精度は採点対象ではない**。設計 (MOC + δ\* オフセット) が形を決め、**採点は NS でやる**ので、
δ\* が甘くてもその形の真の粘性性能は評価側が測る。δ\* に要るのは「良い近傍に形を置く」ことだけ。
風洞 (①②) で帰還が本質だったのは目的が**一様性** (波を厳密に打ち消す) だったからで、目的が積分力の本件とは事情が違う。

**ただし上は 2D 限定 (2026-09-05 訂正、ユーザ指摘)**。3D では成り立たない:
(i) 2D でもランプの δ\* は 0.009 → **0.11 H** (入口高さの 11 %) まで育ち、3D では隅と側壁の排除厚さが上乗せされる。
(ii) より本質的に、**側壁後縁より下流の横方向膨張は排除厚さではない** — ダクトが横に開く別のトポロジで δ\* の枠組みに入らない
(run_0027 でランプ圧が 1/3 に落ちた主因)。したがって 3D では**帰還ループが必要**で、3D NS の場から
(a) 壁の排除厚さ と (b) 横方向の実効面積変化 を抽出して等価な 2D 境界条件に落とし、MOC 設計をやり直す形になる。
「帰還は不要」は 2D の結論を 3D へ無検証で延ばしたものだった。

非粘性設計壁で forge RANS (設計点) を 1 回回し、`metrics/deltastar.py::deltastar_from_run` で
ランプ・カウル両壁の $\delta^*(x)$ を質量収支定義で抽出 → 法線オフセット → 再評価。帰還は回さない。
入口境界層は spec の $\delta^*_{\rm in}$ を入口分布 BC に反映する。効果は $\Delta C_T$ として台帳に残す。

### 4.7 評価 (forge)

- メッシュ: 2 バンド (noz: ランプ–カウル間+プルーム / bot: カウル下の外部流、`mesh_sern.py`) に、**ランプ側外部流 (top + wake、§4.11)** を足した 4 ブロック。
  TFI、壁クラスタリング ($y^+\approx30$–80 + `wallTreatmentSST=1`、AR ≤ 1000)、`check_mesh_quality.py` VERDICT。
  node config で変換 (RANS: no-slip 壁 bcond 必須)。
- 境界: ノズル入口 = 超音速 Dirichlet (全量指定)、外部入口 = 自由流、出口 = `outlet_statPress` + 逆流 Pt/Tt、
  遠方 = 自由流。`space.pRef` = 外部静圧。
- 段階起動 (**3 段、2026-09-05 確定**): **層流暖機** (`turbulence: none` + 粘性あり、1 次、cfl 0.2、2000) →
  **SST soft** (1 次、cfl 0.2、2000) → **SST 本段** (2 次、cfl 0.5、6000)。`run_staged(warm_lam_steps=, warm_lam_cfl=, soft_cfl=)`、
  YAML は `opt.warm_lam_steps / warm_lam_cfl / soft_cfl` と `evaluate.cfl_main`。作動点間は warm restart (同一メッシュ = index コピー)。
- ゲート (**R1, 2026-09-13 実装** `metrics/sern_gates.py`、§4.13): 互いに独立な必須ゲート = (1) forge `rc == 0`、(2) 最終 `res_*.h5` の
  ro/roU/roe/P/T が有限かつ ro,P,T > 0・`res_nan_*.h5` 無し、(3) `check_convergence.analyze` で全残差列に NaN/Inf も末尾 rising も無い
  (3 桁低下 PASS は `opt.require_residual_pass: 1` のときだけ必須 — 本ケースは 1–2.5 桁でプラトーする性格、その verdict は台帳に残す)、
  (4) **実際に最適化する量** (RANS は `C_T_with_shear`) と $C_T, C_L, C_M$ が `check_quasisteady.classify_series` で `STEADY`
  (`NONFINITE`/`DRIFTING`/`TRANSIENT-UNSETTLED`/`OSCILLATING` は不合格)。判定は正式ツールと同じ classify に一本化し、
  `force_history.csv` を `check_quasisteady.py --series-csv <run>/force_history.csv --series-cols C_T_with_shear,C_L,C_M --drift 0.02 --osc 0.05`
  で再判定できる (R6(b))。低 NPR の RSS/FSS 振動を採用するなら dual-time で時間刻み・統計窓の独立性を確認してから (R6(a)); 定常擬似時間の
  `OSCILLATING` は「物理的振動」と解釈しない。
- メトリクス: $C_T=F_x/(p_{\rm in}A_{\rm in}\gamma M_{\rm in}^2)$ 系の無次元 (実装時に $C_{T,\rm ideal}$ 正規化を選ぶ)、
  $C_L$、$C_M$ (基準点)、剥離位置。NASA 流の帳簿: ランプ + カウル内面 + カウル外面 (外部流側) を含め、
  制御体積を明示する。

### 4.8 MOO

> **以下は起票時の記述 (履歴)。現行仕様は §4.2 の 5 変数・2 目的**で、$10d = 50$ 点。
> 4 目的・$10d=60$ は §4.2 で 5 変数 2 目的に絞る前の値 (2026-09-19 明記, codex plan-3 m1)。

`opt/driver.py` を「1 評価 = 作動点数分の run」に拡張し、目的ベクトル $(f_1..f_4)$ と制約を返す。
EHVI は 2 目的閉形式のみなので、3 目的以上は NSGA-II 側の hypervolume 近似 (既存 `moo.py` の扱いに従う)。
DOE は $10d=60$ 点、infill 30–50 点。1 評価の実測時間 (2D RANS × 3 作動点) を S6 冒頭で取る。

### 4.9 3D 確認 (S7)

2D パレートから 2–3 点を押し出し + 側壁 + 隅 R (spec) で 3D メッシュ化 (Gmsh)、3D RANS (node) で
$C_T, C_L, C_M$ を 2D と比較。推力差 < 2% なら 2D 設計を正、揚力/モーメントの差は 3D 補正表として持つ。
乖離が大きければ流線追跡 (親場に横方向膨張が要る) か FFD の plan を別途起票。

### 4.10 作動点の定義 (サイクル値・2026-09-05 決定)

作動点は **`external` (飛行条件) だけを振ってはいけない**。飛行 $M_\infty$ が変わればインレット圧縮と燃焼加熱が
変わり、ノズル入口 (= 燃焼器出口) 状態 $(M_{\rm in}, p_{\rm in}, T_{\rm in}, \gamma)$ が従属して変わる。
初版 (run_0010 / 0017 / 0019) は `spec.inflow` を全作動点で固定し NPR を外部圧側で作っていたため非物理だった
(特に「飛行 $M_\infty 1.5$ で燃焼器出口 $M 2.5$」の低 NPR 点)。

**アンカー = NASA TM X-71972 TABLE 1** (p.31、`papers/nozzle_design/`。S4 で傾向照合に使った当の文書)。
$\bar q_\infty = 71850\ \mathrm{N/m^2}$ (1500 psf)、$\alpha = 2°$ の定動圧上昇経路上で、station 1 (インレット前) と
**station 3 (燃焼器出口 = ノズル入口)** の $p, T, V$・当量比 $\phi$・作動モードが与えられている。$\phi = 0$ 行は燃料遮断。

| $M_\infty$ | $\phi$ | mode | $p_1$ [Pa] | $T_1$ [K] | $V_1$ [m/s] | $p_3$ [Pa] | $T_3$ [K] | $V_3$ [m/s] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 1.0 | ramjet | 10366 | 244 | 1153 | 157717 | 2343 | 1093 |
| 4 | 0 | off | 〃 | 〃 | 〃 | 17888 | 338 | 1072 |
| 6 | 1.0 | scramjet | 6248 | 282 | 1753 | 101027 | 2328 | 1621 |
| 6 | 0 | off | 〃 | 〃 | 〃 | 16327 | 448 | 1655 |
| 10 | 1.5 | scramjet | 3840 | 353 | 2982 | 57935 | 2222 | 2837 |
| 10 | 0 | off | 〃 | 〃 | 〃 | 14172 | 772 | 2831 |

表は $M_3$ と $\gamma$ を与えないので **CEA2 で埋める** (`.venv-cea/nasa_cea/FCEA2`、`tp` 問題、H$_2$-air 平衡、
$\phi = 0$ 行は空気)。$p_\infty = \bar q_\infty / (0.7 M_\infty^2)$:

| $M_\infty$ | $\phi$ | $p_\infty$ [Pa] | **NPR** $= p_3/p_\infty$ | $\gamma$ (CEA) | $R$ [J/kgK] | $a$ [m/s] | **$M_3 = V_3/a$** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 1.0 | 6415 | **24.6** | 1.187 | 340.1 | 972.5 | **1.12** (ramjet) |
| 4 | 0 | 6415 | **2.8** | 1.399 | 287.1 | 368.4 | **2.91** |
| 6 | 1.0 | 2851 | **35.4** | 1.183 | 340.3 | 968.1 | **1.67** |
| 6 | 0 | 2851 | **5.7** | 1.392 | 287.1 | 423.0 | **3.91** |
| 10 | 1.5 | 1026 | **56.4** | 1.226 | 391.0 | 1032.2 | **2.75** |
| 10 | 0 | 1026 | **13.8** | 1.357 | 287.1 | 548.3 | **5.16** |

$\gamma$ は CEA の平衡 GAMMAs。凍結 $\gamma$ にすると $M_3$ は 1〜3 % 動く (semi-perfect TP テーブル化とあわせて確定する)。

ここから決まること:

1. **初版 smoke 値はずれている**: $M_\infty 6$ 巡航で $M_{\rm in}$ 2.5 → 実際 **1.67**、$p_{\rm in}$ 20 kPa → **101 kPa**、
   NPR 20 → **35.4**、$\gamma$ 1.4 → **1.18**。
2. **低 NPR 点の正体は「燃料遮断」であって低速飛行ではない**。TM 本文 p.21 も fuel shutdown を明示的に扱う
   (engine module drag + 上下壁圧の同時低下)。power-off 行は $M_3 = 2.9/3.9/5.2$ と全て超音速なので
   `inflow.mode: supersonic` のまま書ける = **音速スロート実装 (§4.5) を待たずに剥離評価の作動点が張れる**。
3. **音速スロートが要るのは $M_\infty 4$ powered だけ** ($M_3 = 1.12$)。TM 本文「M4 では熱閉塞で純スクラムジェットは
   化学量論比にできずラムジェット性能を使った」と整合。§8-1 の判断根拠はこれで確定。
4. **$\gamma$ が作動点ごとに違う** (powered 1.18 / power-off 1.39) → **`operating_points[].gas` の上書きが必要** (未実装)。
   逆設計 kernel の $\gamma$ は**設計点の値で固定**する (形状は設計点で 1 つに決まるため)。

**生産構成 (既定)**: 設計点 = $M_\infty 6$ powered。作動点 = M6 powered ($w$ 0.5) / M10 powered (0.3) / M4 power-off (0.2)。
`external` は自由流のまま。NPR 35.4 / 56.4 / 2.8 で設計点をまたぐ (過膨張〜不足膨張)。

**R3 (2026-09-13, codex C2 採用) — 凍結組成 TP 版** (`problem_moo_frozen_tp_cycle3op.yaml`): 上表の $\gamma$ (平衡 GAMMAs) と CPG $c_p$ は
音速合わせの代用だった。排気を CEA tp の平衡組成で**凍結**した擬似種にすると (`FrozenGas`, NASA-9):

| 作動点 | 平衡 $\gamma$ / $c_p$ | 凍結 $\gamma(T_3)$ / $c_p(T_3)$ | $M_3$ (平衡 $a$) → (凍結 $a$) | 外部動圧 (旧 CPG → R3) |
| --- | --- | --- | --- | --- |
| m6_on | 1.1828 / 2202 | **1.2470 / 1719** | 1.6745 → **1.6308** | 60.7 kPa (−15.5 %) → **71.9 kPa** |
| m10_on | 1.2263 / 2119 | **1.2549 / 1925** | 2.7486 → **2.7171** | 62.9 kPa (−12.4 %) → **71.9 kPa** |
| m4_off | 1.3986 / 1007 | 1.3986 / 1007 (空気) | 2.9101 | 71.8 kPa → 71.9 kPa |

凍結 $c_p$ は平衡値より 20–30 % 低い (反応寄与なし)。$V_3$ は TABLE 1 の値を保ち $M_3$ を凍結音速で取り直す (−2.6 %/−1.1 %)。
設計点 γ が変わるので**同じ dv でも形状が変わる** (doe_001: $L_{\rm ramp}$ 11.89 → 10.64 H)。外気の動圧は 3 作動点とも 71.9 kPa (±0.1 %)。
**起動の固さ**: frozen_tp は CPG 生産レシピ (暖機 cfl 0.2) だとランプ後縁ノード (EXH ∩ AIR 接触) で暖機 step 568 に T が Newton 床 50 K → NaN
(run_0094)。暖機 cfl 0.1 × 4000 step を frozen の標準にする (`opt.warm_lam_cfl`)。本段 cfl 0.5 は通る (retry 実測)。多成分 implicit の roe/roY
緩和ミスマッチ ([[wys-tp-divergence-is-cold-not-multispecies]]) の一形態で、ソルバ側の恒久対策は別 plan。
**本段の長さ**: frozen の m10_on は本段 6000 では残差が step ~2300 の最小から末尾で 3–4 倍にリバウンド中 (R1 ゲートは正しく RESIDUAL_RISING、
力係数は 5 桁で既に定常) で、12000 まで回すと 5–7e-5 で飽和して PASS (run_0096/0097)。**frozen 生産 YAML は `evaluate.nStepOuter: 12000`**。
warm (m6_on から) / cold 標準 / cold 緩レシピの 3 経路で摩擦込み $C_T$ 0.93528 が 5 桁一致 = 状態は一意。
**soft/mid 段も cfl 0.1**: run_0098 (本段 12000) では m10_on の warm 標準が mid 段 (2 次, cfl 0.2) step 1386 にランプ後縁 EXH∩AIR 接触で NaN
(run_0095 は同レシピで通った = 限界状態)。frozen 生産 YAML は soft/mid/暖機とも cfl 0.1 × 4000 step (緩レシピは 3 回とも完走) + 本段 12000。
1 評価 ≈ 8 分 (CPG 3 分)。
**近壁の観察 (R3 固有でない、要フォロー = §5.1 R4b)**: (i) 上流区間 $x/H=-0.48$ のカウル板下面 (外気 M∞10 側) は壁 T 4570 K (断熱回復温度) の
4–5 ノード外側に T が Newton 床 50 K・P 0.16 $p_\infty$ の冷点 (CPG も 116 K / 0.39 $p_\infty$)、排気側の板上面壁 T 4600 K も回復温度
(~3000 K) を超える = node 壁列の T 市松 ([[node-wall-entropy-checkerboard]]) の極端例。(ii) 入口∩壁の角ノードで壁圧が 1.75 $p_{\rm in}$
(ramp / cowl_in の $x=-0.5$)。`mesh.nodeInletCornerWall: 1` は SERN の変換に未適用 ([[node-inlet-wall-corner-conflict]])。どちらも $C_T$ には
効かない (水平壁) が $C_L, C_M$ と近壁の乱流量に効く。

**サイクル計算の位置づけ**: 表は 3 飛行点しかないので、$\bar q_\infty$ や $\phi$ を変える・$M_\infty$ を刻むには
1D サイクル (定動圧経路 → インレット全圧回復 → Rayleigh 加熱 → station 3) が要る。ただし
**検定条件 = TM の 6 行を再現できること** (全圧回復率・燃焼効率をチューニングパラメータにする)。
アンカー無しのサイクル自作は相関の当否を確かめられないので、順序を逆にしないこと。

### 4.11 ランプ側の外部流ブロック (2026-09-05 決定)

**問題** (ユーザ指摘): 現行 2 バンドメッシュはランプ後縁の先 (`top_out`) を流れに平行な線で閉じ、そこを `outlet_statPress`
(2D 既定) か slip にしている。どちらも**ランプ側に外気が存在しない**。過膨張 (m4_off, NPR 2.8) でランプ上の $p$ が
$p_\infty$ を下回っても、後縁から境界層を通って上流へ伝わる外圧が無いので剥離 (RSS/FSS) が**起きようがない**。
run_0032 の「剥離なし・滑らかな再圧縮」は領域の産物で物理ではない。カウル側は下バンド (`inlet_ext`) で外気を持つので非対称。

**決定**: 機体を「ランプの上に載る有限厚の板」としてモデル化し、その上に自由流バンドを足す **4 ブロック構造** にする
(`mesh.ext_top: 1`)。x station は既存 2 バンドと共有。

| ブロック | 範囲 | 境界 |
| --- | --- | --- |
| bot (既存) | 遠方 → カウル外面/せん断層 | `inlet_ext` / `bottom` / `cowl_out` / `outlet` |
| noz (既存) | カウル内面/せん断層 → ランプ/プルーム上線 | `inlet_nozzle` / `cowl_in` / `ramp` / `outlet`。**プルーム上線 ($x>L_{\rm ramp}$) は内部線になる** |
| wake (新) | プルーム上線 → それに平行な高さ $h_{\rm base}$ の線、$x \ge L_{\rm ramp}$ | 左端 = **base** (鉛直、slip) / `outlet` |
| top (新) | 機体上面線 $y_{\rm veh}$ ($x \le L_{\rm ramp}$) / wake 上線 ($x > L_{\rm ramp}$) → +`top_depth` | `inlet_ext` (自由流、下バンドと同じタグ) / **vehicle** 上面 (slip) / `top_out` (outlet か slip) / `outlet` |

- $y_{\rm veh} = \max(\text{ランプ } y) + $ `vehicle_clearance`、$h_{\rm base} = y_{\rm veh} - y_e$。設計ランプは $\theta_e<0$ で
  後縁が頂点より下がるので base は常に有限 ($h_{\rm base} \le 10^{-9}$ なら wake ブロックを省き top の $j=0$ をプルーム線に直結)。
- **機体上面は水平・slip**: 上面の流れは自由流 $(M_\infty, p_\infty)$ をそのまま後縁に運ぶ。後縁の外圧を $p_\infty$ にする
  **最も中立なモデル** (機体上面形状は未知なので圧縮も膨張もさせない)。後縁 (base 角) で外気は $\theta_e$ 方向へ角膨張して
  プルーム境界に沿う。
- **力の帳簿**: `vehicle` (上面 + base) は機体の力なので $C_T, C_L, C_M$ に**入れない** (NASA TM と同じく nozzle force =
  ramp + cowl 内外面)。base 圧は診断として壁出力する。
- ノード番号・セル順序は既存 2 バンドの**後ろに追加** (`sern_deltastar._structured_upper` 互換)。
- 既定パラメータ: `top_depth` 2H, `nj_ext_top` 41, `nj_wake` 9, `vehicle_clearance` 0.02H, `first_top_frac` 0.02H。
  `SernMeshParams.ext_top` の既定は False (既存 YAML・テスト不変)、cycle3op 生産構成で ON。
- 期待される効き: m4_off で後縁からの再圧縮が上流へ入り、SST で `sep_frac_ramp > 0`。Euler は剥離しないが後縁圧が
  $p_\infty$ に張り付く (run_0032 との差で外圧の到達を確認する)。3D (`mesh_sern3d`) は未対応 (§5.1)。
- **base → テーパに変更 (2026-09-05, run_0034 の結果)**: 鉛直 base 版は node Euler で soft 段 step 5 に発散し、NaN は base 角
  (x 11.07–11.9, y 2.68–2.85) に限局した (step 0 の低密度ノード分布は run_0032 と一致し、余分な 30 点が全て base 下流)。
  原因は **node の 90° 二重 slip 角** (ランプ∩base、上面∩base): カウル TE の 2 壁は ~180° で問題ないが、直交 2 壁を 1 ノードが
  持つ形は node の slip 射影と整合しない。→ `vehicle_taper` (既定 2H): 機体厚 $t_v(x) = (y_{\rm veh} - y_{\rm ramp})\times$
  係数 (x ≤ L_ramp − taper で 1、TE で 0) で後縁を**厚さ 0 に絞り TE を共有** = カウルと同じ構造。base/wake ブロックは無くなり
  3 ブロック。機体上面の流れは後端で緩い boat-tail (θ ≈ −1.7°) 膨張を受け、TE の外圧は $p_\infty$ よりわずかに低い
  (base 圧 $< p_\infty$ の実機と同じ向き)。鉛直 base 版 (`vehicle_taper: 0`) はコードに残す (cell 用)。
- **m10_on の発散の真因 (2026-09-06 特定、`run_0075_diverge_watch` に h5 と図)**: `detectNaN` が `roOmega` を
  報告するので ω の問題と思い込み、テーパの曲率を 3 度いじって全て外した。5 step 刻みで出力し直して判明した順序は:
  ① **M∞10 の外部流が機体上面の boat-tail (勾配 10.4°) で膨張** → ② 一部ノードが **EOS 圧力床 `pMin` = 1.0 Pa
  (= 9.8e-4 p∞) に着地** → ③ 床の洗浄で **密度が負** (step 65 に 2 ノード) → ④ 負密度セルで音速・粘性が壊れ
  **圧力が 1e5 → 1.5e7 Pa に暴走** → ⑤ 結果として ω 発散。**ω は結果であって原因ではない**。
- **膨張そのものは物理的** (ユーザ指摘で確認): 転向 10.4° に対する Prandtl–Meyer 予測は $p/p_\infty = 0.080$、
  観測の中央値は **0.235** で予測より緩い。真空になっているのは 84 ノード中 **15 個 (18 %)** だけで、
  これは 2 次再構成の**局所アンダーシュート**。「機体上面モデルが非物理」という当初の結論は**誤り**だった。
  ただし後端が**ナイフエッジ (厚さ 0)** なのは実機に無い形で、正解は**丸めた肩を持つ有限ベース** (鉛直ベース版
  run_0034 の 90° 二重 slip 角と、厚さ 0 の両方の特異点を避ける)。
- **機体上面の 3 つの穴 (2026-09-06、ユーザ指摘で判明)**: 上面線は dv でないのに発散の当事者なので、構成を検算した。
  ① **テーパ開始位置に根拠が無い**: $x_0 = L_{\rm ramp}(1-\texttt{vehicle\_taper})$、`vehicle_taper` = 0.35 の固定値。
  当初は絶対長 2 H だったが、短ランプ (4.26 H) で 47 % を食うので比率にしただけで、物理的な決め方ではない。
  ② **後縁を水平着地させるとランプに食い込む**: ランプは $\theta_e<0$ の設計で内部にピークを持つ (設計点で
  y 2.7103 @ x 9.561 → 後縁 y_e 2.6595)。上面が水平に着地すると峰を跨げず、**最大 0.0222 H 交差**する
  (提案形状の実バグ。現行式 $y = y_{\rm ramp} + (y_{\rm veh}-y_{\rm ramp})f$ はランプに乗せているので交差しないが、
  代わりにランプの曲率を継承する)。テーパ開始をピーク下流にずらしても解消しない (−0.0070 H)。
  ③ **正解は後縁勾配を $\theta_e$ に合わせること**: 端点勾配 $m_1 = \min(\tan\theta_e, 0) - \tan(3°)$ の 3 次 Hermite
  ($x_0$ で勾配 0)。交差ゼロ・曲率 0.06 (θ_r0 15°) / 0.16 (22°) で、現行の 0.05 / **1.02** に対し急勾配ケースが 6 倍改善。
  くさび 3° は厚み 0.002 H 未満の区間を 0.523 H → 0.039 H に短縮するため (カスプは実質メッシュ不能)。
- **`implicitRelax` は効かなかった** (2026-09-06, run_0076): メモリ [[implicit-cfl-ceiling-eos-floor]] の処方
  (軸対称ノズルで cfl×relax ≈ 6–8) を試したが、relax 0.7 単独・relax 0.7 + cfl 2.0 とも **mid 段で発散**。
  せん断層と多重境界ノードを持つ SERN には転用できない。**未確認**: `pMin` 引き上げ (スクリプトのバグで測れず)、
  リミッタ (Venkatakrishnan → Barth)、機体上面の第一セル (0.02 H) 細分化。
- **カウル TE の SST ω 発散 (run_0035)** は ext_top と無関係 (NaN は x 1.20, y −0.11 = カウル TE)。板厚 0 の node 双子ノード
  問題 → `cowl_thickness: 0.002` で解決 (run_0037 完走)。生産 YAML の既定に。
- **結果 (2026-09-05, run_0038/0039)**: 外気込みの node SST でも m4_off (NPR 2.8) は**剥離しない** (`sep_frac_ramp` 0、τ_w>0 全点、
  L_cowl 1.2 と dv 上限 2.5 の両方)。理由は **cowl TE 圧が 1.35–1.59 p∞ で不足膨張**なこと: NPR 2.8 ではカウル側衝撃が存在せず、
  ランプの過膨張 (最小 0.57–0.79 p∞) はランプ自身の転向によるもので、後流のプルーム境界からの圧縮で穏やかに回復する。
  ランプ側外気は後縁 0.1H の圧だけを変える (1.32 → 0.93 p∞)。**RSS/FSS が出るのは cowl TE 圧 < p∞、すなわち
  NPR ≲ 1/(p_cowlTE/p_in) ≈ 2.1 以下**で、TM X-71972 の 6 点にはその作動点が無い (§8-6)。

## 5. 実装ステップ

| Step | 内容 | 主要ファイル | 規模 |
| --- | --- | --- | --- |
| **S0** | `sern_2d` 問題定義 (spec/derived/dv/作動点スキーマ、過拘束検査)・幾何コンテナ (ランプ/カウル折れ線、包絡検査)・CAPABILITIES 更新 | `probdef.py`, `geometry/sern_geometry.py`, `design/CAPABILITIES.md` | 小 |
| **S1** | 平面最大推力理論の導出と key point 逆設計: 非対称 2 壁前進 MOC (角部扇・カウル後縁自由境界)、乗数関係・縁条件、$c$ 探索、目標 C⁻ 張り、壁流線抽出。解析検算 (§6) | `geometry/moc_sern.py`, `geometry/rao_planar.py`, `design/tests/run_sern_moc_tests.py` | 中 |
| **S2** | 3 ブロック TFI メッシュ生成 (msh4.1 直書き)・品質ゲート・node 変換の配管 | `meshing/mesh_sern.py` | 中 |
| **S3** | 評価 runner (作動点ごとの run 生成、段階起動、warm restart)・力積分メトリクス・ゲート | `evaluate/runner_sern.py`, `metrics/sern_forces.py`, `evaluate/health.py` | 中 |
| **S4** | 検証: (a) 設計点 Euler で MOC の $C_T$ と forge の一致、(b) NASA TM X-71972 基準形 (ランプ 20°/18.5H、カウル 6°/3.12H) の M10/M4 傾向再現、(c) 外部流ブロックの妥当性 (case/23 プリューム資産と突合) | `case/46.sern_design/` (新設、README run 一覧) | 中 |
| **S5** | δ\* 一発補正の配管 (抽出 → オフセット → 再評価) と $\Delta C_T$ 台帳。非一様入口の評価側 BC | `feedback/deltastar_sern.py` (薄いラッパ), `runner_sern.py` | 小 |
| **S6** | MOO: 多作動点束ね・目的ベクトル・制約・DOE → パレート。1 評価時間の実測と予算更新 | `opt/driver.py`, `opt/moo.py`, `design/problem_sern_*.yaml` | 中 |
| **S7** | 3D 確認 run (側壁・隅 R)・2D との差の表・3D 設計要否の判定 | `case/46.sern_design/mesh3d/`, README | 中 |
| S8 (任意) | 二段膨張オプション: 基部の壁圧プラトー規定 (④延長部の壁圧帰還と共通機構) + 延長部は最大推力 | 別 plan に切り出す | 中 |

S0→S1 は CFD 不要で先行できる。S2–S3 は S1 と並行可 (固定形状で配管)。
- **機体上面の作り替えを実装・検証 (2026-09-06, run_0078/0079)**: 上式を後縁 θ_e 接線 Hermite + くさび 3°
  (`mesh.vehicle_wedge_deg`, 既定 3.0) に変更。**力係数は動かない** (C_T 0.9184 vs 旧 0.9187 = 0.03 % 差、
  run_0079 vs run_0077) ので、これは形状の健全化であって物理の変更ではない。
  **ただし発散は直らない** (run_0078, `pMin` 1.0 Pa で soft 段 step 61 に NaN)。NaN 位置は
  **x 0.42, y 0.223 = ノズル内部**で機体上面ではなく、上面の凸角を消しても床は別の場所で割れる。
  → **§4.11 の根因診断 (床が原因、曲率は結果でない) の追加証拠**。生産設定は `p_min: 20.0` のまま。
- **単点 CLI が生産レシピを使っていなかった** (2026-09-06): `runner_sern` / `runner_sern3d` の `main()` が
  `run_staged` を既定引数 (層流暖機なし・soft CFL 0.5・mid なし) で呼んでいた。`opt:` を読むよう修正
  (driver と同じ正本)。この罠で run_0078 を 1 本捨てた。
- **3D node SST が通った (2026-09-06, run_0082)** — §5.1-3 の決着。**真因は三重点そのもの**だった。
  2D の段階起動 (層流暖機 2000 + soft 1 次 cfl 0.2 + mid 2 次 cfl 0.2 + 本段 cfl 0.5) と `p_min: 20` を
  移植しただけでは足りず (run_0080 は暖機段こそ完走したが SST 投入直後の soft 段 **step 3** に同じ三重点で ω 発散、
  soft 段 step 0 の **rms_roK が 1.1e7** = 受け渡しの時点で乱流量が壊れている)。
  **`mesh3d.L_sw` を 0.8 (= 0.8 L_cowl) にしてカウル後縁と側壁後縁を x 方向に離す**と全段完走した。
  すなわち `L_sw = L_cowl` は**2 本のせん断層が 1 本の交線を共有する退化した配置**であり、避ければよい。
  実機の SERN も側壁とカウルを同じ station で終わらせる必然性は無く、§5.1-4 で $L_{\rm sw}$ は 3D dv の第一候補なので、
  設計空間から $L_{\rm sw} = L_{\rm cowl}$ を外すことに設計上の損失は無い。
  結果 (accel 点, `top_out` slip, 524628 セル, メッシュ PASS AR 535.8 / skew 0.365, 本段 4000 step 262 s):
  **C_T 0.9385 (圧力) / 0.9262 (摩擦込み) / C_T,friction −0.0123 / C_L −0.1787 / C_M +0.9892**、
  力係数は 8 スナップショットで **3 つとも STEADY**。同条件 3D Euler (run_0027: C_T 0.932 / C_L −0.204 / C_M +1.07) と
  整合し、差は摩擦の大きさ (−1.3 %) に見合う。面別の内訳では**ランプ幅外 (機体下面) が揚力の主役** (C_L,ramp_outside −0.0903 vs
  幅内 −0.0014) で、Euler で見えていた横方向膨張の描像は粘性込みでも変わらない。
  → **隅 R と 3D の δ\* を測る前提が揃った**。残差は plateau (NOT CONVERGED) だが §4.13 の受入方針どおり力の定常性で判定する。
- **3D 生産作動点はまだ通らない (2026-09-06, run_0083–0088)**。加速点 (run_0082, $p_{\rm in}/p_{\rm ext}$ = 5) は
  通るが、生産作動点 m6_on (同 **35**) は mid 段 (2 次) で落ちる。落ちる場所は毎回**カウル後縁の刃先**で、
  スパン端から 0.014–0.028 H 内側 (x/H 1.187–1.200, z/H 0.972–0.986)。潰した仮説:
  | 試行 | 結果 |
  | --- | --- |
  | `p_min` 20 → 60 Pa (run_0084) | **同じ step 48・同じ 22 ノード** → 圧力床は無関係 (2D の m10_on とは別機構) |
  | リミッタ Venkatakrishnan → Barth (run_0085) | step 34 に**悪化** → リミッタでもない |
  | カウル板厚 0.005 (run_0087) | step 48 → 57。z 一様に効かせた版 (run_0086) は側壁スリットの内外ずれで暖機段 step 4 発散 |
  | mid 段 CFL 0.2 → 0.1 (run_0088) | step 57 → 109 = ほぼ倍 ~~= 同じ物理時刻~~ → CFL はつまみでない。**撤回 (2026-09-13, R6(a), codex M4)**: 定常局所時間刻み (`unsteady: 0`) の step 数は物理時刻でない。「破綻までの擬似時間積分量が同じ」と読み替える |
  残る手は幾何: 板厚法則は TE で 0 に絞るので刃先の半角が **0.6°** しかない。2D では同じ刃先が持つので、
  効いているのは**側端との近接**。→ **カウル後縁を有限ベース (鈍頭) にする** (2D のランプ角部丸めと同じ、
  特異点を幾何で消す手当て)。base 面の追加が要るのでメッシャ側の作業。
### 4.14 / 4.15 3D 計算領域とトポロジ → [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md)

本節 (§4.14 3D 計算領域のトポロジ / §4.14.1・§4.15.1 codex 採用結果 / §4.15 3D の残り 3 件 /
§4.15.2 R4e の案 (d) / §4.15.3 R4e の確定設計 / §4.15.4 2D 有限ベース診断) は**切り出し済み**。
正本は [`plans/active/tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) の同番号の節。
3D の残作業 (R4c/R4d/R4e/R4f/R5/R5b) も同 plan の §5.1 が正本。

### 4.16 CFL・収束判定・起動レシピ → [`tooling-nozzle-sern-startup.md`](tooling-nozzle-sern-startup.md)

本節 (§4.16 CFL と `implicit_relax` の方針 / §4.16.1 収束の判定基準 / §4.16.2 CFL ramp と起動レシピ) は
**切り出し済み**。正本は [`plans/active/tooling-nozzle-sern-startup.md`](tooling-nozzle-sern-startup.md) の同番号の節。
ここに重複を置くと片方だけ更新されるので本文は持たない (AGENTS.md「詳細手順は本ファイルに重複記載せず、該当文書を参照する」)。
起動レシピの残作業 R-a〜R-h も同 plan の §5.1 が正本。

### 5.1 残作業 (優先順・2026-09-05 時点)

**本節が残作業の正本**。作業ログ側 ([`notes/sessions/2026-09-05-handover-sern-design.md`](../../notes/sessions/2026-09-05-handover-sern-design.md))
は使い方・罠の記録に限り、優先順はここを見る。

| # | 項目 | 内容 |
| --- | --- | --- |
| R1 | ~~**評価ゲート修復**~~ **完了 (2026-09-13, §4.7/§4.13, case/46 run_0090 再判定・run_0091 実機スモーク)** (codex C1 採用): 発散評価の採用を撤回 (§4.13-3 廃止)、`steadiness` の NaN→`STEADY` バグ修正、保存場の有限値・全残差・**実際に最適化する量 (`C_T_with_shear`)** の定常性を独立した必須ゲートに、`degraded` を `pareto.json` にも残す、数値失敗を物理的 `INFEASIBLE` と混同しない | `design/forge_design/metrics/sern_forces.py` (`steadiness` L108: `[1,1,1,NaN]` が STEADY を返す = 再現済)、`opt/driver_sern.py` (採用条件 L118 / 書き出し L139, L204)。accepted `tooling-nozzle-moo-loop.md` L33 (ゲート不合格は学習対象外) と揃える |
| R2 | ~~**力の集計範囲の分離**~~ **完了 (2026-09-13)** (codex M5 採用): メッシャ `vehicle` タグ (W/2 < z ≤ W_vehicle/2、外は top_out)、`forces3d` の $C_T,C_L,C_M$ はノズル面のみ・機体力は別枠 (`mesh3d.W_vehicle` で `Z_ext` から独立)、`metrics/sern_momentum.py` で BCONDS 全面の運動量収支を検算 (2D 閉じ残差 0.1–1 %, 3D 1.6 % of $F_{\rm ideal}$; 壁力は帳簿と 1e-5 [no-slip] / 0.1 % [slip 弱 BC] 一致)。再計算 run_0092 (3D, 6000 step) / run_0093 (2D, 8000 step): **ノズル $C_T$ 3D/2D = −1.9 %** (run_0027 再集計と同値 = タグ分離は帳簿だけを変える)、$C_L$ −0.10 vs 0.00、$C_M$ +0.72 vs 0.00 (§6)。3D の frozen_tp 配管 (bcond Y・IC・species_db) は書いたが未実走 | `mesh_sern3d.py` / `runner_sern3d.py` / `tests/run_sern_mesh3d_tests.py` (閉性・タグ分割 20 項目) |
| R3 | ~~**作動点の凍結 TP 化 + 外気 = 空気**~~ **完了 (2026-09-13, §4.5/§4.10)** (codex C2 採用 = §8-7 の (c)): `gas.model: frozen_tp` (排気 EXH = CEA 凍結組成の擬似種, 外気 AIR, thermalMethod 2)、`gas/frozen.py` (NASA-9 に H2/OH/H/NO/O/CO 追加, `FrozenGas`, 等エントロピー理想推力)、`gas_states` を排気/外気で分離、IC/BC/理想推力を同一物性で。外部動圧 3 作動点とも 71.9 kPa。検証: 単体 `run_sern_frozen_gas_tests.py` 40 項目 (NIST N2, 空気 R/γ/a, Ar で CPG 一致, MW = CEA 24.430, IC の T 反転厳密, q∞) + 実機 run_0094/0095 (§10)。**残**: 3D 実走、暖機 cfl 0.1 が全 dv 箱で足りるかは MOO 再取得 (R7) で確認、ソルバ側の多成分 implicit 安定化は別 plan | `case/46/problem_moo_frozen_tp_cycle3op.yaml` (生産), `cea/tmx_operating_points.py` |
| R4 | **2D/3D 外部領域と格子・領域独立性** (codex M6 採用) — **実装済 (2026-09-13)**: `mesh_sern3d` に `ext_top` (2D の機体上面テーパ + 自由流バンドを z 一様に押し出し, タグ `vehicle_top`)、`ramp_fillet`、`underside_far` を追加 (閉性・タグのテスト 30 項目)。領域独立性は `r4_domain_study.py` (case/46 run_0102: base / Z_ext 3 / x_out 4 / bot_depth 1.5 / top_depth 4 / 格子 1.25 倍, 3D Euler 加速点) — 結果は §10: 3D メッシャにランプ後縁以降の上側外部領域 (2D の `ext_top` 相当) を追加、固定した機体形状で遠方境界・出口距離・格子の感度を確認、メッシュ品質 PASS と領域独立性を別ゲートに | `meshing/mesh_sern3d.py` L192 (`top_out` で閉じている)。§4.11 の「作動点の性質なので形状によらず剥離しない」は 2 形状の一般化なので「設計箱全域」の根拠には使わない (剥離制約を外す判断 §8-6 は維持) |
| R4b | **近壁の 2 点を片付ける** (2026-09-13 起票, R3 の実機で顕在化): (i) ~~入口∩壁角の壁圧 1.75 $p_{\rm in}$ → SERN の変換で `mesh.nodeInletCornerWall: 1` を有効にし $C_L/C_M$ の差を測る~~ **測定済・却下 (2026-09-13, run_0100/0101)**: 角ノードは 1.39 → 0.91 $p_{\rm in}$ になるが上流ダクト壁全体が 0.88–0.96 $p_{\rm in}$ に沈む (無しは 1.005)、$C_T$ −0.3 %、m10_on は暖機で NaN。既定 0 のまま (`mesh.node_inlet_corner_wall` は配管だけ残す)。角 1 ノードの 1.39 $p_{\rm in}$ は水平壁なので $C_T$ に効かず、$C_L/C_M$ への寄与は 1 ノード分 (<0.1 %)、(ii) ~~M∞10 外気側カウル板下面の冷点 (T 床 50 K) と壁 T の回復温度超え → 壁の熱境界 (`spec.wall_thermal` 等温) か node 壁列の処方を決める~~ **決着 (2026-09-13, run_0103): 等温壁 1000 K を生産既定に** (残差 3.5 桁低下 = 断熱の 100 倍改善、T max 2320 K、冷点 143 K に緩和; 摩擦込み $C_T$ −0.4 %, $C_M$ +0.11)。残る冷点はカウル板 LE が入口面に露出する上流延長の産物 (水平壁で $C_T$ 不変)。Tw の値は §8-11 でユーザ確認 | §4.10 末尾の観察。`nodeInletCornerWall` は converter オプション |
| R6 (**2026-10-05 整理**: (a)(b)(c)(e) 完了。残 = (d) のうち **数値許容と検出目標** [R7a で事前固定済み: 許容は §8、検出目標 0.002] と **作動点別 C_M 窓の実機値** [MOO 本番 R7b の前、ユーザ判断]。以下の「残 = (a)(c)(e)…」は 2026-09-13 時点の古い記述) | **検証・問題定義の整合** (codex M4/M7/M8/M9 採用) — **(a) 完了** (§4.11 の「同じ物理時刻」を撤回・§4.7 で擬似時間の OSCILLATING を物理振動と解釈しないと明記; 2026-09-13)、**(c) 完了** (`run_sern_moc_tests.py` 6b: 縁条件残差 0 の点が等長拘束の下で C_T 最大であることを assert 付きで検算)、**(e) 完了** (§8-10 で帰還必須論を撤回済、§5.1-4 は L_sw を足した小規模探索を先にする記述)、**(b) 完了** (`check_quasisteady.py --series-csv` + `force_history.csv`、判定器一本化)、**(d) 一部完了** (§4.2 を 5 変数・2 目的に訂正、作動点別 C_M 窓 `opt.cm_window`、`L_ramp_max` 最終輪郭再検査 = §5.1-8b 決着)。残 = (a)(c)(e) と (d) の許容値明記: (a) 定常擬似時間の履歴を物理時間と解釈しない (§4.11 の「CFL 半減で破綻 step 倍 = 同じ物理時刻」と §4.7 の RSS/FSS 統計は撤回。振動を採用するなら dual-time で時間刻み・内部反復・統計窓の独立性を確認、P 床到達は診断指標)、(b) `check_quasisteady.py --quantity C_T,C_L,C_M` を接続し `steadiness` と判定器を一本化、(c) Rao 検証は同一ガス・作動点・長さ拘束で独立に行い Shyne 式 11–15 との対応を assert 付きで (`run_sern_moc_tests.py` L117 の掃引は assert なし)、§4.3 の「Pareto 端点に Rao 点が出なければ実装誤り」は撤回、格子・領域誤差の許容 (例 $\|\Delta C_T\| < 0.002$) を明記、(d) 問題定義を **5 変数・推力効率と長さの 2 目的**と明文化 (§4.2 の 6 変数・4 目的は `driver_sern.py` L33 と不一致、`theta_c` は `moc_sern.py` L387 で場から決まる)、$C_M$ は加重平均でなく作動点別の許容窓、`L_ramp_max` を最終輪郭で再検査 (§5.1-8b と統合)、(e) 「3D 最適化には MOC への帰還が必須」(§8-10-4) は撤回: MOC は形状パラメータ化として維持し、3D 評価器の成立後に `L_sw` を足した小規模探索で改善を測る。帰還は既存族の不足を実測してから別 plan | `methods/design/overview.md` L799 (無帰還・3D 確認) と本 plan を一致させる |
| **3D** | **3D 計算領域の残作業 (R4c/R4d/R4e/R4f/R5/R5b) は [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) §5.1 が正本** | 本体からは重複記載しない |
| ~~**R8**~~ (決定 2026-09-27、**完了 2026-10-01**: 段 (i)・(ii) と生産 YAML の輸送切替) | **種 DB をソルバ所有の lump 記法へ切り替える** (依頼元: 別セッション、plan `thermophysics-solver-owned-species-db.md` @ feature/gap-heating-precision `085e00e1`〜`5ba303c6`) | **ユーザ決定 (2026-09-27)**: (1) SERN の lump 名 `AIR` はソルバ内蔵の擬似種 `AIR` (cp/R 3.5) と衝突するので **SERN 側を `AMB` に改名** (内蔵は触らない)。(2) 取り込みは **main 経由で両ブランチを合わせる** (cherry-pick しない)。(3) 着手は **3D の R5h → R5j → R4d の後**。作業: runner_sern の config 生成を `physProp.species: [{name, lump, basis: mole}]` に (参考 `runner_axismach._apply_gas_to_config`・`composition.solver_species_config`)、EXH の構成種 (CO, H2, OH, H, NO, O) の生エントリだけ `species_db_external.yaml`、meta・後処理・`_species_signature` の `AIR`→`AMB`、AWS の clone 再ビルド (`FORGE_BIN`・変換器も同ビルド)。旧場からの継続は `restart_field.py --force-species` で 1 回だけ移行。**合格条件 (事前)**: 切り替え前後で代表作動点の推力・モーメント差がノイズ床以内 (ノイズ床 = 同一設定 3 反復の差の 3 倍と下限、case/44 の V0 方式)。予告: 輸送物性を CEA trans.inp に寄せる変更 (#5) で NS/SST の結果が変わる見込み (H2O μ +23 %・λ +47 % @600 K)。**更新 (2026-09-30、依頼元の正本 = feature/species-transport `plans/active/thermophysics-solver-owned-species-db.md` §5.2 @ `45d8d523`)**: (a) EXH の構成種 (CO/H2/OH/H/NO/O) も内蔵になったので `species_db_external.yaml` も生成 `species_db.yaml` も不要、`write_species_db` (runner_sern.py:528) による擬似種の書き出しをやめる。改名の対象は SPECIES_ORDER・tp_species の既定 (:88, :110)・`FrozenGas.from_mole(..., "AIR", ...)` (:108)・段間継承 (:729-732)・meta・後処理 (plot_sern*.py 等)・problem YAML。順序 (EXH, AMB) は変えない。(b) 輸送物性の指定が必須: viscMethod 2 は `physProp.transport` 無しで起動エラー (混合則 CEA frozen)。problem YAML に `gas.transport` (H2O は `custom:h2o_iapws_cea_v1`、N2/O2/AR) を足す。現行 viscMethod 1 (空気 Sutherland) からの切替で NS/SST の結果は変わる。(c) 種の照合が厳密: 属性の無い・未検証の場はソルバも Python ツールも停止、許可はその実行だけ `FORGE_ALLOW_UNVERIFIED_SPECIES=1` / `--force-species` (`FORGE_REQUIRE_VERIFIED_SPECIES` は撤去)。(d) 旧バイナリは新 config を読めない → AWS の clone も取り込んで再ビルド、`FORGE_BIN` で指定。(e) **回帰は 2 段** (依頼元の指定): (i) 改名 + lump 記法だけ (輸送は旧のまま) で代表作動点の推力・モーメントがノイズ床以内、(ii) 輸送物性の切替は変化量を記録 (不変を合格にしない)。凝縮関係の変更は凝縮 OFF なら影響なし。**着手時期 (ユーザ決定 2026-09-30)**: farfield plan §5.1 #4c の g4 の 2 本 (`case/46.sern_design/run_0999`・`run_1000`) が終わってから。取り込み方は 2026-09-27 の決定どおり main 経由 → **変更 (ユーザ決定 2026-09-30)**: species-transport が main に未マージ (plan in_progress、main との差 604 commit) のため、`origin/feature/species-transport` を **`feature/sern-design` に直接マージ**する (事前の merge-tree で衝突 5 ファイル: `output.cpp`・`restart_field.py` [add/add]・`plans/README.md`・case/48 README・削除済み plan 1 本)。main への統合は後で各ブランチから。**実装 (2026-09-30)**: マージ `237357fd` (衝突 5 件解消)、ランナー `5eee3ff2` (AMB 改名・lump 記法・`species_db.yaml` の書き出し停止、旧名 `[EXH, AIR]` は理由付きで拒否、試験 `run_sern_frozen_gas_tests`・`run_gas_tests` ALL PASS)。**依頼文との食い違い**: EXH の構成種 H2/OH/NO/H/O/CO は共通データにあるが `legacy_builtin: design` だけで、ソルバの内蔵 (`speciesDB.cpp:118`、`legacy_builtin: solver` のみ読む) には入らない → これら 6 種の**生の**係数は `species_db_external.yaml` で渡す (axis-Mach ランナーと同じ形、合成済み擬似種は書かない)。**回帰 (i) の事前登録 (2026-09-30)**: 問題 = `case/46.sern_design/problem_moo_frozen_tp_cycle3op.yaml`、作動点 m6_on、2D node SST、ランナー既定の段階起動 (暖機・soft・mid・本段 12000)。(0) 熱物性の照合: 旧ランナーが書く `species_db.yaml` の EXH/AIR と、新バイナリの `forge --resolve-species` が合成する EXH/AMB の MW・NASA-9 係数が相対 1e-15 以内。(1) A = 旧 (AWS `~/forge-pgrad-new`、ランナー c758e59d・バイナリ build-ff = a6ceee0b) × 3 反復、B = 新 (AWS `~/forge-r8`、5eee3ff2 のビルド) × 3 反復、同じ設定。比較量 = ランナー collect の C_T・C_T_with_shear・C_L・C_M (本段の末尾平均)。ノイズ床 N = max(3 × 群内の最大差 [A・B それぞれ], 1e-6)。**合格 = 全 run GATES PASS かつ 4 量とも \|平均_B − 平均_A\| ≤ N**。未達なら原因を切り分ける (合格条件は変えない)。(ii) 輸送物性の切替は (i) の合格後、変化量を記録 (合否なし)。**経過 (2026-09-30)**: AWS に作業ツリー `~/forge-r8` (8b87c0e7) を作りビルド — AWS の CUDA 13.2 は thrust が `include/cccl` にあるので `-DCMAKE_CXX_FLAGS="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl"` が要る (旧 build-ff も同じ指定。`cuda_forge` への CUDA include 追加 8b87c0e7 はその指定があれば不要だが無害)。**(0) の結果: FAIL (僅差)** — AMB vs 旧 AIR は MW・係数とも完全一致、EXH は MW 相対 2.8e-16 だが係数の最大相対差 **1.02e-15** (`nasa9_low[0]` 1560.6275407330086 vs …07、1〜5 ulp の差が 10 係数) で基準 1e-15 を超えた。合成の計算順序 (Python と C++) の丸めと見るが、解釈は (1) の結果とまとめて諮る。(1) A = `~/forge-pgrad-new/case/46.sern_design/run_1002_r8A_1`–`run_1004_r8A_3`、B = `~/forge-r8/case/46.sern_design/run_1005_r8B_1`–`run_1007_r8B_3`。**(1) の結果: PASS** — 6 本とも GATES PASS、本段末尾 6000 step 平均: C_T |Δ| 1.5e-6 (N 1.3e-5)、C_T_with_shear 2.0e-6 (1.2e-5)、C_L 8.6e-6 (4.1e-5)、C_M 1.8e-4 (9.4e-4)。4 量とも `check_quasisteady` STEADY、NaN・床 0、`check_convergence` は 6 本とも全列プラトー (本段で 0.1–0.6 桁、上昇なし、A・B 同じ = SERN の既知の性質)。**判断: 2026-09-30 codex (diagnose) [記録](../../notes/reviews/2026-09-30-sern-r8-stage0-thermo-diagnose.md) — 全件採用**: (0) は FAIL のまま「丸め差」と記録 (差は最大 **7 ulp**、5 ulp は誤り)、(1) の PASS だけで合格にせず実 run の記録で演算経路 A/B を追加、1e-15 を観測に合わせて緩めない。**演算経路 A/B (`case/46.sern_design/r8_stage0_path_ab.py`、CFD なし)**: A = 旧コードの `lump_entry` で旧 run の `species_db.yaml` を**ビット一致**で再現、B = C++ `synthesizeLump` の逐次演算を新 run の config・構成種生データで再計算し新 run の `resolved_species_*.yaml` を**ビット一致**で再現 → 保存値の差 (MW 2.8e-16、係数 10/18 個・最大 7 ulp) は演算経路 (質量↔モル換算・正規化・加算順と係数の相殺 [codex: 相殺倍率 17.8]) だけで説明できる。物性差は cp 最大相対 9.0e-16、顕内部エネルギー差/(cp·max(T,298)) 2.6e-15 (≤ 1e-12)。構成種の温度区切りは全て 200/1000/6000。→ **段 (i) を合格とする** ((0) は FAIL・丸めとして記録、(1) PASS)。段 (ii) へ進む。**段 (ii) 実装と結果 (2026-09-30、記録のみ)**: frozen_tp も `gas.transport` を受け付けるようにし (`probdef.py` の制限を解除)、runner_sern は `gas.transport` があるときだけ viscMethod 2 + `physProp.transport` (無ければ viscMethod 1 のまま = 既存 config と同一) — 9f95c764、試験 `run_transport_spec_tests` に frozen_tp の正例を追加 (ALL PASS)。問題 = `problem_moo_frozen_tp_cycle3op_transport.yaml` (生産 YAML の複製 + 全 11 構成種の `gas.transport`、H2O = custom:h2o_iapws_cea_v1、他 cea)。C = `~/forge-r8/case/46.sern_design/run_1008_r8C_tr_1`–`run_1010_r8C_tr_3` (新バイナリ、3 反復): 3 本とも GATES PASS・4 量 STEADY・NaN/床 0、残差は B と同じくプラトー。**B (Sutherland) → C (種ごと) の変化 (本段末尾 6000 step 平均、群内最大差は B/C とも ≤ 3.6e-4)**: C_T −2.7e-5 (−0.003 %)、C_T_with_shear +1.72e-4 (+0.019 %)、C_T_friction +1.99e-4 (+0.64 %、摩擦抵抗が減る)、C_L −8.0e-5 (−0.029 %)、C_M +1.65e-3 (+0.023 %)。→ 変化は主に壁摩擦 (排気 H2O 分の粘性の違い) で、§8 許容 (C_T・C_L 2e-3、C_M 5e-2) より 1 桁以上小さい。~~**未決 (ユーザ判断)**: SERN の生産 YAML を種ごとの輸送物性 (gas.transport) に切り替えるか~~ **決着 (2026-10-01 ユーザ「はい」)**: 生産 YAML (`problem_moo_frozen_tp_cycle3op.yaml`・`problem_3d_prod_m6on_*.yaml` 8 件) に `gas.transport` を入れる。作動点で実種が違う (m4_off は N2/O2/AR/CO2 の 4 種) ので、YAML には**全作動点の実種の和集合**を書き、runner_sern は作動点ごとにその実種へ絞ってから照合する (書き漏れはエラーのまま)。切替後の生産 run は旧 run と比べて力係数が約 0.03 % 以内で変わる (段 (ii))。**実施 (b96f4fa1)**: 生産 YAML 8 件に `gas.transport` (11 種の和集合)、runner_sern の `frozen_transport` が作動点の実種に絞る (m4_off は 4 種、試験 `run_transport_spec_tests` に正例・書き漏れの負例を追加、ALL PASS)。実ソルバの起動確認: m4_off を `--prepare-only` → 5 step (`run_case.sh`、AWS `~/forge-r8`) で viscMethod 2 + 4 種の transport を読み込み、NaN なしで完走。**R8 はこれで完了** (旧名の生産 run を継続するときは `restart_field.py --force-species` で移行)。**訂正 (codex diagnose 2026-10-01、[記録](../../notes/reviews/2026-10-01-farfield-r8-status-review-diagnose.md))**: 移行した場の出力は `species_input_unverified=1` のままで、標準の restart は毎回 `--force-species` (と起動時の `FORGE_ALLOW_UNVERIFIED_SPECIES=1`) が要る — 「1 回だけ」は誤り。仕様変更の要否は種 DB plan #3d (依頼元) の判断。段 (i) の合格は「(0) FAIL・丸め差として例外受理」と「(1) 準定常回帰 PASS (残差はプラトー、収束解の一致ではない)」の組み合わせに限定。「EXH 構成種も内蔵」は誤りと依頼元へ伝える |
| **R9** (予告受領 2026-10-01、**完了 2026-10-02**: 段 3 bfd8c3c8・#14 f28ca2fa、回帰 R9/R9b) | **種 DB 段 3 (CEA 熱物性の全面採用) の取り込み** (依頼元: feature/species-transport、plan `thermophysics-solver-owned-species-db.md` §5.1 #13-3、HEAD 561ba4d1 以降) | 内容: ソルバ内蔵の NASA-9 を CEA thermo.inp そのものに (N2/O2/CO2/Ar/He は 3 区間、H2O/He は CEA の分子量、CO/H2/OH/H/NO/O も内蔵に) → SERN の `species_db_external.yaml` には実種が出なくなる (R8 で記録した「内蔵種」の食い違いはこれで解消)。6000 K 未満の数値影響: SERN の構成種は変化 0、H2O は分子量の丸めで cp 相対 1.1e-6・h 15 J/kg、μ・λ ≤ 6.5e-7。**互換性ハッシュが全部変わる** → 進行中の restart 連鎖は係数不一致で止まるので、明示許可 (`--force-species` / `FORGE_ALLOW_UNVERIFIED_SPECIES=1`) で移行が要る。**取り込み条件**: (1) 依頼元から段 3 完了 (case/44 回帰 PASS) の連絡が来るまで取り込まない、(2) SERN の restart 連鎖の切れ目 (farfield plan #4e の 3D 輸送対照 A/B が終わった後) に合わせる、(3) 取り込み後は R8 と同じく回帰 (代表作動点の推力・モーメントがノイズ床内) を取る。`composition.py`・`semiperfect.py`・`evaluate/ic.py` に SERN 側の未 commit の変更は無い (2026-10-01 確認)。**依頼元へ伝えること**: 予告の「以後は新しい記録で継承される」は、こちらの観測 (codex diagnose 2026-10-01: `restart_field.py` は `species_input_unverified=1` の場を毎回拒否し、許可すると属性を消す → 継続のたびに許可が要る) と食い違う — 種 DB plan #3d の確認を依頼。**依頼元の返信 (2026-10-01)**: 指摘どおりで予告は誤り。現状はソルバが印付きの場をハッシュ一致なら通して印を継承する一方、Python ツール (restart_field / interp_field / convert_species_field / runner の段間継承) は毎回拒否・許可で印を消す。対処: ツールをソルバと同じ規則にそろえる (印付きかつ元のハッシュ = 宛先のハッシュなら許可なしで通し印を継承、属性なし・ハッシュ不一致は従来どおり停止; 種 DB plan §5.1 #3d、e5fdf2fd)。#3d は段 3 と一緒に入り、完了連絡に含まれる → **取り込み後は明示許可 1 回で移行、以後の継続は許可なし**。3D 輸送対照 A/B を先に済ませる順序と、取り込み条件 (完了連絡待ち・連鎖の切れ目・R8 と同じ回帰) は了承済み。**完了連絡 (2026-10-01)**: species-transport HEAD 9a048e6c で段 3 と #3d (5302b404) が入った。SERN 代表 config のハッシュは ff57c6ef1fdf15c2 → 28834c8ae07252e8、`species_db_external.yaml` は出なくなる。継続の規則: 取り込み後の最初の 1 回だけソルバを `FORGE_ALLOW_UNVERIFIED_SPECIES=1` で回す → その出力 (未検証の印付き) は以後 restart_field / interp_field / 段間継承が許可なしで通り印を継承、`--force-species` で写した場は属性なしになり次の起動だけ許可が要る。case/44 の再計算 (run_0532) は反復ばらつきの内側、SERN 構成種の 6000 K 未満の変化は H2O 分子量の分 (1e-6 程度) の見込み。既知の残り: `run_sern_frozen_gas_tests` の restart fixture (SERN 側は 2026-09-30 に試験区間だけ許可済み)。取り込み条件 2 (連鎖の切れ目) も満たす (3D 輸送対照 A/B 完了、farfield #4f は未投入のまま保留) **取り込み (2026-10-01)**: species-transport 9a048e6c を `feature/sern-design` に直接マージ (bfd8c3c8)。衝突は試験 2 件 (`run_gas_tests.py`・`run_sern_frozen_gas_tests.py`) だけで、AMB 名を残し Y_H2O を段 3 の値 0.2411089155 に。試験 `run_gas_tests`・`run_sern_frozen_gas_tests`・`run_transport_spec_tests` ALL PASS、CO/H/H2/NO/O/OH は `solver_builtin_names` に入り `write_species_db` は `species_meta.yaml` だけを書く (確認済み)。 **回帰の事前登録 (2026-10-01、run 前)**: 問題 = `problem_moo_frozen_tp_cycle3op_transport.yaml` (生産と同じ種ごとの輸送)、作動点 m6_on、ランナー既定の段階起動 (MOC IC から、restart なし → 許可は不要)。A = R8 段 (ii) の `~/forge-r8/case/46.sern_design/run_1008_r8C_tr_1`–`run_1010_r8C_tr_3` (バイナリ sha256 dc560910、R9 前)、B = 同じコマンドで `run_1014_r9_1`–`run_1016_r9_3` (bfd8c3c8 の再ビルド)。比較量 = C_T・C_T_with_shear・C_L・C_M の本段末尾 6000 step 平均。ノイズ床 N = max(3 × 群内の最大差 [A・B それぞれ], 1e-6)。**合格 = 6 本とも GATES PASS かつ 4 量とも |平均_B − 平均_A| ≤ N**。併せて記録: B の種ハッシュ (予告 28834c8ae07252e8)、B の run に `species_db_external.yaml` が無いこと。未達なら原因を切り分ける (合格条件は変えない)。 **回帰の結果 (2026-10-01): PASS** — B = `~/forge-r8/case/46.sern_design/run_1014_r9_1`–`run_1016_r9_3`。初回投入は段間継承 (`restart_by_index` → `_species_signature`) が段 3 の解決済み記録の区間可変形式 (`Tbounds`/`nasa9_intervals`) を読めず KeyError で停止 → 2 区間・区間可変の両方を読むよう修正し試験を追加 (dba4b582、`run_sern_frozen_gas_tests` ALL PASS)、同名で再投入。6 本とも GATES PASS、4 量 STEADY、残差は全列プラトー (A と同じ、上昇なし)。末尾 6000 step 平均 |B−A| / N: C_T 5.5e-7 / 5.0e-5、C_T_with_shear 1.4e-6 / 4.8e-5、C_L 2.0e-6 / 2.0e-4、C_M 4.3e-5 / 4.5e-3 (相対 ≤ 0.0007 %)。B の群内ばらつきは A の 4–5 倍 (C_M 1.5e-3 vs 3.6e-4; 原因は未調査、N はこれで決まっている)。種ハッシュ 3e96b116 → f0a9035e (EXH の区間が 200/1000/6000/20000 の 3 区間に。予告の 28834c8a は依頼元の代表 config の値で、比較の対象外)、`species_db_external.yaml` は書かれない。→ #4f (farfield plan) を新バイナリで投入してよい。 **追記の受領 (2026-10-01、取り込み・回帰・#4f 投入の後)**: 依頼元から「段 3 の取り込みは LJ 出典の既定変更 (種 DB plan §4.10・§5.1 #14、`physProp.ljSource`、既定 [gri30, svehla1962]) が入るまで待つ」(互換性ハッシュの変更を 1 回にまとめるため)。こちらは連絡の前に段 3 を取り込み済み (bfd8c3c8) なので、ハッシュは #14 で**もう 1 回**変わる。SERN への影響: `speciesDiffusionMethod` の既定は 1 (kinetic 混合平均) で SERN の config は指定していない → EXH/AMB の拡散係数が LJ (EXH は構成種の質量分率平均) で決まり、#14 で H2/OH/H/O/CO の値が変わると変わる。μ・λ は CEA モードなので影響なし。**扱い (2026-10-01)**: 走行中の farfield #4f は止めずに段 3 バイナリ (LJ = 現行値 = 今後の legacy_v1) で完走させる。#14 の完了連絡を受けたら取り込み、R8 型の回帰 (3 反復ずつ) を取る。その差は拡散係数の変化を含むので不変を合格にせず変化量を記録し (R8 段 (ii) と同じ扱い)、§8 許容・ε と比べて #4f の結論 (G・幅) に効くかを判断する。効く大きさなら #4f を新 LJ で回し直す。依頼元へは「段 3 は連絡の前に取り込み済み、#14 は改めて取り込む」と返す。 **#14 の完了連絡 (2026-10-02)**: species-transport 33520773 (段 3 + #14、SERN 代表 config のハッシュ 20ad6332)。既定 LJ = [gri30, svehla1962]、旧値は `ljSource: [legacy_v1]`。変化量の表 `notes/investigations/2026-10-01-svehla1962-lj/delta_table.md` (拡散係数 最大 +34 %、OH–H)。**取り込み (2026-10-02)**: f28ca2fa (衝突は `procedures/recommended-settings.md` の変更ログ 1 件、両方残す)。段 3 以降のソルバ差分は種データ・DB 読み込み・輸送 DB だけ (数値カーネルの変更なし)。試験 3 本 ALL PASS。**回帰 R9b の事前登録 (2026-10-02、run 前)**: (1) 2D (R8 型): A = `run_1014_r9_1`–`run_1016_r9_3` (段 3 バイナリ 2ac853cf)、B = 同じ問題・コマンドで `run_1024_r9b_1`–`run_1026_r9b_3` (f28ca2fa の再ビルド)。SERN は `speciesDiffusionMethod` 既定 1 (kinetic 混合平均) なので EXH (構成種の LJ の質量分率平均) と AMB の拡散係数が変わる → **合否なし、変化量を記録** (R8 段 (ii) と同じ扱い)。報告量 = C_T・C_T_with_shear・C_T_friction・C_L・C_M の末尾 6000 step 平均の B − A と、ノイズ床 N (R8 と同じ定義) との比。併せて EXH/AMB の解決済み LJ (σ, ε/k) の A・B 値、種ハッシュ、`species_db_external.yaml` 無し。(2) 3D (farfield #4f の結論の持ち越し判定): `run_1027_r9b_3d_g3_2p50` = run_1017 の最終場から新バイナリで 20000 step (同一格子・同一設定、restart は #3d の印継承)。**判定**: run_1017 との D が全 4 量 ≤ τ (0.2ε: C_T・C_T_with_shear・C_L 1e-4、C_M 1e-3) なら #4f の結論 (2.50 H・G+D §8 内) を新バイナリに持ち越す。超えたら #4f の 4 本を新バイナリで回し直す (同じ規則)。窓条件未達なら +20000。 **R9b の結果 (2026-10-02)**: バイナリ f28ca2fa 再ビルド (sha256 bec8db5f)。解決済み LJ: EXH σ 3.36827 → 3.36691 Å・ε/k 212.555 → 212.484 K (構成種の質量分率平均なので H2/OH 等の変化は薄まる)、AMB は不変。種ハッシュ f0a9035e → b72ff462、`species_db_external.yaml` 無し。(1) 2D: B = `run_1024_r9b_1`–`run_1026_r9b_3` 3 本とも GATES PASS・4 量 STEADY。B − A (末尾 6000 step 平均) / N: C_T −6.1e-7 / 5.0e-5、C_T_with_shear −6.4e-8 / 4.8e-5、C_T_friction +5.5e-7 / 2.2e-6 (0.25 N)、C_L −5.0e-6 / 2.0e-4、C_M +1.05e-4 / 4.5e-3 → 全量ノイズ床の 1/4 以下 (記録、合否なし)。(2) 3D: `run_1027_r9b_3d_g3_2p50` GATES PASS・窓条件 OK、restart は**ハッシュ変更のため許可なしは拒否** (想定どおり) → `--force-species`。run_1017 との D: C_T 6.4e-6・C_T_with_shear 6.4e-6・C_L 2.6e-5・C_M 5.8e-4、全量 ≤ τ → **#4f の結論を新バイナリに持ち越す** (事前規則)。**R9 完了** (段 3 + #14)。 **訂正 (codex diagnose 2026-10-02、[記録](../../notes/reviews/2026-10-02-farfield-4f-r9-result-diagnose.md)、全件採用)**: 「R9・R9b は力係数をノイズ床以下でしか動かさない」と一般化しない。N は A・B の群内最大差から作るので B のばらつきが増えると判定も緩む。記述は「登録 N 以下の群平均差」(R9 の C_M 平均差 4.3e-5 は旧 A だけから作る N_A 1.08e-3 よりも小さいので PASS は撤回しない) と「**R9 の B 群は再現性幅が増えた** (C_M 群内差 3.6e-4 → 1.5e-3、約 4.2 倍)、原因未確認」を併記。3D の τ はノイズ床ではない。残作業: 反復ごとの平均・振幅・前後窓差を回収し、過渡混入と反復間差を分ける (担当 O)。生産の限定運用は維持。 **回収結果 (2026-10-03)**: 2D 9 本 (R8C run_1008–1010 / R9 run_1014–1016 / R9b run_1024–1026) の本段末尾 6000 step: 各 run の窓内振幅 a は C_M 2.2–4.2e-3・C_L 1.0–1.9e-4・C_T_with_shear 2.5–4.8e-5、窓内の線形傾き ×6000 step は C_M ±0.3–3.2e-3 (符号は run ごとにばらばら)、前窓 (本段 0–6000 step) との差は全 run 同じ (C_M −2.2〜−2.3e-2 = 本段前半の過渡で、全 run 共通)。群内最大差は C_M 3.6e-4 (R8C)・1.5e-3 (R9)・5.4e-4 (R9b) で、**いずれも各 run の窓内振幅 a より小さい**。→ 群内差の違いは、有界な振れ (a ≈ 3e-3) を 6000 step 平均で切り取る位置の違いで説明でき、R9 だけが大きいのは標本 3 本の偶然と読む (ほぼ同じバイナリの R9b で 5.4e-4 に戻った; バイナリ由来の再現性悪化の証拠は無い)。**含意**: 「3 × 群内最大差」のノイズ床は、群内差が a より小さい限り実際の不確かさ (≈ a) を過小評価しうる。R8/R9/R9b/R10 の平均差 (≤ 1.05e-4) は a よりさらに小さいので結論は変わらないが、今後の回帰のノイズ床は max(3 × 群内最大差, 窓内振幅 a の群平均) とするのが妥当 (次の事前登録から適用、過去の判定は変えない)。 | — | F |
| **R10** (連絡受領 2026-10-03、**完了 2026-10-03**: 取り込み 5dd3a0df・回帰 PASS) | **species-transport の追加分の取り込み** (72 commit、2026-10-03 時点) | 中身: (a) SERN 入口組成の拡散係数の変化予測表 (ede50130、notes のみ)、(b) 二相拡散 (凝縮、opt-in・既定 OFF → SERN は凝縮 OFF なので動かない)、(c) リミッタ調査 (opt-in 診断のみ、既定経路不変)。依頼元: 急ぎ不要。**取り込む場合の条件**: 代表作動点 (2D m6_on) で同じ IC から短い回帰 1 本、推力・モーメントの差が事前に決めたノイズ床以内 (R9 の A 群 run_1024–1026 を基準にし、ノイズ床と区間は投入前にこの行へ書く)。**依頼元からの注意 (2026-10-03)**: (1) 既定リミッタ (limiterScaled 1) では少数節点でリミッタ係数が切り替わり続け残差が下げ止まる例 (case/16、原因未確定、`methods/limiter.md` の新節) — SERN のプラトー (衝撃・せん断層) は別機構の可能性が高い (依頼元の診断役の見立て)。収束は従来どおり報告量の `check_quasisteady` (判定区間と閾値を明記) で判定し、残差の位置を見て収束を合格にしない。(2) `space.limiterEpsConst` は撤去済みで存在しない (SERN の YAML・ランナーは未使用、2026-10-03 確認)。ψ 凍結は入っていない。(3) 参考: `notes/investigations/limiter-unstructured-convergence-survey.md`・`limiter-recommended-and-recent-survey.md` (取り込み後に読める)。取り込みの時期は restart 連鎖の切れ目でユーザと決める。**取り込み (2026-10-03)**: 5dd3a0df (衝突 2 件: `solver_density_cuda/main.cpp` は依頼元の `assembleResidualPre/Post` 分割に帳簿の診断フック (FORGE_DUMP_LEDGER、既定 no-op) を移して両立、`plans/README.md` は両方の行を残す)。試験 3 本 ALL PASS。**回帰の事前登録 (2026-10-03、run 前)**: 問題・コマンドは R9b と同じ (`problem_moo_frozen_tp_cycle3op_transport.yaml` m6_on、ランナー既定の段階起動、MOC IC)。A = R9b の `run_1024_r9b_1`–`run_1026_r9b_3`、B = `run_1030_r10_1` (5dd3a0df の再ビルド、1 本)。比較量 = C_T・C_T_with_shear・C_L・C_M の末尾 6000 step 平均。ノイズ床 N = max(3 × A 群内の最大差, 1e-6)。**合格 = B が GATES PASS・4 量 STEADY かつ 4 量とも |B − 平均_A| ≤ N**。未達なら原因を切り分ける (合格条件は変えない)。予測: 凝縮 OFF・既定リミッタ経路不変なので数値は変わらない (同一 IC でも atomicAdd の非決定性があるのでビット一致は求めない)。**結果 (2026-10-03): PASS** — `~/forge-r8/case/46.sern_design/run_1030_r10_1` (sha256 ee6f0c09)、GATES PASS・4 量 STEADY。B − 平均_A / N: C_T −2.7e-6 / 1.8e-5、C_T_with_shear −3.1e-6 / 1.8e-5、C_L −1.4e-5 / 7.2e-5、C_M +3.1e-4 / 1.6e-3。**R10 完了** | — | O |
| **R11** (連絡受領 2026-10-05、**完了 2026-10-05**: ff 77318d0e・回帰 PASS) | **main 統合 (PR #3、merge 77318d0e) の取り込み** | 依頼元の連絡: sern-design (4703a727) と species-transport を main に統合、sern-design は fast-forward 可。増分 50 commit (ソルバ 23 ファイル): 二相拡散既定化 plan の S0 と診断 G1–G3 (二相拡散は既定 OFF、診断は `FORGE_DIAG_TP_*` 設定時のみ)、凝縮・化学種・受動種カーネルと main.cpp (SLAU・勾配・リミッタ・境界は不変)、衝突解消 3 件 (species_advection_faceY_d の引数 ffY・tpo、stage_manifest の hard キー twophase_diffusion_effective は凝縮 run のみ、test_farfield_flux.cu を NASA-9 区間可変に追随)。依頼元の回帰は case/16・case/44 で PASS、case/46 は未実施。**取り込み (2026-10-05)**: `git merge --ff-only origin/main` → 77318d0e (衝突なし)。設計側試験 3 本 ALL PASS。**回帰の事前登録 (2026-10-05、run 前)**: R10 と同じ問題・コマンド (`problem_moo_frozen_tp_cycle3op_transport.yaml` m6_on、MOC IC、段階起動)。A = `run_1024_r9b_1`–`run_1026_r9b_3`、B = `run_1033_r11_1` (77318d0e の再ビルド、`FORGE_CUDA_BLOCKSIZE=128`、1 本)。比較量 = C_T・C_T_with_shear・C_L・C_M の末尾 6000 step 平均。**ノイズ床 N = max(3 × A 群内の最大差, A 群の窓内振幅 a の平均, 1e-6)** (R9 回収結果 2026-10-03 の含意を本件から適用)。**合格 = B が GATES PASS・4 量 STEADY かつ 4 量とも |B − 平均_A| ≤ N**。未達なら切り分ける (合格条件は変えない)。**結果 (2026-10-05): PASS** — `~/forge-r8/case/46.sern_design/run_1033_r11_1` (sha256 eb0ea331)、GATES PASS・4 量 STEADY。B − 平均_A / N: C_T −8.0e-6 / 4.1e-5、C_T_with_shear −7.9e-6 / 4.2e-5、C_L −2.4e-5 / 1.4e-4、C_M +5.4e-4 / 3.2e-3 (旧定義の 3 × 群内差だけでも全量内: C_M 1.6e-3)。**R11 完了** | — | O |
| ~~R7~~ (**2026-10-05 R7a/R7b に分割**、下の 2 行) | **小規模探索で判別能力を確認してから MOO 再取得** | R1–R6 の後。#1 の run_0054 系は R1 のゲートで再判定 (発散評価を含む Pareto は使わない)。早期停止 (#1c) は R7 の後 |
| **R7a** (**最優先**、判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-next-priority-diagnose.md) — 全件採用) | **3D 判別試験: `L_sw` 2 水準 × 3 反復** | 固定格子 (g3)・単作動点 (m6_on)・MOC 輪郭固定で、3D 評価が設計差を見分けられるかだけを認定する (多作動点 MOO の成立確認ではない)。**事前登録 (2026-10-05、run 前)**: A = `L_sw` 0.8 (現行生産、実位置は既存 station に丸め ≈ 0.8091 H)、B = `L_sw` 1.0 (≈ 1.0043 H)、生成後の実位置を台帳に残す。各 3 反復 (`run_1034`–`1036` = A、`run_1037`–`1039` = B)。固定: MOC 輪郭・幅・側壁以外の形状・z 座標列・`first_z_frac`・`sz` テーパ・物性・BC (生産 YAML `problem_3d_prod_m6on_wallres.yaml`、B は `L_sw` だけ変えた複製)・リミッタ基準値 (run_0986 値)・バイナリ (eb0ea331)・`FORGE_CUDA_BLOCKSIZE=128`、診断介入なし。初期場: A は `run_1017_ff4f_g3_2p50` 最終場から `restart_field.py` (同一格子)、B は runner で格子を作り同じ最終場から `interp_field.py` (壁内外ノードの取り違え・新しい壁距離を確認)。B が直接継続で発散したら段階起動を適応区間とし、その後を判定区間にする (区間は stage_manifest で分ける)。長さ: 判定区間を各 20000 step・500 step 出力。末尾 10000 と直前 10000 の平均差が C_T・C_T_with_shear・C_L ≤ 5e-5、C_M ≤ 5e-4 (未達なら各 +20000、なお未達は判定不能)。必要条件: メッシュ品質・有限値・床ゲート・4 量 `check_quasisteady` STEADY、`check_convergence` は区間付きで併記 (プラトーだけで失格にしない、§4.39)。**判定 (q = C_T_with_shear)**: 群ごとに末尾平均の最大差 r・最大窓内振幅 a・最大前窓差 d、N = max(3r, a, d, 1e-6)、E = N_A + N_B。分岐 A: |平均_B − 平均_A| − E > 0.002 → 二水準を識別できる (符号から優位側を記録し、作動点別制約と格子・領域感度へ)。分岐 B: |平均_B − 平均_A| + E ≤ 0.002 → この二水準では 0.002 規模の利得を示せない (大規模 DOE に進まない)。他は判定不能。0.002 は検出目標で、数値誤差の許容 (§8) とは別の量。やらないこと: R5h・W2 の再実施を前提に戻す、R5h の単一形状の感度を設計箱全域の誤差上限にする、2D driver に `L_sw` を足して 3D 最適化と扱う、g3/g4 差を誤差上限と呼ぶ。**結果 (2026-10-05、事前規則で判定)**: 6 本とも GATES PASS。A (`run_1034`–`1036`) は 20000 step で窓条件 OK。B (`run_1037`–`1039`、interp 起点) は窓条件未達 (C_T_with_shear 前窓差 2.6e-4) → 事前規則で +20000 (`run_1040`–`1042`) → 全量 OK。末尾 10000 step 平均 (B は延長 run): C_T_with_shear A 0.8705560 / B 0.8726419、**B − A = +0.002086**、N_A 3.3e-6・N_B 3.6e-6、E 6.9e-6 → |Δ| − E = 0.002079 > 0.002 → **分岐 A (二水準を識別できる、検出目標をわずかに上回る)**。他の量: C_T +0.00343、C_L +0.00898、**C_M −0.213** (§8 の C_M 許容 0.05 の 4 倍超; 側壁を延ばすとピッチモーメントが大きく動く)。原本 `notes/investigations/2026-10-05-sern-r7a/R7A_VERDICT.txt`。**解釈は未確定** (エスカレーション条件 7、codex 諮問中): 余裕が検出目標の 4 % しかない点、B の起点が双子節点を同値で埋めた補間だった点、C_M の大変化の扱い (作動点別 C_M 窓が未定)。**判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-r7a-result-diagnose.md) — 全件採用**: (1) 分岐 A は「g3・m6_on・今回の初期化で登録基準を満たした」までの限定付きで保持 (E は反復・窓内変動の指標で初期場依存・格子・領域誤差を含まない)。(2) `interp_field` は既存の座標一致の双子節点 (厚さ 0 の板の両側) も内外を区別せず写す。B の 3 反復は同じ入力の複製なので共通の写像誤差は群内差に出ない → 写像だけを変えた対照が要る。(3) ΔC_M −0.213 は形状間の設計応答であって §8 の 0.05 (同一形状の格子・領域許容) と比べる量ではない (基準点 x_ref = −20 H の長い腕で増幅、x_ref = 0 換算で約 −0.034)。(4) 現行 C_M はノズル面の圧力モーメントだけ (摩擦・機体面なし) で、実機トリム窓と比べるには面集合・基準点・正規化・摩擦をそろえる必要。(5) R7a だけで MOO 本番へ進まない: 設計差の格子依存は G_B − G_A で決まり、既存 G は A 形状のもの。**写像対照の事前登録 (2026-10-05、run 前)**: 同じ B 格子・設定・バイナリで新規 1 本 `run_1043_r7a_lsw10_map` → 継続 `run_1044_r7a_lsw10_map_c40k` (計 40000 step、最後の 20000 step を判定区間)。起点は run_1017 res_20000、**初期写像だけ変更**: 座標で対応づけ、A 側で座標の一致する既存の双子は境界タグの組で見分けて保存量を直接コピー、B で新しく分かれた節点 (A では 1 節点) は両側ともその A 節点の値 (従来と同じ値)、`wall_dist` は新格子の値 (`r7a_twin_restart.py`)。投入前に、共通部の保存量一致と、従来 IC との差 (変更節点数・位置・量別最大差) を記録。必要条件: 前後 10000 step の窓条件・4 量 STEADY・有限・床・品質ゲート、残差 VERDICT 併記。未達は判定不能。比較幅 U = 既存 B 群の N + 新 run の max(a, d, 1e-6)。**結果 A**: 全量で |平均差| + U が C_T・C_T_with_shear・C_L ≤ 5e-5、C_M ≤ 5e-4 → 大きな初期場依存を支持しない (短側壁群との比較で分岐 A が維持されるかも確認)。**結果 B**: いずれかで |平均差| − U が閾値超え → 初期化非依存を棄却し R7a の認定範囲を訂正。他は判定不能。やらないこと: 同じ入力の 4 本目で初期場依存を調べたことにする、余裕が小さいことだけで事前判定を撤回する、E を全数値誤差の上限と扱う。**写像対照の結果 (2026-10-05)**: 写像変更で初期場は 5110 節点で ρ が変わった (最大 Δρ 0.33、入口の隅 x −0.05 z 0.1 = 従来の interp が厚さ 0 の板の反対側 [排気/外気] を拾っていた; 他の量は interp の原始量→保存量の組み直しの丸めで広く微差)。`run_1043` (20000、窓条件未達) → `run_1044` (+20000、全量 OK・GATES PASS)。B 群平均との差 / U / |差| + U: C_T −4.9e-7 / 6.9e-6 / 7.3e-6、**C_T_with_shear −4.3e-7 / 6.8e-6 / 7.2e-6**、C_L −1.2e-6 / 2.4e-5 / 2.5e-5 (いずれも結果 A の閾値 5e-5 内)、**C_M +3.0e-5 / 5.1e-4 / 5.40e-4 (閾値 5e-4 をわずかに超え中間)** → 事前規則の全体判定は**判定不能** (C_M で U が窓内振幅 2 本分だけで閾値に達するため。差そのものは 3e-5)。推力 2 量と C_L は初期写像に依存しない (結果 A の条件内)。B を写像対照で置き換えても短側壁群との比較は |Δ| − E 0.002079 > 0.002 で分岐 A は維持。原本 `notes/investigations/2026-10-05-sern-r7a/R7A_MAP_VERDICT.txt`。**次 (R7b の前提、codex 2026-10-05)**: 設計差の格子依存 G_B − G_A (B 形状でも g4 を取る)、C_M の定義を実機トリム窓とそろえる、作動点別 C_M 窓 (ユーザ判断)。**作動点別 C_M 窓 (2026-10-05 ユーザ「特に当てもない」)**: 実機の値・出典は無い → 保留。MOO の C_M 制約はそれまで現行の加重平均 `cm_min` −7.0 のまま (定義も現行 = ノズル面の圧力モーメント・x_ref −20 H)。**g4 での設計差の事前登録 (2026-10-05、run 前)**: A_g4 = `run_1045_r7a_g4_lsw08` (`problem_3d_prod_m6on_g4.yaml`、run_1021 [g4・A 形状] 最終場から restart_field)、B_g4 = `run_1046_r7a_g4_lsw10` (`problem_3d_prod_m6on_g4_lsw10.yaml`、格子を作り run_1021 の場を `r7a_twin_restart.py` で双子の内外を保存して写す)。バイナリ eb0ea331・設定は run_1021 の時間・space、各 20000 step、窓条件 (R7a と同じ) 未達なら +20000、なお未達は判定不能。g4 は各 1 本 (反復なし) なので幅は U = Σ max(a, d, 1e-6) (各 run)。**判定**: Δ_g4 = C_T_with_shear(B_g4) − (A_g4)。|Δ_g4| − U > 0.002 → g4 でも識別できる (R7a の分岐 A が格子で維持)。|Δ_g4| + U ≤ 0.002 → g4 では識別できない (R7a の結論は g3 限定と訂正)。他は判定不能。併せて記録 (合否なし): 4 量の G_A = A_g4 − 平均(A_g3: run_1034–1036)、G_B = B_g4 − 平均(B_g3: run_1040–1042)、**G_B − G_A** (= 設計差の格子依存)、側壁の実位置 (g4 の station に丸め)。**結果 (2026-10-05)**: 側壁の実位置は g3 と同じ (0.080910 / 0.100433 m)。A_g4 `run_1045` は 20000 step で全量 OK・GATES PASS。B_g4 `run_1046` は未達 → +20000 の `run_1048` で推力 2 量・C_L は OK、**C_M だけ前窓差 7.1e-4 > 5e-4** (9.75e-3 → 7.1e-4 と減衰中) → **事前規則で判定不能** (延長 1 回後も未達)。記録 (合否なし、原本 `notes/investigations/2026-10-05-sern-r7a/R7A_G4_VERDICT.txt`): C_T_with_shear Δg3 +0.002086 / **Δg4 +0.002086** (U_g4 6.3e-6、参考 |Δg4| − U 0.002079)、G_A −1.947e-4・G_B −1.949e-4・**G_B − G_A −1.3e-7**; C_T Δ +0.003434 / +0.003442、G_B − G_A +8.0e-6; C_L Δ +0.008975 / +0.009012、G_B − G_A +3.7e-5; C_M Δ −0.2132 / −0.2142、G_A +0.0168・G_B +0.0158 (B は未定常)、G_B − G_A −1.0e-3。→ ~~推力・揚力の設計差は g3/g4 で変わらない~~ (下の判断で訂正)。C_M は B_g4 が未定常のため未確定。扱いは諮問 (エスカレーション条件 3)。**判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-r7a-g4-diagnose.md) — 全件採用**: (1) 旧 g4 試験は **判定不能のまま確定・保存** (結果を書き換えない)。(2) C_M を判定から外す救済は却下 (C_M は最適化の制約にも使う)。(3) 書ける範囲 = 「**g3 で登録基準を満たした (分岐 A)。g4 の参考中心値も約 +0.002086 だが全量判定は未成立**」。「設計差の格子依存は 1e-7 級」は撤回 (G_B − G_A −1.3e-7 は比較幅 E_g3 + U_g4 = 1.3e-5 の 1/100 で分解できていない)。(4) 写像対照は「今回の写像変更による差は推力 2 量・C_L で指定許容内」まで (B を写像対照に置き換えると差の差の中心値が −1.3e-7 → +3.0e-7 と動く)。(5) 実機トリムの認定は保留、`cm_min` −7.0 は探索用の暫定制約として維持 (実機窓を作らない)。C_L が OK で C_M が未達なのは矛盾しない: x_ref = −20 H の腕で C_L の変化 3.6e-5 だけでモーメントに約 7.1e-4 の項が乗る (原因の確定ではない)。**追加診断の事前登録 (2026-10-05、run 前、旧試験とは別)**: `run_1049_r7a_g4_lsw10_c60k` = `run_1048` の最終場から restart_field、変更は反復数だけ (格子・設定・バイナリ・出力 500 step 固定)、**固定長 20000 step** (途中で止めない)、追加区間の前半/後半 10000 で判定。必要条件: C_T・C_T_with_shear・C_L 前窓差 ≤ 5e-5、C_M ≤ 5e-4、4 量 STEADY、有限・床・品質、`check_convergence` VERDICT と区間を保存。**結果 A** (満たす) → 反復不足説を支持し、A_g4 との推力差と U を再計算、|Δ| − U > 0.002 なら追加試験として g4 の識別成立を記録。**結果 B** (満たさない) → 固定長延長で解消する仮説を棄却し判定不能を維持 (再延長・閾値変更はしない)。併せて force_history・residual_history・個別判定・restart 来歴・集約スクリプトを保存。**追加診断の結果 (2026-10-05): 結果 A** — `run_1049_r7a_g4_lsw10_c60k` (累積 60000 step) で 4 量とも窓条件 OK (C_M 前窓差 1.9e-4 ≤ 5e-4)・STEADY・GATES PASS (`check_convergence` は NOT CONVERGED [プラトー] を併記、§4.39)。A_g4 との比較: C_T_with_shear Δg4 +0.002085、U 2.0e-6、|Δ| − U 0.002083 > 0.002 → **追加試験として g4 でも識別成立** (旧 g4 試験の判定不能は書き換えない)。記録: G_A −1.947e-4・G_B −1.956e-4・G_B − G_A −8.8e-7 (C_T_with_shear)、C_T +7.5e-6、C_L +2.6e-5、C_M −6.8e-4 (比較幅 E_g3 + U_g4 の内側で、分解できた格子感度とは言わない)。原本 `notes/investigations/2026-10-05-sern-r7a/R7A_G4X_VERDICT.txt`。**R7a のまとめ**: g3 (3 反復) と g4 (追加試験) の両方で、`L_sw` 0.8 → 1.0 の推力差 +0.0021 (+0.24 %) を登録基準で識別できた (m6_on・固定基準値・この 2 水準に限る)。C_M の設計応答は約 −0.21 (x_ref −20 H) | AWS | O |
| **R7b** (R7a の後。**順序 (codex 2026-10-05)**: ① R7a の g4 追加診断 → ② ~~現行 C_M の仕様を methods と plan に明記~~ **済 (2026-10-05、`methods/design/overview.md` の「$C_M$ の定義」)** → ③ ~~`L_sw` を物理 station として指定・再現できるようにする~~ **実装 (2026-10-05)**: `mesh3d.L_sw_exact: true` でノズル区間 [0 (角丸めありは xf2), L_cowl] の station を区分線形に写し最寄り station を L_sw に置く (両端固定、間隔比の変化 ~1 %、上限 1.25 を超える・区間端と重なる指定は生成失敗)。既定 false は従来どおり丸め (既存 run の形は不変)。試験 `run_sern_mesh3d_tests` に 6 項目 (指定値どおり・0.8 と 0.8001 が別形状・閉性・Jacobian・区間端の拒否) — ALL PASS。生産格子での品質確認 (2026-10-05、AWS で prepare のみ): `L_sw` 0.8 + `L_sw_exact` で g3 は `PASS` (skew max 0.449)、g4 は `SOFT-PASS` (skew max 0.701) — どちらも元の格子と同じ判定・同じ skew、側壁後縁 x = 0.0800 m (指定どおり)。**③ 済** → ④ その形状で他作動点の評価成立と設計差の格子・領域感度 (**判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-r7b-step4-diagnose.md) — 全件採用**: 2D の再診断を前提にしない [`chi_evidence` の run_0968 m10_on・run_0966 m4_off は 2D で 4 量 STEADY・残差プラトー、ただし viscMethod 1・旧種 DB]; `L_sw_exact` の形状は力係数の不変を保証しないので exact 形状で m6_on の二水準を取り直す; 他作動点は 3 作動点とも A/B を対にして格子・領域感度; 0.002 は検出目標で各作動点の必須改善量ではない (加重目的 ΔJ = 0.5Δ₆ + 0.3Δ₁₀ + 0.2Δ₄、E_J = 0.5E₆ + 0.3E₁₀ + 0.2E₄)。**順序**: 4-1 exact 形状・m6_on・g3 の 0.8/1.0 対 → 4-2 m10_on → m4_off の g3・A 形状の起動成立 → 4-3 両作動点の B → 4-4 3 作動点 × A/B の g4 → 4-5 3 作動点 × A/B の領域感度 (一方向ずつ + 同時拡大、共通領域を保持して外側にセルを足す) → ⑤。**4-1 の事前登録 (2026-10-05、run 前)**: A = `run_1050_r7b_x_lsw08` (`L_sw` 0.8 + `L_sw_exact`)、B = `run_1051_r7b_x_lsw10` (1.0 + exact)、他は R7a と同じ (m6_on・g3・生産 YAML の複製・時間/space は元 run を継承 = リミッタ基準値固定・eb0ea331・blocksize 128)。初期場: A は `run_1034` (A・g3 の発達場)、B は `run_1040` (B・g3) から**節点番号で写す** (`r7b_index_restart.py`: 節点数・各節点の境界タグの組が一致し、座標のずれがノズル区間の x だけで半間隔以下であることを検査してから保存量を index コピー、wall_dist は新格子の値)。各 1 本、判定区間 20000 step・500 step 出力、窓条件未達なら事前登録の +20000 を 1 回、なお未達は判定不能。q = C_T_with_shear、U = Σ max(a, d, 1e-6)。**結果 A** (ゲート成立かつ |Δ| − U > 0.002) → exact 形状への持ち越しを支持 (各 1 本の予備判別。反復込みの認定には各 3 反復で取り直す) → 4-2 へ。**結果 B** (|Δ| + U ≤ 0.002) → 持ち越しを棄却。他は判定不能。**4-1 の結果 (2026-10-05): 結果 A** — `run_1050` (A、後縁 0.0800 m)・`run_1051` (B、0.1000 m) とも 20000 step で窓条件 OK・4 量 STEADY・GATES PASS。C_T_with_shear A 0.8704579 / B 0.8725856、Δ +0.002128、U 6.2e-6、|Δ| − U 0.002122 > 0.002 → exact 形状への持ち越しを支持 (各 1 本の予備判別)。他: C_T Δ +0.003514、C_L +0.009323、C_M −0.2209。丸め版との差 (exact − 丸め) は A で C_T_with_shear −9.8e-5・C_L −5.3e-4・C_M +0.012、B で −5.6e-5・−1.9e-4・+0.0045 — 側壁後縁が A で 0.91 mm・B で 0.43 mm 短くなった形状差を含む (A→B の傾き dC_L/dL_sw ≈ 0.47 /m からの見積もり A −4.3e-4・B −2.0e-4 と整合、station 写像自体の寄与は切り分けていない)。原本 `notes/investigations/2026-10-05-sern-r7a/R7B_X_VERDICT.txt`。次は 4-2。**4-2 の事前登録 (2026-10-05、run 前)**: 問題 = `problem_3d_prod_3op_wallres_lswx08.yaml` (exact A 形状の生産 3D に、2D MOO YAML と同じ m10_on・m4_off の定義を追加、`evaluate.limiter_ref` を run_0986 の値に固定、輸送は作動点の実種に絞られる)。m10_on → m4_off の順に各 1 本、g3。起動 = runner_sern3d の既定の段階起動 (生産 YAML の層流暖機 2000 [CFL 0.2] → SST 1 次 2000 → 2 次 2000 → 本段 4000 [CFL 0.25])、MOC の領域別 IC から (m6_on の場からの warm start はしない)。その後、本段と同一設定で restart して 20000 step を**判定区間**とし (段階起動の区間は連結しない)、窓条件未達なら +20000 を 1 回。run 名: m10_on `run_1054_r7b_m10_A` (段階起動) → `run_1055_r7b_m10_A_c` (判定区間) [→ `run_1056` 延長]、m4_off `run_1057_r7b_m4_A` → `run_1058_r7b_m4_A_c` [→ `run_1059`]。**起動成立の条件**: 全段が NaN なく完走 (rc 0、detectNaN 無し)。**評価成立の条件 (判定区間)**: GATES PASS・4 量 STEADY・窓条件 (C_T・C_T_with_shear・C_L ≤ 5e-5、C_M ≤ 5e-4)・床/ω 下限への張り付き 0・化学種の有限・非負・組成和、`check_convergence` は区間付きで併記 (プラトーは失格にしない)。壁解像は `check_wall_resolution.py` で局所 y₁⁺ を記録 (m6_on の緩和を自動適用しない)。成立しなければ症状を記録して止め、諮る (2D に戻るかを含む)。**4-2 の結果 (2026-10-05)**: 段階起動は両作動点とも完走 (NaN なし)。**m4_off** (`run_1057` → 判定区間 `run_1058`): GATES PASS・4 量 STEADY・窓条件 OK (前窓差 ≤ 2.6e-4) → **評価成立** (C_T_with_shear 0.9676590、C_L −0.02494、C_M +0.5065)。壁解像 (局所 y₁⁺) は未計測。**m10_on** (`run_1054` → `run_1055`): 4 量 STEADY・窓条件 OK (C_T_with_shear 0.8970375、C_M −1.7785) だが **GATES FAIL `FLOOR_STUCK` (T ≤ 50 K が 1 節点)** → 事前規則で**評価不成立**。その節点は (x 0.12000, y −0.01053, z 0.0979) = カウル後縁 (厚さ 0) で、R5h の低温点と同じ場所 (m6_on では 95–105 K)。2 番目に低い節点は 211 K で、孤立した 1 点。→ R5h の再開条件 (「新形状・新作動点で床到達が増えたら再開」、3D plan §5.1 R5h) に該当。次の手は諮る (エスカレーション条件 2/3)。**判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-r7b-m10-te-floor-diagnose.md) — 全件採用**: R5h を再開 (恒久修正は保留)、有限厚後縁の先行は却下 (形状・帳簿が同時に変わる)、**床ゲートからの除外は却下** (TP は温度反転後に roe を再構成するので床到達は表示だけの問題でない)、旧 m6_on の感度で無影響とするのは却下、**m10_on 抜きで 4-3 以降を進めるのは却下** (先行は m4_off の壁解像確認だけ)。観測と解釈の分離: 「50 K の 1 点」は観測、「厚さ 0 の壁ノード」「旧機序と同一」「積分力に効かない」は未確認 (旧低温点は後縁直下の**内部節点**と記録)、`FLOOR_STUCK` は最終場の検査で判定区間全体の張り付きは別に見る。**R5h-m10 診断 A/B の事前登録 (2026-10-05、run 前)**: 起点 `run_1055_r7b_m10_A_c` 最終場を restart_field で 2 run に分岐、変更は `FORGE_DIAG_FACE_VEL_CELL` の有無だけ (A = 対照、B = 節点 517160 に接する後流側 1 面の同節点側速度 3 成分だけセル値へ; 面 ID は今回の接続から帳簿ダンプで特定、旧 1535830 は流用しない)。(1) 機序: 最初の共通状態の 1 更新を帳簿で記録 — 対象点・相手点・隣接 CV の面流束和・対流/粘性残差・実保存量増分・EOS 前後の変更を分け、内部エネルギー射影 R_ρE − u·R_ρu + (|u|²/2 − e_int)R_ρ と実更新を別に確認、閉合誤差 ≤ 介入差の 1 % または丸め上限の大きい方。(2) 100 step で早期確認し、機序を支持するときだけ各 20000 step まで継続 (窓条件未達なら +20000 を 1 回)。床離脱 = B の末尾 10000 step の全保存場で温度床到達 0・他の床/化学種ゲート合格、近傍最低温度も追う (冷点の移動を除く)。(3) 力の感度: 両者 4 量 STEADY・窓条件、Dq = |平均_B − 平均_A| + max(a_A, d_A, 1e-6) + max(a_B, d_B, 1e-6) ≤ C_T・C_T_with_shear・C_L 5e-5、C_M 5e-4 (領域許容の 1/10、この介入への感度が小さいとだけ認定)。**分岐**: 床離脱・加熱側の収支・閉合が成立 → 第 1 仮説 (後流向き面の速度外挿が冷却を維持) を支持し恒久対策の別 plan へ。有効な介入でも加熱応答なし・準定常後も床が残る → 棄却。介入無効・閉合不足・過渡未終了は判定不能。A の FLOOR_STUCK は対照として保持し生産評価に受理しない。設計差への無影響はこの 1 組では証明できない (δ₁.₀ − δ₀.₈ が要る)。やらない: 床ゲート除外・床値変更・`bndFirstOrder`・診断の面/節点指定を生産へ入れる・100 step の昇温で影響なしと結論。**m4_off の壁解像 (2026-10-05、`check_wall_resolution.py`)**: run_1058 で FAIL — 側面の壁 (sidewall・vehicle_side、もともと非解像 = R5m) が y₁⁺ 平均 19–35 で主因。y 法線壁は ramp 0.55 (>1 が 6.9 %)・cowl_in 1.03 (49 %)・cowl_out 0.82 (6.1 %)。**比較: 同じ exact A 形状の m6_on (run_1050) も FAIL** で y 法線壁は ramp 2.05 (57 %)・cowl_in 4.69 (97 %) — m4_off のほうが良い。生産 g3 格子は m6_on でも局所 y₁⁺ ≤ 1 を満たしていない (§4.42 の ramp 0.549・cowl_in 1.234 より悪い; 記録の食い違いは未調査)。**R5h-m10 診断 A/B の結果 (2026-10-05)**: 面 ID = 1535833 (517160 の後流側 +x 面、相手 522113 x 0.12255、A の F_roe +11.0 で 6 面中最大)。(1) 1 更新 (`run_1061` A / `run_1062` B、帳簿): 517160 の res_final roe A −0.662 / B +10.53、内部エネルギー射影 A −0.801 / B +10.745 (加熱側へ)、差は対流段で生じ、相手 522113 の差がちょうど逆符号 (roe B − A +11.19 / −11.19 で閉じる)。1 step 後の T[517160] は A・B とも 50.00002 K (1 更新では温度に出ない)。(2) 100 step (`run_1063` A / `run_1064` B): A は T[517160] 50 K のまま。**B は T[517160] 384.6 K に上がったが、隣の 517199 (y −0.01051、同じ x・z、面 1535832 の相手) が 413 → 50 K に落ち、T ≤ 50 は 1 節点のまま** (517200 も 205 → 115 K)。→ 介入した 1 面の速度外挿が 517160 の冷却を担うことは支持されるが、床からは離脱せず冷点が隣へ移った (事前に除くとした「冷点の移動」)。事前登録の「機序を支持するときだけ 20000 step へ」の扱いは諮問 (中間・判定不能の可能性) → ⑤ 3D 最適化へ接続) | **3D 最適化への接続と MOO 再取得** | 3D runner と `L_sw` を最適化 driver に接続 (`L_sw` はメッシャが既存 station に丸めるので、連続変数にする前に物理 station 化)。対象作動点の評価成立・設計差の格子/領域感度・作動点別 C_M 窓 (実機の値と出典はユーザ判断、生産 YAML は加重平均 `cm_min` −7.0 のみで `cm_window` 未指定) を確かめてから MOO 再取得。2D は形状生成・傾向把握に使う (MOO を 2D で回して 3D は代表点確認、という分担は 3D SST をループに入れるユーザ判断 [本表 L442] と矛盾するので採らない) | — | F |
| 1 | ~~作動点のサイクル値化~~ **実装・検証済 (2026-09-05, §10)**。残 = **MOO の再取得** (第 1 回 run_0040 は起動不安定で破棄、確定レシピ §4.12 で run_0054 を投入済)。**2026-09-09: R7 に従属 (R1–R6 が先)** | `problem_moo_sst_node_cycle3op.yaml` (ext_top・板厚 2e-3・cm_min −7.0)。前提の §4.11 / 1b / §8-6 は全て決着。run_0019 の結論は非物理な作動点が駆動したので破棄 |
| 1c | ~~warm start が plan と実装で食い違っている~~ **実装・検証済 (2026-09-06)**。残 = 本段の早期停止 | §4.7 は「作動点間は warm restart (同一メッシュ = index コピー)」と書いているが、`driver_sern._eval_op` は **3 作動点それぞれを一様 IC から 4 段梯子で立ち上げている**。1 評価 = 3 × 12000 step。案: (A) op 間 warm start (同一メッシュなので `restart_by_index` が使える。ただし NPR が 35.4 ↔ 2.8 と 12 倍違うので相似スケーリングが要る、順序は NPR の近い順)、(B) 本段 CFL を固定 0.5 から風洞チェーン (`runner_axismach.run_staged_ns` の `stages="ramp"`) 方式の段階昇圧へ (風洞は cfl_main 2.0–3.5 に到達している)、(C) 候補間 seeding は `interp_field` の座標一致ノード問題があるので後回し。**codex レビュー反映済**: (A) は素の index コピーが熱力学的に非互換 (同じ roe を別 γ で読むと圧力が (γ_dst−1)/(γ_src−1) 倍ずれる。m6→m10 で 1.24 倍、→m4_off で 2.18 倍) → **相似リマップ + 適応段**の A′ に変更し `warm_from_run` として実装。m10_on のみ m6_on から (NPR 35.4→56.4)、m4_off (2.8) は cold のまま。2 形状で cold と C_T が 5 桁一致・**27–28 % 短縮** (run_0072/0073)。(B) の CFL ランプは**撤回** — 風洞の ramp は「収束済み場から同条件へ restart」の文脈で効くもので、case/45 run_0027 の A/B で **warm start では利得なし**と実測されている。短縮の残りは**本段 (6000 step) を力係数の定常判定で早期停止**するほうが効く |
| 1b | ~~`opt.cm_min` の再設定~~ **決着 (2026-09-05): −7.0** | dv 箱 324 点の MOC スイープ (`case/46/sweep_moc_dv_cycle3op.{py,log}`, 成立 246): 設計点 $C_M$ min −13.8 / p10 −10.7 / med −6.9 / p90 −2.5 / max +1.5、$C_T$ 0.81–0.965。重み付き $C_M$ は設計点の ≈0.7 倍 (m10_on −6.1、m4_off +1.0 で緩む) なので −7.0 で最悪 ~15 % を落とす。旧 −2.5 は候補の 9 割を落とす |
| 2 | ~~ランプ後縁側の外部流ブロック~~ **実装・検証済 (2026-09-05, §4.11)**。残 = 剥離が出る作動点の追加可否 (§8-6) | `mesh.ext_top` + `vehicle_taper` で領域は閉じた。m4_off では形状によらず剥離しない (cowl TE 不足膨張)。剥離制約を効かせるなら NPR ≲ 2 の燃料遮断点 (M∞ ≲ 3) を外挿で足す |
| 3 | ~~**3D SST を通す**~~ ~~**解決 (2026-09-06, run_0082)**~~ **撤回 (2026-09-09, codex M3) → R5**。以下は旧記述: `L_sw` を `L_cowl` からずらす (0.8) だけで全段完走。残 = `L_sw = L_cowl` のまま通す必要は無いという判断の追認と、生産解像度での再確認 | NaN は run_0028 の **x 0.98–1.04, y −0.098〜−0.048, z 0.986–0.996** = カウル後縁 ∩ 側壁後縁の交線 (2 本のせん断層の 3 重点)。**2D で確立したレシピ (§4.12/§4.13) がまだ 3D に移植されていない**: (a) 層流暖機段 (2D の決定打)、(b) 本段 CFL 1.0 → 0.5、(c) 後縁板厚を**側壁後縁にも**入れる、(d) 後縁交線の丸め (2D のランプ角部丸めと同じ発想)。新規開発ではなく移植。これが通らないと**隅 R の効果も 3D の δ\* も測れない** (どちらも本質的に粘性) |
| 4 | **3D 形状を dv に入れる** (R5・R6(e) の後。帰還は別 plan) | ~~補正サロゲート ($\Delta(C_T,C_L,C_M)$ = 3D − 2D を Kriging)~~ **撤回 (2026-09-05, ユーザ指摘)**: あれは 2D で作った形を 3D 相当のスコアで**並べ替えるだけ**で、形が 3D 環境に応答しない (「最適化してる感が無い」)。$C_T$ −4.5 %・揚力符号反転という大きさの効果を外挿で扱うのも無理がある。**正しくは 3D SST を in-loop にする** (§5.1-3 の後)。dv 候補: **$L_{\rm sw}$ (側壁の x 到達距離、`SernMesh3DParams.L_sw` に既存)** = 3D 効果の主因、**側壁後縁のカットバック掃引角** (実機は斜めに切る)、**隅 R**。**「側壁の高さ」は自由変数にならない** — 側壁はカウルからランプまでを塞ぐので高さは局所ダクト高さで決まる。隅 R の効きは $L_{\rm sw}$ より一桁小さいと**推測**しているが、3D SST が通るまで検証できない |
| 5 | 亜音速外部流のチャンバー構成 | 静止・地上試験用。case/23 方式。超音速外部流なら不要 |
| 6 | カウル輪郭の設計 | 今は直線固定 (Lv & Xu 2021 はカウルにも最適化理論) |
| 7 | 音速給気のスロート接続 (`sonic_throat`) | §4.10-3 より**対象は $M_\infty 4$ powered のみ**。Hall + kernel MOC で $M$ 1.3–1.5 まで対称内部ノズル |
| 8 | 非一様入口 (回転流 MOC)・凍結 $\gamma$ → 有限速度化学 | 入口分布 BC は評価側に既にある |
| 8b | ~~`L_ramp_max` が probe と設計で不一致~~ **決着 (2026-09-13, R6(d))**: 最終輪郭で再検査 (`driver_sern._check_design`, 許容 `opt.l_ramp_tol`) | 成立性判定 (`_design_probe`) は粗い kernel (`nj_moc_probe` 151, `dx` 4e-3)、実際の設計は細かい kernel (301, 2e-3) を使うため、probe が 19.9 と判定した点が設計で 20.14 になる (run_0071 doe_014)。上限を ~1 % 超える点が通る。probe を細かくするか判定に余裕を持たせる |
| 2b | **m10_on の発散を止める (最優先の技術課題)**。上面線は **後縁接線 Hermite + くさび 3°** に作り替える (§4.11、曲率 1.02 → 0.16)。テーパ開始 0.65 L に根拠が無い件も同時に解消 | §4.11 の真因 (P 床 → 負密度 → ω) に対し `implicitRelax` は無効と判明。次に試す順: (a) `pMin` を p∞ の数 % (50 Pa) へ — **未測定**、(b) リミッタを Barth に / mid 段を 1 次に留める、(c) 機体上面の第一セル細分化、(d) 後端を**丸めた肩 + 有限ベース**に作り替える (ナイフエッジと 90° 角の両方を避ける)。`evaluate.implicit_relax` / `evaluate.p_min` は配管済み (commit 9a53bb9a) |
| 9 | CFD 領域の縮小 | `mesh.bot_depth` 3H → 超音速外部流なら 1H 程度 |
| 10 | **Shyne の式 11–15 と突き合わせ、厳密 Rao 最適の位置を確かめる** | まず NASA TM-100955 / TM-103175 (§4.3) の平面版乗数条件と `rao_planar` を照合する。次に `rao_planar.lip_residual` ($\theta_c-\theta_{\rm lip}$、0 が Rao 最適) は今は診断値で、dv には課していない。残差 0 になる $(M_c,f)$ を箱の中で求め、**MOO のパレートがその点を含むか**を見る。理論の最適点は推力側の端点として出てくるはずで、出てこなければ平面縮約か実装のどこかが違うという判定になる (§4.3 の依拠箇所の検算) |

## 6. 検証

- **単体 (S1)**: (i) 平面単純波則 ($\theta+\nu$ 保存) 1e-10、(ii) 対称極限 (カウル = ランプの鏡像) で
  平面最小長ノズルの既知解 ($\theta_{\max}=\nu(M_e)/2$、Argrow–Emanuel) を再現、(iii) $p_e/p_a=1$ で
  制御面が一様平行流、(iv) 制御面上の推力積分 = 出口面の推力積分 (保存則、<0.1%)、(v) 逆設計で
  張った目標 C⁻ を前進 MOC で再計算した状態が目標と $M$ 0.1% / $\theta$ 0.05° 以内、(vi) 格子収束。
- **CFD (S4)**: 設計点 Euler (`visc: 0.0`, `thermCond: 0.0`) で $C_T$ が MOC 値と 1% 以内、RANS で
  全作動点 `check_convergence.py` PASS (低 NPR は `check_quasisteady.py` の統計)。NASA 基準形で
  「カウル短縮 → $C_T$ 低下・頭下げ $C_M$」「カウル角 6° 近傍で $C_M$ 最良」の傾向再現。
- **メッシュ**: 各 run で `check_mesh_quality.py` VERDICT PASS を README に記録。
- **δ\* (S5)**: 補正前後で $C_T$ の変化が符号・桁で妥当 (壁摩擦込みで −1〜−3% 程度) であること。
- **3D (S7)**: 2D vs 3D の $C_T$ 差 < 2% を「2D 設計で足りる」判定基準とする。
  → **2026-09-05 時点で不成立**: run_0027 (3D Euler, 側壁あり) vs run_0029 (同条件 2D) で
  $C_T$ 0.932 vs 0.976 (**−4.5 %**)、$C_L$ −0.204 vs +0.004 (**符号反転**)、$C_M$ +1.07 vs −0.001。
  **2D チェーンが精密に解いている量 (摩擦 1 %、δ\* が $C_L$ で 5 %) より 3D で落ちている量が大きい**。
  したがって「2D 設計 + 3D 確認」では閉じず、§5.1-4 の手当てが要る。
  → **2026-09-09 (codex M5 採用)**: この −4.5 % には幅外の機体下面の集計が混入している (除くと −1.9 %)。
  両 run とも NOT CONVERGED なので 2 % 判定は §5.1 R2 (集計分離) と R5 の後に再計算する。
  → **2026-09-13 (R2 完了後の再計算, Euler)**: run_0092 (3D, `vehicle` タグ, 6000 step) vs run_0093 (2D, 8000 step)、両 run とも力係数 STEADY・
  残差 2.4–3.2 桁でプラトー (NOT CONVERGED stalled)。**ノズル $C_T$ 0.9578 vs 0.9762 = −1.9 %** (2 % 内)、$C_L$ −0.100 vs +0.004、$C_M$ +0.72 vs −0.00。
  推力は「2D 設計 + 3D 確認」で足りるが揚力・モーメントは 3D 効果が支配的 = 3D 補正表が要る。SST での判定は R5 の後。

### 6.1 レビュー記録 (codex)

[`AGENTS.md`](../../AGENTS.md) 「codex レビュー」の記録表 (`solver_density_cuda/tools/codex_review.py`、手順 [`procedures/codex-review.md`](../../procedures/codex-review.md))。
本 plan は 2026-09-04 起票で S0–S6 実装済みのため `plan` 段は事後レビュー (残作業表の優先 3 件に重点)。

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-09 | [2026-09-09-tooling-nozzle-sern-chain-plan.md](../../notes/reviews/2026-09-09-tooling-nozzle-sern-chain-plan.md) | NO-GO, C2/M7/m0 | **全件採用 (ユーザ決定 2026-09-09) → §5.1 R1–R7 に展開、§4.13-3 撤回、§8-7/§8-10 決着**。当日のスポット検証: C1 の `steadiness([1,1,1,NaN])` → `STEADY` は再現 (実バグ)、M3 の run_0083/0088 `DIVERGED (NaN/Inf)` は `check_convergence.py` で再確認。C2 (外気と排気で同じ γ,R → 外部動圧 −15 %) と M5 (3D 推力 −4.5 % に幅外機体下面の集計混入) は未検証 |
| diagnose | 2026-10-05 | [2026-10-05-sern-next-priority-diagnose.md](../../notes/reviews/2026-10-05-sern-next-priority-diagnose.md) (brief [2026-10-05-sern-next-priority.md](../../notes/reviews/briefs/2026-10-05-sern-next-priority.md)) | Major 5 / Minor 1 | **全件採用**: R7 を R7a (3D `L_sw` 2 水準 × 3 反復の判別試験、最優先) と R7b (3D 最適化接続 → MOO) に分割、R6 (d) の数値許容・検出目標を R7a 前に固定し実機 C_M 窓は R7b 前、R5h は診断完了・恒久対策保留、W2 は実施済み、farfield の限定受理を R4 全体の完了と読まない、`L_sw` は既存 station に丸められるので二水準は実位置で記録 |
| diagnose | 2026-10-05 | [2026-10-05-sern-r7a-result-diagnose.md](../../notes/reviews/2026-10-05-sern-r7a-result-diagnose.md) (brief [2026-10-05-sern-r7a-result.md](../../notes/reviews/briefs/2026-10-05-sern-r7a-result.md)) | Major 5 | **全件採用**: 分岐 A を限定付きで保持、写像だけを変えた対照 1 本 (双子の内外を保存)、ΔC_M は設計応答で §8 と比べない、C_M は圧力・ノズル面のみで実機窓と比べるには定義をそろえる、MOO 本番の前に G_B − G_A (§5.1 R7a) |
| diagnose | 2026-10-05 | [2026-10-05-sern-r7a-g4-diagnose.md](../../notes/reviews/2026-10-05-sern-r7a-g4-diagnose.md) (brief [2026-10-05-sern-r7a-g4.md](../../notes/reviews/briefs/2026-10-05-sern-r7a-g4.md)) | Major 5 | **全件採用**: 旧 g4 は判定不能で確定、B_g4 の固定長 +20000 を別試験として事前登録 (§5.1 R7a)、C_M を外す救済は却下、「格子依存 1e-7 級」撤回、写像対照の表現を限定、実機トリムは保留で `cm_min` は暫定、R7b の順序 ①–⑤ |
| diagnose | 2026-10-05 | [2026-10-05-sern-r7b-m10-te-floor-diagnose.md](../../notes/reviews/2026-10-05-sern-r7b-m10-te-floor-diagnose.md) (brief [2026-10-05-sern-r7b-m10-te-floor.md](../../notes/reviews/briefs/2026-10-05-sern-r7b-m10-te-floor.md)) | Major 5 | **全件採用**: R5h 再開・診断 A/B (`FORGE_DIAG_FACE_VEL_CELL` の有無だけ)、床ゲート除外・有限厚後縁の先行・m10_on 抜きの先行は却下、先行は m4_off の壁解像だけ (§5.1 R7b 4-2) |
| diagnose | 2026-10-05 | [2026-10-05-sern-r7b-step4-diagnose.md](../../notes/reviews/2026-10-05-sern-r7b-step4-diagnose.md) (brief [2026-10-05-sern-r7b-step4.md](../../notes/reviews/briefs/2026-10-05-sern-r7b-step4.md)) | Major 5 | **全件採用**: 2D 再診断を前提にしない、exact 形状で m6_on 二水準を取り直す (4-1)、生産 YAML 名だけでは同条件にならない (時間/space・基準値を継承)、3 作動点 × A/B で感度、0.002 は各作動点の必須量ではなく加重目的で判断 (§5.1 R7b ④) |

## 7. 影響範囲

- 新規: `design/forge_design/geometry/{sern_geometry,moc_sern,rao_planar}.py`, `meshing/mesh_sern.py`,
  `evaluate/runner_sern.py`, `metrics/sern_forces.py`, `design/tests/run_sern_moc_tests.py`, `case/46.sern_design/`
- 変更: `probdef.py` (`KNOWN_TYPES` に `sern_2d`), `opt/driver.py` / `opt/moo.py` (多作動点・多目的), `design/CAPABILITIES.md`
- forge 本体: 変更なし (2D 平面 RANS・超音速入口・外部流は既存機能)。不足が出たら別 plan
- 文書: `methods/design/overview.md` 「SERN チェーン」節を実装に同期、親 plan §4.6 ⑤ / §5 Phase 4

### 4.12 SST の起動安定性 (2026-09-05 確定)

サイクル値作動点では**カウル後縁のせん断層**で `roOmega` が発散する。排気と外部流の密度・温度比が大きい
(m10_on: ρ 0.067 vs 0.012、T 2222 vs 228 K) ためで、3 作動点すべてで同じ場所 ($x \approx L_{\rm cowl}$)。効かなかった対策と効いた対策:

| 対策 | 結果 |
| --- | --- |
| カウル板厚 | **不可**。作動点で要求が逆転する (m4_off は 2e-3 で通り 5e-3 で落ち、m10_on は 5e-3 で通り 2e-3 で落ちる)。MOO では成立しない |
| soft 段 CFL 0.5 → 0.2 | 入口面∩カウル外面の角 (M∞10 で顕在化、run_0041) には有効。後縁には無効 (step 1719 に遅れるだけ、run_0042) |
| **層流暖機段**を soft の前に入れる | **有効**。平均場が立ってから乱流方程式を入れる。3 作動点とも通り、力係数は暖機なしの成功例と一致 (run_0048) |
| 本段 CFL 1.0 → 0.5 | **必要**。1.0 は限界的で、同一メッシュ・同一 config でも run により発散する (run_0048 成功 / run_0049 失敗)。0.5 で m10_on 2 回とも同一値 (run_0050/0051) |
| plan | 2026-09-19 | [2026-09-19-tooling-nozzle-sern-chain-plan.md](../../notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan.md) | NO-GO, C1/M4/m0 | **全件採用** (§4.14.1)。C1 → R4c、M2 → R5b、M4 → R4d を残作業表に起票。M1 は §4.14-5 を撤回し自分で再検証 (T=6000 K は 1 ノードで既に NaN、床は T=50 K 側)。M3 は即修正 (`runner_sern3d` の壁を `wall_bcond_line` に接続) したが run_0119 は step 76 で発散し単体では直らない |
| plan | 2026-09-19 | [2026-09-19-tooling-nozzle-sern-chain-plan-2.md](../../notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan-2.md) | NO-GO, C2/M5/m1 | **全件採用** (§4.15.1)。C1 → §4.15-1 を撤回し R5b に統合、C2 → R4e、M1 → R4d、M3 → R4f を起票。M2 は取り下げ、M4/M5 は `mirror_symmetry.py` に実装済。run_0122 の結果は C2 の壁 (5.72 MPa) の影響下にあるので**評価に使わない** |
| plan | 2026-09-19 | [2026-09-19-tooling-nozzle-sern-chain-plan-3.md](../../notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan-3.md) | GO-with-changes, C0/M5/m1 | **全件採用** (§4.15.3)。M1 → 後縁でブロックを切り替えない形に設計変更、M2 → `run_0034` の診断を撤回 (実測は 2 次・`cfl_pseudo=4` で「soft 段・90° 二重 slip」ではない、根因未確定に修正)、M3 → 物理寸法を格子から独立させる、M4 → BC と帳簿を案 (d) の生産仕様に内包 (R4f を前倒し)、M5 → 実装前に固定する検証 5 段を記載、m1 → §4.8 の旧 4 目的記述を履歴と明示 |

確定レシピ = 板厚 5e-3 + 層流暖機 (cfl 0.2, 2000) + SST soft (cfl 0.2, 2000) + SST 本段 (cfl 0.5, 6000)。
検証値 (dv 既定): m6_on C_T 0.9089 (摩擦込み 0.9018) / m10_on 0.9246 (0.9175) / m4_off 0.9759 (0.9701)、全て STEADY。

### 4.13 評価の受理方針と再試行梯子 (2026-09-05 決定)

dv 箱を動かすと**作動点 × 形状の組み合わせごとに固い場所が変わる** (ランプ膨張角部 → §4.12 の丸めで解決、
カウル後縁 → m4_off × 短カウルで残存)。全点に効く単一設定は見つからなかったので、**梯子**にする。

1. **標準レシピ** (§4.12 + `mesh.ramp_fillet: 0.1`): 層流暖機 2000 @cfl 0.2 → SST soft 2000 @0.2 → SST mid (2 次) 2000 @0.2 → 本段 6000 @0.5
2. 落ちたら **緩レシピ**で 1 回だけ再試行: 各段 step 倍・cfl 半分 (4000 @0.1 ×3 → 本段 8000 @0.25)。別 run ディレクトリ (`*_retry`) に作る
3. それでも `rc != 0` の場合でも、**力係数の履歴が 4 スナップショット以上あり 3 成分とも `STEADY`** なら採用し、
   台帳に `degraded: true` を記録する。カウル後縁の破綻は**力積分が収束しきった後に起きる**ため
   (run_0067: C_T 0.95766→0.95741, trend 3.6e-5 の 6 点で `STEADY` の後に発散)、積分値は使える。
4. 上を満たさなければ FAIL (最適化側には INFEASIBLE として返る)

**撤回 (2026-09-09, codex C1 採用 → §5.1 R1)**: 上の 3. (「`rc != 0` でも力係数が `STEADY` なら採用」) は**廃止**する。
`steadiness` が NaN 系列を `STEADY` と判定する欠陥があり、標準レシピで条件を満たすと再試行も省略されていたため、
発散・物性誤差が設計性能として取り込まれていた。発散・非有限値を含む評価はサロゲート学習と Pareto から除外し、
保存場の有限値・全残差・実目的量 (`C_T_with_shear`) の定常性を独立ゲートにする。以下の `degraded` の記述は R1 で置換される。

**R1 実装 (2026-09-13)**: 採用条件は **`rc == 0` かつ §4.7 の 4 ゲート全 PASS** のみ (`driver_sern._eval_op`)。標準レシピがゲート不合格なら
緩レシピで 1 回再試行し、それで通れば `degraded: true` (= 再試行で通った評価) を**台帳と `pareto.json` の両方**に残す。
分類: `INFEASIBLE` (物理: `DESIGN` 逆設計不成立 / `L_RAMP_MAX` 最終輪郭超過 / `CM_WINDOW`) と `FAIL` (数値: `DIVERGED` / `RESIDUAL_RISING` /
`NOT_CONVERGED` / `UNSTEADY` / `NO_FORCES` / `ERROR`) を分け、いずれもサロゲート学習・Pareto から外す。既存キャンペーンは
`driver_sern --rejudge <out>` で CFD を回さず再判定できる (元 run は読むだけ)。**旧 `steadiness` の NaN→STEADY は再現・修正済**
(`[1,1,1,NaN]` → `NONFINITE`、テスト `design/tests/run_sern_gates_tests.py`)。なお handover の「node の `rms_roOmega` が本段で 3e18 一定
(壁ノードピン留めの診断値)」は現行バイナリでは解消している (`nodeWallDirichlet_d.cu` が Dirichlet ノード残差を除外; run_0069 の m4_off は
1e1 台で推移)。run_0069 doe_014 m4_off の 5.3e18 は step 5461 からの**本物の ω 発散**で、力係数が STEADY のまま起きた = 旧ゲートの穴。

## 8. 未確定事項 (ユーザ確認)

1. ~~**入口の既定**~~ **決着 (2026-09-05, §4.10-3)**: 超音速給気 (燃焼器出口 $M_{\rm in}$ 指定) を既定で続ける。TM X-71972 の作動点では音速スロートが要るのは $M_\infty 4$ powered ($M_3 = 1.12$、ラムジェット) の 1 点だけで、残り 5 点は $M_3 \ge 1.67$ の超音速。`sonic_throat` は §5.1-7 に後回しでよい。
2. ~~**作動点セット**~~ **決着 (2026-09-05, §4.10)**: NASA TM X-71972 TABLE 1 (定動圧 1500 psf 経路) をアンカーにサイクル値で定義する。既定 = 設計点 $M_\infty 6$ powered、作動点 M6 powered ($w$ 0.5) / M10 powered (0.3) / M4 power-off (0.2)。低 NPR は「低速飛行」ではなく「燃料遮断」で作る。
3. **モーメント基準点と目標**: $C_M$ を最小化するのか目標値に合わせるのか。**部分決着 (2026-09-13, R6(d))**: 目的には入れず制約 (加重平均窓 + 作動点別窓 `opt.cm_window`) とする。許容値の根拠 (トリム能力) は案件で決める。
4. **隅 R・側壁**: S7 で spec として評価するだけ。設計変数にするかは S7 の結果で判断。
5. **外部流を設計段階に入れるか**: 入れない (NASA 流) を既定。
6. ~~**剥離が出る作動点をセットに入れるか**~~ **決着 (2026-09-05, ユーザ判断): (b) 入れない**。TM X-71972 の 3 作動点で
   MOO を回し、剥離制約は持たず `sep_frac_ramp` を台帳記録に留める。根拠: 剥離しないのは**作動点の性質** (最低 NPR 2.8 でも
   cowl TE が 1.35–1.59 p∞ の不足膨張) であって形状の選び方ではなく (L_cowl 1.2 と dv 上限 2.5 の両方で `sep_frac` 0、
   §4.11 結果)、制約として設計を弁別する情報量が無い。NPR ≲ 2 の点を入れるには TM アンカー外の外挿が要り、
   そこまでして得るものが無いと判断した。降下・アボート側の要件が案件で立ったときに再検討する。

7. **外部流のガスが排気ガスのまま** (2026-09-05 発見、未決): `gas_states` は単一 CPG の $(\gamma, R)$ を排気にも外部流にも使う。
   §4.10 で作動点ごとに排気ガス (燃焼生成物) を入れた結果、**外部流 (空気) が燃焼生成物の $\gamma, R$ で計算される**。
   $M_\infty, p_\infty, T_\infty$ は指定どおりだが $\rho_\infty, u_\infty$ がずれる: m10_on で $\gamma R$ = 479.5 (空気 401.8) →
   $u_\infty$ +9.3 %、$\rho_\infty$ −27 %。m6_on は $\gamma R$ = 402.5 で偶然ほぼ一致 (誤差 0.1 %)、m4_off は空気そのもので厳密。
   影響: 外部流の斜め衝撃角・Prandtl–Meyer が変わり、カウル TE のせん断層が実際より強くなる (m10_on の数値不安定の一因の疑い)。
   選択肢: (a) 現状維持 (ノズル性能は排気ガスで決まるので設計目的には妥当) + 制限として明記、(b) 外部流を空気に合わせる
   (排気の膨張が狂うので不可)、(c) `gas.model: semiperfect` の多成分化 (§4.5) — 別 plan 規模。~~**推奨 = (a)** だが、
   m10_on の安定性が (c) を要求するなら再考する。~~ **決着 (2026-09-09, codex C2 採用 → §5.1 R3): (c)**。
   外部動圧が m6_on でも −15.5 % (密度 −15.7 %) で (a) は外力・せん断層の評価を正当化しない。既存の多成分 TP を再利用する。

8. **`L_cowl` の上限 2.5 が狭すぎる** (2026-09-05 発見、要判断): 現行 dv は 0.6–2.5 H だが、**NASA TM X-71972 の基準形は
   カウル長 3.12 H** で箱の外にある。C_T は上限で飽和しておらず単調増加で、最適化が上限に張り付く:

   | 出典 | カウル長 [H] | $C_T$ |
   | --- | --- | --- |
   | S4 run_0003 (NASA 平板, M∞10) | 2.0 | 0.9605 |
   | S4 run_0004 (**NASA 基準形**) | 3.12 | 0.9704 |
   | S4 run_0005 | 4.5 | 0.9720 (ほぼ飽和) |
   | MOC dv スイープ (中央値) | 0.60 / 1.55 / 2.50 | 0.8749 / 0.9060 / 0.9161 |

   スイープの $C_T$ 上位 8 点は**全て `L_cowl` = 2.5 (上限)**、キャンペーン run_0069 の高 $C_T$ 点も
   doe_004 (2.30) / doe_014 (2.06)。→ **上限を 4.5 H まで広げるべき** (S4 で飽和が確認できている値)。
   副次的に、数値的に固い短カウル域 (`L_cowl` ≲ 0.7 で m4_off / m10_on が落ちる、§4.13) から離れる利点もある。
   **決着 (2026-09-05, ユーザ判断): 打ち切って上限 4.5 H で回し直す**。run_0069 は 22 評価 / PASS 14 で打ち切り
   (最良 doe_004 C_T_w 0.9434 @ Lc 2.30)。広い箱の MOC スイープでは C_M が Lc とともに緩む (中央値 −7.25 @0.6 →
   −2.02 @4.5) ので `cm_min` −7.0 は据え置き。`x_max_kernel` は 14 → 20。

9. **評価を Euler + 摩擦相関で代用できるか** (2026-09-05 提起): 剥離が出ないと確定した (§4.11 結果) ので、
   NS を回す実質的な理由は**摩擦が目的関数と結合していること**に絞られる。run_0071 の m6_on 13 点で
   摩擦の寄与は $L_{\rm ramp}$ に明確に相関する:

   | $L_{\rm ramp}$ | 5.77 | 9.25 | 13.92 | 18.26 |
   | --- | --- | --- | --- | --- |
   | 摩擦 / $C_T$ | 0.53 % | 0.97 % | 1.04 % | **1.13 %** |

   ランプを 5.8 → 18.3 H で摩擦は **2.3 倍**。Euler だけだと長いランプの罰則が幾何長さだけになり、パレートが
   長いランプ側へ偏る。→ Euler (1 評価 ≈ 3 分、**10 倍速**) + 摩擦の相関モデルで代用できる可能性がある。
   ただし相関を作るには NS を数十本回す必要があり、現行規模 (1 評価 3 分半) では割に合わない。
   **3D を in-loop に入れる (§5.1-4 第 2 段 (a)) ときに再検討する**。

10. **ロードマップを組み替えるか** (2026-09-05 提起、**要ユーザ確認 / plan 未反映**): 上の 3 点 (§4.6 の 3D では帰還必須、
    §4.3 の Rao 最適性は 3D で失われる、§5.1-4 の補正サロゲート撤回) から、本筋は次の順が正しいと考えている:

    1. **3D SST を通す** (§5.1-3。2D レシピの移植)
    2. **3D で δ\* と横方向実効面積を測る** (帰還の入力を作る)
    3. **dv を (MOC 5 個 + $L_{\rm sw}$ + 側壁カットバック掃引 + 隅 R) に拡張**し 3D SST を in-loop
    4. **設計への帰還を入れる** (3D の排除効果を MOC 境界条件へ)。ここまでで初めて「3D で最適化した」と言える

    これに伴い、走行中の 2D キャンペーン (run_0071) は**パラメータ化の素性と 2D 傾向を掴む役**に格下げして完走させる。
    ~~§5.1 の優先順の組み替え (3D SST を最優先へ) はユーザ判断待ちで、**まだ §5.1 の順序には反映していない**。~~
    **決着 (2026-09-09, codex レビュー採用 → §5.1 R1–R7)**: 順序は「評価器の成立確認 (ゲート・集計・物性・領域・3D SST 再現)
    → 小規模探索 → MOO 再取得」に組み替えた。上の 4. (MOC への帰還必須) は撤回 (R6(e))。

11. **壁の熱境界 (2026-09-13 決定・要ユーザ確認)**: 生産 YAML (`problem_moo_frozen_tp_cycle3op.yaml`) の既定を **等温壁 $T_w$ = 1000 K** にした。根拠: 断熱だと M∞10 の外気側カウル板下面が
    回復温度 4600 K (実機では成立しない) になり、node 壁列の T 市松で隣接ノードが T 床 50 K に落ちて残差が 5e-5 でプラトーする (run_0095–0099)。
    等温 1000 K では残差が 3.5 桁落ち (run_0103)、摩擦込み $C_T$ は −0.4 % (摩擦 +30 %)、$C_M$ +0.11。$T_w$ = 1000 K は冷却構造の代表値として置いた暫定値で、
    案件の壁温 (材料・冷却) が決まれば差し替える。作動点ごとに $T_w$ を変える必要があれば `operating_points[].wall_thermal` を追加する (未実装)。

## 9. 完了条件

- [ ] `methods/design/overview.md` の SERN 節を実装と同期
- [ ] S0–S7 完了、§6 の判定を満たす (S8 は別 plan)
- [ ] `design/CAPABILITIES.md` の `sern_2d` を ✅ に更新 (検証ケース = case/46)
- [ ] `status: done` にして §10 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` 同期

## 10. 変更ログ

- `2026-09-13` — **R2 完了・R3 実装** (§4.5/§4.10/§5.1)。R2: `mesh_sern3d` に `vehicle` タグ + `W_vehicle`、`forces3d` はノズル面のみ (vehicle 別枠)、
  `metrics/sern_momentum.py` (BCONDS の運動量収支; 2D 閉じ残差 0.1 % [Euler] / 1.1 % [SST], 3D 1.6 %; 壁力の帳簿一致)、run_0092/0093 で
  3D/2D ノズル $C_T$ −1.9 % (§6)。R3: `gas/frozen.py` + `runner_sern` の frozen_tp 配管 (擬似種 EXH/AIR, thermoHrefTemp 298.15, IC の roe を同じ datum で、
  理想推力は凍結等エントロピー膨張)、`tmx_operating_points.py` が組成と凍結 $M_3$ を出力、生産 YAML `problem_moo_frozen_tp_cycle3op.yaml`。
  単体 40 項目 ALL PASS (Ar で CPG 一致・MW = CEA・q∞ 71.9 kPa)。実機 run_0094 (doe_001 dv): m6_on は CPG 生産レシピ (暖機 cfl 0.2) で
  ランプ後縁ノードの EXH∩AIR 接触で暖機 step 568 に NaN → 緩レシピで完走 GATES PASS ($C_T(p)$ 0.9520 = MOC 0.9470 +0.5 %, 摩擦込み 0.9447,
  $C_M$ −5.11)。m10_on は `warm_from_run` frozen 経路のバグ (作動点未適用の YAML から組成) → 修正。frozen の標準レシピは暖機 cfl 0.1 × 4000
  (run_0095 で 3 作動点を再実行、結果は下)。
- `2026-09-13` — **R3 の実機確認 (run_0095–0099)**: m6_on は暖機 cfl 0.1 で標準レシピ PASS (摩擦込み $C_T$ 0.9447, $C_M$ −5.11)。m10_on は本段 6000 では
  RESIDUAL_RISING (残差が最小から 3–4 倍にリバウンド中、力係数は定常) → 12000 で飽和して PASS (run_0097, 0.93528; warm/cold/緩の 3 経路で 5 桁一致)。
  run_0098 (本段 12000, 3 作動点): PASS (degraded) $C_{T,w}$ 0.9415 / $C_{M,w}$ −3.00; m10_on warm 標準は mid 段 cfl 0.2 で NaN (限界状態) →
  frozen 生産レシピを soft/mid/暖機 cfl 0.1 × 4000 + 本段 12000 に確定 (**run_0099: retry 無しで 3 作動点 PASS**, $C_{T,w}$ 0.94152, 463 s/評価)。m4_off 0.9430 / $C_M$ +2.96。近壁の 2 観察 (入口∩壁角 1.75 $p_{\rm in}$、
  M∞10 板下面の冷点) は §5.1 R4b に起票。R3 の**設計への影響**: 凍結 γ 1.247 / $M_{\rm in}$ 1.631 で同じ dv の形状が変わる (doe_001: $L_{\rm ramp}$ 11.89 → 10.64 H)
  ので、run_0069 系の CPG 台帳と frozen 台帳は形状も違う = 比較不可、MOO は frozen YAML で再取得 (R7)。

- `2026-09-13` — **R1 評価ゲート修復を実装・検証** (§4.7/§4.13, codex C1 採用): `metrics/sern_gates.py` (rc / 保存場有限・正値 / 全残差 NaN・rising /
  実目的量 + $C_T,C_L,C_M$ の STEADY を独立ゲート化)、`sern_forces.steadiness` を正式ツール `check_quasisteady.classify_series` に委譲
  (NaN→STEADY バグ修正、`NONFINITE` 語彙を追加)、`check_quasisteady.py --series-csv/--series-cols` (CSV 系列モード) と `force_history.csv`
  (R6(b))、`driver_sern` は発散評価の採用 (旧 §4.13-3) を撤回し PASS/INFEASIBLE/FAIL + fail_class に分類、`degraded` を `pareto.json` にも保存、
  `--rejudge` で既存キャンペーンを CFD 無しで再判定、`L_ramp_max` を最終輪郭で再検査、作動点別 C_M 窓 `opt.cm_window` (R6(d))。
  R2 の集計側: `forces3d` をノズル力 / 機体力に分離 (`mesh3d.W_vehicle`)。単体テスト `design/tests/run_sern_gates_tests.py` 62 項目 ALL PASS、
  既存 `run_sern_moc_tests` 19 / `run_sern_mesh_tests` 29 ALL PASS。**再判定** (`case/46/run_0090_rejudge_r1_gates/`): run_0069 は
  PASS 14 → **10** (HV 2.6835 → 2.6450): 除外 4 件 = 残差 rising 3 (doe_004/015 m6_on は本段末尾で $\rho,\rho u,\rho e$ 残差が 2 倍リバウンド、
  doe_014 m4_off は step 5461 から $\rho\omega$ 残差 1e1 → 5e18 の**本物の発散**、いずれも力係数は STEADY だった) + $C_M$ TRANSIENT-UNSETTLED 1 (doe_016);
  run_0054 は変わらず 1/1 DIVERGED。3D run_0027 再集計はノズル 0.9578 / 機体 −0.0253 / 総計 0.9325 で codex 検算と一致 (2D 0.9762 比 −1.9 %、両 run とも
  NOT CONVERGED のまま = 2 % 判定は R5 の後)。実機スモーク run_0091 (run_0069 doe_001 の再評価、現行 driver): 結果は case README。
  **未実施**: R2 のメッシャタグ分離・運動量収支、R3–R5、R6(a)(c)(e)。

- `2026-09-09` — codex plan 段レビュー (NO-GO, C2/M7) を §6.1 に記録し、**全件採用 (ユーザ決定)**: §5.1 を R1–R7 (評価ゲート修復 / 力の集計分離 / 凍結 TP 化 + 外気空気 / 外部領域と領域独立性 / 3D SST 再現 / 検証・問題定義の整合 / 小規模探索後に MOO) に組み替え、§4.13-3 (発散評価の採用) 撤回、§5.1-3 「3D SST 解決」撤回、§8-7 は (c)、§8-10 は帰還必須論を撤回して決着。

- `2026-09-04` — 初稿。調査ノート [`sern-design-method-survey.md`](../../notes/investigations/sern-design-method-survey.md) の推奨 (§4) を計画化。親 plan §4.6 ⑤ の壁圧 Bézier dv・局所帰還・3D FFD in-loop を撤回し本計画に差し替え。branch `feature/sern-design`。
- `2026-09-04` — **S0–S1 実装** (`geometry/{sern_geometry,rao_planar,moc_sern}.py`, `probdef.KNOWN_TYPES` に `sern_2d`,
  `tests/run_sern_moc_tests.py` 全 PASS)。平面 MOC は逆 (格子) 法 + 角部扇の解析閉包 + **TE 扇と自由境界反射の
  C⁺ レイ束** (格子だけでは扇が最初の station で 1 セルに潰れ伝播しない・楔域が古い値のまま残る、の 2 件を実測して
  導入)。閉包判定の罠: 「K⁻ が入口値のまま」は TE 扇の下流でも真になるので TE 先頭レイより上流に限定した。
  §6 解析検算: 一様流機械精度 / ランプ扇 M,θ 誤差 0 (解析閉包) / カウル反射後 M(ν_in+2θ_r) 一致 /
  **対称 MLN 極限 (M_in 1.5→M_e 3)**: θ_c=0.0000°, 出口高さ = 等エントロピー面積比 (0.01%), c–e 流量 = f (0.01%),
  K⁻ 一定 1e-12, 上半分総推力 = 出口運動量流束 (0.03%) / カウル付き SERN (M_in 2.5, 15°/5°, L_cowl 1, p_ext 0.05):
  格子倍で C_T 変化 0.001%・L_ramp 0.04%。dv 掃引 33 点 (f 0.35–0.6, M_c 3.2–4.4): 等長候補間で C_T 最大点が
  縁条件残差最小側 — 粗いので最適性の確証は別途 (fine sweep か等長拘束の直接最適化)。
  **速度**: kernel march 14 s (nj 301, dx 2e-3, x 9H)、逆設計 1 本 <1 s。図: scratchpad `sern/` (artifact 化)。
  **未対応**: semi-perfect ガス (圧力比写像)、非一様入口、p_ext > p_TE の衝撃。
- `2026-09-04` — **S2–S3 実装 + S4(a) 合格**。`meshing/mesh_sern.py` (2 バンド構造 quad、カウルはスリット =
  中間線ノードを上下 2 重、**TE 点は共有** — 重複させると TE 直後に幅 0 の隙間 = 未タグ境界辺が出る、
  `tests/run_sern_mesh_tests.py` で検出・修正)、`evaluate/runner_sern.py` (逆設計→メッシュ→品質ゲート→領域別一様 IC→
  段階起動 soft 1 次 cfl0.5 3000 step → 本段 2 次 cfl4。staging の場移植は `interp_field.py`)、
  `metrics/sern_forces.py` (壁面出力 `res_<name>_<id>_<step>.h5` の CONNE 面ごとに (p−p_a)·n を積分。cell = 面値 / node = 節点平均)。
  **S4(a)**: case/46 run_0002 (cell Euler slip, 56k セル, 品質 PASS) で C_T 0.9660 vs MOC 0.9666 (−0.05 %)、
  C_L 0.1427 vs 0.1392、C_M −0.939 vs −0.946。壁圧分布はランプ・カウル内面とも MOC に重なる
  (x/H 6–7 の後縁扇反射位置に小差)。残差は 1.4 桁プラトー (`NOT CONVERGED`) だが力係数は STEADY (1e-5)。
  `geometry.mode: straight` (平板ランプを切るだけ、NASA 照合用) を追加。S4(b) は run_0003–0007 (カウル長 2/3.12/4.5H、
  カウル角 3/6/12°, M∞10, γ1.3) を投入。
- `2026-09-04` — **S4(b) 合格** (case/46 run_0003–0007, `geometry.mode: straight`, cell Euler, M∞10/γ1.3): 内面の力は
  forge と MOC が 5 形状すべてで C_T +0.001 以内・C_M 0.05 以内で一致。カウル長 2.0/3.12/4.5H で C_T 0.9605/0.9704/0.9720
  (短縮で大きく減、延長の利得は小 = NASA)、前方基準 (−20H) の C_M −1.31/−0.68/−0.36 (短いカウルで大きな頭下げ = NASA)。
  カウル角 3/6/12° で C_T 0.9786/0.9704/0.9314 (単調減 = NASA)。**カウル外面の衝撃圧は MOC に無い寄与** (12° で C_T −0.035) で、
  評価器で必ず取る。C_M の傾向は基準点に依存するため、問題定義の `moment_ref` は機体 CG を必ず与える (§8-3 の回答)。
  詳細は case/46 README。
- `2026-09-04` — **RANS 化・S5・S6 (Euler) 完了** (case/46 run_0008–0013, run_0010_moo_euler_2op):
  (i) node + SST は本段でカウル角部直下流の内面側から ω が発散 (cfl 4/1 とも、run_0008/0009) → **未解決**。評価器は
  cell + SST (壁関数, cfl 2) で成立 (run_0011: C_T(p) 0.9685 / 摩擦 −0.0085 / C_L 0.154 / C_M −0.980、力 STEADY)。
  twall の符号は「流体に働く traction」(viscousFlux_d.cu L527/L679) で確認し `sern_forces` に反映。
  (ii) S5: 単独抽出の δ* は膨張扇を欠損と誤認 (0.24 H) → **Euler 基準の質量流束欠損** `deltastar_sern_vs_euler` を採用
  (ランプ 0.009→0.11 H)。オフセット壁 RANS (run_0013) で壁圧が MOC に全域一致、C_L 0.146 / C_M −0.930 (Euler 0.143 / −0.939)、
  C_T は不変。(iii) S6: `opt/driver_sern.py` で 2 作動点 (cruise w0.6 / accel w0.4) Euler MOO、18 評価 17 PASS、HV 1.179、
  パレート 8 点 (L 3.9–12.5 H で C_T,w 0.960–0.977)。1 評価 ≈ 135 s。RANS 化した MOO は評価器の切替 (`evaluate.model: sst`,
  `discretization: cell`) だけで回るが、1 評価 ≈ 3× のコスト。作動点 accel (p_ext/p_in 0.15) は Euler では剥離が出ないため
  RANS で再評価が要る。
- `2026-09-04` — **node + SST 発散の真因を解決** (case/46 run_0014–0016): 真因は solver ではなく、段階起動の stage 間移植
  `interp_field.py` (座標最近傍) が**スリットの座標一致双子壁ノードを同じ元ノードに写す**こと (合成場で 134/134 station 誤写像を確認)。
  排気側壁ノードが外部流の圧力を持ち 2 次で爆発していた。`runner_sern.restart_by_index` (index コピー) に替えて板厚 0 のまま完走
  (run_0016: C_T(p) 0.9691 / 摩擦 −0.0065 / C_L 0.155 / C_M −0.981、cell と一致)。板厚オプション `mesh.cowl_thickness` は任意
  (入口側で 0 に絞ると 4 境界ノードができて発散するので入口から一定にする)。twall の符号規約は node (壁に働く力) と cell (流体に働く力)
  で逆 → `sern_forces` で離散化ごとに切替。**残観察**: node の rms_roOmega が本段で 3e18 一定 (壁ノード ω ピン留め残差の混入、場は健全) と
  境界角の単ノード圧力外れ (入口角・TE)。**教訓**: 座標一致ノードを持つメッシュ (スリット/薄板) では最近傍補間の restart を使わない。
- `2026-09-04` — **S6 RANS 版 (node SST) 完了** (case/46 run_0017): 2 作動点 (cruise / accel p_ext 0.2) + C_M ≥ −2.5 制約、
  目的は摩擦込み C_T。14 評価 10 PASS、HV 0.953、パレート 5 点 (L 5.5–9.7 H, C_T,w 0.960–0.967)。1 点 76 s。
  剥離指標 (`sep_frac_ramp`: 壁に働く接線力が逆向きの長さ割合) を `sern_forces` に追加したが、この作動点範囲では全点 0 —
  剥離を評価するには p_ext/p_in ≳ 0.5 の作動点が要る (未実施)。RANS でも前線の形は Euler と同じ。
- `2026-09-05` — **設計点固定バグ修正**: 作動点ごとに kernel (自由境界圧) を再設計していた (run_0010/0017 は作動点間で形状不一致)。
  `design_from_problem(p, design_external)` で spec.external に固定。**低 NPR 作動点** (M∞1.5, p_ext/p_in 0.6) を単点評価 (run_0018:
  C_T 0.92, C_L −0.32, C_M +1.8 頭上げ、剥離なし・滑らかな再圧縮) → 第 3 作動点 (w 0.2) にした **3 作動点 node SST MOO** (run_0019,
  10 PASS): 低 NPR が支配的で長いランプは C_T 0.83–0.90 に落ち、前線は最短ランプ (L 3.8H) に退化、剥離割合 sep_frac が最大 0.10 で発火。
  チャンバー: 超音速外部流 (M∞>1) では不要、亜音速/静止 (地上試験) には未実装 (case/23 方式が要る)。
  **S7 3D**: `meshing/mesh_sern3d.py` (2 バンド押し出し hex、カウル/側壁スリット、TE/側壁後縁共有、ランプ線も内外 2 重)、
  `evaluate/runner_sern3d.py` (index IC、quad 力積分 [CONNE = 5 整数/面])。**外側空間なし基準 (run_0023) は 2D と 1 % 以内で一致**。
  外側空間あり (run_0020) は横端トポロジで soft 段 NaN が未解決 (3 バグ修正済み: 入口タグ / 2 種入口の共有ノード / IC)。
- `2026-09-05` — **kernel の打ち切り** (ユーザ指摘: 解析領域が必要以上に広い): `PlanarMOC.march(stop_at=(f, M_c))` で f 流線を
  同時積分し key point 到達の 0.3H 先で march を止める (x_max は上限)。設計はビット同一、領域 9H → x_c+0.3。CFD 側の
  `bot_depth` 3H も超音速外部流なら 1H 程度で足りる (次キャンペーンで縮める)。
  **3D 外側空間あり**: Euler は 1 次 soft 段を完走、2 次本段でノズル幅外のランプ角部 (M∞6 が 15° 凸角で p/10 に膨張) から発散。
  加速作動点 (M∞3.5) で Euler/SST を再試行中 (run_0024/0025)。
- `2026-09-05` — **S7 3D: Euler で外側空間あり構成が成立** (run_0027, 加速点, 52.5 万セル, top_out = slip):
  幅内ランプ T −0.005/L −0.019 (2D +0.014/+0.085)、幅外ランプ T −0.025/L −0.104 → C_T 0.932 (2D 0.976), C_L −0.20 (2D 0.00),
  C_M +1.07。側壁がカウル後縁で終わると排気が横に逃げてランプ圧が 1/3 に落ちる = 3D 効果は推力 −4.5 %・揚力反転で大きい。
  3D の落とし穴 (修正済み): 入口面の幅外は inlet_ext / ランプ∩側壁の共有ノードは 2 種入口 → 2 重化 / index IC のカウル外面 /
  暖機コピーは状態量のみ (wall_dist を潰さない) / top_out 静圧出口は流れ平行で不安定 → slip。M∞6 では幅外ランプ角部の真空膨張で
  2 次が落ちる (加速点 M∞3.5 で回避)。**SST 3D は後縁 3 重点で ω → inf (run_0028) が未解決**。
  壁面出力 CONNE (3D) = [5, n0..n3]。kernel は key point で打ち切り。
- `2026-09-05` — **作動点をサイクル値で定義し直す方針を決定** (§4.10 新設、§8-1 / §8-2 決着、§5.1 に残作業の正本を新設)。
  初版の作動点は `spec.inflow` を全点で固定し NPR を外部圧側で作っていて非物理だった。アンカーに
  **NASA TM X-71972 TABLE 1** (定動圧 1500 psf 経路の station 1 / station 3 状態、$\phi$、作動モード) を採り、
  CEA2 (`tp`, H$_2$-air 平衡) で $\gamma, R, a$ を補って $M_3$ と NPR を確定した。結果: $M_\infty 6$ 巡航は
  $M_{\rm in} = 1.67$ / $p_{\rm in} = 101$ kPa / NPR 35.4 / $\gamma = 1.18$ (初版は 2.5 / 20 kPa / 20 / 1.4)。
  **低 NPR 点の正体は燃料遮断** (power-off, $M_3 = 2.9$–5.2 で全て超音速) であり低速飛行ではない。
  実装は未着手 (`operating_points[].gas` 上書きが必要)。既存 MOO (run_0010/0017/0019) の結論は作動点が
  非物理なので破棄・再取得とする。
- `2026-09-05` — **作動点サイクル値化を実装・検証** (§4.10)。`runner_sern.design_snapshot` で設計点 (入口・外部流・ガス) を
  作動点適用前に控え、`design_from_problem(p, design=)` が**それで逆設計を固定**する (従来は `p.spec["inflow"]["M_in"]` を
  読んでいたので、作動点が inflow を上書きした瞬間に形状が作動点依存になっていた)。`select_operating_point` が
  `operating_points[].gas` (gamma/cp) を上書き。`collect` の `cfd_vs_moc` は設計点 run 限定 (`on_design_point`)。
  `runner_sern3d` も同形。`driver_sern` は無変更で通る。導出スクリプト `case/46/cea/tmx_operating_points.py`。
  検証: 単体テスト ALL PASS / run_0030 (prepare ×3) で 3 作動点の**輪郭 sha 一致** (L_ramp 11.36 H) / run_0031 (m6_on, node Euler)
  **C_T 0.9059 = MOC 0.9073 (−0.15 %)**, C_L 0.3135 (0.3105), C_M −8.23 (−8.17), 力係数 STEADY, `check_quasisteady` ALL STEADY,
  残差は 1.0–1.7 桁プラトーで NOT CONVERGED (run_0002 と同性格) / run_0032 (m4_off) C_T 0.974, C_M +0.96 / run_0033 (m10_on)
  C_T 0.923, C_M −6.10。**dv 範囲を張り直した**: NPR 35.4・γ 1.183 では完全膨張 M_e 3.59 で key point は M_c ≲ 3.4 まで、
  M_c 3.2 で既に L_ramp 16–20 H (旧 3.2–4.4 は不成立) → M_c 2.4–3.4, L_cowl 上限 2.5, ni_noz 160。
  **`opt.cm_min: -2.5` は使えない** (§5.1-1b)。
- `2026-09-05` — **ランプ側外部流ブロック (§4.11) を実装・Euler で検証**。`mesh_sern.py` に `ext_top` (top バンド + `vehicle` 壁、
  `_add_ext_top`)、`vehicle_taper` (後端を厚さ 0 に絞り TE 共有; 鉛直 base + wake 版は node で step 5 発散 run_0034)、
  `_geom_start`。`runner_sern`: `vehicle` bcond (slip, 壁出力)、`paste_region_ic` が `y_top` で燃焼器領域を閉じる、YAML キー
  `mesh.ext_top/top_depth/nj_ext_top/nj_wake/vehicle_clearance/first_top_frac/vehicle_taper`。mesh テスト 27 項目 ALL PASS
  (既存 9 + ext_top 12 + taper 6; 順序互換・閉境界・辺数・CCW)。**run_0036** (node Euler m4_off, ext_top テーパ): MESH PASS
  (AR 443, skew 0.464, 60k cells)、C_T 0.9736 / C_M +0.973 (run_0032 と同値)、STEADY、残差 0.8 桁 still converging。
  **後縁圧 p/p∞ 1.317 (run_0032, 領域の産物) → 0.932 (機体上面 boat-tail 膨張後の外気)** で外気の到達を確認。
  ただし m4_off (NPR 2.8) ではランプ圧の主因は**カウル側**プルーム境界 (0.5–7H で 0.84 p∞ の軽い過膨張 → 7H 以降カウル側からの
  圧縮で 1.3 p∞ へ) で、ランプ側外気は後縁 0.1H だけを変える。剥離が出るかは cowl TE 衝撃の強さ (= L_cowl 依存) の問題で、
  L_cowl 1.2 の smoke 形状は m4_off で剥離しない可能性が高い (SST run_0038 で確認中)。
  **node SST のカウル TE ω 発散** (run_0035, 板厚 0) は ext_top と無関係で、`cowl_thickness: 0.002` で解決 (run_0037 完走、
  sep_frac 0) → 生産 YAML の既定に。
- `2026-09-05` — **`opt.cm_min` を −7.0 に決定** (§5.1-1b)。dv 箱の MOC スイープで設計点 $C_M$ の分布を取り (CFD と 0.7 % 一致するので代理可)、重み付き比 ≈0.7 を掛けて最悪 ~15 % を落とす値にした。$C_T$ 上位は M_c 3.0 / θ_r0 22° / L_cowl 2.5 (L_ramp 12–16 H, C_T 0.957–0.965) に集中 —
  θ_r0 と L_cowl が上限に張り付くので、MOO 前に範囲の妥当性 (θ_r0 上限 22°、L_cowl 2.5) を見直す余地あり。
- `2026-09-05` — **S6 MOO 再取得を投入** (case/46 `run_0040_moo_cycle3op`)。§8-6 をユーザ判断で (b) に決着させ、
  剥離制約は持たず `sep_frac_ramp` は台帳記録のみとした。構成: `problem_moo_sst_node_cycle3op.yaml` (設計点 M∞6 powered、
  3 作動点 m6_on 0.5 / m10_on 0.3 / m4_off 0.2、node SST、ext_top テーパ、`cowl_thickness` 2e-3、`cm_min` −7.0)、
  dv 5 変数、DOE 50 + infill 10×2 = 70 評価。実測 1 評価 ≈ 5 分 (SST 1 run ≈ 80 s × 3 作動点) → 約 6 時間、ディスク ≈ 6 GB。
  HV 基準点 = (−0.75, 20.5) (旧既定 −0.90 は新しい C_T 帯の半分を切り捨てる)。
- `2026-09-05` — **SST 起動レシピを確定** (§4.12) し **MOO 第 2 回を投入** (case/46 `run_0054_moo_cycle3op`)。第 1 回 (run_0040) は
  m10_on が soft 段 step 3 で `roOmega` 発散して全滅 → 切り分け (run_0041–0053): 発散は**カウル後縁のせん断層**が本体で、
  板厚は作動点で要求が逆転して使えず、CFL は場所を移すだけ。**層流暖機段**を `run_staged(warm_lam_steps=)` に追加し、
  本段 CFL も 1.0 → 0.5 に下げた (1.0 は同一メッシュ・同一 config で結果が割れる限界状態)。3 作動点の検証値は暖機なしの
  成功例と一致。DOE 40 + infill 8×2 = 56 評価、1 評価 ≈ 8 分。
- `2026-09-05` — **`L_cowl` 上限を 2.5 → 4.5 H に拡大し、キャンペーンを AWS g5 (A10G) へ移して再投入** (§8-8 決着)。
  旧箱は NASA 基準形 (3.12 H) を含まず C_T が上限で張り付いていた。run_0069 は 22 評価 (PASS 14) で打ち切り。
  併せて driver の容量節約が `.xmf` を消さず「中身の無い xmf」を残していたのを修正。
- `2026-09-05` — **平面版 Rao の一次資料を発見** (§4.3)。当初「Rao 本文未入手のため自前導出」としていたが、
  **Shyne 1988 NASA TM-100955** と **Shyne & Keith 1990 NASA TM-103175 (AIAA-90-2222)** が Rao 法の 2 次元修正 +
  **カウル切り詰め (scarf) + 外部流**を扱っており、本チェーンとほぼ同一の問題設定だった。両方 NTRS から取得して
  `papers/nozzle_design/` に配置。乗数条件は式 11–15 の形で書かれており、`rao_planar` との照合を §5.1-10 に積んだ。
  Shyne は**カウル切り詰め位置に最適値がある**と結論しており §8-8 (L_cowl 上限拡大) の直接の比較対象になる。
- `2026-09-05` — **評価戦略の判断材料を記録** (ユーザ質問への回答を反映)。(i) §4.6 に δ\* 一発で足りる理由 (C_T は不変・
  精度が採点対象でない) を実測付きで追記。(ii) §6 の「2D vs 3D で $C_T$ 差 < 2 %」判定が **不成立** (−4.5 %、揚力符号反転) と明記。
  (iii) §5.1-4 を「側壁長を dv に」から **2 段構えの具体案** (3D Euler 3 分での順位確認 → 補正サロゲート) に差し替え。
  (iv) §8-9 に「NS を Euler + 摩擦相関で代用できるか」を、摩擦が $L_{\rm ramp}$ に相関する実測 (0.53 % → 1.13 %) 付きで起票。
  (v) §5.1-8b に `L_ramp_max` の probe/設計 解像度不一致を起票。
- `2026-09-05` — **3D についての誤りを訂正** (ユーザ指摘)。(i) §4.6 の「δ\* 一発で足りる」は **2D 限定**と明記 —
  3D では横方向膨張が排除厚さの枠組みに入らないので帰還が必要。(ii) §4.3 に **3D・粘性下で Rao 最適性は保たれない**
  ことと、主張できるのは「良いパラメータ化」までであることを明記。(iii) §5.1-4 の**補正サロゲート案を撤回** —
  形が 3D 環境に応答しないので最適化になっていない。(iv) §5.1-3 を「3D SST を通す」に格上げし、2D レシピ
  (層流暖機・CFL 0.5・後縁板厚・交線丸め) の移植という具体策と NaN 位置を明記。(v) §8-10 にロードマップ
  組み替え案を起票 (**優先順の組み替え自体はユーザ判断待ちで未反映**)。
- `2026-09-06` — **作動点間 warm start を実装・検証** (§5.1-1c、ユーザ指示の最優先項目、codex レビュー A′ 準拠)。
  `warm_from_run` = 入口状態比の相似スケール + **目標 γ での roe 再構成** (素の index コピーは γ が変わると圧力が
  1.24〜2.18 倍ずれるので不可) + 適応段 500 step。YAML の `operating_points[].warm_from` で指定する。
  2 形状 × m10_on で **cold と C_T・摩擦込み C_T が 5 桁一致**、C_M 差 0.02–0.04 %、**27–28 % 短縮** (run_0072/0073)。
  本段 CFL ランプ案は撤回 (風洞の実測で warm start では利得なし)。**残る短縮余地は本段の早期停止**。
- `2026-09-06` — **m10_on 発散の真因を特定し、ω 説とテーパ形状説を棄却** (§4.11)。5 step 刻みの出力
  (`run_0075_diverge_watch`) で「P 床 → 負密度 → 圧力暴走 → ω 発散」の順序を確定。膨張自体は Prandtl–Meyer より
  緩く物理的で、床に落ちるのは 18 % のノードだけ = 局所アンダーシュート。**「形状が非物理」という結論は誤りだった**
  (ユーザ指摘)。`implicitRelax` 0.7 (+ cfl 2.0) は無効 (run_0076)。`evaluate.implicit_relax` / `p_min` を配管。
  MOO を ext_top なしで回す案は**保留** — 診断を誤ったまま逃げる判断だったため。残作業は §5.1-2b。
```

## 参考: `plans/active/tooling-nozzle-sern-3d.md`

```
# ⑤ SERN の 3D 計算領域 (本体 plan §4.14/§4.15 系の切り出し)

## メタ

- **area**: `tooling / optimization`
- **status**: `in_progress`  <!-- 2026-09-19 切り出し (本体 plan が codex の引数上限 128 KB を超えたため)。R4c/R4f/R5b 完了、R4e はメッシュ側完了・CFD 側が残 -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「SERN チェーン」節
- **related_plans**:
  - 本体: [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) (⑤ SERN チェーン全体。本 plan はその §4.14/§4.15 系と 3D の残作業を持つ)
  - 併走: [`tooling-nozzle-sern-startup.md`](tooling-nozzle-sern-startup.md) (起動レシピと収束判定。2D が対象)
- **created**: `2026-09-19`
- **owner**: `sano`

## 1. 目的

⑤ SERN の **3D 計算領域を、評価に使える形にする**。2D で確立した逆設計・起動レシピを 3D
(側壁・機体側面・機体ベース・幅外の外部流) に持ち上げ、3D SST を生産 3 作動点で成立させる。

## 2. スコープ

**対象**: `mesh_sern3d.py` の領域トポロジ、`runner_sern3d.py` の BC と力の帳簿、3D の起動と受理ゲート。
**対象外**: 2D の起動レシピと収束判定 ([`tooling-nozzle-sern-startup.md`](tooling-nozzle-sern-startup.md))、
逆設計カーネル・MOO (本体 plan)。

## 3. 関連 docs と前提

- 本体 plan の §4.9 (3D 確認 S7)・§4.11 (ランプ側外部流ブロック) が前提。
- `case/46.sern_design/README.md` の run 一覧が run の一次情報。

## 4. 設計方針

### 4.14 3D 計算領域のトポロジ (2026-09-19 起票・未決)

R5 (`run_0118`, 生産構成 m6_on, 現行バイナリ) は**層流暖機段 step 180** で `ro` NaN。位置は
**x/H 0.00–0.05, y/H 0.978–1.018, z/H 1.16–2.50** = ランプ膨張角部の直下、側壁 (z/H 1.0) の外側から
遠方境界まで。P は床 20.00 Pa、T は反転クランプ上限 6000 K。**1 次・層流・CFL 0.1 で落ちる**ので
スキームでも乱流でもない。

診断で分かったこと。

1. **2D 断面の z 一様押し出しなので、ノズルの形が全幅に複製されている**。`y_top(x)` に z が入らない
   (`mesh_sern3d.py:96`)。x/H 0.49 の断面でノード配置は z によらず同一 (すき間 1.180 → 3.841 = 固体 2.661 H)。
2. **機体幅の外 (z > W_vehicle/2 = 2.0) にも厚さ 2.66 H の固体が残る**。両面 slip (`underside_far` / `top_out`)。
   自由流しか無いはずの場所に幻の物体がある。**ここは無条件で fluid にすべき** (ユーザ指摘 2026-09-19)。
3. **同じ 20.2° の角を、内側は排気 (101 kPa)、外側は自由流 (2.85 kPa) が回る**。Prandtl–Meyer:
   排気 → 36.8 kPa (床の 1800 倍)、自由流 m6_on → 50.6 Pa (床の 2.5 倍)、m10_on → **0.2 Pa = 床の下**。
4. **node 壁列の市松が重畳**。角直後の壁法線方向に P を並べると 20 Pa と 1200 Pa が 1 ノードおき
   (`y/H` 1.0051 で 1171, 1.0052 で 20.0, 1.0062 で 26.1, 1.0063 で 1309)。滑らかな膨張なら 50 Pa 前後で揃うはず。
   既知の [[node-slip-spurious-flow]] / [[node-wall-entropy-checkerboard]]。
5. **解離の有無は無関係** (ユーザ指摘で訂正): NASA-9 は `Thi` 6000 K 超を端点 cp で線形外挿する
   (`thermo_d.cuh:193`) ので数値的に健全。~~効いているのは `DEPVAR_TMAX 6000` ハードクランプ~~
   → **撤回 (2026-09-19, codex M1、自分で再検証済)**: `run_0118` で T=6000 K は **1 ノードだけで既に NaN**。
   有限のまま床に張り付いているのは **P=20 Pa の 140 ノードのうち 139 が T=50 K** =
   **低温側 `DEPVAR_TMIN`**。膨張は冷える方向なので高温側を疑ったのは方向を誤っていた。
   codex の独立積分 (凍結組成・等エントロピー): m6 20.2° → 50.4 Pa / 69.5 K / 2.53e-3、
   **m10 20.2° → 0.147 Pa / 18.0 K / 2.83e-5** で、`pMin` 20 Pa・`DEPVAR_TMIN` 50 K・`roMin` 1e-4 の
   **三つの床に同時に抵触**する。さらに `dependentVariables_d.cu:209` は**クランプ後の T から `roe` を再構成**
   するので、床は「不整合」ではなく**保存エネルギーの非保存的な書き換え**である。
6. **側壁は slip ではない** (`sidewall_in/out` = `wall`)。ただし **`vehicle_top` は slip** で、
   §4.11 の 2D の扱い (機体上面は中立モデル) をそのまま 3D に持ち込んでいる。

**未決 = 機体幅内 (W/2 < z ≤ W_vehicle/2) の下面をどう与えるか。** 候補:
(A) 前体の延長で水平 — 側壁外側に最大 2.78 H の段差、側壁は x ≤ L_sw にしか無いので下流を塞ぐ面が要る。
(B) z 方向に滑らかにブレンド — 段差・面の追加なしだが z/H 1.16 はブレンド初期で効果が薄い。
(C) 緩い boat-tail 角 (6° 程度) に固定 — 膨張が m6_on 51 → ~700 Pa, m10_on 0.2 → ~60 Pa に緩む。
(D) `W_vehicle` を W まで詰め、幅外は全て fluid — 最も単純。横方向レリーフは側壁 (no-slip, x ≤ L_sw) が決める。


#### 4.14.1 codex plan レビュー (2026-09-19) の採用結果

判定 **NO-GO** (C1/M4/m0)。記録 [`notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan.md`](../../notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan.md)。

| # | 指摘 | 採否 |
| --- | --- | --- |
| C1 | **D 案は単体で成立しない**。`W_vehicle` は境界タグを選ぶだけで、4→2 で座標・接続は同一 (`mesh_sern3d.py:228`、上面は `:288` で全幅押し出し)。幅外を fluid にするには旧下面〜旧上面を流体で埋め、`z=W/2` に**機体側面** (ダクト側壁とは別物、`L_sw` より下流にも在る) を新設し、交線のノード定義と新領域の組成・初期値を決める必要がある。既存メッシュテスト 30 項目は**幻の固体を含む領域を合格させる** | **採用** → §5.1 R4c |
| M1 | 「1 次だからスキームではない」「6000 K クランプが原因」は断定不可。T=6000 K は 1 ノードで既に NaN、床は T=50 K 側 | **採用** (§4.14-5 を撤回・書き換え、自分で再検証済) |
| M2 | m10 は `pMin` だけでなく `DEPVAR_TMIN` 50 K・`roMin` 1e-4 にも抵触。幾何修正だけでは保証できない。ゲート (`sern_gates.py:39`) は**床に張り付いた有限解を拒否しない** | **採用** → §5.1 R5b |
| M3 | **3D runner が生産の等温壁を無視**。`runner_sern3d.py:44` が `kind: wall` を直書きし、2D の `p.wall_bcond_line()` (`runner_sern.py:271`) を使っていない。run_0118 は断熱で回っていた | **採用・修正済 (2026-09-19)**。ただし run_0119 (等温 1000 K) は **step 76 で発散**し、これ単体では直らない (codex の予告どおり) |
| M4 | `side_far` が slip 固定で横方向流出を遮り反射する。`r4_domain_study.py` は加速点 Euler で許容値は $\|\Delta C_T\|<0.002$ のみ、生産 TP・SST と $C_L/C_M$ の独立性は未保証。`run_0102_r4_domain/` の成果物は本 checkout に無い | **採用** → §5.1 R4d |

**run_0119 (等温壁を効かせた R5)**: step 76 で `ro` NaN。NaN は **z/H 1.000–1.108 = 側壁のすぐ外**に移り (断熱の run_0118 は 1.16–2.50)、
T=50 K 床の 138 ノードは**全て z/H 2.500 の遠方境界列**で NaN 位置とは別。床は同時多発しており単一の連鎖ではない。

### 4.15 3D の残り 3 件 (2026-09-19 起票)

R4c (領域トポロジ) と出口境界を直した結果、**3D SST が初めて全段完走**した (`run_0122_r5_outflow`, 12000 step,
`check_convergence.py` = NOT CONVERGED (stalled) で **DIVERGED ではない**)。到達までの分解:

| 問題 | 原因 | 対処 | run |
| --- | --- | --- | --- |
| ランプ角部の幅外で床に張り付く | 2D 断面の z 一様押し出しでノズル形状が全幅に複製 | 機体をノズル幅に閉じ幅外を流体化 (`vehicle_side`) | 0120 |
| 後縁の三重点 | カウル後縁と側壁後縁が同位置 | `L_sw` 0.8 で分離 | 0121 |
| **出口から逆流が遡る** | `outlet_statPress` の Ps 2851 Pa が実出口圧 ~6.5 kPa より桁違いに低く、node の亜音速壁列から unstart | **`evaluate.outlet_kind: outflow`** (全量外挿) | 0122 |

出口の効果は明確: 出口プルーム圧 7.5→128 kPa (発散) が **6.38→6.52 kPa で静定**、Ux 1714→153 が **1975→1881**、
逆流ノードは全 10 スナップショットで **0**。

**残り 3 件。**

1. **`rms_roOmega` が 8.95e+15 で張り付く (最優先)**。
   ~~位置は入口面 ∩ 側壁の稜線~~ **撤回 (2026-09-19, codex plan レビュー 2 の C1)**。
   `roOmega` の**最大値の位置**と**残差の発生場所**は別物で、私は前者を見て原因と断じた。真因:

   - 残差は最初から張り付いていない。**step 1450 で 0.224 → 1451 で 1.235 → 1452 で 8.29e15** と跳ぶ
   - **ノード ID 646932 の 1 点**。位置 (x/H, y/H, z/H) = **(10.5397, 3.76611, 0)**、所属は **`sym` のみ**
   - そこで `roOmega = 1e-20`、`k = 0` に張り付く。実装と同じ面補間で独立再計算すると
     `∇k·∇ω = −1.07048e16`、`k=0` より `F1=0`、ω の分母下限 `1e-12` を使った
     **交差拡散項 × 体積 = −9.58791e18**。全 1,146,866 ノードに対するこの 1 点の RMS 寄与 = **8.95298e15**
   - **保存場が有限に見えるのは、負の更新を下限で切っているため** (`update_d.cu:343`)。
     「値が動かない = ピンされている」という私の解釈は逆だった

   根拠: `ransSource_d.cu:136`, `同:212`, `update_d.cu:343`。
   **切り分け手順**: `res_1000.h5` から診断 run を作り step 1452 相当まで細かく保存し、
   `res_roOmega` / `omg_trans,prod,dest,cross` / `dK*` / `dOmega*` / `sstF1` / `src_jac_omega` /
   `transport_diag_omega` / `dt_local` / 更新前後の `roK,roOmega` を取る。閾値超過**点数ではなく二乗和寄与率**で
   局在化し、下限到達**前**の輸送・負の交差拡散・更新量を追う。最初に下限へ落ちる原因は未確定。

2. ~~**`vehicle_top` を probe として出力**~~ **取り下げ (2026-09-19 ユーザ判断 + codex M2)**: 既に `outputHDFflg: 1` で面出力済み。`probe/point_probes.cu:64` は `points` しか読まず `surfaces` の実装が無い。以下は旧記述。現状 `probe.yaml` は `surfaces:` が空で、壁面出力は
   `outputHDFflg` 経由のみ。機体上面は **slip** なので摩擦は無いが、そこを流れる外気の状態
   (圧力・温度・マッハ数) は外圧の境界条件として効く。probe に載せて時系列で見られるようにする。
   **未決**: 機体上面を slip のままにするか粘性壁にするか (codex M3 の「別途明示」)。
   帳簿 (ノズル $C_T$) には入らないので力には効かないが、外気の状態には効く。

3. ~~**対称面の鏡像フィルタ**~~ **実装済 (2026-09-19, `solver_density_cuda/tools/mirror_symmetry.py`)**。codex M4/M5 を全件反映: 名前の末尾判定 → **型の登録表** (`limiter_Uz` の誤反転を修正)、`ids[::-1]` → **要素種別ごとの向き置換** (単位 hex で Jacobian +1 → −1 → +1 を確認)、未知フィールドは**警告**、h5 属性を保持、対称面の `Uz` (実測 −239.85〜+329.22 m/s) は黙ってゼロ化せず元解の対称条件誤差として報告。複数面を `patch` id つきで束ねる `--merge` 相当も実装。以下は旧記述。計算は半スパン (z ≥ 0) なので、可視化・確認のたびに片側しか見えない。
   z = 0 の対称面で形状と値を反転して全幅を復元するツールを作る。
   **ベクトルの扱いが肝**: スカラー (ro, P, T, k, omega, wall_dist) はそのまま複製、
   **z 成分を持つベクトルは符号反転** (`roUz`, `Uz`, および勾配の `d*dz` 系)。
   テンソル成分 (`dUxdz`, `dUzdx` など) は z が奇数回現れる成分だけ反転する。
   対称面上のノードは重複させない (座標一致ノードを作ると [[interp-field-coincident-nodes-trap]])。


#### 4.15.1 codex plan レビュー 2 (2026-09-19) の採用結果

判定 **NO-GO** (C2/M5/m1)。記録 [`notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan-2.md`](../../notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan-2.md)。

| # | 指摘 | 採否 |
| --- | --- | --- |
| C1 | **ω 診断が誤り**。真因は `sym` 上 1 ノードの下限張り付き + 交差拡散 (上記)。かつ **ゲートが異常解を受理している** (`sern_gates.py:104` は有限かつ非上昇なら PASS、`require_residual_pass=false`) | **採用**。§4.15-1 を書き換え、R5b に統合 |
| C2 | **R4c 未完了**。`mesh_sern3d.py` の側面バンド末端を `vehicle_side` の壁で閉じており、その壁は z = W/2 に限らず**遠方境界 (z/H 2.5) まで延びる 224 面**。幅外の流体が下流へ抜けられず `side_far` との交線で最大 **5.72 MPa** (外気 2851 Pa)。**run_0122 の結果はこの影響下にある** | **採用**。楔で縮退させる案は閉性が壊れた (漏れ 32・余り 10) ので保留、欠陥を docstring に明記 → R4e |
| M1 | **`outflow` の効果を「出口だけ」と断定できない**。`runner_sern3d.py:47` の共通関数が `outlet` と `bottom` の両方に効き、同時に変わっている。`bottom` の面法線 Mach は **−5.2e-6〜4.3e-6** で超音速流出面ではない。2851 と 6.5 kPa の差は 2.3 倍で「桁違い」でもない | **採用** → R4d に統合 (出口と遠方境界を分離して比較) |
| M2 | probe の `surfaces` は実装が無い。面出力は既にある | **採用** (取り下げ) |
| M3 | 「帳簿外だから力に効かない」は**誤り**。上面の境界層・熱伝達・後縁流れが変わればランプ/カウルの圧力も変わる。**生産の `vehicle_top` は機体側面と同じ等温粘性壁に統一**を推奨 | **採用** → R4f |
| M4/M5 | 鏡像ツールの符号規則は**型**で定義すべき (`limiter_Uz` 誤反転)、`ids[::-1]` は hex の向きを戻さない、未知フィールドの黙認、h5 属性喪失、対称面の反対称成分 | **採用・実装済** |
| m1 | plan / README / 検証条件の不同期 (§6 の「全作動点 PASS」と §4.7 の「PASS 任意」、run_0122 が README に無い 等) | **採用** (本コミットで同期) |


#### 4.15.2 R4e の案 (d): 機体後縁に**有限の厚み**を持たせる (2026-09-19 ユーザ提案)

3 案 (楔縮退 / 鈍頭化 / prism) が境界書き出しの連鎖で止まったのを受けての提案。
**機体をテーパで厚さ 0 まで絞らず、有限の厚みのままベース面で終わらせる**。

**「カット (手前で切る)」ではなく「厚みを用意する (全長 + ベース)」を採る。**
SERN はランプそのものが機体下面なので、ランプの上には後縁まで機体構造が要る。
後縁より手前でカットするとランプの最後の区間だけ上に何も無くなり、実機として成立しない。
したがって機体は $x = L_{\rm ramp}$ まで全長あり、そこで厚み $t_{\rm base}$ のベース面で終わる。

利点:
- 側面バンドが**厚さ 0 に潰れない**ので末端の退化 (楔・鈍頭・prism のどれも不要) が消える
- 末端を閉じる面が「実在する機体ベース」になり物理的な意味を持つ (今の閉じ面は人工物)

レビューで見てほしい論点:
1. $L_{\rm ramp}$ より下流では上バンドの下端が旧ランプ線 (プルーム上線) に落ちる。幅外では
   上流がバンド (yt..y3)、下流が上バンド (yt..天井) となり**同じ領域を別の節点分布が担当**するので
   節点が合わない。合わせるにはバンドの分布と上バンド下部の分布を一致させる必要があり、
   案 (e) (バンド上端を直線 $Y_{\rm out}$ にして出口まで通し上バンド下端を z 依存にする) の制約と同型ではないか
2. ベースの下流に**後流域**ができる。2D の鉛直 base + wake は **node で step 5 発散**した
   (`run_0034`, ランプ∩base と上面∩base の 90° 二重 slip 角)。3D の機体ベースで同じ病が出ないか。
   肩を丸める (2D のランプ角部丸めと同じ発想) で回避できるか
3. $t_{\rm base}$ をどう決めるか。現行の `vehicle_clearance` 0.02 H / `first_wall_frac` 0.004 H との関係
4. ベース面の境界条件 (等温粘性壁か slip か) と帳簿の扱い (ノズル力に入れない)

案 (e) も比較対象として残す。実装量・数値的頑健性・物理的妥当性でどちらが優れるかを判定してほしい。

#### 4.15.3 R4e の確定設計 (2026-09-19, codex plan-3 GO-with-changes → **案 (d) 採用**)

**採用: 案 (d) 「全長機体 + 有限ベース + 適合した側面・後流バンド」。** 案 (e) は薄い後縁モデルのままで
「有限厚の機体」という要件を満たさず、z 依存の分割線の設計が余分に要るため不採用。

**後縁でブロック構成を切り替えない** (codex M1)。旧案の論点 1 (後縁下流で上バンド下端をランプ線に落とす) は
不要な前提だった。幅外バンドを**出口まで維持**する:

- 下端 $a(x) = y_t(x)$ (ランプ/プルーム上線)、上端 $b(x)$、後縁で $b(L_{\rm ramp}) = y_e + t_{\rm base}$。
  下流の $b(x)$ は**後流ブロックの上端として連続に延長**する。
- 幅外 ($z > W/2$): $a..b$ の流体バンドを全長に置く。
- 幅内 ($z \le W/2$): 同じ $a..b$ のブロックを**後縁下流だけ**に置く。
- 上バンド: 全幅で $b..$天井 を担当する (担当の切り替えが無くなる)。
- ベース壁: $x = L_{\rm ramp}$、**幅内の $a..b$ だけ**に置く。

節点は共通分布 $y_j = a + \eta_j (b-a)$ と**同じ節点 ID** を使う。これで側面バンドと後流の接続が適合し、
幅外を塞ぐ面 (現行 224 面) が消える。codex は接続骨格 (92 hex) で非多様体面 0・流体域連結・幅外の閉じ面 0 を確認済。
**これは接続の成立確認であって CFD 安定性の実証ではない。**

**物理寸法を格子から独立させる** (codex M3)。現行は `clearance = max(vehicle_clearance, 3*first_top_frac)`
([`mesh_sern3d.py:274`](../../design/forge_design/meshing/mesh_sern3d.py))、機体終端は `first_wall_frac`
([同:288](../../design/forge_design/meshing/mesh_sern3d.py)) で決まるので、**格子を変えると形状が動く**
(実測: `first_top_frac` 0.02 → 0.01 で機体上面基準高さが 3.840894 H → 3.810894 H = 0.03 H 移動)。
これを残すと格子独立性試験が成立しない。$t_{\rm base}/H$・機体上面輪郭・丸め半径を**物理入力として固定**し、
格子が足りなければ**メッシュ生成を失敗させる**。`vehicle_clearance` 0.02 H は $t_{\rm base}$ の根拠にならないので、
$t_{\rm base}/H = 0.02$ は**暫定モデル値**と明記し、0.01 / 0.02 / 0.04 の形状感度を別試験にする。

**BC と帳簿を先に決める** (codex M4)。`vehicle_top` / `vehicle_side` / `vehicle_base` を
**同一の `wall_isothermal`・同一 $T_w$ に統一**する (R4f を後回しにせず案 (d) の生産仕様に含める)。slip は診断用に限定。
集計は明示的な面集合で三分する: **ノズル力** = ランプ・カウル・ダクト側壁 / **機体力** = 上面・機体側面・ベース /
**全収支** = 両者を含む全境界。ベースはノズル力から除くが**機体力として保存する** (現行 `forces3d` は
集計対象に `vehicle_side` が無く、特別扱いしない面をノズル力に加算してしまう)。

**~~`run_0034` の base+wake 発散は soft 段・90° 二重 slip 角が原因~~ 撤回 (2026-09-19, codex M2)**。
実行ログは `cfl_pseudo=4`・`convMethod=1` (**2 次**) で「soft 段」ではない。`res_nan_5.h5` の密度非有限 79 点が
ベース近傍に集中することは確認できるが、それだけでは角の処理を真因と特定できない。現行の
[`slip_d`](../../solver_density_cuda/cuda_forge/boundaryCond_d.cu) は面ごとに ghost を作るので
「共有ノードを二方向へ順次射影するから必ず破綻する」という説明は実装に裏付けられていない。
**正しい記述: ベース近傍で発散、根因未確定。** 同じ有限ベース形状で現行バイナリ・1 次・低 `cfl_pseudo` の
2D 診断を先に置く。生産の等温 no-slip は $u=0$ を共通に課すので旧 Euler/slip と同一条件ではない。
**肩の丸めは検証候補であって発散回避の保証ではない** — 採るなら半径・接点・変更するランプ区間まで定義する。

**検証を実装前に固定する** (codex M5)。現行の `run_sern_mesh3d_tests.py` は今回も ALL PASS だが、
機体上面の z を**中央値**で見ており、角の行列式も**絶対値**しか見ないので人工壁・負向き要素を排除できない。順序:

1. **接続**: 内部面は同じ ID で 2 セル共有、境界面は 1 セル所有、重複タグ 0、**幅外の横断壁 0**、ベース面積を設計値と照合
2. **幾何**: **符号付き** Jacobian、双対体積・閉性、float32 変換後の節点衝突、`check_mesh_quality` PASS
3. **起動**: 2D 有限ベース診断 → 3D 層流 → SST。追加流体には領域別の組成・状態を与え、**変更前メッシュからの無条件 index コピーを禁止**
4. **受入**: R5b を先に有効化し、全残差 PASS・目的量 STEADY・持続する床/更新クリップなし (壁の正当な `k=0` ピンは区別する)
5. **独立性**: 生産 3 作動点で固定形状の格子・領域感度。暫定基準 $|\Delta C_T| \le 0.002$, $|\Delta C_L| \le 0.002$, $|\Delta C_M| \le 0.02$

ベース圧・再循環長も時系列監視に入れる。定常擬似時間で振動する場合、**その振動を物理的な後流変動と解釈しない**。

#### 4.15.4 R4e 着手順 ③: 2D 有限ベース診断の結果 (2026-09-19)

3D に進む前に、**同じ形状 (テーパ + 厚み `t_base` のベース) を 2D で**回した (codex plan-3 M2/M5 の指定)。
2D メッシャに「テーパ + 薄いベース」モードを追加 (`mesh.t_base`)。従来の `vehicle_taper = 0` は**全高の鉛直ベース**
($h_{\rm base} \approx 1.9$ H) で、`run_0034` が落ちたのはそちら — 今回のベースは**その約 1/100**。

| run | 形状・設定 | 結果 |
| --- | --- | --- |
| `run_0202_r4e_2d_base` | `t_base` 0.02、ベース = **slip** (2D 既定)、soft cfl 1.0 | 暖機 (層流) 3 段**完走** (残差 −2.5 桁) → **soft (SST 1 次) step 7 で NaN**、42 節点 |
| `run_0203_r4e_2d_base_cfl01` | 同上、soft cfl **0.1** | **soft step 57 で NaN**、同じ 42 節点・同じ位置 |
| `run_0204_r4e_2d_base_t005` | `t_base` **0.005**、soft cfl 1.0 | **soft step 20 で NaN**、同じ 42 節点 |
| `run_0205_r4e_2d_base_wall` | `t_base` 0.02、ベース = **等温粘性壁** (`evaluate.vehicle_kind: wall`) | **soft 1500 step 完走** (`rms_ro` 2.9e-5 → 1.0e-5) → **mid (2 次) step 28 で NaN**、34 節点 |
| `run_0206_r4e_2d_base_wall_barth` | 同上 + `limiter: 1` (Barth) | **mid step 27 で NaN** (Venkatakrishnan と同じ) |

**読み取れたこと**:

1. **CFL 律速ではない**。cfl を 1/10 にすると step が約 8 倍 (7 → 57) = **同じ擬似時間で落ちる**。構造的な問題。
2. **厚みを 1/4 にしても消えない** (step 7 → 20、節点数は同じ 42)。ベースを薄くするだけでは解決しない。
3. **slip のベースが SST と相性が悪い**。鋭い 90° 角の slip 面で `wall_dist` が 0 に落ちる
   ([[sst-mesh-walldist-gotcha]] と同型)。**等温粘性壁に変えると SST 1 次段を完走する** —
   codex M4 の「上面・側面・ベースを同一の壁条件に統一せよ」を裏付ける実測。
4. **残る壁は 2 次再構成**。壁版の NaN は x/H 10.773 (= L_ramp) から下流・y/H 2.706–2.728 (ベース高さ帯) に
   34 節点で、有限 `P` の最小が **0** = ベース角で 2 次再構成が圧力を潰している。
   **リミッタでは直らない** (Barth も Venkatakrishnan も同じ step)。

**上の 4 点のうち 2 つは撤回した (2026-09-19, codex plan-4)**:

- ~~「ベース角の 2 次再構成が `P` を 0 に潰している」~~ → **`P=0` は症状で、しかも角ではない**。
  `run_0205/res_nan_28.h5` の `P=0` は**節点 62103 の 1 点**で、位置は (x/H, y/H) = (10.7732115, 2.7181568) =
  **ベース下端 2.7081569 と上端 2.7281569 のちょうど中央**。その点の `ro`/`roe` は**既に NaN**。
  `ro`/`roe` が有限かつ `ro>0` の点に限れば**最小圧力は 584.09 Pa** (床 20 Pa の遥か上)。
  正しい記述: **2 次切替後、ベース近傍で発散。最初に破綻する演算は未特定** ([[detectnan-reports-symptom-not-cause]] と同じ罠)。
- ~~「リミッタでは直らない (Barth も Venkatakrishnan も同じ)」~~ → **両者は同じ不整合を共有しているので切り分けになっていない**。
  §4.15.5 を参照。

**切り分けが不足していた点** (codex M3/M4/M5):

- `run_0202`–`0204` の CFL・厚み掃引は **slip ベース・SST 1 次**でしかやっていない。**等温壁・2 次段の低 CFL 試験は無い**。
- `evaluate.vehicle_kind` は**機体上面とベースをまとめて**変える。「ベースだけを壁にした効果」ではない。
  slip 版のベース 9 節点の `wall_dist` も 0 ではなく 0 → 0.002 m に増える。
- 2D 診断は 3D 生産メッシュと**同じ形状ではない** (`L_ramp` 10.773 vs 10.643、`theta_e` −1.49° vs +5.86°、
  `L_cowl` 1.20 vs 1.60、`y_veh` 2.78 vs 3.84)。さらに `t_base` の作り方も違う
  (2D は既存上面に `t_base·s³` を加算 = 後縁勾配も変わる、3D はエルミート終点を `y_e + t_base` にする)。
  **「同種トポロジの診断」としてのみ扱う**。
- `nj_wake` はベース高さ方向 = 鉛直ベースに対して**接線方向**。ベース直後の壁法線間隔 (Δx/H ≈ 0.082 = ベース厚の 4.1 倍)
  は `nj_wake` を増やしても変わらない。`t_base` 0.02 → 0.005 は Δy も 0.0025 → 0.000625 に変えるので純粋な形状感度でもない。
- soft 終端の `rms_ro` は **7.00e-7** (上表の 1.0e-5 は途中 step の読み違い)。

#### 4.15.5 node の再構成点とリミッタ評価点の不一致 (2026-09-19, codex plan-4 M2)

**forge の node 2 次計算は、リミッタが制限する点と流束が実際に再構成する点が違う。**

- 再構成先: `convectiveFlux_{slau,roe,hlle}_d.inc.cuh` は `g_reconEdgeMid == 1` のとき目標点を**エッジ中点**にする。
  `convectiveFlux_d.cu:69` は `const int rem = (cfg.discretization == "node") ? 1 : 0;` = **node では常に有効**で config から切れない。
- リミッタ: `limiter_d.cu` の pass2 は試行値を `Qt = qc + g·(pc[ip] − c0)` = **双対面重心**で評価する。
  `g_reconEdgeMid` の分岐が**無い**。Barth (`SCHEME==1`) と Venkatakrishnan の両方が同じ経路。

したがって**リミッタの単調性保証は、実際に使われる面値には掛かっていない**。
実測 (`run_0205` の 2D SERN メッシュ、内部面 123586):

| ズレ / エッジ長 | 面数 | 割合 |
| --- | --- | --- |
| 中央値 | — | 0.0090 |
| > 1 % | 58899 | **47.7 %** |
| > 10 % | 2852 | 2.3 % |
| 最大 | — | 0.478 |

**これは SERN 固有でなく node 2 次の一般的な欠陥**なので、**別のソルバ plan** で扱う (codex の指定)。
本 plan では「ベース発散の原因候補」としてのみ扱い、**丸めより先に検証する**。
共有勾配のゼロ化は不要 (codex)。

**次の試験順 (codex 指定、上から)**:

1. **現行失敗の再現と演算段階の特定**。`run_0205/res_0.h5` を共通初期場にし、1 次継続と 2 次切替を比較。
   最初の異常まで毎 step 保存し、面再構成 → 流束 → block-DPLUR → 化学種更新 → EOS → 壁ピンのどこで
   非有限・非許容が出るかを捕捉する。**NaN ダンプから時間順序を逆算しない**。
2. **再構成点の整合試験** (§4.15.5)。現行と、リミッタ評価点をエッジ中点に揃えた実装を比較する。
3. **同じ等温壁・2 次条件で時間積分を分離**。`cfl_pseudo` 1 → 0.1、次に `implicitRelax` 1 → 0.7 を個別に。
   step 数 × CFL を「同じ擬似時間」の証明にしない。単なる延命を成功扱いしない。
4. **固定形状の局所格子試験**。`t_base/H = 0.02` のまま、まず後縁前後の Δx を縮める → 上下接続の間隔比 → `nj_wake` 9 → 17。
5. **最後に肩の丸めと厚みの形状感度**。丸めは仕様を先に固定する:

   | 項目 | 仕様 |
   | --- | --- |
   | 形 | 原輪郭と鉛直ベースに接する円弧を**固体側**に置く。上下の半径・接点・下側ランプの変更区間を明記 |
   | 半径 | 暫定比較 `r/H` = 0 / 0.0025 / 0.005。**物理入力として固定**し、円弧は最低 8 区間で解像して再細分化 |
   | ベース残存高さ | 各円弧がベース上で消費する長さを `d_lower, d_upper` として **`t_flat = t_base − d_lower − d_upper > 0` を必須**。非直角では `t_base − 2r` としない |
   | `t_base` | 0.01 / 0.02 / 0.04 H は**別の**形状感度。**薄くして安定化する方針は採らない** |

   丸めた下側肩は**ノズル輪郭そのものを変える**ので、変更区間の力をノズル帳簿に含め、
   「ベース面積 = `t_base` × 幅」の試験を新しい面集合に合わせて更新する。

#### 4.15.6 R4e ①: 破綻する演算の特定 (2026-09-19, 毎 step 診断)

`run_0205/res_0.h5` (soft 完走後 = mid 段の初期場) を**共通初期場**に、`output.level: 2` (勾配・リミッタ込み)・
`outStepInterval: 1` で回した。run はすべて `case/46.sern_design/` 配下:

| run | 条件 | 結果 |
| --- | --- | --- |
| `run_0207_r4e_diag_1st` | 1 次 (`convMethod: 0`), cfl 1.0 | **40 step 完走・場は不動** (`P` min 1064.85 → 1066.58、`T` min 220.57 一定) |
| `run_0208_r4e_diag_2nd` | 2 次, cfl 1.0 | **step 28 で NaN** (run_0205 と同じ step = 再現) |
| `run_0209_r4e_diag_nobase` | **ベース無し**の生産形状 (`run_0198` から), 2 次, cfl 1.0 | **40 step 完走・場は不動** (`P` min 1382.05 一定) |
| `run_0210_r4e_diag_2nd_nowd` | 2 次, cfl 1.0, **`nodeWallDirichlet: 0`** | **40 step 完走** |
| `run_0211_r4e_diag_2nd_cfl01` | 2 次, **cfl 0.1** | **step 184 で NaN** |

**① 最初に壊れるのは「機体ベースの壁節点そのもの」**。`P` 最小は step 1 から一貫して
**x/H = 10.77321 (= `L_ramp` 厳密)** の節点 62100–62106 (`vehicle_base` 面の壁節点) にあり、そこで

| step | 0 | 5 | 10 | 20 | 27 |
| --- | --- | --- | --- | --- | --- |
| `P` [Pa] (節点 62104) | 1159.3 | 964.7 | 736.7 | 302.8 | **31.6** |
| `ro` [kg/m³] | 0.00358 | 0.00284 | 0.00215 | 0.00082 | **0.00003** |

と **ρ が 120 倍減**する。一方**ベース帯 41 節点の平均 `P` は 4694 → 4346 (7 % 減)** にすぎない。
つまり帯全体の排出ではなく**壁節点の双対 CV だけが単調に空になる**。
その節点の `wall_dist` = 0・`Uy` = 0 はどちらも正常 (壁節点として正しい)。双対体積は 1.024e-06。

**② CFL 律速ではない** (codex M3 が指摘した「等温壁・2 次での低 CFL 試験が無い」を実施):
cfl 1.0 → step 28、cfl 0.1 → step 184。**10 分の 1 の CFL で 6.6 倍の step** = ほぼ同じ擬似時間で落ちる。

**③ 壁速度 Dirichlet が排出に必要**: `nodeWallDirichlet: 0` にすると同じ 2 次 run が 40 step 完走し、
同じ節点の `P` は 1159 → 3364 → 1397 と**むしろ上がる**。ただしその節点は `Ux` 1093 / `Uy` 1041 m/s を持ち
**壁になっていない**ので、これは**診断であって対策ではない** ([[node-wall-velocity-needs-dirichlet-flag]])。<br>
**機構の仮説**: ベース壁節点の双対 CV は、境界半割面が壁 (対流質量流束なし) で、自分側の再構成も `u = 0` を運ぶため、
**質量の流入経路がない**。内部面の正味が流出なら、補充されずに単調に空になる。2 次は隣接側の再構成を強めるので排出が立ち、
1 次では立たない。`nodeWallDirichlet: 0` で止まるのは、節点が有限速度を持てば流入できるため。

**④ codex M2 (リミッタと再構成点の不一致) は実在するが、ベース発散の原因ではない**。
level 2 ダンプから実装と同じ式 (`phif = phiC + ψ·g·cpd`) で両点の再構成 `P` を計算すると:

| step | リミッタの評価点 (双対面重心) 最小 | SLAU が使う点 (エッジ中点) 最小 | 中点で P≤0 の面 |
| --- | --- | --- | --- |
| 1 | 974.0 Pa | **−10252.6 Pa** | 1 |
| 5 | 452.6 Pa | **−20690.7 Pa** | 1 |
| 27 | −280.3 Pa | −16717.5 Pa | 5 |

**step 1 から流束に負圧が入る**。場所は**カウル後縁 (x/H 1.2005)** — ベースではない。
しかし**ベース無しの生産形状 (`run_0209`) でも同じ違反が同じ場所で毎 step 起きており** (中点 −13747 Pa、P≤0 の面 1)、
その run は 40 step 完全に安定。**したがってこの不整合は常設の欠陥ではあるが、ベースを殺しているものではない**。
§4.15.5 の別ソルバ plan 化は維持する (それ自体は直すべき) が、R4e の原因からは外す。

**次**: codex ③ の残り (`implicitRelax` の分離) と ④ (固定形状の局所格子)。
④ は仮説と整合する — ベース壁 CV は双対体積 1.02e-06 と極小で、その直後の流れ方向間隔はベース厚の 4.1 倍 (codex M5)。
**飢えた小 CV が巨大な隣接へ吐く**構図なので、後縁前後の Δx を縮める試験を先に置く。

#### 4.15.7 R4e ③④: 時間積分と局所格子の分離、そして機構の特定 (2026-09-19)

**③ 時間積分は機構ではない** (codex M3 が求めた「等温壁・2 次での分離」):

| 条件 | NaN step |
| --- | --- |
| cfl 1.0 / `implicitRelax` 1 | 28 |
| **cfl 0.1** / relax 1 (`run_0211`) | 184 (10 分の 1 の CFL で 6.6 倍 = ほぼ同じ擬似時間) |
| cfl 1.0 / **relax 0.7** (`run_0212`) | 39 (**1.4 倍の延命のみ**) |

乱流は先行しない: ベース帯の `k` が跳ねるのは step 25 で、`P` は step 2 から落ちている。
したがって乱流更新凍結の対照は不要 (codex の条件付き指定を満たさない)。

**④ 局所格子を細かくすると悪化する**。`mesh.split_plume_at_te` を実装し (プルーム区間を後縁で分割し
前後にクラスタを置く)、ベース直後の Δx を **0.08194 → 0.02451** (ベース厚の 4.1 倍 → 1.2 倍、`ni` は同じ 360) に縮めた
(`run_0213_r4e_texd`)。結果は **mid 段 step 10 で NaN** (従来 28)。
ベース壁節点の双対体積は 1.024e-06 → 3.063e-07 (**3.3 倍減**) で、落ちる step は **2.8 倍早い** = ほぼ体積比例。
→ **「極小 CV が巨大な隣接へ吐く」仮説は棄却**。正しくは **一定の流出率が固定体積の CV を空にしている**。

**機構 (特定)**: **速度を Dirichlet でピンした壁節点から、2 次 MUSCL 再構成が面に非ゼロの速度を置く。**

実装と同じ式 (`phif = phiC + ψ·g·cpd`、node は `cpd` = エッジ中点) で `run_0208` step 1 を再現すると、
ベース壁節点 62100 は**節点値 `Ux = Uy = 0`** なのに面への再構成が

| 面 | 再構成 `uL` [m/s] | `vL` [m/s] |
| --- | --- | --- |
| 121053 | **457.42** | 143.96 |
| 123578 | −72.16 | −22.71 |

となる。境界半割面は壁なので**対流質量流入が無く**、内部面の正味が流出ならその CV は補充されずに空になる。
1 次では `uL = u_node = 0` なので立たない。

**なぜベースだけか** (同 step、速度ピン壁節点 1028 個を集計。再構成状態から見積もった正味外向き質量流束を
その節点の質量で割った値。SLAU の風上化は含まないので**見積り**):

| 壁節点 | 数 | 正味流出の節点 | 流出率/質量 中央値 |
| --- | --- | --- | --- |
| **ベース壁** | 9 | 8 | **1.24e+05 1/s** |
| その他の壁 | 1019 | 566 | 6.02e-01 1/s |

**ベース壁は他の壁の約 2e5 倍**。接線方向の壁では再構成の接線速度が片側で流入・片側で流出になって釣り合うが、
**流れに正対する背面壁 (法線 −x) では全ての内部面が外向き**になり、純粋な吸い出しになる。
これが `nodeWallDirichlet: 0` で止まり (`run_0210`)、1 次で立たず、CV 体積に比例する理由をすべて説明する。

**~~⑤ (肩の丸め) に進む前の判断~~ 撤回 (2026-09-19、SU2 のソースを読んだ結果)**。
当初は「速度ピン壁ノードからの 2 次再構成が異常なので丸めても直らない、ソルバ側の課題」と書いたが、
**再構成は異常ではなかった**。詳細は [`convection-node-wall-reconstruction.md`](convection-node-wall-reconstruction.md) §4.4/§4.5:

- SU2 v8.5.0 も **no-slip は強制 Dirichlet のみ・壁ノードの連続の式は自由・MUSCL に壁の特別扱いなし**で、
  forge と同型。壁ノードの扱いに差は無い。
- forge の再構成値は **隣接ノード速度の半分の 0.70 倍** (62100: 457.4 / (0.5 x 1303.5))。
  片側勾配の理論値 `0.5 x u_隣接` より**小さく**、SU2 が出す値より保守的。`psi` も 1.0 で SU2 の Venkatakrishnan と同じ。
- SU2 の汎用セーフガード (非物理再構成のエッジ 1 次化 / 更新率 20 % 上限 / 点ロールバック) も**この排出を止めない**
  (1 step あたり密度変化率は 4.2 % から育ち、20 % を超えるのは 27 step 中 5 step だけ。再構成は非物理にならない)。

**したがって「ソルバの欠陥だから丸めは筋が悪い」は誤り。codex が指定した順序 (④ 局所格子 → ⑤ 丸めと厚み) が正しい。**

**原因は未確定のまま残す (2026-09-19 再訂正, codex plan レビュー C1/M7)**。
下の「未解像が主因」説は**まだ証明されていない**。測定器の座標を直したところ
(`CELLS/centCoords` は実行時にノード座標へ置換される)、**ベース無しの対照 `run_0209` は 40 step 通して
非物理な再構成がゼロ、ベース run `run_0208` は step 6 から発生** (−66.78 Pa → 3 面) と分かった。
**SU2 の `bad_recon` 検査 (該当エッジを 20 反復 1 次化) はベース run で発火し対照では発火しない**ので、
ソルバ側の対策が効く可能性が残っている。
[`convection-node-wall-reconstruction.md`](convection-node-wall-reconstruction.md) の W2 (面単位フォールバック) の
A/B を**局所格子試験より先に**行うこと。以下は候補として残す仮説:

ベースを横切るノードは 9 個、直後の Δx は**ベース厚の 4.1 倍**で、
ベース圧を決める近傍後流がまったく解像されていない。初期値のベース圧 1159 Pa は外気 2851 Pa の 0.41 倍で
**超音速外流の base pressure として妥当**であり、そこから 32 Pa (真空) へ落ちるのは解像不足の帰結という読み。
~~**1 次が安定なのは一方弁だから**~~ **撤回 (codex Major 2)**: SLAU の質量流束には圧力差の散逸項
`-chi/c_diss*(P_R-P_L)` が入るので、左右の速度が 0 でも質量は動く。1 次が 40 step 不動だったのは
観測事実として残すが、理由づけは取り下げる。
さらに §4.15.7 の格子試験 (`split_plume_at_te`) は **Δx を縮めたが壁 CV も 3.3 倍小さくした**ので、
**後流を解像する試験になっていない**。ベースを横切る点数と後流の解像を**分けて振り直す**こと。

ソルバ側の課題は **リミッタ評価点の不整合 (§4.15.5)** と **面単位の非物理フォールバックの不在**の 2 件で、
[`convection-node-wall-reconstruction.md`](convection-node-wall-reconstruction.md) が引き取った。


**受入条件に足りないもの** (codex M6): `run_0205` は `wallTreatmentSST: 1` = **壁関数**で、等温 no-slip にしただけでは
低 Re 壁解像にならない。`floor_gate` は最終保存場しか見ないので「**持続する**床・更新クリップ無し」を確認できない。
3D の `require_residual_pass` も既定 `False`。受入表に「床・更新クリップの時系列」「ベース圧・再循環長」
「局所壁解像 (`ypls` でなく第一内部点距離と接線壁応力)」を足し、2 次切替の前後で履歴を分けて保存する。




**⚠ 本節の A/B 判定は §4.20 で撤回 (生存 step は再現しない)**

### 4.16 R4e 有限ベースの 2D 起動 — 排出は消え、残るのは**肩の成長振動** (2026-09-20)

**前段**: ベース面からの質量排出は解像度不足が真因で、メッシャの `first_wake_frac` (ベース直後の第一 station を
絶対値で `t_base/5` 以下に強制) で根治した。経緯と根拠は
[`convection-node-wall-reconstruction.md`](convection-node-wall-reconstruction.md) §4.28 / §4.29。

**残った現象は別物**。細分メッシュ + cfl 5 (`run_0292_wake_cfl5`) の NaN 前を 1 step 刻みで追うと、
**単調排出ではなく振動が成長**している (ベース節点 id 62103 の ρ_min):

| step | 8 | 14 | 16 | 20 | 22 | 26 | 28 | 30 | 34 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ρ_min ×1e3 | 3.390 | 3.417 | 3.341 | **3.678** | **3.214** | **3.118** | 2.176 | 2.362 | — |
| P_min [Pa] | 1122 | 1122 | 1122 | 1122 | 1122 | 1122 | 970 | **51.8** | **20 (床)** |
| \|U\|max [m/s] | 2529 | 2529 | 2529 | 2530 | 2530 | 2530 | 2530 | 2530 | **15806** |

振れ幅が step 8 の ±0.5 % から step 26 の ±7 % へ**指数的に成長**し、そこで P が床まで落ちて速度が発散する。
NaN 節点 75 個の分布は x 1.0702–1.0808 / y 0.2694–0.2791 で、**ベース面 (x=1.07732, y=0.27082–0.27282) の
上流側・上下の肩にまたがる**。§4.28 の排出 (ベース面から下流だけ) とは場所も指紋も違う。

**`lowMachPrecond: 2` は効かないどころか悪化させる** (低マッハ市松の指紋に見えるので最初に試した):

| 構成 | 結果 |
| --- | --- |
| `lmp 0`, cfl 5 (`run_0292`) | NaN step 34 |
| **`lmp 2`, cfl 5** (`run_0294`) | NaN step **5** |
| `lmp 0`, cfl 1 (`run_0295`) | **500 step 完走** |
| **`lmp 2`, cfl 1** (`run_0296`) | NaN step **7** |

→ **[[backstep-lowmach-checkerboard-precond2]] の処方はここでは使えない**。別機序である。

**CFL 上限 (細分メッシュ, 500 step, `lmp 0`)**:

| cfl | 1.0 | 1.5 | 2.0 | 3.0 | **4.0** |
| --- | --- | --- | --- | --- | --- |
| 完走 | ○ | ○ | ○ | ○ | **NaN step 79** |
| `rms_roe` (500 step 時点) | **3.61** | 25.9 | 118 | 208 | — |

**粗メッシュは cfl 1 でも 5 でも step 28 だった**ので、細分が初めて CFL 上限を作った。
ただし **CFL を上げるほど残差が悪化**する (roe が 3.61 → 208) ので、上限ぎりぎりを狙う意味は無い。

**暫定の生産設定**: 有限ベース (`t_base > 0`) の 2D は **`opt.cfl_main: 1.0`**
(`problem_r4e_2d_base_wake_cfl1.yaml`)。ベース無しの生産値 5.0 は使えない。
`run_0295` の 500 step は `check_convergence.py` で **NOT CONVERGED (still converging / plateau)** だが
`rms_ro` `rms_roUx` `rms_roe` `rms_roY1` は falling で、**発散ではなく収束途上**である
(粗メッシュはここで既に NaN だった)。長時間 run で収束するかは `run_0301_wake_cfl1_full` で確認中。

**機序 → §4.17 で判明**: 陰解法の反復不安定。起点は肩ではなく**ベース帯**で、(b) の角 BC 競合は無かった。

**⚠ 本節の機序断定は §4.20 で撤回 (implicitRelax の A/B は未実施・生存 step は再現しない)**

### 4.17 肩の振動は~~**陰解法の反復不安定**で、sweep を増やすと悪化する~~ (撤回) (2026-09-20)

§4.16 の振動を 3 節点で 1 step 刻みに追った (`run_0292_wake_cfl5`, Δρ/ρ [%]):

| step | 12 | 18 | 21 | 24 | 27 | 30 | 33 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ベース中央 (62103) | +0.02 | +0.26 | −0.14 | +1.11 | **−4.39** | **+8.84** | **−1832** (ρ<0) |
| 下肩 ramp∩base (44679) | +0.05 | +0.44 | −0.91 | +3.41 | −0.18 | −7.64 | −4.22 |
| 上肩 vehicle (47668) | −0.00 | −0.01 | +0.07 | +0.09 | +0.43 | +1.33 | +2.14 |

**step 15 までは平坦**で、そこから**符号を交互に変えながら約 3 step で倍増**する。
**ベース帯で始まり肩は巻き込まれるだけ** — §4.16 で「肩の現象」と書いたが、**起点はベース帯**である (記述を訂正)。
空間的な市松ではなく**時間方向の奇偶振動** = 反復のスペクトル半径 > 1 の指紋。

**角の BC 競合ではない**: 下肩ノード 44679 は physID 4 (`ramp`) と 10 (`vehicle`) の**両方に属する**が、
どちらも `wall_isothermal` で値も同一 (Ux=Uy=Uz=0, Ts=1000) なので [[node-inlet-wall-corner-conflict]] 型の競合は無い。
上肩 47668 は physID 10 のみだが同じように振動する。

**陰解法ノブの A/B (cfl 5, 300 step)**:

| 変更 | 結果 | `rms_roe` |
| --- | --- | --- |
| 基準 `nStepInner 5, implicitRelax 1.0` (`run_0292`) | NaN step 34 | — |
| **`nStepInner: 10`** (`run_0302`) | NaN step **26** (悪化) | — |
| **`nStepInner: 20`** (`run_0303`) | NaN step **26** (悪化) | — |
| `implicitRelax: 0.7` (`run_0304`) | **300 step 完走** | 4.01e+02 |
| `implicitRelax: 0.5` (`run_0305`) | **300 step 完走** | 4.01e+02 |
| (参考) cfl 1.0, relax 1.0 (`run_0295`, 500 step) | 完走 | **3.61e+00** |

**DPLUR の sweep を増やすと悪化する**のが決定的: 線形解の収束不足ではなく、**近似 Jacobian そのものが
不安定な更新へ収束している**。[[implicit-cfl-ceiling-eos-floor]] の「陰的 defect-correction の反復不安定」と同型。
`implicitRelax` は止められるが、**残差水準は cfl 1 + relax 1.0 の 100 倍悪い** (401 vs 3.61) ので使う価値が無い
([[sern-cfl-no-implicit-relax]] のユーザ方針とも整合)。

**生産設定 (2D 有限ベース)**: `opt.cfl_main: 1.0`, `implicitRelax: 1.0`, `nStepInner: 5`, `lowMachPrecond: 0`。

### 4.18 R4e 有限ベース 2D の現状 — **ベースは解決、起動が残る** (2026-09-20)

`run_0306_wake_cfl1_long` (細分メッシュ + `cfl 1.0` + `implicitRelax 1.0` + `lmp 0`, `run_0295` の場から 12000 step):

- **NaN 無し・12000 step 完走**。粗メッシュは cfl を問わず step 28 で死んでいたので、これが初めての通し。
- **ベース領域は完全に定常**: ベース中央 ρ 0.00515636 → 0.00515634 (有効数字 6 桁で不変)、
  下肩 0.0184513、上肩 0.00620986、`P_min` 1122 Pa、`ro_min` 0.003384 がいずれも step 2000 以降**不変**。
- `check_quasisteady.py` **`OVERALL: ALL STEADY`** (shock / machmax 8.559 / pmax 1.206e5、drift 0.0 %/tail・fluct 0.0 %)。
- `check_convergence.py` **`NOT CONVERGED (stalled/plateau)`**。全列がプラトーで、末尾 2000 step の変動係数は
  `rms_ro` 20.9 % / `rms_roUx` 21.8 % / `rms_roe` 17.3 % = **残差のリミットサイクル**。
  ベース領域が完全に静止している以上、この振動は**別の場所** (プルーム/せん断層) 由来である。

**まとめ**: §4.28/§4.29 の後流細分で**ベース面の排出は根治**し、§4.17 の cfl 1.0 で**陰解法の反復不安定も回避**できた。
**「ベースを付けると数十 step で壊れる」という ⑤ の詰まりは解消**している。

**残る 2 件**:

1. **起動が不安定**。MOC IC からの通し (`run_0301_wake_cfl1_full`) は **soft 段 (1 次 + SST + cfl 1.0) の step 229 で NaN**。
   warm 梯子 (層流 cfl 0.1→0.3→1.0) は通る。収束した場からは 12000 step 安定なので、**壊れているのは起動経路**。
   問題 YAML の既存コメントも「`cfl_main` 1.0 は限界的で、同一メッシュ・同一 config でも run により発散した」と
   記録しており、起動の当たり外れがある。**次はここ** (soft 段の梯子を細かくする / warm を延ばす)。
2. **残差プラトー**。積分量は STEADY なので設計チェーンの目的量は取れるが、
   `require_residual_pass` を立てると落ちる。プルーム側の振動源を特定するまではゲートを積分量で見る。

### 4.19 R4e 有限ベース 2D が**通しで完走** — 生産レシピ確定 (2026-09-20)

§4.18 の残課題「起動が不安定」を潰した。soft 段の CFL が高すぎたのが原因:

| soft の設定 | 結果 |
| --- | --- |
| `soft_cfl: 1.0` 固定 (`run_0301`) | **soft 段 step 229 で NaN** |
| `soft_ramp: [0.2, 0.5, 1.0]` (`run_0307`) | 0.2 と 0.5 の脚は通り、**1.0 の脚で step 87 NaN** |
| **`soft_ramp: [0.2, 0.35, 0.5]`, `soft_cfl: 0.5`** (`run_0308`) | **全段通過** |

**`run_0308_wake_soft05` (`problem_r4e_2d_base_wake_soft05.yaml`) — 通しで完走**:

- 段: warm 0.1 → 0.3 → 1.0 (層流) / soft 0.2 → 0.35 → 0.5 (SST, 1 次) / mid 0.2 → 0.35 → 0.5 (2 次) / 本段 8000 step cfl 1.0。
- 設計チェーンの **`GATES: PASS`** (`fail_class=None`, `reasons=[]`)、目的量 **`C_T_with_shear` = 0.9154138**
  (`C_T` 0.9244435, `C_M` −7.368420)。
- `check_quasisteady.py` **`OVERALL: ALL STEADY`** (81 スナップショット。shock / machmax 8.559 / pmax 1.206e5、drift 0.0 %/tail・fluct 0.0 %)。
- `check_convergence.py` **`NOT CONVERGED (stalled/plateau)`**。全列 1.5〜2.3 桁低下してプラトー。
  §4.18 と同じ**残差のリミットサイクル**で、ベース領域は静止している。

**生産レシピ (2D 有限ベース `t_base > 0`)** — ベース無しの値とは別物なので取り違えないこと:

| 項目 | 有限ベース | ベース無し (従来) |
| --- | --- | --- |
| `mesh.first_wake_frac` | **必須, ≤ `t_base/5`** (無いと生成が失敗する) | 不要 |
| `opt.soft_ramp` / `soft_cfl` | **[0.2, 0.35, 0.5] / 0.5** | 1.0 |
| `opt.cfl_main` | **1.0** (上限 3、4 で発散) | 5.0 |
| `implicitRelax` / `nStepInner` / `lowMachPrecond` | **1.0 / 5 / 0** (変えると悪化: §4.17) | 同じ |

**残件**: 残差プラトーの振動源 (プルーム側) の特定。積分量は STEADY なので目的量は取れるが、
`opt.require_residual_pass` は立てられない。

### 4.20 codex result レビューによる**大幅な撤回** (2026-09-20)

[`notes/reviews/2026-09-20-tooling-nozzle-sern-3d-result.md`](../../notes/reviews/2026-09-20-tooling-nozzle-sern-3d-result.md) (**NO-GO**, Major 7 / Minor 2)。
**根拠は全件自分で確認した。§4.16–§4.19 の機序断定はほぼ全部撤回する。**

#### (1) 「粗細 × CFL の 2×2」は成立していなかった (Major 1)

`run_0269_v1_base` を「粗・cfl 5」として扱ったが、**実際の config は `cfl: 1.0, cfl_pseudo: 1.0`**
(ログも `cfl_pseudo=1`)。比較先の `run_0293_coarse_cfl1` も 1。**同じ条件を比較していた**ので、
排出率がビット一致したのは当然であり、**「粗メッシュでは CFL を変えても変わらない = 時間刻みの不安定ではない」は撤回**する。
`setDT_d.cu:216` は `cfl_pseudo` に比例した局所刻みを設定しており、CFL 非依存という結論自体が誤り。

#### (2) **NaN が出る step は再現しない** — 機序断定の土台が無い

`run_0292_wake_cfl5` と同一入力 (mesh の md5 一致、config 差は `nStepOuter`/`outStepInterval` のみ) で繰り返すと:

| run | `nStepOuter` / `outStepInterval` | NaN |
| --- | --- | --- |
| `run_0292_wake_cfl5` | 40 / 1 | step **34** |
| `run_0323_repro_rep` (同一) | 40 / 1 | step **38** |
| `run_0325_repro_rep2` | 300 / 1 | step **53** |
| `run_0324_repro_out50` | 40 / **50** | **40 step 完走** |

**同じ設定で 34 / 38 / 53 とばらつく**。§4.15 で「生存 step は指標にならない」と自分で棄却したのに、
§4.16–§4.17 の A/B (lowMachPrecond 2 が悪化 / nStepInner 10・20 が悪化 / CFL 上限は 3) を**全部その指標で判定していた**。
**これらは全部保留**とする。`outStepInterval` だけで完走/発散が変わることも未解明 (出力が解を摂動しているのか、
単に上のばらつきの範囲内なのか切り分けていない)。

#### (3) `implicitRelax` の A/B は**実施されていなかった** (Major 2)

`run_0304_sh_relax07` / `run_0305_sh_relax05` の**ログはどちらも `implicitRelax=1`**。
元 config に `implicitRelax` キーが無く、`sed` が空振りしていた。
**「relax 0.7/0.5 なら止まるが残差が 100 倍悪い」は撤回**する (300 step 完走は relax 1 での結果)。

#### (4) `first_wake_frac` は 3D でベース直後を細分していなかった (Major 3)

`mesh_sern3d.py` は `_wake_stations(L_cowl, ...)` と**カウル後縁**を起点に渡しており、
ベース直後 (ランプ後縁 `L_ramp`) は粗いままだった。codex の再計算で **要求 0.004 に対し実測 0.171 = 42.75 倍**。
さらに共有する `_geom_start` は**公比の探索上限 3 で打ち切ってから正規化**するので、
点数が足りないと要求間隔を実現しないまま返る (長さ 2・5 点・要求 0.004 で実測 0.05)。
**入力値の検査では保証にならない。**

**修正 (実装済み)**: 3D も `L_ramp` で区間を分割して**ベース直後から**細分する。
さらに 2D/3D 共通で **生成後の実座標を検査する `check_wake_first_spacing()`** を追加し、
実間隔が要求を超えたら `ValueError` にした。検証: 2D **0.000400 m (t_base の 0.20 倍)**、3D **0.004000 (要求 ≤ 0.004)**。

#### (5) 必須化が既存 3D 機能を壊していた (Major 4)

3D の既定は `t_base=0.02`, `first_wake_frac=0` なので、私のハードエラーで
`design/tests/run_sern_mesh3d_tests.py` が**最初のケースで落ちていた**。
**修正**: 未設定は失敗させず **`t_base/5` を既定**にする (保証は生成後検査で行う)。**テストは `ALL PASS` に復帰**。
`t_base/5`・soft ramp・CFL 上限は **`m6_on` の一形状・一厚みでしか確かめていない**ので、
**「全有限ベースの一般則」ではなく暫定レシピ**と書き直す。

#### 生き残る結果

- **`run_0306_wake_cfl1_long`**: 12000 step 完走、ベース領域の ρ が step 2000 以降**有効数字 6 桁で不変**、
  `check_quasisteady` ALL STEADY。これは生存 step でなく**長時間の定常性**なので有効。
- **`run_0308_wake_soft05`**: 全段通過・`GATES: PASS`・`C_T_with_shear` 0.9154138・ALL STEADY。有効。
- **排出率の測定** (粗 −3.94〜−4.53 %/step vs 細 −0.00 %/step): 連続量の測定なので有効。
  ただし**粗と細で初期密度が 0.00446167 vs 0.00515834 と違う**ので、格子だけの介入ではない (Major 1)。
  **「解像度不足が真因・根治」は「細分を含む条件変更で排出が抑制された」に限定する。**

#### 残作業 (§5.1 へ)

1. 共通の場を正しく移植し、**実効 CFL をログで検査した本当の 2×2** を再実施する。
2. **判定を生存 step から外す** — 発散までの step でなく、同一 step 数での排出率・残差・定常性で比べる。
3. `implicitRelax` の A/B を**実効設定を確認して**やり直す。
4. `outStepInterval` が結果を変えるのかを切り分ける。

### 4.21 撤回後のやり直し設計 (2026-09-20)

§4.20 で A/B の土台が崩れたので、**判定の設計から作り直す**。

#### O1 出力が解を変えるのか (最優先)

`run_0292` (`outStepInterval: 1`) は step 34 で NaN、`run_0324` (同一入力・`outStepInterval: 50`) は 40 step 完走。
**出力頻度が解に影響してはならない**ので、これは (a) 出力経路が状態を触っている実バグ、
(b) 発散が指数成長していて微小差が数 step の違いに増幅されているだけ、のどちらか。

**試験**: 同一 config を 2 回走らせ**最終場をビット比較**する (`outStepInterval` 50 同士、1 同士)。
次に `outStepInterval` 1 と 50 を**同じ step 数**で走らせ、共通の出力 step で場を比較する。
(a) なら「2 回走のビット差 0、かつ 1 vs 50 に差」になる。(b) なら 2 回走にも差が出る。

#### O2 判定指標を生存 step から外す

発散までの step は再現しない (34/38/53)。代わりに**同一 step 数で止めて連続量で比べる**:

| 指標 | 取り方 |
| --- | --- |
| ベース節点の排出率 | Δρ/ρ [%/step] の初期 10 step 平均 (単調排出は再現性が高い) |
| 振動の成長率 | ρ の step 間差分の符号反転回数と包絡の倍加 step |
| 残差 | 同一 step での `rms_*` 全列 |
| 定常性 | `check_quasisteady.py` |

**A/B はすべて「同一 step 数で完走させて上の量を比べる」**。発散するまで回してはいけない。

#### O3 本当の 2×2

粗/細 × cfl 1/5 を、**同じ物理場から**出発させ、**ログで実効 CFL を確認**して回す。
粗と細で初期密度が違った (0.00446167 vs 0.00515834) 問題は、**細メッシュの場を粗メッシュへ
`interp_field.py` で移植**して揃える ([[cross-mesh-restart-interp]])。

#### O4 `implicitRelax` の A/B をやり直す

元 config にキーが無く `sed` が空振りしていた。**キーを明示的に追記**し、**ログの `implicitRelax=` で確認**してから比べる。

### 4.22 やり直した A/B — **機序の結論は戻ってきた。ただし土台は別物** (2026-09-20)

#### O1 出力は無罪。**run が本当に非決定的**だった

| 比較 (30 step 後, `run_0327`–`run_0330_det_*`) | 最大相対差 | 差のあるノード |
| --- | --- | --- |
| 同一 config 2 回走 (出力 30 step ごと) | **1.577e+01** (`Uy`) | 60062 |
| 同一 config 2 回走 (出力 毎 step) | 6.162e-04 (`Uy`) | 60103 |
| 出力 30 ごと vs 毎 step | 2.161e-01 (`Uy`) | 60112 |

**同一 config の 2 回走だけで最大の差が出る**ので `outStepInterval` は原因でない。

**1 step で切り分けた** (`run_0331`–`run_0333_det1`): step 0 (入力) は**ビット一致**、
**step 1 で 1.481e-07〜1.851e-07 の差が 16288〜16403 ノード** (最大は `roY1`)。
これは node の既知のノイズ床 ([[cell-atomicadd-nondeterminism]]: node は ~1e-6) で、
**発散が指数成長してこれを増幅するので生存 step が散る**。§4.20 の 34/38/53 はこれで説明が付く。

#### O2/O4 **同一 30 step・3 反復・実効設定をログで確認**した A/B

`run_0334`–`run_0351_ab_*_r{1,2,3}` (細メッシュ, `run_0292` と同じ場から)。
`implicitRelax` は元 config にキーが無かったので**明示追記し、ログの `implicitRelax=` で確認**した。

| 構成 | ベース節点 \|Δρ/ρ\| 初期 10 step 平均 [%/step] | **振動成長 (末 5 / 初 5)** | `rms_ro` @30 | 30 step 到達 |
| --- | --- | --- | --- | --- |
| base (cfl 5) | 0.0119 ± 4.6e-6 | **2.25e+03 ± 3e+02** | 5.61e-5 | 3/3 完走 |
| **cfl 1** | 0.000475 ± 4.3e-7 | **1.44 ± 0.004** | 4.14e-6 | 3/3 完走 |
| `lowMachPrecond: 2` | — | — | — | **3/3 とも step 5 で NaN** |
| `nStepInner: 20` | 0.0116 ± 1.2e-5 | — | — | **3/3 とも step 26 で NaN** |
| `implicitRelax: 0.7` | 0.000723 ± 8.5e-7 | 2.71 ± 0.010 | 1.46e-5 | 3/3 完走 |
| **`implicitRelax: 0.5`** | 0.000295 ± 8.1e-13 | **0.867 ± 0.008** (**減衰**) | 3.96e-6 | 3/3 完走 |

**反復ばらつきが極めて小さい** (成長率で ±0.4 %) = **これは再現する指標**である。生存 step とは違う。

**§4.16/§4.17 の結論は戻ってきた**が、根拠は入れ替わった:

- **CFL は効く**。cfl 5 は振動が 30 step で **2250 倍**、cfl 1 は **1.44 倍**。
  §4.20 で撤回した「CFL 非依存」は**誤りのまま**で正しい (そもそも `run_0269` が cfl 1 だった)。
- **`lowMachPrecond: 2` は決定的に悪い** — 3 反復とも **step 5**。強い不安定では生存 step も再現する。
- **`nStepInner: 20` も悪い** — 3 反復とも **step 26** (base は 30 step 完走)。
  sweep を増やすほど悪化するので、線形解の収束不足ではないという読みは**支持される**。
- **`implicitRelax` の A/B は今度こそ実施した**。0.7 はまだ成長 (2.71)、**0.5 は減衰 (0.867)** し、
  `rms_ro` 3.96e-6 は **cfl 1 の 4.14e-6 と同等**。
  §4.20 で撤回した「relax で止まる」は**結果としては正しかった**が、当時の run では効いていなかった。

**`cfl 5 + implicitRelax 0.5` は棄却** (`run_0352_long_cfl1ref` vs `run_0353_long_relax05`, 各 4000 step):

| | `rms_ro` | `rms_roUx` | `rms_roe` | `rms_roOmega` |
| --- | --- | --- | --- | --- |
| **cfl 1, relax 1** | **1.37e-6** | **2.23e-3** | **2.44** | **0.136** |
| cfl 5, relax 0.5 | 5.42e-6 | 1.62e-2 | **81.7** | 1.93 |

両者 `check_quasisteady` **STEADY**、`check_convergence` は両者 **NOT CONVERGED (plateau)**。
**30 step の指標では同等に見えたが、4000 step では cfl 1 が全列で 4〜33 倍良い**。
→ **短窓の指標は再現するが、生産設定の順位付けには使えない**。`cfl_main 1.0` を維持する。
(relax は解自体も動かす [[sern-cfl-no-implicit-relax]]。)

### 4.23 codex 相談 (plan 段) の反映と、指標の作り直し (2026-09-20)

[`notes/reviews/2026-09-20-tooling-nozzle-sern-3d-plan.md`](../../notes/reviews/2026-09-20-tooling-nozzle-sern-3d-plan.md) (**GO-with-changes**, Major 6 / Minor 1)。

#### 非決定性の**発生源が特定された** (Major 1)

codex が独立に `run_0292` と `run_0323` を比較し (入力の SHA-256 一致、`res_0` ビット一致、
`res_1` の `ro` が 18811 点・最大 3.73e-8、`roe` が 13613 点・最大 0.0625 で相違)、
**SLAU の残差集積に浮動小数点 `atomicAdd`** があることを示した:

- `convectiveFlux_slau_d.inc.cuh:563` — `atomicAdd(&res_ro[ic0], ...)` ほか流れ 5 変数。**1 スレッド = 1 面**なので
  1 つの CV に複数面が非決定的な順序で加算される。**node でも通る主ループ**である。
- `ransTransport_d.cu:55` — SST 勾配も同じ構造。

**[[cell-atomicadd-nondeterminism]] の「node は決定的」は正確でない** — node も同じ `atomicAdd` を通り、
振幅が cell より小さい (~1e-7/step vs ~6e-4) だけである。

また **「2 回走に差があるから出力バグではない」という私の判別論理は成立しない** (両者は併存しうる)。
結論自体は次で裏づけられた。

#### 出力経路はソース読解でも無罪 (Major 2)

codex がソースで確認: `output.cpp:65` の `outputH5_XDMF` は **D2H コピーのみで EOS を呼ばない**、
`main.cpp:1374` の `dependentVariables` は残差組立側で呼ばれる。
**通常の HDF5 出力から EOS 再適用へ至る直接経路は無い**。
ただし `dependentVariables_d.cu:76,209` は**密度床を適用し TP 経路で `roe` を再構成する**ので、
**診断のために EOS を追加で呼ぶと診断自体が介入になる**。今後の診断は複製した状態の上だけで行う。

#### 指標の作り直し (Major 3) — `%/step` をやめて**質量収支**にする

`dt_local` は CFL と局所状態に依存する (`setDT_d.cu:216`) ので、`%/step` は粗細・CFL 間で**同じ尺度ではない**。
固定ノード集合の **ΣρV** で測り直した (ベース帯 8 節点):

| | 初期 ΣρV [kg/m] | 経過 |
| --- | --- | --- |
| **粗 (cfl 1, `run_0269`)** | 7.972854e-08 | step 4 **94.0 %** → 12 **81.5 %** → 24 **65.4 %** → step 28 **0 % (NaN)** |
| **細 (cfl 1, `run_0306`)** | 1.036875e-08 | step 2000 **99.974 %** → 12000 **99.973 %** (**12000 step 平坦**) |

**正規化不要で結論が言える**: 粗メッシュは 24 step でベース帯の質量を **1/3 失って**消滅し、
細メッシュは **12000 step で 0.027 % しか動かず横ばい**。
(粗/細で開始場が違う点は残る [Major 1/5] が、この差はその感度をはるかに超える。)

#### 未反映のまま残す指摘 (§5.1 へ)

- **Major 3 の残り**: 面 `massflux` の流入/流出別収支と圧力差項の寄与、`dt_local`・DPLUR 補正・EOS 床/壁ピンによる
  保存量変更を**時相ごとに分けて**記録する。固定線形系に対する sweep ごとの線形残差で、線形解法と非線形更新を分ける。
- **Major 4**: `run_0324` は 40 step 走ったが `res_0.h5` しか書いていない (`output.cpp:264` は単純な剰余判定で
  **最終 step の強制保存が無い**)。短窓 A/B では**最終 step の保存を明示する**こと。
- **Major 5**: `interp_field.py` は共通初期場を保証しない。粗/細の比較をやり直すときは移植後の場を検査する。
- **Major 6**: 設定キーだけでなく**バイナリと実効作用素**も固定して記録する。

### 4.24 R5 3D SST を現行バイナリで取り直した — **入口∩側壁の角の格子**が壁 (2026-09-20)

2D の有限ベースが通った (§4.19) ので ④ に進み、**生産 3D YAML (`problem_3d_sst_cycle_m6on.yaml`) を
現行バイナリで**回した (`run_0400_3d_base`, 1161650 節点, **21 ms/step**)。
新しいもの: R4e トポロジ + `first_wake_frac` 自動 (`t_base/5`) + **リミッタ新既定 (`limiterScaled: 1`, K=0.05)** + `inlet_Pressure` 修正。

**結果**: 層流暖機 (cfl 0.2, 2000 step) は通り、**soft 段 (SST, 1 次, cfl 0.2) の step 229 で `roOmega` が非有限**。

**壊れ方が極めて局所的**:

| 量 | 値 |
| --- | --- |
| 非有限ノード数 | **`roOmega` が 1 節点のみ**。`ro`/`P`/`k`/`vis_turb` は**全て有限** |
| その節点 | id 218659 **(0.0809, 0.1208, 0.1004)**、`wall_dist` **4.000e-4** (= `first_wall_frac`)、**`k = 0`** |
| 境界所属 | **無し** (内部節点)。壁第 1 層 |
| `omega` | 最大 **8.4e25** (1 節点のみ > 1e20) |

**初期場 (層流暖機の出力) が既に壊れている**:

| 量 | 値 |
| --- | --- |
| `k` | min **0** / 中央値 480 / p99 1.40e4 / **max 2.06e9** |
| `k > 1e6` の節点 | **ちょうど 1 つ** = **(0.081, 0.122, 0.100)** = **後に ω が発散する節点そのもの** |
| `k ≤ 0` の節点 | **429** (step 229 時点で 617)。**全て `wall_dist ≈ 4.000e-4` = 壁第 1 層**、壁上は 0 |
| その帯 | x 0.054–1.136 / y −0.007–0.272 / **z 0.000–0.100** |

**z = 0.1000 は側壁面**である (`sidewall_in` / `sidewall_out` / `vehicle_side` の z がいずれも 0.1000 固定)。
つまり **k≤0 帯は側壁までの範囲**にあり、発散節点は**側壁のすぐ内側 (z=0.1004)**。

**メッシュ最小セルも同じ場所**: 体積最小の 5 節点はいずれも **z = 0.1000 ちょうど、x 0.0002–0.0011、y ≈ 0** で
体積 **8.4777e-12** (全て同値 = 退化帯)。中央値 4.47e-8 の **5000 分の 1**、発散節点 (4.82e-10) のさらに **57 分の 1**。
= **ノズル入口 ∩ 側壁の角**に退化したセル帯がある。

**発散節点自体の双対 CV は歪んでいない** (隣接との体積比 0.70–1.69) ので、
**その節点が特異なのではなく、上流の角の格子が k の初期場を壊している**とみるのが自然。

**位置づけ**: 3D YAML のコメントにある「`soft_cfl` 既定 0.5 では m10_on が**入口面∩カウル外面の角**で step 3 発散 (run_0041)」
と**同じ場所**である。したがってこれは **R4c (3D 領域トポロジの作り直し) の管轄**であり、
2D で効いた後流細分やリミッタ修正では届かない。

**次の切り分け**: (a) 入口∩側壁の角の格子を直す (退化帯 8.48e-12 の解消)、
(b) SST の初期場を暖機の出力でなく**明示的に自由流値で patch** する (`k` 0〜2e9 は初期値として論外)、
(c) 壁第 1 層で `k ≤ 0` が 429 節点出ること自体の切り分け (SST の壁境界か初期化か)。
**(b) が最も安く、(a) より先に試すべき**。

### 4.25 3D SST の壁は **SST 壁関数の壁モデル渦粘性に上限が無い**こと (2026-09-20)

§4.24 で「入口∩側壁の角の格子」と書いたが、**もっと手前で特定できた**。

#### (1) 保存量 `roK` は **1 step も進む前に**壊れている

暖機の最終場 (`sern.h5`) と soft 段の初期出力 (`res_0.h5`) を比べると **`ro` はビット同一**なのに:

| | `roK` の範囲 | `k≤0` |
| --- | --- | --- |
| 暖機出力 (soft の入力) | **18.21 〜 50.26** (k 54.75〜1.70e4) | **0 節点** |
| soft 段 `res_0` | **0 〜 1.789e7** (k 0〜2.06e9) | 429 節点 |

**31153 節点で `roK` が書き換わる** (壁上は 210 だけ、**内部 25550**)。
`roK > 1e3` の 6994 節点は**全て `wall_dist = 4.000e-4` = 壁第 1 層**。

#### (2) 壁関数を外すと直る

同じ場・同じ設定で `wallTreatmentSST` だけ変えた:

| 設定 | 結果 | `res_0` の `k` |
| --- | --- | --- |
| **`wallTreatmentSST: 0`** (低 Re) | **300 step 完走** | [0, **1.23e4**] (暖機出力の範囲内) |
| `wallTreatmentSST: 1` (壁関数) | **step 37 で NaN** | [0, **2.06e9**] |

**2D 生産 (`run_0308`) も `wallTreatmentSST: 1` で通る**ので、設定ではなく**3D 固有の条件**で壁関数が壊れる。

#### (3) 機構: `ν_t,wall = ν(1/g − 1)` に上限が無い

`ransWallFunction_d.cu:288-294`:

```
omega_w  = max(sqrt(omega_vis² + omega_log²), ωmin),  omega_vis = 6ν/(β₁y²)
mut_wall = ν · max(1/max(g, kSmall) − 1, 0)          g = du⁺/dy⁺ (Reichardt)
k_wf     = omega_w · mut_wall                         ← **上限なし**
```

**失敗節点 (0.0809, 0.1217, 0.1004) の実測**: ρ **8.686e-3**、μ 6.823e-5 → **ν = 7.855e-3**
(低密度域なので海面値の ~400 倍)、`wall_dist` 4.0e-4、ω 441895、k **2.059e9**。
逆算すると **ν_t,wall = k/ω_w ≈ 4660 m²/s**、したがって **g ≈ 1.7e-6** (ほぼゼロ)。

$\nu_t \approx 4660\,\mathrm{m^2/s}$ は物理的にあり得ない (通常 1e-3〜1e-1)。
**`g → 0` で `ν(1/g − 1)` が発散し、低密度で `ν` が大きいほど悪化する。**

#### (4) 提案する対策 (別 plan)

**混合長による上限を掛ける**: 壁モデルの渦粘性は本来 $\nu_t \le \kappa u_\tau y$ なので、

$$\nu_{t,\rm wall} = \min\!\left[\nu\left(\frac{1}{g}-1\right),\ \kappa\,u_\tau\,y\right]$$

とすれば `g → 0` でも発散しない。**任意のクリップではなく壁法則そのものの上限**なので筋が通る。
あわせて **`k_wf` が異常に大きい節点数を診断で数える** (今回は 1 節点でも run を殺した)。

**これは乱流モデルの変更なので別 plan を立てる** → §5.1 R5c。
**当面の回避**: `wallTreatmentSST: 0` は回るが、このメッシュは低 Re 解像ではない (第一層 `first_wall_frac` 0.004) ので
**物理は正しくない**。回避であって解ではない。

### 4.26 壁解像 (y⁺≈1) メッシュの見積もり — **実現可能** (2026-09-20)

**方針変更**: `wallTreatmentSST: 1` (SST 壁関数) は**使わない**ことになった (ユーザ指示 2026-09-20)。
ソルバ既定も 0 に変え、指定すると警告が出る ([`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) の壁処理の節)。
したがって **§4.25 の R5c (壁モデル渦粘性の上限) は優先度が下がり**、代わりに**壁解像メッシュ**が要る。

**現状の壁解像** (`check_wall_resolution.py`, `run_0401_3d_wf0` step 300):

| 壁 | y₁ [m] | y₁⁺ 平均 | p99 | 最大 |
| --- | --- | --- | --- | --- |
| `sidewall_in` | 4.00e-4 | **38.3** | 48.2 | 48.4 |
| `sidewall_out` | 4.00e-4 | 11.2 | 32.5 | 35.5 |
| `vehicle_base` | 4.00e-4 | 16.4 | **241** | 391 |
| `vehicle_side` | 4.00e-4 | 12.0 | 31.3 | 33.6 |
| **`vehicle_top`** | 2.00e-3 | **87.0** | **589** | 908 |

**`VERDICT: FAIL`** — 目標 y⁺ 1 を超える面積が **100 %** (許容 2 %)。

**必要な細分と、その代償**:

| 指標 | 必要な y₁ | 細分率 | 伸長比 1.15 で必要な層数 | **AR 中央値** |
| --- | --- | --- | --- | --- |
| 現行 | 4.00e-4 | — | — | **5.3** |
| 平均 y⁺ を 1 に (38×) | 1.04e-5 | 38× | **+24 層** | **203** |
| `vehicle_top` 平均 (87×) | 2.30e-5 | 87× | **+32 層** | **464** |
| `vehicle_top` p99 (589×) | — | 589× | **+46 層** | **3145** |

(接線方向の刻みは実測で中央値 2.136e-3 m。AR = 接線/法線)

**結論: y⁺≈1 は実現可能**。最悪ケース (p99) でも **AR 中央値 3145** で、
壁法線構造格子に認められている **AR ≤ 5000** ([[mesh-quality-gate]], 2026-09-12 ユーザ決定) の内側に収まる。
**壁関数に逃げる必要はない。**

**規模の見積もり**: 現行 1161650 節点。壁法線 (j) に +30 層、側壁も壁なので z 方向 (`nz_in` 25 / `nz_out` 15) にも同等の層が要る。
j 方向が 80 → ~110 (×1.4)、z 方向が 40 → ~65 (×1.6) で **概算 2.3M 節点**。
[[heavy-postproc-on-aws]] のとおりローカルは 2.7M 節点の変換で OOM するので、**変換と計算は AWS で行う**。

**残作業** → §5.1 R5d。

### 4.27 3D のカウル板厚が効いていなかった (厚さ 0 スリット) — **修正 (2026-09-20)**

§4.26 の壁解像メッシュを組む前に、`check_wall_resolution.py --geometry-only` (本日追加) で
**生成後の実座標**から第一層厚を測ったところ、`cowl_out` だけが `first_wall_frac` の **2.24 倍**
(8.966e-04 m = 4.0e-04 + 5.0e-04) になっていた。内訳の 5.0e-04 m は `cowl_thickness` (0.005 H) そのもの。

**原因**: [`mesh_sern3d.py`](../../design/forge_design/meshing/mesh_sern3d.py) の y 分布で

```python
Y[i, :njb] = y_bot + s_bot * (lo - y_bot)     # j = 0..jm、最後 (j=jm) は lo = ym - t/2
Y[i, jm:]  = up    + s_top * (yt[i] - up)     # j = jm..NJ-1、先頭 (j=jm) を up = ym + t/2 で**上書き**
```

`j = jm` は**下バンドの壁ノード (`cowl_out`)** なのに、2 行目が `up` で上書きしていた。
結果、`base(i,jm,k)` (cowl_out) も `dup1` (cowl_in) も **`ym + t/2`** に置かれ、
板厚 `cowl_thickness` は**まったく効かず**、上下の壁ノードが同一座標の**厚さ 0 スリット**になっていた。

**実測 (生産メッシュ `problem_3d_sst_cycle_m6on.yaml`, 1161650 節点)**:

| 量 | 修正前 | 修正後 |
| --- | --- | --- |
| `cowl_in` y − `cowl_out` y (中央値) | **0.0000e+00** | **5.0000e-04** (= `cowl_thickness`×H) |
| 同 (最小) | 0 | 0 (TE と z=W/2 で設計どおり閉じる) |
| `cowl_out` の第一層厚 | 8.966e-04 | 4.0e-04 (= `first_wall_frac`) の見込み |
| `run_sern_mesh3d_tests.py` | ALL PASS | ALL PASS |

**2D は無罪**: [`mesh_sern.py`](../../design/forge_design/meshing/mesh_sern.py) は下バンドを `ym_lo` で終え、
上側は `dup` ノードに `ym_up` を置くので板厚は正しく出ている。**3D だけの欠陥**。

**修正の中身** (2 点):

1. `Y[i, jm] = lo` を最後に**明示代入**する。順序を入れ替えるだけでは
   `y_bot + 1.0*(lo - y_bot)` の丸めで `lo` と一致せず、厚さ 0 の設計区間で双子ノードが
   ~1e-16 だけ離れる。変換器は座標を **float32** で書くので、この隙間は解けずに
   「衝突するかどうかが丸め次第」になり、`run_sern_mesh3d_tests.py` の
   「float32 化で新たな節点衝突が起きない」が落ちる (実際に落ちた: 27882 vs 28008)。
2. 板の半厚が float32 で解けない列 (`0.5·t·sz < 1e-5·max(1,|y_m|)`) は**厳密に 0 に丸める**。
   設計どおり閉じるところは厳密に閉じる ([[axisym-rweight-closure-fp32]] と同じ方針)。

**効いた範囲**: `cowl_thickness > 0` を指定した**すべての 3D run** (R5 系の SST 試験、`run_0400`/`run_0401` を含む)。

**ただし R5 の真因ではなかった (§4.28 の A/B)**。2D で厚さ 0 スリットが発散原因だった実績
(`run_0035`) から「3D SST が立たない件の候補要因」と見たが、**同一バイナリ・同一 config で
板厚バグのあるメッシュの対照 (`run_0404`) も完走した**。板厚が答えに効く量は `C_T_with_shear` で
**0.003 %** (0.8563184 vs 0.8563473)。**幾何のバグとしては直すべきだが、R5 の説明にはならない**。

### 4.28 R5 が立たなかった真因は **SST 壁関数** — 板厚ではない (A/B, 2026-09-20)

`run_0403` (板厚修正) と `run_0404` (板厚バグのまま、**同一バイナリ・同一 config・同一 20000 step**) の A/B:

| | `run_0403_3d_cowlfix` | `run_0404_3d_cowlcollapse_ctrl` (対照) |
| --- | --- | --- |
| カウル板厚 | 5.0000e-04 m (設計値) | **0** (双子ノードが同一座標) |
| 完走 | **20000 step 完走** | **20000 step 完走** |
| `GATES` | **PASS** (fail_class=None) | **PASS** (fail_class=None) |
| `C_T_with_shear` | 0.8563184 | 0.8563473 (**差 0.003 %**) |
| `check_convergence.py` | **NOT CONVERGED (stalled/plateau)** 全列 0.1–0.3 dec 横ばい | **NOT CONVERGED (stalled/plateau)** 同左 |
| `check_quasisteady.py` | **ALL STEADY** (shock/machmax/pmax とも drift 0.0 %) | **ALL STEADY** |
| NaN/Inf | 残差全列 0・最終場の全変数 0 (ρ min 1.66e-3, P min 863 Pa, T 158–3060 K) | — |

**したがって R5 (soft 段 step 229 で `roOmega` 非有限、`k=2.06e9`) を解いたのは板厚ではない**。
config を突き合わせると:

| run | `wallTreatmentSST` | 結果 |
| --- | --- | --- |
| `run_0400_3d_base` | **1** (生成器が直書きしていた) | soft 段 step 229 で発散 (§4.24) |
| `run_0403` / `run_0404` | **0** | 20000 step 完走 |

**= §4.25 で名指しした「壁モデル渦粘性 $\nu_{t,w}=\nu(1/g-1)$ に上限が無い」がそのまま R5 の壁だった**。
壁関数を使わない方針 (2026-09-20 ユーザ指示、生成器・ソルバ既定とも 0 に変更) が、結果として R5 を通した。
[[sst-wall-function-banned]]。

**残る 2 つ** (どちらも壁関数とは別問題):

1. **未収束**: 全残差が 0.1–0.3 dec で横ばい。本段は 2 次 + リミッタ有効なので
   **リミッタのチャタリング**が候補 (2D で `limiter: 0` なら 5.11e-7 まで落ちることを確認済)。
   段階起動なので低下桁数が小さく出るのは当然だが、**横ばいであること自体**は別途詰める。
2. **壁解像 FAIL**: y₁⁺ 平均 12.7–82.3 (§4.26.1 を step 20000 の実測で更新)。→ R5d。

**`GATES: PASS` は収束の根拠にならない** (この YAML は `require_residual_pass` 未設定 = false)。
R5 の受入では `require_residual_pass: true` を要求する (R5b の残作業)。

### 4.26.1 壁解像の実測 (全 8 壁) と必要な細分 — §4.26 の表の補完

§4.26 の表は `cowl_in` / `cowl_out` / `ramp` を落としていた。全 8 壁は次のとおり
(`run_0401_3d_wf0` step 300、`VERDICT: FAIL`、目標超過面積 最大 100 %):

| 壁 | y₁ [m] | y₁⁺ 平均 | p99 | 最大 | y₁⁺/y₁ [1/m] | 支配する knob |
| --- | --- | --- | --- | --- | --- | --- |
| `cowl_in` | 3.985e-04 | 38.1 | 46.2 | 144.2 | 9.56e4 | `first_wall_frac` |
| `cowl_out` | 8.966e-04 | 36.9 | 54.2 | 144.2 | 4.11e4 | 同 (§4.27 で 4.0e-04 になる) |
| `ramp` | 3.864e-04 | 17.7 | 38.9 | 45.7 | 4.58e4 | 同 |
| `sidewall_in` | 4.000e-04 | 38.3 | 48.2 | 48.4 | 9.58e4 | `first_z_frac` |
| `sidewall_out` | 4.000e-04 | 11.2 | 32.5 | 35.5 | 2.79e4 | 同 |
| `vehicle_side` | 4.000e-04 | 12.0 | 31.3 | 33.6 | 3.00e4 | 同 |
| `vehicle_top` | 2.000e-03 | 87.0 | 588.8 | 907.9 | 4.35e4 | `first_top_frac` |
| `vehicle_base` | 3.999e-04 | 16.4 | 241.1 | 391.3 | 4.11e4 | `first_wake_frac` |

**step 300 の場で測った上の値は過渡**である。`run_0403` の **step 20000 (ALL STEADY)** で測り直すと:

| 壁 | y₁ [m] | y₁⁺ 平均 | p99 | y₁⁺/y₁ [1/m] |
| --- | --- | --- | --- | --- |
| **`cowl_in`** | 3.985e-04 | **57.7** | 84.3 | **1.449e5** (y 方向の律速) |
| `cowl_out` | 3.985e-04 | 22.4 | 31.6 | 5.62e4 |
| `ramp` | 3.864e-04 | 25.0 | 80.3 | 6.48e4 |
| **`sidewall_in`** | 4.000e-04 | **64.2** | 86.3 | **1.605e5** (z 方向の律速) |
| `sidewall_out` | 4.000e-04 | 14.0 | 39.4 | 3.50e4 |
| `vehicle_base` | 3.999e-04 | 12.7 | 214.5 | 3.19e4 |
| `vehicle_side` | 4.000e-04 | 15.4 | 36.7 | 3.86e4 |
| `vehicle_top` | 2.000e-03 | 82.3 | 706.1 | 4.11e4 |

平均 y₁⁺ を 0.7 に落とすのに要る第一層厚と、伸長比 1.15 で要る層数 (**step 20000 実測ベース**):

| knob | 現行 (/H) | 目標 (/H) | 比 | 追加層数 | 効く格子数 |
| --- | --- | --- | --- | --- | --- |
| `first_wall_frac` | 4.0e-03 | **4.8e-05** | 82.8 | +32 | `nj_bot` 31→63, `nj_top` 49→81 |
| `first_z_frac` | 4.0e-03 | **4.4e-05** | 91.7 | +32 | `nz_in` 25→57, `nz_out` 15→47 |
| `first_top_frac` | 2.0e-02 | **1.7e-04** | 117.6 | +34 | `nj_ext_top` 41→75 |
| `first_wake_frac` | 4.0e-03 | **2.2e-04** | 18.3 | +21 | `ni_plume` 130→151 |

**規模**: 主ブロック 255×143×103 = 3.76M、`ext_top` ~2.0M、`vehicle_side` ~0.2M で **概算 6.0M 節点**
(§4.26 の 2.3M は `first_top_frac` と `first_wake_frac` を勘定に入れていなかったので**過小**)。
接線刻み中央値 2.136e-03 m に対し AR は 30〜300 程度で、壁法線構造層に認められた 5000 の内側。
**`.msh` は 1.5 GB 級**になるので生成・変換とも AWS で行い、ホスト RAM (g5.xlarge 16 GB) が足りるかを
先に確かめる。足りなければ `vehicle_top` を slip に戻して `ext_top` を現行のままにする案があるが、
R4f (slip と等温壁の差を同一幾何で測る) と交絡するので**最後の手段**とする。

### 4.29 R5f 本段プラトーの切り分け — 第一便は**設計が悪く**空振り (2026-09-20)

`run_0403` の step 20000 の場から**本段だけ**を 10000 step 回す A/B (`run_staged(run_dir, "main")`)。
第一便 (`run_0405`–`run_0407`) は**全滅**した。いずれも step 6–7 で `roOmega`→`roe` が非有限:

| run | 変更点 | 結果 |
| --- | --- | --- |
| `run_0405_main_nolim` | `limiter: 2` → **0** | step 7 で `ro` 非有限 |
| `run_0406_main_cfl2` | `cfl_pseudo` 0.5 → **2.0** | 同上 |
| `run_0407_main_cfl5` | `cfl_pseudo` 0.5 → **5.0** | 同上 |

**restart は無罪**: `run_0405/res_0.h5` と `run_0403/res_20000.h5` を突き合わせると
`ro`/`roUx`/`wall_dist` は**完全一致**、`roe` 5.4e-8、`roK` 5.8e-5、`roOmega` 1.2e-3 の相対差
(= res 書き出し前の従属変数再評価ぶん)。`vis_turb` だけ相対 1.0 でずれるのは既知の
restart 非忠実性 ([[forge-sst-restart-nonfidelity]])。

**私の設計ミス**: `limiter: 0` は [`procedures/solver-settings.md`](../../procedures/solver-settings.md) で
**原則使用しないと明記されている**設定で、M6 外部流 + 衝撃の 2 次では発散して当然。
2D で `limiter: 0` が 5.11e-7 まで落ちたのは別ケースの診断であって、3D に持ち込める対照ではない。
**これら 3 run は破棄**した (`case/46.sern_design/run_0405_main_nolim` / `run_0406_main_cfl2` /
`run_0407_main_cfl5` — AWS 上で削除済)。

**ついでに分かったこと**: `cfl_main: 0.5` は壁関数を切った後も**上げられない** (2.0 で 7 step)。
YAML のコメントは「1.0 は限界的」と壁関数込みの経験で書かれていたが、壁関数を外しても
上限は動いていない。**CFL の壁は壁関数とは別要因**。

**第二便の設計** (`run_0408`–`run_0411`、同じ場から 10000 step):

| run | 変更点 | 見たいこと |
| --- | --- | --- |
| `run_0408_main_ref` | なし (cfl 0.5 / K 0.05) | restart 込みの基準プラトー |
| `run_0409_main_cfl025` | `cfl_pseudo` → 0.25 | **プラトーが dt に比例するか** (比例ならリミッタ/限界周期のチャタリング) |
| `run_0410_main_k05` | `venkatK` 0.05 → 0.5 | 緩いリミッタでプラトーが下がるか |
| `run_0411_main_k0005` | `venkatK` 0.05 → 0.005 | 厳しいリミッタで上がるか |

### 4.32 R5d: **側壁を y⁺≈1 にすると AR が破綻する** — 段階 1 は y 法線壁のみ (2026-09-20)

§4.26.1 の knob (`first_wall_frac` / `first_z_frac` とも 7.0e-05) でメッシュを作ると
**生成は通るが品質ゲートで落ちる**:

| | `run_0412_3d_wallres` (全方向) |
| --- | --- |
| 節点 / セル | 5,511,173 / 5,382,864 |
| AR | **max 32704 / p99 3956 / >1000 が 7.70 %** |
| skew | max **0.963** / >0.90 が 0.33 % |
| **VERDICT** | **FAIL** (`prepare` が例外で停止) |

**原因は z 分布がテンソル積であること**。`first_z_frac` は `zs` を 1 本だけ作って**全 (i, j) で共有**するので、
側壁 (`sidewall_in` は x ≤ `L_sw` = 0.8 H にしか無い) のための細層が**領域全体を貫く**。
そこに遠方の粗いセル (下バンドの y 刻みは最大 2.3 H) が隣接して AR が跳ねる。
**局所的に測って確認**:

| 変種 | `first_wall_frac` | `first_z_frac` | 節点 | AR max | >5000 |
| --- | --- | --- | --- | --- | --- |
| probA | 1.5e-04 | 4.0e-03 (現行) | 1.92M | **5125** | 0.004 % (76 セル、全て x > `L_ramp`) |
| **stage1** | **1.6e-04** | 4.0e-03 (現行) | **1.92M** | **4819** | **0 %** (p99 995) |
| probB | 1.5e-04 | 1.5e-04 | 4.29M | 15473 | 0.211 % |
| probC | 1.6e-04 | 1.6e-04 + 遠方を細分 (`nj_bot` 65) | 4.52M | 12083 | 0.139 % |

AR > 5000 のセルの**最短辺は probB/probC では z**で、位置は **z = 0.0998–0.1001 = 側壁**。
`nj_bot` を増やして下バンドの粗さを抑えても 12083 までしか下がらない
(遠方の最粗セルが 1.93 H 残るため)。**試した分布・規模 (probA–probC) では側壁の y⁺≈1 と AR ≤ 5000 は両立しなかった**。**「現行トポロジでは不可能」とまでは言えない** (codex plan-2 M6): AR は最長辺/最短辺なので、遠方の長辺をさらに割れば下がる余地がある。その費用を測っていない。

**したがって段階を分ける**:

- **段階 1 (`run_0413_3d_wallres_y`, `problem_3d_sst_cycle_m6on_wallres_y.yaml`)**:
  **y 法線壁だけ**を解像する。`first_wall_frac` 4.0e-03 → **1.6e-04**、`first_top_frac` 0.02 → **6.0e-04**、
  `first_wake_frac` → 1.7e-04、`nj_bot` 31→55 / `nj_top` 49→73 / `nj_ext_top` 41→57 / `ni_plume` 130→153。
  **1,919,057 節点・AR max 4819 (>5000 が 0 セル)**。予想 y₁⁺: `cowl_in` 2.3 / `ramp` 1.0 / `cowl_out` 0.9 /
  `vehicle_top` 2.5、**側壁は 14–64 のまま**。
  → **プラトーが下がるかを測る**。下がれば壁解像説が裏づけられ、残る分が側壁の寄与になる。
- **段階 2 (未着手)**: 側壁の解像は**受理の必須条件に残す** (codex plan-2 M6: 「場の変化が側壁に集中しない」は側壁の離散化誤差が小さい証拠ではない)。手は (i) z のクラスタを側壁の範囲に局所化する構造変更、(ii) 遠方の長辺を割って AR を下げる、の 2 つ
  (テンソル積の `zs` を (i, j) 依存にする、または側壁近傍だけ別ブロックにする)。
  費用対効果は段階 1 の結果を見てから決める。

**段階 1 の投入は一度止まった — SERN チェーンに `mesh.ar_max` が無かった** (2026-09-20)。
AR max 4819 / p99 995 / **skew max 0.466 (>0.90 は 0 セル)** = AGENTS.md が認めた
「壁法線に沿った構造格子の境界層セルは AR ≤ 5000 まで可」(2026-09-12 ユーザ決定) そのものなのに、
`runner_sern3d.py` / `runner_sern.py` は `check_mesh_quality.py` を**既定の `--ar-max 1000` で呼んでいた**
(`runner_axismach.py` / `runner_wt.py` は既に `p.mesh.get("ar_max", 1000)` を渡していた)。
両 runner に同じ knob を足し、`problem_3d_sst_cycle_m6on_wallres_y.yaml` に
**`mesh.ar_max: 5000` と「AR 緩和 (≤5000)」の根拠**を明記した。既定 1000 は変えていないので既存 run に影響は無い。

**`run_0412_3d_wallres` は破棄**(`sern.msh` / `sern_qc.h5` を削除し `MESH_QUALITY.txt` だけ残した)。

### 4.31 R5f プラトーの正体は **近壁** — リミッタでも CFL でもない (2026-09-20)

`run_0403` の step 20000 の場から本段だけ 10000 step、4 通り (`run_0408`–`run_0411`):

| run | `cfl_pseudo` | `venkatK` | `rms_ro` 終値 | 基準比 | `check_convergence` |
| --- | --- | --- | --- | --- | --- |
| `run_0408_main_ref` | 0.5 | 0.05 | 3.26e-06 | 1.00 | NOT CONVERGED (stalled/plateau) |
| `run_0409_main_cfl025` | **0.25** | 0.05 | 2.41e-06 | **0.74** | 同左 |
| `run_0410_main_k05` | 0.5 | **0.5** | 3.15e-06 | **0.97** | 同左 |
| `run_0411_main_k0005` | 0.5 | **0.005** | 3.30e-06 | **1.01** | 同左 |

**リミッタ説は否定された**: `venkatK` を **100 倍 (0.005 → 0.5)** 振ってもプラトーは **4 %** しか動かない。
**単純な dt 比例のチャタリングでもない**: CFL を半分にしても 0.74 倍 (比例なら 0.5 倍)。
**「2D で `limiter: 0` が 5.11e-7 まで落ちた」からの類推は 3D では成り立たなかった** (§4.29 でその対照自体が発散)。

**場の動きは壁に集中している**。`run_0408` の step 5000 → 10000 (5000 step) の ρ 相対変化を
`wall_dist` の帯ごとに見ると:

| 帯 | 節点数 | rel 中央値 | p99 |
| --- | --- | --- | --- |
| **壁第 1 層** (`wall_dist` < 6e-4) | 47566 | **4.78e-05** | 4.15e-03 |
| 近壁 (6e-4 – 5e-3) | 121051 | 4.64e-05 | 1.75e-03 |
| 中間 (5e-3 – 1e-1) | 484392 | 1.44e-05 | 6.61e-04 |
| 遠方 (> 1e-1) | 508641 | **1.28e-06** | 1.85e-04 |

遠方から壁第 1 層で **37 倍**。変化上位 0.1 % の節点の **52.9 % が壁第 1 層**にある
(壁第 1 層は全節点の **4.1 %** しかない = **13 倍の濃縮**)。最大変化点は
(0.1924, 0.1515, **0.0996**) — **z = 0.0996 は側壁 (W/2·H = 0.1) のすぐ内側**で、
`sidewall_in` は y₁⁺ 平均 **64.2** = 8 壁で最悪の解像度 (§4.26.1)。

**結論**: プラトーは**壁解像の問題**とみるのが最も自然。低 Re SST の壁 ω は
$\omega_w = 60\nu/(\beta_1 y_1^2)$ で**第一節点が粘性低層にある**ことを前提にしており、
y₁⁺ 13–82 ではその前提が成り立たない。**R5d (壁解像メッシュ) は精度の話だけでなく、
収束の律速そのもの**であり、優先度を上げる。

**残る確認**: 壁解像メッシュで同じ設定を回してプラトーが下がるか (R5d の受入条件に入れる)。
下がらなければ、次は近壁の非定常性 (側壁角の渦) と float32 ノイズ床を分ける。

### 4.30 R5d 壁解像メッシュ: 生成は通り、**変換がホスト RAM で落ちた** (2026-09-20)

§4.26.1 の knob で `problem_3d_sst_cycle_m6on_wallres.yaml` を組み、
`case/46.sern_design/run_0412_3d_wallres` で `--prepare-only`:

| 段 | 結果 |
| --- | --- |
| 形状・格子生成 | **成功**。`sern.msh` 649 MB、**5,511,173 節点 / 5,638,114 要素** (見積もり 6.0M に対し 5.5M)。約 4 分 |
| `convertGmshToForge` (QC 用 cell 変換) | **OOM kill** (`anon-rss` 8.58 GB / `total-vm` 17.1 GB) |

**箱**: AWS g5.xlarge = **RAM 15 GB**。落ちた時点で R5f のスイープ (1.16M 節点) と
**別セッションの `~/forge49` の計算**が同居しており、空きは 8 GB だった。

**GPU は足りる**: 1.16M 節点の 3D SST が **1830 MiB** なので 5.5M で**約 9.5 GB**、A10G 23 GB の内側。
**律速はホスト RAM (変換) だけ**。

**注意 (node は変換が 2 回走る)**: `prepare` は QC 用に `discretization: "cell"` で 1 回、
本番の node で 1 回、**計 2 回** `convertGmshToForge` を呼ぶ (`runner_sern3d.py` の
`convert_mesh(run_dir, "sern.msh", "sern_qc.h5")` → `check_mesh_quality` → `convert_mesh(..., MESH)`)。
8.6 GB のピークが 2 回来る。

**対処の順** (安い順):

1. **直列化して単独で変換する** (追加コスト 0)。他の計算を止めてから `--prepare-only` をやり直す。
   単独なら空きが ~13 GB になるので 8.6 GB は通る見込み。
2. スワップを足す (ディスク残 15 GB なので 6–8 GB が上限。`sern_qc.h5` + `sern.h5` で ~4 GB 要る)。
3. インスタンスを 32 GB 級に上げる (**ユーザ判断**)。
4. メッシュを削る (`vehicle_top` を粗いままにすると `ext_top` の ~2M 節点が浮くが、
   R4f (slip と等温壁の差を同一幾何で測る) と交絡するので最後の手段)。

### 4.33 段階 1 の結果: 壁解像はプラトーを**改善するが解消しない** (2026-09-20)

**メッシュ**: `run_0413_3d_wallres_y` 系、1,919,057 節点、**`check_mesh_quality` VERDICT: PASS
(AR max 4819 / p99 995 / skew max 0.466, `--ar-max 5000`)**。

**達成した壁解像** (`check_wall_resolution.py`, `run_0414` step 20000):

| 壁 | y₁ [m] | y₁⁺ 平均 | 粗メッシュ (`run_0403`) |
| --- | --- | --- | --- |
| `ramp` | 1.546e-05 | **0.700** | 25.0 |
| `cowl_out` | 1.594e-05 | **0.846** | 22.4 |
| `cowl_in` | 1.594e-05 | **1.608** | 57.7 |
| `vehicle_top` | 6.002e-05 | **2.354** | 82.3 |
| `sidewall_in` | 4.000e-04 | 50.9 | 64.2 (**段階 1 では未解像**) |
| `vehicle_side` | 4.000e-04 | 16.1 | 15.4 (同上) |
| `vehicle_base` | 1.693e-05 | 59.4 (p99 1013) | 12.7 |

`vehicle_base` はベース角の**幾何的特異点**で traction が格子収束しない (ツールが最大値で判定しない理由そのもの)。
全体 VERDICT は側壁が残るので FAIL。

**壁解像メッシュは本段の CFL 上限を下げる**: `cfl_main 0.5` は **step 97 で発散**
(`rms_roOmega` が step 94 に 5.28e+08 → 96 に 7.01e+19 と先行、既知の ω 先行指紋)。
2 次自体は mid 段 (cfl 0.2) で 2000 step 通っているので、次数でなく CFL。**0.25 で完走**。
壁第 1 層を 4.0e-04 → 1.6e-05 m (25 倍薄く) した剛性増と整合する。

**公平な比較 (同じ cfl 0.25・同じ 40000 step・どちらも落ち着いた場から継続)**:

| 残差列 | 壁解像 `run_0415_wallres_y_cont40k` | 粗 `run_0416_coarse_cfl025_cont40k` |
| --- | --- | --- |
| `rms_ro` | **0.4 dec** (2.11e-06 → 1.04e-06) | 0.1 dec (2.16e-06 → 2.55e-06) |
| `rms_roUx` | **0.4 dec** | 0.1 dec |
| `rms_roUy` | **0.8 dec** | 0.1 dec |
| `rms_roUz` | 0.2 dec | 0.2 dec |
| `rms_roe` | **0.5 dec** | 0.1 dec |
| `rms_roK` | **1.2 dec (falling)** | 0.3 dec |
| `rms_roOmega` | **0.5 dec** | 0.1 dec |
| `check_convergence` | **NOT CONVERGED (stalled/plateau)** | **NOT CONVERGED (stalled/plateau)** |
| `check_quasisteady` | **ALL STEADY** | **ALL STEADY** |

**メッシュ間で残差の絶対値は比較できない** (DOF 数も dt も違う) ので、比べたのは**低下桁数とトレンド**。

**結論**: **壁解像は効くが足りない**。y 法線壁を y⁺ 0.7–2.4 にすると低下桁数が 0.1 → 0.4–1.2 dec に増え、
`rms_roK` はまだ落ちている。しかし**両者とも `NOT CONVERGED (stalled/plateau)` のまま**で、
「プラトー = 壁解像」という §4.31 の読みは**部分的にしか正しくなかった**。

**側壁主犯説も支持されない**。段階 1 後に場の動き (step 12000 → 20000 の ρ 相対変化) を測ると:

| 帯 (`wall_dist`) | 節点数 | rel 中央値 |
| --- | --- | --- |
| 壁第 1 層 (< 3e-5) | 42104 | 6.29e-04 |
| 近壁 (3e-5 – 6e-4) | 163663 | 6.41e-04 |
| 中間 | 960410 | 7.7e-05 / 2.3e-05 |
| 遠方 (> 1e-1) | 752880 | 5.89e-07 |

壁近傍と遠方の比は **1068 倍**で依然「近壁が律速」だが、**変化上位 0.1 % のうち側壁帯
(z = 0.0995–0.1005、全節点の 8.1 %) は 0.3 % しかない** (むしろ枯れている)。
中心は `wall_dist` ≈ 4.7e-03 m (0.047 H)、x 0.12–1.34 m = **ランプ上の境界層〜プルームのせん断層**。

**局所ピーク量には実際の格子感度がある** (どちらも `ALL STEADY` だが値が違う):
`machmax` **6.09** (壁解像) vs **7.08** (粗)、`pmax` **1.27e+05** vs **1.42e+05**。
= 粗メッシュの局所ピークは格子収束していない。

**残る候補と次の一手** (§5.1 R5g):

1. **せん断層の真の非定常性**。定常擬似時間では減衰できない小スケール変動なら、残差は原理的に下げ止まる。
   `check_quasisteady` が積分・派生量で STEADY なのと矛盾しない。
2. **float32 / `atomicAdd` のノイズ床** ([[cell-atomicadd-nondeterminism]]: node は ~1e-7/step)。
   プラトーは `rms_ro` 1e-06 前後で、桁が近い。**最も安い判別**は同一入力の反復で、
   [[noise-floor-both-sides]] のとおり**片側 3 反復では足りない**ので両側から取る。
3. 側壁の解像 (段階 2)。上の局在から**優先度は下げる**。

### 4.34 R5g: **ノイズ床ではない**。ただし「リミットサイクル」は**撤回** (2026-09-20)

> **撤回 (同日, codex plan-2 レビュー M3 を採用)**: 本節の指標は **4000 step 離れた場の差の L2 ノルム**で、
> **符号と位相を落とすので単調ドリフトでも一定値になる**。周期軌道への到達も振幅も示せない。
> 「リミットサイクル」「到達限界は残り 1 桁」「受理条件を残差から力の平均±振幅へ置き換える」は
> **いずれも撤回**する。残るのは **(i) 場の変化は 5 反復の再現性幅の 3〜62 倍ある**、
> **(ii) `ro`/`roUy`/`roOmega` はこの指標では減衰が止まって見える**、**(iii) 力は 7 桁で定常**の 3 点だけ。
> 5 反復の終点間差は「4000 step 後の再現性幅」であって**丸め誤差床そのものではない**。
> **`NOT CONVERGED` は保持**し、受理条件は置き換えない (§5.1 R5g を未完了へ戻した)。

**(1) ノイズ床の測定** (`run_0417_noise_r1`–`r5`): `run_0415` の step 40000 の場から
**同一入力・同一設定で 4000 step を 5 本**。node の `atomicAdd` は集積順が走るたびに変わるので同じ入力でも一致しない。
床は**全ペアの最大**、[[noise-floor-both-sides]] のとおり 3 本では桁で足りないので 5 本。
正規化は最終場の L2 ノルム。

| 量 | ノイズ床 (全ペア最大, rel L2) | 同じ 4000 step の変化 | **変化/床** |
| --- | --- | --- | --- |
| `ro` | 7.38e-05 | 9.30e-04 | **12.6** |
| `roUx` | 7.81e-05 | 6.47e-04 | 8.3 |
| `roUy` | 1.06e-04 | 3.70e-03 | **34.8** |
| `roUz` | 2.28e-04 | 2.85e-03 | 12.5 |
| `roe` | 5.57e-05 | 4.49e-04 | 8.1 |
| `roK` | 4.82e-05 | 2.96e-03 | **61.5** |
| `roOmega` | 5.28e-05 | 1.66e-04 | 3.1 |

**場は床の 3〜62 倍動いている = プラトーはノイズ床ではない**。候補 (2) は却下。
ただし余裕は大きくない: 残差があと 1 桁下がると変化が床に届く。**この離散化で到達しうる収束は残り 1 桁程度**。

**(2) それは収束の途中か、リミットサイクルか** — `run_0415` の窓ごと (4000 step) の相対 L2 変化:

| 窓 | `ro` | `roUy` | `roK` | `roOmega` |
| --- | --- | --- | --- | --- |
| 0→4000 | 3.07e-03 | 7.07e-03 | 2.41e-02 | 4.34e-03 |
| 8000→12000 | 1.22e-03 | 3.82e-03 | 1.33e-02 | 1.83e-03 |
| 20000→24000 | 1.10e-03 | 3.78e-03 | 6.93e-03 | 6.26e-04 |
| 28000→32000 | 6.24e-04 | 2.55e-03 | 5.24e-03 | 4.67e-04 |
| 36000→40000 | 7.17e-04 | 4.07e-03 | 3.44e-03 | 8.08e-04 |

`roK` だけは単調に減衰し続けている (2.4e-2 → 3.4e-3)。**`ro` / `roUy` / `roOmega` は減衰を止め、
7e-4 / 4e-3 / 8e-4 のあたりで振動している** = **リミットサイクル**。振幅はノイズ床の **8〜35 倍**。

**(3) それでも力は定常** (`collect`, 末尾 10 スナップショット):

| 量 | 平均 | 振幅 |
| --- | --- | --- |
| `C_T_with_shear` | **0.8595940** | **8.01e-07** |
| `C_T` | 0.8722673 | 6.37e-07 |
| `C_L` | 0.0783774 | 3.52e-06 |
| `C_M` | −2.049706 | 8.80e-05 |

**場のリミットサイクルは積分量をまったく動かさない** (相対振幅 1e-6 台)。
`check_quasisteady` も `ALL STEADY`。

**結論**: 本段の残差プラトーは **小スケールの真のリミットサイクル**であって、リミッタでも CFL でも
ノイズ床でもない。**「残差を下げる」方向は筋が悪い**。受理は
**(a) 派生量の `check_quasisteady`、(b) 力の平均±振幅、(c) リミットサイクル振幅とノイズ床の比**
で行い、残差プラトーはその振幅として**数値で記録する**のが正しい形。

**(4) 壁解像は推力を受理許容を超えて動かす**:

| | 粗 (`run_0403`) | 壁解像 (`run_0415`) | 差 |
| --- | --- | --- | --- |
| `C_T_with_shear` | 0.8563184 | **0.8595940** | **+0.0033** |

§8 の受理許容は $|\Delta C_T| \le 0.002$。**差 0.0033 はこれを超える**ので、
**粗メッシュは要求精度で格子収束していない**。生産は壁解像メッシュ側を使う。
局所ピークの差も同方向 (`machmax` 7.08 → 6.09、`pmax` 1.42e+05 → 1.27e+05)。

**(5) 壁解像が炙り出した欠陥: カウル後縁 × 側壁の三重線** — `GATES: FAIL fail_class=FLOOR_STUCK`:

| 節点 | 位置 | T [K] | P [Pa] | ρ | `wall_dist` |
| --- | --- | --- | --- | --- | --- |
| 517197 | (0.1200, −0.0105, 0.0963) | **37.6** | 8586 | 0.671 | 1.60e-05 |
| 517198 | (0.1200, −0.0105, 0.0972) | **49.8** | 8402 | 0.496 | 1.60e-05 |

x = 0.1200 m は **`L_cowl` ちょうど**、z 0.096–0.097 は側壁 (W/2 = 0.1) の直内、
`wall_dist` 1.60e-05 = **壁第 1 層**。つまり**カウル後縁と側壁が交わる三重線**。
全体では T < 80 K が **4 節点**、< 120 K が 8 節点だけで、p0.01 は 213 K (中央値 225 K) なので
**完全に局所**。粗メッシュでは第 1 層が 25 倍厚く、この特異点がならされて見えなかった。
**受理を止めているのはこれ 1 つ** → §5.1 R5h。

### 4.35 R5h: カウル後縁の低温スポット (**診断を 2 点訂正**, 2026-09-20)

> **訂正 1 (codex plan-2 M2)**: **`FLOOR_STUCK` は誤判定だった**。温度の実効下限は経路で違う —
> 単相 CPG (`thermalMethod: 0`, 凝縮なし) は `dependentVariables_d.cu` L290 の
> `T = max(intE/(cp/gamma), tMin)` で床は **config の `tMin` (既定 1e-4 K)**、
> **50 K の `DEPVAR_TMIN` は TP/多成分/凝縮の温度反転経路だけ** (同 L148/151/160)。
> `sern_gates.py::_solver_floors` が経路を見ずに一律 50 K を下限にしていたのが原因。
> **修正済み** (`thermalMethod`/`condensation` を見るようにした)。修正後 `run_0415` は **`GATES: PASS`**。
> **37.6 K は床に張り付いてはいない** — 非物理な低温スポットであることは変わらないが、分類は誤りだった。
>
> **訂正 2 (codex plan-2 M4)**: **後縁は「厚さ 0 のスリット」ではない**。`dup1` は
> [`mesh_sern3d.py`](../../design/forge_design/meshing/mesh_sern3d.py) の `for i in range(i_te)` で
> **`i < i_te` にしか作られない**ので、後縁 station では `cowl_in` と `cowl_out` が**同一ノード ID**
> (実測: 両タグの共有ノード **25 個** = z-station 数)。これは**正しく閉じた後縁**であって、
> 上流の「別 ID で座標が一致するスリット」とは別物。したがって「`cowl_thickness` で避けたはずの形が残っていた」
> という断定は誤り。残る事実は「**1 本の壁ノードが片側 2772 K・反対側 ~1000 K を受け持ち、
> 第 1 層 1.6e-05 m でその状態が潰れる**」こと。

`run_0415` で `tMin` に張り付いた 2 節点 (517197 / 517198) を追った。

**位置**: x = 0.12000 m = **`L_cowl` ちょうど (カウル最終 station)**、z 0.0963 / 0.0972 (側壁 W/2 = 0.1 の直内)、
`wall_dist` 1.60e-05 = **壁第 1 層**。どちらも**境界に属さない内部節点**で、座標一致の双子も無く、
双対 CV 体積は壁第 1 層の中央値の 0.31–0.36 倍 = **退化していない**。

**壁法線に沿った状態** (x = `L_cowl`, z = 0.09631、y 昇順):

| y | T [K] | ρ | P [Pa] | \|U\| | `wall_dist` |
| --- | --- | --- | --- | --- | --- |
| −0.010657 | 160.9 | 0.182 | 9938 | 539 | 1.58e-04 |
| −0.010557 | 86.8 | 0.351 | 10364 | 396 | 5.81e-05 |
| −0.010534 | 57.9 | 0.581 | 11432 | 268 | 3.52e-05 |
| **−0.010515** | **37.6** | **0.671** | 8586 | 72 | **1.60e-05** |
| **−0.010499** | 94.7 | 0.261 | 8417 | **0** | **0 (壁)** |
| **−0.010483** | **370.7** | 0.056 | 7020 | 160 | **1.60e-05** |
| −0.010438 | 1101.9 | 0.020 | 7344 | 662 | 6.03e-05 |
| −0.010324 | 1854.0 | 0.014 | 8530 | 1344 | 1.74e-04 |

壁の上 (ノズル側) は 2772 K の排気、下 (外部流側) は回復温度 ~1000 K。
**その 2 つが厚さ 0 の板を挟んで 3.2e-05 m で向き合っている**。

**1 station 手前まではまったく正常** (z 固定・x 昇順、壁第 1 層の上下ペア):

| x | 下側 T / ρ / P | 上側 T / ρ / P |
| --- | --- | --- |
| 0.11871 | 995.7 / 0.0176 / 5961 | 2769.5 / 0.0281 / 26490 |
| 0.11936 | 801.9 / 0.0202 / 5512 | 2864.8 / 0.0285 / 27828 |
| **0.12000 (`L_cowl`)** | **37.6 / 0.671 / 8586** | **370.7 / 0.056 / 7020** |

**最終 station だけが崩れている**。

**原因**: 板厚の法則が `np.interp(xs, [−L_up, 0.8·L_cowl, L_cowl], [t, t, 0])` で、
**x = `L_cowl` で厳密に 0** になる。実測でも x = `L_cowl` の 25 z-station すべてで
`cowl_in` と `cowl_out` の y 差が **0.000e+00** (板全体の中央値は 5.000e-04)。
= **後縁 1 station だけが厚さ 0 のスリットとして残っている**。これは
`cowl_thickness` を導入して避けたはずの形そのもの (2D `run_0035`: 厚さ 0 は m6_on/m10_on とも発散)。

**壁解像が作った新しい欠陥であって、前からあったものが鋭くなったのではない**:

| run | T 最小 | 位置 | `wall_dist` | T < 120 K |
| --- | --- | --- | --- | --- |
| `run_0403` (粗) | 157.95 | (1.174, 0.469, 0.250) | 2.51e-01 | **0 節点** |
| `run_0416` (粗・継続) | 157.10 | 同上 | 2.51e-01 | **0 節点** |
| `run_0415` (壁解像) | **37.60** | **(0.120, −0.0105, 0.0963)** | **1.60e-05** | 8 節点 |

粗メッシュの T 最小は**プルームの過膨張 (遠方, `wall_dist` 0.25)** であってカウル後縁ではない。
第 1 層が 4.0e-04 のときは厚さ 0 スリットの両側の状態がならされていたが、
**1.6e-05 (25 倍薄く) にすると 2772 K と 1000 K が 3.2e-05 m で向き合い、状態が潰れる**。

**対処案** (R5h):

| 案 | 中身 | 費用 | 副作用 |
| --- | --- | --- | --- |
| **(a) 後縁に有限の厚みを残す** | 板厚のテーパを 0 でなく有限値で止め、`cowl_base` タグ (25 z-cell × 第 1 層) で閉じる。**機体ベース (R4e 案 (d), `t_base` + `vehicle_base`) と同じ型**で、実績がある | 中 (新タグ + 力の帳簿 + 試験) | 形状が変わる (物理的にはむしろ正しい: 実物のカウル後縁は有限厚) |
| (b) 後縁近傍だけ壁層を粗いままにする | `first_wall_frac` を x でブレンド (§4.32 で必要性を確認した仕組み) | 小 | **後縁を意図的に解像しない** = 逃げ |
| (c) 受理して床でクリップ | 何もしない | 0 | `FLOOR_STUCK` で落ちるので受理できない |

**(a) を推す** — (b) は後縁という一番効く場所の解像を落とすので、`C_T` の格子収束を主張できなくなる。
ただし新タグを足すので **codex の plan レビューを先に通す**。

### 4.36 codex plan-2 が見つけた 2 つの前提崩れ (2026-09-20)

**(1) R5 系の run は生産条件ではない — CPG + 断熱壁** (M1)。
`problem_3d_sst_cycle_m6on{,_wallres,_wallres_y}.yaml` はいずれも `gas` 節と `spec.wall_thermal` を
省略しており、`probdef.py` の既定 (ガスは CPG、壁は断熱) が効く。実際の run の config で確認:
**`thermalMethod: 0`**、`ramp`/`cowl_in`/… はすべて **`kind: wall`・`floats:` 空 = 断熱**。
したがって §4.28–§4.35 の数値は**すべて CPG・断熱の診断値**であり、
R4f の「等温壁 (Ts 1000 K) に統一済み」は**この問題 YAML には効いていない**。
外部動圧も生成値 **60698.9 Pa** で YAML コメントの 71850 Pa と違う。
**生産の受理はこの条件で取り直す必要がある** → §5.1 R5i。

**(2) カウル側端のテーパ幅が格子に追随する** (M5)。
`mesh_sern3d.py` の `kt = max(k_sw - 2, 0)` は板厚を落とす範囲を **「最後の 2 z-セル」**で決めている。
物理幅は格子で変わる:

| `nz_in` / `first_z_frac` | 側端テーパ幅 |
| --- | --- |
| 25 / 4.0e-03 | 0.868 mm |
| 57 / 4.4e-05 | **0.00949 mm** |

**91 倍**。側壁を解像する変更が**カウル形状そのものを変える**ので、粗細の差を格子誤差に帰属できない。
[[geometry-must-not-follow-mesh-spacing]] (2026-09-19 に `vehicle_clearance` で同じ直しをした) の
**同型の違反が残っていた**。→ §5.1 R5j。

### 4.37 R5i: 生産条件を 3D に入れ直した — **断熱壁とプラトーは 2D で既知だった** (2026-09-20)

§4.36-1 を受けて、2D 生産 YAML (`problem_cflA_base.yaml`) から生産条件を確定させた:

```yaml
gas:  {model: frozen_tp, thermo_href_temp: 298.15, gamma: 1.2470, cp: 1718.5}
spec: {wall_thermal: {mode: isothermal, Tw: 1000.0}}
      inflow.M_in: 1.6308          # 3D は 1.6745 だった (凍結 γ に対応する値)
      m6_on.gas: γ 1.2470 / cp 1718.5 / CEA 凍結組成 11 種
```

**2D 生産 YAML の `wall_thermal` のコメントが、私が 3D で追っていたものをそのまま書いている**:

> 断熱だと M∞10 板下面が回復温度 **4600 K (非物理)** で**残差プラトー 5e-5**;
> **等温 1000 K で残差 3.5 桁低下** (run_0103)

§4.31–§4.34 でリミッタ・CFL・壁解像・ノイズ床と順に潰してきたが、**その全期間、条件が生産と違っていた**。
R5h のカウル後縁 37.6 K も「**断熱壁の回復温度がそこだけ崩れている**」で筋が通る。
**断熱 CPG と生産 (等温 TP) の差を同一メッシュ・同一 CFL で測る**のが先決 → `run_0418_3d_prod_wallres`。

**新しい問題 YAML**: `problem_3d_prod_m6on_wallres.yaml`。壁解像の knob (`first_wall_frac` 1.6e-04 等・
`ar_max` 5000) は据え置き、gas / wall / 作動点だけ生産へ。`cfl_main` は **0.5 → 0.25**
(壁解像メッシュで 0.5 は本段 step 97 で発散する、§4.33)。

**生成された config を投入前に検証**した (`_prodcheck`, 破棄済):

| 項目 | 値 |
| --- | --- |
| `thermalMethod` | **2** (TP) |
| `ramp` / `cowl_in` / … | **`wall_isothermal`, `Ts: 1000.0`** |
| メッシュ品質 | **`VERDICT: PASS`** (1,838,592 セル、AR max **4212.9** / p99 992 / skew max 0.449、>5000 が 0 セル) |

**形状が変わる点に注意**: γ が 1.1828 → 1.2470 になるので MOC 設計の輪郭が動く
(セル数 1,837,440 → 1,838,592)。**`run_0403`–`run_0417` とは別形状**であり、直接比較はできない。

**YAML の罠 2 件** ([[config-rewrite-yaml-hazards]] に同型): 化学種キー `NO` は**真偽値に化ける**のでクォートが要る、
`5e-05` は**文字列に化ける**ので `0.00005` と書く。どちらも `load_problem` / `sum()` で落ちて気づいた。

### 4.37.1 生産条件の結果: **プラトーは消えない**。ただし後縁の低温スポットは 3 倍改善し、y⁺ は 3 倍悪化 (2026-09-20)

`run_0418_3d_prod_wallres` (TP + 等温壁 1000 K、壁解像メッシュ、`cfl_main` 0.25、本段 20000 step)。
比較対象は**同じ位置の**断熱 CPG run `run_0414_wallres_y_cfl025` (同じ壁解像 knob・同じ CFL・同じ 20000 step)。

| 残差列 | 生産 (TP + 等温) `run_0418` | 断熱 CPG `run_0414` |
| --- | --- | --- |
| `rms_ro` | 1.0 dec | 1.1 dec |
| `rms_roUx` | 1.0 dec | 0.9 dec |
| `rms_roUy` | **1.5 dec** | 0.8 dec |
| `rms_roUz` | 0.9 dec | 0.9 dec |
| `rms_roe` | 0.9 dec | 0.7 dec |
| `rms_roK` | **1.1 dec (falling)** | 0.4 dec |
| `rms_roOmega` | 0.6 dec | 0.5 dec |
| `check_convergence` | **NOT CONVERGED (stalled/plateau)** | **NOT CONVERGED (stalled/plateau)** |
| `check_quasisteady` | **ALL STEADY** | **ALL STEADY** |
| `GATES` | **PASS** | PASS (ゲート修正後) |

**等温壁でも 3D のプラトーは消えなかった**。2D のコメントが言う「3.5 桁低下」は起きていない
(改善は `roUy` と `roK` で 0.7 桁程度)。**「断熱条件だったからプラトーしていた」という見立ては外れ**。
ただし条件は生産に戻ったので、以後の議論はこの run を基準にする。

**後縁の低温スポットは大きく改善した** (R5h):

| | 断熱 CPG `run_0415` | 生産 `run_0418` |
| --- | --- | --- |
| T 最小 | **37.60 K** | **109.63 K** |
| 位置 | (0.1200, −0.01051, 0.09631) | **(0.1200, −0.01056, 0.09719)** — **同じカウル後縁 × 側壁近傍** |
| `wall_dist` | 1.60e-05 (第 1 層) | 5.80e-05 (第 2 層) |
| T < 120 K | 8 節点 | **1 節点** |
| T の p0.01 | 213.1 K | 195.6 K |

**3 倍近く上がり、影響節点は 8 → 1 に減ったが消えてはいない**。場所も同じ。
= ~~断熱壁の回復温度は**寄与していたが主因ではない**~~ **訂正 (2026-09-27, codex diagnose Major)**: 比較は断熱 CPG → 等温 TP で熱物性も同時に変わり、最小点の位置・壁距離も違うので、回復温度への帰属はできない。**「条件一式の変更で最低温度が上昇した」**までが言えること。R5h は残る。

**冷却壁は y⁺ を 3 倍悪くする** — これが R5d に直接効く:

| 壁 | 断熱 CPG `run_0414` | 生産 `run_0418` | 倍率 |
| --- | --- | --- | --- |
| `cowl_in` | 1.608 | **4.919** | ×3.1 |
| `ramp` | 0.700 | **2.167** | ×3.1 |
| `cowl_out` | 0.846 | 0.718 | ×0.85 |
| `vehicle_top` | 2.354 | 2.254 | ×0.96 |
| `sidewall_in` | 50.9 | **213.2** | ×4.2 |
| `vehicle_side` | 16.1 | 21.8 | ×1.4 |

[[isothermal-wall-chain-plate-validation]] の「冷却で y1+ ×5–6」と同じ向き
(壁が冷えると近傍の密度が上がり粘性が下がる)。
**断熱条件で設計した壁解像メッシュは、生産条件では解像不足**。
`cowl_in` を y⁺ 0.7 にするには `first_wall_frac` をさらに **約 7 倍**細かく (1.6e-04 → 2.3e-05) する必要があり、
AR は現在の 4213 から **~29000** に跳ねる (§4.32 と同じ壁)。
**→ 遠方の長辺を割って AR を下げる費用を測る** (codex plan-2 M6 が指摘した未評価の道) → R5m。

### 4.38 壁解像ツールを輸送モデル非依存にした (2026-09-20)

`check_wall_resolution.py` の `mu_of` は Sutherland を再現していて **CPG 専用**、
多成分 (`viscMethod` 2 以降) では `None` を返して**判定不能**になっていた。
生産が TP になると**壁解像が測れない**ので、**体積出力 `res_<step>.h5` の `vis_lam` を優先して使う**ようにした
(DOF ごとの分子粘性で、輸送モデルに依らず正しい)。無ければ従来の `mu_of` にフォールバックする。

回帰 (`run_0415`, CPG): `cowl_in` 1.608 / `ramp` 0.700 / `sidewall_in` 50.833 (旧 50.874) で**一致**。

### 4.39 **ユーザ決定 (2026-09-20): 残差プラトーは受理の阻害要因としない**

> 「残渣が落ち切らないからと言って、正直問題とは思っていないんです。もとから変動が小さいのかもしれませんから。」

**この判断を採る**。裏づけとして揃っているもの:

| 根拠 | 値 |
| --- | --- |
| 力の定常性 (`collect`, 末尾 10 スナップショット) | `C_T_with_shear` **0.8595940 ± 8.01e-07**、`C_T` ± 6.37e-07、`C_L` ± 3.52e-06、`C_M` ± 8.80e-05 |
| 派生量 (`check_quasisteady`) | **ALL STEADY** (shock / machmax / pmax とも drift 0.0 %) |
| 場の動き vs 再現性幅 (§4.34) | 3〜62 倍 (= 桁違いに大きくはない) |

**ただし判定ツールの出力は消さない**。`check_convergence.py` の **`NOT CONVERGED (stalled/plateau)` は
そのまま併記**し、「未収束を見落とした」のか「小振幅と判断した」のかが後から区別できるようにする
(codex plan-2 M3 は「**撤回した論拠で受理条件を置き換えるな**」という指摘であって、
ユーザが振幅を見て許容すること自体を否定していない)。

**したがって受理の本体は次の 2 つに移る**:

1. **生産条件での力係数の格子収束** (§8 の $|\Delta C_T| \le 0.002$)。断熱 CPG では粗↔壁解像で **+0.0033** と
   許容を超えていた (§4.34-4)。生産条件 (TP + 等温) で粗・中・細を揃えて示す → **R5n**。
2. **壁解像の扱い**。生産の冷却壁で y⁺ は ×3.1 悪化し (`cowl_in` 4.92)、y⁺≈1 は AR ~29000 を呼ぶ。
   AGENTS.md の規定は「y₁⁺≤1 を目標とし、**緩和するときは熱流束・摩擦の格子感度で裏付ける**」なので、
   **1 を追うより 1. の格子列で裏づけて緩和する**。R5m の AR 費用測定はその上限を知るために続ける。

**§5.1 の R5g は「追わない」に変更**し、残差プラトーの原因究明は打ち切る。

### 4.40 R5n: 生産条件の格子収束 — **未収束なのは摩擦だけ** (2026-09-20)

生産条件 (frozen_tp + 等温壁 1000 K)、**形状・BC は同一**、壁法線解像だけを 4 倍ずつ振った 3 点
(`cfl_main` 0.25、本段 20000 step、3 本とも `GATES: PASS`):

| | g1 `run_0419` (2.56e-03) | g2 `run_0420` (6.4e-04) | g3 `run_0418` (1.6e-04) | g1→g2 | **g2→g3** | §8 許容 |
| --- | --- | --- | --- | --- | --- | --- |
| **`C_T` (圧力のみ)** | 0.8921423 | 0.8901614 | 0.8892752 | −0.00198 | **−0.00089** | 0.002 ✅ |
| **`C_L`** | 0.0611365 | 0.0613020 | 0.0624612 | +0.00017 | **+0.00116** | 0.002 ✅ |
| `C_T_with_shear` | 0.8595643 | 0.8673550 | 0.8705395 | +0.00779 | **+0.00318** | 0.002 ❌ |
| `C_M` | −1.622503 | −1.619331 | −1.641517 | +0.00317 | **−0.02219** | 0.02 ❌ (僅か) |
| `C_T_friction` | −0.0325780 | −0.0228064 | −0.0187356 | +0.00977 | **+0.00407** | — |

**未収束の正体は摩擦である**。`C_T` (圧力のみ) は **−0.00198 → −0.00089** と半減して許容内に入り、
`C_L` も許容内。`C_T_with_shear` が許容を外れるのは `C_T_friction` が
**−0.03258 → −0.02281 → −0.01874** と動き続けているからで、差の比 2.40
(= 収束次数 $p=\ln 2.40/\ln 4 = 0.63$)。Richardson 外挿すると
`C_T_friction` → **−0.0158**、g3 に残る誤差は **~0.0029** で、そのまま `C_T_with_shear` の残差になる。

**これは AGENTS.md が壁解像の緩和に要求する「熱流束・摩擦の格子感度」そのもの**であり、
**緩和できないことを示している**: 摩擦はまだ収束していない。

**あと 1 段 (`first_wall_frac` 4e-05) で残差 ~0.0013 まで落ちる見込みだが、AR が ~16900 になる**。
§4.32 / R5m と同じ壁。**R5m の測定によれば x ブレンド (壁が無い下流に細層を置かない) を入れれば
AR は 1/2.8 に下がる**ので、4e-05 でも ~6000 と射程に入る。

**したがって次の分岐点は「x ブレンドを実装するか」**:

| 案 | 中身 | 得られるもの | 費用 |
| --- | --- | --- | --- |
| **(A) x ブレンドを入れて 1 段細かくする** | `first_wall_frac` を x でブレンド (壁の無い下流は粗いまま) + `first_wall_frac` 4e-05 | **摩擦が収束し `C_T_with_shear` と `C_M` が許容内に入る見込み**。節点も**減る** (無駄な細層が消える) | 生成器 1 箇所 + 試験 |
| (B) 摩擦の非収束を明記して受理 | 現状 g3 を採用し、`C_T_with_shear` の格子誤差 **±0.0029** を成績に併記 | すぐ受理できる | 許容 0.002 を満たさないと明記が要る |

**(A) を推す** — x ブレンドは §4.32 (AR)・§4.35 (後縁) ・R5m のいずれでも必要性が出ており、
**節点を減らしながら**精度を上げる唯一の手。ただし生成器変更なので実装前にユーザ判断を仰ぐ。

### 4.41 R5o 設計: 壁第 1 層の **x ブレンド** (実装前に固定, 2026-09-20 ユーザ承認)

**問題**: `first_wall_frac` は全 station に一律に効く。`s_bot` の細端は**中間線** (カウル = 壁なのは
x ≤ `L_cowl` まで)、`s_top` の細端は**中間線と上線** (ランプ = 壁なのは x ≤ `L_ramp` まで) なので、
**それより下流は壁が無いのに壁用の細層を敷いている**。R5m の実測で AR > 5000 セルの **57.8 %** がそこにあり、
外側 z を細分しても **max AR は 1 も下がらなかった** (29712 のまま)。

**法則** (物理長で決める。格子間隔から決めない = [[geometry-must-not-follow-mesh-spacing]]):

$$ f(x) = \exp\big[(1-w)\ln f_{\rm wall} + w \ln f_{\rm far}\big], \qquad
   w = s\!\left(\mathrm{clip}\frac{x - x_{\rm end}}{L_{\rm blend}}\right),\quad s(t)=t^2(3-2t) $$

- $f_{\rm wall}$ = `first_wall_frac` (現行)、$f_{\rm far}$ = **新規 `first_wall_frac_far`**
- $x_{\rm end}$ = **下バンドは `L_cowl`、上バンドは `L_ramp`** (それぞれの壁の終端)
- $L_{\rm blend}$ = **新規 `wall_frac_blend_len`** (/H, 既定 0.5)。**物理長**であって「最後の N セル」ではない
- 間隔は 25 倍も跨ぐので**対数補間** (線形だと遷移の後半でほぼ $f_{\rm far}$ に張り付く)

**既定は挙動不変**: `first_wall_frac_far` 未設定 (= 0) なら $f_{\rm far} = f_{\rm wall}$ として
**ブレンドを一切かけない**。既存の 3D メッシュはビット一致で再現される。

**上バンドは片側しか選べない**: `_tanh_two_sided(n, first)` は対称で両端が同じ `first` になる。
station ごとに 1 つしか選べないので、**どちらかの端が壁なら細い方** (= x ≤ `L_ramp` で $f_{\rm wall}$) を採る。
`L_cowl` < x ≤ `L_ramp` の区間は下端 (プルーム界面) にも不要な細層が付くが、そこの接線刻みは ~2e-3 m なので
**AR は ~400 で無害**。

**適用範囲は主ブロックの `s_bot`/`s_top` のみ**。`ext_top` の `first_top_frac` は機体上面が全長にわたって壁なので触らない。
2D (`mesh_sern.py`) は同じ無駄があるが AR が問題になっていないので**今回は変えない** (差分を最小にする)。

**検証** (`run_sern_mesh3d_tests.py` に追加):

1. **既定で挙動不変**: `first_wall_frac_far` 未設定なら座標が**ビット一致**すること。
2. **ブレンドが実際に効く**: `first_wall_frac_far` を与えたとき、**生成後の実座標**で
   x ≤ `x_end` の第 1 層厚が `first_wall_frac`、x ≥ `x_end` + `L_blend` の第 1 層厚が
   `first_wall_frac_far` に一致すること (§4.30 と同じく**入力値でなく生成結果**を見る)。
3. **形状が動かないこと**: `first_wall_frac_far` を変えても壁の輪郭 (`y_veh`, ランプ・カウル座標) が不変。
4. **境界の閉性・符号付き Jacobian・float32 衝突**は既存試験がそのまま掛かる。

**受け入れ**: `first_wall_frac` 4e-05 + `first_wall_frac_far` 4e-03 で
**AR ≤ 5000 かつ節点数が現行 (1.92M) 以下**になること。そのうえで R5n の格子列に 4 点目として足し、
`C_T_with_shear` と `C_M` が §8 の許容に入るかを見る。

### 4.42 R5o の結果: 摩擦は収束した。残るのは **C_M の許容と側壁** (2026-09-21)

**メッシュ**: `run_0421` / `run_0422` (3,122,688 セル)。x ブレンド + `nz_out` 32 で
**`VERDICT: SOFT-PASS`** (AR max 12866 / **>5000 が 0.02 %**、**skew max 0.701 / >0.90 が 0**)。

**壁解像は y 法線壁で目標到達**:

| 壁 | g3 (1.6e-04) | **g4 (4.0e-05)** |
| --- | --- | --- |
| `cowl_in` | 4.919 | **1.234** |
| `ramp` | 2.167 | **0.549** |
| `cowl_out` | 0.718 | **0.185** |
| `sidewall_in` | 213.2 | 188.3 (**未解像のまま**) |

**力係数** (`run_0422` = g4 を +40000 step 継続して定常化。`check_quasisteady` **ALL STEADY**、
力の振幅は `C_T_with_shear` ±1.95e-06 / `C_M` ±3.83e-04):

| | g3 `run_0418` | **g4 `run_0422`** | **Δ (g3→g4)** | §8 許容 | 判定 |
| --- | --- | --- | --- | --- | --- |
| `C_T` (圧力) | 0.8892752 | 0.8891410 | **−0.00013** | 0.002 | ✅ |
| **`C_T_with_shear`** | 0.8705395 | 0.8702979 | **−0.00024** | 0.002 | ✅ |
| `C_L` | 0.0624612 | 0.0639124 | **+0.00145** | 0.002 | ✅ |
| `C_M` | −1.641517 | −1.683246 | **−0.04173** | 0.02 | ❌ |
| `C_T_friction` | −0.0187356 | −0.0188430 | −0.00011 | — | **収束** |

**摩擦は収束した** — g2→g3 の +0.00407 が g3→g4 で **−0.00011** になり、
`C_T_with_shear` も +0.00318 → **−0.00024** と許容内に入った。**R5o の目的は達成**。

**C_M だけが外れる理由は許容値の不整合である**。モーメント基準は $x_{\rm ref} = -20H$ で、
力の作用点までの腕は 20〜29 H ある。したがって $\Delta C_L$ はそのまま
$\Delta C_M \approx \Delta C_L \times (d/L_{\rm ref})$ に増幅される:

$$ 0.00145 \times 29 = 0.042 \;\simeq\; |\Delta C_M| = 0.0417 $$

**観測値はこの増幅でほぼ全部説明できる**。つまり §8 の対
($|\Delta C_L| \le 0.002$, $|\Delta C_M| \le 0.02$) は**内部矛盾している**:
腕が 20〜29 H なら $\Delta C_L = 0.002$ だけで $\Delta C_M$ は 0.04〜0.06 になる。
**どちらかを直す必要がある** (C_L を ~0.0007 まで締めるか、C_M を ~0.05 に緩めるか) → **ユーザ判断**。

**残る阻害要因は側壁**: `GATES: FAIL fail_class=FLOOR_STUCK` — **ρ が床 (1e-4) に張り付いた壁ノードが 8 個**
(20000 step 時点では 13 個)。位置は **z = 0.10000 = 側壁の壁ノードそのもの** (`wall_dist` 0)、
x 0.025–0.035、T は等温壁の 1000 K 固定、P ≈ 34 Pa (圧力床 20 Pa の直上)。
**側壁は y⁺ 188 で未解像**のままなので整合的。

**交絡の明示**: g4 は g3 に対し **3 つ同時に変えている** (`first_wall_frac` 1/4・`nz_out` 15→32・x ブレンド ON)。
純粋な壁法線 1 段細分ではないので、**Δ(g3→g4) を「格子収束の 1 段」と読むのは厳密には誤り**。
摩擦の収束は壁量なので壁間隔で説明できるが、`C_M` の差に `nz_out`/ブレンドがどれだけ効いたかは分けていない。

### 4.43 R5q の正体: **厚さ 0 が 2 つ交わる線** (2026-09-21)

`run_0422` で ρ が床に張り付いた 12 ノードを追った。**所属は `sidewall_out` (physID 12) のみ**、
位置は x 0.0272–0.0481 m・y −0.00415–−0.00232 m・**z = 0.10000 (= W/2)**。
**同じ座標に双子が乗っており、片方だけが枯れている**:

| 双子 | ρ | P [Pa] | T | \|U\| |
| --- | --- | --- | --- | --- |
| 内側 (ノズル) | 2.498e-01 | **85018.6** | 1000.0 | 0 |
| **外側** | **1.000e-04 (床)** | **28.7** | 1000.0 | 0 |

外部静圧 2851 Pa の **1/100**。**格子を細かくするほど悪化する**:

| run | `first_wall_frac` | ρ min | P min [Pa] | 床ノード |
| --- | --- | --- | --- | --- |
| g1 `run_0419` | 2.56e-03 | 1.95e-03 | 811.5 | 0 |
| g2 `run_0420` | 6.4e-04 | 1.61e-03 | 719.5 | 0 |
| g3 `run_0418` | 1.6e-04 | 9.51e-04 | 273.1 | 0 |
| **g4 `run_0422`** | **4.0e-05** | **1.00e-04 (床)** | **28.7** | **12** |

**格子収束していない = 床の副作用ではなく真の欠陥**。§4.35 のカウル後縁と同じ指紋。

**幾何の実測 — 厚さ 0 が 2 つある**:

| | 実測 |
| --- | --- |
| 側壁 (`sidewall_in` vs `out`) の z 差 | **min = 中央値 = max = 0.000e+00** = **厚さ 0 のスリット** (双子 5893 対) |
| カウル板厚 | z=0.0000 で 5.000e-04、z=0.0980 で 5.000e-04、**z=0.1000 (側壁) で 0.000e+00** |

枯れるのは**この 2 つの交線**。`sz` の法則が「最後の 2 z-セルで板厚を 0 に絞る」ので、カウル板は
**側壁ちょうどで刃先に潰れる**。そこに厚さ 0 の側壁スリットが重なり、4 面が 1 点で出会う。

**「カウルを潰さない」だけでは直らない**。コードに理由が書いてある:
「側壁 (z = W/2) は厚さ 0 のスリットなので、そこまで板厚を効かせるとスリットの内外で中間線の高さがずれ、
側壁後縁が破綻する (`run_0086`: 暖機段 step 4 で `ro` NaN)」。
**側壁が厚さ 0 である限り、内側は上下 2 本・外側は 1 本の中間線を要求して矛盾する**。

### 4.44 R5q 設計: 側壁に厚みを与え、カウルの潰しをやめる (実装前に固定)

**2 つは 1 組**。側壁に厚み $t_{sw}$ を与えると内面 $W/2 - t_{sw}/2$ と外面 $W/2 + t_{sw}/2$ が
別座標になり、内面が厚いカウルの上下 2 本を、外面が単一の中間線を受けられる。

**(1) 側壁厚さは物理入力**、x でテーパして `L_sw` で 0 に閉じる ([[geometry-must-not-follow-mesh-spacing]]):

$$ t_{sw}(x) = t_{sw}\cdot\mathrm{interp}\big(x;\; [-L_{up},\, 0.8L_{sw},\, L_{sw}],\; [1,1,0]\big) $$

新規パラメータ `sidewall_thickness` (/H, 既定 **0.0 = 従来挙動**)。カウルの `cowl_thickness` と同じ法則の形。

**(2) z 分布を station 依存にする**。現在 `zs` は 1 次元で全 (i,j) 共有。両端だけを動かす**相似縮小**にする
(内部の分布比 `s_in`/`s_out` は不変なので $t_{sw}=0$ で**ビット一致**):

$$ z^{\rm in}_i = \big(\tfrac{W}{2} - \tfrac{t_{sw}(x_i)}{2}\big)\, s_{\rm in},\qquad
   z^{\rm out}_i = \big(\tfrac{W}{2} + \tfrac{t_{sw}(x_i)}{2}\big) + \big(Z_{\rm ext} - \tfrac{t_{sw}(x_i)}{2}\big)\,(1 - s_{\rm out}[::-1]) $$

触る箇所は **4 つだけ** (調査済):
`coords[b:b+nz, 2] = zs` (L262) → `zsx[i]`、`dup1` の z (L269) → `zsx[i][k]`、
`dup2` の z (L273) → $W/2 + t_{sw}(x_i)/2$、面分類の `zm` (L321) → `zsx[i]`
(分類は $z \lessgtr W/2$ の比較なので $t_{sw}/2$ のずれでは変わらない)。

**(3) カウルの側端テーパ `sz` を廃す**。`sz` は「$k \le k_{sw}$ で 1、$k > k_{sw}$ で 0」の**ステップ**になる
(側壁の外に板は無いので中間線は単一のまま)。**これが R5j (側端テーパ幅が格子に追随、`nz_in` 25→57 で 91 倍) の是正も兼ねる**。

**(4) 共有ノード ID は変えない**。`dup2` の生成条件 ($i < i_{sw}$, $j > j_m$) も `k_sw` の索引も不変。
トポロジは同一で、**動くのは座標だけ**。

**検証** (`run_sern_mesh3d_tests.py` に追加):

1. **既定 (`sidewall_thickness: 0`) で座標がビット一致** — 既存メッシュを一切動かさない。
2. **側壁厚が設計値どおり**: 生成後の実座標で `sidewall_in` と `sidewall_out` の z 差が
   $x \le 0.8L_{sw}$ で $t_{sw}$、$x \ge L_{sw}$ で 0。
3. **カウル板厚が側壁まで保たれる**: `cowl_in`−`cowl_out` の y 差が z=W/2 近傍でも `cowl_thickness`。
4. **形状が格子から独立**: `nz_in` を変えても側壁厚・カウル板厚・輪郭が不変 (R5j の是正確認)。
5. 境界の閉性・符号付き Jacobian・float32 新規衝突 0 は既存試験がそのまま掛かる。

**受け入れ**: メッシュ品質が g4 と同等以上 (AR >5000 が ≤0.1 %、skew >0.90 が 0) で、
**ρ 床ノードが 0**、かつ格子列を取り直して 4 係数が §8 の許容内。

**形状が変わる**: ノズル幅が内側へ $t_{sw}/2$ 狭まるので **g1–g4 は全部取り直し**。

### 4.45 R5r: カウル板の**自由な側端**の双対 CV が開いていた (旧メッシャ。2026-09-22 修正)

新メッシャの受入用に作った **CV ごとの閉性検査** (`solver_density_cuda/tools/check_dual_closure.py`。各 CV の $|\Sigma S|/\Sigma|S|$) を
旧メッシャの格子に掛けて見つけた。変換器のログの閉性 (全域の最大面積で正規化した 1 つの数) には埋もれていた。

| 格子 | 開いた CV ($>10^{-2}$) | 最大 | 場所 |
| --- | --- | --- | --- |
| `case/46.sern_design/run_0402_3d_wf1/sern.h5` (116 万節点) | **68** | **0.304** | 全部 $z = W/2$・$x \in [0.081, 0.119]$ m = $[L_{sw}, L_{cowl}]$。cowl_in 34 + cowl_out 34 = $2\,(i_{te} - i_{sw})$ |
| 小さい再現格子・修正前 (3.7 万節点) | **18** = $2 \times 9$ | 0.321 | 同上 |
| 同・**修正後** | **0** | 6.2e-06 (`VERDICT: PASS`) | — |

**原因**: `mesh_sern3d.py` の `dup1` (カウル上面側のコピー) が $k = k_{sw}$ を含んでいた。側壁がある $x < L_{sw}$ ではそこは T 字の交線で、上下を分ける必要がある。
しかし**側壁が終わった $L_{sw} \le x < L_{cowl}$ では、そこは厚さ 0 の板の自由な側端**で、上下の流体はつながっている。端の節点まで二重化すると、
双子の双対 CV は互いの間の面を持たず、30 % 開く (primal の面の閉性試験では見えない。上下とも壁面として数えられるため)。板の後縁 ($i_{te}$) は最初から共有だった。
**修正**: `share_cowl_free_edge` (既定 True)。$i \ge i_{sw}$ の $k_{sw}$ 列は二重化しない。変換器の全域閉性も 1.2e-06 → 3.8e-09 に下がった。
回帰試験 4 件を `design/tests/run_sern_mesh3d_tests.py` に追加 (`ALL PASS`)。

**R5q との関係 — 別の場所**。R5q の床ノード 12 個は $x$ 0.027–0.048 m = **$x < L_{sw}$ の T 字部**の `sidewall_out` で、今回の自由端 ($x \ge L_{sw}$) ではない。
`run_0402` の開いた CV のうち側壁タグに属するのは 2 個だけ。**この修正で R5q が直るとは言えない**。ただし自由端の 68 CV は質量・運動量の偽の湧き出しなので、
力の係数への影響は測り直す必要がある (R5s)。g4 (`run_0422`) の格子で床ノードと開いた CV が重なるかも、同じツールで確かめる (格子は AWS 上)。

### 4.46 R5o-chi + R5s: 現行コードで g3/g4 を取り直した (2026-09-27)

**条件**: 現行コード (chi 既定 auto = 1・スカラー勾配 node 既定 lsq・R5r のメッシャ修正込み、バイナリ `2fa3826c`、
`FORGE_CUDA_BLOCKSIZE=128`)。**出口は `outflow`** (runner の既定。[`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) §5.1 B1f で 2026-09-23 に変更)。**旧列 `run_0418`–`run_0422` は `statPress`** なので、旧列との比較には出口 BC の変更も混ざる (B1f は「直接比較しない」としている。2026-09-27 に記録後に気づいて訂正)。起動エコーで両 run とも `slauWallNormalChi` 1 (auto)・`scalarGradient` lsq (default) を確認。
g3 は旧 `run_0418` と同じ本段 20000 step、g4 は旧 `run_0421`→`run_0422` と同じ本段 20000 + 継続 40000 step。
継続は `restart_field.py` (9 量ビット一致) で同一設定のまま。

| run | 内容 |
| --- | --- |
| `case/46.sern_design/run_0969_3d_g3_chidef` → `run_0971_3d_g3_chidef_cont16k` | g3 段階起動 + 本段 4000 → +16000 (本段計 20000) |
| `case/46.sern_design/run_0970_3d_g4_chidef` → `run_0972_3d_g4_chidef_cont40k` | g4 段階起動 + 本段 20000 → +40000 |

**力係数** (最終 step。`check_quasisteady` は両 run とも 4 量 **STEADY**。振幅: g3 `C_T_with_shear` 2.2e-6 / `C_M` 9.0e-4、
g4 `C_T_with_shear` 1.9e-6 / `C_M` 3.4e-4):

| | 新 g3 `run_0971` | 新 g4 `run_0972` | **Δ (g3→g4)** | §8 許容 | 判定 | 旧 Δ (`run_0418`→`run_0422`) |
| --- | --- | --- | --- | --- | --- | --- |
| `C_T` | 0.8892245 | 0.8891100 | **−0.00011** | 0.002 | ✅ | −0.00013 |
| `C_T_with_shear` | 0.8705591 | 0.8703353 | **−0.00022** | 0.002 | ✅ | −0.00024 |
| `C_L` | 0.0625997 | 0.0639807 | **+0.00138** | 0.002 | ✅ | +0.00145 |
| `C_M` | −1.6444370 | −1.6849318 | **−0.04049** | 0.05 | ✅ | −0.04173 |
| `C_T_friction` | −0.0186653 | −0.0187747 | −0.00011 | — | — | −0.00011 |

**GATES**: 新 g3 `PASS`、新 g4 `PASS` (**`floors` 全項目 0**: ρ ≤ 1e-4 が 0 ノード、ρ min 2.45e-3、P min 801 Pa、T min 94.7 K)。
旧 g4 `run_0422` は同じ段数で `FAIL FLOOR_STUCK` (側壁の壁ノード 8 個が ρ 床 1e-4) だった。
新 g4 は本段 20000 step の時点でも床 0 (旧 g4 は同時点で 13 ノード)。
**床が消えた理由は切り分けていない** (chi 既定化・lsq 化・R5r のメッシャ修正を同時に入れている)。R5s(1) (旧格子で床ノードと開いた CV の重なり) は
旧 run の場が削除済みのため、この床張り付き数の比較で代えた。

**旧列との差** (同じ格子・同じ段数、**chi・lsq・メッシャ修正・出口 BC (`statPress`→`outflow`) の 4 つの合算**。参考値): g3 は `C_T` −0.00005 / `C_T_with_shear` +0.00002 / `C_L` +0.00014 / `C_M` −0.0029、
g4 は `C_T` −0.00003 / `C_T_with_shear` +0.00004 / `C_L` +0.00007 / `C_M` −0.0017。いずれも §8 の許容の 1/10 以下なので、
R5o-chi の手順どおり切り分けはしない。**新しい列 (g1 `run_0973`・g3 `run_0971`・g4 `run_0972`、すべて同条件) を以後の 3D 生産の基準列とする** (g2 は未取得)。

**残差**: `check_convergence` は両 run とも **NOT CONVERGED** (プラトー。g4 継続区間の低下 0.2–0.7 dec、rising なし)。旧 run と同じ性質で、
**3D の固定点 (前 plan V7) は未解決のまま**。上の判定は §8 の許容と準定常 (STEADY) によるもので、残差の収束は主張しない。

**壁解像は測り直していない** (格子は旧列と同じ。y⁺ は §4.42 の値が目安)。

**カウル後縁の低温スポット (R5h) は新しい列にも残る**: g1 なし (T 最小 176.6 K はプルーム遠方) / g3 **105.4 K** (1 節点 < 120 K) / g4 **94.7 K** (1 節点 < 120 K、3 節点 < 150 K)、位置はいずれも (0.12000, −0.01053, 0.09719) = 厚さ 0 の後縁ノードの直下 (外部流側) の内部節点。**格子を細かくするほど冷たくなる**。**上の力係数の格子差が許容内であることは、低温点の影響の上限にはならない** (g3・g4 とも低温点を含み、共通の誤差は差分に出ない。2026-09-27 codex diagnose Major)。

**R5s(2) g1 の取り直し** (`case/46.sern_design/run_0973_3d_g1_chidef`、同じ現行コード・段階起動 + 本段 20000 step、旧 `run_0419` と同条件。
メッシュ SOFT-PASS、`GATES: PASS` (床 0)、4 量 STEADY、残差はプラトー):

| | 旧 g1 `run_0419` | 新 g1 `run_0973` | 差 | 事前の閾値 (§5.1 R5s) | 判定 |
| --- | --- | --- | --- | --- | --- |
| `C_T` | 0.8921423 | 0.8918832 | −0.00026 | 0.002 | ✅ |
| `C_L` | 0.0611365 | 0.0617238 | +0.00059 | 0.002 | ✅ |
| `C_M` | −1.622503 | −1.633609 | −0.01111 | 0.05 | ✅ |
| `C_T_with_shear` | 0.8595643 | 0.8597132 | +0.00015 | (0.002) | ✅ |
| `C_T_friction` | −0.0325780 | −0.0321700 | +0.00041 | — | — |

3 係数とも閾値内。ただし**この差は chi・lsq・メッシャ修正・出口 BC の 4 つの合算**で (上の注記)、B1f の「直接比較しない」に当たるので参考値に留める。
結論は比較でなく**新しい列そのもの**で立てる: 同条件の g1・g3・g4 がそろい、g3→g4 は §8 許容内 (上表) なので、**g 系列の取り直しは g2 を除いて済んだ**扱いにする。
新しい列の g1→g3 は `C_T` −0.00266 / `C_T_with_shear` +0.01085 / `C_L` +0.00088 / `C_M` −0.01083 (2 段分、摩擦の未収束が主) で、旧列 (§4.40) と同じ傾向。

## 5. 実装ステップ

本体 plan §5 の R4 系を引き継ぐ。着手順は §4.15.3 末尾 (codex plan-3 の指定):
**① 接続・物理寸法の確定 → ② ベース込みの BC と帳簿の確定 → ③ `run_0034` の診断訂正と起動試験の定義 → ④ R5b・定量検証条件の確定**。

### 5.1 残作業 (優先順)

**本節が 3D の残作業の正本**。本体 plan §5.1 からはここへのポインタだけを置く。

| # | 項目 | 内容 |
| --- | --- | --- |
| **R4e-次** | **順序を再訂正 (2026-09-19, codex C1/M7)**: ① [`convection-node-wall-reconstruction.md`](convection-node-wall-reconstruction.md) の **W2 (面単位の非物理フォールバック) の A/B** を先に。測定器を直すと**ベース無しの対照は 40 step クリーン、ベース run は step 6 から非物理な再構成**が出るので、SU2 の `bad_recon` 相当が効く可能性がある。→ ② 効かなければ **局所格子** (前回は Δx と壁 CV を同時に変えたのでやり直し。**(i) ベースを横切る点数と (ii) 後流の流れ方向解像を分けて振る**) → ③ 肩の丸めと `t_base` の形状感度。判定量はベース圧が外気の 0.3〜0.5 倍に収まるか (初期値 0.41 倍は妥当だった)。**原因は未確定**として扱う **整理 (2026-10-05、codex diagnose [記録](../../notes/reviews/2026-10-05-sern-next-priority-diagnose.md))**: W2 の A/B は [convection-node-wall-reconstruction](convection-node-wall-reconstruction.md) で実施済み (両設定とも step 28 NaN、「発散を防げなかった」までに限定)。再投入を前提にしない。後流 station の対策は R4e で実装・検証済み。**本行は閉じる** | §4.15.4 / §4.15.7 |
| R5 | **現行バイナリでの再取得 完了 (2026-09-20) → §4.24**: soft 段 step 229 で `roOmega` が**1 節点だけ**非有限。発散節点は壁第 1 層の内部節点で **`k=0`**、初期場で既に **`k=2.06e9` (単独最大)**。`k≤0` が**壁第 1 層に 429 節点**。退化セル帯 (体積 8.48e-12 = 中央値の 1/5000) が**入口∩側壁の角**にある → **R4c の管轄**。次は (b) SST 初期場の明示 patch が最安。旧記述: **3D SST を現行バイナリで再現** (codex M3 採用): 「解決済 (run_0082)」を**撤回**し状態を「加速点で完走、生産点は未成立」に戻す。バイナリ・実効設定を固定して生産 3 作動点を同一幾何で再試験、`L_sw` 分離・鈍頭化はその後の比較対象。通らなければソルバ修正を別 plan 化 **整理 (2026-10-05)**: 「現形状・m6_on の成立 (§4.46 の基準列) と側方境界の限定受理 (farfield plan)」と「未検証の方向 (上方・下方・出口距離、同時拡大)・他作動点・新形状」を分ける。2D 回帰を 3D 全作動点へ転用しない | run_0082 NOT CONVERGED (stalled, 終端 `rms_roOmega` 7.6e16 = 完走の証拠であって収束ではない)、run_0083/0084/0087/0088 DIVERGED (`check_convergence.py` 再確認済)。9/8 の SST 既定 (`solverConfig.hpp` L202) と `scalarTransport_d.cu` L104 の相対ガード変更後なので 9/6 の結果から現行の限界は言えない |
| R4c | **3D 領域トポロジを作り直す (codex C1、最優先)**: `z > W/2` の旧ランプ線〜旧機体上面線を**流体で埋める**、`z = W/2` に**機体側面**を専用タグで新設 (`L_sw` より下流にも及ぶ・ダクト側壁と帳簿を分ける)、交線のノード定義、新領域の組成/初期値を領域情報から与える (`runner_sern3d.py:72` の index 算術を流用しない)。**テストに断面検査を足す** (固体/流体の連結性・内部面の共有・正体積・双対閉性。現行 30 項目は幻の固体を通す) **整理 (2026-10-05、codex diagnose)**: 旧トポロジの修正は R4e (§4.15.3) の実装・検証に対応付けて閉じる。有限厚・隅フィレットの接合部は R5q ([tooling-sern-mesh-blocking](tooling-sern-mesh-blocking.md)、B1c で 245 万節点案まで進行) と分けて扱う | §4.14 / §4.14.1 |
| R4e | ~~**側面バンド末端の壁を除去する (codex C2)**~~ **メッシュ側 完了 (2026-09-19)**。案 (d) を実装 (§4.15.3)。`_vehicle_top_line` はバンド上端 `b(x)` を全 x で返し後縁で厚み `t_base`、`b ≥ yt + t_base` の下限で退化区間を作らない (旧 `i_end` は廃止)。上バンドは j=0 を常に自前で持ち**後縁で担当を切り替えない**。バンドは幅外で全長・幅内で後縁下流のみ。新タグ `vehicle_base` (physID 18)。<br>**生産メッシュ (`case/46.sern_design/run_0200_r4e_mesh/`, 1105796 cells) の実測**:<br><br>| 量 | 旧 | 新 |<br>| --- | --- | --- |<br>| 幅外を横断する壁面 | 224 | **0** |<br>| `vehicle` / `underside_far` 面 | 176 / 44 | **0 / 0** |<br>| `vehicle_base` 面積 | — | 0.020000 H² (= `t_base`×W/2 と厳密一致) |<br>| `y_veh` (形状) | 3.8408936552 | 3.8408936552 (**不変**) |<br>| メッシュ品質 | AR 789.2 / skew 0.401 PASS | AR 789.2 / skew 0.401 PASS |<br><br>**形状を格子から独立させた** (M3): `clearance = max(vehicle_clearance, 3*first_top_frac)` を廃し物理値のみに。粗すぎる格子は `ValueError` で**生成を失敗させる**。旧実効値 0.06 を既定・既存 config 78 件に明示 (2D も同じ扱いにし、2D メッシュがビット一致することを確認)。<br>**帳簿を三分割した** (M4): `vehicle_top`/`vehicle_side`/`vehicle_base` を同一の等温壁にし、`_VEHICLE_FACES` でノズル力から外して機体力に保存 (`C_T_vehicle_side` / `C_T_vehicle_base` を追加)。<br>**試験を強化した** (M5): 幅外横断壁 0・ベース面積照合・**符号付き** Jacobian・float32 の新規衝突 0・格子を変えても形状が動かないこと・粗い格子で生成失敗すること。`run_sern_mesh3d_tests.py` ALL PASS。<br>**着手順 ③ = 2D 有限ベース診断 完了 (§4.15.4、`run_0202`–`run_0206`)**: ベースを **slip から等温粘性壁に変えると SST 1 次段を完走**する (codex M4 の裏付け)。残る壁は**ベース角の 2 次再構成で `P` が 0 に潰れる**ことで、CFL・厚み・リミッタのいずれでも直らない。次の候補は肩の丸め (半径・接点を先に定義) とベース近傍の格子分布。<br>**③ 完了 (2026-09-20)**: 2 次の突破は**メッシュ側**だった — ベース直後の第一 station が厚さの 4.1 倍で後流を 1 セルで跨いでおり、壁 CV から質量が抜けていた。`first_wake_frac ≤ t_base/5` を必須化して根治 (§4.20–§4.23)。起動は `soft_ramp [0.2,0.35,0.5]` + `soft_cfl 0.5` + `cfl_main 1.0`。**`run_0308_wake_soft05` が全段通過・`GATES: PASS`・`C_T_with_shear` 0.9154138・`check_quasisteady` ALL STEADY** (§4.19)。「ベース角の 2 次再構成で P が潰れる」「肩の丸めが次の候補」は**どちらも外れ**だった。<br>**残り = ④ のみ**: 3D 層流 → 3D SST | `mesh_sern3d.py` / `mesh_sern.py` / `runner_sern3d.py` / `tests/run_sern_mesh3d_tests.py` |
| R4f | ~~**`vehicle_top` を等温粘性壁に統一する (codex M3)**~~ **実装済 (2026-09-19)**: 既定を `wall_isothermal` (Ts 1000 K) に変更し機体側面と揃えた。旧 slip は `evaluate.vehicle_top_kind: slip` で比較用に残す。残 = R4e/R5b の後に同一幾何で差を測る。以下は旧記述: 現状 slip。機体側面は等温粘性壁で不整合。slip は比較用に残し、R4e/R5b の後に同一幾何で差を測る。上面の壁距離再生成と近壁解像度も検証対象 | `runner_sern3d.py:75` |
| R5b | ~~**床に張り付いた解を受理しないゲート**~~ **実装済 (2026-09-19)**: `sern_gates.py` に `floor_gate` (EOS 床 `pMin`/`tMin`/`roMin` と乱流下限 `roOmega` 1e-20・`k`≤0 の張り付きノード数) と `residual_scale_gate` (残差列の桁の揃い、中央値の 1e6 倍超を検出) を追加し `evaluate_gates` に接続 (`fail_class` = `FLOOR_STUCK` / `RESIDUAL_UNBALANCED`)。**run_0122 を正しく落とす**: T≤50K 70 ノード / roOmega≤1e-20 **1 ノード** (codex が名指しした ID 646932) / k≤0 28 ノード、残差中央値 3.42e-05 に対し rms_roOmega 8.95e+15。残 = R5 の受入で `require_residual_pass` を必須にする。以下は旧記述:  (codex M2): `sern_gates.py:39` は有限・正値しか見ない。`pMin`/`DEPVAR_TMIN`/`roMin` 到達ノード数、床による保存量変更、床感度を受理条件に足す。**さらに `roOmega` 下限到達と更新クリップも追加し、R5 の受入では `require_residual_pass` を必須にする** (codex C1)。m10 の解析値は 0.147 Pa / 18.0 K / 2.83e-5 で三つの床に抵触 | §4.14-5 |
| R5c-旧 | (§4.26 で優先度低下) **SST 壁関数の壁モデル渦粘性に上限を入れる** (§4.25)。`ν_t,wall = ν(1/g − 1)` が `g→0` で発散し、低密度域 (ν が大きい) で 4660 m²/s に達して `k=2.06e9` を作り 3D SST を殺す。混合長 `κ u_τ y` で上限を掛ける案。**乱流モデルの変更なので別 plan**。診断 (`k_wf` 異常節点の計数) も併せて |
| ~~**R5e**~~ **完了 (2026-09-20, §4.28)**: 板厚を直した `run_0403_3d_cowlfix` と、板厚バグのままの対照 `run_0404_3d_cowlcollapse_ctrl` が**どちらも 20000 step 完走・GATES PASS・ALL STEADY・NOT CONVERGED**。**板厚は R5 の真因ではなかった** (`C_T_with_shear` 差 0.003 %)。真因は `wallTreatmentSST: 1` (= §4.25 の壁モデル渦粘性に上限が無い件) で、壁関数を使わない方針がそのまま R5 を通した | `case/46.sern_design/run_0403_3d_cowlfix` / `run_0404_3d_cowlcollapse_ctrl` (AWS `~/forge-sern`) |
| **R5f** | ~~リミッタのチャタリング説~~ **否定 (2026-09-20, §4.31)**: `venkatK` 100 倍でプラトー 4 % 変化、CFL 半分で 0.74 倍。**場の動きは壁第 1 層に 13 倍濃縮**しており、最大変化点は最悪解像の `sidewall_in` 直近。**プラトーは壁解像の問題** → **R5d に統合**。壁解像メッシュでプラトーが下がることを R5d の受入条件にする。下がらなければ近壁の非定常性と float32 ノイズ床を分ける | `run_0408_main_ref` / `run_0409_main_cfl025` / `run_0410_main_k05` / `run_0411_main_k0005` |
| R5d | ~~3D を壁解像 (y⁺≈1) メッシュに作り直す~~ **段階 1 完了 (2026-09-20, §4.32/§4.33)**: y 法線壁は y⁺ 0.70–2.35 を達成 (1.92M 節点・AR max 4819 PASS)。本段 CFL は 0.5→**0.25** に下げる必要あり。**低下桁数は 0.1 → 0.4–1.2 dec に改善したが `NOT CONVERGED (stalled/plateau)` のまま**。**段階 2 (側壁の z 局所クラスタ化) は優先度低** — 残った動きは側壁でなくランプ上のせん断層 | `run_0413`–`run_0416` |
| ~~**R5g**~~ | **打ち切り (2026-09-20, ユーザ決定 §4.39)**: 残差プラトーは受理の阻害要因としない。力は 7 桁で定常・派生量は ALL STEADY・場の動きは再現性幅の 3–62 倍で、「もとから変動が小さい」と判断。**`check_convergence` の `NOT CONVERGED` は消さず併記する** | 記録のみ |
| ~~**R5n**~~ | **完了 (2026-09-20, §4.40)**: 生産条件・同一形状で壁法線解像を 4 倍ずつ振った 3 点。**`C_T` (圧力) −0.00089 と `C_L` +0.00116 は許容内**、`C_T_with_shear` +0.00318 と `C_M` −0.02219 が外れる。**原因は摩擦** (`C_T_friction` が比 2.40 で動き続ける、収束次数 0.63、外挿残差 ~0.0029)。= **壁解像の緩和は根拠づけられない** | `run_0418` / `run_0419` / `run_0420` |
| ~~**R5o**~~ | **完了 (2026-09-21, §4.41/§4.42)**: x ブレンドを実装 (既定は挙動不変・座標ビット一致、急なブレンドは `df/dx > 4e-3` で生成失敗)。メッシュ **SOFT-PASS** (AR>5000 が 0.02 %・skew 0.701)。**摩擦が収束し `C_T_with_shear` が −0.00024 で許容内**、`C_T`/`C_L` も許容内。y 法線壁の y⁺ は `ramp` 0.549 / `cowl_in` 1.234 | `run_0421` / `run_0422` |
| ~~**R5p**~~ | **決着 (2026-09-21, ユーザ決定)**: **`C_M` の許容を 0.02 → 0.05 に緩めた** (§8 に反映)。腕が 20–29 H あるので $|\Delta C_L| \le 0.002$ と $|\Delta C_M| \le 0.02$ は同時に満たせなかった。**厳密な整合は 0.06** なので、$\Delta C_L$ が 0.002 に近い case が出たら見直す。これにより **g3→g4 の `C_M` −0.04173 は許容内**になり、**4 係数すべてが許容内**に入った | plan §8 |

| **R5q** | **別 plan へ移管 (2026-09-21, codex plan NO-GO + 方式相談を受けて)**。§4.44 の「座標だけ動かす」設計は**撤回**。接合部のトポロジ (露出カウル側端の壁・側壁下端のノード定義・`ext_top`/`vehicle_side` の座標写像・`L_sw` の物理 station 化・厚み 0 の互換分岐・`half_W_m` と入口項の帳簿・新モード試験一式) を明示設計する**新 plan** を立てる。方式は ~~(a) 自作メッシャ改修で確定~~ **保留に戻した (2026-09-21)**: ユーザから「フィレットを付けてメッシングすることも想定」との要件が後出しで出た。codex 相談はフィレット無しの前提だった。**輪郭面内 (x–y) のフィレットは自作で可** (`ramp_fillet` 実績) だが、**断面内 (y–z) の隅フィレット (側壁×ランプ/カウルの接合隅) はテンソル積構造が壊れ O グリッド等のマルチブロックが要る** = 自作の優位が消える。第三案 **gmsh-OCC + Python API** (msh4.1→`convertGmshToForge` の既存 node 経路に乗る) も候補に追加。論点は y⁺≈1 の 3D 境界層を接合隅で積めるか (case/49 で tet+VL は系統細分不可の実績) と node の混合要素可否。~~フィレットの種類をユーザに確認中~~ → **ユーザ回答 (2026-09-21): 側壁×ランプ/カウルの隅の丸めを含む**。codex 再相談の結果 **(d) gmsh Python API 全ヘキサ・マルチブロック**を推奨 (記録 `2026-09-21-codex-mesh-approach-2.md`)。主方向 x の transfinite ブロック列、側壁終端は専用 3D 接続ブロック、露出カウル側端は独立した物理壁、壁が終わっても格子帯を潰さない。**最初の実証 = 側壁終端を含む 3D 接続模型 (第一層 4 µm・約 40 層、3 解像度、無人)**。ユーザの方式決定待ち。**本 plan の 3D 受理は「厚さ 0 交線に床張り付き 12 ノード / 320 万節点、格子細分で悪化」を既知の限界として明記して締める** **整理 (2026-10-05)**: R7a と並行・後続 (新メッシャ完成まで判別試験を止めない)。移管先 plan は B1c で 245 万節点案まで進んでおり、「1250 万節点で載らない」だけでは現況を表さない | 新 plan (未作成) |
| ~~R5r~~ (済 2026-09-22: §4.45。格子の段階で確認。開いた CV 18 → 0、`VERDICT: PASS`) | **旧メッシャ: カウル板の自由端の二重節点を共有にする** | `share_cowl_free_edge` (既定 True)。回帰試験 4 件追加 |
| ~~**R5s**~~ (**完了 2026-09-27 §4.46**: (1) 代替で済 = 新 g4 の床張り付き 0。(2) 新 g1 `run_0973` を取得。旧 `run_0419` との差 C_T −0.00026 / C_L +0.00059 / C_M −0.01111 は閾値内だが出口 BC の変更を含む参考値。**新しい列 g1/g3/g4 (同条件) を基準列とし、g2 以外の取り直しは済**) | **R5r の修正の CFD への効き目を測る (AWS)** | (1) `check_dual_closure.py` を g4 `run_0422/sern.h5` に掛け、R5q の床ノード 12 個が開いた CV と重なるかを見る。(2) g1 (`problem_3d_prod_m6on_g1.yaml`) を修正後の格子で回し直し、C_T・C_L・C_M と ρ min の差を出す。差が許容 (0.002 / 0.002 / 0.05) を超えるなら g 系列を取り直す |

| ~~R5j~~ | **R5q に統合 (2026-09-21)**: カウル側端テーパ `sz` を廃すことで、「最後の 2 z-セル」で物理幅が決まる問題 (`nz_in` 25→57 で 0.868 → 0.00949 mm、91 倍) が消える | §4.44-3 |


| **R5h** (**再開 2026-10-05**: exact A 形状の m10_on 3D で後縁の 1 節点が T 50 K 床に張り付き評価不成立 → 診断 A/B を本体 plan [tooling-nozzle-sern-chain](tooling-nozzle-sern-chain.md) §5.1 R7b 4-2 に事前登録。以下は旧記述: **整理 2026-10-05: 診断と現形状の力感度試験 [run_0984/0985、16000 step] は完了、恒久対策は別件で保留。新形状で低温域・床到達が増えたら再開** (codex diagnose 2026-10-05)。**F**。**判断: 2026-09-27 codex (diagnose) [`notes/reviews/2026-09-27-sern3d-r5h-te-coldspot-diagnose.md`](../../notes/reviews/2026-09-27-sern3d-r5h-te-coldspot-diagnose.md) — 現形状で診断を続ける。有限厚後縁 plan への移管・監視項目への格下げは不採用**。**次の A/B (codex 指定)**: g3 `run_0971` の最終場から `restart_field.py` でビット一致の 2 run、A = `convMethod` 1 (現行)・B = 0、他は固定、各 100 step。主判定は最初の共通状態での残差: 低温点・接する後縁壁点・隣接 CV の保存量別の更新量、EOS 前後の保存量と T。分岐 A = 再構成差が低温点を加熱側へ動かす更新差を説明 (EOS では説明されない) → 再構成依存を支持。分岐 B = 説明しない → 第 1 仮説棄却。収支が閉じなければ判定不能。**一次化で最低温度が上がっただけでは原因確定にしない**。やらないこと: 温度床を上げる・後縁だけ粗くして受理・新旧メッシャの力差を影響上限と呼ぶ。2026-09-27 時点の現況は §4.46。**A/B 実施 (2026-09-27, run_0974–0977)**: 共通状態 (両 res_0 ビット一致) で低温点 (節点 517159) の 1 更新は A (conv 1) で T 不変 (res_roe −2.4e-4)、B (conv 0) で **+10.6 K** (res_roe +6.36)。100 step で A は 105.7 K 不動、B は 535 K まで単調増加・未飽和。EOS 検算の差 0.32 K、T・P は床より十分上。**判断 (結果): 2026-09-27 codex (diagnose) [`notes/reviews/2026-09-27-sern3d-r5h-ab-result-diagnose.md`](../../notes/reviews/2026-09-27-sern3d-r5h-ab-result-diagnose.md) — 「再構成への感度を確認」までは記録してよい。EOS の除外と「A は偽の定常」の認定は保留** (保存量から T を再現できても EOS 処理で保存量が変わっていない証明にならない / `limiter_T` 0 は ρ 再構成経路では面温度の 1 次化を意味しない / B は未飽和の過渡。また (dt/vol)·res_ro = +8.2e-4 に対し実 Δro = −1.6e-3 で残差と更新を直結した収支は閉じていない)。**次の手 (codex 指定)**: 同じ A/B (各 1 step・`nStepInner` 1) に**局所の帳簿ダンプ**を足す — 低温点・接する後縁壁点・隣接 CV について、各接続面の両側状態 (ρ・P・速度・組成・面温度・全エンタルピー・χ)・質量/運動量/エネルギー流束、対流・粘性・ソース・境界を分けた残差、陰解法の保存量増分と EOS・床・境界処理による追加変更。**判定基準 (事前)**: 収支の不一致が「A/B 差の 1 % または丸め誤差上限の大きい方」以下で閉じたとき分岐 A/B を判定、閉じなければ判定不能)。**局所帳簿 実施 (2026-09-27, run_0978/0979、診断ダンプ `FORGE_DUMP_LEDGER` を追加 — 出力専用で massflux・状態が旧バイナリとビット一致 [run_0980/0981])**: (1) 面流束の和が対流残差を再現 (ρ 差 ≤ 3e-13、ρe 差 ≤ 2.3e-7 = A/B 差 6.36 の 4e-8) = **収支は閉じた**。(2) 段別: EOS・境界処理による保存量の変更は 0 (両 call)、ソース 0、粘性 +3.874 は A/B 共通、**A/B 差は対流だけ** (conv A −3.875 / B +2.483)。更新は A Δe_int −0.17 J/kg、B +7866 J/kg (EOS 後 T 116.1 K)。→ **事前の規則で分岐 A 成立 (再構成依存を支持、EOS 単独原因は除外)**。(3) 面別: B−A 差 ρe +6.36 のうち **+5.76 が +x 面 1 枚** (相手 522112 = 後縁の後流側、|U| ≈ 1500 m/s)。その面で A は低温点側の再構成速度が **セル 10 m/s → 1333 m/s** (limiter_Ux 0.28)、ρ は 1 次 (limiter_ro 0) で 0.305 のまま、P_L 14023 (セル 9367)。運び出す全エンタルピー h_p 7.98e5 J/kg (セル −1.78e5)、質量流束は B の 5.2 倍。A では対流 −3.875 と粘性 +3.874 が釣り合う。**これは観測であって原因の確定ではない**。**判断: 2026-09-27 codex (diagnose) [`notes/reviews/2026-09-27-sern3d-r5h-ledger-result-diagnose.md`](../../notes/reviews/2026-09-27-sern3d-r5h-ledger-result-diagnose.md) — 分岐 A は今回の初回応答に限って採用 (形成の履歴まで EOS 無関係とは言わない。面和→残差の閉合と陰解法の更新収支は別)。「速度外挿が原因」「粘性加熱と釣り合う 105 K の定常」は要再検証 (全域 1 次化は両側の ρ・P・U を同時に変える / res_roe の粘性寄与は内部エネルギー加熱量でない / Venkat は全近傍の成分別 min/max で制限するので「相手値を超えないから正常」は除外根拠にならない)。****次の A/B (事前固定)**: 面 1535830 の低温点側の速度 3 成分だけをセル値に戻す (α=0) 診断介入 vs 現行 (α=1)、同一 restart・各 1 更新・`nStepInner` 1、相手節点 522112 も記録。合格 = 当該面の差 5.76 の **80 % 以上**を同方向に再現し、内部エネルギー射影 ρV·de/dt = R_ρE − u·R_ρu + (|u|²/2 − e)R_ρ と実更新も加熱側へ動く → 「局所残差を維持する速度再構成の機序」を採用 (正しい温度・恒久対策は認定しない)。それ以外 → 第 1 仮説棄却、収支が 1 % で閉じなければ判定不能。やらないこと: 有限厚化を原因の証明に使う・h0 クリップや共通リミッタを直ちに全ケースへ・(c) 監視への格下げ (力への影響は未測定)。**1 変数介入 実施 (2026-09-27, run_0982 α=1 / run_0983 α=0、診断介入 `FORGE_DIAG_FACE_VEL_CELL`)**: massflux が変わったのは面 1535830 の 1 面だけ。当該面の ρe 寄与 −5.52 → +0.15、差 +5.67 = **既報 5.76 の 98.5 %** (基準 ≥ 80 %)、内部エネルギー射影 −1.4e-4 → **+6.80** (全域 1 次化 +6.70)、実更新 Δe_int −0.17 → **+7945 J/kg**、EOS 後 T 105.37 → 116.23 K (全域 1 次化 116.12)、収支の閉じ ≤ 2.3e-7。→ **事前の規則で「局所残差を維持する速度再構成の機序」を採用** (後流側の面で低温点側の速度 3 成分の外挿が、全域 1 次化の効果をほぼすべて再現する)。正しい温度・恒久対策は認定しない。**力の感度 (codex 指定の最小評価) 走行中**: run_0971 最終場から同一設定で 16000 step、介入維持 `run_0984_r5h_sens_alpha0` / 対照 `run_0985_r5h_sens_alpha1` (同じ診断ビルド)。判定は両者 `check_quasisteady` STEADY の上で平均差に双方の変動幅を加えた幅 (影響上限ではなく感度)。**結果 (2026-09-27)**: 両者 `GATES: PASS`、4 量 STEADY。介入維持では**低温スポットが消える** (T 最小 176.3 K = プルーム遠方、T<150 K は 0 節点。対照は 105.4 K のまま)。感度幅 (|平均差| + 双方の振幅): `C_T` 6.0e-6 / `C_T_with_shear` 5.7e-6 / `C_L` 2.0e-5 / `C_M` 4.4e-4 (平均差そのものは 3e-7 / 3e-7 / 2e-7 / 2.9e-6) = **§8 許容の 1/100 以下**。これは「この 1 面の介入に対する力の感度」であって誤差の上限ではない (codex の注記どおり)。**残り (未着手・判断待ち)**: 恒久対策は再構成の変更 (全ケースに効く数値カーネル) なので別 plan + 上位の判断 + codex plan 段が要る。R5h はここで止め、ユーザ指示の順 (R5h → R5j → R4d) で R5j へ進む。診断 run (run_0974–0983) の場は削除済み (帳簿 CSV は残す) | **カウル後縁の非物理な低温スポット** (§4.35)。x = `L_cowl`・側壁直内・壁第 1 層で T 37.6 K / ρ 0.671 (1 station 手前は 801.9 K で正常)。**床張り付きではない** (訂正済) が非物理。後縁は**同一ノード ID の正しく閉じた後縁**なので、原因は「1 本の壁ノードが 2772 K と ~1000 K を同時に受け持つこと」。まず **(a) 当該節点の保存量・T の時系列と EOS 前後の保存量変更**を見て、**(b) 断熱壁の回復温度がそこだけ崩れているのか、対流流束の不整合か**を分ける。**有限厚後縁 (案 a) の実装はその後** — codex は推奨しつつ、接続設計 (後縁厚・上下輪郭・後流ブロック・側端の閉じ方・共有ノード ID) を**図と式で先に固定**し、`cowl_base` をノズル力に含め、BC・壁距離・初期場・圧力/摩擦/モーメント帳簿を定義してからにせよ、と指定 | `run_0415` の節点 517197/517198 |
| ~~**R5i**~~ | **完了 (2026-09-20, §4.37/§4.37.1)**: 生産条件 (TP + 等温壁 1000 K) を 3D に入れ直し `run_0418_3d_prod_wallres` で取得。**プラトーは消えない** (低下 0.6–1.5 dec、断熱 CPG と同程度)。後縁の低温スポットは **37.6 → 109.6 K**・影響 8 → 1 節点と改善するが**残る**。**冷却壁は y⁺ を ×3.1 悪化**させるので、断熱で設計した壁解像メッシュは生産条件では解像不足 | `run_0418_3d_prod_wallres` |
| **R5m** | **AR を下げる道の費用を測る (codex plan-2 M6)**。生産条件で `cowl_in` を y⁺ 0.7 にするには `first_wall_frac` 1.6e-04 → **2.3e-05** が要り、AR は 4213 → **~29000**。AR は最長辺/最短辺なので、**遠方の長辺 (下バンドの最粗 y 刻み 1.9–2.3 H、外側 z 0.77 H) を割れば下がる**。節点数と RAM の費用を測り、y⁺≈1 が現実的かを数値で決める。測るまで「不可能」とは言わない **整理 (2026-10-05)**: R7a (本体 plan の 3D 判別試験) と並行・後続。ただし側壁解像不足のまま最終順位を認定しない。z を細分するときは `sz` テーパによる形状変化を分離する | ローカルの AR プローブ (`probe_ar.py` 相当) で先に見積もる |
| ~~R5l~~ | **完了 (2026-09-20, §4.38)**: `check_wall_resolution.py` が多成分で判定不能になる件を、体積出力の `vis_lam` を優先使用する形で解消。CPG 回帰は一致 | `solver_density_cuda/tools/check_wall_resolution.py` |

| ~~**R5j**~~ (**重複行を閉じる 2026-09-27**: 上の行のとおり 2026-09-21 に R5q へ統合済みで、R5q は全ヘキサ接続模型の plan [`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md) に移管されている。**現行メッシャ `mesh_sern3d.py` には `sz` テーパが残る**が、生産の格子列 g1–g4 は `nz_in` がすべて 25 なのでテーパ幅は列の中で一定 (格子収束の判定 §4.46 には交絡しない)。`nz_in` を変える比較をするときだけ効く) | **カウル側端テーパを格子から独立させる (M5)**。`kt = max(k_sw - 2, 0)` が「最後の 2 z-セル」で物理幅を決めており、`nz_in` 25→57 で **0.868 mm → 0.00949 mm (91 倍)**。物理入力に移し、**粗細で輪郭と面積が一致する試験**を足す。有限厚後縁を入れるとベース面積にも及ぶ | `mesh_sern3d.py` L198–206 / `run_sern_mesh3d_tests.py` |
| ~~R5k~~ | **完了 (2026-09-20)**: `sern_gates.py::_solver_floors` が経路を見ずに温度下限を一律 50 K にしていたのを、`thermalMethod`/`condensation` を見て CPG は `tMin` (既定 1e-4 K)、TP/凝縮だけ `DEPVAR_TMIN` 50 K にする形へ修正。`run_0415` が `FLOOR_STUCK` → **`GATES: PASS`** | `design/forge_design/metrics/sern_gates.py` |

| R5c | ~~SST 壁関数の壁モデル渦粘性に上限~~ **優先度低下 (2026-09-20)**: 壁関数を使わない方針になったため。欠陥の記録は §4.25 に残す。将来 `wallTreatmentSST: 1` を復活させるなら必須 |
| R4d (**F**。**判断: 2026-09-27 codex (diagnose) [`notes/reviews/2026-09-27-sern3d-r4d-domain-design-diagnose.md`](../../notes/reviews/2026-09-27-sern3d-r4d-domain-design-diagnose.md) — 6 変種の一括投入は不可、まず g3 で `Z_ext` 1.5 → 2.25 の 1 組だけ**。**事前固定**: (i) 認定は g3 (g1 は予備のみ: g1→g3 の圧力 C_T 差 0.0027 が領域許容を超える)。(ii) 共通領域の座標・接続を保持し、追加領域にだけセルを足す (現行メッシャは寸法から分布を再生成するので、そのままの比較は不可)。(iii) BC は両方とも現行のまま。A = 同一格子 restart、B = cross-mesh restart (共通部の補間差を記録)。(iv) 各 20000 step・500 step 間隔、末尾 10000 step で評価し、その前 10000 step との平均差が各量 0.1ε 以下 (未達なら延長)。(v) **領域感度の許容 ε**: `C_T`・`C_T_with_shear`・`C_L` 0.0005、**`C_M` 0.005** (§8 の 1/4 だと C_M は格子差 0.0405 と合わせて 0.05 を超える)。比較量 D = |平均_B − 平均_A| + a_A + a_B、a = max|値 − 平均|。最終認定は格子差の幅 G と合わせて G + D ≤ §8 の総許容。(vi) 両者 `check_quasisteady` STEADY・床/NaN・メッシュ品質が必要条件、`check_convergence` 併記。(vii) 分岐: 1 量でも D > ε → 現側方領域は不合格 / 全量 D ≤ ε → この拡大に対する感度は許容内 (無限遠への独立性ではない)。(viii) 補助: 各遠方境界 (side_far・top_out・bottom・outlet) の法線 Mach Mₙ = U·n/a 分布・流入面積率・組成。**やらないこと**: `side_far`/`top_out` の slip → outflow を「非反射化」として生産採用 (outflow は全量外挿で非反射の保証ではない。しかも runner は `top_out_kind: outflow` を slip に落とし、`side_far` は slip 固定 = 生成器の対応が先)、g1 の無感度を g3 へ転用、変種ごとに独立に誤差予算を配る。V3 (Z_ext) の後に V4–V6 (top_depth・bot_depth・x_out_extra) を同条件で、最後に最終 BC を固定した同時拡大領域の確認と m6 以外の作動点)。**V3 実施中 (2026-09-27)**: `mesh3d.z_append` を追加 (既定 0 でビット不変、共通領域の座標・hex を保つ試験つき)。A `run_0986_r4d_g3_base` (run_0971 最終場を同一格子 restart、20000 step) / B `run_0987_r4d_g3_zapp075` (`z_append` 0.75 → 最外セル 0.92 H を 1 層追加、遠方面 z 2.5 → 3.42 H、1,939,968 セル、品質 PASS)。B の IC: 共通 1,920,103 節点は `case/46.sern_design/r4d_common_restart.py` で保存量を index コピー (ビット一致)、追加 102,286 節点は `interp_field` の最近傍。**途中で `interp_field.py` の欠陥を発見・修正**: 照合座標を常に x,y だけで取っていたので 3D では z を無視した別節点の値を貼っていた (初回 B は共通節点の ρ が最大 55 倍違い異常終了)。3D→3D は x,y,z、2D→3D (押し出し) は x,y、3D→2D は拒否に変更、試験 `solver_density_cuda/tools/test_interp_field_3d.py` (旧コードで FAIL・新コードで PASS)。過去に 3D で `interp_field` を使った run (case/16 run_0219 は index コピーで回避済み、case/49 は使用しない注記あり) は、この欠陥を含むので読み直しが要る)。**V3 結果 (2026-09-27、事前の規則で判定)**: B は C_M の窓条件未達 (前後窓の平均差 1.2e-3 > 0.1ε) のため規則どおり `run_0988_r4d_g3_zapp075_cont20k` で +20000 延長し、窓条件は全量 OK。A `run_0986` / B `run_0988` の末尾 10000 step: `C_T` 0.8892228 / 0.8892285 (D 1.1e-5 ✅)、`C_T_with_shear` 0.8705567 / 0.8705795 (D 2.8e-5 ✅)、**`C_L` 0.0625905 / 0.0618810 (D 7.4e-4 > ε 5e-4 ❌)**、**`C_M` −1.6441445 / −1.6228715 (D 0.0218 > ε 0.005 ❌)**。両者 GATES PASS・4 量 STEADY・床 0、`check_convergence` NOT CONVERGED (プラトー)。→ **分岐 A: 現側方領域 (遠方面 z = 2.5 H) は C_L・C_M で領域独立でない** (この拡大への感度が許容超。無限遠への距離は未知)。格子差と合わせると C_M 0.040 + 0.022 = 0.062、C_L 0.0014 + 0.0007 = 0.0021 で §8 の総許容 (0.05 / 0.002) も超える。補助: 差はほぼランプ面 (ΔFy −3.77 N / 1735 N)、カウルは 0.04 N。`side_far` (slip) 上で |Mn| > 0.05 の節点が 18 % (他タグと共有しない面内節点、最大 Mn +0.61 @ x/H 4.1 = プルーム域)、`top_out` (slip) も同様 (最大 0.44)。**原因 (側方反射か等) は未判定**。**判断 (結果): 2026-09-27 codex (diagnose) [`notes/reviews/2026-09-27-sern3d-r4d-v3-result-diagnose.md`](../../notes/reviews/2026-09-27-sern3d-r4d-v3-result-diagnose.md) — 分岐 A を採用 (「現 BC・g3 でこの側方拡大への感度が許容超」まで。拡大側の十分性・反射が真因とは認定しない)。****3D の C_L・C_M は当面、探索用の参考値に留め、設計の合否・最終順位の確定には使わない** (格子差との単純和でも C_L 0.00212・C_M 0.0623 で総許容超)。節点 Mₙ ≠ 0 は slip を通る流入の証拠ではない (slip は境界状態で法線成分を除いた速度から流束を組む: `boundaryCond_d.cu:68`・`convectiveFlux_boundary_d.inc.cuh:192`)。「流入 51.7 %」は内部節点の速度の符号の割合であって流入率ではない (訂正)。**次の手 (codex 指定)**: 幅の系列より先に、g3 同一格子・遠方面 z 2.5 H のまま `side_far` だけ slip → **検証済みの特性型遠方境界** (流入特性に外気状態、流出特性は内部) の A/B を 1 組 (top_out・bottom・outlet は維持、同じ新バイナリで A を対照、各 20000 step、同じ ε・窓条件)。結果 A (C_L か C_M で D > ε) → 側方 BC の影響を支持、BC を確定してから領域系列へ。結果 B (両方 D ≤ ε) → slip 閉塞だけでは説明しない (旧境界節点が内部になる双対 CV の変化などを候補に残す)。順序: 側方 BC の診断 → 遠方 BC の確定 → 最終 BC で V3–V6 と同時拡大確認。V4–V6 を現 BC で回すのは保留。**前提: forge に特性型遠方境界は無い** (bcondKind は slip / outflow / inlet_* / outlet_statPress / wall* / periodic / axis)。新規 BC の開発 (TP・組成・SST 整合、自由流保持と波の出入りの検証) = 開発フロー (methods → 新 plan → codex plan 段) が要る。~~着手はユーザ判断待ち~~ **決着 (2026-09-27 ユーザ決定「広げる」)**: 特性型 BC の開発より先に、**側方を広げる系列で C_L・C_M が落ち着く幅を実測する** (codex の推奨順 [BC 先行] とは異なるユーザ判断)。**事前固定**: `z_append` 1.8 (最外 0.92 H を 2 層、遠方面 z 4.34 H) と 3.5 (4 層、6.18 H) を追加し、既存の 2.5 H (run_0986) / 3.42 H (run_0988) と合わせて 4 点。BC は現行 (slip) のまま、g3 同一条件、各 20000 step・同じ ε (C_T・C_T_with_shear・C_L 5e-4、C_M 5e-3)・窓条件 (前後 10000 step の平均差 ≤ 0.1ε、未達なら延長)。IC は run_0988 最終場から共通領域を `r4d_common_restart.py` で index コピー + 追加層は最近傍。**判定**: 隣り合う幅の組で D ≤ ε (C_L・C_M 両方) となる最初の組の**狭い側**を必要幅とする。6.18 H まで満たさなければ「slip のままでは 6 H でも収束しない」と記録し、特性型 BC の開発を再検討する。slip のまま広げても反射そのものは無くならない (遠くなるだけ) ので、収束しても「無限遠と一致」とは言わない)。**幅系列の結果 (2026-09-27、事前規則で判定)**: 4 点とも GATES PASS・4 量 STEADY・床 0・窓条件 OK (`check_convergence` はいずれも NOT CONVERGED = プラトー)。末尾 10000 step 平均: C_L 0.0625905 (2.50 H, run_0986) / 0.0618810 (3.42 H, run_0988) / 0.0618842 (4.35 H, run_0989) / 0.0618870 (6.19 H, run_0990)、C_M −1.6441445 / −1.6228715 / −1.6229336 / −1.6229828。隣り合う組の D: 2.50→3.42 **不合格** (C_L 7.4e-4、C_M 0.0218)、**3.42→4.35 合格** (C_L 2.5e-5、C_M 5.5e-4)、4.35→6.19 合格 (C_L 2.3e-5、C_M 4.8e-4)。C_T・C_T_with_shear は全組 D ≤ 3e-5。→ **事前規則で必要な側方幅 = 遠方面 z 3.42 H** (`Z_ext` 1.5 + `z_append` 0.75 = 最外セル 1 層追加)。3.42 H から先の変化は許容の 1/9 以下で頭打ち。**言っていないこと**: slip のままなので無限遠との一致ではない / 2.50→3.42 だけが大きく動く理由 (旧境界付近の離散化の変化 [codex 第 2 仮説] か、プルームが境界に当たるかどうかの閾値か) は未判定 / 生産の格子列 (g1/g3/g4) は 2.50 H で取ったので、C_L・C_M の格子差 G を 3.42 H で取り直すかは未決 (C_T・C_T_with_shear は領域の影響が 3e-5 以下で格子列の結論は不変)) **追記 (2026-09-29、farfield plan §6 V3)**: side_far を特性型遠方境界 `farfield` にすると、2.50 H で 3.42/4.35 H の farfield と全 4 量 D ≤ ε (試験系列内の限定付き、初期場履歴 A/B 走行中)。slip の必要幅 3.42 H は当時の設定系列 (自動リミッタ基準値・旧バイナリ) での結果として残す。生産の側方 BC 切替はユーザ判断で、前提に farfield での格子差 G の取り直しが要る ([farfield plan](../accepted/boundary-node-farfield-characteristic.md) §5.1 #4c)。**決着 (2026-10-01)**: G を取り直し (farfield で G + D が §8 内) のうえ、ユーザ決定で生産 3D YAML の `side_far_kind` を farfield に | **幅外を開いた後の遠方境界・領域独立性** (codex M4): `side_far` の slip 固定 (反射) を見直し、生産 TP・SST で $C_L/C_M$ まで含めた領域独立性を測る。許容値を係数ごとに数値で固定 | `r4_domain_study.py` は加速点 Euler・$C_T$ のみ |
| ~~**R5o-chi**~~ (**完了 2026-09-27 (§4.46)**: 新 g3 `run_0971`/新 g4 `run_0972` とも `GATES: PASS`・4 量 STEADY、Δ(g3→g4) `C_T` −0.00011 / `C_T_with_shear` −0.00022 / `C_L` +0.00138 / `C_M` −0.04049 で**全て §8 許容内**、旧 g4 の FLOOR_STUCK は**床 0 に** (理由は未切り分け)。旧列との差は許容の 1/10 以下。残差はプラトー (V7 は未解決)。当初の範囲: 現行コード (chi 既定 auto = 1 [plan convection-slau-wall-normal-chi-default accepted]・スカラー勾配の node 既定 lsq [plan gradient-scalar-lsq-unification accepted]・R5r のメッシャ修正 `share_cowl_free_edge` 込み) で g3 `problem_3d_prod_m6on_wallres.yaml` と g4 `problem_3d_prod_m6on_g4.yaml` を runner_sern3d で回し直し、g4 は旧 run_0422 と同じく定常化まで継続。判定は §8 の許容 (Δ(g3→g4) C_T・C_L 0.002、C_M 0.05)、GATES (特に旧 g4 の FLOOR_STUCK = 側壁 8 節点の床が消えるか)、`check_quasisteady`。旧列 (run_0418/0422) との差は chi・lsq・メッシャ修正が混ざるので、差が許容を超えたときだけ切り分ける。R5s(1) の旧格子での照合は、旧 run の場が削除済みで床ノード ID が得られないため、新 g4 の床張り付き数で代える) | **R5n 型の格子収束列を `space.slauWallNormalChi: 1` で再取得** (2026-09-25 委譲、元: [`convection-slau-wall-normal-chi-usage-rule.md`](../accepted/convection-slau-wall-normal-chi-usage-rule.md) §5.1 #7)。3D 生産は flag 0 で解を持たず常に flag 1 なので、決めるべきは flag 1 の生産解の格子収束。メッシュ生成コマンドと `FORGE_CUDA_BLOCKSIZE` を run に記録。**3D の固定点 (前 plan V7 のドリフト減衰) も未解決のまま**ここに属する | 許容は §8 (R5n)。前 plan [`convection-slau-wall-normal-chi.md`](../accepted/convection-slau-wall-normal-chi.md) §5.1 #9b・V7 |

## 6. 検証

§4.15.3 の「実装前に固定する検証 5 段」(接続 / 幾何 / 起動 / 受入 / 独立性) に従う。
メッシュ側の機械的検査は `design/tests/run_sern_mesh3d_tests.py`。

### 6.1 レビュー記録 (codex)

| stage | 日付 | 記録 | 判定 | 採否 |
| --- | --- | --- | --- | --- |
| 相談 (再) | 2026-09-21 | [2026-09-21-codex-mesh-approach-2.md](../../notes/reviews/2026-09-21-codex-mesh-approach-2.md) | **推奨を (a) → (d) gmsh Python API 全ヘキサ・マルチブロックへ変更** | **採用を推す (最終決定はユーザ)**。前提変更 2 点: ①ユーザ要件に**断面 (y–z) の隅フィレット**が入った ②当方の「Salome→forge node 経路は無い」は**誤り** (case/37 に実在)。隅フィレットで (a) も本格的なマルチブロック設計になり自作の優位が消える。(d) は case/49 に生産実績 (O グリッド/バタフライ・`--scale` 細分列・AWS 生成) があり msh4.1 で既存 node 経路に乗る。(b) Salome VL は case/49 の実測 (VL 総厚と接線サイズの結合で細分列が作れない) で棄却、4 層の成功は 40 層の裏付けにならない。M1–M4 は共有幾何・接続設計で消えるが **M5/M6/M7 は残る**。**方式確定の条件 = 「側壁終端を含む 3D 接続模型」の無人生成が品質ゲートを通ること** (未実証なので鵜呑みにしない) |
| plan | 2026-09-21 | [2026-09-21-tooling-nozzle-sern-3d-plan.md](../../notes/reviews/2026-09-21-tooling-nozzle-sern-3d-plan.md) | **NO-GO**, C0/M7/m1 | **全件採用 → §4.44 の設計を撤回**。M1 `L_sw`(0.8)<`L_cowl`(1.2) で露出カウル側端を閉じる面が無い (未分類面 12)。M2 `dup2` は `j>jm` のみで側壁下端に厚みを与えられない。M3 `_add_ext_top3d`/`_add_vehicle_side3d` も `zs[k]` を使う。M4 `i_sw=argmin(|xs−L_sw|)` で `L_sw` が station に無い。M5 厚み 0 でも `sz` ステップ化で 1108 節点が動きビット一致しない。M6 `half_W_m` 固定で入口項と理想推力に系統誤差。M7 既存試験は新モード未覆。**規模が 3〜5 倍**なので R5q を別 plan へ移管 |
| 相談 | 2026-09-21 | [2026-09-21-codex-mesh-approach.md](../../notes/reviews/2026-09-21-codex-mesh-approach.md) | **(a) 自作メッシャ改修を推奨** | **全面採用**。CAD+Salome に移っても **M5/M6/M7 は消えない** (とくに入口帳簿はメッシャ非依存)。**Salome→forge の node 経路が存在せず** (`fluent_h5_to_forge.py:442,498` は cell 中心で代用不可)、MED/UNV→msh4.1→`convertGmshToForge` の新設と受入試験が要る。無人生成の成功率も未検証。**(1)(2)(4) の本質は「物理形状と接続を格子索引から独立させる」ことで (a) でも可能**。→ (a) で進め、接合部のトポロジを明示設計する別 plan を立てる |
| plan | 2026-09-20 | [2026-09-20-tooling-nozzle-sern-3d-plan-2.md](../../notes/reviews/2026-09-20-tooling-nozzle-sern-3d-plan-2.md) | **NO-GO**, C0/M6/m1 | **全件採用**。**M1** → R5 系は生産条件でなく **CPG + 断熱** (`thermalMethod: 0`、壁の `floats:` 空) と確認 → §4.36-1 / R5i。**M2 (実バグ)** → `FLOOR_STUCK` は誤判定。CPG の温度床は `tMin` (1e-4 K) で、50 K の `DEPVAR_TMIN` は TP/凝縮の反転経路だけ。`sern_gates.py::_solver_floors` を経路依存に修正し、`run_0415` は **`GATES: PASS`** になった → §4.35 訂正 1。**M3** → L2 差は符号・位相を落とすので「リミットサイクル」「残り 1 桁が限界」「受理条件の置き換え」を**撤回**、`NOT CONVERGED` を保持 → §4.34 / R5g を未完了へ。**M4** → `dup1` は `i < i_te` のみなので後縁は**同一ノード ID = 正しく閉じた後縁**。「厚さ 0 スリット」の断定を撤回 → §4.35 訂正 2。**M5** → カウル側端テーパが「最後の 2 z-セル」で決まり物理幅が **91 倍**変わる ([[geometry-must-not-follow-mesh-spacing]] の同型違反) → §4.36-2 / R5j。**M6** → §4.32 の「不可能」を「試した分布では不成立」に限定し、側壁解像を受理条件に戻した。**m1** → §5.1 を現状同期 |
| plan | 2026-09-20 | [2026-09-20-tooling-nozzle-sern-3d-plan.md](../../notes/reviews/2026-09-20-tooling-nozzle-sern-3d-plan.md) | GO-with-changes, M6/m1 | **全件採用 → §4.23**。M1 → 非決定性の発生源が `convectiveFlux_slau_d.inc.cuh:563` と `ransTransport_d.cu:55` の float `atomicAdd` と特定 (node も通る)。私の「2 回走に差があるから出力バグでない」という判別論理は不成立と認め、出力の無罪は M2 のソース読解で裏づけ。M2 → 通常出力から EOS 再適用への経路は無いが `dependentVariables` は床適用と `roe` 再構成をするので診断で EOS を呼ばないこと。M3 → `%/step` を廃し**固定ノード集合の ΣρV** に変更 (粗 24 step で 65.4 % → 0、細 12000 step で 99.973 %)。残りの機序分離指標は §5.1。M4/M5/M6 → 最終 step の強制保存・`interp_field` 後の検査・バイナリ固定を §5.1 に起票 |
| **result** | 2026-09-20 | [2026-09-20-tooling-nozzle-sern-3d-result.md](../../notes/reviews/2026-09-20-tooling-nozzle-sern-3d-result.md) | NO-GO, M7/m2 | **全件採用 → §4.20 で大幅撤回**。**M1** → `run_0269` は config もログも `cfl_pseudo=1` で「粗・cfl 5」ではなかった。2×2 は成立せず「CFL 非依存」を撤回。**M2** → `run_0304`/`0305` のログは `implicitRelax=1` で A/B 未実施。「relax で止まる」を撤回。**M3 (実バグ)** → 3D は `_wake_stations(L_cowl, ...)` とカウル後縁から細分しており**ベース直後は 42.75 倍粗いまま**だった。`L_ramp` 分割に修正し、入力値でなく**生成後の実座標**を見る `check_wake_first_spacing()` を追加 (`_geom_start` は公比上限 3 で打ち切るので入力検査では保証できない)。検証: 2D 0.000400 m / 3D 0.004000。**M4 (私の退行)** → ハードエラーが `run_sern_mesh3d_tests.py` を壊していた。未設定は `t_base/5` を既定にして **ALL PASS に復帰**。`t_base/5`・soft ramp・CFL は暫定レシピと明記。**M5/M6/M7** → `run_0308` を「通し完走」までに限定、段階の実効設定と途中場の保存を残作業に、§5.1 と完了条件を同期。**m1/m2** → 現在仕様の旧トポロジ記述と、phase を跨いだ変動係数の計算を訂正。**自分で追加確認**: 同一入力の繰り返しで NaN が step **34/38/53** とばらつく = §4.16–§4.17 の A/B は全部保留 |
| plan | 2026-09-19 | [2026-09-19-tooling-nozzle-sern-chain-plan-2.md](../../notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan-2.md) | NO-GO, C2/M5/m1 | **全件採用** (§4.15.1)。C2 → R4e、M1 → R4d、M3 → R4f を起票 |
| plan | 2026-09-19 | [2026-09-19-tooling-nozzle-sern-chain-plan-3.md](../../notes/reviews/2026-09-19-tooling-nozzle-sern-chain-plan-3.md) | GO-with-changes, C0/M5/m1 | **全件採用** (§4.15.3)。案 (d) 採用、`run_0034` の診断撤回、形状を格子から独立、帳簿三分割、検証 5 段 |
| plan | 2026-09-19 | [2026-09-19-tooling-nozzle-sern-3d-plan.md](../../notes/reviews/2026-09-19-tooling-nozzle-sern-3d-plan.md) | GO-with-changes, C0/M6/m1 | **全件採用** (§4.15.4 訂正 + §4.15.5)。M1 → `P=0` は症状・位置もベース中央と確認し原因認定を撤回、M2 → **リミッタ評価点と再構成点の不一致**を実装で確認し定量 (内部面の 47.7 % が 1 % 以上ズレ)・別ソルバ plan へ、M3/M4/M5 → 切り分け不足を明記 (等温壁 2 次の低 CFL 試験なし / 2D と 3D は別形状 / `nj_wake` は接線方向)、M6 → 受入条件に壁解像・床の時系列を追加、m1 → 数値と実装の同期 |

## 7. 影響範囲

`design/forge_design/meshing/mesh_sern3d.py` / `mesh_sern.py` (2D にも同じ形状モードと格子非依存化)、
`design/forge_design/evaluate/runner_sern3d.py` / `runner_sern.py`、`design/tests/run_sern_mesh3d_tests.py`、
`case/46.sern_design/problem_*.yaml` (`vehicle_clearance` を実効値 0.06 に明示)。

## 8. 完了条件

3D SST が生産 3 作動点で受理ゲートを通り、格子・領域独立性が次の許容値内に入ること:

$$ |\Delta C_T| \le 0.002,\qquad |\Delta C_L| \le 0.002,\qquad |\Delta C_M| \le 0.05 $$

**`C_M` の許容を 0.02 → 0.05 に緩めた (2026-09-21, ユーザ決定。§4.42 / R5p)**。
理由: モーメント基準が $x_{\rm ref} = -20H$ で力の作用点までの腕が 20〜29 H あるため、
$\Delta C_L$ はそのまま $\Delta C_M \approx \Delta C_L \times (d/L_{\rm ref})$ に増幅される。
**$|\Delta C_L| \le 0.002$ を満たしても $|\Delta C_M|$ は 0.04〜0.06 になる**ので、
旧 0.02 は $C_L$ 側の許容と**同時には満たせない数字**だった。
実測 (§4.42) でも $\Delta C_L$ 0.00145 × 29 = 0.042 ≒ $|\Delta C_M|$ 0.0417 と、増幅でほぼ全部説明できる。

**残る不整合 (小)**: 厳密に整合させるなら $0.002 \times 29 = 0.058$ なので **0.06** が要る。
0.05 は $\Delta C_L$ が許容一杯 (0.002) のときにわずかに binding する。
実運用では $\Delta C_L$ は 0.00145 で $\Delta C_M$ 0.0417 < 0.05 なので通る。
**将来 $\Delta C_L$ が 0.002 に近い case が出たら 0.06 へ見直すこと**。

## 9. 変更ログ

- `2026-09-27` — R5o-chi と R5s を合わせて着手 (現行コードで g3/g4 を再取得)。

- `2026-09-22` — **R5r**: CV ごとの閉性検査で、旧メッシャのカウル板の自由な側端に開いた CV (閉性 0.30) を 68 個発見。`dup1` が自由端まで節点を二重化していた。共有に修正 (§4.45)、再現格子で 18 → 0。R5q とは場所が別。CFD での確認は R5s。

- `2026-09-19` — 本体 plan §4.14/§4.15 系と 3D の残作業を切り出し (本体が codex の 128 KB 引数上限を超えたため)。
  内容は移設のみで変更なし。
```

## 出力形式 (この形のまま)

```
結論: <次にやる一手を 1 文で>
第 1 仮説: <内容>  確度: <高/中/低>
  根拠: <ファイル:行 / run パスと数値>
  反証条件: <何が観測されたらこの仮説は誤りか>
第 2・第 3 仮説: <あれば 1 行ずつ>
判別 A/B: <変える設定 1 点、回す長さ、見る量>  → A なら … / B なら …
やらない方がよいこと: <呼び出し側が取りそうな誤った一手>
呼び出し側の前提への異議: <ブリーフの枠組み・除外判断・指標の定義で受け入れなかったものと理由。無ければ「無し」>
不足情報: <あれば>
```
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
