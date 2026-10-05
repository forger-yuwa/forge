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

## ブリーフ (`notes/reviews/briefs/2026-10-05-sern-cowl-blunt-te-route.md`)

# 諮問: SERN 3D のカウル後縁を有限厚にする — 作り方の選択 (2026-10-05、ユーザ決定「有限の厚み」)

関連 plan: `plans/active/tooling-nozzle-sern-chain.md` §5.1 R7b 4-2 (m10_on の後縁低温点と帳簿)、`plans/active/tooling-nozzle-sern-3d.md` (§4.43–4.46、R5q・R5h)、
`plans/active/tooling-sern-mesh-blocking.md` (全ヘキサ gmsh メッシャ、接続模型 PASS・B1d CFD 成立確認と B1e は F で未了、B4–B6 未着手)。
コード: `design/forge_design/meshing/mesh_sern3d.py` (L200–330 カウル板厚の実装)。

## 事実
- m10_on の 3D (exact A・g3) でカウル後縁 (x = L_cowl、厚さ 0 の後縁) の第 2 層節点が T 50 K 床に張り付き評価不成立。後流側 1 面の速度外挿を止めると第 1 層へ冷点が移る (帳簿: 冷却は対流を含む更新、EOS は床到達後の足し戻しだけ)。ユーザ決定: 後縁を有限厚にする (カーネル修正・m10_on 除外は採らない)。
- 現行メッシャ: カウル板は上下 2 列の節点 (双子) を板厚だけ離した「スリット」で、板の中にセルは無い。板厚 `cowl_thickness` (生産 0.005 H) は入口〜0.8 L_cowl で一定、0.8 L_cowl → L_cowl で 0 に絞り、後縁 i_te で上下の節点を共有、下流は 1 列。z 方向は側壁 (z = W/2 = k_sw) で板厚を 0 に絞る (`sz`、最後の 2 セル、側壁スリットとの整合のため。run_0086 で必要と判明)。側壁が終わった後 (x > L_sw) の板の自由な側端は上下共有 (R5r)。
- ランプ側の機体後縁は R4e で有限厚ベース (`t_base` 0.02 H、ベース直後の第一 station ≤ t/5) を実装済み (上バンドの構造の中で `b ≥ yt + t_base` を保つ方式)。
- 全ヘキサ新メッシャ plan: 側壁・カウルを有限厚の固体にし端面を実在の壁にする方針、接続模型 (scale 0.5・0.71) で全ゲート PASS、規模は第一層 16 µm・端面 0.25 mm で 245 万節点 (B1c)、CFD 成立確認 B1d と EOS 床の件 B1e は未了 (F)、SERN 全体への展開 B4・runner/帳簿 B5・CFD 受理の取り直し B6 は未着手。

## 選択肢
- **T1 (現行メッシャを拡張)**: 下バンドと上バンドの間に「ギャップ層」(j 方向の数層) を全 x・全 z に挿入。板の内側 (z ≤ W/2) で x < L_cowl のギャップセルは固体 (除外)、x ≥ L_cowl は流体 (ベース後流)、x = L_cowl のギャップ面がベース壁 (新タグ)。z > W/2 のギャップセルは常に流体。課題: カウル板の側端 (z = W/2) に厚みのある壁面 (側壁がある x ≤ L_sw では側壁との接合、x > L_sw では露出した側端面) が要る = R5q と同じ接合部トポロジ。`sz` テーパ (側端で厚み 0) を残すと、ベース面の高さが z で 0 に潰れ退化セル。
- **T2 (全ヘキサ新メッシャへ)**: tooling-sern-mesh-blocking の B1d → B4–B6 を進め、有限厚カウル (後縁含む) の SERN 全体メッシュを作る。規模・期間大、R7b はその間止まる。
- **T3 (限定)**: 側端は現行どおり厚み 0 のスリットのまま、後縁の厚みを z の内側だけで持たせる (側端近傍の数セルで厚みを 0 に絞る現行 `sz` と同じ扱い) — 退化を避けられるかは未検討。

## 問い
1. T1/T2/T3 のどれ (または他) を選ぶべきか。判断の根拠 (退化・接合部の検証コスト・R7b を止める期間・後の全ヘキサ移行との二重投資)。
2. 有限厚の後縁の物理パラメータ (後縁厚 t_te、板厚分布、ベース後流の格子規則) の決め方。t_te は物理入力か (`t_base` は「暫定モデル値」扱い)。
3. 受入試験の最小設計 (メッシュゲート、m10_on の床張り付き解消、m6_on・m4_off の力係数への影響の測り方と事前閾値)。新形状は基準列 (g3/g4) を取り直す必要があるか。

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
| **R7b** (R7a の後。**順序 (codex 2026-10-05)**: ① R7a の g4 追加診断 → ② ~~現行 C_M の仕様を methods と plan に明記~~ **済 (2026-10-05、`methods/design/overview.md` の「$C_M$ の定義」)** → ③ ~~`L_sw` を物理 station として指定・再現できるようにする~~ **実装 (2026-10-05)**: `mesh3d.L_sw_exact: true` でノズル区間 [0 (角丸めありは xf2), L_cowl] の station を区分線形に写し最寄り station を L_sw に置く (両端固定、間隔比の変化 ~1 %、上限 1.25 を超える・区間端と重なる指定は生成失敗)。既定 false は従来どおり丸め (既存 run の形は不変)。試験 `run_sern_mesh3d_tests` に 6 項目 (指定値どおり・0.8 と 0.8001 が別形状・閉性・Jacobian・区間端の拒否) — ALL PASS。生産格子での品質確認 (2026-10-05、AWS で prepare のみ): `L_sw` 0.8 + `L_sw_exact` で g3 は `PASS` (skew max 0.449)、g4 は `SOFT-PASS` (skew max 0.701) — どちらも元の格子と同じ判定・同じ skew、側壁後縁 x = 0.0800 m (指定どおり)。**③ 済** → ④ その形状で他作動点の評価成立と設計差の格子・領域感度 (**判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-r7b-step4-diagnose.md) — 全件採用**: 2D の再診断を前提にしない [`chi_evidence` の run_0968 m10_on・run_0966 m4_off は 2D で 4 量 STEADY・残差プラトー、ただし viscMethod 1・旧種 DB]; `L_sw_exact` の形状は力係数の不変を保証しないので exact 形状で m6_on の二水準を取り直す; 他作動点は 3 作動点とも A/B を対にして格子・領域感度; 0.002 は検出目標で各作動点の必須改善量ではない (加重目的 ΔJ = 0.5Δ₆ + 0.3Δ₁₀ + 0.2Δ₄、E_J = 0.5E₆ + 0.3E₁₀ + 0.2E₄)。**順序**: 4-1 exact 形状・m6_on・g3 の 0.8/1.0 対 → 4-2 m10_on → m4_off の g3・A 形状の起動成立 → 4-3 両作動点の B → 4-4 3 作動点 × A/B の g4 → 4-5 3 作動点 × A/B の領域感度 (一方向ずつ + 同時拡大、共通領域を保持して外側にセルを足す) → ⑤。**4-1 の事前登録 (2026-10-05、run 前)**: A = `run_1050_r7b_x_lsw08` (`L_sw` 0.8 + `L_sw_exact`)、B = `run_1051_r7b_x_lsw10` (1.0 + exact)、他は R7a と同じ (m6_on・g3・生産 YAML の複製・時間/space は元 run を継承 = リミッタ基準値固定・eb0ea331・blocksize 128)。初期場: A は `run_1034` (A・g3 の発達場)、B は `run_1040` (B・g3) から**節点番号で写す** (`r7b_index_restart.py`: 節点数・各節点の境界タグの組が一致し、座標のずれがノズル区間の x だけで半間隔以下であることを検査してから保存量を index コピー、wall_dist は新格子の値)。各 1 本、判定区間 20000 step・500 step 出力、窓条件未達なら事前登録の +20000 を 1 回、なお未達は判定不能。q = C_T_with_shear、U = Σ max(a, d, 1e-6)。**結果 A** (ゲート成立かつ |Δ| − U > 0.002) → exact 形状への持ち越しを支持 (各 1 本の予備判別。反復込みの認定には各 3 反復で取り直す) → 4-2 へ。**結果 B** (|Δ| + U ≤ 0.002) → 持ち越しを棄却。他は判定不能。**4-1 の結果 (2026-10-05): 結果 A** — `run_1050` (A、後縁 0.0800 m)・`run_1051` (B、0.1000 m) とも 20000 step で窓条件 OK・4 量 STEADY・GATES PASS。C_T_with_shear A 0.8704579 / B 0.8725856、Δ +0.002128、U 6.2e-6、|Δ| − U 0.002122 > 0.002 → exact 形状への持ち越しを支持 (各 1 本の予備判別)。他: C_T Δ +0.003514、C_L +0.009323、C_M −0.2209。丸め版との差 (exact − 丸め) は A で C_T_with_shear −9.8e-5・C_L −5.3e-4・C_M +0.012、B で −5.6e-5・−1.9e-4・+0.0045 — 側壁後縁が A で 0.91 mm・B で 0.43 mm 短くなった形状差を含む (A→B の傾き dC_L/dL_sw ≈ 0.47 /m からの見積もり A −4.3e-4・B −2.0e-4 と整合、station 写像自体の寄与は切り分けていない)。原本 `notes/investigations/2026-10-05-sern-r7a/R7B_X_VERDICT.txt`。次は 4-2。**4-2 の事前登録 (2026-10-05、run 前)**: 問題 = `problem_3d_prod_3op_wallres_lswx08.yaml` (exact A 形状の生産 3D に、2D MOO YAML と同じ m10_on・m4_off の定義を追加、`evaluate.limiter_ref` を run_0986 の値に固定、輸送は作動点の実種に絞られる)。m10_on → m4_off の順に各 1 本、g3。起動 = runner_sern3d の既定の段階起動 (生産 YAML の層流暖機 2000 [CFL 0.2] → SST 1 次 2000 → 2 次 2000 → 本段 4000 [CFL 0.25])、MOC の領域別 IC から (m6_on の場からの warm start はしない)。その後、本段と同一設定で restart して 20000 step を**判定区間**とし (段階起動の区間は連結しない)、窓条件未達なら +20000 を 1 回。run 名: m10_on `run_1054_r7b_m10_A` (段階起動) → `run_1055_r7b_m10_A_c` (判定区間) [→ `run_1056` 延長]、m4_off `run_1057_r7b_m4_A` → `run_1058_r7b_m4_A_c` [→ `run_1059`]。**起動成立の条件**: 全段が NaN なく完走 (rc 0、detectNaN 無し)。**評価成立の条件 (判定区間)**: GATES PASS・4 量 STEADY・窓条件 (C_T・C_T_with_shear・C_L ≤ 5e-5、C_M ≤ 5e-4)・床/ω 下限への張り付き 0・化学種の有限・非負・組成和、`check_convergence` は区間付きで併記 (プラトーは失格にしない)。壁解像は `check_wall_resolution.py` で局所 y₁⁺ を記録 (m6_on の緩和を自動適用しない)。成立しなければ症状を記録して止め、諮る (2D に戻るかを含む)。**4-2 の結果 (2026-10-05)**: 段階起動は両作動点とも完走 (NaN なし)。**m4_off** (`run_1057` → 判定区間 `run_1058`): GATES PASS・4 量 STEADY・窓条件 OK (前窓差 ≤ 2.6e-4) → **評価成立** (C_T_with_shear 0.9676590、C_L −0.02494、C_M +0.5065)。壁解像 (局所 y₁⁺) は未計測。**m10_on** (`run_1054` → `run_1055`): 4 量 STEADY・窓条件 OK (C_T_with_shear 0.8970375、C_M −1.7785) だが **GATES FAIL `FLOOR_STUCK` (T ≤ 50 K が 1 節点)** → 事前規則で**評価不成立**。その節点は (x 0.12000, y −0.01053, z 0.0979) = カウル後縁 (厚さ 0) で、R5h の低温点と同じ場所 (m6_on では 95–105 K)。2 番目に低い節点は 211 K で、孤立した 1 点。→ R5h の再開条件 (「新形状・新作動点で床到達が増えたら再開」、3D plan §5.1 R5h) に該当。次の手は諮る (エスカレーション条件 2/3)。**判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-r7b-m10-te-floor-diagnose.md) — 全件採用**: R5h を再開 (恒久修正は保留)、有限厚後縁の先行は却下 (形状・帳簿が同時に変わる)、**床ゲートからの除外は却下** (TP は温度反転後に roe を再構成するので床到達は表示だけの問題でない)、旧 m6_on の感度で無影響とするのは却下、**m10_on 抜きで 4-3 以降を進めるのは却下** (先行は m4_off の壁解像確認だけ)。観測と解釈の分離: 「50 K の 1 点」は観測、「厚さ 0 の壁ノード」「旧機序と同一」「積分力に効かない」は未確認 (旧低温点は後縁直下の**内部節点**と記録)、`FLOOR_STUCK` は最終場の検査で判定区間全体の張り付きは別に見る。**R5h-m10 診断 A/B の事前登録 (2026-10-05、run 前)**: 起点 `run_1055_r7b_m10_A_c` 最終場を restart_field で 2 run に分岐、変更は `FORGE_DIAG_FACE_VEL_CELL` の有無だけ (A = 対照、B = 節点 517160 に接する後流側 1 面の同節点側速度 3 成分だけセル値へ; 面 ID は今回の接続から帳簿ダンプで特定、旧 1535830 は流用しない)。(1) 機序: 最初の共通状態の 1 更新を帳簿で記録 — 対象点・相手点・隣接 CV の面流束和・対流/粘性残差・実保存量増分・EOS 前後の変更を分け、内部エネルギー射影 R_ρE − u·R_ρu + (|u|²/2 − e_int)R_ρ と実更新を別に確認、閉合誤差 ≤ 介入差の 1 % または丸め上限の大きい方。(2) 100 step で早期確認し、機序を支持するときだけ各 20000 step まで継続 (窓条件未達なら +20000 を 1 回)。床離脱 = B の末尾 10000 step の全保存場で温度床到達 0・他の床/化学種ゲート合格、近傍最低温度も追う (冷点の移動を除く)。(3) 力の感度: 両者 4 量 STEADY・窓条件、Dq = |平均_B − 平均_A| + max(a_A, d_A, 1e-6) + max(a_B, d_B, 1e-6) ≤ C_T・C_T_with_shear・C_L 5e-5、C_M 5e-4 (領域許容の 1/10、この介入への感度が小さいとだけ認定)。**分岐**: 床離脱・加熱側の収支・閉合が成立 → 第 1 仮説 (後流向き面の速度外挿が冷却を維持) を支持し恒久対策の別 plan へ。有効な介入でも加熱応答なし・準定常後も床が残る → 棄却。介入無効・閉合不足・過渡未終了は判定不能。A の FLOOR_STUCK は対照として保持し生産評価に受理しない。設計差への無影響はこの 1 組では証明できない (δ₁.₀ − δ₀.₈ が要る)。やらない: 床ゲート除外・床値変更・`bndFirstOrder`・診断の面/節点指定を生産へ入れる・100 step の昇温で影響なしと結論。**m4_off の壁解像 (2026-10-05、`check_wall_resolution.py`)**: run_1058 で FAIL — 側面の壁 (sidewall・vehicle_side、もともと非解像 = R5m) が y₁⁺ 平均 19–35 で主因。y 法線壁は ramp 0.55 (>1 が 6.9 %)・cowl_in 1.03 (49 %)・cowl_out 0.82 (6.1 %)。**比較: 同じ exact A 形状の m6_on (run_1050) も FAIL** で y 法線壁は ramp 2.05 (57 %)・cowl_in 4.69 (97 %) — m4_off のほうが良い。生産 g3 格子は m6_on でも局所 y₁⁺ ≤ 1 を満たしていない (§4.42 の ramp 0.549・cowl_in 1.234 より悪い; 記録の食い違いは未調査)。**R5h-m10 診断 A/B の結果 (2026-10-05)**: 面 ID = 1535833 (517160 の後流側 +x 面、相手 522113 x 0.12255、A の F_roe +11.0 で 6 面中最大)。(1) 1 更新 (`run_1061` A / `run_1062` B、帳簿): 517160 の res_final roe A −0.662 / B +10.53、内部エネルギー射影 A −0.801 / B +10.745 (加熱側へ)、差は対流段で生じ、相手 522113 の差がちょうど逆符号 (roe B − A +11.19 / −11.19 で閉じる)。1 step 後の T[517160] は A・B とも 50.00002 K (1 更新では温度に出ない)。(2) 100 step (`run_1063` A / `run_1064` B): A は T[517160] 50 K のまま。**B は T[517160] 384.6 K に上がったが、隣の 517199 (y −0.01051、同じ x・z、面 1535832 の相手) が 413 → 50 K に落ち、T ≤ 50 は 1 節点のまま** (517200 も 205 → 115 K)。→ 介入した 1 面の速度外挿が 517160 の冷却を担うことは支持されるが、床からは離脱せず冷点が隣へ移った (事前に除くとした「冷点の移動」)。事前登録の「機序を支持するときだけ 20000 step へ」の扱いは諮問 (中間・判定不能の可能性)。**判断: 2026-10-05 codex (diagnose) [記録](../../notes/reviews/2026-10-05-sern-r5h-m10-ab-diagnose.md) — 全件採用 (訂正 3 件含む)**: (1) 事前登録の総合判定は**判定不能** (100 step は過渡未終了、局所残差への効果だけ採用)。(2) 「1 step 後 T 同じ」を無応答と読むのは誤り — 順序は EOS → 残差 → commit → 出力で、出力の T は更新前の値。call k の残差は call k+1 の entry 保存量と結ぶ。(3) **517199 を「板の反対側の双子」と読んだのは誤り**: 番号式 (i·NJ + j)·nz + k (NJ 127、nz 39、カウル線 j 54) で 517160 = (104, 52, 20)、517199 = (104, 53, 20)、後縁壁点 517238 = (104, 54, 20) → 両点ともカウル下側で、冷点は第 2 層から第 1 層 (壁隣) へ移った。(4) **壁解像の「悪化」は g3/g4 の取り違え**: §4.42 の ramp 0.549・cowl_in 1.234 は g4、g3 は 2.167・4.919 で今回 g3 の 2.05・4.69 と近い → R7b ④ の順序は変えない (局所解像不足は保持)。(5) `check_wall_resolution.py` の「面積割合」は実装上は**個数割合** (面積重みなし) → 上の 57 %・97 % 等は個数割合と読む (ツールの修正は残作業)。**次 (事前登録 2026-10-05、run 前)**: 同じ起点・同じ 1 面 A/B を各 100 step、帳簿を毎 call (`FORGE_DUMP_LEDGER_CALLS=100`、節点 517160・517199・517200・522113・517238) で取り、call k の残差 → call k+1 の entry (commit 後の保存量) → call k+1 の after_eos_bc (EOS・境界の上書き) を対応づける。517199 が初めて床へ届く call を特定し、その直前の更新が対流を含む残差で説明できるか (結果 A: 近壁層間の連成を支持) / EOS・境界の上書きで生じるか (結果 B: 対流再構成だけの説明を棄却) を見る。閉合許容は介入差の 1 % または丸め上限。20000 step 継続・2 面介入・有限厚後縁への直行はしない。**結果 (2026-10-05、`run_1065` A / `run_1066` B、帳簿 毎 call 100 回、原本 `notes/investigations/2026-10-05-sern-r7a/R5H_M10_LEDGER100.txt`): 結果 A** — A: 517160 は全 call で床、517199 は 413 → 413 K、517200 218 → 205 K。B: 517160 は 50 → 384.6 K、**517199 は call 92 で初めて床**。517199 の call ごとの内訳: call 1–5 は更新で加熱 (413 → 421 K、内部エネルギー射影 +0.02〜+0.15)、call 6 以降は射影が負に転じ更新で単調に冷却 (Δe_int の更新分 −1.7e3〜−7.8e3 J/kg/step、EOS の上書きは ≤ 1e-2 で実質 0) → call 91 で 50.7 → 50.0 K。床到達後は更新が負のエネルギーを入れ、EOS がそれと同量を足し戻して T を 50 K に保つ (打ち消す量は −3.8e3 → −3.4e4 J/kg/step と毎 step 増える)。→ 517199 の冷却は EOS 前の更新 (対流を含む残差) で説明でき、EOS の上書きで冷点が生まれる説明は棄却。**近壁層間 (第 2 層 517160 ↔ 第 1 層 517199) の対流の連成を支持**。床到達後は EOS が毎 step エネルギーを足す床として働く (保存量を変える; 床ゲート除外が不可である根拠を補強)。**次の手 (ユーザ決定 2026-10-05「有限の厚みですね」)**: カウル後縁を**有限厚** (鈍い後縁 + ベース面) にする。再構成のカーネル修正・m10_on を除いて進める案は採らない。作り方 (現行テンソル積メッシャの拡張か、全ヘキサ新メッシャ [tooling-sern-mesh-blocking] か) は諮問中: 有限厚にするとカウル板の側端 (z = W/2) で板厚を 0 に絞る現行処理と両立せず、R5q の接合部トポロジと同じ問題になる → ⑤ 3D 最適化へ接続) | **3D 最適化への接続と MOO 再取得** | 3D runner と `L_sw` を最適化 driver に接続 (`L_sw` はメッシャが既存 station に丸めるので、連続変数にする前に物理 station 化)。対象作動点の評価成立・設計差の格子/領域感度・作動点別 C_M 窓 (実機の値と出典はユーザ判断、生産 YAML は加重平均 `cm_min` −7.0 のみで `cm_window` 未指定) を確かめてから MOO 再取得。2D は形状生成・傾向把握に使う (MOO を 2D で回して 3D は代表点確認、という分担は 3D SST をループに入れるユーザ判断 [本表 L442] と矛盾するので採らない) | — | F |
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

## 参考: `plans/active/tooling-sern-mesh-blocking.md`

```
# ⑤ SERN 3D メッシュの CAD ベース化 — **非構造 (Salome) を先行**、全ヘキサ・ブロッキングは形状確定後

⑤ SERN の 3D メッシュ生成を、自作テンソル積メッシャ (`mesh_sern3d.py`) から
**gmsh Python API による全ヘキサ・マルチブロック** (case/49 と同じ方式) へ移す。
動機は 2 つ: (1) 厚さ 0 スリットの交線 (ρ が床に張り付く欠陥の現場。格子細分で悪化) を**有限厚の形状で除去し、床張り付きの解消を検証できる形状を作る** (根治したかどうかの判定は後続 B6 の CFD に委ねる。codex plan m1)、
(2) ユーザ要件の**断面隅フィレット** (側壁×ランプ・側壁×カウル) はテンソル積では作れない。

## メタ

- **area**: `tooling / optimization`
- **status**: `draft`
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) の「メッシュ」節 (SERN 3D ブロックメッシャの項を追加済み)
- **related_plans**:
  - 親: [`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (§4.43 欠陥の特定 / §4.44 撤回した設計 / R5q の移管元)
  - 先例: case/49 の全ヘキサ・ブロッキング (`case/49.plate_annular_cavity_m5/cad/build_hex_mesh.py`。**別セッションの資産なので読むだけ**)
- **created**: `2026-09-21`
- **owner**: `CFD Dev (Claude セッション: SERN 3D)`

## 1. 目的

1. SERN 3D の形状から**厚さ 0 の面を無くす** (側壁・カウルを有限厚の固体にし、端面を実在の壁にする)。
2. **断面 (y–z) の隅フィレット**を O グリッドで受ける。
3. 受理に必要な性質を保つ: **系統的な格子細分列**・**冷却壁 y⁺≈1 (第一層 4 µm・約 40 層)**・
   **dv を振っても無人で生成**・**位相固定で寸法と輪郭だけ可変**・msh4.1 → `convertGmshToForge` の既存 node 経路。

## 2. スコープ

- **含む (本 plan)**: 新バックエンドの設計、**方式確定のための実証模型「側壁終端を含む 3D 接続模型」** (メッシュのみ・CFD なし)、
  その受入試験、gmsh 実装方針。
- **含まない (方式確定後に別フェーズ)**: SERN 全体 (外部流ブロック全域・機体上面/ベース・プルーム全長) への展開、
  runner 統合、力の帳簿の改修、CFD の格子収束列の再取得。§5.1 に後続フェーズとして置く。
- **旧バックエンドは消さない**: `mesh_sern3d.py` は互換用に保持し、既定のままにする (親 plan codex M5)。

## 3. 関連 docs と前提

- 方式選定の経緯: [`notes/reviews/2026-09-21-codex-mesh-approach.md`](../../notes/reviews/2026-09-21-codex-mesh-approach.md) (フィレット無し前提 → 自作)、
  [`…-approach-2.md`](../../notes/reviews/2026-09-21-codex-mesh-approach-2.md) (隅フィレット要件 + 前提訂正 → **本方式**)。
- 撤回した設計と NO-GO の Major 7: [`notes/reviews/2026-09-21-tooling-nozzle-sern-3d-plan.md`](../../notes/reviews/2026-09-21-tooling-nozzle-sern-3d-plan.md)。
- Salome ViscousLayers (tet+prism) は採らない: case/49 の実測で **VL 総厚と接線サイズが結合し系統的な細分列が作れなかった**。
- 品質ゲート: AR ≤ 5000 (壁法線の構造層) / skew ≤ 0.90 (`check_mesh_quality.py`)。形状は格子間隔から決めない
  ([[geometry-must-not-follow-mesh-spacing]])。壁が終わった下流に壁用の細層を敷かない ([[mesh-wall-layer-x-blend]])。

## 4. 設計方針

### 4.1 形状モデル — 固体を明示し、流体はその補集合

半スパン模型 (z = 0 対称)。長さは H で無次元化。**形状パラメータはすべて物理入力**で、格子パラメータと独立。

| 固体 | 範囲 | パラメータ |
| --- | --- | --- |
| ランプ (上壁) | 輪郭 $y_r(x)$ (逆設計 MOC の点列をスプラインで保持)、全 z | 設計から |
| カウル板 | $0 \le x \le L_{cowl}$、$y \in [y_c(x) - t_c,\; y_c(x)]$、$z \in [0,\; z_{out}]$ | `cowl_thickness` $t_c$ |
| 側壁板 | $0 \le x \le L_{sw}$、$z \in [W/2,\; z_{out}]$、$y$ はカウル上面〜ランプ | **`sidewall_thickness` $t_{sw}$ (新規)**、$z_{out} = W/2 + t_{sw}$ |
| 隅フィレット | 側壁内面×ランプ、側壁内面×カウル上面。半径 $r_f(x)$ | **`corner_fillet` $r_f$ (新規)**、終端テーパ長 `fillet_taper_len` |

**側壁の厚みは外向きに付ける** ($W/2 \to W/2 + t_{sw}$)。**ノズル幅 W は変わらない**。~~入口面積も変わらないので親 plan M6 を形状側で回避できる~~ **← 誤り (codex plan M4)**: 入口から存在する内隅フィレットが流路面積を減らす (半スパン断面で $\Delta A = 2(1-\pi/4) r_f^2$、$r_f = 0.10H$ で 0.00429 H²)。現行の入口項は矩形面積 (`half_W` 倍) のままなので、**入口帳簿の実面積化 (B5) は必須**のまま。

**端面はすべて実在の壁**: 側壁後端面 `sidewall_end` (x = $L_{sw}$)、露出カウル側端 `cowl_side`
($L_{sw} < x \le L_{cowl}$ の z = $z_{out}$ 面)、カウル後端面 `cowl_base` (x = $L_{cowl}$。現行の「後縁で厚さ 0 に絞る」をやめ
有限厚ベースにする = 親 plan §4.35 の低温スポットの形状対策も兼ねる)。

**フィレットの終端**: $r_f(x)$ を $[L_{sw} - L_{ft},\; L_{sw}]$ で **0 へ滑らかに絞る** (smoothstep、$L_{ft}$ は物理長)。
鈍頭のまま終えると、後端面の後ろの後流ブロックの断面が「円弧がランプに**接する**尖点」を持ち、内角 0° でヘキサが切れないため。
(→ §未確定事項 U1: 終端形状はユーザ確認)

### 4.2 断面 (y–z) のブロッキング

**ダクト内 (N 領域)**: 3 面が壁 (上 = ランプ、右 = 側壁内面、下 = カウル上面)、左 = 対称面。
**バタフライ型の C リング**: コア 1 + リング 3 (上・側・下)。コアの右上/右下の角から**フィレット円弧の中点**へ対角線を引いて
リングを分割する。リング厚 $\delta_r$ は**物理入力** (境界層を収める厚さ。既定 0.08 H)。
壁法線方向 = リングを横切る方向で、ここに第一層 $h_1$ からの幾何級数を置く。
$r_f \to 0$ でも位相は変わらない (円弧半分 + 直線で 1 辺なので辺長は正のまま。角部は 45°/135° のセルになり skew ≈ 0.5)。

**外部 (U = カウル下、S = 側壁の外)**: 隅は尖ったままなので **H 型の直交ブロック**。z 線 {0, $z_{out}$, $z_{out}+\delta_r$, $Z_{far}$}、
y 線 {$y_{bot}$, カウル下面 $-\delta_r$, カウル下面, カウル上面, ランプ}。壁に接する帯が壁層、帯は壁の外へ**内部帯としてそのまま延びる**。
S の上は機体下面 (ランプ面の z 方向延長、壁 `vehicle`)。

**共有面の整合**: 側壁内面上では N のリング (側) が面全高を受け持つ (対角線は $r_f=0$ で角に落ちる)。
側壁より下流ではこの面が内部面になり、スロット後流ブロック・S の帯と**同じ y 分割**を共有する。
分割数は「辺の同値類」ごとに 1 つの変数で持ち、全ブロックで一致させる (case/49 の実装と同じ作法)。

### 4.3 x 方向 — 物理端点で区間を切り、位相は全区間で固定

区間: $[0,\; L_{sw}-L_{ft}]$ (フィレット一定) / $[L_{sw}-L_{ft},\; L_{sw}]$ (フィレット絞り) /
$[L_{sw},\; L_{cowl}]$ (露出カウル側端) / $[L_{cowl},\; x_{end}]$ (短い後流)。
**端点は物理入力そのもの**で、最近傍 station への丸めをしない (親 plan M4)。

**1 断面の平行移動 `extrude` は使わない**。各区間の両端に y–z 断面のエンティティを作り、対応する頂点を
**x 方向の輪郭曲線** (ランプはスプライン、その他は直線/テーパ) で結び、曲面パッチ + transfinite 体積で六面体ブロックにする。

**壁が終わっても格子帯を潰さない**:
- 側壁の後ろ: 固体だったスロット ($z \in [W/2, z_{out}]$) を**スロット後流ブロック**で埋める。上流面 = `sidewall_end` (壁)。
  N のリング (側) の外面と S の帯の内面は、ここから内部面としてスロット後流ブロックと共有される。
- カウルの後ろ: 同様に**カウル後流ブロック** (上流面 = `cowl_base`)。N と U はこれを介してつながる。
- $L_{sw} < x \le L_{cowl}$: カウル板は続くので、その z = $z_{out}$ 面が `cowl_side` (壁) になる。

**固定位相の成立条件を事前判定する**: $0 < L_{ft} < L_{sw} < L_{cowl}$、$2 r_f + 2\delta_r <$ ダクト高さ/幅の最小値、
$t_{sw}, t_c > 0$、区間長 ≥ 下限。満たさない dv は**生成前にエラー** (黙って別位相にしない)。

### 4.4 分割と細分列

- 壁法線: リング厚 $\delta_r$ 固定のもとで (層数 $n$, 第一層 $h_1$) を与え、成長率は $(\delta_r, n, h_1)$ から解く。
  生産目標は $h_1 = 4\times10^{-5}$ H (= 4 µm)・$n \approx 40$・成長率 ≈ 1.15。
- 細分列は **`scale`** で全方向の分割数を一律倍・$h_1$ を 1/倍 (case/49 と同じ「真の細分列」)。
  実証では $h_1$ を 4 倍刻みの 3 水準 (6.4e-4 / 1.6e-4 / 4e-5 H) × 層数 (20 / 30 / 40)、接線分割は ×1 / ×1.5 / ×2。
- 形状はブロック境界の幾何エンティティが持つ。**解像度を変えても節点は同じ曲線・曲面の上に乗る** = 形状は格子に依らない。
- 壁が終わった下流では、帯の第一層を物理長で粗い値へ移す (x ブレンド。判定量は df/dx ≤ 4e-3)。

### 4.5 タグと帳簿

既存 physID の意味は保つ (`ramp` 4 / `cowl_in` 5 / `cowl_out` 6 / `sidewall_in` 11 / `sidewall_out` 12 / `vehicle` 14 ほか)。
新設: `sidewall_end` / `cowl_side` / `cowl_base`。~~フィレット面は `sidewall_in` に含める~~ → **フィレット面は 45° 点で二分し、隣接する主壁に帰属させる** (`e1`×x = `ramp` は上フィレットのランプ側半分を含む / `e5`×x = `sidewall_in` は上下フィレットの側壁側半分を含む / `e8`×x = `cowl_in` は下フィレットのカウル側半分を含む。2026-09-22 codex M4: 1 本の壁曲線を掃引した面に辺名でタグを付ける方式では、フィレットだけを別の壁に帰属できない)。
**力への帰属** (ノズル力か機体力か) は後続フェーズで決めるが、実証の段階で「全外部面がちょうど 1 つのタグに属する」ことは検査する。

### 4.6 gmsh 実装方針

`gmsh.model.geo` (組込みカーネル): 点・直線・`addCircleArc`・`addSpline` → 4 辺 (または 3 辺) の曲面パッチ
(`addSurfaceFilling`) → 6 面の体積 → `setTransfiniteCurve/Surface/Volume` + `setRecombine`。msh4.1 で書き出し既存変換器へ。
コードは `design/forge_design/meshing/mesh_sern3d_blocks.py` (新規)。gmsh は pip パッケージ (AWS でも導入可)。
**未実証の核心**: transfinite 体積が「x で変わる断面 + 曲面の側面 + 第一層 4 µm × 40 層」で品質ゲートを通るか。
これを確かめるのが実証模型で、通らなければ方式を再検討する。

### 4.7 実現性スパイクで分かったこと (2026-09-21, scratchpad・未 commit)

gmsh 4.15.2 (`.venv-mesh`) の `geo` + transfinite で、バタフライの上リングブロック 1 個を x 方向 2 区間つないで試した
(ランプは曲線、壁法線は第一層 4 µm × 40 層、成長率は**辺ごとに解く**)。

| 構成 | 結果 |
| --- | --- |
| $r_f$ 一定 (0.10 H) | **100 % ヘキサ (30720)・負 Jacobian 0 (float32 後も 0)・最短辺 4.000e-06 (指定どおり)・AR max 768**。ただし **skew max 0.861 / 平均 0.558**、壁→第一内部ノード距離が **4.0e-06〜2.0e-05 (中央値 8.9e-06)** とばらつく |
| $r_f$ 0.10 → 0.02 H に絞る | **負 Jacobian 72・skew max 1.000 (>0.90 が 277)** — 折れる |
| $r_f$ 0.10 → **0** に絞る | **負 Jacobian 782・skew max 1.000** — 折れる。円弧専用ブロックに分ける案は $r_f < \delta_r$ で裏返る |

**読み**: (1) 「x で変わる断面 + 曲面 + 4 µm × 40 層」の単純ブロックは**切れる** (codex の独立試験でも全ヘキサ・最小 Jacobian 3.04e-08 > 0)。
(2) **フィレットを絞る区間は素朴な作りでは折れる** — 壁断面曲線のパラメータ化が station ごとに違うのが原因とみられる。
(3) **線形 TFI は第一層厚をブロックの局所厚みに比例させてしまう** — 端の辺で $h_1$ を指定するだけでは壁全体で保てない。
リングを**一様厚み** (内側の辺を壁曲線のオフセットにする) にし、壁法線方向の辺を**3D の壁法線**に合わせる必要がある。
いずれも codex plan レビューの M1 / M3 と同じ結論。

**追試 (2026-09-22): 絞り区間の折れは壁曲線のパラメータ化が原因だった**。壁断面曲線 (対称面 → ランプ直線 → 円弧) を
**弧長一様に標本化した 1 本のスプライン**にして全 station で揃えると、同じ絞り ($r_f$ 0.10 H → **0**) で:

| | 点列がばらばら (旧) | **弧長一様の 1 本スプライン** |
| --- | --- | --- |
| 負 Jacobian | 782 | **0 (float32 後も 0)** |
| skew max / 平均 | 1.000 / 0.50 | **0.501 / 0.277** (0.501 は角の 45° セル = 設計どおり) |
| AR max | 2829 | 822 |
| 壁→第一内部ノード (中央値) | 5.6e-06 | **4.06e-06** (指定 4.0e-06) |

絞り区間を半分 (0.3 H) にしても同じ。**これで codex M1 (直線 + 円弧を 1 辺にできない・半径 0 の円弧は作れない) を回避できる**:
円弧エンティティを持たず、1 本のスプラインが $r_f \to 0$ で直線に縮退するだけなので位相は変わらない。
**形状公差** ($r_f$ = 10 mm の四分円 + 直線、最大ずれ): 標本 33 点 81.7 µm / 65 点 20.3 µm / 129 点 5.1 µm / **257 点 1.3 µm (0.013 % of $r_f$)**。
標本点数の 2 乗で減り、制御点は形状側で固定するので**格子解像度に依らない一定の公差**。257 点を既定にする想定 (→ U3)。
未解決: 壁→第一内部ノード距離の**最小値 2.9e-07** が 1 点出ている (中央値・最大は正常)。B0c の点別検査で原因を特定する。

**絞り区間はどこまで短くできるか (2026-09-22)**: $r_f$ = 0.10 H を 0 へ絞る長さを振った (同じ上リングブロック、4 µm × 40 層)。

| 絞り長さ | $r_f$ 比 | 負 Jacobian | skew max |
| --- | --- | --- | --- |
| 0.30 H | 3 倍 | 0 | 0.501 |
| 0.10 H | 1 倍 | 0 | 0.501 |
| **0.05 H** | **0.5 倍** | **0** | **0.506** |

**$r_f$ の半分の長さで絞っても崩れない**。「終端までフィレットを付けたままの鈍頭」は全ヘキサでは尖点 (内角 0°) のため品質ゲート内では切れないが、**終端の手前 0.5–1 $r_f$ で消す形**なら実質的にそれに近い形状を全ヘキサで持てる。
鈍頭との差が問題になるなら、鈍頭をそのまま切れる Salome 非構造 (§4.9) と同条件で回して推力係数への影響を実測する。

**非構造の節点数は壁の接線サイズで決まる** (§4.9 の補足): 接続模型で 1.5 mm → 1,381,415 節点、**3.0 mm → 395,009 節点** (tet 35 万 + prism 64 万)。SERN 全体で 3 mm なら 300–400 万節点級で全ヘキサ (322 万) と同じ桁になるが、流れ方向の解像は落ちる (全ヘキサの利点は流れ方向に長い異方性セルを置けること)。

**全ヘキサの試作形状 (2026-09-22, `case/46.sern_design/cad/hex_junction_demo.py`)**: ダクト半スパン断面を バタフライ C リング 3 (上・側・下) + コア 1 で切り、x 方向 5 区間 = **20 ブロック**を transfinite 体積でつないだ。上隅は入口で角 → 15 mm で r 10 mm へ徐変 → 終端の手前 10 mm (= 1 $r_f$) で消える。下隅は r 6 mm 一定 → 同じ位置で消える。ランプは曲線。
**321,786 節点 / 307,200 要素 (100 % ヘキサ)・負 Jacobian 0 (float32 後も 0)・skew max 0.501 / 平均 0.207 (>0.90 は 0)・AR max 1177 (>5000 は 0)**。最短辺 2.83e-06 は 45° の対角線上で第一層 4e-06 が $\cos 45°$ 倍に投影されたもの。
図は Artifact (全ヘキサと Salome 非構造の比較ページ) でユーザに提示。出力は `cad/demo_out/` (git 追跡外)。

**鈍頭終端の尖点はいつ出るか (2026-09-22, ユーザ指摘で整理)**: 尖点 (内角 0°) ができるのは、フィレットが接する 2 つの壁のうち**片方だけが終わる**場合 (側壁が終わりランプが下流へ続く)。**接する 2 壁が同じ位置で厚みを持って一緒に終わるなら尖点は消える** (ガセットの後端面が相手のベース面と同一平面でつながり、接点は輪郭の滑らかな点 = 内角 180°)。
SERN では: **下隅 (カウル×側壁) は `L_sw = L_cowl` にして両方を有限厚ベースにすれば鈍頭で終えられる**。上隅 (ランプ×側壁) は側壁をランプ後端まで延ばす場合だけで、途中で終わるなら手前で絞る。(旧メッシュで `L_sw = L_cowl` が退化とされたのは板が厚さ 0 だったためで、有限厚なら問題にならない。)

選択肢の整理 (案 1 フェードアウト / 案 2 鈍頭は尖点で全ヘキサ不成立) は Artifact にまとめ、ユーザに 4 点の決定を依頼中:
終端の形・絞り長さ $L_{ft}$・入口側の始め方・$r_f$ の扱い (固定か dv か、上下同値か)。

### 4.8 **ユーザ決定 (2026-09-22): 形状が未確定なので、非構造で切れるようにする。ツールは Salome**

> 「正直どういう形状にするかなんて決まってないのだから、非構造できれるようにしといてくれないか。
> 出発点が徐変で、終端はフィレットつけたままの可能性もあるし」
> 「下側についているかもしれないし、上側についているかもしれない」「salome の方がよさそうに思います」

**受け止め**: 全ヘキサ・ブロッキング (§4.1–§4.7) は形状ごとに位相を手で設計する方式で、
「終端にフィレットを付けたままの鈍頭」は尖点のため**全ヘキサでは切れない** (§4.7 追試の案 2)。
形状が固まる前にブロッキングを先行させるのは順序が逆だった。**優先を入れ替える**:

1. **先行: CAD (Salome GEOM) + 非構造 (SMESH: NETGEN tet + ViscousLayers prism)**。形状の自由度を最優先。
   隅フィレットは**隅ごとに独立な半径プロファイル r(x)** (上 = ランプ×側壁、下 = カウル×側壁。
   始端の徐変・途中一定・終端は「付けたまま」か「消す」を選べる) を入力にする。
2. **保留: 全ヘキサ・ブロッキング (B0–B3)**。形状が確定し、系統的な格子細分列が必要になった段階で再開する。
   §4.7 の知見 (弧長一様の 1 本スプラインで絞り区間が切れる等) はそのとき使う。

経路は repo に実績がある (case/37): GEOM → SMESH → MED → `med_to_msh41.py` → msh4.1 → `convertGmshToForge` (node で tet/prism/pyramid 対応済み)。

### 4.8.1 **ユーザ決定 (2026-09-22, 再): 全ヘキサを本線にする。Salome は保険**

§4.8 の後、非構造の代償 (層数が板厚で縛られる・系統的細分列が作れない・節点数) を見たユーザが
「メッシュ数多くなるなら、全ヘキサでやる方針にしようかな」と述べ、全ヘキサの試作形状 (§4.7, 20 ブロック・負 Jacobian 0・skew 0.501) の
図を見て「いい感じ」と評価した。**優先を次のとおり確定する**:

1. **本線 = 全ヘキサ・ブロッキング (B0 → B1 → …)**。フィレットは上隅・下隅で**独立な r(x)** (付けない側は 0、始端は徐変可)。
   終端は**接する 2 壁が一緒に終わる隅は鈍頭可** (下隅は `L_sw = L_cowl` のとき)、**片方だけ終わる隅は手前 0.5–1 $r_f$ で絞って消す**
   (絞り長さは物理入力)。
2. **保険 = Salome 非構造** (`case/46.sern_design/cad/salome_junction_demo.py`)。用途は (i) 全ヘキサで切れない形の試行、
   (ii) 鈍頭終端と「絞って消す」終端の差を同条件で実測する対照。S2 以降 (受入試験の整備・全体展開) は必要になるまで保留。

### 4.9 Salome スパイク: 形状の自由度は得られた。**粘性層の総厚は最薄の板厚より薄くする必要がある** (2026-09-22, scratchpad)

**形状**: 接続模型 (ランプ体 ∪ 側壁 5 mm ∪ カウル板 2 mm、H = 0.1 m)。フィレットは
**「隅の断面 (角 + 円弧の曲線三角形) を x の各位置で r(x) に合わせて描き `MakeThruSections` でロフトしたガセット固体」**
を固体に融合して作る。上隅 = **徐変で始まり (1 → 10 mm) 終端まで付けたまま (鈍頭)**、下隅 = **8 mm 一定 → 終端で消す (→ 1 mm)**。
流体 = 箱 − 固体。**`CheckShape` valid・流体 solid 1 個・面 15 枚**。
(gmsh-OCC の可変半径 `occ.fillet` は同じ形状で「Could not compute fillet」で失敗した。)

**メッシュ** (NETGEN 1D-2D-3D + ViscousLayers、層は壁 9 面のみ・外枠 6 面は除外):

| 仕様 | 粘性層 | 結果 |
| --- | --- | --- |
| a | 20 µm × 20 層 (1.2) = 3.7 mm | プリズム 20 万のみ、**tet 0** (`Compute: True` なのに充填が黙って失敗) |
| b | 4 µm × 40 層 (1.15) = 7.1 mm | 同上 (プリズム 40 万、tet 0)。到達層厚は平均 3.3 mm |
| **c** | **4 µm × 22 層 (1.2) = 1.08 mm**、壁の接線サイズ 1.5 mm | **tet 1,085,035 + prism 2,304,698、1,381,415 節点。層厚未達の警告なし** |

**tet 0 の無言失敗は case/49 の記録と同じ症状**。今回の条件では
**粘性層の総厚 > 最薄の板厚 (カウル 2 mm)** で起き、総厚を板厚未満にすると通った。
→ **受入試験には「tet 数 > 0」「層厚未達の警告なし」を必ず入れる** (`Compute: True` は合格の根拠にならない)。

**forge への取り込みと品質** (仕様 c。`med_to_msh41.py` の必須グループを汎用化したコピー → `convertGmshToForge` cell 変換):

| 項目 | 値 |
| --- | --- |
| セル | 3,389,733 (wedge 2,304,698 + tetra 1,085,035)、境界面 151,508 |
| AR | **max 686.1** / p99 393.7 (>5000 は 0) |
| skew | max 0.936 / p99 0.460、**>0.90 が 46 セル (0.00 %)** |
| **`check_mesh_quality.py --ar-max 5000`** | **VERDICT: SOFT-PASS (<0.1% outliers)** |
| 壁第一層 (`--geometry-only`, cell 重心基準) | 中央値 **4.39e-06** / 平均 4.67e-06 / 最小 1.03e-06 / **最大 6.68e-04** |

**この方式の代償 (隠さず書く)**:

1. **層数と成長率が板厚で縛られる**: 40 層 × 1.15 (総厚 7.1 mm) は張れず、22 層 × 1.2 (1.08 mm) になった。
   境界層の外側 (数 mm) は等方的な tet (壁近傍 1.5 mm) で受けることになる。
2. **節点数が大きい**: この小さな接続模型 (0.24 × 0.2 × 0.2 m) だけで 138 万節点。壁面の三角形が等方 (1.5 mm) なので、
   SERN 全体 (ランプ ~1 m) では **1000 万節点級**になる見込み (全ヘキサは 322 万節点だった)。AWS の GPU (23 GB) には載るが、変換のホスト RAM が要る。
3. **系統的な格子細分列が作れない** (case/49 の実測: VL 総厚と接線サイズが結合する)。§8 の格子収束の受理は
   「入れ子でない 3 水準」で代替することになり、判定力は落ちる。
4. skew > 0.90 が 46 セル、壁第一層の最大 6.7e-04 が残る (フィレットの徐変始端/消える終端の微小面が疑わしい)。要調査。

### 4.10 B0 接続設計 — 接続模型のブロック構成と共有面 (2026-09-22, 実装前に固定)

§4.7 の試作 (ダクト内 N の 4 ブロック) を、**固体を有限厚で持つ接続模型**へ拡張する。codex plan レビュー M1–M4・M7 への回答を兼ねる。
§4.1–§4.3 の記述と食い違う点は**本節が優先**する (食い違い: 外部の y 線、フィレットの持ち方、上バンドの 2 分割をやめた点)。

**(a) 形状パラメータ (すべて物理入力, /H)**

| 記号 | 意味 | 試作値 |
| --- | --- | --- |
| $y_r(x)$, $y_c(x)$ | ランプ輪郭 / カウル上面 | 曲線 / 0 |
| $t_c$, $t_{sw}$ | カウル板厚 / 側壁板厚 (**外向き**: $z \in [z_w, z_o]$, $z_w = W/2$, $z_o = z_w + t_{sw}$) | 0.02 / 0.05 |
| $L_{sw} < L_{cowl} < x_{end}$ | 側壁後端 / カウル後端 / 模型の出口 | 1.2 / 1.6 / 2.4 |
| $r_u(x)$, $r_l(x)$ | 上隅・下隅のフィレット半径 (独立。$x \ge L_{sw}$ では 0) | §4.7 と同じ |
| $\delta_r$ | 壁層帯の厚さ (N のリング厚、外部の壁帯の厚さ) | 0.08 |
| $y_{bot}$, $Z_{far}$ | 外部領域の下端 / 側方端 | −1.0 / 2.0 |

カウル板は $z \in [0, z_o]$ (側壁の下まで延びる)。$y_{cl} = y_c - t_c$ をカウル下面とする。
**事前判定** (満たさなければ生成前にエラー): $0 < L_{sw} < L_{cowl} < x_{end}$、$r_u, r_l \ge 0$ かつ $r_u + r_l + 2\delta_r < y_r - y_c$、
$\max(r_u, r_l) + \delta_r < z_w$、$t_c, t_{sw} > 0$、$x \ge L_{sw}$ で $r_u = r_l = 0$、各区間長 ≥ 下限、$y_{bot} < y_{cl} - \delta_r$、$z_o + \delta_r < Z_{far}$。
**コアが流体内に残る条件** (2026-09-22 codex m6): コア角 `CUR` $= (y_r - \delta_r, z_w - \delta_r)$ と円弧中点 `WU` は $r = \delta_r/(1 - 1/\sqrt2) = 3.41\delta_r$ で一致し e2 が潰れる。
各 station と区間中点で **全辺長 > 0・e2/e6 の長さ ≥ $\delta_r/2$・輪郭の非交差** を検査する (試作値 $\delta_r = 0.08H$ で $r < 0.20H$)。

**(b) 断面 (y–z) の胞体複体** — 各 station で同じ名前の点・辺・ブロックを作る (位相は全 station で同一)

| 領域 | ブロック | 範囲 | 壁法線の細分 |
| --- | --- | --- | --- |
| **N** ダクト内 | `N_TOP` / `N_SIDE` / `N_BOT` / `N_CORE` (バタフライ C リング + コア。§4.7) | $z \in [0, z_w]$, $y \in [y_c, y_r]$ | リング横断 (辺 e2/e4/e6/e9) に $h_1$ から幾何級数 |
| **SW** 側壁の跡 | `SW` | $z \in [z_w, z_o]$, $y \in [y_c, y_r]$ | y は e5 と同一分布、z は $n_{sw}$ 分割 |
| **CW** カウルの跡 | `CW1` ($z \in [0, z_w]$) / `CW2` ($z \in [z_w, z_o]$) | $y \in [y_{cl}, y_c]$ | y は $n_{tc}$ 分割 |
| **U** カウル下 | `U1a`/`U1b` ($z \in [0,z_w]$), `U2a`/`U2b` ($z \in [z_w,z_o]$)。a = $[y_{bot}, y_{cl}-\delta_r]$、b = $[y_{cl}-\delta_r, y_{cl}]$ | | b 帯の上端 (カウル下面) へ $h_1$ |
| **S** 側壁の外 | `Sb*`/`Sc*` × {a, b, c, de}。b = $z \in [z_o, z_o+\delta_r]$、c = $[z_o+\delta_r, Z_{far}]$。y は a, b (U と同じ), c = $[y_{cl}, y_c]$, de = $[y_c, y_r]$ | | b 列の左端 (側壁外面・カウル側端) へ $h_1$、de の上端 (機体下面) へ $h_1$ |

計 **4 + 1 + 2 + 4 + 8 = 19 ブロック/断面**。**`z = z_w` の分割を U (U1/U2) と CW (CW1/CW2) に通す** (codex M2)。
S の上バンドは分けず `de` 1 本にして、**e5 (側壁内面の壁曲線) と同じ節点数・同じ両端細分**を使う (下記 (d))。

**(c) x 区間とブロックの有無** — 端点は物理入力そのもの (最近傍丸めをしない。M4)

| 区間 | N (4) | SW | CW1, CW2 | U (4) | S (8) | 計 |
| --- | --- | --- | --- | --- | --- | --- |
| A = $[0, L_{sw}]$ (フィレットの節点でさらに分割) | ✓ | — (固体) | — (固体) | ✓ | ✓ | 16 |
| B = $[L_{sw}, L_{cowl}]$ | ✓ | ✓ | — (固体) | ✓ | ✓ | 17 |
| C = $[L_{cowl}, x_{end}]$ | ✓ | ✓ | ✓ | ✓ | ✓ | 19 |

**壁が終わっても格子帯を潰さない**: N のリングは B・C でも同じ形で続き、壁だった面が内部の共有面になる。
フィレットは A の末尾で $r \to 0$ に絞る (上下独立)。$L_{sw} = L_{cowl}$ で下隅を鈍頭にする形 (§4.7 末尾) は本模型の後に扱う。

**(d) 共有面の一覧** — 「同じ幾何エンティティを 2 つの体積が参照する」ことで節点 ID・分布・向きを一致させる

| 面 | A | B | C |
| --- | --- | --- | --- |
| e5 × x (N_SIDE の外面, $z = z_w$, $y \in [y_c,y_r]$) | 壁 `sidewall_in` | **N_SIDE ↔ SW** | N_SIDE ↔ SW |
| SW の右面 ($z = z_o$) | (壁 `sidewall_out` = Sb_de の左面) | **SW ↔ Sb_de** | SW ↔ Sb_de |
| SW の下面 ($y = y_c$, $z \in [z_w,z_o]$) | — | 壁 `cowl_in` (側壁の跡に露出したカウル上面) | **SW ↔ CW2** |
| SW の上面 ($y = y_r$) | — | 壁 `vehicle` (露出したランプ面) | 同左 |
| e8 × x (N_BOT の壁面, $y = y_c$, $z \in [0,z_w]$) | 壁 `cowl_in` | 壁 `cowl_in` | **N_BOT ↔ CW1** |
| CW1 ↔ CW2 ($z = z_w$) / CW2 ↔ Sb_c ($z = z_o$) | — / 壁 `cowl_side` | — / 壁 `cowl_side` | 共有 / **共有** |
| CW1・CW2 の下面 ($y = y_{cl}$) | 壁 `cowl_out` (= U1b・U2b の上面) | 同左 | **CW1 ↔ U1b、CW2 ↔ U2b** |
| U2* ↔ Sb_a/Sb_b ($z = z_o$, $y < y_{cl}$)、U1* ↔ U2*、Sb* ↔ Sc*、a ↔ b ↔ c ↔ de の上下 | 共有 | 共有 | 共有 |
| 断面 (x = const) の面 | 区間の境で隣の区間と共有。**片側しか体積が無い面が端面**: x = $L_{sw}$ の SW 断面 = 壁 `sidewall_end`、x = $L_{cowl}$ の CW1/CW2 断面 = 壁 `cowl_base`、x = 0 = 入口、x = $x_{end}$ = 出口 |

**N 内部の共有面** (全区間で共有。2026-09-22 codex m5 で明示): e2×x (TOP↔SIDE)、e6×x (SIDE↔BOT)、e3×x (TOP↔CORE)、e7×x (SIDE↔CORE)、e10×x (BOT↔CORE)。
断面の内部共有辺の数は **A 21 / B 23 / C 29** (SW で +2、CW1/CW2 で +6。codex が独立に数えて整合を確認)。

**外部面のタグ表** (参照する体積が 1 つの面。これ以外に参照数 1 の面が出たら失敗):

| 面 | タグ | 区間 |
| --- | --- | --- |
| e1×x (上フィレットのランプ側半分を含む) / SW 上面 ($y = y_r$, $z \in [z_w, z_o]$) | `ramp` | A–C / B–C |
| $y = y_r$, $z > z_o$ (S の de 上面) | `vehicle` | A–C |
| e5×x (上下フィレットの側壁側半分を含む) | `sidewall_in` | A |
| S の b 列 de の左面 ($z = z_o$) | `sidewall_out` | A |
| e8×x (下フィレットのカウル側半分を含む) / SW 下面 | `cowl_in` | A–B / B |
| U1b・U2b の上面 ($y = y_{cl}$) | `cowl_out` | A–B |
| S の b 列 c の左面 ($z = z_o$, $y \in [y_{cl}, y_c]$) | `cowl_side` | A–B |
| x = $L_{sw}$ の SW 断面 / x = $L_{cowl}$ の CW1・CW2 断面 | `sidewall_end` / `cowl_base` | — |
| $z = 0$ (e4, e11, e9、U1a・U1b の左面、C の CW1 左面) | `symmetry` | A–C |
| x = 0 の N 断面 / x = 0 の U・S 断面 | `inlet` / `ext_in` | — |
| x = $x_{end}$ の全断面 | `outlet` | — |
| $y = y_{bot}$ / $z = Z_{far}$ | `far_bottom` / `far_side` | A–C |

共有面は 1 つの正準向きで持ち、各体積からは符号付きで参照する。**参照数に加えて、体積ごとに頂点 Jacobian の符号が揃っていること**を検査する
(向きの混在は msh 出力の前に節点順を直す)。

**節点数の拘束** (向かい合う辺は同数): z 方向 $N_Z$ (e1, e3, e8, e10, U1 列, CW1)、側壁の跡 $n_{sw}$ (SW, CW2, U2 列)、
y 方向 $N_Y$ (e5, e7, e11, SW, S の de)、カウル厚 $n_{tc}$ (CW, S の c)、壁帯 $N_L$ (リング横断、U の b、S の b 列、S の de は $N_Y$ に含む)、
外部の粗い帯 $n_{far,y}$ (a)、$n_{far,z}$ (S の c 列)。**境界タグは「1 つの体積しか参照しない面」を機械的に拾い、(辺名 or 断面ブロック名, 区間) の表でタグを決める**
(全外部面がちょうど 1 タグに属することを試験する)。

**(e) 分布**

- リング横断・U の b・S の b 列: 壁側第一層 $h_1^{eff}$、成長率は**辺ごとに**辺長から解く (§4.7)。
- (**2026-09-22 追記: 下の $s_1$ の決め方と隅の station は §4.11 の 1・2・8 で置き換えた**。目標が法線距離であることは変わらない)
- **第一層の目標は「辺に沿った長さ」ではなく、実際に生成する壁面の 3D 法線に対する距離 $|\hat n\cdot\Delta x| = h_1$** (2026-09-22 codex M1。
  試作は辺長を 4 µm にしていたので、$r = 0$ の隅では対角辺 e2/e6 の第一節点の壁距離が $4/\sqrt2 = 2.83$ µm = **−29.3 %** だった。
  ランプ角の $1/\cos\theta_r$ だけではこれも、フィレット遷移面の $r'(x)$ の寄与 (水平ランプでも $0.10H \to 0$ を $0.05H$ で絞ると円弧中点で 0.770) も直らない)。
  壁法線辺 (e2, e4, e6, e9 と外部の壁帯) の第一間隔は station ごとに $s_1 = h_1 / |\hat n \cdot \hat d|$ で決める
  ($\hat d$ = 辺の断面内単位方向、$\hat n$ = 壁面パラメータ化 $X(x,\alpha)$ の数値微分 $X_x \times X_\alpha$ から得る 3D 単位法線)。
  **鋭い隅 ($r < h_1$) では接する両壁の条件を同時に課す**: $s_1 = h_1 / \min_w |\hat n_w\cdot\hat d|$ (水平の直角隅で $\sqrt2 h_1$)。
  $h_1 \le r$ では円弧の法線 1 本で足りる (対角辺は円弧の半径方向)。
- **検査は station の辺だけでなく壁面の全節点で行う**: 壁の四角形から節点法線を作り、その壁に接するヘキサの対向節点までの $|\hat n\cdot\Delta x|$ を
  **壁タグ別・節点別**に出す。許容は「隅・フィレット遷移・端面の稜から $2\delta_r$ 以上離れた壁で ±5 %、それ以外で ±10 %」。
  station 間は transfinite 補間で $s_1$ が線形に変わるので、絞り区間で外れたら **station を足す** (station は物理位置なので解像度に依らない)。
- 両端細分の辺 (Bump) は両端の間隔が同じになる。ランプ側とカウル側で $1/|\hat n\cdot\hat d|$ が違うので、**両端の平均** $h_1 (1 + 1/\cos\theta_r)/2$ を端の間隔に使う
  (試作のランプ角 16.7° で ±2.1 %、22° で ±3.8 %。±5 % に入らない角度では辺を分割する)。
- **e5 (側壁内面に沿う y 分布) は両端細分** ($h_1$ から成長率 ≤ 1.2)。B・C では SW がこの分布を引き継ぎ、
  **側壁の跡に露出したランプ面・カウル上面の壁層**になる (一様だと露出面の第一層が mm 級になる)。
  **節点数は同値類 (e5/e7/e11・SW・S の de) の全 station の最長辺から決める** (2026-09-22 codex M2: 模型出口で e5 は 0.160 m、90 節点・成長率 1.2 では 0.134 m しか覆えない)。
  辺は分割せず gmsh の `Bump` 分布を使う: 間隔が弧長の放物線 $h(s) \propto b - a(s - L/2)^2$、端/中央の比 = 係数、**最大成長率 ≈ $1 + 4h_{mid}/L$**。
  実測 (L 0.10 m・端 4 µm): 91 節点で 1.210、111 節点で 1.163、131 節点で 1.133。係数は station ごとに端の間隔から二分法で解き、**向かい合う辺には同じ係数**を使う (格子線がまっすぐ通る)。
  同じ理由で **e1・e8 は隅側へ片側細分**する (鋭い隅では TOP/BOT の接線間隔が側壁の法線解像になるため。e3・e10・U1 列・CW1 は同じ比で追従)。
  SW の z ($t_{sw}$)・CW の y ($t_c$) も両隣が壁層なので `Bump` で両端を $h_1$ に合わせる (界面で間隔が 100 倍跳ぶのを避ける)。
- **端面の x 分布**: `sidewall_end`・`cowl_base` は壁なので、**端面の第一層も $h_1$ (法線 = x)** にする (2026-09-22 codex M3:
  当初の「第一間隔 $\le t/5$」は後流を板厚で分解する条件で、壁法線解像の条件ではない。$t/5$ は 4 µm の 250 倍/100 倍を許してしまう)。
  条件は $\Delta x_1 = h_1$ (±5 %) **かつ** 端面から $t$ の範囲で $\Delta x \le t/5$、成長率 ≤ 1.2。上流側の区間末尾も同じ間隔に寄せて段差を作らない。
  区間 B は両端が端面なので `Bump`、A の最終区間は終端へ・C は始端から片側細分。**同じ区間の全ブロックは x 方向の節点数と分布を共有**し、
  実際の x 間隔の成長率と区間境界での連続 (隣接比 ≤ 1.2) を検査する。
- 壁が無くなった下流の壁帯は物理長で粗い値へ移す (df/dx ≤ 4e-3。[[mesh-wall-layer-x-blend]])。**本模型は短いので未適用、全体展開 (B4) で入れる**。

**(f) 実装** — `case/46.sern_design/cad/hex_junction_demo.py` を汎用化する: station ごとに {点, 辺, ブロック} の辞書を作る関数、
区間ごとの「存在するブロック」集合、辺ごとに 1 枚の x 面・ブロックごとに 1 つの体積を作る共通ビルダ (辞書でキャッシュして二重生成しない)、
参照数 1 の面からの境界タグ付け、msh4.1 出力。**方式確定後に `design/forge_design/meshing/` へ移す**。

### 4.11 B1 の結果 — 接続模型は 2 解像度で受入ゲートを通った。**生産解像度は手元に載らない** (2026-09-22)

実装は `case/46.sern_design/cad/hex_junction_model.py` (§4.10 の 19 ブロック × 3 区間。検査も同じスクリプトに入れた)。
出力は `case/46.sern_design/cad/demo_out/jm_s05*`・`jm_s071*` (git 追跡外・再生成可)。

| `--scale` | $h_1$ | 節点 / ヘキサ | ブロック | 負 Jacobian (float64 / float32) | skew max | AR max (>1000 のセルの skew) | 位相・タグ | 壁第一層 (9 タグ) | VERDICT | gmsh 最大 RSS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 8 µm | 1,932,448 / 1,873,258 | 216 | 0 / 0 | 0.530 | 2575 (≤ 0.251) | 境界面 117,496 = タグ付き四角形 117,496、タグ無し 0 | 全タグ合格 | **PASS** | 2.6 GB |
| 0.71 | 5.6 µm | 4,966,050 / 4,853,569 | 216 | 0 / 0 | 0.534 | 3253 (≤ 0.229) | 223,742 = 223,742、タグ無し 0 | 全タグ合格 | **PASS** | 6.4 GB |
| 1.0 | 4 µm | (約 1250 万と外挿) | — | — | — | — | — | — | **未実行** (手元 11 GB では gmsh が 90 秒で 8.9 GB に達し中断) | 約 16 GB と外挿 |

壁第一層の実測 (法線距離 / $h_1$。scale 0.71): 平面部 ramp 0.987–1.001・sidewall_in 1.003–1.023・cowl_in 1.000・外部と端面の 6 タグは 1.000 (vehicle のみ 0.979–0.991)、
フィレット遷移区間 0.962–1.061、フィレット円弧上 0.935–1.288。

**設計の変更点 (§4.10 からの差分。すべて実測で決めた)**

1. **第一層は「辺長の比」で決まる**。gmsh の transfinite 補間は層 1 を $X = (1-\eta)\,$壁$\, + \eta\,$内側線、$\eta = s_1/L$ (両側辺の間で線形) に置く。
   内側線が平面壁から $\delta_r$ の平行線なら、**リングの全側辺で $\eta = h_1/\delta_r$ に揃えると平面部の法線距離は厳密に $h_1$**。
   対角辺は $s_1 = h_1 L_{e2}/\delta_r$ ($r=0$ で $\sqrt2 h_1$ = codex M1 の両壁条件と一致し、$r$ について連続なので、$r=0$ の手前の特別な station は不要になった)。
   §4.10 (e) に書いた $s_1 = h_1/|\hat n\cdot\hat d|$ を辺ごとに課す方式は、平面部が対角辺に近づくほど薄くなり (−26 %)、逆に両壁条件だけだと厚くなる (+50 %) ので**不採用**。
2. **内側線は壁のオフセット曲線にする** ($r > \delta_r$ で同心の円弧 $r - \delta_r$)。直線のままだと $r > \delta_r$ の円弧上でリングが $\delta_r$ より薄く、第一層が −20〜−30 % になる。
   $r$ が $\delta_r$ を横切る位置は内側線の角が円弧に変わる折れ点なので **station にする**。
3. **隅へ向けた接線方向の細分 (e1・e8) はやめた**。バタフライでは細分がコアの内部 (壁から $\delta_r$ 離れた線) に伝播する。e5 の両端細分も同様に伝播するが、
   x 間隔を $h_1$ の 750 倍に抑えればコア内の細い行は AR < 1000 で直交なので許容した (e7/e11 は e5 と同じ係数 = 扇形なし)。
4. **分布を変えるのは長いブロックをまたぐ所だけ**。薄い帯 ($\delta_r$, $t_c$, $t_{sw}$) をまたいで分布を変えると skew 0.76。遠方の最大格子幅に上限 (`HFAR`) を置く (端面層が断面全体に伝播し、遠方の幅が AR を決める)。
5. **e5 の一様 → 両端細分の移行は側壁の全長で、間隔について線形に行う**。対数で混ぜると節点の移動が上流に偏り、壁面の格子線が x–y 面で 45° 傾く。
6. **壁層セルの接線幅に上限 ($900\,h_1$)**。リングは最大 45° 傾くので、壁層でも AR ≤ 1000 に収める (N_Y が 101 → 133 級に増える)。
7. **半径を変える区間の長さは $2r$ 以上** (物理入力への制約)。短いと $r'(x)$ とランプ勾配で壁法線が x に傾き、円弧上の第一層が薄くなる (長さ $r$ で −19 %、$2r$ で −9 %)。
   対角辺にはこの傾きの補正を 3 割だけ掛ける (全部掛けると平面部が厚くなる)。
8. **判定区分**: 平面部 ±5 %、稜から $2\delta_r$ 以内と半径が変わる区間は ±10 %、**フィレット円弧上は −10 % / +45 %** (バタフライ位相では 45° 点でリングが最大 $\sqrt2$ 倍厚い。薄くはしない)。
   節点法線は、タグの境が滑らか (50° 未満) なら隣のタグの面も入れて作る (45° 点は ramp と sidewall_in の境)。
9. AR の判定は「> 5000 が 0、かつ AR > 1000 のセルの skew ≤ 0.30」(近直交の層だけを緩和の対象にする、の実装)。AR > 1000 は全て端面層が遠方へ伝播したセル。

**node 変換後の検査** (scale 0.5、`case/46.sern_design/cad/demo_out/conv_s05/`。`convertGmshToForge` 42 秒・最大 RSS 7.1 GB = **3.7 kB/節点**):

- 変換器ログ: 体積和の相対差 9.5e-6、全域正規化の閉性 4.9e-7。
- **CV ごとの検査 (新規ツール `solver_density_cuda/tools/check_dual_closure.py`)**: 双対体積は全 CV で正。閉性 $|\Sigma S|/\Sigma|S|$ は median 2.9e-9・p99 2.0e-6・**max 3.9e-4** (> 1e-5 が 13,050 CV、> 1e-3 は 0)。
  → §6-3 の 1e-5 には**不合格**。最悪は端面直前の壁 CV (Δx 8 µm × Δy 8 µm × Δz 5 mm)。誤差の大きさは「座標の float32 の 1 ulp (1.5e-8) × 辺長 5 mm = 7.5e-11 m²」に一致し、
  変換器が座標を `stof` で float32 に読んで双対面を単精度で組んでいること (`mesh/gmshReader.hpp:504`) が原因と見ている (未検証)。
- **同じツールを旧メッシャの格子 (`case/46.sern_design/run_0402_3d_wf1/sern.h5`) に掛けると、閉性 0.30 の「開いた CV」が 68 個出た**。全部 $z = W/2$・$x \in [L_{sw}, L_{cowl}]$ の
  カウル板の露出した側端 (cowl_in 34 + cowl_out 34)。旧メッシャは板の自由端でも節点を二重化している (`mesh_sern3d.py` の `dup1` が $k = k_{sw}$ を含む)。
  **新しい接続模型には開いた CV は無い**。旧メッシャは同日修正した (親 plan §4.45 / R5r: 再現格子で 18 → 0)。**親 plan R5q の床ノードとは場所が別** (R5q は $x < L_{sw}$ の T 字部) なので、R5q の原因とは言えない。
- `check_mesh_quality.py` は **FATAL (体積の符号が混在: 負 151)** を返した。該当セルは頂点 Jacobian が 8 つとも正の妥当なヘキサで、上フィレットの壁層にある。
  ツールの 5 四面体分割は向かい合う面を別の対角線で切るので、**面の反り (弦の矢高 312 µm) が層厚 (8 µm) より大きい薄いセル**では分割した四面体どうしが交差して和が負になる。
  codex が指摘した「反例を通す」の逆向きの誤判定。フィレット円弧を接線方向 1.5 セルで切っているのも粗い (B2 で接線分割を円弧長から決める)。

**規模の見積もり**: 節点数は scale の約 2.7 乗で増え、scale 1.0 は約 1250 万。gmsh 1.3 kB/節点 → 16 GB、変換器 3.7 kB/節点 → **46 GB**。
長さ 2.4 H の接続模型だけでこの規模になる理由は、(i) 全ての壁 (側壁の内外・カウルの上下と側端・機体下面・端面 2 枚) に $h_1$ の層を置き、(ii) 構造格子なので層が断面と x の全域に伝播するため。
x 方向 231 面のうち 190 面は端面 2 枚の両側の層。**SERN 全体にそのまま展開すると現在の計算環境 (手元 11 GB、AWS g5.xlarge 16 GB) には載らない**。削る候補と効果は未確定事項 U4。

### 4.12 **ユーザ決定 (2026-09-22): 第一層を粗くして 200〜300 万節点に収める。検証できる規模を優先する**

ユーザ: 「4 µ も小さくしないといけない理由？ 粗くして、200 か 300 万くらいにしないと検証もできなそうなので、そのように」。

**4 µm の出どころ**: AGENTS.md の既定目標「局所 $y_1^+ \le 1$」。g4 (4 µm) の実測 $y_1^+$ は ramp 0.549 / cowl_in 1.234 / cowl_out 0.185 で、cowl_in を 1 付近にするための値だった。
**緩める裏付け** (同じ規則が求める格子感度): 親 plan §4.42 の g3 (16 µm) → g4 (4 µm) で ΔC_T −0.00013・ΔC_T(摩擦込み) −0.00024・ΔC_L +0.00145・ΔC_M −0.0417、
いずれも許容 (0.002 / 0.002 / 0.05) の内側。**力の係数を見る限り 16 µm で足りている** ($y_1^+$ は約 4 倍: ramp 約 2.2、cowl_in 約 5)。壁熱流束の感度は未確認なので、熱を評価に使うときは測り直す。

**採用した設定** (`hex_junction_model.py --set H1=1.6e-4 H1_END=2.5e-3 HX=0.05 HX_FAR=0.08 NZ=31`):
第一層 16 µm (成長率 1.2 は据え置き)、**端面の第一層は 0.25 mm** ($t_c/8$。後流を板厚で分解する条件 $t/5$ は満たす。codex M3 の「端面にも壁層」はユーザ決定で見送り。端面の上は剥離域で壁せん断が小さい)、
x 間隔 5 mm、z 方向 31 節点。効いたのは端面の緩和 (x 面数 189 → 108)。第一層を 4 倍にしても層数は 37 → 30 にしか減らない (成長率 1.2 の等比なので対数でしか効かない)。

| 設定 | 節点 / ヘキサ | skew max | AR max | 負 Jacobian (f64/f32) | 位相・タグ | 壁第一層 | 模型の VERDICT | gmsh RSS | 変換 RSS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16 µm・端面 0.25 mm (`demo_out/jm_c16*`) | **2,446,572** / 2,378,448 | 0.533 | **676** (> 1000 は 0) | 0 / 0 | 一致・タグ無し 0 | 全 9 タグ合格 (平面部 0.979–1.022、円弧上 0.934–1.232) | **PASS** | 3.4 GB | 9.0 GB |

node 変換後 (`case/46.sern_design/cad/demo_out/conv_c16/junction.h5`):
`check_mesh_quality.py` → **`VERDICT: PASS (AR<= 1000, skew<= 0.90)`** (AR 緩和なし。max 676 / skew max 0.533)。
`check_dual_closure.py` → 体積は全 CV 正、閉性 max 1.7e-4 (> 1e-5 が 6,410 CV) = **FAIL** (変換器の単精度。B1b)。変換器ログは体積和の相対差 5.0e-6・全域閉性 6.0e-7。

**`check_mesh_quality.py` の体積判定を直した** (B0c の残り)。六面体は 5 四面体分割の和をやめ、頂点 Jacobian の平均を体積、**符号の混在を「局所反転」として FATAL** にする。
旧判定は (i) フィレット壁層の妥当なセル 111 個を「向きの不整合」と誤判定し、(ii) 逆に 1 頂点を押し込んだ局所反転を見逃していた (両方とも回帰試験に追加。旧ツールで 2 件とも NG、新ツールで OK)。
station の置き方も 1 点直した: $r = \delta_r$ の折れ点 station がある区間では中点 station を足さない (station が詰まると x 間隔が 1.35 倍飛ぶ)。scale 0.5 は再実行して PASS のまま。

### 4.13 接続模型の初回 CFD (AWS, 2026-09-22 着手)

目的は力の係数ではなく (形状が MOC 設計でない試作なので生産値と比べられない)、**有限厚の板と接続部で計算が成立するか**: 発散しないこと、密度・圧力の床に張り付く節点が 0 であること、壁解像。
ドライバは `case/46.sern_design/cad/run_junction_model.py` (気体状態・ソルバ設定・段階起動は生産 runner のものを使い、格子と境界タグだけ差し替える)。
設定は `problem_3d_prod_m6on_wallres.yaml` (frozen_tp・等温壁 1000 K・SST・本段 CFL 0.25・層流暖機 2000 → soft 2000 → mid 2000 → 本段 4000)。

| run (AWS `~/forge-sern/case/46.sern_design/`) | 差分 | 結果 |
| --- | --- | --- |
| `run_0423_3d_junction_c16` | IC = 領域別の一様流 (生産 runner と同じ作り) | **発散**。層流暖機 (1 次・CFL 0.2) の step 454 で NaN。NaN 254 節点は**カウル後端面の直後** (x/H 1.6006, y/H −0.003〜−0.006, z/H 1.025〜1.03 = CW2)、ρ min 3.95e-05 は**側壁後端面の直後** (x/H 1.2006, SW)。暖機中の残差は下がっていなかった (`rms_roe` 81 → 44 で横ばい、`rms_ro` は微増) |
| `run_0424_3d_junction_c16_wakerest` | + **後流ブロック (SW・CW) を静止で始める** | **発散**: step 2 で NaN 41,014 節点。後流ブロックの壁節点 (露出したカウル上面・ランプ面) と、動く流体との境界 (z = W/2 と z_o の 16 µm 壁層) に 1600 m/s の速度不連続を置いたため。ピストン仮説の A/B にはならなかった |
| `run_0425_3d_junction_c16_wakeblend` | + **端面からの 3D 距離 d で速度を smoothstep ($d < 10t$) で立ち上げる** (不連続を作らない) | **発散**: step 2 で NaN。速度の不連続は原因でなかった |
| (診断, 一時 run・破棄済) | run_0425 の IC で暖機を 3 step、1 step ごとに出力 | **step 1 で ρ = −0.28** at (x/H 1.6124, y/H −0.010, z/H 1.009) = カウル後端面の直後・**板厚の中央**。そこは IC の排気 (101 kPa, 2328 K) / 外気 (2851 Pa, 221 K) の境目で、35 倍の圧力段差が 60 µm 隣の節点間にある。静止させると受け側 (冷・c 300 m/s) の局所 dt が、高圧側から来る波 (c 900 m/s + 流速) に対して 3 倍以上大きく、1 step で CV の質量を超える流入になる。run_0423 は流速 1788 m/s が dt を 6 倍小さくしていたので 454 step 生きた。**ピストン仮説は未検証のまま** (段差の方が先に効いた) |
| `run_0426_3d_junction_c16_smoothic` | + **端面の後ろでは排気/外気の境目を幅 5 mm の smoothstep で混ぜる** (保存量の線形混合) | **暖機 1720 step で NaN** (残差は低下中)。IC の問題は解消し、**側壁後端面の壁 CV の排出**が残った (下の診断) |

**run_0423 の見立て (仮説。run_0424 が A/B)**: 一様流の IC だと、no-slip の後端面から流体が 1620〜1790 m/s で遠ざかる。これはピストンを引き抜くのと同じで、外部流 (音速 約 300 m/s) では
逃げ速度 $2c/(\gamma-1) \approx 1500$ m/s を超えるので、1 次元では端面に真空ができる。旧メッシャのカウル後縁は厚さ 0 で端面が無かった (機体ベースはテーパ付き) ので、この起動問題は出ていなかった。


**run_0426 (境目を滑らかにした IC) の結果と診断 (2026-09-22)**: 暖機 1720 step で NaN。残差は下がっていた (`rms_roe` 100 → 9.2)。
同じ IC・設定で 200 step ごとに場を出した診断 run (`run_0427_3d_junction_diag200`, 一時 run) で順序が取れた:

| step | ρ min | P min | 場所 |
| --- | --- | --- | --- |
| 600 | 2.8e-3 | 736 Pa | 外部流の膨張 (機体下面。試作ランプの 16.7° をそのまま外気が回る) = 想定内 |
| **800** | **4.9e-4** | **164 Pa** | **側壁後端面の壁節点** (x = L_sw。`centCoords` は双対重心なので 0.06 mm 後ろに写る)、y = 0.9 mm (露出したカウル上面から)、z = 0.25 mm (内側の縁から) |
| 1000 | 1.8e-4 | 59 Pa | 同じ節点。周囲の壁節点へ広がる (ρ < 5e-4 が 42 節点) |
| 1400 | 3.0e-5 | 10 Pa | 同じ場所。**床 (pMin 20 Pa, roMin 1e-4) を割る** |
| 1600 | 4.2e-5 | 14 Pa | ランプ側の隅 (y = 1.31 H) も同じ経路で落ちる |
| 1720 | NaN 207 節点 | | |

**機序** (壁節点の CV の排出。[[base-wake-resolution-rule]] と同じ型で、第一 station を $t/20$ にしても起きる):
後端面の壁節点は `nodeWallDirichlet: 1` で $u = 0$・$T = T_w$ に固定され、**密度だけが自由**。その 0.26 mm 後ろの内点は排気が内側の縁を回り込む流れで
$|U| \approx 1550$ m/s (面から遠ざかる $u_x$ 850 + 縁を回る $u_z$ 1290)。壁 CV と内点の間の面の質量流束は外向き一方で、壁節点の運動量が固定なので
**圧力勾配 (壁 67 Pa vs 内点 10 kPa) が流体を壁 CV へ戻せない**。排出率 = $u_{face}\,\Delta t/\Delta x_{half}$ ≈ 425 × 4.6e-9 / 0.125 mm ≈ 1.5 %/step (実測 P 2000 → 164 Pa / 200 step と一致)。
後端面の中央 (y 0.5〜1.0 H) は P 8〜13 kPa で健全。落ちるのは**後端面 × 続く壁 (カウル上面・ランプ) × 内側の縁の 3 重の隅**で、
壁境界層の低運動量流体が縁を回るとき最も強く膨張する所。物理的にはそこに再循環が立って質量が戻るが、局所時間刻みの擬似過渡ではその前に壁 CV が抜ける。

旧メッシャでカウル後縁が厚さ 0 (端面なし)・機体ベースがテーパ (親 plan: 「鉛直 base は node で発散」→ `vehicle_taper`) だったのは、この型を避けていたことになる。
**有限厚の板をそのまま鈍頭で終える形は、node 方式の起動でこの排出に当たる**。

**選択肢** (§5.1 B1d):

| 案 | 内容 | 見込み | 代償 |
| --- | --- | --- | --- |
| (a) 壁節点の Dirichlet 固定を外す (`nodeWallDirichlet: 0`) | 壁 CV の運動量が圧力勾配で応答でき、質量が戻る | **A/B 済 (`run_0428_3d_junction_diag_nowd`, 同じ IC・同じ 1800 step): 完走、後端面は健全 (P 2.9〜25 kPa) = 機序を確定**。ただし縁の壁節点が 1754 m/s で滑り、**入口近くの機体下面 × 側壁外面の凹角 (x 4 mm, z = z_o) で別の排出** (ρ 6.7e-4 → 8.3e-5、P は床 20 Pa に張り付き) | 全壁に効く。壁速度が厳密 0 でなくなる ([[node-wall-velocity-needs-dirichlet-flag]])。生産の格子列は 1 で取っている |
| (b) 後端をテーパにする (板厚を後端の手前で薄く絞る。機体ベースと同じ) | 面が流れに直交しなくなり、壁 CV は接線流を見る | 機体ベースで実績 (t_base 2 mm + テーパ 0.35 で全段通過) | 形状の決定 (どこから・どこまで絞るか)。メッシャは station ごとの $z_o(x)$・$y_{cl}(x)$ で対応可 |
| (c) 起動だけ後端面を slip にする | 壁 CV に接線速度を許す | 未検証 | 段を 1 つ増やす。等温壁との整合 |
| **(e) 段階起動: 弱い壁 (0) で再循環を立ててから Dirichlet (1) に戻す** | 後端面に向かう流れができれば固定 CV も満たされる | **不成立** (`run_0429_3d_junction_diag_wd_after_nowd`: run_0428 の res_1800 から 1 に戻すと **55 step で NaN**、場所はランプ側の隅の後端面節点)。弱い壁の場でもその隅は P 54 Pa・ρ 1.6e-4 で、再循環が立った場ではなかった | — |
| (d) 床を上げる (`roMin`/`pMin`) | 過渡を床で持たせる (run_0428 で P は床 20 Pa に張り付いたまま ρ が抜け続けた = 床は効かない) | 床は保存量の排出を止めない (ρ 3.6e-5 < roMin 1e-4 まで落ちた) | 負密度 → NaN の経路 ([[detectnan-reports-symptom-not-cause]]) |


~~**判断 (2026-09-22, Fable)**: (a)(d)(e) は却下、(c) は (e) と同型なので見送り。**(b) 後端のテーパを推奨する**~~ (**撤回 2026-09-22, codex plan レビュー 3 回目 NO-GO M1–M5**。下の「訂正」を参照) — 根拠は機体ベースの実績 (鉛直 base は node で発散、テーパ + 薄い残りベースで全段通過) と、
本模型の排出が「流れに直交する面 × 縁を回る高速流」の組合せで起きること。テーパなら面が流れに沿い、壁 CV は接線流を見る。
(旧記述: **ユーザ判断待ち** (形状の決定): 側壁後端・カウル後端の板厚を、後端の手前 $L_t$ で $t \to t_{base}$ に絞る。案: $L_t = 5t$、$t_{base} = t/5$ (側壁 5 → 1 mm、カウル 2 → 0.4 mm)、
残るベースの直後は第一 station $t_{base}/5$。メッシャは station ごとの $z_o(x)$・$y_{cl}(x)$ で表せる (位相不変)。**T 字接合部の厚みはそのまま** (R5q の対策は保つ)。
テーパでも残りベースの壁 CV は同じ機序を持つので、$t_{base}$ をどこまで薄くできるかは B1d の A/B で決める。
)

**訂正 (2026-09-22, codex レビュー `2026-09-22-tooling-sern-mesh-blocking-plan-2.md` を全件採用)**:

- 「壁 CV の排出」は**有力な仮説**に格下げ。「機序確定」「段階起動は不成立」「テーパ推奨」は撤回する。
- (M1) 「圧力勾配では戻れない」は実装と不整合: SLAU の質量流束は圧力差項 $-\chi(P_R-P_L)/\hat c$ を持ち、低圧の壁 CV へ質量を押し戻せる。
  ただし $\chi$ を決める Mach 数は接線速度を含む $|U|$ で計算されるので、内点 $(850, 0, 1290)$ m/s では $\hat c \le 1092$ m/s で $\chi = 0$
  (`convectiveFlux_slau_d.inc.cuh:533`)。**検証すべき仮説は「高速の接線流で SLAU の圧力差項が消えている」**。排出率の見積もり 1.5 %/step も
  $u_n/2$ の粗い式で、実装式は $\dot m/A = \rho_w\rho_i/(\rho_w+\rho_i)\,u_n$ + 陰解法の補正。実測に対応する率は 1.24 %/step で「一致」とは言えない。
- (M2) `nodeWallDirichlet: 0` は温度ピン・エネルギー残差ゼロ化・陰解法のエネルギー行固定も外す (`nodeWallDirichlet_d.cu:147`, `timeIntegration_d.cu:1511`)。
  run_0428 が示すのは「壁の強制処理一式を変えると発散経路が変わる」まで。run_0429 は再循環の立っていない場からの復帰なので (e) の反証にならない。
- (M3) IC に不連続が残っていた: `dist_face` が上流側を ∞ に切るので、側壁内面のすぐ内側の流体で $x = L_{sw}$ を跨ぐと速度係数が 1 → 3e-7 に跳ぶ。
  混合重み $w_z, w_y$ も $x > L - \delta$ で二値から混合へ突然切り替わる (流体内で 1 → 0.85 / 0.65)。「速度不連続は原因でなかった」は維持できない。**修正済** (x 方向にも smoothstep)。
- (M4) TP の EOS 床は `ro = max(ro, roMin)` を**保存量に書き戻す** (`dependentVariables_d.cu:76,203`) が、その後の温度ピンは床なしの $P = \rho R T_w$ を書く (`nodeWallDirichlet_d.cu:80`)。
  出力が床未満でも「床が効いていない」証拠にならない (床補正 → 排出の繰り返しと両立)。局所 dt は受け側の音速だけでなく面速度・両 CV の V/A・粘性項で決まる (`setDT_d.cu:64,181`)。
  run_0425 の「受け側の流入過大で ρ < 0」も、正味流入なら密度は増えるので流束の符号と陰解法補正の確認が要る (事実として残すのは「step 1 で ρ −0.28、隣が ρ max 0.17」まで)。
  閉性検査 FAIL (max 1.7e-4, 6,410 CV) をドライバが無視して続行していた → 警告を出すよう修正、解釈時に考慮する。
- (M5) テーパの A/B で $\Delta x_1 = t_{base}/5$ とすると CV 体積・局所 dt・解像度が同時に変わり、形状効果と格子効果を分離できない。残るベース面は依然 x に垂直。
  試すなら鈍頭対照を残し $t_{base}$ と $\Delta x_1$ を独立に振る。合否は「完走」でなく、全壁タグの床到達数・固定 CV 群の密度履歴・ベース周囲の符号付き質量収支・逆流域・壁圧/熱流束・床補正量 +
  `check_convergence.py` (同一設定区間) と `check_quasisteady.py` の VERDICT。「連続式だけ弱形式に」は連続式が既に残っているので変更内容が定まらない。
- (m6) 比較対象は node ID・壁タグ・時刻で固定し、AWS の commit・バイナリ・実効 config・IC 生成条件を診断記録に添える。

**次の手順 (codex 推奨、鈍頭のまま)**: ① IC の不連続を除去 (済) → ② 問題 CV の毎 step (実際は出力刻み) の面別流束・更新・床補正を記録 → ③ 運動量固定と熱的固定を分離した A/B (断熱壁 = 運動量固定のみ)。

#### 4.13.1 ①②③ の実施 (2026-09-22, run_0430 / run_0431)

codex 3 回目の指示どおり **鈍頭のまま** ① IC 修正 → ② 面別 SLAU 流束 → ③ 断熱壁 A/B を実施した。
以下は**観測事実**と**解釈 (仮説)** を分けて書く (前回「機序確定」と書いて撤回になったため)。

**共通の再現条件**: AWS `~/forge-sern`、`feature/sern-design` HEAD `5964312c`。
**forge バイナリは `c21fdd19` ビルドのものを流用** (sha256 `eec7f767…`)。run_0427 と同一バイナリにして IC 以外を揃えるため
(`c21fdd19`→`5964312c` の solver 差分は `cuda_forge/transition_d.cu` のみ、本 run は `turbulence: none` で不活性)。
メッシュ `cad/demo_out/jm_c16.msh` (2,446,572 節点)。`check_mesh_quality.py` **PASS** (AR max 676, skew max 0.533)、
`check_dual_closure.py` は **FAIL** (max 1.7e-4, 6,410 CV。B1b の既知欠陥)。
`solverConfig.yaml` は **run_0427 とバイト一致** (node / `nodeWallDirichlet: 1` / SLAU / `convMethod: 0` / blockDPLUR /
`unsteady: 0` / `cfl = cfl_pseudo = 0.2` / `pMin 20` / frozen_tp / `turbulence: none` / 1800 step / 200 step 出力)。

| run (AWS `case/46.sern_design/`) | 差分 | 結果 |
| --- | --- | --- |
| `run_0430_3d_junction_diag_icfix` | **IC の不連続を修正** (`13ea76ca`)。数値設定は run_0427 と同一 | 1800 step 完走・NaN 0・床到達 0。ただし**排出は続いている** (下) |
| `run_0431_3d_junction_diag_adiabatic` | run_0430 の bcond だけ `wall_isothermal`(Ts 1000 K) → `wall` (断熱)。他は同一 | 1800 step 完走。**同じ節点が同じ率で抜ける** |

**事実 1 (排出は止まっていない)**: run_0430 の ρ_min は単調に下がり続ける。
`res_1200` 5.83e-4 (x/H 1.6) → `res_1400` 3.99e-4 (1.6) → `res_1600` 2.86e-4 (1.2) → `res_1800` **1.71e-4** (1.2)。
側壁後端面の壁節点 CV 153797 では 200 step ごとにほぼ半減し、相対排出率 $(\dot\rho/\rho)$ は
`res_1200` 以降 −4.5e5 → −6.7e5 → −7.9e5 → −7.4e5 /s と**ほぼ一定 = 指数減衰**。
→ **「1800 step 完走・床到達 0」は打ち切りが早かっただけ**で、この率のまま外挿すると
床 `roMin` 1e-4 への到達は **step ≈ 1950**。完走を合格の根拠にしない (codex M5)。
**IC の修正は発散時刻を動かしたが、排出そのものは消していない** (run_0427 は同設定・同バイナリで step 1720 NaN)。

**事実 2 (面別 SLAU 収支。B1d ②)**: `case/46.sern_design/cad/diag_wall_cv_budget.py` で
`convectiveFlux_slau_d.inc.cuh:533-564` の mdot を 1 次 (`convMethod 0`) で場から再計算した
(`slauContactFloor` 既定 0 なので欠落項なし。境界半割面は `convectiveFlux_boundary_d.inc.cuh:114,135,192` で
$\dot m = A\rho_R(\mathbf U_b\!\cdot\!\mathbf n)$、no-slip 壁は $\mathbf U_b = 0$ なので厳密 0)。
CV 153797 (node (x,y,z)/H = (1.2, 1.2675, 1.0046)、$V$ = 1.1723e-11 m³) の `res_1800`、$\dot m>0$ = 流出:

| face | 相手 | $A$ [m²] | $\mathbf n$ (外向き) | $\dot m$ | 移流項 | 圧力差項 | $\chi$ | $\hat c$ | $V_{n,\text{self}}$ | $V_{n,\text{oth}}$ | $P_{\text{self}}$ | $P_{\text{oth}}$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 575140 | 1135684 (内点) | 9.381e-08 | (1, 0, 0) | **+4.345e-09** | +4.345e-09 | **0** | **0.0000** | 700.5 | 0.0 | **276.2** | 57.8 | **3545** |
| 574891 | 153714 (壁) | 1.350e-07 | (0,0,−1) | −1.095e-09 | 0 | −1.095e-09 | 1.0000 | 665.1 | 0.0 | 0.0 | 57.8 | 68.6 |
| 575139 | 153880 (壁) | 1.350e-07 | (0,0,1) | −1.746e-09 | 0 | −1.746e-09 | 1.0000 | 665.1 | 0.0 | 0.0 | 57.8 | 75.1 |
| 575135 / 575138 | 153796 / 153798 (壁) | 1.117e-08 | ±(0.24,∓0.97,0) | −1.45e-11 / −2.45e-12 | 0 | 同左 | 1.0000 | 665.1 | 0.0 | 0.0 | 57.8 | 59.6 / 58.1 |
| 合計 | | | | **+1.488e-09 kg/s** | | | | | | | | |

$(\dot\rho/\rho) = -7.4\times10^5$ /s。dump 間で 200 step あたり約 1/2 という実測から逆算した実効 $\Delta t$ は 4.7e-9 s で、
記録済みの局所 dt 4.6e-9 s と **2 %** で整合する。
式の手計算も一致する: $L$ = 壁節点 ($\rho$ 1.711e-4, $V_n=0$)、$R$ = 内点 ($\rho$ 8.436e-3, $V_n$ = +276.2)、$g=0$、
$\widehat{|V_n|} = \rho_R\,276.2/(\rho_L+\rho_R) = 270.7$、移流項 $= \tfrac12(\rho_L(0+270.7) + \rho_R(276.2-270.7)) = 0.0464$ kg/m²s、$\times A$ = 4.35e-9 ✓。

$\chi = (1-\hat M)^2$、$\hat M = \min(1, \sqrt{(|U_L|^2+|U_R|^2)/2}/\hat c)$ は**速度ベクトルの大きさ**で決まる (`:538`) ので、
内点の $|U| \ge \sqrt2\,\hat c$ = 991 m/s のこの面では $\hat M = 1$、**圧力差項 $-\chi(P_R-P_L)/\hat c$ が恒等的に 0**。
$P$ が 3545 Pa 対 57.8 Pa と 61 倍あってもこの面からは質量が戻らない。
**これは Shima–Kitamura の SLAU 原式どおりであって実装の欠陥ではない** (「χ が消えている」をバグのように書かないこと)。

**事実 3 (釣り合う壁節点もある)**: カウル後端面の CV 189814 (x/H 1.6) は同じ指紋 ($V_{n,\text{oth}}$ = 381.5, $\chi$ = 0、流出 +4.449e-08) を持つが、
壁↔壁面 713819 ($A$ 4.166e-07 > 流出面 $A$ 2.267e-07、$\chi=1$、$P$ 335 対 180 Pa) から −4.884e-08 の流入があり、
正味 +7.7e-10 でほぼ均衡する。実際 ρ は `res_1400` 3.99e-4 で**反転**し `res_1800` 5.39e-4 まで戻った。
→ 「$V_{n,\text{oth}}>0$ かつ $\chi=0$」は**必要条件であって十分条件ではない**。

**事実 4 (断熱壁 A/B。B1d ③)**: run_0431 では壁節点の音速が 665 → 937 m/s と明確に変わり ($T$ が自由に動いた)、
それでも CV 153797 の ρ は `res_1800` で **1.836e-4** (等温 1.711e-4)、相対排出率 −6.9e5 /s (等温 −7.4e5 /s)、
面別の内訳も同型 (流出面 $\chi$ = 0.0000 で移流項のみ +5.105e-09、壁↔壁面は $\chi$ = 1 で流入)。
ρ_min の場所も dump ごとに等温と同一。→ **排出は運動量拘束の側であって熱的拘束ではない** (codex M2 の交絡を分離できた)。

**事実 5 (閉性は効いていない)**: 当該 CV の $|\Sigma \mathbf S|/\Sigma|\mathbf S|$ は
153797 **4.19e-6** / 153796 3.94e-6 / 153798 4.94e-6 / 189814 **6.35e-8** / 1135684 3.27e-9 (全域 max 1.705e-4)。
収支は実際の面ベクトルで計算して実測減衰率と 2 % で合うので、閉性誤差が支配するには当該 CV で $\gtrsim 0.3$ が要る。
→ **B1b (変換器の単精度) はこの排出の原因ではない** (B1b 自体は別途直す)。

**収束ゲート** (AGENTS.md 必須。どちらも合格ではない):

```
run_0430: NOT CONVERGED (stalled/plateau) rms_ro 0.6dec flat / rms_roUx 0.9dec flat / rms_roUy 1.1dec flat /
          rms_roUz 0.6dec falling / rms_roe 1.0dec falling / rms_roY1 0.3dec flat
run_0431: NOT CONVERGED (stalled/plateau) 全列 flat、rms_roY1 は peak 1.15e-5 > fin 8.58e-6
```

**仮説 H1 (機序ではなく仮説。確度: 高 = `res_1800` の離散収支として / 中 = 「床 → NaN に至る」「起動で解消しない」として)**:
壁節点 ($u$ ピン) の CV は、W↔I 面の隣接内点が壁から遠ざかり ($V_{n,i}>0$) かつ $|U_i| \ge \sqrt2\,\hat c$ ($\chi=0$) のとき
**移流の流入経路を持たず**、壁列の $\chi=1$ 面 (両側 $V_n=0$ なので圧力差項だけが残る) からの圧力拡散だけが補充源になる。
**壁節点 1 個でなく壁列 (row) の収支**で決まり、列のどこかに $V_{n,i}\le0$ か $\chi>0$ の節点があれば列全体が補充され、
無ければ列ごと抜ける。側壁後端 × カウル上面 × 内側の縁の 3 重隅では補充が足りず、正味の時定数
$\tau_{\text{net}} = \rho/|\dot\rho| \approx 1.35$ µs (局所擬似時間) の指数減衰になる
(単面の排出だけなら 2.2e6 /s だが、壁↔壁の流入で正味 7.4e5 /s。**数字は正味で書く**)。

- **反証条件**: 延長 run で内点 1135684 の $V_{n,i}$ が減衰し ρ_w が床の手前で反転する (= 起動の競争であって構造的でない)。
- **まだ言えないこと**: 「床 → NaN」の経路。EOS 床は `ro` を保存量に書き戻す (`dependentVariables_d.cu:207,271,282,296`) ので、
  §4.13 に記録した run_0427 の「ρ 3.0e-5 < roMin」は **1 step で 70 % 落ちた**ことになり、0.35 %/step の排出では説明できない。
  床以降は別の機構 (陰解法の対角が床密度で崩れる等) の可能性があり、5 step 刻みで見ないと分からない。
- **第 3 の可能性 (未確認)**: 内点 1135684 の状態 ($P$ 3545 Pa, $|U|\ge991$) は内側の縁を回る Prandtl–Meyer 膨張の**過渡**で、
  後流の再循環が立つ前の姿かもしれない。隅の粘性層は粗い見積りで $\delta \approx \sqrt{\nu s/U} \approx 50$ µm、
  第一 station $\Delta x_1$ = 250 µm はその 5 倍外にある。**`t_base/5` 則はここでは尺度が違う** ([[base-wake-resolution-rule]] は候補に戻す)。

**判別 A/B (`run_0432_3d_junction_diag_ext4k`)**: run_0430 の `res_1800` から `restart_by_index` で継続、
**config は `nStepOuter` 1800→4000 以外 run_0430 と一致**、加えて 4 点の**毎 step プローブ**
(CV 153797 壁 / 1135684 内点 / 189814 壁 / 1785247 内点、`point_probe_[0-3].out` に `T,P,Ux,Uy,Uz`)。

- **分岐 A** (ρ_w が床の手前で反転し $V_{n,i}$ が減少): **起動の競争**であって構造的でない。対策は起動側
  (暖機の刻み方・run_0428 の場で再循環が立つまで待ってから `nodeWallDirichlet` を 1 に戻す等)。**幾何は触らない**。
- **分岐 B** ($V_{n,i}\approx276$ のまま床到達): このメッシュ + スキームで**構造的**。次の二択を起票する —
  (i) 隅の $\Delta x_1$ を**粘性層の尺度**に振る ($t_{base}$ 固定。codex M5 の「独立に振れ」)、
  (ii) ソルバ側: 壁隣接面の $\chi$ を**面法線 Mach** で評価する ($M_n$ = 0.39 → $\chi$ = 0.52 → 補充 ≈ 1.2e-7 ≫ 排出 4.3e-9)。
  **(ii) は `cuda_forge/` の数値カーネル変更なので、着手前に plan 起票と `diagnostician` (条件 6) を通す。**
  NaN に至る step は指標にしない。

**やらないと決めたこと (2026-09-22, `diagnostician` 判断)**: テーパ実装 (3 重隅が残るので codex M5 と同じ結論で候補外) /
`slauVariant: 2` (mdot の $\chi$ 項は SLAU/SLAU2 共通で無意味) / 床を上げる / `nodeWallDirichlet` をタグ別に効かせる
(`wall_flag` が bcond 横断なので共有節点にピンが残る) / A/B 前に $\chi$ の定義をいじる / `check_quasisteady` で STEADY を取りに行く
(減衰中は `DRIFTING` が正しい答え)。

**判別 A/B の結果 (2026-09-22, `run_0432_3d_junction_diag_ext4k`) → 分岐 B**:
run_0430 の `res_1800` から `restart_by_index` で継続。**config は `nStepOuter` 1800→4000 以外 run_0430 と一致** (diff で確認)。
加えて 4 点の**毎 step プローブ** (`point_probe_[0-3].out`, `T,P,Ux,Uy,Uz`。CV 153797 壁 / 1135684 内点 / 189814 壁 / 1785247 内点。
probe は `centCoords` 最近傍で、起動ログに与えた座標がそのまま出るので取り違えなし)。**step 815 (通算 2615) で `ro` に NaN**。

- **事実 6 (駆動条件は緩まない = 分岐 B)**: 内点 1135684 の `Ux` は 276.3 → **286.6** m/s (微増)、`|U|` 1316 → 1303、
  面 575140 の $\chi$ は全 dump で 0.0000。815 step (通算 2615 step) かけて**減衰しない**。
  → 「起動の競争」ではない。**起動側の工夫 (暖機の刻み方) でこの排出は消えない**。
- **事実 7 (壁列を伝う)**: CV 153797 の面別内訳の時間変化 ($\dot m>0$ = 流出):

  | dump | 575140 (内点, $\chi$=0) | 575139 (壁 153880, $\chi$=1) | 574891 (壁 153714, $\chi$=1) | 正味 |
  | --- | --- | --- | --- | --- |
  | res_0 | +4.346e-09 | −1.746e-09 | −1.095e-09 | +1.488e-09 (流出) |
  | res_200 | +2.882e-09 | −4.852e-10 | −1.728e-09 | +6.591e-10 (流出) |
  | res_400 | +2.143e-09 | −1.385e-10 | −2.109e-09 | −1.184e-10 |
  | res_600 | +1.365e-09 | **+2.304e-10** | −2.674e-09 | −1.099e-09 |
  | res_800 | +7.105e-10 | **+7.825e-10** | −3.099e-09 | −1.621e-09 |

  補充源だった $\chi$=1 面が、隣の壁節点 153880 が先に抜けた結果 (res_800 で $\rho$ 3.735e-06 < self 2.654e-05)
  **res_600 以降は流出面に転じる**。H1 の「壁列 (row) の収支で決まる」定式化と整合する。

**事実 8 (EOS 床は保存量に入らない。定常陰解法の commit が床前の基準から足すため)** — コードで確認:

1. `dependentVariables_d.cu:76` で `ro_temp = max(ro[ic], roMin)`、`:207` で `ro[ic] = ro_temp` と**作業配列には書き戻す**
   (流束・`setDT`・block-DPLUR の Jacobian はこの床後の状態で組まれる)。
2. しかし定常の commit は `update_d.cu:248` の **`ro[ic] = roN[ic] + d0`** で、`roN` = **床前**の基準値。
   → 床の分は捨てられ、保存量は床を割って減り続ける。**dump/probe に出る $\rho$ が `roMin` を下回るのはこのため**
   (codex M4「床は保存量に書き戻す」と §4.13「床は保存量の排出を止めない」は、**書き戻す先が作業配列**であることで両立する)。
3. 帰結: $\rho_{raw} < \rho_{Min}$ の節点では流束が $\rho_f = \rho_{Min}$ (定数) で評価されるので、
   排出が $\rho_{raw}$ に比例しない**一定のシンク**になる。**指数減衰 (0 に漸近するだけ) が線形減衰 (有限 step で 0 を切る) に変わる**。

**事実 8 の定量検証 (パラメータ・$\Delta t$ に依存しない形)**: `point_probe_0.out` の $\rho = P/(R T_w)$ ($R$ = 338.3, $T_w$ = 1000 K) を
区間に分けて調べた。予測される折れ点は $\rho_{raw} = \rho_{Min}$ ⇔ **$P$ = 33.83 Pa** と $P = p_{Min}$ = **20 Pa**。

| 区間 | step | 形 (決定係数) | 傾き |
| --- | --- | --- | --- |
| A: $\rho_{raw} > \rho_{Min}$ | 1–275 | **指数** ($R^2$ 0.9888 対 線形 0.9700) | 率 −1.93e-3 /step |
| B: $\rho_{Min} > \rho_{raw}$, $P>p_{Min}$ | 277–549 | **線形** ($R^2$ **0.999967** 対 指数 0.9948) | −1.504e-7 /step |
| C: $P < p_{Min}$ | 551–814 | ほぼ線形 (線形 0.9909 / 指数 0.9990) | −1.308e-7 /step |

$P$ が 33.83 Pa を切るのは step **276**、20 Pa を切るのは step **550** で、**形の変化は予測どおりの場所で起きている**。
さらに、床直前 (step 200–275) の局所指数率 $C$ = 1.399e-3 /step から予測される**床後の一定シンクの傾き** $-C\rho_{Min}$ = **−1.399e-7 /step** は、
実測の −1.504e-7 (B) / −1.308e-7 (C) と **予測/実測 = 0.93 / 1.07** で一致する ($\Delta t$ も自由パラメータも使わない)。
→ **床は「排出を止めない」だけでなく、自己抑制的な指数減衰を一定シンクに変えて負密度への到達を保証している**。

**付随して見つかったこと (別欠陥)**: dual-time 経路の commit は `update_d.cu` の `applyBlockImplicitCorrectionInPlace_d`
(`ro += dq`) なので、**床の意味が定常と非定常で違う**。また既存の opt-in `time.deltaT.updateGuardAlpha`
(`solverConfig.cpp:337`, 既定 0) が commit 時に `dq` を正値性でスケールする。§5.1 **B1e** に起票した
(既定の変更は全ケースに効くので、まず `updateGuardAlpha` で足りるかを確認する)。

**まだ書かないこと**: 「機序確定」。事実 2・6・7 は**離散収支として確認済み** ($\rho \ge \rho_{Min}$ の区間)、
事実 8 は**床の実装事実 + 定量一致**だが、**dump の収支がソルバの収支と一致するのは $\rho\ge\rho_{Min}$ かつ $\rho R T_w \ge p_{Min}$ の間だけ**
(それ以降は「commit 後・pin 後」の dump 状態と「pin → 床・クランプ後」のソルバ状態が別物になる。
`diag_wall_cv_budget.py` にソルバ前処理を足した再収支は**未実施**)。

**未了 (AWS 停止のため)**: `res_nan_815.h5` で最初に負・NaN になった**節点 ID** の確認。一定シンクなら先に 0 を切るのは
$\rho_{raw}$ が最小の **153880** (res_800 で 3.735e-6) であって 153797 ではないはず。**153797 を NaN の出所と書かないこと**。

**次手 (2026-09-22 `diagnostician` 判断。(i)(ii) の前に config だけで済む判別を置く)**:
層流暖機の `solver` を **SLAU → ROE** にした同一 restart (`res_1800` 起点、他は同一)。
Roe の質量流束は $|u\pm c|\,\Delta p/c^2$ の音響波散逸を**全 Mach で**持ち、SLAU の $\chi$ による切替が無い。

- **A** ($\rho_w$ の減衰が止まる/大幅に鈍る。**床前区間**で判定): 「壁隣接面で圧力–質量結合が $\chi$=0 で消える」が本質
  → **(ii) 壁隣接面の $\chi$ を面法線 Mach で評価する**を plan §4 に起票する根拠になる ($M_n$ = 0.39 → $\chi$ = 0.52 → 補充 ≈1.2e-7 ≫ 排出 4.3e-9)。
  **(ii) は `cuda_forge/` の流束カーネル変更なので、plan §4 起票 + codex plan 段 + 条件 6 が先** (実装から入らない)。
- **B** (同じ率で抜ける): 圧力項の有無でなく「$u$ ピン節点に移流の流入経路が無い」構造が本質
  → **(i) 隅の $\Delta x_1$ を粘性層の尺度に振る** ($t_{base}$ 固定) か、Dirichlet 側の設計を起票。

判定は**床前区間の $\rho_w$ 減衰率 (連続量)** で行い、**NaN の step は指標にしない**。Roe は TP 2 種に対応 (`convectiveFlux_roe_d.inc.cuh:159,192,242`)、HLLE は TP 未対応。

**Roe の判別 A/B (2026-09-23, `run_0433_3d_junction_diag_roe`) → 分岐 A (狭い主張に限って確定)**:
run_0430 の `res_1800` 起点 (run_0432 と同一)、**`solver: "SLAU"`→`"ROE"` のみ変更** (+ `outStepInterval` 200→1000、数値に無関係)。
probe は run_0432 と同一 4 点・毎 step。**step 15 で NaN, rc=1**。

| 節点 | SLAU (run_0432) | Roe (run_0433) |
| --- | --- | --- |
| CV 153797 (側壁後端面 壁) | 1.711e-4 から **−0.19 %/step** | 1.951e-4 → 8.114e-4 (**+7.5 %/step**) |
| CV 189814 (カウル後端面 壁) | 減衰 | 5.738e-4 → 1.331e-3 |
| CV 1135684 (内点) | `Ux` 276.3→286.6 | `Ux` 276.3→280.1 (ほぼ同じ) |

**step 1 で既に符号が反転している** (1.711e-4 → 1.951e-4, +14 %)。大きさも Roe の音響散逸項の手計算と整合する:
面 575140 ($\Delta P$=+3487 Pa, $\Delta u_n$=276, $\tilde\rho\approx$1.2e-3, $\tilde c\approx$700) で
$\alpha_\pm\approx$3.8e-3/3.3e-3, $|\lambda_\pm|\approx$941/459 → 散逸の質量成分 ≈ −2.55 kg/m²s、中心項 +1.16 →
**正味 −1.39 kg/m²s = 壁へ 1.3e-7 kg/s の流入** (SLAU の排出 4.3e-9 の約 30 倍)。陽的なら +30 %/step 相当で、
観測の +14 %→+7.5 % は陰解法の減衰込みで整合する。→ **過補正ではなく Roe の項がそのまま出ている**。

- **確定してよい (分岐 A)**: 「**$\chi=0$ が壁 CV への補充を消している**」。
  根拠: step 1 の符号反転 + 大きさの一致 + 同起点の SLAU が −0.19 %/step。
- **確定してはいけない**: 「圧力–質量結合があれば壁 CV は**定常に**満たされる」「振動しない」。
  到達点は $P_w \to P_i$ (≈3.5 kPa, $\rho_w\approx$1e-2) で、+7.5 %/step なら 1.7e-4 から約 55 step 要る。
  **15 step は 1/4 の道のりで飽和は未観測**。これは対策 plan の §6 で測る。
- **Roe は診断専用**。run_0433 の NaN は **x/H 2.101–2.400, y/H −1.000〜−0.565, z/H 1.049–1.173 の外部流膨張域**で
  144 節点、step 14 に `rms_roY1` 3.0e+23 / `rms_roe` inf と**化学種から**爆発した (後端面は健全、153797 は 8.117e-4)。
  Roe の**低密度非正値性** (強膨張で中間状態が負になりうる。AUSM⁺ と違い正値保存でない) という既知の性質で、**本件と独立**。
  生産・段階起動には使わない。

**NaN の順序 (B1d の未了項目) — 測れないが推定で閉じる**:

- `run_0434_3d_junction_diag_nanorder` (run_0432 の `res_800` から 60 step, 5 step 出力, 5 点毎 step probe) は
  **方法として不成立**。restart した時点で床下の $\rho$ が `roMin` へ戻り (153797 2.654e-5 → 9.962e-5、
  **153880 3.735e-6 → 9.953e-5**)、元は 15 step 後に落ちたところが 60 step 完走した。
- **理由 (コードで確認)**: 初期化は `dependentVariables` (`main.cpp:1225`, 床を作業配列へ) → `updateVariablesOuter`
  (`main.cpp:1252`, `roN = ro`) の順なので**床後の値が基準に入る**。ループ中は `updateVariablesOuter` が commit の後
  (`main.cpp:1805`) なので入らない。→ **床下の節点では「restart ≠ 継続」**。§5.1 B1e に追記。
- **推定**: 153880 は `res_800` で $\rho$ 3.735e-6、床後の一定シンク −1.5e-7 (run_0432 区間 B) 〜 −4.5e-7/step (run_0434) で
  **0 を切るのは step 808–825**。NaN の step 815 はこの窓の中で、153797 (残り約 190 step) より先。
  207 節点が一斉に NaN になるのは block-DPLUR の `nStepInner 5` sweep が 1 step で dq の NaN を 5 リング先まで運ぶため
  (壁近傍 5 リングの一部 ≈ 200 節点)。**順序は `res_nan_*.h5` からは原理的に読めない**。
- 厳密に測る唯一の方法は run_0432 を **`res_200`** ($\rho$ 1.11e-4 > `roMin` でリセットされない) から `outStepInterval 5` で 615 step
  回すこと。**結果が推定と違っても対策は変わらない**ので、価値が低いと判断して実施しない。

**対策の起票**: [`convection-slau-wall-normal-chi.md`](convection-slau-wall-normal-chi.md) (2026-09-23)。
**文献調査で先行研究が見つかったため本命は再検討中** (SLAU の $|\mathbf u|$ ベース $\widehat M$ は原典の意図的な選択だが経験則。
AUSM 系には超音速での圧力流束散逸不足を補う改良と、質量流束の pressure diffusion 項を作り直して正値性を保つ改良がある。
AUSM⁺ は CFL 条件下で密度正値保存性が示されているが Roe には無い)。案は A 自作 $\chi_n$ / B 壁隣接面だけ ES 行列散逸 (KEEP のカーネル流用) /
C 既存の改良 AUSM / D 正値性ガード (単独では不可)。**codex plan 段が先**。

**一般性 (未確定)**: この指紋は node + `nodeWallDirichlet: 1` + 圧力拡散を $|U|$ で切る流束 (SLAU/SLAU2, AUSM⁺-up も同型) の
組合せで、**超音速の剥離縁で第一内点が壁面粘性層の外にある**ときに出る。接続模型固有とは限らないが、
B3 (メッシュ方式の受入) の中で先に A/B を終え、分岐 B の結果が出てからソルバ側の plan を起こす。


## 5. 実装ステップ

1. 実証模型のジェネレータ (N の C リング + スロット後流 + 露出カウル側端 + カウル後流 + U/S の H ブロック)。
2. 受入試験スクリプト (§6)。3 解像度 × dv の厳しい組合せで無人生成。
3. 結果を codex result レビューに掛け、方式を確定 (または再検討)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| ~~S1~~ (済: スクリプト 3 本を `case/46.sern_design/cad/` に配置。S2–S5 は保険として保留) | **Salome パイプラインを repo に入れる (§4.8/§4.9)** | scratchpad のスパイクを `case/46.sern_design/cad/` (GEOM+SMESH スクリプト) と 汎用の `med_to_msh41` (必須グループを引数化) に整理。隅ごとの半径プロファイル r(x) を JSON/YAML 入力にする |
| **S2** | 受入試験 (非構造版) | **tet 数 > 0・層厚未達の警告なし**・全外部面がタグに属する・`check_mesh_quality` の VERDICT・壁別の第一層 (B0c の点別検査を流用)・node 変換後の双対閉性。`Compute: True` を合格の根拠にしない |
| **S3** | 残課題の調査 | skew>0.90 の 46 セルと壁第一層最大 6.7e-04 の位置特定 (徐変始端・消える終端の微小面)。最小半径 r_min の下限を決める |
| **S4** | SERN 全体形状への展開と規模見積もり | ランプ輪郭 (MOC 点列) を GEOM のスプライン面で持つ。節点数・変換 RAM・AWS での生成 (Salome は AWS 未導入) |
| **S5** | 粘性層の制約の緩和策 | 板厚で総厚が縛られる件。板を厚くする/端を丸める/板の周囲だけ層を薄くする (面別 VL) の比較 |
| ~~B0~~ (済 2026-09-22: §4.10。codex plan レビュー 2 回目 GO-with-changes の全件を §4.10/§4.5/§6 に反映) | **接続設計を完成させる (B1 の前提。codex plan M1/M2)** | 全ブロックの**頂点・曲線・6 面の対応表**と、`L_sw`・`L_cowl` の前後の**全共有面一覧** (同一 surface ID・節点 ID・分布・向きまで規定)。`z = W/2` の分割をカウル後流と U に伝播させる。**フィレット終端は「消える円弧をブロック辺にしない」専用の 3D 遷移として設計** (`addSurfaceFilling` は境界 3–4 曲線のみ・半径 0 の円弧は不可)。近似曲面を使うなら U3 の形状誤差・接線誤差の許容値を先に決める |
| ~~B0b~~ (済 2026-09-22: §4.11 の 1・2・7。第一層は η = h1/dr で揃える方式に変更) | **壁層の配置を 3D 法線で定義する (M3)** | 壁面ごとに 3D 法線に対する目標距離を定義 (傾斜 22° のランプで y 間隔 $h_1$ は法線距離 0.927 $h_1$ = −7.3 %)。`sidewall_end`・`cowl_base` 直後の x 分布を追加。リングは一様厚み (§4.7)。全壁一律の $h_1$ が不要なら壁別の目標値を明示 |
| **B0c** (一部済 2026-09-22: 頂点 Jacobian・位相/タグ・壁第一層は `hex_junction_model.py` に、CV ごとの閉性は `check_dual_closure.py` に実装。`check_mesh_quality.py` の体積判定は §4.12 で修正済。**残: 隣接間隔比の実測検査** (x は設計値で検査済、断面内は未)) | **実行可能な受入ゲートを作る (M5/M6)** | (i) 物理座標を float32 化した後の**局所 (頂点) Jacobian** 検査を独立実装 (`check_mesh_quality.py` は四面体分割の体積和しか見ず、最小頂点 Jacobian −0.175 の反例を通す)。(ii) 双対閉性は**各 CV の $|\Sigma S|/\Sigma|S|$** で判定 (変換器ログは全域最大面積で正規化、双対体積は絶対値加算なので反転を排除できない)。(iii) 壁層は**壁別・点別**に指定値との誤差・評価不能数・選択点が内部点であること・層数を検査 (`check_wall_resolution.py --geometry-only` は統計を出すだけで `MEASURED` を返し、別の壁上の点を除外しない)。**`SOFT-PASS` は拒否**、AR 5000 の緩和は近直交の壁法線層に限定 (45° の隅セルは対象外) |
| ~~B1~~ (済 2026-09-22: §4.11。scale 0.5・0.71 で PASS、scale 1.0 は手元のメモリに載らず未実行) | **実証模型「側壁終端を含む 3D 接続模型」** | §4.1–§4.4 を 1 つの小規模模型に入れる。壁法線は第一層 4 µm・約 40 層。メッシュのみ (CFD なし) |
| B2 | 受入試験 | §6 の 4 群を B0c のゲートで判定。**方式受入は実際の `scale` 3 段階**で行い、第一層 4 倍刻みの系列は別の感度試験にする (M7)。dv の具体値・許容範囲・拒否境界を表にし、**短いフィレットテーパ・短い露出カウル区間 (`L_cowl − L_sw` 下限近傍)** を含める。入口面積は「CAD 形状・端点の不変」と「離散境界の面積誤差 (弦近似。$r_f$ 0.10 H の四分円 8/16/32 分割で 1.0e-4 / 2.5e-5 / 6.3e-6 H²)」を**別ゲート**にする (M4) |
| **B1b** | **変換器の双対面を倍精度で組む** | CV ごとの閉性が max 3.9e-4 (§4.11)。座標を double で読み、面ベクトル・重心・体積を double で計算して float32 で書く。別セッションの変換器作業と衝突しないよう着手前に確認 |
| ~~B1c~~ (済 2026-09-22: §4.12 ユーザ決定。第一層 16 µm・端面 0.25 mm で 245 万節点、全ゲート PASS) | **規模を決める (U4)** | scale 1.0 は 1250 万節点・変換 46 GB の見積もり。削る候補ごとに節点数を実測し、ユーザと決める |
| **B1d** (**F**) | **接続模型の CFD 成立確認 (§4.13 / §4.13.1)** | ①②③ + 判別 A/B 済 (run_0430–0434)。**判定: 分岐 A** — 「$\chi=0$ が壁 CV への補充を消している」まで確定 (run_0433 の Roe で step 1 から符号反転、大きさも Roe の音響項の手計算と一致)。**「結合が戻れば定常に満たされる」は未確定** (15 step は飽和点の 1/4)。**判断: 2026-09-23 `diagnostician` — 対策は [`convection-slau-wall-normal-chi.md`](convection-slau-wall-normal-chi.md) へ切り出し、本 plan からは実装しない**。NaN 順序は**測れないが推定で閉じた** (153880 が step 808–825 で先。§4.13.1)。**残**: 対策 plan の V1 が通ったら本 plan の B3 (方式確定) へ戻る |
| **B1e** (**F**) | **EOS 床が保存量に入らない (定常陰解法。§4.13.1 事実 8)** | `dependentVariables_d.cu:76,207` は床を**作業配列**に書き戻すが、定常の commit は `update_d.cu:248` の `ro = roN + d0` で**床前の `roN` 起点**なので床が捨てられる。結果、$\rho_{raw}<\rho_{Min}$ では流束が $\rho_f=\rho_{Min}$ 固定で評価され、**指数減衰が一定シンクに変わり有限 step で負密度に達する** (定量検証済み: 予測/実測 0.93–1.07)。dual-time は `applyBlockImplicitCorrectionInPlace_d` の `ro += dq` で**床の意味が定常と違う**。**再起動時だけは床が基準に入る** (`main.cpp:1225` → `:1252` の順。ループ中は `:1805` が commit の後) ので、**床下の節点では「restart ≠ 継続」** (run_0434 で実測: 153880 3.735e-6 → 9.953e-5)。**既定の変更は全ケースに効くので、まず既存 opt-in `time.deltaT.updateGuardAlpha` (`solverConfig.cpp:337`, 既定 0) で足りるかを確認する**。接続模型の (i)(ii) とは独立。着手は plan §4 起票 + 条件 6 + codex plan 段を経てから |
| **B1f** (済 2026-09-23) | **出口 BC を `outflow` に (同じ罠を 2 回踏んだ)** | `outlet_statPress` は node の壁列・後流の**亜音速ノード**に Ps ≪ 実出口圧を課し、そこから圧力が育つ。**run_0121 で一度踏み** (対策が `problem_r5_3d_sst_outflow.yaml` の個別 YAML に留まった)、**接続模型で再発** (run_0430–0436: 出口の亜音速率 3.5 → 99.5 %、`far_bottom` 22 MPa、残差 +0.8 桁)。**A/B**: 同一起点・同一設定で出口と `far_bottom` を `outflow` にしただけの `run_0437_3d_junction_chi_n_outflow` は 亜音速率 4.4 % / 逆流 0.1 % / 全域 P>2e5 が **0 節点** / 残差 **−1.3 桁**。**処置**: `run_junction_model.py`・`runner_sern3d.py`・`runner_sern.py` の `outlet_kind` 既定を `outflow` に変更 (2D は切替を新設)。`procedures/recommended-settings.md` の「出口」に追記。**影響は case/46 のみ** (`PHYS_SERN` 専用・他 runner は自前の bcond)。**残**: この既定変更より前の生産値 (格子収束列 run_0418–0421 等) は `statPress` 系列なので直接比較しない |
| B3 | 方式確定の判定 | B1/B2 が 3 解像度 × dv 組合せで通れば確定。codex result レビュー |
| B4 | (後続) SERN 全体への展開 | 上流外部流・機体上面/側面/ベース・プルーム全長。`mesh3d.backend: blocks` で切替 (既定 legacy) |
| B5 | (後続) runner・帳簿 | 新タグの力の帰属、入口の実面積・流束積分、`prepare_info` |
| B6 | (後続) CFD 受理の取り直し | 生産条件で格子収束列 3〜4 点、密度床ノード 0 の確認、起動レシピの再検証 |

## 6. 検証

**実証模型の受入条件** (メッシュのみ。全部通って方式確定):

1. **位相**: 100 % ヘキサ。体積の欠落・重複なし。内部面は 2 要素が共有、外部面は 1 要素 + ちょうど 1 つの境界タグ (§4.10 (d) のタグ表どおり)。
   体積ごとに頂点 Jacobian の符号が揃う。
2. **形状の格子非依存**: `scale` 3 段階で、(i) **CAD 量** = station の点座標・端点位置 ($L_{sw}$, $L_{cowl}$, フィレット終端)・壁曲線の標本点が一致 (相対 1e-9)、
   (ii) **離散量** = 入口面積は解析値 ($-\,(1-\pi/4)(r_u^2 + r_l^2)$ を含む) との差が弦近似の見積もり以内で、細分で単調に減る。(i) と (ii) は別ゲート。
3. **数値的健全性**: float32 化後の**頂点 Jacobian** が全ヘキサで正。node 変換後の双対体積が正、**CV ごとの** $|\Sigma S|/\Sigma|S| \le 1\times10^{-5}$。
   skew ≤ 0.90 (**`SOFT-PASS` は不合格**)。AR ≤ 1000、壁法線の近直交層 (壁面・端面の層) に限り ≤ 5000。
4. **壁層の保証**: §4.10 (e) の節点別検査で、全壁タグ (端面を含む) の $|\hat n\cdot\Delta x|$ が $h_1$ に対し **平面部 ±5 %、稜から $2\delta_r$ 以内と半径が変わる区間 ±10 %、フィレット円弧上 −10 % / +45 %** (§4.11 の 8)。
   層数は指定どおり (黙って減らさない)。壁法線・端面法線・両端細分の各方向で**隣接間隔比 ≤ 1.2** (区間境界を含む)。
5. **評価不能は失敗**: 対向節点が取れない壁節点・タグの無い外部面・判定に必要な量の欠落は、件数が 1 でも不合格。

dv の厳しい組合せ: ランプ初期角 (10° / 22°)、$L_{cowl}$ (短/長) と $L_{sw} = 0.67 L_{cowl}$、$r_f$ (0.02 / 0.10 H)、
$t_{sw}$・$t_c$ (0.005 / 0.02 H)。事前判定で弾かれるべき組合せが**エラーで止まる**ことも試験に含める。

### 6.1 レビュー記録 (codex)

| stage | 日付 | 記録 | 判定 | 採否 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-21 | [2026-09-21-tooling-sern-mesh-blocking-plan.md](../../notes/reviews/2026-09-21-tooling-sern-mesh-blocking-plan.md) | **NO-GO**, C0/M7/m1 (**方式は維持、B1 の前に接続設計を完成させよ**) | **全件採用**。M1 (円弧+直線を 1 辺にはできない・$r_f \to 0$ は固定位相で不成立) → B0、当方のスパイクも同じ結果 (§4.7)。M2 (カウル後流で N/スロット後流/U の面分割が不一致、`z=W/2` を伝播) → B0。M3 (断面リングの第一層は 3D 法線距離と違う、端面直後の x 分布が無い) → B0b。M4 (入口面積不変は誤り: フィレットが面積を減らす。離散面積は弦近似で解像度依存) → §4.1 訂正 + B2。M5 (`--geometry-only` は受入試験にならない) / M6 (体積和では局所反転を見逃す、閉性は CV ごとに) → B0c。M7 (3 解像度が宣言した `scale` 列でない、dv の数値が無い) → B2。m1 (「根治」は本フェーズで検証不能) → 目的を訂正 |
| plan (3 回目, §4.13 限定) | 2026-09-22 | [2026-09-22-tooling-sern-mesh-blocking-plan-2.md](../../notes/reviews/2026-09-22-tooling-sern-mesh-blocking-plan-2.md) | **NO-GO**, C0/M5/m1 (**「真因確定」と「テーパ推奨」を撤回、鈍頭のまま診断を先行**) | **全件採用**。M1 (SLAU の圧力差項と χ=0 の仮説) → §4.13 訂正 + B1d ②。M2 (Dirichlet 0 は熱的固定も外す = 交絡) → 撤回 + B1d ③ (断熱壁 A/B)。M3 (IC の不連続 2 箇所) → `run_junction_model.py` 修正済。M4 (床は保存量に書き戻す・局所 dt・閉性 FAIL 無視) → 記述訂正 + 警告追加。M5 (テーパ A/B の交絡・合否基準) → B1d の設計。m6 (比較対象の固定) → 診断記録の様式 |
| plan (2 回目, §4.10 限定) | 2026-09-22 | [2026-09-22-tooling-sern-mesh-blocking-plan.md](../../notes/reviews/2026-09-22-tooling-sern-mesh-blocking-plan.md) | **GO-with-changes**, C0/M4/m3 (**19 ブロック構成は成立、全ヘキサを維持**) | **全件採用**。M1 (第一層は辺長でなく 3D 法線距離。$r=0$ の隅で −29.3 %、絞り面で 0.770) → §4.10 (e) を $s_1 = h_1/|\hat n\cdot\hat d|$ + 節点別検査に書換。M2 (e5 の節点数は同値類の最長辺から。90 節点では 0.134 m まで) → `Bump` + 111 節点級、当方の実測を §4.10 (e) に。M3 (端面は $t/5$ でなく $h_1$) → §4.10 (e)。M4 (辺名タグとフィレット全部 `sidewall_in` は両立しない) → §4.5 を「45° 点で二分し主壁へ帰属」に。m5 (N 内部 5 共有面と外部タグ表) → §4.10 (d)。m6 (コアの流体内包含) → §4.10 (a)。m7 (§6 と不一致) → §6 を書換 |
| 相談 | 2026-09-21 | [2026-09-21-codex-mesh-approach-2.md](../../notes/reviews/2026-09-21-codex-mesh-approach-2.md) | (a)→**(d) 本方式を推奨** | 採用。本 plan の起点。「最初の実証は側壁終端を含む 3D 接続模型」を B1 に採用 |

## 7. 影響範囲

新規: `design/forge_design/meshing/mesh_sern3d_blocks.py`、`design/tests/run_sern_blockmesh_tests.py`。
既存の `mesh_sern3d.py`・runner・問題 YAML は**本 plan の範囲では変更しない**。依存に gmsh (pip) が加わる。

## 8. 完了条件

実証模型が §6 の 4 群を 3 解像度 × dv 組合せで通り、codex result レビューを経て方式が確定すること。
(SERN 全体の受理は B4–B6 の後続フェーズで、親 plan §8 の許容値に従う。)

## 未確定事項

- **U1 フィレットの終端形状**: 側壁後端の手前で $r_f \to 0$ に絞る案を既定にした (鈍頭終端は後流ブロックに尖点が出てヘキサが切れない)。
  実機形状の意図と違えばユーザ確認のうえ見直す。
- **U2 外側の隅** (側壁外面×機体下面など) にフィレットが要るか。現状は尖ったまま H ブロック。
- ~~**U4 規模と壁解像の範囲**~~ → **決着 (2026-09-22, §4.12)**: (a) 端面を緩める + (c) 第一層 16 µm を採用。外部の壁の層は残した。(元の記述) (2026-09-22): 全壁 $h_1$・成長率 1.2・端面にも層、では SERN 全体が現環境に載らない (§4.11)。候補: (a) 端面の第一層を緩める (`H1_END`。x 面数の 8 割が端面層)、
  (b) 外部の壁 (機体下面・側壁外面・カウル下面) を壁解像の対象から外す、(c) $h_1$ を 5.6–8 µm にする (ランプの実測 $y^+$ は 0.55)、(d) 層の伝播を止める非構造の遷移帯、(e) 大きいホスト。
- **U3 フィレット面の厳密さ**: station 間は曲面パッチ (Coons 補間) で、一定半径の厳密なフィレット面とは O(Δx²) の差がある。
  station を物理的に固定するので解像度には依らない。許容できるか。

## 9. 変更ログ

- `2026-09-23` — **出口 BC の欠陥を特定** (B1f)。ユーザが場を見て「謎の壁」を指摘 → 出口から上流へ育つ圧力域。`outlet_statPress` → `outflow` の A/B (run_0435 対 run_0437) で残差が +0.8 桁 → −1.3 桁。runner 3 本の既定を変更。**壁排出の対策 (χ_n) とは独立の欠陥**で、step を揃えた比較 (800 対 1000) では対策の有無で出口の劣化に差が無い。
- `2026-09-23` — 判別 A/B (run_0433 Roe) で **分岐 A** を確定 (狭い主張: $\chi=0$ が補充を消している)。NaN 順序は測れないと判明し推定で閉じた (run_0434: 床下の節点は restart で `roMin` に戻るため元の破綻を再現しない)。対策は [`convection-slau-wall-normal-chi.md`](convection-slau-wall-normal-chi.md) へ切り出し (文献調査で先行研究が見つかり本命は再検討中)。
- `2026-09-22` — B1d ①②③ を実施 (§4.13.1, run_0430/0431/0432)。**IC 修正は発散時刻を動かすだけで排出は残る**、**断熱壁 A/B で熱的拘束は無関係と分離**、**面別 SLAU 収支で $\chi$=0 の流出面を特定**、**判別 A/B は分岐 B (構造的)**。あわせて **EOS 床が定常 commit で捨てられる**ことを確認し B1e に起票 (`diagnostician` 判断)。
- `2026-09-22` — codex plan レビュー 3 回目 (§4.13) **NO-GO M5/m1 を全件採用**: 「真因確定」「テーパ推奨」を撤回、IC の不連続 2 箇所を修正、鈍頭のまま面別収支と断熱壁 A/B へ (B1d)。
- `2026-09-22` — 接続模型の CFD: 7 run (run_0423–0429) で発散の真因を**鈍頭後端面の壁 CV の排出** (Dirichlet 固定壁 CV は圧力勾配で満たされない) に特定。`nodeWallDirichlet: 0` の A/B で機序確定、段階起動は不成立。**後端のテーパを推奨、ユーザ判断待ち** (§4.13)。
- `2026-09-22` — 接続模型の初回 CFD に着手 (§4.13)。run_0423 は層流暖機 step 454 で発散 (後端面の直後)、run_0424 を後流静止 IC で投入。
- `2026-09-22` — **ユーザ決定: 粗くして 200〜300 万節点** (§4.12)。第一層 16 µm・端面 0.25 mm で 2,446,572 節点、模型の全ゲート PASS、`check_mesh_quality.py` PASS (AR max 676)。同ツールの六面体の体積判定を頂点 Jacobian に変更 (誤判定と見逃しの両方を回帰試験に追加)。
- `2026-09-22` — **B1 完了** (§4.11): 接続模型が scale 0.5 (193 万節点)・0.71 (497 万節点) で受入ゲート PASS。第一層の決め方を η = h1/dr に変更、内側線を壁のオフセットに。
  CV ごとの閉性検査ツールを新設 (新模型 max 3.9e-4 = 変換器の単精度、**旧メッシャの格子に閉性 0.30 の開いた CV 68 個** → 親 plan R5r で修正)。scale 1.0 は 1250 万節点の見積もりで手元に載らない → U4。
- `2026-09-22` — codex plan レビュー 2 回目 (§4.10 限定) **GO-with-changes** C0/M4/m3 を全件採用: 第一層を 3D 法線距離で定義・e5 は `Bump` で節点数を最長辺から決定・端面も $h_1$・フィレットは 45° 点で主壁へ帰属・§6 を実行可能な合否仕様に統一。B0 完了、B1 へ。
- `2026-09-22` — §4.10 に B0 接続設計 (19 ブロック/断面・3 区間・共有面一覧・分布規則・境界タグ規則)。codex plan レビューへ。
- `2026-09-22` — **ユーザ決定 (再): 全ヘキサを本線、Salome は保険** (§4.8.1)。B0 を再開、S2–S5 は保留。
- `2026-09-22` — 全ヘキサの試作形状 (20 ブロック、32 万節点、負 Jacobian 0・skew 0.501) を作りユーザに図で提示。鈍頭終端の尖点の成立条件を整理 (接する 2 壁が一緒に終われば消える)。試作スクリプトを `case/46.sern_design/cad/` に配置。
- `2026-09-22` — 絞り区間の短さの試験 (0.5 $r_f$ でも成立) と非構造の節点数の接線サイズ依存 (1.5 mm 138 万 / 3 mm 39.5 万) を §4.7 に追記。ユーザは「メッシュ数が多くなるなら全ヘキサ」に傾いており、方式の最終確認待ち。
- `2026-09-22` — **ユーザ決定: 非構造 (Salome) を先行、全ヘキサは形状確定後** (§4.8)。Salome スパイクが forge の品質ゲートまで通った (§4.9: SOFT-PASS、AR max 686)。§5.1 に S1–S5 を追加し B0–B3 を保留に。
- `2026-09-22` — §4.7 に追試: 弧長一様の 1 本スプラインで絞り区間の折れが消えた (負 Jacobian 782 → 0)。形状公差の表。フィレット終端の選択肢をユーザに提示。
- `2026-09-21` — codex plan レビュー NO-GO (M7/m1) を全件採用。§5.1 に B0/B0b/B0c を B1 の前提として追加、§4.1 の入口面積の主張を訂正、§4.7 にスパイク結果。
- `2026-09-21` — 起票。親 plan R5q (厚さ 0 交線の欠陥) の移管先。方式は codex 再相談の推奨 (d) を採用。
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
