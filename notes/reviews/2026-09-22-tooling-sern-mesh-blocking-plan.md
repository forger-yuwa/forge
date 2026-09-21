# codex レビュー: tooling-sern-mesh-blocking (plan)

- **plan**: [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md)
- **stage**: `plan`
- **date**: 2026-09-22
- **commit**: `645f25a4` (feature/sern-design)
- **codex**: effort `high`, 7.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m3
- **focus**: §4.10 の B0 接続設計だけを見てほしい (§4.8 の Salome 先行は §4.8.1 で撤回済み、本線は全ヘキサ)。19 ブロック/断面 × 3 区間の胞体複体が全ヘキサで閉じるか: 共有面一覧 (d) に漏れ・不一致がないか、節点数の拘束が矛盾しないか、壁が終わった後の面の扱い、3D 法線補正と許容 (±5%/±10%)、e5 の両端細分、端面直後の x 分布。試作 (case/46.sern_design/cad/hex_junction_demo.py, N の 4 ブロックは負 Jacobian 0・skew 0.501 で成立済み) を拡張する実装方針 (f) の落とし穴。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

19 ブロックの接続構成は成立し、全ヘキサ方式を維持することを推奨します。  
ただし、現行の法線補正では第一層の ±5%／±10% を満たせず、両端細分・端面解像・タグ割当にも実装前の修正が必要です。  
以下は §4.10 に限定した評価です。CFD の収束・定常性は判定していません。

1. **Major — `h₁/cosθr` だけでは、対角辺とフィレット遷移面の法線距離を補正できない。**

   根拠: [plan:315](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:315)、[試作コード:64](/home/sano/work/forge/case/46.sern_design/cad/hex_junction_demo.py:64)、[同:117](/home/sano/work/forge/case/46.sern_design/cad/hex_junction_demo.py:117)。

   試作は**辺に沿った距離**を 4 µm にしています。保存済み `cad/demo_out/hex_junction.msh` を再計測すると、`r=0` の `e2/e6` は辺長 4.000 µm、y・z への投影がそれぞれ **2.828 µm**。水平壁でも **−29.3%** です。`e2` にランプ角補正を掛けても、この断面内の 45° 投影は残ります。

   また、フィレット遷移面の法線には `r′(x)` が入ります。試作の上隅形状から導くと、円弧角を α とした断面法線方向の変位の投影係数は
   
   $$
   \left[1+\{\sin\alpha\,y_r'+(1-\sin\alpha-\cos\alpha)r_u'\}^2\right]^{-1/2}
   $$
   
   です。水平ランプでも、半径を `0.10H → 0`、長さ `0.05H` で線形に絞る場合、円弧中点では **0.770**。ランプ角だけでは補正されません。

   **対案:** 実際に生成する曲面の法線と格子変位について、`|n·Δx|` を目標にする設計へ変更してください。鋭い隅では接する両壁の条件を同時に課します。水平の直角隅なら対角辺の第一間隔は `√2h₁` が必要です。補正は station の辺だけでなく、壁面内部と station 間も検査対象にしてください。

   なお、保存済み 307,200 ヘキサの float32 座標を独立再検査し、非正の頂点 `cornerJ` は **0 個**でした。これは第一層距離の合格を意味しません。

2. **Major — `e5` の両端細分には、全 station を通した節点数の決定手順が必要。**

   根拠: [plan:308](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:308)、[同:320](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:320)。試作は現在、`e5/e7/e11` に一様分布を指定しています（[コード:120](/home/sano/work/forge/case/46.sern_design/cad/hex_junction_demo.py:120)）。

   試作の `ramp` を模型出口 `x/H=2.4` まで延長すると、`e5` の長さは **0.16048 m**。90 節点＝89 間隔、両端 4 µm、隣接間隔の成長率 ≤1.2 で覆える最大長は
   
   $$
   \sum_{i=0}^{88}4\times10^{-6}\,1.2^{\min(i,88-i)}
   =0.13406\ {\rm m}
   $$
   
   なので、90 節点では成立しません。

   **対案:** 両端の目標間隔・辺長・成長率上限から必要節点数を求め、同値類全体の最大値を `N_Y` に採用してください。`e5/e7/e11`、`SW`、`S_de` の分布をまとめて決める必要があります。

   辺を分割せず両端細分する手段はあります。Gmsh の [`Bump`](https://gmsh.info/doc/texinfo/#Structured-grids) をメモリ上で試したところ、この直線・同じ長さでは **101 節点、両端 4 µm、最大成長率 1.1965** が得られました。ただし、この一次元試験は変更後の `N_SIDE` の品質保証にはならないため、非対称半径・短いテーパを含む再検証が必要です。

3. **Major — 端面直後の `Δx≤t/5` は、端面の壁第一層を保証しない。**

   根拠: [plan:322](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:322)、全壁の第一層を要求する [同:364](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:364)。

   `H=0.1 m` の試作値では、`sidewall_end` は `t/5=1 mm`、`cowl_base` は `0.4 mm`。4 µm の **250 倍／100 倍**まで許します。これは後流を板厚に対して分解する条件であり、端面に対する壁法線解像の条件とは別です。

   **対案:** 端面も第一層 4 µm を目標とし、`Δx₁` にその許容範囲と `t/5` 上限を同時に課してください。区間 B は始端・終端の両側に細分が必要です。同じ区間の全ブロックで x 方向節点数を共有し、実際の x 間隔について成長率と区間境界の接続を検査してください。

4. **Major — 「辺名による面タグ」と「フィレット全部を `sidewall_in`」は両立しない。**

   根拠: [plan:110](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:110)、[同:310](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:310)、[試作コード:61](/home/sano/work/forge/case/46.sern_design/cad/hex_junction_demo.py:61)。

   `e1` はランプ直線＋上フィレット半分、`e8` は下フィレット半分＋カウル直線を含む一本の曲線です。これを掃引した一枚の面に、辺名からタグを一つ付ける方式では、フィレット部分だけを `sidewall_in` にできません。「全外部面に一つのタグ」は満たしても、指定したタグの意味が破れます。

   **対案:** 本模型では**複合面を隣接する主壁へ帰属**させ、`e1×x=ramp`、`e5×x=sidewall_in`、`e8×x=cowl_in` と明記してください。フィレットの半分をそれぞれ含む定義に §4.5 を改めれば、19 ブロックを維持できます。

5. **Minor — 共有面一覧には省略があるが、接続そのものの矛盾は見つからない。**

   根拠: [plan:294](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:294)、既存の N の接続を定義する [試作コード:74](/home/sano/work/forge/case/46.sern_design/cad/hex_junction_demo.py:74)。

   断面の接続を独立に数えると、内部共有辺は **A:21、B:23、C:29**。`SW` 追加で 2 本、`CW1/CW2` 追加で 6 本増え、接続は整合します。壁の終了後に `N_SIDE↔SW`、`N_BOT↔CW1` とする扱いも妥当です。

   ただし (d) は N 内部の **`e2/e3/e6/e7/e10` の五面**を明示していません。`ramp`、対称面、S 上端、遠方境界のタグ表も未完成です。さらに A の `e5/e8` はフィレットを含むため、表中の平面座標の説明は厳密には成立しません。

   **対案:** 既存 N の五共有面と全外部面を表へ追加してください。共有面は一つの正準向きで保持し、各体積からの参照符号は逆向きにします。参照数だけでなく、各体積の面の向きが閉じることも検査してください。

6. **Minor — 事前判定は、試作コードの固定コアが流体内に残る条件を欠く。**

   根拠: [plan:267](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:267)、[試作コード:64](/home/sano/work/forge/case/46.sern_design/cad/hex_junction_demo.py:64)。

   固定コア角 `CUR=(y_r−δr,z_w−δr)` と円弧中点 `WU` は、
   `r_u=δr/(1−1/√2)=0.273137H` で一致し、`e2` が潰れます。例えば `r_u=0.30H, r_l=0, 高さ=幅=H` は現在の事前条件を通りますが、このコア配置には使えません。

   **対案:** 各 station と遷移区間について、コアの流体内包含・辺長の正値・輪郭の非交差を事前条件に追加してください。外部帯についても `y_bot<y_cl−δr`、`z_o+δr<Z_far` を明示してください。

7. **Minor — §4.10 と §6 の受入条件を同期する必要がある。**

   根拠: [plan:319](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:319)、[同:347](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:347)、[同:359](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:359)。

   §6 は依然として全壁 ±5%、離散入口面積の解像度間一致、`--geometry-only` による保証を要求しています。§4.10 の隅 ±10%、B0c/B2 の修正と一致していません。

   **対案:** §6 を実行可能な合否仕様に統一してください。壁の許容区分、端面目標、成長率、CAD 誤差と離散化誤差の区別、CV ごとの閉性閾値を明記し、評価不能を失敗扱いにします。

**推奨は、19 ブロック構成を維持して全ヘキサで進めることです。** 実装前の優先順は **1 → 2 → 3 → 4 → 5〜7**。その後、B0c のゲートで既存 N を変更後の分布まで検証し、外部ブロックへ拡張してください。節点数の同値類や壁終了後の接続を理由に、方式を撤回する必要はありません。

ファイルは変更していません。上記の提案は **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 3
