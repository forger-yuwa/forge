#!/usr/bin/env python3
"""接続模型 (全ヘキサ): 側壁終端・カウル終端を含む 3D 接続を gmsh の transfinite ブロックで切る。
設計は plans/active/tooling-sern-mesh-blocking.md §4.10 (19 ブロック/断面 × 3 区間)。

  N  : ダクト内のバタフライ (TOP/SIDE/BOT/CORE)。壁曲線は弧長一様の 1 本スプライン (§4.7)
  SW : 側壁の跡 (x >= L_sw)          CW1/CW2 : カウルの跡 (x >= L_cowl)
  U  : カウル下 (2 列 × 2 帯)         S : 側壁の外 (2 列 × 4 帯)

第一層は「辺長」でなく壁面の 3D 法線に対する距離で決める (§4.10 (e))。端面 (sidewall_end / cowl_base) にも壁層を置く。
境界タグは「参照する体積が 1 つの面」を機械的に拾い、(面の名前, 区間) の表で決める。

B4a (plan §4.14-2〜5・§6.2、§5.1 B4-1) で MOC 輪郭に置き換えた (旧: 放物線ランプ・平坦カウル・断面隅フィレット。後方互換なし):
  - ランプ y_r(x) = 参照輪郭 (`ramp_fillet` 適用後の 2D 輪郭)、カウル内壁 (cowl_in) y_c(x) = cowl_xy。どちらも点列の 3 次 B-spline 補間
    (`Contour`、点列の外は端の接線で直線延長)。入力は `--contours DIR` の ramp_contour.csv / cowl_contour.csv (列 x_m,y_m、単位 m。
    export_contours.py で作る)。
  - カウル外壁 (cowl_out) = 内壁を内壁法線方向 (外 = 下) に板厚 TC だけオフセットした曲線 (テーパなし)。後端は内外とも平面 x = L_cowl で切る
    (オフセット曲線を延長して平面と交わらせる)。x > L_cowl のカウル跡 (CW1/CW2) の上下は後縁から接線方向に直線延長。
  - 側壁: 内面 z = W/2、厚さ TSW を +z 側へ、後端は平面 x = L_sw。断面隅フィレットは全て 0。
  - x 方向は**格子の全 x 位置を station にする** (各 station で断面を物理位置から作り、x 区間は 1 セル幅)。全節点が station 面上にあり、
    壁節点は参照輪郭そのものの上に乗る (transfinite の x 方向補間で輪郭から外れることがない)。
  - 後流区間 [L, L + LWAKE] (L = L_sw, L_cowl) の全体で Δx <= DXW、以降は比 <= 1.2 で HX / HX_FAR へ。成長率は生成 1.2・実測の隣接比 <= 1.2。
  - 検査 (VERDICT): 位相・タグ・頂点 Jacobian (float32)・skew/AR・壁第一層 (first_layer_check) に加えて、
    出力実座標の隣接間隔比 (adjacent_spacing_check、全ヘキサの面隣接・継ぎ目を含む)、後流 Δx (wake_dx_check)、形状ゲート (shape_check)。

B4-2 (plan §5.1、§6.3 の局所 A′/B′ = seam_coupon.py で B′ 合格を確認してから全継ぎ目へ): 断面の分布を継ぎ目で連動させた。
  - 共有辺を挟む 2 ブロックの「共有辺の端から出る辺」の第一間隔を比 1.2 以内にする (match_spacing: 隣の全部と比 g 以内の幾何平均、解けなければ停止)。
  - z ∈ [0, W/2] の同値類は継ぎ目 z = W/2 側 (内側の隅 CUR/CLR) へ片側等比、第一間隔は全辺で F。節点数は NZ = max(P.NZ, prog_n(W/2, h1, 1.15) + 1)。
  - リング横断 (e2/e4/e6/e9) は両端細分 (Bump、NR 区間、増加側 >= NL)。リング外縁が h1 級になりコアの端 (e7/e11/e3/e10) と 1.2 以内でつながる。
  - e5 (と SW・S の de) の端は E5K 倍、側壁外面帯 hz32/hz42 の第一間隔は KV 倍 (P0 の注記)。e5 の x 方向ブレンドと e7/e11 の max クランプは撤去。
  - cowl_side は一般壁 (h1・NL) として検査 (END_TAGS から外した)。
  - 検査は msh を書き出してから gmsh を解放して numpy だけで行う (ピーク RSS 約 1.4 kB/節点)。
  - --set SEC2D=x1,x2,... は指定 station の断面だけを 2D で切って継ぎ目の比と skew を出す (分布の調整用の診断。判定はしない)。

B4-4 (plan §5.1、§6.4 の局所 A/B): リング対角辺 e2/e6 の拘束解除。
  - RING_CORE_RATIO (= コア側端間隔 / 壁側端間隔) が 1 なら従来の対称 Bump (A)。1 以外なら壁側端間隔を保持したまま、コア側端だけを
    その比にした非対称分布 (B、同じ区間数・全長。ring_diag_spacing)。gmsh の transfinite 曲線は 2 端独立の端間隔を持てないので、
    generate(1) の後に直線辺の内部節点を明示位置へ移してから 2D/3D を張る (apply_node_laws)。
  - 壁側の F (e1/e8/hz?0 等) は保持し、コア側 e3/e7/e10 の端 Fc は新しいコア側端間隔から連動 (match_spacing)。比 1 では Fc = F (A と同一)。
  - --set COUPON_X=x COUPON_N=k: 全模型と同じ分布・同じ station のうち x/H に最も近い station を含む k 個だけを押し出した 3D 接続試験片
    (本体と同じ検査関数。区間の端面が無いので後流 Δx・形状ゲートは対象外、壁層は試験片にある壁タグだけ)。

B4-4 (plan §6.5): SW の z 分布 hz31/hz41 の端間隔 a3/a4 の選定則 SW_END_RULE (0 = match_spacing = A・既定、1 = s6/2^(1/4)・s2/2^(1/4) = B。sw_end_spacing)。
  --set SECX=x1,x2,... (/H) は全長診断: 模型の終端を参照ランプ輪郭の終端 (または SECX_END) まで現行の延長則で延ばし、station と分割数を
    その全域から現行規則で決め直したうえで、**指定 x そのもの** (station への丸めなし) の断面を 2D で切って全継ぎ目・ブロック内の比と skew を出す。
    要求 x と評価 x を並べて出す。XEND (1.5 H) を超える x の延長則は暫定 (判定はしない)。

B4b-4 (plan §5.1、§4.15 Minor・§4.16 (6)): 1D 分布の桁あふれと最大間隔の判定。
  - prog_r: r^n があふれる所 (n ln r > 700) だけ等比和を対数で比べ、二分区間は旧と同じ区間が根を挟まないときだけ広げる (旧は prog_r(1.224, 220, 0.001) で
    OverflowError、区間の外の根には端 50 / 1e-3 を黙って返していた)。解の無い入力は ValueError。march の等比和も geo_sum で桁あふれしない。
    全域を対数にしないのは、既定の生成で旧と最大 4 ulp ずれて生成結果がビット同一でなくなるため (旧が解けていた入力では旧と同じビット)。
  - 遠方帯 (vy0k・hz?3) の区間数 far_n: 実際に張る全組 (帯長, 第一間隔) の実間隔列 (prog_spacing) の最大が HFAR 以下。旧の近似式 L (1 - 1/r) は
    既定の模型で 0.0991 H と判定していたが実際の最大は 0.109 H だった → 既定の NFY・NFZ は 20 -> 22 に変わる (それ以外の分布はビット同一)。
B4b-5 (plan §5.1、§4.16): 壁タグごとの第一内部点距離を入力にする。
  - --set H1_<tag>=v (/H) で壁ごとに上書き (既定は一般壁 H1・端面 H1_END)。cowl_side は sidewall_out と同一平面で z 分布を共有するので同じ値
    (指定しなければ従い、違えば停止)。--preset design|final は plan §4.16 の暫定の試験値の組 (PRESETS)。期待層数も壁別 (prog_n(DR, h1_w, 1.2) + 4)。
  - 隅の対角 e2/e6 の壁側端は接する 2 壁の h1 から (diag_first)。リングの補間は隅の比を壁の全幅へ運ぶので、隅の比が ±5 % の外の組合せは止める。
    同一平面の 2 タグの継ぎ目 (x >= L_sw の ramp | vehicle) は節点が両タグの検査を受けるので、1 つの間隔で両目標の ±10 % に入らない組合せは止める。
  - 両端の第一間隔が違う同値類: SW 列 hz?1 (sidewall_in 側 -> sidewall_out 側) と CW 列 vy2k (cowl_out 側 -> cowl_in 側) は両端分布 two_end_spacing、
    側壁帯の y (e5・vy32・vy33・vy34) は同じ Bump を滑らかに曲げた分布 (e5 = vy32 の逆順。薄い SW ブロックの両側で節点の y を揃える)。全壁が等しい既定は旧のまま。
  - --set ESTIMATE=1: gmsh の格子を張らずに節点数・ヘキサ数 (count_mesh、全節点が station 面上にあるので断面の和) と解けない組合せの一覧を返す。

usage (mesh venv):  .venv-mesh/bin/python case/46.sern_design/cad/hex_junction_model.py OUT_PREFIX --contours DIR [--preset design|final] [--scale S] [--set k=v ...]
出力: OUT.msh / OUT_report.json / OUT_sections.npz (描画用)。ESTIMATE=1 は標準出力に JSON だけ
"""
import sys, math, json, argparse, hashlib
from pathlib import Path
import numpy as np
import gmsh

# ---------------------------------------------------------------- 物理入力 (すべて /H)
# 数値仕様は plan §6.2 の表: h1 1.6e-4 H、端面 h1e 1.0e-3 H、後流 [L, L + 0.02 H] で Δx <= 1.0e-3 H、成長率 1.2。
# HX / HX_FAR / NZ は §4.12 の採用値 (0.05 / 0.08 / 31)。TC (カウル板厚 = 後縁厚) は §6.2 の A/B で 0.02 / 0.005 を与える (既定は生産値 0.005)
P0 = dict(H=0.1, ZW=1.0, TSW=0.005, TC=0.005, LSW=0.8, LCOWL=1.2, XEND=1.5, DR=0.08, YBOT=-1.0, ZFAR=2.0,
          ZONE=0.3,                             # 端面の上下流に置く物理 station までの距離
          LWAKE=0.02, DXW=1.0e-3,               # 後流区間の長さと、その区間全体の x 間隔の上限
          H1=1.6e-4, H1_END=1.0e-3, G=1.2, NZ=31, HX=0.05, HX_FAR=0.08, M_SPL=129,
          E5K=1.03,                             # e5 (と SW・S の de の y) の両端 = h1 (c + cc)/2 × E5K。N_SIDE で e7 の端 F (>= 対角の s2/1.2) との差を縮め、
                                                # リングの格子線の傾きを抑える (N_SIDE|SW の比)。B・C の露出壁の第一層はランプ +1.4 %・カウル上面 +4.6 % (稜の近傍 ±10 %)
          KV=1.02,                              # 側壁外面帯 (hz32・hz42) の第一間隔 / h1。SW の z (hz41) の両端を e2 の第一間隔と hz42 の両方に 1.2 以内で
                                                # つなぐための余裕 (KV = 1 だと窓 [s2/1.2, 1.2 h1] の幅が 0.06 %)。sidewall_out の第一層は +2 % (±5 % の内側)
          RING_CORE_RATIO=0.7071067811865476,                   # リング対角辺 e2/e6 のコア側端間隔 / 壁側端間隔 (plan §6.4。1.0 = 旧の対称 Bump = A、B は 1/√2)
          SW_END_RULE=1,                         # SW の z 分布 hz31/hz41 の端間隔 a3/a4 の選定則 (plan §6.5。0 = match_spacing = A、1 = s6/2^(1/4)・s2/2^(1/4) = B)
          HFAR=0.10, AR_TAN=900.0)                            # 外部遠方の最大格子幅 (端面層が断面全体に伝播するので、遠方の幅が AR を決める)
# B4-2 (plan §5.1、§6.3): NZ は下限。z 方向 (e1/e3/e8/e10/hz?0) の節点数は「継ぎ目側の第一間隔 (>= h1) から成長率 1.15 で 1 H を覆う」数との大きい方
# (x/H, 上隅 r/H, 下隅 r/H)。区間内は smoothstep (両端で r' = 0)
# 半径を変える区間の長さは 2 r 以上にする: 短いと r'(x) で壁法線が x に傾き、円弧上の第一層が薄くなる (長さ r で -19 %、2 r で -5 %)
# B4a (plan §4.14-5): 断面隅フィレットは初版ゼロ。半径を明示的に 0 にする (旧 PROFILE の非ゼロ半径は流用しない)
PROFILE = [(0.00, 0.00, 0.00)]


def smooth(t):
    t = min(max(t, 0.0), 1.0); return t * t * (3.0 - 2.0 * t)


# ---------------------------------------------------------------- 1D 分布
# B4b-4 (plan §5.1、§4.15 Minor・§4.16 (6)): 等比和 Σ_{i<n} r^i の評価で r^n が桁あふれしない (旧 prog_r(1.224, 220, 0.001) は OverflowError)。
# n ln r が LN_BIG を超える所だけ対数 ln Σ = n ln r + ln(1 - r^-n) - ln(r - 1) で比べ、それ以外は旧と同じ算術のまま
# (全域を対数にすると既定の生成で 4218 回中 130 回が最大 4 ulp ずれ、生成結果がビット同一でなくなるため)
LN_BIG = 700.0


def geo_sum(r, n):
    """等比和 Σ_{i<n} r^i (r != 1)。桁あふれする r (> 1) では inf を返す (旧式 (r^n - 1)/(r - 1) と同じ算術)"""
    if r > 1.0 and n * math.log(r) > LN_BIG: return math.inf
    return (r ** n - 1) / (r - 1)


def _geo_lt(h, r, n, L):
    """h Σ_{i<n} r^i < L の判定 (桁あふれなし)。r^n があふれる所だけ対数で比べる"""
    if r > 1.0 and n * math.log(r) > LN_BIG:
        t = n * math.log(r)
        return math.log(h) + t + math.log(-math.expm1(-t)) - math.log(r - 1.0) < math.log(L)
    return h * (r ** n - 1) / (r - 1) < L


def prog_r(L, n, h):
    """n 区間・第一間隔 h・全長 L の公比 (r < 1 も許す)。解が無い (n < 2 で n h != L、または h >= L > n h でない r < 1 側) なら ValueError。
    二分区間は旧と同じ [1 + 1e-12, 50] (r > 1) / [1e-3, 1 - 1e-12] (r < 1) を使い、根を挟まないときだけ広げる (r > 1 は上端を
    2 (L/h)^(1/(n-1)) へ: h r^(n-1) <= L なので根はその下。r < 1 は下端を 0 へ)。旧は区間の外の根に区間の端 (50 / 1e-3) を黙って返していた"""
    if abs(n * h - L) < 1e-12 * L: return 1.0
    if n < 2 or not (h > 0.0 and L > 0.0): raise ValueError("等比分布の公比が決まらない (L=%.6g, n=%d, h=%.6g)" % (L, n, h))
    if n * h < L:
        lo, hi = 1.0 + 1e-12, 50.0
        if _geo_lt(h, hi, n, L): hi = 2.0 * (L / h) ** (1.0 / (n - 1))
    else:
        if h >= L: raise ValueError("第一間隔が全長以上で等比分布にならない (L=%.6g, n=%d, h=%.6g)" % (L, n, h))
        lo, hi = 1e-3, 1.0 - 1e-12
        if not _geo_lt(h, lo, n, L): lo = 0.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if _geo_lt(h, r, n, L) else (lo, r)
    return 0.5 * (lo + hi)


def prog_n(L, h, g):
    """公比 <= g で覆える最小の区間数"""
    if h >= L: return 1
    return max(int(math.ceil(math.log(1.0 + L * (g - 1.0) / h) / math.log(g) - 1e-9)), 1)


def prog_spacing(L, n, h):
    """n 区間・全長 L・第一間隔 h の等比分布の実間隔列 (始点側から。gmsh の Progression と同じ則)。h r^i <= L なので r^i はあふれない"""
    r = prog_r(L, n, h); d = h * r ** np.arange(n); d *= L / d.sum()
    return d


def prog_nodes(L, n, h, at_end=False):
    """n 区間・全長 L の等比分布の節点位置 (0..L)。第一間隔 h を始点側 (at_end=True なら終点側) に置く。gmsh の Progression と同じ則"""
    d = prog_spacing(L, n, h)
    if at_end: d = d[::-1]
    return np.r_[0.0, np.cumsum(d)]


def far_n(pairs, hmax, g_t, nmax=400):
    """遠方帯 (壁帯の外の vy0k・hz?3) の区間数 (B4b-4、plan §4.15 Minor)。pairs は実際に張る (帯長, 第一間隔) の組の全部。
    公比 <= g_t で覆える数から始め、**実間隔列の最大** (prog_spacing の最大) が全組で hmax 以下になるまで増やす。
    旧は 1 組の概算 (第一間隔 = 壁帯の最外間隔の概算 × g_t) と末尾間隔の近似式 L (1 - 1/r) で判定しており、既定 (h1 16 µm) の
    接続模型で実際の最大は 0.109 H (> HFAR 0.10 H) だった。nmax で止まったら最大間隔は hmax を超えたまま (戻り値の 2 番目で判る)"""
    n = max(prog_n(L, h0, g_t) + 1 for L, h0 in pairs)
    mx = lambda n_: max(float(prog_spacing(L, n_, h0).max()) for L, h0 in pairs)
    while mx(n) > hmax and n < nmax: n += 1
    return n, mx(n)


def match_spacing(nbr, g):
    """継ぎ目の間隔の連動 (plan §5.1 B4-2): 共有辺の端で隣り合う間隔 nbr (複数) の**全部と比 g 以内**になる間隔を返す。
    可能な範囲 [max(nbr)/g, g min(nbr)] の幾何平均 (両側の比が等しく最小)。範囲が空なら ValueError (分布だけでは解けない継ぎ目)"""
    nbr = [float(v) for v in nbr]; lo, hi = max(nbr) / g, g * min(nbr)
    if lo > hi * (1 + 1e-12): raise ValueError("継ぎ目の間隔を比 %.3g 以内で連動できない: 隣の間隔 %s (比 %.4f > %.4f)" % (g, nbr, max(nbr) / min(nbr), g * g))
    return math.sqrt(lo * hi)


def ring_diag_spacing(base, ratio, g=1.2):
    """リング対角辺の非対称分布 (plan §6.4 の B、§5.1 B4-4)。base は対称分布の間隔列 (先頭 = 壁側、末尾 = コア側、合計 = 辺長)。
    壁側の端間隔 base[0] を保持し、コア側の端間隔を ratio × base[0] にする。区間数と全長も保持する。
    d_i = base_i exp(p(τ_i))、τ_i = i/(n-1)、p(τ) = ln(ratio base[0]/base[-1]) τ + c τ(1-τ) で、c は全長が base と等しくなるよう二分法で決める
    (ratio = base[-1]/base[0] なら c = 0 で base そのもの)。内部の隣接比が g を超えたら ValueError (分布だけでは満たせない)"""
    b = np.asarray(base, float); n = len(b)
    if n < 3: raise ValueError("区間数が少なすぎる (n=%d)" % n)
    t = np.arange(n) / (n - 1.0); a = math.log(ratio * b[0] / b[-1]); L = b.sum()
    f = lambda c: float((b * np.exp(a * t + c * t * (1 - t))).sum()) - L
    lo, hi = -60.0, 60.0
    if f(lo) > 0 or f(hi) < 0: raise ValueError("非対称分布の全長を合わせられない")
    for _ in range(200):
        c = 0.5 * (lo + hi); lo, hi = (c, hi) if f(c) < 0 else (lo, c)
    d = b * np.exp(a * t + 0.5 * (lo + hi) * t * (1 - t)); d *= L / d.sum()      # 二分法の残差 (1e-16 相対) を全長に吸収
    q = np.maximum(d[1:] / d[:-1], d[:-1] / d[1:]).max()
    if q > g * (1 + 1e-12): raise ValueError("非対称分布の内部隣接比 %.4f > %.3g (ratio %.4g, n %d)" % (q, g, ratio, n))
    return d


SW_R0 = 2.0 ** 0.25


def sw_end_spacing(s2, s6, F, kv, g, rule=0):
    """SW の z 分布の端間隔 (a3 = hz31 の WL 側、a4 = hz41 の WU 側) の選定則 (plan §6.5、§5.1 B4-4)。
    rule 0 (A、既定): a3 = match_spacing([s6, F, kv])、a4 = match_spacing([s2, kv]) (隅の対角端・e8 の F・側壁外面帯 kv と比 g 以内の幾何平均)。
    rule 1 (B): a3 = s6/2^(1/4)、a4 = s2/2^(1/4)。両隅の s/a と、e5 中央の N_SIDE 第一層 t_N = (s2 + s6)/(2√2) に対する
    平均間隔 (a3 + a4)/2 の比がどれも 2^(1/4) になる (codex diagnose 2026-10-06 の端点・中央の拘束)。B は F・kv との比を拘束しない (3D の検査で見る)"""
    if rule == 0: return match_spacing([s6, F, kv], g), match_spacing([s2, kv], g)
    if rule == 1: return s6 / SW_R0, s2 / SW_R0
    raise ValueError("SW_END_RULE は 0 か 1 (%r)" % rule)


# ---------------------------------------------------------------- 壁別の第一内部点距離 (B4b-5、plan §4.16・§5.1 B4b-5)
# 壁タグ -> 第一内部点距離 h1 (/H)。既定は一般壁 H1・端面 H1_END (旧の共通値)。P["H1_<tag>"] で壁ごとに上書きする。
# 全壁が等しいときの生成は旧と同じ演算 (ビット同一) になるよう、各式は「旧の h1 を壁別の値に置き換えた形」にしてある。
WALL_GENERAL = ("ramp", "vehicle", "sidewall_in", "sidewall_out", "cowl_in", "cowl_out", "cowl_side")
WALL_ENDS = ("sidewall_end", "cowl_base")
# 設定のプリセット (plan §4.16 の「暫定の試験値」の表。受理値ではない。単位 H = 0.1 m なので 16 µm = 1.6e-4 H)。
#   design: 設計ループ用の粗い設定。細かくする壁 16 µm、緩める壁 64 µm、端面 h1e 0.1 mm、後流 Δx <= 0.1 mm を 0.02 H、HX 0.10・HX_FAR 0.16 H
#   final : 最終 (暫定)。細かくする壁は y1+ <= 1 の比例の目安 (cowl_in 3・ramp 7・sidewall_in 2 µm)、cowl_out は呼び出し側の暫定 16 µm
#           (plan の表は「2〜7 µm、壁ごと」で cowl_out の個別値なし。新しい形状の壁応力で決め直す)。緩める壁・端面・後流は design と同じ、HX/HX_FAR は現行 scale 1
_RELAXED = dict(H1_sidewall_out=6.4e-4, H1_cowl_side=6.4e-4)   # 側壁跡上の帯 vehicle は ramp と同一平面なので ramp に従う (wall_h1)
_ENDS_WAKE = dict(H1_sidewall_end=1.0e-3, H1_cowl_base=1.0e-3, LWAKE=0.02, DXW=1.0e-3, HFAR=0.10, TC=0.005, TSW=0.005)
PRESETS = {
    "design": dict(H1_ramp=1.6e-4, H1_cowl_in=1.6e-4, H1_cowl_out=1.6e-4, H1_sidewall_in=1.6e-4, **_RELAXED, **_ENDS_WAKE, HX=0.10, HX_FAR=0.16),
    "final": dict(H1_ramp=7.0e-5, H1_cowl_in=3.0e-5, H1_cowl_out=1.6e-4, H1_sidewall_in=2.0e-5, **_RELAXED, **_ENDS_WAKE, HX=0.05, HX_FAR=0.08),
}


def wall_h1(P):
    """壁タグ -> h1 (/H)。cowl_side は sidewall_out と同一平面 (z = z_o) で z 分布 (hz?2) を共有するので同じ値にする:
    H1_cowl_side を与えなければ sidewall_out に従い、違う値を与えたら ValueError (plan §6.2 の仕様訂正・B4b-5)。未知の H1_* キーも ValueError"""
    tags = WALL_GENERAL + WALL_ENDS
    bad = [k for k in P if k.startswith("H1_") and k != "H1_END" and k[3:] not in tags]
    if bad: raise ValueError("未知の壁タグの第一層指定 %s (壁タグ: %s)" % (bad, ", ".join(tags)))
    out = {t: float(P.get("H1_" + t, P.get("H1_END", P["H1"]) if t in WALL_ENDS else P["H1"])) for t in tags}
    if "H1_cowl_side" not in P: out["cowl_side"] = out["sidewall_out"]
    elif not math.isclose(out["cowl_side"], out["sidewall_out"], rel_tol=1e-12, abs_tol=0.0):
        raise ValueError("cowl_side (%.6g H) と sidewall_out (%.6g H) の第一層が違う: 同一平面 z = z_o で z 分布 hz?2 を共有するので同じ値にする"
                         % (out["cowl_side"], out["sidewall_out"]))
    # 側壁跡上の帯 vehicle (y = y_r、z ∈ [z_w, z_o]) は ramp と同一平面で第一内部点の列 (vy32 の第一節点) を共有する。
    # 違う第一層は両立しない (16/64 µm で比 2.03/0.51、稜の許容 ±10 % の外) ので ramp に従う (plan §4.16 の仕様訂正、2026-10-07)
    if "H1_vehicle" not in P: out["vehicle"] = out["ramp"]
    elif not math.isclose(out["vehicle"], out["ramp"], rel_tol=1e-12, abs_tol=0.0):
        raise ValueError("vehicle (%.6g H) と ramp (%.6g H) の第一層が違う: 同一平面 y = y_r で第一内部点の列を共有するので同じ値にする"
                         % (out["vehicle"], out["ramp"]))
    return out


def diag_first(h_side, h_w, cw, k, elen, dr):
    """隅の対角辺 e2 (ランプ × sidewall_in) / e6 (cowl_in × sidewall_in) の壁側端の間隔。接する 2 壁の目標 (側壁は z 方向の h_side、
    ランプ/カウル内壁は鉛直 cw h_w) の算術平均を、対角の辺長に対する比 (elen/dr 倍) で置く。h_w == h_side なら旧式 h1 0.5 (1 + cw) k elen/dr と同じ演算"""
    return h_side * 0.5 * (1 + cw * (h_w / h_side)) * k * elen / dr


def corner_ratios(s, elen, dr, h_side, h_w, cw):
    """隅の第一内部点 (対角の第一節点) の、接する 2 壁の目標に対する比 (sidewall_in の比, ランプ/カウル内壁の比)。
    対角の第一節点の壁法線成分は s dr/elen (r = 0 の 45° 対角で y・z とも s/√2)。リングの transfinite 補間は層 1 の位置を
    両側辺 (e4/e2、e2/e6、e6/e9) の値から線形に補うので、隅の比はそのまま壁の全幅に運ばれる (平面部の許容 ±5 % に収める必要がある)"""
    d = s * dr / elen
    return d / h_side, d / (cw * h_w)


CORNER_TOL = 0.05       # 隅の比の許容 (平面部の ±5 %。上の理由で稜の区分 ±10 % ではない)
RIDGE_TOL = 0.10        # 同一平面の継ぎ目 (2 つの壁タグに属する節点 = 稜の区分 ±10 %)


def station_conflicts(x_h, c, cc, s, elen, dr, hw, a_top, a_v, sw_present):
    """壁別 h1 で解けない組合せの検出 (B4b-5。黙って h1 を変えない)。1 station 分の問題点のリスト (空なら解ける)。
      - 隅 WU (ramp × sidewall_in)・WL (cowl_in × sidewall_in): 対角の第一節点 1 つで 2 壁の目標を同時に満たす。比が ±CORNER_TOL の外なら不可
      - 同一平面の ramp | vehicle (x >= L_sw の y = y_r, z = z_o): 継ぎ目の節点は両タグに属し、第一内部点 (vy32 の第一節点) が 1 つなので、
        1 つの鉛直間隔 d で両方の目標の ±RIDGE_TOL に入る必要がある (d は 2 つの目標の幾何平均)"""
    out = []
    for nm, e, hw_, cw in (("WU (ramp × sidewall_in)", "e2", hw["ramp"], c), ("WL (cowl_in × sidewall_in)", "e6", hw["cowl_in"], cc)):
        rs, rw = corner_ratios(s[e], elen[e], dr, hw["sidewall_in"], hw_, cw)
        if abs(rs - 1) > CORNER_TOL or abs(rw - 1) > CORNER_TOL:
            out.append(dict(kind="corner", where=nm, x=x_h, ratio_sidewall_in=rs, ratio_wall=rw, tol=CORNER_TOL,
                            h1_um=[hw["sidewall_in"] * 1e6, hw_ * 1e6],
                            reason="隅の対角の第一節点 1 つで 2 壁の目標 (第一内部点距離) を同時に満たせない (リングの補間で隅の比が壁全体へ運ばれる)"))
    if sw_present and a_top != a_v:
        d = math.sqrt(a_top * a_v); rr, rv = d / (c * hw["ramp"]), d / (c * hw["vehicle"])
        if abs(rr - 1) > RIDGE_TOL or abs(rv - 1) > RIDGE_TOL:
            out.append(dict(kind="coplanar", where="ramp | vehicle (y = y_r, z = z_o, x >= L_sw)", x=x_h, ratio_ramp=rr, ratio_vehicle=rv, tol=RIDGE_TOL,
                            h1_um=[hw["ramp"] * 1e6, hw["vehicle"] * 1e6],
                            reason="同一平面の 2 タグの継ぎ目の節点は両タグの検査を受け、第一内部点が 1 つなので 1 つの間隔で両目標の ±10 % に入れない"))
    return out


def two_end_spacing(L, n, a, b, g):
    """n 区間・全長 L、始端の間隔 a・終端の間隔 b の分布 (B4b-5: 両端の壁・継ぎ目の第一間隔が違う辺)。
    d_i = min(a q^i, b q^(n-1-i), M): 両端から公比 q (<= g) で増やし、中央は M で頭打ち。端の間隔は a・b、隣接比は <= g。
    M = max(a, b)・q = g で長さが足りなければ q = g で M を二分法、余れば M = max(a, b) で q を [(max/min)^(1/(n-1)), g] の二分法。
    区間数が少なすぎ・多すぎで合わなければ ValueError"""
    if not (L > 0 and a > 0 and b > 0 and n >= 2): raise ValueError("両端分布の入力が不正 (L=%.6g, n=%d, a=%.6g, b=%.6g)" % (L, n, a, b))
    i = np.arange(n, dtype=float); lo_, hi_ = min(a, b), max(a, b); la, lb = math.log(a), math.log(b)
    qmin = (hi_ / lo_) ** (1.0 / (n - 1))
    if qmin > g * (1 + 1e-12): raise ValueError("両端分布: 区間が少なすぎる (端の比 %.4g を公比 %.3g 以下で %d 区間に収められない)" % (hi_ / lo_, g, n))
    def dist(q, M):
        lq = math.log(q); return np.exp(np.minimum(np.minimum(la + i * lq, lb + (n - 1 - i) * lq), math.log(M)))
    S = lambda q, M: float(dist(q, M).sum())
    if S(qmin, hi_) > L * (1 + 1e-12): raise ValueError("両端分布: 区間が多すぎる (n=%d、端 %.4g・%.4g の等比だけで全長 %.4g を超える)" % (n, a, b, L))
    if S(g, math.inf) < L * (1 - 1e-12): raise ValueError("両端分布: 区間が少なすぎる (n=%d、公比 %.3g でも全長 %.4g に届かない)" % (n, g, L))
    if S(g, hi_) >= L:
        lo, hi = qmin, g
        for _ in range(200):
            q = 0.5 * (lo + hi); lo, hi = (q, hi) if S(q, hi_) < L else (lo, q)
        d = dist(0.5 * (lo + hi), hi_)
    else:
        lo, hi = hi_, float(dist(g, math.inf).max())
        for _ in range(200):
            M = 0.5 * (lo + hi); lo, hi = (M, hi) if S(g, M) < L else (lo, M)
        d = dist(g, 0.5 * (lo + hi))
    d = d * (L / d.sum())
    q_ = float(np.maximum(d[1:] / d[:-1], d[:-1] / d[1:]).max())
    if q_ > g * (1 + 1e-9): raise ValueError("両端分布の隣接比 %.6f > %.3g" % (q_, g))
    return d


def two_end_n_range(L, a, b, g):
    """two_end_spacing が解ける区間数の範囲 (n_lo, n_hi)。無ければ None。下限 = 公比 g で両端から増やして全長に届く最小、
    上限 = 端から端への等比 (公比 (max/min)^(1/(n-1))) だけで全長を超えない最大"""
    lo_, hi_ = min(a, b), max(a, b); rho = hi_ / lo_
    smin = lambda n: n * lo_ if rho == 1.0 else lo_ * (rho * rho ** (1.0 / (n - 1)) - 1) / (rho ** (1.0 / (n - 1)) - 1)
    def smax(n):
        i = np.arange(n, dtype=float); lg = math.log(g)
        return float(np.exp(np.minimum(math.log(a) + i * lg, math.log(b) + (n - 1 - i) * lg)).sum())
    n = max(2, int(math.ceil(math.log(rho) / math.log(g) - 1e-12)) + 1)
    while smax(n) < L * (1 - 1e-12):
        n += 1
        if n > 10 * L / lo_ + 10: return None
    n_lo = n
    if smin(n_lo) > L * (1 + 1e-12): return None
    lo, hi = n_lo, n_lo + 1
    while smin(hi) <= L * (1 + 1e-12): lo, hi = hi, 2 * hi
    while hi - lo > 1:
        m_ = (lo + hi) // 2; lo, hi = (m_, hi) if smin(m_) <= L * (1 + 1e-12) else (lo, m_)
    return n_lo, lo


def two_end_law(L, n, a, b, g):
    """辺の分布則 (節点数, "Nodes", 相対位置の列) を two_end_spacing から作る (set_curve_law / apply_node_laws で明示配置)"""
    d = two_end_spacing(L, n - 1, a, b, g); return (n, "Nodes", tuple(np.r_[0.0, np.cumsum(d)] / d.sum()))


class Bump:
    """gmsh の Bump 分布を単位長の直線で実測して較正する (間隔は弧長の放物線、端/中央 = coef)"""
    def __init__(self): self.c = {}

    def spacing(self, n, coef):
        k = (n, float("%.10e" % coef))
        if k not in self.c:
            gmsh.model.add("cal"); g = gmsh.model.geo
            l = g.addLine(g.addPoint(0, 0, 0), g.addPoint(1, 0, 0)); g.mesh.setTransfiniteCurve(l, n, "Bump", k[1]); g.synchronize()
            gmsh.model.mesh.generate(1); _, x, _ = gmsh.model.mesh.getNodes(1, l, includeBoundary=True)
            self.c[k] = np.diff(np.sort(x.reshape(-1, 3)[:, 0])); gmsh.model.remove()
        return self.c[k]

    def coef(self, n, h_rel):
        """端の間隔/全長 = h_rel になる係数。一様以上なら None (一様)"""
        if h_rel >= 0.98 / (n - 1): return None
        mk = ("c", n, float("%.6e" % h_rel))
        if mk in self.c: return self.c[mk]
        lo, hi = 1e-8, 0.999
        for _ in range(60):
            c = math.sqrt(lo * hi); lo, hi = (c, hi) if self.spacing(n, c)[0] < h_rel else (lo, c)
        self.c[mk] = math.sqrt(lo * hi); return self.c[mk]

    def fit_n(self, h_rel, g, n0=11, hmax_rel=1.0):
        n = n0
        while True:
            c = self.coef(n, h_rel)
            if c is None: return n
            h = self.spacing(n, c)
            if np.maximum(h[1:] / h[:-1], h[:-1] / h[1:]).max() <= g and h.max() <= hmax_rel: return n
            n += 2 if n < 60 else 4


def set_curve_law(tg, law, moves):
    """辺の分布則を gmsh に渡す。"Nodes" (節点の相対位置の列、B4-4 の非対称対角) は節点数だけ一様で張り、moves に積んで apply_node_laws で移す"""
    n, typ, cf = law
    if typ == "Nodes": gmsh.model.geo.mesh.setTransfiniteCurve(tg, n, "Progression", 1.0); moves.append((tg, np.asarray(cf, float)))
    else: gmsh.model.geo.mesh.setTransfiniteCurve(tg, n, typ, cf)


def apply_node_laws(moves):
    """1D を張り、明示分布の直線辺 (geo の Line、媒介変数は始点 -> 終点で線形) の内部節点を所定の相対位置へ移す (座標と媒介変数の両方)。
    この後の generate(2/3) は 1D を作り直さない (gmsh は既存の次元より上だけ張る) ので、面・体の transfinite はこの節点を使う。
    移した後の節点を読み直して、位置が指定どおりか (相対 1e-12) を確かめる"""
    if not moves: return
    M = gmsh.model.mesh; M.generate(1)
    for tg, fr in moves:
        tags, _, par = M.getNodes(1, tg, includeBoundary=False)
        if len(tags) != len(fr) - 2: raise RuntimeError("明示分布の節点数が合わない (curve %d: %d != %d)" % (tg, len(tags), len(fr) - 2))
        lo, hi = gmsh.model.getParametrizationBounds(1, tg); lo, hi = float(lo[0]), float(hi[0])
        A = np.array(gmsh.model.getValue(1, tg, [lo])); B = np.array(gmsh.model.getValue(1, tg, [hi]))
        for k, i in enumerate(np.argsort(par)):
            f = float(fr[k + 1]); M.setNode(int(tags[i]), list(A + f * (B - A)), [lo + f * (hi - lo)])
        _, xyz, _ = M.getNodes(1, tg, includeBoundary=False); xyz = xyz.reshape(-1, 3)
        f_ = np.sort((xyz - A) @ (B - A) / ((B - A) @ (B - A)))
        if np.abs(f_ - fr[1:-1]).max() > 1e-12: raise RuntimeError("明示分布の節点が移っていない (curve %d)" % tg)


def march(lengths, h0, g_t, g, hmax, caps=None):
    """端面から離れる向きに区間を順にたどり、区間ごとの (区間数, 公比) を決める。最初の区間は第一間隔 = h0。
    以降は「第一間隔 = 前の最終間隔 × 公比」で継ぎ目の間隔比も公比に揃える。なるべく粗く (公比 <= g_t、最終間隔 <= 1.15 hmax)。
    caps[i] (None 以外) は区間 i の**全体**の間隔上限 (後流区間の Δx <= DXW。plan §4.14-4): 区間内の全間隔 (= 単調なので始端と終端) が
    caps[i] 以下の分割だけを採る (1.15 の緩和なし。丸めの 1e-9 相対だけ許す)"""
    out, last = [], None
    caps = caps or [None] * len(lengths)
    for L, cap in zip(lengths, caps):
        if cap is not None:
            n_cap = int(math.ceil(L / cap - 1e-9))
            for n in range(n_cap, n_cap + 400):
                if last is None:
                    try: r = prog_r(L, n, h0)
                    except ValueError: continue         # n = 1 等で公比が決まらない (旧は区間の端 50 / 1e-3 を返し下の公比の判定で棄却していた)
                    first = h0
                else:                                   # 前区間の最終間隔から公比で続ける。合わなければ一様
                    r = 1.0; first = L / n
                end = first * r ** (n - 1)
                if max(first, end) <= cap * (1 + 1e-9) and 1.0 / g - 1e-9 <= r <= g + 1e-9 and (last is None or max(first / last, last / first) <= g + 1e-9):
                    break
            else:
                raise ValueError("後流区間に間隔上限 %.3g を課す分割が見つからない (L=%.4g)" % (cap, L))
            out.append((n, 1.0 if n == 1 else r)); last = end; continue
        found = None
        for lim in (g_t, g):
            for n in range(1, int(L / (h0 if last is None else last / g)) + 4):
                if last is None:
                    try: r = prog_r(L, n, h0)
                    except ValueError: continue         # 同上 (n = 1)
                    first = h0
                elif n == 1:
                    r = L / last; first = L
                else:
                    f = lambda q: last * q * (geo_sum(q, n) if abs(q - 1) > 1e-12 else n)      # geo_sum: q^n があふれる所は inf (B4b-4)
                    if f(1.0 / g) > L or f(g) < L: continue
                    lo, hi = 1.0 / g, g
                    for _ in range(100):
                        q = 0.5 * (lo + hi); lo, hi = (q, hi) if f(q) < L else (lo, q)
                    r = 0.5 * (lo + hi); first = last * r
                if not (1.0 / g - 1e-9 <= r <= lim + 1e-9): continue
                end = first * (1.0 if n == 1 else r ** (n - 1))
                if end > 1.15 * hmax: continue
                found = (n, 1.0 if n == 1 else r, end); break
            if found: break
        if not found:                       # 物理 station の間隔が格子間隔より狭い等。一様で入れ、継ぎ目の比は報告側で検査する
            n = max(int(round(L / min(last, hmax))), 1); found = (n, 1.0, L / n)
        out.append(found[:2]); last = found[2]
    return out


# ---------------------------------------------------------------- 形状
def _out(x, v):
    return float(v) if np.ndim(x) == 0 else v


class Contour:
    """MOC 輪郭 (点列, 単位 m) の 3 次 B-spline 補間 y(x) (plan §4.14-3)。点列の範囲外は端点の接線で直線延長する
    (カウル後縁の外壁交点・カウル跡の上下境界。§6.2 接続規則)。x は狭義単調増加であること"""
    def __init__(s, xy):
        from scipy.interpolate import make_interp_spline
        xy = np.asarray(xy, float)
        if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 4: raise ValueError("輪郭は (N>=4, 2) の点列")
        if np.any(np.diff(xy[:, 0]) <= 0): raise ValueError("輪郭の x が狭義単調増加でない")
        s.xy = xy; s.xa, s.xb = float(xy[0, 0]), float(xy[-1, 0])
        s.sp = make_interp_spline(xy[:, 0], xy[:, 1], k=3); s.d1 = s.sp.derivative(); s.d2 = s.sp.derivative(2)
        s.ya, s.yb = float(s.sp(s.xa)), float(s.sp(s.xb)); s.sa, s.sb = float(s.d1(s.xa)), float(s.d1(s.xb))

    def y(s, x):
        xv = np.asarray(x, float); v = s.sp(np.clip(xv, s.xa, s.xb))
        v = np.where(xv < s.xa, s.ya + s.sa * (xv - s.xa), np.where(xv > s.xb, s.yb + s.sb * (xv - s.xb), v))
        return _out(x, v)

    def dy(s, x):
        xv = np.asarray(x, float); v = s.d1(np.clip(xv, s.xa, s.xb))
        v = np.where(xv < s.xa, s.sa, np.where(xv > s.xb, s.sb, v)); return _out(x, v)

    def ddy(s, x):
        xv = np.asarray(x, float); v = s.d2(np.clip(xv, s.xa, s.xb))
        v = np.where((xv < s.xa) | (xv > s.xb), 0.0, v); return _out(x, v)

    def project(s, pts, it=12):
        """点 (x, y) から曲線への最近点の x 座標 p と距離 (符号: 曲線より上 = 正)。延長部を含む"""
        P = np.asarray(pts, float).reshape(-1, 2); p = P[:, 0].copy()
        for _ in range(it):                     # Newton: f(p) = (p - x0) + (y(p) - y0) y'(p) = 0
            yy, d1, d2 = s.y(p), s.dy(p), s.ddy(p)
            f = (p - P[:, 0]) + (yy - P[:, 1]) * d1; df = 1.0 + d1 * d1 + (yy - P[:, 1]) * d2
            p = p - f / df
        d = np.hypot(P[:, 0] - p, P[:, 1] - s.y(p)) * np.sign(P[:, 1] - s.y(p))
        return p, d


def load_contours(d):
    """DIR/ramp_contour.csv と DIR/cowl_contour.csv (列 x_m,y_m、単位 m) を読む。runner_sern3d.prepare と同じ書式"""
    d = Path(d); out = {}
    for k in ("ramp", "cowl"):
        f = d / f"{k}_contour.csv"; raw = f.read_bytes()
        out[k] = np.loadtxt(f, delimiter=",", skiprows=1); out[k + "_src"] = dict(path=str(f.resolve()), sha256=hashlib.sha256(raw).hexdigest(), n=int(len(out[k])))
    return out


class Geom:
    def __init__(s, P, ramp_xy, cowl_xy):
        s.P = P; H = P["H"]; s.H = H
        for k in ("ZW", "TSW", "TC", "LSW", "LCOWL", "XEND", "DR", "YBOT", "ZFAR", "ZONE", "LWAKE"): setattr(s, k, P[k] * H)
        s.ZO = s.ZW + s.TSW
        s.prof = [(a * H, b * H, c * H) for a, b, c in PROFILE]
        s.ramp, s.cowl = Contour(ramp_xy), Contour(cowl_xy)
        # 参照輪郭の範囲と模型の寸法の整合 (黙って外挿・切り詰めをしない)
        if s.ramp.xa > 1e-12 * H or s.ramp.xb < s.XEND: raise ValueError("ランプ輪郭が模型の x 範囲 [0, XEND] を覆わない (%.4g..%.4g m)" % (s.ramp.xa, s.ramp.xb))
        if abs(s.cowl.xa) > 1e-12 * H: raise ValueError("カウル輪郭の始点 x が 0 でない (%.3g m)" % s.cowl.xa)
        if abs(s.cowl.xb - s.LCOWL) > 1e-9 * H: raise ValueError("カウル輪郭の終点 x %.9g m が L_cowl %.9g m と一致しない" % (s.cowl.xb, s.LCOWL))
        xx = np.linspace(0.0, s.LCOWL, 4001); kap = np.abs(s.cowl.ddy(xx)) / (1 + s.cowl.dy(xx) ** 2) ** 1.5
        if kap.max() * s.TC >= 0.5: raise ValueError("板厚がカウル内壁の曲率半径に対して厚すぎる (κ t = %.3g)" % (kap.max() * s.TC))
        s.p_te = float(s._po(np.array([s.LCOWL]))[0]); s.s_te = float(s.cowl.dy(s.p_te))
        s.yo_te = float(s.cowl.y(s.p_te) - s.TC / math.sqrt(1 + s.s_te ** 2))     # 外壁オフセット曲線と平面 x = L_cowl の交点の y

    def yr(s, x): return s.ramp.y(x)
    def dyr(s, x): return s.ramp.dy(x)
    def yc(s, x): return s.cowl.y(x)                 # カウル内壁 (流れ側) = cowl_xy
    def dyc(s, x): return s.cowl.dy(x)

    def _po(s, x):
        """外壁オフセット曲線上で x 座標が x になる点の、内壁側の媒介変数 p (内壁点 (p, y_c(p)) を法線方向に t ずらすと x)"""
        x = np.asarray(x, float); p = x.copy()
        for _ in range(60):
            sl = s.cowl.dy(p); p = x - s.TC * sl / np.sqrt(1 + sl * sl)
        return p

    def yo(s, x):
        """カウル外壁 (cowl_out) の y。x > L_cowl は外壁交点から接線方向の直線延長"""
        xv = np.asarray(x, float); p = s._po(np.minimum(xv, s.LCOWL)); sl = s.cowl.dy(p)
        v = s.cowl.y(p) - s.TC / np.sqrt(1 + sl * sl)
        v = np.where(xv > s.LCOWL, s.yo_te + s.s_te * (xv - s.LCOWL), v); return _out(x, v)

    def dyo(s, x):
        xv = np.asarray(x, float); v = s.cowl.dy(s._po(np.minimum(xv, s.LCOWL))); v = np.where(xv > s.LCOWL, s.s_te, v); return _out(x, v)

    def r(s, x):
        pr = s.prof
        if x <= pr[0][0]: return pr[0][1], pr[0][2]
        for (xa, ua, la), (xb, ub, lb) in zip(pr[:-1], pr[1:]):
            if x <= xb:
                w = smooth((x - xa) / (xb - xa)); return ua + w * (ub - ua), la + w * (lb - la)
        return pr[-1][1], pr[-1][2]

    def section(s, x):
        """断面の点 (y,z) と 3 本の壁曲線 (密な折れ線)。角度は z 軸から y 軸へ測る"""
        ru, rl = s.r(x); yr, yc, yo, ZW, DR = s.yr(x), s.yc(x), s.yo(x), s.ZW, s.DR
        arc = lambda cy, cz, r, a0, a1: np.stack([cy + r * np.sin(np.linspace(a0, a1, 200)), cz + r * np.cos(np.linspace(a0, a1, 200))], 1)
        if ru > 0: up1, up2 = arc(yr - ru, ZW - ru, ru, math.pi / 2, math.pi / 4), arc(yr - ru, ZW - ru, ru, math.pi / 4, 0.0)
        else: up1 = up2 = np.array([(yr, ZW)])
        if rl > 0: lo1, lo2 = arc(yc + rl, ZW - rl, rl, 0.0, -math.pi / 4), arc(yc + rl, ZW - rl, rl, -math.pi / 4, -math.pi / 2)
        else: lo1 = lo2 = np.array([(yc, ZW)])
        walls = dict(e1=np.vstack([[(yr, 0.0)], [(yr, ZW - ru)], up1]),
                     e5=np.vstack([up2, [(yr - ru, ZW)], [(yc + rl, ZW)], lo1]),
                     e8=np.vstack([lo2, [(yc, ZW - rl)], [(yc, 0.0)]]))
        # 内側線は壁を dr だけ内側へずらした線 (r > dr なら同心の円弧 r - dr、r <= dr なら角)。リングの法線方向の厚みが円弧上でも dr になる
        iu, il = max(ru - DR, 0.0), max(rl - DR, 0.0)
        if iu > 0: iu1, iu2 = arc(yr - ru, ZW - ru, iu, math.pi / 2, math.pi / 4), arc(yr - ru, ZW - ru, iu, math.pi / 4, 0.0)
        else: iu1 = iu2 = np.array([(yr - DR, ZW - DR)])
        if il > 0: il1, il2 = arc(yc + rl, ZW - rl, il, 0.0, -math.pi / 4), arc(yc + rl, ZW - rl, il, -math.pi / 4, -math.pi / 2)
        else: il1 = il2 = np.array([(yc + DR, ZW - DR)])
        walls.update(e3=np.vstack([[(yr - DR, 0.0)], [(yr - DR, ZW - max(ru, DR))], iu1])[::-1],
                     e7=np.vstack([iu2, [(yr - max(ru, DR), ZW - DR)], [(yc + max(rl, DR), ZW - DR)], il1])[::-1],
                     e10=np.vstack([il2, [(yc + DR, ZW - max(rl, DR))], [(yc + DR, 0.0)]])[::-1])
        pt = dict(PT0=(yr, 0.0), WU=tuple(up1[-1]), CUR=tuple(iu1[-1]), CUL=(yr - DR, 0.0),
                  WL=tuple(lo1[-1]), CLR=tuple(il1[-1]), CLL=(yc + DR, 0.0), PB0=(yc, 0.0))
        ys = [s.YBOT, yo - DR, yo, yc, yr]; zs = [0.0, ZW, s.ZO, s.ZO + DR, s.ZFAR]      # y 位置は x 依存 (外壁 yo・内壁 yc・ランプ yr)
        for j in range(5):
            for k in range(5):
                if (j, k) in ((3, 0), (3, 1), (4, 1), (4, 0)): continue
                pt["G%d%d" % (j, k)] = (ys[j], zs[k])
        return pt, walls, (ru, rl)


# 断面の辺: 名前 -> (始点, 終点)。hz30 = -e8、vy31 = -e5 は別名
N_EDGE = dict(e1=("PT0", "WU"), e2=("WU", "CUR"), e3=("CUR", "CUL"), e4=("CUL", "PT0"), e5=("WU", "WL"), e6=("WL", "CLR"),
              e7=("CLR", "CUR"), e8=("WL", "PB0"), e9=("PB0", "CLL"), e10=("CLL", "CLR"), e11=("CUL", "CLL"))
ALIAS = {"hz30": ("e8", -1), "vy31": ("e5", -1)}
GNAME = {(3, 0): "PB0", (3, 1): "WL", (4, 1): "WU", (4, 0): "PT0"}
gp = lambda j, k: GNAME.get((j, k), "G%d%d" % (j, k))
EDGE = dict(N_EDGE)
for j in range(5):
    for k in range(4):
        if (j, k) != (4, 0) and "hz%d%d" % (j, k) not in ALIAS: EDGE["hz%d%d" % (j, k)] = (gp(j, k), gp(j, k + 1))
for j in range(4):
    for k in range(5):
        if (j, k) != (3, 0) and "vy%d%d" % (j, k) not in ALIAS: EDGE["vy%d%d" % (j, k)] = (gp(j, k), gp(j + 1, k))


def ref(name, sign=1):
    if name.startswith("-"): return ref(name[1:], -sign)
    if name in ALIAS: b, s_ = ALIAS[name]; return b, s_ * sign
    return name, sign


def grid_block(j, k): return ["hz%d%d" % (j, k), "vy%d%d" % (j, k + 1), "-hz%d%d" % (j + 1, k), "-vy%d%d" % (j, k)]


BLK = {"N_TOP": ["e1", "e2", "e3", "e4"], "N_SIDE": ["e5", "e6", "e7", "-e2"], "N_BOT": ["e8", "e9", "e10", "-e6"], "N_CORE": ["-e3", "-e7", "-e10", "-e11"],
       "SW": grid_block(3, 1), "CW1": grid_block(2, 0), "CW2": grid_block(2, 1),
       "U1a": grid_block(0, 0), "U1b": grid_block(1, 0), "U2a": grid_block(0, 1), "U2b": grid_block(1, 1)}
for k, c in ((2, "b"), (3, "c")):
    for j, b in ((0, "a"), (1, "b"), (2, "c"), (3, "de")): BLK["S%s_%s" % (c, b)] = grid_block(j, k)
BLK = {b: [ref(e) for e in loop] for b, loop in BLK.items()}
ALWAYS = [b for b in BLK if b not in ("SW", "CW1", "CW2")]
PRESENT = {"A": set(ALWAYS), "B": set(ALWAYS) | {"SW"}, "C": set(ALWAYS) | {"SW", "CW1", "CW2"}}
WALL_TAGS = ("ramp", "vehicle", "sidewall_in", "sidewall_out", "sidewall_end", "cowl_in", "cowl_out", "cowl_side", "cowl_base")


def edge_counts(cnt):
    """辺 -> 節点数 (節点数は辺の同値類ごとに 1 つ。laws の割り当てと同じ表)。cnt は counts の dict (NL_co・NL_so は壁帯の区間数)"""
    ec = {}
    for e in ("e1", "e3", "e8", "e10", "hz00", "hz10", "hz20"): ec[e] = cnt["NZ"]
    for e in ("e2", "e4", "e6", "e9"): ec[e] = cnt["NR"] + 1
    for e in ("e5", "e7", "e11", "vy32", "vy33", "vy34"): ec[e] = cnt["NY"]
    for j in range(5):
        if j <= 2: ec["hz%d1" % j] = cnt["NSW"]
        ec["hz%d2" % j] = cnt["NL_so"] + 1; ec["hz%d3" % j] = cnt["NFZ"] + 1
        ec["vy2%d" % j] = cnt["NTC"]; ec["vy1%d" % j] = cnt["NL_co"] + 1; ec["vy0%d" % j] = cnt["NFY"] + 1
    ec["hz31"] = ec["hz41"] = cnt["NSW"]
    return ec


def count_mesh(regions, ec):
    """接続模型の節点数・ヘキサ数を gmsh なしで数える (B4b-5 の事前見積もり)。全節点が station 面上にあるので
    節点数 = Σ_station (断面の点 + 辺の内部節点 + 面の内部節点)。regions は区間ごとの領域 (A/B/C) の列 (station 数 - 1 個)、ec は edge_counts"""
    nb_n = lambda b: (ec[BLK[b][0][0]], ec[BLK[b][1][0]])
    nodes = hexes = 0; ns = len(regions) + 1
    for si in range(ns):
        nb = (PRESENT[regions[si - 1]] if si > 0 else set()) | (PRESENT[regions[si]] if si < ns - 1 else set())
        ed = {e for b in nb for e, _ in BLK[b]}
        nodes += len({p for e in ed for p in EDGE[e]}) + sum(ec[e] - 2 for e in ed) + sum((ni - 2) * (nj - 2) for ni, nj in map(nb_n, nb))
    for rg in regions: hexes += sum((ni - 1) * (nj - 1) for ni, nj in map(nb_n, PRESENT[rg]))
    return nodes, hexes


def xface_tag(e, reg):
    """参照数 1 の x 面 (辺 e を x に掃引した面) のタグ。表に無ければ None = 失敗"""
    if e == "e1" or e == "hz41": return "ramp"
    if e in ("hz42", "hz43"): return "vehicle"
    if e == "e5": return "sidewall_in" if reg == "A" else None
    if e == "vy32": return "sidewall_out" if reg == "A" else None
    if e == "e8": return "cowl_in" if reg in "AB" else None
    if e == "hz31": return "cowl_in" if reg == "B" else None
    if e in ("hz20", "hz21"): return "cowl_out" if reg in "AB" else None
    if e == "vy22": return "cowl_side" if reg in "AB" else None
    if e in ("e4", "e11", "e9", "vy00", "vy10") or (e == "vy20" and reg == "C"): return "symmetry"
    if e.startswith("hz0"): return "far_bottom"
    if e.startswith("vy") and e.endswith("4"): return "far_side"
    return None


# ---------------------------------------------------------------- 本体
def build(out, P, contours, quiet=False):
    import time
    t0 = time.time()
    if P.get("SECX") is not None:           # 全長診断: 模型の終端 XEND を参照ランプ輪郭の終端 (既定) か SECX_END まで延ばし、station と分割数をその全域から現行規則で決め直す
        P = dict(P); xe = float(contours["ramp"][-1, 0]) / P["H"]
        while xe * P["H"] > float(contours["ramp"][-1, 0]): xe = float(np.nextafter(xe, 0.0))      # 丸めで輪郭の外に出ない (Geom の範囲検査)
        P["XEND"] = min(float(P["SECX_END"]), xe) if P.get("SECX_END") is not None else xe
    G = Geom(P, contours["ramp"], contours["cowl"]); H = G.H; g = P["G"]; DR = G.DR
    # B4b-5: 壁別の第一内部点距離 (m)。全壁が等しければ旧の h1 / h1e と同じ値・同じ演算になる
    hw = {t: v * H for t, v in wall_h1(P).items()}
    h_r, h_v, h_si, h_so, h_ci, h_co = (hw[t] for t in ("ramp", "vehicle", "sidewall_in", "sidewall_out", "cowl_in", "cowl_out"))
    h_ring = min(h_r, h_si, h_ci)           # リングの 3 壁 (ramp・sidewall_in・cowl_in) の最小: リング横断 NR・z の NZ の下限を決める
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    bump = Bump()
    # ---- 物理 station (すべて物理位置) と物理区間。後流区間の終端 L + LWAKE も物理 station にする
    xm = 0.5 * (G.LSW + G.LCOWL)
    if not (G.LSW + G.LWAKE < xm and G.LCOWL + G.LWAKE < G.XEND): raise ValueError("後流区間 LWAKE が区間に収まらない")
    xs = {p[0] for p in G.prof} | {G.LSW, G.LCOWL, G.XEND, xm, G.LSW - G.ZONE, G.LCOWL + G.ZONE, G.LSW + G.LWAKE, G.LCOWL + G.LWAKE}
    for (xa, ua, la), (xb, ub, lb) in zip(G.prof[:-1], G.prof[1:]):
        if (ua, la) != (ub, lb):
            kink = False
            for ra, rb in ((ua, ub), (la, lb)):      # r が dr を横切る位置: 内側線の角が円弧に変わり CUR(x) が折れるので station にする
                if (ra - G.DR) * (rb - G.DR) < 0:
                    lo, hi = xa, xb
                    for _ in range(80):
                        xm_ = 0.5 * (lo + hi); rm = ra + smooth((xm_ - xa) / (xb - xa)) * (rb - ra)
                        lo, hi = (xm_, hi) if (rm - G.DR) * (ra - G.DR) > 0 else (lo, xm_)
                    xs.add(0.5 * (lo + hi)); kink = True
            if not kink: xs.add(0.5 * (xa + xb))   # 折れ点の station があれば中点は足さない (station が詰まると x 間隔が飛ぶ)          # 対角辺の第一間隔 s1 = h1 L/dr は r について連続なので、r = 0 の手前に特別な station は要らない
    xp = sorted(x for x in xs if -1e-12 <= x <= G.XEND + 1e-12)
    xp = [x for i, x in enumerate(xp) if i == 0 or x - xp[i - 1] > 1e-9]
    reg = lambda xa, xb: "A" if xb <= G.LSW + 1e-12 else ("B" if xb <= G.LCOWL + 1e-12 else "C")
    ivp = [(xp[i], xp[i + 1], reg(xp[i], xp[i + 1])) for i in range(len(xp) - 1)]
    # ---- x 分布: 端面から離れる向きに march。A は L_sw から上流へ、B は両端から中央へ、C は L_cowl から下流へ。
    # 後流区間 [L, L + LWAKE] は区間全体で Δx <= DXW (plan §6.2)、端面の第一間隔は h1e (B4b-5: L_sw 側は sidewall_end、L_cowl 側は cowl_base の値)
    hx = P["HX"] * H; g_t = 1.0 + 0.75 * (g - 1.0); dxw = P["DXW"] * H; h1e_sw, h1e_cb = hw["sidewall_end"], hw["cowl_base"]
    near = lambda a, b: abs(a - b) <= 1e-9 * H
    capf = lambda i: dxw if any(near(ivp[i][0], L_) and near(ivp[i][1], L_ + G.LWAKE) for L_ in (G.LSW, G.LCOWL)) else None
    idx = {r_: [i for i, v in enumerate(ivp) if v[2] == r_] for r_ in "ABC"}
    xlaw = {}
    def run(ids, toward_plus, hmax, h0):
        res = march([ivp[i][1] - ivp[i][0] for i in ids], h0, g_t, g, hmax, caps=[capf(i) for i in ids])
        for i, (n, r) in zip(ids, res): xlaw[i] = (n, r if toward_plus else -r)       # 負 = 逆向き (第一間隔が終端側)
    run(idx["A"][::-1], False, hx, h1e_sw)
    run([i for i in idx["B"] if ivp[i][1] <= xm + 1e-12], True, hx, h1e_sw); run([i for i in idx["B"] if ivp[i][0] >= xm - 1e-12][::-1], False, hx, h1e_cb)
    run(idx["C"], True, P["HX_FAR"] * H, h1e_cb)
    # ---- 格子の全 x 位置を station にする (各物理区間を march の分布で分割。物理 station の位置はそのまま)
    xs = [xp[0]]
    for i, (xa, xb, _) in enumerate(ivp):
        n, r = xlaw[i]; q = abs(r) if r > 0 else 1.0 / abs(r); w = q ** np.arange(n); cw = np.cumsum(w) / w.sum()
        xs += [xa + (xb - xa) * float(c_) for c_ in cw[:-1]] + [xb]
    dxs = np.diff(xs); x_ratio_design = float(np.maximum(dxs[1:] / dxs[:-1], dxs[:-1] / dxs[1:]).max())
    ivs = [(xs[i], xs[i + 1], reg(xs[i], xs[i + 1])) for i in range(len(xs) - 1)]
    # ---- 事前判定 (全 station と区間中点)
    assert 0 < G.LSW < G.LCOWL < G.XEND and G.TC > 0 and G.TSW > 0 and G.ZO + DR < G.ZFAR
    for x in xs + [0.5 * (a + b) for a, b, _ in ivs]:
        pt, _, (ru, rl) = G.section(x)
        if x >= G.LSW - 1e-12 and (ru > 0 or rl > 0): raise ValueError("x >= L_sw で r != 0")
        if ru + rl + 2 * DR >= G.yr(x) - G.yc(x) or max(ru, rl) + DR >= G.ZW: raise ValueError("フィレットが入らない x=%.4g" % x)
        if not G.YBOT < G.yo(x) - DR: raise ValueError("外部の下端 YBOT がカウル外壁の壁帯と交わる x=%.4g" % x)
        for e in ("e2", "e6"):
            if math.dist(pt[N_EDGE[e][0]], pt[N_EDGE[e][1]]) < 0.5 * DR: raise ValueError("コアが壁に近すぎる (%s, x=%.4g)" % (e, x))
    # ---- 断面の分布 (節点数は同値類ごとに 1 つ)
    # B4b-5: 期待層数は壁別 (plan §6.2 の式 prog_n(DR, h1, 1.2) + 4 を壁ごとの h1 で。端面も同式)。リング横断 NR は 3 壁の最大の層数 (= 最小の h1) から、
    # カウル外壁帯 vy1k は cowl_out、側壁外面帯 hz?2 は sidewall_out (= cowl_side) の層数で張る
    NLw = {t: prog_n(DR, hw[t], g) + 4 for t in hw}
    NL = prog_n(DR, h_ring, g) + 4                                       # リングの壁法線の区間数の基準 (plan §6.2: prog_n(DR, h1, 1.2) + 4)
    NL_co, NL_so = NLw["cowl_out"], NLw["sidewall_out"]
    xa_ = np.array(xs); Lmax = float(np.max(G.yr(xa_) - G.yc(xa_)))
    ccs = np.sqrt(1 + G.dyc(xa_) ** 2); LTmax = float(np.max(G.yc(xa_) - G.yo(xa_)))
    # 壁層セルの接線方向の幅の上限 (リングは最大 45° 傾くので、壁層でも AR <= 1000 に収める)。y の同値類 (e5・vy32) を接線の辺に持つ側壁の小さい方の h1
    tan_max = P["AR_TAN"] * min(h_si, h_so)
    NY = bump.fit_n(min(h_r, h_ci, h_v) / Lmax, g, hmax_rel=tan_max / Lmax); NSW = bump.fit_n(min(h_si, h_so) / G.TSW, g)
    NTC = bump.fit_n(min(h_ci, h_co) * float(ccs.min()) / LTmax, g)   # 板の y 方向 (鉛直厚 t/cos θ_c、両端 h1 cc)
    # z 方向 (z ∈ [0, W/2] の同値類 e1/e3/e8/e10/hz00/hz10/hz20) は継ぎ目 z = W/2 側へ片側細分 (plan §6.3 の B′ を全継ぎ目へ、§5.1 B4-2)。
    # 第一間隔はどれも h1 以上・長さは W/2 以下なので、h1 から成長率 g_t で W/2 を覆う区間数があれば全辺の公比 <= g_t (分割数の固定 NZ=31 は下限に)
    NZ = max(int(P["NZ"]), prog_n(G.ZW, h_ring, g_t) + 1)
    # リングの横断 (同値類 e2/e4/e6/e9): 対角辺 e2/e6 は両端細分 (Bump)。対角の内側端 CUR/CLR の間隔をコア (e3/e7/e10) と 1.2 以内でつなぐため
    # (等比のままだと内側端が 0.015 H になり、壁側の隅で h1 級に細分した e1/e5/e8 と同じリングの中で格子線が 75° 以上傾く)。
    # 増加側の区間数が壁の期待層数 NL 以上になるよう 2 NL 区間から始め、公比 <= g まで増やす
    cmax = float(np.max(np.sqrt(1 + G.dyr(xa_) ** 2))); s2_rel = h_ring * 0.5 * (1 + cmax) / DR      # 対角辺の第一間隔 / 辺長 (r = 0 で L = sqrt2 dr)
    NR = 2 * NL
    while True:
        cf_ = bump.coef(NR + 1, s2_rel); sp_ = bump.spacing(NR + 1, cf_)
        if np.maximum(sp_[1:] / sp_[:-1], sp_[:-1] / sp_[1:]).max() <= g_t: break
        NR += 2
    kv = P["KV"] * h_so
    # 遠方帯 vy0k (カウル外壁帯の外) と hz?3 (側壁外面帯の外) の区間数 (B4b-4): laws で実際に張る (帯長, 第一間隔 = 壁帯の最外間隔 × 公比) の
    # 全組について、実間隔列の最大が HFAR 以下 (far_n)。旧は概算の第一間隔と近似式 L (1 - 1/r) で、既定の実最大 0.109 H を合格にしていた
    def band_out(hh, nl):                   # 壁帯 (nl 区間・第一間隔 hh・全長 DR) の外の第一間隔 = 最外間隔 × 公比 (laws と同じ式)
        hh = min(hh, 0.999 * DR / nl); return hh * prog_r(DR, nl, hh) ** nl
    NFY, far_y = far_n([(float(G.yo(x_)) - DR - G.YBOT, band_out(h_co * math.sqrt(1 + G.dyo(x_) ** 2), NL_co)) for x_ in xs], P["HFAR"] * H, g_t)
    NFZ, far_z = far_n([(G.ZFAR - G.ZO - DR, band_out(kv if j >= 3 else h_so, NL_so)) for j in range(5)], P["HFAR"] * H, g_t)
    seams = {}                                                           # 継ぎ目の間隔の連動の記録 (station ごとの値。report に最小・最大を出す)
    # B4b-5: 両端の壁 (継ぎ目) の第一間隔が違う同値類は、両端の間隔を別々に持つ分布にする。全壁が等しいとき (既定) は旧の対称な Bump のまま (ビット同一)。
    # SW 列 (z ∈ [z_w, z_o] の hz?1)・CW 列 (板厚の y の vy2k) は両端分布 two_end_spacing (公比 <= g_two)。
    # 側壁帯の y (e5・vy32・vy33・vy34) は同じ対称 Bump を滑らかに曲げる y_class (薄い SW ブロックの両側 e5・vy32 の節点の y を揃えるため。
    # 別々の両端分布にすると頭打ちの位置がずれて SW の格子線が傾き、接続の試験片で継ぎ目の比 2.9〜3.7 になった)
    asym = dict(SW=not (h_so == h_si == h_r == h_ci), TC=h_ci != h_co, Y=not (h_r == h_v == h_ci))
    # 両端分布の公比の上限: 既存の Bump の当てはめ (公比 <= g) と同等の粗さにし、1.2 ちょうどの丸めで隣接比の判定 (> 1.2) を割らない余裕 0.8 % を残す
    g_two = 1.0 + 0.95 * (g - 1.0)

    def bl(n, Led, h):
        cf = bump.coef(n, h / Led); return (n, "Progression", 1.0) if cf is None else (n, "Bump", cf)

    def y_class(Ly, v, n):
        """側壁帯の y の同値類 (e5・vy32・vy33・vy34) の、両端の間隔が違うときの分布 (B4b-5)。全員を同じ基準分布 (端 a_bot の対称 Bump) から
        ring_diag_spacing で滑らかに曲げて上端だけを変える (下端 a_bot・区間数・全長は同じ)。e5 は vy32 の逆順そのもの (薄い SW ブロックの両側で節点の y を揃える)。
        戻り値: 辺 -> 間隔列 (各辺の始点から)。内部隣接比が g を超えれば ValueError"""
        n_, typ_, cf_ = bl(n, Ly, v["a_bot"])
        base = bump.spacing(n_, cf_) * Ly if typ_ == "Bump" else np.full(n_ - 1, Ly / (n_ - 1))
        d32 = ring_diag_spacing(base, v["top32"] / base[0], g); dv = ring_diag_spacing(base, v["a_v"] / base[0], g)
        return dict(vy32=d32, e5=d32[::-1].copy(), vy33=dv, vy34=dv)

    def seam_vals(x, Q):
        """station x の継ぎ目の間隔 (壁別 h1 から) と、壁別 h1 で解けない組合せ (station_conflicts) のリスト。laws と事前検査が共有する"""
        pt, elen = Q["pt"], Q["elen"]
        c = math.sqrt(1 + G.dyr(x) ** 2); ru, rl = G.r(x); rc = float(P.get("RING_CORE_RATIO", 1.0))
        cc = math.sqrt(1 + G.dyc(x) ** 2); co = math.sqrt(1 + G.dyo(x) ** 2)       # カウル内壁・外壁の傾き (鉛直の間隔 = 法線距離 × c)
        # gmsh の transfinite 補間は、層 1 を「辺長に対する比 eta = s1/L」で置く (X = (1-eta) 壁 + eta 内側線、eta は両側辺の間で線形)。
        # 内側線は平面壁から dr だけ離れた平行線なので、**リングの全側辺で eta = h1/dr に揃えれば平面部の法線距離は厳密に h1**。
        # 対角辺 e2/e6 は長さが sqrt2 dr - (sqrt2-1) r なので s1 = h1 L/dr (r=0 で sqrt2 h1 = 接する両壁の条件と同じ)。
        # フィレット円弧上は内側線までの法線距離が dr より長いぶん厚くなる (45° 点で L_e2/dr 倍。薄くはならない)
        def tilt(e, wp):                      # 45° 点の壁面法線の x 方向の傾き a_n = d(壁点)/dx . d (r'(x) とランプ勾配で生じる)。法線距離は 1/sqrt(1+a_n^2) 倍に縮む
            a_, b_ = N_EDGE[e]; d = np.array(pt[b_]) - np.array(pt[a_]); d /= np.linalg.norm(d); dx = 1e-5 * H
            x0, x1 = max(x - dx, 0.0), min(x + dx, G.XEND); w_ = lambda xx: np.array(G.section(xx)[0][wp])
            return math.sqrt(1 + (float((w_(x1) - w_(x0)) @ d) / (x1 - x0)) ** 2)
        # 傾きの補正は 3 割だけ掛ける: 全部掛けると e2 側の eta が増えて平面部が厚くなる (円弧上 -9 % -> -6 %、平面部 +3 % 以内)
        k2 = (1 + 0.3 * (tilt("e2", "WU") - 1)) if ru > 0 else 1.0; k6 = (1 + 0.3 * (tilt("e6", "WL") - 1)) if rl > 0 else 1.0
        # B4b-5: 隅の対角 e2/e6 の壁側端は接する 2 壁 (sidewall_in と ramp / cowl_in) の h1 から (diag_first)、対称面側 e4/e9 はランプ / カウル内壁の h1
        s = dict(e2=diag_first(h_si, h_r, c, k2, elen["e2"], DR), e6=diag_first(h_si, h_ci, cc, k6, elen["e6"], DR), e4=h_r * c, e9=h_ci * cc)
        # e5 (と SW・S の de の y) の端: 旧はランプ側とカウル側の傾きの平均 0.5 h1 (c + cc) E5K (B・C で露出壁の第一層)。壁別には端ごとにその壁の h1:
        # 上端 (y_r) はランプ a_top・機体下面 a_v、下端 (y_c) はカウル内壁 a_bot (全壁が等しければ 3 つとも旧の hend と同じ値)
        a_top, a_bot, a_v = (0.5 * h_ * (c + cc) * P["E5K"] for h_ in (h_r, h_ci, h_v))
        sw_present = "SW" in Q["nb"]
        sv = dict(c=c, cc=cc, co=co, s=s, a_top=a_top, a_bot=a_bot, a_v=a_v, sw=sw_present, LT=G.yc(x) - G.yo(x), LFY=G.yo(x) - DR - G.YBOT)
        issues = station_conflicts(x / H, c, cc, s, elen, DR, hw, a_top, a_v, sw_present)
        # vy32 の上端 (z = z_o, y = y_r): B・C (SW あり) はランプ (SW) と機体下面 (Sb_de) の継ぎ目で、1 つの間隔を両方の目標に合わせる
        # (幾何平均。station_conflicts で ±10 % を判定。全壁が等しければ a_top そのもの)。e5 の上端も同じ値にする: SW は幅 t_sw の薄いブロックなので
        # 両側の辺 e5・vy32 の節点の y がずれると格子線が大きく傾く。x の連続性のため A (SW なし) でも同じ規則
        sv["top32"] = a_top if a_top == a_v else math.sqrt(a_top * a_v)
        # リング横断: 対角 e2/e6 は両端 s (Bump)、平面の対称面側 e4/e9 は壁から等比 (NR 区間)
        Lr = {e: bl(NR + 1, elen[e], s[e]) for e in ("e2", "e4", "e6", "e9")}
        # B4-4 (plan §6.4): RING_CORE_RATIO != 1 なら e2/e6 は壁側端 (WU/WL) の間隔を保持し、コア側端 (CUR/CLR) だけを比 rc にする非対称分布。
        # 区間数 NR・全長は同じ。gmsh の則では表せないので "Nodes" (節点の相対位置の列) で渡し、generate(1) の後に節点を移す (apply_node_laws)
        sc = dict(e2=s["e2"], e6=s["e6"])                                # コア側端の間隔
        if abs(rc - 1.0) > 1e-12:
            for e in ("e2", "e6"):
                n_, typ_, cf_ = Lr[e]
                base = bump.spacing(n_, cf_) * elen[e] if typ_ == "Bump" else np.full(n_ - 1, elen[e] / (n_ - 1))
                d_ = ring_diag_spacing(base, rc, g); Lr[e] = (n_, "Nodes", tuple(np.r_[0.0, np.cumsum(d_)] / d_.sum())); sc[e] = float(d_[-1])
        sv.update(Lr=Lr, sc=sc)
        try:
            # 隅 WU/WL と内側の隅 CUR/CLR: e1@WU・e5@WU (TOP|SIDE)、e3@CUR・e7@CUR・e2@CUR (TOP|SIDE・TOP|CORE・SIDE|CORE)、WL/CLR も同様。
            # リングの中で格子線を 45° 以内に保つため e1/e3 (e8/e10, e7) の第一間隔を 1 つの値 F にする: F は e5 の両端 (a_top・a_bot) と対角の端 s2/s6 と 1.2 以内
            # N_SIDE の中ほどのリング外縁の間隔 (対角 e2/e6 の端の法線成分の平均 = (s2 + s6)/(2 sqrt2)) も N_CORE の z 間隔 (e3/e10 の端 F) と接する
            F = match_spacing([a_top, a_bot, s["e2"], s["e6"], (s["e2"] + s["e6"]) / (2 * math.sqrt(2))], g)
            # コア側 (e3@CUR・e7・e10@CLR) の端 Fc: B4-4 では対角のコア側端 sc から同じ式で連動させる (rc = 1 なら Fc = F で A と同一)
            Fc = F if abs(rc - 1.0) <= 1e-12 else match_spacing([sc["e2"], sc["e6"], (sc["e2"] + sc["e6"]) / (2 * math.sqrt(2))], g)
            o11 = match_spacing([s["e4"], s["e9"]], g)                   # コアの対称面側 e11: リング外縁 o4 (CUL)・o9 (CLL) に合わせる
            a3, a4 = sw_end_spacing(s["e2"], s["e6"], F, kv, g, int(P.get("SW_END_RULE", 0)))
            # z ∈ [W/2, W/2 + t_sw]: U2/CW2 の列 hz01/hz11/hz21 は z = W/2 側 F と z = z_o 側 hz?2 (h_so) に連動 (対称なら両端 b21)
            b21 = None if asym["SW"] else match_spacing([F, h_so], g)
            # 板厚 vy21 (C の CW1|CW2) の上端は e6@WL・e5@WL に、下端は vy11@G21 (h_co co) に連動 (対称なら両端 v21)
            v21 = match_spacing([s["e6"], a_bot, h_co * co], g) if not asym["TC"] else None
            v21t = match_spacing([s["e6"], a_bot], g) if asym["TC"] else None
        except ValueError as e_:
            issues.append(dict(kind="seam", where="継ぎ目の間隔の連動 (match_spacing)", x=x / H, reason=str(e_))); return sv, issues
        sv.update(F=F, Fc=Fc, o11=o11, a3=a3, a4=a4, b21=b21, v21=v21, v21t=v21t)
        return sv, issues

    def conflict_text(iss):
        """解けない組合せの一覧を、種類・位置ごとに x の範囲と比の範囲へまとめた文字列にする"""
        grp = {}
        for it in iss: grp.setdefault((it["kind"], it["where"]), []).append(it)
        out = []
        for (kd, wh), v in grp.items():
            xr = (min(i["x"] for i in v), max(i["x"] for i in v))
            rk = [k for k in v[0] if k.startswith("ratio")]
            rr = ", ".join("%s %.4g..%.4g" % (k, min(i[k] for i in v), max(i[k] for i in v)) for k in rk)
            hv = (" h1 %s µm" % "/".join("%.3g" % u for u in v[0]["h1_um"])) if "h1_um" in v[0] else ""
            out.append("[%s] %s: x/H %.4g..%.4g (%d station)%s %s (許容 ±%s) — %s" % (kd, wh, xr[0], xr[1], len(v), hv, rr, v[0].get("tol", "-"), v[0]["reason"]))
        return "\n  ".join(out)

    def laws(x, Q):
        """station x での各辺の (節点数, 種別, 係数)。
        継ぎ目の規則 (plan §5.1 B4-2): 共有辺を挟む 2 ブロックの「共有辺の端から出る辺」の第一間隔を比 g 以内に揃える。
        端点ごとの拘束は match_spacing で解き (両隣との比が等しくなる幾何平均)、解けなければ ValueError で止める。
        gmsh の transfinite は同値類ごとに節点数 1 つ・辺ごとに分布 1 つなので、継ぎ目を挟む両側の第一/最終間隔を同じ値から作る。
        B4b-5: 壁別 h1 で解けない組合せ (station_conflicts) があれば止める (黙って h1 を変えない)"""
        sv, issues = Q["sv"] if "sv" in Q else seam_vals(x, Q)
        if issues: raise ValueError("壁別 h1 で解けない組合せ:\n  " + conflict_text(issues))
        elen = Q["elen"]; L = dict(sv["Lr"]); s, sc, c, cc, co = sv["s"], sv["sc"], sv["c"], sv["cc"], sv["co"]
        a_top, a_bot, a_v, F, Fc, a3, a4, LT, LFY = (sv[k] for k in ("a_top", "a_bot", "a_v", "F", "Fc", "a3", "a4", "LT", "LFY"))
        L["e7"] = bl(NY, elen["e7"], Fc)                                 # コアの側壁側: 内側の隅 CUR/CLR で e3/e10 と同じ Fc
        L["e11"] = bl(NY, elen["e11"], sv["o11"])
        if not asym["Y"]:                       # 全壁が等しい: 旧と同じ対称な Bump (両端 hend = a_top = a_bot = a_v)
            L["e5"] = bl(NY, elen["e5"], a_top)
            for k in (2, 3, 4): L["vy3%d" % k] = bl(NY, elen["vy32"], a_top)
        else:                                   # e5 は WU (ランプ) -> WL (カウル内壁)、vy3k は y_c (カウル内壁側) -> y_r (ランプ / 機体下面)
            for e, d_ in y_class(elen["vy32"], sv, NY).items(): L[e] = (NY, "Nodes", tuple(np.r_[0.0, np.cumsum(d_)] / d_.sum()))
        # z ∈ [0, W/2]: 継ぎ目 z = W/2 側 (e3/e10 は内側の隅) へ片側細分。**同値類の全辺で継ぎ目側の第一間隔を F に揃える**
        # (CW1 は板厚 0.005 H の薄いブロックなので、上下の辺 e8・hz20 の分布が違うと格子線が大きく傾く。U1b・U1a も同じ分布で追従)
        zl = lambda e: elen.get(e, G.ZW)
        L["e1"] = (NZ, "Progression", -prog_r(zl("e1"), NZ - 1, F)); L["e3"] = (NZ, "Progression", prog_r(zl("e3"), NZ - 1, Fc))
        L["e8"] = (NZ, "Progression", prog_r(zl("e8"), NZ - 1, F)); L["e10"] = (NZ, "Progression", -prog_r(zl("e10"), NZ - 1, Fc))
        for j in (0, 1, 2): L["hz%d0" % j] = (NZ, "Progression", -prog_r(G.ZW, NZ - 1, F))
        # z ∈ [W/2, W/2 + t_sw] (SW 列): hz01/hz11/hz21 は z = W/2 側 F・z = z_o 側 h_so (sidewall_out = cowl_side)。
        # SW の hz31 (WL)・hz41 (WU) は z = W/2 側が e6/e2 の第一間隔からの a3/a4 (sw_end_spacing)、z = z_o 側が側壁外面帯 hz32/hz42 の第一間隔 kv。
        # 全壁が等しければ旧の対称な Bump (b21・a3・a4)、sidewall_in と sidewall_out の h1 が違えば両端の間隔を別に持つ分布
        if not asym["SW"]:
            for j in (0, 1, 2): L["hz%d1" % j] = bl(NSW, G.TSW, sv["b21"])
            L["hz31"] = bl(NSW, G.TSW, a3); L["hz41"] = bl(NSW, G.TSW, a4)
        else:
            for j in (0, 1, 2): L["hz%d1" % j] = two_end_law(G.TSW, NSW, F, h_so, g_two)
            L["hz31"] = two_end_law(G.TSW, NSW, a3, kv, g_two); L["hz41"] = two_end_law(G.TSW, NSW, a4, kv, g_two)
        # y: カウル外壁帯 vy1k (NL_co 区間、壁 = yo 側 h_co co)、板厚 vy2k (CW 列: 上 = cowl_in、下 = cowl_out)。vy21 (C の CW1|CW2) は e6@WL・e5@WL・vy11@G21 に連動
        cap_co = lambda hh: min(hh, 0.999 * DR / NL_co)
        for k in range(5):
            if not asym["TC"]: L["vy2%d" % k] = bl(NTC, LT, sv["v21"] if k == 1 else h_ci * cc)
            else: L["vy2%d" % k] = two_end_law(LT, NTC, h_co * co, sv["v21t"] if k == 1 else h_ci * cc, g_two)      # 始点 = y_o (cowl_out 側)
            hh = cap_co(h_co * co); r_ = prog_r(DR, NL_co, hh); L["vy1%d" % k] = (NL_co + 1, "Progression", -r_)
            L["vy0%d" % k] = (NFY + 1, "Progression", -prog_r(LFY, NFY, hh * r_ ** NL_co))
        cap_so = lambda hh: min(hh, 0.999 * DR / NL_so)
        for j in range(5):                       # z = z_o から外へ: 側壁外面帯 hz{j}2 (NL_so 区間、第一間隔 h_so)。hz32・hz42 は kv (SW の分布との連動)
            hh = cap_so(kv if j >= 3 else h_so); r_ = prog_r(DR, NL_so, hh); L["hz%d2" % j] = (NL_so + 1, "Progression", r_)
            L["hz%d3" % j] = (NFZ + 1, "Progression", prog_r(G.ZFAR - G.ZO - DR, NFZ, hh * r_ ** NL_so))
        hr_ = h_ring                            # 記録の基準 (旧の h1。全壁が等しければ同じ値)
        rec = dict(F_over_h1=F / hr_, Fc_over_h1=Fc / hr_, s2c_over_h1=sc["e2"] / hr_, s6c_over_h1=sc["e6"] / hr_, a3_over_h1=a3 / hr_, a4_over_h1=a4 / hr_,
                   s2_over_h1=s["e2"] / hr_, s6_over_h1=s["e6"] / hr_, hend_over_h1=a_top / hr_, o4_over_h1=s["e4"] / hr_, o9_over_h1=s["e9"] / hr_)
        if sv["b21"] is not None: rec["b21_over_h1"] = sv["b21"] / hr_
        if sv["v21"] is not None: rec["v21_over_h1"] = sv["v21"] / hr_
        if asym["Y"]: rec.update(a_bot_over_h1=a_bot / hr_, a_v_over_h1=a_v / hr_, top32_over_h1=sv["top32"] / hr_)
        for k_, v_ in rec.items():
            q = seams.setdefault(k_, [math.inf, -math.inf]); q[0] = min(q[0], v_); q[1] = max(q[1], v_)
        return L

    # ---- 断面の前計算 (gmsh の較正モデルは本体モデルを作る前に使い切る)
    need_blk = lambda si: (PRESENT[ivs[si - 1][2]] if si > 0 else set()) | (PRESENT[ivs[si][2]] if si < len(ivs) else set())
    def sec_geo(x, nb):
        """位置 x の断面の幾何 (点・辺・壁曲線・辺長)。nb はその断面に面を張るブロックの集合"""
        pt, walls, _ = G.section(x)
        ed = sorted({e for b in nb for e, _ in BLK[b]}); elen, wr = {}, {}
        for e in ed:
            a_, b_ = EDGE[e]
            if e in walls:
                W = np.asarray(walls[e]); dd = np.r_[0, np.cumsum(np.linalg.norm(np.diff(W, axis=0), axis=1))]
                W = W[np.r_[True, np.diff(dd) > 1e-14]]; dd = np.r_[0, np.cumsum(np.linalg.norm(np.diff(W, axis=0), axis=1))]
                ch = W[-1] - W[0]; off = np.abs(ch[0] * (W[:, 1] - W[0, 1]) - ch[1] * (W[:, 0] - W[0, 0])) / max(np.linalg.norm(ch), 1e-300)
                if off.max() <= 1e-12 * H: elen[e] = math.dist(pt[a_], pt[b_]); continue       # 直線の壁 (半径 0) は Line で張る (節点が壁の直線上に厳密に乗る)
                sv_ = np.linspace(0, dd[-1], int(P["M_SPL"])); wr[e] = np.stack([np.interp(sv_, dd, W[:, q]) for q in (0, 1)], 1); elen[e] = dd[-1]
            else: elen[e] = math.dist(pt[a_], pt[b_])
        for e in ("e2", "e5", "e6", "e7", "e11", "vy32"): elen.setdefault(e, math.dist(pt[EDGE[e][0]], pt[EDGE[e][1]]))
        return dict(pt=pt, nb=nb, ed=ed, elen=elen, wr=wr)
    def sec_at(x, nb):
        """位置 x の断面 (点・辺・壁曲線・各辺の分布則)"""
        Q = sec_geo(x, nb); Q["law"] = laws(x, Q); return Q
    # ---- B4b-5: 壁別 h1 の事前検査 (全 station。解けない組合せは位置・値・理由を付けて止める) と、両端の間隔が違う同値類の節点数
    SECG = [sec_geo(x, need_blk(si)) for si, x in enumerate(xs)]
    conflicts = []
    for x, Q in zip(xs, SECG):
        Q["sv"] = seam_vals(x, Q); conflicts += Q["sv"][1]
    asym_ranges = {}
    def class_n(cls, items, n_min):
        """両端の間隔が違う同値類の区間数: 全 station・全辺の (長さ, 始端, 終端) で two_end_n_range の共通部分の最小 (n_min 以上)。無ければ conflicts へ"""
        rng = [two_end_n_range(L_, a_, b_, g_two) for L_, a_, b_ in items]
        if any(r_ is None for r_ in rng) or not rng:
            conflicts.append(dict(kind="count", where=cls, x=float("nan"), reason="両端の間隔が違う分布の区間数が取れない (端の比が大きすぎる、または長さに対して端の間隔が大きすぎる)")); return n_min
        lo, hi = max(max(r_[0] for r_ in rng), n_min), min(r_[1] for r_ in rng)
        asym_ranges[cls] = (lo, hi)
        if lo > hi: conflicts.append(dict(kind="count", where=cls, x=float("nan"), reason="両端の間隔が違う分布の区間数の範囲が全 station で交わらない (下限 %d > 上限 %d)" % (lo, hi)))
        return lo
    ok_sv = [(x, Q["sv"][0]) for x, Q in zip(xs, SECG) if "F" in Q["sv"][0]]
    if asym["SW"] and ok_sv:
        NSW = class_n("SW 列 (hz?1)", [(G.TSW, v["F"], h_so) for _, v in ok_sv] + [(G.TSW, v[k], kv) for _, v in ok_sv for k in ("a3", "a4")], 1) + 1
    if asym["TC"] and ok_sv:
        NTC = class_n("CW 列 (vy2k)", [(v["LT"], h_co * v["co"], h_ci * v["cc"]) for _, v in ok_sv] + [(v["LT"], h_co * v["co"], v["v21t"]) for _, v in ok_sv], 1) + 1
    if asym["Y"] and ok_sv:                 # 側壁帯の y: 全 station で y_class が解け (内部隣接比 <= g)、側壁の線の最大間隔が AR の上限 tan_max 以下になる最小の NY
        ly = lambda x_: float(G.yr(x_) - G.yc(x_)); ny0 = NY
        while True:
            try:
                if max(float(d_.max()) for x_, v in ok_sv for d_ in y_class(ly(x_), v, NY).values()) <= tan_max: break
            except ValueError:
                pass
            NY += 1
            if NY > 3 * ny0:
                conflicts.append(dict(kind="count", where="側壁帯の y (e5・vy3k)", x=float("nan"), reason="両端の間隔が違う分布 (Bump を曲げる) が NY %d..%d で解けない" % (ny0, NY))); break
        asym_ranges["側壁帯の y (e5・vy3k)"] = (ny0 - 1, NY - 1)
    counts = dict(NL=NL, NL_co=NL_co, NL_so=NL_so, NL_END={t: NLw[t] for t in WALL_ENDS}, NR=NR, NY=NY, NZ=NZ, NSW=NSW, NTC=NTC, NFY=NFY, NFZ=NFZ)
    wall_info = dict(wall_h1={t: hw[t] / H for t in hw}, wall_h1_um={t: hw[t] * 1e6 for t in hw}, wall_layers=NLw, asymmetric_classes=asym,
                     asym_ranges={k: list(v) for k, v in asym_ranges.items()}, preset=P.get("PRESET"))
    if P.get("ESTIMATE"):                   # B4b-5: 節点数・ヘキサ数の事前見積もり (gmsh の格子を張らない)。解けない組合せも一覧で返す
        nodes, hexes = count_mesh([v[2] for v in ivs], edge_counts(counts))
        gmsh.finalize()
        return dict(VERDICT="ESTIMATE", nodes=nodes, hexes=hexes, stations=len(xs), counts=counts, **wall_info,
                    far_spacing_max=dict(y=far_y / H, z=far_z / H, limit=P["HFAR"]), conflicts=conflicts,
                    conflicts_text=conflict_text([c_ for c_ in conflicts if c_["kind"] != "count"]) if conflicts else "",
                    nx=[abs(xlaw[i][0]) for i in range(len(ivp))], phys_regions=[v[2] for v in ivp], time_s=time.time() - t0)
    if conflicts:
        gmsh.finalize()
        cnt_ = [c_ for c_ in conflicts if c_["kind"] == "count"]
        raise ValueError("壁別 h1 で解けない組合せ (生成しない):\n  " + conflict_text([c_ for c_ in conflicts if c_["kind"] != "count"])
                         + "".join("\n  [count] %s: %s" % (c_["where"], c_["reason"]) for c_ in cnt_))
    if P.get("SECX") is not None:           # 全長診断 (plan §6.5 末尾): 指定 x (/H) の断面を station に丸めずその位置で作る (分割数は延長した全域の station から決めたもの)
        items = [(xv, x, sec_at(x, set().union(*[PRESENT[r_] for r_ in rg]))) for xv, x, rg in section_targets(np.atleast_1d(P["SECX"]), ivp, H)]
        res = section_2d(G, items, seams, quiet)
        js = np.array(xs) / H; ii = np.searchsorted(js, [it[0] for it in items])
        res.update(counts=counts, **wall_info, domain_x_end=G.XEND / H, n_stations=len(xs),
                   neighbor_stations={"%.6g" % it[0]: [float(js[max(i_ - 1, 0)]), float(js[min(i_, len(js) - 1)])] for it, i_ in zip(items, ii)},
                   note="XEND (模型 1.5 H) を超える x は現行の延長則 (カウル跡の上下 = 後縁から接線延長、ランプは参照輪郭) で作った断面。延長則は暫定 "
                        "(B4a(1) の接続表で確定させる)。分割数 (NY・NTC・NR・NFY・NFZ 等) は x ∈ [0, domain_x_end] の全 station から現行規則で決めた値")
        return res
    for x, Q in zip(xs, SECG): Q["law"] = laws(x, Q)
    SEC = SECG
    coupon = P.get("COUPON_N") is not None
    if coupon:                              # B4-4 (plan §6.4): 全模型と同じ station・同じ分布のうち、COUPON_X に最も近い station を含む COUPON_N 個だけを押し出す
        k_ = int(P["COUPON_N"]); i_ = int(np.argmin(np.abs(np.array(xs) - float(P["COUPON_X"]) * H)))
        i0 = max(0, min(i_ - k_ + 1 if i_ == len(xs) - 1 else i_ - k_ // 2, len(xs) - k_)); sel = list(range(i0, i0 + k_))
        xs = [xs[i] for i in sel]; SEC = [SEC[i] for i in sel]
        ivs = [(xs[i], xs[i + 1], reg(xs[i], xs[i + 1])) for i in range(len(xs) - 1)]
        for si, Q in enumerate(SEC):        # 試験片の両端は片側の区間のブロックだけ (余分な面を作らない)
            Q["nb"] = need_blk(si); Q["ed"] = sorted({e for b in Q["nb"] for e, _ in BLK[b]})
            Q["wr"] = {e: v for e, v in Q["wr"].items() if e in Q["ed"]}
    if P.get("SEC2D") is not None:          # 分布の調整用: 指定 x (/H) に最も近い station の断面だけを 2D で切り、継ぎ目の比と skew を出す (判定には使わない)
        items = []
        for xv in [float(v) for v in np.atleast_1d(P["SEC2D"])]:
            si = int(np.argmin(np.abs(np.array(xs) - xv * H))); items.append((xv, xs[si], SEC[si]))
        return section_2d(G, items, seams, quiet)
    # ---- gmsh モデル
    gmsh.model.add("junction"); ge = gmsh.model.geo; gm = ge.mesh
    S = []
    for si, x in enumerate(xs):
        Q = SEC[si]; pt, ed = Q["pt"], Q["ed"]; pn = sorted({p for e in ed for p in EDGE[e]})
        pid = {p: ge.addPoint(x, pt[p][0], pt[p][1]) for p in pn}; eid = {}
        for e in ed:
            a_, b_ = EDGE[e]
            if e in Q["wr"]: eid[e] = ge.addSpline([pid[a_]] + [ge.addPoint(x, p_[0], p_[1]) for p_ in Q["wr"][e][1:-1]] + [pid[b_]])
            else: eid[e] = ge.addLine(pid[a_], pid[b_])
        face = {b: ge.addSurfaceFilling([ge.addCurveLoop([sg * eid[e] for e, sg in BLK[b]])]) for b in sorted(Q["nb"])}
        S.append(dict(x=x, pt=pt, pid=pid, eid=eid, face=face, law=Q["law"]))
    vols, xfaces, refs = {}, {}, {}
    for ii, (xa, xb, rg) in enumerate(ivs):
        A, B = S[ii], S[ii + 1]; blks = sorted(PRESENT[rg])
        ed = sorted({e for b in blks for e, _ in BLK[b]}); pn = sorted({p for e in ed for p in EDGE[e]}); xl = {}
        for p in pn:                    # x 区間は 1 セル幅: 両端の station 点を結ぶ直線 (内部節点なし。全節点が station 面上)
            xl[p] = ge.addLine(A["pid"][p], B["pid"][p]); gm.setTransfiniteCurve(xl[p], 2, "Progression", 1.0)
        for e in ed:
            a, b = EDGE[e]
            xfaces[(ii, e)] = ge.addSurfaceFilling([ge.addCurveLoop([A["eid"][e], xl[b], -B["eid"][e], -xl[a]])])
        for b in blks:
            fs = [("s", ii, b), ("s", ii + 1, b)] + [("x", ii, e) for e, _ in BLK[b]]
            for f in fs: refs[f] = refs.get(f, 0) + 1
            sl = ge.addSurfaceLoop([A["face"][b], B["face"][b]] + [xfaces[(ii, e)] for e, _ in BLK[b]]); vols[(ii, b)] = ge.addVolume([sl])
    moves = []
    for St in S:
        for e, tg in St["eid"].items(): set_curve_law(tg, St["law"][e], moves)
        for f in St["face"].values(): gm.setTransfiniteSurface(f); gm.setRecombine(2, f)
    for f in xfaces.values(): gm.setTransfiniteSurface(f); gm.setRecombine(2, f)
    for v in vols.values(): gm.setTransfiniteVolume(v)
    # ---- 境界タグ: 参照数 1 の面だけ
    groups, untag = {}, []
    for f, c in refs.items():
        assert c in (1, 2), (f, c)
        if c == 2: continue
        if f[0] == "s":
            _, si, b = f; x = xs[si]
            if si == 0: t = "inlet" if b.startswith("N_") else "ext_in"
            elif si == len(xs) - 1: t = "outlet"
            elif b == "SW" and abs(x - G.LSW) < 1e-12: t = "sidewall_end"
            elif b in ("CW1", "CW2") and abs(x - G.LCOWL) < 1e-12: t = "cowl_base"
            else: t = None
            tg = S[si]["face"][b]
        else:
            _, ii, e = f; t = xface_tag(e, ivs[ii][2]); tg = xfaces[(ii, e)]
        if t is None: untag.append(f)
        else: groups.setdefault(t, []).append(tg)
    # 表にあるのに共有面だった、の逆向きの取りこぼしも見る (参照数 2 の面にタグは付かない)
    ge.synchronize()
    for t in sorted(groups): gmsh.model.addPhysicalGroup(2, groups[t], name=t)
    gmsh.model.addPhysicalGroup(3, list(vols.values()), name="fluid")
    apply_node_laws(moves)                  # B4-4: 明示分布の辺 (非対称の対角 e2/e6) の節点を移す (1D を張ってから)
    if P.get("DIM2"):                       # 断面だけ切って 2D の歪みをブロック別・station 別に見る (分布の調整用。速い)
        gmsh.model.mesh.generate(2); ntag, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3)
        idx = np.zeros(int(ntag.max()) + 1, np.int64); idx[ntag.astype(np.int64)] = np.arange(len(ntag)); res = {}
        for si, St in enumerate(S):
            for b, f in St["face"].items():
                q = xyz[idx[np.asarray(gmsh.model.mesh.getElements(2, f)[2][0], np.int64)].reshape(-1, 4)]; sk = np.zeros(len(q))
                for c_ in range(4):
                    u, v = q[:, (c_ - 1) % 4] - q[:, c_], q[:, (c_ + 1) % 4] - q[:, c_]
                    th = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', u, v) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)), -1, 1)))
                    sk = np.maximum(sk, np.abs(th - 90) / 90)
                res.setdefault(b, []).append((round(xs[si] / H, 3), round(float(sk.max()), 3), int((sk > 0.7).sum())))
        for b in res: print("%-7s" % b, " ".join("%.3f:%.2f(%d)" % t for t in res[b] if t[1] > 0.3))
        n2 = int(len(ntag)); cnt = counts
        print("counts", cnt, "nx", [xlaw[i][0] for i in range(len(ivp))], "stations", len(xs)); print("seams", {k: [round(v[0], 4), round(v[1], 4)] for k, v in seams.items()})
        print("nodes (全節点が station 面上なので 2D の節点数 = 3D の節点数):", n2, " 生成時間 %.1f s" % (time.time() - t0))
        gmsh.finalize(); return dict(VERDICT="DIM2", nodes=n2, stations=len(xs), counts=cnt, seam_spacings=seams)
    gmsh.model.mesh.generate(3)
    t_gen = time.time() - t0
    rep = dict(params={k: P[k] for k in sorted(P)}, contours={k: contours[k + "_src"] for k in ("ramp", "cowl")},
               geometry=dict(L_sw=G.LSW / H, L_cowl=G.LCOWL / H, x_end=G.XEND / H, t_cowl=G.TC / H, t_sw=G.TSW / H, W_half=G.ZW / H,
                             y_c_te=float(G.yc(G.LCOWL)) / H, y_o_te=G.yo_te / H, cowl_base_height=(float(G.yc(G.LCOWL)) - G.yo_te) / H,
                             cowl_base_height_expected=G.TC / math.sqrt(1 + G.s_te ** 2) / H * (1 + G.s_te ** 2), y_r_inlet=float(G.yr(0.0)) / H),
               phys_stations=[x / H for x in xp], phys_regions=[v[2] for v in ivp], nx=[abs(xlaw[i][0]) for i in range(len(ivp))],
               x_ratio=[abs(xlaw[i][1]) for i in range(len(ivp))], wake_caps=[capf(i) / H if capf(i) else None for i in range(len(ivp))],
               n_stations=len(xs), x_adjacent_ratio_design=x_ratio_design, dx_min=float(dxs.min()), dx_max=float(dxs.max()),
               counts=counts, seam_spacings=seams, seam_spacings_ref_h1=h_ring / H, **wall_info,
               far_spacing_max=dict(y=far_y / H, z=far_z / H, limit=P["HFAR"]),
               layer_rules=dict(wall="NL_w = prog_n(DR, h1_w, 1.2) + 4 (壁別の h1_w で。端面 sidewall_end・cowl_base も同式。cowl_side は sidewall_out と同じ h1。plan §6.2・B4b-5)",
                                ring="NR: 対角 e2/e6 の Bump の増加側が 3 壁 (ramp・sidewall_in・cowl_in) の最大の層数 NL 以上・公比 <= 1.15 (B4-2)",
                                z="NZ = max(P.NZ, prog_n(W/2, min h1 (ramp・sidewall_in・cowl_in), 1.15) + 1) (B4-2)",
                                h1=wall_info["wall_h1"], dx_wake_max=dxw / H, wake_length=G.LWAKE / H, growth=g),
               n_blocks=len(vols), blocks_per_region={r_: len(PRESENT[r_]) for r_ in "ABC"}, untagged_boundary_faces=[str(u) for u in untag],
               tags={t: len(v) for t, v in groups.items()}, time_generate_s=t_gen)
    gmsh.write(out + ".msh")                        # 検査の前に書き出す (check が gmsh を finalize して格子を解放する)
    if coupon: rep["coupon"] = dict(x_requested=float(P["COUPON_X"]), stations=[x / H for x in xs], regions=sorted({v[2] for v in ivs}))
    rep.update(check(G, P, vols, groups, hw, None, out, NLw, None, [x for x, _, _ in ivs], coupon=coupon))     # B4b-5: 壁別の目標 h1・期待層数
    rep["time_total_s"] = time.time() - t0
    gates = dict(untagged=not untag, other_elems=rep["other_elems"] == {}, neg_jac_f32=rep["neg_jac_f32"] == 0, mixed_sign_volumes=rep["mixed_sign_volumes"] == 0,
                 skew=rep["skew_max"] <= 0.90, ar=rep["ar_gt5000"] == 0 and rep["ar1000_skew_max"] <= 0.30, topology=rep["topology"]["ok"],
                 first_layer=all(v["ok"] for v in rep["first_layer"].values()), adjacent_spacing=rep["adjacent_spacing"]["ok"],
                 wake_dx=rep["wake_dx"]["ok"], shape=rep["shape"]["ok"])
    if coupon:                              # 試験片には端面・後流区間・カウル内壁が無い (対象外。本体の生成で判定する)
        for k_ in ("wake_dx", "shape"): gates.pop(k_)
    rep["gates"] = gates
    rep["VERDICT"] = "PASS" if all(gates.values()) else "FAIL"
    json.dump(rep, open(out + "_report.json", "w"), indent=1, ensure_ascii=False)
    if not quiet: print(json.dumps({k: v for k, v in rep.items() if k not in ("params", "phys_stations")}, indent=1, ensure_ascii=False))
    return rep


def section_targets(xq, ivp, H):
    """全長診断 (SECX) の評価位置: 要求 x/H をそのまま評価位置 x = x/H × H にし (station への丸めなし)、その位置を含む物理区間の領域 (A/B/C) の集合を返す。
    物理 station 上の x は両側の区間の領域 (station に張る面と同じ規則)。物理区間の範囲外は ValueError。戻り値は (要求 x/H, 評価 x, 領域の集合) の列"""
    out = []; x0, x1 = ivp[0][0], ivp[-1][1]
    for xv in [float(v) for v in xq]:
        x = xv * H
        if not (x0 - 1e-12 * H <= x <= x1 + 1e-12 * H): raise ValueError("SECX の x/H = %.6g が診断の範囲 [%.6g, %.6g] の外" % (xv, x0 / H, x1 / H))
        out.append((xv, x, {r_ for a_, b_, r_ in ivp if a_ - 1e-12 * H <= x <= b_ + 1e-12 * H}))
    return out


def section_2d(G, items, seams, quiet=False):
    """断面 1 枚の 2D 四角形格子で、辺を共有する 2 つの四角形の「共有辺の端から出る辺」の長さの比 (3D の隣接間隔比の断面方向と同じ定義) と
    skew をブロックの組ごとに出す (分布の調整用の診断。受入の判定は 3D の adjacent_spacing_check)。
    items は (要求 x/H, 評価 x [m], 断面 Q) の列。SEC2D は最近傍 station の断面、SECX (全長診断) は要求位置そのものの断面を渡す。
    各組の最大比には位置と、比を作った 2 本の辺 (共有辺の端の節点から両側のブロックへ出る辺) のベクトル (/H) を付ける"""
    out = {}
    for xv, x, Q in items:
        gmsh.model.add("sec"); ge = gmsh.model.geo; pt, ed = Q["pt"], Q["ed"]
        pid = {p: ge.addPoint(x, pt[p][0], pt[p][1]) for p in sorted({p for e in ed for p in EDGE[e]})}; eid = {}
        for e in ed:
            a_, b_ = EDGE[e]
            eid[e] = ge.addSpline([pid[a_]] + [ge.addPoint(x, p_[0], p_[1]) for p_ in Q["wr"][e][1:-1]] + [pid[b_]]) if e in Q["wr"] else ge.addLine(pid[a_], pid[b_])
        face = {b: ge.addSurfaceFilling([ge.addCurveLoop([sg * eid[e] for e, sg in BLK[b]])]) for b in sorted(Q["nb"])}
        moves = []
        for e, tg in eid.items(): set_curve_law(tg, Q["law"][e], moves)
        for f in face.values(): ge.mesh.setTransfiniteSurface(f); ge.mesh.setRecombine(2, f)
        ge.synchronize(); apply_node_laws(moves); gmsh.model.mesh.generate(2)
        ntag, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3)
        idx = np.zeros(int(ntag.max()) + 1, np.int64); idx[ntag.astype(np.int64)] = np.arange(len(ntag))
        qs, lab = [], []; names = sorted(face)
        for b in names:
            q_ = idx[np.asarray(gmsh.model.mesh.getElements(2, face[b])[2][0], np.int64)].reshape(-1, 4); qs.append(q_); lab.append(np.full(len(q_), names.index(b)))
        q = np.vstack(qs); lab = np.concatenate(lab)
        # 辺 (k, k+1) を共有する四角形の組。端点 a から出る辺は (a, 前の頂点)、端点 b から出る辺は (b, 次の頂点)
        rec = {}
        for k in range(4):
            a, b = q[:, k], q[:, (k + 1) % 4]; pa, nb_ = q[:, (k + 3) % 4], q[:, (k + 2) % 4]
            for i in range(len(q)):
                rec.setdefault((min(a[i], b[i]), max(a[i], b[i])), []).append((i, {a[i]: pa[i], b[i]: nb_[i]}))
        pair, skw = {}, {}
        for key, v in rec.items():
            if len(v) != 2: continue
            (i1, d1), (i2, d2) = v
            for n_ in key:
                v1, v2 = xyz[d1[n_]] - xyz[n_], xyz[d2[n_]] - xyz[n_]; l1, l2 = np.linalg.norm(v1), np.linalg.norm(v2)
                r = max(l1 / l2, l2 / l1); b1, b2 = names[lab[i1]], names[lab[i2]]
                if b1 > b2: b1, b2, v1, v2, l1, l2 = b2, b1, v2, v1, l2, l1
                k_ = "%s|%s" % (b1, b2)
                if r > pair.get(k_, (0,))[0]:
                    pair[k_] = (float(r), [float(c) / G.H for c in xyz[n_][1:]],
                                {b1: [float(c) / G.H for c in v1[1:]], b2 + ("'" if b1 == b2 else ""): [float(c) / G.H for c in v2[1:]]})
        for b in names:
            qq = xyz[q[lab == names.index(b)]]; sk = np.zeros(len(qq))
            for c_ in range(4):
                u, v = qq[:, (c_ - 1) % 4] - qq[:, c_], qq[:, (c_ + 1) % 4] - qq[:, c_]
                th = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', u, v) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)), -1, 1)))
                sk = np.maximum(sk, np.abs(th - 90) / 90)
            skw[b] = float(sk.max())
        gmsh.model.remove()
        key_ = "%.4f" % (x / G.H)
        out[key_] = dict(x_requested=xv, x_evaluated=x / G.H, pairs=dict(sorted(pair.items(), key=lambda kv: -kv[1][0])), skew=skw, nodes=int(len(xyz)))
        if not quiet:
            print("x/H 要求 %.6g 評価 %s  nodes %d  max ratio %.4f  max skew %.3f (%s)" % (xv, key_, len(xyz), max(v[0] for v in pair.values()), max(skw.values()), max(skw, key=skw.get)))
            for k_, v in list(out[key_]["pairs"].items())[:10]: print("   %-14s %.4f at (y,z)/H = (%.5f, %.5f)" % (k_, v[0], v[1][0], v[1][1]))
    gmsh.finalize(); return dict(VERDICT="DIM2", sections=out, seam_spacings=seams)


# ---------------------------------------------------------------- 検査
NBR = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7), 4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}
HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
HE = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]


def cornerJ(Pn):
    return np.stack([np.einsum('ij,ij->i', np.cross(Pn[:, a] - Pn[:, c], Pn[:, b] - Pn[:, c]), Pn[:, d] - Pn[:, c]) for c, (a, b, d) in NBR.items()], 1)


def hex_quality(xyz, hx, chunk=400000):
    """頂点 Jacobian (float64 / float32 座標)・skew (面の内角の 90° からのずれ / 90)・AR (最長辺 / 最短辺) をヘキサごとに返す (gmsh 非依存の純関数)。
    戻り値: (J, J32, sk, AR, Led)。符号の基準は呼び出し側で J の中央値の符号に揃える"""
    x32 = xyz.astype(np.float32).astype(float)
    J = np.vstack([cornerJ(xyz[hx[i:i + chunk]]) for i in range(0, len(hx), chunk)])
    J32 = np.vstack([cornerJ(x32[hx[i:i + chunk]]) for i in range(0, len(hx), chunk)]); del x32
    sk = np.zeros(len(hx))
    for f in HF:
        Pf = [xyz[hx[:, q]] for q in f]
        for c in range(4):
            a, b = Pf[(c - 1) % 4] - Pf[c], Pf[(c + 1) % 4] - Pf[c]
            th = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)), -1, 1)))
            sk = np.maximum(sk, np.abs(th - 90) / 90)
    Led = np.stack([np.linalg.norm(xyz[hx[:, a]] - xyz[hx[:, b]], axis=1) for a, b in HE], 1); AR = Led.max(1) / Led.min(1)
    return J, J32, sk, AR, Led


# 端面 (第一内部点距離 h1e・期待層数 NL_END で判定。plan §6.2 の表: cowl_base・sidewall_end)。
# cowl_side は一般壁 (h1・NL。2026-10-06 仕様訂正: z = 側壁外面の平面上で sidewall_out と節点・z 分布を共有するので h1 と h1e を同時に要求できない)
END_TAGS = ("sidewall_end", "cowl_base")


def model_classifier(G):
    """接続模型の判定区分 (フィレット円弧上 / 半径が変わる区間) を返す関数を作る。first_layer_check の classify に渡す"""
    def classify(t, pts):
        arc = np.zeros(len(pts), bool); edge = np.zeros(len(pts), bool)
        if t in ("ramp", "sidewall_in", "cowl_in"):
            dl = 0.05 * G.H      # r -> 0 の端点 (フィレットの生え際) も円弧の区分に入れる
            rr = np.array([np.max([G.r(x_ + q) for q in (-dl, 0.0, dl)], axis=0) for x_ in pts[:, 0]]).reshape(-1, 2)
            yrn = np.array([G.yr(x_) for x_ in pts[:, 0]]); mg = 0.5 * G.DR
            arc |= (rr[:, 0] > 0) & (pts[:, 2] > G.ZW - rr[:, 0] - mg) & (pts[:, 1] > yrn - rr[:, 0] - mg)
            arc |= (rr[:, 1] > 0) & (pts[:, 2] > G.ZW - rr[:, 1] - mg) & (pts[:, 1] < G.yc(pts[:, 0]) + rr[:, 1] + mg)
            for (xa, ua, la), (xb, ub, lb) in zip(G.prof[:-1], G.prof[1:]):       # 半径が変わる区間 (フィレット遷移) はダクト内の壁全周を ±10 % の区分に入れる (plan §4.10 (e))
                if (ua, la) != (ub, lb): edge |= (pts[:, 0] >= xa - 1e-12) & (pts[:, 0] <= xb + 1e-12)
        return arc, edge
    return classify


def first_layer_check(xyz, hx, bq, h1, h1e, nl, nl_end, dr, classify=None, wall_tags=WALL_TAGS, end_tags=END_TAGS,
                      expected_tags=None, gmax=1.2, unit=1.0):
    """壁第一層の検査 (gmsh に依存しない純関数。plan tooling-sern-mesh-blocking §5.1 B4-0)。
    入力: 節点座標 xyz (N,3)、ヘキサ hx (M,8, gmsh 順)、壁タグ -> 境界四角形 (K,4) の dict bq、第一層の指定 h1 (端面は h1e)、
    期待層数 nl (端面は nl_end)。B4b-5: h1 と nl に壁タグ -> 値の dict を渡すと壁ごとの目標・期待層数で判定する (h1e・nl_end は使わない)、稜の近傍幅の基準 dr、区分関数 classify(t, pts) -> (arc, edge) (None なら円弧区分なし)。
    壁タグ別に次を見て、1 つでも外れたら ok=False:
      (1) 第一内部点の資格: 壁節点 w に接する全ヘキサの節点のうち、**どの壁タグにも属さない**節点で、壁法線の内側 (n.d > 0) にあるものが候補。
          (n.d が |d| の 1e-6 以下の節点 = 壁と同一平面で丸めだけ内側に見える節点は候補にしない。2026-10-06 B4-1 で端面 station 面上の
          隣接節点 (n.d ~ 1e-18 m) が選ばれ比 3e-14 で偽 FAIL したため)。
          候補のうち w と辺でつながるものがあればその中から、法線からの横ずれ |d - (n.d) n| が最小のものを第一内部点 p1 とする
          (隅の近くで格子線が扇形に傾いても、同じヘキサの対角の点でなく格子線上の点を取る)。辺でつながる候補が無いときは
          接する全ヘキサの候補から選ぶ (直交格子の凹の稜では p1 は対角の節点になる)。これを許すのは稜 (複数の壁タグに属する節点) だけで、
          非稜の節点で辺でつながる候補が無い・候補が 1 つも無い壁節点は「評価不能」。
      (2) 全数被覆: タグの全節点を評価する。法線が作れない節点 (ヘキサの面でない四角形だけに属する) と (1) の評価不能を数え、1 以上なら不合格。
          expected_tags に挙げたタグが bq に無ければそのタグは不合格 (missing)。
      (3) 実層数: w -> p1 の向きから、辺でつながる節点を「直前の向きとの cos が最大 (> 0.5)・壁節点でない・来た節点でない」で 1 つずつたどり、
          間隔 d_k = |p_k - p_{k-1}| が減らない (d_{k+1}/d_k >= 0.98) 間を層として数える (上限 2 × 期待層数)。期待層数 (nl / 端面は nl_end) より少なければ不合格。
          稜で p1 が辺でつながらない (対角) 節点は、たどる向きが格子線に乗らないので層数・間隔比を数えない (layers_skipped_diag に件数を出す。
          非稜の節点で p1 が辺でつながらなければ評価不能)。
      (4) 隣接間隔比: (3) でたどった連続 3 節点の間隔比 d_{k+1}/d_k を壁から期待層数ぶん実測し、> gmax (許容の上乗せなし) の節点があれば不合格。
    第一層の比 |n.(p1 - w)| / h1 の区分 (平面部 ±5 %、稜から 2 dr 以内と半径が変わる区間 ±10 %、フィレット円弧上 −10〜+45 %) は従来どおり。
    """
    from scipy.spatial import cKDTree
    N = len(xyz); wall_tags = tuple(wall_tags)
    wallnode = np.zeros(N, bool)
    for t in wall_tags:
        if t in bq: wallnode[np.asarray(bq[t]).ravel()] = True
    multi = np.zeros(N, int)
    for t in wall_tags:
        if t in bq: multi[np.unique(bq[t])] += 1
    ridge = xyz[multi >= 2]; tree = cKDTree(ridge) if len(ridge) else None
    # 壁四角形 -> それを面に持つヘキサ (法線の向きを内側に揃えるのに使う)
    cand = np.where(wallnode[hx].sum(1) >= 4)[0]; key = {}
    for hi in cand:
        for f in HF:
            if wallnode[hx[hi, list(f)]].all(): key[tuple(sorted(hx[hi, list(f)]))] = (hi, f)
    # 辺でつながる節点 (全ヘキサの 12 辺)。行ごとに -1 で詰めた表
    E = np.vstack([hx[:, [a, b]] for a, b in HE]).astype(np.int64); E = np.vstack([E, E[:, ::-1]])
    E = np.unique(E[:, 0] * N + E[:, 1]); src, dst = E // N, E % N; del E
    deg = np.bincount(src, minlength=N); D = int(deg.max()) if N else 0; st = np.r_[0, np.cumsum(deg)[:-1]]
    NB = -np.ones((N, max(D, 1)), np.int64); NB[src, np.arange(len(src)) - st[src]] = dst; del src, dst
    # 節点 -> 接するヘキサ (壁節点だけで使う)
    hflat = hx.ravel().astype(np.int64); od = np.argsort(hflat, kind="stable"); hcnt = np.bincount(hflat, minlength=N); hst = np.r_[0, np.cumsum(hcnt)[:-1]]
    def hexes_of(nd):
        Dh = int(hcnt[nd].max()) if len(nd) else 1; Hm = -np.ones((len(nd), Dh), np.int64)
        for j in range(Dh):
            m = hcnt[nd] > j; Hm[m, j] = od[hst[nd[m]] + j] // 8
        return Hm
    acc = {}
    for t in wall_tags:                      # 1 巡目: 壁タグ別の節点法線 (隣接する壁四角形の面積ベクトルの和)
        if t not in bq: continue
        q = np.asarray(bq[t]); nrm = np.zeros((N, 3)); uneval = 0
        for row in q:
            hit = key.get(tuple(sorted(row)))
            if hit is None: uneval += 1; continue
            hi, f = hit; loc = list(f); Pq = xyz[hx[hi, loc]]; n_ = np.cross(Pq[2] - Pq[0], Pq[3] - Pq[1])
            cen_in = xyz[hx[hi]].mean(0) - Pq.mean(0); n_ = n_ * np.sign(n_ @ cen_in)
            nrm[hx[hi, loc]] += n_
        acc[t] = (nrm, uneval)
    fl = {}
    for t in (expected_tags or ()):
        if t not in bq: fl[t] = dict(missing=True, ok=False)
    for t in wall_tags:
        if t not in bq: continue
        q = np.asarray(bq[t]); nrm, uneval = acc[t]
        nd_all = np.unique(q); has_n = np.linalg.norm(nrm[nd_all], axis=1) > 0; no_normal = int((~has_n).sum()); nd = nd_all[has_n]
        n_ = nrm[nd] / np.linalg.norm(nrm[nd], axis=1)[:, None]
        ns = nrm[nd].copy()                                                   # タグの境が滑らか (フィレットの 45° 点など) なら隣のタグの面も法線に入れる。稜 (50° 超) は入れない
        for t2 in acc:
            if t2 == t: continue
            v2 = acc[t2][0][nd]; l2 = np.linalg.norm(v2, axis=1); okk = l2 > 0
            cs = np.zeros(len(nd)); cs[okk] = np.einsum('ij,ij->i', n_[okk], v2[okk]) / l2[okk]
            ns += np.where((cs > math.cos(math.radians(50.0)))[:, None], v2, 0.0)
        n_ = ns / np.linalg.norm(ns, axis=1)[:, None]
        # (1) 第一内部点: 接する全ヘキサの節点から、壁でない・内側・横ずれ最小
        p1 = -np.ones(len(nd), np.int64)
        for c0 in range(0, len(nd), 20000):                 # メモリを抑えるため節点を分けて処理
            sl = slice(c0, c0 + 20000); ndc, nc = nd[sl], n_[sl]
            Hm = hexes_of(ndc); C = np.where(Hm[:, :, None] >= 0, hx[np.maximum(Hm, 0)], -1).reshape(len(ndc), -1)
            dv = xyz[np.maximum(C, 0)] - xyz[ndc][:, None, :]; dn = np.einsum('ijk,ik->ij', dv, nc)
            tang = np.linalg.norm(dv - dn[:, :, None] * nc[:, None, :], axis=2)
            # 内側 = n.d > 0。ただし壁と同一平面の節点 (端面の station 面上の隣など) が座標の丸め (1 ulp) で n.d > 0 になるのは除く
            valid = (C >= 0) & ~wallnode[np.maximum(C, 0)] & (dn > 1e-6 * np.linalg.norm(dv, axis=2))
            ise = (C[:, :, None] == NB[ndc][:, None, :]).any(2)          # 辺でつながる候補を優先 (扇形に傾いた格子線でも線上の点を取る)
            valid &= np.where(ise.__and__(valid).any(1)[:, None], ise, True)
            tang = np.where(valid, tang, np.inf); j = np.argmin(tang, axis=1); r_ = np.arange(len(ndc))
            p1[sl] = np.where(np.isfinite(tang[r_, j]), C[r_, j], -1)
        qual = p1 >= 0
        edgecon = qual & np.any(NB[nd] == p1[:, None], axis=1)
        isridge = multi[nd] >= 2; diag = qual & ~edgecon & isridge
        noqual = ~qual | (qual & ~edgecon & ~isridge)               # 候補なし / 非稜で辺につながらない = 評価不能
        if isinstance(h1, dict): target, nexp = h1[t], int(nl[t])          # 壁別 (B4b-5)
        else: is_end = t in end_tags; target = h1e if is_end else h1; nexp = int(nl_end if is_end else nl)
        ev = ~noqual; nde, ne, pe = nd[ev], n_[ev], p1[ev]
        ratio = np.abs(np.einsum('ij,ij->i', ne, xyz[pe] - xyz[nde])) / target
        # (3)(4) 実層数と隣接間隔比: 辺でつながる p1 から格子線をたどる
        walk = ~diag[ev]; lay = np.zeros(len(nde), int); lay[walk] = 1; rmax = np.full(len(nde), np.nan)
        wi = np.flatnonzero(walk); prev, cur = nde[wi], pe[wi]; dprev = np.linalg.norm(xyz[cur] - xyz[prev], axis=1)
        alive = np.ones(len(wi), bool); rcur = np.full(len(wi), np.nan)
        for k in range(1, 2 * max(nexp, 1)):
            if not alive.any(): break
            a = np.flatnonzero(alive); Cn = NB[cur[a]]; u = xyz[cur[a]] - xyz[prev[a]]; u /= np.linalg.norm(u, axis=1)[:, None]
            vv = xyz[np.maximum(Cn, 0)] - xyz[cur[a]][:, None, :]; ln = np.linalg.norm(vv, axis=2)
            cs = np.einsum('ijk,ik->ij', vv, u) / np.where(ln > 0, ln, 1.0)
            ok_ = (Cn >= 0) & (Cn != prev[a][:, None]) & ~wallnode[np.maximum(Cn, 0)] & (cs > 0.5)
            cs = np.where(ok_, cs, -np.inf); jj = np.argmax(cs, axis=1); go = np.isfinite(cs[np.arange(len(a)), jj])
            nx_ = Cn[np.arange(len(a)), jj]; dnext = np.where(go, ln[np.arange(len(a)), jj], 0.0); rr = dnext / dprev[a]
            go &= rr >= 0.98
            if k < nexp: rcur[a[go]] = np.fmax(rcur[a[go]], rr[go])     # 期待層数の範囲の間隔比だけを判定に使う
            alive[a[~go]] = False; g_ = a[go]; lay[wi[g_]] += 1
            prev[g_], cur[g_], dprev[g_] = cur[g_], nx_[go], dnext[go]
        rmax[wi] = rcur
        # 第一層の比の区分: (i) フィレット円弧の上は [-10 %, +45 %] (バタフライ位相では 45° 点で最大 sqrt2 倍厚い。薄くはしない)、
        #       (ii) 稜 (複数の壁タグに属する節点) から 2 dr 以内は ±10 %、(iii) それ以外の平面部は ±5 %
        arc = np.zeros(len(nde), bool); edge = np.zeros(len(nde), bool)
        if tree is not None and len(nde): edge |= tree.query(xyz[nde])[0] < 2 * dr
        if classify is not None and len(nde):
            a_, e_ = classify(t, xyz[nde]); arc |= a_; edge |= e_
        edge &= ~arc; core = ~(edge | arc); dev = ratio - 1
        bad = (core & (np.abs(dev) > 0.05)) | (edge & (np.abs(dev) > 0.10)) | (arc & ((dev < -0.10) | (dev > 0.45)))
        short = walk & (lay < nexp); steep = walk & (rmax > gmax)
        n_uneval = no_normal + int(noqual.sum())
        mm = lambda v: (float(v.min()), float(v.max())) if len(v) else (float("nan"), float("nan"))
        fl[t] = dict(target_h1=float(target) / unit, nodes=int(len(nd_all)), evaluated=int(len(nde)), unevaluable_nodes=n_uneval, unevaluable_quads=int(uneval),
                     ratio_min=mm(ratio)[0], ratio_max=mm(ratio)[1],
                     layers_expected=nexp, layers_min=int(lay[walk].min()) if walk.any() else 0, layers_short=int(short.sum()),
                     layers_skipped_diag=int(diag.sum()), adjacent_ratio_max=float(np.nanmax(rmax)) if np.isfinite(rmax).any() else float("nan"),
                     adjacent_ratio_gt=int(steep.sum()), adjacent_ratio_limit=gmax)
        for nm, mk in (("core", core), ("edge", edge), ("arc", arc)):
            fl[t][nm] = dict(nodes=int(mk.sum()), out=int((bad & mk).sum()), min=float(ratio[mk].min()) if mk.any() else 1.0, max=float(ratio[mk].max()) if mk.any() else 1.0)
        fl[t]["ok"] = bool(uneval == 0 and n_uneval == 0 and len(nde) > 0 and not bad.any() and not short.any() and not steep.any()); dev = np.abs(dev)
        fl[t]["_bad"] = [[float(v) for v in xyz[nde[i]] / unit] + [float(ratio[i]), "core" if core[i] else ("edge" if edge[i] else "arc")] for i in np.flatnonzero(bad)[:8]]
        fl[t]["_short"] = [[float(v) for v in xyz[nde[i]] / unit] + [int(lay[i])] for i in np.flatnonzero(short)[:8]]
        fl[t]["_steep"] = [[float(v) for v in xyz[nde[i]] / unit] + [float(rmax[i])] for i in np.flatnonzero(steep)[:8]]
        fl[t]["_uneval"] = [[float(v) for v in xyz[i] / unit] for i in np.r_[nd_all[~has_n], nd[noqual]][:8]]
        fl[t]["_worst"] = [[float(v) for v in xyz[nde[i]]] + [float(ratio[i])] for i in np.argsort(-dev)[:5]]
    return fl


# 各面 HF[i] の節点に対応する、向かいの面の節点 (面に直交する辺の他端)
HF_OPP = [(4, 7, 6, 5), (0, 1, 2, 3), (3, 2, 6, 7), (0, 3, 7, 4), (1, 0, 4, 5), (2, 1, 5, 6)]


def adjacent_spacing_check(xyz, hx, gmax=1.2, labels=None, label_names=None, unit=1.0, xdir=0.5, nworst=12, chunk=2000000):
    """出力実座標の隣接間隔比 (plan §6.2「実測の隣接間隔比 <= 1.2、全方向・継ぎ目を含む」。§5.1 B4-1)。
    面を共有する 2 つのヘキサについて、共有面の各節点から**面に直交する辺** (向かいの面への辺) の長さを両側で取り、
    比 max(L1/L2, L2/L1) を節点ごとに出す。ブロック内の隣接も、ブロック・区間の継ぎ目を跨ぐ隣接も同じ式で見る (gmsh 非依存の純関数)。
    分類: 辺の向きの |dx|/|d| > xdir を x 方向、それ以外を断面内とする。比 > gmax (許容の上乗せなし) が 1 つでもあれば ok=False。
    labels (ヘキサごとの整数、例えばブロック名の番号) を渡すと、超過をブロックの組ごとに数える"""
    xyz = np.asarray(xyz, float); hx = np.asarray(hx, np.int64); M = len(hx)
    if M == 0: return dict(ok=False, reason="ヘキサなし")
    def keys(Q):
        Q = np.sort(Q, axis=1); return (Q[:, 0] << 32) | Q[:, 1], (Q[:, 2] << 32) | Q[:, 3]
    k1 = np.empty(6 * M, np.int64); k2 = np.empty(6 * M, np.int64)
    for fi, f in enumerate(HF):
        k1[fi * M:(fi + 1) * M], k2[fi * M:(fi + 1) * M] = keys(hx[:, list(f)])
    od = np.lexsort((k2, k1)); k1s, k2s = k1[od], k2[od]; del k1, k2
    same = (k1s[1:] == k1s[:-1]) & (k2s[1:] == k2s[:-1]); del k1s, k2s
    ia = od[:-1][same]; ib = od[1:][same]; del od, same
    n_pairs = len(ia)
    def side(e):                    # 面 entry -> (面の節点 (節点番号順), 直交辺の長さ, 直交辺のベクトル)
        h = e % M; fi = e // M; F = np.array(HF)[fi]; O = np.array(HF_OPP)[fi]
        nf = np.take_along_axis(hx[h], F, 1); no = np.take_along_axis(hx[h], O, 1)
        o = np.argsort(nf, axis=1); nf = np.take_along_axis(nf, o, 1); no = np.take_along_axis(no, o, 1)
        v = xyz[no] - xyz[nf]; return nf, np.linalg.norm(v, axis=2), v, h
    rmax_all = np.zeros(n_pairs); isx = np.zeros(n_pairs, bool); hpair = np.zeros((n_pairs, 2), np.int64); worst_node = np.zeros(n_pairs, np.int64)
    for c0 in range(0, n_pairs, chunk):
        sl = slice(c0, c0 + chunk)
        nfa, La, va, ha = side(ia[sl]); nfb, Lb, vb, hb = side(ib[sl])
        assert np.array_equal(nfa, nfb)
        r = np.maximum(La / Lb, Lb / La); j = np.argmax(r, axis=1); rr = np.arange(len(j))
        rmax_all[sl] = r[rr, j]; worst_node[sl] = nfa[rr, j]
        vm = va.mean(1); isx[sl] = np.abs(vm[:, 0]) > xdir * np.linalg.norm(vm, axis=1)
        hpair[sl, 0], hpair[sl, 1] = ha, hb
    bad = rmax_all > gmax
    out = dict(gmax=gmax, pairs=int(n_pairs), ratio_max=float(rmax_all.max()) if n_pairs else float("nan"), n_gt=int(bad.sum()),
               x=dict(pairs=int(isx.sum()), ratio_max=float(rmax_all[isx].max()) if isx.any() else 1.0, n_gt=int((bad & isx).sum())),
               section=dict(pairs=int((~isx).sum()), ratio_max=float(rmax_all[~isx].max()) if (~isx).any() else 1.0, n_gt=int((bad & ~isx).sum())))
    nm = (lambda i: label_names[i]) if label_names is not None else (lambda i: int(i))
    item = lambda i: dict(ratio=float(rmax_all[i]), dir="x" if isx[i] else "section", at=[float(v) / unit for v in xyz[worst_node[i]]],
                          blocks=[nm(labels[hpair[i, 0]]), nm(labels[hpair[i, 1]])] if labels is not None else None)
    out["worst"] = [item(i) for i in np.argsort(-rmax_all)[:nworst]]
    for cls, m in (("x", isx), ("section", ~isx)):
        ii_ = np.flatnonzero(m); out[cls]["worst"] = [item(i) for i in ii_[np.argsort(-rmax_all[ii_])[:nworst]]]
    if labels is not None and n_pairs:              # 参考: 全ブロック組の最大比 (超過の有無によらず。判定には使わない)
        la0, lb0 = np.sort(np.stack([labels[hpair[:, 0]], labels[hpair[:, 1]]], 1), axis=1).T
        kk = (la0.astype(np.int64) * 65536 + lb0) * 2 + isx; uk, inv = np.unique(kk, return_inverse=True)
        mx = np.zeros(len(uk)); np.maximum.at(mx, inv, rmax_all)
        out["max_by_block_pair"] = {"%s|%s|%s" % (nm(int(k // 2 // 65536)), nm(int(k // 2 % 65536)), "x" if k % 2 else "section"): float(v)
                                    for k, v in zip(uk, mx)}
        # 各組の最大比の位置と、比を作った 2 本の辺 (共有面の節点から面に直交して両側のヘキサへ出る辺) のベクトル (plan §6.5「失敗節点の接続ベクトル」)
        arg = np.full(len(uk), -1, np.int64); hit = rmax_all >= mx[inv]
        arg[inv[hit][::-1]] = np.flatnonzero(hit)[::-1]                # 組ごとに最大比を取る最初の面対
        nfa, La, va, ha = side(ia[arg]); _, Lb, vb, hb = side(ib[arg]); jj = np.argmax(np.maximum(La / Lb, Lb / La), axis=1)
        out["at_by_block_pair"] = {"%s|%s|%s" % (nm(int(k // 2 // 65536)), nm(int(k // 2 % 65536)), "x" if k % 2 else "section"):
                                   dict(ratio=float(mx[q]), at=[float(v) / unit for v in xyz[nfa[q, jj[q]]]],
                                        vec={nm(labels[ha[q]]): [float(v) / unit for v in va[q, jj[q]]], nm(labels[hb[q]]) + ("'" if labels[ha[q]] == labels[hb[q]] else ""): [float(v) / unit for v in vb[q, jj[q]]]})
                                   for q, k in enumerate(uk)}
    if labels is not None and bad.any():
        la, lb = np.sort(np.stack([labels[hpair[bad, 0]], labels[hpair[bad, 1]]], 1), axis=1).T; bi = np.flatnonzero(bad); rb = rmax_all[bad]
        grp = {}
        for a_, b_, r_, i in zip(la, lb, rb, bi):
            k = "%s|%s|%s" % (nm(a_), nm(b_), "x" if isx[i] else "section"); q = grp.setdefault(k, [0, 0.0])
            q[0] += 1; q[1] = max(q[1], float(r_))
        out["gt_by_block_pair"] = {k: dict(n=v[0], ratio_max=v[1]) for k, v in sorted(grp.items(), key=lambda kv: -kv[1][1])[:30]}
    out["ok"] = bool(n_pairs > 0 and not bad.any())
    return out


def wake_dx_check(xyz, hx, ranges, cap, unit=1.0, xdir=0.5, tol=1e-9):
    """後流区間の実測 Δx (plan §6.2: x ∈ [L, L + 0.02 H] で Δx <= 1.0e-3 H)。ヘキサの 12 辺のうち x 方向 (|dx|/|d| > xdir) で
    両端の x が区間 [a, b] に入る辺の |dx| の最大を区間ごとに出す。丸めの許容は相対 tol (1e-9) だけ。辺が 1 本も無い区間は評価不能 = 不合格"""
    xyz = np.asarray(xyz, float); hx = np.asarray(hx, np.int64); res = {}; ok = True
    span = max(float(np.ptp(xyz[:, 0])), 1e-300)
    for name, (a, b) in ranges.items():
        e_ = tol * span; mx, n = 0.0, 0
        for p_, q_ in HE:
            A, B = xyz[hx[:, p_]], xyz[hx[:, q_]]; d = B - A; dx = np.abs(d[:, 0])
            m = (dx > xdir * np.linalg.norm(d, axis=1)) & (np.minimum(A[:, 0], B[:, 0]) >= a - e_) & (np.maximum(A[:, 0], B[:, 0]) <= b + e_)
            if m.any(): mx = max(mx, float(dx[m].max())); n += int(m.sum())
        good = n > 0 and mx <= cap * (1 + tol); ok &= good
        res[name] = dict(range=[a / unit, b / unit], edges=n, dx_max=mx / unit, cap=cap / unit, ok=bool(good))
    res["ok"] = bool(ok); return res


def shape_check(xyz, bq, G, tol_c=1e-6, tol_x=1e-9):
    """形状ゲート (plan §6.2 の形状の表。品質・閉性とは独立)。単位 H で判定する。
      cowl_in 節点と参照輪郭 cowl_xy の距離 <= tol_c、ramp 節点と参照輪郭 (ramp_fillet 適用後) の距離 <= tol_c、
      cowl_out 節点から内壁輪郭 (後縁の先は接線延長) までの法線距離 = t ± tol_c、
      sidewall_end / cowl_base の全節点 |x - L| <= tol_x、cowl_base の y 範囲 = [外壁交点 y_o(L), 内壁後縁 y_c(L)] (tol_x)。
    G は Geom (ramp / cowl の Contour、TC、LSW、LCOWL、yo_te を使う)。必要なタグが無ければ不合格"""
    H = G.H; out = {}; ok = True
    def need(t):
        if t not in bq: out[t] = dict(missing=True, ok=False); return None
        return xyz[np.unique(np.asarray(bq[t]))]
    for t, C in (("cowl_in", G.cowl), ("ramp", G.ramp)):
        P_ = need(t)
        if P_ is None: ok = False; continue
        _, d = C.project(P_[:, :2]); dmax = float(np.abs(d).max()) / H; i = int(np.argmax(np.abs(d)))
        g_ = dmax <= tol_c; ok &= g_
        out[t] = dict(nodes=int(len(P_)), dist_max=dmax, at=[float(v) / H for v in P_[i]], tol=tol_c, ok=bool(g_))
    P_ = need("cowl_out")
    if P_ is None: ok = False
    else:
        _, d = G.cowl.project(P_[:, :2]); e = (-d - G.TC) / H; i = int(np.argmax(np.abs(e)))     # 外壁は内壁の下 (d < 0)
        g_ = float(np.abs(e).max()) <= tol_c; ok &= g_
        out["cowl_out_thickness"] = dict(nodes=int(len(P_)), t=G.TC / H, err_min=float(e.min()), err_max=float(e.max()), at=[float(v) / H for v in P_[i]], tol=tol_c, ok=bool(g_))
    for t, L_ in (("sidewall_end", G.LSW), ("cowl_base", G.LCOWL)):
        P_ = need(t)
        if P_ is None: ok = False; continue
        dx = float(np.abs(P_[:, 0] - L_).max()) / H; g_ = dx <= tol_x
        out[t + "_x"] = dict(nodes=int(len(P_)), dx_max=dx, tol=tol_x, ok=bool(g_)); ok &= g_
        if t == "cowl_base":
            ylo, yhi = float(P_[:, 1].min()), float(P_[:, 1].max()); yc_te = float(G.yc(G.LCOWL))
            e_lo, e_hi = abs(ylo - G.yo_te) / H, abs(yhi - yc_te) / H; g_ = e_lo <= tol_x and e_hi <= tol_x; ok &= g_
            out["cowl_base_y"] = dict(y_min=ylo / H, y_out_expected=G.yo_te / H, err_out=e_lo, y_max=yhi / H, y_in_expected=yc_te / H, err_in=e_hi,
                                      height=(yhi - ylo) / H, tol=tol_x, ok=bool(g_))
    out["ok"] = bool(ok); return out


def check(G, P, vols, groups, h1, h1e, out, nl, nl_end, xstart, coupon=False):
    from scipy.spatial import cKDTree
    M = gmsh.model.mesh; ntag, xyz, _ = M.getNodes(); xyz = xyz.reshape(-1, 3)
    idx = np.zeros(int(ntag.max()) + 1, np.int64); idx[ntag.astype(np.int64)] = np.arange(len(ntag))
    other, hx_l, mixed = {}, [], 0
    for v in vols.values():
        et, _, cn = M.getElements(3, v)
        for t, c in zip(et, cn):
            if int(t) != 5: other[int(t)] = other.get(int(t), 0) + len(c); continue
            hv = idx[np.asarray(c, np.int64)].reshape(-1, 8); hx_l.append(hv)
    et2, _, cn2 = M.getElements(2)
    for t, c in zip(et2, cn2):
        if int(t) != 3: other[int(t)] = other.get(int(t), 0) + len(c)
    bq = {}
    for t, surfs in groups.items():
        q = []
        for s_ in surfs:
            _, _, c2 = M.getElements(2, s_); q.append(idx[np.asarray(c2[0], np.int64)].reshape(-1, 4))
        bq[t] = np.vstack(q)
    del idx, ntag
    gmsh.finalize()                                 # msh は書き出し済み。検査は numpy 配列だけで行い、gmsh の格子を先に解放する (ピーク RSS を下げる)
    hx = np.vstack(hx_l); used = np.unique(hx)
    J, J32, sk, AR, Led = hex_quality(xyz, hx); min_edge = float(Led.min()); del Led
    sg = np.sign(np.median(J)); o = 0
    for hv in hx_l:
        s_ = np.sign(J[o:o + len(hv)]); mixed += int(not (np.all(s_ > 0) or np.all(s_ < 0)) or np.sign(s_.flat[0]) != sg); o += len(hv)
    neg_jac, neg_jac_f32 = int((J * sg <= 0).sum()), int((J32 * sg <= 0).sum()); del J, J32
    o = 0; per = {}
    for (ii, b), hv in zip(vols.keys(), hx_l):
        q = per.setdefault(b, dict(skew_max=0.0, ar_max=0.0, n=0)); m = slice(o, o + len(hv)); o += len(hv)
        if sk[m].max() > q["skew_max"]: q["skew_max"] = float(sk[m].max()); q["skew_at_x"] = xstart[ii] / G.H
        if AR[m].max() > q["ar_max"]: q["ar_max"] = float(AR[m].max()); q["ar_at_x"] = xstart[ii] / G.H
        q["n"] += len(hv)
        mm = (AR[m] > 1000) & (sk[m] > 0.30)
        if mm.any():
            q["ar1000_skewed"] = q.get("ar1000_skewed", 0) + int(mm.sum()); i_ = int(np.flatnonzero(mm)[np.argmax(AR[m][mm])])
            cand = (float(AR[m][i_]), float(sk[m][i_]), [float(v) for v in xyz[hv[i_]].mean(0) / G.H], [float(v) for v in np.ptp(xyz[hv[i_]], axis=0)])
            if cand[0] > q.get("ar1000_skewed_worst", (0,))[0]: q["ar1000_skewed_worst"] = cand
    # ---- 位相: 全ヘキサの面を数え、1 回だけ現れる面 = 境界面がタグ付き四角形とちょうど一致すること
    def keys(Q):                                   # 面 (4 節点) -> 並べ替えに使う 2 本の 64 bit キー
        Q = np.sort(Q.astype(np.int64), axis=1); return (Q[:, 0] << 32) | Q[:, 1], (Q[:, 2] << 32) | Q[:, 3]
    h32 = hx.astype(np.int32); k1, k2 = keys(np.vstack([h32[:, f] for f in HF])); od = np.lexsort((k2, k1)); k1, k2 = k1[od], k2[od]; del od
    new = np.r_[True, (k1[1:] != k1[:-1]) | (k2[1:] != k2[:-1])]; st = np.flatnonzero(new); cnt = np.diff(np.r_[st, len(k1)])
    b1, b2 = k1[st[cnt == 1]], k2[st[cnt == 1]]; del k1, k2
    a1, a2 = keys(np.vstack(list(bq.values()))); od = np.lexsort((a2, a1)); a1, a2 = a1[od], a2[od]
    dup = int(((a1[1:] == a1[:-1]) & (a2[1:] == a2[:-1])).sum())
    same = len(b1) == len(a1) and bool(np.all(b1 == a1) and np.all(b2 == a2))
    topo = dict(faces_gt2=int((cnt > 2).sum()), boundary_faces=int(len(b1)), tagged_quads=int(len(a1)), quads_in_two_tags=dup,
                unused_nodes=int(len(xyz) - len(used)), ok=bool(same and (cnt > 2).sum() == 0 and dup == 0))
    # ---- 壁第一層: 壁タグ別・節点別の |n . dx| / h1、第一内部点の資格、全数被覆、実層数、隣接間隔比 (検査の中核は first_layer_check)
    # 期待層数: 壁別に prog_n(DR, h1_w, 1.2) + 4 (端面も同式。plan §6.2 の正式規則・B4b-5 で壁別。build は h1 と nl に壁タグ -> 値の dict を渡す)。判定は接続から数えた実層数
    exp_tags = tuple(t for t in WALL_TAGS if t in bq) if coupon else WALL_TAGS        # 試験片は含む壁タグだけ (本体は全壁タグの存在も要求)
    fl = first_layer_check(xyz, hx, bq, h1, h1e, nl, nl_end, G.DR, classify=model_classifier(G), expected_tags=exp_tags, unit=G.H)
    # ---- 出力実座標の隣接間隔比・後流 Δx・形状ゲート
    names = sorted({b for (_, b) in vols}); bidx = {b: i for i, b in enumerate(names)}
    lab = np.concatenate([np.full(len(hv), bidx[b], np.int16) for (_, b), hv in zip(vols.keys(), hx_l)])
    adj = adjacent_spacing_check(xyz, hx, gmax=P["G"], labels=lab, label_names=names, unit=G.H); del lab
    if coupon: wk = shp = dict(ok=None, note="試験片では対象外")
    else:
        wk = wake_dx_check(xyz, hx, {"cowl": (G.LCOWL, G.LCOWL + G.LWAKE), "sidewall": (G.LSW, G.LSW + G.LWAKE)}, P["DXW"] * G.H, unit=G.H)
        shp = shape_check(xyz, bq, G)
    # ---- 描画用
    dump = {"wall_" + t: xyz[bq[t]] for t in bq if t in WALL_TAGS}
    for xq in P.get("DUMP_X", [0.0, 0.5, 0.8, 1.0, 1.2, 1.5]):
        xv = xq * G.H; tol = 1e-9
        for f in HF:
            Pf = xyz[hx[:, f]]; m = np.all(np.abs(Pf[:, :, 0] - xv) < tol, axis=1)
            if m.any(): dump.setdefault("sec_%.3f" % xq, []).append(Pf[m])
    dump = {k: (np.vstack(v) if isinstance(v, list) else v) for k, v in dump.items()}
    np.savez_compressed(out + "_sections.npz", **dump)
    return dict(nodes=int(len(used)), hexes=int(len(hx)), other_elems=other, neg_jac=neg_jac, neg_jac_f32=neg_jac_f32,
                mixed_sign_volumes=mixed, skew_max=float(sk.max()), skew_p99=float(np.percentile(sk, 99)), skew_gt090=int((sk > 0.9).sum()),
                skew_gt070=int((sk > 0.7).sum()), ar_max=float(AR.max()), ar_gt1000=int((AR > 1000).sum()), ar_gt5000=int((AR > 5000).sum()), ar1000_skew_max=float(sk[AR > 1000].max()) if (AR > 1000).any() else 0.0,
                min_edge=min_edge, per_block=per, topology=topo, first_layer=fl, adjacent_spacing=adj, wake_dx=wk, shape=shp)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--contours", required=True, help="ramp_contour.csv / cowl_contour.csv のあるディレクトリ (export_contours.py で作る)")
    ap.add_argument("--preset", choices=sorted(PRESETS), help="壁別 h1 などの設定の組 (B4b-5、plan §4.16 の暫定の試験値)。--set はその後に上書き")
    ap.add_argument("--set", nargs="*", default=[]); a = ap.parse_args(); P = dict(P0)
    if a.preset: P.update(PRESETS[a.preset]); P["PRESET"] = a.preset
    for kv in a.set: k, v = kv.split("="); P[k] = [float(u) for u in v.split(",")] if "," in v else float(v)
    if a.scale != 1.0:          # 細分列 (plan §6.2): h1 (壁別の H1_* も)・h1e・後流 Δx 上限・HX・HX_FAR を 1/scale、NZ を scale 倍。成長率 G は 1.2 固定 (層数で吸収)
        for k in ["H1", "H1_END", "DXW", "HX", "HX_FAR"] + [k_ for k_ in P if k_.startswith("H1_") and k_ != "H1_END"]: P[k] /= a.scale
        P["NZ"] = int(round((P["NZ"] - 1) * a.scale)) + 1
    r = build(a.out, P, load_contours(a.contours), quiet=False)
    if r["VERDICT"] == "ESTIMATE": print(json.dumps({k: v for k, v in r.items() if k != "conflicts"}, indent=1, ensure_ascii=False))
    sys.exit(0 if r["VERDICT"] in ("PASS", "DIM2", "ESTIMATE") else 1)
