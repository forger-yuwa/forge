"""⑤ SERN の 3D 構造 hex メッシュ (S7: 側壁・有限スパンの確認用)。

2D の 2 バンド構造格子 (mesh_sern) を z 方向に押し出す。z=0 は対称面 (半スパン)、z = W/2 に側壁 (スリット面)、
その外側 z ∈ (W/2, Z_far] は外部流が回り込む空間。ランプ (上境界) は全スパン (機体下面)、カウルは z ≤ W/2 の
有限幅の板 (z > W/2 では中間線ノードを共有 = 板なし)。側壁は x ≤ L_sw の上バンド (カウル〜ランプ) にある。

ノード重複 (スリット):
  D1 カウル: i < i_te かつ z_k ≤ W/2 の中間線ノードを上下 2 重 (2D と同じ)。
  D2 側壁: k = k_sw (z = W/2), i < i_sw (後縁 station は共有), 上バンド j (jm < j ≤ NJ−1、ランプ線含む) を内外 2 重。j = jm では内側 = cowl_in 側の
     上コピー、外側 = 中間線の元ノード (cowl_out 兼)。ランプ線を共有すると入口面で 2 種の入口に属する矛盾ノードになる。
境界面 (quad): inlet_nozzle / inlet_ext / outlet / ramp / top_out / cowl_in / cowl_out / bottom / sym (z=0) /
              side_far (z = Z_far) / sidewall_in / sidewall_out / **vehicle** (幅外の機体下面: x ≤ L_ramp, W/2 < z ≤ W_vehicle/2)。
  R2 (codex M5, 2026-09-13): `ramp` タグはノズル幅内 (z ≤ W/2) だけにし、幅外の機体下面は `vehicle` に分ける。
  `W_vehicle` (全幅/H, None = 遠方境界まで) の外は `top_out` (遠方境界の産物で力の帳簿に入れない)。
hex は gmsh 型 5 (底面 4 点 CCW → 上面 4 点)。座標一致ノードを持つので stage 間 restart は index コピーにすること。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .mesh2d import _radial_fracs
from .mesh_sern import _cluster_stations, _tanh_two_sided, _wake_stations, check_wake_first_spacing

PHYS_SERN3D = {"inlet_nozzle": 1, "inlet_ext": 2, "outlet": 3, "ramp": 4, "cowl_in": 5, "cowl_out": 6, "bottom": 7,
               "top_out": 8, "sym": 9, "side_far": 10, "sidewall_in": 11, "sidewall_out": 12, "fluid": 13, "vehicle": 14, "vehicle_top": 15, "underside_far": 16, "vehicle_side": 17, "vehicle_base": 18}


@dataclass
class SernMesh3DParams:
    ni_up: int = 10
    ni_noz: int = 60
    ni_plume: int = 110
    nj_top: int = 49
    nj_bot: int = 31
    nz_in: int = 25          # z ∈ [0, W/2]
    nz_out: int = 17         # z ∈ (W/2, Z_far]
    W: float = 2.0           # ノズル幅 / H (全幅)
    Z_ext: float = 1.5       # 側壁外側の空間 / H
    L_sw: float | None = None  # 側壁の x 範囲 (None → L_cowl)
    W_vehicle: float | None = None  # 機体の物理幅 / H (全幅)。None = 遠方境界まで機体下面 (旧挙動)。W/2 < z ≤ W_vehicle/2 が vehicle タグ
    # R4 (codex M6, 2026-09-13): ランプ側外部流ブロック (2D の ext_top を z 一様に押し出し)。機体上面 (vehicle_top, 後端テーパ) +
    # 自由流バンド top_depth。ランプ後縁より下流はプルーム上線とノード共有 (2D と同じ)。vehicle_taper > 0 必須 (鉛直 base は node で発散)
    ext_top: bool = False
    top_depth: float = 2.0
    nj_ext_top: int = 41
    vehicle_side: bool = True    # True (R4c, 2026-09-19): z > W/2 の旧ランプ線〜旧機体上面線を**流体で埋め**、
                                 # z = W/2 を機体側面 (vehicle_side, no-slip) にする。False = 旧挙動 (幅外にも
                                 # ノズル形状の固体が残る = 自由流の中に架空の物体。codex plan レビュー C1)
    nj_vside: int = 17           # 機体側面バンドの y 方向点数
    vehicle_clearance: float = 0.06   # 機体上面 = max(ランプ y) + これ (/H)。**物理入力**。
                                      # 旧実装は max(0.02, 3*first_top_frac) で**実効 0.06** だったので、
                                      # 同じ形状を保つためその値を明示既定にした (R4e / codex plan-3 M3)
    t_base: float = 0.02         # 機体後縁ベースの厚み /H。**物理入力** (格子から独立)。R4e 案 (d)。
                                 # `vehicle_clearance` は t_base の根拠にならないので暫定モデル値として扱い、
                                 # 0.01/0.02/0.04 の形状感度は別試験 (codex plan-3 M3)
    first_top_frac: float = 0.02
    vehicle_taper: float = 0.0
    vehicle_wedge_deg: float = 3.0
    ramp_fillet: float = 0.0     # ランプ膨張角部 (x=0, y=1) の丸め半径 /H (2D mesh_sern と同じ; §4.12 の SST 起動対策)。0 = 鋭角
    L_up: float = 0.5
    x_out_extra: float = 2.0
    bot_depth: float = 3.0
    first_wake_frac: float = 0.0  # ベース直後の第一 station 間隔 /H (**絶対値**)。t_base > 0 のとき必須で
                                  # t_base/5 以下 (plan convection-node-wall-reconstruction §4.28)
    first_wall_frac: float = 4.0e-3
    first_wall_frac_far: float = 0.0   # 壁が終わった下流の第一層厚 /H。**0 = ブレンドしない (既定・挙動不変)**。
                                       # 正値なら `first_wall_frac` からここへ滑らかに移す (plan sern-3d §4.41)。
                                       # 下バンドの壁終端は `L_cowl`、上バンドは `L_ramp`
    wall_frac_blend_len: float = 0.5   # ブレンドの**物理長** /H。格子間隔から決めない
    first_z_frac: float = 4.0e-3
    interface_angle: float = 0.0
    top_ext_angle: float = 0.0
    scale: float = 1.0
    cowl_thickness: float = 0.0  # カウル板厚 /H (0 = 厚さ 0 のスリット)。2D と同じ x 分布で TE に向け 0 に絞る。
                                 # 側端 (z = W/2) の外は板が無いので、最後の z セルで厚さ 0 に閉じる
    x_cluster_w: float = 0.15
    x_cluster_a: float = 3.0


def _sern3d_plume(L_cowl, L_ramp, x_out, prm):
    """プルーム区間の station。ベース厚さ `t_base` があるなら、後縁直後の第一間隔を**絶対値**で押さえる。

    3D の R4e は `t_base` が既定 > 0 (有限ベース) なので、ここを外すと
    ベース後流の剪断層を 1 セルで跨いで発散する (plan convection-node-wall-reconstruction §4.28)。
    形状でなく解像度の問題なので、粗い指定は黙って通さず**生成を失敗させる**。"""
    tb = float(getattr(prm, "t_base", 0.0))
    fw = float(getattr(prm, "first_wake_frac", 0.0))
    if tb > 0.0:
        if not (fw > 0.0):
            fw = tb / 5.0   # 未設定は既定 t_base/5 (codex Major 4: 既存 config/テストを壊さない)
            print(f"[mesh_sern3d] first_wake_frac 未設定 → t_base/5 = {fw:g} を使う")
        if fw > tb / 5.0:
            raise ValueError(
                f"mesh_sern3d: first_wake_frac {fw:g} がベース厚さ {tb:g} に対して粗すぎる "
                f"(t_base/5 = {tb/5.0:g} 以下にすること。plan convection-node-wall-reconstruction §4.28)")
    # **ベース直後 (= ランプ後縁 L_ramp) から**細分する。前版は L_cowl を渡しており、
    # カウル後縁が細かくなるだけでベース直後は粗いままだった
    # (codex 2026-09-20 result レビュー Major 3: 要求 0.004 に対し実測 0.171 = 42.75 倍)。
    n = int(prm.ni_plume)
    if not (L_cowl < L_ramp < x_out) or not (tb > 0.0):
        return _wake_stations(L_cowl, x_out, n, 0.0, prm.x_cluster_w, prm.x_cluster_a)
    fa = (L_ramp - L_cowl) / (x_out - L_cowl)
    na = max(int(round(n * fa)), 5)
    nb = max(n - na + 1, 5)
    return np.concatenate([
        _cluster_stations(L_cowl, L_ramp, na, (True, True), prm.x_cluster_w, prm.x_cluster_a),
        _wake_stations(L_ramp, x_out, nb, fw, prm.x_cluster_w, prm.x_cluster_a)[1:]])


def generate_sern_mesh3d(design, prm: SernMesh3DParams):
    L_cowl = float(design.cowl_xy[-1, 0]); y_te = float(design.cowl_xy[-1, 1])
    tan_c = -y_te / L_cowl if L_cowl > 0 else 0.0
    L_ramp = float(design.L_ramp); y_e = float(design.ramp_xy[-1, 1])
    L_sw = L_cowl if prm.L_sw is None else float(prm.L_sw)
    x_out = L_ramp + prm.x_out_extra
    rx, ry = design.ramp_xy[:, 0], design.ramp_xy[:, 1]
    # ランプ角部の丸め (2D mesh_sern と同じ式): 接点 x1 = −t, x2 = t cosθ, t = R tan(θ/2)。丸め区間はブロック境界にする
    R_f = float(prm.ramp_fillet)
    if R_f > 0.0 and len(rx) > 1:
        th_r0 = float(np.arctan2(ry[1] - ry[0], rx[1] - rx[0])); t_f = R_f * np.tan(0.5 * th_r0)
        xf1, xf2, nf = -t_f, t_f * np.cos(th_r0), 7
        xs = np.concatenate([
            _cluster_stations(-prm.L_up, xf1, max(prm.ni_up - 6, 4), (False, True), prm.x_cluster_w, prm.x_cluster_a),
            np.linspace(xf1, 0.0, nf)[1:], np.linspace(0.0, xf2, nf)[1:],
            _cluster_stations(xf2, L_cowl, prm.ni_noz, (True, True), prm.x_cluster_w, prm.x_cluster_a)[1:],
            _sern3d_plume(L_cowl, L_ramp, x_out, prm)[1:],
        ])
    else:
        R_f = 0.0; t_f = xf1 = xf2 = 0.0
        xs = np.concatenate([
            _cluster_stations(-prm.L_up, 0.0, prm.ni_up, (False, True), prm.x_cluster_w, prm.x_cluster_a),
            _cluster_stations(0.0, L_cowl, prm.ni_noz, (True, True), prm.x_cluster_w, prm.x_cluster_a)[1:],
            _sern3d_plume(L_cowl, L_ramp, x_out, prm)[1:],
        ])
    xs[int(np.argmin(np.abs(xs - L_ramp)))] = L_ramp
    check_wake_first_spacing(xs, L_ramp, float(getattr(prm, "t_base", 0.0)),
                             float(getattr(prm, "first_wake_frac", 0.0)))
    i_te = int(np.argmin(np.abs(xs - L_cowl))); assert abs(xs[i_te] - L_cowl) < 1e-12
    i_sw = int(np.argmin(np.abs(xs - L_sw)))

    def y_top(x):
        x = np.asarray(x, dtype=float)
        y = np.where(x < 0.0, 1.0, np.where(x <= L_ramp, np.interp(x, rx, ry), y_e + (x - L_ramp) * np.tan(prm.top_ext_angle)))
        if R_f > 0.0:
            u = x + t_f
            arc = 1.0 + R_f - np.sqrt(np.clip(R_f * R_f - u * u, 0.0, None))
            y = np.where((x >= xf1) & (x <= xf2), np.maximum(y, arc), y)
        return y

    def y_mid(x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0.0, 0.0, np.where(x <= L_cowl, -x * tan_c, y_te + (x - L_cowl) * np.tan(prm.interface_angle)))
    y_bot = float(y_mid(x_out)) - prm.bot_depth
    yt, ym = y_top(xs), y_mid(xs)
    ni, njt, njb = len(xs), prm.nj_top, prm.nj_bot
    NJ = njb + njt - 1; jm = njb - 1
    # z 分布: [0, W/2] は側壁側 (z=W/2) にクラスタ、(W/2, Z_far] は側壁側にクラスタ
    hw = 0.5 * prm.W
    s_in = _radial_fracs(prm.nz_in, min(prm.first_z_frac / hw, 0.5 / (prm.nz_in - 1)))
    z_in = hw * s_in
    if prm.nz_out > 0:
        s_out = _radial_fracs(prm.nz_out, min(prm.first_z_frac / prm.Z_ext, 0.5 / (prm.nz_out - 1)))
        z_out = hw + prm.Z_ext * (1.0 - s_out[::-1])   # 側壁側が細かい
        zs = np.concatenate([z_in, z_out[1:]])
    else:
        zs = z_in                                        # 外側空間なし: z = W/2 が遠方境界 (側壁 = 境界壁)
    nz = len(zs); k_sw = prm.nz_in - 1
    no_outer = prm.nz_out == 0
    # --- 2D 断面 (各 station の y 列) ---
    # カウル板厚 (2D の `cowl_thickness` と同じ法則): 入口から 0.8 L_cowl まで t、TE で 0。
    # 厚さ 0 のスリットは node で双子ノードになり、2D では m6_on/m10_on とも発散した (case/46 run_0035)。
    # 板は z <= W/2 にしか無いので、内側 (k <= k_sw) だけ ym±t/2 に割り、外側は単一の ym にする。
    t_c = float(prm.cowl_thickness)
    tk = (np.interp(xs, [-prm.L_up, 0.8 * L_cowl, L_cowl], [t_c, t_c, 0.0], left=t_c, right=0.0)
          if t_c > 0.0 else np.zeros_like(xs))
    Y2 = np.zeros((ni, NJ))       # 板の外側 (z > W/2) 用: 中間線は単一
    Yin = np.zeros((ni, NJ))      # 板の内側 (z <= W/2) 用: 中間線は ym ± t/2
    # --- 壁第 1 層の x ブレンド (plan sern-3d §4.41) ---
    # `first_wall_frac` は全 station に効くので、**壁が終わった下流にも壁用の細層**を敷いていた。
    # R5m の実測では AR > 5000 セルの 57.8 % がそこで、外側 z を細分しても max AR は下がらなかった。
    # 壁の終端 (下バンド = `L_cowl`、上バンド = `L_ramp`) から**物理長** `wall_frac_blend_len` かけて
    # `first_wall_frac_far` へ smoothstep で移す。間隔は 25 倍も跨ぐので**対数補間**。
    _fw = float(prm.first_wall_frac)
    _ff = float(getattr(prm, "first_wall_frac_far", 0.0) or 0.0)
    _bl = max(float(getattr(prm, "wall_frac_blend_len", 0.5)), 1.0e-12)

    def _first_at(x, x_end):
        if not (_ff > 0.0) or _ff == _fw:
            return _fw                      # 既定: ブレンドしない (既存メッシュはビット一致)
        t = min(max((x - x_end) / _bl, 0.0), 1.0)
        w = t * t * (3.0 - 2.0 * t)
        return float(np.exp((1.0 - w) * np.log(_fw) + w * np.log(_ff)))

    # **ブレンドが急すぎると skew が出る**ので生成を失敗させる (2026-09-21)。
    # 判定量は **第一層厚の x 方向勾配** df/dx (無次元)。これは第 1 j 線が壁からずれる角の tan で、
    # skew を直接押し上げる。**station あたりの比では駄目**: 遠方は Δx が大きいので比が勝手に上がり、
    # 実際には skew が出ていないところで落ちる (生産 config, blend 3.0 で上バンド比 3.13 でも skew 0.701)。
    # 較正 (生産 config, `first_wall_frac` 4e-05 → far 4e-03):
    #   blend 0.5 下バンド df/dx **0.0111** → skew max 0.933 (>0.90 が 0.10 %) = 不可
    #   blend 3.0 下バンド 0.00098 / 上バンド 0.00162 → skew max 0.701 (>0.90 が 0) = 可
    # 閾値 0.004 は両者の間 (可の 2.5 倍・不可の 1/2.8)。最終判定は `check_mesh_quality` が行う。
    _flo = np.array([_first_at(x, L_cowl) for x in xs])
    _fup = np.array([_first_at(x, L_ramp) for x in xs])
    if _ff > 0.0 and _ff != _fw:
        _dx = np.diff(xs)
        _g = max(float(np.max(np.abs(np.diff(_flo)) / np.maximum(_dx, 1e-30))),
                 float(np.max(np.abs(np.diff(_fup)) / np.maximum(_dx, 1e-30))))
        if _g > 4.0e-3:
            raise ValueError(
                f"mesh_sern3d: 壁第 1 層の x ブレンドが急すぎる (df/dx = {_g:.4g} > 4.0e-3)。"
                f"`wall_frac_blend_len` を大きくするか station 数を増やすこと。"
                f"急なブレンドは skew を生む (df/dx 0.0111 で skew max 0.933; plan sern-3d §4.41)")
    for i in range(ni):
        f_lo = _flo[i]     # 下バンドの細端 = 中間線 (カウル)
        f_up = _fup[i]     # 上バンドの細端 = 上線 (ランプ)。両側 tanh なので片側で決める
        for Y, lo, up in ((Y2, ym[i], ym[i]), (Yin, ym[i] - 0.5 * tk[i], ym[i] + 0.5 * tk[i])):
            h_lo = max(lo - y_bot, 1e-12); h_up = max(yt[i] - up, 1e-12)
            s_bot = _radial_fracs(njb, min(f_lo / h_lo, 0.5 / (njb - 1)))
            s_top = _tanh_two_sided(njt, min(f_up / h_up, 0.5 / (njt - 1)))
            Y[i, :njb] = y_bot + s_bot * (lo - y_bot)
            Y[i, jm:] = up + s_top * (yt[i] - up)
            # j = jm は**下バンドの壁ノード (cowl_out)**。上の 2 行では `up` で上書きされ、
            # `Yin[i,jm]` が上面になって板厚が消えていた (cowl_in と cowl_out が同一座標の
            # **厚さ 0 スリット**。`cowl_thickness` が 3D で効いていなかった。2026-09-20 修正)。
            # 代入は `lo` **そのもの** (`y_bot + 1.0*(lo-y_bot)` は丸めで一致しない)。
            Y[i, jm] = lo
    # --- ノード番号 ---
    N_base = ni * NJ * nz
    def base(i, j, k): return (i * NJ + j) * nz + k
    dup1 = {}   # (i,k) -> id  (カウル上コピー)
    dup2 = {}   # (i,j) -> id  (側壁外コピー, k=k_sw)
    nid = N_base
    for i in range(i_te):
        for k in range(k_sw + 1):
            dup1[(i, k)] = nid; nid += 1
    for i in (range(0) if no_outer else range(i_sw)):   # 側壁後縁 (i_sw) は共有 (カウル TE と同じ)。外側空間なしなら重複なし
        for j in range(jm + 1, NJ):            # ランプ線 (j = NJ−1) も内外 2 重: 共有すると入口面で
            dup2[(i, j)] = nid; nid += 1       # inlet_nozzle と inlet_ext の両方に属し step 3 で発散した
    # 板厚の z 分布: 側壁 (z = W/2) は**厚さ 0 のスリット**なので、そこまで板厚を効かせると
    # スリットの内外で中間線の高さがずれ、側壁後縁が破綻する (run_0086: 暖機段 step 4 で ro NaN)。
    # 最後の 2 セルで厚さを 0 に絞り、k_sw では内外が一致するようにする。
    sz = np.ones(nz)
    if t_c > 0.0:
        kt = max(k_sw - 2, 0)
        if k_sw > kt:
            u = (zs[kt:k_sw + 1] - zs[kt]) / max(zs[k_sw] - zs[kt], 1e-30)
            sz[kt:k_sw + 1] = 1.0 - u * u * (3.0 - 2.0 * u)
        sz[k_sw + 1:] = 0.0
    # 板の半厚 (x station, z ごと)。**float32 で解けない厚みは 0 に丸める**: 変換器は座標を float32 で
    # 書くので、解けない隙間を残すと双子ノード (cowl_in / cowl_out) が「ほぼ一致」になり、
    # 衝突するかどうかが丸めに依存する。**設計どおり閉じるところは厳密に閉じる**
    # ([[axisym-rweight-closure-fp32]] と同じ方針)。閉じた列は上下とも同じ量だけ動かさない。
    _half = 0.5 * tk[:, None] * sz[None, :]                       # (ni, nz)
    _closed = _half < (1.0e-5 * np.maximum(1.0, np.abs(ym))[:, None])
    _szw = np.where(_closed, 0.0, 1.0) * sz[None, :]              # (ni, nz) 実効 sz
    coords = np.zeros((nid, 3))
    for i in range(ni):
        for j in range(NJ):
            b = base(i, j, 0)
            coords[b:b + nz, 0] = xs[i]; coords[b:b + nz, 2] = zs
            # base(i, jm, k) は「下側 (cowl_out)」。板がある内側 (k <= k_sw) だけ Yin、外側は Y2 (= 単一の中間線)。
            coords[b:b + nz, 1] = Y2[i, j]
            if t_c > 0.0 and tk[i] > 0.0:
                coords[b:b + nz, 1] += (Yin[i, j] - Y2[i, j]) * _szw[i]
    n_dup_cowl_closed = 0
    for (i, k), n in dup1.items():
        coords[n] = (xs[i], ym[i] + 0.5 * tk[i] * _szw[i, k], zs[k])   # 上側 (cowl_in)
        if t_c <= 0.0 or tk[i] <= 0.0 or _closed[i, k]:
            n_dup_cowl_closed += 1                 # 座標一致の双子 (設計どおり閉じた列)
    for (i, j), n in dup2.items():
        coords[n] = (xs[i], Y2[i, j], zs[k_sw])

    def node(i, j, k, side: str):
        """side: 'lo' (下バンド), 'up_in' (上バンド, 内側 z<=W/2), 'up_out' (上バンド, 外側 z>W/2)。"""
        if j == jm:
            if side == "up_in" and (i, k) in dup1:
                return dup1[(i, k)]
            return base(i, j, k)
        if side == "up_out" and k == k_sw and (i, j) in dup2:
            return dup2[(i, j)]
        return base(i, j, k)

    hexes = []
    for i in range(ni - 1):
        for k in range(nz - 1):
            for j in range(NJ - 1):
                if j < jm:
                    sd = ("lo", "lo")
                elif k < k_sw:
                    sd = ("up_in", "up_in")
                elif k == k_sw:
                    sd = ("up_in", "up_out")    # k=k_sw 側は内側面 (境界 k) と外側 (k+1)
                else:
                    sd = ("up_out", "up_out")
                # k = k_sw の cell (k_sw → k_sw+1) は外側 cell: 面 k_sw のノードは outer コピー
                s0, s1 = ("up_out", "up_out") if (j >= jm and k == k_sw) else sd
                n = lambda ii, jj, kk, ss: node(ii, jj, kk, ss)
                hexes.append((n(i, j, k, s0), n(i + 1, j, k, s0), n(i + 1, j + 1, k, s0), n(i, j + 1, k, s0),
                              n(i, j, k + 1, s1), n(i + 1, j, k + 1, s1), n(i + 1, j + 1, k + 1, s1), n(i, j + 1, k + 1, s1)))
    hexes = np.asarray(hexes, dtype=np.int64)
    # --- 境界 quad ---
    B = {n: [] for n in PHYS_SERN3D if n != "fluid"}
    for k in range(nz - 1):
        side_k = "up_in" if k < k_sw else "up_out"
        side_k1 = "up_in" if k + 1 <= k_sw and k + 1 != k_sw else ("up_in" if k + 1 < k_sw else "up_out")
        for j in range(NJ - 1):
            sd = "lo" if j < jm else ("up_in" if k < k_sw else "up_out")
            sd1 = "lo" if j < jm else ("up_in" if k + 1 < k_sw or (k + 1 == k_sw and k < k_sw) else "up_out")
            # k+1 == k_sw の面ノードは内側 (inner) — 外側 cell は k = k_sw から
            # 入口面: ノズル幅 (z ≤ W/2) の上バンドだけが燃焼器出口。側壁の外側 (k ≥ k_sw) は上バンドでも外部流
            nm = "inlet_ext" if (j < jm or k >= k_sw) else "inlet_nozzle"
            B[nm].append((node(0, j, k, sd), node(0, j + 1, k, sd), node(0, j + 1, k + 1, sd1), node(0, j, k + 1, sd1)))
            B["outlet"].append((node(ni - 1, j, k, sd), node(ni - 1, j, k + 1, sd1), node(ni - 1, j + 1, k + 1, sd1), node(ni - 1, j + 1, k, sd)))
    vs_on = bool(prm.ext_top and prm.vehicle_side and not no_outer)
    for i in range(ni - 1):
        for k in range(nz - 1):
            sdk = "up_in" if k < k_sw else "up_out"; sdk1 = "up_in" if k + 1 < k_sw or k + 1 == k_sw and k < k_sw else "up_out"
            B["bottom"].append((base(i, 0, k), base(i + 1, 0, k), base(i + 1, 0, k + 1), base(i, 0, k + 1)))
            xm = 0.5 * (xs[i] + xs[i + 1]); zm = 0.5 * (zs[k] + zs[k + 1])
            top = (node(i, NJ - 1, k, sdk), node(i, NJ - 1, k + 1, sdk1), node(i + 1, NJ - 1, k + 1, sdk1), node(i + 1, NJ - 1, k, sdk))
            if xm > L_ramp:
                if not prm.ext_top:
                    B["top_out"].append(top)      # ext_top ではプルーム上線は上のブロックとの内部面
            elif k < k_sw:                       # ノズル幅内 (z ≤ W/2) = ramp
                B["ramp"].append(top)
            elif vs_on:
                pass                             # R4e: 幅外はバンドが**全長**を覆うので常に内部面
            elif prm.W_vehicle is None or zm <= 0.5 * float(prm.W_vehicle):
                B["vehicle"].append(top)         # 幅外の機体下面 (R2、旧挙動)
            else:
                B["underside_far"].append(top)   # 機体幅の外 (z > W_vehicle/2) = 遠方境界の産物。slip 壁・帳簿外
            if i + 1 <= i_te and k + 1 <= k_sw:
                B["cowl_in"].append((node(i, jm, k, "up_in"), node(i + 1, jm, k, "up_in"), node(i + 1, jm, k + 1, "up_in"), node(i, jm, k + 1, "up_in")))
                B["cowl_out"].append((base(i, jm, k), base(i, jm, k + 1), base(i + 1, jm, k + 1), base(i + 1, jm, k)))
        for j in range(NJ - 1):
            sd = "lo" if j < jm else "up_in"
            B["sym"].append((node(i, j, 0, sd), node(i, j + 1, 0, sd), node(i + 1, j + 1, 0, sd), node(i + 1, j, 0, sd)))
            sdo = "lo" if j < jm else ("up_in" if no_outer else "up_out")
            far_q = (node(i, j, nz - 1, sdo), node(i + 1, j, nz - 1, sdo), node(i + 1, j + 1, nz - 1, sdo), node(i, j + 1, nz - 1, sdo))
            if no_outer and j >= jm and i + 1 <= i_sw:
                B["sidewall_in"].append(far_q)          # 外側空間なし: 遠方境界のノズル区間が側壁 (境界壁)
            else:
                B["side_far"].append(far_q)
            if (not no_outer) and i + 1 <= i_sw and j >= jm:
                B["sidewall_in"].append((node(i, j, k_sw, "up_in"), node(i + 1, j, k_sw, "up_in"), node(i + 1, j + 1, k_sw, "up_in"), node(i, j + 1, k_sw, "up_in")))
                B["sidewall_out"].append((node(i, j, k_sw, "up_out"), node(i, j + 1, k_sw, "up_out"), node(i + 1, j + 1, k_sw, "up_out"), node(i + 1, j, k_sw, "up_out")))
    ext = {"ext_top": bool(prm.ext_top)}
    if prm.ext_top:
        coords, hexes, top_fn = _add_ext_top3d(coords, hexes, B, xs, yt, zs, L_ramp, y_e, node, NJ, prm, ext,
                                               k_sw=k_sw, vs_on=bool(prm.vehicle_side and not no_outer))
        if prm.vehicle_side and not no_outer:
            coords, hexes = _add_vehicle_side3d(coords, hexes, B, xs, yt, zs, k_sw, node, top_fn, NJ, prm, ext, L_ramp, y_e)
    coords *= prm.scale
    info = {"ni": ni, "NJ": NJ, "nz": nz, "jm": jm, "k_sw": k_sw, "i_te": i_te, "i_sw": i_sw, "cells": int(hexes.shape[0]),
            "n_dup_cowl_closed": int(n_dup_cowl_closed),
            "first_wall_frac_far": float(_ff), "wall_frac_blend_len": float(_bl),
            "nodes": int(coords.shape[0]), "W": prm.W, "Z_far": float(zs[-1]), "L_sw": L_sw, "cowl_thickness": t_c, "x_out": x_out, "y_bot": y_bot,
            "L_cowl": L_cowl, "L_ramp": L_ramp, "n_dup_cowl": len(dup1), "n_dup_side": len(dup2),
            "W_vehicle": prm.W_vehicle, "n_vehicle_faces": len(B["vehicle"]), "ramp_fillet": R_f, **ext}
    return coords, hexes, B, info, y_mid


def _vehicle_top_line(xs, yt, L_ramp, y_e, prm):
    """バンド上端線 b(x) と、機体後縁 station k_r を返す (R4e 案 (d), codex plan-3 M1/M3)。

    `b(x)` は**全 x で定義**され、後縁で機体ベースの厚み `t_base` を持つ:

    - `x < x0` … `y_veh` = ランプ最大 y + `vehicle_clearance` (平らな機体上面)
    - テーパ区間 … 3 次エルミートで `y_e + t_base` に着地
    - `x > L_ramp` … `yt(x) + t_base` (後流ブロック上端として**連続に延長**)

    どこでも `b ≥ yt + t_base` に下限を掛けるので、**厚さ 0 に潰れる区間が無い** —
    旧実装の `i_end` (厚さが `first_wall_frac` を切る station) は不要になった。

    **形状は格子間隔に依存させない** (codex plan-3 M3): 旧実装は
    `clearance = max(vehicle_clearance, 3*first_top_frac)` で、`first_top_frac` を
    0.02 → 0.01 にするだけで機体上面が 0.03 H 動いた = 格子独立性試験が成立しない。
    いまはクリアランスを物理値のみで決め、**格子が足りなければ生成を失敗させる**。
    """
    clr = float(prm.vehicle_clearance)
    tb = float(prm.t_base)
    if tb <= 0.0:
        raise ValueError("mesh_sern3d: t_base > 0 が要る (R4e 案 (d) は有限厚のベースで終わる)")
    if float(prm.first_top_frac) > clr / 3.0:
        raise ValueError(
            f"mesh_sern3d: first_top_frac {prm.first_top_frac:g} が vehicle_clearance {clr:g} に対して粗すぎる "
            f"(clr/3 = {clr / 3.0:g} 以下が要る)。**形状を格子に合わせて動かさない**ので、格子側を細かくすること")
    ni = len(xs)
    k_r = int(np.argmin(np.abs(xs - L_ramp)))
    y_veh = float(yt[:k_r + 1].max()) + clr
    taper_len = float(prm.vehicle_taper) * L_ramp
    x0 = L_ramp - taper_len
    m_e = (yt[k_r] - yt[k_r - 1]) / max(xs[k_r] - xs[k_r - 1], 1e-12)
    m1 = min(m_e, 0.0) - np.tan(np.radians(float(prm.vehicle_wedge_deg)))
    y_end = y_e + tb                       # 後縁は厚み t_base で終わる (= yt[k_r] + tb)
    _s = np.clip((xs - x0) / max(taper_len, 1e-12), 0.0, 1.0)
    b = ((2.0 * _s ** 3 - 3.0 * _s ** 2 + 1.0) * y_veh + (-2.0 * _s ** 3 + 3.0 * _s ** 2) * y_end
         + (_s ** 3 - _s ** 2) * taper_len * m1)
    b = np.where(xs < x0, y_veh, b)
    b = np.maximum(b, yt + tb)             # どこでも最低 t_base の厚み (退化区間を作らない)
    b[k_r] = y_end
    b = np.where(np.arange(ni) <= k_r, b, yt + tb)   # 後縁より下流 = プルーム線 + t_base (x=L_ramp で連続)
    return b, y_veh, k_r


def _add_ext_top3d(coords, hexes, B, xs, yt, zs, L_ramp, y_e, node, NJ, prm, ext, k_sw=None, vs_on=False):
    """ランプ側外部流ブロック: バンド上端線 b(x) から天井までを z 一様に押し出す (R4e 案 (d))。

    **後縁でブロックの担当を切り替えない** (codex plan-3 M1)。旧実装は `i_end` より下流で j=0 を
    プルーム上線と共有していたため、幅内と幅外で同じ領域を別の節点分布が担当し、側面バンドの末端を
    人工的な壁で閉じる必要があった (幅外の流体が抜けられず 5.72 MPa)。いまは j=0 が常に b(x) を持ち、
    その下 (yt..b) は幅外では全長バンド、幅内では後縁下流の後流ブロックが担当する。

    境界: `vehicle_top` (x ≤ L_ramp かつ z ≤ W/2 の下面), `top_out` (上面), `inlet_ext` / `outlet` / `sym` / `side_far`。"""
    from .mesh2d import _radial_fracs  # noqa: F401  (既存 import と同じ経路)
    ni, nz, njT = len(xs), len(zs), int(prm.nj_ext_top)
    if not prm.vehicle_taper > 0.0:
        raise ValueError("mesh_sern3d ext_top は vehicle_taper > 0 (テーパ版) のみ (鉛直 base + wake は node で発散するので未実装)")
    y3, y_veh, k_r = _vehicle_top_line(xs, yt, L_ramp, y_e, prm)
    assert abs(xs[k_r] - L_ramp) < 1e-12
    N0 = coords.shape[0]
    n_own0 = ni * nz                                      # j=0 (= b(x)) は **全 station で自前**
    N_topj = ni * (njT - 1) * nz

    def top(i, j, k):
        if j == 0:
            return N0 + i * nz + k
        return N0 + n_own0 + ((i * (njT - 1)) + (j - 1)) * nz + k
    new = np.zeros((n_own0 + N_topj, 3))
    s_t = _geom_start3(njT, min(prm.first_top_frac / prm.top_depth, 0.5))
    for i in range(ni):
        yy = y3[i] + s_t * prm.top_depth
        for k in range(nz):
            new[top(i, 0, k) - N0] = (xs[i], y3[i], zs[k])
            for j in range(1, njT):
                new[top(i, j, k) - N0] = (xs[i], yy[j], zs[k])
    coords = np.vstack([coords, new])
    extra = []
    for i in range(ni - 1):
        for k in range(nz - 1):
            for j in range(njT - 1):
                extra.append((top(i, j, k), top(i + 1, j, k), top(i + 1, j + 1, k), top(i, j + 1, k),
                              top(i, j, k + 1), top(i + 1, j, k + 1), top(i + 1, j + 1, k + 1), top(i, j + 1, k + 1)))
    hexes = np.vstack([hexes, np.asarray(extra, dtype=np.int64)])
    B.setdefault("vehicle_top", [])
    for i in range(ni - 1):
        for k in range(nz - 1):
            # 機体上面 = 幅内 (z ≤ W/2) かつ後縁より上流だけ。幅外はバンド、後縁より下流は後流ブロックとの
            # **内部面**になる (R4e 案 (d): 担当の切り替えが無いので条件が 1 本になった)
            if i + 1 <= k_r and (k < k_sw or not vs_on):
                B["vehicle_top"].append((top(i, 0, k), top(i, 0, k + 1), top(i + 1, 0, k + 1), top(i + 1, 0, k)))
            B["top_out"].append((top(i, njT - 1, k), top(i + 1, njT - 1, k), top(i + 1, njT - 1, k + 1), top(i, njT - 1, k + 1)))
        for j in range(njT - 1):
            B["sym"].append((top(i, j, 0), top(i, j + 1, 0), top(i + 1, j + 1, 0), top(i + 1, j, 0)))
            B["side_far"].append((top(i, j, nz - 1), top(i + 1, j, nz - 1), top(i + 1, j + 1, nz - 1), top(i, j + 1, nz - 1)))
    for k in range(nz - 1):
        for j in range(njT - 1):
            B["inlet_ext"].append((top(0, j, k), top(0, j + 1, k), top(0, j + 1, k + 1), top(0, j, k + 1)))
            B["outlet"].append((top(ni - 1, j, k), top(ni - 1, j, k + 1), top(ni - 1, j + 1, k + 1), top(ni - 1, j + 1, k)))
    ext.update({"nj_ext_top": njT, "y_veh": y_veh, "i_ramp_te": int(k_r), "t_base": float(prm.t_base), "vehicle_taper": float(prm.vehicle_taper),
                "vehicle_wedge_deg": float(prm.vehicle_wedge_deg), "top_depth": float(prm.top_depth), "n_top_nodes": int(n_own0 + N_topj)})
    return coords, hexes, top


def _add_vehicle_side3d(coords, hexes, B, xs, yt, zs, k_sw, node, top_fn, NJ, prm, ext, L_ramp, y_e):
    """機体まわりのバンド a(x)=yt .. b(x) を張る (R4e 案 (d), codex plan-3 M1)。

    **後縁でブロックの担当を切り替えない**のが要点。旧実装は幅外バンドを機体終端で打ち切り、
    そこを `vehicle_side` の壁で閉じていた。その壁は z=W/2 に限らず遠方境界まで延びる 224 面で、
    **幅外の流体が下流へ抜けられず 5.72 MPa** (外気 2851 Pa) まで溜まっていた (旧 run_0122 はこの影響下)。

    いまの割り当て:

    - **幅外** (z > W/2): `a..b` の流体バンドを**全長**に置く。末端を閉じる面は要らない。
    - **幅内** (z ≤ W/2): 同じ `a..b` のブロックを**後縁より下流だけ**に置く (後流ブロック)。
    - 機体 (幅内・後縁より上流の `a..b`) は**流体を張らない** = そこが機体構造。
    - `vehicle_base`: x = L_ramp の幅内 `a..b` (厚み `t_base` の物理的なベース面)。
    - `vehicle_side`: z = W/2 の幅内側面 (後縁より上流)。

    節点は共通分布 `y_j = a + η_j (b−a)` と**同じ節点 ID** を使うので、幅外バンドと後流ブロックは適合する。
    """
    b_line, _yv, k_r = _vehicle_top_line(xs, yt, L_ramp, y_e, prm)
    ni, nz = len(xs), len(zs)
    njv = int(prm.nj_vside)
    if njv < 3:
        raise ValueError("nj_vside は 3 以上")

    def k0(i):
        """station i でバンドが始まる z index。後縁より上流は幅外のみ、下流は全幅。"""
        return 0 if i >= k_r else k_sw

    # 内部行 (j=0 は noz バンド上端、j=njv-1 は top バンド j=0 と共有) の節点を i ごとに詰めて確保する
    starts = np.zeros(ni + 1, dtype=np.int64)
    for i in range(ni):
        starts[i + 1] = starts[i] + (nz - k0(i)) * (njv - 2)
    N0 = coords.shape[0]
    n_int = int(starts[-1])

    def vs(i, j, k):
        if j == 0:
            return node(i, NJ - 1, k, "up_out" if k >= k_sw else "up_in")
        if j == njv - 1:
            return top_fn(i, 0, k)
        return int(N0 + starts[i] + (k - k0(i)) * (njv - 2) + (j - 1))

    new_xyz = np.zeros((n_int, 3))
    for i in range(ni):
        h = max(b_line[i] - yt[i], 0.0)
        # 両端とも内部面なので壁クラスタは不要。下端 (noz バンド上端) のセル厚に合わせて伸ばす
        frac = (_geom_start3(njv, min(prm.first_wall_frac / h, 0.5)) if h > 1e-12
                else np.linspace(0.0, 1.0, njv))
        yy = yt[i] + frac * h
        for k in range(k0(i), nz):
            for j in range(1, njv - 1):
                new_xyz[vs(i, j, k) - N0] = (xs[i], yy[j], zs[k])
    coords = np.vstack([coords, new_xyz])

    extra = []
    for i in range(ni - 1):
        kb = max(k0(i), k0(i + 1))          # 両 station に在る z 範囲 (= k0(i)。k0 は非増加)
        for k in range(kb, nz - 1):
            for j in range(njv - 1):
                extra.append((vs(i, j, k), vs(i + 1, j, k), vs(i + 1, j + 1, k), vs(i, j + 1, k),
                              vs(i, j, k + 1), vs(i + 1, j, k + 1), vs(i + 1, j + 1, k + 1), vs(i, j + 1, k + 1)))
    if extra:
        hexes = np.vstack([hexes, np.asarray(extra, dtype=np.int64)])

    B.setdefault("vehicle_side", [])
    B.setdefault("vehicle_base", [])
    for j in range(njv - 1):
        # 機体ベース (x = L_ramp、幅内)。上流側に機体が在るので -x 向きの境界面
        for k in range(0, k_sw):
            B["vehicle_base"].append((vs(k_r, j, k), vs(k_r, j + 1, k), vs(k_r, j + 1, k + 1), vs(k_r, j, k + 1)))
        # 機体側面 (z = W/2、後縁より上流)。バンドの cell は +z 側にある
        for i in range(k_r):
            B["vehicle_side"].append((vs(i, j, k_sw), vs(i, j + 1, k_sw), vs(i + 1, j + 1, k_sw), vs(i + 1, j, k_sw)))
    for i in range(ni - 1):
        kb = max(k0(i), k0(i + 1))
        for j in range(njv - 1):
            B["side_far"].append((vs(i, j, nz - 1), vs(i + 1, j, nz - 1), vs(i + 1, j + 1, nz - 1), vs(i, j + 1, nz - 1)))
            if kb == 0:                      # 後縁より下流は対称面まで在る
                B["sym"].append((vs(i, j, 0), vs(i, j + 1, 0), vs(i + 1, j + 1, 0), vs(i + 1, j, 0)))
    for k in range(k0(0), nz - 1):
        for j in range(njv - 1):
            B["inlet_ext"].append((vs(0, j, k), vs(0, j + 1, k), vs(0, j + 1, k + 1), vs(0, j, k + 1)))
    for k in range(k0(ni - 1), nz - 1):
        for j in range(njv - 1):
            B["outlet"].append((vs(ni - 1, j, k), vs(ni - 1, j, k + 1), vs(ni - 1, j + 1, k + 1), vs(ni - 1, j + 1, k)))
    ext.update({"nj_vside": njv, "n_vside_nodes": int(n_int), "n_vside_faces": len(B["vehicle_side"]),
                "n_vbase_faces": len(B["vehicle_base"]), "i_ramp_te": int(k_r)})
    return coords, hexes


def _geom_start3(n, first):
    """[0,1] を n 点で、最初の間隔が first になる等比分布 (mesh_sern._geom_start と同じ)。"""
    from .mesh_sern import _geom_start
    return _geom_start(n, first)


def write_msh41_3d(path, coords, hexes, bquads: dict, phys: dict) -> None:
    names = [n for n in phys if n != "fluid" and len(bquads.get(n, [])) > 0]
    n_nodes = coords.shape[0]
    n_elems = hexes.shape[0] + sum(len(bquads[g]) for g in names)
    mn, mx = coords.min(0), coords.max(0)
    L = []; ap = L.append
    ap("$MeshFormat\n4.1 0 8\n$EndMeshFormat")
    ap(f"$PhysicalNames\n{len(names) + 1}")
    for nm in names:
        ap(f'2 {phys[nm]} "{nm}"')
    ap(f'3 {phys["fluid"]} "fluid"')
    ap("$EndPhysicalNames")
    ap("$Entities"); ap(f"0 0 {len(names)} 1")
    bb = f"{mn[0]:.9g} {mn[1]:.9g} {mn[2]:.9g} {mx[0]:.9g} {mx[1]:.9g} {mx[2]:.9g}"
    for si, nm in enumerate(names, start=1):
        ap(f"{si} {bb} 1 {phys[nm]} 0")
    ap(f"1 {bb} 1 {phys['fluid']} 0")
    ap("$EndEntities")
    ap("$Nodes"); ap(f"1 {n_nodes} 1 {n_nodes}"); ap(f"3 1 0 {n_nodes}")
    ap("\n".join(str(i + 1) for i in range(n_nodes)))
    ap("\n".join(f"{c[0]:.10g} {c[1]:.10g} {c[2]:.10g}" for c in coords))
    ap("$EndNodes")
    ap("$Elements"); ap(f"{1 + len(names)} {n_elems} 1 {n_elems}")
    et = 1
    ap(f"3 1 5 {hexes.shape[0]}")
    ap("\n".join(f"{et + k} " + " ".join(str(v + 1) for v in h) for k, h in enumerate(hexes)))
    et += hexes.shape[0]
    for si, nm in enumerate(names, start=1):
        q = bquads[nm]
        ap(f"2 {si} 3 {len(q)}")
        ap("\n".join(f"{et + k} {a+1} {b+1} {c+1} {d+1}" for k, (a, b, c, d) in enumerate(q)))
        et += len(q)
    ap("$EndElements")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
