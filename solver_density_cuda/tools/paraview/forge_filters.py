"""forge の res_*.xmf / res_*.h5 用 ParaView Python プラグイン。

提供フィルタ:
  Forge Derived Quantities  … マッハ数 / シュリーレン (無次元密度勾配) /
                              Q 値 / ヘリシティ を一括で計算して出力に追加する。
  Forge Saturation          … 凝縮種 (H2O / N2) の蒸気分圧・飽和蒸気圧・飽和温度・
                              過冷却度・過飽和度を計算する (forge の condS/condTsat と同じ式)。

読み込み方法 (ParaView GUI):
  Tools > Manage Plugins > Load New... で本ファイルを選び、
  必要なら "Auto Load" を有効にする。以後 Filters > Alphabetical に
  "Forge Derived Quantities" が現れる。
  "Failed to import paraview.detail.pythonalgorithm" で失敗する環境
  (ParaView 5.11 + Python 3.12 の組み合わせ。ParaView 側が Python 3.11 で削除された
  inspect.getargspec を import するため) では、代わりに同ディレクトリの
  macro_load_forge_filters.py を Macros > Add new macro... で登録し、そのボタンから読み込む。

pvpython / pvbatch から使う場合:
  from paraview.simple import *
  LoadPlugin("<repo>/solver_density_cuda/tools/paraview/forge_filters.py", ns=globals())
  src = XDMFReader(FileNames=["res_12000.xmf"])
  drv = ForgeDerivedQuantities(Input=src)

入力に期待する配列 (cell モードは Cell Data、node モードは Point Data):
  必須: ro, Ux, Uy, Uz          (速度ベクトル U と各種微分量の元)
  任意: P, sonic                (マッハ数)
  任意: dUxdx…dUzdz, drodx…drodz (ソルバ自身の勾配。あれば既定で優先使用)
勾配配列が無い場合は vtkGradientFilter でメッシュから勾配を計算する。

出力配列:
  U, Mach, sound_speed
  grad_ro, schlieren_mag, schlieren_dir, schlieren
  vorticity, vorticity_mag, Q, Q_norm, (lambda2)
  helicity, helicity_norm

Forge Saturation (凝縮 ON/OFF どちらの run にも使える後処理):
  入力: T (必須), ro, P, 蒸気質量分率配列, 液相質量分率配列 (既定 g_0, 無ければ 0)
  蒸気配列の決め方 (Y1 を既定採用しない: 5 種順序では Y1=N2 を水蒸気と誤計算するため):
    - `Run Config (solverConfig.yaml path)` に run の solverConfig.yaml を指定すると、`physProp.species` と
      `condensation.condensationSpecies` / `condGasSpecies` (無ければ H2O) から `Y{index}` を名前で解決する
      (tools/forge_species.py)。
    - 指定しないときは `Vapor Mass Fraction Array` を明示する (空ならエラー)。空気凝縮 (CPG carrier, 配列なし) は
      Vapor Mass Fraction Array を "none" にして Vapor Mass Fraction Constant を使う。
  蒸気分圧 p_v = ro (Y_v − g) R_v T (carrier 形, forge cond_vapor_state と同一)。蒸気配列が無く定数も 0 なら
  純蒸気 (p_v = P)。空気凝縮 (CPG carrier) は Vapor Mass Fraction Constant に condVaporMassFraction (0.7671) を入れる。
  出力: p_vapor, p_sat, T_sat, subcooling (= T_sat − T; 正が過冷却), supersaturation (S = p_v / p_sat(T)), log10_S
        (+ Compute Ice で氷基準 p_sat_ice, S_ice, T_sat_ice, subcooling_ice; H2O のみ)
  p_sat: H2O = Murphy & Koop (2005) 過冷却液 (forge h2o_psat)、氷 = 同 ice 式。
         N2 = Jacobsen 液 + 50 K 未満の Clausius–Clapeyron 外挿 (forge n2_psat_ex; 低温整合オプション付き)。
  T_sat: p_sat(T_sat) = p_v の Newton 反転 (forge cond_Tsat と同じ初期値・クリップ)。p_v ≤ 1e-6 Pa は 0。
"""

import inspect

# ParaView 5.11 系の paraview.detail.pythonalgorithm は Python 3.11+ で
# 削除された inspect.getargspec を import するため、ここで補っておく
# (使われ方は `.args` の参照のみなので getfullargspec で代替できる)。
if not hasattr(inspect, "getargspec"):
    inspect.getargspec = inspect.getfullargspec

import numpy as np

from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa  # noqa: F401  (ParaView 側の初期化用)
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkDataObject
from vtkmodules.vtkFiltersGeneral import vtkGradientFilter

POINTS = vtkDataObject.FIELD_ASSOCIATION_POINTS
CELLS = vtkDataObject.FIELD_ASSOCIATION_CELLS

VEL_GRAD_NAMES = (
    "dUxdx", "dUxdy", "dUxdz",
    "dUydx", "dUydy", "dUydz",
    "dUzdx", "dUzdy", "dUzdz",
)
RO_GRAD_NAMES = ("drodx", "drody", "drodz")


# ---------------------------------------------------------------- 補助関数

def _attr(ds, assoc):
    return ds.GetPointData() if assoc == POINTS else ds.GetCellData()


def _has(ds, assoc, name):
    return _attr(ds, assoc).GetArray(name) is not None


def _get(ds, assoc, name):
    arr = _attr(ds, assoc).GetArray(name)
    if arr is None:
        return None
    return vtk_to_numpy(arr).astype(np.float64)


def _add(ds, assoc, name, values):
    values = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    varr = numpy_to_vtk(values, deep=1)
    varr.SetName(name)
    _attr(ds, assoc).AddArray(varr)


def _detect_assoc(ds):
    """ro / Ux が乗っている側を自動判定する (cell モード=Cell, node モード=Point)。"""
    for assoc in (POINTS, CELLS):
        if _has(ds, assoc, "ro") or _has(ds, assoc, "Ux"):
            return assoc
    # 判定不能ならデータ数の多い方
    return POINTS if ds.GetPointData().GetNumberOfArrays() >= ds.GetCellData().GetNumberOfArrays() else CELLS


def _velocity(ds, assoc):
    ux, uy, uz = (_get(ds, assoc, n) for n in ("Ux", "Uy", "Uz"))
    if ux is None:
        vec = _get(ds, assoc, "U")
        if vec is not None and vec.ndim == 2 and vec.shape[1] == 3:
            return vec
        return None
    n = ux.shape[0]
    uy = np.zeros(n) if uy is None else uy
    uz = np.zeros(n) if uz is None else uz
    return np.column_stack((ux, uy, uz))


def _mesh_gradient(ds, assoc, values, name="__forge_tmp"):
    """vtkGradientFilter でメッシュから勾配を計算する。

    スカラー入力なら (N,3)、3 成分ベクトル入力なら (N,9) を返す。
    (N,9) の並びは行優先 = (du/dx, du/dy, du/dz, dv/dx, …) で
    forge の dUxdx… と同じ規約。
    """
    tmp = ds.NewInstance()
    tmp.ShallowCopy(ds)
    _add(tmp, assoc, name, values)

    gf = vtkGradientFilter()
    gf.SetInputData(tmp)
    gf.SetInputArrayToProcess(0, 0, 0, assoc, name)
    gf.SetResultArrayName("__forge_grad")
    gf.Update()

    out = gf.GetOutput()
    grad = vtk_to_numpy(_attr(out, assoc).GetArray("__forge_grad")).astype(np.float64)
    return grad


def _velocity_gradient(ds, assoc, vel, use_solver):
    """速度勾配テンソル g[:, i, j] = dU_i/dx_j を返す。"""
    if use_solver and all(_has(ds, assoc, n) for n in VEL_GRAD_NAMES):
        cols = [_get(ds, assoc, n) for n in VEL_GRAD_NAMES]
        return np.column_stack(cols).reshape(-1, 3, 3), "solver"
    return _mesh_gradient(ds, assoc, vel).reshape(-1, 3, 3), "mesh"


def _density_gradient(ds, assoc, use_solver):
    if use_solver and all(_has(ds, assoc, n) for n in RO_GRAD_NAMES):
        return np.column_stack([_get(ds, assoc, n) for n in RO_GRAD_NAMES]), "solver"
    ro = _get(ds, assoc, "ro")
    if ro is None:
        return None, None
    return _mesh_gradient(ds, assoc, ro), "mesh"


def _leaf_pairs(inp, out):
    """(入力 leaf, 出力 leaf) の組を返す。composite / 単一データセット両対応。"""
    if inp.IsA("vtkCompositeDataSet"):
        out.CopyStructure(inp)
        pairs = []
        it = inp.NewIterator()
        it.InitTraversal()
        while not it.IsDoneWithTraversal():
            blk = it.GetCurrentDataObject()
            if blk is not None:
                new_blk = blk.NewInstance()
                new_blk.ShallowCopy(blk)
                out.SetDataSet(it, new_blk)
                pairs.append((blk, new_blk))
            it.GoToNextItem()
        return pairs
    out.ShallowCopy(inp)
    return [(inp, out)]


# ---------------------------------------------------------------- フィルタ本体

@smproxy.filter(label="Forge Derived Quantities")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet", "vtkCompositeDataSet"], composite_data_supported=True)
class ForgeDerivedQuantities(VTKPythonAlgorithmBase):
    """forge の結果からマッハ数・シュリーレン・Q 値・ヘリシティを計算する。"""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=1,
            inputType="vtkDataObject", outputType="vtkDataObject")
        self._gamma = 1.4
        self._use_solver_grad = True
        self._use_solver_sonic = True
        self._mach = True
        self._schlieren = True
        self._qcrit = True
        self._helicity = True
        self._lambda2 = False
        self._sch_dir = [1.0, 0.0, 0.0]
        self._sch_exp = 15.0

    # -- properties -------------------------------------------------

    @smproperty.doublevector(name="Gamma", default_values=1.4)
    @smdomain.doublerange(min=1.0, max=2.0)
    def SetGamma(self, value):
        """比熱比 (音速 a=sqrt(gamma*P/ro) 用)。"""
        self._gamma = float(value)
        self.Modified()

    @smproperty.intvector(name="UseSolverGradients", label="Use Solver Gradients", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetUseSolverGradients(self, value):
        """ON: ソルバが出力した勾配 (dUxdx…, drodx…) をそのまま使う。
        OFF (または配列が無い場合): vtkGradientFilter でメッシュから勾配を計算する。"""
        self._use_solver_grad = bool(value)
        self.Modified()

    @smproperty.intvector(name="UseSolverSonic", label="Use Solver Sound Speed", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetUseSolverSonic(self, value):
        """ON: 配列 sonic をそのまま音速に使う。OFF/無い場合は a=sqrt(gamma*P/ro)。"""
        self._use_solver_sonic = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeMach", label="Compute Mach", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeMach(self, value):
        self._mach = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeSchlieren", label="Compute Schlieren", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeSchlieren(self, value):
        self._schlieren = bool(value)
        self.Modified()

    @smproperty.doublevector(name="SchlierenDirection", label="Schlieren Direction",
                             default_values=[1.0, 0.0, 0.0], number_of_elements=3)
    def SetSchlierenDirection(self, x, y, z):
        """方向シュリーレン schlieren_dir = (n^*grad ro)/max|n^*grad ro| の方向 n。"""
        self._sch_dir = [float(x), float(y), float(z)]
        self.Modified()

    @smproperty.doublevector(name="SchlierenExponent", label="Schlieren Exponent", default_values=15.0)
    @smdomain.doublerange(min=0.0, max=200.0)
    def SetSchlierenExponent(self, value):
        """schlieren = exp(-k*|grad ro|/max|grad ro|) の k。大きいほど弱い勾配を強調。"""
        self._sch_exp = float(value)
        self.Modified()

    @smproperty.intvector(name="ComputeQ", label="Compute Q", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeQ(self, value):
        self._qcrit = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeLambda2", label="Compute Lambda2", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeLambda2(self, value):
        self._lambda2 = bool(value)
        self.Modified()

    @smproperty.intvector(name="ComputeHelicity", label="Compute Helicity", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeHelicity(self, value):
        self._helicity = bool(value)
        self.Modified()

    # -- pipeline ---------------------------------------------------

    def RequestDataObject(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0])
        if inp is None:
            return 0
        out = vtkDataObject.GetData(outInfo)
        if out is None or not out.IsA(inp.GetClassName()):
            out = inp.NewInstance()
            outInfo.GetInformationObject(0).Set(vtkDataObject.DATA_OBJECT(), out)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0])
        out = vtkDataObject.GetData(outInfo)

        pairs = _leaf_pairs(inp, out)

        # 1 パス目: 各 leaf の派生量を計算 (シュリーレン正規化用に全体最大も集める)
        results = []
        gmax_mag = 0.0
        gmax_dir = 0.0
        for _, leaf in pairs:
            if not leaf.IsA("vtkDataSet"):
                results.append(None)
                continue
            res = self._compute_leaf(leaf)
            results.append(res)
            if res is not None and "grad_ro_mag" in res:
                gmax_mag = max(gmax_mag, float(res["grad_ro_mag"].max(initial=0.0)))
                gmax_dir = max(gmax_dir, float(np.abs(res["grad_ro_dir"]).max(initial=0.0)))

        # 2 パス目: 正規化してから書き出す
        for (_, leaf), res in zip(pairs, results):
            if res is None:
                continue
            assoc = res["assoc"]
            for name, values in res["arrays"].items():
                _add(leaf, assoc, name, values)
            if "grad_ro_mag" in res:
                mag = res["grad_ro_mag"]
                dir_ = res["grad_ro_dir"]
                norm_mag = mag / gmax_mag if gmax_mag > 0.0 else np.zeros_like(mag)
                norm_dir = dir_ / gmax_dir if gmax_dir > 0.0 else np.zeros_like(dir_)
                _add(leaf, assoc, "schlieren_mag", norm_mag)
                _add(leaf, assoc, "schlieren_dir", norm_dir)
                _add(leaf, assoc, "schlieren", np.exp(-self._sch_exp * norm_mag))
        return 1

    # -- 中身 --------------------------------------------------------

    def _compute_leaf(self, ds):
        n = ds.GetNumberOfPoints()
        if n == 0:
            return None
        assoc = _detect_assoc(ds)
        if assoc == CELLS and ds.GetNumberOfCells() == 0:
            return None

        res = {"assoc": assoc, "arrays": {}}
        arrays = res["arrays"]

        vel = _velocity(ds, assoc)
        if vel is not None:
            arrays["U"] = vel

        ro = _get(ds, assoc, "ro")

        # --- マッハ数
        if self._mach and vel is not None:
            sonic = _get(ds, assoc, "sonic") if self._use_solver_sonic else None
            if sonic is None:
                pres = _get(ds, assoc, "P")
                if pres is not None and ro is not None:
                    sonic = np.sqrt(np.maximum(self._gamma * pres / np.maximum(ro, 1e-30), 0.0))
            if sonic is not None:
                umag = np.linalg.norm(vel, axis=1)
                arrays["sound_speed"] = sonic
                arrays["Mach"] = umag / np.maximum(sonic, 1e-30)

        # --- シュリーレン (正規化は RequestData 側の全体最大で行う)
        if self._schlieren and ro is not None:
            grad_ro, _ = _density_gradient(ds, assoc, self._use_solver_grad)
            if grad_ro is not None:
                arrays["grad_ro"] = grad_ro
                d = np.asarray(self._sch_dir, dtype=np.float64)
                dnorm = np.linalg.norm(d)
                d = d / dnorm if dnorm > 0.0 else np.array([1.0, 0.0, 0.0])
                res["grad_ro_mag"] = np.linalg.norm(grad_ro, axis=1)
                res["grad_ro_dir"] = grad_ro @ d

        # --- 速度勾配ベースの量 (Q / lambda2 / ヘリシティ)
        need_grad = (self._qcrit or self._helicity or self._lambda2) and vel is not None
        if need_grad:
            g, _src = _velocity_gradient(ds, assoc, vel, self._use_solver_grad)

            # 渦度 omega = rot U
            omega = np.column_stack((
                g[:, 2, 1] - g[:, 1, 2],
                g[:, 0, 2] - g[:, 2, 0],
                g[:, 1, 0] - g[:, 0, 1],
            ))
            omega_mag = np.linalg.norm(omega, axis=1)
            arrays["vorticity"] = omega
            arrays["vorticity_mag"] = omega_mag

            if self._qcrit or self._lambda2:
                gT = np.transpose(g, (0, 2, 1))
                S = 0.5 * (g + gT)
                W = 0.5 * (g - gT)
                s2 = np.einsum("nij,nij->n", S, S)
                w2 = np.einsum("nij,nij->n", W, W)
                if self._qcrit:
                    arrays["Q"] = 0.5 * (w2 - s2)
                    denom = 0.5 * (w2 + s2)
                    arrays["Q_norm"] = np.where(denom > 0.0, 0.5 * (w2 - s2) / np.maximum(denom, 1e-30), 0.0)
                if self._lambda2:
                    M = np.matmul(S, S) + np.matmul(W, W)
                    M = 0.5 * (M + np.transpose(M, (0, 2, 1)))  # 数値的な非対称を除去
                    # 発散した run では場に NaN が混ざる。固有値計算が落ちないよう
                    # 有限な要素だけ解き、残りは NaN のままにする。
                    ok = np.isfinite(M).all(axis=(1, 2))
                    lam2 = np.full(M.shape[0], np.nan)
                    if ok.any():
                        lam2[ok] = np.linalg.eigvalsh(M[ok])[:, 1]  # 昇順の中間固有値
                    arrays["lambda2"] = lam2

            if self._helicity:
                hel = np.einsum("ni,ni->n", vel, omega)
                arrays["helicity"] = hel
                denom = np.linalg.norm(vel, axis=1) * omega_mag
                arrays["helicity_norm"] = np.where(denom > 0.0, hel / np.maximum(denom, 1e-30), 0.0)

        return res


# ---------------------------------------------------------------- 飽和量 (凝縮後処理)
#
# forge の condensationProperties_d.cuh / condensationSource_d.cuh と同じ式を numpy で写したもの。
# 凝縮 OFF の run (condS_0 / condTsat_0 が無い) でも、T / ro / Y から過飽和・過冷却を評価できる。

H2O_RV = 461.5          # J/(kg K)  (forge condProps_H2O().R)
N2_RV = 296.8           # J/(kg K)  (forge condProps_N2().R)
N2_TC = 126.192
N2_PSAT_TSWITCH = 50.0  # COND_PSAT_TSWITCH
N2_LATENT_TA = 70.0     # COND_N2_LATENT_TA
N2_CPV = 1038.8         # COND_N2_CPV
COND_T_PROP_FLOOR = 45.0


def _h2o_psat_liquid(T):
    """Murphy & Koop (2005) 過冷却液の飽和蒸気圧 [Pa] (forge h2o_psat と同一, 120 K 下限クランプ)。"""
    Tc = np.maximum(T, 120.0)
    lnp = (54.842763 - 6763.22 / Tc - 4.210 * np.log(Tc) + 0.000367 * Tc
           + np.tanh(0.0415 * (Tc - 218.8)) * (53.878 - 1331.22 / Tc - 9.44523 * np.log(Tc) + 0.014025 * Tc))
    return np.exp(lnp)


def _h2o_psat_ice(T):
    """Murphy & Koop (2005) 氷の飽和蒸気圧 [Pa] (T > 110 K)。"""
    Tc = np.maximum(T, 110.0)
    return np.exp(9.550426 - 5723.265 / Tc + 3.53068 * np.log(Tc) - 0.00728332 * Tc)


def _n2_latent_poly(T):
    Tcl = np.clip(T, COND_T_PROP_FLOOR, N2_TC - 0.5)
    L = (-2.137e-8 * Tcl ** 4 + 7.18e-6 * Tcl ** 3 - 9.142e-4 * Tcl ** 2 + 0.05069 * Tcl - 0.809) * 1.0e6
    return np.maximum(L, 0.0)


def _n2_psat_jacobsen(Tcl):
    n1, n2, n3 = 8394.409444, -1890.045259, -7.282229165
    n4, n5, n6 = 0.01022850966, 5.556063825e-4, -5.944544662e-6
    n7, n8, n9 = 2.715433932e-8, -4.879535904e-11, 509.5360824
    dTc = N2_TC - Tcl
    lnP = (n1 / Tcl + n2 + n3 * Tcl + n4 * np.power(np.maximum(dTc, 0.0), 1.95)
           + n5 * Tcl ** 3 + n6 * Tcl ** 4 + n7 * Tcl ** 5 + n8 * Tcl ** 6 + n9 * np.log(Tcl))
    return np.exp(lnP) * 101325.0


def _n2_psat(T, psat_lowT=1, latent_lowT=1, cl=2000.0):
    """N2 過冷却液の飽和蒸気圧 [Pa] (forge n2_psat_ex と同一)。50 K 未満は Clausius–Clapeyron 外挿。"""
    T = np.asarray(T, dtype=np.float64)
    Tcl = np.minimum(np.maximum(T, N2_PSAT_TSWITCH), N2_TC - 0.5)
    hi = _n2_psat_jacobsen(Tcl)
    Tsw = N2_PSAT_TSWITCH
    psw = _n2_psat_jacobsen(np.array(Tsw))
    Tlo = np.maximum(T, 5.0)
    if psat_lowT and latent_lowT:
        Ta = N2_LATENT_TA; La = float(_n2_latent_poly(np.array(Ta))); cpr = N2_CPV - cl
        lnr = ((La - cpr * Ta) * (1.0 / Tsw - 1.0 / Tlo) + cpr * np.log(Tlo / Tsw)) / N2_RV
        lo = psw * np.exp(lnr)
    else:
        Lref = float(_n2_latent_poly(np.array(Tsw)))
        lo = psw * np.exp(-(Lref / N2_RV) * (1.0 / Tlo - 1.0 / Tsw))
    return np.where(T >= Tsw, hi, lo)


def _tsat_newton(psat_fn, pv, T_guess):
    """p_sat(T_sat) = p_v を Newton で解く (forge cond_Tsat と同じ初期値・クリップ・反復数)。
    導関数は同じ p_sat の数値微分 d ln p_sat/dT (forge は L/(R T²) = Clausius–Clapeyron; 根は同じ)。"""
    pv = np.asarray(pv, dtype=np.float64)
    valid = pv > 1.0e-6
    lnpv = np.log(np.where(valid, pv, 1.0))
    T = np.where((T_guess > 50.0) & (T_guess < 1000.0), T_guess, 250.0).astype(np.float64)
    h = 0.01
    for _ in range(25):
        ps = psat_fn(T)
        f = np.log(np.maximum(ps, 1.0e-300)) - lnpv
        dfdT = (np.log(np.maximum(psat_fn(T + h), 1.0e-300)) - np.log(np.maximum(psat_fn(T - h), 1.0e-300))) / (2 * h)
        dfdT = np.maximum(dfdT, 1.0e-6)
        dT = np.clip(f / dfdT, -0.3 * T, 0.3 * T)
        T = np.clip(T - dT, 50.0, 1000.0)
    return np.where(valid, T, 0.0)


@smproxy.filter(label="Forge Saturation")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet", "vtkCompositeDataSet"], composite_data_supported=True)
class ForgeSaturation(VTKPythonAlgorithmBase):
    """凝縮種の蒸気分圧・飽和蒸気圧・飽和温度・過冷却度・過飽和度を計算する。"""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=1,
            inputType="vtkDataObject", outputType="vtkDataObject")
        self._species = 0
        self._yv_array = ""          # 既定は空: Run Config か明示指定が必須 (Y1 を自動採用しない)
        self._run_config = ""
        self._g_array = "g_0"
        self._yv_const = 0.0
        self._ice = False
        self._n2_psat_lowT = 1
        self._n2_latent_lowT = 1
        self._n2_liquid_cp = 2000.0

    # -- properties -------------------------------------------------

    @smproperty.intvector(name="Species", default_values=0)
    @smdomain.xml('<EnumerationDomain name="enum">'
                  '<Entry text="H2O" value="0"/>'
                  '<Entry text="N2" value="1"/>'
                  '</EnumerationDomain>')
    def SetSpecies(self, value):
        """凝縮種。H2O = Murphy & Koop 過冷却液 (R_v 461.5)、N2 = Jacobsen + C–C 外挿 (R_v 296.8)。"""
        self._species = int(value)
        self.Modified()

    @smproperty.stringvector(name="RunConfig", label="Run Config (solverConfig.yaml path)", default_values="")
    def SetRunConfig(self, value):
        """run の solverConfig.yaml のパス。指定すると凝縮種の配列 Y{index} を physProp.species / condensation から
        名前で解決する (Vapor Mass Fraction Array より優先)。"""
        self._run_config = str(value).strip()
        self.Modified()

    @smproperty.stringvector(name="VaporMassFractionArray", label="Vapor Mass Fraction Array", default_values="")
    def SetVaporMassFractionArray(self, value):
        """蒸気 (凝縮種) の質量分率配列名 (例 Y1)。Run Config 未指定なら必須 (空はエラー; Y1 を自動採用しない)。
        "none" で配列を使わず定数 (Vapor Mass Fraction Constant) を使う。"""
        self._yv_array = str(value).strip()
        self.Modified()

    @smproperty.stringvector(name="LiquidMassFractionArray", label="Liquid Mass Fraction Array", default_values="g_0")
    def SetLiquidMassFractionArray(self, value):
        """液相質量分率配列名 (凝縮 ON の run の g_0)。無ければ 0 として扱う。"""
        self._g_array = str(value).strip()
        self.Modified()

    @smproperty.doublevector(name="VaporMassFractionConstant", label="Vapor Mass Fraction Constant", default_values=0.0)
    @smdomain.doublerange(min=0.0, max=1.0)
    def SetVaporMassFractionConstant(self, value):
        """蒸気配列が無いときの一様質量分率 (空気凝縮 CPG carrier は condVaporMassFraction = 0.7671)。
        0 なら純蒸気 (p_v = P) として扱う。"""
        self._yv_const = float(value)
        self.Modified()

    @smproperty.intvector(name="ComputeIce", label="Compute Ice (H2O)", default_values=0)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetComputeIce(self, value):
        """H2O で氷基準の p_sat_ice / S_ice / T_sat_ice / subcooling_ice も出す。"""
        self._ice = bool(value)
        self.Modified()

    @smproperty.intvector(name="N2PsatLowT", label="N2 psat low-T consistent", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetN2PsatLowT(self, value):
        """N2: 50 K 未満の外挿を低温整合潜熱で行う (forge condN2PsatLowT)。"""
        self._n2_psat_lowT = int(bool(value))
        self.Modified()

    @smproperty.intvector(name="N2LatentLowT", label="N2 latent low-T consistent", default_values=1)
    @smdomain.xml('<BooleanDomain name="bool"/>')
    def SetN2LatentLowT(self, value):
        """N2: 潜熱の低温整合 (forge condN2LatentLowT)。"""
        self._n2_latent_lowT = int(bool(value))
        self.Modified()

    @smproperty.doublevector(name="N2LiquidCp", label="N2 liquid cp [J/kg/K]", default_values=2000.0)
    def SetN2LiquidCp(self, value):
        """N2 液比熱 (forge condN2LiquidCp)。"""
        self._n2_liquid_cp = float(value)
        self.Modified()

    # -- pipeline ---------------------------------------------------

    def RequestDataObject(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0])
        if inp is None:
            return 0
        out = vtkDataObject.GetData(outInfo)
        if out is None or not out.IsA(inp.GetClassName()):
            out = inp.NewInstance()
            outInfo.GetInformationObject(0).Set(vtkDataObject.DATA_OBJECT(), out)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        inp = vtkDataObject.GetData(inInfo[0])
        out = vtkDataObject.GetData(outInfo)
        for _, leaf in _leaf_pairs(inp, out):
            if leaf.IsA("vtkDataSet"):
                self._compute_leaf(leaf)
        return 1

    def _resolve_vapor_array(self):
        """蒸気配列名を決める: Run Config があれば名前解決、無ければ明示配列 (空はエラー)。"none" は配列なし。"""
        if self._run_config:
            import os, sys
            run_dir = self._run_config
            if os.path.isfile(run_dir):
                run_dir = os.path.dirname(os.path.abspath(run_dir))
            tools_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if tools_dir not in sys.path:
                sys.path.insert(0, tools_dir)
            from forge_species import species_info
            info = species_info(run_dir)
            if info["vapor_array"] is not None:
                return info["vapor_array"]
            # 凝縮 OFF の run: H2O があればそれ、無ければ明示配列へ
            for cand in ("H2O", "WATER"):
                if cand in info["index"]:
                    return f"Y{info['index'][cand]}"
            if self._yv_array:
                return None if self._yv_array.lower() == "none" else self._yv_array
            raise RuntimeError(f"Forge Saturation: {run_dir} に凝縮種/H2O が無い。Vapor Mass Fraction Array を明示すること")
        if not self._yv_array:
            raise RuntimeError("Forge Saturation: Vapor Mass Fraction Array が空。Run Config (solverConfig.yaml) を指定するか、"
                               "配列名 (例 Y1) を明示する (Y1 を既定採用しない: 種順序で H2O の index は変わる)。"
                               " 配列を使わない (純蒸気/定数) なら \"none\"")
        return None if self._yv_array.lower() == "none" else self._yv_array

    def _psat_fn(self):
        if self._species == 1:
            return lambda T: _n2_psat(T, self._n2_psat_lowT, self._n2_latent_lowT, self._n2_liquid_cp)
        return _h2o_psat_liquid

    def _compute_leaf(self, ds):
        if ds.GetNumberOfPoints() == 0:
            return
        assoc = _detect_assoc(ds)
        if not _has(ds, assoc, "T"):
            assoc = POINTS if _has(ds, POINTS, "T") else CELLS
        T = _get(ds, assoc, "T")
        if T is None:
            return
        ro = _get(ds, assoc, "ro")
        P = _get(ds, assoc, "P")
        Rv = N2_RV if self._species == 1 else H2O_RV

        g = _get(ds, assoc, self._g_array) if self._g_array else None
        g = np.zeros_like(T) if g is None else np.maximum(g, 0.0)
        yv_name = self._resolve_vapor_array()
        yv = _get(ds, assoc, yv_name) if yv_name else None
        if yv_name and yv is None:
            raise RuntimeError(f"Forge Saturation: 蒸気配列 '{yv_name}' が入力に無い")
        if yv is None and self._yv_const > 0.0:
            yv = np.full_like(T, self._yv_const)

        # 蒸気分圧 (forge cond_vapor_state と同一)
        if yv is not None and ro is not None:
            pv = ro * np.maximum(yv - g, 0.0) * Rv * T
        elif P is not None:
            pv = P.copy()   # 純蒸気
        else:
            return

        psat = self._psat_fn()
        ps = psat(T)
        S = pv / np.maximum(ps, 1.0e-300)
        Tsat = _tsat_newton(psat, pv, T)
        _add(ds, assoc, "p_vapor", pv)
        _add(ds, assoc, "p_sat", ps)
        _add(ds, assoc, "supersaturation", S)
        _add(ds, assoc, "log10_S", np.log10(np.maximum(S, 1.0e-300)))
        _add(ds, assoc, "T_sat", Tsat)
        _add(ds, assoc, "subcooling", np.where(Tsat > 0.0, Tsat - T, 0.0))
        if self._ice and self._species == 0:
            psi = _h2o_psat_ice(T)
            Tsi = _tsat_newton(_h2o_psat_ice, pv, T)
            _add(ds, assoc, "p_sat_ice", psi)
            _add(ds, assoc, "S_ice", pv / np.maximum(psi, 1.0e-300))
            _add(ds, assoc, "T_sat_ice", Tsi)
            _add(ds, assoc, "subcooling_ice", np.where(Tsi > 0.0, Tsi - T, 0.0))
