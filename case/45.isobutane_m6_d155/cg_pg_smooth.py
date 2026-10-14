"""plan tooling-nozzle-core-grid §4.94 の材料 (探索的・登録外): Q0・Q0+δ★/2・Q0+δ★ の原始量 (各 run の res_1 = step の始めの状態) の比例の誤差と、壁法線の速度の符号の反転。usage: python3 cg_pg_smooth.py (case/45 で)"""
import json, h5py, numpy as np
# 状態 Q0・Q0+δ★/2・Q0+δ★ の原始量 (各 run の res_1 = その step の始めの状態) で、場の変化が比例しているかを見る (GI と GB)。
BASE = "run_0693_pdt_off"; info = json.load(open(BASE + "/prepare_info.json")); ni = int(info["mesh"]["ni"]); S = float(info["scale_m"])
def rd(run, keys):
    with h5py.File(run + "/res_1.h5", "r") as h: return {k: np.asarray(h["VALUE/" + k][:], float) for k in keys}
K = ("Ux", "Uy", "P", "T", "ro", "vis_turb", "volume")
with h5py.File(BASE + "/res_20000.h5", "r") as h: xy = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)[:, :2]
x = xy[:, 0] / S; nj = x.size // ni; jw = nj - 1 - np.arange(x.size) % nj
# 壁法線: 各列で jw = 0 の節点から jw = 1 の節点への向き
X2 = xy.reshape(ni, nj, 2)[:, ::-1, :]; nvec = X2[:, 1, :] - X2[:, 0, :]; nvec /= np.linalg.norm(nvec, axis=1, keepdims=True)
nn = np.repeat(nvec[:, None, :], nj, axis=1)[:, ::-1, :].reshape(-1, 2)
for arm, (r0, rh, r1) in {"GI": ("run_0696_pg_gi5", "run_0702_pg_gih", "run_0704_pg_gib"), "GB": ("run_0699_pg_gb5", "run_0703_pg_gbh", "run_0705_pg_gbb")}.items():
    A, H, B = rd(r0, K), rd(rh, K), rd(r1, K); V = A["volume"]
    m = (x >= -8) & (x <= -2) & (jw >= 1) & (jw <= 10)
    print(arm, "縮流部の壁際 (壁から 1〜10) の比例の誤差 ‖Δ(½) − ½Δ(1)‖/‖½Δ(1)‖:")
    for k in ("Ux", "Uy", "P", "T", "ro", "vis_turb"):
        d1 = B[k] - A[k]; dh = H[k] - A[k]; num = np.sqrt(np.sum(V[m] * (dh - 0.5 * d1)[m] ** 2)); den = np.sqrt(np.sum(V[m] * (0.5 * d1)[m] ** 2))
        print("   %-8s %.3e   (‖½Δ(1)‖ %.3e)" % (k, num / max(den, 1e-300), den / np.sqrt(np.sum(V[m]))))
    un0 = A["Ux"] * nn[:, 0] + A["Uy"] * nn[:, 1]; un1 = B["Ux"] * nn[:, 0] + B["Uy"] * nn[:, 1]; unh = H["Ux"] * nn[:, 0] + H["Uy"] * nn[:, 1]
    flip = (np.sign(un0) != np.sign(un1)) | (np.sign(un0) != np.sign(unh))
    mm = (x >= -8) & (x <= -2) & (jw <= 10)
    print("   壁法線の速度の符号が変わる節点: 縮流部の壁から 0〜10 で %d / %d (jw の内訳 %s)" % (int(np.sum(flip & mm)), int(np.sum(mm)), np.bincount(jw[flip & mm], minlength=11)[:11].tolist()))
    print("   壁法線の速度の大きさ (縮流部の壁から 1〜10 の中央値): %.3e m/s、|Δ(1)| の中央値 %.3e" % (np.median(np.abs(un0[m])), np.median(np.abs((un1 - un0)[m]))))
