"""runner.run_converter の試験 (2026-10-13、plan tooling-nozzle-isothermal-wall-chain §5.1 #31 の準備が変換器の終了時の既知の
GPU エラーで止まった件)。偽の変換器 (bash) で次を確かめる:
(a) rc 0 はそのまま通る。
(b) rc 1 + stdout `Write Input HDF5` + stderr `GPUassert: invalid argument` + 出力が MESH/COORD・VALUE を持つ → 通し、記録を書く。
(c) 同じ症状でも出力が無い・壊れている・VALUE が無い → CalledProcessError。
(d) 別のエラー (stderr が違う、stdout に書き出しの行が無い) → CalledProcessError。
使い方: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_converter_exit_tests.py"""
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.evaluate.runner import run_converter  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(("ok   " if ok else "FAIL ") + name + (f"  [{info}]" if info else ""))
    if not ok:
        FAIL.append(name)


def fake(tmp: Path, rc: int, stdout: str, stderr: str, write: str) -> Path:
    """偽の変換器: write = 'good' (MESH/COORD と VALUE)、'novalue' (VALUE 無し)、'junk' (HDF5 でない)、'none' (書かない)。"""
    py = tmp / "fake_write.py"
    py.write_text(
        "import sys, h5py, numpy as np\n"
        "kind, out = sys.argv[1], sys.argv[2]\n"
        "if kind == 'junk':\n    open(out, 'w').write('not hdf5')\n"
        "elif kind in ('good', 'novalue'):\n"
        "    with h5py.File(out, 'w') as f:\n"
        "        f['MESH/COORD'] = np.zeros(6)\n"
        "        if kind == 'good': f['VALUE/ro'] = np.ones(2)\n")
    sh = tmp / f"conv_{rc}_{write}.sh"
    sh.write_text("#!/bin/bash\n"
                  f"[ {write!r} != none ] && {sys.executable} {py} {write} \"$2\"\n"
                  f"printf '%s\\n' {stdout!r}\n"
                  f"printf '%s\\n' {stderr!r} >&2\n"
                  f"exit {rc}\n")
    sh.chmod(0o755)
    return sh


def run(tmp, rc, so, se, write):
    d = Path(tempfile.mkdtemp(dir=tmp))
    conv = fake(d, rc, so, se, write)
    try:
        r = run_converter([conv, "nozzle.msh", "nozzle.h5"], cwd=d)
        return "ok", r, d
    except subprocess.CalledProcessError as e:
        return "raise", e, d


def main():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        W, G = "*** Write Input HDF5 ***\nbaseName=nozzle", "GPUassert: invalid argument cudaWrapper.cu 81"
        st, r, d = run(tmp, 0, "done", "", "good")
        check("(a) rc 0 は通る", st == "ok" and r.returncode == 0 and not list(d.glob("*_converter_exit.log")))
        st, r, d = run(tmp, 1, W, G, "good")
        check("(b) 既知の終了時のエラー + 書き終えた HDF5 は通す", st == "ok" and r.returncode == 1)
        check("(b) 記録を書く", (d / "nozzle_converter_exit.log").is_file())
        for kind in ("none", "junk", "novalue"):
            st, e, _ = run(tmp, 1, W, G, kind)
            check(f"(c) 既知の症状でも出力が {kind} なら止める", st == "raise")
        st, e, _ = run(tmp, 1, W, "Segmentation fault", "good")
        check("(d) stderr が違うエラーは止める", st == "raise")
        st, e, _ = run(tmp, 1, "*** Set Initial Values ***", G, "good")
        check("(d) 書き出しの行の無い GPU エラーは止める", st == "raise")
        st, e, _ = run(tmp, 2, W, G, "good")
        check("(b') rc 2 でも同じ症状なら通す (rc の値では区別しない)", st == "ok")
    print(f"FAIL 件数: {len(FAIL)}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
