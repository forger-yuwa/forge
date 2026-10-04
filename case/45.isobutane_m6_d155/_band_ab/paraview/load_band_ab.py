"""ParaView GUI 用: 帯修正 A/B の結果を読み込み、Turbo で色付けする。
起動: ~/opt/ParaView-6.1.1-MPI-Linux-Python3.12-x86_64/bin/paraview --script=load_band_ab.py   (このディレクトリで)
"""
import os
from paraview.simple import *  # noqa: F401,F403

here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
specs = [("B1", "M_vs_Euler_pct", 0.5, 0.0), ("B0", "M_vs_Euler_pct", 0.5, -12.0), ("A0p", "M_vs_Euler_pct", 0.5, -24.0),
         ("B1_minus_B0", "dM_pct", 0.5, -36.0)]
view = GetActiveViewOrCreate("RenderView")
for name, arr, lim, dy in specs:
    r = XDMFReader(registrationName=name, FileNames=[os.path.join(here, name + ".xmf")])
    t = Transform(registrationName=name + "_shift", Input=r); t.Transform.Translate = [0.0, dy, 0.0]
    d = Show(t, view); ColorBy(d, ("POINTS", arr))
    lut = GetColorTransferFunction(arr); lut.ApplyPreset("Turbo", True); lut.RescaleTransferFunction(-lim, lim)
    d.SetScalarBarVisibility(view, True)
    Hide(r, view)
    txt = Text(registrationName=name + "_label"); txt.Text = name
for n in ("Euler",):
    XDMFReader(registrationName=n, FileNames=[os.path.join(here, n + ".xmf")])
view.InteractionMode = "2D"; ResetCamera(view); Render(view)
