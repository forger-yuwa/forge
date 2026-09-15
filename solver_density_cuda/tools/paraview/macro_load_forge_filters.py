"""ParaView マクロ: forge_filters.py プラグインをワンクリックで読み込む。

使い方:
  ParaView GUI で Macros > Add new macro... から本ファイルを登録する。
  以後ツールバーの `macro_load_forge_filters` を押すと
  Filters > Alphabetical > "Forge Derived Quantities" と "Forge Saturation" が使えるようになる。

Ubuntu の ParaView 5.11 + Python 3.12 では paraview.detail.pythonalgorithm が
Python 3.11 で削除された inspect.getargspec を import するため、
Tools > Manage Plugins からの直接ロードが失敗する。本マクロは
その shim を入れてからプラグインを読むので、その環境でも動く。

プラグインの場所は次の順に探す。別の場所に置いた場合は環境変数
FORGE_PARAVIEW_FILTERS に絶対パスを設定する。
"""

import inspect
import os

if not hasattr(inspect, "getargspec"):
    inspect.getargspec = inspect.getfullargspec

# shim を効かせた状態で先に import してキャッシュさせる (以降のプラグイン読込で再利用される)
import paraview.detail.pythonalgorithm  # noqa: F401,E402

from paraview.simple import LoadPlugin  # noqa: E402

_REL = os.path.join("solver_density_cuda", "tools", "paraview", "forge_filters.py")

_CANDIDATES = [
    os.environ.get("FORGE_PARAVIEW_FILTERS", ""),
    # マクロを repo から直接開いた場合 (GUI 登録時は ~/.config/ParaView/Macros に複製される)
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "forge_filters.py"),
    os.path.join("/workspace", _REL),                                  # run_paraview_gui.sh (Docker)
    os.path.join(os.path.expanduser("~"), "work", "forge", _REL),      # WSL native
]

for _path in _CANDIDATES:
    if _path and os.path.isfile(_path):
        LoadPlugin(_path, remote=False, ns=globals())
        print("[forge] loaded ParaView plugin: %s" % _path)
        break
else:
    raise RuntimeError(
        "forge_filters.py が見つかりません。環境変数 FORGE_PARAVIEW_FILTERS に絶対パスを設定してください。"
        " 探した場所: %s" % [p for p in _CANDIDATES if p])
