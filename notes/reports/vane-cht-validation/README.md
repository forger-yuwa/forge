# Cooled Vane CHT Validation — 報告の原稿

Artifact「Cooled Vane CHT Validation」(https://claude.ai/artifact/LqUXeaKXTSvNbAd2Sz46bw) の**原稿の控え**。
2026-09-21 に `/tmp` の作業領域が再起動で消えて原稿を Artifact から復元する羽目になったので、公開のたびにここへ写す。

- `index.html` … 公開した本文 (Artifact の skeleton に包まれる前のもの)。V36 = 2026-09-22。
- `img/` … 図。**リポジトリの `.gitignore` が `*.png` を除外しているので追跡されない** (手元の控えのみ。消えたら下の作図元から作り直す)。作図元: 線グラフ `case/53.c3x_vane_cht/tools/report_figures.py`、§04 の 2 枚 `fx_ab_figure.py` / `resid_split.py`、
  コンタ `solver_density_cuda/tools/plot_field_contours.py`、形状 `check_geometry.py`、連成壁温 `plot_tw.py`、固体 `plot_solid.py`。
  `domain.png` / `mesh_forcing.png` / `shed.png` は作図スクリプトを失っており、**Artifact 上の画像が唯一の原本** (Artifact tool の read で `img/<name>.png` を取得できる)。
- 数字は `case/53.c3x_vane_cht/tools/report_numbers.py` (領域別偏差・帯内率・準定常判定) と `plot_tw.py` (連成壁温) の出力から転記する。
