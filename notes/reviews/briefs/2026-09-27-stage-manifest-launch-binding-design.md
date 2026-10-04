# 諮問ブリーフ: 段の区間判定を起動記録に結び付ける plan の §4・§6 確定 (2026-09-27)

plan: `plans/active/tooling-stage-manifest-launch-binding.md` (draft、§4・§6 は案)。親: `plans/active/gradient-scalar-lsq-unification.md` §4.4・§5.1 #2g・#6。

## 決めてほしいこと

§4 (設計 1–5) と §6 (T1–T6) をこのまま確定してよいか。直すならどこをどう直すか。特に:
1. バイナリ id: CMake 埋め込みの git rev と `/proc/self/exe` の sha256 のどちらを正本にするか (両方か)。起動ごとに実行ファイルを読むコストは許容してよいか。
2. 段と起動の単調対応 (設定ハッシュ + 起動順) で十分か。manifest の段数と起動の回数がずれる (途中で落ちて再起動した等) 場合の扱い。
3. バイナリが違えば常に別区間にする (同一ソースの再ビルドでも分かれる) 保守側でよいか。
4. T1 用に「既定だけ lsq に変えた試験ビルド」を作る方法でよいか。

## 読んでよいファイル

plan 2 本、`solver_density_cuda/tools/stage_manifest.py` 全体 (341 行)、`solver_density_cuda/main.cpp:90-125` (`appendLaunchRecord`)、`solver_density_cuda/tools/test_stage_manifest_*.py`、`design/forge_design/evaluate/runner_sern.py:30-60`。

## 観測事実

- `load_launches` は `cfg_fnv` → 最後の起動の chi 値 (`stage_manifest.py:152-163`)、`_normalize` がそれで同じ `cfg_fnv` の全段を上書き (`:279-295`)。
- 起動記録には `exe_size`・`exe_mtime` はあるがバイナリ内容の id は無い。`RUN_PROVENANCE.txt` の `forge_sha256` は run_case.sh 側で、起動記録と結び付いていない。
- 2026-09-26 に暫定で `mesh.scalarGradient` を YAML hard キーに追加し、PyYAML 欠如時は停止にした (回帰試験 5 例 PASS)。

## 仮説

- 単調対応 + 実効値 + バイナリ id で、codex M7 の「旧既定と新既定を誤連結」は防げる。
