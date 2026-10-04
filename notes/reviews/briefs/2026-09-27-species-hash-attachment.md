# 諮問ブリーフ: 場の species ハッシュを「誰が付けるか」(2026-09-27)

関連 plan: `plans/active/thermophysics-solver-owned-species-db.md` §4.3・§5.1 #3・§6 V1、仕様 `methods/thermophysics.md` §1b.4。

## 問い (1 つ)

保存場・初期場に付ける species ハッシュ属性について、**誰がいつ付けるか (案 A/B/C またはそれ以外) を 1 つに絞って推奨**し、
その案で「属性が無い場」をどう扱うかと、ハッシュに含める項目 (`source` を含めるか、datum) を決めてほしい。
推奨案で V1 の各試験 (下) が成立するかも示してほしい。

## 観測事実

- #3 の実装を委譲した implementer が着手前に停止して返した (編集なし)。理由: 「属性の無い場は照合不能として終了、`physProp.speciesAllowUnverified: 1` のときだけ通す」とすると、
  **旧い保存場だけでなく新規 TP run の初期場もすべて止まる**。
  - TP の初期場はソルバも変換器も作らない: `solver_density_cuda/input/setInitial.hpp:32-216` は CPG 式のみ。
    実際は Python (`design/forge_design/evaluate/ic.py:51-85` `paste_isentropic_ic` 等、case 内の `gen_*_ic.py`) が変換器出力 `nozzle.h5` の `/VALUE` を上書きして作る。
    変換器 (`solver_density_cuda/mesh/convertGmshToForge.cpp:64`) も Python も属性を書かない。
  - SERN の段間引き継ぎ `design/forge_design/evaluate/runner_sern.py:626-632` (`restart_by_index`) は index コピーで属性を移さない。
  - 設計 runner の config は `valueFileName: nozzle.h5` を直接読む (`runner.py:87`, `runner_wt.py:178/226/256`, `runner_sern.py:252`)。
  - 回避にフラグを生成 config へ書くと、段・restart の複製で引き継がれて照合が恒常的に切れる。
- 変換器は `convertGmshToForge.cpp:28` で `speciesDB_init` を呼ぶので、ソルバと同じ resolver でハッシュを計算できる。
- plan §6 V1 の試験: A 同一設定 → 許可 / B 内蔵 (or 外部 DB) 係数 +0.001 → 拒否 / (a) 保存場を別 run にコピー → 拒否 / (b) 記録の取り違え → 拒否 /
  (c) 外部 DB 変更 → 拒否 / (d) 起動前の宛先解決 / (e) 液相エントリだけ変更 → 拒否 (#10 で) / 記録の無い旧い場 → 照合不能で止まり明示フラグでのみ通る。

## implementer が出した案

- **案 A**: ソルバに「記録だけ書いて終わる」モード (`forge --resolve-species`: GPU を使わず `resolved_species_<h16>.yaml` を書きハッシュを表示して終了) を足し、
  Python の IC 作成 (`ic.py`, `case/*/gen_*_ic.py`) と runner の段間処理がその結果で属性を付ける。V1 (d) を満たせる。範囲: `ic.py`・各 runner。
- **案 B**: 変換器が自分の解決結果のハッシュを出力 h5 に付ける。`restart_field.py` / `interp_field.py` / `restart_by_index` は DST の属性を SRC のもので置換し、SRC に無ければ DST の属性を**削除**
  (削除しないと旧い場が変換器のハッシュで素通り)。弱点: Python が上書きした初期場も変換器のハッシュのまま通る (初期場の熱物性不整合は起動過渡の問題で restart の忠実性とは別)。
  既存 run の旧 `nozzle.h5` は 1 回だけフラグが要る。範囲: `convertGmshToForge.cpp`・`runner_sern.py`。
- **案 C**: 属性が無いものは初期場とみなして通し、ツールが属性の無い旧 res を移すときだけ `species_hash: "unverified-legacy"` を書いてソルバがそれを止める。
  弱点: 旧 `res_*.h5` を直接 `valueFileName` にした場合と既に restart 済みの既存 `nozzle.h5` は素通り。

小さな判断:
- `source` (builtin/file) をハッシュに含めると、係数がビット同一でも出所の表記だけで restart が拒否され通す手段がない (#4 の共通データ化で内蔵種の出所表記が変わると全 run 不一致)。記録には残しハッシュから外す案。
- datum: 係数は datum 前の絶対基準で記録し `thermoHrefTemp` をハッシュに含める (datum が違えば roe の意味が変わるので拒否すべき)。
- CPG (`thermalMethod != 2`): 記録・属性・照合をしない。
- 凝縮の液相物性: 今回は記録に含めない (#10 で)。

## 当方の見立て (棄却してよい)

- 案 B を基本に、案 A の「resolve-only モード」を V1 (d) 用に後から足すのが範囲と効果の釣り合いがよいのではないか。
  変換器が付けたハッシュは「このメッシュ h5 はこの物性の config で作られた」を意味し、IC を上書きする Python は物性を変えないので、属性の意味は保たれる。
- `source` はハッシュから外す、datum (`thermoHrefTemp`) は含める、に賛成。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの

- `sed -n '70,140p' plans/active/thermophysics-solver-owned-species-db.md`
- `sed -n '186,240p' methods/thermophysics.md`
- `sed -n '20,70p' solver_density_cuda/mesh/convertGmshToForge.cpp`
- `sed -n '40,90p' design/forge_design/evaluate/ic.py`
- `sed -n '615,640p' design/forge_design/evaluate/runner_sern.py`
- `sed -n '170,230p' solver_density_cuda/tools/forge_species.py`
- `sed -n '1,60p' solver_density_cuda/tools/restart_field.py`
- `sed -n '55,80p' solver_density_cuda/tools/interp_field.py`
