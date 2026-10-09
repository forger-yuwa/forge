# codex 諮問 (diagnose): line-viscous-resume

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-resume.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-resume.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `5f0b10cf` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.8 min, rc=0
- **結論**: **次はノズルの同一初期場・同一新バイナリで、値3のマスク7／5を比較し、熱伝導Kが早期破綻に必要かを判別する。**
- **extra**: `case/48.flat_plate_cooled_m4/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| E1で熱伝導Dを残し、壁拘束を統一する | **条件付き採用／Major** | [薄層D:129](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:129)と[キー5:1009](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1009)は、正の物性では代数的に同じ熱伝導D。ただし値3は[965行](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:965)でスカラーを既に加える。**キー5の処理全体を追加するとスカラーを二重計上する**。熱伝導Dを共通経路で常に一度だけ加え、ビット2では熱伝導Kだけを切り替える実装を推奨する。 |
| マスク1・2・4を旧runと比較する | **却下／Major** | [plan:244](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:244)では、対照だけマスク追加前のバイナリ。現コードにはマスク自体が未実装。**同じ新バイナリの7／5を先に比較**し、既知の破綻の再現と熱伝導Kの必要性を一組で調べる。 |
| E2の準備スクリプトをそのまま実行する | **却下／Major** | [e2_plate.sh:27](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/e2_plate.sh:27)はrestartの終了コードを確認せず、最後の`rm`が成功すると`prep`も成功する。**初期場の移送失敗後でも起動できる**。restart失敗時は即returnし、保存量の移送検査とメッシュ品質確認を起動条件にする。 |
| E2のA／Bで粘性結合だけを判別する | **要再検証／Major** | [e2_plate.sh:32](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/e2_plate.sh:32)は両側キー5だが、値2だけ[壁拘束:1146](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1146)が変わる。ノズルでの除外結果を平板へ移せない。将来E2を行うなら両側キー7で拘束を統一する。 |
| 平板で非再現なら「ノズル固有条件が必要」 | **却下／Major** | [plan:249](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:249)ではCFL・初期場・格子・境界条件も同時に違う。結論は**「この平板条件では2000 step以内に再現しなかった」まで**。ARだけでなく実際のΔτとV/Δτを比較する。 |
| E1・E2の6本同時投入 | **今回は却下／Minor** | [E2:31](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/e2_plate.sh:31)は既に3本並行起動する構成だが、上記の準備・帰属問題が未解決。**まずE1一組を逐次実行**する。並行実行そのものを数値的不正とはしないが、速度評価には使えない。 |

結論: **次はノズルの同一初期場・同一新バイナリで、値3のマスク7／5を比較し、熱伝導Kが早期破綻に必要かを判別する。**

第 1 仮説: **熱伝導Kの追加が、運動量・仕事の結合および壁拘束との相互作用を通じ、今回の早期破綻に必要になっている。** 確度: **低、未確認**。  
  根拠: [熱伝導K:138](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:138)は、既存の熱伝導Dに対して新たに加わる近傍結合で、独立に除去できる。planの記録では、`case/45.isobutane_m6_d155/run_0301_vn1b_lvc3/`は122 stepで非有限、一方`run_0304_wallB_tj7/`は2000 step有限。ただし両者の差には運動量・仕事も含まれ、**熱伝導を他の項より疑う実測上の根拠はまだない**。この順序を選ぶ理由は、一組で相互作用を含めた必要性を検査できるため。  
  反証条件: **マスク7で破綻が再現し、熱伝導Kを除いたマスク5も2000 step以内に破綻すること**。「熱伝導Kが破綻に必要」を棄却する。ただし、別の機構による破綻かは発生場所・順序で確認する。

第 2 仮説: **運動量・仕事のD/Kだけでも破綻し、熱伝導Kは発生時刻を変えるにとどまる。** 確度: 低、未確認。

判別 A/B:

- **A＝マスク7、B＝マスク5（1＋4）**。熱伝導D・全行スカラー・壁拘束・運動量と仕事のD/Kは同一にし、**熱伝導Kの有無だけ**を変える。
- 両側とも`run_0183_ns_coldmesh_tw300_ext/res_100000.h5`、値3、キー5、方向別、上限0、`cfl_pseudo: 4`、緩和0.7、5 sweep、最大2000 step。同一SHAのバイナリと同一入力を使う。
- 初回の同一状態で行列差を確認する。**D・rhs・拘束行は同一、Kの差は熱伝導成分だけ**でなければ、この比較は無効。基準合わせにはマスク0が値0＋キー7の行列を再現する確認も必要。
- 全残差を毎step、序盤200 stepは既存の局所帳簿を継続する。全域の最初の非有限・非正のρ/P/T、EOS床到達、δρ/ρ・δu・δTの増幅場所を記録する。200 stepごとの場だけでは、既知の122 step破綻の前後関係を追えない。

**→ Aだけ破綻、Bは2000 step有限・非物理値なし:** この条件・期間で、熱伝導Kを除くことが破綻回避に十分。第1仮説を支持するが、熱伝導単独の誤りや長期安定は証明しない。  
**→ A・Bとも破綻:** 第1仮説の必要性を棄却し、運動量・仕事側を候補に残す。  
**→ Aが破綻を再現しない:** 帰属不能。旧`run_0301`を代替対照にしない。

両側の`check_convergence.py`のVERDICTと実際の判定区間を保存する。これは短期破綻の試験であり、有限で完走しても収束・準定常とは呼ばない。

やらない方がよいこと: **最初からマスク1・2・4と平板3本を投入すること。** 単独項の十分性探索より、まず既知の破綻から熱伝導Kだけを引く。マスク5はこの目的で最初から入れる。値2・3の本線採用、既定変更、コード削除にも進まない。

呼び出し側の前提への異議:

- **［Major］「FP32場をFP64格子に載せるだけ」は無影響ではない。** 有限のFP32保存量のFP64への拡張自体は正確だが、幾何の再計算・現行離散化・物性再評価による残差は変わり得る。[restart_field.py:74](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/restart_field.py:74)は座標の完全一致を要求するため、再変換で拒否される可能性もある。対案は節点対応・接続・座標差を確認し、保存量を再構成せず移すこと。`centCoords`等の置換で検査を通してはいけない。
- **［Major］AR≈900／4000はΔτの伸びそのものではない。** [setDT_d.cu:82](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:82)には速度・音速・粘性が入り、[167行](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:167)でライン面を除外し、軸対称では別の寄与も加わる。平板とノズルではMach分布、壁の曲率、入口条件、初期残差も違う。将来の単因子候補は**実測Δτに基づく`lineDtDirectionalCap`**とし、先に格子ARを変えない。
- **［Minor］平板起点を「収束場」と呼ぶ根拠が不足している。** [README:34](/home/sano/work/forge-integ-1005/case/48.flat_plate_cooled_m4/README.md:34)は残差床と派生量の定常性を記録している。対案は「既存の発達場」とし、残差収束と目的量の定常性を別々に示すこと。

不足情報: 対象`run_0183`・`run_0300/0301/0303/0304`・平板`run_0005`の原本がこの作業ツリーにないため、残差CSV・HDF5・VERDICTを独立確認できなかった。plan記載の`run_0303/0304`の判定は、区間0〜1999でともに **NOT CONVERGED**。発散時刻・有限性・除外判断は現時点では記録に依存する。また、FP64ビルドでも[線形組立精度の既定はfloat](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1681)なので、実効`implicitSolvePrecision`も必要。

ファイル変更・forge起動なし。**plan未反映**。呼び出し側で§6.8に採否と上記分岐を反映すること。
