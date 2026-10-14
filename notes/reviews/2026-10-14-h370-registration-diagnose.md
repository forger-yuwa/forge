# codex 諮問 (diagnose): h370-registration

- **brief**: [`notes/reviews/briefs/2026-10-14-h370-registration.md`](../../notes/reviews/briefs/2026-10-14-h370-registration.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-14
- **commit**: `48951bf5` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 7.0 min, rc=0
- **結論**: **Euler 参照を案 A・`throat_refine: 5` で成立させ、抽出座標と以下の判定条件を固定してから、V-h370 の2本を投入する。**
- **extra**: `case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml`, `case/45.isobutane_m6_d155/_band_ab/half_inlet/h370_prep_check.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 論点 | 重大度・採否 | 根拠と推奨 |
|---|---|---|
| ① Euler 壁のゲート | **Major・案 A を採用** | 設計壁の `monotone_gate` と物理壁の厳密な単調性判定が衝突している。読み取り計算で `min_r1 = −4.0164e−8`、累積減少 `6.1293e−9 r_t`、設計壁ゲート合格と Euler 構築例外を再現した。根拠：[wall_axismach.py:235](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:235)、[runner_axismach.py:1378](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1378)。**H 上流付き設計壁として扱い、下流には設計壁ゲートを適用する。** 上流の単調収縮、正の半径、一意な大域最小、C² 接続は引き続き検査する。壁の保存・復元も対応させる。共有ゲートの緩和、δ_r=0 特例の物理壁クラスへの追加、S の修正は採らない。 |
| ② δ_in の許容超過 | **Major・`throat_refine: 5` を採用** | **今回、生成座標で再計算した。** `4 → 5` だけで max｜δ_in−δ_r｜は **0.001136651 → 0.000796404 mm**、全列被覆を維持し、既存許容 **0.001 mm** を満たした。2000×97 は不変。照合式は [delta_in_check.py:70](/home/sano/work/forge-integ-1005/design/forge_design/metrics/delta_in_check.py:70)、補間は [deltastar.py:475](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:475)。**最初の pass 前に修正する。** 旧 FAIL は残し、修正版を別入力として登録する。許容の事後緩和・抽出器全体の変更・誤差予算への振替は不要。変換後の再照合はまだ必要。 |
| ③ リミッタ基準値 | **Major・今回は C 由来の固定値を維持** | 基準値は再構成増分の無次元化に入り、変更すると数値作用素が変わる：[limiter_d.cu:394](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/limiter_d.cu:394)。実装も restart ごとの再計算による変化を避けるため明示固定を尊重する：[main.cpp:1792](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1792)。準一次元平均との差だけでは不適切と判定できない。**ρ・P・a は現値を今回の全本段で固定**し、pass ごとの場平均への更新はしない。H の値への変更は別の数値感度確認として扱う。 |
| ④ 収束ゲート | **Major・通常判定を維持、横ばいによる必然的 FAIL という前提は却下** | 判定条件は「同一区間のピークから3桁以上低下、末尾が rising でない」であり、**flat 自体は禁止されていない**：[check_convergence.py:193](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_convergence.py:193)。`8.6e−8` という一つの水準だけでは判定できない。**`--segment --drop 3` の全列判定を登録する。** 不合格なら壁更新を止める。`--from-floor` は参照自身の通常判定 PASS が必要で、300 K の未確認床は使えない：[同:359](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_convergence.py:359)。現行 handoff に床指定の配管もない：[deltastar_loop.py:594](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:594)。傾きだけの代替合格や `--allow-unconverged` の常用は採らない。 |
| ⑤ 継続の意味 | **Minor・外れたら数え直す方式を採用** | 「予算内に連続4万 step 保持できる区間を得る」という停止条件として妥当。[stop_rule.py:404](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/stop_rule.py:404) の再計数を合成系列でも確認した。最初の到達、各中断、最終保持開始を保存し、**中断しても上限20万 step は延ばさない**。NaN・データ異常・実行異常は再計数で救済しない。 |
| ⑥ point の追加 | **Major・今回の最小登録には追加しない** | 300 K の比較は未確認で、370 K の根拠にはできない。一方、point も同じ格子・空間離散化を共有するため、共通の空間誤差まで検出する独立参照にはならない。plan が限定した「H・B0 による暫定の設計」の資格確認にとどめる：[plan:206](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:206)。B0 の収束ゲートを免除する代わりに point を足す、という扱いもしない。 |
| 追加：抽出座標 | **Major・絶対寸法に対する除外済み扱いを却下** | 現行抽出は `remap_constant_x=False`。壁法線格子の曲がった列を一定 x として扱う：[deltastar.py:74](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:74)。既往結果は冷却／断熱の**比**への影響0.060%を棄却したもので、δ_E 自体は最大0.70%変わったと記録されている：[plan:267](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:267)。今回の **0.01 mm の絶対寸法評価では一定 x 断面への補間を明示採用**し、停止帯の生成と後段抽出で共通化する。現在の呼び出しには指定がない：[stop_rule.py:667](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/stop_rule.py:667)、[deltastar_loop.py:674](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:674)。端の補間範囲不足を定数外挿で隠して全域合格にしない。 |
| ⑦ 事前登録 | **Major・下記条件で採用** | **Euler 1本＋NS 1本、固定壁、壁更新なし**。Euler を先に資格確認して参照スナップショットを固定する。抽出器は参照ディレクトリの最新出力を読むため、計算中の Euler を参照すると停止帯と後段抽出の参照が変わり得る：[deltastar.py:89](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:89)。 |

結論: **Euler 参照を案 A・`throat_refine: 5` で成立させ、抽出座標と以下の判定条件を固定してから、V-h370 の2本を投入する。**

第 1 仮説: 現在の先行障害は流れの未整定ではなく、Euler 参照の幾何・抽出契約の不整合である。確度: **高**。  
根拠: CFD 0 step で壁ゲート例外と δ_in の許容超過を再現した。NS 生成器の SHA は記録どおり `3980ac54…`。診断用 `HRefWall` による Euler 配点変更では既存許容を満たした。  
反証条件: 案 A の正式な保存・復元経路と変換後 HDF5 で、同じ S・δ_r・尺度・被覆を確認しても許容超過が残ること。その場合、今回の生成座標での説明だけでは不十分となる。

第 2 仮説: 旧 `run_0174` と H の Euler 参照の違いが、共通使用域でも δ_E を0.01 mm超動かす。確度: **中・場では未確認**。

判別 A/B: **同一 NS 保存系列に対し、Euler 参照だけを変える。追加の NS は回さない。**

- **腕 A:** 旧 `run_0174`。実際に覆う共通使用域のみ。
- **腕 B:** 資格確認済みの H の Euler 参照。新しい出口まで含む全域判定はこちらで行う。
- 両腕で NS スナップショット、一定 x への補間、`band_select="edge"`、ω=0.5、平滑化条件を固定する。未緩和 δ_E、帯位置・感度、`held`、次壁の位置・真のスロートを比較する。
- **結果 A:** 共通使用域の全断面・判定窓で max｜δ_E,B−δ_E,A｜≤0.01 mm、抽出の合否も不変 → 第2仮説を退ける。ただし旧参照の出口不足は解消しない。
- **結果 B:** 0.01 mm超、または抽出の合否が変わる → 旧参照の流用を退ける。時間変動や参照の未整定で判別できなければ判定不能。

この A/B を含む **V-h370 の登録内容**は次に絞る。

1. **Euler:** 本段上限60,000 step、保存1,000 step、`cfl_main: 2`。現在の YAML は12,000なので明示変更する。全残差の `PASS`、NaN/Infなし、出口 M・流量・参照場の時系列確認後、使う `res_*.h5` と SHA を固定する。参照の時間変動が抽出差の判定を左右するなら資格未達とする。
2. **NS:** 370 K・H・Gref・既定の B0 固定値。S1/S2 各3,000、本段上限200,000、保存2,500 step。初期 δ_r は `24251cb6…` の表。自動延長しない。
3. **停止設定:** `stop_rule: true`、`stop_rule_cap: 200000`、`stop_rule_hold_steps: 40000`。`x_fracs` は **[0.419738616669543, 0.734542579171700, 0.986385749173426]**。窓20,000、最低9点、r の平均上限 `5e−4`・振幅上限 `5e−3`、θ_r の drift 上限 `5e−4`・振幅上限 `5e−3` を明記する。要求位置と実際に選ばれた断面位置を記録する。保持4万 step は20万の内数。
4. **判定区間:** `stage_manifest` が示す同一設定の本段全体で収束判定。設計量は停止時刻 s に対する **[s−40000, s−20000] と [s−20000, s]**。各窓の全使用域で δ_E・δ_target の最大−最小≤0.01 mm、窓平均同士の差も全域で≤0.01 mm。`check_quasisteady.py` の VERDICT を別途残し、既定の相対閾値だけで寸法予算の合格を代用しない。
5. **受け渡し:** `REACHED` と全残差 `PASS` に加え、保存則、局所 y₁⁺、全域被覆、使用域の `held`・欠測なし、次壁の幾何健全性を要求する。平滑化の λ と昇格上限も投入前に固定する。未収束時に診断抽出を行っても、壁更新の許可とは分ける。
6. **投入前:** 正式な Euler 壁の保存・復元照合、変換後 `delta_in_check`、両格子の `check_mesh_quality.py` を通す。今回の生成座標では Euler の AR **4.984**・skew **0.29568**だが、**変換後品質の VERDICT は未取得**。新 run の番号・相対パス、初期場・バイナリ・実効設定の SHA を登録する。

やらない方がよいこと: **0.001 mm の許容緩和、共有壁ゲートの緩和、未確認床や上書きによる収束合格化、`REACHED` だけでの壁更新、旧 Euler の出口外挿。** 0.05 mm の固定点条件は設計終了用であり、最初の更新開始条件にはしない。

呼び出し側の前提への異議: 「残差が横ばいなら通常判定は必ず FAIL」は実装と異なる。また「固定点のずれ e・残差2e」は、e=δ_in−δ_r、ω=0.5、再平滑化を無視した近似としては成立するが、実装は更新後にも平滑化するため、実際の誤差床の保証にはならない：[deltastar_loop.py:116](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:116)。

不足情報: `run_0174`・`run_0917`・`run_0918` の場・全残差・VERDICT、採用初期場の資格と SHA、変換後の品質・照合結果。300 K の「差≤0.002%」は今回独立検証できていない。**ファイル変更・forge 起動なし。plan 未反映**。反映先は `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.9・§5.1 #36a〜#36c・§6 V-h370。
