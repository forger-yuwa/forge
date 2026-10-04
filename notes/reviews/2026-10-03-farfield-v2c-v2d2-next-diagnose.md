# codex 諮問 (diagnose): farfield-v2c-v2d2-next

- **brief**: [`notes/reviews/briefs/2026-10-03-farfield-v2c-v2d2-next.md`](../../notes/reviews/briefs/2026-10-03-farfield-v2c-v2d2-next.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-10-03
- **commit**: `d70d897f` (feature/sern-design)
- **codex**: effort `high`, 5.8 min, rc=0
- **結論**: **直近の一手は V2d-2 の最小 dt における `nSub` 20／40比較とし、V2c は形状を再変更せず一次化への感度を調べる。**
- **extra**: `case/58.farfield_verification/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0／Major 4／Minor 1）

対象 `run_*` は手元に **0件**でした。以下の実測値はブリーフ・台帳の報告値として扱い、独立には再認定していません。コードの確認と、評価器への合成入力による検査を行いました。`forge` の起動・ファイル変更はありません。

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | V2c を「凸角の欠陥、farfield 無関係」とする解釈：**却下** | 凸角を除いた v2 でも発散・停滞が報告されています。また A/C も入口は `farfield` です。[setup_v2c.py:76](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2c.py:76)。上面 slip 対照は「farfield を含まない対照」ではありません。**形状と BC を固定し、まず再構成次数への依存を調べる**方針を採用します。候補(d)の陰解法は既に使用中です（同:65、68）。 |
| **Major** | `limiterScaled: 0` でリミッタ係数の切替だけを診断：**却下** | この設定は無次元化だけでなく、再構成との評価点整合も切り替えます。[solverConfig.cpp:685](/home/sano/work/forge-sern-design/solver_density_cuda/input/solverConfig.cpp:685)。旧式の K も固定値です。**`convMethod: 1→0` を最初の診断変更に限定**し、成功しても「二次再構成への感度」とだけ解釈します。係数切替が真因とは認定しません。 |
| **Major** | V2d-2 の誤差増大を float32 に確定：**要再検証** | BDF 項と状態への加算は実際に `flow_float` です。[update_d.cu:459](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/update_d.cu:459)、[同:517](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/update_d.cu:517)。ただし TP 条件で `nSub` 感度が未確認です。**先に固定 dt で20／40を比較**します。CPG の V2a 合格は代用できません。 |
| **Major** | 現行評価器の PASS を必要条件の確認込みとする扱い：**却下** | `same` は終了時刻の短い方まで比較するだけです。[eval_v2d.py:98](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2d.py:98)。合成入力で、入射約0.615 ms・反射約1.495 msに対し両短領域を0.1 msで打ち切っても **`VERDICT: PASS`** を再現しました。**固定評価区間の完全被覆、時刻の単調性、全入力の有限性、有限かつ正の入射振幅を必須にする**必要があります。この欠陥が既報 FAIL の原因だったとは主張しません。 |
| **Minor** | 「1e−3 P∞＝2.2 Pa」：**訂正を採用** | 入力振幅は **2.851 Pa** です。[setup_v2d.py:129](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2d.py:129)。2.2 Pa が評価点で減衰した入射振幅なら、その実測値として別記してください。時間精度の分母は長領域の入射窓振幅です（[eval_v2d.py:96](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2d.py:96)）。 |

結論: **直近の一手は V2d-2 の最小 dt における `nSub` 20／40比較とし、V2c は形状を再変更せず一次化への感度を調べる。**

第 1 仮説: **V2d-2 は、微小擾乱に対する保存量・BDF 項・更新の float32 丸めが時間誤差を支配し始めている。** 確度: **中**

根拠: `flow_float` は float（[flowFormat.hpp:6](/home/sano/work/forge-sern-design/solver_density_cuda/flowFormat.hpp:6)）。BDF2 は `aQ−bQn+cQnn` をその型で計算し、補正も現在値へ直接加算します。報告された dt 細分での差の増大とは整合しますが、**BDF の減算、更新消失、EOS、流束積算のどれが支配的かは未確認**です。圧力1 ulp＝0.000244140625 Paだけから、累積誤差0.01〜0.05 Paの発生源は特定できません。

反証条件: 演算精度だけを十分に上げた比較でも、同じ評価区間の dt 感度がほぼ変わらなければ「丸め支配」を棄却する。なお、typedef の変更だけでは `sqrtf`／`cbrtf` 等が残るため、直ちに「全域 FP64」とは呼べません（[limiter_d.cu:379](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiter_d.cu:379)）。

第 2 仮説: **V2d-2 の `nSub=20` は要求精度に不足している。** 確度: **中、未確認**。同条件での倍増比較がなく、除外できません。

第 3 仮説: **V2c の停滞は、二次再構成を含む空間離散化と反復の組合せに強く依存する。** 確度: **低〜中**。現在設定は二次風上＋Venkatakrishnanですが、case/16 の実績だけでは本ケースの機構を同定できません。

判別 A/B: **症状ごとに変更は1点に限定する。実施優先は①。**

**① V2d-2：`time.nSubIterDualTime` だけ20→40。**

- `run_0142` 相当の tp2・短領域・`--left-hot`・dt≈2.2e−7 を使い、同じ元のパルス初期場から比較する。`nStepInner` は変更しない。
- 終了時刻は両者とも、登録済みの反射窓末尾まで。同じバイナリ・格子・初期保存量・物性・プローブ節点を照合できれば、20側は既存データを再利用してよい。
- 固定区間全体で `D_N＝max|P₂₀−P₄₀|/A_inc` を測る。`A_inc` は同じ長領域の入射窓から一度だけ決め、差の最大時刻も記録する。
- **D_N≤0.002なら**、この dt における20→40の感度は入射振幅の0.2%以下。反復数依存による説明は弱まり、精度変更の診断へ進む。ただし内部反復誤差そのものの上限証明ではない。
- **D_N>0.01なら**、「20回で十分」を棄却する。反復不足と丸めによる反復依存は、まだ分離できない。
- 中間は判定保留。既存の **dt 半減差≤1%** は別条件として維持する。

**② V2c：`space.convMethod` だけ1→0。**

- 旧形状 B の `case/58.farfield_verification/run_0064_v2c_B/` 最終場を、同一格子の新規2 runへ保存量のままコピーする。A＝二次、B＝一次。形状・BC・CFL・リミッタ基準値・バイナリを固定する。
- 各6000 step。壁圧を500 step間隔で保存し、初期・最終の全節点／全5成分の局所残差を射影なしで評価する。
- 事前指標は全域最大の規格化残差 `Rmax` と、末尾2000 stepの壁圧変動幅 `W＝max_x(max_t p−min_t p)/Δp`。
- **一次側だけ Rmax≤1e−5、W<0.001を満たし、二次側が停滞を再現するなら**、一次化を安定起動に使う方針を支持する。ただし、散逸増加でも改善し得るため「リミッタ係数の切替が真因」とは言わない。
- **一次側でも Rmax・W が二次側の1/2以上で必要条件未達なら**、「一次化だけで停滞を解消できる」を棄却する。それ以外は判定保留。
- これは診断であり、一次側の改善を V2c 合格にしない。収束認定には別途 `check_convergence.py` の PASS、準定常性には `check_quasisteady.py` と登録済み壁圧条件が必要。次数変更前後の残差系列は連結しない。

やらない方がよいこと: **形状・角度・次数・CFLを同時に変えること、全 slip 閉鎖箱を元の通過流の対照にすること、振幅を10倍にして元の微小音響試験を合格扱いすること、ulp比へ合格基準を緩めること。** dt の追加細分だけを続けるのも、現状では診断になりません。

呼び出し側の前提への異議:

- 作用素差と追加射影の影響について、既存診断で棄却した事項を最初からやり直す必要はありません。ただし、その結論を全域・全状態へ一般化しないでください。
- README の `run_0094`–`0096` に残る「離散定常解が2つ」は、後続の絶対残差診断と矛盾します。**「初期場依存の停滞状態」へ訂正**すべきです。
- **accepted への経路は変更しません。** V2c の必要条件と2系列の全線判定、V2d-2 の時間精度、独立参照の精度確認・比較記録、未完の帳簿対応・docs、result レビューが必要です。今回追加するのは評価器の被覆・有限性検査です。V3 や限定運用の承認では代替できません。

不足情報: AWS の実 config、バイナリ識別、プローブ時系列、全残差履歴、局所残差帳簿、壁圧時系列、各 VERDICT 原本が不足しています。中間全場の削除後も、壁圧の全評価節点の系列が残っているか確認が必要です。**今回、既報の収束・反射・時間精度を独立再認定していません。**

**plan 未反映（依頼どおり変更なし）。** 呼び出し側の反映先は `plans/active/boundary-node-farfield-characteristic.md` §5.1 #3c／#3f、§6、§6.1です。
