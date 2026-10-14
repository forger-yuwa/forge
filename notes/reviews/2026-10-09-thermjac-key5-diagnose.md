# codex 諮問 (diagnose): thermjac-key5

- **brief**: [`notes/reviews/briefs/2026-10-09-thermjac-key5.md`](../../notes/reviews/briefs/2026-10-09-thermjac-key5.md)
- **plan**: [`plans/active/time_integration-implicit-thermal-jacobian.md`](../../plans/active/time_integration-implicit-thermal-jacobian.md)
- **date**: 2026-10-09
- **commit**: `831d0b41` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.4 min, rc=0
- **結論**: **次の計算は、同一バイナリ・同一起点の point、`cfl_pseudo: 4` で `implicitThermalJacobian: 1 → 5` だけを変える最大2000 stepのA/Bに絞り、下流の温度補正と運動量残差の成長を計測する。**
- **extra**: `solver_density_cuda/cuda_forge/timeIntegration_d.cu`, `methods/time_integration/implementation.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 判定 | 根拠と対案 |
|---|---|---|
| 行ごとの減衰の不整合という見立て | **採用。ただし主因認定は要再検証［Major］** | 局所モデルの式は正しい。しかし、実装には対流 A⁺ と近傍補正がある。[timeIntegration_d.cu:875](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:875)。実際の行列・残差に対する温度補正で検証する。 |
| キー 5 の限定試験 | **採用** | スカラー A を全行に戻す実装は案 (b) と一致する。[timeIntegration_d.cu:969](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:969)。既定化・一般的な安定性の保証には進めない。 |
| 過去のキー 0 との比較による帰属 | **要再検証［Major］** | [plan:124](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:124) 自身が比較の前提にした V0 が未実施。直近の診断は同一バイナリのキー 1／5 で行い、`run_0191` を使う後退判定は V0 合格後に行う。 |
| Vb-point・Vb-dir の判定 | **修正して採用［Major］** | [plan:158](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:158) の5000 stepごとの保存では、2000 stepの早期判定中に場を追えない。既存527節点の帳簿も、今回壊れた下流79〜89層目を含まない。序盤出力と監視領域を追加する。 |
| U1 による Jacobian 検証の代替 | **却下［Major］** | 静止CPGの純伝導では、高速TPの交差項を検査できない。[plan:118](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:118)。過去レビューで求めた非零速度・TP・異なるエネルギー基準の温度微分／行列作用の検査は未消化として残す。 |

結論: **次の計算は、同一バイナリ・同一起点の point、`cfl_pseudo: 4` で `implicitThermalJacobian: 1 → 5` だけを変える最大2000 stepのA/Bに絞り、下流の温度補正と運動量残差の成長を計測する。**

第 1 仮説: **キー 1 がエネルギー行だけからスカラー減衰 A を除いたため、保存量補正の比率が崩れ、高速域で温度・圧力への余計な結合を生み、残差更新とのフィードバックで成長する。** 確度: **中**。

  根拠:

`t = V/Δτ`、`K = ½|u|²`、`z = (0,0,0,0,1)ᵀ`、`w = (K−e,−u,−v,−w_velocity,1)` と置く。温度補正は `δT = w·δQ/(ρc_v)`。対流・近傍結合を省いた内部点では、実装は次の行列になる。

```text
キー 0: D₀ = (t + A)I
キー 1: D₁ = tI + A·diag(1,1,1,1,0) + Bzw
キー 5: D₅ = (t + A)I + Bzw
```

`R = R₀(1,u,v,w_velocity,E)ᵀ` に対して、キー 1 は確かに

```text
δT = (E/c_v)·(δρ/ρ)·A/(t+B)
```

を生む。提供された値から独立に5×5行列を組んで解くと、`δρ/ρ = 0.01` に対し、キー 0 は丸め誤差程度、キー 1 は **8.0075 K**、キー 5 は **0 K** となった。これは**局所モデルの計算結果であり、run の実測ではない**。

ただし、同じモデルのキー 1 の固有値は `1.71` が4個、`1.54` が1個で、すべて正である。**8 Kの交差応答だけでは、発散や行列の特異化を証明しない。** 実際の残差・EOS・空間結合を通じて、その補正が増幅されることを確認する必要がある。

提供記録の「高速・高 μ_t の層で壊れる」はこの機構と整合する。一方、922 stepの非有限位置は**最初に不安定モードが発生した位置とは限らない**。また、係数の見積もりは `run_0208/res_200000`、発散試験の起点は `run_0183/res_100000` であり、異なる状態である。[plan:164](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:164)

**対流 A⁺を入れた場合:** 不整合の存在は変わらないが、8 Kという大きさや符号は保証されない。対流込みのキー 0 の局所行列を `D₀`、同じ右辺を `b` とすると、キー 5 は `D₅ = D₀+Bzw` なので、正則な場合は

```text
w·D₅⁻¹b = (w·D₀⁻¹b) / (1 + B·w·D₀⁻¹z)
```

となる。したがって、**キー 0 の解が温度を変えないなら、同じ右辺のキー 5 も温度を変えない**。しかし、`w·b=0` だけでは、対流を解いた後も温度不変とはいえない。実際には近傍右辺もsweepごとに変わるため、局所式だけで全反復の安定性は決まらない。

**二重減衰について:** 対流を省いたモデルでは `H=zw` は `H²=H` を満たし、キー 5 の温度不変な4方向は `t+A`、温度を変える1方向は `t+A+B` で抑えられる。この範囲では、温度補正をさらに小さくする合理的な前処理である。RHSを変えないため根も変えない。ただし、対流・ライン・有限sweepを含む系について、**「遅くなるだけで、安定性への別の影響はない」とまでは断定できない**。

  反証条件: 同じ状態・右辺で対流を含む補正を比較しても、キー 1 に特有の温度応答が成長箇所に現れない、またはキー 5 がその応答を除いても同じモードが成長するなら、**この不整合を今回の主因・十分な説明とする仮説を退ける**。キー 5 が安定化しただけでは、以下の第2仮説との区別はつかない。

第 2 仮説: **キー 1 は熱伝導を加える一方、エネルギー行の粘性仕事に対する従来の安定化も失った。** 確度: 中〜低。RHSには `τ·u_face` があるが、新しい行は温度微分だけである。[viscousFlux_d.cu:336](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:336)。高速・高 μ_t 域との対応はこの説でも説明できる。対案は、熱伝導と粘性仕事を分離し、それぞれが補正へ与える作用を測ること。

第 3 仮説: **directionalでは、時間対角の縮小によって、残る空間結合の近似誤差が別の成長モードを許している。** 確度: 中〜低、未確認。ラインのKは対流FVSから構成され、熱伝導の近傍結合は追加されていない。[timeIntegration_d.cu:906](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:906)。対案は、対流段の増大を原因と決めず、固定状態で内側線形残差と外側の増幅を分けて調べること。

判別 A/B: **A＝キー 1、B＝キー 5。point・`cfl_pseudo: 4`、最大2000 step。変更する数値設定はこのキーだけ。**

- 同じ新バイナリ、同じ `run_0183/res_100000` の保存量、同じ精度・リミッタ基準・BC・sweep数・緩和を使う。既存runを上書きせず、実行先で空き番号を確認して新規runを作る。
- 全残差を毎step、場を少なくとも100 stepごとに保存する。局所帳簿には既存の縮流部に加え、`x_w≈70`・79〜89層目と周囲を含める。
- `after_eos_bc` に位相を揃え、`δρ/ρ`、`δT/T`、`δP/P`、速度補正を初期の固定尺度で測る。序盤には同じ状態・右辺で両行列の補正も比較し、熱伝導・粘性仕事・対流の寄与を分離する。
- **Aが成長し、Bだけが成長を止める場合:** 「エネルギー行のAを戻しても救えない」を棄却する。局所温度応答の差が成長に先行すれば第1仮説を支持するが、第2仮説は残る。
- **Bでも同じモードが成長する場合:** 「Aの欠落を戻せば今回の不安定性を解消できる」という単独説を棄却する。
- Aが元の成長を再現しない場合は判別不能。バイナリ・実効設定・restartの同一性へ戻る。

**Vbの合否には追加が必要。** 「全残差が開始時の10倍未満」は停止条件としては使えるが、9倍まで増大し続けても通るため、安定判定として不十分である。末尾500 stepのトレンドと局所振幅の成長率、全域の非有限・非物理値、線形solve失敗件数を併記する。開始値がゼロに近い残差列には固定した絶対尺度を使う。`check_convergence.py` の区間・VERDICTも保存し、短期合格は「2000 stepで成長を検出しなかった」に限定する。

**pointで安定し、directionalで対流段が育つ場合に次に見るものは、成長箇所でのsweep別の補正と線形残差である。** 内側反復自体が増幅しているのか、内側は解けているのにcommit後の残差が増えるのかを分ける。後者なら、sweep追加ではなく近似LHSと実残差の作用の差が問題になる。キー 0 の近壁熱モードに対する過去の精度・sweep試験を、キー 5 の別モードの除外証拠にはしない。

やらない方がよいこと:

- キー 5 の短期安定を、提案機構の確定・収束・高速化成功と呼ぶこと。
- 対流段が大きくなったという理由だけで、流束スキーム変更へ進むこと。
- 初期200 stepの残差低下や、欠損の総和だけで採用すること。
- 対流込みの確認をせず、8 Kの見積もりを実際の補正量として使うこと。

呼び出し側の前提への異議:

- **［Major］観測と解釈が混在している。** 「粘性段の成長が弱まった」は提供された観測だが、「熱伝導が駆動するモードを抑えた」は解釈である。[plan:151](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:151)。対案は `wi_eheat` と `wi_ework` の分離記録を根拠にすること。
- **［Minor］TPで壁際の `E/c_v≈T` は一般に成立しない。** エネルギー基準と温度依存比熱があるためである。また `ρc_vδT` と `δ(ρc_vT)` は異なる。[plan:89](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:89)。対案は、前者を実際の `e(T_w)/c_v(T_w)` に置き換え、温度行を `Bρc_vδT` と記述すること。
- **［Major］NaN位置と μ_t 最大位置の対応だけでは発生源を同定できない。** 対案は、発散run自身の初期成長区間で係数と補正を測ること。

不足情報: 対象の `case/45.isobutane_m6_d155/run_0183*`、`run_0191*`、`run_0203〜0219*` はこの作業ツリーに存在しない。実効YAML、バイナリ照合、残差・HDF5・帳簿、各VERDICTを直接再検証できていない。発散時刻・位置は提供記録として扱った。現在のHEADは `831d0b41` だが、指定プロンプトに収録されたplan・`timeIntegration_d.cu`・methods本文は現行ファイルと一致することを確認した。

ファイル変更・forge起動なし。**plan未反映。呼び出し側が §4.3・§6.0・§6.1 に反映する診断結果である。**
