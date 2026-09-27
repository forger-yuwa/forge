# V-ax4 再開契約の負例 (2026-09-27、codex result M1 の修正後。ローカル float ビルド、各 1 step)

状態ファイルは平面 = case/58 run_0012_v6p_df5_i50_nlog/conjugate_state_5.h5 (属性なしの旧状態。p_full は 3 属性を付与)、軸対称 = case/61 run_0005_annulus_r32s16_w20k/conjugate_state_4.h5。1 つずつ属性を消して再開。

| ケース | C++ (forge の再開) | 評価器 check_cht_balance._state_geometry |
| --- | --- | --- |
| p_legacy | 受理 (復元) | (下の表) |
| p_full | 受理 (復元) | (下の表) |
| p_nogeo | 拒否: ERROR: conjugate_state_5.h5 の識別属性が一部だけある (geometry 無 / | (下の表) |
| p_nounit | 拒否: ERROR: conjugate_state_5.h5 の識別属性が一部だけある (geometry 有 / | (下の表) |
| p_nocon | 拒否: ERROR: conjugate_state_5.h5 の識別属性が一部だけある (geometry 有 / | (下の表) |
| a_full | 受理 (復元) | (下の表) |
| a_nogeo | 拒否: ERROR: conjugate_state_4.h5 の識別属性が一部だけある (geometry 無 / | (下の表) |
| a_nounit | 拒否: ERROR: conjugate_state_4.h5 の識別属性が一部だけある (geometry 有 / | (下の表) |
| a_nocon | 拒否: ERROR: conjugate_state_4.h5 の識別属性が一部だけある (geometry 有 / | (下の表) |

平面の p_legacy / p_full は復元後 step 0 で固体温度の安全停止 (流体場を引き継がない再開による既知の起動の罠、旧 plan #99 ⑤)。契約の判定とは無関係。

```
p_legacy  -> (False, None)
p_full    -> (False, None)
p_nogeo   -> (False, "conjugate_state_5.h5 の識別属性が一部だけある (geometry=None, load_unit='W/m', state_contract=2)。壊れた状態として判定しない")
p_nounit  -> (False, "conjugate_state_5.h5 の識別属性が一部だけある (geometry='planar', load_unit=None, state_contract=2)。壊れた状態として判定しない")
p_nocon   -> (False, "conjugate_state_5.h5 の識別属性が一部だけある (geometry='planar', load_unit='W/m', state_contract=None)。壊れた状態として判定しない")
a_full    -> (True, None)
a_nogeo   -> (False, "conjugate_state_4.h5 の識別属性が一部だけある (geometry=None, load_unit='W/rad', state_contract=2)。壊れた状態として判定しない")
a_nounit  -> (False, "conjugate_state_4.h5 の識別属性が一部だけある (geometry='axisym0', load_unit=None, state_contract=2)。壊れた状態として判定しない")
a_nocon   -> (False, "conjugate_state_4.h5 の識別属性が一部だけある (geometry='axisym0', load_unit='W/rad', state_contract=None)。壊れた状態として判定しない")
```
