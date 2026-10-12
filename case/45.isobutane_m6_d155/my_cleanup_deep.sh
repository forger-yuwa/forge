#!/bin/bash
# my_cleanup.sh の 2 回目 (2026-10-13、ユーザ指示「可能な限り容量空けて。要らなそうなものは消していい」)。使い方: bash my_cleanup_deep.sh list|delete
# 1 回目は RUN_RC の無い古い run (2026-10-08 頃まで) を飛ばしていた。ここでは「走行中でなく、最後の書き込みが 12 時間より前」の run も判定済みとして扱う。
# 対象は 1 回目と同じ自分の run 番号 (case/45: 0184〜0499 (0483 を除く)・0590〜0599・0620〜0624・0700〜0803 (走行中の試験 G 0760〜0762 を除く)、case/48: 0057〜0082)
# と、自分の作業ディレクトリ _jprobe・_conv614・_conv9・_stage・_g0c_d10 (10 桁で書いた壊れた G0c、§5.1 #25)。
# 各 run で残す: 最終の res・m9_watch.json の到達の step の res・壁/出口の時系列 (res_wall_*・res_outlet_*)・残差の履歴・判定の json・VERDICT・設定。
# 消す: 上記以外の res_<N>.h5/.xmf (res_0 を含む: 初期場は _tr/q32_init.h5 や引き継ぎ元の最終の場から作り直せる)、入力・診断の h5 の写し
#   (nozzle.h5・mesh.h5・stage.h5・hoop.h5・geomab*.h5・qacc_adopt_*.h5・res_nan_*.h5)、メッシュの写し *.msh。シンボリックリンクは触らない。
# 50 MB を超える ledger*・*.csv・*.log は gzip する。消したもの・gzip したものは _cleanup_<日付>.tsv に記録する。
set -u
MODE=${1:-list}
busy=$(for p in $(pgrep -x forge); do readlink /proc/$p/cwd; done)
C45=~/forge-wallfit/case/45.isobutane_m6_d155
LOG=$C45/_cleanup_$(date +%Y-%m-%d).tsv
now=$(date +%s)
myruns() {
  for C in $C45 ~/forge-wallfit/case/48.flat_plate_cooled_m4; do
    for r in $C/run_0[0-9][0-9][0-9]_*; do
      [ -d $r ] || continue
      n=$(basename $r | cut -c5-8); n=$((10#$n)); c=$(basename $C)
      if [ $c = 45.isobutane_m6_d155 ]; then
        ok=0; { [ $n -ge 184 ] && [ $n -le 499 ] && [ $n -ne 483 ]; } && ok=1
        { [ $n -ge 590 ] && [ $n -le 599 ]; } && ok=1; { [ $n -ge 620 ] && [ $n -le 624 ]; } && ok=1; { [ $n -ge 700 ] && [ $n -le 803 ]; } && ok=1
        { [ $n -ge 760 ] && [ $n -le 762 ]; } && ok=0
      else ok=0; { [ $n -ge 57 ] && [ $n -le 82 ]; } && ok=1; fi
      [ $ok = 1 ] || continue
      echo "$busy" | grep -qx "$r" && continue
      if [ ! -e $r/RUN_RC ]; then
        t=$(find $r -maxdepth 1 -type f -printf '%T@\n' | sort -n | tail -1); t=${t%.*}
        [ -n "$t" ] && [ $((now - t)) -gt 43200 ] || continue          # 印が無く 12 時間以内に書き込みがある run は飛ばす
      fi
      echo $r
    done
  done
}
cands() {
  myruns | while read r; do
    last=$(ls $r 2>/dev/null | sed -n 's/^res_\([0-9]\+\)\.h5$/\1/p' | sort -n | tail -1)
    reach=$(python3 -c "import json; print(json.load(open('$r/m9_watch.json')).get('reach_step') or '')" 2>/dev/null)
    for f in $r/res_*.h5 $r/res_*.xmf; do
      [ -e $f ] && [ ! -L $f ] || continue; s=$(basename $f | sed -n 's/^res_\([0-9]\+\)\.\(h5\|xmf\)$/\1/p')
      [ -n "$s" ] || continue
      { [ "$s" = "$last" ] && [ "$s" != 0 ]; } && continue; [ "$s" = "$reach" ] && continue
      echo $f
    done
    for f in $r/nozzle.h5 $r/mesh.h5 $r/stage.h5 $r/hoop.h5 $r/geomab*.h5 $r/qacc_adopt_*.h5 $r/res_nan_*.h5 $r/*.msh; do [ -e $f ] && [ ! -L $f ] && echo $f; done
  done
  find $C45/_jprobe $C45/_conv614 $C45/_conv9 $C45/_stage -type f \( -name '*.h5' -o -name '*.npz' -o -name '*.msh' -o -name '*.npy' \) -size +1M 2>/dev/null
  [ -d $C45/_g0c_d10 ] && find $C45/_g0c_d10 -type f
}
gzc() { myruns | while read r; do find $r -maxdepth 1 -type f \( -name 'ledger*' -o -name '*.csv' -o -name '*.log' -o -name '*.faces' \) ! -name '*.gz' -size +50M; done; }
if [ $MODE = list ]; then
  cands | xargs -r du -s --apparent-size -B1 | awk '{s+=$1} END {printf "消す: %.1f GB\n", s/1e9}'
  gzc | xargs -r du -s --apparent-size -B1 | awk '{s+=$1} END {printf "gzip する: %.1f GB\n", s/1e9}'
  echo "対象 run: $(myruns | wc -l) 本"
else
  echo -e "# 2 回目に消したもの ($(hostname)、$(date -Is)、M6ノズルの設計の判定済み run と作業ディレクトリ)\tbytes" >> $LOG
  cands | while read f; do b=$(stat -c %s $f); rm -f -- "$f" && echo -e "$f\t$b" >> $LOG; done
  gzc | while read f; do gzip -f "$f" && echo -e "$f\tgzip" >> $LOG; done
  echo "記録: $LOG"; df -h ~ | tail -1
fi
