#!/bin/bash
# M6ノズルの設計セッションの判定済みの run の大きなファイルを消す (2026-10-13、ユーザ指示「可能な限り容量空けて。要らなそうなものは消していい」)。
# 使い方: bash my_cleanup.sh list|delete。~/forge-wallfit/case/45 と case/48 の、自分の run 番号だけを対象にする。
#   case/45: run_0181〜0499 (保護: 0181〜0183 の参照・0483 は forge-25 との約束)、run_0590〜0599、0620〜0624、0700〜0803 (走行中の run は飛ばす)。
#   case/48: run_0057〜0082。forge-8a (05xx の 0500〜0589) と forge-25 (06xx、0620〜0624 を除く) の番号は対象外。
# 各 run で残す: res_0・最終の res (番号が最大)・m9_watch.json の到達の step の res・小さいファイル。消す: それ以外の res_<N>.h5/.xmf、run の中の入力の写し
#   (nozzle.h5・mesh.h5)、stage.h5・hoop.h5 以外の 5 MB を超える診断の h5 は一覧に出すだけ。ログ・残差の履歴は 50 MB を超えたら gzip。消したものは _cleanup_<日付>.tsv に記録。
set -u
MODE=${1:-list}
busy=$(for p in $(pgrep -x forge); do readlink /proc/$p/cwd; done)
LOG=~/forge-wallfit/case/45.isobutane_m6_d155/_cleanup_$(date +%Y-%m-%d).tsv
cands() {
  for C in ~/forge-wallfit/case/45.isobutane_m6_d155 ~/forge-wallfit/case/48.flat_plate_cooled_m4; do
    for r in $C/run_0[0-9][0-9][0-9]_*; do
      [ -d $r ] || continue
      n=$(basename $r | cut -c5-8); n=$((10#$n)); c=$(basename $C)
      if [ $c = 45.isobutane_m6_d155 ]; then
        ok=0; { [ $n -ge 184 ] && [ $n -le 499 ] && [ $n -ne 483 ]; } && ok=1
        { [ $n -ge 590 ] && [ $n -le 599 ]; } && ok=1; { [ $n -ge 620 ] && [ $n -le 624 ]; } && ok=1; { [ $n -ge 700 ] && [ $n -le 803 ]; } && ok=1
      else ok=0; { [ $n -ge 57 ] && [ $n -le 82 ]; } && ok=1; fi
      [ $ok = 1 ] || continue
      echo "$busy" | grep -qx "$r" && continue
      [ -e $r/RUN_RC ] || continue                                      # 終わった印の無い run は飛ばす (走行中・記録の無い古い run)
      { [ $c = 45.isobutane_m6_d155 ] && [ $n -ge 760 ] && [ $n -le 762 ]; } && continue   # 走行中の試験 G (§6.58) ははっきり除く
      last=$(ls $r 2>/dev/null | sed -n 's/^res_\([0-9]\+\)\.h5$/\1/p' | sort -n | tail -1)
      reach=$(python3 -c "import json,sys; print(json.load(open('$r/m9_watch.json')).get('reach_step') or '')" 2>/dev/null)
      for f in $r/res_*.h5 $r/res_*.xmf; do
        [ -e $f ] || continue; b=$(basename $f); s=$(echo $b | sed -n 's/^res_\([0-9]\+\)\.\(h5\|xmf\)$/\1/p')
        [ -n "$s" ] || continue
        { [ "$s" = 0 ] || [ "$s" = "$last" ] || [ "$s" = "$reach" ]; } && continue
        echo $f
      done
      for f in $r/nozzle.h5 $r/mesh.h5 $r/stage.h5; do [ -e $f ] && [ ! -L $f ] && echo $f; done
    done
  done
}
if [ $MODE = list ]; then
  cands | xargs -r du -s --apparent-size -B1 | awk '{s+=$1; n=split($2,a,"/"); t[a[n-1]]+=$1} END {for (k in t) printf "%7.0f MB %s\n", t[k]/1e6, k; printf "合計 %.1f GB (%d run)\n", s/1e9, length(t)}' | sort -rn | head -70
else
  echo -e "# 消したもの ($(hostname)、$(date -Is)、M6ノズルの設計の判定済み run)\tbytes" >> $LOG
  cands | while read f; do b=$(stat -c %s $f); rm -f -- "$f" && echo -e "$f\t$b" >> $LOG; done
  for C in ~/forge-wallfit/case/45.isobutane_m6_d155 ~/forge-wallfit/case/48.flat_plate_cooled_m4; do
    find $C -maxdepth 2 -path "*/run_0*" -type f \( -name forge_run.log -o -name residual_history.csv \) -size +50M 2>/dev/null | while read f; do
      d=$(dirname $f); echo "$busy" | grep -qx "$d" && continue; gzip -f "$f" && echo -e "$f\tgzip" >> $LOG; done
  done
  echo "記録: $LOG"; df -h ~ | tail -1
fi
