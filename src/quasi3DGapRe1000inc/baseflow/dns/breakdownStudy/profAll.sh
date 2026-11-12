#!/bin/bash
cd ~/Desktop/PhD/runs/quasi3DGapRe1000inc/baseflow/dns
for d in d1.5_w65_Lz6.0_k32 d1.5_w66_Lz6.0_k32 d1.5_w70_Lz6.0_k32 d1.5_w66_Lz6.0_k32_seed1e-{9,5,3} d1.5_w70_Lz6.0_k32_seed1e-{9,5,3}; do
  for f in $d/data/avg/pointsPO_n600_chk*.dat; do
    n=$(basename $f .dat); n=${n#pointsPO_n600_chk}
    o=breakdownStudy/profiles/${d}_chk$n.txt
    [ -s $o ] || awk -f breakdownStudy/prof.awk $f > $o
  done
done
echo done
