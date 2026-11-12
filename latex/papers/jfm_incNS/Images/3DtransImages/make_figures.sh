#!/bin/sh
# Regenerate the transition-sequence panels of fig:3dtransAtGap and the panel of
# fig:3dtransShallowGap.
# Set RENDER=1 to also re-render the PNGs from the DNS output (slow).
set -e
cd "$(dirname "$0")"
DNS=/home/victor/hosts/hpc/PhD/runs/quasi3DGapRe1000inc/baseflow/dns/d4_w19_Lz20_k32_BFwith0.01Times3dMode

# index  fade strength of the downstream plate  render flags
while read -r i s flags; do
	n=d4w19_3dMode_L2crit_$i
	[ "${RENDER:-0}" = 1 ] && pvbatch render.py "$DNS/mesh_${i}_pert_L2Crit.vtu" "pngs/$n.png" --left 0.07 $flags
	python3 gen_tex.py "pngs/$n.png" "pdfs/$n.tex" --fade-strength "$s"
	(cd pdfs && pdflatex -interaction=batchmode "$n.tex" >/dev/null && rm -f "$n.aux" "$n.log")
	echo "$n"
done <<END
11 1.00
16 0.88
17 0.76
19 0.64 --expand
20 0.52 --expand
23 0.40 --expand
END

# shallow-gap case of fig:3dtransShallowGap
DNS=/home/victor/hosts/hpc/PhD/runs/quasi3DGapRe1000inc/baseflow/dns/d1.5_w80_Lz5_k32
n=d1.5w80Lz5_3d_L2crit_31
export GAP_CASE=d1.5w80
[ "${RENDER:-0}" = 1 ] && pvbatch render.py "$DNS/mesh_31.vtu" "pngs/$n.png" --expand
python3 gen_tex.py "pngs/$n.png" "pdfs/$n.tex"
(cd pdfs && pdflatex -interaction=batchmode "$n.tex" >/dev/null && rm -f "$n.aux" "$n.log")
echo "$n"

# forced region-III case now shown in fig:3dtransShallowGap (seed amplitude 1e-3, chk 225;
# the .vtu is made with FieldConvert -m L2Criterion mesh.xml session.xml mesh_225.chk mesh_225.vtu)
DNS=/home/victor/hosts/hpc/PhD/runs/quasi3DGapRe1000inc/baseflow/dns/d1.5_w70_Lz6.0_k32_seed1e-3
n=d1.5w70Lz6_seed1e-3_L2crit_225
export GAP_CASE=d1.5w70
[ "${RENDER:-0}" = 1 ] && pvbatch render.py "$DNS/mesh_225.vtu" "pngs/$n.png" --expand
python3 gen_tex.py "pngs/$n.png" "pdfs/$n.tex" --front-faces --fade-start 0.85
(cd pdfs && pdflatex -interaction=batchmode "$n.tex" >/dev/null && rm -f "$n.aux" "$n.log")
echo "$n"
