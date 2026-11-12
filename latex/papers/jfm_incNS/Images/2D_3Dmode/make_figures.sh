#!/bin/sh
# Regenerate the eigenmode panels of fig:modes and fig:3dfloquetMode with the
# camera and domain styling of Images/3DtransImages.
# Set RENDER=1 to also re-render the PNGs from the eigenmode files (slow).
set -e
cd "$(dirname "$0")"
TOOLS=../3DtransImages
RUNS=/home/victor/hosts/hpc/PhD/runs
RENDER_FLAGS="--mode --thr -2 --expand --left 0.15"

# name  downstream fade start  eigenmode file
while read -r n fs vtu; do
	[ "${RENDER:-0}" = 1 ] && pvbatch $TOOLS/render.py "$RUNS/$vtu" "pngs/$n.png" $RENDER_FLAGS
	python3 $TOOLS/gen_tex.py "pngs/$n.png" "pdfs/$n.tex" --fade-start "$fs"
	(cd pdfs && pdflatex -interaction=batchmode "$n.tex" >/dev/null && rm -f "$n.aux" "$n.log")
	echo "$n"
done <<END
d4w19_2dmode 0.5 incGapRe1000/arnoldiLinearSolver/d4_w19/NumSteps1e3/mesh_eig_0.vtu
d4w19_3dmode 0.2 quasi3DGapRe1000inc/arnoldiLinearSolver/multipleModes/d4_w19_Lz20_k32/mesh_eig_0.vtu
d4w19_3dFloquetmode 0.5 quasi3DGapRe1000inc/floquetSolver/multipleModes/d4_w19_Lz20_k32/NumSteps10e3/mesh_eig_0.vtu
END
