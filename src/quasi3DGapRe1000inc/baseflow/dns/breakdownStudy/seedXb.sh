#!/bin/bash
# <|u'|^2>_z and <w'^2>_z of selected checkpoints of a run, interpolated onto the
# pointsPO_n600 points (data/avg/pointsPO_n600_chk<n>.dat), to measure x_b vs seed.
# Usage: seedXb.sh <case> <chk> [<chk> ...]   (run from the dns directory)
source /etc/profile
module load tools/eb-dev CMake/3.24.3-GCCcore-11.3.0 SCOTCH/7.0.1-gompi-2022a Boost/1.79.0-GCC-11.3.0 OpenBLAS/0.3.20-GCC-11.3.0 FlexiBLAS/3.2.0-GCC-11.3.0 FFTW/3.3.10-GCC-11.3.0 ScaLAPACK/2.2.0-gompi-2022a-fb HDF5/1.12.2-gompi-2022a VTK/9.2.2-foss-2022a
export PATH=$HOME/nektar-v5.9.0/build/dist/bin:$PATH
c=$1; shift
cd $c
mkdir -p data/avg
[ -f data/pointsPO_n600.pts ] || cp ../${c%_seed*}/data/pointsPO_n600.pts data/
for n in "$@"; do
    out=data/avg/pointsPO_n600_chk$n.dat
    [ -s $out ] && continue
    t=tmp_$n
    FieldConvert -f -m addfld:fromfld=baseflow2DwithW.fld:scale=-1 mesh.xml session.xml mesh_$n.chk ${t}_pert.fld > /dev/null &&
    FieldConvert -f -m fieldfromstring:fieldstr="u^2+v^2+w^2":fieldname="uvw" mesh.xml session.xml ${t}_pert.fld ${t}_a.fld > /dev/null &&
    FieldConvert -f -m fieldfromstring:fieldstr="w^2":fieldname="ww" mesh.xml session.xml ${t}_a.fld ${t}_b.fld > /dev/null &&
    FieldConvert -f -m meanmode mesh.xml session.xml ${t}_b.fld ${t}_Zmean.fld > /dev/null &&
    FieldConvert -f -m interppoints:fromxml=mesh.xml:fromfld=${t}_Zmean.fld:topts=data/pointsPO_n600.pts $out > /dev/null ||
    echo "FAILED $c chk $n"
    rm -rf ${t}_pert.fld ${t}_a.fld ${t}_b.fld ${t}_Zmean.fld
done
echo "DONE $c"
