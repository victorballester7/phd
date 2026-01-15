#!/bin/bash

# Description
# 1. generate a .pts file with N points in a line for multiple x values
# 2. This file is then sent to the nodes, which uses fieldconvert of a .fld file to interpolate the field to the points in the .pts file.
# 3. The output of fieldconvert is copied back to the local machine.
# Usage: ./createPointsFileToInterpolateFieldTo.sh

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

constVar="x" #constVar="x"
# constValue=(-5.2 -4.333333333333334 -3.4666666666666677 -2.6000000000000014 -1.7333333333333352 -0.8666666666666689 -2.6645352591003757e-15 0.8666666666666636 1.7333333333333298 2.599999999999996 3.4666666666666623 4.333333333333328 5.199999999999995 6.066666666666662 6.933333333333327 7.799999999999993 8.66666666666666 9.533333333333328 10.399999999999991 11.266666666666659 12.133333333333326 12.999999999999993 13.866666666666656 14.733333333333324 15.59999999999999 16.466666666666658 17.333333333333325 18.19999999999999 19.066666666666656 19.933333333333323 20.799999999999986 21.666666666666654 22.53333333333332 23.399999999999988 24.266666666666655 25.13333333333332 25.999999999999986 26.86666666666665 27.733333333333317 28.599999999999984 29.46666666666665 30.333333333333318 31.199999999999985) #constValue=(-100 -75)
constValue=(200)
varValueMIN=0 #varValueMIN=0
varValueMAX=150 #varValueMAX=30
N=100 #N=800

datadir="data"
HOST="hpc"
# CASE="d3_w26" #CASE="d1.5_w50" 
CASE="flat"
# DIR="incNSboeingGapRe1000/directLinearSolver/blowingSuction/${CASE}"
# DIR="incNSboeingGapRe1000/baseflow/dns/d3_w26_stableIC" #DIR="flatSurfaceRe1000Ma0.6ComNS/dns"
# DIR="flatSurfaceRe1000Ma0.6ComNS/dns_shortDomain"
DIR="incGapRe1000/baseflow/dns/d3_w26_coarsecoarsecoarseMesh"
MESH_REMOTE="mesh.xml"
FLD_REMOTE="mesh_99.chk" #FLD_REMOTE="mesh_${CASE}_13.chk"
FLD_REMOTE_NOEXTENSION="${FLD_REMOTE%.*}"

# automatic variables
USER="vb824"
DIR_REMOTE_PRE="/rds/general/user/${USER}/home"
DIR_REMOTE="${DIR_REMOTE_PRE}/Desktop/PhD/runs/${DIR}"
DIR_LOCAL="${HOME}/Desktop/PhD/src/${DIR}"
mkdir -p "${DIR_LOCAL}/${datadir}"

generate_pts_file_for_chkfile() {
    local output_file="${datadir}/points${FLD_REMOTE_NOEXTENSION}_n${N}"
    local file_path="${DIR_LOCAL}/${output_file}.pts"

    {
        echo '<?xml version="1.0" encoding="utf-8" ?>'
        echo '<NEKTAR>'
        echo '<POINTS DIM="2" FIELDS="">'

        for varA in "${constValue[@]}"; do
            for ((i = 0; i < N; i++)); do
                step=$(echo "scale=6; $i / ($N - 1)" | bc)
                varB=$(echo "scale=6; $varValueMIN + ($varValueMAX - $varValueMIN) * $step * $step" | bc)
                if [ "$constVar" = "x" ]; then
                    x=$varA
                    y=$varB
                else
                    x=$varB
                    y=$varA
                fi
                printf "%.6f %.6f\n" "$x" "$y"
            done
            echo -e "${CYAN}Generated points for constValue=${varA}${RESET}" >&2
        done

        echo '</POINTS>'
        echo '</NEKTAR>'
    } > "$file_path"

    echo -e "${GREEN}Generated combined .pts file: ${output_file}.pts${RESET}"
}

# Create .pts file for first CHKFILE and then cp all the others to the rest CHKFILES
generate_pts_file_for_chkfile

echo -e "${GREEN}All .pts files generated for case ${CASE}!${RESET}"


# Transfer all .pts files
rsync -avz --progress "${DIR_LOCAL}/" "${USER}@${HOST}:${DIR_REMOTE}/"

# Create array string for remote use
constValue_str=$(printf "%s " "${constValue[@]}")

# Do all interpolations remotely in one SSH session
ssh "${USER}@${HOST}" /bin/bash << EOF
    source /etc/profile
    source ~/.bashrc  # Ensure modules are available
    
    # Try to load nektar++ module, continue if it fails
    module load nektar++ || echo "Warning: Could not load nektar++ module, continuing anyway..."
    
    cd "${DIR_REMOTE}"
    
    # Convert string back to array
    constValue_array=(${constValue_str})
    
    FLD="${FLD_REMOTE}"
    output_file="${datadir}/points${FLD_REMOTE_NOEXTENSION}_n${N}"
    echo "\${FLD}"
    echo "\${output_file}.pts"
    echo "Processing file ${FLD_REMOTE}..."
    rm -f \${output_file}.dat
    FieldConvert -m interppoints:fromxml=${MESH_REMOTE}:fromfld=\${FLD}:topts=\${output_file}.pts \${output_file}.dat
    echo "Completed CHKFILE ${FLD_REMOTE} (case: ${CASE})"
EOF

# Retrieve all .dat files
rsync -avz --progress "${USER}@${HOST}:${DIR_REMOTE}/${datadir}/" "${DIR_LOCAL}/${datadir}/"

echo -e "${GREEN}All interpolations completed and files retrieved for case ${CASE}!${RESET}"

