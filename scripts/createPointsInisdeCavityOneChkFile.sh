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

SCRIPTS_DIR=$HOME/Desktop/PhD/scripts

# generate equispaced data in a box of inside the cavity
# the box goes from x = 0 to x = width, and from y = depth/2 to y = -depth

dx=0.5
dy=0.05


datadir="data"
HOST="hpc"
# DIR="incNSboeingGapRe1000/directLinearSolver/blowingSuction/${CASE}"
# DIR="incNSboeingGapRe1000/baseflow/dns/d3_w26_stableIC" #DIR="flatSurfaceRe1000Ma0.6ComNS/dns"
# DIR="flatSurfaceRe1000Ma0.6ComNS/dns_shortDomain"
DIR="incGapRe1000/baseflow/dns/d4_w17_stableIC"
MESH_REMOTE="mesh.xml"
FLD_REMOTE="mesh_39_d2Udy2.fld" 
FLD_REMOTE_NOEXTENSION="${FLD_REMOTE%.*}"

# automatic variables
USER="vb824"
DIR_REMOTE_PRE="/rds/general/user/${USER}/home"
DIR_REMOTE="${DIR_REMOTE_PRE}/Desktop/PhD/runs/${DIR}"
DIR_LOCAL="${HOME}/Desktop/PhD/src/${DIR}"
mkdir -p "${DIR_LOCAL}/${datadir}"

generate_pts_file_for_chkfile() {
    local output_file="${datadir}/points${FLD_REMOTE_NOEXTENSION}"
    local file_path="${DIR_LOCAL}/${output_file}.pts"

    {
        echo '<?xml version="1.0" encoding="utf-8" ?>'
        echo '<NEKTAR>'
        echo '<POINTS DIM="2" FIELDS="">'

        for x in $(seq 0 $dx $width); do
            for y in $(seq $(echo "$depth / 2" | bc -l) -$dy $(echo "-$depth" | bc -l)); do
                printf "%.6f %.6f\n" "$x" "$y"
            done
        done

        echo '</POINTS>'
        echo '</NEKTAR>'
    } > "$file_path"

    echo -e "${GREEN}Generated combined .pts file: ${output_file}.pts${RESET}"
}

source ${SCRIPTS_DIR}/bashFunctions/getDepthANDWidth.sh

# Create .pts file for first CHKFILE and then cp all the others to the rest CHKFILES
getDepthANDWidth ${DIR_LOCAL}
generate_pts_file_for_chkfile

echo -e "${GREEN}All .pts files generated!${RESET}"


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
    
    FLD="${FLD_REMOTE}"
    output_file="${datadir}/points${FLD_REMOTE_NOEXTENSION}"
    echo "\${FLD}"
    echo "\${output_file}.pts"
    echo "Processing file ${FLD_REMOTE}..."
    rm -f \${output_file}.dat
    FieldConvert -m interppoints:fromxml=${MESH_REMOTE}:fromfld=\${FLD}:topts=\${output_file}.pts \${output_file}.dat
EOF

# Retrieve all .dat files
rsync -avz --progress "${USER}@${HOST}:${DIR_REMOTE}/${datadir}/" "${DIR_LOCAL}/${datadir}/"

echo -e "${GREEN}All interpolations completed and files retrieved!${RESET}"

