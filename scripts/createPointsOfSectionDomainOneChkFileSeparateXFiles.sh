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

# Get user input
# X=(15 20 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000)
# X=(-70 -50 -25 0 20 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000)

constVar="x"
# constValue=(-70 -50 -25 0 20 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000)
constValue=(-70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000)
varValueMIN=0
varValueMAX=150

N=600
datadir="data"

# For remote access
# HOST="typhoon"
HOST="hpc"

# DIR="incNSboeingGapRe1000/directLinearSolver/blowingSuction/d1.5_w45/wgn"
# MESH_REMOTE="mesh.xml"
# FLD_REMOTE="baseflow.fld"

CASE="d2_w36"
# CASE="flat"
# DIR="boeingGapRe1000IncNS/baseflow/dns/${CASE}"
DIR="incNSboeingGapRe1000/baseflow/dns/${CASE}/"
# DIR="flatSurfaceRe1000IncNS/directLinearSolver/blowingSuction/wgn"
MESH_REMOTE="mesh.xml"
CHKFILE="90"
CHKFILE_INIT=8
CHKFILE_FINAL=91
FLD_REMOTE="mesh_${CHKFILE}.chk"


# do not edit
USER="vb824"
if [ "$HOST" = "typhoon" ]; then
    DIR_REMOTE_PRE="/home/${USER}"
else
    DIR_REMOTE_PRE="/rds/general/user/${USER}/home"
fi
DIR_REMOTE="Desktop/PhD/runs/${DIR}"
DIR_REMOTE="${DIR_REMOTE_PRE}/${DIR_REMOTE}"
DIR_LOCAL="${HOME}/Desktop/PhD/src/${DIR}"
overwrite_all=false
skip_all=false

function genPointsFile {
    # Check if .dat file already exists
    if [ -f "${DIR_LOCAL}/${output_file}.dat" ]; then
        if [ "$skip_all" = true ]; then
            return 1
        fi
        if [ "$overwrite_all" = false ]; then
            while true; do
                echo -e "${YELLOW}File '${output_file}.dat' already exists. Would you like to overwrite it? (y/n/Y (yes all)/N (no all)) (default = N)${RESET}"
                read -n 1 -r choice
                echo ""
                choice=${choice:-N} # Default to 'A' if no input is given

                case $choice in
                    y)
                        break
                        ;;
                    n) 
                        return 1
                        ;;
                    Y)
                        overwrite_all=true;
                        break
                        ;;
                    N)
                        skip_all=true;
                        return 1
                        ;;
                    *)
                    echo -e "${RED}Invalid option. Introduce 'y' (yes), 'n' (no), 'Y' (yes all) or 'N' (no all).${RESET}"
                    ;;
                esac
            done
        fi
    fi

    echo -e "${YELLOW}Generating .pts file for ${constVar} = $varA...${RESET}"
    mkdir -p "${DIR_LOCAL}/${datadir}"
    
    # Create and write to the .pts file
    {
        echo '<?xml version="1.0" encoding="utf-8" ?>'
        echo '<NEKTAR>'
        echo '<POINTS DIM="2" FIELDS="">'
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
        echo '</POINTS>'
        echo '</NEKTAR>'
    } > "${DIR_LOCAL}/${output_file}.pts"
    
    echo -e "${GREEN}File '${output_file}'.pts generated successfully!${RESET}"

    return 0
}

function interpolateData {
    echo -e "${CYAN}Interpolating field data to points for ${constVar}=$varA...${RESET}"
    
    ssh "${USER}@${HOST}" /bin/bash << EOF
        cd "${DIR_REMOTE}"
        mkdir -p "${datadir}"
EOF
    
    scp -r "${DIR_LOCAL}/${output_file}.pts" "${USER}@${HOST}:${DIR_REMOTE}/$datadir"
    
    ssh "${USER}@${HOST}" /bin/bash << EOF
        source /etc/profile
        source ~/.bashrc  # Ensure modules are available
        module load nektar++ # will prompt an error when we are in HPC, but doesn't matter
        cd "${DIR_REMOTE}"
        rm -rf ${datadir}/*.dat # remove old data to prevent not overwriting
        FieldConvert -m interppoints:fromxml="${MESH_REMOTE}":fromfld="${FLD_REMOTE}":topts="${output_file}".pts ${output_file}.dat
EOF
    
    scp "${USER}@${HOST}:${DIR_REMOTE}/${output_file}.dat" "${DIR_LOCAL}/$datadir"
    echo -e "${GREEN}Field data interpolated to points for ${constVar} = $varA successfully!${RESET}"
}

# Create a directory to store all the data
mkdir -p "${DIR_LOCAL}/${datadir}"

# Process each X value
# constantchkfile
for varA in "${constValue[@]}"; do
    echo -e "${CYAN}Processing ${constVar} = $varA${RESET}"

    output_file="${datadir}/points${CHKFILE}_${constVar}${varA}_n${N}"

    genPointsFile 
    interpolateData
    echo -e "${GREEN}Completed processing for ${constVar} = ${varA}${RESET}"
    echo "----------------------------------------"
done

# varying chkfile
# for ((CHKFILE = ${CHKFILE_INIT}; CHKFILE <= ${CHKFILE_FINAL}; CHKFILE++)); do
#     FLD_REMOTE="mesh_${CHKFILE}.chk"
#     for varA in "${constValue[@]}"; do
#         echo -e "${CYAN}Processing ${constVar} = $varA${RESET}"

#         output_file="${datadir}/points${CHKFILE}_${constVar}${varA}_n${N}"

#         genPointsFile 
#         interpolateData
#         echo -e "${GREEN}Completed processing for ${constVar} = ${varA}${RESET}"
#         echo "----------------------------------------"
#     done
# done


echo -e "${GREEN}All x values have been processed successfully!${RESET}"
