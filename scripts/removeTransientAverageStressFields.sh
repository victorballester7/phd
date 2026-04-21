#!/bin/bash

# Description
# Finds the *_stress.fld average folders and subtracts the transient part.
# Specifically: selects the highest numbered folder and the first folder where
# <FinalTime> >= 3000, then averages between them with FieldConvert.

# Usage: ./removeTransientAverageStressFields.sh

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

HOST="hpc"
CASE="d2_w24/omega0.16" #CASE="d1.25_w70"
# CASE="flat"
# DIR="bfsRe1000inc/directLinearSolver/blowingSuction/${CASE}"
# DIR="deepGapRe1000inc/directLinearSolver/blowingSuction/${CASE}/"
DIR="incGapRe1000/directLinearSolver/omegaBlowSuct/${CASE}"
MESH_REMOTE="mesh.xml"
CHKFILE="_stress.fld"
COMBINEAVG="true" # whether to do combineAvg operation
MINTIME=5000 # minimum FinalTime to consider for combineAvg

# automatic variables
USER="vb824"
DIR_REMOTE_PRE="/rds/general/user/${USER}/home"
DIR_REMOTE="${DIR_REMOTE_PRE}/Desktop/PhD/runs/${DIR}"
DIR_LOCAL="${HOME}/Desktop/PhD/src/${DIR}"

mkdir -p "${DIR_LOCAL}"

echo -e "${CYAN}Connecting to ${HOST} to process stress averages...${RESET}"


# Do all interpolations remotely in one SSH session
ssh "${USER}@${HOST}" /bin/bash << EOF
    # set -e
    
    source /etc/profile
    source ~/.bashrc  # Ensure modules are available

    cd "${DIR_REMOTE}"

    echo -e "${CYAN} Scanning for *_stress.fld folders... ${RESET}"

    # list stress folders sorted numerically
    
    # sort by the number in the file name
    # stressFolders=(\$(ls -d *${CHKFILE} 2>/dev/null | sed -E 's/.*_([0-9]+)_stress\.fld/\1 \0/'  | sort -n -k1,1  | cut -d' ' -f2-))

    # sort by time creation
    # stressFolders=(\$(ls -dtr *${CHKFILE} 2>/dev/null)) # this sorts by modification time, and just looking at the file with vim (without modification) changes the modification time, which is not ideal.
    stressFolders=(\$(for d in *${CHKFILE}; do printf "%s %s\n" "\$(stat -c %Y "\$d/P0000000.fld" 2>/dev/null)" "\$d"; done | sort -n | awk '{print \$2}'))


    if [ \${#stressFolders[@]} -eq 0 ]; then
        echo -e "${RED} No stress average folders found! ${RESET}"
        exit 1
    fi

    stressFINAL=\${stressFolders[-1]}
    echo -e "${GREEN} Final stress folder: \${stressFINAL} ${RESET}"

    rm -rf mesh_avg.fld

    echo -e "${CYAN} CombineAvg is set to: ${COMBINEAVG} ${RESET}"
    
    if [ "${COMBINEAVG}" = "true" ]; then
        echo -e "${CYAN} Proceeding with combineAvg operation... ${RESET}"
        stressBEG=""
        for folder in "\${stressFolders[@]}"; do
            xmlfile="\${folder}/Info.xml"
            if [ ! -f "\$xmlfile" ]; then
                continue
            fi
            finalTime=\$(grep -oPm1 '(?<=<FinalTime>)[^<]+' "\$xmlfile" || echo 0)
            cmp=\$(awk -v t="\$finalTime" 'BEGIN{if (t >= ${MINTIME}) print 1; else print 0}')
            if [ "\$cmp" -eq 1 ]; then
                stressBEG="\$folder"
                break
            fi
        done

        if [ -z "\$stressBEG" ]; then
            echo -e "${RED} Could not find any folder with FinalTime >= ${MINTIME}! ${RESET}"
            exit 1
        fi

        # Patch NumberOfFieldDumps in Info.xml of stressBEG (make it negative if not already)
        xmlfile="\${stressBEG}/Info.xml"
        if [ -f "\$xmlfile" ]; then
            current=\$(grep -oPm1 '(?<=<NumberOfFieldDumps>)[^<]+' "\$xmlfile")
            echo -e "${YELLOW} Current NumberOfFieldDumps in \$xmlfile → \$current ${RESET}"
            if [[ "\$current" =~ ^[0-9]+$ ]]; then
                sed -i "s#<NumberOfFieldDumps>\$current</NumberOfFieldDumps>#<NumberOfFieldDumps>-\$current</NumberOfFieldDumps>#" "\$xmlfile"
                echo -e "${YELLOW} Patched NumberOfFieldDumps in \$xmlfile → -\$current ${RESET}"
            fi
        fi

        echo -e "${GREEN} Starting average from: \${stressBEG} to \${stressFINAL} ${RESET}"

        FieldConvert -m combineAvg:fromfld=\${stressBEG} mesh.xml \${stressFINAL} mesh_avg.fld

        echo -e "${CYAN} Completed combineAvg operation on ${CASE} ${RESET}"
    else
        echo -e "${YELLOW} Skipping combineAvg operation as COMBINEAVG is set to false. Copying final folder to _avg.fld instead... ${RESET}"
        cp -r "\${stressFINAL}" "mesh_avg.fld"
    fi

EOF
