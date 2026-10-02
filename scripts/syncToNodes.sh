#!/bin/bash

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

# Define local and remote paths
# LOCAL_DIR is derived from this script's location (scripts/ lives inside the PhD repo),
# so it works both on the desktop and on Android/Termux, where $HOME is
# /data/data/com.termux/files/home and the repo may live elsewhere (e.g. ~/storage/shared/...).
# It can be overridden with PHD_DIR=/path/to/PhD syncToNodes.sh
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
LOCAL_DIR="${PHD_DIR:-$(dirname "$SCRIPT_DIR")}"
REMOTE_USER="vb824"  
REMOTE_DIR="~/Desktop/PhD"

# List of remote hosts
# HOSTS=("typhoon" "hpc")
HOSTS=("hpc" "typhoon")
SYNC_DIRS=("${LOCAL_DIR}/scripts")

# I set the src directory separately as it is synced to the runs directory on the remote hosts (which does not exist locally)
SRC_DIR="${LOCAL_DIR}/src"
RUNS_DIR="${REMOTE_DIR}/runs"

# list of directories to exclude from syncing following the rsync syntax
EXCLUSIONS=(
    "--exclude=*.venv*"
    "--exclude=*.egg-info*"
    "--exclude=*__pycache__*"
    "--exclude=*uv.lock*"
    "--exclude=*.dat"
    "--exclude=*.pts"
    "--exclude=*.chk"
    "--exclude=*.vtu"
    "--exclude=*.npz"
    "--exclude=*scripts/pp*" 
    "--exclude=*scripts/data*" 
    "--exclude=*scripts/images*" 
    "--exclude=*scripts/blasius*" 
    "--exclude=*spod*/*.dat"
)

if [ ! -d "$LOCAL_DIR/src" ]; then
    echo -e "${RED}Could not find PhD directory at $LOCAL_DIR (set PHD_DIR to override)${RESET}"
    exit 1
fi
echo -e "${CYAN}Local PhD directory: $LOCAL_DIR${RESET}"

# Loop through each host and sync the directories
for HOST in "${HOSTS[@]}"; do
    echo -e "${YELLOW}Syncing with $HOST...${RESET}"
    STATUS=0
    for DIR in "${SYNC_DIRS[@]}"; do
        rsync -avz --progress "${EXCLUSIONS[@]}" "$DIR" "$REMOTE_USER@$HOST:$REMOTE_DIR/" || STATUS=1
    done
    rsync -avz --progress "${EXCLUSIONS[@]}" "$SRC_DIR/" "$REMOTE_USER@$HOST:$RUNS_DIR" || STATUS=1

    if [ $STATUS -eq 0 ]; then
        echo -e "${GREEN}Sync with $HOST completed successfully!${RESET}"
    else
        echo -e "${RED}Sync with $HOST failed!${RESET}"
    fi
done

# CFDcourse26 is expected to sit next to the PhD directory; skip it if it is not there (e.g. on the phone)
CFDcourseDIR="$(dirname "$LOCAL_DIR")/CFDcourse26/src"
if [ -d "$CFDcourseDIR" ]; then
    rsync -avz --progress "${EXCLUSIONS[@]}" "$CFDcourseDIR/" "$REMOTE_USER@hpc:~/Desktop/CFDcourse26/runs"
else
    echo -e "${YELLOW}Skipping CFDcourse26 sync: $CFDcourseDIR not found${RESET}"
fi
