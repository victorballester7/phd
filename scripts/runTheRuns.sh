#!/bin/bash

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'  # No Color

function determineJobType() {
    if command -v qstat &> /dev/null; then
        echo 'p'  # PBS Pro detected
    elif command -v scontrol &> /dev/null; then
        echo 's'  # Slurm detected
    else
        echo 'n'  # Neither detected
    fi
}


# Function to read user input
function inputRead() {
    # Prompt the user to select the execution mode: local or cluster
    # echo -e "${CYAN}Select execution mode: Local (n), PBS (p), or Slurm (s):${NC}"
    # read -n 1 -p "Run in nodes (n) or submit to cluster (p for pbspro or s for slurm)? [n/p/s]: " mode
    # echo -e " " # New line
    
    mode=$(determineJobType)

    if [[ "$mode" == "n" ]]; then
        echo -e "${RED}No job scheduler detected. Exiting.${NC}"
        exit 1
    fi

    echo -e "${CYAN}Execution mode: $mode${NC}"

    # Cluster mode (PBS/Slurm) should not expect mesh/session files, just folders
    if [[ "$#" -lt 1 ]]; then
        echo -e "${YELLOW}Usage: $0 <folder1> <folder2> ...${NC}"
        echo -e "${YELLOW}Example: $0 NumSteps500 NumSteps1000 NumSteps2000${NC}"
        exit 1
    fi
}

function getFolders() {
    for arg in "$@"; do
        if [[ "$arg" == *"*"* ]]; then
            folders+=( $(find . -maxdepth 1 -type d \( -name "$arg" \) | sort) )
        else
            folders+=("$arg")
        fi
    done
}

function runJobs() {
    # Iterate through folders
    folder0="$(pwd)"
    for folder in "${folders[@]}"; do
        if [ -d "$folder" ]; then
            if [[ "$mode" == "p" ]]; then
                if [[ -f "$folder/pbspro.job" ]]; then
                    echo -e "${GREEN}Submitting PBS job in $folder${NC}"
                    cd "$folder"
                    qsub $PBSJOB
                    cd "$folder0"
                else
                    echo -e "${RED}No file pbspro.job found in $folder${NC}"
                fi
            else  # Slurm mode
                if [[ -f "$folder/slurm.job" ]]; then
                    echo -e "${GREEN}Submitting Slurm job in $folder${NC}"
                    cd "$folder"
                    sbatch $SLURMJOB
                    cd "$folder0"
                else
                    echo -e "${RED}No file slurm.job found in $folder${NC}"
                fi
            fi
        else
            echo -e "${RED}$folder is not a directory${NC}"
        fi
    done
}

SLURMJOB="slurm.job"
PBSJOB="pbspro.job"

inputRead "$@"

getFolders "$@"

runJobs
