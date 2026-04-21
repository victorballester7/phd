#!/bin/bash

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

SCRIPTS_DIR=$HOME/Desktop/PhD/scripts
# HOST="typhoon"
HOST="hpc"
USER="vb824"
if [ "$HOST" = "typhoon" ]; then
    DIR_REMOTE_PRE="/home/${USER}"
else
    DIR_REMOTE_PRE="/rds/general/user/${USER}/home"
fi

function readInput {

  # Minimum required args: depth and width
  if [ "$#" -lt 8 ]; then
    echo -e "${RED}Usage: $0 <depth> <width> [--Lz LENGTHS (as function of <depth>)] [--NumSteps STEPS] [--mesh MESH_FOLDER]${RESET}"
    echo -e "${YELLOW}Example: $0 3 15 --Lz 1 2 3 --NumSteps 10e3 20e3 30e3 --mesh ../../incNSboeingGapRe1000/baseflow/dns/d3_w15${RESET}"
    exit 1
  fi

  depth="$1"
  width="$2"
  codename="d${depth}_w${width}"
  shift 2   # Move past depth and width

  Lz_values=()
  NumSteps_values=()
  mesh_folder=""

  # Parse optional flags
  while [[ "$#" -gt 0 ]]; do
    case "$1" in
      --Lz)
        shift
        while [[ "$#" -gt 0 && ! "$1" =~ ^-- ]]; do
          Lz_values+=("$1")
          shift
        done
        continue
        ;;
      --NumSteps)
        shift
        while [[ "$#" -gt 0 && ! "$1" =~ ^-- ]]; do
          NumSteps_values+=("$1")
          shift
        done
        continue
        ;;
      --mesh)
        mesh_folder="$2"
        shift 2
        continue
        ;;
      *)
        echo -e "${RED}Unknown option: $1${RESET}"
        echo -e "${RED}Usage: $0 <depth> <width> [--Lz LENGTHS (as function of d)] [--NumSteps STEPS] [--mesh MESH_FOLDER]${RESET}"
        echo -e "${YELLOW}Example: $0 3 15 --Lz 1 2 3 --NumSteps 10e3 20e3 30e3 --mesh ../../incNSboeingGapRe1000/baseflow/dns/d3_w15${RESET}"

        exit 1
        ;;
    esac
  done


  if [ -z "$mesh_folder" ]; then
    echo -e "${RED}Error: mesh folder not provided. Use --mesh MESH_FOLDER to specify it.${RESET}"
    exit 1
  fi

  localDIR_LS=$(pwd)
  localDIR_baseflow="${localDIR_LS}/${mesh_folder}"

  localDIRtmp="${localDIR_LS##*/Desktop/}"
  localDIRtmp="Desktop/${localDIRtmp/src/runs}"
  remoteDIR_LS="${DIR_REMOTE_PRE}/${localDIRtmp}/${codename}"

  localDIRtmp="${localDIR_baseflow##*/Desktop/}"
  localDIRtmp="Desktop/${localDIRtmp/src/runs}"
  remoteDIR_baseflow="${DIR_REMOTE_PRE}/${localDIRtmp}/"
  
  echo -e "${CYAN}Local LS directory: $localDIR_LS${RESET}"
  echo -e "${CYAN}Local baseflow directory: $localDIR_baseflow${RESET}"
  echo -e "${CYAN}Remote LS directory: $remoteDIR_LS${RESET}"
  echo -e "${CYAN}Remote baseflow directory: $remoteDIR_baseflow${RESET}"
}

source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh
source $SCRIPTS_DIR/bashFunctions/createBaseflowFile.sh
source $SCRIPTS_DIR/bashFunctions/updateTimeStepLS.sh

readInput "$@"
# createBaseflowFile
# updateTimeStepLS
# createGapFolder.sh $depth $width 

if [ $? -ne 0 ]; then
  echo -e "${RED}Error creating base gap folder.${RESET}"
  exit 1
fi
cd ${codename}
createParameterRuns.sh Lz ${Lz_values[@]}
if [ $? -ne 0 ]; then
  echo -e "${RED}Error creating Lz parameter runs.${RESET}"
  exit 1
fi
# for each Lz_values, cd into the folde "Lz{value}d" and run createParameterRuns.sh NumSteps ...
for Lz in "${Lz_values[@]}"; do
  cd "Lz${Lz}d"
  createParameterRuns.sh NumSteps ${NumSteps_values[@]}
  if [ $? -ne 0 ]; then
    echo -e "${RED}Error creating NumSteps parameter runs for Lz=${Lz}.${RESET}"
    cd ..
    exit 1
  fi
  cd ..
done
cd ..
echo -e "${GREEN}Quasi-3D gap folder structure created successfully!${RESET}"

