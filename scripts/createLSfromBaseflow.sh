#!/bin/bash

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

SCRIPTS_DIR=$HOME/Desktop/PhD/scripts

# For remote access
# HOST="typhoon"
HOST="hpc"

# do not edit
USER="vb824"
if [ "$HOST" = "typhoon" ]; then
    DIR_REMOTE_PRE="/home/${USER}"
else
    DIR_REMOTE_PRE="/rds/general/user/${USER}/home"
fi

function readInput {
  # prompt a message if there are less than 2 arguments
  if [ "$#" -ne 1 ]; then 
    echo -e "${RED}Usage: $0 <folderDNS>${RESET}"
    echo -e "${YELLOW}For example: $0 d4_w15 ${RESET}"
    echo -e "${YELLOW}             $0 d1.5_w45longInflow ${RESET}"
    exit 1
  fi

  folder=$1
  
  localDIR_LS=$(pwd)
  # localDIR_baseflow="/home/victor/Desktop/PhD/src/ffsRe1000inc/baseflow/dns/${folder}"
  # localDIR_baseflow="/home/victor/Desktop/PhD/src/bfsRe1000inc/baseflow/dns/${folder}"
  # localDIR_baseflow="/home/victor/Desktop/PhD/src/deepGapRe1000inc/baseflow/dns/${folder}"
  localDIR_baseflow="/home/victor/Desktop/PhD/src/incGapRe1000/baseflow/dns/${folder}"
  
  localDIRtmp="${localDIR_LS##*/Desktop/}"
  localDIRtmp="Desktop/${localDIRtmp/src/runs}"
  remoteDIR_LS="${DIR_REMOTE_PRE}/${localDIRtmp}/${folder}"

  localDIRtmp="${localDIR_baseflow##*/Desktop/}"
  localDIRtmp="Desktop/${localDIRtmp/src/runs}"
  remoteDIR_baseflow="${DIR_REMOTE_PRE}/${localDIRtmp}"


  cd "$localDIR_baseflow" || exit 1

  echo -e "${CYAN}Session file: $session_file${RESET}"
  echo -e "${CYAN}Mesh file: $mesh_file${RESET}"
  echo -e "${CYAN}Folder: $folder${RESET}"
  echo -e "${CYAN}Local LS directory: $localDIR_LS${RESET}"
  echo -e "${CYAN}Local baseflow directory: $localDIR_baseflow${RESET}"
  echo -e "${CYAN}Remote LS directory: $remoteDIR_LS${RESET}"
  echo -e "${CYAN}Remote baseflow directory: $remoteDIR_baseflow${RESET}"
  echo ""

  # Check if the input files exist (mesh file and session file)
  if [[ ! -f "$session_file" ]]; then
    echo -e "${RED}Error opening the input files.${RESET}"
    exit 1
  fi
}



function modifyFile {
  cd "$localDIR_LS" || exit 1
  mkdir -p "$folder"
  cd "$folder" || exit 1
  cp "../$session_file" .
  cp "$localDIR_baseflow/$geo_file" .
  cp "$localDIR_baseflow/$mesh_file" .
  cp "$localDIR_baseflow/pbspro.job" .
  cp "$localDIR_baseflow/slurm.job" .

  updateTimeStepLS

  ma=0

  if [[ "$PWD" =~ Ma([0-9]+(\.[0-9]+)?) ]]; then
    ma="${BASH_REMATCH[1]}"
  fi


  prepareRun.sh ${depth} ${width} ${ma}
  cd ..

  echo -e "${GREEN}Session file changed properly.${RESET}"

}

echo -e "${YELLOW} THIS SCRIPT SHOULD BE RUN FROM INSIDE THE LINEARSOLVER DIRECTORY (WHICH HAS TO CONTAIN A REFERENCE SESSION FILE. THE .GEO, .XML AND PBS FILES ARE TAKEN FROM THE BASEFLOW/DNS DIRERCTORY)\n ${RESET}"

source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh
source $SCRIPTS_DIR/bashFunctions/createBaseflowFile.sh
source $SCRIPTS_DIR/bashFunctions/updateTimeStepLS.sh
source $SCRIPTS_DIR/bashFunctions/getDepthANDWidth.sh

readInput "$@"
createBaseflowFile
depthANDwidth=$(getDepthANDWidth $folder)
depth=$(echo $depthANDwidth | awk '{print $1}')
width=$(echo $depthANDwidth | awk '{print $2}')
modifyFile

