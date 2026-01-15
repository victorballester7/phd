#!/bin/bash
# Colors
RED='\e[31m'
GREEN='\e[32m'
YELLOW='\e[33m'
CYAN='\e[36m'
NC='\e[0m' # No Color


echo -ne "${YELLOW}THIS FILE NEEDS TO BE CALLED FROM THE DIRECTORY WHERE THE FILES OF THE EDGE STATE ARE LOCATED.${NC}\n"

function createICs {
    SOURCE_DIR_LOW="../d${depth}_w${width}_edgeStateLow"
    SOURCE_DIR_HIGH="../d${depth}_w${width}_edgeStateHigh"

    # Find the last index matching the pattern and compute +1
    last_indexconv=$(ls | grep -oP "convIC\K\d+" | sort -n | tail -1)
    last_indexdiv=$(ls | grep -oP "divIC\K\d+" | sort -n | tail -1)

    lowfile="convIC${last_indexconv}.fld"
    highfile="divIC${last_indexdiv}.fld"

    currentdir=$(pwd)

    # print energy of both current states


    echo -e "${CYAN}\nLow state Energy:${NC}"
    FieldConvert -m printfldnorms $mesh_file "convIC${last_indexconv}_pert.fld" out.stdout
    echo -e "${CYAN}High state Energy:${NC}"
    FieldConvert -m printfldnorms $mesh_file "divIC${last_indexdiv}_pert.fld" out.stdout

    
    scaleAddFields.sh 0.65 0.35 $lowfile $highfile
    mv initialCond.fld $SOURCE_DIR_LOW
    cd $SOURCE_DIR_LOW
    rm -r *.chk *.his *.err
    runTheRuns.sh .
    cd $currentdir

    scaleAddFields.sh 0.35 0.65 $lowfile $highfile
    mv initialCond.fld $SOURCE_DIR_HIGH
    cd $SOURCE_DIR_HIGH
    rm -r *.chk *.his *.err
    runTheRuns.sh .
    cd $currentdir
}

# Get the script's directory
DIR_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $DIR_SCRIPT/bashFunctions/getMeshSessionFiles.sh
source $DIR_SCRIPT/bashFunctions/getDepthANDWidth.sh

depthANDwidth=$(getDepthANDWidth $mesh_file)
depth=$(echo $depthANDwidth | awk '{print $1}')
width=$(echo $depthANDwidth | awk '{print $2}')

echo -e "${CYAN}Depth: $depth${NC}"
echo -e "${CYAN}Width: $width${NC}"

createICs
