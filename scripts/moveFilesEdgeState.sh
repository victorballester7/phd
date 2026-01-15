#!/bin/bash
# Colors
RED='\e[31m'
GREEN='\e[32m'
YELLOW='\e[33m'
CYAN='\e[36m'
NC='\e[0m' # No Color



function movefiles {
    while true; do
        echo -ne "${YELLOW}Would you like to copy from Low or High folder? (l/h): ${NC}"
        read -n 1 -r choiceLH
        echo "" # new line
        case "$choiceLH" in
            l|L) break ;;  # overwrite
            h|H) break;;
            *) echo -e "${RED}Invalid option. Introduce 'l' or 'h'.${NC}" ;;
        esac
    done

    while true; do
        echo -ne "${YELLOW}Has the initial condition been converged or diverged? (c/d): ${NC}"
        read -n 1 -r choiceCD
        echo "" # new line
        case "$choiceCD" in
            c|C) break ;;  # overwrite
            d|D) break;;
            *) echo -e "${RED}Invalid option. Introduce 'c' or 'd'.${NC}" ;;
        esac
    done


    if [ "$choiceLH" == "l" ] || [ "$choiceLH" == "L" ]; then
        SOURCE_DIR="../d${depth}_w${width}_edgeStateLow"
    elif [ "$choiceLH" == "h" ] || [ "$choiceLH" == "H" ]; then
        SOURCE_DIR="../d${depth}_w${width}_edgeStateHigh"
    fi

    if [ "$choiceCD" == "c" ] || [ "$choiceCD" == "C" ]; then
        PATTERN="convIC"
    elif [ "$choiceCD" == "d" ] || [ "$choiceCD" == "D" ]; then
        PATTERN="divIC" 
    fi

    # Find the last index matching the pattern and compute +1
    last_index=$(ls | grep -oP "${PATTERN}\K\d+" | sort -n | tail -1)
    if [ -z "$last_index" ]; then
        next_index=0
    else
        next_index=$((last_index + 1))
    fi

    prefix="${PATTERN}${next_index}"

    echo "Next index: $prefix"

    mv $SOURCE_DIR/initialCond.fld ${prefix}.fld
    mv $SOURCE_DIR/HistoryPoints.his ${prefix}_historyPoints.his
    mv $SOURCE_DIR/EnergyError.err ${prefix}_energyerror.err

    FieldConvert -m addfld:fromfld=baseflowSFD.fld/:scale=-1 $mesh_file ${prefix}.fld ${prefix}_pert.fld
}

echo -ne "${YELLOW}THIS FILE NEEDS TO BE CALLED FROM THE DIRECTORY WHERE THE FILES OF THE EDGE STATE ARE LOCATED.${NC}\n"

# Get the script's directory
DIR_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $DIR_SCRIPT/bashFunctions/getMeshSessionFiles.sh
source $DIR_SCRIPT/bashFunctions/getDepthANDWidth.sh

depthANDwidth=$(getDepthANDWidth $mesh_file)
depth=$(echo $depthANDwidth | awk '{print $1}')
width=$(echo $depthANDwidth | awk '{print $2}')

echo -e "${CYAN}Depth: $depth${NC}"
echo -e "${CYAN}Width: $width${NC}"

movefiles
