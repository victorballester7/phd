#!/bin/bash

# Description: This script produces the perturbation field u' such that u(x,y,z,t) = E_z(u(x,y,z,t)) + u'(x,y,z,t), where E_z(u) is the z-mean of u.

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

function getFiles (){
    files=()

    if [ "$#" -gt 0 ]; then
      for arg in "$@"; do
          if [[ "$arg" =~ ^[0-9]+$ ]]; then
              # we remove extension from mesh_file and add the argument (number) to it + .chk
              files+=( $(echo "${mesh_file}" | sed "s/\.[^.]*$//")_"$arg".chk )
          # if arg contains *, then we assume it is a pattern and we add all files that match it
          elif [[ "$arg" == *"*"* ]]; then
              files+=( $(find "${directory}" -maxdepth 1 -type d \( -name "$arg" \) | sort) )
          else  
              files+=("$arg")
          fi
      done
    else
      files=($(find "${directory}" -maxdepth 1 -type d \( -name "*.chk" -o -name "*.fld" \) | sort))
    fi

    if [ ${#files[@]} -eq 0 ]; then
      echo -e "${RED}No .chk or .fld files found in the directory.${NC}"
      exit 1
    fi
}

function run() {
    field_file=$1
    echo -e "${CYAN}Field file: $field_file${RESET}"
    field_name_no_ext="${field_file%.*}"
    Zmean_field_file="${field_name_no_ext}_Zmean.fld"
    perturbation_field_file="${field_name_no_ext}_pertZmean.fld"

    # Interpolate the field file
    echo -e "${CYAN}Computing Z-mean field...${RESET}"
    FieldConvert -m meanmode $mesh_file $session_file $field_file $Zmean_field_file

    echo -e "${CYAN}Computing perturbation field...${RESET}"
    FieldConvert -m addfld:fromfld=${Zmean_field_file}:scale=-1 $mesh_file $session_file $field_file $perturbation_field_file
}


if [ "$#" -eq 0 ]; then
    echo "Usage: $0 field_file.fld field_file_2.fld ..."
    echo "Example: $0 field_file.fld field_file_2.fld"
    exit 1
fi

# Get the script's directory
DIR_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $DIR_SCRIPT/bashFunctions/getMeshSessionFiles.sh

echo -e "${CYAN}Mesh file: $mesh_file${RESET}"
echo -e "${CYAN}Session file: $session_file${RESET}"

# get the files to convert
getFiles "$@"

for file in "${files[@]}"; do
    if [[ ! -d "${file}" && ! -f "${file}" ]]; then
        echo -e "${YELLOW}File ${file} does not exist. Skipping...${RESET}"
        continue
    fi

    echo -e "${GREEN}Processing ${file}...${RESET}"
    run "$file"
done

