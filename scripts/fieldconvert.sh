#!/bin/bash

# Description: This script converts the output (checkpoint) files (*.chk or *.fld) to .vtu format 
# in order to visualize the results in Paraview or similar.

# Different kinds of usage: Run from the directory where the .chk or .fld files are stored:
# $directory_of_scripts/fieldconvert.sh 01 02 03  # Converts mesh_01.chk, mesh_02.chk, mesh_03.chk
# $directory_of_scripts/fieldconvert.sh -m vorticity myfield.fld
# $directory_of_scripts/fieldconvert.sh -m vorticity -m CFL 48

# IMPORTANT: The mesh file should be of the form mesh*.xml.

# Colors
RED='\e[31m'
GREEN='\e[32m'
YELLOW='\e[33m'
CYAN='\e[36m'
NC='\e[0m' # No Color

SCRIPTS_DIR=$HOME/Desktop/PhD/scripts

function getFiles (){
    files=()
    options=()

    # if the first argument is '-m' or '-n', then read the second argument and store it to later pass it to FieldConvert
    while [[ "$1" == "-m" || "$1" == "-n" ]]; do
        options+=("$1" "$2")
        shift 2
    done

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

function convertFiles () {
    for file in "${files[@]}"; do
        if [[ ! -d "${file}" && ! -f "${file}" ]]; then
            echo -e "${YELLOW}File ${file} does not exist. Skipping...${NC}"
            continue
        fi

        name=$(basename $(basename "${file}" .chk) .fld)
        output_file="${name}.vtu"

        if [ -f "${output_file}" ]; then
            if [ "$overwrite_all" = false ]; then
                while true; do
                    echo -ne "${YELLOW}File ${output_file} already exists. Do you want to overwrite it? (y/n/A): ${NC}"
                    read -n 1 -r choice
                    echo "" # new line
                    choice=${choice:-A}  # if the user presses enter, the default option is 'A'
                    case "$choice" in
                        y|Y) break ;;  # overwrite
                        n|N) echo -e "${GREEN}Omitting ${output_file}.${NC}"; continue 2;;
                        a|A) overwrite_all=true; break ;;  # overwrite all
                        *) echo -e "${RED}Invalid option. Introduce 'y' (yes), 'n' (no) o 'a' (all).${NC}" ;;
                    esac
                done
            fi
        fi

        echo -e "${GREEN}Converting ${file} to ${output_file}...${NC}"
        FieldConvert -f ${options[@]} "$mesh_file" "$session_file" "${file}" "${output_file}"
    done
}

source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh

# Directory to search (change '.' to the target directory if needed)
directory=$(pwd)

overwrite_all=false  # flag to overwrite all files

echo -e "${CYAN}Mesh file: ${mesh_file}${NC}"
echo -e "${CYAN}Session file: ${session_file}${NC}"

# get the files to convert
getFiles "$@"

# convert the files
convertFiles

echo -e "${CYAN}Conversion completed.${NC}"

