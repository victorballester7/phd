#!/bin/bash

# Description: This script interpolates a field file (of the form .fld .chk, entered as an argument) that has been computed using an old mesh, to a new mesh. So it outputs the a field .fld file for the new mesh.

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 field_file.fld old_mesh_file.xml"
    echo "Example: $0 field_file.fld ../d4_w16.35/mesh.xml"
    exit 1
fi

field_file=$1
old_mesh_file=$2
new_field_file="initialCond.fld"

# Get the script's directory
DIR_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $DIR_SCRIPT/bashFunctions/getMeshSessionFiles.sh

echo -e "${CYAN}Old mesh file: $old_mesh_file${RESET}"
echo -e "${CYAN}Interpolating field file: $field_file${RESET}"
echo -e "${CYAN}New mesh file: $mesh_file${RESET}"
echo -e "${CYAN}Interpolated field file: $new_field_file${RESET}"

# Interpolate the field file

FieldConvert -m interpfield:fromxml=${old_mesh_file}:fromfld=${field_file} ${mesh_file} ${new_field_file}

