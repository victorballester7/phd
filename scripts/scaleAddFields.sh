#!/bin/bash

# Description: This script scales a field file by alpha, scales a second field file by beta, and adds them together to produce a new field file.

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"


if [ "$#" -ne 4 ]; then
    echo "Usage: $0 alpha beta first_field.fld second_field.fld"
    echo "Example: $0 0.3 0.7 field1.fld field2.fld"
    exit 1
fi

alpha=$1
beta=$2
first_field=$3
second_field=$4
output_field="initialCond.fld"

# Get the script's directory
DIR_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $DIR_SCRIPT/bashFunctions/getMeshSessionFiles.sh

echo -e "${CYAN}alpha: $alpha${RESET}"
echo -e "${CYAN}beta: $beta${RESET}"
echo -e "${CYAN}First field file (F1): $first_field${RESET}"
echo -e "${CYAN}Second field file (F2): $second_field${RESET}"
echo -e "${CYAN}Output field file (alpha * F1 + beta * F2): $output_field${RESET}"
echo -e "${CYAN}mesh file: $mesh_file${RESET}"
echo -e "${CYAN}session file: $session_file${RESET}"

FieldConvert -m scaleinputfld:scale="${alpha}" ${mesh_file} ${session_file} ${first_field} tmp.fld
FieldConvert -m addfld:fromfld=${second_field}:scale="${beta}" ${mesh_file} ${session_file} tmp.fld ${output_field}
rm tmp.fld

