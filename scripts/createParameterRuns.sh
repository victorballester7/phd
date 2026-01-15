#!/bin/bash

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

SCRIPTS_DIR=$HOME/Desktop/PhD/scripts

function readInput {
  # prompt a message if there are less than 2 arguments
  if [ "$#" -lt 2 ]; then 
    echo -e "${RED}Usage: $0 <parameter_to_change> <values>${RESET}"
    echo -e "${YELLOW}For example: $0 mode_num 1 2 3 4 5${RESET}"
    exit 1
  fi

  parameter=$1
  parameterlabel=$parameter
  parameter_values=("${@:2}") # all the values after the 3rd argument
  parameter_valueslabel=("${@:2}")
  parent_dir=$(basename "$(dirname "$(realpath $session_file)")")

  if [[ "$parameter" == "Lz" || "$parameter" == "LZ" ]]; then
      parameter="LZ"
      parameterlabel="Lz"
      new_values=()
      new_valueslabel=()
      for v in "${parameter_values[@]}"; do
          new_values+=("${v}*depthGap")
          new_valueslabel+=("${v}d")
      done

      # Reassign arrays correctly
      parameter_values=("${new_values[@]}")
      parameter_valueslabel=("${new_valueslabel[@]}")
  fi

  echo -e "${CYAN}Input file: $session_file${RESET}"
  echo -e "${CYAN}Mesh file: $mesh_file${RESET}"
  echo -e "${CYAN}Parameter to change: $parameter${RESET}"
  echo -e "${CYAN}Values: $parameter_values${RESET}"
  echo -e "${CYAN}Values labels: $parameter_valueslabel${RESET}"
  echo -e "${CYAN}Parent directory: $parent_dir${RESET}"
  echo ""

  # Check if the input files exist (mesh file and session file)
  if [[ ! -f "$mesh_file" || ! -f "$session_file" ]]; then
    echo -e "${RED}Error opening the input files.${RESET}"
    exit 1
  fi

  # check if there exists such a parameter in the session file
  if grep -q "<p>[[:space:]]*${parameter}[[:space:]]*=.*</p>" "$session_file"; then
    pCapital=0
  elif grep -q "<P>[[:space:]]*${parameter}[[:space:]]*=.*</P>" "$session_file"; then
    pCapital=1
  else
    echo -e "${RED}Parameter $parameter not found in the session file.${RESET}"
    exit 1
  fi
}

function generateFolders {
  local n=${#parameter_values[@]}

  # Loop over the length of parameter_values
  for ((i=0; i<n; i++)); do
    p="${parameter_values[$i]}"
    plabel="${parameter_valueslabel[$i]}"    # Create directory for the current p
    dir_name="${parameterlabel}${plabel}"
    mkdir -p "$dir_name"
    cp "$mesh_file" "${dir_name}/${mesh_file}"
    
    # Copy the input file and modify it
    output_file="${dir_name}/${session_file}"
    
    if [ $pCapital -eq 1 ]; then
      sed "s/<P>[[:space:]]*${parameter}[[:space:]]*=.*<\/P>/<P> ${parameter} = ${p} <\/P>/g" "$session_file" > "$output_file"
    else
      sed "s/<p>[[:space:]]*${parameter}[[:space:]]*=.*<\/p>/<p> ${parameter} = ${p} <\/p>/g" "$session_file" > "$output_file"
    fi
    
    # Check if there are any .job files in the current directory
    if ls *.job 1> /dev/null 2>&1; then
      for job_file in *.job; do
        # if the job_file contains the pattern 'pbspro'
        if [[ "$job_file" == *pbspro* ]]; then
          # copy the job file to the new directory and change the job name
          cp "$job_file" "${dir_name}/${job_file}"
          sed -i "s/^#PBS -N .*/#PBS -N ${parent_dir}_${dir_name}/" "${dir_name}/${job_file}"
        else 
          # if the job_file contains the pattern 'slurm'
          cp "$job_file" "${dir_name}/${job_file}"
          sed -i "s/^#SBATCH --job-name=.*/#SBATCH --job-name=${parent_dir}_${dir_name}/" "${dir_name}/${job_file}"
        fi
      done
    else
      echo -e "${YELLOW}No .job files found in the current directory.${RESET}"
    fi

    echo -e "${GREEN}Generated file for ${parameter} $p in $dir_name${RESET}"
  done

  echo -e "${GREEN}All files have been generated.${RESET}"

}

source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh

readInput "$@"

generateFolders

