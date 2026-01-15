#!/bin/bash

# Description: This script converts all the mesh files .msh (from gmsh) in the current directory to .xml format for later Nektar use.
# On top of that, it adds the history points to the session file (otherwise I forget to execute that script before submitting simulations).

# Usage: Run from the directory where the .msh files are stored:
# $directory_of_scripts/nekmesh.sh  

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

# execute historyPoints.sh
function hisPoints {
    # Execute the historyPoints.sh script
    echo -e "${CYAN}Executing historyPoints.sh...${RESET}"
    historyPoints $1 $2 $3
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error executing historyPoints.sh${RESET}"
        exit 1
    fi
}


# Function to convert .geo files to .msh
function convert_geo2msh {
    # Find .msh files in the directory
    files=($(find "${directory}" -maxdepth 1 -type f -name "*.geo" | sort))

    # If no files are found, exit with a message
    if [ ${#files[@]} -eq 0 ]; then
        echo -e "${YELLOW}No .geo files found in the directory.${RESET}"
        exit 1
    fi

    # Loop through the files and convert them
    for file in "${files[@]}"; do
        # Extract the basename of the file
        name=$(basename "${file}" .msh)

        echo -e "${CYAN}Processing: ${file}${RESET}"

        # Convert the file to .msh format
        gmsh -2 "${file}"

        if [ $? -eq 0 ]; then
            echo -e "${GREEN}File ${name}.msh created successfully.${RESET}"
        else
            echo -e "${RED}Error converting ${file} to msh.${RESET}"
        fi
    done
}

# Function to convert .msh files to .xml
function convert_msh2xml {
    # Find .msh files in the directory
    files=($(find "${directory}" -maxdepth 1 -type f -name "*.msh" | sort))

    # If no files are found, exit with a message
    if [ ${#files[@]} -eq 0 ]; then
        echo -e "${YELLOW}No .msh files found in the directory.${RESET}"
        exit 1
    fi

    # Loop through the files and convert them
    for file in "${files[@]}"; do
        # Extract the basename of the file
        name=$(basename "${file}" .msh)

        echo -e "${CYAN}Processing: ${file}${RESET}"
        
        # Convert the file to .xml format
        NekMesh "${file}" "${name}.xml"

        if [ $? -eq 0 ]; then
            echo -e "${GREEN}File ${name}.xml created successfully.${RESET}"
        else
            echo -e "${RED}Error converting ${name}.msh to XML.${RESET}"
        fi
    done
}

function editJobs {
  codename="d${depth%.}_w${width%.}"

  # remove leading 0 if ma is of the form 0.x
  if [[ "$ma" =~ ^0\. ]]; then
    ma=".${ma#0.}"
  fi
  jobname="d${depth%.}w${width%.}M${ma}"

  # Check if there are any .job files in the current directory
  if ls *.job 1> /dev/null 2>&1; then
    for job_file in *.job; do
      # If the job_file contains the pattern 'pbspro'
      if [[ "$job_file" == *pbspro* ]]; then
        sed -i "s/^#PBS -N .*/#PBS -N ${jobname}/" "${job_file}"
      # If the job_file contains the pattern 'slurm'
      else 
        sed -i "s/^#SBATCH --job-name=.*/#SBATCH --job-name=${jobname}/" "${job_file}"
      fi
    done
  else
    echo -e "${YELLOW}No .job files found in the current directory.${RESET}"
  fi   
}

# Get the script's directory
DIR_SCRIPT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source $DIR_SCRIPT/bashFunctions/getMeshSessionFiles.sh
source $DIR_SCRIPT/bashFunctions/historyPoints.sh

# Directory to search (current directory)
directory=$(pwd)

convert_geo2msh
convert_msh2xml

depth=$1
width=$2
ma=$3
sed -i -E "s#(<P> *depthGap *= *)[0-9.+-]+( *</P>)#\1${depth}\2#" "$session_file"
sed -i -E "s#(<P> *widthGap *= *)[0-9.+-]+( *</P>)#\1${width}\2#" "$session_file"

# execute historyPoints.sh
hisPoints $session_file $depth $width

# Edit the job files
editJobs
