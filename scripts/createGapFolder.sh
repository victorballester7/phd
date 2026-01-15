#!/bin/bash

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
RESET="\e[0m"

SCRIPTS_DIR=$HOME/Desktop/PhD/scripts

function readInput {

  # Minimum required args: depth and width
  if [ "$#" -lt 2 ]; then
    echo -e "${RED}Usage: $0 <depth> <width> [--foldername NAME] [--mesh PATH]${RESET}"
    exit 1
  fi

  depth="$1"
  width="$2"
  shift 2   # Move past depth and width

  folder_name=""
  mesh_folder=""

  # Parse optional flags
  while [[ "$#" -gt 0 ]]; do
    case "$1" in
      --foldername)
        folder_name="$2"
        shift 2
        ;;
      --mesh)
        mesh_folder="$2"
        shift 2
        ;;
      *)
        echo -e "${RED}Unknown option: $1${RESET}"
        echo -e "${RED}Usage: $0 <depth> <width> [--foldername NAME] [--mesh PATH]${RESET}"
        exit 1
        ;;
    esac
  done

  # Default folder name if user didn’t provide one
  codename="d${depth}_w${width}"
  if [ -z "$folder_name" ]; then
    folder_name="$codename"
  fi

  # Default geo file unless overridden
  if [ -n "$mesh_folder" ]; then
    if [ ! -e "$mesh_folder" ]; then
      echo -e "${RED}Error: mesh path '$mesh_folder' does not exist.${RESET}"
      exit 1
    fi
    geo_file="$mesh_folder/mesh.geo"
  fi

  parent_dir=$(basename "$(dirname "$(realpath "$session_file")")")

  echo -e "${CYAN}Input file: $session_file${RESET}"
  echo -e "${CYAN}Geo file: $geo_file${RESET}"
  echo -e "${CYAN}Codename: $codename${RESET}"
  echo -e "${CYAN}Folder: $folder_name${RESET}"
  echo -e "${CYAN}Parent directory: $parent_dir${RESET}"
  if [ -n "$mesh_folder" ]; then
    echo -e "${CYAN}Custom mesh path: $mesh_folder${RESET}"
  fi
  echo ""

  # Check files exist
  if [[ ! -f "$session_file" ]]; then
    echo -e "${RED}Session file not found: $session_file${RESET}"
    exit 1
  fi

  if [[ ! -e "$geo_file" ]]; then
    echo -e "${RED}Geo file not found: $geo_file${RESET}"
    exit 1
  fi
}


function generateFolders {
  pbspro_file="pbspro.job"
  slurm_file="slurm.job"
  mkdir -p "$folder_name"

  cp $session_file "$folder_name/$session_file"
  cp $pbspro_file "$folder_name/$pbspro_file"
  cp $slurm_file "$folder_name/$slurm_file"
  
  # update the depth and width in the new geo file
  sed -e "s/^D = .*deltaStar;/D = ${depth} * deltaStar;/" \
    -e "s/^W = .*deltaStar;/W = ${width} * deltaStar;/" "$geo_file" > "${folder_name}/mesh.geo"

  if [[ -n "$mesh_folder" ]]; then
    cp "${mesh_folder}/mesh.msh" "${folder_name}/mesh.msh"
    cp "${mesh_folder}/mesh.xml" "${folder_name}/mesh.xml"
  fi

  cd "$folder_name"
  ma=0

  if [[ "$PWD" =~ Ma([0-9]+(\.[0-9]+)?) ]]; then
    ma="${BASH_REMATCH[1]}"
  fi

  prepareRun.sh ${depth} ${width} ${ma}
  cd ..

  echo -e "${GREEN}All files have been generated.${RESET}"

}

source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh

readInput "$@"
generateFolders

