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
    echo -e "${RED}Usage: $0 <folder>${RESET}"
    echo -e "${YELLOW}For example: $0 d4_w16.5 ${RESET}"
    echo -e "${YELLOW}             $0 NumSteps200e3 ${RESET}"
    exit 1
  fi

  folder=$1
  localDIR=$(pwd)

  cd "$folder" || {
    echo -e "${RED}Error: Folder '$folder' does not exist.${RESET}"
    exit 1
  }
}

function plotInfo {
  localDIRtmp="${localDIR##*/Desktop/}"
  localDIRtmp="Desktop/${localDIRtmp/src/runs}"
  remoteDIR="${DIR_REMOTE_PRE}/${localDIRtmp}"

  # check if the folder name contains "inc" (indicating incompressible Navier-Stokes) or not
      
  if [[ "$localDIR" == *inc* ]]; then
      incNS=1
  else
      incNS=0
  fi

  if [[ "$localDIR" == *3d* || "$localDIR" == *3D* ]]; then
      is3D=1
  else
      is3D=0
  fi

  echo -e "${CYAN}Session file: $session_file${RESET}"
  echo -e "${CYAN}isIncNS: $isIncNS${RESET}"
  echo -e "${CYAN}Mesh file: $mesh_file${RESET}"
  echo -e "${CYAN}Folder: $folder${RESET}"
  echo -e "${CYAN}Local directory: $localDIR${RESET}"
  echo -e "${CYAN}Remote directory: $remoteDIR${RESET}"
  echo ""

  # Check if the input files exist (mesh file and session file)
  if [[ ! -f "$session_file" ]]; then
    echo -e "${RED}Error opening the input files.${RESET}"
    exit 1
  fi
}


function make_ic_block() {
    if [[ "$incNS" -eq 1 ]]; then
        if [[ "$is3D" -eq 1 ]]; then
            printf "<FUNCTION NAME=\"InitialConditions\">\n<F VAR=\"u,v,w,p\" FILE=\"initialCond.fld\" />\n</FUNCTION>\n"
        else
            printf "<FUNCTION NAME=\"InitialConditions\">\n<F VAR=\"u,v,p\" FILE=\"initialCond.fld\" />\n</FUNCTION>\n"
        fi
    else
        if [[ "$is3D" -eq 1 ]]; then
            printf "<FUNCTION NAME=\"InitialConditions\">\n<F VAR=\"rho,rhou,rhov,rhow,E\" FILE=\"initialCond.fld\" />\n</FUNCTION>\n"
        else
            printf "<FUNCTION NAME=\"InitialConditions\">\n<F VAR=\"rho,rhou,rhov,E\" FILE=\"initialCond.fld\" />\n</FUNCTION>\n"
        fi
    fi
}

function mesh_ic_state() {
    # prints: 0 = no mesh IC block
    #         1 = exists and uncommented
    #         2 = exists but commented

    awk '
        # detect any InitialConditions block
        /<FUNCTION NAME="InitialConditions">/ {inside=1; next}

        inside && /<F VAR=.*FILE=".*"/ {
            print ($0 ~ /<!--/ ? 2 : 1)
            found=1
            exit
        }

        inside && /<\/FUNCTION>/ {inside=0}
        END { if (!found) print 0 }
    ' "$session_file"
}

function is_block_commented_or_noblock() {
    awk '
        BEGIN { inside=0; found=0 }

        # Detect start of IC block
        /<FUNCTION NAME="InitialConditions">/ {
            inside=1
            next
        }

        # Inside IC block → look for Blasius line (<E VAR=...>)
        inside && /<E VAR=/ {
            found=1
            if ($0 ~ /<!--/) print 1; else print 0
            exit
        }

        # End of any IC block → stop searching inside
        inside && /<\/FUNCTION>/ {
            inside=0
            next
        }

        END {
            if (found==0)
                print 1    # No Blasius IC found
        }
    ' "$session_file"
}

function comment_block() {
    awk '
        BEGIN {inside=0}
        /<FUNCTION NAME="InitialConditions">/ {
            inside=1
            print "<!-- " $0 " -->"
            next
        }
        inside && /<\/FUNCTION>/ {
            print " <!-- " $0 " -->"
            inside=0
            next
        }
        inside {
            print " <!-- " $0 " -->"
            next
        }
        {print}
    ' "$session_file" > "$tmp"
    mv "$tmp" "$session_file"
}

function uncomment_block() {
    awk '
        {
            line=$0
            sub(/^\s*<!--\s*/, "", line)
            sub(/\s*-->\s*$/, "", line)
            print line
        }
    ' "$session_file" > "$tmp"
    mv "$tmp" "$session_file"
}

function insert_mesh_ic() {
    tmp2="$(mktemp)"
    awk -v new_block="$(make_ic_block | sed 's/"/\\"/g')" '
        /<FUNCTION NAME="InitialConditions">/ && !done {
            print new_block
            done=1
        }
        {print}
    ' "$session_file" > "$tmp2"
    mv "$tmp2" "$session_file"
}

function substitute_ic_file() {
    perl -0777 -pi -e 's#(<FUNCTION NAME="InitialConditions">\s*<F VAR="[^"]*" FILE=")[^"]*(" />\s*</FUNCTION>)#${1}'"$new_chk_file"'${2}#s' "$session_file"
}

function uncomment_mesh_ic() {
    tmp="$(mktemp)"

    awk '
    BEGIN {
        inside = 0
        is_mesh = 0
        block = ""
    }

    # Start of any InitialConditions block: begin buffering
    /<FUNCTION NAME="InitialConditions">/ {
        inside = 1
        is_mesh = 0
        block = $0 "\n"
        next
    }

    # While inside a block, keep buffering every line and check for mesh F-line
    inside {
        block = block $0 "\n"
        # detect an F line that references a file (mesh IC)
        if ($0 ~ /<F[^>]*FILE[[:space:]]*=[[:space:]]*".*"/) {
            is_mesh = 1
        }
        # continue buffering until we hit the closing tag (handled below)
        if ($0 ~ /<\/FUNCTION>/) {
            # block is complete; process it
            if (is_mesh) {
                # split the buffered block into lines and strip ALL comment wrappers
                n = split(block, lines, /\n/)
                for (i=1; i<=n; i++) {
                    line = lines[i]
                    if (line == "") { print ""; continue }
                    # remove every occurrence of "<!--" and "-->" (and surrounding spaces)
                    # do this repeatedly until none remain (to handle nested wrappers)
                    while (match(line, /<!--/)) {
                        gsub(/<!--[[:space:]]*/, "", line)
                    }
                    while (match(line, /-->/)) {
                        gsub(/[[:space:]]*-->/, "", line)
                    }
                    # trim trailing and leading spaces introduced by removals (optional)
                    # but preserve original leading indentation if present
                    # We will remove only a single leading space if it is result of comment removal
                    # (leave other whitespace alone)
                    sub(/^[[:space:]]+/, "&", line)
                    print line
                }
            } else {
                # Not a mesh block — print it exactly as buffered
                printf "%s", block
            }
            # reset state
            inside = 0
            is_mesh = 0
            block = ""
            next
        }
        next
    }

    # Outside any block: print line unchanged
    { print }
    ' "$session_file" > "$tmp" && mv "$tmp" "$session_file"
}


function process_initial_conditions() {

    tmp="$(mktemp)"

    block_commented=$(is_block_commented_or_noblock)
    printf "block_commented: %s\n" "$block_commented"

    # 1. If the Blasius IC block is already commented → nothing to do
    if [[ "$block_commented" == "1" ]]; then
      echo "[OK] Blasius IC already commented. Nothing else to do."
      substitute_ic_file
      return
    fi

    # 2. Comment all 
    echo "[INFO] Commenting IC block..."
    comment_block

    mesh_state=$(mesh_ic_state)

    echo "mesh_state: $mesh_state"

    # 3. Handle the mesh block
    case "$mesh_state" in
        0)
            echo "[INFO] No mesh IC found → creating it."
            insert_mesh_ic
            ;;
        1)
            echo "[OK] Mesh IC already exists and is uncommented → nothing to do."
            ;;
        2)
            echo "[INFO] Mesh IC exists but is commented → uncommenting it."
            uncomment_mesh_ic
            ;;
    esac

    echo "mesh_state: $mesh_state"
    substitute_ic_file
}



function getTimeStepAndChkFile {
  # Capture SSH output into local variables
  read -r chk_number cfl1 cfl2 cfl3 cfl4 <<< "$(ssh "${USER}@${HOST}" /bin/bash << EOF
    cd "${remoteDIR}/$folder" || exit 1

    # Extract last chk number
    chk_number=\$(ls -td mesh_*.chk | head -1 |  sed -E 's/.*_([0-9]+)\.chk/\1/')

    # Extract last 4 CFL values
    mapfile -t cfls < <(grep 'CFL:' output.txt | tail -n 4 | awk '{print \$2}')

    echo "\$chk_number \${cfls[0]} \${cfls[1]} \${cfls[2]} \${cfls[3]}"
EOF
)"

  # Store CFLs in an array
  cfls=("$cfl1" "$cfl2" "$cfl3" "$cfl4")

  # Print results
  echo -e "${CYAN}Last chk file number: $chk_number${RESET}"
  echo -e "${CYAN}Last 4 CFL values:${RESET}"
  for i in "${!cfls[@]}"; do
    echo -e "${CYAN}CFL[-$((4 - i))]: ${cfls[$i]}${RESET}"
  done
}


function modifyFile {
  new_chk_file="mesh_${chk_number}.chk"

  update_just_restart_file=true
  update_timestep=true
  while true; do
    echo -e "${YELLOW}Would you like to update ONLY the restart file, or change as well the number of modes to 8-9 (which changes timestep as well) and remove SVV? (Y (only restart file)/n (change everything)) (default = y)${RESET}"
    read -n 1 -r choice
    echo ""
    choice=${choice:-Y} # Default to 'Y' if no input is given

    case $choice in
        y|Y)
            update_just_restart_file=true
            break
            ;;
        n|N)
            update_just_restart_file=false
            break
            ;;
        *)
        echo -e "${RED}Invalid option. Introduce 'y' (to change the restart file only), 'n' (to change everything).${RESET}"
        ;;
    esac
  done
  if [[ "$update_just_restart_file" == "true" ]]; then
    while true; do
      echo -e "${YELLOW}Would you like to update the timestep or keep the old one? (Y (yes, update the timestep)/n (keep the old timestep)) (default = y)${RESET}"
      read -n 1 -r choice
      echo ""
      choice=${choice:-Y} # Default to 'Y' if no input is given

      case $choice in
          y|Y)
              update_timestep=true
              break
              ;;
          n|N)
              update_timestep=false
              break
              ;;
          *)
          echo -e "${RED}Invalid option. Introduce 'y' (to change the restart file only), 'n' (to change everything).${RESET}"
          ;;
      esac
    done
  fi

  echo -e "${YELLOW}Update JUST the restart file: $update_just_restart_file${RESET}"
  echo -e "${YELLOW}Update timestep: $update_timestep${RESET}"

  process_initial_conditions

  # Use proper syntax for boolean logic
  if [[ "$update_just_restart_file" == "true" ]] || \
     (grep -q 'NUMMODES="9"' "$session_file" && grep -q 'NUMMODES="8"' "$session_file"); then
    mult="1"
  else
    mult="0.5625"
  fi

  # 2. Replace NUMMODES
  if [[ "$update_just_restart_file" != "true" ]]; then
    sed -i 's/NUMMODES="7"/NUMMODES="9"/' "$session_file"
    sed -i 's/NUMMODES="6"/NUMMODES="8"/' "$session_file"
  fi


  # 3. Extract timestep
  old_timestep=$(awk '/<P> *TimeStep *=/ {
    for(i=1;i<=NF;i++) {
      if ($i ~ /^[0-9.]+$/) {
        print $i; exit
      }
    }
  }' "$session_file")

  if [[ -z "$old_timestep" ]]; then
    echo -e "${RED}Could not extract TimeStep from $session_file.${RESET}"
    return 1
  fi

  # 4. Compute new timestep
  newCFL="0.3"
  new_timestep=$(echo "$old_timestep * $mult * $newCFL / ${cfls[0]}" | bc -l)


  if [[ "$update_timestep" == "false" ]]; then
    new_timestep=$old_timestep
  fi

  echo -e "${CYAN}Old TimeStep: $old_timestep${RESET}"
  echo -e "${CYAN}New TimeStep: $new_timestep${RESET}"

  # 5. Replace TimeStep line
  sed -i "s|<P> *TimeStep *= *${old_timestep} *<\/P>|<P> TimeStep = ${new_timestep} </P>|" "$session_file"

  # 6. Delete the line containing SpectralVanishingViscosity
  if [[ "$update_just_restart_file" != "true" ]]; then
    sed -i '/SpectralVanishingViscosity/d' "$session_file"
  fi
  
  echo -e "${GREEN}Session file changed properly.${RESET}"

}

source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh


readInput "$@"
plotInfo
getTimeStepAndChkFile
modifyFile
