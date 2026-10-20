#!/bin/bash
# Point InitialConditions in <folder>/session.xml to the last .chk file of the run on the remote host.
# Compressible: the analytic (Blasius) IC block is kept as SpongeRefFlow, the reference flow of the sponge.
# Incompressible: the analytic IC block is commented out.

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
  if [ "$#" -ne 1 ]; then
    echo -e "${RED}Usage: $0 <folder>${RESET}"
    echo -e "${YELLOW}For example: $0 d4_w16.5 ${RESET}"
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

  [[ "$localDIR" == *inc* ]] && incNS=1 || incNS=0
  [[ "$localDIR" == *3d* || "$localDIR" == *3D* ]] && is3D=1 || is3D=0

  echo -e "${CYAN}Session file: $session_file${RESET}"
  echo -e "${CYAN}isIncNS: $incNS${RESET}"
  echo -e "${CYAN}Folder: $folder${RESET}"
  echo -e "${CYAN}Local directory: $localDIR${RESET}"
  echo -e "${CYAN}Remote directory: $remoteDIR${RESET}"
  echo ""

  if [[ ! -f "$session_file" ]]; then
    echo -e "${RED}Error: $session_file not found.${RESET}"
    exit 1
  fi
}

function getChkFile {
  chk_number=$(ssh "${USER}@${HOST}" "cd '${remoteDIR}/$folder' && ls -td mesh_*.chk | head -1" | sed -E 's/.*_([0-9]+)\.chk/\1/')

  if [[ -z "$chk_number" ]]; then
    echo -e "${RED}Error: no mesh_*.chk found in ${remoteDIR}/$folder${RESET}"
    exit 1
  fi

  new_chk_file="mesh_${chk_number}.chk"
  echo -e "${CYAN}Last chk file: $new_chk_file${RESET}"
}

function make_ic_block() {
    if [[ "$incNS" -eq 1 ]]; then
        [[ "$is3D" -eq 1 ]] && vars="u,v,w,p" || vars="u,v,p"
    else
        [[ "$is3D" -eq 1 ]] && vars="rho,rhou,rhov,rhow,E" || vars="rho,rhou,rhov,E"
    fi
    printf "<FUNCTION NAME=\"InitialConditions\">\n<F VAR=\"%s\" FILE=\"initialCond.fld\" />\n</FUNCTION>\n" "$vars"
}

function restore_commented_analytic_ic() {
    # Older versions of this script commented out the analytic (Blasius) IC block, which also killed
    # the sponge (REFFLOW = InitialConditions then pointed to the restart file). Uncomment such a
    # block and rename it to SpongeRefFlow so it keeps acting as the sponge reference flow.
    tmp="$(mktemp)"
    awk '
        function flush(   n, i, line) {
            n = split(block, lines, /\n/)
            for (i = 1; i < n; i++) {
                line = lines[i]
                if (is_analytic) {
                    sub(/^[[:space:]]*<!--[[:space:]]?/, "", line)
                    sub(/[[:space:]]?-->[[:space:]]*$/, "", line)
                    if (i == 1) sub(/NAME="InitialConditions"/, "NAME=\"SpongeRefFlow\"", line)
                }
                print line
            }
            inside = 0; is_analytic = 0; block = ""
        }

        /^[[:space:]]*<!--[[:space:]]*<FUNCTION NAME="InitialConditions">/ {
            inside = 1; block = $0 "\n"; next
        }
        inside {
            block = block $0 "\n"
            if ($0 ~ /<E VAR=/) is_analytic = 1
            if ($0 ~ /<\/FUNCTION>/) flush()
            next
        }
        { print }
        END { if (inside) flush() }
    ' "$session_file" > "$tmp" && mv "$tmp" "$session_file"
}

function rename_analytic_ic() {
    # Rename the live analytic IC block (one containing <E VAR=...>) to SpongeRefFlow
    perl -0777 -pi -e 's#^(\s*)<FUNCTION NAME="InitialConditions">(?=(?:(?!</FUNCTION>).)*?<E VAR=)#${1}<FUNCTION NAME="SpongeRefFlow">#ms' "$session_file"

    if grep -q 'NAME="SpongeRefFlow"' "$session_file"; then
        sed -i 's#<REFFLOW> *InitialConditions *</REFFLOW>#<REFFLOW> SpongeRefFlow </REFFLOW>#' "$session_file"
    fi
}

function comment_analytic_ic() {
    # Incompressible runs have no sponge: comment out the live analytic (Blasius) IC block, also undoing
    # an earlier rename to SpongeRefFlow. Blocks reading from a file (<F ... FILE=...>) are left alone.
    tmp="$(mktemp)"
    awk '
        function flush(   n, i, line) {
            n = split(block, lines, /\n/)
            for (i = 1; i < n; i++) {
                line = lines[i]
                if (is_analytic) {
                    if (i == 1) sub(/NAME="SpongeRefFlow"/, "NAME=\"InitialConditions\"", line)
                    line = "<!-- " line " -->"
                }
                print line
            }
            inside = 0; is_analytic = 0; block = ""
        }

        /^[[:space:]]*<FUNCTION NAME="(InitialConditions|SpongeRefFlow)">/ {
            inside = 1; block = $0 "\n"; next
        }
        inside {
            block = block $0 "\n"
            if ($0 ~ /<E VAR=/) has_expr = 1
            if ($0 ~ /<F[^>]*FILE[[:space:]]*=/) has_file = 1
            if ($0 ~ /<\/FUNCTION>/) { is_analytic = has_expr && !has_file; has_expr = has_file = 0; flush() }
            next
        }
        { print }
        END { if (inside) flush() }
    ' "$session_file" > "$tmp" && mv "$tmp" "$session_file"
}

function mesh_ic_state() {
    # prints: 0 = no file IC block
    #         1 = exists and uncommented
    #         2 = exists but commented
    awk '
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

function insert_mesh_ic() {
    tmp="$(mktemp)"
    awk -v new_block="$(make_ic_block | sed 's/"/\\"/g')" '
        /<FUNCTION NAME="(InitialConditions|SpongeRefFlow)">/ && !done {
            print new_block
            done=1
        }
        {print}
    ' "$session_file" > "$tmp" && mv "$tmp" "$session_file"
}

function uncomment_mesh_ic() {
    # Strip the comment wrappers of the InitialConditions block that reads from a file
    tmp="$(mktemp)"
    awk '
        /<FUNCTION NAME="InitialConditions">/ { inside = 1; is_mesh = 0; block = $0 "\n"; next }
        inside {
            block = block $0 "\n"
            if ($0 ~ /<F[^>]*FILE[[:space:]]*=[[:space:]]*".*"/) is_mesh = 1
            if ($0 ~ /<\/FUNCTION>/) {
                if (is_mesh) {
                    n = split(block, lines, /\n/)
                    for (i = 1; i < n; i++) {
                        line = lines[i]
                        gsub(/<!--[[:space:]]*/, "", line)
                        gsub(/[[:space:]]*-->/, "", line)
                        print line
                    }
                } else {
                    printf "%s", block
                }
                inside = 0
            }
            next
        }
        { print }
    ' "$session_file" > "$tmp" && mv "$tmp" "$session_file"
}

function substitute_ic_file() {
    perl -0777 -pi -e 's#(<FUNCTION NAME="InitialConditions">\s*<F VAR="[^"]*" FILE=")[^"]*(" />\s*</FUNCTION>)#${1}'"$new_chk_file"'${2}#s' "$session_file"
}

function updateInitialConditions() {
    if [[ "$incNS" -eq 1 ]]; then
        comment_analytic_ic
    else
        restore_commented_analytic_ic
        rename_analytic_ic
    fi

    case "$(mesh_ic_state)" in
        0)
            echo "[INFO] No file IC found → creating it."
            insert_mesh_ic
            ;;
        2)
            echo "[INFO] File IC exists but is commented → uncommenting it."
            uncomment_mesh_ic
            ;;
    esac

    substitute_ic_file

    if grep -q "FILE=\"$new_chk_file\"" "$session_file"; then
      echo -e "${GREEN}InitialConditions now read from $new_chk_file.${RESET}"
    else
      echo -e "${RED}Error: could not set $new_chk_file as InitialConditions in $session_file.${RESET}"
      exit 1
    fi
}

function showDiff() {
    echo ""
    if cmp -s "$session_backup" "$session_file"; then
      echo -e "${YELLOW}No changes to $session_file.${RESET}"
      return
    fi
    echo -e "${CYAN}Changes to $session_file:${RESET}"
    diff -u --color=always --label "old/$session_file" --label "new/$session_file" \
        "$session_backup" "$session_file"
}

source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh

readInput "$@"
plotInfo
getChkFile

session_backup="$(mktemp)"
trap 'rm -f "$session_backup"' EXIT
cp "$session_file" "$session_backup"

updateInitialConditions
showDiff
