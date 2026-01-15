function createBaseflowFile {
  # Capture SSH output into local variables
  ssh "${USER}@${HOST}" /bin/bash << EOF
    cd "${remoteDIR_baseflow}" || exit 1

    # Extract last chk number
    chk_number=\$(ls -td mesh_*.chk | head -1 |  sed -E 's/.*_([0-9]+)\.chk/\1/')

    mkdir -p "${remoteDIR_LS}"

    # Safely remove if baseflow.fld already exists (whether file or dir)
    rm -rf "${remoteDIR_LS}/baseflow.fld"
    
    echo -e "${CYAN}Copying mesh_\${chk_number}.chk to ${remoteDIR_LS}/baseflow.fld${RESET}"
    cp -r mesh_\${chk_number}.chk ${remoteDIR_LS}/baseflow.fld
EOF
}

