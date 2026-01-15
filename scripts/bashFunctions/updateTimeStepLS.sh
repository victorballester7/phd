function updateTimeStepLS() {
  # 3. Extract timestep
  timestep=$(awk '/<P> *TimeStep *=/ {
    for(i=1;i<=NF;i++) {
      if ($i ~ /^[0-9.]+$/) {
        print $i; exit
      }
    }
  }' "$localDIR_baseflow/$session_file")

  if [[ -z "$timestep" ]]; then
    echo -e "${RED}Could not extract TimeStep from $session_file${RESET}"
    return 1
  fi

  echo -e "${CYAN}TimeStep: $timestep${RESET}"

  # 5. Replace TimeStep line
  sed -i "s|<P> *TimeStep *= *[^<]*</P>|<P> TimeStep = ${timestep} </P>|" "$session_file"

  echo -e "${GREEN}Session file changed properly.${RESET}"

}
