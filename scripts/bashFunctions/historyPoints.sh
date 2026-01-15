function historyPoints {
    # Description: This script inserts history points into a Nektar session .xml file 
    # based on the depth and width of the domain using float division.

    # Check argument count
    if [ "$#" -ne 3 ]; then
        echo -e "${RED}Usage: $0 sessionFile.xml depth width${RESET}"
        echo -e "${YELLOW}Example: $0 sessionFile.xml 4 15${RESET}"
        exit 1
    fi

    sessionFile=$1
    depth=$2
    width=$3

    # Validate the existence of the session file
    if [ ! -f "$sessionFile" ]; then
        echo -e "${RED}Error: File $sessionFile not found.${RESET}"
        exit 1
    fi

    geoFile=$(ls mesh*.geo 2>/dev/null | head -n 1)

    if [[ -z "$geoFile" ]]; then
        echo -e "${RED}No mesh*.geo file found!${RESET}"
        exit 1
    fi

    # Extract the numeric factor XXX (handles integer or float)
    lengthOutflowFactor=$(grep -oP 'lengthOutflow\s*=\s*\K[0-9]*\.?[0-9]+' "$geoFile")

    if [[ -z "$lengthOutflowFactor" ]]; then
        echo -e "${RED}Could not extract lengthOutflow factor from ${geoFile}!${RESET}"
        exit 1
    fi

    echo -e "${CYAN}Detected lengthOutflow factor: ${lengthOutflowFactor}${RESET}"

    echo -e "${CYAN}Processing session file: $sessionFile${RESET}"

    # Use sed to remove only the lines between <PARAM NAME="Points"> and </PARAM>
    sed '/<PARAM NAME="Points">/,/<\/PARAM>/ {
        /<PARAM NAME="Points">/b
        /<\/PARAM>/b
        d
    }' "$sessionFile" > "file.tmp"

    mv "file.tmp" "$sessionFile"
    echo -e "${GREEN}Previous history points removed.${RESET}"

    # Create the list of history points
    historyPoints=()

    # upstream points
    for offset in 100 75 50 35 25 15 10 5; do
        # historyPoints+=("-$offset 0.25 0")
        historyPoints+=("-$offset 0.5 0")
        # historyPoints+=("-$offset 0.75 0")
        # historyPoints+=("-$offset 1 0")
        # historyPoints+=("-$offset 1.25 0")
        # historyPoints+=("-$offset 1.5 0")
        # historyPoints+=("-$offset 1.75 0")
        # historyPoints+=("-$offset 2 0")
    done

    historyPoints+=("0 0.5 0")

    # points inside the gap
    step=5
    # we remove the -l option to use integer division
    numPointsInsideGap=$(echo "($width / $step)" | bc) 
    distanceBetweenPoints=$(echo "$width / ($numPointsInsideGap + 1)" | bc -l)

    echo -e "${CYAN}Adding $numPointsInsideGap history points inside the gap...${RESET}"

    for i in $(seq 1 $numPointsInsideGap); do
        x=$(echo "$i * $distanceBetweenPoints" | bc -l)
        historyPoints+=("$x 0.5 0")
        historyPoints+=("$x $(echo "-$depth / 2" | bc -l) 0")
    done

    # x=5
    # while (( $(echo "$x < $width - 0.0001" | bc -l) )); do
    #     historyPoints+=("$x 0.5 0")
    #     historyPoints+=("$x $(echo "-$depth / 2" | bc -l) 0")
    #     x=$(echo "$x + $step" | bc -l)
    # done


    for offset in 0 5 10 15 25 35 $(seq 50 25 $lengthOutflowFactor); do
        x=$(echo "$width + $offset" | bc -l)
        historyPoints+=("$x 0.5 0")
    done

    echo -e "${CYAN}Inserting history points into $sessionFile...${RESET}"

    # Insert the history points into the session file
    awk -v points="${historyPoints[*]}" '
        BEGIN { split(points, arr, " ") }
        /<PARAM NAME="Points">/ {
            print
            for (i = 1; i <= length(arr); i++) {
                if (i % 3 == 1) x = arr[i]
                else if (i % 3 == 2) y = arr[i]
                else print "    " x, y, arr[i]
            }
            next
        }
        { print }
    ' "$sessionFile" > temp.xml && mv temp.xml "$sessionFile"
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error inserting history points into $sessionFile.${RESET}"
        exit 1
    fi

    echo -e "${GREEN}History points successfully added to $sessionFile.${RESET}"
}
