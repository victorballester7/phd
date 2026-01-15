function checkOldFileExists {
    local filename="$1"
    if [[ -f "${filename}.old" ]]; then
        grep -qv '^#' "$filename" >> "${filename}.old"
        rm -f "$filename"
    else
        mv "$filename" "${filename}.old"
    fi
}


function updateHistoryEnergyFiles {
    filenameHis="HistoryPoints.his"
    filenameE="EnergyError.err"
    filenameFE="EnergyFourierModes.mdl"


    checkOldFileExists "$filenameHis"
    checkOldFileExists "$filenameE"
    checkOldFileExists "$filenameFE"
}

