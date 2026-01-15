# get the depth and width of a file path
function getDepthANDWidth() {
    # the file is of the form: /path/to/file/mesh_d4_w15.5.geo
    # Extract the filename
    filename=$(basename $1)

    # Use regex to extract depth and width
    if [[ $filename =~ d([0-9.]+)_w([0-9.]+) ]]; then
        depth="${BASH_REMATCH[1]}"
        width="${BASH_REMATCH[2]}"
    else
        echo "0 0" # Return 0 0 if the pattern does not match
        return 0
    fi

    # removes trailing zeros, if they exist
    echo "${depth%.} ${width%.}"
}
