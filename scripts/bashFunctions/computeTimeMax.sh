# Remove leading zeros from hours, minutes, and seconds; like 01 -> 1
function computeTimeMax() {
    hours=$(echo $hours | sed 's/^0*//')
    minutes=$(echo $minutes | sed 's/^0*//')
    seconds=$(echo $seconds | sed 's/^0*//')

    # Convert days, hours, minutes, seconds to total seconds
    TIMEMAX=$((days * 86400 + hours * 3600 + minutes * 60 + seconds))

    # substract 30 seconds for safe termination
    TIMEMAX=$((TIMEMAX - 60))
}
