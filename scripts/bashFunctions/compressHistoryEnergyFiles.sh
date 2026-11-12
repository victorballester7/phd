# Keep the history/energy outputs of a run directory as ONE gzip file each (2026-11-10, replaces the .old scheme of
# updateHistoryEnergyFiles.sh). For every X in HISTORY_ENERGY_FILES:
#   X.gz   all previous runs (one header at the top, then data only; gzip members appended run after run)
#   X      the run in progress (plain text, written by Nektar++)
# compressHistoryEnergyFiles folds X.old (old scheme) and X into X.gz and removes them. Called by pbspro.sh before the
# solver (leftovers of a run that was killed before folding) and after it. Header lines (starting with # or a letter,
# e.g. "Time u_L2 ..." in EnergyError.err) are kept only at the top of X.gz. If the header of X differs from the one in
# X.gz (e.g. other history points), X.gz is kept as X.<date>.gz and a new X.gz is started, so readers never mix them.
# Readers: pp.fileManagement (X.gz + X), steadyBaseflows.py (histail).
HISTORY_ENERGY_FILES="HistoryPoints.his EnergyError.err EnergyFourierModes.mdl EnergyFile.mdl"

function _gz {
    if command -v pigz > /dev/null; then pigz -p ${NP:-4} -c; else gzip -c; fi
}

function _header {  # leading header block of a plain or gzip file
    if [[ "$1" == *.gz ]]; then zcat "$1" 2> /dev/null | awk '!/^[#A-Za-z]/{exit} {print}'
    else awk '!/^[#A-Za-z]/{exit} {print}' "$1"; fi
}

function compressOneHistoryFile {
    local f="$1" gz="$1.gz" src
    # an append to X.gz that was interrupted (walltime): cut X.gz back to its size before the append
    if [[ -f "$gz.size" ]]; then truncate -s "$(cat "$gz.size")" "$gz"; rm -f "$gz.size"; fi
    rm -f "$gz.part"
    for src in "$f.old" "$f"; do
        [[ -f "$src" ]] || continue
        if ! grep -qm1 '^[^#A-Za-z]' "$src"; then rm -f "$src"; continue; fi  # header only, no data
        if [[ -f "$gz" ]] && ! cmp -s <(_header "$gz") <(_header "$src"); then
            mv "$gz" "${f}.$(date +%Y%m%d_%H%M%S).gz"
            echo "compressHistoryEnergyFiles: header of $src differs from $gz; old one kept as ${f}.<date>.gz"
        fi
        if [[ ! -f "$gz" ]]; then
            # first member: header block once, later header lines (old .old files repeated them) dropped
            awk '/^[#A-Za-z]/{if (!d) print; next} {d=1; print}' "$src" | _gz > "$gz.part" &&
                mv "$gz.part" "$gz" && rm -f "$src" || { echo "compressHistoryEnergyFiles: failed on $src"; return 1; }
        else
            grep -v '^[#A-Za-z]' "$src" | _gz > "$gz.part" || { echo "compressHistoryEnergyFiles: failed on $src"; return 1; }
            stat -c %s "$gz" > "$gz.size"
            cat "$gz.part" >> "$gz" && rm -f "$gz.size" "$gz.part" "$src"
        fi
    done
}

function compressHistoryEnergyFiles {
    local x
    for x in $HISTORY_ENERGY_FILES; do
        compressOneHistoryFile "$x"
    done
}
