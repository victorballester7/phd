# Paths and solvers
INC_SOLVER=IncNavierStokesSolver
SCRIPTS_DIR=$HOME/Desktop/PhD/scripts

# Calculate the total number of processors
NP=$SLURM_NTASKS
JOB_ID=$SLURM_JOB_ID
slurmjob="$SLURM_SUBMIT_DIR/slurm.job"
datefile=$(date +"%Y%m%d_%H%M%S".date)

# source necessary functions
source $SCRIPTS_DIR/bashFunctions/computeTimeMax.sh
source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh
source $SCRIPTS_DIR/bashFunctions/uploadDateFile.sh
source $SCRIPTS_DIR/bashFunctions/updateHistoryEnergyFiles.sh
source $SCRIPTS_DIR/bashFunctions/runIncNS.sh

function getTime {
    # Extract the time line
    time_str=$(grep "^#SBATCH --time=" "$slurmjob" | awk '{print $2}')

    # Split into days and HH:MM:SS
    if [[ "$time_str" =~ ([0-9]+)-([0-9]+):([0-9]+):([0-9]+) ]]; then
        days="${BASH_REMATCH[1]}"
        hours="${BASH_REMATCH[2]}"
        minutes="${BASH_REMATCH[3]}"
        seconds="${BASH_REMATCH[4]}"
    elif [[ "$time_str" =~ ([0-9]+):([0-9]+):([0-9]+) ]]; then
        days=0
        hours="${BASH_REMATCH[1]}"
        minutes="${BASH_REMATCH[2]}"
        seconds="${BASH_REMATCH[3]}"
    else
        echo "Invalid time format"
        exit 1
    fi

}


function getOutputFiles {
    # Extract log and output file names from slurm.job
    log_file=$(grep "^#SBATCH --error=" "$slurmjob" | awk '{print $2}' | cut -d= -f2)
    output_file=$(grep "^#SBATCH --output=" "$slurmjob" | awk '{print $2}' | cut -d= -f2)
}

# create a file with the date and time of the job submission
function setup {
    getTime
    computeTimeMax
    getOutputFiles
    uploadDateFile
    updateHistoryEnergyFiles

    rm -f $output_file $log_file 
}

setup
runIncNS
