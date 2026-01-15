function getTime { 
    # get the time limit
    time_str=$(grep "^#PBS -l walltime=" "$pbsjob" | awk '{print $3}')
    if [[ "$time_str" =~ ([0-9]+):([0-9]+):([0-9]+) ]]; then
        hours="${BASH_REMATCH[1]}"
        minutes="${BASH_REMATCH[2]}"
        seconds="${BASH_REMATCH[3]}"
    else
        echo "Invalid time format"
        exit 1
    fi
}

function getOutputFiles {
    # Extract log and output file names from pbspro.job
    log_file=$(grep "^#PBS -e" "$pbsjob" | awk '{print $3}')
    output_file=$(grep "^#PBS -o" "$pbsjob" | awk '{print $3}')
    log_file_nek="log.txt"
    output_file_nek="output.txt"

}

function setup {
    getTime
    computeTimeMax
    getOutputFiles
    uploadDateFile
    updateHistoryEnergyFiles

    rm -f $output_file_nek $log_file_nek $output_file $log_file

    touch $output_file_nek $log_file_nek $output_file $log_file

    # Load required modules
    module load tools/eb-dev CMake/3.24.3-GCCcore-11.3.0 SCOTCH/7.0.1-gompi-2022a Boost/1.79.0-GCC-11.3.0 OpenBLAS/0.3.20-GCC-11.3.0 FlexiBLAS/3.2.0-GCC-11.3.0 FFTW/3.3.10-GCC-11.3.0 ScaLAPACK/2.2.0-gompi-2022a-fb HDF5/1.12.2-gompi-2022a VTK/9.2.2-foss-2022a
}

function run {
    if [[ $PBS_O_WORKDIR == *inc* ]]; then
        runNS inc
    else
        runNS com
    fi
}

# Paths and solvers
NEK_PATH=$HOME/nektar-v5.9.0
INC_SOLVER=$NEK_PATH/build/dist/bin/IncNavierStokesSolver 
COM_SOLVER=$NEK_PATH/build/dist/bin/CompressibleFlowSolver 
SCRIPTS_DIR=$HOME/Desktop/PhD/scripts

JOB_ID=$PBS_JOBID
pbsjob="$PBS_O_WORKDIR/pbspro.job"
datefile=$(date +"%Y%m%d_%H%M%S".date)

# Find the number of CPUs
NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')

source $SCRIPTS_DIR/bashFunctions/computeTimeMax.sh
source $SCRIPTS_DIR/bashFunctions/getMeshSessionFiles.sh
source $SCRIPTS_DIR/bashFunctions/uploadDateFile.sh
source $SCRIPTS_DIR/bashFunctions/updateHistoryEnergyFiles.sh
source $SCRIPTS_DIR/bashFunctions/runNS.sh

cd $PBS_O_WORKDIR

setup
run
