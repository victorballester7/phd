function runNS {
    # Run the job start script
    python3 $SCRIPTS_DIR/jobStart.py $JOB_ID

    if [ "$1" == "inc" ]; then
        SOLVER="$INC_SOLVER -v" 
    else
        SOLVER="$COM_SOLVER"
    fi
    
    # check if variable nforuier is set
    if [ -n "$nfourier" ]; then
        SOLVER="$SOLVER --npz $nfourier"
    fi

    # if NP > 1, use mpirun
    if [ $NP -ge 1 ]; then
        SOLVER="mpirun --timeout $TIMEMAX -np $NP $SOLVER"
    fi

    echo "Running solver with command: $SOLVER $mesh_file $session_file"

    # Run the solver
    $SOLVER $mesh_file $session_file > $output_file_nek 2> $log_file_nek
    
    python3 $SCRIPTS_DIR/jobFinish.py $JOB_ID

}
