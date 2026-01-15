function onCancel {
    # if cat log_file is blank 
    echo "Job what killed intentionally by the user" >> $log_file
    python3 $SCRIPTS_DIR/jobFinish.py $JOB_ID
    exit 1
}

function trapCancel {
    trap "onCancel" SIGTERM
}
