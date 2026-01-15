# Kill process that are running even after scanceling a job from other nodes without direct access to them

## Problem

We have access to typhoon, from where we submit out slurm jobs. Now, sometimes (initially because the timeout command), programs keep running even after the job is canceled wioth scancel. I noticed this problem because the machines were running so slowly after that. In order to solve this just follow the steps below. We assume that those processes in the background are running on the node 'dragon'

```bash
ssh typhoon
srun --nodelist=dragon --ntasks=1 --pty top  # This will show the processes running on dragon
```

Now once here, select a process (hit enter) with your username, and then press 15 for SIGTERM or 9 for SIGKILL. I usually go for 9, since it's more aggressive.
