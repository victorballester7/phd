from pp.colors import colors
import subprocess


def run_cmd(cmd: str):
    """
    Executes a shell command and prints it in blue for better visibility.
    Parameters:
    - cmd: The command to be executed as a string.
    """
    print(colors.OKBLUE + f"Running command: {cmd}" + colors.ENDC)
    subprocess.run(cmd, shell=True, check=True)


def get_remote_dir(dir_local: str) -> str:
    """
    Converts a local directory path to a remote directory path.
    Assumes the remote base path needs to be change to "/rds/general/user/vb824/home/Desktop/PhD/runs" instead of "/home/victor/Desktop/PhD/src".
    Parameters:
    - dir_local: The local directory path as a string.
    Returns:
    - The corresponding remote directory path as a string.
    """
    remote_base = "/rds/general/user/vb824/home/Desktop/PhD/runs"
    local_base = "/home/victor/Desktop/PhD/src"
    if not dir_local.startswith(local_base):
        raise ValueError(
            colors.FAIL + f"Local directory must start with {local_base}" + colors.ENDC
        )
    dir_remote = dir_local.replace(local_base, remote_base)
    return dir_remote
