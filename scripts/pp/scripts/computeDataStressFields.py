import os
import subprocess
from pp.fileManagement import editFile
import numpy as np
from pp.colors import colors
from pp.codenamesNfactor import code_names
import glob

def main():
    """
        for each folder in code_names, edit the files removeTransientAverage.sh and createPointsSectionDomainOneChkFile.sh and execute them from here
    """
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    bash_dir = os.path.join(script_dir, "../..")
    
    removeTransientScript = os.path.join(bash_dir, "removeTransientAverageStressFields.sh")
    createPointsScript = os.path.join(bash_dir, "createPointsOfSectionDomainOneChkFile.sh")

    blowSuction_dir = os.path.join(script_dir, "../../../src/incNSboeingGapRe1000/directLinearSolver/blowingSuction")

    combineavg = False 
    combineavg_str = "true" if combineavg else "false"
    line_startswithRTS = np.array(["CASE=", "COMBINEAVG="])
    line_startswithCPS = np.array(["CASE="])


    # just to do a quick check
    # code_names2 =np.array(["d3.75_w10",])

    for code in code_names:
        print(colors.HEADER + f"Processing case: {code}" + colors.ENDC)

        # check if there is already data in the folder

        case_dir = os.path.join(blowSuction_dir, code, "data", "*.dat")
        if glob.glob(case_dir):  # returns a list of matching files
            print(colors.WARNING + f"Data already exists for case: {code}, would you like to reprocess this case? (Y/n): " + colors.ENDC, end="")
            user_input = input().strip().lower()
            if user_input != 'y' and user_input != '':
                print(colors.WARNING + f"Skipping case: {code}" + colors.ENDC)
                continue
            else:
                print(colors.WARNING + f"Reprocessing case: {code}" + colors.ENDC)

        # edit files
        # change the line CASE="d4_w15" in createPointsScript and removeTransientScript
        replacement_lineRTS = np.array([f'CASE="{code}"', f'COMBINEAVG="{combineavg_str}"'])
        replacement_lineCPS = np.array([f'CASE="{code}"'])
        editFile(removeTransientScript, replacement_lineRTS, line_startswithRTS)
        editFile(createPointsScript, replacement_lineCPS, line_startswithCPS)
        
        print(colors.OKBLUE + f"Edited scripts for case: {code}" + colors.ENDC)

        print(colors.OKBLUE + f"Running remove transient script for case: {code}" + colors.ENDC)

        # run the files
        subprocess.run([removeTransientScript], cwd=bash_dir, check=True)

        print(colors.OKGREEN + f"Completed remove transient for case: {code}" + colors.ENDC)
        print(colors.OKBLUE + f"Running create points script for case: {code}" + colors.ENDC)

        subprocess.run([createPointsScript], cwd=bash_dir, check=True)

        print(colors.OKGREEN + f"Completed create points for case: {code}" + colors.ENDC)
        print(colors.OKGREEN + f"Completed processing for case: {code}" + colors.ENDC)
    


if __name__ == "__main__":
    main()

