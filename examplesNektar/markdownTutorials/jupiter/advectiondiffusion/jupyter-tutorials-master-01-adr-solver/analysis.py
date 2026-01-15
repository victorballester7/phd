import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt

############ h-type convergence tests ################################
meshList = [
    "mesh-1elmnt.xml",
    "mesh-4elmnt.xml",
    "mesh-25elmnt.xml",
    "mesh-36elmnt.xml",
    "mesh-49elmnt.xml",
    "mesh-81elmnt.xml",
    "mesh-100elmnt.xml",
]

# we will use Polynomial order 4 (NUMMODES=5) for h-type study. You can change this to any polynomial order you wish by adding the
# associate exp sessions at the end of command line. For example, if you want to use NUMMODES=8 modify the commands as follwos
#
# FieldConvert -m dof mesh.xml session.xml exp8.xml a.stdout
#
# ADRSolver mesh.xml session.xml exp8.xml
#
# note: the name of mesh and session in the above commands should be replaced with appropriate ones you are using.


# first we need to set up the array of Ndof (number of degree of freedom)
mypath = "../../../../nektar-v5.6.0/build/dist/bin/"
Ndof_htype = []
for mesh in meshList:
    cmd = mypath + "FieldConvert -m dof " + mesh + " errorAnalysisSession.xml a.stdout"
    res = subprocess.check_output(cmd, shell=True, text=True)
    ndof = res.split("Total number of DOF:", maxsplit=1)[-1].split(maxsplit=1)[0]
    Ndof_htype.append(ndof)

Ndof_htype = np.array(Ndof_htype, dtype=int)
print("Ndof array for h type analysis:")
print(Ndof_htype)
print("Dof count are finished")
print(" ")
print("Starting simulations for h type convergence")
# now run the simulations and store the errors

L2Err_htype = []
LinfErr_htype = []
# now solving the equations using several polynomial orders and save L2 and Linf errors
for mesh in meshList:
    cmd = mypath + "ADRSolver " + mesh + " errorAnalysisSession.xml "
    # grab the std out
    stdout = subprocess.check_output(cmd, shell=True, text=True)

    # get L2 error
    l2 = stdout.split("L 2 error (variable u) :", maxsplit=1)[-1].split(maxsplit=1)[0]
    L2Err_htype.append(l2)

    # get the Linf error
    lifn = stdout.split("L inf error (variable u) :", maxsplit=1)[-1].split(maxsplit=1)[
        0
    ]
    LinfErr_htype.append(l2)

# convert the results to numeric arrays
L2Err_htype = np.array(L2Err_htype, dtype=np.float64)
LinfErr_htype = np.array(LinfErr_htype, dtype=np.float64)


print("h type convergence simulatins are finished")
print("h type L2 errors: ")
print(L2Err_htype)
print("")
print("h type Linf erros:")
print(LinfErr_htype)
print("\n")
print("Starting The P type convergence analysis")
print("\n")
############ P-type convergence tests ################################
expOrders = [
    "exp2.xml",
    "exp3.xml",
    "exp4.xml",
    "exp5.xml",
    "exp6.xml",
    "exp7.xml",
    "exp8.xml",
    "exp9.xml",
    "exp10.xml",
    "exp11.xml",
]

# For this experiment, we will use the mesh with 25 elements. you can change the mesh to any mesh you desire

# first set up the array of Ndof (number of degree of freedom)
Ndof_ptype = []
for exp in expOrders:
    cmd = (
        mypath
        + "FieldConvert -m dof mesh-25elmnt.xml errorAnalysisSession.xml "
        + exp
        + " a.stdout"
    )
    res = subprocess.check_output(cmd, shell=True, text=True)
    ndof = res.split("Total number of DOF:", maxsplit=1)[-1].split(maxsplit=1)[0]
    Ndof_ptype.append(ndof)

Ndof_ptype = np.array(Ndof_ptype, dtype=int)
print("Ndof array for P type analysis:")
print(Ndof_ptype)
print("Dof count finished")
print(" ")
print("Starting simulations for P type convergence")

L2Err_ptype = []
LinfErr_ptype = []
# now solving the equations using several polynomial orders and save L2 and Linf errors
for exp in expOrders:
    cmd = mypath + "ADRSolver mesh-25elmnt.xml errorAnalysisSession.xml " + exp
    # grab the std out
    stdout = subprocess.check_output(cmd, shell=True, text=True)

    # get L2 error
    l2 = stdout.split("L 2 error (variable u) :", maxsplit=1)[-1].split(maxsplit=1)[0]
    L2Err_ptype.append(l2)

    # get the Linf error
    lifn = stdout.split("L inf error (variable u) :", maxsplit=1)[-1].split(maxsplit=1)[
        0
    ]
    LinfErr_ptype.append(l2)

# convert the results to numeric arrays
L2Err_ptype = np.array(L2Err_ptype, dtype=np.float64)
LinfErr_ptype = np.array(LinfErr_ptype, dtype=np.float64)

print("P type simulations are finished.")
print("P type L2 errors: ")
print(L2Err_ptype)
print("")
print("P type Linf erros:")
print(LinfErr_ptype)
print("\n\n")

# Creating figure and axes
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[7, 11])

# plot h_type L2 error
ax1.loglog(
    Ndof_htype, L2Err_htype, "-or", linewidth=2, markersize=10, label="h Convergence"
)
# Plot p_type L2 error
ax1.loglog(
    Ndof_ptype, L2Err_ptype, "-sb", linewidth=2, markersize=10, label="P Convergence"
)

ax1.set_title("L2 Error", fontsize=15)
ax1.set_xlabel("Number of Dof", fontsize=13)
ax1.set_ylabel("||E||2", fontsize=13)
ax1.legend()
ax1.grid(True, which="both", ls="-")


# plot h_type Linf error
ax2.loglog(
    Ndof_htype, LinfErr_htype, "-or", linewidth=2, markersize=10, label="h Convergence"
)
# Plot p_type Linf error
ax2.loglog(
    Ndof_ptype, LinfErr_ptype, "-sb", linewidth=2, markersize=10, label="P Convergence"
)

ax2.set_title("Linf Error", fontsize=15)
ax2.set_xlabel("Number of Dof", fontsize=13)
ax2.set_ylabel("||E||inf", fontsize=13)
ax2.legend()
ax2.grid(True, which="both", ls="-")

plt.tight_layout()
plt.show()

