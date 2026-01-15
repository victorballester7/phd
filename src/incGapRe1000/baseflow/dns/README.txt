This is the history of the Re1000_DNS.

1. We started running at fixed Reynolds (1000), d = 4 and varying the width of the gap w to see where the instability happens. So we ran for w = 10, 14, 18, 22, 26, 30 ('_old'-labelled folders).
2. I observed that the cases w = 10, 14 were stable and all the other ones look unstable. For the w = 14 case, we even ran the eigenvalue solver and we obtained all eigenvalues with negative real part.
3. With Yongyun we talk about how the different lengthsscales in the domain affect the solution. Thus, we did a sensitivity analysis that you can find inside this folder. In there we basically vary the length before and after the gap and the height of the domain. Have a look at the README file inside that folder to understand the coding of the folders.
4. Since the transition was between w=14 and w =18, we ran w=15,16,17 with a good mesh (a bigger one). 
5. Actually it was extremlly big the mesh. Enormous. So those examples are labelled as '_dense'. After a talk with my two supervisors, we agreed to reduce the size of the mesh in order to run fast now, and eventually in the future in 3d and compressible simulations.

