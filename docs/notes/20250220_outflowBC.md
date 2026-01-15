# Outflow boundary conditions

The outflow boundary conditions implemented as

<N VAR='u' USERDEFINED="HOutflow" VALUE="0" />
<N VAR='v' USERDEFINED="HOutflow" VALUE="0" />
<N VAR='p' USERDEFINED="HOutflow" VALUE="0" />

is giving errors in the preassure and v-component of the velocity field at the outflow boundary for the baseflow and also for the EV solver. Initially I thought there were only problems in the EV solver but no! So from now on, we use normal Neumann and Dirichlet boundary conditions for the outflow boundary.
