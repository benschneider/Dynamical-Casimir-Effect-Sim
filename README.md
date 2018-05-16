# DCE simulation script

## Why: Two things: 

- 1. To verify our current understanding on the DCE  and compare it with the measured data.
- 2. To quantify the non-linearity in the system.
- 3. Approximate response of the SQUID to the second order.

![Critical Current in a SQUID](https://github.com/benschneider/Sim_DCE/files/2008820/squirIc.pdf)
![Hamiltonian of the System](https://github.com/benschneider/Sim_DCE/files/2008813/hamiltonian.pdf)
(J. R. Johansson et al.: Phys. Rev. A 82, 052509 (2010))

## What this script does:
	
- 1. Takes the SQUIDs non-linear phase response as a function of magnetic flux.
- 2. Creates a magnetic flux offset plus a sinusoidal change emulating a Flux signal.
- 3. Obtains phase response values corresponding to that Flux signal from (1)
- 4. (Similarly Obtains loss response values from the magnitude response.)
- 5. Calculates Frequency components of this phase response.
- 6. Multiplies these frequency components with a parabola expected by DCE.
- 7. Adds calculated losses due to frequencies above the superconducting gap frequency and the loss response values from the SQUID's magnitude response.

## How:
- [phase_res.py](https://github.com/benschneider/Sim_DCE/blob/master/phase_res.py), runs the simulation.
-  The Phase and Loss responses where obtained from measurement (presented in the paper in Fig.2c). The corresponding data files included are: **SQUID_phase.dat** and **SQUID_magnitude.dat**

- Both responses are then fitted using the python script.

![Flow_Chart](https://github.com/benschneider/Sim_DCE/files/2008815/Flow-hart.pdf)
