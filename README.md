# DCE simulation script

## Why: Two things: 

- 1. To verify our current understanding on the DCE  and compare it with the measured data.
- 2. To quantify the non-linearity in the system.
- 3. Approximate response of the SQUID to the second order.

![Critical Current in a SQUID](https://user-images.githubusercontent.com/4573907/40114209-50e3e604-590c-11e8-9805-73b1c692887b.png)
![Hamiltonian of the System](https://user-images.githubusercontent.com/4573907/40114212-544f9f36-590c-11e8-9134-163e7d6c2706.png)
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

![Flow_Chart](https://user-images.githubusercontent.com/4573907/40114197-4aeff86e-590c-11e8-917f-b8e18ca9b6f4.png)
