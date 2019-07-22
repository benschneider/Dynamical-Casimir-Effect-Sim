# DCE simulation script

## Why:

- 1. To verify our current understanding on the DCE  and compare it with the measured data.
- 2. To quantify the non-linearity in the system.
- 3. Approximate response of a SQUID at the end of a Transmission line.

![Hamiltonian of the System](https://user-images.githubusercontent.com/4573907/40114212-544f9f36-590c-11e8-9134-163e7d6c2706.png)
(J. R. Johansson et al.: Phys. Rev. A 82, 052509 (2010))

## What this script does:

- 1. Takes the SQUIDs non-linear phase response as a function of magnetic flux.
- 2. Creates a magnetic flux offset plus a sinusoidal change emulating a Flux signal.
- 3. Obtains phase response values corresponding to that Flux signal from (1)
- 4. (Similarly Obtains loss response values from the magnitude response.)
- 5. Calculates Frequency components of this phase response.
- 6 changed (Multiplies these frequency components with a parabola expected by DCE.) 
Now these frequency components enter the following model:
```python
    part1 = F0*D*2j*pi #  F0*D is the output of these frequency components resulting in effective phase shifts
    part2 = np.sqrt(w*(wd-w))/(wp*Z0)  # DCE Parabola
    part3 = 1-(w**2)/(wp**2)+1j*Ldc*w/Z0  # SQUID reflection eqn r(w)(negative w)
    part4 = 1-(wd-w)**2/(wp**2)-1j*Ldc*(wd-w)/Z0  # SQUID reflection eqn r(w)(positive w)
    beta = part1*part2/(part3*part4)  # beta^2 is the photon rate
```
 
- 7. Adds calculated losses due to frequencies above the superconducting gap frequency and the loss response values from the SQUID's magnitude response.

![Flow_Chart](https://user-images.githubusercontent.com/4573907/40114197-4aeff86e-590c-11e8-917f-b8e18ca9b6f4.png)

## How:

- [start.py](https://github.com/benschneider/Sim_DCE/blob/master/start.py), runs the simulation.
-  The Phase and Loss responses where obtained from measurement (presented in the paper in Fig.2c). The corresponding data files included are: **SQUID_phase.dat** and **SQUID_magnitude.dat**
- Both responses are then fitted using the python script.

- Output MTX files generated by the script can then be visualized using [Spyview](http://nsweb.tn.tudelft.nl/~gsteele/spyview/).
