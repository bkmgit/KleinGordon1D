# KleinGordon1D

Programs for solving the one dimensional cubic Klein Gordon equation

<img src="https://latex.codecogs.com/svg.latex?\Large&space;u_{tt}=\Delta{u}-u+u^3" title="\Large u_{tt}=\Delta u - u + u^3" />

The numerical approximations are compared to the exact traveling wave solution

<img src="https://latex.codecogs.com/svg.latex?\Large&space;u=\sqrt{2}\textup{sech}\left(\frac{x-ct}{\sqrt{1-c^2}}\right)" title="u=\sqrt{2}\textup{sech}\left(\frac{x-ct}{\sqrt{1-c^2}}\right)" />

for c = 0.5, t ∈ [0, 5] and x ∈ [−9π, 9π) with periodic boundary conditions.

## Using the programs

The programs have been tested on Linux systems and are written in [Fortran](https://wg5-fortran.org/). 
They use [GNU make](https://www.gnu.org/software/make/) or a compatible build system and require a 
Fortran 90 compiler. To compile the codes type make in the appropriate directory.

Example batch queue system submission scripts are included for use on a cluster. For the Fourier spectral 
codes, a Fast Fourier transform routine is required, at present either [FFTW 3](http://fftw.org/) or [MathKeisan FFT](http://mathkeisan.com/). Paths in the makefiles will need to be changed appropriately for your system.
