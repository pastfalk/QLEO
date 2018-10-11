QLEO - A full-f quasilinear simulation code 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a manual for the Fortran-90 code QLEO containing a short description of the program, an explanation of the input parameters, and some advice for the correct usage of the code. Note that this code version is adapted to the GNU Fortran Compiler and was tested with the compiler version gcc 6.1. The code was partly parallelized with OpenMP. For executing QLEO, first run 'make' and than execute './qsolve'.

General remarks
---------------
The QLEO ('Quasilinear Electromagnetic Oscillations') code is based on the quasilinear kinetic equation which self-consistently describes how a given velocity distribution function of a homogeneous gyrotropic plasma evolves in time in the presence of an unstable spectrum of parallel propagating, electromagnetic linear eigenmodes. The implemented formalism is based on the quasilinear kinetic equations given in, e.g., 'Methods in Nonlinear Plasma Theory' by R. C. Davidson (1972). The code allows for the inclusion of an arbitrary number of particle species where each species can be described either by an arbitrary gyrotropic velocity distribution function or as a static (bi-)Maxwellian background species. The velocity distribution function is provided to the code as a data set sampled on the 2d (vpara,vperp) velocity space and is advanced in time using an explicit Euler method.
The QLEO code is intended to enable a systematic study of the quasilinear stabilization of parallel propagating velocity space instabilities, assuming realistic velocity distribution functions.

You are welcome to use and distribute this code, and to adapt QLEO to your own purposes. If you publish some work which is based on results obtained with QLEO, please cite the corresponding paper Astfalk & Jenko, JGR 2018 and make sure that you mark any substantial changes in the code with respect to the original source code.


The velocity distribution
-------------------------
If you choose to include a particle species with a prescribed arbitrary gyrotropic distribution, switch the input parameter 'mode' to '1'. The code then reads the distribution from the file 'distribution/distributionX.dat' where 'X' is the index ('1', '2', '3', ...)  of the corresponding particle species where the numbering is according to all included species with arbitrary velocity distributions. The file is expected to sample the distribution F in parallel and perpendicular velocity space, v_para and v_perp,  with an equidistant grid in each direction. The data has to be arranged in three columns where column 1 = v_para, column 2 = v_perp, column 3 = F(vpara,vperp). The velocities are expected to be in ascending order and for each v_para, F is scanned for all v_perp, before proceeding with the next v_para. For reference, see the sample distribution file which is produced by 'distribution/print_bimax.py'.
The velocity components and the distribution values have to be normalized with respect to the Alfvén velocity of the first particle species. Furthermore, the velocity distribution has to be normalized with respect to the species density such that \int F(v_para,v_perp)*v_perp dv_perp dv_para = 1.0.


The program structure
---------------------
The QLEO code advances a given velocity distribution function and spectrum of eigenmodes in time by solving the quasilinear kinetic equations using an explicit Euler method. This requires (1) the evaluation of an integral over the considered spectrum of wavenumbers ('get_deltaf()') and, at each time step, (2) a computation of the instantaneous dispersion relation of the eigenmodes ('disp_rel()'). Furthermore, (3) an estimation of the velocity derivatives of the distribution function is required ('deriv_dist()').

(1) For computing the dispersion relation, the formalism is coupled to the LEOPARD code described in Astfalk & Jenko, JGR 2017. At each time step, the code evaluates the dispersion relation twice. First, it covers the full spectrum of eigenmodes requested by the user. Then, it adjusts the spectrum to the unstable modes only ('adapt_krange()') and recomputes the dispersion relation. 

(2) The integration is carried out piecewise-analytically after interpolating the frequency, growth rate, and magnetic energy spectrum with cubic splines. After interpolation, the denominator of the integrand turns into a cubic polynomial which may have roots in the considered wavenumber interval. The roots may introduce singularities in the integration interval and are computed with 'cubic_sol()'. The contributions of these poles to the integral have to be added accordingly.

(3) The velocity derivatives of the distribution function are computed based on local exponential fits. For each grid point, the distribution values of the two adjacent grid points are taken into account to fit an exponential to the velocity distribution. Then, The derivatives are applied to the exponential. This method is more accurate than standard central difference methods because the exponential better reflects the local shape of the velocity distribution.


The MPFUN package
-----------------
For evaluating the dispersion relation based on arbitrary velocity distributions, the computation of the parallel velocity integrals and the determination of the generalized hypergeometric functions required for the computation of the perpendicular velocity integrals sometimes demands higher accuracy than double or quadruple precision could provide. Therefore, the MPFUN2015 package by David H. Bailey was included which allows thread-safe Fortran computations with arbitrary precision.
The package was found to sometimes cause a breakdown of the code. The problem seems to be related to issues with type conversion but the origin of this has not been located yet. If you come across this issue and find a fix for it, feel free to share your solution.


The input parameters
--------------------

&wavenumber

kstart - The lower border of the chosen wavenumber interval.

kend   - The upper border of the chosen wavenumber interval.

nk     - This determines at how many points the code will evaluate the dispersion relation within the chosen wavenumber interval.


Note:
All wavenumbers are given in units of the inertial length of the first particle species.


&initial_guess

omega_r     - The initial guess for the real frequency from which the Muller method starts to iterate a root, omega(k), of the dispersion relation at the wavenumber k(1)=kstart.

omega_i     - The initial guess for the growth or damping rate from which the Muller method starts to iterate a root, omega(k), of the dispersion relation at the wavenumber k(1)=kstart.

increment_r - The frequency value by which the previously found root, omega(k), is incremented to provide the starting value for the next Muller iteration at the subsequent wavenumber, k(i+1)=k(i)+dk.

increment_i - The growth rate value by which the previously found root, omega(k), is incremented to provide the starting value for the next Muller iteration at the subsequent wavenumber, k(i+1)=k(i)+dk.


Note:
A proper initial guess is crucial for a successful root finding. If the guess lies too far away from the dispersion branch of interest, you may land on another branch. In general, the algorithm always converges to a certain root. It has to be figured out by the user, whether this is the root of interest.

The increments are only necessary for the initial guesses for the root finding at k(2) and k(3). For subsequent wavenumbers a quadratic polynomial approximation determines all following initial guesses. If dk is not too high, this usually works well.

Both frequencies and growth rates are always given in units of the gyro frequency of the first particle species.


&setup

Nspecies - The number of particle species the user wants to include.

theta 	 - The propagation angle of the waves, i.e. the angle between the wave vector k and the background magnetic field (which is aligned with the z-axis in the chosen coordinate system). So far, QLEO can only simulate a spectrum of parallel propagating eigenmodes, i.e. theta=0.

delta 	 - The squared ratio of gyro frequency and plasma frequency of the first particle species.

Bstart	 - The amplitude of the initial magnetic energy perturbation delta B_k ^2.

dt	 - The time step in units of the first particle species' gyro frequency.

tmax	 - The maximum time, up to which simulation is carried out, in units of the first particle species' gyro frequency.

sgn	 - The sign of the chosen wave polarization. For omega_r > 0, sgn=1 corresponds to left-hand and sgn=-1 corresponds to right-hand polarization.

restart	 - When restarting a simulation, set restart to the desired time step. Otherwise set restart=0.


Note:

The parallel and perpendicular wavenumbers are given as k_para=k*cos(theta) and k_perp=k*sin(theta).

Delta gives a measure for the magnetization of the plasma. Low delta corresponds to weak, high delta corresponds to strong magnetization.


&accuracy

rf_error    - The 'root finding error' gives the exit-condition for the Muller iteration. An error of 1.0d-2 or 1.0d-3 generally gives good results. But - of course - the choice depends on the accuracy requested by the user.

eps_error   - The 'epsilon error' gives the exit condition for the sum over the Bessel index n. Once the relative contribution of the computed dielectric tensor components for a given n gets smaller than the given eps_error, the code exits the loop.


Note:
If a solution seems fishy, play with these parameters and check whether the solution is numerically converged.

Choose the rf_error to be not too demanding, otherwise the algorithm may run into convergence problems.



&species

mode_in      - Choose '0' for a bi-Maxwellian plasma or '1' for a plasma with arbitrary gyrotropic velocity distribution. For mode_in=0, the particle species will be considered as a static Maxwellian background. Only for mode_in=1, the particle velocity distribution will be evolved in time.

q_in         - Charge of the particles in units of the charge of the first particle species.

mu_in 	     - Mass of the particles in units of the mass of the first particle species.

dens_in	     - Density of the particles in units of the density of the first particle species

drift_in     - This introduces a drift velocity to the bi-Maxwellian distribution (mode '0' only). The drift is normalized with respect to the Alfvén velocity.

sym_in       - Choose '1' if provided velocity distribution function (mode '1' only) is symetric with respect to the parallel velocity, choose '0' if it is asymmetric. Please note that sym_in = 0 has not been tested thoroughly yet. 

beta_para_in - Beta parameter parallel to the background magnetic field (mode '0' only).

beta_perp_in - Beta parameter perpendicular to the background magnetic field (mode '0' only).


Note:
If you need more than the two default particle species, just add additional parameter blocks below the two default &species blocks. The choice, which particle species is declared in the first &species block, is of major importance since the normalization of all output data depends on this choice. E.g., if you choose protons to be the first particle species, then all frequencies and growth rates will be given in units of the proton gyrofrequency and the wavenumbers will be in units of the proton inertial length.

When including particle species with arbitrary velocity distribution, the size of the provided velocity grid will significantly affect the performance of the code. Good accuracy at sufficiently fast run times was found for distribution grids with 200 points in v_para and 100 points in v_perp - but in general the performance is highly dependent on the detailed velocity space structure of the distribution and the considered dispersion branch.



The output files
----------------
anis_X.dat:
	column 1 - simulation time
	column 2 - total magnetic energy
	column 3 - parallel beta of species X
	column 4 - perpendicular beta of species X
	column 5 - maximum growth rate

omega.dat:
	column 1 - wavenumber
	column 2 - frequency
	column 3 - growth/damping rate

distribution/evolution/dist/X_Y.dat:
	column 1 - parallel velocity of particle species Y
	column 2 - perpendicular velocity of particle species Y
	column 3 - distribution value of particle species Y at time step X

distribution/evolution/delfdelpa/X_Y.dat:
	column 1 - parallel velocity of particle species Y
	column 2 - perpendicular velocity of particle species Y
	column 3 - derivative df / dv_pa  of particle species Y at time step X
	column 4 - derivative d^2 f / dv_para^2  of particle species Y at time step X
	column 5 - derivative d^2 f / dv_para dv_perp  of particle species Y at time step X

distribution/evolution/delfdelpe/X_Y.dat:
	column 1 - parallel velocity of particle species Y
	column 2 - perpendicular velocity of particle species Y
	column 3 - derivative df / dv_perp  of particle species Y at time step X
	column 4 - derivative d^2 f / dv_perp^2  of particle species Y at time step X

distribution/evolution/df/X_Y.dat
	column 1 - parallel velocity of particle species Y
	column 2 - perpendicular velocity of particle species Y
	column 3 - change of velocity distribution function df computed from kinetic quasilinear equation for particle species Y at time step X

distribution/evolution/T1/X_Y.dat
	column 1 - parallel velocity of particle species Y
	column 2 - solution of the wavenumber integration over term T11 = -0.25*dt*(q*mu)**2 * dB_k^2 * |omega|^2 /k^2 /(omega-k*vpa-sgn*mu*q) at time step X
	column 3 - solution of the wavenumber integration over term T12 = -0.25*dt*(q*mu)**2 * dB_k^2 * omega^* /k /(omega-k*vpa-sgn*mu*q) at time step X
	column 4 - solution of the wavenumber integration over term T13 = -0.25*dt*(q*mu)**2 * dB_k^2 * omega/k /(omega-k*vpa-sgn*mu*q) at time step X
	column 5 - solution of the wavenumber integration over term T14 = -0.25*dt*(q*mu)**2 * dB_k^2 /(omega-k*vpa-sgn*mu*q) at time step X

distribution/evolution/T2/X_Y.dat
	column 1 - parallel velocity of particle species Y
	column 2 - solution of the wavenumber integration over term T21 = -0.25*dt*(q*mu)**2 * dB_k^2 * omega /(omega-k*vpa-sgn*mu*q)^2 at time step X
	column 3 - solution of the wavenumber integration over term T22 = -0.25*dt*(q*mu)**2 * dB_k^2 * k /(omega-k*vpa-sgn*mu*q)^2 at time step X

distribution/evolution/omega/X.dat
	simulation time
	column 1 - wavenumber
	column 2 - magnetic energy
	column 3 - frequency
	column 4 - growth rate


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For more information on QLEO, see Astfalk & Jenko, JGR 2018. If you have further questions concerning the usage of the code or if you like to discuss some general issues, feel free to contact patrick.astfalk@ipp.mpg.de.
 QLEO
