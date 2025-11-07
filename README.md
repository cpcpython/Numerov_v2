# Numerov Package (v2)

**Numerov (v2) Contains a set of codes for variety of 1D Schrodinger equation problems.**

The Numerov package mainly consists of Python3 subroutines of Numerov method for 1D problems in time-independent quantum systems, and its application in molecular physics, including quantum simple harmonic oscillators, Particle in 1D Box problems, α particle decay, Franck-Condon vibronic transitions of diatomic molecules, quantum tunneling in various models (including quantum dots), Quantum wells, Hydrogen atom solutions, Calculation of radial functions from Pseudopotentias etc.

Originally Numerov is developed for calculating 1D Franck-Condon factors for various diatomic molecules using ab initio data points and later it is extended to other applications too.

*A key aspect of the Numerov package is that it includes a spline-based (linear/cubic) Numerov algorithm to support various types of potential representations, which are otherwise difficult to represent analytically.*

For example, if one has only few, but highly accurate ab inito potential energy values, using Spline-Numerov one can easily calculate highly accurate wavefunction and energies by avoiding analytically fitted potential functions (like Morse Potential). See our H2-H2+ Photoelectron Spectrum simulation which use Morse potential and Spline potentials (Refer: ).

We also added interactive-Matplotlib based graphical support to visualize data like wavefunctions. See Manual: Numerov_2025.pdf for more details.
completed / Tested codes are:

1.	P1Dbox - for the standard 1D 'particle in a box' problems. Here the particle is an electron.


This free software (GPLv3) is authored by: 

**Krishna Mohan G P [a], K R Aiswarya [b], Anusree P [b], Twinkle A R [b], Kiran Biju [a] and Roy K B [c]** 
 [a] Mar Baselios Engineering Colege and Technology (autonomous), Trivandrum, Kerala, India. 
 [b] Mar Ivanios College (autonomous), Trivandrum, Kerala India [c] SNGS Govt. College, Pattambi, Kerala, India.


