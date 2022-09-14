# Maass Form Trace Formula Code

This code is complementary to the paper [arxiv:2201.08760](https://arxiv.org/abs/2201.08760)

This code is for computing and verifying the Laplace and Hecke eigenvalues for squarefree level and trivial character. This code is in a rather unreadable state but I am creating this repository in case it may help anyone in the future.

The files are as follows:
- testfunc.c - This file contains all the functions needed to generate and compute the test function used in the computation. It also contains some other smaller functions needed for various computations.
- tracecomp.c - This file computes all the trace formula values for a range of Hecke operators T_n with -MAXN <= n <  MAXN. The inputs, which are hard-coded at the start of the file, are: NRPOCS = number of processors, MAXLEVEL = max level to compute, MAXN = maximum Hecke operator to compute, NUMSAMPLEPOINTS = number of sample points for elliptic term integrals, NUMTAYLORTERMS = number of terms in Taylor series approximation for elliptic integrals.
- tracecomp.sh - Shell script to organise the data from tracecomp.c into the individual levels.
- intetaylorterms.c - This files precomputes the integrals needed for the Taylor series approximation of the elliptic terms. 
- maasscomp.c - This file uses the data generated from tracecomp.c to compute and verify the Laplace and Hecke eigenvalues for the Maass forms we can compute. 
- completeness.c - This file computes the number of Maass forms we can prove completeness for and give a bound of the Laplace eigenvalue R for how far we can prove completeness to.
- L.gp - Pari code to compute the L function values for the hyperbolic terms.
- L.sh - Script to parallelise L.gp.
- Makefile - This is just to illustrate how to compile each of the files.

The only non-hard coded input to each of these files is the value d, which is the power to raise the test function.

tracecomp.c and maasscomp.c both rely on the prime factorisation library by Andrew Booker: https://github.com/arbooker/factor64 .
They also rely on the C-library [ARB](https://arblib.org/).

The code in this state is not really meant to be easily useable, and more to use as a guidance if you wish to implement this algorithm yourself.
