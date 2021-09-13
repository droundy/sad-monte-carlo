Monte Carlo
===========

A library for simulations in statistical mechanics using Monte Carlo algorithms.
This code  is architected for quickly testing changes to algorithms, and for
comparing different algorithms accross a variety of different physical systems.

This code features a wide variety of Monte Carlo algorithms, including ordinary
canonical Monte Carlo, Wang Landau (WL), $1/t$-Wang Landau ($1/t$-WL),
Stochastic Approximation Monte Carlo (SAMC), and Statistical Association with
Dynamic update factor (SAD), as well as a new Zeno's Monte Carlo algorithm.

This code also features a relatively broad set of relatively simple
systems upon which algorithms can be tested.  The code supports three simple
materials with periodic boundary conditions:  the Ising model, the hard-sphere
fluid, and the purely-repulsive Weeks-Chandler-Andersen fluid.  It supports
simulation of isolated clusters of Lennard-Jones atoms.  Finally, the code
supports a few artificial test systems for which we have analytic densities of
states.