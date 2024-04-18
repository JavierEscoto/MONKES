# MONKES
 Code for calculating neoclassical transport in stellarators.

The MONoenergetic Kinetic Equation Solver (MONKES) is a freely available, open-source code that calculates neoclassical transport in stellarators by solving a radially local, monoenergetic drift-kinetic equation. MONKES main feature is its speed of calculation: The combination of spectral discretization in spatial and velocity coordinates with block sparsity allows MONKES to compute monoenergetic coefficients at low collisionality, in a single core, in approximately one minute. MONKES is sufficiently fast to be integrated into stellarator optimization codes for direct optimization of the bootstrap current (and/or neoclassical transport in general) and to be included in predictive transport suites.

Details about the drift-kinetic equation and the method used to solve it can be found at [Escoto, Velasco, Calvo, Landreman and Parra, Nuclear Fusion (2024)](https://iopscience.iop.org/article/10.1088/1741-4326/ad3fc9).

# License and usage of the code
In order to obtain permission to use MONKES and to receive an updated user manual and working examples, contact CIEMAT at fjavier.escoto@ciemat.es
