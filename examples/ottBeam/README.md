
Deflection Example
==================

This simulation demonstrates how FDTD can be used for simulating scattering
of a beam represented by a vector spherical wave function expansion
by a spherical particle.
The simulation writes out a two text files for the beam axial slice
when the sphere is a two different positions in the beam, similar
to the middle two panels of Figure 3 from

> Isaac C. D. Lenton, Alexander B. Stilgoe, Halina Rubinsztein-Dunlop
> and Timo A. Nieminen., "Visual Guide to Optical Tweezers",
> *European Journal of Physics* **38**(3), 034009 (2017)
> https://doi.org/10.1088/1361-6404/aa6271

To run the simulation you will need an ab file with a column of
complex values representing the power in each vector spherical mode,
and a corresponding nm file with a column of integers specifying the
ab weights.
We use the convention used in version 1 of the optical tweezers toolbox
for the ab and nm vectors.
The default nm and ab file is for a Gaussian beam, almost any beam can
be substituted instead as long as you have an appropriate VSWF representation.
These files can be generated using a more recent version of the optical
tweezers toolbox, see the `generateBeam.m` script for an example.

