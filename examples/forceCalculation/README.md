Force Calculation
=================

This example shows how to calculate force on a spherical particle
in a plane wave beam.
There are a couple of methods for calculating (or estimating) force,
these typically include:

  * Integral over a surface surrounding the particle

  * Integral over the particle volume

Where the bounds of these volumes/surfaces are chosen and the property
that is integrated varies in the literature.
This example uses the Maxwell stress tensor over a cube containing
the particle inside the TFSF source layer for the volume and surface
integral.
This gives the flux of the total scattered field.
A second surface outside the TFSF source layer can also be used for
the TFSF integral to give the flux of the scattered-only field.

The accuracy of the results depends on various simulation parameters
including resolution and number of integration points.
For a sphere, we can fairly accurately calculate the scattering and
corresponding force using Mie theory.
This is a good comparison for the accuracy of the FDTD implementation.

This directory contains the simulation source file `sim.cpp`
as well as a Matlab script `mie.m` to calculate the force using the
[optical tweezers toolbox](https://github.com/ilent2/ott).

