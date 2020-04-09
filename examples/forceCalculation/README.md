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
A similar integral on a surface outside the TFSF source layer
doesn't give the force but does given an estimate for how much
of the field is scattered or delayed by the object.

The accuracy of the results depends on various simulation parameters
including resolution and number of integration points.
For a sphere, we can fairly accurately calculate the scattering and
corresponding force using Mie theory.
This is a good comparison for the accuracy of the FDTD implementation.

This directory contains the simulation source file `sim.cpp`
as well as a Matlab script `mie.m` to calculate the force using the
[optical tweezers toolbox](https://github.com/ilent2/ott).

