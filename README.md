ot-cpp-fdtd
===========

This package is a prototype for a C++ template library for
doing optical tweezers simulations using the finite difference
time domain (FDTD) method.
The code was initially developed as part of an honours year project
and parts of the code were used to generate the figures in

>  Isaac C. D. Lenton, Alexander B. Stilgoe, Halina Rubinsztein-Dunlop
>  and Timo A. Nieminen., "Visual Guide to Optical Tweezers",
>  *European Journal of Physics* **38**(3), 034009 (2017)
>  https://doi.org/10.1088/1361-6404/aa6271

The code is release in the hope that it will be useful, however it
is far from complete and has very little documentation.
The code is part an experiment with template meta programming and
part a framework for testing different features and implementations
of FDTD.
That said, parts of the code could be a useful starting point for
writing memory efficient FDTD or adding optical force/torque
calculations to an existing FDTD package.
Verification is still needed for the different force and torque
calculation methods, and there are a few other features
that could make the package more widely usable.

This repository only includes the more completed parts of the project.
Feel free to make suggestions for additional features or submit
pull requests.

Usage
-----

To get started, it is probably best to take a look at one of
the example codes and modify it to achieve the desired result.
The package is a header-only C++ template library.
To describe a particular scattering problem you need to
construct the layers of a `Fdtd::Simulation` object to specify how
the simulation should be run, the problem dimensionality and coordinate
system, and any outputs such as force and torque.
This is done with a command similar to

  using Simulation = Fdtd::Simulation<Timing, Indices, Coordinates,
      Permittivity, Permeability, SourceLayer, BoundaryLayers,
      Initialization, Output>;

The only required parts are the `Timing` and `Indices`, but the other
parts are typically required in order to do useful work.
The other layers can be provided in almost any order, as long as
all the required features for a layer have been provided by
previous layers.
In the above example, the parameters are

  * `Timing` a class specifying how time should evolve through
    the simulation.  For example, for linear time steps you could
    use the `Fdtd::Timing::Linear` class

      const unsigned step_count = 50;
      const double step_size = 1.0e-4;
      using Timing = Fdtd::Timing::Linear<num_steps, step_size>;

  * `Indices` a class specifying the dimensionality of the problem.
    For most problems this will be a 3-D grid, for example

      const unsigned grid_size = 10;
      using Indices = Fdtd::Indices::Simple<grid_size, grid_size, grid_size>;

  * `Coordinates` specifies the coordinate system to use.  For
    Cartesian coordinate this could be

      const double spacing = 1.0e-6;
      using Coordinates = Fdtd::Coordinate::Cartesian<
        spacing, spacing, spacing>;

  * `Permittivity` and `Permeability` specify the material properties.
    These can be homogeneous or inhomogeneous, for example

      using Permittivity = Fdtd::Materials::SimplePermittivity;
      using Permeability = Fdtd::Materials::HomogeneousPermeability<
          Vacuum::permeability>;

  * `SourceLayer` can be a source such as a total field scattered field (TFSF)
    source, for example we could use a vector spherical wave function
    beam generated using the optical tweezers toolbox

      using Beam = Fdtd::Sources::ToolboxBeam<nm_file, ab_file, Water::speed,
        Fdtd::Offset::Centre, Fdtd::Offset::Location<beam_offset>>;
      using Source = Fdtd::Sources::ContinuousWave<Beam,
          Vacuum::angular_frequency>;
      const unsigned tfsf_depth = 5;
      using SourceLayer = Fdtd::Tfsf::Simple<Source, tfsf_depth>;

  * `BoundaryLayers` can be any absorbing or perfectly matched
    boundary layer, for example

      const unsigned cpml_depth = 5;
      using Cpml = Fdtd::Cpml::Simple<cpml_depth, Timing::step_size,
        Vacuum::permittivity, Vacuum::permeability>;

  * `Initialisation` this class is used to initialise any properties
    such as the permittivity for inhomogenous materials.
    For an example initialisation layer, see the example directory.

  * `Output` output layers can be specified to write field values,
    calculate properties such as forces or torques or simply output
    the progress of the simulation.  For example, to write out a
    report of the FDTD setup/progress

      using Output = Fdtd::Output::Report;

Then it is simply a matter of constructing the `Simulation` class
and running the simulation.
For example, a program's main function might look like

```c++
    int main(int argc, char** argv) {
      Simulation sim;
      sim.evaluate();
      return 0;
    }
```

Requirements
------------
* C++14 compatible compiler (tested with GNU GCC 7.3.1)

For the example Makefiles you will also need GNU Make.

How to cite
-----------

There is not official publication for this software yet,
however if you would like to cite something you can either
cite this repository

> Isaac C. D. Lenton, "ot-cpp-fdtd"
> https://github.com/ilent2/ot-cpp-fdtd

or cite the paper that first used this code

>  Isaac C. D. Lenton, Alexander B. Stilgoe, Halina Rubinsztein-Dunlop
>  and Timo A. Nieminen., "Visual Guide to Optical Tweezers",
>  *European Journal of Physics* **38**(3), 034009 (2017)
>  https://doi.org/10.1088/1361-6404/aa6271

If you use any or all of this software in a publication it would be
nice to hear from you, please get in contact.
If the software is useful to enough people it might be eligible for
submission to
[Software-Impacts](https://www.journals.elsevier.com/software-impacts/)
or a similar journal.

License
-------

This version of the code is licensed under the GNU GPLv3.
The package includes an alpha-release of the Optical Tweezers
Toolbox C++ implementation (from 2012), this is not part of
ot-cpp-fdtd, refer to Alexander Stilgoe or the [Optical Tweezers
Toolbox](https://github.com/ilent2/ott) for more information about
licensing these parts.

> Copyright (C) 2016 Isaac Lenton
>
> This program is free software: you can redistribute it and/or modify
> it under the terms of the GNU General Public License as published by
> the Free Software Foundation, either version 3 of the License, or
> (at your option) any later version.
>
> This program is distributed in the hope that it will be useful,
> but WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
> GNU General Public License for more details.

Further details can be found in the ``LICENSE`` file.
If you would like to use the toolbox for something not covered by
the license, please get in contact.

Suggestions/Improvements/Bug-fixes
----------------------------------

Feel free to send pull-requests, issue/feature requests or
send emails directly to me ([Isaac Lenton](mailto:uqilento@uq.edu.au)).
At the moment my goal is completing my PhD, so this project is
a low priority.  Feel free to get in contact if you would like to
work on the project or if you have any questions/feedback.

Acknowledgement
---------------

The original code was completed by Isaac Lenton as part of
his honours year project (supervised by Alexander Stilgoe and
Timo Nieminen).
This version was released by Isaac Lenton as part of his PhD
(suppervised by Timo Nieminen, Alexander Stilgoe and
Halina Rubinsztein-Dunlop) with the support of the Australian
Government Research Training Program (RTP) scholarship.

