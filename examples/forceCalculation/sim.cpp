/** Example force calculation script
 *
 * Part of ot-cpp-fdtd, Copyright 2020 Isaac Lenton
 * See LICENSE for details about using/distributing this file.
 */

#include <string>
#include "FdtdIncludeAll.hpp"

class ExampleForceCalculation
{

  // Materials and wavelengths used in the simulation
  static constexpr double wavelength = 1064e-9;
  using Vacuum = Common::SiVacuum<wavelength>;
  using Sphere = Common::Dielectric<Common::Indices::polystyrene, Vacuum>;

  //
  // Indices and coordinate layer
  //

  // Sphere radius
  static constexpr double sphere_radius = 0.3*Vacuum::wavelength;

  // Calculate grid spacing
  static constexpr double resolution = 20.0;   // pts per wavelength
  static constexpr double wavelength_min
      = Common::wavelength_min_max<Vacuum, Sphere>::min;
  static constexpr double spacing = wavelength_min/resolution;

  // Calculate the size of the simulation space
  //    Particle size + output layer + TFSF layer + CPML layer
  // This probably doesn't need as larger padding, try experimenting with it
  static constexpr unsigned PAD = 3;
  static constexpr unsigned cpml_depth = 10;
  static constexpr unsigned output1_depth = cpml_depth + PAD;
  static constexpr unsigned tfsf_depth = output1_depth + PAD;
  static constexpr unsigned output2_depth = tfsf_depth + PAD;
  static constexpr unsigned grid_size = 2*(output2_depth + PAD)
      + make_even(ceil(2*sphere_radius/spacing));

  // Using a 3-D coordinate system, so need 3-D indices
  using Indices = Fdtd::Indices::Simple<grid_size, grid_size, grid_size>;

  // Cartesian coordinate system
  using Coordinates = Fdtd::Coordinate::Cartesian<spacing, spacing, spacing>;

  //
  // Timing layer
  //

  // Simulation length (number of optical cycles)
  static constexpr unsigned num_cycles = 20;

  // Use linear timing with a multiple of 3 steps per optical cycle
  // A multiple of three is often more accurate for calculating
  // the force integrals
  static constexpr unsigned num_dimensions = 3;
  static constexpr unsigned time_multiple = 3;
  static constexpr double max_time_step_size = Common::MaxTimeStepSize<
      spacing, num_dimensions, Vacuum, Sphere>::value();
  using Timing = Common::Timing::Linear<num_cycles, max_time_step_size,
        Vacuum, time_multiple>;

  //
  // Specify material layers and sphere shape
  //

  // Homogeneous permeability, inhomogeneous permittivity
  using Permittivity = Fdtd::Materials::SimplePermittivity;
  using Permeability = Fdtd::Materials::HomogeneousPermeability<
      Vacuum::permeability>;

  // Specify the class for initializing the permittivity
  struct Init : public Fdtd::EmptyModule
  {

    template <typename Base> class EvaluateMethods : public Base
    {
    public:
      void initialize(void)
      {
        Base::initialize();

        // Initialize E-field, H-field and Permittivity
        for (auto it = Base::dataBegin(); it != Base::dataEnd(); ++it)
        {
          it->efield() = Utilities::Vec3d();
          it->hfield() = Utilities::Vec3d();

          Utilities::Vec3d loc = Base::location(it);
          Utilities::Vec3d centre = Base::centre();
          if ((loc - centre).length() <= sphere_radius)
          {
            it->permittivity() = Sphere::permittivity;
          }
          else
          {
            it->permittivity() = Vacuum::permittivity;
          }
        }
      }
    };

  };

  //
  // TFSF Layer (beam and source)
  //

  static constexpr Utilities::Vec3d direction
      = Utilities::Vec3d(0, 0, 1);
  static constexpr Utilities::Vec3d polarisation
      = Utilities::Vec3d(0, 1, 0);
  static constexpr double intensity = 1.0;
  using Beam = Fdtd::Sources::PlaneWaveBeam<direction, polarisation,
        intensity, Vacuum::speed>;
  using Source = Fdtd::Sources::ContinuousWave<Beam, Vacuum::angular_frequency>;
  using Tfsf = Fdtd::Tfsf::Simple<Source, tfsf_depth>;

  //
  // Boundary layer (CPML)
  //

  // Using CPML surrounding the simulation space
  using Cpml = Fdtd::Cpml::Simple<cpml_depth, Timing::step_size,
      Vacuum::permittivity, Vacuum::permeability>;

  //
  // Setup force calculations and output files
  //

  // This section is a bit messy and probably shouldn't be implemented
  // using so many template classes.

  // Exterior and interior output volume/surfaces
  using force_offset1 = Utilities::OffsetBox<output1_depth>;
  using force_offset2 = Utilities::OffsetBox<output2_depth>;

  // Iterators for volume and surface
  using SurfaceIter1 = Fdtd::Output::Integrate::
      Iterators::Surface<force_offset1>;
  using SurfaceIter2 = Fdtd::Output::Integrate::
      Iterators::Surface<force_offset2>;
  using VolumeIter = Fdtd::Output::Integrate::
      Iterators::Volume<force_offset2>;

  // Output time window
  using StrideWindow = Fdtd::Output::StrideWindow<
      Timing::output_stride, Timing::multiple, -1, -1>;

public:

  // File names for output
  // Hmm, this feels like KLUDGE, they shouldn't be static!!!!
  // We had a way with extern strings that worked, but was bad too!
  static std::string output_volume;
  static std::string output_surface1;
  static std::string output_surface2;

private:

  // Only output the last values calculated to the output files
  static constexpr double offset = 0;  // <-- ignore this
  using SurfaceOutput1 = Fdtd::Output::Integrate::Output::FileLast<
    output_surface1, double, offset, false, false, true>;
  using SurfaceOutput2 = Fdtd::Output::Integrate::Output::FileLast<
    output_surface2, double, offset, false, false, true>;
  using VolumeOutput = Fdtd::Output::Integrate::Output::FileLast<
    output_volume, double, offset, false, false, true>;

  // Force calculation method (Maxwell stress tensor)
  using ForceMethod = Fdtd::Output::Integrate::Parameters::StressForce;

  // Assemble integrators
  using VolumeModule = Fdtd::Output::Integrate::Integrate<
      VolumeOutput, ForceMethod, VolumeIter, StrideWindow>;
  using SurfaceModule1 = Fdtd::Output::Integrate::Integrate<
      SurfaceOutput1, ForceMethod, SurfaceIter1, StrideWindow>;
  using SurfaceModule2 = Fdtd::Output::Integrate::Integrate<
      SurfaceOutput2, ForceMethod, SurfaceIter2, StrideWindow>;

public:

  //
  // Assemble simulation class
  //

  // Combine modules into a simulation object
  using Simulation = Fdtd::Simulation<Timing, Indices, Coordinates,
      Permittivity, Permeability, Tfsf, Cpml, Init,
      Fdtd::Materials::TimeAverageE<0>, Fdtd::Materials::TimeAverageH<0>,
      Fdtd::Output::Integrate::MaxwellStressTensor<true>,
      VolumeModule, SurfaceModule1, SurfaceModule2, Fdtd::Output::Report,
      Fdtd::Output::Progress<1*Timing::cycle_length>>;
};

std::string ExampleForceCalculation::output_volume;
std::string ExampleForceCalculation::output_surface1;
std::string ExampleForceCalculation::output_surface2;


int main(int argc, char** argv)
{
  // Configure simulation outputs
  ExampleForceCalculation::output_volume = "volume.dat";
  ExampleForceCalculation::output_surface1 = "surface_scat.dat";
  ExampleForceCalculation::output_surface2 = "surface_total.dat";

  // Setup and run the simulation
  ExampleForceCalculation::Simulation sim;
  sim.evaluate();
}

