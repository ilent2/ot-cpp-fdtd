/** Demonstrate the deflection of an optical tweezers beam.
 *
 * Based on the simulation code for Figure 3 from
 *
 *   Isaac C. D. Lenton, Alexander B. Stilgoe, Halina Rubinsztein-Dunlop
 *   and Timo A. Nieminen., "Visual Guide to Optical Tweezers",
 *   *European Journal of Physics* **38**(3), 034009 (2017)
 *   https://doi.org/10.1088/1361-6404/aa6271
 *
 * Part of ot-cpp-fdtd, Copyright 2020 Isaac Lenton
 * See LICENSE for details about using/distributing this file.
 */

#include "FdtdIncludeAll.hpp"

struct DeflectionSim
{
  // Materials and wavelengths used in the simulation
  static constexpr double wavelength = 1064e-9;
  using Vacuum = Common::SiVacuum<wavelength>;
  using Sphere = Common::Dielectric<Common::Indices::water, Vacuum>;

  // Sphere properties
  static double sphere_offset;
  static constexpr double sphere_radius = 1.0*Vacuum::wavelength;

  //
  // Indices and coordinate layer
  //

  static constexpr double resolution = 20;
  static constexpr double spacing = Vacuum::wavelength/resolution;

  // Calculate size of simulation space
  static constexpr unsigned cpml_depth = 10;
  static constexpr unsigned PAD = 3;
  static constexpr unsigned tfsf_depth = cpml_depth + PAD;

  static constexpr unsigned sim_width = 2*tfsf_depth + 8*resolution;
  static constexpr unsigned sim_height = 2*tfsf_depth + 2*PAD
      + 2*sphere_radius/spacing;
  static constexpr unsigned sim_depth = 2*tfsf_depth + 16*resolution;

  using Indices = Fdtd::Indices::Simple<sim_height, sim_width, sim_depth>;
  using Coordinates = Fdtd::Coordinate::Cartesian<spacing, spacing, spacing>;

  //
  // Timing layer
  //

  // Simulation length (number of optical cycles)
  static constexpr unsigned num_cycles = 25;

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
  // Permittivity/Permeability layers
  //

  using Permeability = Fdtd::Materials::HomogeneousPermeability<
      Vacuum::permeability>;
  using Permittivity = Fdtd::Materials::SimplePermittivity;

  struct Init : public Fdtd::EmptyModule
  {
    template <typename Base> class EvaluateMethods : public Base
    {
      bool isSphere(double x, double y, double z)
      {
        z -= sim_depth/2.0;
        y -= sim_width/2.0;
        x -= sim_height/2.0;

        y += sphere_offset/spacing;

        if (sqrt(x*x + y*y + z*z) < sphere_radius/spacing) return true;
        return false;
      }

    public:
      void initialize(void)
      {
        Base::initialize();

        for (auto it = Base::dataBegin(); it != Base::dataEnd(); ++it)
        {
          unsigned x, y, z;
          std::tie(x, y, z) = Base::getOffset(it);

          if (isSphere(x, y, z))
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

  static const std::string nm_file;
  static const std::string ab_file;

  using Beam = Fdtd::Sources::ToolboxBeam<nm_file, ab_file, Vacuum::speed,
      Fdtd::Offset::Centre>;
  using Source = Fdtd::Sources::ContinuousWave<Beam, Vacuum::angular_frequency>;
  using Tfsf = Fdtd::Tfsf::Simple<Source, tfsf_depth>;

  //
  // Boundary layer (CPML)
  //

  // Using CPML surrounding the simulation space
  using Cpml = Fdtd::Cpml::Simple<cpml_depth, Timing::step_size,
      Vacuum::permittivity, Vacuum::permeability>;

  //
  // Output module
  //

  // Only output frames from the last optical cycle
  using StrideWindow = Fdtd::Output::StrideWindow<
      Timing::output_stride, Timing::multiple,
      Timing::total_steps-Timing::cycle_length, -1>;

  // Declare the module for writing the E-field (along beam axis)
  static std::string output_file;
  using Write = Fdtd::Output::WritePlane<output_file,
      Fdtd::Output::Parameters::EfieldMinusTimeAverage,
      2, sim_height/2, StrideWindow>;

  //
  // Assemble simulation object
  //

  using Simulation = Fdtd::Simulation<Timing, Indices, Coordinates,
    Permittivity, Fdtd::Materials::TimeAverageE<0>,
    Permeability, Fdtd::Materials::TimeAverageH<0>,
    Tfsf, Cpml, Init, Write,
    Fdtd::Output::Report, Fdtd::Output::Progress<1*Timing::cycle_length>>;
};

const std::string DeflectionSim::nm_file = "lgmn.dat";
const std::string DeflectionSim::ab_file = "lgab.dat";
std::string DeflectionSim::output_file;
double DeflectionSim::sphere_offset;

int main(int argc, char** argv)
{
  DeflectionSim::Simulation sim;

  // Simulate no deflection (just focussing)
  DeflectionSim::sphere_offset = 0.1*DeflectionSim::sphere_radius;
  DeflectionSim::output_file = "plane_0.dat";
  sim.evaluate();

  // Simulate deflection (almost near maximum force location)
  DeflectionSim::sphere_offset = DeflectionSim::sphere_radius;
  DeflectionSim::output_file = "plane_1.dat";
  sim.evaluate();

  return 0;
}

