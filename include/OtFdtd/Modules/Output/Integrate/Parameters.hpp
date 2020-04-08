/* Integrate/Parameters.hpp - Description of the integration kernels.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_INTEGRATE_PARAMETERS_HPP
#define FDTD_OUTPUT_INTEGRATE_PARAMETERS_HPP

#include "Utilities/Vec.hpp"

namespace Fdtd {
namespace Output {
namespace Integrate {
namespace Parameters {

/** A kernel for the Poynting vector flux. */
struct PoyntingFlux
{
  template <typename Base>
  struct Parameter : public Base
  {
    Utilities::Vec3d SurfaceKernel(typename Base::TimingIterator it,
        typename Base::DataIterator loc, Utilities::Vec3d norm, int Axis)
    {
      return norm * Base::surEatE(loc, Axis).cross(Base::surHatE(loc, Axis));
    }
  };

  static std::string repr(void)
  {
    return "PoyntingFlux";
  }
};

/** Calculate the spin angular momentum assuming a single frequency.
 *
 * This method calculates the degree of circular polarization of the
 * scattered light.  An extra negation can be included when the method
 * is used outside one TFSF plane.
 */
template <bool Negate=false>
struct SpinAngularMomentum
{
  template <typename Base>
  class Parameter : public Base
  {
    Utilities::Vec3d thisEfield(typename Base::DataIterator loc)
    {
      // TODO: Needs cleaning
      return Base::efieldLast(loc);
    }

    /** Calculate the effective E field at the E-origin.
     * Averages the surrounding E field values. */
    Utilities::Vec3d delayedEfieldAvg(typename Base::DataIterator loc)
    {
      // TODO: Needs cleaning
      return Utilities::Vec3d(
          (thisEfield(Base::next(loc, 2))
              + thisEfield(loc)).axis(2) / 2.0,
          (thisEfield(Base::next(loc, 1))
              + thisEfield(loc)).axis(1) / 2.0,
          (thisEfield(Base::next(loc, 0))
              + thisEfield(loc)).axis(0) / 2.0);
    }

  public:

    Utilities::Vec3d SurfaceKernel(typename Base::TimingIterator it,
        typename Base::DataIterator loc, Utilities::Vec3d norm, int Axis)
    {
      auto efield = Base::surEatE(loc, Axis);
      auto delayedEfield = 0.5*(delayedEfieldAvg(loc)
          + delayedEfieldAvg(Base::prev(loc, Axis)));
      auto eps = pow(loc->isoPermittivity()/loc->isoPermeability(), 0.5);
      auto D1 = eps * Utilities::Vec3d(efield.y*delayedEfield.z,
              efield.x*delayedEfield.z, efield.x*delayedEfield.y) * norm;

      if (Negate && norm.dot(Utilities::Vec3d(1.0, 1.0, 1.0)) < 0.0)
      {
        return -2.0 * D1;
      }

      return 2.0 * D1;
    }
  };

  static std::string repr(void)
  {
    return "SpinAngularMomentum";
  }
};

/** Calculate the spin angular momentum assuming a single frequency. */
template <typename Origin>
struct StressTorque
{
  template <typename Base>
  struct Parameter : public Base
  {
    /** Calculate the torque using the surface integral method.
     *
     * This might omit the spin angular momentum of the beam.
     */
    Utilities::Vec3d SurfaceKernel(typename Base::TimingIterator it,
        typename Base::DataIterator loc, Utilities::Vec3d norm, int Axis)
    {
      // Calculate the maxwell stress tensor for each point
      auto stress = Base::stressTensorE(loc, Axis);

      // Calculate the radial vectors
      auto co = Base::location(loc) - Origin::offset(*this);

      // Cross the stress tensor with the location (from Benito08)
      //    This doesn't appear to include the spin angular momentum
      return norm * stress.cross(co);
    }

    /** Calculate the torque using the volume method.
     *
     * This only works for explicit permittivity/permeability,
     * and will probably not work for ADE conductivity.
     *
     * This will be very expensive, so we might as well calculate
     * 8 stress tensors instead of working out a more optimal way.
     * I know there is at least a method with 6, and I suspect we
     * could achieve better efficiency by expanding it out further.
     */
    Utilities::Vec3d VolumeKernel(typename Base::TimingIterator it,
        typename Base::DataIterator loc, double volume)
    {
      // Calculate 8 adjacent stress tensors
      auto a1 = Base::stressTensorH(loc);
      auto a2 = Base::stressTensorH(Base::next(loc, 0));
      auto a3 = Base::stressTensorH(Base::next(Base::next(loc, 0), 1));
      auto a4 = Base::stressTensorH(Base::next(loc, 1));

      auto b1 = Base::stressTensorH(Base::next(loc, 2));
      auto b2 = Base::stressTensorH(Base::next(Base::next(loc, 2), 0));
      auto b3 = Base::stressTensorH(
          Base::next(Base::next(Base::next(loc, 2), 0), 1));
      auto b4 = Base::stressTensorH(Base::next(Base::next(loc, 2), 1));

      // Calculate 6 adjacent stress tensors (averaged)
      auto x0 = (a1 + a2 + a3 + a4) / 4.0;
      auto x1 = (b1 + b2 + b3 + b4) / 4.0;

      auto y0 = (a1 + a2 + b1 + b2) / 4.0;
      auto y1 = (a3 + a4 + b3 + b4) / 4.0;

      auto z0 = (a1 + a4 + b1 + b4) / 4.0;
      auto z1 = (a2 + a3 + b2 + b3) / 4.0;

      // Calculate derivatives
      auto dz = (z1 - z0) / Base::hstep(loc, 0);
      auto dy = (y1 - y0) / Base::hstep(loc, 1);
      auto dx = (x1 - x0) / Base::hstep(loc, 2);

      auto coord = Base::location(loc) - Origin::offset(*this);

      // Calculate force and hope the compiler optimises
      return coord.cross((dz.row(0) + dy.row(1) + dx.row(2))) * volume;
    }
  };

  static std::string repr(void)
  {
    return "StressTorque";
  }
};

/** Calculate the force using the stress tensor method. */
struct StressForce
{
  template <typename Base>
  struct Parameter : public Base
  {
    /** Calculate the force using the surface integral method. */
    Utilities::Vec3d SurfaceKernel(typename Base::TimingIterator it,
        typename Base::DataIterator loc, Utilities::Vec3d norm, int Axis)
    {
      // Calculate the maxwell stress tensor for each point
      auto stress = Base::stressTensorE(loc, Axis);

      // Calculate the momentum flux
      return norm * stress;
    }

    /** Calculate the torque using the volume method.
     *
     * This only works for explicit permittivity/permeability,
     * and will probably not work for ADE conductivity.
     *
     * This will be very expensive, so we might as well calculate
     * 8 stress tensors instead of working out a more optimal way.
     * I know there is at least a method with 6, and I suspect we
     * could achieve better efficiency by expanding it out further.
     */
    Utilities::Vec3d VolumeKernel(typename Base::TimingIterator it,
        typename Base::DataIterator loc, double volume)
    {
      // Calculate 8 adjacent stress tensors
      auto a1 = Base::stressTensorH(loc);
      auto a2 = Base::stressTensorH(Base::next(loc, 0));
      auto a3 = Base::stressTensorH(Base::next(Base::next(loc, 0), 1));
      auto a4 = Base::stressTensorH(Base::next(loc, 1));

      auto b1 = Base::stressTensorH(Base::next(loc, 2));
      auto b2 = Base::stressTensorH(Base::next(Base::next(loc, 2), 0));
      auto b3 = Base::stressTensorH(
          Base::next(Base::next(Base::next(loc, 2), 0), 1));
      auto b4 = Base::stressTensorH(Base::next(Base::next(loc, 2), 1));

      // Calculate 6 adjacent stress tensors (averaged)
      auto x0 = (a1 + a2 + a3 + a4) / 4.0;
      auto x1 = (b1 + b2 + b3 + b4) / 4.0;

      auto y0 = (a1 + a2 + b1 + b2) / 4.0;
      auto y1 = (a3 + a4 + b3 + b4) / 4.0;

      auto z0 = (a1 + a4 + b1 + b4) / 4.0;
      auto z1 = (a2 + a3 + b2 + b3) / 4.0;

      // Calculate derivatives
      auto dz = (z1 - z0) / Base::hstep(loc, 0);
      auto dy = (y1 - y0) / Base::hstep(loc, 1);
      auto dx = (x1 - x0) / Base::hstep(loc, 2);

      // Calculate force and hope the compiler optimises
      return (dz.row(0) + dy.row(1) + dx.row(2)) * volume;
    }
  };

  static std::string repr(void)
  {
    return "StressForce";
  }
};

/** Calculate the spin angular momentum assuming a single frequency. */
struct EdTorque
{
  template <typename Base>
  struct Parameter : public Base
  {
    /** Calculate the torque from the E and D fields. */
    Utilities::Vec3d VolumeKernel(typename Base::TimingIterator it,
        typename Base::DataIterator loc, double volume)
    {
      auto E = Base::avgEatE(loc);
      auto D = loc->permittivity() * E;
      return E.cross(D) * volume;
    }
  };

  static std::string repr(void)
  {
    return "EdTorque";
  }
};

} // namespace Parameters
} // namespace Integrate
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_INTEGRATE_PARAMETERS_HPP

