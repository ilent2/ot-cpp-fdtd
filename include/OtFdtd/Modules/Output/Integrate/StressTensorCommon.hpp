/* Output/Integrate/StressTensorCommon.hpp - Common parts for stress.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_INTEGRATE_STRESS_TENSOR_COMMON_HPP
#define FDTD_OUTPUT_INTEGRATE_STRESS_TENSOR_COMMON_HPP

#include "Modules/EmptyModule.hpp"
#include "Utilities/Vec.hpp"
#include "Utilities/Mat.hpp"

namespace Fdtd {
namespace Output {
namespace Integrate {

/** Common parts for stress tensor implementations.
 *
 * @tparam SubtractDc     : Subtract the DC field before calculation.
 */
template <bool SubtractDc=true>
struct StressTensorCommon : public EmptyModule
{
  /** Methods for calculating the stress tensor.
   *
   * TODO: Some of these methods could perhaps be merged elsewhere?
   */
  template <typename StressTensor, typename Base>
  class AdvanceMethods : public Base
  {
    static Utilities::Vec3d thisEfield(typename Base::DataIterator loc)
    {
      if (SubtractDc) return Base::acEfield(loc);
      else return loc->efield();
    }

    static Utilities::Vec3d thisHfield(typename Base::DataIterator loc)
    {
      if (SubtractDc) return Base::acHfield(loc);
      else return loc->hfield();
    }

  public:

    /** Calculate the effective E field at the E-origin.
     * Averages the surrounding E field values. */
    Utilities::Vec3d avgEatE(typename Base::DataIterator loc)
    {
      return Utilities::Vec3d(
          (thisEfield(Base::next(loc, 2))
              + thisEfield(loc)).axis(2) / 2.0,
          (thisEfield(Base::next(loc, 1))
              + thisEfield(loc)).axis(1) / 2.0,
          (thisEfield(Base::next(loc, 0))
              + thisEfield(loc)).axis(0) / 2.0);
    }

    /** Calculate the effective E field at the H-origin.
     * Averages the surrounding E field values. */
    Utilities::Vec3d avgEatH(typename Base::DataIterator loc)
    {
      return Utilities::Vec3d(
        0.25*(thisEfield(loc)
            + thisEfield(Base::prev(loc, 0))
            + thisEfield(Base::prev(loc, 1))
            + thisEfield(Base::prev(Base::prev(loc, 0), 1))).axis(2),
        0.25*(thisEfield(loc)
            + thisEfield(Base::prev(loc, 0))
            + thisEfield(Base::prev(loc, 2))
            + thisEfield(Base::prev(Base::prev(loc, 0), 2))).axis(1),
        0.25*(thisEfield(loc)
            + thisEfield(Base::prev(loc, 1))
            + thisEfield(Base::prev(loc, 2))
            + thisEfield(Base::prev(Base::prev(loc, 1), 2))).axis(0));
    }

    /** Calculate the effective H field at the E-origin.
     * Averages the surrounding H field values. */
    Utilities::Vec3d avgHatE(typename Base::DataIterator loc)
    {
      return Utilities::Vec3d(
        0.25*(thisHfield(loc)
            + thisHfield(Base::next(loc, 0))
            + thisHfield(Base::next(loc, 1))
            + thisHfield(Base::next(Base::next(loc, 0), 1))).axis(2),
        0.25*(thisHfield(loc)
            + thisHfield(Base::next(loc, 0))
            + thisHfield(Base::next(loc, 2))
            + thisHfield(Base::next(Base::next(loc, 0), 2))).axis(1),
        0.25*(thisHfield(loc)
            + thisHfield(Base::next(loc, 1))
            + thisHfield(Base::next(loc, 2))
            + thisHfield(Base::next(Base::next(loc, 1), 2))).axis(0));
    }

    /** Calculate the effective H field at the H-origin.
     * Averages the surrounding H field values. */
    Utilities::Vec3d avgHatH(typename Base::DataIterator loc)
    {
      return Utilities::Vec3d(
          (thisHfield(Base::prev(loc, 2))
              + thisHfield(loc)).axis(2) / 2.0,
          (thisHfield(Base::prev(loc, 1))
              + thisHfield(loc)).axis(1) / 2.0,
          (thisHfield(Base::prev(loc, 0))
              + thisHfield(loc)).axis(0) / 2.0);
    }

    Utilities::Vec3d surHatH(typename Base::DataIterator loc, int Axis)
    {
      // TODO: This is bad, but thesis is due
      return 0.5 * (avgHatH(loc) + avgHatH(Base::prev(loc, Axis)));
    }

    Utilities::Vec3d surHatE(typename Base::DataIterator loc, int Axis)
    {
      // TODO: This is bad, but thesis is due
      return 0.5 * (avgHatE(loc) + avgHatE(Base::prev(loc, Axis)));
    }

    Utilities::Vec3d surEatH(typename Base::DataIterator loc, int Axis)
    {
      // TODO: This is bad, but thesis is due
      return 0.5 * (avgEatH(loc) + avgEatH(Base::prev(loc, Axis)));
    }

    Utilities::Vec3d surEatE(typename Base::DataIterator loc, int Axis)
    {
      // TODO: This is bad, but thesis is due
      return 0.5 * (avgEatE(loc) + avgEatE(Base::prev(loc, Axis)));
    }

    /** Calculate the stress tensor centred at the H-origin. */
    Utilities::Mat3d stressTensorH(typename Base::DataIterator loc, int Axis)
    {
      auto avgE = surEatH(loc, Axis);
      auto avgH = surHatH(loc, Axis);
      // TODO: Should this use averaged permittivity?
      return StressTensor::stressTensor(avgE, avgH,
          loc->permittivity(), loc->permeability());
    }

    /** Calculate the stress tensor centred at the E-origin. */
    Utilities::Mat3d stressTensorE(typename Base::DataIterator loc, int Axis)
    {
      auto avgE = surEatE(loc, Axis);
      auto avgH = surHatE(loc, Axis);
      // TODO: Should this use averaged permeability?
      return StressTensor::stressTensor(avgE, avgH,
          loc->permittivity(), loc->permeability());
    }

    /** Calculate the stress tensor centred at the H-origin. */
    Utilities::Mat3d stressTensorH(typename Base::DataIterator loc)
    {
      auto avgE = avgEatH(loc);
      auto avgH = avgHatH(loc);
      // TODO: Should this use averaged permittivity?
      return StressTensor::stressTensor(avgE, avgH,
          loc->permittivity(), loc->permeability());
    }

    /** Calculate the stress tensor centred at the E-origin. */
    Utilities::Mat3d stressTensorE(typename Base::DataIterator loc)
    {
      auto avgE = avgEatE(loc);
      auto avgH = avgHatE(loc);
      // TODO: Should this use averaged permeability?
      return StressTensor::stressTensor(avgE, avgH,
          loc->permittivity(), loc->permeability());
    }
  };
};

} // namespace Integrate
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_INTEGRATE_STRESS_TENSOR_COMMON_HPP

