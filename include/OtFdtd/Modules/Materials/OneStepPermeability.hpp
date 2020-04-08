/* Fdtd/Materials/OneStepPermeability.hpp - Base for simple permeability.
 *
 * Purpose: Provides classes and methods common to homogeneous
 *    and inhomogeneous isotropic permeability.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_ONE_STEP_PERMEABILITY_HPP
#define FDTD_MATERIALS_ONE_STEP_PERMEABILITY_HPP

#include "../EmptyModule.hpp"
#include "../../Utilities/Vec.hpp"

namespace Fdtd {
namespace Materials {
namespace _Permeability {

struct OneStep : public EmptyModule
{

  template <typename Base> class ElementData : public Base
  {
    Utilities::Vec3d m_Hfield;

  public:

    //@{
    /** Methods to get/set the H field. */
    const Utilities::Vec3d& hfield(void) const { return m_Hfield; }
    Utilities::Vec3d& hfield(void) { return m_Hfield; }
    //@}
  };

  template <typename Base> class AdvanceMethods : public Base
  {
  public:

    /** Calculate the curl of the H-field around each E-field axis. */
    Utilities::Vec3d hfieldCurl(typename Base::DataIterator it) const
    {
      const auto& h0 = it->hfield();
      const auto hx1 = Base::nextr(it, 2)->hfield();
      const auto hy1 = Base::nextr(it, 1)->hfield();
      const auto hz1 = Base::nextr(it, 0)->hfield();

      // Approximate the derivative as the slope between the two
      // adjacent H-field locations, separated by the distance hstep.
      return Utilities::Vec3d(
        (hy1.z - h0.z)/Base::hstep(it, 1) - (hz1.y - h0.y)/Base::hstep(it, 0),
        (hz1.x - h0.x)/Base::hstep(it, 0) - (hx1.z - h0.z)/Base::hstep(it, 2),
        (hx1.y - h0.y)/Base::hstep(it, 2) - (hy1.x - h0.x)/Base::hstep(it, 1));
    }

    /** Add a value to the H-field. */
    template <typename T>
    void bfieldAdd(typename Base::DataIterator loc, Utilities::Vec3t<T> val)
    {
      loc->hfield() += val / loc->permeability();
    }

    /** Add a value to the H-field. */
    void bfieldAdd(typename Base::DataIterator loc, double val, int axis)
    {
      loc->hfield().axis(axis) += val / loc->permeability();
    }
  };

  template <typename Base> class EvolveMethods : public Base
  {
  public:

    /** Clear the hfield. */
    void initialize(void)
    {
      Base::initialize();

      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        loc->hfield() = Utilities::Vec3d();
      }
    }

    /** The 1-step update equation for the H-field. */
    void advanceH(typename Base::TimingIterator it)
    {
      //#pragma omp parallel for
      //for (unsigned i = 0; i < Base::rlength(); ++i)
      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        //auto loc = Base::atIndex(i);
        Base::bfieldAdd(loc, -Base::efieldCurl(loc) * it->hfieldDt());
      }

      Base::advanceH(it);
    }
  };

};

} // namespace _Permeability
} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_ONE_STEP_PERMEABILITY_HPP

