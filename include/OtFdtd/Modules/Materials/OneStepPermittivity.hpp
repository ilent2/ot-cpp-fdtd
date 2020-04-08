/* Fdtd/Materials/OneStepPermittivity.hpp - Base for simple permittivity.
 *
 * Purpose: Provides classes and methods common to homogeneous
 *    and inhomogeneous isotropic permittivity.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_ONE_STEP_PERMITTIVITY_HPP
#define FDTD_MATERIALS_ONE_STEP_PERMITTIVITY_HPP

#include "../EmptyModule.hpp"
#include "../../Utilities/Vec.hpp"

namespace Fdtd {
namespace Materials {
namespace _Permittivity {

struct OneStep : public EmptyModule
{

  template <typename Base> class ElementData : public Base
  {
    Utilities::Vec3d m_Efield;

  public:

    //@{
    /** Methods to get/set the E field. */
    const Utilities::Vec3d& efield(void) const { return m_Efield; }
    Utilities::Vec3d& efield(void) { return m_Efield; }
    //@}
  };

  template <typename Base> class AdvanceMethods : public Base
  {
  public:

    /** Calculate the curl of the E-field around each H-field axis. */
    Utilities::Vec3d efieldCurl(typename Base::DataIterator it) const
    {
      const auto& e1 = it->efield();
      const auto ex0 = Base::prevr(it, 2)->efield();
      const auto ey0 = Base::prevr(it, 1)->efield();
      const auto ez0 = Base::prevr(it, 0)->efield();

      // Approximate the derivative as the slope between the two
      // adjacent E-field locations, separated by the distance estep.
      return Utilities::Vec3d(
        (e1.z - ey0.z)/Base::estep(it, 1) - (e1.y - ez0.y)/Base::estep(it, 0),
        (e1.x - ez0.x)/Base::estep(it, 0) - (e1.z - ex0.z)/Base::estep(it, 2),
        (e1.y - ex0.y)/Base::estep(it, 2) - (e1.x - ey0.x)/Base::estep(it, 1));
    }

    /** Add a value to the E-field. */
    template <typename T>
    void dfieldAdd(typename Base::DataIterator loc, Utilities::Vec3t<T> val)
    {
      loc->efield() += val / loc->permittivity();
    }

    /** Add a value to the E-field. */
    void dfieldAdd(typename Base::DataIterator loc, double val, int axis)
    {
      loc->efield().axis(axis) += val / loc->permittivity();
    }
  };

  template <typename Base> class EvolveMethods : public Base
  {
  public:

    /** Clear the efield. */
    void initialize(void)
    {
      Base::initialize();

      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        loc->efield() = Utilities::Vec3d();
      }
    }

    /** The 1-step update equation for the E-field. */
    void advanceE(typename Base::TimingIterator it)
    {
      //#pragma omp parallel for
      //for (unsigned i = 0; i < Base::rlength(); ++i)
      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        //auto loc = Base::atIndex(i);
        Base::dfieldAdd(loc, Base::hfieldCurl(loc) * it->efieldDt());
      }

      Base::advanceE(it);
    }
  };

};

} // namespace _Permittivity
} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_ONE_STEP_PERMITTIVITY_HPP

