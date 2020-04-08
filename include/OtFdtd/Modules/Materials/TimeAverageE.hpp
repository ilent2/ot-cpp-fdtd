/* Materials/TimeAverageE.hpp - Calculate the time average of the E field.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_TIME_AVERAGE_E_HPP
#define FDTD_MATERIALS_TIME_AVERAGE_E_HPP

#include "../EmptyModule.hpp"

namespace Fdtd {
namespace Materials {

/** Calculate the time average since the simulation start.
 *
 * @tparam zeroBias = The initial weighting towards zero DC field.
 */
template <unsigned zeroBias=0, unsigned Start=0>
struct TimeAverageE : public EmptyModule
{

  template <typename Base> class ElementData : public Base
  {
    Utilities::Vec3d m_EfieldAverage;

  public:

    //@{
    /** Methods to get/set the E field. */
    const Utilities::Vec3d& efieldAverage(void) const
        { return m_EfieldAverage; }
    Utilities::Vec3d& efieldAverage(void)
        { return m_EfieldAverage; }
    //@}
  };

  template <typename Base> struct AdvanceMethods : public Base
  {
    static Utilities::Vec3d acEfield(typename Base::DataIterator loc)
    {
      return loc->efield() - loc->efieldAverage();
    }
  };

  template <typename Base> class EvolveMethods : public Base
  {
    double dT;

    void add_value(void)
    {
      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        loc->efieldAverage() =
            loc->efieldAverage()*dT/(dT + 1.0) + loc->efield()/(dT + 1.0);
      }
      dT += 1.0;
    }

  public:

    /** Clear the efieldAverage. */
    void initialize(void)
    {
      Base::initialize();

      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        loc->efieldAverage() = Utilities::Vec3d();
      }
      dT = double(zeroBias);
    }

    /** Accumulate the E-field average. */
    void inspectE(typename Base::TimingIterator it)
    {
      if (it->index() >= Start) add_value();
      Base::inspectE(it);
    }
  };

};

} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_MATERIALS_TIME_AVERAGE_E_HPP

