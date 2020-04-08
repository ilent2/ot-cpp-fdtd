/* Materials/TimeAverageH.hpp - Calculate the time average of the H-field.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_MATERIALS_TIME_AVERAGE_H_HPP
#define FDTD_MATERIALS_TIME_AVERAGE_H_HPP

#include "../EmptyModule.hpp"

namespace Fdtd {
namespace Materials {

/** Calculate the time average since the simulation start.
 *
 * @tparam zeroBias = The initial weighting towards zero DC field.
 */
template <unsigned zeroBias=0, unsigned Start=0>
struct TimeAverageH : public EmptyModule
{

  template <typename Base> class ElementData : public Base
  {
    Utilities::Vec3d m_HfieldAverage;

  public:

    //@{
    /** Methods to get/set the H field. */
    const Utilities::Vec3d& hfieldAverage(void) const
        { return m_HfieldAverage; }
    Utilities::Vec3d& hfieldAverage(void)
        { return m_HfieldAverage; }
    //@}
  };

  template <typename Base> struct AdvanceMethods : public Base
  {
    static Utilities::Vec3d acHfield(typename Base::DataIterator loc)
    {
      return loc->hfield() - loc->hfieldAverage();
    }
  };

  template <typename Base> class EvolveMethods : public Base
  {
    double dT;

    void add_value(void)
    {
      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        loc->hfieldAverage() =
            loc->hfieldAverage()*dT/(dT + 1.0) + loc->hfield()/(dT + 1.0);
      }
      dT += 1.0;
    }

  public:

    /** Clear the hfieldAverage. */
    void initialize(void)
    {
      Base::initialize();

      for (auto loc = Base::dataBegin(); loc != Base::dataEnd(); ++loc)
      {
        loc->hfieldAverage() = Utilities::Vec3d();
      }
      dT = double(zeroBias);
    }

    /** Accumulate the H-field average. */
    void inspectH(typename Base::TimingIterator it)
    {
      if (it->index() >= Start) add_value();
      Base::inspectH(it);
    }
  };

};

} // namespace Materials
} // namespace Fdtd

#endif // #ifndef FDTD_MATERIALS_TIME_AVERAGE_H_HPP

