/* Fdtd/Timing/Linear.hpp - Linear timing description.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_TIMING_LINEAR_HPP
#define FDTD_TIMING_LINEAR_HPP

#include "../EmptyModule.hpp"
#include "TimingIterator.hpp"

namespace Fdtd {
namespace Timing {

/** Default offset for Timing/Linear. */
class _Linear { public: static const double Offset; };
const double _Linear::Offset = 0.0;

template <unsigned StepCount, const double& StepSize,
    const double& Offset=_Linear::Offset>
struct Linear : public EmptyModule
{

  /** Data for a particular time step. */
  class _Data
  {
    unsigned m_unIndex;

  public:
    _Data(unsigned Index) : m_unIndex(Index) {}
    constexpr const double& efieldDt(void) const { return StepSize; }
    constexpr const double& hfieldDt(void) const { return StepSize; }

    /** Time of the E-field step currently being computed. */
    constexpr double efieldT(void) const
        { return StepSize*(m_unIndex+1.0) + Offset; }

    /** Time of the last E-field step computed. */
    constexpr double efieldTlast(void) const
        { return StepSize*(m_unIndex) + Offset; }

    /** Time of the H-field step currently being computed. */
    constexpr double hfieldT(void) const
        { return StepSize*(m_unIndex+0.5) + Offset; }

    /** Time of the last H-field step computed. */
    constexpr double hfieldTlast(void) const
        { return StepSize*(m_unIndex-0.5) + Offset; }

    constexpr unsigned index(void) const { return m_unIndex; }

    _Data& operator++(void)
    {
      m_unIndex = std::min(StepCount, m_unIndex+1);
      return *this;
    }
  };

  template <typename Base> class Timing : public Base
  {
  public:

    /** The timing type to use for the simulation. */
    using TimingIterator = Iterator<_Data>;

    template <typename T>
    void report_construct(T& fp) const
    {
      fp << "Timing:Linear with " << StepCount << " steps "
          << StepSize << " long." << std::endl;
      Base::report_construct(fp);
    }

    /** Return iterator for first index. */
    TimingIterator TimingBegin(void) const
    {
      return TimingIterator(0);
    }

    /** Return iterator for last index. */
    TimingIterator TimingEnd(void) const
    {
      return TimingIterator(StepCount);
    }

  };

};

} // namespace Timing
} // namespace Fdtd

#endif // #ifndef FDTD_TIMING_LINEAR_HPP

