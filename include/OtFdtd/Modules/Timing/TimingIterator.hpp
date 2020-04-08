/* Fdtd/Timing/TimingIterator.hpp - A common timing iterator type.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_TIMING_TIMING_ITERATOR_HPP
#define FDTD_TIMING_TIMING_ITERATOR_HPP

namespace Fdtd {
namespace Timing {

/** Time step iterator. */
template <typename TimingData>
class Iterator
    : public std::iterator<std::forward_iterator_tag, TimingData>
{
public:
  Iterator(unsigned index) : m_Value(index) {}
  Iterator& operator++(void) { ++m_Value; return *this; }
  bool operator==(const Iterator& o) const
      { return o.m_Value.index() == m_Value.index(); }
  bool operator!=(const Iterator& o) const
      { return o.m_Value.index() != m_Value.index(); }
  TimingData& operator*(void) { return m_Value; }
  const TimingData& operator*(void) const { return m_Value; }
  TimingData* operator->(void) { return &m_Value; } 
  const TimingData* operator->(void) const { return &m_Value; }

private:
  TimingData m_Value;
};

} // namespace Timing
} // namespace Fdtd

#endif // #ifndef FDTD_TIMING_TIMING_ITERATOR_HPP

