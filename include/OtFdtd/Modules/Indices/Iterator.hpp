/* Indices/Iterator.hpp - Base class for indices iterators.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_INDICES_ITERATOR_HPP
#define FDTD_INDICES_ITERATOR_HPP

#include <iterator>

namespace Fdtd {
namespace Indices {

/** Iterator for element data type. */
template <typename ElementData>
class Iterator
    : public std::iterator<std::bidirectional_iterator_tag, ElementData>
{
public:

  Iterator(ElementData* start, unsigned length, unsigned offset)
      : m_pValue(start+offset), m_pStart(start), m_unLength(length),
        m_bValid(true) {}

  Iterator(void) : m_bValid(false) {}

  Iterator& operator++(void) { ++m_pValue; return *this; }
  Iterator& operator--(void) { --m_pValue; return *this; }

  Iterator operator+(int i) const
  {
    Iterator it(*this);
    it.m_pValue += i;
    return it;
  }
  Iterator operator-(int i) const
  {
    Iterator it(*this);
    it.m_pValue -= i;
    return it;
  }

  Iterator operator+(unsigned i) const
  {
    Iterator it(*this);
    it.m_pValue += i;
    return it;
  }
  Iterator operator-(unsigned i) const
  {
    Iterator it(*this);
    it.m_pValue -= i;
    return it;
  }

  unsigned getOffset(void) const
  {
    return m_pValue - m_pStart;
  }

  operator bool() const
  {
    // Note, this only checks for a valid memory location
    //    We also need to check if the axis is valid
    return m_bValid && m_pValue >= m_pStart && m_pValue < m_pStart+m_unLength;
  }

  bool operator! (void) const
  {
    return !((bool) *this);
  }

  bool operator==(const Iterator& o) const { return m_pValue == o.m_pValue; }
  bool operator!=(const Iterator& o) const { return m_pValue != o.m_pValue; }

  ElementData& operator*(void) { return *m_pValue; }
  const ElementData& operator*(void) const { return *m_pValue; }

  ElementData* operator->(void) { return m_pValue; }
  const ElementData* operator->(void) const { return m_pValue; }

protected:
  ElementData* m_pValue;
  ElementData* m_pStart;
  unsigned m_unLength;
  bool m_bValid;
};

template <typename ElementData>
class ConstIterator : public Iterator<ElementData>
{
  using Base = Iterator<ElementData>;

public:

  // TODO: This feels like kludge, using a const cast...
  ConstIterator(const ElementData* start, unsigned length, unsigned offset)
    : Base(const_cast<ElementData*>(start), length, offset) {}

  ConstIterator(const Base& base) : Base(base) {}

  ConstIterator(void) : Base() {}

  const ElementData& operator*(void) { return *Base::m_pValue; }
  const ElementData* operator->(void) { return Base::m_pValue; }

};

} // namespace Indices
} // namespace Fdtd

#endif // #ifndef FDTD_INDICES_ITERATOR_HPP

