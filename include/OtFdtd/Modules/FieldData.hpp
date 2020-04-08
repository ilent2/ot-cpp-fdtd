/* Fdtd/FieldData.hpp - Data provider for field objects.
 *
 * This class allocates the memory at each location of the field.
 * The type (and amount of memory) allocated depends the template type.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_FIELD_DATA_HPP
#define FDTD_FIELD_DATA_HPP

#include "EmptyModule.hpp"
#include "Indices/Iterator.hpp"
#include <tuple>

namespace Fdtd {

/** Allocates the field data for field objects.
 *
 * @tparam ElementData : The resolved field element to use for the field.
 */
template <typename ElmData>
struct FieldDataModule : public EmptyModule
{
  template <typename Base> class Fallback : public Base
  {
  public:
    /** Base class for data iterator type. */
    template <typename B>
    using DataIterator = Indices::Iterator<B>;

    /** Base class for const data iterator type. */
    template <typename B>
    using ConstDataIterator = Indices::ConstIterator<B>;
  };

  template <typename Base> class FieldData : public Base
  {
  public:

    /** Allocate memory. */
    FieldData(void)
    {
      m_Data = new ElmData[Base::rlength()];
    }

    /** Free memory. */
    ~FieldData(void)
    {
      delete[] m_Data;
    }

    /** Finalize the DataIterator type. */
    using DataIterator = typename Base::template DataIterator<ElmData>;

    /** Finalize the ConstDataIterator type. */
    using ConstDataIterator =
        typename Base::template ConstDataIterator<ElmData>;

    //@{
    /** Return an iterator to the start of the data structure. */
    DataIterator dataBegin(void)
    {
      return DataIterator(m_Data, Base::rlength(), 0);
    }
    ConstDataIterator cdataBegin(void) const
    {
      return ConstDataIterator(m_Data, Base::rlength(), 0);
    }
    //@}

    //@{
    /** Return iterator to the element past the end of the data structure. */
    DataIterator dataEnd(void)
    {
      return DataIterator(m_Data, Base::rlength(), Base::rlength());
    }
    ConstDataIterator cdataEnd(void) const
    {
      return ConstDataIterator(m_Data, Base::rlength(), Base::rlength());
    }
    //@}

    //@{
    /** Get an iterator for a given coordinate. */
    template <typename... Indices>
    DataIterator at(Indices... args)
    {
      return DataIterator(m_Data, Base::rlength(), Base::toIndex(args...));
    }
    template <typename... Indices>
    ConstDataIterator at(Indices... args) const
    {
      return ConstDataIterator(m_Data,
          Base::rlength(), Base::toIndex(args...));
    }
    //}

    //@{
    /** Get an iterator for a given index. */
    DataIterator atIndex(unsigned index)
    {
      return DataIterator(m_Data, Base::rlength(), index);
    }
    ConstDataIterator at(unsigned index) const
    {
      return ConstDataIterator(m_Data, Base::rlength(), index);
    }
    //}

    /** Get the memory offset from an iterator. */
    unsigned getIndex(ConstDataIterator it) const
    {
      // TODO: We were using iterator difference, but this didn't
      //    seem to work correctly, don't know why.  Use caution!
      return it.getOffset();
    }

    /** Get the offset along an axis from an iterator. */
    unsigned getOffset(ConstDataIterator it, int Axis) const
    {
      auto lengths = Base::fromIndex(getIndex(it));
      if (Axis == 2) return std::get<0>(lengths);
      if (Axis == 1) return std::get<1>(lengths);
      if (Axis == 0) return std::get<2>(lengths);
      throw std::logic_error("Invalid axis in getOffset");
    }

    /** Get offset for each axis. */
    auto getOffset(ConstDataIterator it) const
        -> decltype(Base::fromIndex(0))
    {
      return Base::fromIndex(getIndex(it));
    }

    /** Advance to the next element location along axis. */
    template <typename ItType>
    ItType next(const ItType& loc, int Axis, unsigned count=1) const
    {
      if (!loc) return ItType();

      unsigned offset = getOffset(loc, Axis);
      int nextIndex = Base::nextIndex(offset, Axis, count);
      if (nextIndex < 0) return ItType();

      int newIndex = nextIndex - (int) offset;
      newIndex *= Base::rlength(Axis-1);

      return loc + newIndex;
    }

    /** Advance to the previous element location along axis. */
    template <typename ItType>
    ItType prev(const ItType& loc, int Axis, unsigned count=1) const
    {
      if (!loc) return ItType();

      unsigned offset = getOffset(loc, Axis);
      int prevIndex = Base::prevIndex(offset, Axis, count);
      if (prevIndex < 0) return ItType();

      int newIndex = prevIndex - (int) offset;
      newIndex *= Base::rlength(Axis-1);

      return loc + newIndex;
    }

    // TODO: There are multiple values we could return for nextr/prevr,
    //    at the moment we are returning the input so we
    //    calculate a zero derivative in h/efieldCurl.

    /** Return the (read-only) next element to use for calculations. */
    template <typename ItType>
    ConstDataIterator nextr(const ItType& loc,
        int Axis, unsigned count=1) const
    {
      auto newloc = next(loc, Axis, count);
      // TODO: return newloc ? newloc : Utilities::Vec3d();
      return newloc ? newloc : loc;
    }

    /** Return the (read-only) previous element to use for calculations. */
    template <typename ItType>
    ConstDataIterator prevr(const ItType& loc,
        int Axis, unsigned count=1) const
    {
      auto newloc = prev(loc, Axis, count);
      // TODO: return newloc ? newloc : Utilities::Vec3d();
      return newloc ? newloc : loc;
    }

  private:
    ElmData* m_Data;        ///< Data stored at field locations
  };
};

} // namespace Fdtd

#endif // #ifndef FDTD_FIELD_DATA_HPP

