/* Indices/Linear.hpp - Recursive linear indices definition.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_INDICES_LINEAR_HPP
#define FDTD_INDICES_LINEAR_HPP

#include "../EmptyModule.hpp"
#include <tuple>

namespace Fdtd {
namespace Indices {

/** A linear indices axis.
 *
 * @tparam Index        : The axis we are describing.
 * @tparam Length       : The length of the axis.
 */
template <int Index, unsigned Length>
struct Linear : public EmptyModule
{
  template <typename Base> class Fallback : public Base
  {
  public:

    /** Declare the starting index for resolving AxisIndex. */
    constexpr static int AxisIndex = 0;

    /** Return the length of an axis (0 indexed). */
    constexpr static unsigned length(int Axis)
    {
      return 0;
    }

    /** Returns the index volume of this and lower axes. */
    constexpr static unsigned rlength(int Axis=AxisIndex-1)
    {
      return 1;
    }

    /** Convert from a (i, j, k) index to a memory index. */
    template <typename... Args>
    constexpr static int toIndex(Args... args)
    {
      return 0;
    }

    /** Convert from a memory index to a (i, j, k) index. */
    template <typename... T>
    constexpr static std::tuple<T...> fromIndex(unsigned i)
    {
      return std::make_tuple();
    }

    /** Advance to the next element. */
    int nextIndex(int loc, int Axis, unsigned count=1) const
    {
      return -1;
    }

    /** Advance to the previous element. */
    int prevIndex(int loc, int Axis, unsigned count=1) const
    {
      return -1;
    }

  };

  template <typename Base> class Indices : public Base
  {
  public:

    /** Calculate our number of axes. */
    constexpr static int AxisIndex =
        std::max(Base::AxisIndex, Index + 1);

    /** Report we are using linear indices on this axis. */
    template <typename T> void report_construct(T& fp) const
    {
      fp << "Using Indices::Linear on axis " << Index
          << " with " << Length << " elements" << std::endl;
      Base::report_construct(fp);
    }

    /** Return the length of an axis (0 indexed). */
    constexpr static unsigned length(int Axis)
    {
      return Axis == Index ? Length : Base::length(Axis);
    }

    /** Returns the index volume of this and lower axes. */
    constexpr static unsigned rlength(int Axis=AxisIndex-1)
    {
      return Axis == Index ?
          Length * Base::rlength(Index-1) : Base::rlength(Axis);
    }

		/** Convert from a (i, j, k) index to a memory index. */
    template <typename... Args>
    constexpr static int toIndex(unsigned N, Args... args)
    {
      static_assert(sizeof...(Args) == Index,
          "Too many indices passed to toIndex");

      if (N >= length(Index)) return -1;
      int base = Base::toIndex(args...);
      if (base < 0) return -1;
      return base + N*Base::rlength(Index-1);
    }

    // TODO: Do we need to reintroduce an AxisCap (for fromIndex, etc.)

    /** Convert from a memory index to a (i, j, k) index. */
    constexpr static auto fromIndex(unsigned i)
        -> decltype(std::tuple_cat(std::make_tuple(1u), Base::fromIndex(0)))
    {
      unsigned n = i / Base::rlength(Index-1);
      unsigned m = i % Base::rlength(Index-1);
      return std::tuple_cat(std::make_tuple(n), Base::fromIndex(m));
    }

    /** Advance to the next element. */
    int nextIndex(unsigned loc, int Axis, unsigned count=1) const
    {
      if (Axis != Index) return Base::nextIndex(loc, Axis, count);
      return loc + count < length(Axis) ? loc + count : -1;
    }

    /** Advance to the previous element. */
    int prevIndex(unsigned loc, int Axis, unsigned count=1) const
    {
      if (Axis != Index) return Base::prevIndex(loc, Axis, count);
      return loc >= count ? loc - count : -1;
    }

  };
};

} // namespace Indices
} // namespace Fdtd

#endif // #ifndef FDTD_INDICES_LINEAR_HPP

