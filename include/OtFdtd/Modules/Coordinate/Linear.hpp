/* Coordinate/Linear.hpp - Linear coordinate axis.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_COORDINATE_LINEAR_HPP
#define FDTD_COORDINATE_LINEAR_HPP

#include "../EmptyModule.hpp"

namespace Fdtd {
namespace Coordinate {

/** A linear coordinate axis.
 *
 * @tparam Index   : The axis to use linear indices for.
 * @tparam Spacing : The linear step size.
 */
template <int Index, const double& Spacing>
struct Linear : public EmptyModule
{

  template <typename Base> class Fallback : public Base
  {
  public:

    /** Declare the starting index for resolving CoordinateIndex. */
    constexpr static int CoordinateIndex = 0;

    /** Return the distance between two adjacent E-field locations. */
    template <typename T>
    double hstep(T it, int Axis=Index) const
    {
      throw std::logic_error("Invalid Axis for hstep.");
    }

    /** Return the distance between two adjacent H-field locations. */
    template <typename T>
    double estep(T it, int Axis=Index) const
    {
      throw std::logic_error("Invalid Axis for estep.");
    }

    /** Compute the offset of the component vector from the origin. */
    template <typename T>
    double elength(T it, int Axis=Index, int Component=-1) const
    {
      return 0.0;
    }

    /** Compute the offset of the component vector from the origin. */
    template <typename T>
    double hlength(T it, int Axis=Index, int Component=-1) const
    {
      return 0.0;
    }
  };

  template <typename Base> class CoordinateMethods : public Base
  {
  public:

    /** Calculate the number of coordinate axes. */
    constexpr static int CoordinateIndex =
        std::max(Base::CoordinateIndex, Index + 1);

    /** Report we are using linear coordinates on this axis. */
    template <typename T> void report_construct(T& fp) const
    {
      fp << "Using Coordinates::Linear on axis " << Index
          << " with spacing " << Spacing << std::endl;
      Base::report_construct(fp);
    }

    /** Return the distance between two adjacent E-field locations. */
    double hstep(typename Base::DataIterator it,
        int Axis=Index) const
    {
      return Axis == Index ? Spacing : Base::hstep(it, Axis);
    }

    /** Return the distance between two adjacent H-field locations. */
    double estep(typename Base::DataIterator it,
        int Axis=Index) const
    {
      return Axis == Index ? Spacing : Base::estep(it, Axis);
    }

    /** Compute the offset of the component vector from the origin. */
    double elength(typename Base::DataIterator it,
        int Axis=Index, int Component=-1) const
    {
      if (Axis != Index)
        return Base::elength(it, Axis, Component);

      unsigned i = Base::getOffset(it, Axis);
      return Component == Axis ? Spacing * i : Spacing * (i+0.5);
    }

    /** Compute the offset of the component vector from the origin. */
    double hlength(typename Base::DataIterator it,
        int Axis=Index, int Component=-1) const
    {
      if (Axis != Index)
        return Base::hlength(it, Axis, Component);

      unsigned i = Base::getOffset(it, Axis);
      return Component == Axis ? Spacing * (i+0.5) : Spacing * i;
    }
  };

};

} // namespace Coordinate
} // namespace Fdtd

#endif // #ifndef FDTD_COORDINATE_LINEAR_HPP

