/* Fdtd/Coordinate/Cartesian.hpp - Cartesian coordinates.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_COORDINATE_CARTESIAN_HPP
#define FDTD_COORDINATE_CARTESIAN_HPP

#include "../EmptyModule.hpp"
#include "Linear.hpp"
#include "Utilities/DecoratedChild.hpp"
#include "Utilities/Vec.hpp"
#include "../ModuleSelectors.hpp"

namespace Fdtd {
namespace Coordinate {

/** Cartesian coordinates for 1, 2 and 3 dimensions. */
template <const double&... Args> struct Cartesian;

/** Methods for building Cartesian specialisations. */
namespace _Cartesian {

template <const double& first, const double&... others>
class Linear : public EmptyModule
{
  using OurLinear = Coordinate::Linear<sizeof...(others), first>;
  using BaseLinear = _Cartesian::Linear<others...>;

public:

  template <typename Base> using Fallback =
      typename Utilities::DecoratedChild<FallbackSelector,
      BaseLinear, OurLinear>::template Type<Base>;

  template <typename Base> using CoordinateMethods =
      typename Utilities::DecoratedChild<CoordinateMethodsSelector,
      BaseLinear, OurLinear>::template Type<Base>;
};

template <const double& last>
class Linear<last> : public EmptyModule
{
  using OurLinear = Coordinate::Linear<0, last>;

public:

  template <typename Base> using Fallback =
      typename OurLinear::template Fallback<Base>;

  template <typename Base> using CoordinateMethods =
      typename OurLinear::template CoordinateMethods<Base>;
};

template <typename Base>
struct Common : public Base
{
  /** Return the coordinate of a E-field axis. */
  Utilities::Vec3d ecoord(typename Base::DataIterator it,
      int Component=-1) const
  {
    return Utilities::Vec3d(Base::elength(it, 2, Component),
                            Base::elength(it, 1, Component),
                            Base::elength(it, 0, Component));
  }

  /** Return the coordinate of a H-field axis. */
  Utilities::Vec3d hcoord(typename Base::DataIterator it,
      int Component=-1) const
  {
    return Utilities::Vec3d(Base::hlength(it, 2, Component),
                            Base::hlength(it, 1, Component),
                            Base::hlength(it, 0, Component));
  }

  /** Return location at the centre of the given element. */
  Utilities::Vec3d location(typename Base::DataIterator it) const
  {
    // For now we set the centre location of an element to
    //    the location of the E-field coordinate origin, but we
    //    might experiment with using the H-field or average.
    return ecoord(it, -1);
  }

  /** Calculate the area of loc perpendicular to axis0 and axis1. */
  double earea(typename Base::DataIterator loc, int axis0, int axis1) const
  {
    return Base::estep(loc, axis0) * Base::estep(loc, axis1);
  }

  /** Calculate the area of loc towards axis. */
  double earea(typename Base::DataIterator loc, int axis) const
  {
    int axis0 = axis < 1 ? 1 : 0;
    int axis1 = axis < 2 ? 2 : 1;
    return earea(loc, axis0, axis1);
  }

  /** Calculate the area of loc perpendicular to axis0 and axis1. */
  double harea(typename Base::DataIterator loc, int axis0, int axis1) const
  {
    return Base::hstep(loc, axis0) * Base::hstep(loc, axis1);
  }

  /** Calculate the area of loc towards axis. */
  double harea(typename Base::DataIterator loc, int axis) const
  {
    int axis0 = axis < 1 ? 1 : 0;
    int axis1 = axis < 2 ? 2 : 1;
    return harea(loc, axis0, axis1);
  }

};

} // namespace _Cartesian

template <const double& z>
struct Cartesian<z> : public _Cartesian::Linear<z>
{
  template <typename B> class Coordinates : public _Cartesian::Common<B>
  {
    using Base = _Cartesian::Common<B>;

  public:

    /** Return location at the centre of the centre element. */
    Utilities::Vec3d centre(void) const
    {
      auto it = Base::at(Base::length(0)/2);
      return Base::location(it);
    }
  };
};

template <const double& y, const double& z>
struct Cartesian<y, z> : public _Cartesian::Linear<y, z>
{
  template <typename B> class Coordinates : public _Cartesian::Common<B>
  {
    using Base = _Cartesian::Common<B>;

  public:

    /** Return location at the centre of the centre element. */
    Utilities::Vec3d centre(void) const
    {
      auto it = Base::at(Base::length(1)/2, Base::length(0)/2);
      return Base::location(it);
    }
  };
};

template <const double& x, const double& y, const double& z>
struct Cartesian<x, y, z> : public _Cartesian::Linear<x, y, z>
{
  template <typename B> class Coordinates : public _Cartesian::Common<B>
  {
    using Base = _Cartesian::Common<B>;

  public:

    /** Return location at the centre of the centre element. */
    Utilities::Vec3d centre(void) const
    {
      auto it = Base::at(Base::length(2)/2,
          Base::length(1)/2, Base::length(0)/2);
      return Base::location(it);
    }
  };
};

} // namespace Coordinate
} // namespace Fdtd

#endif // #ifndef FDTD_COORDINATE_CARTESIAN_HPP

