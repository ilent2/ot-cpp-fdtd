/* Cpml/CoordinateFromStretching.hpp - Cpml stretched coordinate helper.
 *
 * Provides a Grading object suitable for using in Coordinate/Stretching.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_CPML_COORDINATE_FROM_STRETCHING_HPP
#define FDTD_CPML_COORDINATE_FROM_STRETCHING_HPP

namespace Fdtd {
namespace Cpml {

/** Convert from a stretching to a object suitable for Stretched. */
template <typename Stretching, unsigned Depth, bool Forward>
class CoordinateFromStretching
{
public:
  /** Calculate the e-step stretching for a particular axis index. */
  template <typename Base>
  static double estep(const Base& base, unsigned i, int Axis)
  {
    if (Forward)
    {
      if (i >= Depth) return 1.0;
      return Stretching::kappa((double) (Depth - i)/Depth);
    }
    else
    {
      unsigned offset = base.length(Axis) - Depth;
      if (i < offset) return 1.0;
      return Stretching::kappa((double) (i - offset + 0.5)/Depth);
    }
  }

  /** Calculate the h-step stretching for a particular axis index. */
  template <typename Base>
  static double hstep(const Base& base, unsigned i, unsigned Axis)
  {
    if (Forward)
    {
      if (i >= Depth) return 1.0;
      return Stretching::kappa((double) (Depth - i - 0.5)/Depth);
    }
    else
    {
      unsigned offset = base.length(Axis) - Depth;
      if (i < offset) return 1.0;
      return Stretching::kappa((double) (i - offset + 1.0)/Depth);
    }
  }

  template <typename Base>
  static double offset(const Base& base, unsigned Axis, double old)
  {
    if (Forward)
    {
      double dh = 0.0, de = 0.0;
      for (unsigned i = 0; i < Depth; ++i)
      {
        if (i > 0) de += estep(base, i, Axis);
        dh += hstep(base, i, Axis);
      }

      return (dh - 0.5 - de) * (2.0 * old);
    }
    else
    {
      return old;
    }
  }
};

} // namespace Cpml
} // namespace Fdtd

#endif // #ifndef FDTD_CPML_COORDINATE_FROM_STRETCHING_HPP

