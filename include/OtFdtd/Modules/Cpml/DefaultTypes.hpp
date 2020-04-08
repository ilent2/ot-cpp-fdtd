/* DefaultTypes.hpp - Default stretching and grading types for CPML.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_CPML_DEFAULT_TYPES_HPP
#define FDTD_CPML_DEFAULT_TYPES_HPP

namespace Fdtd {
namespace Cpml {

/** Provides default graded parameters for stretched CPML. */
class DefaultStretching
{
  constexpr static const double m = 4.0;
  constexpr static const double ma = 1.0;
  constexpr static const double smax = -(m+1.0)*log(1.0e-6)/2.0;
  constexpr static const double amax = 0.05;
  constexpr static const double kmax = 5.0;

public:
  /** Grading function 1. (d in range 0, 1) */
  static double sigma(double d)
  {
    return smax*pow(d, m);
  }

  /** Grading function 2. (d in range 0, 1) */
  static double alpha(double d)
  {
    return amax*pow(1.0 - d, ma);
  }

  /** Stretched coordinate grading function. (d in range 0, 1) */
  static double kappa(double d)
  {
    return 1.0 + pow(d, m)*(kmax - 1.0);
  }
};

/** Provides default graded parameters for non-stretched CPML. */
class DefaultGrading
{
  constexpr static const double m = 4.0;
  constexpr static const double ma = 1.0;
  constexpr static const double smax = -(m+1.0)*log(1.0e-6)/2.0;
  constexpr static const double amax = 0.05;
  constexpr static const double kmax = 5.0;

public:
  /** Grading function 1. (d in range 0, 1) */
  static double sigma(double d)
  {
    return smax*pow(d, m);
  }

  /** Grading function 1.  (d in range 0, 1) */
  static double alpha(double d)
  {
    return amax*pow(1.0 - d, ma);
  }

  static double kappa(double d)
  {
    return 1.0;
  }
};

} // namespace Cpml
} // namespace Fdtd

#endif // #ifndef FDTD_CPML_DEFAULT_TYPES

