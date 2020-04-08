/* Indices/Polar.hpp - Polar coordinate definitions.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_INDICES_POLAR_HPP
#define FDTD_INDICES_POLAR_HPP

#include "Pbc.hpp"

namespace Fdtd {
namespace Indices {

// TODO: Copy polar implementation from version 2

/** Determine if base has a spherical fold on the given axis. */
template <typename Base, int Axis>
struct has_spherical_fold
{
  // TODO
  static constexpr bool value = false;
};

/** Is base a polar coordinate system. */
template <typename Base>
struct is_polar
{
  static constexpr bool value = Base::AxisIndex == 2
      && has_pbc<Base, 0>::value;
};

/** Is base a cylindrical coordinate system. */
template <typename Base>
struct is_cylindrical
{
  static constexpr bool value = Base::AxisIndex == 3
      && has_pbc<Base, 0>::value;
};

/** Is base a spherical coordinate system. */
template <typename Base>
struct is_spherical
{
  static constexpr bool value = Base::AxisIndex == 3
      && has_pbc<Base, 0>::value && has_spherical_fold<Base, 1>::value;
};

} // namespace Indices
} // namespace Fdtd

#endif // #ifndef FDTD_INDICES_POLAR_HPP

