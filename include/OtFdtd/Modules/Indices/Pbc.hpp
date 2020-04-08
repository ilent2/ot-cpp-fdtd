/* Indices/Pbc.hpp - Periodic boundary condition definitions.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_INDICES_PBC_HPP
#define FDTD_INDICES_PBC_HPP

namespace Fdtd {
namespace Indices {

// TODO: Copy Pbc implementation from version2

/** Check if base has a periodic boundary on axis. */
template <typename Base, int Axis>
struct has_pbc
{
  // TODO
  static constexpr bool value = false;
};

} // namespace Indices
} // namespace Fdtd

#endif // #ifndef FDTD_INDICES_PBC_HPP

