/* Output/Integrate/Offset.hpp - Classes for describing the torque origin.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_OUTPUT_INTEGRATE_OFFSET_HPP
#define FDTD_OUTPUT_INTEGRATE_OFFSET_HPP

#include "Utilities/Vec.hpp"

namespace Fdtd {
namespace Output {
namespace Integrate {
namespace Offset {

/** Use the centre of the simulation grid. */
struct Centre
{
  template <typename T>
  static Utilities::Vec3d offset(const T& base)
  {
    return base.centre();
  }
};

} // namespace Offset
} // namespace Integrate
} // namespace Output
} // namespace Fdtd

#endif // #ifndef FDTD_OUTPUT_INTEGRATE_OFFSET_HPP

