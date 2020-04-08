/* Utilities/OffsetBox.hpp - Provides an object for describing 6 offsets.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef UTILITIES_OFFSET_BOX_HPP
#define UTILITIES_OFFSET_BOX_HPP

#include <string>
#include <cstdio>

namespace Utilities {

/** A simple grouping of a few useful input parameters. */
template <unsigned X0, unsigned X1=X0,
    unsigned Y0=X0, unsigned Y1=X0, unsigned Z0=X0, unsigned Z1=X0>
struct OffsetBox
{
  static constexpr unsigned x0 = X0;
  static constexpr unsigned x1 = X1;
  static constexpr unsigned y0 = Y0;
  static constexpr unsigned y1 = Y1;
  static constexpr unsigned z0 = Z0;
  static constexpr unsigned z1 = Z1;

  static std::string repr(void)
  {
    char szBuf[128];
    snprintf(szBuf, 128, "(%d, %d, %d, %d, %d, %d)", x0, x1, y0, y1, z0, z1);
    return std::string(szBuf);
  }
};

} // namespace Utilities

#endif // #ifndef UTILITIES_OFFSET_BOX_HPP

